# ── Stage-5 CRC suite: native constructive primitive surfaces ───────────────────
#
# Correctness  : each primitive is a closed manifold oriented surface (is_meshable)
#                and fills (via mesh_volume/tetrahedralize) to its exact analytic
#                volume — box, cylinder/cone (N-gon solids), geodesic sphere,
#                box-with-tunnel (genus-1).
# Robustness   : arbitrary axis cylinder, off-origin boxes; bad args rejected.
# Completeness : surfaces feed straight into the volume pipeline.

using Test
using Tessella
using Tessella.MeshTypes
using Tessella.Geometry
using Tessella.Heal

gvol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
              for t in 1:ntets(m); init=0.0)

# Signed volume enclosed by a triangle surface via the divergence theorem:
# V = (1/6) Σ_faces (v1 × v2) · v3.  For an OUTWARD-oriented closed surface this is
# +(enclosed volume); it is sign-sensitive to triangle winding, so it is the oracle
# that distinguishes a correctly-oriented cavity (inner normals into the cavity)
# from an un-reversed one.
function signed_surface_volume(s)
    v = 0.0
    for k in 1:ntris(s)
        a = node(s, s.tris[1,k]); b = node(s, s.tris[2,k]); c = node(s, s.tris[3,k])
        cx = a[2]*b[3]-a[3]*b[2]; cy = a[3]*b[1]-a[1]*b[3]; cz = a[1]*b[2]-a[2]*b[1]
        v += (cx*c[1] + cy*c[2] + cz*c[3])
    end
    v/6
end

function sphere_surface_allocations(subdivisions)
    GC.gc()
    return @allocated sphere_surface((0.,0.,0.),1.;subdivisions=subdivisions)
end

function cone_surface_allocations(sectors)
    GC.gc()
    return @allocated cone_surface((0.,0.,0.),(0.,0.,1.),1.,.5,2.;
                                   nθ=sectors,nz=3)
end

@testset "Geometry primitives (Stage 5)" begin

    @testset "finite geometry and resource contracts" begin
        @test_throws ArgumentError box_surface(0,1,0,1,0,NaN)
        @test_throws ArgumentError cylinder_surface((0.,0.,0.),(Inf,0.,0.),1.,1.)
        @test_throws ArgumentError cylinder_surface((0.,0.,0.),(0.,0.,1.),Inf,1.)
        @test_throws ArgumentError cylinder_surface((0.,0.,0.),(0.,0.,1.),1.,1.;
                                                    nθ=typemax(Int), nz=2)
        @test_throws ArgumentError sphere_surface((0.,0.,0.),0.)
        @test_throws ArgumentError sphere_surface((0.,0.,Inf),1.)
        @test_throws ArgumentError sphere_surface((0.,0.,0.),1.;subdivisions=-1)
        @test_throws ArgumentError sphere_surface((0.,0.,0.),1.;subdivisions=true)
        @test_throws ArgumentError sphere_surface((0.,0.,0.),1.;subdivisions=2,
                                                  max_nodes=65)
        @test_throws ArgumentError sphere_surface((floatmax(Float64),0.,0.),
                                                  floatmax(Float64))
        @test_throws ArgumentError cone_surface((0.,0.,0.),(0.,0.,1.),0.,0.,1.)
        @test_throws ArgumentError cone_surface((0.,0.,0.),(0.,0.,1.),-1.,0.,1.)
        @test_throws ArgumentError cone_surface((0.,0.,0.),(0.,0.,0.),1.,0.,1.)
        @test_throws ArgumentError cone_surface((0.,0.,0.),(0.,0.,1.),1.,0.,0.)
        @test_throws ArgumentError cone_surface((0.,0.,0.),(0.,0.,1.),1.,0.,1.;nθ=2)
        @test_throws ArgumentError cone_surface((0.,0.,0.),(0.,0.,1.),1.,0.,1.;nz=1)
        @test_throws ArgumentError cone_surface((0.,0.,0.),(0.,0.,1.),1.,0.,1.;
                                                nθ=16,nz=3,max_triangles=63)
        @test_throws ArgumentError box_tunnel_surface(0,4,0,4,0,2,1,3,1,NaN)
        @test_throws ArgumentError box_shell_surface(0,4,0,4,0,4,1,3,1,3,1,Inf)
    end

    @testset "sphere_surface: projected geodesic surface and convergent volume" begin
        coarse=sphere_surface((0.,0.,0.),1.;subdivisions=0)
        medium=sphere_surface((0.,0.,0.),1.;subdivisions=1)
        surface=sphere_surface((0.,0.,0.),1.;subdivisions=2)
        @test (nnodes(coarse),ntris(coarse))==(6,8)
        @test (nnodes(medium),ntris(medium))==(18,32)
        @test (nnodes(surface),ntris(surface))==(66,128)
        @test mesh_crc(surface).sha ==
              "2c7bf12222ab5796df858b3ef015349be3fc8acecff5d445963f41215377bd54"
        @test all(isapprox(hypot(node(surface,i)...),1.;atol=8eps(Float64),rtol=0)
                  for i in 1:nnodes(surface))
        volumes=signed_surface_volume.((coarse,medium,surface))
        @test 0 < volumes[1] < volumes[2] < volumes[3] < 4pi/3
        @test is_meshable(surface)[1]
        volume=tetrahedralize(surface)
        @test validate(volume).ok
        @test gvol(volume) ≈ volumes[3] rtol=2e-13

        translated=sphere_surface((1.,-2.,.5),2.;subdivisions=1)
        @test all(isapprox(hypot((node(translated,i)[1]-1.,
                                 node(translated,i)[2]+2.,
                                 node(translated,i)[3]-.5)...),2.;
                           atol=32eps(Float64),rtol=8eps(Float64))
                  for i in 1:nnodes(translated))
        @test is_meshable(translated)[1]
    end

    @testset "cone_surface: cones and frusta are closed and volume-correct" begin
        sectors=24;height=2.;r1=1.;r2=.5
        frustum=cone_surface((0.,0.,0.),(0.,0.,1.),r1,r2,height;
                             nθ=sectors,nz=3)
        @test (nnodes(frustum),ntris(frustum))==(74,144)
        @test mesh_crc(frustum).sha ==
              "8afc4d9f3bf9740313d9ce099302acb32335a9de392d8940d60cdb42b9465115"
        @test is_meshable(frustum)[1]
        polygon_factor=0.5*sectors*sinpi(2/sectors)
        expected=polygon_factor*height*(r1^2+r1*r2+r2^2)/3
        @test signed_surface_volume(frustum) ≈ expected rtol=2e-14
        filled=tetrahedralize(frustum)
        @test validate(filled).ok
        @test gvol(filled) ≈ expected rtol=1e-12

        lower_apex=cone_surface((1.,-1.,.5),(1.,2.,3.),0.,1.,2.;
                                nθ=16,nz=3)
        upper_apex=cone_surface((1.,-1.,.5),(1.,2.,3.),1.,0.,2.;
                                nθ=16,nz=3)
        @test (nnodes(lower_apex),ntris(lower_apex))==(34,64)
        @test (nnodes(upper_apex),ntris(upper_apex))==(34,64)
        @test signed_surface_volume(lower_apex) > 0
        @test signed_surface_volume(upper_apex) > 0
        @test validate(tetrahedralize(lower_apex)).ok
        @test validate(tetrahedralize(upper_apex)).ok
    end

    @testset "sphere/cone construction has linear allocation growth" begin
        # Warm each specialized path before measuring. A sphere subdivision grows
        # its output by 4x; doubling cone sectors grows its output by 2x.
        sphere_surface((0.,0.,0.),1.;subdivisions=3)
        cone_surface((0.,0.,0.),(0.,0.,1.),1.,.5,2.;nθ=96,nz=3)
        sphere2=sphere_surface_allocations(2)
        sphere3=sphere_surface_allocations(3)
        cone48=cone_surface_allocations(48)
        cone96=cone_surface_allocations(96)
        @test sphere3 <= 5sphere2
        @test cone96 <= 3cone48
    end

    @testset "box_surface: closed, meshable, exact volume" begin
        s = box_surface(-1, 2, 0, 3, 1, 5)          # 3×3×4 = 36
        @test is_meshable(s)[1]
        m = tetrahedralize(s)
        @test gvol(m) ≈ 36.0 rtol=1e-6
        @test validate(m).ok
        @test_throws ArgumentError box_surface(1, 0, 0, 1, 0, 1)
    end

    @testset "cylinder_surface: watertight, N-gon prism volume" begin
        s = cylinder_surface((0.0,0.0,0.0), (0.0,0.0,1.0), 2.0, 5.0; nθ=24, nz=3)
        @test is_meshable(s)[1]
        m = tetrahedralize(s)
        @test validate(m).ok
        # exact 24-gon prism volume = ½·N·R²·sin(2π/N)·H
        @test gvol(m) ≈ 0.5*24*2.0^2*sinpi(2/24)*5.0 rtol=1e-6
        # tilted axis still watertight + meshable
        s2 = cylinder_surface((1.0,1.0,1.0), (1.0,1.0,1.0), 1.0, 3.0; nθ=16, nz=2)
        @test is_meshable(s2)[1]
        @test validate(tetrahedralize(s2)).ok
    end

    @testset "box_tunnel_surface: genus-1 bore, exact volume" begin
        s = box_tunnel_surface(1,5, 1,5, 1,3, 2,4, 2,4)   # 32 − 8 = 24
        @test is_meshable(s)[1]
        m = tetrahedralize(s)
        @test gvol(m) ≈ 24.0 rtol=1e-6
        @test validate(m).ok
        @test boundary_faces(m.tets)[2] == 2               # watertight fill
        @test_throws ArgumentError box_tunnel_surface(1,5,1,5,1,3, 0,6,2,4)  # inner not inside
    end

    @testset "box_shell_surface: hollow box (Boolean difference), oriented, exact volume" begin
        # outer [-1,2]×[0,3]×[1,5] = 3·3·4 = 36 ; inner [0,1]×[1,2]×[2,3] = 1 ⇒ 35
        s = box_shell_surface(-1,2, 0,3, 1,5, 0,1, 1,2, 2,3)
        @test is_meshable(s)[1]                             # closed & manifold (2 components, 0 open edges)
        # ORIENTATION ORACLE (mutation-sensitive): outer outward (+36) + inner into
        # cavity (−1) ⇒ +35. An un-reversed inner winding would read +37, not 35.
        @test signed_surface_volume(s) ≈ 35.0 rtol=1e-12
        m = tetrahedralize(s)
        @test gvol(m) ≈ 35.0 rtol=1e-6                     # fill volume = outer − inner (cavity carved)
        @test validate(m).ok
        @test boundary_faces(m.tets)[2] == 2               # manifold fill (max face incidence 2)
        # non-strictly-inside inner boxes are rejected
        @test_throws ArgumentError box_shell_surface(-1,2,0,3,1,5, -1,1, 1,2, 2,3)  # inner touches outer face
        @test_throws ArgumentError box_shell_surface(-1,2,0,3,1,5, 0,1, -1,4, 2,3)  # inner spans past outer
    end

    @testset "primitives feed mesh_volume directly" begin
        @test validate(mesh_volume(box_surface(0,1,0,1,0,1))).ok
        @test validate(mesh_volume(cylinder_surface((0.,0,0),(0.,0,1),1.,2.; nθ=20))).ok
        @test validate(mesh_volume(sphere_surface((0.,0.,0.),1.;subdivisions=1))).ok
        @test validate(mesh_volume(cone_surface((0.,0.,0.),(0.,0.,1.),1.,0.,2.;
                                                nθ=16,nz=2))).ok
    end
end
