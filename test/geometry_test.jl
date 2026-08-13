# ── Stage-5 CRC suite: native constructive primitive surfaces ───────────────────
#
# Correctness  : each primitive is a closed manifold oriented surface (is_meshable)
#                and fills (via mesh_volume/tetrahedralize) to its exact analytic
#                volume — box, cylinder (N-gon prism), box-with-tunnel (genus-1).
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

@testset "Geometry primitives (Stage 5)" begin

    @testset "finite geometry and resource contracts" begin
        @test_throws ArgumentError box_surface(0,1,0,1,0,NaN)
        @test_throws ArgumentError cylinder_surface((0.,0.,0.),(Inf,0.,0.),1.,1.)
        @test_throws ArgumentError cylinder_surface((0.,0.,0.),(0.,0.,1.),Inf,1.)
        @test_throws ArgumentError cylinder_surface((0.,0.,0.),(0.,0.,1.),1.,1.;
                                                    nθ=typemax(Int), nz=2)
        @test_throws ArgumentError box_tunnel_surface(0,4,0,4,0,2,1,3,1,NaN)
        @test_throws ArgumentError box_shell_surface(0,4,0,4,0,4,1,3,1,3,1,Inf)
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
    end
end
