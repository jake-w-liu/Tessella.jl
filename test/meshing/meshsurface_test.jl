# ── Stage-2 CRC suite: surface meshing (planar, cylinder, parametric) ───────────
#
# Correctness  : planar face area is exact (isometric projection); cylinder and
#                parametric areas converge to the analytic value as h→0; min-angle
#                quality bound holds on the surface.
# Robustness   : tilted planes, holes, welded cylinder seam (watertight), poles-
#                avoided parametric patch.
# Completeness : lifted meshes validate(); size field honored.

using Test
using Tessella.MeshSurface
using Tessella.MeshTypes
using Tessella.SizeField

surf_area(m) = sum(triangle_area(node(m,m.tris[1,t]),node(m,m.tris[2,t]),node(m,m.tris[3,t]))
                   for t in 1:ntris(m); init=0.0)
function surf_min_angle(m)
    mn = 180.0
    for t in 1:ntris(m)
        p = (node(m,m.tris[1,t]), node(m,m.tris[2,t]), node(m,m.tris[3,t]))
        for i in 1:3
            a=p[i]; b=p[mod1(i+1,3)]; c=p[mod1(i+2,3)]
            v1=(b[1]-a[1],b[2]-a[2],b[3]-a[3]); v2=(c[1]-a[1],c[2]-a[2],c[3]-a[3])
            d=(v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3])/(sqrt(sum(v1.^2))*sqrt(sum(v2.^2)))
            mn=min(mn, acosd(clamp(d,-1.0,1.0)))
        end
    end
    return mn
end

@testset "MeshSurface (Stage 2)" begin

    @testset "surface input contracts" begin
        sf = ConstantSize(1.0)
        @test_throws ArgumentError plane_frame([(0.0,0.0,0.0),(1.0,0.0,0.0),(NaN,1.0,0.0)])
        @test_throws ArgumentError mesh_planar_face(
            [[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,1.0),(0.0,1.0,0.0)]], sf)
        @test_throws ArgumentError mesh_cylinder_face(
            (0.0,0.0,0.0),(0.0,0.0,1.0),1.0,Inf,sf)
        @test_throws ArgumentError mesh_cylinder_face(
            (0.0,0.0,0.0),(0.0,0.0,0.0),1.0,1.0,sf)
        @test_throws ArgumentError mesh_parametric_face(
            (u,v)->(u,v,0.0),1.0,0.0,0.0,1.0,sf)
        @test_throws ArgumentError mesh_parametric_face(
            (u,v)->(NaN,v,0.0),0.0,1.0,0.0,1.0,sf)
        @test_throws ArgumentError mesh_parametric_face(
            (u,v)->(u,0.0,0.0),0.0,1.0,0.0,1.0,sf)
    end

    @testset "tilted planar rectangle: exact area + quality" begin
        o=(1.0,2.0,3.0); u=(1/sqrt(2),1/sqrt(2),0.0); v=(0.0,0.0,1.0)
        c1=o; c2=o.+6.0.*u; c3=o.+6.0.*u.+4.0.*v; c4=o.+4.0.*v
        m=mesh_planar_face([[c1,c2,c3,c4]], ConstantSize(0.5); min_angle_deg=25.0)
        @test surf_area(m) ≈ 24.0 atol=1e-8        # 6×4, isometric ⇒ exact
        @test surf_min_angle(m) >= 25.0 - 1e-6
        @test validate(m).ok
    end

    @testset "planar face with a hole" begin
        outer=[(0.0,0.0,0.0),(10.0,0.0,0.0),(10.0,10.0,0.0),(0.0,10.0,0.0)]
        hole =[(4.0,4.0,0.0),(4.0,6.0,0.0),(6.0,6.0,0.0),(6.0,4.0,0.0)]
        m=mesh_planar_face([outer,hole], ConstantSize(1.0); min_angle_deg=25.0)
        @test surf_area(m) ≈ 96.0 atol=1e-8        # 100 − 4 hole
        @test surf_min_angle(m) >= 25.0 - 1e-6
        @test validate(m).ok
    end

    @testset "cylinder lateral surface: watertight + area convergence" begin
        exact = 2π*2.0*5.0
        mc = mesh_cylinder_face((0.0,0.0,0.0),(0.0,0.0,1.0), 2.0, 5.0, ConstantSize(0.8))
        mf = mesh_cylinder_face((0.0,0.0,0.0),(0.0,0.0,1.0), 2.0, 5.0, ConstantSize(0.3))
        @test surf_area(mc) < exact                # chordal under-estimate
        @test surf_area(mf) < exact
        @test abs(surf_area(mf) - exact) < abs(surf_area(mc) - exact)   # converges
        @test surf_area(mf) ≈ exact rtol=0.02      # fine mesh within 2%
        @test validate(mc).ok
        @test surf_min_angle(mc) >= 25.0 - 2.0     # quality (small seam tolerance)
        # watertight: no non-manifold edges; only the two rims are boundary
        be, maxinc = boundary_edges(mc.tris)
        @test maxinc == 2                          # manifold
        # even under a strongly varying size field, the seam has NO gap edges
        # (structured wrap ⇒ watertight; a Ruppert-unrolled mesh would gap)
        sf2 = FunctionSize((x,y,z) -> 0.2 + 0.15*abs(z) + 0.1*abs(x))
        mv = mesh_cylinder_face((0.0,0.0,0.0),(0.0,0.0,1.0), 2.0, 5.0, sf2)
        bev, _ = boundary_edges(mv.tris)
        seamgap = count(e -> !((abs(node(mv,e[1])[3])<1e-6 && abs(node(mv,e[2])[3])<1e-6) ||
                               (abs(node(mv,e[1])[3]-5)<1e-6 && abs(node(mv,e[2])[3]-5)<1e-6)), bev)
        @test seamgap == 0
        @test validate(mv).ok
        controlled=mesh_cylinder_face((0.,0.,0.),(0.,0.,1.),2.0,5.0,ConstantSize(0.8);
                                      min_angle_deg=25.0,max_area=0.05)
        @test surf_min_angle(controlled)>=25.0-1e-9
        @test maximum(triangle_area(node(controlled,controlled.tris[1,t]),
                                    node(controlled,controlled.tris[2,t]),
                                    node(controlled,controlled.tris[3,t])) for t in 1:ntris(controlled))<=0.05*(1+1e-12)
        @test_throws ArgumentError mesh_cylinder_face((0.,0.,0.),(0.,0.,1.),2.0,5.0,
                                                       ConstantSize(0.8);min_angle_deg=40.0)
    end

    @testset "parametric planar patch equals 2-D meshing" begin
        s(u,v) = (u, v, 0.0)                       # flat parametrization
        m=mesh_parametric_face(s, 0.0, 4.0, 0.0, 3.0, ConstantSize(0.5); min_angle_deg=25.0)
        @test surf_area(m) ≈ 12.0 atol=1e-8        # 4×3, isometric
        @test validate(m).ok
    end

    @testset "parametric curved patch: valid + area convergence" begin
        # paraboloid z = 0.1(u²+v²) over [-2,2]²; analytic area computed by fine quadrature
        s(u,v) = (u, v, 0.1*(u^2+v^2))
        analytic = let A=0.0, N=400, du=4.0/N
            for i in 0:N-1, j in 0:N-1
                u=-2+ (i+0.5)*du; v=-2+(j+0.5)*du
                zu=0.2u; zv=0.2v; A += sqrt(1+zu^2+zv^2)*du*du
            end; A
        end
        mc=mesh_parametric_face(s, -2.0,2.0,-2.0,2.0, ConstantSize(0.5))
        mf=mesh_parametric_face(s, -2.0,2.0,-2.0,2.0, ConstantSize(0.2))
        @test validate(mc).ok && validate(mf).ok
        @test surf_area(mc) < analytic + 1e-6      # chordal ≤ analytic (convex-ish)
        @test abs(surf_area(mf)-analytic) < abs(surf_area(mc)-analytic)  # converges
        @test surf_area(mf) ≈ analytic rtol=0.01
    end
end
