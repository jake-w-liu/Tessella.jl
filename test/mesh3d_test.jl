# ── Stage-3 CRC suite: 3-D Delaunay kernel ──────────────────────────────────────
#
# Correctness  : exact empty-circumsphere (is_delaunay3 via insphere_sos) is the
#                defining oracle; convex-region volume conservation (Σ tet vol =
#                bounding-box volume for a box-filling point set); Euler χ=1 (ball);
#                boundary faces form a closed manifold (2-sphere), χ=2.
# Robustness   : random clouds (multiple seeds), the maximally-degenerate cube
#                (8 cospherical corners) resolved by symbolic perturbation, and
#                a structured box grid.
# Completeness : exported tet meshes validate() (positive volumes, manifold).

using Test
using Tessella.Mesh3D
using Tessella.MeshTypes

mesh_vol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
                  for t in 1:ntets(m); init=0.0)

mutable struct _R3; s::UInt64; end
_nf(r::_R3) = (r.s ⊻= r.s<<13; r.s ⊻= r.s>>7; r.s ⊻= r.s<<17; (r.s>>11)/Float64(2^53))

@testset "Mesh3D Delaunay (Stage 3)" begin

    @testset "random clouds: exact empty-circumsphere + valid + manifold" begin
        for seed in (1, 7, 42)
            r = _R3(UInt64(seed)*0x9E3779B97F4A7C15 + 1)
            n = 90
            xs=Float64[_nf(r) for _ in 1:n]; ys=Float64[_nf(r) for _ in 1:n]; zs=Float64[_nf(r) for _ in 1:n]
            T = delaunay3d(xs, ys, zs; rng_seed=seed)
            @test check_consistency3(T)[1]
            dok, nv = is_delaunay3(T)
            @test dok
            @test nv == 0
            m = to_mesh3(T)
            @test validate(m).ok                       # positive volumes, manifold
            @test euler_characteristic(m) == 1         # tetrahedralized ball
            _, maxinc = boundary_faces(m.tets)
            @test maxinc == 2                          # closed manifold boundary
            @test boundary_euler(m) == 2               # boundary is a 2-sphere
        end
    end

    @testset "convex box: Σ tet volume = box volume" begin
        # a grid of points filling [0,2]×[0,3]×[0,1]; convex hull = the box.
        xs=Float64[]; ys=Float64[]; zs=Float64[]
        for i in 0:3, j in 0:4, k in 0:2
            push!(xs, 2*i/3); push!(ys, 3*j/4); push!(zs, 1*k/2)
        end
        T = delaunay3d(xs, ys, zs; rng_seed=3)
        @test check_consistency3(T)[1]
        m = to_mesh3(T)
        @test validate(m).ok
        @test mesh_vol(m) ≈ 6.0 rtol=1e-6              # 2·3·1, perturbation-tolerant
    end

    @testset "degenerate unit cube (8 cospherical corners) resolved" begin
        cx=Float64[0,1,1,0,0,1,1,0]; cy=Float64[0,0,1,1,0,0,1,1]; cz=Float64[0,0,0,0,1,1,1,1]
        T = delaunay3d(cx, cy, cz)
        @test check_consistency3(T)[1]
        @test is_delaunay3(T)[1]
        m = to_mesh3(T)
        @test validate(m).ok                           # NO flat tets
        @test mesh_vol(m) ≈ 1.0 rtol=1e-6
        @test euler_characteristic(m) == 1
    end

    @testset "single tet (4 points)" begin
        T = delaunay3d([0.0,1,0,0],[0.0,0,1,0],[0.0,0,0,1])
        m = to_mesh3(T)
        @test ntets(m) == 1
        @test mesh_vol(m) ≈ 1/6 rtol=1e-6
        @test validate(m).ok
    end

    @testset "seed-independence (deterministic perturbation ⇒ same mesh)" begin
        r = _R3(0xDEAD_BEEF_1234_5678)
        n = 80
        xs=Float64[_nf(r) for _ in 1:n]; ys=Float64[_nf(r) for _ in 1:n]; zs=Float64[_nf(r) for _ in 1:n]
        shas = String[]
        for seed in (1, 5, 99)
            m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=seed))
            push!(shas, mesh_crc(m).sha)
        end
        @test all(==(shas[1]), shas)
    end

    @testset "tetrahedralize: fill a domain from its boundary surface" begin
        # closed cube surface (12 outward triangles) → filled volume 1
        C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
        F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
        ct=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); ct[:,k]=Int32[f...]; end
        cs=Mesh(C; tris=ct)
        @test boundary_edges(cs.tris)[1] |> isempty       # closed surface
        m=tetrahedralize(cs)
        @test mesh_vol(m) ≈ 1.0 rtol=1e-6
        @test validate(m).ok
        _, mi = boundary_faces(m.tets); @test mi == 2      # watertight fill

        # non-convex L-prism (2×2×1 minus a 1×1×1 corner) → volume 3
        base=[(0.0,0.0),(2.0,0.0),(2.0,1.0),(1.0,1.0),(1.0,2.0),(0.0,2.0)]; nb=length(base)
        LC=Matrix{Float64}(undef,3,2nb)
        for (i,(x,y)) in enumerate(base); LC[:,i]=[x,y,0.0]; LC[:,i+nb]=[x,y,1.0]; end
        bt=[(1,2,3),(1,3,4),(1,4,5),(1,5,6)]
        lt=NTuple{3,Int32}[]
        for (a,b,c) in bt; push!(lt,(Int32(a),Int32(c),Int32(b))); end
        for (a,b,c) in bt; push!(lt,(Int32(a+nb),Int32(b+nb),Int32(c+nb))); end
        for i in 1:nb; j=i%nb+1; push!(lt,(Int32(i),Int32(j),Int32(j+nb))); push!(lt,(Int32(i),Int32(j+nb),Int32(i+nb))); end
        ltm=Matrix{Int32}(undef,3,length(lt)); for (k,f) in enumerate(lt); ltm[:,k]=Int32[f...]; end
        ls=Mesh(LC; tris=ltm)
        ml=tetrahedralize(ls)
        @test mesh_vol(ml) ≈ 3.0 rtol=1e-6                 # concave corner excluded
        @test validate(ml).ok

        # convex octahedron (vertices ±1 on axes) → volume 4/3
        OC=Float64[1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1]
        OF=[(1,3,5),(3,2,5),(2,4,5),(4,1,5),(3,1,6),(2,3,6),(4,2,6),(1,4,6)]
        ot=Matrix{Int32}(undef,3,length(OF)); for (k,f) in enumerate(OF); ot[:,k]=Int32[f...]; end
        os=Mesh(OC; tris=ot)
        mo=tetrahedralize(os)
        @test mesh_vol(mo) ≈ 4/3 rtol=1e-6
        @test validate(mo).ok
    end

    @testset "degenerate input handled gracefully" begin
        # exact-degenerate (perturb=false): no non-coplanar 4-tuple ⇒ empty mesh
        @test ntets(to_mesh3(delaunay3d([0.0,1,2],[0.0,1,2],[0.0,1,2]; perturb=false))) == 0     # collinear
        @test ntets(to_mesh3(delaunay3d([0.0,1,0,1],[0.0,0,1,1],[0.0,0,0,0]; perturb=false))) == 0 # coplanar
        @test ntets(to_mesh3(delaunay3d([0.0],[0.0],[0.0]))) == 0                                 # 1 point
    end
end
