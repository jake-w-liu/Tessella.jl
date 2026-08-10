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

    @testset "degenerate input handled gracefully" begin
        # exact-degenerate (perturb=false): no non-coplanar 4-tuple ⇒ empty mesh
        @test ntets(to_mesh3(delaunay3d([0.0,1,2],[0.0,1,2],[0.0,1,2]; perturb=false))) == 0     # collinear
        @test ntets(to_mesh3(delaunay3d([0.0,1,0,1],[0.0,0,1,1],[0.0,0,0,0]; perturb=false))) == 0 # coplanar
        @test ntets(to_mesh3(delaunay3d([0.0],[0.0],[0.0]))) == 0                                 # 1 point
    end
end
