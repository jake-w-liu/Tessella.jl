# The 2-D Delaunay drivers start each point-location walk from a nearby inserted
# vertex (bucket grid) instead of the previous insertion; with random insertion
# order the old walk was O(sqrt(n)) per point (constrained_delaunay of 1e5 random
# points: 0.86 s -> 0.10 s).  The triangulation must be bit-identical to the
# unhinted insertion sequence, including storage order (refinement scans slots).
using Test
using Tessella
const M2 = Tessella.Mesh2D

function _unhinted_delaunay(xs, ys, seed)
    ux, uy, _ = M2.dedup_points(xs, ys)
    T = M2.Triangulation(ux, uy)
    placed = M2._init_triangulation!(T)
    order = Int32[i for i in 1:T.nreal if !(Int32(i) in placed)]
    M2._shuffle_det!(order, UInt64(seed))
    for v in order; M2._insert_point!(T, v); end
    return T
end

@testset "hinted point location leaves the triangulation unchanged" begin
    s = UInt64(0x9E3779B97F4A7C15)
    rnd() = (s ⊻= s << 13; s ⊻= s >> 7; s ⊻= s << 17; Float64(s >> 11) / 2.0^53)
    clouds = Any[]
    # Sizes for the interpreted default suite.
    push!(clouds, ([rnd() for _ in 1:500], [rnd() for _ in 1:500]))
    push!(clouds, (Float64[i for i in 0:15 for j in 0:15], Float64[j for i in 0:15 for j in 0:15]))
    push!(clouds, ([1e9 + rnd() for _ in 1:300], [rnd() * 1e-3 for _ in 1:300]))
    push!(clouds, ([cospi(k / 100) for k in 0:199], [sinpi(k / 100) for k in 0:199]))
    for (xs, ys) in clouds, seed in (1, 7)
        T0 = _unhinted_delaunay(xs, ys, seed)
        T1 = M2.delaunay2d(xs, ys; rng_seed=seed)
        T2 = M2.constrained_delaunay(xs, ys, Tuple{Int,Int}[]; rng_seed=seed)
        @test T1.tv == T0.tv && T1.tn == T0.tn && T1.alive == T0.alive
        @test T2.tv == T0.tv && T2.tn == T0.tn && T2.alive == T0.alive
        @test M2.check_consistency(T1)[1]
        length(xs) <= 900 && @test M2.is_delaunay(T1)[1]
    end
end
