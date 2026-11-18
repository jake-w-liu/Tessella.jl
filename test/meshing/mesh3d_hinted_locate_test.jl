# delaunay3d starts each point-location walk from a nearby inserted vertex (bucket
# grid) and matches cavity spokes in an open-addressing table instead of a Dict.
# The tetrahedralization, including tet storage order, must be bit-identical to
# the unhinted insertion sequence; points on a face plane of the found tet fall
# back to the default walk (their ghost/SoS location is path dependent).
using Test
using Tessella
const M3 = Tessella.Mesh3D

function _unhinted_delaunay3(xs, ys, zs, seed)
    ux, uy, uz, _ = M3._dedup3(xs, ys, zs)
    T = M3.Triangulation3(ux, uy, uz)
    placed = M3._init3!(T)
    order = Int32[i for i in 1:T.nreal if !(Int32(i) in placed)]
    M3._shuffle3!(order, UInt64(seed))
    for v in order; M3.insert_point3!(T, v); end
    return T
end

@testset "hinted 3-D point location leaves the tetrahedralization unchanged" begin
    s = UInt64(0x2545F4914F6CDD1D)
    rnd() = (s ⊻= s << 13; s ⊻= s >> 7; s ⊻= s << 17; Float64(s >> 11) / 2.0^53)
    clouds = Any[]
    push!(clouds, ([rnd() for _ in 1:1500], [rnd() for _ in 1:1500], [rnd() for _ in 1:1500]))
    g = [(i, j, k) for i in 0:6 for j in 0:6 for k in 0:6]            # cospherical lattice
    push!(clouds, (Float64[p[1] for p in g], Float64[p[2] for p in g], Float64[p[3] for p in g]))
    # coplanar layers (hull faces with many coplanar points)
    push!(clouds, ([rnd() for _ in 1:600], [rnd() for _ in 1:600], Float64[i % 4 for i in 1:600]))
    for (xs, ys, zs) in clouds, seed in (1, 5)
        T0 = _unhinted_delaunay3(xs, ys, zs, seed)
        T1 = M3.delaunay3d(copy(xs), copy(ys), copy(zs); rng_seed=seed, perturb=false)
        @test T1.tv == T0.tv && T1.tn == T0.tn && T1.alive == T0.alive
        @test M3.check_consistency3(T1)[1]
    end
end
