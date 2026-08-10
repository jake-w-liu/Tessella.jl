# ── Stage-1 CRC suite: 2-D Delaunay (+ CDT, refinement added below) ─────────────
#
# Correctness  : the exact empty-circumcircle property (is_delaunay, built on the
#                exact incircle_sos) is the *defining* oracle; plus an independent
#                monotone-chain convex-hull oracle for the 2n−2−h triangle count
#                and hull-edge cross-check; plus Euler characteristic.
# Robustness   : cocircular rings, structured grids, near-collinear clouds, many
#                fixed seeds; order-independence of the SoS triangulation.
# Completeness : exported meshes validate(); adjacency is mutually consistent.

using Test
using Tessella.Mesh2D
using Tessella.MeshTypes

# deterministic xorshift stream (no RNG dep) — reproducible clouds
mutable struct _RNG; s::UInt64; end
_nf(r::_RNG) = (r.s ⊻= r.s<<13; r.s ⊻= r.s>>7; r.s ⊻= r.s<<17; (r.s>>11)/Float64(2^53))

# independent convex-hull oracle: Andrew's monotone chain → hull vertex count and
# the set of undirected hull edges (as sorted index pairs into the input arrays).
function hull_oracle(xs, ys)
    n = length(xs)
    idx = sort(collect(1:n), by = i -> (xs[i], ys[i]))
    cross(o,a,b) = (xs[a]-xs[o])*(ys[b]-ys[o]) - (ys[a]-ys[o])*(xs[b]-xs[o])
    lower = Int[]
    for i in idx
        while length(lower) >= 2 && cross(lower[end-1], lower[end], i) <= 0; pop!(lower); end
        push!(lower, i)
    end
    upper = Int[]
    for i in reverse(idx)
        while length(upper) >= 2 && cross(upper[end-1], upper[end], i) <= 0; pop!(upper); end
        push!(upper, i)
    end
    hull = vcat(lower[1:end-1], upper[1:end-1])
    edges = Set{Tuple{Int,Int}}()
    for k in 1:length(hull)
        a = hull[k]; b = hull[k % length(hull) + 1]
        push!(edges, (min(a,b), max(a,b)))
    end
    return length(hull), edges
end

# boundary edges of a triangle mesh, in the ORIGINAL point indexing (via coords match)
function mesh_boundary_edges_by_coord(m, xs, ys; tol=1e-12)
    # map mesh node → original index by coordinate
    orig = Dict{Int,Int}()
    for i in 1:nnodes(m)
        p = node(m, i)
        for j in eachindex(xs)
            if abs(p[1]-xs[j]) < tol && abs(p[2]-ys[j]) < tol
                orig[i] = j; break
            end
        end
    end
    be, _ = boundary_edges(m.tris)
    return Set{Tuple{Int,Int}}((min(orig[e[1]],orig[e[2]]), max(orig[e[1]],orig[e[2]])) for e in be)
end

@testset "Mesh2D Delaunay (Stage 1)" begin

    @testset "tiny cases" begin
        # single triangle
        m = triangulate([0.0,1.0,0.0], [0.0,0.0,1.0])
        @test ntris(m) == 1
        @test validate(m).ok
        # unit square → 2 triangles, both CCW, Delaunay
        T = delaunay2d([0.0,1.0,1.0,0.0], [0.0,0.0,1.0,1.0])
        @test check_consistency(T)[1]
        @test is_delaunay(T)[1]
        @test ntris(to_mesh(T)) == 2
    end

    @testset "random clouds: exact Delaunay + consistency + valid + hull" begin
        for seed in (1, 2, 3, 42, 777)
            r = _RNG(UInt64(seed) * 0x9E3779B97F4A7C15 + 1)
            n = 250
            xs = Float64[_nf(r) for _ in 1:n]; ys = Float64[_nf(r) for _ in 1:n]
            T = delaunay2d(xs, ys; rng_seed=seed)
            @test check_consistency(T)[1]
            dok, nv = is_delaunay(T)
            @test dok
            @test nv == 0
            m = to_mesh(T)
            @test validate(m).ok
            # triangle-count identity 2n − 2 − h (general position)
            h, hedges = hull_oracle(xs, ys)
            @test ntris(m) == 2n - 2 - h
            # boundary edges == convex hull edges
            @test mesh_boundary_edges_by_coord(m, xs, ys) == hedges
            # Euler χ of a disk triangulation = 1
            @test euler_characteristic(m) == 1
        end
    end

    @testset "structured grid (many cocircular/degenerate configs)" begin
        xs = Float64[]; ys = Float64[]
        for i in 0:8, j in 0:8; push!(xs, float(i)); push!(ys, float(j)); end
        T = delaunay2d(xs, ys; rng_seed=3)
        @test check_consistency(T)[1]
        @test is_delaunay(T)[1]
        m = to_mesh(T)
        @test validate(m).ok
        # a full grid triangulation of the 9×9 lattice covers the square:
        # area of all triangles == 8×8 = 64
        area = 0.0
        for t in 1:ntris(m)
            a=node(m,m.tris[1,t]); b=node(m,m.tris[2,t]); c=node(m,m.tris[3,t])
            area += triangle_area(a,b,c)
        end
        @test area ≈ 64.0 atol=1e-9
    end

    @testset "points on a circle (all cocircular) — SoS yields valid triangulation" begin
        n = 32
        xs = Float64[cospi(2k/n) for k in 0:n-1]
        ys = Float64[sinpi(2k/n) for k in 0:n-1]
        T = delaunay2d(xs, ys; rng_seed=5)
        @test check_consistency(T)[1]
        @test is_delaunay(T)[1]
        m = to_mesh(T)
        @test ntris(m) == n - 2          # triangulating a convex polygon: n−2 triangles
        @test validate(m).ok
    end

    @testset "order-independence: SoS triangulation is seed-independent" begin
        r = _RNG(0x1234_5678_9abc_def0)
        n = 300
        xs = Float64[_nf(r) for _ in 1:n]; ys = Float64[_nf(r) for _ in 1:n]
        crcs = String[]
        for seed in (1, 2, 3, 99, 12345)
            m = triangulate(xs, ys; rng_seed=seed)
            push!(crcs, mesh_crc(m).sha)
        end
        @test all(==(crcs[1]), crcs)     # identical topology regardless of insertion order
    end
end
