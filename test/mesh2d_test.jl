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

    @testset "public input and state contracts" begin
        @test_throws ArgumentError dedup_points([0.0], [0.0, 1.0])
        @test_throws ArgumentError delaunay2d([0.0, 1.0, NaN], [0.0, 0.0, 1.0])
        @test_throws ArgumentError delaunay2d([0.0, 1.0, 0.0], [0.0, Inf, 1.0])
        @test_throws ArgumentError delaunay2d([0.0, 1.0, 0.0], [0.0, 0.0, 1.0]; rng_seed=-1)
        @test_throws ArgumentError constrained_delaunay(
            [0.0,1.0,0.0], [0.0,0.0,1.0], Tuple{Int,Int}[]; rng_seed=big(2)^80)

        ux, uy, remap = dedup_points([-0.0, 0.0], [0.0, -0.0])
        @test ux == [0.0] && uy == [0.0]
        @test !signbit(ux[1]) && !signbit(uy[1])
        @test remap == Int32[1, 1]

        T = delaunay2d([0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
        @test_throws ArgumentError insert_point!(T, 0)
        @test_throws ArgumentError insert_point!(T, 4)
        @test_throws ArgumentError insert_point!(T, 1)       # already inserted
        @test_throws ArgumentError insert_segment!(T, 0, 1)
        @test_throws ArgumentError insert_segment!(T, 1, 1)
        @test_throws ArgumentError to_mesh(T; interior=Bool[true])

        @test_throws ArgumentError constrained_delaunay(
            [0.0,1.0,0.0], [0.0,0.0,1.0], [(1,4)])
        @test_throws ArgumentError constrained_delaunay(
            [0.0,1.0,0.0], [0.0,0.0,1.0], [(2,2)])
        @test_throws ArgumentError constrained_delaunay(
            [0.0,-0.0,1.0,0.0], [0.0,0.0,0.0,1.0], [(1,2)])

        @test_throws ArgumentError Tessella.mesh_planar(
            Float64[0,1,0], Float64[0,0,1], Tuple{Int,Int}[])
        @test_throws ArgumentError Tessella.mesh_planar(
            Float64[0,1,1,0], Float64[0,0,1,1], [(1,2)])
    end

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

# ── independent oracle: total mesh area over a triangle Mesh ─────────────────────
mesh_area(m) = sum(triangle_area(node(m,m.tris[1,t]),node(m,m.tris[2,t]),node(m,m.tris[3,t]))
                   for t in 1:ntris(m); init=0.0)
edge_in(T, vi, vj) = Tessella.Mesh2D._edge_between(T, Int32(vi), Int32(vj))[1] != 0

@testset "Mesh2D CDT (Stage 1)" begin

    @testset "interior classification requires closed boundary loops" begin
        xs=Float64[0,1,1,0]; ys=Float64[0,0,1,1]
        openT = constrained_delaunay(xs, ys, [(1,2)])
        @test_throws ArgumentError classify_interior(openT)
    end

    @testset "skew constraint crossing a grid" begin
        xs=Float64[]; ys=Float64[]
        for i in 0:4, j in 0:4; push!(xs,float(i)); push!(ys,float(j)); end
        p00=findfirst(k->xs[k]==0&&ys[k]==0,1:25); p41=findfirst(k->xs[k]==4&&ys[k]==1,1:25)
        T=constrained_delaunay(xs,ys,[(p00,p41)];rng_seed=1)
        @test check_consistency(T)[1]
        ok,miss,viol = is_constrained_delaunay(T)
        @test ok && miss==0 && viol==0
        @test edge_in(T, p00, p41)                    # the constraint is an edge
    end

    @testset "through-vertex constraint splits into collinear sub-edges" begin
        xs=Float64[]; ys=Float64[]
        for i in 0:4, j in 0:4; push!(xs,float(i)); push!(ys,float(j)); end
        p00=findfirst(k->xs[k]==0&&ys[k]==0,1:25); p44=findfirst(k->xs[k]==4&&ys[k]==4,1:25)
        T=constrained_delaunay(xs,ys,[(p00,p44)];rng_seed=1)
        @test check_consistency(T)[1]
        @test is_constrained_delaunay(T)[1]
        # the collinear sub-edges (0,0)-(1,1)-(2,2)-(3,3)-(4,4) are all present
        pt(a,b)=findfirst(k->xs[k]==a&&ys[k]==b,1:25)
        diag = [(0,0),(1,1),(2,2),(3,3),(4,4)]
        for i in 1:length(diag)-1
            @test edge_in(T, pt(diag[i]...), pt(diag[i+1]...))
        end
    end

    @testset "already-an-edge constraint is a no-op (just marked)" begin
        xs=Float64[0,1,0,1]; ys=Float64[0,0,1,1]
        T0=delaunay2d(xs,ys)
        # find any existing edge and constrain it
        # (0,0)-(1,0) is a hull edge → an edge
        T=constrained_delaunay(xs,ys,[(1,2)])
        @test check_consistency(T)[1]
        @test is_constrained_delaunay(T)[1]
        @test edge_in(T,1,2)
    end

    @testset "square domain: interior area == domain area" begin
        r = _RNG(0xBEEF)
        xs=Float64[0,10,10,0]; ys=Float64[0,0,10,10]
        for _ in 1:150; push!(xs, _nf(r)*10); push!(ys, _nf(r)*10); end
        segs=[(1,2),(2,3),(3,4),(4,1)]
        T=constrained_delaunay(xs,ys,segs;rng_seed=3)
        @test check_consistency(T)[1]
        @test is_constrained_delaunay(T)[1]
        m=to_mesh(T; interior=classify_interior(T))
        @test mesh_area(m) ≈ 100.0 atol=1e-9
        @test validate(m).ok
    end

    @testset "domain with a hole: interior excludes the hole" begin
        r = _RNG(0xC0FFEE)
        xs=Float64[0,10,10,0, 4,6,6,4]; ys=Float64[0,0,10,10, 4,4,6,6]
        for _ in 1:200
            x=_nf(r)*10; y=_nf(r)*10
            (4<x<6 && 4<y<6) && continue
            push!(xs,x); push!(ys,y)
        end
        segs=[(1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5)]
        T=constrained_delaunay(xs,ys,segs;rng_seed=2)
        @test check_consistency(T)[1]
        @test is_constrained_delaunay(T)[1]
        m=to_mesh(T; interior=classify_interior(T))
        @test mesh_area(m) ≈ 96.0 atol=1e-9    # 100 − 4 (hole)
        @test validate(m).ok
    end

    @testset "many random non-crossing constraints stay Delaunay-constrained" begin
        # constrain a fan of chords from one hull-ish point that don't cross
        r = _RNG(0x5EED)
        n=120
        xs=Float64[_nf(r) for _ in 1:n]; ys=Float64[_nf(r) for _ in 1:n]
        # constrain edges of the point with smallest x to a few others (star, no crossings)
        c = argmin(xs)
        others = sort(setdiff(1:n, [c]); by=k->atan(ys[k]-ys[c], xs[k]-xs[c]))
        segs = [(c, others[k]) for k in 1:10:length(others)]
        T=constrained_delaunay(xs,ys,segs;rng_seed=7)
        @test check_consistency(T)[1]
        ok,miss,viol = is_constrained_delaunay(T)
        @test ok
        @test miss==0
        for s in segs; @test edge_in(T, s[1], s[2]); end
    end

    @testset "crossing constraints are rejected" begin
        xs=Float64[0,2,0,2]; ys=Float64[0,2,2,0]
        # segments (1,2)=(0,0)-(2,2) and (3,4)=(0,2)-(2,0) cross at center
        @test_throws ArgumentError constrained_delaunay(xs,ys,[(1,2),(3,4)])
    end
end

# minimum interior angle (degrees) over a triangle Mesh — independent quality oracle
function mesh_min_angle(m)
    mn = 180.0
    for t in 1:ntris(m)
        p = (node(m,m.tris[1,t]), node(m,m.tris[2,t]), node(m,m.tris[3,t]))
        for i in 1:3
            a=p[i]; b=p[mod1(i+1,3)]; c=p[mod1(i+2,3)]
            v1=(b[1]-a[1],b[2]-a[2]); v2=(c[1]-a[1],c[2]-a[2])
            cang = clamp((v1[1]*v2[1]+v1[2]*v2[2])/(hypot(v1...)*hypot(v2...)), -1.0, 1.0)
            mn = min(mn, acosd(cang))
        end
    end
    return mn
end
mesh_max_tri_area(m) = maximum(triangle_area(node(m,m.tris[1,t]),node(m,m.tris[2,t]),node(m,m.tris[3,t]))
                               for t in 1:ntris(m); init=0.0)

@testset "Mesh2D Ruppert refinement (Stage 1)" begin

    @testset "parameter validation and explicit step-limit blocker" begin
        makeT() = constrained_delaunay(Float64[0,1,1,0], Float64[0,0,1,1],
                                       [(1,2),(2,3),(3,4),(4,1)])
        @test_throws ArgumentError refine!(makeT(); min_angle_deg=NaN)
        @test_throws ArgumentError refine!(makeT(); min_angle_deg=-1)
        @test_throws ArgumentError refine!(makeT(); min_angle_deg=60)
        @test_throws ArgumentError refine!(makeT(); max_area=0)
        @test_throws ArgumentError refine!(makeT(); max_area=NaN)
        @test_throws ArgumentError refine!(makeT(); maxsteps=0)
        @test_throws ArgumentError refine!(makeT(); min_angle_deg=0,
                                           size=(x,y)->NaN, maxsteps=1)
        @test_throws ArgumentError refine!(makeT(); min_angle_deg=0,
                                           size=(x,y)->"bad", maxsteps=1)
        @test_throws ErrorException refine!(makeT(); min_angle_deg=0,
                                            max_area=0.01, maxsteps=1)
    end

    @testset "angle bound achieved, area & constraints preserved" begin
        r = _RNG(0xBEEF01)
        xs=Float64[0,10,10,0]; ys=Float64[0,0,10,10]
        for _ in 1:20; push!(xs, 1+_nf(r)*8); push!(ys, 1+_nf(r)*8); end
        segs=[(1,2),(2,3),(3,4),(4,1)]
        T=constrained_delaunay(xs,ys,segs;rng_seed=3)
        m0=to_mesh(T; interior=classify_interior(T))
        @test mesh_min_angle(m0) < 20.0            # starts skinny
        interior=refine!(T; min_angle_deg=20.0)
        @test check_consistency(T)[1]
        @test is_constrained_delaunay(T)[1]
        m=to_mesh(T; interior=interior)
        @test validate(m).ok
        @test mesh_min_angle(m) >= 20.0 - 1e-6     # angle guarantee met
        @test mesh_area(m) ≈ 100.0 atol=1e-9       # domain area preserved
    end

    @testset "area-bounded refinement" begin
        xs=Float64[0,10,10,0]; ys=Float64[0,0,10,10]; segs=[(1,2),(2,3),(3,4),(4,1)]
        T=constrained_delaunay(xs,ys,segs)
        interior=refine!(T; min_angle_deg=20.0, max_area=2.0)
        m=to_mesh(T; interior=interior)
        @test mesh_max_tri_area(m) <= 2.0 + 1e-9
        @test mesh_min_angle(m) >= 20.0 - 1e-6
        @test mesh_area(m) ≈ 100.0 atol=1e-9
        @test is_constrained_delaunay(T)[1]
    end

    @testset "refinement respects a hole" begin
        # interior points kept well clear of both the outer boundary and the hole
        # (Ruppert's angle guarantee requires input to respect local feature size).
        r = _RNG(0xABC)
        xs=Float64[0,10,10,0, 4,6,6,4]; ys=Float64[0,0,10,10, 4,4,6,6]
        for _ in 1:30
            x=1.0+_nf(r)*8.0; y=1.0+_nf(r)*8.0
            (3.0 < x < 7.0 && 3.0 < y < 7.0) && continue   # ≥1 from the hole
            push!(xs,x); push!(ys,y)
        end
        segs=[(1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5)]
        T=constrained_delaunay(xs,ys,segs;rng_seed=1)
        interior=refine!(T; min_angle_deg=20.0)
        m=to_mesh(T; interior=interior)
        @test validate(m).ok
        @test mesh_min_angle(m) >= 20.0 - 1e-6
        @test mesh_area(m) ≈ 96.0 atol=1e-9        # hole preserved through refinement
        @test is_constrained_delaunay(T)[1]
    end
end
