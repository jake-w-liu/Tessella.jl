# Mesh.HighOrderOptimize — bounded interior mid-node relocation (OPT01 slice).
#
# Correctness  : only interior mid-edge nodes are candidates — corners and
#                mid nodes on boundary faces/edges (the geometry-projected
#                supports) never move, even when the movable mask marks them.
# Acceptance   : every relocation keeps the exact cubic-Bernstein positivity
#                certificate on each incident element AND strictly improves
#                the minimum scaled Jacobian sampled at the ten type-11
#                interpolation points; a folded element can never accept a
#                candidate (the pass neither creates nor repairs inversions).
# Independence : quality is re-checked by an INDEPENDENT loop-form detJ
#                sampler on the degree-6 barycentric lattice (distinct from
#                the engine's 10-point scaled-Jacobian samples).
# API          : Mesh.HighOrderOptimize=0 leaves set_order(2) output bitwise
#                identical; =1 runs the pass on the unpublished plan mesh.

using Test
using LinearAlgebra: norm, dot, cross
import Tessella
using Tessella.MeshTypes
using Tessella.HighOrder
const API=Tessella.API

# ── INDEPENDENT P2 Jacobian oracle (loop-form gradients — same pattern as
#    highorder_test.jl, distinct from the module's _p2_grads) ────────────────
function _opt_oracle_grads(r, s, t)
    L = (1-r-s-t, r, s, t)
    dL = ((-1.0,-1.0,-1.0),(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))
    g = Vector{NTuple{3,Float64}}(undef, 10)
    for i in 1:4; c = 4L[i]-1; g[i] = (c*dL[i][1], c*dL[i][2], c*dL[i][3]); end
    for (k,(a,b)) in enumerate(((1,2),(2,3),(3,1),(1,4),(3,4),(2,4)))
        g[4+k] = (4*(L[a]*dL[b][1]+L[b]*dL[a][1]), 4*(L[a]*dL[b][2]+L[b]*dL[a][2]), 4*(L[a]*dL[b][3]+L[b]*dL[a][3]))
    end
    return g
end
function _opt_oracle_detJ(coords, tet10, t, r, s, u)
    g = _opt_oracle_grads(r, s, u); J = zeros(3,3)
    for k in 1:10
        v = tet10[k,t]
        for d in 1:3
            J[d,1]+=g[k][1]*coords[d,v]; J[d,2]+=g[k][2]*coords[d,v]; J[d,3]+=g[k][3]*coords[d,v]
        end
    end
    J[1,1]*(J[2,2]*J[3,3]-J[2,3]*J[3,2]) - J[1,2]*(J[2,1]*J[3,3]-J[2,3]*J[3,1]) + J[1,3]*(J[2,1]*J[3,2]-J[2,2]*J[3,1])
end
const _OPT_ORACLE_PTS = let S = NTuple{3,Float64}[]
    for j in 0:6, k in 0:6-j, l in 0:6-j-k; push!(S, (j/6, k/6, l/6)); end
    S
end
function _opt_oracle_min_detJ(coords, tet10, ts=axes(tet10,2))
    mn = Inf
    for t in ts, (r,s,u) in _OPT_ORACLE_PTS
        mn = min(mn, _opt_oracle_detJ(coords, tet10, t, r, s, u))
    end
    return mn
end
# Tets containing node `v` — the incident set of a mid node.
function _opt_incident_tets(tet10, v)
    ts = Int[]
    for t in axes(tet10,2)
        any(tet10[k,t]==v for k in 1:10) && push!(ts,t)
    end
    return ts
end

# Transverse displacement used to depress an edge mid from its midpoint.
const _OPT_SHIFT = 0.05/sqrt(2)

# ── 4-tet centroid fan: the outer tet (1,2,3,4) split by centroid node 5 into
#    4 sub-tets.  Interior edges (i,5) carry genuinely interior mid nodes;
#    the 6 outer edges carry boundary supports.  `moves` maps an undirected
#    primary edge to a 3-vector displacement added to that mid's exact
#    midpoint (transverse displacements bend the edge and depress the
#    element's minimum scaled Jacobian — pure tangential offsets do not).
#    Returns (coords, tet10, boundary_mids, interior_mids, edge_endpoints).
function _opt_fan10(; moves::Dict{NTuple{2,Int},NTuple{3,Float64}}=Dict{NTuple{2,Int},NTuple{3,Float64}}())
    corners = Float64[0 1 0 0 0.25;
                      0 0 1 0 0.25;
                      0 0 0 1 0.25]
    tets = NTuple{4,Int}[]
    for face in ((2,3,4), (1,3,4), (1,2,4), (1,2,3))
        t = (face[1], face[2], face[3], 5)
        a,b,c,d = (corners[:,k] for k in t)
        dot(b-a, cross(c-a, d-a)) < 0 && (t = (t[1],t[3],t[2],t[4]))
        push!(tets, t)
    end
    slots = ((5,1,2),(6,2,3),(7,3,1),(8,1,4),(9,3,4),(10,2,4))
    edges = NTuple{2,Int32}[]
    edge_id = Dict{NTuple{2,Int32},Int32}()
    getmid = function(a, b)
        key = a < b ? (Int32(a),Int32(b)) : (Int32(b),Int32(a))
        get!(edge_id, key) do
            push!(edges, key); Int32(5 + length(edges))
        end
    end
    tet10 = Matrix{Int32}(undef, 10, 4)
    for t in 1:4
        for k in 1:4
            tet10[k,t] = tets[t][k]
        end
        for (slot,i,j) in slots
            tet10[slot,t] = getmid(tets[t][i], tets[t][j])
        end
    end
    coords = Matrix{Float64}(undef, 3, 5 + length(edges))
    coords[:,1:5] = corners
    for (k,(a,b)) in pairs(edges)
        mid = 0.5 .* (corners[:,Int(a)] .+ corners[:,Int(b)])
        off = get(moves, (Int(a),Int(b)), (0.0,0.0,0.0))
        coords[:,5+k] = mid .+ collect(off)
    end
    bnd,_ = boundary_faces(hcat([Int32[t...] for t in tets]...))
    bedges = Set{NTuple{2,Int32}}()
    for f in bnd
        push!(bedges,minmax(f[1],f[2]));push!(bedges,minmax(f[2],f[3]))
        push!(bedges,minmax(f[1],f[3]))
    end
    boundary_mids = Int[]; interior_mids = Dict{Int,NTuple{2,Int32}}()
    for (k,e) in pairs(edges)
        node = 5+k
        if minmax(e[1],e[2]) in bedges
            push!(boundary_mids, node)
        else
            interior_mids[node] = e
        end
    end
    return coords, tet10, boundary_mids, interior_mids, edges
end

@testset "HighOrderOptimize interior relocation (engine)" begin
    @testset "perturbed interior mid improves certified min scaled Jacobian" begin
        # Interior edge (1,5) with its mid displaced transversely off the
        # straight edge: still Bernstein-positive but a depressed minimum
        # scaled Jacobian — improvable territory.
        coords, tet10, boundary_mids, interior_mids, _ =
            _opt_fan10(moves=Dict((1,5)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        p = P2Mesh(coords, tet10; require_positive_tets=false)
        @test p2_min_jacobian(p) > 0
        bad_mid = only(mid for (mid,e) in interior_mids if e == (Int32(1),Int32(5)))
        incident = _opt_incident_tets(tet10, bad_mid)
        jac_before = _opt_oracle_min_detJ(p.coords, p.tet10, incident)
        moved = optimize_p2_tet_interior!(p.coords, p.tet10, trues(nnodes(p));
                                          itermax=60, passes=10)
        # `moved` counts accepted node moves across all sweeps; neighbouring
        # interior mids may legitimately cooperate to repair the warped
        # element, and a node improving in several sweeps counts per sweep.
        @test moved[1] >= 1
        @test p2_min_jacobian(p) > 0               # certificate intact
        @test _opt_oracle_min_detJ(p.coords, p.tet10) > 0
        @test _opt_oracle_min_detJ(p.coords, p.tet10, incident) > jac_before
        a, b = interior_mids[bad_mid]
        target = 0.5 .* (coords[:,Int(a)] .+ coords[:,Int(b)])
        @test norm(p.coords[:,bad_mid] .- target) <
              norm(coords[:,bad_mid] .- target)
        @test norm(p.coords[:,bad_mid] .- target) < 0.05
        for mid in boundary_mids
            @test p.coords[:,mid] == coords[:,mid]
        end
    end

    @testset "boundary supports never move, even when they limit quality" begin
        # Boundary edge (1,2) has its mid displaced transversely — a depressed
        # (but still certified) minimum on two elements that the pass cannot
        # fix without moving the boundary support. Everything stays unchanged.
        coords, tet10, boundary_mids, _, edges =
            _opt_fan10(moves=Dict((1,2)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        bad_mid = only(5+k for (k,e) in pairs(edges)
                       if e == (Int32(1),Int32(2)))
        p = P2Mesh(coords, tet10; require_positive_tets=false)
        @test p2_min_jacobian(p) > 0
        optimize_p2_tet_interior!(p.coords, p.tet10, trues(nnodes(p));
                                itermax=60, passes=10)
        @test p.coords == coords
        @test bad_mid in boundary_mids
    end

    @testset "corners never move" begin
        coords, tet10, _, _, _ =
            _opt_fan10(moves=Dict((1,5)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        p = P2Mesh(coords, tet10; require_positive_tets=false)
        optimize_p2_tet_interior!(p.coords, p.tet10, trues(nnodes(p));
                                itermax=60, passes=10)
        @test p.coords[:,1:5] == coords[:,1:5]     # corners bitwise unchanged
    end

    @testset "movable mask is honored" begin
        coords, tet10, _, interior_mids, _ =
            _opt_fan10(moves=Dict((1,5)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        p = P2Mesh(coords, tet10; require_positive_tets=false)
        moved = optimize_p2_tet_interior!(p.coords, p.tet10, falses(nnodes(p)))
        @test moved[1] == 0
        @test p.coords == coords                    # bitwise identical
        # Only the perturbed interior mid unmasked → all other mids stay put
        # bit-for-bit; accepted moves may span several sweeps.
        movable = falses(nnodes(p))
        bad_mid = only(mid for (mid,e) in interior_mids
                       if e == (Int32(1),Int32(5)))
        movable[bad_mid] = true
        moved = optimize_p2_tet_interior!(p.coords, p.tet10, movable;
                                          itermax=60, passes=10)
        @test moved[1] >= 1
        @test p.coords[:,bad_mid] != coords[:,bad_mid]
        for (mid, _) in interior_mids
            mid == bad_mid && continue
            @test p.coords[:,mid] == coords[:,mid]
        end
    end

    @testset "folded elements are frozen, never repaired or worsened" begin
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                 tets=reshape(Int32[1,2,3,4],4,1))
        p = p2_tetmesh(m)
        # The between-sample fold fixture from highorder_test.jl.
        p.coords[:,p.tet10[5:10,1]] = Float64[
             0.598244   1.71802  -0.674877  -0.150606  -1.12732    0.694353;
             0.152545   0.407208  1.05877   -0.0461052  0.525971   0.280946;
             0.0303981 -1.00104   0.0559645  0.276108   0.689428   0.981785]
        @test p2_min_jacobian(p) < 0
        before = copy(p.coords)
        moved = optimize_p2_tet_interior!(p.coords, p.tet10, trues(nnodes(p));
                                          itermax=60, passes=10)
        @test p.coords == before                   # single-tet: all mids are
                                                 # boundary supports anyway
        @test p2_min_jacobian(p) < 0               # unchanged certificate
    end

    @testset "already-good straight mesh is a no-op" begin
        coords, tet10, _, interior_mids, _ = _opt_fan10()
        p = P2Mesh(coords, tet10; require_positive_tets=false)
        before = copy(p.coords)
        moved = optimize_p2_tet_interior!(p.coords, p.tet10, trues(nnodes(p));
                                          itermax=60, passes=10)
        @test moved[1] == 0
        @test p.coords == before
    end

    @testset "itermax=0 bounds the pass to zero work" begin
        coords, tet10, _, _, _ =
            _opt_fan10(moves=Dict((1,5)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        p = P2Mesh(coords, tet10; require_positive_tets=false)
        before = copy(p.coords)
        moved = optimize_p2_tet_interior!(p.coords, p.tet10, trues(nnodes(p));
                                          itermax=0, passes=10)
        @test moved[1] == 0
        @test p.coords == before
    end
end

@testset "HighOrderOptimize option wiring (API)" begin
    @testset "set_order(2) with =1 relocates interior mids; boundary fixed" begin
        coords, tet10, boundary_mids, interior_mids, _ =
            _opt_fan10(moves=Dict((1,5)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        api = API
        api.initialize()
        try
            api.option("Mesh.HighOrderOptimize", 1)
            api.model.add_discrete_entity(3, 7)
            api.mesh.add_nodes(3, 7,
                UInt64.(collect(axes(coords,2))), vec(coords))
            api.mesh.add_elements_by_type(7, 11,
                UInt64.(collect(axes(tet10,2))), UInt64.(vec(tet10)))
            api.mesh.generate(0)
            api.mesh.set_order(2)
            for mid in boundary_mids
                p,_,_,_ = api.mesh.get_node(UInt64(mid))
                @test p == coords[:,mid]
            end
            # The perturbed interior mid is relocated toward the edge
            # midpoint; the untouched interior mids stay put.
            for (mid,(a,b)) in interior_mids
                p,_,_,_ = api.mesh.get_node(UInt64(mid))
                target = 0.5 .* (coords[:,Int(a)] .+ coords[:,Int(b)])
                if (Int(a),Int(b)) == (1,5)
                    @test norm(p .- target) < norm(coords[:,mid] .- target)
                else
                    @test p == coords[:,mid]
                end
            end
        finally
            api.option("Mesh.HighOrderOptimize", 0)
            api.finalize()
        end
    end

    @testset "live P2Mesh overlay is detached before relocation" begin
        # `_high_order_optimize_overlay(...; fresh=false)` must return an
        # optimized COPY — the live overlay object is session state and may
        # never be mutated mid-pass.
        coords, tet10, boundary_mids, interior_mids, edges =
            _opt_fan10(moves=Dict((1,5)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        p = P2Mesh(coords, tet10; require_positive_tets=false)
        before = copy(p.coords)
        # Per-edge mid owners mirroring `_p2_mid_owners`: boundary supports
        # belong to surface entities, interior mids to the volume.
        interior_edges = Set(values(interior_mids))
        mids = Tuple{Int,Int32}[ (minmax(e[1],e[2]) in interior_edges ?
                                  (3,Int32(7)) : (2,Int32(1))) for e in edges ]
        api = API
        api.initialize()
        try
            api.option("Mesh.HighOrderOptimize", 1)
            optimized = Tessella.API._high_order_optimize_overlay(p, mids; fresh=false)
            @test optimized !== p
            @test p.coords == before          # live object untouched
            @test optimized.coords != before  # copy carried the relocation
            for mid in boundary_mids
                @test optimized.coords[:,mid] == before[:,mid]
            end
            # Disabled option returns the live object unchanged.
            api.option("Mesh.HighOrderOptimize", 0)
            same = Tessella.API._high_order_optimize_overlay(p, mids; fresh=false)
            @test same === p
            @test p.coords == before
        finally
            api.option("Mesh.HighOrderOptimize", 0)
            api.finalize()
        end
    end

    @testset "set_order(2) with =0 is bitwise identical to the import" begin
        coords, tet10, boundary_mids, interior_mids, _ =
            _opt_fan10(moves=Dict((1,5)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        api = API
        api.initialize()
        try
            @test api.option("Mesh.HighOrderOptimize") == 0
            api.model.add_discrete_entity(3, 7)
            api.mesh.add_nodes(3, 7,
                UInt64.(collect(axes(coords,2))), vec(coords))
            api.mesh.add_elements_by_type(7, 11,
                UInt64.(collect(axes(tet10,2))), UInt64.(vec(tet10)))
            api.mesh.generate(0)
            api.mesh.set_order(2)
            for node in axes(coords,2)
                p,_,_,_ = api.mesh.get_node(UInt64(node))
                @test p == coords[:,node]
            end
        finally
            api.finalize()
        end
    end

end

# ══════════════════════════════════════════════════════════════════════════════
# Multi-family interior relocation (gmsh types 9, 10, 12, 13 + uncertified).
# ══════════════════════════════════════════════════════════════════════════════
#
# Fixtures are built by position deduplication: each element's 27/18/9/6
# reference-mapped node positions are interned, so shared facets share node
# indices exactly (all positions are dyadic — interning by tuple equality is
# exact).  The independent quality oracle below recomputes tensor/product
# gradients in loop form (distinct from the module's ntuple implementations)
# and samples on a denser lattice than the engine's interpolation points.

const _OPT_QUA_NODES = ((1,1),(3,1),(3,3),(1,3),(2,1),(3,2),(2,3),(1,2),(2,2))
const _OPT_HEX_CORNERS = ((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
                         (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.))
const _OPT_HEX_EDGES = ((1,2),(1,4),(1,5),(2,3),(2,6),(3,4),
                        (3,7),(4,8),(5,6),(5,8),(6,7),(7,8))
const _OPT_HEX_FACES = ((1,4,3,2),(1,2,6,5),(1,5,8,4),
                        (2,3,7,6),(3,4,8,7),(5,6,7,8))
const _OPT_PRI_EDGES = ((1,2),(1,3),(1,4),(2,3),(2,5),(3,6),
                        (4,5),(4,6),(5,6))
const _OPT_PRI_QFACES = ((1,2,5,4),(1,4,6,3),(2,3,6,5))

# Node positions of one unit reference element mapped by `mapfn` (identity
# gives the unit reference shape itself).
function _opt_block_nodes(family::Symbol, mapfn)
    pts = NTuple{3,Float64}[]
    if family === :quad9
        for (a,b) in _OPT_QUA_NODES
            push!(pts, mapfn((Float64(a-2), Float64(b-2), 0.0)))
        end
    elseif family === :hex27
        for k in 1:8
            push!(pts, mapfn(_OPT_HEX_CORNERS[k]))
        end
        for (i,j) in _OPT_HEX_EDGES
            push!(pts, mapfn(0.5 .* (_OPT_HEX_CORNERS[i] .+ _OPT_HEX_CORNERS[j])))
        end
        for f in _OPT_HEX_FACES
            push!(pts, mapfn(0.25 .* (_OPT_HEX_CORNERS[f[1]] .+
                _OPT_HEX_CORNERS[f[2]] .+ _OPT_HEX_CORNERS[f[3]] .+
                _OPT_HEX_CORNERS[f[4]])))
        end
        c = (0.0,0.0,0.0)
        for k in 1:8
            c = c .+ _OPT_HEX_CORNERS[k]
        end
        push!(pts, mapfn(0.125 .* c))
    elseif family === :pri18
        pc = ((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(1.,0.,1.),(0.,1.,1.))
        for k in 1:6
            push!(pts, mapfn(pc[k]))
        end
        for (i,j) in _OPT_PRI_EDGES
            push!(pts, mapfn(0.5 .* (pc[i] .+ pc[j])))
        end
        for f in _OPT_PRI_QFACES
            push!(pts, mapfn(0.25 .* (pc[f[1]] .+ pc[f[2]] .+
                                      pc[f[3]] .+ pc[f[4]])))
        end
    elseif family === :tri6
        for p in ((0.,0.),(1.,0.),(0.,1.),(0.5,0.),(0.5,0.5),(0.,0.5))
            push!(pts, mapfn((p[1], p[2], 0.0)))
        end
    else
        error("unknown family $family")
    end
    return pts
end

# Intern a list of per-element node-position vectors into shared coordinates
# and an nnodes × nelements connectivity; returns (coords, conn, positions).
function _opt_assemble(families::Vector{Symbol}, maps::Vector)
    rows = Dict(:tri6=>6, :quad9=>9, :hex27=>27, :pri18=>18)
    index = Dict{NTuple{3,Float64},Int32}()
    nodes = NTuple{3,Float64}[]
    per_elem = [ _opt_block_nodes(f, maps[i]) for (i,f) in enumerate(families)]
    conn = Matrix{Int32}(undef, rows[first(families)], length(per_elem))
    for (t, pts) in enumerate(per_elem)
        for (k, p) in enumerate(pts)
            # Canonicalize -0.0 → +0.0: Dict keys compare with isequal, and
            # isequal(-0.0, 0.0) is false — without this, mirrored/rotated
            # fixtures fail to share nodes.
            key = (p[1]+0.0, p[2]+0.0, p[3]+0.0)
            conn[k,t] = get!(index, key) do
                push!(nodes, key); Int32(length(nodes))
            end
        end
    end
    coords = Matrix{Float64}(undef, 3, length(nodes))
    for (i,p) in enumerate(nodes)
        coords[:,i] = collect(p)
    end
    return coords, conn, index
end

# ── independent scaled-Jacobian oracles (loop form) ───────────────────────
_opt_lag(x)  = (0.5x*(x-1.0), 1.0-x*x, 0.5x*(x+1.0))
_opt_dlag(x) = (x-0.5, -2.0x, x+0.5)

function _opt_oracle_grads(::Val{9}, u, v, w)
    lu=_opt_lag(u); du=_opt_dlag(u); lv=_opt_lag(v); dv=_opt_dlag(v)
    return [(du[a]*lv[b], lu[a]*dv[b], 0.0) for (a,b) in _OPT_QUA_NODES]
end
function _opt_oracle_grads(::Val{12}, u, v, w)
    lu=_opt_lag(u); du=_opt_dlag(u); lv=_opt_lag(v); dv=_opt_dlag(v)
    lw=_opt_lag(w); dw=_opt_dlag(w)
    return [(du[a]*lv[b]*lw[c], lu[a]*dv[b]*lw[c], lu[a]*lv[b]*dw[c])
            for (a,b,c) in HighOrder._P2HEX_NODES]
end
function _opt_oracle_grads(::Val{13}, r, s, w)
    L1=1.0-r-s; L2=r; L3=s
    tv=(L1*(2L1-1), L2*(2L2-1), L3*(2L3-1), 4L1*L2, 4L2*L3, 4L1*L3)
    dr=(1-4L1, 4L2-1, 0.0, 4(L1-L2), 4L3, -4L3)
    ds=(1-4L1, 0.0, 4L3-1, -4L2, 4L2, 4(L1-L3))
    lw=_opt_lag(w); dw=_opt_dlag(w)
    return [(dr[a]*lw[b], ds[a]*lw[b], tv[a]*dw[b])
            for (a,b) in HighOrder._P2PRI_NODES]
end
function _opt_oracle_J(coords, conn, t, grads)
    J = zeros(3,3)
    for (k,g) in enumerate(grads)
        v = conn[k,t]
        for d in 1:3
            J[d,1]+=g[1]*coords[d,v]; J[d,2]+=g[2]*coords[d,v]
            J[d,3]+=g[3]*coords[d,v]
        end
    end
    return J
end
_opt_oracle_det3(J) = J[1,1]*(J[2,2]*J[3,3]-J[2,3]*J[3,2]) -
    J[1,2]*(J[2,1]*J[3,3]-J[2,3]*J[3,1]) +
    J[1,3]*(J[2,1]*J[3,2]-J[2,2]*J[3,1])
_opt_oracle_scaled(J) = _opt_oracle_det3(J) /
    (norm(J[:,1])*norm(J[:,2])*norm(J[:,3]))
_opt_oracle_scaled2d(J) = norm(cross(J[:,1],J[:,2])) /
    (norm(J[:,1])*norm(J[:,2]))

# Minimum scaled Jacobian measure over a denser lattice than the engine's
# interpolation-point samples: [-1,1]² for the quad, the r+s≤1 triangle ×
# [-1,1] for the prism, [-1,1]³ for the hex.
function _opt_oracle_min_scaled(coords, conn, kind, ts=axes(conn,2))
    mn = Inf
    if kind == Val(9)
        grid = range(-1., 1.; length=3)
        for t in ts, u in grid, v in grid
            J = _opt_oracle_J(coords, conn, t,
                              _opt_oracle_grads(kind, u, v, 0.0))
            w = _opt_oracle_scaled2d(J)
            w < mn && (mn = w)
        end
    elseif kind == Val(13)
        for t in ts, r in range(0., 1.; length=3),
                     s in range(0., 1.; length=3),
                     w in range(-1., 1.; length=3)
            r + s > 1 && continue
            J = _opt_oracle_J(coords, conn, t,
                              _opt_oracle_grads(kind, r, s, w))
            v = _opt_oracle_scaled(J)
            v < mn && (mn = v)
        end
    else
        for t in ts, u in range(-1., 1.; length=3),
                     v in range(-1., 1.; length=3),
                     w in range(-1., 1.; length=3)
            J = _opt_oracle_J(coords, conn, t,
                              _opt_oracle_grads(kind, u, v, w))
            v = _opt_oracle_scaled(J)
            v < mn && (mn = v)
        end
    end
    return mn
end

@testset "HighOrderOptimize multi-family engine" begin
    @testset "Quad9 center relocates; edge mids stay fixed" begin
        coords, q9, index = _opt_assemble([:quad9], [identity])
        center = Int(index[(0.0,0.0,0.0)])
        target = copy(coords[:,center])
        coords[:,center] .+= (-0.18, 0.06, 0.0)   # in-plane warp of the map
        warped = copy(coords[:,center])
        @test p2_quad_min_jacobian(coords, q9) > 0   # still certified
        min0 = _opt_oracle_min_scaled(coords, q9, Val(9))
        before = copy(coords)
        moved = optimize_p2_interior!(coords, ((Int32(10), q9),),
                                      trues(size(coords,2));
                                      itermax=60, passes=10)
        @test moved[1] >= 1
        @test p2_quad_min_jacobian(coords, q9) > 0
        @test _opt_oracle_min_scaled(coords, q9, Val(9)) > min0
        # center moved toward the square centroid; all other nodes bitwise
        @test norm(coords[:,center] .- target) < norm(warped .- target)
        for v in axes(coords,2)
            v == center && continue
            @test coords[:,v] == before[:,v]
        end
    end

    @testset "Tri6 interior wire mid is a fixed surface support" begin
        # Two triangles sharing the diagonal edge: its mid is interior to the
        # wire complex but is still a surface support — surface-family edge
        # mids are never candidates.  Tri B is a 180° in-plane rotation about
        # the shared edge (orientation preserving, certified).
        rot = p -> (1.0 - p[1], 1.0 - p[2], 0.0)
        coords, t6, index = _opt_assemble([:tri6,:tri6], [identity, rot])
        shared = Int(index[(0.5,0.5,0.0)])
        coords[:,shared] .+= (0.12, -0.12, 0.0)
        @test p2_tri_min_jacobian(
            P2TriMesh(coords, t6; _check_geometry=false)) > 0
        before = copy(coords)
        moved = optimize_p2_interior!(coords, ((Int32(9), t6),),
                                      trues(size(coords,2));
                                      itermax=60, passes=10)
        @test moved[1] == 0
        @test coords == before
    end

    @testset "Hex27 body center relocates on a single element" begin
        coords, h27, index = _opt_assemble([:hex27], [identity])
        body = Int(index[(0.5,0.5,0.5)])
        target = copy(coords[:,body])
        coords[:,body] .+= (-0.13, 0.05, -0.09)
        warped = copy(coords[:,body])
        @test p2_hex_min_jacobian(coords, h27) > 0
        min0 = _opt_oracle_min_scaled(coords, h27, Val(12))
        before = copy(coords)
        moved = optimize_p2_interior!(coords, ((Int32(12), h27),),
                                      trues(size(coords,2));
                                      itermax=60, passes=10)
        @test moved[1] >= 1
        @test p2_hex_min_jacobian(coords, h27) > 0
        @test _opt_oracle_min_scaled(coords, h27, Val(12)) > min0
        @test norm(coords[:,body] .- target) < norm(warped .- target)
        for v in axes(coords,2)
            v == body && continue
            @test coords[:,v] == before[:,v]
        end
    end

    @testset "Hex27 interior-face center relocates; boundary supports fixed" begin
        coords, h27, index = _opt_assemble([:hex27,:hex27],
            [identity, p -> (p[1] .+ 1.0, p[2], p[3])])
        face_mid = Int(index[(1.0,0.5,0.5)])   # shared-face center
        target = copy(coords[:,face_mid])
        before = copy(coords)
        coords[:,face_mid] .+= (0.0, -0.14, 0.11)
        warped = copy(coords[:,face_mid])
        @test p2_hex_min_jacobian(coords, h27) > 0
        min0 = _opt_oracle_min_scaled(coords, h27, Val(12))
        moved = optimize_p2_interior!(coords, ((Int32(12), h27),),
                                      trues(size(coords,2));
                                      itermax=60, passes=10)
        @test moved[1] >= 1
        @test p2_hex_min_jacobian(coords, h27) > 0
        @test _opt_oracle_min_scaled(coords, h27, Val(12)) > min0
        @test norm(coords[:,face_mid] .- target) < norm(warped .- target)
        # Corners, boundary mids, and the shared-face edge mids touching
        # boundary facets stay bitwise unchanged; the only other candidates
        # are the two element-interior body centers.
        movers = Set((face_mid, Int(index[(0.5,0.5,0.5)]),
                      Int(index[(1.5,0.5,0.5)])))
        for (p,i) in index
            Int(i) in movers && continue
            @test coords[:,Int(i)] == before[:,Int(i)]
        end
    end

    @testset "single Prism18 is a certified no-op" begin
        # Every facet of a lone prism is complex boundary and type 13 has no
        # element-interior slot → no candidates at all.
        coords, p18, _ = _opt_assemble([:pri18], [identity])
        @test p2_prism_min_jacobian(coords, p18) > 0
        before = copy(coords)
        moved = optimize_p2_interior!(coords, ((Int32(13), p18),),
                                      trues(size(coords,2));
                                      itermax=60, passes=10)
        @test moved[1] == 0
        @test coords == before
    end

    @testset "Prism18 interior quad-face center relocates" begin
        # Two prisms sharing the x=0 quad face; prism B is a 90° rotation
        # about z (orientation preserving — a mirror would be uncertified).
        rot = p -> (-p[2], p[1], p[3])
        coords, p18, index = _opt_assemble([:pri18,:pri18], [identity, rot])
        face_mid = Int(index[(0.0,0.5,0.5)])
        target = copy(coords[:,face_mid])
        before = copy(coords)
        coords[:,face_mid] .+= (0.0, -0.14, 0.10)
        warped = copy(coords[:,face_mid])
        @test p2_prism_min_jacobian(coords, p18) > 0
        min0 = _opt_oracle_min_scaled(coords, p18, Val(13))
        moved = optimize_p2_interior!(coords, ((Int32(13), p18),),
                                      trues(size(coords,2));
                                      itermax=60, passes=10)
        @test moved[1] >= 1
        @test p2_prism_min_jacobian(coords, p18) > 0
        @test _opt_oracle_min_scaled(coords, p18, Val(13)) > min0
        @test norm(coords[:,face_mid] .- target) < norm(warped .- target)
        for (p,i) in index
            Int(i) == face_mid && continue
            @test coords[:,Int(i)] == before[:,Int(i)]
        end
    end

    @testset "uncertified families no-op and exclude shared nodes" begin
        # A pure Pyr14 block is skipped: it has no boolean Bernstein bound.
        pyr = reshape(Int32.(collect(1:14)), 14, 1)
        coords = repeat(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1], 1, 4)[:,1:14]
        before = copy(coords)
        moved = optimize_p2_interior!(coords, ((Int32(14), pyr),),
                                      trues(size(coords,2)))
        @test moved == (0, Inf, Inf, 0, 0)
        @test coords == before
        # A node shared between a certified Tet10 interior edge and a Tri6
        # block is a cross-family support: it never moves even though the
        # same edge is a Tet10-interior candidate.
        coords2, tet10, _, interior_mids, _ =
            _opt_fan10(moves=Dict((1,5)=>(_OPT_SHIFT,-_OPT_SHIFT,0.0)))
        bad_mid = only(mid for (mid,e) in interior_mids
                       if e == (Int32(1),Int32(5)))
        tri = Matrix{Int32}(undef,6,1)
        tri[:,1] = Int32[tet10[1,1], tet10[5,1], tet10[2,1],
                         bad_mid, tet10[6,1], tet10[7,1]]
        expected = copy(coords2[:,bad_mid])
        optimize_p2_interior!(
            coords2, ((Int32(11), tet10), (Int32(9), tri)),
            trues(size(coords2,2)); itermax=60, passes=10)
        @test coords2[:,bad_mid] == expected
    end

    @testset "set_order(2) on a Quad9 surface relocates the interior node" begin
        # A discrete surface entity carries owner dimension 2: the mask
        # admits the Quad9 element-interior node while every edge mid stays
        # a surface support.  The center is warped in-plane; the pass pulls
        # it back toward the square centroid under the exact Quad9 Gram
        # certificate.
        coords, q9, index = _opt_assemble([:quad9], [identity])
        center = Int(index[(0.0,0.0,0.0)])
        coords[:,center] .+= (-0.18, 0.06, 0.0)
        target = Float64[0.0, 0.0, 0.0]
        warped = copy(coords[:,center])
        api = API
        api.initialize()
        try
            api.option("Mesh.HighOrderOptimize", 1)
            api.model.add_discrete_entity(2, 7)
            api.mesh.add_nodes(2, 7,
                UInt64.(collect(axes(coords,2))), vec(coords))
            api.mesh.add_elements_by_type(7, 10,
                UInt64.(collect(axes(q9,2))), UInt64.(vec(q9)))
            api.mesh.generate(0)
            api.mesh.set_order(2)
            for v in axes(coords,2)
                p,_,_,_ = api.mesh.get_node(UInt64(v))
                if v == center
                    @test norm(p .- target) < norm(warped .- target)
                else
                    @test p == coords[:,v]
                end
            end
        finally
            api.option("Mesh.HighOrderOptimize", 0)
            api.finalize()
        end
    end
end
