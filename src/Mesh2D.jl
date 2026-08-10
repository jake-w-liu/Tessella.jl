"""
    Mesh2D

Stage-1 planar meshing (PLAN.md §3 "Mesh2D"): incremental **Delaunay**
triangulation (Bowyer–Watson) with a symbolic point at infinity, **constrained
Delaunay** (CDT) for planar straight-line graphs, and **Ruppert** quality
refinement.

Everything rides on the exact predicates in [`Predicates`](@ref) — orientation
and in-circle tests are `*_sos` (Simulation of Simplicity), so the kernel never
sees a genuine tie and its topological decisions stay globally consistent. That
exactness makes the empty-circumcircle property hold *exactly* (the oracle the
tests check).

Robustness note: rather than a finite bounding super-triangle (which is *not*
robust — a point outside the evolving hull can be wrongly excluded from a
super-adjacent triangle's circumcircle, cracking the boundary), the convex hull
is closed by **ghost triangles** sharing a single vertex `GHOST` at infinity.
A ghost triangle's "in-circumcircle" test degenerates to an orientation test on
its real hull edge — the exact, correct behaviour in the limit. Ghost triangles
are stripped at export.

The working triangulation is index-compact: two flat `Int32` arrays (`tv`
vertices, `tn` neighbours, 3 per triangle) plus a free list, so insert/delete is
`O(1)` amortized. [`check_consistency`](@ref) validates the half-edge adjacency
after operations in the test suite.
"""
module Mesh2D

using ..Predicates: orient2, orient2_sos, incircle_sos
using ..MeshTypes: Mesh

export Triangulation, triangulate, delaunay2d, dedup_points
export check_consistency, is_delaunay, to_mesh, ntriangles_live, insert_point!

const GHOST = Int32(0)   # the single vertex at infinity closing the convex hull

# ════════════════════════════════════════════════════════════════════════════════
# Working triangulation
# ════════════════════════════════════════════════════════════════════════════════

mutable struct Triangulation
    x::Vector{Float64}          # real point x, indices 1..nreal (GHOST=0 has no coords)
    y::Vector{Float64}
    tv::Vector{Int32}           # 3 per triangle: vertex ids
    tn::Vector{Int32}           # 3 per triangle: neighbour across edge opposite vertex k (0-slot = none)
    alive::Vector{Bool}
    freelist::Vector{Int32}
    nreal::Int
    last::Int32                 # walk hint
end

@inline _pt(T::Triangulation, i::Integer) = @inbounds (T.x[i], T.y[i])
@inline _vert(T::Triangulation, t::Integer, k::Integer) = @inbounds T.tv[3*(t-1)+k]
@inline _setv!(T::Triangulation, t::Integer, k::Integer, val) = @inbounds (T.tv[3*(t-1)+k] = Int32(val))
@inline _nbr(T::Triangulation, t::Integer, k::Integer) = @inbounds T.tn[3*(t-1)+k]
@inline _setn!(T::Triangulation, t::Integer, k::Integer, val) = @inbounds (T.tn[3*(t-1)+k] = Int32(val))
@inline _is_ghost_v(v) = v == GHOST

@inline function _ghost_slot(T::Triangulation, t::Integer)
    base = 3*(t-1)
    @inbounds (T.tv[base+1]==GHOST) && return 1
    @inbounds (T.tv[base+2]==GHOST) && return 2
    @inbounds (T.tv[base+3]==GHOST) && return 3
    return 0
end
@inline _is_ghost_tri(T::Triangulation, t::Integer) = _ghost_slot(T, t) != 0

ntriangles_live(T::Triangulation) = count(T.alive)

@inline function _newtri!(T::Triangulation, a, b, c)
    if !isempty(T.freelist)
        t = pop!(T.freelist); base = 3*(t-1)
        @inbounds T.tv[base+1]=a; T.tv[base+2]=b; T.tv[base+3]=c
        @inbounds T.tn[base+1]=0; T.tn[base+2]=0; T.tn[base+3]=0
        @inbounds T.alive[t]=true
        return Int32(t)
    else
        push!(T.tv, Int32(a), Int32(b), Int32(c))
        push!(T.tn, Int32(0), Int32(0), Int32(0))
        push!(T.alive, true)
        return Int32(length(T.alive))
    end
end

@inline function _killtri!(T::Triangulation, t)
    @inbounds T.alive[t] = false
    push!(T.freelist, Int32(t))
    return nothing
end

@inline function _nslot(T::Triangulation, t, s)
    base = 3*(t-1)
    @inbounds (T.tn[base+1]==s) && return 1
    @inbounds (T.tn[base+2]==s) && return 2
    @inbounds (T.tn[base+3]==s) && return 3
    return 0
end

# the two vertices of the edge opposite local vertex k, in the triangle's stored order
@inline function _edge(T::Triangulation, t, k)
    base = 3*(t-1)
    if k == 1
        @inbounds return (T.tv[base+2], T.tv[base+3])
    elseif k == 2
        @inbounds return (T.tv[base+3], T.tv[base+1])
    else
        @inbounds return (T.tv[base+1], T.tv[base+2])
    end
end

@inline function _nslot_by_edge(T::Triangulation, t, u, v)
    @inbounds for k in 1:3
        a, b = _edge(T, t, k)
        ((a==u && b==v) || (a==v && b==u)) && return k
    end
    return 0
end

@inline function _link_spoke!(T, spoke, w::Int32, t::Int32, slot::Int)
    if haskey(spoke, w)
        (t2, slot2) = spoke[w]
        _setn!(T, t, slot, t2)
        _setn!(T, t2, Int(slot2), t)
        delete!(spoke, w)
    else
        spoke[w] = (t, Int32(slot))
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# In-circle (ghost-aware) — the heart of robustness
# ════════════════════════════════════════════════════════════════════════════════

# Is real point `vid` inside triangle `t`'s circumdisk?  Real triangle → exact
# incircle_sos. Ghost triangle [.,.,GHOST] → its real hull edge (u,v) is stored so
# the solid is to the RIGHT of u→v; the ghost's circumdisk is the OUTER half-plane
# (left of u→v). The test is *collinear-aware* and uses exact `orient2` (not SoS):
# strictly outside ⇒ in; strictly inside ⇒ out; exactly on the edge's line ⇒ in
# only if `vid` lies strictly *between* u and v. This is what prevents flat
# (zero-area) triangles on collinear convex-hull edges while keeping `locate` and
# the cavity test mutually consistent (see the analysis in the module docs/tests).
@inline function _in_circumcircle(T::Triangulation, t, vid)
    gs = _ghost_slot(T, t)
    if gs == 0
        a = _vert(T,t,1); b = _vert(T,t,2); c = _vert(T,t,3)
        return incircle_sos(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,vid), a,b,c,Int32(vid)) > 0
    else
        u, v = _edge(T, t, gs)
        pu = _pt(T,u); pv = _pt(T,v); pp = _pt(T,vid)
        o = orient2(pu, pv, pp)
        o > 0 && return true
        o < 0 && return false
        return _strictly_between(pu, pv, pp)   # collinear ⇒ only mid-segment points
    end
end

@inline function _strictly_between(u, v, p)
    dx = v[1]-u[1]; dy = v[2]-u[2]
    s = (p[1]-u[1])*dx + (p[2]-u[2])*dy
    L = dx*dx + dy*dy
    return 0.0 < s < L
end

# ════════════════════════════════════════════════════════════════════════════════
# Construction & initialization
# ════════════════════════════════════════════════════════════════════════════════

function Triangulation(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs); @assert length(ys) == n
    return Triangulation(xs, ys, Int32[], Int32[], Bool[], Int32[], n, Int32(0))
end

# Build the seed: first non-collinear triple → 1 real triangle ringed by 3 ghosts.
# Returns the set of the (up to 3) real vertices already placed, or an empty set
# if all points are collinear (no triangle possible).
function _init_triangulation!(T::Triangulation)
    n = T.nreal
    n < 3 && return Int32[]
    a = Int32(1)
    b = Int32(0)
    @inbounds for i in 2:n
        if (T.x[i], T.y[i]) != (T.x[a], T.y[a]); b = Int32(i); break; end
    end
    b == 0 && return Int32[]
    c = Int32(0)
    @inbounds for i in 2:n
        (i == b) && continue
        if orient2(_pt(T,a), _pt(T,b), _pt(T,i)) != 0; c = Int32(i); break; end
    end
    c == 0 && return Int32[]
    # orient CCW
    if orient2(_pt(T,a), _pt(T,b), _pt(T,c)) < 0
        b, c = c, b
    end
    Rt = _newtri!(T, a, b, c)                     # [a,b,c] CCW
    # ghosts across each hull edge, stored reversed so solid is right of the edge
    Ga = _newtri!(T, c, b, GHOST)                 # across edge opp a (b,c) → reversed (c,b)
    Gb = _newtri!(T, a, c, GHOST)                 # across edge opp b (c,a) → reversed (a,c)
    Gc = _newtri!(T, b, a, GHOST)                 # across edge opp c (a,b) → reversed (b,a)
    _setn!(T, Rt, 1, Ga); _setn!(T, Ga, 3, Rt)
    _setn!(T, Rt, 2, Gb); _setn!(T, Gb, 3, Rt)
    _setn!(T, Rt, 3, Gc); _setn!(T, Gc, 3, Rt)
    # link the three ghosts into a ring by their shared REAL hull vertex
    spoke = Dict{Int32, Tuple{Int32,Int32}}()
    # Ga=[c,b,GHOST]: slot1 (opp c)=edge(b,GHOST) shares vertex b; slot2 (opp b)=edge(GHOST,c) shares c
    _link_spoke!(T, spoke, b, Ga, 1); _link_spoke!(T, spoke, c, Ga, 2)
    # Gb=[a,c,GHOST]: slot1 (opp a)=edge(c,GHOST) shares c; slot2 (opp c)=edge(GHOST,a) shares a
    _link_spoke!(T, spoke, c, Gb, 1); _link_spoke!(T, spoke, a, Gb, 2)
    # Gc=[b,a,GHOST]: slot1 (opp b)=edge(a,GHOST) shares a; slot2 (opp a)=edge(GHOST,b) shares b
    _link_spoke!(T, spoke, a, Gc, 1); _link_spoke!(T, spoke, b, Gc, 2)
    isempty(spoke) || error("Mesh2D._init: ghost ring did not close")
    T.last = Rt
    return Int32[a, b, c]
end

# ════════════════════════════════════════════════════════════════════════════════
# Point location (Lawson walk, ghost-aware)
# ════════════════════════════════════════════════════════════════════════════════

# Remembering *stochastic* walk (Devillers–Pion–Teillaud): at each triangle test
# the three edges starting from a pseudo-random offset and cross the first edge
# `p` is outside of, never immediately back to `prev`. Randomizing the edge order
# breaks the cycles a purely deterministic visibility walk can fall into, giving
# guaranteed termination on a Delaunay triangulation. The offset stream is keyed
# on `vid`, so the walk (hence the whole triangulation) stays reproducible.
function locate(T::Triangulation, px, py, vid; start::Integer=T.last)
    t = start
    (t == 0 || !T.alive[t]) && (t = _first_alive(T))
    # step onto a real triangle to start (the real neighbour is opposite GHOST,
    # whose slot varies — insertion-created ghosts don't keep GHOST in slot 3)
    let gs = _ghost_slot(T, t)
        gs != 0 && (t = _nbr(T, t, gs))
    end
    prev = Int32(0)
    rs = UInt64(vid) * 0x9E3779B97F4A7C15 + 0xD1B54A32D192ED03
    guard = 0; maxsteps = 12 * length(T.alive) + 64
    @inbounds while true
        guard += 1
        guard > maxsteps && error("Mesh2D.locate: walk did not terminate")
        rs ⊻= rs << 13; rs ⊻= rs >> 7; rs ⊻= rs << 17
        r = Int(rs % UInt64(3))
        moved = false
        for i in 0:2
            k = (r + i) % 3 + 1
            nb = _nbr(T, t, k)
            (nb == prev || nb == 0) && continue
            if _is_ghost_tri(T, nb)
                # cross into the ghost only if it genuinely contains p (collinear-
                # aware); otherwise this hull edge is not the way out.
                _in_circumcircle(T, nb, vid) && return nb
                continue
            end
            a, b = _edge(T, t, k)
            # p strictly right of CCW edge a→b ⇒ outside t across edge k
            if orient2_sos((T.x[a],T.y[a]),(T.x[b],T.y[b]),(px,py), a, b, Int32(vid)) < 0
                prev = Int32(t); t = nb; moved = true; break
            end
        end
        moved || return t
    end
end

function _first_alive(T::Triangulation)
    @inbounds for t in eachindex(T.alive)
        T.alive[t] && return Int32(t)
    end
    error("Mesh2D: no live triangles")
end

# ════════════════════════════════════════════════════════════════════════════════
# Bowyer–Watson insertion
# ════════════════════════════════════════════════════════════════════════════════

function insert_point!(T::Triangulation, vid::Integer)
    px, py = _pt(T, vid)
    t0 = locate(T, px, py, vid)
    cavity = Int32[t0]
    incav = Set{Int32}(); push!(incav, t0)
    boundary = Tuple{Int32,Int32,Int32}[]     # (edge_a, edge_b, outside neighbour)
    stack = Int32[t0]
    @inbounds while !isempty(stack)
        t = pop!(stack)
        for k in 1:3
            nb = _nbr(T, t, k)
            (nb != 0 && (nb in incav)) && continue
            inside = nb != 0 && _in_circumcircle(T, nb, vid)
            if inside
                push!(incav, nb); push!(cavity, nb); push!(stack, nb)
            else
                a, b = _edge(T, t, k)
                push!(boundary, (a, b, nb))
            end
        end
    end
    for t in cavity; _killtri!(T, t); end
    _retriangulate_cavity!(T, boundary, Int32(vid))
    return nothing
end

function _retriangulate_cavity!(T::Triangulation, boundary, vid::Int32)
    spoke = Dict{Int32, Tuple{Int32,Int32}}()   # boundary vertex w → (new tri, slot)
    anytri = Int32(0)
    for (a, b, nb) in boundary
        t = _newtri!(T, a, b, vid)               # [a, b, vid]; a or b may be GHOST
        anytri = t
        _setn!(T, t, 3, nb)                       # slot3 (opp vid) = edge (a,b) → nb
        if nb != 0
            j = _nslot_by_edge(T, nb, a, b)
            j != 0 && _setn!(T, nb, j, t)
        end
        _link_spoke!(T, spoke, a, t, 2)           # slot2 (opp b) = edge (vid,a), spoke at a
        _link_spoke!(T, spoke, b, t, 1)           # slot1 (opp a) = edge (b,vid), spoke at b
    end
    isempty(spoke) || error("Mesh2D: cavity boundary not a simple polygon (unmatched spoke)")
    anytri != 0 && (T.last = anytri)
    return nothing
end

# ════════════════════════════════════════════════════════════════════════════════
# Driver + export
# ════════════════════════════════════════════════════════════════════════════════

"""
    delaunay2d(xs, ys; rng_seed=1) -> Triangulation

Delaunay triangulation of the (deduplicated) 2-D points. Points are inserted in a
fixed-seed randomized order; the resulting triangulation is order-independent
(SoS makes even degenerate ties deterministic in the vertex indices).
"""
function delaunay2d(xs::Vector{Float64}, ys::Vector{Float64}; rng_seed::Integer=1)
    ux, uy, _ = dedup_points(xs, ys)
    T = Triangulation(ux, uy)
    placed = _init_triangulation!(T)
    isempty(placed) && return T                    # <3 non-collinear points → no triangles
    placedset = Set{Int32}(placed)
    order = Int32[i for i in 1:T.nreal if !(Int32(i) in placedset)]
    _shuffle_det!(order, UInt64(rng_seed))
    for v in order
        insert_point!(T, v)
    end
    return T
end

"""
    dedup_points(xs, ys) -> (ux, uy, remap)

Merge exactly-coincident points (Delaunay of duplicated points is undefined).
`remap[i]` is the unique-point index original point `i` maps to (for CDT segment
relabelling).
"""
function dedup_points(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs); @assert length(ys) == n
    seen = Dict{Tuple{Float64,Float64},Int32}()
    ux = Float64[]; uy = Float64[]; remap = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        key = (xs[i], ys[i])
        id = get(seen, key, Int32(0))
        if id == 0
            push!(ux, xs[i]); push!(uy, ys[i]); id = Int32(length(ux)); seen[key] = id
        end
        remap[i] = id
    end
    return ux, uy, remap
end

# deterministic Fisher–Yates via a splitmix64 stream (no RNG dep in the core)
function _shuffle_det!(a::Vector{Int32}, seed::UInt64)
    s = seed
    @inbounds for i in length(a):-1:2
        s += 0x9E3779B97F4A7C15
        z = s
        z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
        z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
        z = z ⊻ (z >> 31)
        j = Int(z % UInt64(i)) + 1
        a[i], a[j] = a[j], a[i]
    end
    return a
end

"""
    to_mesh(T; interior=nothing) -> Mesh

Export live *real* triangles (ghosts dropped) to a [`Mesh`](@ref), compacting the
referenced nodes to a canonical `1:N` ordered by coordinate (so identical
geometry yields an identical `mesh_crc` regardless of insertion order). If
`interior` (a per-triangle `Bool` vector) is given, only interior triangles are
kept — used by CDT/refinement to drop exterior/hole regions.
"""
function to_mesh(T::Triangulation; interior::Union{Nothing,Vector{Bool}}=nothing)
    tris = NTuple{3,Int32}[]
    @inbounds for t in eachindex(T.alive)
        T.alive[t] || continue
        _is_ghost_tri(T, t) && continue
        interior !== nothing && !interior[t] && continue
        push!(tris, (_vert(T,t,1), _vert(T,t,2), _vert(T,t,3)))
    end
    # canonical node ordering by coordinate
    usedset = Set{Int32}()
    for tr in tris, v in tr; push!(usedset, v); end
    used = sort(collect(usedset); by = v -> (T.x[v], T.y[v]))
    newid = Dict{Int32,Int32}()
    for (k, v) in enumerate(used); newid[v] = Int32(k); end
    coords = Matrix{Float64}(undef, 3, length(used))
    @inbounds for (k, v) in enumerate(used)
        coords[1,k]=T.x[v]; coords[2,k]=T.y[v]; coords[3,k]=0.0
    end
    triM = Matrix{Int32}(undef, 3, length(tris))
    @inbounds for (j, tr) in enumerate(tris)
        triM[1,j]=newid[tr[1]]; triM[2,j]=newid[tr[2]]; triM[3,j]=newid[tr[3]]
    end
    return Mesh(coords; tris=triM)
end

triangulate(xs::Vector{Float64}, ys::Vector{Float64}; rng_seed::Integer=1) =
    to_mesh(delaunay2d(xs, ys; rng_seed=rng_seed))

# ════════════════════════════════════════════════════════════════════════════════
# Verification helpers (CRC test suite)
# ════════════════════════════════════════════════════════════════════════════════

"""
    check_consistency(T) -> (ok::Bool, msg::String)

Structural invariants: real triangles are CCW; ghost triangles carry exactly one
`GHOST`; every neighbour link is mutual and matches on a shared *reversed* edge;
no dangling references to dead triangles.
"""
function check_consistency(T::Triangulation)
    @inbounds for t in eachindex(T.alive)
        T.alive[t] || continue
        gs = _ghost_slot(T, t)
        a = _vert(T,t,1); b = _vert(T,t,2); c = _vert(T,t,3)
        (a==b || b==c || a==c) && return (false, "degenerate triangle $t ($a,$b,$c)")
        if gs == 0
            orient2(_pt(T,a),_pt(T,b),_pt(T,c)) > 0 || return (false, "real triangle $t not CCW")
        end
        for k in 1:3
            s = _nbr(T, t, k)
            s == 0 && return (false, "triangle $t has an open edge (slot $k = 0)")
            T.alive[s] || return (false, "triangle $t neighbour $s dead")
            j = _nslot(T, s, t)
            j == 0 && return (false, "triangle $t→$s not mutual")
            u, w = _edge(T, t, k); u2, w2 = _edge(T, s, j)
            (u==w2 && w==u2) || return (false, "triangle $t/$s edge mismatch")
        end
    end
    return (true, "ok")
end

"""
    is_delaunay(T) -> (ok::Bool, nviol::Int)

The defining oracle: for every live *real* triangle, no other real vertex lies
strictly inside its circumcircle (exact `incircle_sos`). Exact ⇒ a true
independent verification of the empty-circumcircle property.
"""
function is_delaunay(T::Triangulation)
    nviol = 0
    @inbounds for t in eachindex(T.alive)
        T.alive[t] || continue
        _is_ghost_tri(T, t) && continue
        a = _vert(T,t,1); b = _vert(T,t,2); c = _vert(T,t,3)
        for v in 1:T.nreal
            (v==a || v==b || v==c) && continue
            if incircle_sos(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,v), a,b,c,Int32(v)) > 0
                nviol += 1
            end
        end
    end
    return (nviol == 0, nviol)
end

end # module Mesh2D
