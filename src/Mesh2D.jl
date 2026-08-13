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
using ..MeshTypes: Mesh, triangle_area

export Triangulation, triangulate, delaunay2d, dedup_points
export check_consistency, is_delaunay, to_mesh, ntriangles_live, insert_point!
export insert_segment!, constrained_delaunay, classify_interior, is_constrained_delaunay
export refine!

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
    vtri::Vector{Int32}         # per real vertex: one incident triangle (O(1) lookup hint)
    seg::Set{NTuple{2,Int32}}   # constrained edges (sorted vertex pairs) for CDT
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
        _touch_vtri!(T, Int32(t), a, b, c)
        return Int32(t)
    else
        length(T.alive)<typemax(Int32) ||
            throw(ArgumentError("Mesh2D: triangle storage exceeds the Int32 indexing limit"))
        push!(T.tv, Int32(a), Int32(b), Int32(c))
        push!(T.tn, Int32(0), Int32(0), Int32(0))
        push!(T.alive, true)
        t = Int32(length(T.alive))
        _touch_vtri!(T, t, a, b, c)
        return t
    end
end

# record `t` as an incident triangle for each of its real vertices (lookup hint)
@inline function _touch_vtri!(T::Triangulation, t::Int32, a, b, c)
    @inbounds begin
        (1 <= a <= T.nreal) && (T.vtri[a] = t)
        (1 <= b <= T.nreal) && (T.vtri[b] = t)
        (1 <= c <= T.nreal) && (T.vtri[c] = t)
    end
    return nothing
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

function _strictly_between(u, v, p)
    # The caller has already certified exact collinearity.  Comparing along any
    # coordinate on which the endpoints differ is exact for Float64 and avoids
    # overflowing the dot-product form on opposite-sign extreme coordinates.
    u[1]!=v[1] && return min(u[1],v[1])<p[1]<max(u[1],v[1])
    u[2]!=v[2] && return min(u[2],v[2])<p[2]<max(u[2],v[2])
    return false
end

# ════════════════════════════════════════════════════════════════════════════════
# Construction & initialization
# ════════════════════════════════════════════════════════════════════════════════

function Triangulation(xs::Vector{Float64}, ys::Vector{Float64})
    n = _validate_points(xs, ys, "Triangulation")
    return Triangulation(xs, ys, Int32[], Int32[], Bool[], Int32[], n, Int32(0),
                         zeros(Int32, n), Set{NTuple{2,Int32}}())
end

function _validate_points(xs::Vector{Float64}, ys::Vector{Float64}, caller::AbstractString)
    n = length(xs)
    length(ys) == n ||
        throw(ArgumentError("$caller: xs and ys must have the same length (got $n and $(length(ys)))"))
    n <= typemax(Int32) ||
        throw(ArgumentError("$caller: $n points exceed the Int32 indexing limit"))
    @inbounds for i in 1:n
        (isfinite(xs[i]) && isfinite(ys[i])) ||
            throw(ArgumentError("$caller: point $i has non-finite coordinates ($(xs[i]), $(ys[i]))"))
    end
    return n
end

function _seed64(seed::Integer, caller::AbstractString)
    0 <= seed <= typemax(UInt64) ||
        throw(ArgumentError("$caller: rng_seed must be in 0:$(typemax(UInt64)) (got $seed)"))
    return UInt64(seed)
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

function insert_point!(T::Triangulation, vid::Integer; constrained::Bool=false,
                      newtris::Union{Nothing,Vector{Int32}}=nothing)
    1 <= vid <= T.nreal ||
        throw(ArgumentError("insert_point!: vertex $vid is outside 1:$(T.nreal)"))
    @inbounds T.vtri[vid] == 0 ||
        throw(ArgumentError("insert_point!: vertex $vid is already present in the triangulation"))
    vid = Int32(vid)
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
            a, b = _edge(T, t, k)
            # a constrained edge is never crossed by the growing cavity (CDT-preserving)
            blocked = constrained && (_skey(a, b) in T.seg)
            inside = !blocked && nb != 0 && _in_circumcircle(T, nb, vid)
            if inside
                push!(incav, nb); push!(cavity, nb); push!(stack, nb)
            else
                push!(boundary, (a, b, nb))
            end
        end
    end
    for t in cavity; _killtri!(T, t); end
    _retriangulate_cavity!(T, boundary, Int32(vid), newtris)
    return nothing
end

function _retriangulate_cavity!(T::Triangulation, boundary, vid::Int32,
                                newtris::Union{Nothing,Vector{Int32}}=nothing)
    spoke = Dict{Int32, Tuple{Int32,Int32}}()   # boundary vertex w → (new tri, slot)
    anytri = Int32(0)
    for (a, b, nb) in boundary
        t = _newtri!(T, a, b, vid)               # [a, b, vid]; a or b may be GHOST
        anytri = t
        newtris !== nothing && push!(newtris, t)
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
# Edge flip (foundation of constrained insertion & Lawson legalization)
# ════════════════════════════════════════════════════════════════════════════════

# Flip the edge opposite local vertex `k` in triangle `t` (shared with its
# neighbour). The two triangles t=[p,q,r] (apex p) and t2 (apex s) forming the
# quad p→q→s→r are replaced in place by [p,q,s] and [p,s,r] (the other diagonal
# p–s). All four outer neighbours and both apex spokes are relinked. Returns the
# two new triangle ids. Neither triangle may be a ghost, and the quad must be
# strictly convex (checked by callers via orient2).
function _flip!(T::Triangulation, t::Integer, k::Integer)
    p = _vert(T, t, k)
    q, r = _edge(T, t, k)                     # CCW edge q→r, p to its left
    t2 = _nbr(T, t, k)
    j = _nslot(T, t2, t)
    s = _vert(T, t2, j)                       # apex of t2
    # outer neighbours by undirected edge
    Npq = _nbr(T, t, _nslot_by_edge(T, t, p, q))
    Nrp = _nbr(T, t, _nslot_by_edge(T, t, r, p))
    Nqs = _nbr(T, t2, _nslot_by_edge(T, t2, q, s))
    Nsr = _nbr(T, t2, _nslot_by_edge(T, t2, s, r))
    # rebuild in place: A reuses t = [p,q,s]; B reuses t2 = [p,s,r]
    A = Int32(t); B = Int32(t2)
    _setv!(T, A, 1, p); _setv!(T, A, 2, q); _setv!(T, A, 3, s)
    _setv!(T, B, 1, p); _setv!(T, B, 2, s); _setv!(T, B, 3, r)
    # A neighbours: opp p(slot1)=(q,s)→Nqs; opp q(slot2)=(s,p)→B; opp s(slot3)=(p,q)→Npq
    _setn!(T, A, 1, Nqs); _setn!(T, A, 2, B); _setn!(T, A, 3, Npq)
    # B neighbours: opp p(slot1)=(s,r)→Nsr; opp s(slot2)=(r,p)→Nrp; opp r(slot3)=(p,s)→A
    _setn!(T, B, 1, Nsr); _setn!(T, B, 2, Nrp); _setn!(T, B, 3, A)
    # relink outer neighbours' back-pointers (match by undirected edge)
    _relink!(T, Npq, A, p, q)
    _relink!(T, Nqs, A, q, s)
    _relink!(T, Nsr, B, s, r)
    _relink!(T, Nrp, B, r, p)
    _touch_vtri!(T, A, p, q, s); _touch_vtri!(T, B, p, s, r)
    T.last = A
    return A, B
end

@inline function _relink!(T::Triangulation, nb, newt, u, w)
    nb == 0 && return
    slot = _nslot_by_edge(T, nb, u, w)
    slot != 0 && _setn!(T, nb, slot, newt)
    return
end

# ════════════════════════════════════════════════════════════════════════════════
# Constrained Delaunay: segment insertion (flip-based, Sloan/Anglada)
# ════════════════════════════════════════════════════════════════════════════════

@inline _skey(a, b) = a < b ? (Int32(a), Int32(b)) : (Int32(b), Int32(a))
@inline _vslot(T::Triangulation, t, v) =
    (@inbounds T.tv[3*(t-1)+1]==v) ? 1 : ((@inbounds T.tv[3*(t-1)+2]==v) ? 2 :
     ((@inbounds T.tv[3*(t-1)+3]==v) ? 3 : 0))

# a live real triangle incident to real vertex vi (hint, with scan fallback)
function _find_incident(T::Triangulation, vi::Integer)
    h = @inbounds T.vtri[vi]
    (h != 0 && T.alive[h] && _vslot(T, h, vi) != 0) && return Int32(h)
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && _vslot(T, t, vi) != 0) && (T.vtri[vi] = Int32(t); return Int32(t))
    end
    error("Mesh2D: vertex $vi not incident to any live triangle")
end

# rotate CCW around vi: from triangle t (vi at slot s), the next incident triangle
# is across the edge (vi, w) where w is the CCW-previous vertex — that edge is
# opposite the CCW-next vertex u (slot s%3+1).
@inline function _rotate_ccw(T::Triangulation, t, s)
    return _nbr(T, t, s % 3 + 1)
end

"""Is `(vi,vj)` already an edge of the triangulation? Return `(t, k)` or `(0,0)`."""
function _edge_between(T::Triangulation, vi::Integer, vj::Integer)
    t0 = _find_incident(T, vi); t = t0
    while true
        s = _vslot(T, t, vi)
        u = _vert(T, t, s % 3 + 1)          # CCW-next
        if u == vj
            k = _nslot_by_edge(T, t, vi, vj)
            return (Int32(t), Int32(k))
        end
        nt = _rotate_ccw(T, t, s)
        (nt == 0 || nt == t0) && break
        t = nt
    end
    return (Int32(0), Int32(0))
end

# proper crossing of open segments (e1,e2) and (vi,vj); shared endpoints ⇒ false
function _crosses_segment(T::Triangulation, e1, e2, vi, vj)
    (e1==vi||e1==vj||e2==vi||e2==vj) && return false
    o1 = orient2(_pt(T,vi),_pt(T,vj),_pt(T,e1))
    o2 = orient2(_pt(T,vi),_pt(T,vj),_pt(T,e2))
    o3 = orient2(_pt(T,e1),_pt(T,e2),_pt(T,vi))
    o4 = orient2(_pt(T,e1),_pt(T,e2),_pt(T,vj))
    return (o1>0) != (o2>0) && (o3>0) != (o4>0) && o1!=0 && o2!=0 && o3!=0 && o4!=0
end

# entry into the segment walk: rotate around vi to the triangle whose opposite
# edge vi→vj crosses. Returns (:cross, t, u, w) with the entry triangle t and its
# opposite edge (u,w) — u left, w right of vi→vj — or (:vertex, w, t) if vj is
# collinear beyond an incident vertex w (split there).
function _segment_entry(T::Triangulation, vi::Integer, vj::Integer)
    t0 = _find_incident(T, vi); t = t0
    while true
        s = _vslot(T, t, vi)
        u = _vert(T, t, s % 3 + 1)          # CCW-next
        w = _vert(T, t, (s+1) % 3 + 1)      # CCW-prev
        if !_is_ghost_v(u) && !_is_ghost_v(w)
            ou = orient2(_pt(T,vi), _pt(T,u), _pt(T,vj))
            ow = orient2(_pt(T,vi), _pt(T,w), _pt(T,vj))
            (ou == 0 && _forward(T, vi, u, vj)) && return (:vertex, Int32(t), Int32(u), Int32(0))
            (ow == 0 && _forward(T, vi, w, vj)) && return (:vertex, Int32(t), Int32(w), Int32(0))
            (ou > 0 && ow < 0) && return (:cross, Int32(t), Int32(u), Int32(w))
        end
        nt = _rotate_ccw(T, t, s)
        (nt == 0 || nt == t0) && break
        t = nt
    end
    error("Mesh2D._segment_entry: no entry triangle for ($vi,$vj)")
end

# is w strictly between vi and vj on their common line (assumes collinear)?
@inline function _forward(T, vi, w, vj)
    pv=_pt(T,vi); pw=_pt(T,w); pj=_pt(T,vj)
    dx=pj[1]-pv[1]; dy=pj[2]-pv[2]
    s=(pw[1]-pv[1])*dx+(pw[2]-pv[2])*dy
    L=dx*dx+dy*dy
    return 0.0 < s < L
end

# triangles sharing undirected edge (a,b): returns (t1, k1, t2) with k1 the slot
# in t1 opposite (a,b); (0,0,0) if the edge doesn't exist.
function _edge_triangles(T::Triangulation, a, b)
    t0 = _find_incident(T, a); t = t0
    while true
        s = _vslot(T, t, a)
        u = _vert(T, t, s % 3 + 1)
        if u == b
            k = _nslot_by_edge(T, t, a, b)
            return (Int32(t), Int32(k), _nbr(T, t, k))
        end
        nt = _rotate_ccw(T, t, s)
        (nt == 0 || nt == t0) && break
        t = nt
    end
    return (Int32(0), Int32(0), Int32(0))
end

"""
    insert_segment!(T, vi, vj)

Force segment `(vi, vj)` to appear as an edge (constrained Delaunay), then mark
it constrained. Crossed edges are removed by flips; if a vertex lies on the
segment the segment is split there. Finally the affected non-constrained edges
are re-legalized to restore the constrained-Delaunay property.
"""
function insert_segment!(T::Triangulation, vi::Integer, vj::Integer)
    1 <= vi <= T.nreal ||
        throw(ArgumentError("insert_segment!: vertex $vi is outside 1:$(T.nreal)"))
    1 <= vj <= T.nreal ||
        throw(ArgumentError("insert_segment!: vertex $vj is outside 1:$(T.nreal)"))
    vi == vj && throw(ArgumentError("insert_segment!: a constraint must have two distinct endpoints"))
    (@inbounds T.vtri[vi] != 0) ||
        throw(ArgumentError("insert_segment!: vertex $vi has not been inserted"))
    (@inbounds T.vtri[vj] != 0) ||
        throw(ArgumentError("insert_segment!: vertex $vj has not been inserted"))
    vi = Int32(vi); vj = Int32(vj)
    et, _ = _edge_between(T, vi, vj)
    if et != 0
        push!(T.seg, _skey(vi, vj)); return
    end
    kind, t, u, w = _segment_entry(T, vi, vj)
    if kind === :vertex
        insert_segment!(T, vi, u); insert_segment!(T, u, vj); return
    end
    # kind === :cross : entry triangle t, opposite (crossed) edge = {u,w}. Assign
    # left/right by side of the segment vi→vj (the wedge order u,w is unrelated).
    crossed = Tuple{Int32,Int32}[]
    la, ra = orient2(_pt(T,vi), _pt(T,vj), _pt(T,u)) > 0 ? (u, w) : (w, u)
    push!(crossed, _skey(la, ra))
    guard = 0; maxstep = 8*length(T.alive) + 64
    while true
        guard += 1; guard > maxstep && error("Mesh2D: segment walk stalled")
        k = _nslot_by_edge(T, t, la, ra)
        nt = _nbr(T, t, k)
        _is_ghost_tri(T, nt) && error("Mesh2D: segment walk left the hull (bug)")
        c = _vert(T, nt, _nslot(T, nt, t))           # apex of the triangle ahead
        c == vj && break
        o = orient2(_pt(T,vi), _pt(T,vj), _pt(T,c))
        if o == 0
            insert_segment!(T, vi, c); insert_segment!(T, c, vj); return
        elseif o > 0                                  # c left → next edge (c, ra)
            la = Int32(c)
        else                                          # c right → next edge (la, c)
            ra = Int32(c)
        end
        push!(crossed, _skey(la, ra))
        t = nt
    end
    # flip crossed edges until the segment exists
    newedges = Tuple{Int32,Int32}[]
    # Amortized-O(1) FIFO. `popfirst!` shifts the full vector on every flip and made
    # long constraint walks quadratic in both copied bytes and wall time.
    front = reverse!(crossed)
    back = Tuple{Int32,Int32}[]
    guard = 0; maxflip = 200*length(crossed) + 1000
    while !isempty(front) || !isempty(back)
        guard += 1; guard > maxflip && error("Mesh2D: segment flip loop stalled")
        if isempty(front)
            reverse!(back)
            front, back = back, front
        end
        (a, b) = pop!(front)
        _skey(a, b) in T.seg &&
            throw(ArgumentError("Mesh2D: constraint ($vi,$vj) crosses existing constraint ($a,$b); PSLG segments must not cross (split them at the intersection first)"))
        t1, k1, t2 = _edge_triangles(T, a, b)
        t1 == 0 && continue                          # edge already gone
        (_is_ghost_tri(T, t1) || _is_ghost_tri(T, t2)) && continue
        p = _vert(T, t1, k1)                          # apex of t1
        x, y = _edge(T, t1, k1)                        # the shared edge, CCW in t1
        s = _vert(T, t2, _nslot(T, t2, t1))           # apex of t2
        # convex quad ⇔ x, y on opposite sides of line p→s (orientation-agnostic)
        ox = orient2(_pt(T,p),_pt(T,s),_pt(T,x)); oy = orient2(_pt(T,p),_pt(T,s),_pt(T,y))
        if (ox > 0) != (oy > 0) && ox != 0 && oy != 0
            _flip!(T, t1, k1)                          # (x,y) → (p,s)
            if _crosses_segment(T, p, s, vi, vj)
                push!(back, _skey(p, s))
            else
                push!(newedges, _skey(p, s))
            end
        else
            push!(back, (a, b))                        # non-convex: defer
        end
    end
    push!(T.seg, _skey(vi, vj))
    _legalize!(T, newedges)
    return nothing
end

# Constrained Lawson legalization: flip any non-constrained edge in `edges` (and
# edges exposed by those flips) that violates the empty-circumcircle test.
function _legalize!(T::Triangulation, edges::Vector{Tuple{Int32,Int32}})
    stack = copy(edges)
    guard = 0; maxg = 500*length(edges) + 2000
    while !isempty(stack)
        guard += 1
        guard > maxg &&
            error("Mesh2D: constrained-edge legalization did not converge after $maxg operations")
        (a, b) = pop!(stack)
        (_skey(a,b) in T.seg) && continue             # never flip a constraint
        t1, k1, t2 = _edge_triangles(T, a, b)
        t1 == 0 && continue
        (_is_ghost_tri(T, t1) || _is_ghost_tri(T, t2)) && continue
        p = _vert(T, t1, k1)
        x, y = _edge(T, t1, k1)                        # shared edge, CCW in t1 ⇒ (p,x,y) CCW
        s = _vert(T, t2, _nslot(T, t2, t1))
        ox = orient2(_pt(T,p),_pt(T,s),_pt(T,x)); oy = orient2(_pt(T,p),_pt(T,s),_pt(T,y))
        ((ox > 0) != (oy > 0) && ox != 0 && oy != 0) || continue   # convex only
        # s inside circumcircle of t1=(p,x,y) ⇒ edge illegal, flip it
        if incircle_sos(_pt(T,p),_pt(T,x),_pt(T,y),_pt(T,s), p,x,y,s) > 0
            _flip!(T, t1, k1)                          # (x,y)→(p,s)
            push!(stack, _skey(p,x)); push!(stack, _skey(x,s))
            push!(stack, _skey(s,y)); push!(stack, _skey(y,p))
        end
    end
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
    seed = _seed64(rng_seed, "delaunay2d")
    ux, uy, _ = dedup_points(xs, ys)
    T = Triangulation(ux, uy)
    placed = _init_triangulation!(T)
    isempty(placed) && return T                    # <3 non-collinear points → no triangles
    placedset = Set{Int32}(placed)
    order = Int32[i for i in 1:T.nreal if !(Int32(i) in placedset)]
    _shuffle_det!(order, seed)
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
    n = _validate_points(xs, ys, "dedup_points")
    seen = Dict{Tuple{Float64,Float64},Int32}()
    ux = Float64[]; uy = Float64[]; remap = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        # Dict uses `isequal`, which distinguishes -0.0 from +0.0 even though they
        # are the same geometric coordinate. Canonicalize zeros before keying.
        x = xs[i] == 0.0 ? 0.0 : xs[i]
        y = ys[i] == 0.0 ? 0.0 : ys[i]
        key = (x, y)
        id = get(seen, key, Int32(0))
        if id == 0
            push!(ux, x); push!(uy, y); id = Int32(length(ux)); seen[key] = id
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
function to_mesh(T::Triangulation; interior::Union{Nothing,AbstractVector{Bool}}=nothing)
    interior !== nothing && length(interior) != length(T.alive) &&
        throw(ArgumentError("to_mesh: interior mask length $(length(interior)) does not match triangle storage length $(length(T.alive))"))
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

"""
    constrained_delaunay(xs, ys, segments; rng_seed=1) -> Triangulation

Constrained Delaunay triangulation of the PSLG: points `(xs, ys)` plus `segments`
(a vector of `(i, j)` index pairs into the original points). Every segment is
forced to appear as an edge and marked constrained. Coincident points are merged
and segment indices relabelled consistently.
"""
function constrained_delaunay(xs::Vector{Float64}, ys::Vector{Float64},
                              segments::AbstractVector{<:Tuple{Integer,Integer}};
                              rng_seed::Integer=1)
    seed = _seed64(rng_seed, "constrained_delaunay")
    ux, uy, remap = dedup_points(xs, ys)
    n = length(xs)
    segseen=Set{NTuple{2,Int32}}()
    for (s, (i, j)) in enumerate(segments)
        1 <= i <= n ||
            throw(ArgumentError("constrained_delaunay: segment $s endpoint $i is outside 1:$n"))
        1 <= j <= n ||
            throw(ArgumentError("constrained_delaunay: segment $s endpoint $j is outside 1:$n"))
        i != j ||
            throw(ArgumentError("constrained_delaunay: segment $s has identical endpoints $i"))
        remap[i] != remap[j] ||
            throw(ArgumentError("constrained_delaunay: segment $s collapses after coincident-point deduplication"))
        key=_skey(remap[i],remap[j])
        key in segseen && throw(ArgumentError("constrained_delaunay: duplicate segment $s after point deduplication"))
        push!(segseen,key)
    end
    T = Triangulation(ux, uy)
    placed = _init_triangulation!(T)
    if isempty(placed)
        isempty(segments) && return T          # no non-collinear triple → empty (valid)
        throw(ArgumentError("constrained_delaunay: the points are degenerate (fewer than 3 " *
              "non-coincident, non-collinear points); a 2-D triangulation with constraints " *
              "does not exist"))
    end
    placedset = Set{Int32}(placed)
    order = Int32[i for i in 1:T.nreal if !(Int32(i) in placedset)]
    _shuffle_det!(order, seed)
    for v in order; insert_point!(T, v); end
    for (i, j) in segments
        vi = remap[i]; vj = remap[j]
        insert_segment!(T, vi, vj)
    end
    return T
end

"""
    classify_interior(T) -> Vector{Bool}

Per-triangle interior flag by a **parity** flood-fill from the ghost region:
crossing a constrained edge toggles inside/outside, other edges keep it. A real
triangle is interior iff its parity is odd. This correctly excludes nested holes
(inside the outer loop but inside a hole loop → even parity → exterior). Requires
the domain boundary to be closed constrained loops (as a meshing domain must be).
"""
function classify_interior(T::Triangulation)
    ntri = length(T.alive)
    if !isempty(T.seg)
        degree = zeros(UInt8, T.nreal)
        for (a, b) in T.seg
            @inbounds begin
                degree[a] < 2 ||
                    throw(ArgumentError("classify_interior: constrained boundary branches at vertex $a"))
                degree[b] < 2 ||
                    throw(ArgumentError("classify_interior: constrained boundary branches at vertex $b"))
                degree[a] += 1; degree[b] += 1
            end
        end
        @inbounds for v in eachindex(degree)
            degree[v] == 0 && continue
            degree[v] == 2 ||
                throw(ArgumentError("classify_interior: constrained boundary is open at vertex $v (degree $(degree[v]))"))
        end
    end
    visited = falses(ntri)
    par = falses(ntri)                 # false = exterior side, true = interior
    stack = Tuple{Int32,Bool}[]
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && _is_ghost_tri(T, t)) || continue
        visited[t] = true; par[t] = false; push!(stack, (Int32(t), false))
    end
    @inbounds while !isempty(stack)
        t, p = pop!(stack)
        for k in 1:3
            nb = _nbr(T, t, k)
            (nb == 0 || !T.alive[nb]) && continue
            a, b = _edge(T, t, k)
            np = (_skey(a, b) in T.seg) ? !p : p    # crossing a constraint toggles
            if visited[nb]
                par[nb] == np ||
                    throw(ArgumentError("classify_interior: inconsistent boundary parity across edge $(_skey(a,b))"))
                continue
            end
            visited[nb] = true; par[nb] = np; push!(stack, (nb, np))
        end
    end
    interior = falses(ntri)
    @inbounds for t in eachindex(T.alive)
        interior[t] = T.alive[t] && !_is_ghost_tri(T, t) && par[t]
    end
    return interior
end

# ════════════════════════════════════════════════════════════════════════════════
# Ruppert quality refinement
# ════════════════════════════════════════════════════════════════════════════════

# append a real vertex; returns its id (nreal grows; vtri extended)
function _add_vertex!(T::Triangulation, x::Float64, y::Float64)
    (isfinite(x)&&isfinite(y)) || throw(ArgumentError("refine!: attempted to add a non-finite vertex"))
    T.nreal<typemax(Int32) || throw(ArgumentError("refine!: node count exceeds the Int32 indexing limit"))
    push!(T.x, x); push!(T.y, y); T.nreal += 1
    push!(T.vtri, Int32(0))
    return Int32(T.nreal)
end

@inline _dist(a, b) = hypot(a[1]-b[1],a[2]-b[2])

@inline function _tri_area2(a, b, c)
    return triangle_area((a[1],a[2],0.0),(b[1],b[2],0.0),(c[1],c[2],0.0))
end

function _circumcenter2(a, b, c)
    ax,ay=a;bx=b[1]-ax;by=b[2]-ay;cx=c[1]-ax;cy=c[2]-ay
    scale=max(abs(bx),abs(by),abs(cx),abs(cy))
    (isfinite(scale)&&scale>0)||return (NaN,NaN)
    bx/=scale;by/=scale;cx/=scale;cy/=scale
    d=2.0*(bx*cy-by*cx);d==0&&return (NaN,NaN)
    b2=bx*bx+by*by;c2=cx*cx+cy*cy
    ux=ax+scale*(cy*b2-by*c2)/d
    uy=ay+scale*(bx*c2-cx*b2)/d
    return (ux,uy)
end

# radius-edge ratio = circumradius / shortest edge (≈ 1/(2 sin θ_min))
function _radius_edge2(a, b, c)
    lab=_dist(a,b); lbc=_dist(b,c); lca=_dist(c,a)
    cc=_circumcenter2(a,b,c);R=_dist(a,cc)
    isfinite(R)||return Inf
    return R / min(lab,lbc,lca)
end

@inline function _needs_refine(T, t, B, maxarea, sizefn)
    a=_pt(T,_vert(T,t,1)); b=_pt(T,_vert(T,t,2)); c=_pt(T,_vert(T,t,3))
    (_tri_area2(a,b,c) > maxarea || _radius_edge2(a,b,c) > B) && return true
    if sizefn !== nothing
        cx = a[1]/3+b[1]/3+c[1]/3; cy = a[2]/3+b[2]/3+c[2]/3
        rawh = sizefn(cx, cy)
        rawh isa Real ||
            throw(ArgumentError("refine!: size callback must return a real size (got $(typeof(rawh))) at ($cx,$cy)"))
        h=_mesh2_float(rawh,"size callback result")
        (isfinite(h) && h > 0) ||
            throw(ArgumentError("refine!: size callback returned invalid size $rawh at ($cx,$cy)"))
        emax = max(_dist(a,b), _dist(b,c), _dist(c,a))
        emax > h && return true
    end
    return false
end

# is subsegment (a,b) encroached by point p? (p strictly inside the diametral disk)
function _encroaches(pa,pb,p)
    ax=pa[1]-p[1];ay=pa[2]-p[2];bx=pb[1]-p[1];by=pb[2]-p[2]
    scale=max(abs(ax),abs(ay),abs(bx),abs(by))
    if isfinite(scale)&&scale>0
        ax/=scale;ay/=scale;bx/=scale;by/=scale
        dot=ax*bx+ay*by;perm=abs(ax*bx)+abs(ay*by)
        abs(dot)>16eps(Float64)*perm && return dot<0
    end
    R=Rational{BigInt}
    axb=R(pa[1])-R(p[1]);ayb=R(pa[2])-R(p[2])
    bxb=R(pb[1])-R(p[1]);byb=R(pb[2])-R(p[2])
    return axb*bxb+ayb*byb<0
end

# encroached by an adjacent apex? (in a CDT this suffices to detect any encroachment)
function _subseg_encroached(T::Triangulation, a, b)
    t1, k1, t2 = _edge_triangles(T, a, b)
    t1 == 0 && return false
    pa = _pt(T,a); pb = _pt(T,b)
    p = _vert(T, t1, k1)
    !_is_ghost_v(p) && _encroaches(pa, pb, _pt(T,p)) && return true
    j = _nslot(T, t2, t1); s = _vert(T, t2, j)
    !_is_ghost_v(s) && _encroaches(pa, pb, _pt(T,s)) && return true
    return false
end

@inline function _setflag!(v::Vector{Bool}, i::Integer, val::Bool)
    while length(v) < i; push!(v, false); end
    @inbounds v[i] = val
    return nothing
end

# split subsegment (a,b) at its midpoint; update constraints and interior flags
function _split_subsegment!(T::Triangulation, a::Int32, b::Int32, interior::Vector{Bool},pointids)
    pa = _pt(T,a); pb = _pt(T,b)
    mx=pa[1]/2+pb[1]/2;my=pa[2]/2+pb[2]/2
    mp=(mx==0 ? 0.0 : mx,my==0 ? 0.0 : my)
    (isfinite(mp[1])&&isfinite(mp[2])&&mp!=pa&&mp!=pb) ||
        throw(ErrorException("refine!: constrained-segment midpoint is below Float64 coordinate resolution"))
    existing=get(pointids,mp,Int32(0))
    existing==0 || throw(ErrorException("refine!: constrained-segment midpoint already exists as vertex $existing"))
    mid = _add_vertex!(T,mp[1],mp[2])
    pointids[mp]=mid
    delete!(T.seg, _skey(a,b))
    push!(T.seg, _skey(a,mid)); push!(T.seg, _skey(mid,b))
    insert_point!(T, mid; constrained=true)
    # constraints changed → recompute interior (segment splits are comparatively rare)
    ni = classify_interior(T)
    resize!(interior, length(ni)); @inbounds for i in eachindex(ni); interior[i]=ni[i]; end
    return mid
end

"""
    refine!(T; min_angle_deg=25.0, max_area=Inf, maxsteps=300_000) -> Vector{Bool}

Ruppert refinement of the constrained Delaunay triangulation `T`: split
encroached subsegments and insert circumcenters of skinny interior triangles
(radius-edge ratio > 1/(2 sin θ), or area > `max_area`) until the quality bound
holds. Constraints are preserved (constrained point insertion). Returns the final
per-triangle interior mask. `min_angle_deg ≲ 20.7°` guarantees termination.
If `maxsteps` is exhausted while any criterion remains unmet, refinement throws an
explicit blocker rather than returning a mesh that violates the requested bounds.
"""
function refine!(T::Triangulation; min_angle_deg::Real=25.0, max_area::Real=Inf,
                 size=nothing, maxsteps::Integer=300_000)
    angle = _mesh2_float(min_angle_deg,"min_angle_deg")
    (isfinite(angle) && 0 <= angle < 60) ||
        throw(ArgumentError("refine!: min_angle_deg must be finite and in [0, 60) (got $min_angle_deg)"))
    area = _mesh2_float(max_area,"max_area")
    (!isnan(area) && area > 0) ||
        throw(ArgumentError("refine!: max_area must be positive or Inf (got $max_area)"))
    maxsteps > 0 || throw(ArgumentError("refine!: maxsteps must be positive (got $maxsteps)"))
    maxsteps<=typemax(Int) || throw(ArgumentError("refine!: maxsteps exceeds the platform Int limit"))
    nsteps=Int(maxsteps)
    B = angle == 0 ? Inf : 1.0 / (2.0*sind(angle))
    interior = Vector{Bool}(collect(classify_interior(T)))
    pointids=Dict{NTuple{2,Float64},Int32}()
    sizehint!(pointids,T.nreal)
    @inbounds for i in 1:T.nreal
        key=(T.x[i]==0 ? 0.0 : T.x[i],T.y[i]==0 ? 0.0 : T.y[i])
        haskey(pointids,key) && throw(ArgumentError("refine!: triangulation contains duplicate point coordinates"))
        pointids[key]=Int32(i)
    end
    for _ in 1:nsteps
        # (1) split an encroached subsegment, if any
        enc = _find_encroached(T)
        if enc !== nothing
            _split_subsegment!(T, enc[1], enc[2], interior,pointids)
            continue
        end
        # (2) pick a triangle that violates the angle/area/size criteria
        t = _find_bad(T, interior, B, area, size)
        t == 0 && return interior
        a=_vert(T,t,1); b=_vert(T,t,2); c=_vert(T,t,3)
        cc = _circumcenter2(_pt(T,a),_pt(T,b),_pt(T,c))
        # (3) if the circumcenter encroaches subsegments, split those instead
        sp = _encroached_by_point(T, cc)
        if sp !== nothing
            _split_subsegment!(T, sp[1], sp[2], interior,pointids)
            continue
        end
        # (4) insert the circumcenter (constrained); skip if it falls outside the domain
        pa=_pt(T,a);pb=_pt(T,b);pc=_pt(T,c)
        fallback=(pa[1]/3+pb[1]/3+pc[1]/3,pa[2]/3+pb[2]/3+pc[2]/3)
        _insert_steiner!(T, cc, interior,pointids;fallback=fallback)
    end
    enc = _find_encroached(T)
    bad = _find_bad(T, interior, B, area, size)
    enc === nothing && bad == 0 && return interior
    remaining = enc === nothing ? "a triangle still violates a quality/area/size bound" :
                                  "constrained segment $(enc) remains encroached"
    throw(ErrorException("refine!: reached maxsteps=$maxsteps before convergence; $remaining"))
end

# Deterministic iteration order over the constraint set (a Set has none), so the
# refinement sequence — hence the output mesh — is reproducible across runs.
_sorted_segs(T::Triangulation) = sort!(collect(T.seg))

function _find_encroached(T::Triangulation)
    for (a, b) in _sorted_segs(T)
        _subseg_encroached(T, a, b) && return (a, b)
    end
    return nothing
end

function _find_bad(T::Triangulation, interior::Vector{Bool}, B, maxarea, sizefn)
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && t <= length(interior) && interior[t] && !_is_ghost_tri(T,t)) || continue
        _needs_refine(T, t, B, maxarea, sizefn) && return Int32(t)
    end
    return Int32(0)
end

# first subsegment (deterministic order) encroached by coordinate p, or nothing
function _encroached_by_point(T::Triangulation, p)
    for (a, b) in _sorted_segs(T)
        _encroaches(_pt(T,a), _pt(T,b), p) && return (a, b)
    end
    return nothing
end

# insert Steiner point p (a coordinate) if it lands in an interior triangle.
# The vertex is added first (so locate's SoS/ghost tests can read its coordinates),
# then undone if the point falls outside the domain — locate never mutates the
# triangulation, so the trailing vertex is safe to pop.
function _insert_steiner!(T::Triangulation,p,interior::Vector{Bool},pointids;fallback=nothing)
    if !(isfinite(p[1])&&isfinite(p[2]))
        fallback===nothing &&
            throw(ErrorException("refine!: triangle circumcenter is not representable as a finite Float64 point"))
        return _insert_steiner!(T,fallback,interior,pointids)
    end
    key=(p[1]==0 ? 0.0 : p[1],p[2]==0 ? 0.0 : p[2])
    existing=get(pointids,key,Int32(0))
    if existing!=0
        fallback===nothing &&
            throw(ErrorException("refine!: required Steiner point coincides with existing vertex $existing at Float64 resolution"))
        return _insert_steiner!(T,fallback,interior,pointids)
    end
    vid = _add_vertex!(T, p[1], p[2])
    pointids[key]=vid
    tc = locate(T, p[1], p[2], vid)
    if _is_ghost_tri(T, tc) || tc > length(interior) || !interior[tc]
        pop!(T.x); pop!(T.y); pop!(T.vtri); T.nreal -= 1;delete!(pointids,key)     # undo
        fallback===nothing &&
            throw(ErrorException("refine!: no representable interior Steiner point for the selected triangle"))
        return _insert_steiner!(T,fallback,interior,pointids)
    end
    newt = Int32[]
    insert_point!(T, vid; constrained=true, newtris=newt)
    for nt in newt; _setflag!(interior, nt, true); end
    return true
end

function _mesh2_float(x::Real,name::AbstractString)
    try
        return Float64(x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("refine!: $name must be representable as Float64"))
    end
end

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

"""
    is_constrained_delaunay(T) -> (ok, n_missing_constraints, n_local_violations)

Constrained-Delaunay oracle: every constrained segment is present as an edge, and
every **non-constrained** interior edge is locally Delaunay (the opposite vertex
of its neighbour is not strictly inside the triangle's circumcircle).
"""
function is_constrained_delaunay(T::Triangulation)
    missing = 0
    for (a, b) in T.seg
        _edge_between(T, a, b)[1] == 0 && (missing += 1)
    end
    viol = 0
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tri(T, t)) || continue
        for k in 1:3
            nb = _nbr(T, t, k)
            (nb == 0 || _is_ghost_tri(T, nb)) && continue
            a, b = _edge(T, t, k)
            _skey(a, b) in T.seg && continue
            s = _vert(T, nb, _nslot(T, nb, t))
            p = _vert(T, t, k)                        # t = (p,a,b) CCW
            if incircle_sos(_pt(T,p),_pt(T,a),_pt(T,b),_pt(T,s), p,a,b,s) > 0
                viol += 1
            end
        end
    end
    return (missing == 0 && viol == 0, missing, viol)
end

end # module Mesh2D
