"""
    Mesh3D

Stage-3 three-dimensional meshing (PLAN.md §3 "Mesh3D") — the robustness
milestone. This module owns the 3-D Delaunay kernel, constrained
tetrahedralization / **boundary recovery** (recovering a given closed surface
triangulation as facets of the tet mesh), and sliver handling. It is exactly the
layer where gmsh 4.13/4.15 fail on the enclosure coax junction.

Foundation, mirroring the robust 2-D kernel:

* Incremental Bowyer–Watson with a single `GHOST` vertex at infinity closing the
  convex hull (ghost tetrahedra), not a finite bounding simplex.
* All decisions use the exact `orient3_sos` / `insphere_sos` predicates.
* A ghost tet's in-sphere test degenerates to an orientation test on its real
  hull face; the coplanar case (pervasive for box geometry) is resolved by an
  in-plane in-circle test — the 3-D analogue of the collinear-hull rule that made
  the 2-D grid triangulation exact.

Index-compact storage: flat `Int32` `tv` (4 vertices/tet) and `tn` (4
neighbours/tet, one opposite each vertex) plus a free list.
"""
module Mesh3D

using ..Predicates: orient3, orient3_sos, insphere_sos, incircle_sos, orient2
using ..MeshTypes: Mesh, tet_dihedral_extrema

export Triangulation3, delaunay3d, tetrahedralize, tetrahedralize_multi, tets_per_region
export check_consistency3, is_delaunay3, to_mesh3, ntets_live, present_faces
export flip23!, flip32!, tets_around_edge, optimize_flips!

const GHOST3 = Int32(0)

mutable struct Triangulation3
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    tv::Vector{Int32}          # 4 per tet
    tn::Vector{Int32}          # 4 per tet: neighbour opposite local vertex k
    alive::Vector{Bool}
    freelist::Vector{Int32}
    nreal::Int
    last::Int32
    vtet::Vector{Int32}        # per real vertex: an incident tet (hint)
end

@inline _pt(T::Triangulation3, i) = @inbounds (T.x[i], T.y[i], T.z[i])
@inline _vert(T::Triangulation3, t, k) = @inbounds T.tv[4*(t-1)+k]
@inline _setv!(T::Triangulation3, t, k, v) = @inbounds (T.tv[4*(t-1)+k] = Int32(v))
@inline _nbr(T::Triangulation3, t, k) = @inbounds T.tn[4*(t-1)+k]
@inline _setn!(T::Triangulation3, t, k, v) = @inbounds (T.tn[4*(t-1)+k] = Int32(v))
@inline _is_ghost_v(v) = v == GHOST3

@inline function _ghost_slot(T::Triangulation3, t)
    b = 4*(t-1)
    @inbounds (T.tv[b+1]==GHOST3) && return 1
    @inbounds (T.tv[b+2]==GHOST3) && return 2
    @inbounds (T.tv[b+3]==GHOST3) && return 3
    @inbounds (T.tv[b+4]==GHOST3) && return 4
    return 0
end
@inline _is_ghost_tet(T::Triangulation3, t) = _ghost_slot(T, t) != 0
ntets_live(T::Triangulation3) = count(T.alive)

# The three vertices of the face opposite local vertex k, ordered so that
# orient3(face..., v_k) < 0 (the opposite vertex is on the face's negative side).
# Consequently a point p is *outside* tet t across face k iff orient3(face..., p) > 0.
@inline function _face(T::Triangulation3, t, k)
    b = 4*(t-1)
    @inbounds v1=T.tv[b+1]; @inbounds v2=T.tv[b+2]; @inbounds v3=T.tv[b+3]; @inbounds v4=T.tv[b+4]
    if k == 1;     return (v2, v3, v4)
    elseif k == 2; return (v1, v4, v3)
    elseif k == 3; return (v1, v2, v4)
    else;          return (v1, v3, v2)
    end
end

@inline function _newtet!(T::Triangulation3, a, b, c, d)
    if !isempty(T.freelist)
        t = pop!(T.freelist); base = 4*(t-1)
        @inbounds T.tv[base+1]=a; T.tv[base+2]=b; T.tv[base+3]=c; T.tv[base+4]=d
        @inbounds T.tn[base+1]=0; T.tn[base+2]=0; T.tn[base+3]=0; T.tn[base+4]=0
        @inbounds T.alive[t]=true
        _touch_vtet!(T, Int32(t), a, b, c, d)
        return Int32(t)
    else
        push!(T.tv, Int32(a),Int32(b),Int32(c),Int32(d))
        push!(T.tn, Int32(0),Int32(0),Int32(0),Int32(0))
        push!(T.alive, true)
        t = Int32(length(T.alive)); _touch_vtet!(T, t, a, b, c, d)
        return t
    end
end

@inline function _touch_vtet!(T::Triangulation3, t, a, b, c, d)
    @inbounds begin
        (1<=a<=T.nreal) && (T.vtet[a]=t); (1<=b<=T.nreal) && (T.vtet[b]=t)
        (1<=c<=T.nreal) && (T.vtet[c]=t); (1<=d<=T.nreal) && (T.vtet[d]=t)
    end
end

@inline function _killtet!(T::Triangulation3, t)
    @inbounds T.alive[t]=false; push!(T.freelist, Int32(t)); nothing
end

@inline function _nslot(T::Triangulation3, t, s)
    b = 4*(t-1)
    @inbounds (T.tn[b+1]==s) && return 1
    @inbounds (T.tn[b+2]==s) && return 2
    @inbounds (T.tn[b+3]==s) && return 3
    @inbounds (T.tn[b+4]==s) && return 4
    return 0
end

@inline _sort3t(a,b,c) = begin
    a,b = a<=b ? (a,b) : (b,a); b,c = b<=c ? (b,c) : (c,b); a,b = a<=b ? (a,b) : (b,a); (a,b,c)
end
# slot in tet t whose opposite face has the vertex set {u,v,w}
@inline function _nslot_by_face(T::Triangulation3, t, u, v, w)
    key = _sort3t(u,v,w)
    @inbounds for k in 1:4
        f = _face(T, t, k)
        _sort3t(f[1],f[2],f[3]) == key && return k
    end
    return 0
end

# ── in-sphere (ghost-aware, coplanar-aware) ─────────────────────────────────────
# Is real vertex vid inside tet t's circumsphere?
@inline function _in_sphere(T::Triangulation3, t, vid)
    gs = _ghost_slot(T, t)
    if gs == 0
        a=_vert(T,t,1); b=_vert(T,t,2); c=_vert(T,t,3); d=_vert(T,t,4)
        return insphere_sos(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,d),_pt(T,vid), a,b,c,d,Int32(vid)) > 0
    else
        # ghost tet: real face = the 3 non-GHOST vertices (= _face opposite GHOST),
        # oriented so the SOLID interior is on its POSITIVE side. A point outside
        # the hull is therefore on the NEGATIVE side ⇒ in the ghost's (half-space)
        # circumsphere. Coplanar ⇒ resolve with the in-plane in-circle.
        f = _face(T, t, gs)
        u,v,w = f[1],f[2],f[3]
        o = orient3(_pt(T,u), _pt(T,v), _pt(T,w), _pt(T,vid))
        o < 0 && return true            # strictly outside the hull face → in ghost sphere
        o > 0 && return false
        return _in_face_circle(T, u, v, w, vid)   # coplanar: in-plane in-circle
    end
end

# p coplanar with (u,v,w): is it inside their circumcircle (in the shared plane)?
function _in_face_circle(T::Triangulation3, u, v, w, p)
    pu=_pt(T,u); pv=_pt(T,v); pw=_pt(T,w)
    # in-plane orthonormal axes
    e1 = _subn(pv, pu); e2 = _subn(pw, pu)
    n = _cross(e1, e2)
    nl = sqrt(_dot(n,n)); nl == 0 && return false     # degenerate face
    ax = _unitn(e1)
    ay = _cross((n[1]/nl, n[2]/nl, n[3]/nl), ax)
    pr(q) = (_dot(_subn(q,pu), ax), _dot(_subn(q,pu), ay))
    a2=pr(pu); b2=pr(pv); c2=pr(pw); p2=pr(_pt(T,p))
    ic = incircle_sos(a2,b2,c2,p2, u,v,w,Int32(p))
    o2 = orient2(a2,b2,c2)                            # face triangle orientation in-plane
    # inside ⇔ in-circle sign matches the CCW/CW orientation of (u,v,w)
    return (ic > 0) == (o2 > 0)
end

@inline _subn(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _dot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline function _unitn(a); l=sqrt(_dot(a,a)); (a[1]/l,a[2]/l,a[3]/l); end

# ── construction ────────────────────────────────────────────────────────────────
function Triangulation3(xs::Vector{Float64}, ys::Vector{Float64}, zs::Vector{Float64})
    n = length(xs); @assert length(ys)==n && length(zs)==n
    return Triangulation3(xs, ys, zs, Int32[], Int32[], Bool[], Int32[], n, Int32(0), zeros(Int32,n))
end

# first non-coplanar 4 real vertices → 1 real tet + 4 ghost tets; returns placed set
function _init3!(T::Triangulation3)
    n = T.nreal; n < 4 && return Int32[]
    a = Int32(1); b = Int32(0)
    @inbounds for i in 2:n
        (T.x[i],T.y[i],T.z[i]) != (T.x[a],T.y[a],T.z[a]) && (b=Int32(i); break)
    end
    b == 0 && return Int32[]
    c = Int32(0)
    @inbounds for i in 2:n
        i==b && continue
        # not collinear with a,b
        if _tri_area_sq(_pt(T,a),_pt(T,b),_pt(T,i)) > 0; c=Int32(i); break; end
    end
    c == 0 && return Int32[]
    d = Int32(0)
    @inbounds for i in 2:n
        (i==b||i==c) && continue
        if orient3(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,i)) != 0; d=Int32(i); break; end
    end
    d == 0 && return Int32[]
    # positive orientation
    if orient3(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,d)) < 0
        c, d = d, c
    end
    R = _newtet!(T, a, b, c, d)                       # positive real tet
    # four ghost tets, one per face; ghost stored [f1,f2,f3,GHOST] with GHOST slot 4
    _make_hull!(T, R)
    T.last = R
    return Int32[a,b,c,d]
end

@inline function _tri_area_sq(a,b,c)
    cr = _cross(_subn(b,a), _subn(c,a)); _dot(cr,cr)
end

# create ghost tets on all current boundary faces (faces whose neighbour is 0) and
# link them (to the real neighbour across the real face, and to each other across
# shared ghost edges). Used at init and whenever new hull faces appear.
function _make_hull!(T::Triangulation3, seed::Int32)
    # For the seed tet, each of its 4 faces is a hull face (neighbour 0).
    spoke = Dict{NTuple{2,Int32}, Tuple{Int32,Int32}}()   # ghost edge → (ghost tet, slot)
    @inbounds for k in 1:4
        f = _face(T, seed, k)                # (f1,f2,f3) with orient3(f,vk)<0
        # store ghost as [f1,f2,f3,GHOST] ⇒ real face _face(g,4)=(f1,f3,f2) = the
        # reverse of the seed's face (satisfies the reversed-shared-face invariant).
        g = _newtet!(T, f[1], f[2], f[3], GHOST3)
        _setn!(T, seed, k, g)
        _setn!(T, g, 4, seed)                # GHOST at slot 4 ⇒ real face opposite = slot 4
        _reg_ghost_spokes!(T, spoke, g)
    end
    isempty(spoke) || error("Mesh3D._make_hull: ghost shell did not close")
    return nothing
end

# For a ghost tet [p,q,r,GHOST] register its 3 faces that contain GHOST (the spokes).
function _reg_ghost_spokes!(T, spoke, g::Int32)
    @inbounds p=_vert(T,g,1); q=_vert(T,g,2); r=_vert(T,g,3)
    # faces containing GHOST are opposite p (slot1), q (slot2), r (slot3);
    # each shares the "ghost edge" = the two real vertices of that face.
    _link_ghost_spoke!(T, spoke, _edgekey(q, r), g, 1)   # face opp p = (q,r,GHOST)
    _link_ghost_spoke!(T, spoke, _edgekey(r, p), g, 2)   # face opp q = (r,p,GHOST)? see _face
    _link_ghost_spoke!(T, spoke, _edgekey(p, q), g, 3)   # face opp r
    return nothing
end

@inline _edgekey(a,b) = a<=b ? (Int32(a),Int32(b)) : (Int32(b),Int32(a))

@inline function _link_ghost_spoke!(T, spoke, key, g::Int32, slot::Int)
    if haskey(spoke, key)
        (g2, s2) = spoke[key]
        _setn!(T, g, slot, g2); _setn!(T, g2, Int(s2), g)
        delete!(spoke, key)
    else
        spoke[key] = (g, Int32(slot))
    end
end

# ── point location (stochastic walk) ────────────────────────────────────────────
function locate3(T::Triangulation3, px, py, pz, vid; start::Integer=T.last)
    t = start
    (t==0 || !T.alive[t]) && (t=_first_alive3(T))
    gs = _ghost_slot(T, t); gs != 0 && (t = _nbr(T, t, gs))
    prev = Int32(0)
    rs = UInt64(vid)*0x9E3779B97F4A7C15 + 0xD1B54A32D192ED03
    guard = 0; maxstep = 16*length(T.alive) + 64
    @inbounds while true
        guard += 1
        # a stochastic walk terminates for a Delaunay triangulation, but severely
        # thin tets (perturbed regular grids) can make it wander; fall back to an
        # exhaustive, guaranteed-correct scan rather than fail.
        guard > maxstep && return _locate_scan(T, px, py, pz, vid)
        rs ⊻= rs<<13; rs ⊻= rs>>7; rs ⊻= rs<<17
        r = Int(rs % UInt64(4))
        moved = false
        for i in 0:3
            k = (r+i)%4 + 1
            nb = _nbr(T, t, k)
            (nb == prev || nb == 0) && continue
            if _is_ghost_tet(T, nb)
                _in_sphere(T, nb, vid) && return nb
                continue
            end
            f = _face(T, t, k)
            if orient3_sos(_pt(T,f[1]),_pt(T,f[2]),_pt(T,f[3]),(px,py,pz), f[1],f[2],f[3],Int32(vid)) > 0
                prev = Int32(t); t = nb; moved = true; break
            end
        end
        moved || return t
    end
end

function _first_alive3(T::Triangulation3)
    @inbounds for t in eachindex(T.alive); T.alive[t] && return Int32(t); end
    error("Mesh3D: no live tets")
end

# Exhaustive location fallback (guaranteed correct): return a real tet containing
# p (all faces test interior), else a ghost tet whose half-space circumsphere
# contains p (p outside the hull). Uses SoS to break the on-face/on-hull ties.
function _locate_scan(T::Triangulation3, px, py, pz, vid)
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T, t)) || continue
        inside = true
        for k in 1:4
            f = _face(T, t, k)
            if orient3_sos(_pt(T,f[1]),_pt(T,f[2]),_pt(T,f[3]),(px,py,pz), f[1],f[2],f[3],Int32(vid)) > 0
                inside = false; break
            end
        end
        inside && return Int32(t)
    end
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && _is_ghost_tet(T, t)) || continue
        _in_sphere(T, t, vid) && return Int32(t)
    end
    error("Mesh3D._locate_scan: point not located (corrupt triangulation?)")
end

# ── Bowyer–Watson insertion ─────────────────────────────────────────────────────
function insert_point3!(T::Triangulation3, vid::Integer; newtets::Union{Nothing,Vector{Int32}}=nothing)
    px,py,pz = _pt(T, vid)
    t0 = locate3(T, px, py, pz, vid)
    cavity = Int32[t0]; incav = Set{Int32}(); push!(incav, t0)
    boundary = Tuple{Int32,Int32,Int32,Int32}[]   # (f1,f2,f3, outside neighbour)
    stack = Int32[t0]
    @inbounds while !isempty(stack)
        t = pop!(stack)
        for k in 1:4
            nb = _nbr(T, t, k)
            (nb != 0 && (nb in incav)) && continue
            inside = nb != 0 && _in_sphere(T, nb, vid)
            if inside
                push!(incav, nb); push!(cavity, nb); push!(stack, nb)
            else
                f = _face(T, t, k)
                push!(boundary, (f[1], f[2], f[3], nb))
            end
        end
    end
    for t in cavity; _killtet!(T, t); end
    _retriangulate3!(T, boundary, Int32(vid), newtets)
    return nothing
end

function _retriangulate3!(T::Triangulation3, boundary, vid::Int32, newtets)
    # One new tet per cavity boundary face, connecting it to vid. A real face
    # (all real) → a real tet [a,b,c,vid] made positive (orient3>0). A ghost face
    # (contains GHOST) → a new ghost tet [u,v,vid,GHOST] (a new hull facet). Every
    # face containing vid is a "spoke" shared with a sibling new tet; siblings are
    # matched by their non-vid edge (which may include GHOST). Uniform handling ⇒
    # correct hull extension when vid lands outside the current hull.
    spoke = Dict{NTuple{2,Int32}, Tuple{Int32,Int32}}()
    anyt = Int32(0)
    for (f1,f2,f3,nb) in boundary
        hasG = _is_ghost_v(f1) || _is_ghost_v(f2) || _is_ghost_v(f3)
        if !hasG
            # store [f1,f3,f2,vid] so _face(t,4) (opposite vid) == the boundary face
            # (f1,f2,f3) — i.e. REVERSED w.r.t. the outside neighbour's copy, as the
            # shared-face invariant requires. Automatically positive: vid is inside
            # the star-shaped cavity ⇒ orient3(f1,f2,f3,vid) < 0 ⇒ orient3(f1,f3,f2,vid) > 0.
            t = _newtet!(T, f1, f3, f2, vid)     # vid at slot 4
            vslot = 4
        else
            u, v = _reals_cyclic(f1, f2, f3)     # 2 reals in face-cyclic order
            t = _newtet!(T, u, v, vid, GHOST3)   # vid at slot 3, GHOST slot 4
            vslot = 3
        end
        anyt = t
        newtets !== nothing && push!(newtets, t)
        _setn!(T, t, vslot, nb)                  # face opposite vid → outside nb
        if nb != 0
            j = _nslot_by_face(T, nb, f1, f2, f3)
            j != 0 && _setn!(T, nb, j, t)
        end
        # register every face that CONTAINS vid (all slots but vslot) by its non-vid edge
        for k in 1:4
            k == vslot && continue
            fk = _face(T, t, k)
            e1, e2 = _non(fk, vid)
            _link_ghost_spoke!(T, spoke, _edgekey(e1, e2), t, k)
        end
    end
    isempty(spoke) || error("Mesh3D: cavity boundary not a topological sphere (unmatched face)")
    anyt != 0 && (T.last = anyt)
    return nothing
end

# the two non-GHOST vertices of a face, in the face's cyclic order
@inline function _reals_cyclic(f1, f2, f3)
    _is_ghost_v(f1) && return (f2, f3)
    _is_ghost_v(f2) && return (f3, f1)
    return (f1, f2)
end

# the two vertices of a 3-tuple face that are not `x`
@inline function _non(f, x)
    f[1] == x ? (f[2], f[3]) : (f[2] == x ? (f[1], f[3]) : (f[1], f[2]))
end

# ── driver + export ─────────────────────────────────────────────────────────────
"""
    delaunay3d(xs, ys, zs; rng_seed=1) -> Triangulation3

3-D Delaunay tetrahedralization of the (deduplicated) points.
"""
function delaunay3d(xs::Vector{Float64}, ys::Vector{Float64}, zs::Vector{Float64};
                    rng_seed::Integer=1, perturb::Bool=true)
    ux,uy,uz,_ = _dedup3(xs,ys,zs)
    perturb && _perturb3!(ux,uy,uz)
    T = Triangulation3(ux,uy,uz)
    placed = _init3!(T)
    isempty(placed) && return T
    ps = Set{Int32}(placed)
    order = Int32[i for i in 1:T.nreal if !(Int32(i) in ps)]
    _shuffle3!(order, UInt64(rng_seed))
    for v in order; insert_point3!(T, v); end
    return T
end

# Break exact coplanar/cospherical degeneracies with a deterministic, per-vertex
# symbolic perturbation of magnitude ~1e-8·(bbox diagonal). This is the concrete
# form of SoS: with the points in general position the exact predicates return
# definite signs and the tetrahedralization has NO flat (zero-volume) tets. The
# shift is far below any meshing tolerance (e.g. the enclosure's 1.2e-6 m), and
# coincident points were already merged, so shared vertices stay shared.
function _perturb3!(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})
    n = length(x); n == 0 && return
    xmin,xmax = extrema(x); ymin,ymax = extrema(y); zmin,zmax = extrema(z)
    diag = sqrt((xmax-xmin)^2 + (ymax-ymin)^2 + (zmax-zmin)^2)
    diag == 0 && (diag = 1.0)
    eps = 1e-8 * diag
    @inbounds for i in 1:n
        s = UInt64(i)*0x9E3779B97F4A7C15 + 0xD1B54A32D192ED03
        r() = (s ⊻= s<<13; s ⊻= s>>7; s ⊻= s<<17; (Float64(s >> 11)/Float64(1<<53)) - 0.5)
        x[i] += eps*r(); y[i] += eps*r(); z[i] += eps*r()
    end
    return
end

function _dedup3(xs,ys,zs)
    n=length(xs); seen=Dict{NTuple{3,Float64},Int32}()
    ux=Float64[]; uy=Float64[]; uz=Float64[]; remap=Vector{Int32}(undef,n)
    @inbounds for i in 1:n
        key=(xs[i],ys[i],zs[i]); id=get(seen,key,Int32(0))
        if id==0; push!(ux,xs[i]);push!(uy,ys[i]);push!(uz,zs[i]); id=Int32(length(ux)); seen[key]=id; end
        remap[i]=id
    end
    return ux,uy,uz,remap
end

function _shuffle3!(a::Vector{Int32}, seed::UInt64)
    s=seed
    @inbounds for i in length(a):-1:2
        s += 0x9E3779B97F4A7C15; z=s
        z=(z ⊻ (z>>30))*0xBF58476D1CE4E5B9; z=(z ⊻ (z>>27))*0x94D049BB133111EB; z=z ⊻ (z>>31)
        j = Int(z % UInt64(i)) + 1; a[i],a[j]=a[j],a[i]
    end
    a
end

"""
    to_mesh3(T; keep=nothing) -> Mesh

Export live real tets (ghosts dropped), nodes compacted canonically by coordinate.
If `keep` (a per-tet `Bool` mask indexed by tet id) is given, only kept tets are
exported — used to drop exterior tets after domain classification.
"""
function to_mesh3(T::Triangulation3; keep::Union{Nothing,AbstractVector{Bool}}=nothing)
    # internal orientation is orient3(v1,v2,v3,v4) > 0, which is tet_signed_volume
    # < 0 (opposite convention). Swap v3,v4 on export so MeshTypes sees positive
    # signed volumes (the geometric standard `validate` checks).
    tets = NTuple{4,Int32}[]
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t)) || continue
        (keep !== nothing && !(t <= length(keep) && keep[t])) && continue
        push!(tets, (_vert(T,t,1),_vert(T,t,2),_vert(T,t,4),_vert(T,t,3)))
    end
    used = Set{Int32}(); for te in tets, v in te; push!(used, v); end
    order = sort(collect(used); by=v->(T.x[v],T.y[v],T.z[v]))
    nid = Dict{Int32,Int32}(); for (k,v) in enumerate(order); nid[v]=Int32(k); end
    coords = Matrix{Float64}(undef,3,length(order))
    @inbounds for (k,v) in enumerate(order); coords[1,k]=T.x[v]; coords[2,k]=T.y[v]; coords[3,k]=T.z[v]; end
    M = Matrix{Int32}(undef,4,length(tets))
    @inbounds for (j,te) in enumerate(tets); for i in 1:4; M[i,j]=nid[te[i]]; end; end
    return Mesh(coords; tets=M)
end

# ════════════════════════════════════════════════════════════════════════════════
# Local mesh modification: 2-3 / 3-2 flips (Stage-3/4 kernel primitives)
# ════════════════════════════════════════════════════════════════════════════════

# Replace a small region: create each tet in `newverts` (oriented positive),
# link its faces to the recorded outer neighbours (`outer`: sorted-face → outside
# tet) and, for faces shared between the new tets, to each other. The reversed-
# shared-face invariant holds automatically for positively-oriented tets. Returns
# the new tet ids.
function _rebuild_region!(T::Triangulation3, newverts, outer::Dict{NTuple{3,Int32},Int32})
    internal = Dict{NTuple{3,Int32}, Tuple{Int32,Int32}}()
    ids = Int32[]
    for v in newverts
        a,b,c,d = v
        orient3(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,d)) < 0 && ((c,d)=(d,c))   # → positive
        t = _newtet!(T, a, b, c, d); push!(ids, t)
        for k in 1:4
            f = _face(T, t, k); key = _sort3t(f[1],f[2],f[3])
            if haskey(outer, key)
                nb = outer[key]; _setn!(T, t, k, nb)
                if nb != 0
                    j = _nslot_by_face(T, nb, f[1], f[2], f[3]); j != 0 && _setn!(T, nb, j, t)
                end
            elseif haskey(internal, key)
                (t2, s2) = internal[key]; _setn!(T, t, k, t2); _setn!(T, t2, Int(s2), t); delete!(internal, key)
            else
                internal[key] = (t, Int32(k))
            end
        end
    end
    isempty(internal) || error("Mesh3D._rebuild_region!: unmatched internal face")
    return ids
end

"""
    flip23!(T, t1, k1) -> Bool

2→3 flip across the face opposite local vertex `k1` of tet `t1`: the two tets
sharing that face (apexes `d`, `e`) become three tets sharing the new edge `d–e`.
Performed only if the bipyramid is convex (else it would create inverted tets);
returns `false` if not flippable. Neither tet may be a ghost.
"""
function flip23!(T::Triangulation3, t1::Integer, k1::Integer)
    t2 = _nbr(T, t1, k1)
    (t2 == 0 || _is_ghost_tet(T, t1) || _is_ghost_tet(T, t2)) && return false
    d = _vert(T, t1, k1)
    f = _face(T, t1, k1); a,b,c = f[1], f[2], f[3]
    j2 = _nslot(T, t2, t1); e = _vert(T, t2, j2)
    # convexity: the three candidate tets must share one orientation sign
    o1 = orient3(_pt(T,a),_pt(T,b),_pt(T,d),_pt(T,e))
    o2 = orient3(_pt(T,b),_pt(T,c),_pt(T,d),_pt(T,e))
    o3 = orient3(_pt(T,c),_pt(T,a),_pt(T,d),_pt(T,e))
    (o1 != 0 && o2 != 0 && o3 != 0 && (o1>0)==(o2>0) && (o2>0)==(o3>0)) || return false
    outer = Dict{NTuple{3,Int32},Int32}()
    for k in 1:4
        k == k1 && continue
        ff = _face(T, t1, k); outer[_sort3t(ff[1],ff[2],ff[3])] = _nbr(T, t1, k)
    end
    for k in 1:4
        k == j2 && continue
        ff = _face(T, t2, k); outer[_sort3t(ff[1],ff[2],ff[3])] = _nbr(T, t2, k)
    end
    _killtet!(T, t1); _killtet!(T, t2)
    ids = _rebuild_region!(T, ((a,b,d,e),(b,c,d,e),(c,a,d,e)), outer)
    T.last = ids[1]
    return true
end

@inline function _tet_mindihedral(T::Triangulation3, t)
    a=_pt(T,_vert(T,t,1)); b=_pt(T,_vert(T,t,2)); c=_pt(T,_vert(T,t,3)); d=_pt(T,_vert(T,t,4))
    mn, _ = tet_dihedral_extrema(a,b,c,d); mn
end

"""
    optimize_flips!(T; passes=4, tol=1e-9) -> n_flips

Local quality optimization by hill-climbing **2-3 flips**: a flip is kept only if
it strictly increases the minimum dihedral angle of the tets it touches, otherwise
it is reverted by its exact inverse (3-2 flip). Volume and validity are preserved
and the global minimum dihedral is **non-decreasing** (a safe optimizer). The
result is generally not Delaunay. Effect is *limited* on a Delaunay mesh — 2-3
flips rarely improve the local minimum, and slivers are not removed by flips alone
(that needs 3-2-driven collapse + exudation, remaining Stage-4 work). Returns the
number of flips applied.
"""
function optimize_flips!(T::Triangulation3; passes::Integer=4, tol::Real=1e-9)
    nflips = 0
    for _ in 1:passes
        changed = false
        # 2→3 flips over interior faces
        @inbounds for t1 in 1:length(T.alive)
            (T.alive[t1] && !_is_ghost_tet(T,t1)) || continue
            for k1 in 1:4
                t2 = _nbr(T, t1, k1)
                (t2 != 0 && !_is_ghost_tet(T,t2)) || continue
                before = min(_tet_mindihedral(T,t1), _tet_mindihedral(T,t2))
                d = _vert(T,t1,k1); e = _vert(T,t2,_nslot(T,t2,t1))
                if flip23!(T, t1, k1)
                    ring = tets_around_edge(T, d, e)
                    after = minimum(_tet_mindihedral(T,t) for t in ring)
                    if after > before + tol
                        nflips += 1; changed = true
                    else
                        flip32!(T, d, e, ring)                   # revert
                    end
                    break     # t1 was killed by the flip (kept or reverted) — stop its face loop
                end
            end
        end
        changed || break
    end
    return nflips
end

"""
    tets_around_edge(T, u, v) -> Vector{Int32}

All live (real) tets incident to the undirected edge `(u,v)`.
"""
function tets_around_edge(T::Triangulation3, u::Integer, v::Integer)
    out = Int32[]
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t)) || continue
        has_u=false; has_v=false
        for k in 1:4
            vv=_vert(T,t,k); vv==u && (has_u=true); vv==v && (has_v=true)
        end
        (has_u && has_v) && push!(out, Int32(t))
    end
    return out
end

"""
    flip32!(T, edge_de, tets3) -> Bool

3→2 flip: the three tets in `tets3` sharing edge `(d,e)` become two tets sharing
the triangle `(a,b,c)` of the other three vertices. `tets3` must be exactly the
three tets around the interior edge `(d,e)`. Returns `false` if not applicable.
"""
function flip32!(T::Triangulation3, d::Integer, e::Integer, tets3)
    length(tets3) == 3 || return false
    any(t -> _is_ghost_tet(T, t), tets3) && return false
    # the "ring" vertices a,b,c = the vertices other than d,e across the three tets
    ring = Int32[]
    for t in tets3, k in 1:4
        v = _vert(T, t, k)
        (v == d || v == e) && continue
        v in ring || push!(ring, v)
    end
    length(ring) == 3 || return false
    a,b,c = ring[1], ring[2], ring[3]
    # collect outer faces: the faces of the three tets that do NOT contain edge (d,e)
    outer = Dict{NTuple{3,Int32},Int32}()
    for t in tets3, k in 1:4
        ff = _face(T, t, k)
        (d in ff && e in ff) && continue          # internal (d,e)-faces vanish
        nb = _nbr(T, t, k)
        (nb in tets3) && continue                  # face between two ring tets — also internal
        outer[_sort3t(ff[1],ff[2],ff[3])] = nb
    end
    for t in tets3; _killtet!(T, t); end
    ids = _rebuild_region!(T, ((a,b,c,d),(a,b,c,e)), outer)
    T.last = ids[1]
    return true
end

# ════════════════════════════════════════════════════════════════════════════════
# Boundary / interface recovery (conforming Delaunay via Steiner points)
# ════════════════════════════════════════════════════════════════════════════════

function _add_vertex3!(T::Triangulation3, x::Float64, y::Float64, z::Float64)
    push!(T.x, x); push!(T.y, y); push!(T.z, z); T.nreal += 1; push!(T.vtet, Int32(0))
    return Int32(T.nreal)
end

@inline _sortface(a,b,c) = _sort3t(Int32(a),Int32(b),Int32(c))

"""Set of all triangular faces present in the current tetrahedralization (sorted)."""
function present_faces(T::Triangulation3)
    S = Set{NTuple{3,Int32}}()
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t)) || continue
        for k in 1:4
            f = _face(T,t,k)
            (_is_ghost_v(f[1])||_is_ghost_v(f[2])||_is_ghost_v(f[3])) && continue
            push!(S, _sortface(f[1],f[2],f[3]))
        end
    end
    return S
end

@inline _face_present(T, present, a, b, c) = _sortface(a,b,c) in present

# ════════════════════════════════════════════════════════════════════════════════
# Domain filling from a boundary surface (Stage-3 volume meshing)
# ════════════════════════════════════════════════════════════════════════════════

"""
    tetrahedralize(surface::Mesh; rng_seed=1) -> Mesh

Fill the volume enclosed by a closed triangulated `surface` with tetrahedra:
Delaunay-tetrahedralize the surface vertices, then keep the tets whose centroid
lies inside the surface (robust point-in-polyhedron by ray casting against the
original surface). Exact for convex domains; for non-convex domains the boundary
is the Delaunay restriction (conforming boundary recovery refines this — see
[`tetrahedralize_recover`](@ref)). Returns the interior tet [`Mesh`](@ref).
"""
function tetrahedralize(surface::Mesh; rng_seed::Integer=1)
    nn = size(surface.coords, 2)
    xs = Vector{Float64}(undef, nn); ys = similar(xs); zs = similar(xs)
    @inbounds for i in 1:nn; xs[i]=surface.coords[1,i]; ys[i]=surface.coords[2,i]; zs[i]=surface.coords[3,i]; end
    T = delaunay3d(xs, ys, zs; rng_seed=rng_seed)
    keep = _classify_by_centroid(T, surface)
    return to_mesh3(T; keep=keep)
end

"""
    tetrahedralize_multi(surfaces; rng_seed=1) -> Mesh

Fill each closed boundary surface in `surfaces` independently and combine the
results into a single multi-region tet [`Mesh`](@ref), tagging each region's tets
with its 1-based index (`tet_tag`). Shared interfaces are meshed per region (not
conforming across regions — that needs constrained interface recovery), but every
volume is genuinely filled: the standing anti-false-positive check is that *each*
region contributes a positive tet count (see `tets_per_region`).
"""
function tetrahedralize_multi(surfaces::AbstractVector{Mesh}; rng_seed::Integer=1)
    coords_cols = Vector{NTuple{3,Float64}}()
    tetcols = Vector{NTuple{4,Int32}}()
    tags = Int32[]
    for (r, surf) in enumerate(surfaces)
        m = tetrahedralize(surf; rng_seed=rng_seed)
        off = Int32(length(coords_cols))
        @inbounds for i in 1:size(m.coords,2)
            push!(coords_cols, (m.coords[1,i], m.coords[2,i], m.coords[3,i]))
        end
        @inbounds for t in 1:size(m.tets,2)
            push!(tetcols, (m.tets[1,t]+off, m.tets[2,t]+off, m.tets[3,t]+off, m.tets[4,t]+off))
            push!(tags, Int32(r))
        end
    end
    coords = Matrix{Float64}(undef, 3, length(coords_cols))
    @inbounds for (i,c) in enumerate(coords_cols); coords[1,i]=c[1]; coords[2,i]=c[2]; coords[3,i]=c[3]; end
    tets = Matrix{Int32}(undef, 4, length(tetcols))
    @inbounds for (j,te) in enumerate(tetcols); for i in 1:4; tets[i,j]=te[i]; end; end
    return Mesh(coords; tets=tets, tet_tag=tags)
end

"""
    tets_per_region(m) -> Dict{Int32,Int}

Tet count per `tet_tag` region — the mandatory anti-false-positive check: a valid
multi-region volume mesh has a *positive* count in every region (DEVELOPMENT.md:
"No empty volumes" is vacuous if a region is empty).
"""
function tets_per_region(m::Mesh)
    d = Dict{Int32,Int}()
    @inbounds for t in 1:size(m.tets,2)
        tag = m.tet_tag[t]; d[tag] = get(d, tag, 0) + 1
    end
    return d
end

# per-tet keep mask: true iff the tet's centroid is inside the surface
function _classify_by_centroid(T::Triangulation3, surface::Mesh)
    keep = falses(length(T.alive))
    dir = _unitn((1.0, 0.3141592653589793, 0.01720209895))   # generic, avoids grazing
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t)) || continue
        a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3);d=_vert(T,t,4)
        cx=(T.x[a]+T.x[b]+T.x[c]+T.x[d])/4
        cy=(T.y[a]+T.y[b]+T.y[c]+T.y[d])/4
        cz=(T.z[a]+T.z[b]+T.z[c]+T.z[d])/4
        keep[t] = _inside_surface((cx,cy,cz), dir, surface)
    end
    return keep
end

# ray-cast parity test: is p inside the closed surface?
function _inside_surface(p, dir, surface::Mesh)
    crossings = 0
    @inbounds for f in 1:size(surface.tris, 2)
        a=surface.tris[1,f]; b=surface.tris[2,f]; c=surface.tris[3,f]
        pa=(surface.coords[1,a],surface.coords[2,a],surface.coords[3,a])
        pb=(surface.coords[1,b],surface.coords[2,b],surface.coords[3,b])
        pc=(surface.coords[1,c],surface.coords[2,c],surface.coords[3,c])
        _ray_hits_tri(p, dir, pa, pb, pc) && (crossings += 1)
    end
    return isodd(crossings)
end

# Möller–Trumbore ray/triangle intersection for t > 0 (positive ray direction).
@inline function _ray_hits_tri(o, d, v0, v1, v2)
    e1 = _subn(v1, v0); e2 = _subn(v2, v0)
    pv = _cross(d, e2); det = _dot(e1, pv)
    abs(det) < 1e-300 && return false           # ray parallel to triangle plane
    inv = 1.0/det
    tv = _subn(o, v0)
    u = _dot(tv, pv)*inv
    (u < 0.0 || u > 1.0) && return false
    qv = _cross(tv, e1)
    v = _dot(d, qv)*inv
    (v < 0.0 || u+v > 1.0) && return false
    t = _dot(e2, qv)*inv
    return t > 1e-12                              # strictly ahead of the origin
end

"""
    check_consistency3(T) -> (ok, msg)

Real tets are positively oriented; every neighbour link is mutual and the shared
face is stored as a reversed (odd) permutation; no dangling/dead references.
"""
function check_consistency3(T::Triangulation3)
    @inbounds for t in eachindex(T.alive)
        T.alive[t] || continue
        gs = _ghost_slot(T, t)
        a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3);d=_vert(T,t,4)
        (a==b||a==c||a==d||b==c||b==d||c==d) && return (false, "degenerate tet $t")
        if gs == 0
            orient3(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,d)) > 0 || return (false, "tet $t not positive")
        end
        for k in 1:4
            s = _nbr(T,t,k)
            s == 0 && return (false, "tet $t open face slot $k")
            T.alive[s] || return (false, "tet $t nbr $s dead")
            j = _nslot(T, s, t); j == 0 && return (false, "tet $t→$s not mutual")
            f = _face(T,t,k); g = _face(T,s,j)
            _sort3t(f...) == _sort3t(g...) || return (false, "tet $t/$s face set mismatch")
            # reversed orientation: f is an odd permutation of g
            _same_cyclic(f, g) && return (false, "tet $t/$s shared face not reversed")
        end
    end
    return (true, "ok")
end

# do (a,b,c) and (d,e,f) represent the same cyclic (even) orientation of a triangle?
@inline function _same_cyclic(f, g)
    (f[1],f[2],f[3]) == (g[1],g[2],g[3]) && return true
    (f[1],f[2],f[3]) == (g[2],g[3],g[1]) && return true
    (f[1],f[2],f[3]) == (g[3],g[1],g[2]) && return true
    return false
end

"""
    is_delaunay3(T) -> (ok, nviol)

Exact empty-circumsphere oracle: no real vertex is strictly inside any real tet's
circumsphere (`insphere_sos`). This *defines* 3-D Delaunay.
"""
function is_delaunay3(T::Triangulation3)
    nviol = 0
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t)) || continue
        a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3);d=_vert(T,t,4)
        for v in 1:T.nreal
            (v==a||v==b||v==c||v==d) && continue
            if insphere_sos(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,d),_pt(T,v), a,b,c,d,Int32(v)) > 0
                nviol += 1
            end
        end
    end
    return (nviol==0, nviol)
end

end # module Mesh3D
