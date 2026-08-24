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

using ..Predicates: orient3, orient3_sos, insphere_sos, incircle3_sos, orient2
using ..MeshTypes: Mesh, tet_dihedral_extrema, validate, is_closed_manifold, boundary_edges, tet_signed_volume
using ..Mesh2D: constrained_delaunay, to_mesh
using ..ExactMesh3D: delaunay3d_exact
using ..SizeField: AbstractSizeField, ConstantSize, metric_edge_length,
                   directional_size

export Triangulation3, delaunay3d, tetrahedralize, tetrahedralize_multi,
       tetrahedralize_conforming, tetrahedralize_conforming_exact, tets_per_region, mesh_box, mesh_box_regions, BoxRegion,
       recover_boundary, mesh_boolean, mesh_sized_conforming, mesh_cylinder
export check_consistency3, is_delaunay3, to_mesh3, ntets_live, present_faces
export flip23!, flip32!, tets_around_edge, optimize_flips!, refine_to_size
export insert_steiner3, recover_segment3, recover_triangle3
export mesh_covers_segment3, mesh_covers_triangle3

# Extension hook installed by the downstream RecoverCDT module after both modules
# have loaded.  Keeping the declaration here avoids an include cycle while letting
# the public volume front door use exact conforming-Delaunay recovery as its final
# robust fallback.
function _recover_boundary_exact end
function _recover_partition_exact end

const GHOST3 = Int32(0)

@inline function _finite3(x::Real, caller::AbstractString, name::AbstractString)
    y = try
        Float64(x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable: $(sprint(showerror, err))"))
    end
    isfinite(y) || throw(ArgumentError("$caller: $name must be finite (got $x)"))
    return y
end

@inline function _mesh3_entity_context(value, caller::AbstractString)
    value===nothing && return nothing
    value isa Tuple && length(value)==2 &&
        value[1] isa Integer && !(value[1] isa Bool) &&
        value[2] isa Integer && !(value[2] isa Bool) ||
        throw(ArgumentError("$caller: an entity context must be nothing or a " *
                            "(dimension, tag) integer tuple (got $value)"))
    dim=try
        Int(value[1])
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError) || rethrow()
        throw(ArgumentError("$caller: entity dimension is outside the platform Int range"))
    end
    tag=try
        Int(value[2])
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError) || rethrow()
        throw(ArgumentError("$caller: entity tag is outside the platform Int range"))
    end
    dim in 0:3 || throw(ArgumentError("$caller: entity dimension must be in 0:3"))
    tag>0 || throw(ArgumentError("$caller: entity tag must be positive"))
    return (dim,tag)
end

function _mesh3_vertex_entities(m::Mesh,resolver,caller::AbstractString)
    resolver===nothing && return nothing
    count=size(m.coords,2)
    resolver isa AbstractVector && length(resolver)!=count && throw(ArgumentError(
        "$caller: vertex_entities must have one entry per mesh node"))
    out=Vector{Union{Nothing,Tuple{Int,Int}}}(undef,count)
    @inbounds for index in 1:count
        point=(m.coords[1,index],m.coords[2,index],m.coords[3,index])
        raw=if resolver isa Tuple
            resolver
        elseif resolver isa AbstractVector
            resolver[index]
        elseif resolver isa AbstractDict
            get(resolver,index,nothing)
        elseif applicable(resolver,index,point)
            resolver(index,point)
        else
            throw(ArgumentError(
                "$caller: vertex_entities must be a point-entity tuple, per-node " *
                "vector/dictionary, or callable (index, point)"))
        end
        context=_mesh3_entity_context(raw,"$caller vertex $index")
        (context===nothing || context[1]==0) || throw(ArgumentError(
            "$caller: vertex_entities entry $index must have dimension 0"))
        out[index]=context
    end
    return out
end

function _mesh3_tag_entities(tags::AbstractVector{Int32}, resolver, dim::Int,
                             caller::AbstractString;required::Bool)
    resolver===nothing && return nothing
    out=Dict{Int32,Union{Nothing,Tuple{Int,Int}}}()
    for tag in tags
        haskey(out,tag) && continue
        raw=if resolver isa AbstractDict
            key=(dim,Int(tag))
            if haskey(resolver,key)
                resolver[key]
            elseif dim==3 && haskey(resolver,tag)
                resolver[tag]
            elseif dim==3 && haskey(resolver,Int(tag))
                resolver[Int(tag)]
            elseif required
                throw(ArgumentError(
                    "$caller: entity_resolver has no entry for dimension-$dim cell tag $tag"))
            else
                nothing
            end
        elseif applicable(resolver,dim,Int(tag))
            resolver(dim,Int(tag))
        else
            throw(ArgumentError("$caller: entity_resolver must be a dictionary keyed by " *
                                "tet_tag/(dimension, cell_tag), or a callable (dimension, tag)"))
        end
        out[tag]=_mesh3_entity_context(raw,caller)
    end
    return out
end

@inline function _seed3(seed::Integer, caller::AbstractString)
    seed isa Bool && throw(ArgumentError("$caller: rng_seed must not be Bool"))
    (0 <= seed <= typemax(Int)) ||
        throw(ArgumentError("$caller: rng_seed must be in 0:$(typemax(Int)) (got $seed)"))
    return Int(seed)
end

@inline function _bounded_count3(value::Integer,caller::AbstractString,
                                 name::AbstractString;minimum::Int=0)
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    converted=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name exceeds the platform Int range"))
    end
    converted>=minimum || throw(ArgumentError(
        "$caller: $name must be ≥ $minimum"))
    return converted
end

@inline function _ceil_count3(x::Float64, caller::AbstractString, name::AbstractString;
                              minimum::Int=1)
    (isfinite(x) && x >= 0) ||
        throw(ArgumentError("$caller: computed $name count is not finite and non-negative"))
    value=try
        ceil(Int,x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: requested $name count exceeds the platform Int limit"))
    end
    return max(minimum,value)
end

@inline function _checked_mul3(caller::AbstractString, name::AbstractString, xs::Int...)
    value = 1
    try
        for x in xs
            value = Base.checked_mul(value, x)
        end
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: requested $name count overflows the platform Int limit"))
    end
    return value
end

@inline function _checked_add3(caller::AbstractString, name::AbstractString, x::Int, y::Int)
    try
        return Base.checked_add(x,y)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: requested $name count overflows the platform Int limit"))
    end
end

@inline function _midpoint3(a::Float64, b::Float64)
    # Same-sign subtraction cannot overflow; opposite-sign half-sums cannot overflow.
    return signbit(a) == signbit(b) ? a + (b-a)/2 : a/2 + b/2
end

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
        length(T.alive)<typemax(Int32) ||
            throw(ArgumentError("Mesh3D: tetrahedron storage exceeds the Int32 indexing limit"))
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
    # Pick a non-degenerate coordinate projection for the orientation sign.  The
    # circle predicate itself retains the full 3-D Euclidean metric exactly; a raw
    # affine projection would turn an oblique circle into an ellipse.
    proj=if orient2((pu[1],pu[2]),(pv[1],pv[2]),(pw[1],pw[2]))!=0
        (1,2)
    elseif orient2((pu[2],pu[3]),(pv[2],pv[3]),(pw[2],pw[3]))!=0
        (2,3)
    elseif orient2((pu[3],pu[1]),(pv[3],pv[1]),(pw[3],pw[1]))!=0
        (3,1)
    else
        return false
    end
    pr(q)=(q[proj[1]],q[proj[2]])
    a2=pr(pu); b2=pr(pv); c2=pr(pw)
    ic = incircle3_sos(pu,pv,pw,_pt(T,p),u,v,w,Int32(p))
    o2 = orient2(a2,b2,c2)                            # face triangle orientation in-plane
    # inside ⇔ in-circle sign matches the CCW/CW orientation of (u,v,w)
    return (ic > 0) == (o2 > 0)
end

@inline _subn(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _dot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline function _unitn(a)
    l = hypot(a[1], a[2], a[3])
    (isfinite(l) && l > 0) || throw(ArgumentError("Mesh3D: direction must have finite positive length"))
    return (a[1]/l,a[2]/l,a[3]/l)
end

# ── construction ────────────────────────────────────────────────────────────────
function Triangulation3(xs::Vector{Float64}, ys::Vector{Float64}, zs::Vector{Float64})
    n = length(xs)
    (length(ys)==n && length(zs)==n) ||
        throw(ArgumentError("Triangulation3: coordinate lengths differ ($(length(xs)), $(length(ys)), $(length(zs)))"))
    n <= typemax(Int32) ||
        throw(ArgumentError("Triangulation3: $n points exceed the Int32 indexing limit"))
    @inbounds for i in 1:n
        (isfinite(xs[i]) && isfinite(ys[i]) && isfinite(zs[i])) ||
            throw(ArgumentError("Triangulation3: point $i has a non-finite coordinate"))
    end
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
        if _noncollinear3(_pt(T,a),_pt(T,b),_pt(T,i)); c=Int32(i); break; end
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

@inline _noncollinear3(a,b,c) =
    orient2((a[1],a[2]),(b[1],b[2]),(c[1],c[2]))!=0 ||
    orient2((a[2],a[3]),(b[2],b[3]),(c[2],c[3]))!=0 ||
    orient2((a[3],a[1]),(b[3],b[1]),(c[3],c[1]))!=0

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

# ── point location (jump-and-walk) ──────────────────────────────────────────────
@inline function _dist2_v(T::Triangulation3, v, px, py, pz)
    @inbounds dx=T.x[v]-px; @inbounds dy=T.y[v]-py; @inbounds dz=T.z[v]-pz
    hypot(dx,dy,dz)
end

# Pick a walk-start tet by jump-and-walk (Mücke–Saias–Zhu): sample ~∛n already-
# inserted landmark vertices and start from the live incident tet of whichever is
# nearest the query. This keeps the subsequent walk short (O(∛n) expected) and —
# the point of it — prevents the walk from wandering the whole mesh (then falling
# to the O(n) exhaustive `_locate_scan`) when the query is far from `T.last`, the
# pathology that made the exact-coordinate (perturb=false) Delaunay O(n²) on
# maximally-cospherical input (e.g. a fine cylinder's ~500 hull vertices).
# Output-identical: the walk's terminal tet is a member of the query's unique
# Bowyer–Watson cavity regardless of start, so the mesh is unchanged — only faster.
# Deterministic: the landmark sample stream is keyed by `vid`.
function _pick_start3(T::Triangulation3, px, py, pz, vid)
    start = T.last
    (start==0 || start>length(T.alive) || !T.alive[start]) && (start=_first_alive3(T))
    n = T.nreal
    n <= 4 && return start
    bestd = Inf
    # seed the "nearest" contest with the current start's real vertices so the
    # sample can only improve on it (never regress below the T.last locality).
    @inbounds for k in 1:4
        v = _vert(T, start, k)
        (v==GHOST3) && continue
        d = _dist2_v(T, v, px, py, pz); d < bestd && (bestd = d)
    end
    k = clamp(round(Int, cbrt(n)) + 3, 4, 64)
    rs = UInt64(vid)*0x9E3779B97F4A7C15 + 0xD1B54A32D192ED03
    @inbounds for _ in 1:k
        rs ⊻= rs<<13; rs ⊻= rs>>7; rs ⊻= rs<<17
        v = Int32(rs % UInt64(n)) + Int32(1)
        vt = T.vtet[v]
        (vt==0 || vt>length(T.alive) || !T.alive[vt]) && continue
        d = _dist2_v(T, v, px, py, pz)
        if d < bestd; bestd = d; start = vt; end
    end
    return start
end

function locate3(T::Triangulation3, px, py, pz, vid; start::Integer=0)
    t = start == 0 ? _pick_start3(T, px, py, pz, vid) : Int32(start)
    (t==0 || t>length(T.alive) || !T.alive[t]) && (t=_first_alive3(T))
    gs = _ghost_slot(T, t); gs != 0 && (t = _nbr(T, t, gs))
    prev = Int32(0)
    rs = UInt64(vid)*0x9E3779B97F4A7C15 + 0xD1B54A32D192ED03
    # Guard for the walk. With the jump-and-walk near start, a walk on a valid
    # Delaunay mesh converges in O(∛n) steps (measured ≤100 even at n=20 000). A
    # walk that exceeds a generous multiple of that is NOT converging — the
    # perturb=false complex is degenerate on maximally-cospherical input (flat
    # tets), where a visibility walk has no termination guarantee — so cut the
    # losses and locate by exhaustive scan rather than wander the whole mesh
    # (16·ntets steps, the old bound, was ~28 000 wasted steps × expensive exact
    # predicates per hard point, the O(n²) cliff on the fine cylinder). Purely a
    # performance knob: the scan returns the same uniquely-determined tet, so the
    # mesh is unchanged (CRC-identical) — only faster.
    guard = 0; maxstep = max(1000, 48*round(Int, cbrt(length(T.alive))))
    @inbounds while true
        guard += 1
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
    1 <= vid <= T.nreal ||
        throw(ArgumentError("insert_point3!: vertex $vid is outside 1:$(T.nreal)"))
    T.vtet[vid] == 0 || throw(ArgumentError("insert_point3!: vertex $vid is already inserted"))
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
    seed = _seed3(rng_seed, "delaunay3d")
    (length(ys)==length(xs) && length(zs)==length(xs)) ||
        throw(ArgumentError("delaunay3d: coordinate lengths differ ($(length(xs)), $(length(ys)), $(length(zs)))"))
    length(xs) <= typemax(Int32) ||
        throw(ArgumentError("delaunay3d: $(length(xs)) points exceed the Int32 indexing limit"))
    @inbounds for i in eachindex(xs)
        (isfinite(xs[i]) && isfinite(ys[i]) && isfinite(zs[i])) ||
            throw(ArgumentError("delaunay3d: point $i has a non-finite coordinate"))
    end
    ux,uy,uz,_ = _dedup3(xs,ys,zs)
    perturb && _perturb3!(ux,uy,uz)
    T = Triangulation3(ux,uy,uz)
    placed = _init3!(T)
    isempty(placed) && return T
    ps = Set{Int32}(placed)
    order = Int32[i for i in 1:T.nreal if !(Int32(i) in ps)]
    _shuffle3!(order, UInt64(seed))
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
    diag = hypot(xmax-xmin, ymax-ymin, zmax-zmin)
    isfinite(diag) || throw(ArgumentError(
        "delaunay3d: point-cloud extent is too large for finite Float64 perturbation; use perturb=false"))
    diag == 0 && (diag = 1.0)
    eps = 1e-8 * diag
    @inbounds for i in 1:n
        s = UInt64(i)*0x9E3779B97F4A7C15 + 0xD1B54A32D192ED03
        r() = (s ⊻= s<<13; s ⊻= s>>7; s ⊻= s<<17; (Float64(s >> 11)/Float64(1<<53)) - 0.5)
        x[i] += eps*r(); y[i] += eps*r(); z[i] += eps*r()
        (isfinite(x[i]) && isfinite(y[i]) && isfinite(z[i])) ||
            throw(ArgumentError("delaunay3d: deterministic perturbation overflowed at point $i"))
    end
    return
end

function _dedup3(xs,ys,zs)
    n=length(xs); seen=Dict{NTuple{3,Float64},Int32}()
    ux=Float64[]; uy=Float64[]; uz=Float64[]; remap=Vector{Int32}(undef,n)
    @inbounds for i in 1:n
        key=(xs[i]==0 ? 0.0 : xs[i], ys[i]==0 ? 0.0 : ys[i], zs[i]==0 ? 0.0 : zs[i])
        id=get(seen,key,Int32(0))
        if id==0; push!(ux,key[1]);push!(uy,key[2]);push!(uz,key[3]); id=Int32(length(ux)); seen[key]=id; end
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
    keep === nothing || length(keep) == length(T.alive) ||
        throw(ArgumentError("to_mesh3: keep length $(length(keep)) does not match tet storage length $(length(T.alive))"))
    # internal orientation is orient3(v1,v2,v3,v4) > 0, which is tet_signed_volume
    # < 0 (opposite convention). Swap v3,v4 on export so MeshTypes sees positive
    # signed volumes (the geometric standard `validate` checks).
    tets = NTuple{4,Int32}[]
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t)) || continue
        (keep !== nothing && !keep[t]) && continue
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
    1 <= t1 <= length(T.alive) ||
        throw(ArgumentError("flip23!: tet $t1 is outside 1:$(length(T.alive))"))
    1 <= k1 <= 4 || throw(ArgumentError("flip23!: local face slot must be in 1:4 (got $k1)"))
    T.alive[t1] || return false
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

Local quality optimization by hill-climbing **2-3 and 3-2 flips**: a flip is kept
only if it strictly increases the minimum dihedral angle of the tets it touches,
otherwise reverted by its exact inverse. Volume and validity are preserved and the
global minimum dihedral is **non-decreasing** (a safe optimizer); the result is
generally not Delaunay. The 3-2 flips collapse many slivers — on a random-cloud
Delaunay this cuts the sliver count substantially and raises the mean dihedral
(e.g. 25.6°→28.8° from flips alone, 25.6°→37.8° combined with smoothing). A few
topology-locked slivers can remain; combine this with
`Optimize.remove_slivers`, which targets poor vertex stars geometrically.
Returns the number of flips applied.
"""
function optimize_flips!(T::Triangulation3; passes::Integer=4, tol::Real=1e-9)
    (0<=passes<=typemax(Int)) || throw(ArgumentError(
        "optimize_flips!: passes must be in 0:$(typemax(Int)) (got $passes)"))
    npasses=Int(passes)
    ftol = _finite3(tol, "optimize_flips!", "tol")
    ftol >= 0 || throw(ArgumentError("optimize_flips!: tol must be non-negative (got $tol)"))
    nflips = 0
    for _ in 1:npasses
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
                    if after > before + ftol
                        nflips += 1; changed = true
                    else
                        flip32!(T, d, e, ring)                   # revert
                    end
                    break     # t1 was killed by the flip (kept or reverted) — stop its face loop
                end
            end
        end
        # 3→2 flips over interior edges shared by exactly 3 tets (can collapse slivers)
        seen = Set{NTuple{2,Int32}}()
        @inbounds for t0 in 1:length(T.alive)
            (T.alive[t0] && !_is_ghost_tet(T,t0)) || continue
            for (i,j) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
                u = _vert(T,t0,i); v = _vert(T,t0,j)
                key = _edgekey(u, v)
                key in seen && continue; push!(seen, key)
                ring = tets_around_edge(T, u, v)
                (length(ring) == 3 && !any(t->_is_ghost_tet(T,t), ring)) || continue
                # ring "apex" triangle a,b,c = the vertices other than u,v
                abc = Int32[]
                for t in ring, k in 1:4
                    w=_vert(T,t,k); (w==u||w==v||w in abc) && continue; push!(abc,w)
                end
                length(abc) == 3 || continue
                before = minimum(_tet_mindihedral(T,t) for t in ring)
                if flip32!(T, u, v, ring)
                    new2 = _face_tets(T, abc[1], abc[2], abc[3])
                    if length(new2) == 2
                        after = min(_tet_mindihedral(T,new2[1]), _tet_mindihedral(T,new2[2]))
                        if after > before + ftol
                            nflips += 1; changed = true
                        else
                            k = _nslot_by_face(T, new2[1], abc[1], abc[2], abc[3])   # revert via 2-3
                            k != 0 && flip23!(T, new2[1], k)
                        end
                    end
                    break     # t0's edges are stale after the flip
                end
            end
        end
        changed || break
    end
    return nflips
end

# does any live ghost tet contain both edge endpoints (⇒ boundary edge)?
function _ghost_touches_edge(T::Triangulation3, u, v)
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && _is_ghost_tet(T,t)) || continue
        hu=false; hv=false
        for k in 1:4
            w=_vert(T,t,k); w==u && (hu=true); w==v && (hv=true)
        end
        (hu && hv) && return true
    end
    return false
end

# the (≤2) live real tets that contain all of a,b,c as vertices
function _face_tets(T::Triangulation3, a, b, c)
    out = Int32[]
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t)) || continue
        na=false; nb=false; nc=false
        for k in 1:4
            v=_vert(T,t,k); v==a && (na=true); v==b && (nb=true); v==c && (nc=true)
        end
        (na && nb && nc) && push!(out, Int32(t))
    end
    return out
end

"""
    tets_around_edge(T, u, v) -> Vector{Int32}

All live (real) tets incident to the undirected edge `(u,v)`.
"""
function tets_around_edge(T::Triangulation3, u::Integer, v::Integer)
    (1 <= u <= T.nreal && 1 <= v <= T.nreal) ||
        throw(ArgumentError("tets_around_edge: endpoints ($u,$v) must lie in 1:$(T.nreal)"))
    u != v || throw(ArgumentError("tets_around_edge: endpoints must be distinct"))
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
    (1 <= d <= T.nreal && 1 <= e <= T.nreal && d != e) ||
        throw(ArgumentError("flip32!: edge endpoints ($d,$e) must be distinct vertices in 1:$(T.nreal)"))
    length(unique(tets3)) == 3 || return false
    all(t -> 1 <= t <= length(T.alive) && T.alive[t], tets3) || return false
    any(t -> _is_ghost_tet(T, t), tets3) && return false
    # the edge must be interior: no ghost tet may touch it (else the 3 real tets do
    # not close the ring and the flip would relink to the wrong neighbours).
    _ghost_touches_edge(T, d, e) && return false
    # the "ring" vertices a,b,c = the vertices other than d,e across the three tets
    ring = Int32[]
    for t in tets3, k in 1:4
        v = _vert(T, t, k)
        (v == d || v == e) && continue
        v in ring || push!(ring, v)
    end
    length(ring) == 3 || return false
    a,b,c = ring[1], ring[2], ring[3]
    # validity: d and e must be on opposite sides of the ring triangle (a,b,c),
    # else the two replacement tets would overlap/invert.
    od = orient3(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,d))
    oe = orient3(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,e))
    (od != 0 && oe != 0 && (od > 0) != (oe > 0)) || return false
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
# Longest-edge bisection — size refinement of a tet mesh to a guaranteed max edge
# ════════════════════════════════════════════════════════════════════════════════

# Minimal binary MAX-heap on (violation ratio, a, b). Vertices never move; queued
# records remain valid until their edge is split or disappears from the incidence map.
struct _EHeap
    d::Vector{Tuple{Float64,Int32,Int32}}
end
_EHeap() = _EHeap(Tuple{Float64,Int32,Int32}[])
@inline function _hpush!(h::_EHeap, x)
    d = h.d; push!(d, x); i = length(d)
    @inbounds while i > 1
        p = i >> 1; d[p][1] >= d[i][1] && break; d[p], d[i] = d[i], d[p]; i = p
    end
end
@inline function _hpop!(h::_EHeap)
    d = h.d; top = d[1]; last = pop!(d)
    if !isempty(d)
        d[1] = last; i = 1; n = length(d)
        @inbounds while true
            l = 2i; r = 2i + 1; big = i
            l <= n && d[l][1] > d[big][1] && (big = l)
            r <= n && d[r][1] > d[big][1] && (big = r)
            big == i && break; d[i], d[big] = d[big], d[i]; i = big
        end
    end
    top
end

"""
    refine_to_size(m::Mesh, hmax) -> Mesh
    refine_to_size(m::Mesh, field::AbstractSizeField) -> Mesh

Refine a valid tet mesh by longest-edge subdivision.  The scalar overload guarantees
**every edge is `≤ hmax`**.  The field overload requires every edge's metric length,
sampled at its two endpoints and midpoint, to be at most one; edges are processed by
descending metric violation ratio. This is `length / local_target` for scalar fields
and the directional `√(dᵀMd)` criterion for anisotropic fields.

`entity` supplies one fixed geometric-entity context. For a classified mesh,
`entity_resolver` may instead be a dictionary keyed by `(dimension, cell_tag)`
(with bare `tet_tag` keys retained for compatibility), or a callable
`(dimension, cell_tag) -> entity`. Every edge is evaluated in each distinct
incident segment, triangle, and live-tetrahedron context and the largest metric
violation is used. Interface and boundary edges are therefore refined whenever
any adjacent geometric entity requires it. Dictionary entries are mandatory for
tetrahedron tags and optional for lower-dimensional tags. The two keywords are
mutually exclusive. This explicit mapping is necessary because legacy cell tags
can be physical tags and are not assumed to be geometric entity tags.
`vertex_entities` independently classifies original mesh nodes as Gmsh point
entities. It accepts one `(0, tag)` tuple, a per-node vector, a sparse dictionary
keyed by node index, or a callable `(index, point) -> entity`; `nothing` entries
are unclassified. Point contexts seed only original classified boundary edges:
they are sampled at the corresponding endpoint on the first split, while child
edges and newly inserted edge-interior nodes use their curve/surface/volume
contexts. This mirrors Gmsh's use of vertex sizes during boundary discretization
without turning a point value into a recursive radial volume-size constraint.

Splitting all tets around one edge keeps the mesh **conforming** (the shared faces are
split identically in every incident tet). The representable midpoint is used normally.
Its three independently rounded coordinates can differ infinitesimally from the exact
represented line; that unique canonical point is accepted only after exact child-
orientation checks. If it coincides with an existing vertex, a checked representable
interior point is selected that makes both children nondegenerate in every incident
tet. Additional ULP offsets are permitted only for an unconstrained, single-region
interior edge star, whose cavity boundary remains unchanged; boundary, interface, and
explicit feature edges require an exactly collinear fallback point or return a
representability blocker. For a constant positive target, every selected point lies
strictly inside the edge's coordinate bounds, so repeated subdivision terminates at
`maxedge ≤ hmax`. A spatial field is evaluated on every new edge and terminates
whenever its sampled targets remain resolvable in `Float64`.
Interior edges are refined too (no interior lattice needed). Boundary vertices are
preserved and boundary edges stay on the boundary to Float64 midpoint resolution.
**Region tags (`tet_tag`) are propagated** — each child inherits its parent
tet's tag — so multi-region meshes keep their partition. This is the size-refinement
terminator behind [`Tessella.mesh_sized`](@ref).
"""
function refine_to_size(m::Mesh,hmax::Real;entity=nothing,entity_resolver=nothing,
                        vertex_entities=nothing)
    target = _finite3(hmax, "refine_to_size", "hmax")
    target > 0 || throw(ArgumentError("refine_to_size: hmax must be positive (got $hmax)"))
    return _refine_to_size(m,ConstantSize(target);target_description="hmax=$target",
                           entity=entity,entity_resolver=entity_resolver,
                           vertex_entities=vertex_entities)
end

function refine_to_size(m::Mesh,field::AbstractSizeField;entity=nothing,
                        entity_resolver=nothing,vertex_entities=nothing)
    return _refine_to_size(m,field;target_description=string(nameof(typeof(field))),
                           entity=entity,entity_resolver=entity_resolver,
                           vertex_entities=vertex_entities)
end

function _refine_to_size(m::Mesh, field::AbstractSizeField;
                         target_description::AbstractString="size field",entity=nothing,
                         entity_resolver=nothing,vertex_entities=nothing)
    nt0 = size(m.tets,2)
    nt0 > 0 || throw(ArgumentError("refine_to_size: input contains no tetrahedra"))
    inputdiag = validate(m)
    inputdiag.ok || throw(ArgumentError("refine_to_size: input mesh is invalid — " *
                                        join(inputdiag.messages, "; ")))
    nv0 = size(m.coords, 2)
    entity!==nothing && entity_resolver!==nothing && throw(ArgumentError(
        "refine_to_size: entity and entity_resolver are mutually exclusive"))
    fixed_entity=_mesh3_entity_context(entity,"refine_to_size")
    vertex_contexts=_mesh3_vertex_entities(m,vertex_entities,"refine_to_size")
    tet_entities=_mesh3_tag_entities(m.tet_tag,entity_resolver,3,"refine_to_size";
                                     required=true)
    cx = Vector{Float64}(undef, nv0); cy = similar(cx); cz = similar(cx)
    @inbounds for i in 1:nv0; cx[i]=m.coords[1,i]; cy[i]=m.coords[2,i]; cz[i]=m.coords[3,i]; end
    @inline coordkey(x,y,z)=(x==0 ? 0.0 : x,y==0 ? 0.0 : y,z==0 ? 0.0 : z)
    coordids=Dict{NTuple{3,Float64},Int32}()
    @inbounds for i in 1:nv0; coordids[coordkey(cx[i],cy[i],cz[i])]=Int32(i); end
    tv = Vector{NTuple{4,Int32}}(); alive = Bool[]; ttag = Int32[]      # ttag: per-tet region tag
    tags0 = m.tet_tag                                                    # input tags (length = ntets)
    inc = Dict{Tuple{Int32,Int32},Vector{Int32}}()
    heap = _EHeap()
    queued = Set{Tuple{Int32,Int32}}()
    deferred = Set{Tuple{Int32,Int32}}()
    deferred_signature = Dict{Tuple{Int32,Int32},Tuple}()
    @inline ek(a,b) = a<=b ? (Int32(a),Int32(b)) : (Int32(b),Int32(a))
    @inline elen(a,b) = hypot(cx[a]-cx[b], cy[a]-cy[b], cz[a]-cz[b])
    # Lower-dimensional cells participate both in conformity updates and in field
    # classification. Resolve their tags once; missing dictionary entries mean
    # "unclassified at this dimension" while tetrahedron entries remain required.
    segv = NTuple{2,Int32}[(m.segs[1,s],m.segs[2,s]) for s in axes(m.segs,2)]
    segalive = trues(length(segv)); segtag = copy(m.seg_tag)
    triv = NTuple{3,Int32}[(m.tris[1,f],m.tris[2,f],m.tris[3,f]) for f in axes(m.tris,2)]
    trialive = trues(length(triv)); tritag = copy(m.tri_tag)
    seg_entities=_mesh3_tag_entities(m.seg_tag,entity_resolver,1,"refine_to_size";
                                     required=false)
    tri_entities=_mesh3_tag_entities(m.tri_tag,entity_resolver,2,"refine_to_size";
                                     required=false)
    seginc = Dict{Tuple{Int32,Int32},Vector{Int}}()
    triinc = Dict{Tuple{Int32,Int32},Vector{Int}}()
    @inbounds for (s,q) in enumerate(segv); push!(get!(() -> Int[],seginc,ek(q...)),s); end
    @inbounds for (f,q) in enumerate(triv)
        for e in (ek(q[1],q[2]),ek(q[2],q[3]),ek(q[3],q[1]))
            push!(get!(() -> Int[],triinc,e),f)
        end
    end
    @inline function vertex_score(a,b)
        vertex_contexts===nothing && return 0.0
        (a<=nv0 && b<=nv0) || return 0.0
        e=ek(a,b)
        (haskey(seginc,e) || haskey(triinc,e)) || return 0.0
        context=vertex_contexts[a]
        context===nothing && return 0.0
        p=(cx[a],cy[a],cz[a]);q=(cx[b],cy[b],cz[b])
        direction=(q[1]-p[1],q[2]-p[2],q[3]-p[3])
        distance=hypot(direction...)
        distance==0 && return 0.0
        return distance/directional_size(field,p,direction;entity=context)
    end
    @inline function edge_score(a,b)
        p=(cx[a],cy[a],cz[a]);q=(cx[b],cy[b],cz[b])
        e=ek(a,b)
        best=if tet_entities===nothing
            metric_edge_length(field,p,q;entity=fixed_entity)
        else
            ts=get(inc,e,nothing)
            local_score=0.0
            if ts!==nothing
                @inbounds for i in eachindex(ts)
                    t=ts[i];alive[t]||continue
                    tag=ttag[t]
                    duplicate=false
                    for j in firstindex(ts):i-1
                        u=ts[j]
                        if alive[u] && ttag[u]==tag
                            duplicate=true;break
                        end
                    end
                    duplicate&&continue
                    local_score=max(local_score,metric_edge_length(
                        field,p,q;entity=tet_entities[tag]))
                end
            end
            local_score
        end
        surface_cells=get(triinc,e,nothing)
        if surface_cells!==nothing && tri_entities!==nothing
            @inbounds for i in eachindex(surface_cells)
                f=surface_cells[i];trialive[f]||continue;tag=tritag[f]
                context=tri_entities[tag];context===nothing&&continue
                duplicate=false
                for j in firstindex(surface_cells):i-1
                    g=surface_cells[j]
                    if trialive[g] && tritag[g]==tag
                        duplicate=true;break
                    end
                end
                duplicate&&continue
                best=max(best,metric_edge_length(field,p,q;entity=context))
            end
        end
        curve_cells=get(seginc,e,nothing)
        if curve_cells!==nothing && seg_entities!==nothing
            @inbounds for i in eachindex(curve_cells)
                s=curve_cells[i];segalive[s]||continue;tag=segtag[s]
                context=seg_entities[tag];context===nothing&&continue
                duplicate=false
                for j in firstindex(curve_cells):i-1
                    g=curve_cells[j]
                    if segalive[g] && segtag[g]==tag
                        duplicate=true;break
                    end
                end
                duplicate&&continue
                best=max(best,metric_edge_length(field,p,q;entity=context))
            end
        end
        return max(best,vertex_score(a,b),vertex_score(b,a))
    end
    function queueedge!(u,v,len)
        (isfinite(len) && len>=0) || throw(ArgumentError(
            "refine_to_size: edge length is not finite"))
        e=ek(u,v)
        (e in queued || e in deferred) && return nothing
        score=edge_score(u,v)
        if score>1
            push!(queued,e);_hpush!(heap,(score,e[1],e[2]))
        end
        return nothing
    end

    # Lower-dimensional cells are part of Mesh's public data contract.  Whenever a
    # tet edge is bisected, split every coincident segment/triangle edge as well so
    # boundary/feature cells and all of their tags remain conforming instead of being
    # silently discarded.
    function addseg!(q::NTuple{2,Int32}, tag::Int32)
        length(segv)<typemax(Int32) ||
            throw(ErrorException("refine_to_size: segment count exceeds the Int32 working limit"))
        push!(segv,q); push!(segalive,true); push!(segtag,tag); s=length(segv)
        push!(get!(() -> Int[],seginc,ek(q...)),s)
        return nothing
    end
    function addtri!(q::NTuple{3,Int32}, tag::Int32)
        length(triv)<typemax(Int32) ||
            throw(ErrorException("refine_to_size: triangle count exceeds the Int32 working limit"))
        push!(triv,q); push!(trialive,true); push!(tritag,tag); f=length(triv)
        for e in (ek(q[1],q[2]),ek(q[2],q[3]),ek(q[3],q[1]))
            push!(get!(() -> Int[],triinc,e),f)
        end
        return nothing
    end
    function split_lower!(a::Int32,b::Int32,mid::Int32)
        e=ek(a,b)
        for s in get(seginc,e,Int[])
            segalive[s] || continue
            q=segv[s]; ek(q...)==e || continue
            segalive[s]=false; g=segtag[s]
            if q[1]==a
                addseg!((a,mid),g); addseg!((mid,b),g)
            else
                addseg!((b,mid),g); addseg!((mid,a),g)
            end
        end
        for f in get(triinc,e,Int[])
            trialive[f] || continue
            q=triv[f]
            hit=0
            for i in 1:3
                j=i==3 ? 1 : i+1
                ek(q[i],q[j])==e && (hit=i; break)
            end
            hit==0 && continue
            j=hit==3 ? 1 : hit+1; k=j==3 ? 1 : j+1
            x=q[hit]; y=q[j]; z=q[k]; g=tritag[f]; trialive[f]=false
            addtri!((x,mid,z),g); addtri!((mid,y,z),g)
        end
        delete!(seginc,e); delete!(triinc,e)
        return nothing
    end

    # add a tet (oriented positive) carrying region tag `tag`; register its 6 edges; if
    # `push_edges` enqueues every edge that violates its sampled local target.
    function addtet!(a, b, c, d, push_edges::Bool, tag::Int32,
                     volume_hint::Union{Nothing,Float64})
        length(tv) < typemax(Int32) ||
            throw(ErrorException("refine_to_size: working tet count exceeds the Int32 incidence limit"))
        sv=volume_hint===nothing ?
           tet_signed_volume((cx[a],cy[a],cz[a]),(cx[b],cy[b],cz[b]),
                             (cx[c],cy[c],cz[c]),(cx[d],cy[d],cz[d])) : volume_hint
        if !(isfinite(sv) && sv != 0)
            pts=((cx[a],cy[a],cz[a]),(cx[b],cy[b],cz[b]),
                 (cx[c],cy[c],cz[c]),(cx[d],cy[d],cz[d]))
            exactsign=-orient3(pts...)
            throw(ErrorException(
                "refine_to_size: bisection produced a non-finite or flat tetrahedron " *
                "at vertices ($a,$b,$c,$d), signed_volume=$sv, exact_sign=$exactsign, points=$pts"))
        end
        sv < 0 && ((c,d)=(d,c))
        push!(tv, (Int32(a),Int32(b),Int32(c),Int32(d))); push!(alive, true); push!(ttag, tag); id = Int32(length(tv))
        q = (a,b,c,d)
        @inbounds for (i,j) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
            u,v = q[i],q[j]; push!(get!(inc, ek(u,v)) do; Int32[] end, id)
            len=elen(u,v)
            isfinite(len) || throw(ErrorException("refine_to_size: an edge length is non-finite"))
            push_edges && queueedge!(u,v,len)
        end
        id
    end
    @inbounds for t in 1:size(m.tets,2); addtet!(Int(m.tets[1,t]),Int(m.tets[2,t]),Int(m.tets[3,t]),Int(m.tets[4,t]), false, tags0[t], nothing); end
    for (e,_) in inc
        len=elen(e[1],e[2]); isfinite(len) || throw(ErrorException("refine_to_size: an edge length is non-finite"))
        queueedge!(e[1],e[2],len)
    end

    # Return the two vertices opposite edge (a,b), ordered so (a,b,u,v) has the
    # same positive orientation as the stored tet. This derives the sign from
    # permutation parity and avoids recomputing the parent's exact volume.
    function edge_opposites(q::NTuple{4,Int32},a::Int32,b::Int32)
        u=Int32(0);v=Int32(0)
        for w in q;(w!=a&&w!=b)&&(u==0 ? (u=w) : (v=w));end
        u!=0&&v!=0||return nothing
        r=(a,b,u,v)
        pos=ntuple(i -> q[1]==r[i] ? 1 : q[2]==r[i] ? 2 : q[3]==r[i] ? 3 : 4,4)
        inversions=0
        for i in 1:3,j in i+1:4;pos[i]>pos[j]&&(inversions+=1);end
        if isodd(inversions)
            u,v=v,u
        end
        return (u,v)
    end

    # Certify p = ((den-num)*a + num*b)/den exactly in represented Float64
    # coordinates without allocating Rational{BigInt}. Float64 significands have at
    # most 53 bits; the conservative exponent-gap limit below leaves ample Int128
    # headroom for coefficients through 32 and the two-term sum.
    function exact_affine_coordinate(a::Float64,b::Float64,p::Float64,
                                     num::Int,den::Int)
        0<num<den<=32 || return false
        ma,ea,sa=Base.decompose(a);mb,eb,sb=Base.decompose(b)
        mp,ep,sp=Base.decompose(p)
        lp=Int128(sp)*Int128(mp)*den
        la=Int128(sa)*Int128(ma)*(den-num)
        lb=Int128(sb)*Int128(mb)*num
        (lp==0&&la==0&&lb==0) && return true
        emin=typemax(Int);emax=typemin(Int)
        if lp!=0;emin=min(emin,ep);emax=max(emax,ep);end
        if la!=0;emin=min(emin,ea);emax=max(emax,ea);end
        if lb!=0;emin=min(emin,eb);emax=max(emax,eb);end
        emax-emin<=60 || return false
        lp=lp==0 ? Int128(0) : lp<<(ep-emin)
        la=la==0 ? Int128(0) : la<<(ea-emin)
        lb=lb==0 ? Int128(0) : lb<<(eb-emin)
        return lp==la+lb
    end
    exact_affine_point(a::Int32,b::Int32,p,num::Int,den::Int)=
        exact_affine_coordinate(cx[a],cx[b],p[1],num,den) &&
        exact_affine_coordinate(cy[a],cy[b],p[2],num,den) &&
        exact_affine_coordinate(cz[a],cz[b],p[3],num,den)

    exactly_on_edge(a::Int32,b::Int32,p)=
        orient2((cx[a],cy[a]),(cx[b],cy[b]),(p[1],p[2]))==0 &&
        orient2((cx[a],cz[a]),(cx[b],cz[b]),(p[1],p[3]))==0 &&
        orient2((cy[a],cz[a]),(cy[b],cz[b]),(p[2],p[3]))==0

    function split_children_valid(a::Int32,b::Int32,ts::Vector{Int32},p,
                                  fraction::Union{Nothing,Tuple{Int,Int}})
        affine=fraction!==nothing && exact_affine_point(a,b,p,fraction...)
        @inbounds for t in ts
            alive[t]||return false
            opposite=edge_opposites(tv[t],a,b);opposite===nothing&&return false
            u,v=opposite
            if affine
                # The stored parent (a,b,u,v) is positive. Exact affine linearity
                # proves both child determinants positive without another predicate.
            else
                c1=-orient3((cx[a],cy[a],cz[a]),p,
                            (cx[u],cy[u],cz[u]),(cx[v],cy[v],cz[v]))
                c2=-orient3(p,(cx[b],cy[b],cz[b]),
                            (cx[u],cy[u],cz[u]),(cx[v],cy[v],cz[v]))
                (c1>0&&c2>0)||return false
            end
        end
        return true
    end

    # An off-edge point changes a piecewise-linear boundary or material interface,
    # even when its displacement is only one ULP. It is safe only in a closed,
    # single-region interior edge star: the cavity boundary then consists solely of
    # faces opposite the edge and is retained verbatim by the two child fans.
    function edge_is_constrained(a::Int32,b::Int32,ts::Vector{Int32})
        e=ek(a,b)
        any(s -> segalive[s],get(seginc,e,Int[])) && return true
        any(f -> trialive[f],get(triinc,e,Int[])) && return true
        faces=Dict{NTuple{3,Int32},Vector{Int32}}()
        @inbounds for t in ts
            alive[t] || return true
            opposite=edge_opposites(tv[t],a,b)
            opposite===nothing && return true
            for u in opposite
                push!(get!(() -> Int32[],faces,_sort3t(a,b,u)),ttag[t])
            end
        end
        return any(tags -> length(tags)!=2 || tags[1]!=tags[2],values(faces))
    end

    # A rounded midpoint can equal an existing vertex or land exactly in a child
    # face plane even when the parent is exactly nondegenerate. Search dyadic
    # positions near 1/2 and retain the most balanced valid candidate, preferring
    # fewer ULP shifts. The fused convex interpolation avoids Rational allocations.
    function fallback_split_point(a::Int32,b::Int32,ts::Vector{Int32},constrained::Bool)
        best=nothing;bestscore=nothing
        lerp(x,y,t)=signbit(x)==signbit(y) ? muladd(t,y-x,x) : (1-t)*x+t*y
        shiftulp(x,k)=begin
            y=x
            if k<0;for _ in 1:-k;y=prevfloat(y);end
            elseif k>0;for _ in 1:k;y=nextfloat(y);end
            end
            y
        end
        for radius in 0:2,n in 1:31
            n==16&&continue
            t=n/32
            base=(lerp(cx[a],cx[b],t),lerp(cy[a],cy[b],t),lerp(cz[a],cz[b],t))
            for ox in -radius:radius,oy in -radius:radius,oz in -radius:radius
                max(abs(ox),abs(oy),abs(oz))==radius||continue
                ((cx[a]==cx[b]&&ox!=0)||(cy[a]==cy[b]&&oy!=0)||
                 (cz[a]==cz[b]&&oz!=0))&&continue
                p=(shiftulp(base[1],ox),shiftulp(base[2],oy),shiftulp(base[3],oz))
                all(min((cx[a],cy[a],cz[a])[i],(cx[b],cy[b],cz[b])[i])<=p[i]<=
                    max((cx[a],cy[a],cz[a])[i],(cx[b],cy[b],cz[b])[i]) for i in 1:3)||continue
            (p!=(cx[a],cy[a],cz[a])&&p!=(cx[b],cy[b],cz[b]))||continue
            haskey(coordids,coordkey(p...))&&continue
            constrained && !exactly_on_edge(a,b,p) && continue
            split_children_valid(a,b,ts,p,(n,32))||continue
            score=(abs(n-16),radius,abs(ox)+abs(oy)+abs(oz),ox,oy,oz)
            if bestscore===nothing || score<bestscore
                best=p;bestscore=score
            end
            end
        end
        best===nothing&&return nothing
        return best::NTuple{3,Float64}
    end

    retriangulation_reason=Ref("")

    # Remove an unsplittable interior edge by retriangulating its closed star. This
    # is a checked 3-2/4-4/general fan flip: boundary faces and boundary-side
    # orientations must be identical, all new internal faces must separate opposite
    # vertices, volumes must agree to Float64 resolution, and region tags must match.
    function retriangulate_edge!(a::Int32,b::Int32,ts::Vector{Int32})
        e=ek(a,b)
        if any(s -> segalive[s],get(seginc,e,Int[]))
            retriangulation_reason[]="edge belongs to a segment cell";return false
        end
        if any(f -> trialive[f],get(triinc,e,Int[]))
            retriangulation_reason[]="edge belongs to a triangle cell";return false
        end
        n=length(ts)
        3<=n<=12||(retriangulation_reason[]="edge-star size $n is outside 3:12";return false)
        tags=unique(Int32[ttag[t] for t in ts])
        length(tags)==1||(retriangulation_reason[]="edge star crosses region tags $tags";return false)

        linkedges=NTuple{2,Int32}[]
        @inbounds for t in ts
            alive[t]||(retriangulation_reason[]="edge star contains a dead tet";return false)
            rem=Int32[v for v in tv[t] if v!=a&&v!=b]
            length(rem)==2||(retriangulation_reason[]="edge incidence is inconsistent";return false)
            push!(linkedges,ek(rem[1],rem[2]))
        end
        length(unique(linkedges))==n||
            (retriangulation_reason[]="edge link contains duplicate edges";return false)
        adj=Dict{Int32,Vector{Int32}}()
        for (u,v) in linkedges
            push!(get!(() -> Int32[],adj,u),v);push!(get!(() -> Int32[],adj,v),u)
        end
        length(adj)==n&&all(length(unique(ns))==2 for ns in values(adj))||
            (retriangulation_reason[]="edge link is not a simple cycle";return false)
        start=minimum(keys(adj));cycle=Int32[start];prev=Int32(0);cur=start
        for k in 2:n
            ns=sort!(unique(adj[cur]));nxt=prev==0 ? ns[1] : (ns[1]==prev ? ns[2] : ns[1])
            (nxt!=start&&!(nxt in cycle))||
                (retriangulation_reason[]="edge link closes early";return false)
            push!(cycle,nxt);prev,cur=cur,nxt
        end
        start in adj[cur]||(retriangulation_reason[]="edge link does not close";return false)

        function face_opposites(cells)
            faces=Dict{NTuple{3,Int32},Vector{Int32}}()
            for q in cells,k in 1:4
                f=ntuple(i -> q[i<k ? i : i+1],3)
                push!(get!(() -> Int32[],faces,_sort3t(f...)),q[k])
            end
            faces
        end
        oldcells=NTuple{4,Int32}[tv[t] for t in ts]
        oldfaces=face_opposites(oldcells)
        all(length(v)<=2 for v in values(oldfaces))||
            (retriangulation_reason[]="old cavity is non-manifold";return false)
        oldboundary=Set(k for (k,v) in oldfaces if length(v)==1)
        oldvolume=sum(abs(tet_signed_volume(
            (cx[q[1]],cy[q[1]],cz[q[1]]),(cx[q[2]],cy[q[2]],cz[q[2]]),
            (cx[q[3]],cy[q[3]],cz[q[3]]),(cx[q[4]],cy[q[4]],cz[q[4]]))) for q in oldcells)
        (isfinite(oldvolume)&&oldvolume>0)||
            (retriangulation_reason[]="old cavity volume is invalid";return false)
        cavityverts=unique(Int32[v for q in oldcells for v in q])
        scale=0.0
        for i in eachindex(cavityverts),j in i+1:length(cavityverts)
            scale=max(scale,elen(cavityverts[i],cavityverts[j]))
        end
        voltol=max(1e-10*oldvolume,256eps(Float64)*scale^3)

        bestcells=nothing;bestvolumes=nothing;bestscore=Inf
        stats=Dict(:flat=>0,:duplicate=>0,:incidence=>0,:boundary=>0,
                   :orientation=>0,:volume=>0,:accepted=>0)
        for root in 1:n
            order=Int32[cycle[mod1(root+j,n)] for j in 0:n-1]
            cells=NTuple{4,Int32}[];volumes=Float64[];ok=true
            for j in 2:n-1
                x,y,z=order[1],order[j],order[j+1]
                for q0 in ((a,x,y,z),(b,x,z,y))
                    q=q0
                    sv=tet_signed_volume((cx[q[1]],cy[q[1]],cz[q[1]]),
                                         (cx[q[2]],cy[q[2]],cz[q[2]]),
                                         (cx[q[3]],cy[q[3]],cz[q[3]]),
                                         (cx[q[4]],cy[q[4]],cz[q[4]]))
                    if !(isfinite(sv)&&sv!=0);ok=false;break end
                    sv<0&&(q=(q[1],q[2],q[4],q[3]);sv=-sv)
                    push!(cells,q);push!(volumes,sv)
                end
                ok||break
            end
            if !ok;stats[:flat]+=1;continue end
            if length(Set(Tuple(sort(collect(q))) for q in cells))!=length(cells)
                stats[:duplicate]+=1;continue
            end
            newfaces=face_opposites(cells)
            if !all(length(v)<=2 for v in values(newfaces))
                stats[:incidence]+=1;continue
            end
            if Set(k for (k,v) in newfaces if length(v)==1)!=oldboundary
                stats[:boundary]+=1;continue
            end
            for (f,opps) in newfaces
                fp=ntuple(i -> (cx[f[i]],cy[f[i]],cz[f[i]]),3)
                if length(opps)==2
                    s1=orient3(fp...,(cx[opps[1]],cy[opps[1]],cz[opps[1]]))
                    s2=orient3(fp...,(cx[opps[2]],cy[opps[2]],cz[opps[2]]))
                    (s1!=0&&s2!=0&&s1==-s2)||(ok=false;break)
                else
                    oldopp=oldfaces[f][1]
                    sn=orient3(fp...,(cx[opps[1]],cy[opps[1]],cz[opps[1]]))
                    so=orient3(fp...,(cx[oldopp],cy[oldopp],cz[oldopp]))
                    (sn!=0&&sn==so)||(ok=false;break)
                end
            end
            if !ok;stats[:orientation]+=1;continue end
            drift=abs(sum(volumes)-oldvolume)
            if drift>voltol;stats[:volume]+=1;continue end
            stats[:accepted]+=1
            score=drift/voltol-min(volumes...)/max(volumes...)
            if score<bestscore
                bestcells=cells;bestvolumes=volumes;bestscore=score
            end
        end
        if bestcells===nothing
            linkpoints=Tuple((v,(cx[v],cy[v],cz[v])) for v in cycle)
            retriangulation_reason[]="no fan passed: stats=$stats, cycle=$linkpoints, " *
                                     "old_volume=$oldvolume, volume_tolerance=$voltol"
            return false
        end
        length(tv)+length(bestcells)<typemax(Int32)||
            (retriangulation_reason[]="replacement exceeds Int32 tet storage";return false)
        for t in ts;alive[t]=false;end
        tag=first(tags)
        for (q,sv) in zip(bestcells,bestvolumes)
            addtet!(q[1],q[2],q[3],q[4],true,tag,sv)
        end
        delete!(inc,e)
        return true
    end

    guard = 0
    while true
      while !isempty(heap.d)
        (_, a, b) = _hpop!(heap)
        delete!(queued,ek(a,b))
        e = ek(a, b); haskey(inc, e) || continue
        score=edge_score(a,b)
        score <= 1 && continue
        ts = Int32[]
        @inbounds for t in inc[e]
            alive[t] || continue
            q = tv[t]; ((a==q[1]||a==q[2]||a==q[3]||a==q[4]) && (b==q[1]||b==q[2]||b==q[3]||b==q[4])) || continue
            push!(ts, t)
        end
        isempty(ts) && continue
        guard += 1
        guard <= typemax(Int32) || throw(ErrorException("refine_to_size: bisection count exceeds Int32"))
        length(cx) < typemax(Int32) ||
            throw(ErrorException("refine_to_size: refined node count exceeds the Int32 indexing limit"))
        mx=_midpoint3(cx[a],cx[b]); my=_midpoint3(cy[a],cy[b]); mz=_midpoint3(cz[a],cz[b])
        ((mx,my,mz)!=(cx[a],cy[a],cz[a]) && (mx,my,mz)!=(cx[b],cy[b],cz[b])) ||
            throw(ErrorException("refine_to_size: $target_description is below Float64 coordinate resolution on edge ($a,$b)"))
        splitpoint=(mx,my,mz)
        midpoint_affine=exact_affine_point(a,b,splitpoint,1,2)
        constrained=!midpoint_affine && edge_is_constrained(a,b,ts)
        # `_midpoint3` is the canonical Float64 representation of the geometric
        # midpoint. On a general 3-D edge its independently rounded coordinates need
        # not be exactly collinear with the represented endpoints. Accept that unique
        # canonical point when both exact child-orientation tests pass. If it collides
        # with an existing vertex, however, the fallback may move by additional ULPs;
        # constrained edges then retain the strict exact-collinearity requirement in
        # `fallback_split_point`.
        splitvalid=!haskey(coordids,coordkey(splitpoint...)) &&
                   split_children_valid(a,b,ts,splitpoint,(1,2))
        if !splitvalid
            constrained=constrained||edge_is_constrained(a,b,ts)
            fallback=fallback_split_point(a,b,ts,constrained)
            if fallback===nothing
                retriangulate_edge!(a,b,ts)&&continue
                reason=
                    "refine_to_size: no conforming representable split or valid local " *
                    "edge-star retriangulation exists for edge ($a,$b) from " *
                    "$((cx[a],cy[a],cz[a])) to $((cx[b],cy[b],cz[b])) across " *
                    "$(length(ts)) incident tetrahedra $(Tuple(tv[t] for t in ts)) " *
                    "with tags $(Tuple(ttag[t] for t in ts)); local retriangulation: " *
                    retriangulation_reason[]
                constrained && throw(ErrorException(reason*"; constrained edge geometry cannot be moved"))
                signature=Tuple(sort(ts))
                if get(deferred_signature,e,nothing)==signature
                    throw(ErrorException(reason*"; neighboring refinement did not change the edge star"))
                end
                deferred_signature[e]=signature;push!(deferred,e)
                continue
            end
            splitpoint=fallback
            mx,my,mz=splitpoint
        end
        push!(cx,mx); push!(cy,my); push!(cz,mz)
        mvid = Int32(length(cx));coordids[coordkey(mx,my,mz)]=mvid;split_lower!(a,b,mvid)
        @inbounds for t in ts
            opposite=edge_opposites(tv[t],a,b)
            opposite===nothing&&throw(ErrorException(
                "refine_to_size: edge incidence changed during subdivision"))
            a1,a2=opposite
            g = ttag[t]; alive[t] = false                    # children inherit the parent tet's region tag
            addtet!(a,mvid,a1,a2,true,g,1.0)
            addtet!(mvid,b,a1,a2,true,g,1.0)
        end
        delete!(inc,e)
      end
      isempty(deferred)&&break
      pending=collect(deferred);empty!(deferred)
      for e in pending
          haskey(inc,e)||continue
          len=elen(e[1],e[2])
          isfinite(len)||throw(ErrorException("refine_to_size: a deferred edge length is non-finite"))
          queueedge!(e[1],e[2],len)
      end
      isempty(heap.d)&&break
    end
    # `inc` contains every generated edge. Split edges are deleted; other stale
    # entries can only reference dead tets. Certify each live output edge once here
    # instead of evaluating an expensive spatial field up to once per incident tet.
    @inbounds for ((u,v),ts) in inc
        live=false
        for t in ts
            if alive[t];live=true;break;end
        end
        live||continue
        score=edge_score(u,v)
        score<=1 || throw(ErrorException(
            "refine_to_size: postcondition failed: an output edge exceeds its local metric target from $target_description"))
    end
    keep = Int32[t for t in 1:length(tv) if alive[t]]
    used = Int32[]; seenv = Set{Int32}()
    for t in keep, v in tv[t]; (v in seenv) || (push!(used, v); push!(seenv, v)); end
    for s in eachindex(segv); segalive[s] || continue; for v in segv[s]; (v in seenv)||(push!(used,v);push!(seenv,v)); end; end
    for f in eachindex(triv); trialive[f] || continue; for v in triv[f]; (v in seenv)||(push!(used,v);push!(seenv,v)); end; end
    sort!(used)
    nid = Dict{Int32,Int32}(); for (i,v) in enumerate(used); nid[v]=Int32(i); end
    C = Matrix{Float64}(undef, 3, length(used))
    @inbounds for (i,v) in enumerate(used); C[1,i]=cx[v]; C[2,i]=cy[v]; C[3,i]=cz[v]; end
    M = Matrix{Int32}(undef, 4, length(keep))
    @inbounds for (j,t) in enumerate(keep); q=tv[t]; M[1,j]=nid[q[1]]; M[2,j]=nid[q[2]]; M[3,j]=nid[q[3]]; M[4,j]=nid[q[4]]; end
    tetags = Int32[ttag[t] for t in keep]
    keeps=Int[s for s in eachindex(segv) if segalive[s]]
    keepf=Int[f for f in eachindex(triv) if trialive[f]]
    S=Matrix{Int32}(undef,2,length(keeps)); st=Vector{Int32}(undef,length(keeps))
    @inbounds for (j,s) in enumerate(keeps); q=segv[s];S[1,j]=nid[q[1]];S[2,j]=nid[q[2]];st[j]=segtag[s];end
    F=Matrix{Int32}(undef,3,length(keepf)); ft=Vector{Int32}(undef,length(keepf))
    @inbounds for (j,f) in enumerate(keepf); q=triv[f];F[1,j]=nid[q[1]];F[2,j]=nid[q[2]];F[3,j]=nid[q[3]];ft[j]=tritag[f];end
    out=Mesh(C;segs=S,tris=F,tets=M,seg_tag=st,tri_tag=ft,tet_tag=tetags)
    diag=validate(out)
    diag.ok || throw(ErrorException("refine_to_size: produced an invalid mesh — " * join(diag.messages,"; ")))
    return out
end

# ════════════════════════════════════════════════════════════════════════════════
# Boundary / interface recovery (conforming Delaunay via Steiner points)
# ════════════════════════════════════════════════════════════════════════════════

function _add_vertex3!(T::Triangulation3, x::Float64, y::Float64, z::Float64)
    T.nreal < typemax(Int32) || throw(ArgumentError("Mesh3D: node count exceeds Int32"))
    (isfinite(x)&&isfinite(y)&&isfinite(z)) ||
        throw(ArgumentError("Mesh3D: cannot add a non-finite vertex"))
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

function _require_surface_input3(surface::Mesh,caller::AbstractString)
    size(surface.coords,2) <= typemax(Int32) ||
        throw(ArgumentError("$caller: surface node count exceeds the Int32 indexing limit"))
    size(surface.tris,2)>0 || throw(ArgumentError("$caller: surface has no triangles"))
    size(surface.tets,2)==0 || throw(ArgumentError(
        "$caller: a boundary surface must not contain tetrahedra"))
    d=validate(surface)
    d.ok || throw(ArgumentError("$caller: input surface is invalid — " * join(d.messages,"; ")))
    referenced=falses(size(surface.coords,2))
    @inbounds for f in axes(surface.tris,2),slot in 1:3
        referenced[Int(surface.tris[slot,f])]=true
    end
    missing=findfirst(!,referenced)
    missing===nothing || throw(ArgumentError(
        "$caller: surface node $missing is not referenced by a triangle"))
    return nothing
end

function _require_surface3(surface::Mesh, caller::AbstractString; oriented::Bool=false)
    _require_surface_input3(surface,caller)
    inc=Dict{NTuple{2,Int32},Int}()
    dirs=Dict{NTuple{2,Int32},Int}()
    faces=Set{NTuple{3,Int32}}()
    @inbounds for f in axes(surface.tris,2)
        a=surface.tris[1,f];b=surface.tris[2,f];c=surface.tris[3,f]
        key=_sort3t(a,b,c)
        key in faces && throw(ArgumentError("$caller: surface contains duplicate triangle $key"))
        push!(faces,key)
        for (u,v) in ((a,b),(b,c),(c,a))
            e=_edgekey(u,v);inc[e]=get(inc,e,0)+1
            dirs[(u,v)]=get(dirs,(u,v),0)+1
        end
    end
    for (e,n) in inc
        n==2 || throw(ArgumentError("$caller: surface edge $e has incidence $n (expected exactly 2)"))
        if oriented
            (get(dirs,(e[1],e[2]),0)==1 && get(dirs,(e[2],e[1]),0)==1) ||
                throw(ArgumentError("$caller: surface has inconsistent winding at edge $e"))
        end
    end
    return nothing
end

function _point3_input(raw,caller::AbstractString,name::AbstractString)
    (raw isa Tuple || raw isa AbstractVector) || throw(ArgumentError(
        "$caller: $name must be a tuple or vector with 3 real components"))
    length(raw)==3 || throw(ArgumentError(
        "$caller: $name must have exactly 3 components"))
    values=Tuple(raw)
    return ntuple(3) do component
        value=values[component]
        value isa Bool && throw(ArgumentError(
            "$caller: $name component $component must not be Bool"))
        value isa Real || throw(ArgumentError(
            "$caller: $name component $component must be real"))
        converted=try
            Float64(value)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: $name component $component must be Float64-representable"))
        end
        isfinite(converted) || throw(ArgumentError(
            "$caller: $name component $component must be finite"))
        converted
    end
end

function _interior_points3(raw,limit::Int,nn::Int,caller::AbstractString)
    raw===nothing && return NTuple{3,Float64}[]
    applicable(iterate,raw) || throw(ArgumentError(
        "$caller: interior_points must be an iterable of three-component points"))
    points=NTuple{3,Float64}[]
    for value in raw
        length(points)<limit || throw(ArgumentError(
            "$caller: interior_points exceed max_interior_points=$limit"))
        length(points)<typemax(Int32)-nn || throw(ArgumentError(
            "$caller: surface and interior points exceed the Int32 node limit"))
        index=length(points)+1
        push!(points,_point3_input(value,caller,"interior point $index"))
    end
    return points
end

function _node_at3(mesh::Mesh, p; atol=0.0)
    @inbounds for i in axes(mesh.coords,2)
        hypot(mesh.coords[1,i]-p[1],mesh.coords[2,i]-p[2],mesh.coords[3,i]-p[3])<=atol && return i
    end
    return 0
end

function _containing_tet(mesh::Mesh, p)
    @inbounds for t in axes(mesh.tets,2)
        a=(mesh.coords[1,mesh.tets[1,t]],mesh.coords[2,mesh.tets[1,t]],mesh.coords[3,mesh.tets[1,t]])
        b=(mesh.coords[1,mesh.tets[2,t]],mesh.coords[2,mesh.tets[2,t]],mesh.coords[3,mesh.tets[2,t]])
        c=(mesh.coords[1,mesh.tets[3,t]],mesh.coords[2,mesh.tets[3,t]],mesh.coords[3,mesh.tets[3,t]])
        d=(mesh.coords[1,mesh.tets[4,t]],mesh.coords[2,mesh.tets[4,t]],mesh.coords[3,mesh.tets[4,t]])
        volume=tet_signed_volume(a,b,c,d)
        volume==0 && continue
        subvolumes=(tet_signed_volume(p,b,c,d),tet_signed_volume(a,p,c,d),
                    tet_signed_volume(a,b,p,d),tet_signed_volume(a,b,c,p))
        # Intersections constructed in Float64 need not be bit-exactly coplanar.
        # Use a tolerance relative to this tet's own volume; the former absolute
        # `max(abs(volume),1)` scale collapsed every point in a sub-micron mesh.
        tolerance=1e-12*abs(volume)
        zeros=0; ok=true
        mask=(false,false,false,false)
        for (k,subvolume) in enumerate(subvolumes)
            if abs(subvolume)<=tolerance
                zeros+=1
                mask=k==1 ? (true,mask[2],mask[3],mask[4]) :
                     k==2 ? (mask[1],true,mask[3],mask[4]) :
                     k==3 ? (mask[1],mask[2],true,mask[4]) :
                            (mask[1],mask[2],mask[3],true)
            elseif signbit(subvolume)!=signbit(volume)
                ok=false; break
            end
        end
        ok && return t, zeros, mask
    end
    return 0, 0, (false,false,false,false)
end

function _finish_interior(mesh::Mesh, extra)
    isempty(extra) && return mesh
    m=insert_interior_points(mesh, extra)
    diag=validate(m)
    diag.ok || throw(ErrorException(
        "tetrahedralize: interior-point insertion produced an invalid mesh — "*
        join(diag.messages,"; ")))
    return m
end

function insert_interior_points(mesh::Mesh, extra)
    isempty(extra) && return mesh
    out=mesh
    for (k,p) in enumerate(extra)
        p isa NTuple{3,Float64} || throw(ArgumentError(
            "insert_interior_points: point $k is not a normalized Float64 triple"))
        out,_=_insert_steiner3(out,p)
    end
    return out
end

@inline function _pt3(mesh::Mesh, i::Integer)
    i=Int(i)
    return (mesh.coords[1,i], mesh.coords[2,i], mesh.coords[3,i])
end

function _on_segment3(x, p, q; atol=1e-12)
    vx,vy,vz=q[1]-p[1], q[2]-p[2], q[3]-p[3]
    wx,wy,wz=x[1]-p[1], x[2]-p[2], x[3]-p[3]
    L2=vx*vx+vy*vy+vz*vz
    L2>0 || return hypot(wx,wy,wz)<=atol
    cross=hypot(vy*wz-vz*wy, vz*wx-vx*wz, vx*wy-vy*wx)
    cross<=atol*sqrt(L2) || return false
    t=(wx*vx+wy*vy+wz*vz)/L2
    return -atol<=t<=1+atol
end

function _tet_edge_set(mesh::Mesh)
    edges=Set{NTuple{2,Int32}}()
    @inbounds for t in axes(mesh.tets,2)
        v=(mesh.tets[1,t],mesh.tets[2,t],mesh.tets[3,t],mesh.tets[4,t])
        for (i,j) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
            a,b=v[i],v[j]
            push!(edges, a<b ? (a,b) : (b,a))
        end
    end
    return edges
end

function mesh_covers_segment3(mesh::Mesh, p, q; atol=1e-12)
    a=_node_at3(mesh,p; atol=atol); b=_node_at3(mesh,q; atol=atol)
    a==0 && return false; b==0 && return false
    a==b && return true
    adj=Dict{Int,Vector{Int}}()
    for (u,v) in _tet_edge_set(mesh)
        pu=_pt3(mesh,u); pv=_pt3(mesh,v)
        (_on_segment3(pu,p,q; atol=atol) && _on_segment3(pv,p,q; atol=atol)) || continue
        push!(get!(Vector{Int}, adj, Int(u)), Int(v))
        push!(get!(Vector{Int}, adj, Int(v)), Int(u))
    end
    nn=size(mesh.coords,2)
    seen=falses(nn); qn=Int[a]; seen[a]=true; head=1
    while head<=length(qn)
        v=qn[head]; head+=1
        v==b && return true
        for u in get(adj, v, Int[])
            seen[u] && continue
            seen[u]=true; push!(qn,u)
        end
    end
    return false
end

function _cross3(u,v)
    return (u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1])
end
@inline _dot3(u,v)=u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
@inline _sub3(u,v)=(u[1]-v[1],u[2]-v[2],u[3]-v[3])
@inline _add3(u,v)=(u[1]+v[1],u[2]+v[2],u[3]+v[3])
@inline _scale3(u,s)=(u[1]*s,u[2]*s,u[3]*s)

function _point_in_triangle3(x, a, b, c; atol=1e-12)
    n=_cross3(_sub3(b,a), _sub3(c,a))
    n2=_dot3(n,n)
    n2>0 || return false
    abs(_dot3(n, _sub3(x,a)))<=atol*sqrt(n2) || return false
    # same-side barycentric via areas
    n1=_cross3(_sub3(b,a), _sub3(x,a))
    n2c=_cross3(_sub3(c,b), _sub3(x,b))
    n3=_cross3(_sub3(a,c), _sub3(x,c))
    return _dot3(n,n1)>=-atol*n2 && _dot3(n,n2c)>=-atol*n2 && _dot3(n,n3)>=-atol*n2
end

function _segment_triangle_hit(p, q, a, b, c; atol=1e-12)
    n=_cross3(_sub3(b,a), _sub3(c,a))
    d=_dot3(n, _sub3(q,p))
    abs(d)<=atol*hypot(n...) && return nothing
    t=_dot3(n, _sub3(a,p))/d
    (t>1e-9 && t<1-1e-9) || return nothing
    x=_add3(p, _scale3(_sub3(q,p), t))
    _point_in_triangle3(x,a,b,c; atol=atol) || return nothing
    return x
end

function _segment_face_hit(mesh::Mesh, p, q; atol=1e-12)
    best=nothing; bestt=Inf
    @inbounds for t in axes(mesh.tets,2)
        ids=(mesh.tets[1,t],mesh.tets[2,t],mesh.tets[3,t],mesh.tets[4,t])
        pts=(_pt3(mesh,ids[1]),_pt3(mesh,ids[2]),_pt3(mesh,ids[3]),_pt3(mesh,ids[4]))
        for (i,j,k) in ((2,3,4),(1,3,4),(1,2,4),(1,2,3))
            hit=_segment_triangle_hit(p,q,pts[i],pts[j],pts[k]; atol=atol)
            hit===nothing && continue
            _node_at3(mesh,hit; atol=1e-9)!=0 && continue
            tt=_dot3(_sub3(hit,p), _sub3(q,p))
            tt<bestt && (bestt=tt; best=hit)
        end
    end
    return best
end

function recover_segment3(mesh::Mesh, p, q; max_inserts::Integer=256)
    out,_=insert_steiner3(mesh,p)
    out,_=insert_steiner3(out,q)
    for _ in 1:max_inserts
        mesh_covers_segment3(out,p,q) && return out
        hit=_segment_face_hit(out,p,q)
        hit===nothing && throw(ErrorException(
            "recover_segment3: segment is not a tet-edge chain and no face crossing was found"))
        out,_=insert_steiner3(out,hit)
    end
    throw(ErrorException("recover_segment3: exceeded max_inserts=$max_inserts"))
end

function _triangle_area3(a,b,c)
    n=_cross3(_sub3(b,a), _sub3(c,a))
    return 0.5*hypot(n...)
end

function mesh_covers_triangle3(mesh::Mesh, a, b, c; rtol=1e-6)
    target=_triangle_area3(a,b,c)
    target>0 || return false
    seen=Set{NTuple{3,Int32}}()
    covered=0.0
    @inbounds for t in axes(mesh.tets,2)
        ids=(mesh.tets[1,t],mesh.tets[2,t],mesh.tets[3,t],mesh.tets[4,t])
        for (i,j,k) in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
            u,v,w=ids[i],ids[j],ids[k]
            s1,s2,s3=u,v,w
            s1>s2 && ((s1,s2)=(s2,s1)); s2>s3 && ((s2,s3)=(s3,s2)); s1>s2 && ((s1,s2)=(s2,s1))
            key=(s1,s2,s3)
            key in seen && continue
            pa=_pt3(mesh,s1); pb=_pt3(mesh,s2); pc=_pt3(mesh,s3)
            (_point_in_triangle3(pa,a,b,c) && _point_in_triangle3(pb,a,b,c) &&
             _point_in_triangle3(pc,a,b,c)) || continue
            push!(seen,key)
            covered+=_triangle_area3(pa,pb,pc)
        end
    end
    return abs(covered-target)<=rtol*target
end

function recover_triangle3(mesh::Mesh, a, b, c; max_inserts::Integer=256)
    out=recover_segment3(mesh,a,b)
    out=recover_segment3(out,b,c)
    out=recover_segment3(out,c,a)
    for _ in 1:max_inserts
        mesh_covers_triangle3(out,a,b,c) && return out
        # split the longest uncovered geometric edge of the sheet by inserting its midpoint
        mid=_add3(_scale3(_add3(a,b),0.5), (0.0,0.0,0.0))
        # try midpoints of the three sides, then the centroid
        candidates=(_scale3(_add3(a,b),0.5), _scale3(_add3(b,c),0.5),
                    _scale3(_add3(c,a),0.5), _scale3(_add3(_add3(a,b),c),1/3))
        inserted=false
        for x in candidates
            _node_at3(out,x; atol=1e-9)==0 || continue
            try
                out,_=insert_steiner3(out,x)
                inserted=true
                break
            catch err
                err isa InterruptException && rethrow()
                (err isa ArgumentError || err isa ErrorException) || rethrow()
            end
        end
        inserted || throw(ErrorException("recover_triangle3: could not insert a Steiner point on the sheet"))
    end
    throw(ErrorException("recover_triangle3: exceeded max_inserts=$max_inserts"))
end

"""
    insert_steiner3(mesh, point) -> (Mesh, node_index)

Insert one finite three-component point into a valid tetrahedral mesh. An exact
coordinate duplicate returns the unchanged mesh and its existing node index.
Otherwise the containing tetrahedron, face star, or edge star is split while
copying lower-dimensional cells unchanged and preserving tetrahedron tags.
Points outside the volume are rejected. The returned mesh is validated and
checked for total-volume conservation before return.
"""
function insert_steiner3(mesh::Mesh,p)
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "insert_steiner3: input mesh is invalid — "*join(diagnostic.messages,"; ")))
    size(mesh.tets,2)>0 || throw(ArgumentError(
        "insert_steiner3: input mesh has no tetrahedra"))
    point=_point3_input(p,"insert_steiner3","point")
    return _insert_steiner3(mesh,point)
end

function _insert_steiner3(mesh::Mesh,p::NTuple{3,Float64})
    existing=_node_at3(mesh,p)
    existing!=0 && return mesh, existing
    size(mesh.coords,2)<typemax(Int32) || throw(ArgumentError(
        "insert_steiner3: node count would exceed the Int32 indexing limit"))
    t,zeros,mask=_containing_tet(mesh,p)
    t==0 && throw(ArgumentError("insert_steiner3: point $p is not inside the volume"))
    out=if zeros==0
        _split_tet_one_to_four(mesh,t,p)
    elseif zeros==1
        _split_tets_on_face(mesh,t,p,mask)
    elseif zeros==2
        _split_tets_on_edge(mesh,t,p,mask)
    else
        throw(ArgumentError("insert_steiner3: point $p coincides with a tet vertex but was not found"))
    end
    diagnostic=validate(out)
    diagnostic.ok || throw(ErrorException(
        "insert_steiner3: insertion produced an invalid mesh — "*
        join(diagnostic.messages,"; ")))
    before=_tet_volume_sum3(mesh,"insert_steiner3 input")
    after=_tet_volume_sum3(out,"insert_steiner3 output")
    tolerance=max(1e-12*before,128eps(Float64)*max(before,after))
    abs(after-before)<=tolerance || throw(ErrorException(
        "insert_steiner3: insertion changed total volume from $before to $after"))
    return out,size(out.coords,2)
end

function _tet_volume_sum3(mesh::Mesh,caller::AbstractString)
    total=0.0
    correction=0.0
    @inbounds for t in axes(mesh.tets,2)
        value=abs(tet_signed_volume(
            _pt3(mesh,mesh.tets[1,t]),_pt3(mesh,mesh.tets[2,t]),
            _pt3(mesh,mesh.tets[3,t]),_pt3(mesh,mesh.tets[4,t])))
        y=value-correction
        next=total+y
        correction=(next-total)-y
        total=next
    end
    (isfinite(total) && total>0) || throw(ArgumentError(
        "$caller: tetrahedron volume sum is not finite and positive"))
    return total
end

function _split_tet_one_to_four(mesh::Mesh, t::Int, p)
    nn=size(mesh.coords,2)+1
    nn<=typemax(Int32) || throw(ArgumentError(
        "insert_steiner3: node count exceeds the Int32 indexing limit"))
    coords=hcat(mesh.coords, [p[1],p[2],p[3]])
    v=Int32(nn)
    a,b,c,d=mesh.tets[1,t],mesh.tets[2,t],mesh.tets[3,t],mesh.tets[4,t]
    nt=size(mesh.tets,2)
    nt<=typemax(Int32)-3 || throw(ArgumentError(
        "insert_steiner3: tetrahedron count exceeds the Int32 topology limit"))
    tets=Matrix{Int32}(undef,4,nt+3)
    @inbounds for j in 1:nt
        if j==t
            tets[1,j]=v; tets[2,j]=b; tets[3,j]=c; tets[4,j]=d
        else
            for i in 1:4; tets[i,j]=mesh.tets[i,j]; end
        end
    end
    tets[1,nt+1]=a; tets[2,nt+1]=v; tets[3,nt+1]=c; tets[4,nt+1]=d
    tets[1,nt+2]=a; tets[2,nt+2]=b; tets[3,nt+2]=v; tets[4,nt+2]=d
    tets[1,nt+3]=a; tets[2,nt+3]=b; tets[3,nt+3]=c; tets[4,nt+3]=v
    function _pos!(col)
        pa=(coords[1,tets[1,col]],coords[2,tets[1,col]],coords[3,tets[1,col]])
        pb=(coords[1,tets[2,col]],coords[2,tets[2,col]],coords[3,tets[2,col]])
        pc=(coords[1,tets[3,col]],coords[2,tets[3,col]],coords[3,tets[3,col]])
        pd=(coords[1,tets[4,col]],coords[2,tets[4,col]],coords[3,tets[4,col]])
        if tet_signed_volume(pa,pb,pc,pd)<0
            tets[2,col],tets[3,col]=tets[3,col],tets[2,col]
        end
    end
    _pos!(t); _pos!(nt+1); _pos!(nt+2); _pos!(nt+3)
    tags=if isempty(mesh.tet_tag)
        Int32[]
    else
        tag=mesh.tet_tag[t]
        vcat(mesh.tet_tag, Int32[tag,tag,tag])
    end
    return Mesh(coords; segs=copy(mesh.segs), tris=copy(mesh.tris), tets=tets,
                seg_tag=copy(mesh.seg_tag), tri_tag=copy(mesh.tri_tag), tet_tag=tags)
end

function _orient_tet!(tets, coords, col)
    pa=(coords[1,tets[1,col]],coords[2,tets[1,col]],coords[3,tets[1,col]])
    pb=(coords[1,tets[2,col]],coords[2,tets[2,col]],coords[3,tets[2,col]])
    pc=(coords[1,tets[3,col]],coords[2,tets[3,col]],coords[3,tets[3,col]])
    pd=(coords[1,tets[4,col]],coords[2,tets[4,col]],coords[3,tets[4,col]])
    if tet_signed_volume(pa,pb,pc,pd)<0
        tets[2,col],tets[3,col]=tets[3,col],tets[2,col]
    end
    return nothing
end

function _rebuild_tets(mesh::Mesh, drop::Set{Int}, news::Vector{NTuple{4,Int32}},
                       new_tags::Vector{Int32}, p)
    nn=size(mesh.coords,2)+1
    nn<=typemax(Int32) || throw(ArgumentError(
        "insert_steiner3: node count exceeds the Int32 indexing limit"))
    length(new_tags)==length(news) || throw(ErrorException(
        "insert_steiner3: replacement connectivity/tag counts differ"))
    nkeep=size(mesh.tets,2)-length(drop)
    nkeep>=0 || throw(ErrorException(
        "insert_steiner3: replacement removes more tetrahedra than exist"))
    nkeep<=typemax(Int32)-length(news) || throw(ArgumentError(
        "insert_steiner3: tetrahedron count exceeds the Int32 topology limit"))
    nout=nkeep+length(news)
    coords=hcat(mesh.coords, Float64[p[1],p[2],p[3]])
    keep=Int[j for j in axes(mesh.tets,2) if !(j in drop)]
    length(keep)==nkeep || throw(ErrorException(
        "insert_steiner3: replacement keep count mismatch"))
    tets=Matrix{Int32}(undef,4,nout)
    tags=isempty(mesh.tet_tag) ? Int32[] : Vector{Int32}(undef,nout)
    @inbounds for (k,j) in enumerate(keep)
        for i in 1:4; tets[i,k]=mesh.tets[i,j]; end
        isempty(tags) || (tags[k]=mesh.tet_tag[j])
    end
    base=length(keep)
    @inbounds for (k,tet) in enumerate(news)
        col=base+k
        tets[1,col]=tet[1]; tets[2,col]=tet[2]; tets[3,col]=tet[3]; tets[4,col]=tet[4]
        _orient_tet!(tets,coords,col)
        isempty(tags) || (tags[col]=new_tags[k])
    end
    return Mesh(coords; segs=copy(mesh.segs), tris=copy(mesh.tris), tets=tets,
                seg_tag=copy(mesh.seg_tag), tri_tag=copy(mesh.tri_tag), tet_tag=tags)
end

function _split_tets_on_face(mesh::Mesh, t::Int, p, mask)
    ids=(mesh.tets[1,t],mesh.tets[2,t],mesh.tets[3,t],mesh.tets[4,t])
    # mask[k] true ⇒ p on the face opposite vertex k
    opp = mask[1] ? 1 : mask[2] ? 2 : mask[3] ? 3 : 4
    face=Int32[ids[i] for i in 1:4 if i!=opp]
    face_set=Set{Int32}(face)
    drop=Set{Int}()
    news=NTuple{4,Int32}[]
    ntags=Int32[]
    v=Int32(size(mesh.coords,2)+1)
    @inbounds for j in axes(mesh.tets,2)
        verts=(mesh.tets[1,j],mesh.tets[2,j],mesh.tets[3,j],mesh.tets[4,j])
        count(in(face_set), verts)==3 || continue
        push!(drop,j)
        w=verts[1]
        for u in verts
            u in face_set || (w=u)
        end
        a,b,c=face[1],face[2],face[3]
        tag=isempty(mesh.tet_tag) ? Int32(0) : mesh.tet_tag[j]
        push!(news,(v,b,c,w)); push!(ntags,tag)
        push!(news,(a,v,c,w)); push!(ntags,tag)
        push!(news,(a,b,v,w)); push!(ntags,tag)
    end
    isempty(drop) && throw(ErrorException("insert_steiner3: face split found no incident tetrahedra"))
    return _rebuild_tets(mesh, drop, news, ntags, p)
end

function _split_tets_on_edge(mesh::Mesh, t::Int, p, mask)
    ids=(mesh.tets[1,t],mesh.tets[2,t],mesh.tets[3,t],mesh.tets[4,t])
    excluded=Int[i for i in 1:4 if mask[i]]
    length(excluded)==2 || throw(ErrorException("insert_steiner3: edge split expected two zero barycentric coordinates"))
    # edge is the two vertices that are NOT the excluded (opposite) vertices
    edge=Int32[ids[i] for i in 1:4 if !(i in excluded)]
    length(edge)==2 || throw(ErrorException("insert_steiner3: could not identify the supporting edge"))
    u,w=edge[1],edge[2]
    drop=Set{Int}()
    news=NTuple{4,Int32}[]
    ntags=Int32[]
    v=Int32(size(mesh.coords,2)+1)
    @inbounds for j in axes(mesh.tets,2)
        verts=(mesh.tets[1,j],mesh.tets[2,j],mesh.tets[3,j],mesh.tets[4,j])
        (u in verts && w in verts) || continue
        push!(drop,j)
        rest=Int32[x for x in verts if x!=u && x!=w]
        length(rest)==2 || continue
        c,d=rest[1],rest[2]
        tag=isempty(mesh.tet_tag) ? Int32(0) : mesh.tet_tag[j]
        push!(news,(u,v,c,d)); push!(ntags,tag)
        push!(news,(v,w,c,d)); push!(ntags,tag)
    end
    isempty(drop) && throw(ErrorException("insert_steiner3: edge split found no incident tetrahedra"))
    return _rebuild_tets(mesh, drop, news, ntags, p)
end

"""
    tetrahedralize(surface; rng_seed=1, optimize=false, check=true,
                   interior_points=nothing,
                   max_interior_points=1_000_000) -> Mesh

Fill the volume enclosed by a triangulated `surface` with tetrahedra. The input
must be a structurally valid triangle surface with no unreferenced nodes or
tetrahedra. The fast path Delaunay-tetrahedralizes its vertices and retains tets
whose centroids lie inside. Before a checked result can escape, the PLC gate used
by boundary recovery certifies positive volume, manifold vertex links, creases,
and the complete input boundary. A failed restriction is rebuilt by the exact
conforming-Delaunay backend, with bounded multi-seed [`recover_boundary`](@ref)
as a final alternative. If neither recovery succeeds, an explicit blocker is
thrown.

`interior_points`, when supplied, must be an iterable of finite three-component
real tuples or vectors. Points are inserted after the certified fill, without
moving the boundary; exact duplicates reuse the existing node. Their count is
bounded by `max_interior_points` and the Int32 mesh-indexing limit.

`check=false` is an expert bypass for the closed-PLC certificate and recovery;
basic input structure and all interior-point/output safety checks still apply.
The caller owns the resulting boundary-conformity assessment.
"""
function tetrahedralize(surface::Mesh; rng_seed::Integer=1, optimize::Bool=false,
                        check::Bool=true, interior_points=nothing,
                        max_interior_points::Integer=1_000_000)
    check ? _require_surface3(surface,"tetrahedralize") :
            _require_surface_input3(surface,"tetrahedralize")
    interior_limit=_bounded_count3(max_interior_points,"tetrahedralize",
                                   "max_interior_points")
    nn = size(surface.coords, 2)
    extra=_interior_points3(interior_points,interior_limit,nn,"tetrahedralize")
    seed=_seed3(rng_seed,"tetrahedralize")
    xs = Vector{Float64}(undef, nn); ys = similar(xs); zs = similar(xs)
    @inbounds for i in 1:nn; xs[i]=surface.coords[1,i]; ys[i]=surface.coords[2,i]; zs[i]=surface.coords[3,i]; end
    T = delaunay3d(xs, ys, zs; rng_seed=seed)
    optimize && optimize_flips!(T; passes=4)     # sliver-reducing flips (safe, volume-preserving)
    keep = _classify_by_centroid(T, surface)
    m = to_mesh3(T; keep=keep)
    check || return insert_interior_points(m, extra)
    ok, reason = _certify_surface_fill(surface, m)
    if ok
        return _finish_interior(m, extra)
    end
    # A kernel fan is the cheapest conforming repair for star-shaped PLCs and avoids
    # repeating equivalent cospherical Delaunay builds (notably polygonal cylinders).
    # It is accepted only through the full PLC certificate; non-star-shaped inputs
    # fall through to general boundary recovery.
    Px, Py, Pz, facets = _rb_dedup_surface(surface)
    if length(Px) >= 4 && !isempty(facets)
        fan = _rb_fan_steiner(Px, Py, Pz, facets)
        if fan !== nothing
            fanok, _ = _certify_surface_fill(surface, fan)
            fanok && return _finish_interior(fan, extra)
        end
    end
    exact_reason = ""
    if applicable(_recover_boundary_exact, surface)
        try
            return _finish_interior(_recover_boundary_exact(surface), extra)
        catch err
            err isa InterruptException && rethrow()
            (err isa ArgumentError || err isa ErrorException) || rethrow()
            exact_reason = sprint(showerror, err)
        end
    else
        exact_reason = "exact recovery backend is unavailable"
    end
    # The exact engine is deterministic and normally settles hard non-star-shaped
    # PLCs directly.  Bound the alternative randomized search here so one public
    # call cannot spend its full 64-seed expert-API budget after exact recovery has
    # already reported a blocker.
    recovery_reason = ""
    try
        return _finish_interior(recover_boundary(surface; rng_seed=seed, max_seeds=8), extra)
    catch err
        err isa InterruptException && rethrow()
        (err isa ArgumentError || err isa ErrorException) || rethrow()
        recovery_reason = sprint(showerror, err)
    end
    throw(ErrorException("tetrahedralize: restricted Delaunay failed the PLC gate ($reason); " *
                         "exact conforming-Delaunay recovery also failed: $exact_reason; " *
                         "boundary recovery also failed: $recovery_reason"))
end

"""
    mesh_box(x0, x1, y0, y1, z0, z1; hmax) -> Mesh

Size-controlled tetrahedral mesh of the axis-aligned box `[x0,x1]×[y0,y1]×[z0,z1]`
with **guaranteed maximum edge length ≤ `hmax`**. Unlike the Delaunay path
([`tetrahedralize`](@ref)), which inherits the input surface's resolution and adds
no interior points, this places a regular grid at spacing `hmax/√3` per axis and
tetrahedralizes each grid cube by the **Kuhn/Freudenthal subdivision** — six
path-tetrahedra sharing the cube's main diagonal.

The subdivision is *explicit connectivity* (no Delaunay of the degenerate lattice,
whose cospherical ties break validity — see `validation/stage4_size_refinement/`),
so the result is provably **valid** (all tets positively oriented, no slivers),
**watertight/conforming** (every cube uses the same diagonal ⇒ shared faces match),
**exact in volume** (`∑vol == box volume`), and **terminating** (finite grid).
The cube main diagonal `= (hmax/√3)·√3 = hmax` bounds the longest edge, so
`maxedge ≤ hmax` holds by construction. This is the Stage-4 size-control primitive
for box regions (e.g. the enclosure case/air cavities).
"""
function mesh_box(x0::Real, x1::Real, y0::Real, y1::Real, z0::Real, z1::Real; hmax::Real)
    x0=_finite3(x0,"mesh_box","x0"); x1=_finite3(x1,"mesh_box","x1")
    y0=_finite3(y0,"mesh_box","y0"); y1=_finite3(y1,"mesh_box","y1")
    z0=_finite3(z0,"mesh_box","z0"); z1=_finite3(z1,"mesh_box","z1")
    hm=_finite3(hmax,"mesh_box","hmax")
    Lx = x1-x0; Ly = y1-y0; Lz = z1-z0
    (isfinite(Lx) && isfinite(Ly) && isfinite(Lz) && Lx > 0 && Ly > 0 && Lz > 0) ||
        throw(ArgumentError("mesh_box: box must have positive extent in every axis (got Lx=$Lx, Ly=$Ly, Lz=$Lz)"))
    hm > 0 || throw(ArgumentError("mesh_box: hmax must be positive (got $hmax)"))
    a  = hm / sqrt(3.0)                                     # cube edge; main diagonal = a√3 = hmax
    (isfinite(a) && a > 0) || throw(ArgumentError("mesh_box: hmax is below Float64 spacing resolution"))
    nx=_ceil_count3(Lx/a,"mesh_box","x intervals")
    ny=_ceil_count3(Ly/a,"mesh_box","y intervals")
    nz=_ceil_count3(Lz/a,"mesh_box","z intervals")
    hx = Lx/nx; hy = Ly/ny; hz = Lz/nz
    nid(i,j,k) = ((k*(ny+1) + j)*(nx+1) + i) + 1            # 1-based grid-node id
    npx=_checked_add3("mesh_box","node",nx,1)
    npy=_checked_add3("mesh_box","node",ny,1)
    npz=_checked_add3("mesh_box","node",nz,1)
    npts=_checked_mul3("mesh_box","node",npx,npy,npz)
    npts <= typemax(Int32) ||
        throw(ArgumentError("mesh_box: $npts nodes exceed the Int32 indexing limit"))
    ntotal=_checked_mul3("mesh_box","tetrahedron",6,nx,ny,nz)
    ntotal <= typemax(Int32) ||
        throw(ArgumentError("mesh_box: $ntotal tetrahedra exceed the Int32 working limit"))
    C = Matrix{Float64}(undef, 3, npts)
    @inbounds for k in 0:nz, j in 0:ny, i in 0:nx
        n = nid(i,j,k); C[1,n]=x0+i*hx; C[2,n]=y0+j*hy; C[3,n]=z0+k*hz
    end
    # six path-tetrahedra per cube: 000 → step1 → step2 → 111 (the 3! monotone
    # lattice paths from the low corner to the high corner). Emitted CCW-positive.
    paths = (((1,0,0),(1,1,0)), ((1,0,0),(1,0,1)), ((0,1,0),(1,1,0)),
             ((0,1,0),(0,1,1)), ((0,0,1),(1,0,1)), ((0,0,1),(0,1,1)))
    tets = Matrix{Int32}(undef, 4, ntotal); nt = 0
    @inbounds for k in 0:nz-1, j in 0:ny-1, i in 0:nx-1
        c(di,dj,dk) = nid(i+di, j+dj, k+dk)
        v0 = c(0,0,0); v7 = c(1,1,1)
        for (s1,s2) in paths
            v1 = c(s1...); v2 = c(s2...)
            # orient so the signed volume is positive (swap the last two if not)
            p0=(C[1,v0],C[2,v0],C[3,v0]); p1=(C[1,v1],C[2,v1],C[3,v1])
            p2=(C[1,v2],C[2,v2],C[3,v2]); p3=(C[1,v7],C[2,v7],C[3,v7])
            nt += 1; tets[1,nt]=v0; tets[2,nt]=v1
            if _signed_vol6(p0,p1,p2,p3) >= 0
                tets[3,nt]=v2; tets[4,nt]=v7
            else
                tets[3,nt]=v7; tets[4,nt]=v2
            end
        end
    end
    out=Mesh(C; tets=tets)
    d=validate(out)
    d.ok || throw(ErrorException("mesh_box: constructed mesh failed validation — " * join(d.messages,"; ")))
    return out
end

# 6× signed tet volume (sign only matters here); positive ⇔ (p1-p0,p2-p0,p3-p0) right-handed
@inline _signed_vol6(p0,p1,p2,p3) = -orient3(p0,p1,p2,p3)

"""
    BoxRegion(x0,x1,y0,y1,z0,z1, tag)

An axis-aligned box `[x0,x1]×[y0,y1]×[z0,z1]` carrying an integer region `tag`.
`tag == 0` is reserved for a **void** (a hole): grid cells classified into a void
region are dropped, so a later void box carves a CSG difference out of an earlier
solid one. Used by [`mesh_box_regions`](@ref).
"""
struct BoxRegion
    x0::Float64; x1::Float64; y0::Float64; y1::Float64; z0::Float64; z1::Float64
    tag::Int
    function BoxRegion(x0::Real,x1::Real,y0::Real,y1::Real,z0::Real,z1::Real,tag::Integer)
        a0=_finite3(x0,"BoxRegion","x0");a1=_finite3(x1,"BoxRegion","x1")
        b0=_finite3(y0,"BoxRegion","y0");b1=_finite3(y1,"BoxRegion","y1")
        c0=_finite3(z0,"BoxRegion","z0");c1=_finite3(z1,"BoxRegion","z1")
        (a0<a1 && b0<b1 && c0<c1) ||
            throw(ArgumentError("BoxRegion: every axis must have finite positive extent"))
        0 <= tag <= typemax(Int32) ||
            throw(ArgumentError("BoxRegion: tag must be in 0:$(typemax(Int32)) (got $tag)"))
        new(a0,a1,b0,b1,c0,c1,Int(tag))
    end
end
@inline _br_contains(b::BoxRegion, x,y,z) =
    (b.x0<=x<=b.x1) && (b.y0<=y<=b.y1) && (b.z0<=z<=b.z1)

# 1-D grid on one axis: every region face coordinate is a breakpoint, and each gap
# is subdivided so no sub-interval exceeds `a` (⇒ every cell main diagonal ≤ a√3).
function _region_axis_grid(faces::Vector{Float64}, a::Float64)
    isempty(faces) && throw(ArgumentError("mesh_box_regions: an axis has no region faces"))
    b = sort!(unique(faces)); g = Float64[b[1]]
    @inbounds for i in 1:length(b)-1
        lo=b[i]; hi=b[i+1]; L=hi-lo
        n = _ceil_count3(max(0.0,L/a-1e-9),"mesh_box_regions","axis intervals"); h = L/n
        _checked_add3("mesh_box_regions","axis grid",length(g),n)
        for k in 1:n; push!(g, lo + k*h); end
    end
    return g
end

"""
    mesh_box_regions(regions::AbstractVector{BoxRegion}; hmax) -> Mesh

Conforming, size-controlled, **multi-region** tetrahedral mesh of an assembly of
axis-aligned boxes — the native-CSG counterpart of [`mesh_box`](@ref). Every
region face becomes a plane of one shared global grid, each gap subdivided so
`maxedge ≤ hmax`; each grid cell is split by the same Kuhn/Freudenthal diagonal, so
tets on either side of a region interface **share faces exactly** (conforming, no
T-junctions). A cell is tagged by the **last** (highest-priority) region whose box
contains its centre; cells in a `tag==0` (void) region or in no region are dropped.

This realizes, for axis-aligned box primitives, **union** (several solid boxes),
**difference** (a `tag==0` box carves a hole), **nesting** (priority order), and
multi-material partitions (distinct tags) — all provably valid (positive tets),
watertight/conforming, with exact per-region volumes, at a guaranteed edge size.
Handles non-convex domains (a hollow shell is `box − void`). `tet_tag` holds the
region tag of each tet. Errors if no cell is kept (empty domain).
"""
function mesh_box_regions(regions::AbstractVector{BoxRegion}; hmax::Real)
    hm=_finite3(hmax,"mesh_box_regions","hmax")
    hm > 0 || throw(ArgumentError("mesh_box_regions: hmax must be positive (got $hmax)"))
    isempty(regions) && throw(ArgumentError("mesh_box_regions: no regions"))
    @inbounds for r in regions
        (r.x1>r.x0 && r.y1>r.y0 && r.z1>r.z0) ||
            throw(ArgumentError("mesh_box_regions: every region needs positive extent (bad box tag=$(r.tag))"))
    end
    a = hm/sqrt(3.0)
    (isfinite(a) && a > 0) ||
        throw(ArgumentError("mesh_box_regions: hmax is below Float64 spacing resolution"))
    fx=Float64[]; fy=Float64[]; fz=Float64[]
    @inbounds for r in regions
        push!(fx,r.x0); push!(fx,r.x1); push!(fy,r.y0); push!(fy,r.y1); push!(fz,r.z0); push!(fz,r.z1)
    end
    X=_region_axis_grid(fx,a); Y=_region_axis_grid(fy,a); Z=_region_axis_grid(fz,a)
    nx=length(X)-1; ny=length(Y)-1; nz=length(Z)-1
    LX=length(X); LY=length(Y)
    globalnodes=_checked_mul3("mesh_box_regions","global grid node",LX,LY,length(Z))
    globalnodes <= typemax(Int32) ||
        throw(ArgumentError("mesh_box_regions: $globalnodes global grid nodes exceed the Int32 indexing limit"))
    ncells=_checked_mul3("mesh_box_regions","grid cell",nx,ny,nz)
    ntmax=_checked_mul3("mesh_box_regions","tetrahedron",6,ncells)
    ntmax <= typemax(Int32) ||
        throw(ArgumentError("mesh_box_regions: the grid can produce $ntmax tetrahedra, exceeding the Int32 working limit"))
    gnid(i,j,k) = ((k*LY + j)*LX + i) + 1
    tetv = NTuple{4,Int32}[]; tett = Int32[]; used = Set{Int32}()
    @inbounds for k in 0:nz-1, j in 0:ny-1, i in 0:nx-1
        cx=0.5*(X[i+1]+X[i+2]); cy=0.5*(Y[j+1]+Y[j+2]); cz=0.5*(Z[k+1]+Z[k+2])
        tag = 0
        for r in regions
            _br_contains(r,cx,cy,cz) && (tag = r.tag)     # last (highest-priority) wins
        end
        tag == 0 && continue                              # void / outside → drop cell
        cn(di,dj,dk) = gnid(i+di, j+dj, k+dk)
        pc(di,dj,dk) = (X[i+di+1], Y[j+dj+1], Z[k+dk+1])
        v0=cn(0,0,0); v7=cn(1,1,1); p0=pc(0,0,0); p7=pc(1,1,1)
        for (s1,s2) in ((( 1,0,0),(1,1,0)), ((1,0,0),(1,0,1)), ((0,1,0),(1,1,0)),
                        (( 0,1,0),(0,1,1)), ((0,0,1),(1,0,1)), ((0,0,1),(0,1,1)))
            v1=cn(s1...); v2=cn(s2...); p1=pc(s1...); p2=pc(s2...)
            if _signed_vol6(p0,p1,p2,p7) >= 0
                push!(tetv,(v0,v1,v2,v7))
            else
                push!(tetv,(v0,v1,v7,v2))
            end
            push!(tett, Int32(tag))
            push!(used,v0); push!(used,v1); push!(used,v2); push!(used,v7)
        end
    end
    isempty(tetv) && throw(ArgumentError("mesh_box_regions: no cells kept (empty domain — all void/outside)"))
    usedv = sort!(collect(used)); remap = Dict{Int32,Int32}()
    @inbounds for (n,v) in enumerate(usedv); remap[v]=Int32(n); end
    C = Matrix{Float64}(undef, 3, length(usedv))
    @inbounds for (n,v) in enumerate(usedv)
        v0=Int(v)-1; i=v0 % LX; j=(v0 ÷ LX) % LY; k=v0 ÷ (LX*LY)
        C[1,n]=X[i+1]; C[2,n]=Y[j+1]; C[3,n]=Z[k+1]
    end
    Tm = Matrix{Int32}(undef, 4, length(tetv)); tags = Vector{Int32}(undef, length(tetv))
    @inbounds for (t,f) in enumerate(tetv)
        Tm[1,t]=remap[f[1]]; Tm[2,t]=remap[f[2]]; Tm[3,t]=remap[f[3]]; Tm[4,t]=remap[f[4]]; tags[t]=tett[t]
    end
    out=Mesh(C; tets=Tm, tet_tag=tags)
    d=validate(out)
    d.ok || throw(ErrorException("mesh_box_regions: constructed mesh failed validation — " * join(d.messages,"; ")))
    return out
end

"""
    mesh_cylinder(center, axis, radius, height; hmax) -> Mesh

Size-controlled tetrahedral mesh of a cylinder (axis start `center`, direction
`axis`, `radius`, `height`) with **guaranteed maximum edge length ≤ `hmax`** — the
uniform-size counterpart of [`mesh_box`](@ref) for the *cylindrical* primitive, and
the robust route for the exact-cospherical case where a Delaunay fill degenerates.

Explicit structured `(r,θ,z)` grid (spacings `≤ hmax/√3`): each index cube is split
by the Kuhn/Freudenthal subdivision, the axis column (`r=0`) collapses to one node
per level (its degenerate tets are dropped), and the angular seam wraps — so the
mesh is **valid** (all-positive), **watertight** (boundary χ=2), **exact in volume**
(the faceted-`nθ`-gon prism volume), and `maxedge ≤ hmax` by construction, with **no
Delaunay** (hence no cospherical degeneracy). *Quality note:* the collapsed axis
gives thinner tets there (min dihedral degrades toward the axis as `nθ` grows — the
inherent singularity of structured cylinder meshes); the size/validity guarantees
hold regardless.
"""
function mesh_cylinder(center, axis, radius::Real, height::Real; hmax::Real)
    length(center)>=3 || throw(ArgumentError("mesh_cylinder: center needs three coordinates"))
    length(axis)>=3 || throw(ArgumentError("mesh_cylinder: axis needs three coordinates"))
    R=_finite3(radius,"mesh_cylinder","radius"); H=_finite3(height,"mesh_cylinder","height")
    hm=_finite3(hmax,"mesh_cylinder","hmax")
    cx=_finite3(center[1],"mesh_cylinder","center[1]")
    cy=_finite3(center[2],"mesh_cylinder","center[2]")
    cz=_finite3(center[3],"mesh_cylinder","center[3]")
    av=(_finite3(axis[1],"mesh_cylinder","axis[1]"),
        _finite3(axis[2],"mesh_cylinder","axis[2]"),
        _finite3(axis[3],"mesh_cylinder","axis[3]"))
    (R > 0 && H > 0) || throw(ArgumentError("mesh_cylinder: radius and height must be positive (got R=$R, H=$H)"))
    hm > 0 || throw(ArgumentError("mesh_cylinder: hmax must be positive (got $hmax)"))
    hypot(av...) > 0 ||
        throw(ArgumentError("mesh_cylinder: axis must be a nonzero direction vector (got $((axis[1],axis[2],axis[3])))"))
    a = hm/sqrt(3.0)
    (isfinite(a) && a > 0) || throw(ArgumentError("mesh_cylinder: hmax is below Float64 spacing resolution"))
    nr=_ceil_count3(R/a,"mesh_cylinder","radial intervals")
    nz=_ceil_count3(H/a,"mesh_cylinder","axial intervals")
    nθ=_ceil_count3(2π*R/a,"mesh_cylinder","angular intervals";minimum=3)
    ringnodes=_checked_mul3("mesh_cylinder","ring node",nr,nθ)
    perlevel=_checked_add3("mesh_cylinder","nodes per level",ringnodes,1)
    levels=_checked_add3("mesh_cylinder","axial level",nz,1)
    nmax=_checked_mul3("mesh_cylinder","node",perlevel,levels)
    nmax <= typemax(Int32) ||
        throw(ArgumentError("mesh_cylinder: $nmax nodes exceed the Int32 indexing limit"))
    ntcandidate=_checked_mul3("mesh_cylinder","candidate tetrahedron",6,nr,nθ,nz)
    ntcandidate <= typemax(Int32) ||
        throw(ArgumentError("mesh_cylinder: $ntcandidate candidate tetrahedra exceed the Int32 working limit"))
    ez=_unitn(av)
    axc = abs(ez[1])<=abs(ez[2]) ? (abs(ez[1])<=abs(ez[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
                                   (abs(ez[2])<=abs(ez[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    d=axc[1]*ez[1]+axc[2]*ez[2]+axc[3]*ez[3]
    ex=_unitn((axc[1]-d*ez[1], axc[2]-d*ez[2], axc[3]-d*ez[3]))
    ey=(ez[2]*ex[3]-ez[3]*ex[2], ez[3]*ex[1]-ez[1]*ex[3], ez[1]*ex[2]-ez[2]*ex[1])
    @inline function pos(i,k,j)
        r=R*i/nr; θ=2π*(mod(k,nθ))/nθ; z=H*j/nz; rc=r*cos(θ); rs=r*sin(θ)
        (cx+rc*ex[1]+rs*ey[1]+z*ez[1], cy+rc*ex[2]+rs*ey[2]+z*ez[2], cz+rc*ex[3]+rs*ey[3]+z*ez[3])
    end
    nid=Dict{NTuple{3,Int},Int32}(); coords=NTuple{3,Float64}[]
    sizehint!(nid,nmax); sizehint!(coords,nmax)
    @inline function id(i,k,j)
        key = i==0 ? (0,0,j) : (i, Int(mod(k,nθ)), j)
        get!(nid,key) do
            p=pos(i,k,j)
            (isfinite(p[1])&&isfinite(p[2])&&isfinite(p[3])) ||
                throw(ArgumentError("mesh_cylinder: generated coordinate overflowed Float64"))
            push!(coords,p); Int32(length(coords))
        end
    end
    tets=NTuple{4,Int32}[]
    sizehint!(tets,ntcandidate)
    @inbounds for j in 0:nz-1, k in 0:nθ-1, i in 0:nr-1
        crn(a1,a2,a3)=id(i+a1, k+a2, j+a3)
        v0=crn(0,0,0); v7=crn(1,1,1)
        for (s1,s2) in ((( 1,0,0),(1,1,0)), ((1,0,0),(1,0,1)), ((0,1,0),(1,1,0)),
                        (( 0,1,0),(0,1,1)), ((0,0,1),(1,0,1)), ((0,0,1),(0,1,1)))
            v1=crn(s1...); v2=crn(s2...)
            (v0==v1||v0==v2||v0==v7||v1==v2||v1==v7||v2==v7) && continue   # drop axis-collapsed degenerates
            p0=coords[v0]; p1=coords[v1]; p2=coords[v2]; p3=coords[v7]
            sv=_signed_vol6(p0,p1,p2,p3)
            sv == 0.0 && continue                                           # drop only exactly-flat tets:
            # the vertex-equality check above already removes every genuine axis-collapse
            # degeneracy (distinct (i,k,j) keys ⇒ distinct coords for r>0), so the remaining
            # tets are all non-degenerate. An ABSOLUTE volume tolerance here would be
            # scale-dependent — at small geometric scale it silently drops legitimate thin
            # near-axis tets (6·vol ~ a³·Δθ), leaving an axial void (χ=0, wrong volume). An
            # exact-zero test keeps the watertight/exact-volume guarantee at every scale.
            push!(tets, sv>=0 ? (v0,v1,v2,v7) : (v0,v1,v7,v2))
        end
    end
    isempty(tets) && throw(ErrorException("mesh_cylinder: produced no tets (degenerate parameters)"))
    C=Matrix{Float64}(undef,3,length(coords))
    @inbounds for (n,p) in enumerate(coords); C[1,n]=p[1]; C[2,n]=p[2]; C[3,n]=p[3]; end
    Tm=Matrix{Int32}(undef,4,length(tets))
    @inbounds for (t,f) in enumerate(tets); Tm[1,t]=f[1]; Tm[2,t]=f[2]; Tm[3,t]=f[3]; Tm[4,t]=f[4]; end
    out=Mesh(C; tets=Tm)
    d=validate(out)
    d.ok || throw(ErrorException("mesh_cylinder: constructed mesh failed validation — " * join(d.messages,"; ")))
    @inbounds for t in axes(Tm,2), (i,j) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
        u=Tm[i,t];v=Tm[j,t]
        hypot(C[1,u]-C[1,v],C[2,u]-C[2,v],C[3,u]-C[3,v]) <= hm*(1+16eps(Float64)) ||
            throw(ErrorException("mesh_cylinder: postcondition failed: an output edge exceeds hmax=$hm"))
    end
    return out
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
    seed=_seed3(rng_seed,"tetrahedralize_multi")
    length(surfaces)<=typemax(Int32) ||
        throw(ArgumentError("tetrahedralize_multi: region count exceeds Int32 tags"))
    coords_cols = Vector{NTuple{3,Float64}}()
    tetcols = Vector{NTuple{4,Int32}}()
    tags = Int32[]
    for (r, surf) in enumerate(surfaces)
        m = tetrahedralize(surf; rng_seed=seed)
        length(coords_cols) <= typemax(Int32)-size(m.coords,2) ||
            throw(ArgumentError("tetrahedralize_multi: combined node count exceeds Int32"))
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
    tetrahedralize_conforming(surfaces; rng_seed=1) -> Mesh

Mesh several region boundary `surfaces` into ONE **conforming** partition. All region
vertices are Delaunay-tetrahedralized *together* on **exact coordinates** — the SoS
predicates break coplanar/cospherical degeneracies symbolically, with **no coordinate
jitter** — then each tet is tagged (`tet_tag`, 1-based) by the region whose surface
contains its centroid (ray-cast point-in-polyhedron). Because it is a single
triangulation, an interface between two regions is a shared set of tet faces —
**conforming across regions**, unlike [`tetrahedralize_multi`](@ref) which fills each
region independently and leaves interfaces non-matching. Tets outside every region
(convex-hull filler) are dropped; it errors if any region receives no tet (the
anti-false-positive contract — a valid partition fills every region).

This is the piece the coordinate-jitter of [`delaunay3d`](@ref)`(perturb=true)` made
impossible: jitter of ~1e-8 moves shared interface vertices off their plane, so no tet
face lies on the interface. Exact coordinates + correct SoS keep the interface a union
of Delaunay faces. Suited to partitions with planar/axis-aligned interfaces (the
enclosure air/case pattern); an interface that is *not* a union of Delaunay faces would
still need constrained boundary recovery.
"""
function tetrahedralize_conforming(surfaces::AbstractVector{Mesh}; rng_seed::Integer=1)
    isempty(surfaces) && return Mesh(Matrix{Float64}(undef, 3, 0))
    seed=_seed3(rng_seed,"tetrahedralize_conforming")
    length(surfaces)<=typemax(Int32) ||
        throw(ArgumentError("tetrahedralize_conforming: region count exceeds Int32 tags"))
    for (r,s) in enumerate(surfaces);_require_surface3(s,"tetrahedralize_conforming region $r");end
    # collect every region vertex once (shared interface vertices are deduplicated)
    seen = Set{NTuple{3,Float64}}(); V = NTuple{3,Float64}[]
    for s in surfaces
        @inbounds for i in 1:size(s.coords, 2)
            p = (s.coords[1,i]==0 ? 0.0 : s.coords[1,i],
                 s.coords[2,i]==0 ? 0.0 : s.coords[2,i],
                 s.coords[3,i]==0 ? 0.0 : s.coords[3,i])
            (p in seen) || (push!(seen, p); push!(V, p))
        end
    end
    n = length(V)
    n<=typemax(Int32) || throw(ArgumentError("tetrahedralize_conforming: $n distinct nodes exceed Int32"))
    xs = Vector{Float64}(undef, n); ys = similar(xs); zs = similar(xs)
    @inbounds for i in 1:n; xs[i]=V[i][1]; ys[i]=V[i][2]; zs[i]=V[i][3]; end
    T = delaunay3d(xs, ys, zs; rng_seed=seed, perturb=false)   # EXACT coords ⇒ conforming
    m0 = to_mesh3(T)
    # per-surface bounding boxes: a centroid OUTSIDE a surface's bbox is definitively
    # outside that surface, so skip its ray-cast — the small coax pin's box excludes
    # almost every air/case tet. The grid index then makes each remaining ray-cast
    # O(1) in the surface's face count (see `_raygrid`), turning O(tets·Σtris) into
    # ~O(tets) with a small constant.
    grids = [_raygrid(s) for s in surfaces]
    bboxes = [begin
        lo=(Inf,Inf,Inf); hi=(-Inf,-Inf,-Inf)
        @inbounds for i in 1:size(s.coords,2)
            lo=(min(lo[1],s.coords[1,i]),min(lo[2],s.coords[2,i]),min(lo[3],s.coords[3,i]))
            hi=(max(hi[1],s.coords[1,i]),max(hi[2],s.coords[2,i]),max(hi[3],s.coords[3,i]))
        end; (lo,hi)
    end for s in surfaces]
    keep = Int32[]; tags = Int32[]
    @inbounds for t in 1:size(m0.tets, 2)
        a=m0.tets[1,t]; b=m0.tets[2,t]; c=m0.tets[3,t]; d=m0.tets[4,t]
        cx=(m0.coords[1,a]+m0.coords[1,b]+m0.coords[1,c]+m0.coords[1,d])/4
        cy=(m0.coords[2,a]+m0.coords[2,b]+m0.coords[2,c]+m0.coords[2,d])/4
        cz=(m0.coords[3,a]+m0.coords[3,b]+m0.coords[3,c]+m0.coords[3,d])/4
        for r in 1:length(surfaces)
            lo, hi = bboxes[r]
            (lo[1] <= cx <= hi[1] && lo[2] <= cy <= hi[2] && lo[3] <= cz <= hi[3]) || continue
            if _inside_grid((cx,cy,cz), grids[r])
                push!(keep, Int32(t)); push!(tags, Int32(r)); break
            end
        end
    end
    tets = Matrix{Int32}(undef, 4, length(keep))
    @inbounds for (j, t) in enumerate(keep), i in 1:4; tets[i,j] = m0.tets[i,t]; end
    mm = Mesh(m0.coords; tets=tets, tet_tag=tags)
    tpr = tets_per_region(mm)
    for r in 1:length(surfaces)
        get(tpr, Int32(r), 0) > 0 || error("tetrahedralize_conforming: region $r received no tet")
    end
    # Always-valid-or-explicit-blocker (PLAN principle #4): the single perturb=false
    # Delaunay can leave exact zero-volume/non-manifold tets on cospherical inputs (e.g.
    # an axis-aligned box assembly with shared corners). Validate before returning —
    # never hand back a silently invalid/non-conforming mesh.
    diag = validate(mm)
    diag.ok || throw(ErrorException("tetrahedralize_conforming: the joint Delaunay of the region " *
        "vertices is not a valid conforming mesh (" * join(diag.messages, "; ") * ") — a cospherical/" *
        "coplanar degeneracy the exact kernel cannot break into positive-volume tets. Use " *
        "mesh_box_regions for axis-aligned box unions, or recover_boundary per region."))
    try
        _certify_regions3(surfaces,mm,"tetrahedralize_conforming")
        return mm
    catch err
        err isa InterruptException && rethrow()
        (err isa ArgumentError || err isa ErrorException) || rethrow()
        if all(_bool_is_axis_aligned,surfaces)
            fallback=_axis_partition3(surfaces,"tetrahedralize_conforming")
            _certify_regions3(surfaces,fallback,"tetrahedralize_conforming")
            return fallback
        end
        if applicable(_recover_partition_exact,surfaces)
            return _recover_partition_exact(surfaces)
        end
        rethrow()
    end
end

"""
    tetrahedralize_conforming_exact(surfaces; rng_seed=1) -> Mesh

Like [`tetrahedralize_conforming`](@ref) but tetrahedralizes the shared region-vertex
set with the **exact-coordinate** Delaunay kernel ([`ExactMesh3D.delaunay3d_exact`](@ref))
instead of the Float64 `perturb=false` kernel. The exact kernel breaks the cospherical/
coplanar degeneracies *deterministically and validly* (no jitter, no degenerate tets),
so it **conformingly meshes multi-region assemblies the Float64 path cannot** — e.g. a
2×2×2 assembly of unit boxes (27 cospherical shared corners), where
`tetrahedralize_conforming` correctly raises its blocker. Each tet is region-tagged by
centroid ray-cast (same shared-Delaunay-face ⇒ conforming interfaces); the result is
validated and every region is checked non-empty, else an explicit blocker — never a
silent bad mesh. Costs exact `Rational{BigInt}` arithmetic (heavier than Float64), so
it is the robust fallback for degenerate inputs, not the default path.
"""
function tetrahedralize_conforming_exact(surfaces::AbstractVector{Mesh}; rng_seed::Integer=1)
    isempty(surfaces) && throw(ArgumentError("tetrahedralize_conforming_exact: no surfaces"))
    _seed3(rng_seed,"tetrahedralize_conforming_exact") # reserved for API parity/determinism
    length(surfaces)<=typemax(Int32) ||
        throw(ArgumentError("tetrahedralize_conforming_exact: region count exceeds Int32 tags"))
    for (r,s) in enumerate(surfaces);_require_surface3(s,"tetrahedralize_conforming_exact region $r");end
    seen = Dict{NTuple{3,Float64},Int32}()
    pts = NTuple{3,Rational{BigInt}}[]; xf = Float64[]; yf = Float64[]; zf = Float64[]
    @inbounds for s in surfaces, i in 1:size(s.coords, 2)
        k = (s.coords[1,i]==0 ? 0.0 : s.coords[1,i],
             s.coords[2,i]==0 ? 0.0 : s.coords[2,i],
             s.coords[3,i]==0 ? 0.0 : s.coords[3,i])
        if !haskey(seen, k)
            seen[k] = Int32(length(pts) + 1)
            push!(pts, (Rational{BigInt}(k[1]), Rational{BigInt}(k[2]), Rational{BigInt}(k[3])))
            push!(xf, k[1]); push!(yf, k[2]); push!(zf, k[3])
        end
    end
    length(pts) >= 4 || throw(ArgumentError("tetrahedralize_conforming_exact: need ≥ 4 distinct vertices"))
    length(pts)<=typemax(Int32) ||
        throw(ArgumentError("tetrahedralize_conforming_exact: distinct node count exceeds Int32"))
    etets = delaunay3d_exact(pts)
    np = length(pts)
    coords = Matrix{Float64}(undef, 3, np)
    @inbounds for i in 1:np; coords[1,i]=xf[i]; coords[2,i]=yf[i]; coords[3,i]=zf[i]; end
    grids = [_raygrid(s) for s in surfaces]
    bboxes = [begin
        lo=(Inf,Inf,Inf); hi=(-Inf,-Inf,-Inf)
        @inbounds for i in 1:size(s.coords,2)
            lo=(min(lo[1],s.coords[1,i]),min(lo[2],s.coords[2,i]),min(lo[3],s.coords[3,i]))
            hi=(max(hi[1],s.coords[1,i]),max(hi[2],s.coords[2,i]),max(hi[3],s.coords[3,i]))
        end; (lo,hi)
    end for s in surfaces]
    keep = NTuple{4,Int32}[]; tags = Int32[]
    @inbounds for t in etets
        cx=(xf[t[1]]+xf[t[2]]+xf[t[3]]+xf[t[4]])/4
        cy=(yf[t[1]]+yf[t[2]]+yf[t[3]]+yf[t[4]])/4
        cz=(zf[t[1]]+zf[t[2]]+zf[t[3]]+zf[t[4]])/4
        for r in 1:length(surfaces)
            lo,hi = bboxes[r]
            (lo[1] <= cx <= hi[1] && lo[2] <= cy <= hi[2] && lo[3] <= cz <= hi[3]) || continue
            if _inside_grid((cx,cy,cz), grids[r])
                push!(keep, (Int32(t[1]),Int32(t[2]),Int32(t[3]),Int32(t[4]))); push!(tags, Int32(r)); break
            end
        end
    end
    tets = Matrix{Int32}(undef, 4, length(keep))
    @inbounds for (j,t) in enumerate(keep); tets[1,j]=t[1]; tets[2,j]=t[2]; tets[3,j]=t[3]; tets[4,j]=t[4]; end
    mm = Mesh(coords; tets=tets, tet_tag=tags)
    tpr = tets_per_region(mm)
    for r in 1:length(surfaces)
        get(tpr, Int32(r), 0) > 0 || error("tetrahedralize_conforming_exact: region $r received no tet")
    end
    diag = validate(mm)
    diag.ok || throw(ErrorException("tetrahedralize_conforming_exact: produced an invalid mesh (" *
        join(diag.messages, "; ") * ")"))
    try
        _certify_regions3(surfaces,mm,"tetrahedralize_conforming_exact")
        return mm
    catch err
        err isa InterruptException && rethrow()
        (err isa ArgumentError || err isa ErrorException) || rethrow()
        if all(_bool_is_axis_aligned,surfaces)
            fallback=_axis_partition3(surfaces,"tetrahedralize_conforming_exact")
            _certify_regions3(surfaces,fallback,"tetrahedralize_conforming_exact")
            return fallback
        end
        applicable(_recover_partition_exact,surfaces) && return _recover_partition_exact(surfaces)
        rethrow()
    end
end

function _certify_regions3(surfaces::AbstractVector{Mesh},m::Mesh,caller::AbstractString)
    combined=_combined_surface3(surfaces)
    Px,Py,Pz,facets=_rb_dedup_surface(combined)
    regions=_rb_build_regions(Px,Py,Pz,facets)
    @inbounds for r in eachindex(surfaces)
        count_r=count(==(Int32(r)),m.tet_tag)
        count_r>0 || throw(ErrorException("$caller: effective region $r received no tetrahedra"))
        tm=Matrix{Int32}(undef,4,count_r);j=0
        for t in axes(m.tets,2)
            m.tet_tag[t]==Int32(r) || continue
            j+=1
            for k in 1:4;tm[k,j]=m.tets[k,t];end
        end
        part=Mesh(m.coords;tets=tm)
        d=validate(part)
        d.ok || throw(ErrorException("$caller: effective region $r is not a tetrahedral manifold — " *
                                     join(d.messages,"; ")))
        coord2pid=Dict{NTuple{3,Float64},Int32}()
        for i in eachindex(Px);coord2pid[(Px[i],Py[i],Pz[i])]=Int32(i);end
        Px2=copy(Px);Py2=copy(Py);Pz2=copy(Pz)
        for i in axes(part.coords,2)
            k=(part.coords[1,i],part.coords[2,i],part.coords[3,i])
            if !haskey(coord2pid,k)
                length(Px2)<typemax(Int32) || throw(ErrorException("$caller: certificate node count exceeds Int32"))
                push!(Px2,k[1]);push!(Py2,k[2]);push!(Pz2,k[3]);coord2pid[k]=Int32(length(Px2))
            end
        end
        _,mid=_rb_present_edges(part,coord2pid)
        ok,bad=_rb_boundary_in_surface(mid,part,regions,Px2,Py2,Pz2)
        ok || throw(ErrorException("$caller: effective region $r has boundary face $bad outside every input PLC"))
    end
    return nothing
end

function _combined_surface3(surfaces::AbstractVector{Mesh})
    ids=Dict{NTuple{3,Float64},Int32}();pts=NTuple{3,Float64}[]
    faces=Dict{NTuple{3,Int32},NTuple{3,Int32}}()
    for s in surfaces
        map=Vector{Int32}(undef,size(s.coords,2))
        @inbounds for i in eachindex(map)
            p=(s.coords[1,i]==0 ? 0.0 : s.coords[1,i],
               s.coords[2,i]==0 ? 0.0 : s.coords[2,i],
               s.coords[3,i]==0 ? 0.0 : s.coords[3,i])
            map[i]=get!(ids,p) do
                length(pts)<typemax(Int32) || throw(ArgumentError("combined PLC node count exceeds Int32"))
                push!(pts,p);Int32(length(pts))
            end
        end
        @inbounds for f in axes(s.tris,2)
            q=(map[s.tris[1,f]],map[s.tris[2,f]],map[s.tris[3,f]])
            faces[_sort3t(q...)]=q
        end
    end
    C=Matrix{Float64}(undef,3,length(pts));F=Matrix{Int32}(undef,3,length(faces))
    @inbounds for (i,p) in enumerate(pts);C[1,i]=p[1];C[2,i]=p[2];C[3,i]=p[3];end
    @inbounds for (j,q) in enumerate(values(faces));F[1,j]=q[1];F[2,j]=q[2];F[3,j]=q[3];end
    return Mesh(C;tris=F)
end

# Deterministic structured fallback for purely axis-aligned PLC partitions.  It
# uses every input face coordinate as a shared grid plane, classifies each open
# cell by the same first-containing-region rule, and applies one global Kuhn split.
# Thus shell/cavity partitions cannot acquire the tag-pinches possible when a
# degenerate joint Delaunay is classified only by tet centroids.
function _axis_partition3(surfaces::AbstractVector{Mesh},caller::AbstractString)
    canon(v)=v==0 ? 0.0 : v
    X=sort!(unique(canon.(vcat([collect(s.coords[1,:]) for s in surfaces]...))))
    Y=sort!(unique(canon.(vcat([collect(s.coords[2,:]) for s in surfaces]...))))
    Z=sort!(unique(canon.(vcat([collect(s.coords[3,:]) for s in surfaces]...))))
    nx=length(X)-1;ny=length(Y)-1;nz=length(Z)-1
    _checked_mul3(caller,"axis partition cell",nx,ny,nz)
    LX=length(X);LY=length(Y);LZ=length(Z)
    ng=_checked_mul3(caller,"axis partition node",LX,LY,LZ)
    ng<=typemax(Int32) || throw(ArgumentError("$caller: axis-partition grid exceeds Int32 nodes"))
    ntmax=_checked_mul3(caller,"axis partition tetrahedron",6,nx,ny,nz)
    ntmax<=typemax(Int32) || throw(ArgumentError("$caller: axis-partition tetrahedra exceed Int32"))
    gid(i,j,k)=((k*LY+j)*LX+i)+1
    grids=[_raygrid(s) for s in surfaces]
    tv=NTuple{4,Int32}[];tags=Int32[];used=Set{Int32}()
    paths=(((1,0,0),(1,1,0)),((1,0,0),(1,0,1)),((0,1,0),(1,1,0)),
           ((0,1,0),(0,1,1)),((0,0,1),(1,0,1)),((0,0,1),(0,1,1)))
    @inbounds for k in 0:nz-1,j in 0:ny-1,i in 0:nx-1
        p=(_midpoint3(X[i+1],X[i+2]),_midpoint3(Y[j+1],Y[j+2]),_midpoint3(Z[k+1],Z[k+2]))
        (X[i+1]<p[1]<X[i+2]&&Y[j+1]<p[2]<Y[j+2]&&Z[k+1]<p[3]<Z[k+2]) ||
            throw(ArgumentError("$caller: adjacent axis planes are below Float64 midpoint resolution"))
        tag=0
        for r in eachindex(grids);_inside_grid(p,grids[r])&&(tag=r;break);end
        tag==0&&continue
        cn(a,b,c)=Int32(gid(i+a,j+b,k+c));pc(a,b,c)=(X[i+a+1],Y[j+b+1],Z[k+c+1])
        v0=cn(0,0,0);v7=cn(1,1,1);p0=pc(0,0,0);p7=pc(1,1,1)
        for (s1,s2) in paths
            v1=cn(s1...);v2=cn(s2...);p1=pc(s1...);p2=pc(s2...)
            q=_signed_vol6(p0,p1,p2,p7)>=0 ? (v0,v1,v2,v7) : (v0,v1,v7,v2)
            push!(tv,q);push!(tags,Int32(tag));union!(used,q)
        end
    end
    isempty(tv)&&throw(ErrorException("$caller: axis partition is empty"))
    uv=sort!(collect(used));nid=Dict{Int32,Int32}(v=>Int32(i) for (i,v) in enumerate(uv))
    C=Matrix{Float64}(undef,3,length(uv))
    @inbounds for (n,v) in enumerate(uv)
        q=Int(v)-1;i=q%LX;j=(q÷LX)%LY;k=q÷(LX*LY)
        C[1,n]=X[i+1];C[2,n]=Y[j+1];C[3,n]=Z[k+1]
    end
    T=Matrix{Int32}(undef,4,length(tv))
    @inbounds for (t,q) in enumerate(tv),i in 1:4;T[i,t]=nid[q[i]];end
    out=Mesh(C;tets=T,tet_tag=tags);d=validate(out)
    d.ok||throw(ErrorException("$caller: structured fallback invalid — "*join(d.messages,"; ")))
    return out
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

# ── projected-plane grid index for point-in-surface parity classification ─────────
# Every parity ray-cast uses ONE fixed generic direction, so a triangle can be hit
# by a ray from p only if p projects (⊥ dir) into that triangle's 2-D bbox. We index
# the surface triangles in a uniform grid over that projection and, per query, test
# only the query cell's candidates with the SAME Möller–Trumbore `_ray_hits_tri`.
# This turns the O(n_tets·n_faces) classifier into O(n_faces build + O(1)/query).
#
# Correctness (why `_inside_grid` == `_inside_surface` bit-for-bit): a ray hit point
# q = p + s·dir satisfies proj(q) = proj(p) (dir is annihilated) and proj(q) ∈
# proj(Tri) (an affine map preserves barycentric membership), so proj(p) ∈ bbox of
# proj(Tri). Any triangle whose 2-D bbox excludes proj(p) is therefore never hit and
# dropping it cannot change the crossing parity. Each triangle is bucketed into every
# cell its bbox overlaps plus a 1-cell halo (which absorbs any float rounding at cell
# boundaries), so the query cell's candidate list is a conservative SUPERSET of the
# truly-hittable triangles ⇒ identical crossing count ⇒ identical parity.
const _CLASSIFY_DIR = _unitn((1.0, 0.3141592653589793, 0.01720209895))

struct _RayGrid
    dir::NTuple{3,Float64}
    u::NTuple{3,Float64}
    v::NTuple{3,Float64}
    origin::NTuple{3,Float64}
    smin::Float64; tmin::Float64
    sspan::Float64; tspan::Float64
    nx::Int; ny::Int
    off::Vector{Int32}      # CSR offsets, length nx*ny+1
    items::Vector{Int32}    # triangle ids grouped by cell
    tris::Matrix{Int32}
    coords::Matrix{Float64}
end

@inline _gridproj(p,o,a)=_dot(_subn(p,o),a)
@inline function _gridcell(q, qmin, qspan, n)
    qspan == 0 && return 0
    r=(q-qmin)/qspan
    isfinite(r) || throw(ArgumentError("Mesh3D classifier: non-finite projected grid coordinate"))
    return clamp(floor(Int,r*n),0,n-1)
end

# 2-D cell span [i0,i1]×[j0,j1] of triangle f's projected bbox, +1-cell halo.
@inline function _tri_cellspan(origin,smin,tmin,sspan,tspan,nx,ny,u,v,tris,coords,f)
    a = tris[1, f]; b = tris[2, f]; c = tris[3, f]
    pa = (coords[1, a], coords[2, a], coords[3, a])
    pb = (coords[1, b], coords[2, b], coords[3, b])
    pc = (coords[1, c], coords[2, c], coords[3, c])
    sa=_gridproj(pa,origin,u);sb=_gridproj(pb,origin,u);sc=_gridproj(pc,origin,u)
    ta=_gridproj(pa,origin,v);tb=_gridproj(pb,origin,v);tc=_gridproj(pc,origin,v)
    smn = min(sa, sb, sc); smx = max(sa, sb, sc)
    tmn = min(ta, tb, tc); tmx = max(ta, tb, tc)
    i0=clamp(_gridcell(smn,smin,sspan,nx)-1,0,nx-1)
    i1=clamp(_gridcell(smx,smin,sspan,nx)+1,0,nx-1)
    j0=clamp(_gridcell(tmn,tmin,tspan,ny)-1,0,ny-1)
    j1=clamp(_gridcell(tmx,tmin,tspan,ny)+1,0,ny-1)
    return (i0, i1, j0, j1)
end

function _raygrid(surface::Mesh, dir::NTuple{3,Float64}=_CLASSIFY_DIR)
    tris = surface.tris; coords = surface.coords
    nf = size(tris, 2); nc = size(coords, 2)
    nf <= typemax(Int32) ||
        throw(ArgumentError("Mesh3D classifier: $nf triangles exceed the Int32 CSR limit"))
    ndir=_unitn(dir)
    a = abs(ndir[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0)
    u = _unitn(_cross(ndir, a)); v = _cross(ndir, u)   # v is unit: dir⊥u, both unit
    origin=nc==0 ? (0.0,0.0,0.0) :
        (coords[1,1],coords[2,1],coords[3,1])
    if nf == 0 || nc == 0
        return _RayGrid(ndir,u,v,origin,0.0,0.0,0.0,0.0,1,1,Int32[1,1],Int32[],tris,coords)
    end
    smin = Inf; smax = -Inf; tmin = Inf; tmax = -Inf
    @inbounds for i in 1:nc
        p = (coords[1, i], coords[2, i], coords[3, i])
        (isfinite(p[1])&&isfinite(p[2])&&isfinite(p[3])) ||
            throw(ArgumentError("Mesh3D classifier: node $i has a non-finite coordinate"))
        s=_gridproj(p,origin,u);t=_gridproj(p,origin,v)
        (isfinite(s)&&isfinite(t)) ||
            throw(ArgumentError("Mesh3D classifier: projected coordinate overflowed Float64"))
        s < smin && (smin = s); s > smax && (smax = s)
        t < tmin && (tmin = t); t > tmax && (tmax = t)
    end
    sspan=smax-smin;tspan=tmax-tmin
    (isfinite(sspan)&&isfinite(tspan)&&sspan>=0&&tspan>=0) ||
        throw(ArgumentError("Mesh3D classifier: projected extent is non-finite"))
    csrmax=Int(typemax(Int32))-1
    ncell=max(1,floor(Int,sqrt(min(csrmax,2nf))))       # initially ~2·nf cells
    maxrefs=min(csrmax,nf>csrmax÷64 ? csrmax : max(nf,64nf))
    nx=ncell;ny=ncell;totalrefs=0
    while true
        nx=ncell;ny=ncell;totalrefs=0;fits=true
        @inbounds for f in 1:nf
            i0,i1,j0,j1=_tri_cellspan(origin,smin,tmin,sspan,tspan,nx,ny,u,v,tris,coords,f)
            add=(i1-i0+1)*(j1-j0+1)
            if add>maxrefs-totalrefs;fits=false;break;end
            totalrefs+=add
        end
        fits && break
        ncell==1 && throw(ArgumentError("Mesh3D classifier: CSR reference count exceeds Int32"))
        ncell=max(1,ncell÷2)
    end
    ncells=_checked_mul3("Mesh3D classifier","grid cell",nx,ny)
    cnt=zeros(Int32,ncells)
    @inbounds for f in 1:nf
        i0,i1,j0,j1=_tri_cellspan(origin,smin,tmin,sspan,tspan,nx,ny,u,v,tris,coords,f)
        for j in j0:j1,i in i0:i1;cnt[j*nx+i+1]+=Int32(1);end
    end
    off = Vector{Int32}(undef, ncells + 1)
    off[1] = Int32(1)
    @inbounds for c in 1:ncells
        off[c+1] = off[c] + cnt[c]
    end
    Int(off[end])-1==totalrefs || error("Mesh3D classifier: internal CSR count mismatch")
    items = Vector{Int32}(undef,totalrefs)
    cur = copy(off)
    @inbounds for f in 1:nf
        i0,i1,j0,j1=_tri_cellspan(origin,smin,tmin,sspan,tspan,nx,ny,u,v,tris,coords,f)
        for j in j0:j1, i in i0:i1
            c = j * nx + i + 1
            items[cur[c]] = Int32(f); cur[c] += Int32(1)
        end
    end
    return _RayGrid(ndir,u,v,origin,smin,tmin,sspan,tspan,nx,ny,off,items,tris,coords)
end

# parity ray-cast via the grid: identical crossing count to `_inside_surface`.
@inline function _inside_grid(p, g::_RayGrid)
    isempty(g.items) && return false
    s=_gridproj(p,g.origin,g.u);t=_gridproj(p,g.origin,g.v)
    ci=_gridcell(s,g.smin,g.sspan,g.nx)
    cj=_gridcell(t,g.tmin,g.tspan,g.ny)
    c = cj * g.nx + ci + 1
    crossings = 0
    tris = g.tris; coords = g.coords
    @inbounds for idx in g.off[c]:(g.off[c+1]-Int32(1))
        f = g.items[idx]
        a = tris[1, f]; b = tris[2, f]; cc = tris[3, f]
        pa = (coords[1, a], coords[2, a], coords[3, a])
        pb = (coords[1, b], coords[2, b], coords[3, b])
        pc = (coords[1, cc], coords[2, cc], coords[3, cc])
        _ray_hits_tri(p, g.dir, pa, pb, pc) && (crossings += 1)
    end
    return isodd(crossings)
end

function _classify_by_centroid(T::Triangulation3, surface::Mesh)
    keep = falses(length(T.alive))
    g = _raygrid(surface)
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t)) || continue
        a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3);d=_vert(T,t,4)
        cx=(T.x[a]+T.x[b]+T.x[c]+T.x[d])/4
        cy=(T.y[a]+T.y[b]+T.y[c]+T.y[d])/4
        cz=(T.z[a]+T.z[b]+T.z[c]+T.z[d])/4
        keep[t] = _inside_grid((cx,cy,cz), g)
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
    det == 0.0 && return false                  # ray parallel to triangle plane
    inv = 1.0/det
    tv = _subn(o, v0)
    u = _dot(tv, pv)*inv
    (u < 0.0 || u > 1.0) && return false
    qv = _cross(tv, e1)
    v = _dot(d, qv)*inv
    (v < 0.0 || u+v > 1.0) && return false
    t = _dot(e2, qv)*inv
    # Centroids/cell centres classified by this module are strictly off the input
    # surface, so the scale-independent geometric condition is simply t>0.  An
    # absolute floor incorrectly classified every solid smaller than that floor as
    # exterior (e.g. a valid 1e-15 box returned no tetrahedra).
    return t > 0.0
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


# ════════════════════════════════════════════════════════════════════════════════
# Boundary recovery — add to `module Mesh3D`.
#
# REQUIRED IMPORT ADDITION (module top): extend the existing MeshTypes import to
#     using ..MeshTypes: Mesh, tet_dihedral_extrema, validate, is_closed_manifold
# Everything else below uses ONLY existing Mesh3D internals (delaunay3d,
# _classify_by_centroid, to_mesh3, _is_ghost_tet, _vert, _pt, _sort3t) and the
# already-imported exact predicate `orient3`. No real src was modified during
# development; this block was verified by Core.eval-ing it (unqualified) into the
# live module and running the cases in the tests field.
#
# STRATEGY (grounded on measurement — see notes): on EXACT coordinates
# (perturb=false) the base incremental Delaunay of the surface vertices already
# yields a *geometrically conforming* boundary for closed axis-aligned/faceted PLCs
# (0 crease-segments missing, 0 tet edges piercing any coplanar region) — the only
# obstruction is rare zero-volume tets from cospherical degeneracy. So recovery =
# retry insertion orders (rng_seed), optionally drop kept zero-volume tets, and
# accept ONLY when an exact geometric-conformity + validity gate passes; otherwise
# throw naming the first unrecovered input facet. perturb=true is unusable (jitter
# moves facet vertices off-plane -> boundary shatters). The gate is triangulation-
# INDEPENDENT (a coplanar region re-triangulated with the SAME vertices still
# conforms), so it does not false-fail the box's opposite-diagonal faces, and it is
# SOUND (a genuinely missing facet — e.g. Schönhardt — is always rejected).
# ════════════════════════════════════════════════════════════════════════════════

# ---- exact rational geometry helpers ----
@inline _rbrat(x::Float64) = Rational{BigInt}(x)
@inline _rbsub(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _rbcross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _rbdot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _rbekey(a,b) = a<=b ? (Int32(a),Int32(b)) : (Int32(b),Int32(a))
@inline function _rbproj2(p, drop::Int)
    drop == 1 && return (p[2], p[3]); drop == 2 && return (p[1], p[3]); return (p[1], p[2])
end
@inline function _rborient2(a,b,c)::Int
    d = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]); d>0 ? 1 : (d<0 ? -1 : 0)
end

struct RBPlane
    A0::NTuple{3,Rational{BigInt}}; N::NTuple{3,Rational{BigInt}}; drop::Int
end
struct RBRegion
    facets::Vector{Int}
    tris::Vector{NTuple{3,Int32}}
    plane::RBPlane
    bnd_edges::Set{NTuple{2,Int32}}
    int_edges::Set{NTuple{2,Int32}}
end

mutable struct _RBUF; p::Vector{Int}; end
_rbuf(n) = _RBUF(collect(1:n))
function _rbfind(u::_RBUF,i); while u.p[i]!=i; u.p[i]=u.p[u.p[i]]; i=u.p[i]; end; i; end
_rbunion(u::_RBUF,i,j) = (u.p[_rbfind(u,i)] = _rbfind(u,j))

@inline _rbptP(Px,Py,Pz,i) = (Px[i],Py[i],Pz[i])
@inline _rbratP(Px,Py,Pz,i) = (_rbrat(Px[i]),_rbrat(Py[i]),_rbrat(Pz[i]))

# dedup surface coords -> (Px,Py,Pz, facets in P-id triples)
function _rb_dedup_surface(surface::Mesh)
    seen = Dict{NTuple{3,Float64},Int32}()
    Px=Float64[]; Py=Float64[]; Pz=Float64[]
    remap = Vector{Int32}(undef, size(surface.coords,2))
    @inbounds for i in 1:size(surface.coords,2)
        key=(surface.coords[1,i],surface.coords[2,i],surface.coords[3,i])
        id=get(seen,key,Int32(0))
        if id==0
            push!(Px,key[1]);push!(Py,key[2]);push!(Pz,key[3]); id=Int32(length(Px)); seen[key]=id
        end
        remap[i]=id
    end
    facets = NTuple{3,Int32}[]
    @inbounds for f in 1:size(surface.tris,2)
        a=remap[surface.tris[1,f]]; b=remap[surface.tris[2,f]]; c=remap[surface.tris[3,f]]
        (a==b||b==c||a==c) && continue
        push!(facets,(a,b,c))
    end
    return Px,Py,Pz,facets
end

# coplanar-region grouping. crease edge = incident facets non-coplanar (orient3!=0).
function _rb_build_regions(Px,Py,Pz, facets::Vector{NTuple{3,Int32}})
    nf = length(facets)
    e2f = Dict{NTuple{2,Int32},Vector{Int}}()
    for (fi,(a,b,c)) in enumerate(facets)
        for e in (_rbekey(a,b),_rbekey(b,c),_rbekey(a,c)); push!(get!(e2f,e,Int[]), fi); end
    end
    uf = _rbuf(nf)
    for (e,fs) in e2f
        length(fs)==2 || continue
        f1=facets[fs[1]]; f2=facets[fs[2]]; u,v=e
        w1=f1[1]; (w1==u||w1==v)&&(w1=f1[2]); (w1==u||w1==v)&&(w1=f1[3])
        w2=f2[1]; (w2==u||w2==v)&&(w2=f2[2]); (w2==u||w2==v)&&(w2=f2[3])
        orient3(_rbptP(Px,Py,Pz,u),_rbptP(Px,Py,Pz,v),_rbptP(Px,Py,Pz,w1),_rbptP(Px,Py,Pz,w2))==0 &&
            _rbunion(uf,fs[1],fs[2])
    end
    comp = Dict{Int,Vector{Int}}()
    for fi in 1:nf; push!(get!(comp,_rbfind(uf,fi),Int[]), fi); end
    regions = RBRegion[]
    for (_,fis) in comp
        tris=[facets[fi] for fi in fis]
        a,b,c = tris[1]
        A0=_rbratP(Px,Py,Pz,a); Br=_rbratP(Px,Py,Pz,b); Cr=_rbratP(Px,Py,Pz,c)
        N=_rbcross(_rbsub(Br,A0),_rbsub(Cr,A0))
        an=abs.(N); drop = an[1]>=an[2] ? (an[1]>=an[3] ? 1 : 3) : (an[2]>=an[3] ? 2 : 3)
        inc=Dict{NTuple{2,Int32},Int}()
        for (x,y,z) in tris, e in (_rbekey(x,y),_rbekey(y,z),_rbekey(x,z)); inc[e]=get(inc,e,0)+1; end
        bnd=Set{NTuple{2,Int32}}(); intl=Set{NTuple{2,Int32}}()
        for (e,ci) in inc; (ci==1 ? push!(bnd,e) : push!(intl,e)); end
        push!(regions, RBRegion(fis,tris,RBPlane(A0,N,drop),bnd,intl))
    end
    return regions
end

function _rb_crease_segments(regions::Vector{RBRegion})
    S=Set{NTuple{2,Int32}}(); for r in regions, e in r.bnd_edges; push!(S,e); end; S
end

@inline function _rbside(pl::RBPlane, x)::Int
    d=_rbdot(pl.N,_rbsub(x,pl.A0)); d>0 ? 1 : (d<0 ? -1 : 0)
end

# is rational in-plane point y strictly interior to region (sound pierce test)?
function _rb_pierces_region(reg::RBRegion, Px,Py,Pz, y)
    drop=reg.plane.drop; y2=_rbproj2(y,drop)
    for (a,b,c) in reg.tris
        pa=_rbproj2(_rbratP(Px,Py,Pz,a),drop); pb=_rbproj2(_rbratP(Px,Py,Pz,b),drop); pc=_rbproj2(_rbratP(Px,Py,Pz,c),drop)
        s1=_rborient2(pa,pb,y2); s2=_rborient2(pb,pc,y2); s3=_rborient2(pc,pa,y2)
        if s1!=0 && s2!=0 && s3!=0 && s1==s2 && s2==s3
            return true
        end
        zc=(s1==0)+(s2==0)+(s3==0)
        if zc==1
            nz=filter(!=(0),(s1,s2,s3))
            if length(nz)==2 && nz[1]==nz[2]
                e = s1==0 ? _rbekey(a,b) : (s2==0 ? _rbekey(b,c) : _rbekey(c,a))
                e in reg.int_edges && return true
            end
        end
    end
    return false
end

function _rb_edge_pierces(reg::RBRegion, Px,Py,Pz, p::Int32, q::Int32)
    pp=_rbratP(Px,Py,Pz,p); qq=_rbratP(Px,Py,Pz,q)
    sp=_rbside(reg.plane,pp); sq=_rbside(reg.plane,qq)
    (sp!=0 && sq!=0 && sp!=sq) || return false
    num=_rbdot(reg.plane.N,_rbsub(reg.plane.A0,pp)); den=_rbdot(reg.plane.N,_rbsub(qq,pp))
    den==0 && return false
    t=num//den
    y=(pp[1]+t*(qq[1]-pp[1]), pp[2]+t*(qq[2]-pp[2]), pp[3]+t*(qq[3]-pp[3]))
    return _rb_pierces_region(reg,Px,Py,Pz,y)
end

# on-segment vertices sorted along u->v (exact)
function _rb_onseg(Px,Py,Pz,u::Int32,v::Int32,nreal::Int)
    pu=_rbratP(Px,Py,Pz,u); pv=_rbratP(Px,Py,Pz,v); d=_rbsub(pv,pu); dd=_rbdot(d,d)
    out=Tuple{Rational{BigInt},Int32}[]
    for w in 1:nreal
        (w==u||w==v) && continue
        pw=_rbratP(Px,Py,Pz,w); cr=_rbcross(_rbsub(pw,pu),d)
        (cr[1]==0&&cr[2]==0&&cr[3]==0) || continue
        s=_rbdot(_rbsub(pw,pu),d); (s>0&&s<dd)||continue
        push!(out,(s//dd,Int32(w)))
    end
    sort!(out,by=x->x[1]); [w for (_,w) in out]
end
function _rb_subsegments(Px,Py,Pz,u::Int32,v::Int32,nreal::Int)
    mids=_rb_onseg(Px,Py,Pz,u,v,nreal); chain=vcat(Int32[u],mids,Int32[v])
    [(chain[i],chain[i+1]) for i in 1:length(chain)-1]
end

# present undirected edges of a Mesh's tets, mapped to P ids (coord match)
function _rb_present_edges(m::Mesh, coord2pid::Dict{NTuple{3,Float64},Int32})
    mid=Vector{Int32}(undef,size(m.coords,2))
    for i in 1:size(m.coords,2)
        mid[i]=coord2pid[(m.coords[1,i],m.coords[2,i],m.coords[3,i])]
    end
    E=Set{NTuple{2,Int32}}()
    for t in 1:size(m.tets,2)
        vs=(mid[m.tets[1,t]],mid[m.tets[2,t]],mid[m.tets[3,t]],mid[m.tets[4,t]])
        for i in 1:4,j in i+1:4; push!(E,_rbekey(vs[i],vs[j])); end
    end
    return E, mid
end

# drop kept flats -> new keep mask (zero-volume tets contribute 0 volume; dropping a
# boundary-plane flat just re-triangulates that facet with the opposite diagonal — the
# gate rejects any drop that removes real coverage).
function _rb_drop_flats(T, keep)
    kk=copy(keep)
    for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T,t) && t<=length(kk) && kk[t]) || continue
        a=_vert(T,t,1);b=_vert(T,t,2);c=_vert(T,t,3);d=_vert(T,t,4)
        orient3(_pt(T,a),_pt(T,b),_pt(T,c),_pt(T,d))==0 && (kk[t]=false)
    end
    kk
end

# closed point-in-region: rational in-plane point y is inside the closed region
# (inside or on the boundary of some region triangle). Certifies a tet-mesh boundary
# face actually lies on the input surface (guards against exterior bulges).
function _rb_point_in_region_closed(reg::RBRegion, Px,Py,Pz, y)
    drop=reg.plane.drop; y2=_rbproj2(y,drop)
    for (a,b,c) in reg.tris
        pa=_rbproj2(_rbratP(Px,Py,Pz,a),drop); pb=_rbproj2(_rbratP(Px,Py,Pz,b),drop); pc=_rbproj2(_rbratP(Px,Py,Pz,c),drop)
        s1=_rborient2(pa,pb,y2); s2=_rborient2(pb,pc,y2); s3=_rborient2(pc,pa,y2)
        haspos = (s1>0)||(s2>0)||(s3>0); hasneg=(s1<0)||(s2<0)||(s3<0)
        (haspos && hasneg) || return true
    end
    return false
end

# every boundary face of the kept tet mesh must lie in some input region (on-plane +
# inside): certifies tet-mesh boundary ⊆ input surface (no exterior bulge). P ids.
function _rb_boundary_in_surface(mid, m::Mesh, regions::Vector{RBRegion}, Px,Py,Pz)
    inc=Dict{NTuple{3,Int32},Int}()
    for t in 1:size(m.tets,2)
        v1=mid[m.tets[1,t]];v2=mid[m.tets[2,t]];v3=mid[m.tets[3,t]];v4=mid[m.tets[4,t]]
        for f in (_sort3t(v2,v3,v4),_sort3t(v1,v3,v4),_sort3t(v1,v2,v4),_sort3t(v1,v2,v3))
            inc[f]=get(inc,f,0)+1
        end
    end
    for (f,c) in inc
        c==1 || continue
        a,b,cc=f
        ratc=((_rbratP(Px,Py,Pz,a)[1]+_rbratP(Px,Py,Pz,b)[1]+_rbratP(Px,Py,Pz,cc)[1])//3,
              (_rbratP(Px,Py,Pz,a)[2]+_rbratP(Px,Py,Pz,b)[2]+_rbratP(Px,Py,Pz,cc)[2])//3,
              (_rbratP(Px,Py,Pz,a)[3]+_rbratP(Px,Py,Pz,b)[3]+_rbratP(Px,Py,Pz,cc)[3])//3)
        found=false
        for r in regions
            (_rbside(r.plane,_rbratP(Px,Py,Pz,a))==0 && _rbside(r.plane,_rbratP(Px,Py,Pz,b))==0 &&
             _rbside(r.plane,_rbratP(Px,Py,Pz,cc))==0) || continue
            _rb_point_in_region_closed(r,Px,Py,Pz,ratc) && (found=true; break)
        end
        found || return (false, f)
    end
    return (true, (Int32(0),Int32(0),Int32(0)))
end

# ---- the exact geometric conformity + validity gate ----
# returns (ok, nrecovered, reason). nrecovered = #input facets in conforming regions.
function _rb_gate(surface::Mesh, m::Mesh, regions::Vector{RBRegion}, S, Px,Py,Pz, facets)
    v = validate(m)
    v.ok || return (false, 0, "invalid tet mesh: $(join(v.messages, "; "))")
    # `validate` already runs the complete vertex-link manifold audit.  Repeating
    # `is_closed_manifold` here would rebuild all topology records and double the
    # peak work of every PLC certificate.
    coord2pid=Dict{NTuple{3,Float64},Int32}()
    for i in 1:length(Px); coord2pid[(Px[i],Py[i],Pz[i])]=Int32(i); end
    Px2=copy(Px);Py2=copy(Py);Pz2=copy(Pz)
    for i in 1:size(m.coords,2)
        k=(m.coords[1,i],m.coords[2,i],m.coords[3,i])
        haskey(coord2pid,k) || (push!(Px2,k[1]);push!(Py2,k[2]);push!(Pz2,k[3]);coord2pid[k]=Int32(length(Px2)))
    end
    E,mid = _rb_present_edges(m, coord2pid)
    nreal=length(Px2)
    facet_ok = trues(length(facets))
    reg_of = Dict{Int,Int}()
    for (ri,r) in enumerate(regions), fi in r.facets; reg_of[fi]=ri; end
    bad_region_reason = ""
    region_bad = falses(length(regions))
    for (ri,r) in enumerate(regions)
        for (u,v2) in r.bnd_edges
            present=true
            for (a,b) in _rb_subsegments(Px2,Py2,Pz2,u,v2,nreal)
                (_rbekey(a,b) in E) || (present=false; break)
            end
            if !present
                region_bad[ri]=true
                bad_region_reason=="" && (bad_region_reason="crease segment ($(u),$(v2)) not recovered")
                break
            end
        end
        if !region_bad[ri]
            for (p,q) in E
                if _rb_edge_pierces(r,Px2,Py2,Pz2,p,q)
                    region_bad[ri]=true
                    bad_region_reason=="" && (bad_region_reason="tet edge ($(p),$(q)) pierces region interior")
                    break
                end
            end
        end
    end
    nrec=0
    for fi in 1:length(facets)
        ri=reg_of[fi]
        if region_bad[ri]; facet_ok[fi]=false; else; nrec+=1; end
    end
    if nrec < length(facets)
        fi = findfirst(==(false), facet_ok)
        f = facets[fi]
        return (false, nrec, "facet $(f) (region $(reg_of[fi])) not conforming: $(bad_region_reason)")
    end
    binok, badface = _rb_boundary_in_surface(mid, m, regions, Px2,Py2,Pz2)
    binok || return (false, nrec, "tet-mesh boundary face $(badface) lies outside the input surface")
    return (true, nrec, "ok")
end

# Certify that an arbitrary tet fill has the same piecewise-linear boundary as
# `surface`.  This is the reusable front door to the exact recovery gate: validity
# and closed-manifold checks alone are insufficient because a restricted Delaunay
# fill can cap a non-convex PLC while remaining a perfectly valid convex-hull mesh.
function _certify_surface_fill(surface::Mesh, m::Mesh)
    Px, Py, Pz, facets = _rb_dedup_surface(surface)
    length(Px) >= 4 || return (false, "input has fewer than four distinct vertices")
    isempty(facets) && return (false, "input has no non-degenerate facets")
    regions = _rb_build_regions(Px, Py, Pz, facets)
    isempty(regions) && return (false, "input has no non-degenerate surface regions")
    gate = _rb_gate(surface, m, regions, _rb_crease_segments(regions),
                    Px, Py, Pz, facets)
    return (gate[1], gate[3])
end

# Consistently orient one connected closed triangle surface without trusting input
# winding.  Adjacent faces must traverse their shared edge in opposite directions.
# Returns a copy with the necessary flips, or `nothing` for an open, non-manifold,
# disconnected, or non-orientable facet complex.  The fan fallback intentionally
# handles one boundary component only; cavities are not star-shaped fan domains.
function _rb_orient_facets(facets::Vector{NTuple{3,Int32}})
    nf = length(facets)
    nf == 0 && return nothing
    incid = Dict{NTuple{2,Int32},Vector{Tuple{Int32,Bool}}}()
    sizehint!(incid, 3nf ÷ 2)
    @inbounds for (fi, (a,b,c)) in enumerate(facets)
        for (u,v) in ((a,b),(b,c),(c,a))
            key = _rbekey(u,v)
            push!(get!(() -> Tuple{Int32,Bool}[], incid, key),
                  (Int32(fi), u < v))
        end
    end
    adjacency = [Tuple{Int32,Bool}[] for _ in 1:nf] # neighbour, same original direction
    for entries in values(incid)
        length(entries) == 2 || return nothing
        (f1,d1), (f2,d2) = entries
        same = d1 == d2
        push!(adjacency[f1], (f2,same)); push!(adjacency[f2], (f1,same))
    end
    seen = falses(nf); flipped = falses(nf)
    stack = Int32[1]; seen[1] = true; nseen = 0
    while !isempty(stack)
        f = pop!(stack); nseen += 1
        @inbounds for (g,same) in adjacency[f]
            want = xor(flipped[f], same)
            if seen[g]
                flipped[g] == want || return nothing
            else
                seen[g] = true; flipped[g] = want; push!(stack, g)
            end
        end
    end
    nseen == nf || return nothing
    oriented = copy(facets)
    @inbounds for i in 1:nf
        if flipped[i]
            a,b,c = oriented[i]; oriented[i] = (a,c,b)
        end
    end
    return oriented
end

# Steiner fallback for a STAR-SHAPED polyhedron (e.g. Schönhardt — not tetrahedraliz-
# able without a Steiner point): orient the closed surface, find an interior kernel
# point that sees every facet (exact orient3 shares one sign across all consistently
# oriented faces), then fan-tetrahedralize — one tet {facet ∪ p} per facet.
# Conforming by construction and all-positive. Returns the fan Mesh, or `nothing` if
# no simple candidate kernel point works. Quality is fan-poor; a caller wanting
# quality should run `optimize_flips!` after.
function _rb_fan_steiner(Px,Py,Pz, facets::Vector{NTuple{3,Int32}})
    oriented = _rb_orient_facets(facets)
    oriented === nothing && return nothing
    nn = length(Px); nf = length(oriented)
    vt(i) = (Px[i], Py[i], Pz[i])
    vc = (sum(Px)/nn, sum(Py)/nn, sum(Pz)/nn)
    fcx=0.0; fcy=0.0; fcz=0.0
    @inbounds for (a,b,c) in oriented
        fcx += (Px[a]+Px[b]+Px[c])/3; fcy += (Py[a]+Py[b]+Py[c])/3; fcz += (Pz[a]+Pz[b]+Pz[c])/3
    end
    fca = (fcx/nf, fcy/nf, fcz/nf)
    function kernel(p)                                   # p strictly interior to every face plane?
        s0 = 0
        @inbounds for (a,b,c) in oriented
            o = orient3(vt(a),vt(b),vt(c),p); o == 0 && return false
            sg = o > 0 ? 1 : -1
            s0 == 0 ? (s0 = sg) : (sg == s0 || return false)
        end
        return true
    end
    p = kernel(vc) ? vc : (kernel(fca) ? fca : nothing)
    p === nothing && return nothing
    C = Matrix{Float64}(undef, 3, nn+1)
    @inbounds for i in 1:nn; C[1,i]=Px[i]; C[2,i]=Py[i]; C[3,i]=Pz[i]; end
    C[1,nn+1]=p[1]; C[2,nn+1]=p[2]; C[3,nn+1]=p[3]; pv=Int32(nn+1)
    tets = Matrix{Int32}(undef, 4, nf)
    @inbounds for (t,(a,b,c)) in enumerate(oriented)
        d0=(C[1,a],C[2,a],C[3,a]); d1=(C[1,b],C[2,b],C[3,b]); d2=(C[1,c],C[2,c],C[3,c]); dp=(p[1],p[2],p[3])
        if _signed_vol6(d0,d1,d2,dp) >= 0
            tets[1,t]=a; tets[2,t]=b; tets[3,t]=c; tets[4,t]=pv
        else
            tets[1,t]=a; tets[2,t]=c; tets[3,t]=b; tets[4,t]=pv
        end
    end
    return Mesh(C; tets=tets)
end

"""
    recover_boundary(surface::Mesh; rng_seed=1, max_seeds=64, steiner=false) -> Mesh

Recover the closed, watertight, manifold PLC `surface` (vertices + triangular
facets, possibly NON-CONVEX) as a CONFORMING interior tet mesh: every input facet
appears — by vertex identity, a coplanar region re-triangulated with the SAME
vertices still conforms — as a face of the tet mesh; the mesh is valid (all tets
positive, no zero-volume, manifold <=2 faces/edge); interior/exterior is correctly
classified (ray-cast, orientation-independent). Returns the interior tet [`Mesh`],
or THROWS an explicit blocker naming the first unrecovered input facet (never a
silently non-conforming mesh).

`steiner=true` enables a **Steiner-point fallback** for genuinely non-tetrahedral-
izable but STAR-SHAPED inputs (e.g. the Schönhardt polyhedron): when Delaunay
recovery cannot conform, fan-tetrahedralize from an interior kernel point (one
Steiner vertex, one tet per facet — conforming and valid, but fan-quality; run
`optimize_flips!` for quality). If the input is not star-shaped from a simple
kernel candidate, it still throws the explicit blocker.
"""
function recover_boundary(surface::Mesh; rng_seed::Integer=1, max_seeds::Integer=64,
                          steiner::Bool=false)
    _require_surface3(surface,"recover_boundary")
    base=_seed3(rng_seed,"recover_boundary")
    (1 <= max_seeds <= typemax(Int)) ||
        throw(ArgumentError("recover_boundary: max_seeds must be in 1:$(typemax(Int)) (got $max_seeds)"))
    nseeds=Int(max_seeds)
    base <= typemax(Int)-(nseeds-1) ||
        throw(ArgumentError("recover_boundary: rng_seed + max_seeds overflows Int"))
    Px,Py,Pz,facets = _rb_dedup_surface(surface)
    length(Px) >= 4 || throw(ArgumentError("recover_boundary: need >= 4 distinct surface vertices"))
    isempty(facets) && throw(ArgumentError("recover_boundary: surface has no facets"))
    regions = _rb_build_regions(Px,Py,Pz,facets)
    S = _rb_crease_segments(regions)
    lastreason = "no seed produced a conforming valid mesh"
    lastnrec = 0
    for k in 0:nseeds-1
        seed = base + k
        try
            T = delaunay3d(copy(Px),copy(Py),copy(Pz); perturb=false, rng_seed=seed)
            keep = _classify_by_centroid(T, surface)
            for drop in (false, true)
                kk = drop ? _rb_drop_flats(T, keep) : keep
                any(kk) || continue
                m = to_mesh3(T; keep=kk)
                ok, nrec, reason = _rb_gate(surface, m, regions, S, Px,Py,Pz, facets)
                ok && return m
                nrec >= lastnrec && (lastnrec = nrec; lastreason = reason)
            end
        catch err
            err isa InterruptException && rethrow()
            (err isa ArgumentError || err isa ErrorException) || rethrow()
            lastnrec == 0 &&
                (lastreason = "seed $seed failed internally: $(sprint(showerror, err))")
        end
    end
    # Steiner fallback for star-shaped non-tetrahedralizable inputs (e.g. Schönhardt):
    # fan-tetrahedralize from an interior kernel point. Conforming by construction
    # (each input facet is a literal tet face); accept only if valid + closed-manifold.
    if steiner
        fm = _rb_fan_steiner(Px,Py,Pz, facets)
        if fm !== nothing && validate(fm).ok && is_closed_manifold(fm)
            return fm
        end
    end
    throw(ErrorException("recover_boundary: recovery FAILED after $nseeds seeds " *
        "($(lastnrec)/$(length(facets)) facets recovered)" *
        (steiner ? " and the Steiner fan fallback did not apply (input not star-shaped from a simple kernel point)" : "") *
        ". First blocker: $lastreason" *
        (steiner ? "" : " (try steiner=true for star-shaped non-tetrahedralizable inputs)")))
end



# ═══════════════════════════════════════════════════════════════════════════════════
# mesh_boolean — native mesh-Boolean CSG for module Mesh3D (no OpenCASCADE).
#
# Add these two lines to the module's existing `using` block (both submodules are
# included before Mesh3D in Tessella.jl, so `..Mesh2D` / `..MeshTypes` resolve):
#     using ..Mesh2D: constrained_delaunay, to_mesh
#     using ..MeshTypes: boundary_edges
# (orient3 is already imported; _inside_surface, _unitn are Mesh3D internals.)
# Add `mesh_boolean` to the module's `export` list.
#
# Two exact engines, dispatched automatically:
#   • AXIS-ALIGNED inputs (boxes / box-unions, incl. every coplanar shared-face case)
#     -> exact plane-arrangement path: cut space by all face planes, classify each
#        grid cell by centroid ray-cast, emit cell faces separating kept/not-kept.
#        Volume is the exact sum of grid-cell volumes; surface watertight by construction.
#   • GENERAL-POSITION inputs (e.g. box × cylinder) -> Cork/libigl-style tri-tri path:
#        exact orient3 tri-tri intersection with Rational{BigInt} endpoints (bit-identical
#        shared seam) -> coplanar-region grouping -> per-region Mesh2D CDT retriangulation
#        with the intersection chords as constraints -> centroid ray-cast classification
#        vs the ORIGINAL opposite closed solid -> per-op assembly.  Any non-general-position
#        degeneracy in this path is an explicit THROW, never a wrong surface.
# The assembled result is verified watertight+manifold before return.
# ═══════════════════════════════════════════════════════════════════════════════════


# ── SECTION A — oracles: divergence-theorem volume + watertight/manifold gate ──────

# Signed volume of a closed, consistently outward-oriented triangle surface by the
# divergence theorem: V = (1/6) Σ_faces v0·(v1×v2). Positive for outward normals.
function _bool_signed_volume(m::Mesh)
    V = 0.0
    @inbounds for f in 1:size(m.tris,2)
        a=m.tris[1,f]; b=m.tris[2,f]; c=m.tris[3,f]
        v0=(m.coords[1,a],m.coords[2,a],m.coords[3,a])
        v1=(m.coords[1,b],m.coords[2,b],m.coords[3,b])
        v2=(m.coords[1,c],m.coords[2,c],m.coords[3,c])
        cx = v1[2]*v2[3]-v1[3]*v2[2]
        cy = v1[3]*v2[1]-v1[1]*v2[3]
        cz = v1[1]*v2[2]-v1[2]*v2[1]
        V += v0[1]*cx + v0[2]*cy + v0[3]*cz
    end
    return V/6.0
end

# (watertight, manifold): watertight ⇔ no boundary edge; manifold ⇔ max edge
# incidence ≤ 2. Together ⇒ every edge shared by exactly two triangles (closed 2-mfd).
function _bool_watertight(m::Mesh)
    bnd, maxinc = boundary_edges(m.tris)
    return (isempty(bnd), isempty(bnd) && maxinc <= 2, length(bnd), maxinc)
end

# ── SECTION A2 — input normalization ───────────────────────────────────────────────

# Dedup coords by exact value; drop degenerate tris; return a clean Mesh (tris only).
function _bool_clean(m::Mesh)
    seen = Dict{NTuple{3,Float64},Int32}()
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    remap = Vector{Int32}(undef, size(m.coords,2))
    @inbounds for i in 1:size(m.coords,2)
        k=(m.coords[1,i]==0 ? 0.0 : m.coords[1,i],
           m.coords[2,i]==0 ? 0.0 : m.coords[2,i],
           m.coords[3,i]==0 ? 0.0 : m.coords[3,i])
        id=get(seen,k,Int32(0))
        if id==0
            push!(xs,k[1]); push!(ys,k[2]); push!(zs,k[3]); id=Int32(length(xs)); seen[k]=id
        end
        remap[i]=id
    end
    tris=NTuple{3,Int32}[]
    @inbounds for f in 1:size(m.tris,2)
        a=remap[m.tris[1,f]]; b=remap[m.tris[2,f]]; c=remap[m.tris[3,f]]
        (a==b||b==c||a==c) && continue
        push!(tris,(a,b,c))
    end
    C=Matrix{Float64}(undef,3,length(xs))
    @inbounds for i in 1:length(xs); C[1,i]=xs[i]; C[2,i]=ys[i]; C[3,i]=zs[i]; end
    tm=Matrix{Int32}(undef,3,length(tris))
    @inbounds for (j,t) in enumerate(tris); tm[1,j]=t[1]; tm[2,j]=t[2]; tm[3,j]=t[3]; end
    return Mesh(C; tris=tm)
end

# Ensure outward orientation: if the divergence-theorem signed volume is negative,
# flip every triangle's winding. Returns an outward-oriented Mesh and its |volume|.
function _bool_orient_outward(m::Mesh)
    v = _bool_signed_volume(m)
    isfinite(v) || throw(ArgumentError("mesh_boolean: input surface signed volume is non-finite"))
    v == 0 && throw(ArgumentError("mesh_boolean: input surface has zero signed volume (not a closed solid)"))
    if v < 0
        tm = copy(m.tris)
        @inbounds for f in 1:size(tm,2); tm[2,f],tm[3,f] = tm[3,f],tm[2,f]; end
        return Mesh(copy(m.coords); tris=tm), -v
    end
    return m, v
end

function _bool_prepare(m::Mesh, name::AbstractString)
    _require_surface3(m,"mesh_boolean input $name";oriented=true)
    c = _bool_clean(m)
    wt, mf, nb, mx = _bool_watertight(c)
    wt || throw(ArgumentError("mesh_boolean: input $name is not watertight ($nb boundary edges)"))
    mf || throw(ArgumentError("mesh_boolean: input $name is not manifold (max edge incidence $mx)"))
    o, _ = _bool_orient_outward(c)
    return o
end

# ── SECTION B — exact plane-arrangement Boolean for AXIS-ALIGNED solids ─────────────
# (both inputs' faces all lie in axis planes ⇒ boxes / box-unions).  Handles ALL
#  coplanar/shared-face degeneracies exactly; volume is the exact sum of grid-cell
#  volumes; surface is watertight+manifold by construction.

# A triangle is axis-aligned iff its 3 vertices share one coordinate (it lies in an
# x=, y=, or z= plane). A solid is axis-aligned iff all its triangles are.
function _bool_is_axis_aligned(m::Mesh)
    @inbounds for f in 1:size(m.tris,2)
        a=m.tris[1,f]; b=m.tris[2,f]; c=m.tris[3,f]
        sx = m.coords[1,a]==m.coords[1,b]==m.coords[1,c]
        sy = m.coords[2,a]==m.coords[2,b]==m.coords[2,c]
        sz = m.coords[3,a]==m.coords[3,b]==m.coords[3,c]
        (sx||sy||sz) || return false
    end
    return true
end

# Emit the 2 outward triangles of an axis-aligned rectangle face. `ax` ∈ 1:3 is the
# constant axis, `cval` its coordinate; (a0,a1),(b0,b1) span axes axA,axB; `outsign`
# = +1/−1 outward direction along `ax`. Corners deduped into `V`; tris appended.
function _bool_emit_face!(V::Dict{NTuple{3,Float64},Int32}, coords::Vector{NTuple{3,Float64}},
                          Tr::Vector{NTuple{3,Int32}}, ax::Int, cval::Float64,
                          a0::Float64,a1::Float64, b0::Float64,b1::Float64,
                          axA::Int, axB::Int, outsign::Int)
    function vid(p::NTuple{3,Float64})
        id=get(V,p,Int32(0))
        if id==0
            length(coords)<typemax(Int32) ||
                throw(ArgumentError("mesh_boolean: result node count exceeds Int32"))
            push!(coords,p); id=Int32(length(coords)); V[p]=id
        end
        return id
    end
    mk(a,b) = begin
        p=zeros(3); p[ax]=cval; p[axA]=a; p[axB]=b; (p[1],p[2],p[3])
    end
    p00=mk(a0,b0); p10=mk(a1,b0); p11=mk(a1,b1); p01=mk(a0,b1)
    i00=vid(p00); i10=vid(p10); i11=vid(p11); i01=vid(p01)
    # normal of (p00->p10->p11); flip winding so its `ax` component matches outsign.
    e1=(p10[1]-p00[1],p10[2]-p00[2],p10[3]-p00[3])
    e2=(p11[1]-p00[1],p11[2]-p00[2],p11[3]-p00[3])
    n=(e1[2]*e2[3]-e1[3]*e2[2], e1[3]*e2[1]-e1[1]*e2[3], e1[1]*e2[2]-e1[2]*e2[1])
    ndir = n[ax] > 0 ? 1 : -1
    if ndir == outsign
        push!(Tr,(i00,i10,i11)); push!(Tr,(i00,i11,i01))
    else
        push!(Tr,(i00,i11,i10)); push!(Tr,(i00,i01,i11))
    end
    return nothing
end

@inline function _bool_keep(op::Symbol, inA::Bool, inB::Bool)
    op === :union        && return inA || inB
    op === :intersection && return inA && inB
    op === :difference   && return inA && !inB
    throw(ArgumentError("mesh_boolean: op must be :union, :intersection, or :difference (got $op)"))
end

function _bool_axis_aligned_boolean(A::Mesh, B::Mesh, op::Symbol)
    canon(v)=v==0 ? 0.0 : v
    xs = sort!(unique(canon.(vcat(A.coords[1,:], B.coords[1,:]))))
    ys = sort!(unique(canon.(vcat(A.coords[2,:], B.coords[2,:]))))
    zs = sort!(unique(canon.(vcat(A.coords[3,:], B.coords[3,:]))))
    nx=length(xs)-1; ny=length(ys)-1; nz=length(zs)-1
    _checked_mul3("mesh_boolean","axis-aligned classification cell",nx,ny,nz)
    dir = _unitn((1.0, 0.3141592653589793, 0.01720209895))       # generic ray, avoids grazing
    keep = falses(nx,ny,nz)
    @inbounds for i in 1:nx, j in 1:ny, k in 1:nz
        cx=_midpoint3(xs[i],xs[i+1]);cy=_midpoint3(ys[j],ys[j+1]);cz=_midpoint3(zs[k],zs[k+1])
        (xs[i]<cx<xs[i+1]&&ys[j]<cy<ys[j+1]&&zs[k]<cz<zs[k+1]) ||
            throw(ArgumentError("mesh_boolean: adjacent axis planes are below Float64 midpoint resolution"))
        inA=_inside_surface((cx,cy,cz),dir,A); inB=_inside_surface((cx,cy,cz),dir,B)
        keep[i,j,k]=_bool_keep(op,inA,inB)
    end
    @inline kept(i,j,k) = (1<=i<=nx && 1<=j<=ny && 1<=k<=nz) ? keep[i,j,k] : false
    V=Dict{NTuple{3,Float64},Int32}(); coords=NTuple{3,Float64}[]; Tr=NTuple{3,Int32}[]
    @inbounds for i in 1:nx, j in 1:ny, k in 1:nz
        keep[i,j,k] || continue
        kept(i-1,j,k) || _bool_emit_face!(V,coords,Tr, 1, xs[i],   ys[j],ys[j+1], zs[k],zs[k+1], 2,3, -1)
        kept(i+1,j,k) || _bool_emit_face!(V,coords,Tr, 1, xs[i+1], ys[j],ys[j+1], zs[k],zs[k+1], 2,3, +1)
        kept(i,j-1,k) || _bool_emit_face!(V,coords,Tr, 2, ys[j],   xs[i],xs[i+1], zs[k],zs[k+1], 1,3, -1)
        kept(i,j+1,k) || _bool_emit_face!(V,coords,Tr, 2, ys[j+1], xs[i],xs[i+1], zs[k],zs[k+1], 1,3, +1)
        kept(i,j,k-1) || _bool_emit_face!(V,coords,Tr, 3, zs[k],   xs[i],xs[i+1], ys[j],ys[j+1], 1,2, -1)
        kept(i,j,k+1) || _bool_emit_face!(V,coords,Tr, 3, zs[k+1], xs[i],xs[i+1], ys[j],ys[j+1], 1,2, +1)
    end
    isempty(Tr) && throw(ErrorException("mesh_boolean: $op of the two solids is empty"))
    C=Matrix{Float64}(undef,3,length(coords))
    @inbounds for (n,p) in enumerate(coords); C[1,n]=p[1]; C[2,n]=p[2]; C[3,n]=p[3]; end
    tm=Matrix{Int32}(undef,3,length(Tr))
    @inbounds for (n,t) in enumerate(Tr); tm[1,n]=t[1]; tm[2,n]=t[2]; tm[3,n]=t[3]; end
    return Mesh(C; tris=tm)
end

# ── SECTION C — general-position exact mesh Boolean (Cork/libigl-style) ─────────────

const _RB3 = NTuple{3,Rational{BigInt}}
@inline _bl_rat(p) = (Rational{BigInt}(p[1]), Rational{BigInt}(p[2]), Rational{BigInt}(p[3]))
@inline _bl_sub(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _bl_cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _bl_dot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

# exact rational point where edge (u,v) crosses plane {x : N·(x−P0)=0}. Caller
# guarantees the edge properly crosses (endpoints strictly on opposite sides).
function _bl_edge_plane(u::_RB3, v::_RB3, N::_RB3, P0::_RB3)
    den = _bl_dot(N, _bl_sub(v,u))
    den == 0 && throw(ErrorException("mesh_boolean: edge parallel to plane in a crossing pair (degeneracy)"))
    t = _bl_dot(N, _bl_sub(P0,u)) // den
    return (u[1]+t*(v[1]-u[1]), u[2]+t*(v[2]-u[2]), u[3]+t*(v[3]-u[3]))
end

# The (at most one) tri-tri intersection segment of fa=(a0,a1,a2), fb=(b0,b1,b2),
# given exact rational vertices. Returns nothing (no proper crossing), :degenerate
# (any orient3 sign zero / coplanar — unsupported), or (p,q)::Tuple{_RB3,_RB3}.
# Branch decisions use the EXACT orient3 predicate on the float points fa*,fb*.
function _bl_tritri(a0f,a1f,a2f, b0f,b1f,b2f, a0,a1,a2, b0,b1,b2)
    sb1 = orient3(a0f,a1f,a2f,b0f); sb2 = orient3(a0f,a1f,a2f,b1f); sb3 = orient3(a0f,a1f,a2f,b2f)
    sa1 = orient3(b0f,b1f,b2f,a0f); sa2 = orient3(b0f,b1f,b2f,a1f); sa3 = orient3(b0f,b1f,b2f,a2f)
    (sb1==0||sb2==0||sb3==0||sa1==0||sa2==0||sa3==0) && return :degenerate
    (sb1==sb2==sb3) && return nothing        # neither triangle crosses the other's plane
    (sa1==sa2==sa3) && return nothing
    Na = _bl_cross(_bl_sub(a1,a0), _bl_sub(a2,a0))
    Nb = _bl_cross(_bl_sub(b1,b0), _bl_sub(b2,b0))
    D  = _bl_cross(Na, Nb)
    (D[1]==0 && D[2]==0 && D[3]==0) && return :degenerate    # parallel but crossing ⇒ coplanar
    apts = _RB3[]; averts=((a0,sa1),(a1,sa2),(a2,sa3))
    for k in 1:3
        (u,su)=averts[k]; (v,sv)=averts[k%3+1]
        (su != sv) && push!(apts, _bl_edge_plane(u,v,Nb,b0))
    end
    bpts = _RB3[]; bverts=((b0,sb1),(b1,sb2),(b2,sb3))
    for k in 1:3
        (u,su)=bverts[k]; (v,sv)=bverts[k%3+1]
        (su != sv) && push!(bpts, _bl_edge_plane(u,v,Na,a0))
    end
    (length(apts)==2 && length(bpts)==2) || return :degenerate
    pa1=_bl_dot(D,apts[1]); pa2=_bl_dot(D,apts[2])
    alo,alopt,ahi,ahipt = pa1<=pa2 ? (pa1,apts[1],pa2,apts[2]) : (pa2,apts[2],pa1,apts[1])
    pb1=_bl_dot(D,bpts[1]); pb2=_bl_dot(D,bpts[2])
    blo,blopt,bhi,bhipt = pb1<=pb2 ? (pb1,bpts[1],pb2,bpts[2]) : (pb2,bpts[2],pb1,bpts[1])
    lo,lopt = alo>=blo ? (alo,alopt) : (blo,blopt)
    hi,hipt = ahi<=bhi ? (ahi,ahipt) : (bhi,bhipt)
    lo >= hi && return nothing               # intervals disjoint or touch at a point
    return (lopt, hipt)
end

# global point registry keyed on EXACT rational coordinate => identical geometric
# points (whatever pairing produced them) collapse to one id => watertight seam.
mutable struct _BLReg
    id::Dict{_RB3,Int32}
    X::Vector{Float64}; Y::Vector{Float64}; Z::Vector{Float64}
    RAT::Vector{_RB3}                     # exact rational coord per id (for on-edge tests)
end
_bl_reg() = _BLReg(Dict{_RB3,Int32}(), Float64[], Float64[], Float64[], _RB3[])
function _bl_gid!(R::_BLReg, rp::_RB3)
    id = get(R.id, rp, Int32(0))
    id != 0 && return id
    length(R.X)<typemax(Int32) ||
        throw(ArgumentError("mesh_boolean: exact point registry exceeds Int32"))
    push!(R.X, Float64(rp[1])); push!(R.Y, Float64(rp[2])); push!(R.Z, Float64(rp[3]))
    push!(R.RAT, rp)
    id = Int32(length(R.X)); R.id[rp] = id; return id
end
@inline _bl_pt(R::_BLReg,i) = (R.X[i], R.Y[i], R.Z[i])

# ── coplanar-region grouping (exact) ────────────────────────────────────────────
# A surface's triangles are grouped into maximal coplanar-connected regions (adjacent
# triangles sharing an edge whose four vertices are coplanar, exact orient3==0). Each
# region is retriangulated AS A WHOLE (its true outer boundary + chords), so a face's
# internal triangulation diagonal — an artifact that would otherwise split the shared
# intersection seam at a spurious point — never becomes a constraint.
struct _BLRegion
    faces::Vector{Int}                 # indices into the surface's face list
    bnd::Vector{NTuple{2,Int32}}       # region boundary edges (global-id pairs)
    N::_RB3                            # rational normal
end
@inline _bl_ekey(a,b) = a<=b ? (Int32(a),Int32(b)) : (Int32(b),Int32(a))

function _bl_regions(R::_BLReg, faces::Vector{NTuple{3,Int32}})
    nf=length(faces)
    e2f=Dict{NTuple{2,Int32},Vector{Int}}()
    for (fi,t) in enumerate(faces)
        for e in (_bl_ekey(t[1],t[2]),_bl_ekey(t[2],t[3]),_bl_ekey(t[1],t[3])); push!(get!(e2f,e,Int[]),fi); end
    end
    par=collect(1:nf)
    find(i)=(while par[i]!=i; par[i]=par[par[i]]; i=par[i]; end; i)
    for (e,fs) in e2f
        length(fs)==2 || continue
        t1=faces[fs[1]]; t2=faces[fs[2]]; u,v=e
        w1=t1[1]; (w1==u||w1==v)&&(w1=t1[2]); (w1==u||w1==v)&&(w1=t1[3])
        w2=t2[1]; (w2==u||w2==v)&&(w2=t2[2]); (w2==u||w2==v)&&(w2=t2[3])
        orient3(_bl_pt(R,u),_bl_pt(R,v),_bl_pt(R,w1),_bl_pt(R,w2))==0 && (par[find(fs[1])]=find(fs[2]))
    end
    comps=Dict{Int,Vector{Int}}()
    for fi in 1:nf; push!(get!(comps,find(fi),Int[]),fi); end
    out=_BLRegion[]
    for (_,fis) in comps
        inc=Dict{NTuple{2,Int32},Int}()
        for fi in fis; t=faces[fi]; for e in (_bl_ekey(t[1],t[2]),_bl_ekey(t[2],t[3]),_bl_ekey(t[1],t[3])); inc[e]=get(inc,e,0)+1; end; end
        bnd=NTuple{2,Int32}[]; for (e,c) in inc; c==1 && push!(bnd,e); end
        t0=faces[fis[1]]
        N=_bl_cross(_bl_sub(_bl_rat(_bl_pt(R,t0[2])),_bl_rat(_bl_pt(R,t0[1]))),
                    _bl_sub(_bl_rat(_bl_pt(R,t0[3])),_bl_rat(_bl_pt(R,t0[1]))))
        push!(out,_BLRegion(fis,bnd,N))
    end
    return out
end

# is global vertex v exactly on region boundary edge (a,b): collinear + between (exact)?
function _bl_on_edge(R::_BLReg, v::Int32, a::Int32, b::Int32)
    (v==a||v==b) && return true
    d=_bl_sub(R.RAT[b],R.RAT[a]); w=_bl_sub(R.RAT[v],R.RAT[a]); cr=_bl_cross(w,d)
    (cr[1]==0&&cr[2]==0&&cr[3]==0) || return false
    s=_bl_dot(w,d); L=_bl_dot(d,d); return 0<s<L
end

# Merge chords that meet at a spurious interior degree-2 collinear vertex (an
# intersection point landing on a face-internal diagonal — the two adjacent triangles
# of one coplanar face each contribute half of one straight chord). Real seam vertices
# are never degree-2-collinear (adjacent facets of the OTHER solid bend), so this only
# removes artifacts and keeps the seam identical on both surfaces.
function _bl_merge_collinear(R::_BLReg, reg::_BLRegion, chords::Vector{NTuple{2,Int32}}, origset::Set{Int32})
    onbnd(v) = any(e->_bl_on_edge(R,v,e[1],e[2]), reg.bnd)
    chords = copy(chords); alive=trues(length(chords))
    while true
        adj=Dict{Int32,Vector{Int}}()
        for ci in 1:length(chords)
            alive[ci] || continue
            a,b=chords[ci]; push!(get!(adj,a,Int[]),ci); push!(get!(adj,b,Int[]),ci)
        end
        merged=false
        for (v,cis) in adj
            (v in origset) && continue
            length(cis)==2 || continue
            onbnd(v) && continue
            (a1,b1)=chords[cis[1]]; (a2,b2)=chords[cis[2]]
            n1 = a1==v ? b1 : a1; n2 = a2==v ? b2 : a2
            (n1==n2||n1==v||n2==v) && continue
            cr=_bl_cross(_bl_sub(R.RAT[v],R.RAT[n1]), _bl_sub(R.RAT[n2],R.RAT[n1]))
            (cr[1]==0&&cr[2]==0&&cr[3]==0) || continue
            alive[cis[1]]=false; alive[cis[2]]=false
            push!(chords,_bl_ekey(n1,n2)); push!(alive,true)
            merged=true; break
        end
        merged || break
    end
    return NTuple{2,Int32}[chords[i] for i in 1:length(chords) if alive[i]]
end

# 2-D crossing-number point-in-polygon over the region boundary loop — keeps only
# sub-triangles whose centroid is inside a (possibly non-convex) region. A float test
# suffices here: sub-triangle centroids sit well away from the boundary edges.
function _bl_pip(cx, cy, bnd2::Vector{NTuple{4,Float64}})
    cnt=0
    for (x1,y1,x2,y2) in bnd2
        ((y1>cy) != (y2>cy)) || continue
        xint = x1 + (cy-y1)/(y2-y1)*(x2-x1)
        xint > cx && (cnt+=1)
    end
    return isodd(cnt)
end

# Retriangulate a whole coplanar region: pre-split boundary edges at on-edge chord
# points (a constructed on-edge point rounds ~1 ulp off the float line through the
# corners, so a single long edge constraint would make Mesh2D's exact-on-float
# predicates read a T-junction as a crossing — the vertex chain removes that), add
# the (collinear-merged) chords, CDT, keep interior triangles, lift + orient to N.
function _bl_retri_region(R::_BLReg, reg::_BLRegion, faces::Vector{NTuple{3,Int32}},
                          chords_in::Vector{NTuple{2,Int32}}, origset::Set{Int32})
    Nrat=reg.N
    chords = _bl_merge_collinear(R, reg, chords_in, origset)
    ids = Int32[]
    for (a,b) in reg.bnd; push!(ids,a); push!(ids,b); end
    for (a,b) in chords;  push!(ids,a); push!(ids,b); end
    ids = unique(ids)
    an = (abs(Nrat[1]),abs(Nrat[2]),abs(Nrat[3]))
    drop = an[1]>=an[2] ? (an[1]>=an[3] ? 1 : 3) : (an[2]>=an[3] ? 2 : 3)
    proj(i) = begin p=_bl_pt(R,i); drop==1 ? (p[2],p[3]) : (drop==2 ? (p[1],p[3]) : (p[1],p[2])) end
    edgeseg = NTuple{2,Int32}[]
    for (ci,cj) in reg.bnd
        pci=R.RAT[ci]; pcj=R.RAT[cj]; d=_bl_sub(pcj,pci); L=_bl_dot(d,d)
        onedge=Tuple{Rational{BigInt},Int32}[]
        for g in ids
            (g==ci||g==cj) && continue
            w=_bl_sub(R.RAT[g],pci); cr=_bl_cross(w,d)
            (cr[1]==0&&cr[2]==0&&cr[3]==0) || continue
            s=_bl_dot(w,d); (s>0 && s<L) || continue
            push!(onedge,(s//L,g))
        end
        sort!(onedge, by=x->x[1])
        chain=vcat(Int32[ci],Int32[g for (_,g) in onedge],Int32[cj])
        for k in 1:length(chain)-1; push!(edgeseg,(chain[k],chain[k+1])); end
    end
    g2local=Dict{Int32,Int}(); xs=Float64[]; ys=Float64[]; back=Dict{NTuple{2,Float64},Int32}()
    for (li,g) in enumerate(ids)
        p2=proj(g); push!(xs,p2[1]); push!(ys,p2[2]); g2local[g]=li; back[(p2[1],p2[2])]=g
    end
    segs=Tuple{Int,Int}[]
    for (p,q) in edgeseg; push!(segs,(g2local[p],g2local[q])); end
    for (p,q) in chords;  push!(segs,(g2local[p],g2local[q])); end
    T=constrained_delaunay(xs,ys,segs); m2=to_mesh(T)
    bnd2=NTuple{4,Float64}[]
    for (a,b) in reg.bnd; pa=proj(a); pb=proj(b); push!(bnd2,(pa[1],pa[2],pb[1],pb[2])); end
    Nf=(Float64(Nrat[1]),Float64(Nrat[2]),Float64(Nrat[3]))
    out=NTuple{3,Int32}[]
    for f in 1:size(m2.tris,2)
        gs=ntuple(k->begin nd=m2.tris[k,f]; back[(m2.coords[1,nd],m2.coords[2,nd])] end, 3)
        pc=((proj(gs[1])[1]+proj(gs[2])[1]+proj(gs[3])[1])/3, (proj(gs[1])[2]+proj(gs[2])[2]+proj(gs[3])[2])/3)
        _bl_pip(pc[1],pc[2],bnd2) || continue        # drop triangles outside a non-convex region
        p0=_bl_pt(R,gs[1]); p1=_bl_pt(R,gs[2]); p2=_bl_pt(R,gs[3])
        n=_bl_cross(_bl_sub(_bl_rat(p1),_bl_rat(p0)), _bl_sub(_bl_rat(p2),_bl_rat(p0)))
        (n[1]==0&&n[2]==0&&n[3]==0) && continue
        d=Float64(n[1])*Nf[1]+Float64(n[2])*Nf[2]+Float64(n[3])*Nf[3]
        push!(out, d>=0 ? (gs[1],gs[2],gs[3]) : (gs[1],gs[3],gs[2]))
    end
    return out
end

# AABB of a triangle (min,max) per axis, from float coords
@inline function _bl_tri_aabb(R,g0,g1,g2)
    p0=_bl_pt(R,g0);p1=_bl_pt(R,g1);p2=_bl_pt(R,g2)
    (min(p0[1],p1[1],p2[1]),min(p0[2],p1[2],p2[2]),min(p0[3],p1[3],p2[3])),
    (max(p0[1],p1[1],p2[1]),max(p0[2],p1[2],p2[2]),max(p0[3],p1[3],p2[3]))
end
@inline _bl_aabb_overlap(la,ha,lb,hb) =
    la[1]<=hb[1]&&lb[1]<=ha[1] && la[2]<=hb[2]&&lb[2]<=ha[2] && la[3]<=hb[3]&&lb[3]<=ha[3]

function _bool_general_boolean(A::Mesh, B::Mesh, op::Symbol)
    R = _bl_reg()
    gA = [ _bl_gid!(R, _bl_rat((A.coords[1,i],A.coords[2,i],A.coords[3,i]))) for i in 1:size(A.coords,2) ]
    gB = [ _bl_gid!(R, _bl_rat((B.coords[1,i],B.coords[2,i],B.coords[3,i]))) for i in 1:size(B.coords,2) ]
    faceA = [ (gA[A.tris[1,f]],gA[A.tris[2,f]],gA[A.tris[3,f]]) for f in 1:size(A.tris,2) ]
    faceB = [ (gB[B.tris[1,f]],gB[B.tris[2,f]],gB[B.tris[3,f]]) for f in 1:size(B.tris,2) ]
    aabbA=[_bl_tri_aabb(R,t...) for t in faceA]; aabbB=[_bl_tri_aabb(R,t...) for t in faceB]
    origset = Set{Int32}(vcat(gA, gB))
    regA = _bl_regions(R, faceA); regB = _bl_regions(R, faceB)
    regOfA = Dict{Int,Int}(); for (ri,rg) in enumerate(regA), fi in rg.faces; regOfA[fi]=ri; end
    regOfB = Dict{Int,Int}(); for (ri,rg) in enumerate(regB), fi in rg.faces; regOfB[fi]=ri; end
    chordsA = [NTuple{2,Int32}[] for _ in regA]; chordsB = [NTuple{2,Int32}[] for _ in regB]
    cutA = falses(length(regA)); cutB = falses(length(regB))
    for fa in 1:length(faceA)
        la,ha = aabbA[fa]; g=faceA[fa]
        a0=_bl_pt(R,g[1]);a1=_bl_pt(R,g[2]);a2=_bl_pt(R,g[3])
        a0r=_bl_rat(a0);a1r=_bl_rat(a1);a2r=_bl_rat(a2)
        for fb in 1:length(faceB)
            _bl_aabb_overlap(la,ha,aabbB[fb][1],aabbB[fb][2]) || continue
            h=faceB[fb]
            b0=_bl_pt(R,h[1]);b1=_bl_pt(R,h[2]);b2=_bl_pt(R,h[3])
            res = _bl_tritri(a0,a1,a2,b0,b1,b2, a0r,a1r,a2r, _bl_rat(b0),_bl_rat(b1),_bl_rat(b2))
            res === nothing && continue
            res === :degenerate && throw(ErrorException(
                "mesh_boolean: non-general-position degeneracy between A-face $fa and B-face $fb " *
                "(vertex/edge on the other's plane, or coplanar overlap) — unsupported by the general " *
                "tri-tri path; use axis-aligned inputs for coplanar cases"))
            p,q = res
            ip=_bl_gid!(R,p); iq=_bl_gid!(R,q)
            ip==iq && continue
            ra=regOfA[fa]; rb=regOfB[fb]
            push!(chordsA[ra],_bl_ekey(ip,iq)); push!(chordsB[rb],_bl_ekey(ip,iq))
            cutA[ra]=true; cutB[rb]=true
        end
    end
    subA = NTuple{3,Int32}[]
    for (ri,rg) in enumerate(regA)
        if !cutA[ri]; for fi in rg.faces; push!(subA, faceA[fi]); end
        else; append!(subA, _bl_retri_region(R, rg, faceA, chordsA[ri], origset)); end
    end
    subB = NTuple{3,Int32}[]
    for (ri,rg) in enumerate(regB)
        if !cutB[ri]; for fi in rg.faces; push!(subB, faceB[fi]); end
        else; append!(subB, _bl_retri_region(R, rg, faceB, chordsB[ri], origset)); end
    end
    dir = _unitn((1.0, 0.3141592653589793, 0.01720209895))
    cen(t) = ((_bl_pt(R,t[1])[1]+_bl_pt(R,t[2])[1]+_bl_pt(R,t[3])[1])/3,
              (_bl_pt(R,t[1])[2]+_bl_pt(R,t[2])[2]+_bl_pt(R,t[3])[2])/3,
              (_bl_pt(R,t[1])[3]+_bl_pt(R,t[2])[3]+_bl_pt(R,t[3])[3])/3)
    inB = [ _inside_surface(cen(t), dir, B) for t in subA ]     # classify vs ORIGINAL solids
    inA = [ _inside_surface(cen(t), dir, A) for t in subB ]
    keptA = NTuple{3,Int32}[]; keptB = NTuple{3,Int32}[]
    if op === :union
        for (i,t) in enumerate(subA); inB[i] || push!(keptA,t); end
        for (i,t) in enumerate(subB); inA[i] || push!(keptB,t); end
    elseif op === :intersection
        for (i,t) in enumerate(subA); inB[i] && push!(keptA,t); end
        for (i,t) in enumerate(subB); inA[i] && push!(keptB,t); end
    else # difference A∖B: A outside B, plus B inside A with reversed winding
        for (i,t) in enumerate(subA); inB[i] || push!(keptA,t); end
        for (i,t) in enumerate(subB); inA[i] && push!(keptB,(t[1],t[3],t[2])); end
    end
    allt = vcat(keptA, keptB)
    isempty(allt) && throw(ErrorException("mesh_boolean: $op of the two solids is empty"))
    used = sort!(unique(vcat([t[1] for t in allt],[t[2] for t in allt],[t[3] for t in allt])))
    length(used)<=typemax(Int32) || throw(ArgumentError("mesh_boolean: result node count exceeds Int32"))
    nid = Dict{Int32,Int32}(); for (k,g) in enumerate(used); nid[g]=Int32(k); end
    C=Matrix{Float64}(undef,3,length(used))
    for (k,g) in enumerate(used); p=_bl_pt(R,g); C[1,k]=p[1]; C[2,k]=p[2]; C[3,k]=p[3]; end
    tm=Matrix{Int32}(undef,3,length(allt))
    for (j,t) in enumerate(allt); tm[1,j]=nid[t[1]]; tm[2,j]=nid[t[2]]; tm[3,j]=nid[t[3]]; end
    return Mesh(C; tris=tm)
end

# ── SECTION D — public entry point ─────────────────────────────────────────────────

"""
    mesh_boolean(A::Mesh, B::Mesh, op::Symbol) -> Mesh

Native mesh-Boolean CSG of two closed, watertight, manifold triangulated solids
`A`, `B`. `op` ∈ (`:union`, `:intersection`, `:difference` = A∖B). Returns the
closed, watertight, outward-oriented surface `Mesh` of the result (fillable by
[`recover_boundary`](@ref)), or THROWS an explicit blocker naming an unsupported
degeneracy — never a silently wrong surface. Axis-aligned inputs (boxes / box-unions,
including every coplanar shared-face case) take an exact plane-arrangement path whose
volume is the exact sum of grid-cell volumes; general-position inputs (e.g.
box × cylinder) take a Cork/libigl-style exact tri-tri + Mesh2D-CDT path.
"""
function mesh_boolean(A::Mesh, B::Mesh, op::Symbol)
    (op === :union || op === :intersection || op === :difference) ||
        throw(ArgumentError("mesh_boolean: op must be :union, :intersection, or :difference (got $op)"))
    Ao = _bool_prepare(A, "A"); Bo = _bool_prepare(B, "B")
    r = (_bool_is_axis_aligned(Ao) && _bool_is_axis_aligned(Bo)) ?
        _bool_axis_aligned_boolean(Ao, Bo, op) : _bool_general_boolean(Ao, Bo, op)
    wt, mf, nb, mx = _bool_watertight(r)
    (wt && mf) || throw(ErrorException(
        "mesh_boolean: assembled $op result is not watertight/manifold " *
        "(boundary edges=$nb, max edge incidence=$mx) — refusing to return a leaky surface"))
    _require_surface3(r,"mesh_boolean result";oriented=true)
    return r
end


"""
    mesh_sized_conforming(surface::Mesh; hmax, inset=hmax, rng_seed=1, max_seeds=32) -> Mesh

CONFORMING tet mesh of the closed PLC `surface` **with interior size control**: a
background lattice of interior Steiner points (spacing `hmax/√3`, kept strictly
inside and inset ≥ `inset` from the boundary so the input facets stay recoverable)
is added, then Delaunay-filled, and accepted only if the exact `Rational{BigInt}`
conformity + validity gate passes — every input facet a tet face, all-positive,
closed-manifold. **Interior** edges are then `≤ hmax`; a boundary transition layer
of thickness ~`inset` stays at the input-surface resolution (pass a finely-
triangulated surface — e.g. from `MeshSurface` — for a fine boundary).

Robust for well-conditioned curved domains (spheres, generic surfaces). For
**maximally-cospherical** inputs (e.g. a fine axis-aligned cylinder's exact rings)
the lattice can break validity — the gate then **throws an explicit blocker**
rather than return an invalid mesh (never silent). For axis-aligned box assemblies
prefer [`mesh_box_regions`](@ref) (uniform `maxedge ≤ hmax`, exact).
"""
function mesh_sized_conforming(surface::Mesh; hmax::Real, inset::Real=hmax,
                               rng_seed::Integer=1, max_seeds::Integer=32)
    hm=_finite3(hmax,"mesh_sized_conforming","hmax")
    ins=_finite3(inset,"mesh_sized_conforming","inset")
    hm > 0 || throw(ArgumentError("mesh_sized_conforming: hmax must be positive (got $hmax)"))
    ins >= 0 || throw(ArgumentError("mesh_sized_conforming: inset must be non-negative (got $inset)"))
    _require_surface3(surface,"mesh_sized_conforming")
    base=_seed3(rng_seed,"mesh_sized_conforming")
    (1 <= max_seeds <= typemax(Int)) ||
        throw(ArgumentError("mesh_sized_conforming: max_seeds must be positive and fit Int (got $max_seeds)"))
    nseeds=Int(max_seeds)
    base <= typemax(Int)-(nseeds-1) ||
        throw(ArgumentError("mesh_sized_conforming: rng_seed + max_seeds overflows Int"))
    Px,Py,Pz,facets = _rb_dedup_surface(surface)
    length(Px) >= 4 || throw(ArgumentError("mesh_sized_conforming: need >= 4 distinct surface vertices"))
    isempty(facets) && throw(ArgumentError("mesh_sized_conforming: surface has no facets"))
    regions = _rb_build_regions(Px,Py,Pz,facets)
    S = _rb_crease_segments(regions)
    xlo=minimum(Px); xhi=maximum(Px); ylo=minimum(Py); yhi=maximum(Py); zlo=minimum(Pz); zhi=maximum(Pz)
    a = hm/sqrt(3.0)
    (isfinite(a)&&a>0) || throw(ArgumentError("mesh_sized_conforming: hmax is below Float64 spacing resolution"))
    sg = _raygrid(surface)
    nx=_ceil_count3((xhi-xlo)/a,"mesh_sized_conforming","x intervals")
    ny=_ceil_count3((yhi-ylo)/a,"mesh_sized_conforming","y intervals")
    nz=_ceil_count3((zhi-zlo)/a,"mesh_sized_conforming","z intervals")
    candidates=_checked_mul3("mesh_sized_conforming","lattice candidate",
                             _checked_add3("mesh_sized_conforming","x lattice",nx,1),
                             _checked_add3("mesh_sized_conforming","y lattice",ny,1),
                             _checked_add3("mesh_sized_conforming","z lattice",nz,1))
    candidates <= typemax(Int32)-length(Px) ||
        throw(ArgumentError("mesh_sized_conforming: lattice candidate count exceeds Int32 capacity"))
    nn=length(Px); Ix=Float64[]; Iy=Float64[]; Iz=Float64[]
    @inbounds for i in 0:nx, j in 0:ny, k in 0:nz
        px=xlo+i*(xhi-xlo)/nx; py=ylo+j*(yhi-ylo)/ny; pz=zlo+k*(zhi-zlo)/nz
        _inside_grid((px,py,pz),sg) || continue
        ok=true
        for v in 1:nn
            if hypot(px-Px[v],py-Py[v],pz-Pz[v]) < ins; ok=false; break; end
        end
        ok && (push!(Ix,px); push!(Iy,py); push!(Iz,pz))
    end
    lastreason="no seed produced a conforming valid sized mesh"; lastnrec=0
    for kk in 0:nseeds-1
        seed=base+kk
        T=delaunay3d(vcat(Px,Ix), vcat(Py,Iy), vcat(Pz,Iz); perturb=false, rng_seed=seed)
        keep=_classify_by_centroid(T,surface)
        for drop in (false,true)
            k2 = drop ? _rb_drop_flats(T,keep) : keep
            any(k2) || continue
            m=to_mesh3(T; keep=k2)
            ok,nrec,reason=_rb_gate(surface,m,regions,S,Px,Py,Pz,facets)
            ok && return m
            nrec>=lastnrec && (lastnrec=nrec; lastreason=reason)
        end
    end
    throw(ErrorException("mesh_sized_conforming: no conforming valid sized mesh after $nseeds seeds " *
        "($(lastnrec)/$(length(facets)) facets recovered — likely a cospherical-degenerate input; use " *
        "recover_boundary without interior size control, or mesh_box_regions for boxes). Blocker: $lastreason"))
end

end # module Mesh3D
