# `gmsh.model.occ.fuse` / `cut` / `intersect` / `fragment` —
# `OCC_Internals::booleanOperator`. The bounded implementation takes
# volume (dim-3) operands only: object/tool dim-tag lists materialize into
# a private scratch `GeoModel` through the closest-entities materializer
# (kernel tags carry verbatim, so sub-shapes shared between operands dedup
# through `seen`), the native analytic BREP engine runs there, and the
# resulting volumes publish back into the pending-shape registry.
# Disposition mirrors upstream `occBooleanPreserveNumbering`: a preserved
# operand keeps its kernel binding, a 1:1-replaced operand unbinds and
# republishes under its tag, deleted operands unbind recursively, and
# fresh pieces allocate tags after the surviving maximum in `outDimTags`
# order. Every record is staged before the first kernel mutation, so a
# publication failure cannot leave a half-applied registry. Engine
# limitation propagates loudly: operands must be materializable analytic
# solids (plane/cylinder/sphere/cone faces, line/circle sections), and
# `tag` <= 0 requests automatic numbering like upstream's negative tag
# (an explicit `tag` only binds single-piece results).

mutable struct _ModelOCCBooleanPub
    caller::String
    kernel::_ModelNativeOCCKernel
    sm::GeoModel
    matmap::Dict{Tuple{Int,Int},UInt64}
    cursor::Vector{Int}
    next_id::UInt64
    staged::Dict{UInt64,_ModelOCCShape}
    planned::Vector{Tuple{Int,Int,UInt64}}
    ids::Dict{Tuple{Int,Int},UInt64}
    ktags::Dict{Tuple{Int,Int},Int}
    # (dim, sm tag) -> operand tag the staged binding must rebind
    # (upstream's "replaced by a single shape" preserve-numbering claim)
    desired::Dict{Tuple{Int,Int},Int}
end

function _model_occ_boolean_id!(pub::_ModelOCCBooleanPub)
    pub.next_id=max(pub.next_id,pub.kernel.next_shape)
    pub.next_id<=typemax(UInt64)-1 || throw(ArgumentError(
        "$(pub.caller): OCC shape identity allocation overflow"))
    id=pub.next_id
    pub.next_id+=1
    return id
end

# Fresh-tag allocation draws on `cursor` alone: staging runs before any
# kernel mutation, so live `kernel.maxima` still holds pre-operation
# values — mixing it in would defeat the surviving-maximum seed the
# caller installs for volumes (index 6).
function _model_occ_boolean_tag!(pub::_ModelOCCBooleanPub,dim::Int)
    index=dim+3
    tag=pub.cursor[index]=pub.cursor[index]+1
    tag<=typemax(Int32) || throw(ArgumentError(
        "$(pub.caller): OCC entity tag allocation exceeds Int32"))
    return tag
end

# Reuse rule mirroring upstream: an sm entity materialized from a kernel
# sub-shape that is still bound reuses the binding outright; anything else
# (modified source, engine-created) gets a fresh tag — modified operands'
# sub-shapes do not keep their tags under preserve-numbering.
function _model_occ_boolean_reuse(pub::_ModelOCCBooleanPub,dim::Int,s::Int)
    kid=get(pub.matmap,(dim,s),UInt64(0))
    kid==0 && return UInt64(0)
    get(pub.kernel.shape_tags,kid,nothing)==(dim,s) || return UInt64(0)
    return kid
end

# A reused shape keeps its existing kernel binding — but staging runs
# before the unbind phase, so a reused sub-shape of an operand that is
# about to be unbound still looks bound here. Plan an idempotent re-bind
# of its existing record so the bind phase restores whatever the operand
# unbind cleared (a no-op when the binding survived).
function _model_occ_boolean_keep!(pub::_ModelOCCBooleanPub,dim::Int,s::Int,
                                  kid::UInt64)
    haskey(pub.staged,kid) && return nothing
    pub.staged[kid]=pub.kernel.shapes[kid]
    push!(pub.planned,(dim,s,kid))
    return nothing
end

# The engine can emit a split corner a few ulps off the position its
# siblings store exactly (typed-input section arithmetic differs from
# the untyped path). A face holding both readings is a few-ulps
# non-planar and the exact `orient3` checks then refuse it on the next
# materialization. Snap a fresh vertex's stored coordinates to a
# deterministic canonical — the lexicographically smallest coincident
# bound vertex, else the smallest coincident staged one — so every
# member of a dust cluster lands on the same bitwise value. Bound
# vertices are canonical candidates too, which keeps dirty engine
# output from ever entering the registry: no dusty bound vertex means
# all later snaps agree.
function _model_occ_boolean_snap_vertex(pub::_ModelOCCBooleanPub,
                                        p::NTuple{3,Float64})
    p=_model_occ_snap_coord.(p)
    atol=1e-11*max(1.0,abs(p[1]),abs(p[2]),abs(p[3]))
    # bound vertices are already committed registry geometry — a dusty
    # staged candidate must never win over them, or both coordinates end
    # up referenced by one face and exact coplanarity checks fail
    best=nothing
    for ((d,_),bid) in pub.kernel.bindings
        d==0 || continue
        c=pub.kernel.shapes[bid].geometry.point
        (abs(p[1]-c[1])<=atol && abs(p[2]-c[2])<=atol &&
         abs(p[3]-c[3])<=atol) || continue
        (best===nothing || c<best) && (best=c)
    end
    best===nothing || return best
    for sh in values(pub.staged)
        sh.dimension==0 || continue
        c=sh.geometry.point
        (abs(p[1]-c[1])<=atol && abs(p[2]-c[2])<=atol &&
         abs(p[3]-c[3])<=atol) || continue
        (best===nothing || c<best) && (best=c)
    end
    return best===nothing ? p : best
end

function _model_occ_boolean_vertex!(pub::_ModelOCCBooleanPub,s::Int,caller)
    skey=(0,s)
    haskey(pub.ids,skey) && return pub.ids[skey]
    kid=_model_occ_boolean_reuse(pub,0,s)
    if kid!=0
        _model_occ_boolean_keep!(pub,0,s,kid)
        pub.ids[skey]=kid
        return kid
    end
    id=_model_occ_boolean_id!(pub)
    tag=_model_occ_boolean_tag!(pub,0)
    kid=get(pub.matmap,skey,UInt64(0))
    mesh_size=kid==0 ? get(pub.sm.point_size,s,0.0) :
        pub.kernel.shapes[kid].geometry.mesh_size
    pub.staged[id]=_ModelOCCShape(0,
        (point=_model_occ_boolean_snap_vertex(pub,pub.sm.points[s]),
        mesh_size=mesh_size),UInt64[])
    push!(pub.planned,(0,tag,id))
    pub.ids[skey]=id
    return id
end

# `CurveOnSurface` slot convention: `fwd` is the edge's native-parameter
# pcurve and is always present; `rev` is the second (PCurve2) curve that
# only closed-surface seam edges carry. The boolean engine stores per-use
# slots instead — an edge used only reversed leaves `fwd` empty — so its
# lone entry normalizes back into `fwd` on publication.
_model_occ_boolean_pc_entry(e)=
    e.fwd===nothing ? (fwd=e.rev,rev=nothing) : e

# Remap a face record's `pcurves` keys from scratch curve tags to kernel
# edge tags (degenerate edges embed the record in their `face` field).
function _model_occ_boolean_remap_face(pub::_ModelOCCBooleanPub,face,caller)
    haskey(face,:pcurves) || return face
    pcs=Dict{Int,NamedTuple}()
    for (ct,pc) in face.pcurves
        _model_occ_boolean_edge!(pub,ct,caller)
        pcs[pub.ktags[(1,ct)]]=_model_occ_boolean_pc_entry(pc)
    end
    return merge(face,(pcurves=pcs,))
end

function _model_occ_boolean_edge!(pub::_ModelOCCBooleanPub,s::Int,caller)
    skey=(1,s)
    haskey(pub.ids,skey) && return pub.ids[skey]
    kid=_model_occ_boolean_reuse(pub,1,s)
    if kid!=0
        _model_occ_boolean_keep!(pub,1,s,kid)
        pub.ids[skey]=kid
        pub.ktags[skey]=s
        return kid
    end
    a,b=pub.sm.curves[s]
    aid=_model_occ_boolean_vertex!(pub,a,caller)
    bid=_model_occ_boolean_vertex!(pub,b,caller)
    # Untyped-default materialization leaves :line records empty; the
    # kernel edge they came from carries the verbatim record to publish.
    g=get(pub.sm.curve_geometry,s,nothing)
    if g===nothing
        kid=get(pub.matmap,(1,s),UInt64(0))
        if kid==0
            pa=pub.sm.points[a];pb=pub.sm.points[b]
            dl=sqrt((pb[1]-pa[1])^2+(pb[2]-pa[2])^2+(pb[3]-pa[3])^2)
            g=(occ=:line,origin=pa,
               dir=((pb[1]-pa[1])/dl,(pb[2]-pa[2])/dl,
                    (pb[3]-pa[3])/dl),t0=0.0,t1=dl)
        else
            g=pub.kernel.shapes[kid].geometry.geometry
        end
    end
    haskey(g,:face) && (g=merge(g,(face=_model_occ_boolean_remap_face(
        pub,g.face,caller),)))
    id=_model_occ_boolean_id!(pub)
    tag=_model_occ_boolean_tag!(pub,1)
    pub.staged[id]=_ModelOCCShape(1,(kind=get(pub.sm.curve_types,s,:line),
        geometry=g),UInt64[aid,bid])
    push!(pub.planned,(1,tag,id))
    pub.ids[skey]=id;pub.ktags[skey]=tag
    return id
end

# A wire's signed loop stores the kernel traversal for outer wires and its
# negated reversal for hole wires (the materializer's hole convention — the
# transform is an involution). Wires and shells are internal (-1/-2) shapes
# and always bind fresh, matching upstream's `maxTag` allocation.
function _model_occ_boolean_wire!(pub::_ModelOCCBooleanPub,loop_tag::Int,
                                  hole::Bool,caller)
    signed=pub.sm.loops[loop_tag]
    seq=hole ? Int[-s for s in Iterators.reverse(signed)] : signed
    uses=Tuple{UInt64,Bool}[]
    for sc in seq
        eid=_model_occ_boolean_edge!(pub,abs(sc),caller)
        push!(uses,(eid,sc<0))
    end
    id=_model_occ_boolean_id!(pub)
    tag=_model_occ_boolean_tag!(pub,-1)
    pub.staged[id]=_ModelOCCShape(-1,(uses=uses,),
        UInt64[u[1] for u in uses])
    push!(pub.planned,(-1,tag,id))
    return id
end

# Newell normal (unnormalized) + sample centroid of a scratch wire's 3D
# traversal. The face's natural normal is the wire's winding normal —
# `_brep_wire_pts`/`_brep_wire_newell` share the traversal sampling; the
# engine's fitted `axis` can disagree with it.
function _model_occ_boolean_wire_normal(sm::GeoModel,loop_tag::Int)
    pts=_brep_wire_pts(sm,loop_tag)
    n=(0.0,0.0,0.0);c=(0.0,0.0,0.0)
    for i in eachindex(pts)
        p=pts[i];q=pts[mod1(i+1,length(pts))]
        c=c.+p
        n=(n[1]+(p[2]-q[2])*(p[3]+q[3]),
           n[2]+(p[3]-q[3])*(p[1]+q[1]),
           n[3]+(p[1]-q[1])*(p[2]+q[2]))
    end
    return n,c./length(pts)
end

# 3-D sample points of a scratch wire's traversal (alias kept for the
# probe helpers' call sites).
_model_occ_boolean_wire_pts(sm::GeoModel,loop_tag::Int) =
    _brep_wire_pts(sm,loop_tag)

# A 3-D point strictly inside a scratch face's outer wire (and outside
# its hole wires) plus the face's unit winding normal. Concave faces put
# the polygon centroid outside the polygon, so candidates are edge
# midpoints nudged along the in-plane perpendicular until the projected
# point-in-polygon test accepts one.
function _model_occ_boolean_face_probe(sm::GeoModel,stag::Int,caller)
    loops=sm.surfaces[stag]
    nw,_=_model_occ_boolean_wire_normal(sm,loops[1])
    L=sqrt(_dot3(nw,nw))
    L==0.0 && throw(ErrorException("$caller: degenerate wire"))
    n=(nw[1]/L,nw[2]/L,nw[3]/L)
    k=argmax(abs.(collect(n)))
    pr=k==3 ? (1,2) : k==2 ? (1,3) : (2,3)
    outer=_model_occ_boolean_wire_pts(sm,loops[1])
    ox=[p[pr[1]] for p in outer];oy=[p[pr[2]] for p in outer]
    oidx=collect(1:length(outer))
    function inside(q)
        _point_in_polygon(q[pr[1]],q[pr[2]],ox,oy,oidx) || return false
        for lp in loops[2:end]
            h=_model_occ_boolean_wire_pts(sm,lp)
            _point_in_polygon(q[pr[1]],q[pr[2]],
                [p[pr[1]] for p in h],[p[pr[2]] for p in h],
                collect(1:length(h))) && return false
        end
        return true
    end
    for (i,p) in enumerate(outer)
        q=outer[mod1(i+1,length(outer))]
        e=(q[1]-p[1],q[2]-p[2],q[3]-p[3])
        el=sqrt(_dot3(e,e));el==0.0 && continue
        perp=_b3cross(n,(e[1]/el,e[2]/el,e[3]/el))
        mid=((p[1]+q[1])/2,(p[2]+q[2])/2,(p[3]+q[3])/2)
        for sgn in (1.0,-1.0),sc in (1e-3,1e-4,1e-5)
            cand=mid.+sgn*sc*el.*perp
            inside(cand) && return cand,n
        end
    end
    _,cw=_model_occ_boolean_wire_normal(sm,loops[1])
    inside(cw) && return cw,n
    throw(ErrorException(
        "$caller: no interior probe on boolean result face"))
end

function _model_occ_boolean_face!(pub::_ModelOCCBooleanPub,s::Int,caller)
    skey=(2,s)
    haskey(pub.ids,skey) && return pub.ids[skey]
    desired=get(pub.desired,skey,0)
    kid=desired==0 ? _model_occ_boolean_reuse(pub,2,s) : UInt64(0)
    if kid!=0
        _model_occ_boolean_keep!(pub,2,s,kid)
        pub.ids[skey]=kid
        pub.ktags[skey]=s
        return kid
    end
    children=UInt64[]
    for (index,loop_tag) in enumerate(pub.sm.surfaces[s])
        push!(children,_model_occ_boolean_wire!(pub,loop_tag,index>1,caller))
    end
    # Untyped-default materialization leaves :plane records empty; the
    # kernel face they came from carries the verbatim record to publish.
    g=get(pub.sm.surface_geometry,s,nothing)
    if g===nothing
        kid=get(pub.matmap,(2,s),UInt64(0))
        kid==0 && throw(ErrorException(
            "$caller: boolean result face lacks a geometry record"))
        g=pub.kernel.shapes[kid].geometry.geometry
    end
    if haskey(g,:pcurves)
        pcs=Dict{Int,NamedTuple}()
        for (ct,pc) in g.pcurves
            _model_occ_boolean_edge!(pub,ct,caller)
            pcs[pub.ktags[(1,ct)]]=_model_occ_boolean_pc_entry(pc)
        end
        g=merge(g,(pcurves=pcs,))
    end
    # The kernel face's `axis` is its natural normal — the wire's winding
    # direction — but the engine's record carries only a plane fit that
    # can disagree with the wire. When it does, flip `axis` to the
    # winding normal and toggle `reversed` so the surface parameter
    # normal (X×Y, unchanged) still round-trips the face normal.
    kind=get(pub.sm.surface_types,s,:plane)
    if kind===:plane && haskey(g,:axis)
        nw,_=_model_occ_boolean_wire_normal(pub.sm,pub.sm.surfaces[s][1])
        _dot3(nw,g.axis)<0 &&
            (g=merge(g,(axis=(.-g.axis),reversed=!get(g,:reversed,false))))
    end
    id=_model_occ_boolean_id!(pub)
    tag=desired==0 ? _model_occ_boolean_tag!(pub,2) : desired
    pub.staged[id]=_ModelOCCShape(2,(kind=get(pub.sm.surface_types,s,:plane),
        geometry=g),children)
    push!(pub.planned,(2,tag,id))
    pub.ids[skey]=id;pub.ktags[skey]=tag
    return id
end

# Whether the scratch face `stag` is used reversed in the result solid
# `sv`. The engine's surface-loop signs track provenance (which operand
# use a piece came from), not a consistent outward orientation — under
# untyped-default materialization the record axes don't reflect it
# either — so the use direction is measured directly: a probe offset
# along the outer wire's winding normal lands inside the solid exactly
# when that normal points into it.
function _model_occ_boolean_face_rev(pub::_ModelOCCBooleanPub,stag::Int,
                                     sv::NamedTuple,diag::Float64,
                                     sm_rev::Bool,caller)
    get(pub.sm.surface_types,stag,:plane)===:plane || return sm_rev
    c,n=_model_occ_boolean_face_probe(pub.sm,stag,caller)
    for scale in (1e-6,1e-7,1e-8)
        ε=scale*diag
        pin=pout=nothing
        try pin=_brep_point_in_solid(sv,c.+ε.*n,_BREP_TOL,caller)
        catch
        end
        try pout=_brep_point_in_solid(sv,c.-ε.*n,_BREP_TOL,caller)
        catch
        end
        pin===true && pout===false && return true
        pin===false && pout===true && return false
    end
    throw(ErrorException("$caller: cannot orient boolean result face"))
end

function _model_occ_boolean_shell!(pub::_ModelOCCBooleanPub,sl::Int,
                                   sv::NamedTuple,diag::Float64,caller)
    uses=Tuple{UInt64,Bool}[]
    for sf in pub.sm.surface_loops[sl]
        stag=abs(sf)
        fid=_model_occ_boolean_face!(pub,stag,caller)
        rev=_model_occ_boolean_face_rev(pub,stag,sv,diag,sf<0,caller)
        push!(uses,(fid,rev))
    end
    id=_model_occ_boolean_id!(pub)
    tag=_model_occ_boolean_tag!(pub,-2)
    pub.staged[id]=_ModelOCCShape(-2,(uses=uses,),
        UInt64[u[1] for u in uses])
    push!(pub.planned,(-2,tag,id))
    return id
end

# Publish one scratch volume's shape tree at kernel tag `ktag`.
function _model_occ_boolean_solid!(pub::_ModelOCCBooleanPub,v::Int,
                                   ktag::Int,caller)
    children=UInt64[]
    sv=_brep_extract(pub.sm,v,caller)
    lo=ntuple(i->minimum(p[i] for p in sv.pts),3)
    hi=ntuple(i->maximum(p[i] for p in sv.pts),3)
    diag=sqrt((hi[1]-lo[1])^2+(hi[2]-lo[2])^2+(hi[3]-lo[3])^2)
    diag==0.0 && (diag=1.0)
    for sl in pub.sm.volumes[v]
        push!(children,_model_occ_boolean_shell!(pub,sl,sv,diag,caller))
    end
    prim=nothing
    rec=get(pub.sm.cylinders,v,nothing)
    rec!==nothing && (prim=(occ=:cylinder,record=rec))
    rec===nothing && (rec=get(pub.sm.spheres,v,nothing))
    rec!==nothing && prim===nothing && (prim=(occ=:sphere,record=rec))
    rec===nothing && (rec=get(pub.sm.cones,v,nothing))
    rec!==nothing && prim===nothing && (prim=(occ=:cone,record=rec))
    geometry=prim===nothing ? (kind=:solid,) : (kind=:solid,primitive=prim)
    id=_model_occ_boolean_id!(pub)
    pub.staged[id]=_ModelOCCShape(3,geometry,children)
    push!(pub.planned,(3,ktag,id))
    pub.ids[(3,v)]=id
    return id
end

# Record the (dim, sm-tag) → kernel-shape-id materialization link for one
# operand's transitive closure (sm tags carry kernel tags verbatim).
function _model_occ_boolean_matmap!(matmap,kernel,id,caller,
                                    visited=Set{UInt64}())
    id in visited && return nothing
    push!(visited,id)
    shape=kernel.shapes[id]
    if shape.dimension>=0
        tag=_model_occ_shape_tag(kernel,id,caller)
        matmap[(shape.dimension,tag)]=id
    end
    for c in shape.children
        _model_occ_boolean_matmap!(matmap,kernel,c,caller,visited)
    end
    return nothing
end

function _model_occ_boolean_dim_tags(entries,caller)
    out=Tuple{Int,Int}[]
    for entry in entries
        d=_model_occ_dimension(first(entry),caller)
        t=_model_occ_requested_tag(last(entry),caller)
        push!(out,(d,t))
    end
    return out
end

function _model_occ_boolean(m::GeoModel,op::Symbol,objects,tools,
                            tag,removeObject,removeTool,caller)
    caller=String(caller)
    kernel=m.occ_kernel
    removeObject isa Bool || throw(ArgumentError(
        "$caller: removeObject must be Bool"))
    removeTool isa Bool || throw(ArgumentError(
        "$caller: removeTool must be Bool"))
    tag isa Integer && !(tag isa Bool) || throw(ArgumentError(
        "$caller: tag must be an integer"))
    requested=Int(tag)
    requested<=typemax(Int32) || throw(ArgumentError(
        "$caller: tag exceeds Int32"))
    # upstream treats a nonpositive tag as automatic numbering
    explicit=requested>0 ? requested : 0
    # an empty object list is a no-op returning success upstream
    isempty(objects) && return Tuple{Int,Int}[],Vector{Tuple{Int,Int}}[]
    odt=_model_occ_boolean_dim_tags(objects,caller)
    tdt=_model_occ_boolean_dim_tags(tools,caller)
    explicit!=0 && haskey(kernel.bindings,(first(odt[1]),explicit)) &&
        throw(ArgumentError("$caller: OpenCASCADE entity with tag " *
                            "$explicit already exists"))
    for (d,t) in Iterators.flatten((odt,tdt))
        (d==3 || d==2) || throw(ArgumentError(
            "$caller: OCC Boolean operands of dimension $d are not yet " *
            "supported (volumes and surfaces only)"))
        haskey(kernel.bindings,(d,t)) || throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension $d with " *
            "tag $t"))
    end
    obj_dt=odt;tool_dt=tdt
    in_dt=vcat(obj_dt,tool_dt)
    length(unique(in_dt))==length(in_dt) || throw(ArgumentError(
        "$caller: Boolean operands must be distinct entities"))
    # a lone operand under union/fragments/difference maps to itself
    # upstream (A∪∅, fragments of one arg, A−∅ all return A); skipping
    # materialization keeps operands the engine cannot yet
    # re-materialize — e.g. curved boolean results — valid in that
    # degenerate call. Intersection's empty-tool result is not the
    # operand, and explicit tags run through the engine's phantom-tag
    # path, so both stay off this shortcut.
    if length(in_dt)==1 && explicit==0 && op!==:intersection
        dt=Tuple{Int,Int}[only(in_dt)]
        return dt,Vector{Tuple{Int,Int}}[dt]
    end
    if any(d->d!=3,first.(in_dt))
        return _model_occ_boolean_mixed(m,op,obj_dt,tool_dt,requested,
                                        removeObject,removeTool,caller)
    end
    obj_tags=Int[t for (_,t) in odt]
    tool_tags=Int[t for (_,t) in tdt]
    in_tags=vcat(obj_tags,tool_tags)
    in_set=Set(in_tags)
    n=length(in_tags)
    # materialize operands into the scratch model — kernel tags carry over
    # verbatim and shared sub-shapes dedup through `seen`
    sm=GeoModel()
    seen=Set{UInt64}()
    matmap=Dict{Tuple{Int,Int},UInt64}()
    snap=NTuple{3,Float64}[]
    for t in in_tags
        id=kernel.bindings[(3,t)]
        _model_occ_dist_materialize!(sm,kernel,3,id,seen,caller;snap=snap,
            untyped_defaults=true)
        _model_occ_boolean_matmap!(matmap,kernel,id,caller)
    end
    # the engine's `_alloc_tag!` draws on `next_tag`; seed it above the
    # materialized kernel tags or fresh entities would collide
    sm.next_tag[1]=max(sm.next_tag[1],maximum(keys(sm.points);init=0))
    sm.next_tag[2]=max(sm.next_tag[2],maximum(keys(sm.curves);init=0))
    sm.next_tag[3]=max(sm.next_tag[3],maximum(keys(sm.surfaces);init=0))
    sm.next_tag[4]=max(sm.next_tag[4],maximum(keys(sm.volumes);init=0))
    image=Vector{Int}[]
    out=boolean_volumes_multi!(sm,op,obj_tags,tool_tags;
        tag=explicit,remove_object=removeObject,remove_tool=removeTool,
        caller=caller,image_out=image)
    # explicit-tag corner: the sole result is an unmodified operand kept in
    # place — upstream's `_multiBind` reports its own tag rather than the
    # requested one; the owning operand is the one mapping to itself
    phantom_t=0
    if explicit!=0 && !haskey(sm.volumes,explicit) && explicit in out
        owners=Int[t for (i,t) in enumerate(in_tags) if t in image[i]]
        length(owners)==1 || throw(ErrorException(
            "$caller: cannot resolve the explicit-tag result operand"))
        phantom_t=only(owners)
    end
    resolve(v)=v==explicit && phantom_t!=0 ? phantom_t : v
    rvmap=Dict{Int,Int}()
    for v in Iterators.flatten((out,Iterators.flatten(image)))
        rvmap[v]=resolve(v)
    end
    needs_publish(sv)=haskey(sm.booleans,sv) ||
        (sv==explicit && explicit!=0 && !(sv in in_set))
    desired(sv)=sv in in_set || sv==explicit ? sv : 0
    # stage every result volume's shape tree before the first kernel
    # mutation — publication order follows `out` (upstream's exploration
    # order), so fresh tags land in `outDimTags` order. Fresh volume tags
    # start after the surviving maximum like upstream's post-removal
    # `maxTag`: surviving volumes are the non-operand bindings plus the
    # operand tags that stay bound (kept or rebound), since a volume's
    # binding has no owning parent to protect it. Sub-entity fresh tags
    # stay on the pre-op maximum — shared sub-shapes can survive an
    # operand unbind through ownership, so that bound is conservative.
    surviving3=0
    for (d,t) in keys(kernel.bindings)
        d==3 && !(t in in_set) && (surviving3=max(surviving3,t))
    end
    for t in in_tags
        haskey(sm.volumes,t) && (surviving3=max(surviving3,t))
    end
    cursor=collect(kernel.maxima);cursor[6]=surviving3
    pub=_ModelOCCBooleanPub(caller,kernel,sm,matmap,
        cursor,kernel.next_shape,
        Dict{UInt64,_ModelOCCShape}(),Tuple{Int,Int,UInt64}[],
        Dict{Tuple{Int,Int},UInt64}(),Dict{Tuple{Int,Int},Int}(),
        Dict{Tuple{Int,Int},Int}())
    ktagmap=Dict{Int,Int}()
    for v in out
        sv=rvmap[v]
        (haskey(sm.volumes,sv) && needs_publish(sv)) || continue
        haskey(ktagmap,sv) && continue
        ktag=desired(sv)
        ktag==0 && (ktag=_model_occ_boolean_tag!(pub,3))
        ktagmap[sv]=ktag
        _model_occ_boolean_solid!(pub,sv,ktag,caller)
    end
    # kernel mutation phase: unbind deleted/rebound operands recursively,
    # then bind the staged shapes — upstream's remove-then-bind order
    for t in in_tags
        (haskey(sm.volumes,t) && !haskey(sm.booleans,t)) && continue
        id=get(kernel.bindings,(3,t),nothing)
        id===nothing || _model_occ_unbind!(kernel,id,true)
    end
    for (dim,ktag,id) in pub.planned
        _model_occ_bind_shape!(kernel,dim,ktag,id,pub.staged[id])
    end
    isempty(pub.planned) ||
        (kernel.next_shape=max(kernel.next_shape,pub.next_id))
    _model_occ_collect_shapes!(kernel)
    out_dt=Tuple{Int,Int}[(3,get(ktagmap,rvmap[v],rvmap[v])) for v in out]
    out_map=Vector{Tuple{Int,Int}}[
        Tuple{Int,Int}[(3,get(ktagmap,rvmap[v],rvmap[v])) for v in im]
        for im in image]
    return out_dt,out_map
end

# `gmsh.model.occ.remove_all_duplicates` — upstream
# `OCC_Internals::removeAllDuplicates`, literally
# `booleanFragments(-1, allMaxDimEntities, {}, remove both)`: every
# bound entity at the maximum dimension fragments against its peers
# (`_getAllDimTags(dimTags, _getMaxDim())`), so overlapping solids split
# into shared cells and exact duplicates collapse into one. Non-volume
# maximum dimensions stay outside the bridge's volume-operand bound and
# reject loudly like a direct `fragment` call would.
function _model_occ_remove_all_duplicates(m::GeoModel,caller)
    kernel=m.occ_kernel
    maxdim=-1
    for (d,t) in keys(kernel.bindings)
        d>=0 && (maxdim=max(maxdim,d))
    end
    maxdim<0 && return nothing      # nothing bound — upstream no-ops
    ents=Tuple{Int,Int}[]
    for (d,t) in keys(kernel.bindings)
        d==maxdim && push!(ents,(d,t))
    end
    sort!(ents)
    length(ents)<2 && return nothing
    _model_occ_boolean(m,:fragments,ents,Int[],-1,true,true,caller)
    return nothing
end

# ════════════════════════════════════════════════════════════════════════════
# Dimension-2 operands — surfaces in OCC Booleans
#
# A surface operand runs through the same analytic arrangement as volumes: it
# is extracted into a `_brep_extract`-shaped rep whose single face has no
# enclosing shell (orev=+1), section/imprint machinery is shared verbatim, and
# its domain regions classify like volume face regions with two differences:
#   * point-in-solid containment is only meaningful against volume operands —
#     a sheet owns no inside cell, so `cont` and the cell masks only ever take
#     volume bits;
#   * a region inside volumes becomes a partition wall (uses on both sides of
#     the containing cell), a region inside nothing stays a standalone result
#     face, and a coincident region follows the same lowest-index-owner rule as
#     coincident volume faces — non-owners are credited geometrically at
#     materialization time (`_occ_face_credits`).
#
# Cells are found by connectivity, not by mask: at each result edge the
# incident faces' uses pair by angular fan so a wall that reaches the
# boundary on a closed loop splits the containing cell while a free rim
# merges its two sides. Cells with no volume membership never materialize as
# solids — faces in them surface as free (standalone) result faces, matching
# upstream `_multiBind`'s "solids then free faces" emission order.
# ════════════════════════════════════════════════════════════════════════════

# `_brep_extract` for a bare surface entity: one face, no shell, orev=+1.
function _occ_sheet_rep(m::GeoModel,tag::Int,caller)
    haskey(m.surfaces,tag) || throw(ArgumentError(
        "$caller: unknown Surface[$tag]"))
    pts=NTuple{3,Float64}[];ptags=Int[]
    edges=NamedTuple[];ctags=Int[]
    pmap=Dict{Int,Int}();emap=Dict{Int,Int}()
    function vertex(p)
        get!(pmap,p) do
            push!(pts,m.points[p]);push!(ptags,p);length(pts)
        end
    end
    stag=tag;orev=1
    sg=get(m.surface_geometry,stag,nothing)
    occ_face=sg!==nothing && hasproperty(sg,:occ)
    occ_face || (sg=_brep_plane_record(m,stag,caller))
    wires=Vector{Tuple{Int,Int,NamedTuple}}[]
    for loop_tag in m.surfaces[stag]
        wire=Tuple{Int,Int,NamedTuple}[]
        for sc in m.loops[loop_tag]
            ctag=abs(sc);sgn=sc>0 ? 1 : -1
            eidx=get!(emap,ctag) do
                curve,t0,t1=_brep_edge_curve(m,ctag)
                a,b=m.curves[ctag]
                push!(edges,(curve=curve,t0=t0,t1=t1,
                             v1=vertex(a),v2=vertex(b)))
                push!(ctags,ctag);length(edges)
            end
            local pc
            if occ_face
                ent=get(sg.pcurves,ctag,nothing)
                ent===nothing && throw(ArgumentError(
                    "$caller: Surface[$stag] lacks a pcurve for " *
                    "Curve[$ctag]"))
                pc=(sgn<0 && ent.rev!==nothing) ? ent.rev : ent.fwd
                pc===nothing && throw(ArgumentError(
                    "$caller: Surface[$stag] lacks a pcurve for " *
                    "Curve[$ctag]"))
                erec=get(m.curve_geometry,ctag,nothing)
                if erec===nothing || !hasproperty(erec,:occ)
                    e=edges[eidx]
                    Δ=e.t1-e.t0
                    pc=_occ_pcurve_reparam(pc,-e.t0/Δ,1.0/Δ)
                end
            else
                e=edges[eidx]
                pc=_brep_synth_pcurve(sg,e.curve,e.t0,e.t1,caller)
            end
            push!(wire,(eidx,sgn,pc))
        end
        push!(wires,wire)
    end
    isempty(wires) && throw(ArgumentError(
        "$caller: Surface[$stag] has no boundary wires"))
    return (pts=pts,ptags=ptags,edges=edges,ctags=ctags,
            faces=[(sg=sg,stag=stag,orev=orev,wires=wires)])
end

# `uses` entry of a kept spec: (membership mask of the adjacent cell, sign).
# Sign convention follows the volume path: cell on the −sign·record-normal
# side of the face — for a solid face orev=+1 the inside cell takes +1.
#
# `_brep_collect_faces_n` over a mixed operand list. Returns
# (kept, sheet_out, touched): `kept` carries cell-boundary specs (uses
# non-empty), `sheet_out` the standalone face pieces (uses empty),
# `touched` the geometrically-modified operand flags.
function _occ_mixed_collect(op::Symbol,reps,isvol::Vector{Bool},nobj::Int,
                            caller)
    tol=_BREP_TOL
    N=length(reps)
    N<=62 || throw(ArgumentError(
        "$caller: at most 62 Boolean operands per operation"))
    coinc=NTuple{4,Int}[]
    sections=NamedTuple[]
    touched=falses(N)
    for i in 1:N, j in i+1:N
        sheets=!isvol[i] && !isvol[j]
        for (fi,fI) in enumerate(reps[i].faces),
            (fj,fJ) in enumerate(reps[j].faces)
            if fI.sg.occ===:plane && fJ.sg.occ===:plane &&
                    _brep_planes_coincident(fI.sg,fJ.sg)
                push!(coinc,(i,fi,j,fj))
                continue
            end
            # non-coplanar sheet×sheet sections split both domains only
            # under fragments/fuse — upstream GF imprints the section edge
            # on each face there, while cut returns the object whole and
            # intersect's bare curve filters out below minDim
            sheets && (op===:difference || op===:intersection) && continue
            for s in _brep_face_sections(fI,fJ,caller)
                spansI=_brep_trim_section(reps[i],fI,s.pcA,s.closed,tol)
                spansJ=_brep_trim_section(reps[j],fJ,s.pcB,s.closed,tol)
                merged=_brep_intersect_spans(spansI,spansJ,s.closed,tol,
                                             caller)
                isempty(merged) && continue
                push!(sections,(curve=s.curve,closed=s.closed,pcA=s.pcA,
                                pcB=s.pcB,oi=i,fa=fi,oj=j,fb=fj,
                                spans=merged))
            end
        end
    end
    # operand edges piercing other operands' faces split those edges
    esplits=[Dict{Int,Vector{Float64}}() for _ in 1:N]
    for i in 1:N
        Si=reps[i];es=esplits[i]
        for (eidx,e) in enumerate(Si.edges)
            e.curve.kind===:degenerate && continue
            ps=Float64[]
            for j in 1:N
                j==i && continue
                for of in reps[j].faces
                    append!(ps,_brep_edge_face_params(reps[j],e,of,tol))
                end
            end
            isempty(ps) || (es[eidx]=sort!(ps))
        end
    end
    sec_cuts=Dict{Int,Vector{Float64}}()
    imps=[[Dict{Int,Vector{NamedTuple}}() for _ in reps[i].faces]
          for i in 1:N]
    for (i,fi,j,fj) in coinc
        _brep_coinc_imprints_n!(reps,i,fi,j,fj,esplits,sec_cuts,
            get!(Vector{NamedTuple},imps[i][fi],j),sections,tol)
        _brep_coinc_imprints_n!(reps,j,fj,i,fi,esplits,sec_cuts,
            get!(Vector{NamedTuple},imps[j][fj],i),sections,tol)
    end
    esec_imps=[[NamedTuple[] for _ in reps[i].faces] for i in 1:N]
    esec_skip=falses(length(sections))
    for (si,sec) in enumerate(sections)
        owners=Int[];owned=Int[]
        for k in (sec.oi,sec.oj)
            S=reps[k]
            efound=0
            for wire in S.faces[k==sec.oi ? sec.fa : sec.fb].wires,
                (eidx,_,_) in wire
                if _brep_same_curve3d(sec.curve,S.edges[eidx].curve,tol)
                    efound=eidx;break
                end
            end
            efound==0 || (push!(owners,k);push!(owned,efound))
        end
        isempty(owners) && continue
        esec_skip[si]=true
        for k in 1:N
            (k==sec.oi || k==sec.oj) && continue
            S=reps[k];efound=0
            for tface in S.faces, wire in tface.wires,(eidx,_,_) in wire
                e=S.edges[eidx]
                _brep_same_curve3d(sec.curve,e.curve,tol) || continue
                covers=any(sec.spans) do (lo,hi)
                    pm=_brep_curve_eval(sec.curve,(lo+hi)/2)
                    tp=_brep_curve_param(e.curve,pm,(e.t0+e.t1)/2)
                    tp!==nothing && tp>=min(e.t0,e.t1)-tol &&
                        tp<=max(e.t0,e.t1)+tol
                end
                covers && (efound=eidx;break)
            end
            efound==0 || (push!(owners,k);push!(owned,efound))
        end
        wpos=argmin(owners);w=owners[wpos];weidx=owned[wpos]
        we=reps[w].edges[weidx]
        targets=Tuple{Int,Int,Any}[]
        for k in (sec.oi,sec.oj)
            k==w && continue
            kf=k==sec.oi ? sec.fa : sec.fb
            kpc=k==sec.oi ? sec.pcA : sec.pcB
            push!(targets,(k,kf,kpc))
            kp=findfirst(==(k),owners)
            if kp!==nothing
                ke=owned[kp]
                for (tf,tface) in enumerate(reps[k].faces)
                    tf==kf && continue
                    any(wr->any(x->x[1]==ke,wr),tface.wires) || continue
                    push!(targets,(k,tf,nothing))
                end
            end
        end
        for (k2,oe) in zip(owners,owned)
            (k2==w || k2==sec.oi || k2==sec.oj) && continue
            for (tf,tface) in enumerate(reps[k2].faces)
                any(wr->any(x->x[1]==oe,wr),tface.wires) || continue
                push!(targets,(k2,tf,nothing))
            end
        end
        for (lo,hi) in sec.spans
            subs=_brep_esec_tspans(we,sec.curve,lo,hi,tol)
            subs===nothing && throw(ArgumentError(
                "$caller: a coincident section endpoint does not project " *
                "onto its operand edge"))
            for (sa,sb,ta,tb) in subs
                p1=_brep_curve_eval(sec.curve,sa)
                p2=_brep_curve_eval(sec.curve,sb)
                wlo,whi=minmax(ta,tb)
                push!(get!(Vector{Float64},esplits[w],weidx),wlo,whi)
                for (k2,oe) in zip(owners,owned)
                    k2==w && continue
                    oer=reps[k2].edges[oe]
                    omid=(oer.t0+oer.t1)/2
                    l1=_brep_curve_param(oer.curve,p1,omid)
                    l2=_brep_curve_param(oer.curve,p2,omid)
                    l1===nothing || l2===nothing || push!(get!(
                        Vector{Float64},esplits[k2],oe),min(l1,l2),
                        max(l1,l2))
                end
                for (tk,tf,spc) in targets
                    tface=reps[tk].faces[tf]
                    pc,neg = if spc!==nothing
                        α=(sb-sa)/(tb-ta);β=sa-α*ta
                        _brep_pc_reparam(spc,α,β,tol,caller)
                    else
                        _brep_pc_on_surface(tface.sg,we.curve,wlo,whi,
                                            tol,caller)
                    end
                    sslo,sshi=neg ? (-whi,-wlo) : (wlo,whi)
                    _brep_imprint_contacts_n!(reps,tk,tf,w,weidx,pc,neg,
                                              sslo,sshi,esplits,sec_cuts,
                                              sections,tol,si)
                    for te in (sslo,sshi)
                        uv=_brep_pc_eval(pc,te)
                        for wire in tface.wires,(e2,_,epc2) in wire
                            e2r=reps[tk].edges[e2]
                            lo2,hi2=minmax(e2r.t0,e2r.t1)
                            tp=_brep_pc_param(epc2,uv,lo2,hi2,10tol)
                            (tp===nothing || tp<=lo2+tol || tp>=hi2-tol) &&
                                continue
                            push!(get!(Vector{Float64},esplits[tk],e2),tp)
                        end
                    end
                    ex_i=findfirst(es->es.wside==w && es.weidx==weidx &&
                                       es.neg==neg && es.pc==pc,
                                   esec_imps[tk][tf])
                    if ex_i===nothing
                        push!(esec_imps[tk][tf],
                              (pc=pc,spans=[(wlo,whi)],neg=neg,
                               wside=w,weidx=weidx,crv=we.curve))
                    else
                        es=esec_imps[tk][tf][ex_i]
                        any(sp->abs(sp[1]-wlo)<tol && abs(sp[2]-whi)<tol,
                            es.spans) || push!(es.spans,(wlo,whi))
                    end
                end
            end
        end
    end
    sec_edges=NamedTuple[]
    secs_pf=[[NamedTuple[] for _ in reps[i].faces] for i in 1:N]
    for (si,sec) in enumerate(sections)
        esec_skip[si] && continue
        for (lo,hi) in sec.spans
            cuts=copy(get(sec_cuts,si,Float64[]))
            for k in (sec.oi,sec.oj)
                S=reps[k]
                face=k==sec.oi ? S.faces[sec.fa] : S.faces[sec.fb]
                spc=k==sec.oi ? sec.pcA : sec.pcB
                uper=_brep_uperiod(face.sg)
                uper>0 || continue
                ulo=_brep_face_ubase(S.edges,face)
                pieces=_brep_pc_split_all_useams(spc,lo,hi,uper,ulo,tol)
                for p in pieces[2:end]
                    push!(cuts,p[1])
                end
            end
            for (s0,s1) in _brep_split_ranges(lo,hi,cuts,tol)
                push!(sec_edges,(curve=sec.curve,t0=s0,t1=s1))
                se=length(sec_edges)
                push!(secs_pf[sec.oi][sec.fa],
                      (se=se,span=(s0,s1),pc=sec.pcA,crv=sec.curve))
                push!(secs_pf[sec.oj][sec.fb],
                      (se=se,span=(s0,s1),pc=sec.pcB,crv=sec.curve))
            end
        end
    end
    for i in 1:N
        es=esplits[i];S=reps[i]
        for (fidx,f) in enumerate(S.faces), sf in secs_pf[i][fidx]
            uper=_brep_uperiod(f.sg)
            ulo=_brep_face_ubase(S.edges,f)
            for tend in (sf.span[1],sf.span[2])
                uv=_brep_wrap_u(_brep_pc_eval(sf.pc,tend),ulo,uper)
                for wire in f.wires,(eidx,sign,epc) in wire
                    e=S.edges[eidx]
                    lo,hi=minmax(e.t0,e.t1)
                    t=_brep_pc_param(epc,uv,lo,hi,10tol)
                    t===nothing && continue
                    (t>lo+tol && t<hi-tol) || continue
                    push!(get!(Vector{Float64},es,eidx),t)
                end
            end
        end
    end
    coinc_by_face=[Dict{Int,Vector{Tuple{Int,Int}}}() for _ in 1:N]
    for (i,fi,j,fj) in coinc
        push!(get!(Vector{Tuple{Int,Int}},coinc_by_face[i],fi),(j,fj))
        push!(get!(Vector{Tuple{Int,Int}},coinc_by_face[j],fj),(i,fi))
    end
    kept=NamedTuple[]
    sheet_out=NamedTuple[]
    # an operand is *modified* (upstream's mapModified) only when a face's
    # resulting domain actually differs: more than one region, a dropped
    # region, or region wires that are not exactly the original full edges.
    # Imprint curves drawn through a face without cutting it — e.g. a
    # partner's rim outside its domain — change nothing and must not set
    # `touched`, or an untouched operand would be needlessly restaged.
    wires_original(S,arcs,cycles,region)=all(
        c->all(pair->begin
                ar=arcs[pair[1]]
                (ar.key isa Tuple && length(ar.key)>=4 &&
                 ar.key[1]===:e) || return false
                ei=findfirst(==(ar.key[2]),S.ctags)
                ei===nothing && return false
                e=S.edges[ei];lo,hi=minmax(e.t0,e.t1)
                ar.key[3]==lo && ar.key[4]==hi
            end, cycles[c]), region)
    for i in 1:N
        S=reps[i];es=esplits[i]
        extras=[begin
                    f=S.faces[fidx]
                    ea=NamedTuple[]
                    by_winner=Dict{Int,Vector{NamedTuple}}()
                    for spec in esec_imps[i][fidx]
                        push!(get!(Vector{NamedTuple},by_winner,spec.wside),
                              spec)
                    end
                    for w in sort!(collect(keys(by_winner)))
                        append!(ea,_brep_esec_arcs(reps[w],S,f,
                                                 by_winner[w],esplits[w],
                                                 w,tol))
                    end
                    arcs=copy(ea)
                    for j in sort!(collect(keys(imps[i][fidx])))
                        append!(arcs,_brep_imprint_arcs(
                            S,reps[j],f,imps[i][fidx][j],esplits[j],
                            secs_pf[i][fidx],j,arcs,tol))
                    end
                    arcs
                end for fidx in eachindex(S.faces)]
        for (fidx,f) in enumerate(S.faces)
            regions,arcs,cycles=_brep_face_regions(
                S,f,i,es,secs_pf[i][fidx],extras[fidx],tol,caller)
            (length(regions)==1 &&
             wires_original(S,arcs,cycles,regions[1])) ||
                (touched[i]=true)
            cofaces=get(coinc_by_face[i],fidx,Tuple{Int,Int}[])
            for region in regions
                uv=_brep_region_sample(arcs,cycles,region,tol,caller)
                P=_brep_surface_eval(f.sg,uv[1],uv[2])
                same_p=Int[];opp_p=Int[]
                for (oj,ofidx) in cofaces
                    of=reps[oj].faces[ofidx]
                    uv2=_brep_surface_uv(of.sg,P)
                    _brep_uv_in_face(reps[oj].edges,of,uv2,tol) || continue
                    sameside=(f.orev*of.orev*_b3dot(f.sg.axis,of.sg.axis))>0
                    push!(sameside ? same_p : opp_p,oj)
                end
                coincident=!isempty(same_p)||!isempty(opp_p)
                # containing volume operands only — sheets own no inside
                # cell, and coincident partners classify through same_p/opp_p
                cont=Int[]
                for k in 1:N
                    k==i && continue
                    !isvol[k] && continue
                    (k in same_p || k in opp_p) && continue
                    _brep_point_in_solid(reps[k],P,tol,caller) || continue
                    push!(cont,k)
                end
                local lin,lout,linL,loutL
                if isvol[i]
                    lin=_brep_bit(i);lout=UInt64(0)
                    linL=lin;loutL=lout
                    for oj in same_p
                        isvol[oj] && (lin|=_brep_bit(oj))
                        linL|=_brep_bit(oj)
                    end
                    for oj in opp_p
                        isvol[oj] && (lout|=_brep_bit(oj))
                        loutL|=_brep_bit(oj)
                    end
                    for k in cont
                        b=_brep_bit(k);lin|=b;lout|=b;linL|=b;loutL|=b
                    end
                else
                    # a sheet's own membership is material, not a cell —
                    # masks take only the adjacent volume cells while the
                    # keep labels see the full operand set
                    lin=UInt64(0);lout=UInt64(0)
                    for oj in same_p
                        isvol[oj] && (lin|=_brep_bit(oj))
                    end
                    for oj in opp_p
                        isvol[oj] && (lout|=_brep_bit(oj))
                    end
                    for k in cont
                        b=_brep_bit(k);lin|=b;lout|=b
                    end
                    linL=_brep_bit(i)|lin;loutL=_brep_bit(i)|lout
                    for oj in vcat(same_p,opp_p)
                        b=_brep_bit(oj);linL|=b;loutL|=b
                    end
                end
                keep_in=_brep_label_keep(op,linL,nobj)
                keep_out=_brep_label_keep(op,loutL,nobj)
                uses=Tuple{UInt64,Int}[]
                if coincident
                    if minimum(vcat(i,same_p,opp_p))!=i
                        touched[i]=true;continue
                    end
                    if !isempty(opp_p) &&
                            (op===:union || op===:intersection)
                        touched[i]=true;continue
                    end
                    keep_in && push!(uses,(lin,f.orev))
                    keep_out && push!(uses,(lout,-f.orev))
                else
                    keep_in && push!(uses,(lin,f.orev))
                    keep_out && push!(uses,(lout,-f.orev))
                end
                wires=[ [ (key=arcs[k].key,fwd=fwd,pc=arcs[k].pc,
                           t0=arcs[k].t0,t1=arcs[k].t1,crv=arcs[k].crv,
                           vt0=arcs[k].vt0,vt1=arcs[k].vt1,neg=arcs[k].neg)
                         for (k,fwd) in cycles[c] ] for c in region ]
                if !isvol[i] && isempty(lin) && isempty(lout)
                    # sheet region touching no volume — a free result face
                    if _brep_label_keep(op,linL|loutL,nobj) ||
                            !isempty(uses)
                        push!(sheet_out,
                              (sg=f.sg,stag=f.stag,wires=wires,
                               uses=Tuple{UInt64,Int}[],oi=i,ruv=uv,sgn=1))
                    else
                        touched[i]=true
                    end
                    continue
                end
                if op===:union && length(uses)==2
                    # a wall between two kept cells is internal to a fuse —
                    # for sheet operands it stays a standalone piece instead
                    if isvol[i]
                        empty!(uses);touched[i]=true;continue
                    else
                        push!(sheet_out,(sg=f.sg,stag=f.stag,wires=wires,
                                         uses=Tuple{UInt64,Int}[],oi=i,
                                         ruv=uv,sgn=1))
                        continue
                    end
                end
                isempty(uses) && (touched[i]=true;continue)
                if !isvol[i] && op===:intersection
                    # a kept sheet region under intersection IS the result
                    # face — upstream binds the section as the sheet's
                    # image; it is not a cell wall of any volume piece
                    push!(sheet_out,(sg=f.sg,stag=f.stag,wires=wires,
                                     uses=Tuple{UInt64,Int}[],oi=i,ruv=uv,
                                     sgn=1))
                    continue
                end
                push!(kept,(sg=f.sg,stag=f.stag,wires=wires,uses=uses,
                            sgn=uses[1][2],oi=i,ruv=uv))
            end
        end
    end
    if op===:union
        sheet_out=_occ_shared_rim_keys(sheet_out,tol)
        # operands whose rims meet a merge partner are *modified* — the
        # union replaces them by a shared face, so leaving `touched`
        # false would republish the operand verbatim beside the merge
        uf=collect(eachindex(sheet_out))
        findr=function(x)
            while uf[x]!=x;uf[x]=uf[uf[x]];x=uf[x] end
            x
        end
        for x in eachindex(sheet_out), y in x+1:length(sheet_out)
            _brep_unifiable(sheet_out[x],sheet_out[y],tol) || continue
            join=any(w->any(a->any(w2->any(b->b.key==a.key,w2),
                                     sheet_out[y].wires),
                            w),sheet_out[x].wires)
            if join
                ra,rb=findr(x),findr(y)
                ra!=rb && (uf[ra]=rb)
            end
        end
        comps=Dict{Int,Vector{Int}}()
        for x in eachindex(sheet_out)
            push!(get!(Vector{Int},comps,findr(x)),x)
        end
        for mem in values(comps)
            length(mem)>1 || continue
            for mi in mem
                oi=get(sheet_out[mi],:oi,0)
                oi>0 && (touched[oi]=true)
            end
        end
        sheet_out=_brep_unify_union(sheet_out,tol,caller)
    end
    op===:intersection &&
        (kept=_brep_unify_union(kept,tol,caller;same_uses=true))
    return kept,sheet_out,touched
end

# Relabel boundary arcs that trace the same 3-D segment in different
# operand faces onto a shared canonical key. `_brep_unify_union` groups
# merge components by `(:e,ctag)` keys, and `_brep_merge_union_faces`
# drops interior arcs only within a key group — both miss distinct
# operand edges that merely coincide in space, so without relabeling an
# adjacent pair of coplanar sheets never fuses. Coincidence is judged on
# evaluated endpoints, not params: partner rims run opposite directions.
function _occ_shared_rim_keys(specs,tol)
    locs=Tuple{Int,Int,Int}[]
    for (si,spec) in enumerate(specs),(wi,w) in enumerate(spec.wires),
        (ai,_) in enumerate(w)
        push!(locs,(si,wi,ai))
    end
    length(locs)<2 && return specs
    idx=Dict{Tuple{Int,Int,Int},Int}(l=>k for (k,l) in enumerate(locs))
    uf=collect(eachindex(locs))
    findr=function(x)
        while uf[x]!=x;uf[x]=uf[uf[x]];x=uf[x] end
        x
    end
    for x in eachindex(locs), y in x+1:length(locs)
        sx,wx,ax=locs[x];sy,wy,ay=locs[y]
        sx==sy && continue
        a1=specs[sx].wires[wx][ax];a2=specs[sy].wires[wy][ay]
        a1.crv.kind===:degenerate && continue
        a2.crv.kind===:degenerate && continue
        _brep_same_curve3d(a1.crv,a2.crv,tol) || continue
        p1a=_brep_vertex_eval(a1.crv,a1.vt0)
        p1b=_brep_vertex_eval(a1.crv,a1.vt1)
        p2a=_brep_vertex_eval(a2.crv,a2.vt0)
        p2b=_brep_vertex_eval(a2.crv,a2.vt1)
        ((_b3dist(p1a,p2a)<tol && _b3dist(p1b,p2b)<tol) ||
         (_b3dist(p1a,p2b)<tol && _b3dist(p1b,p2a)<tol)) || continue
        ra,rb=findr(x),findr(y)
        ra!=rb && (uf[ra]=rb)
    end
    canon=Dict{Int,Any}()
    out=NamedTuple[]
    for (si,spec) in enumerate(specs)
        wires=map(enumerate(spec.wires)) do (wi,w)
            map(enumerate(w)) do (ai,ar)
                r=findr(idx[(si,wi,ai)])
                ck=get!(canon,r,ar.key)
                ck==ar.key ? ar : merge(ar,(key=ck,))
            end
        end
        push!(out,merge(spec,(wires=wires,)))
    end
    return out
end

# One face's 3-D tangent at an edge midpoint — the fan's reference direction.
function _occ_edge_tangent(crv,t,tol)
    crv.kind===:degenerate && return nothing
    if crv.kind===:line
        return crv.d
    end
    dt=0.25*max(tol,abs(crv.t1-crv.t0))
    a=_brep_curve_eval(crv,t-dt);b=_brep_curve_eval(crv,t+dt)
    return _b3(b,a)
end

# Decompose kept face specs into result cells. Two uses land in the same
# cell when every shared edge pairs them through the angular fan: at an edge
# the incident faces sort by in-plane radial direction, consecutive faces'
# adjacent wedges pair the uses on the wedge's side, and a face alone on a
# free edge merges its own sides. Uses never paired into a cell keep their
# volume-membership mask — a cell emits a solid only when its mask passes
# the op's keep rule, so dropped regions emit no solid while their boundary
# faces may still materialize.
function _occ_mixed_cells(kept,op::Symbol,nobj::Int,caller)
    tol=_BREP_TOL
    # result-edge incidence table (one fan entry per spec per edge)
    emap=Dict{Any,Int}();inc=Vector{Tuple{Int,NamedTuple}}[]
    for (fi,spec) in enumerate(kept), wire in spec.wires, ar in wire
        ei=get(emap,ar.key,0)
        if ei==0
            ei=length(inc)+1;emap[ar.key]=ei
            push!(inc,Tuple{Int,NamedTuple}[])
        end
        push!(inc[ei],(fi,ar))
    end
    uflat=Tuple{Int,Int}[]
    for (fi,spec) in enumerate(kept),ui in eachindex(spec.uses)
        push!(uflat,(fi,ui))
    end
    pos=Dict{Tuple{Int,Int},Int}(u=>k for (k,u) in enumerate(uflat))
    uf=collect(eachindex(uflat))
    function findx(x)
        r=x;while uf[r]!=r;uf[r]=uf[uf[r]];r=uf[r];end
        while uf[x]!=r;uf[x],x=r,uf[x];end
        r
    end
    function uni(a,b)
        ra,rb=findx(a),findx(b)
        ra!=rb && (uf[ra]=rb)
        return nothing
    end
    function useid(fi,sgn)
        for (ui,u) in enumerate(kept[fi].uses)
            u[2]==sgn && return pos[(fi,ui)]
        end
        return 0
    end
    for lst in inc
        seen=Dict{Int,NamedTuple}()
        for (fi,ar) in lst
            haskey(seen,fi) || (seen[fi]=ar)
        end
        fs=collect(keys(seen))
        if length(fs)==1
            fi=only(fs)
            a=useid(fi,+1);b=useid(fi,-1)
            a!=0 && b!=0 && uni(a,b)
            continue
        end
        ar0=lst[1][2]
        tm=(ar0.vt0+ar0.vt1)/2
        Pe=_brep_curve_eval(ar0.crv,tm)
        t=_occ_edge_tangent(ar0.crv,tm,tol)
        if t===nothing || _b3norm(t)<tol
            for a in fs,b in fs
                a<b || continue
                for s in (1,-1)
                    x=useid(a,s);y=useid(b,s)
                    x!=0 && y!=0 && uni(x,y)
                end
            end
            continue
        end
        t=_b3unit(t)
        # stable basis of the fan plane
        ref=abs(t[1])<0.9 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
        u1=_b3unit(_b3cross(t,ref));u2=_b3cross(t,u1)
        rows=NamedTuple[]
        for fi in fs
            ar=seen[fi];spec=kept[fi]
            tm2=(ar.t0+ar.t1)/2
            uvm=_brep_pc_eval(ar.pc,tm2)
            n=_brep_surface_normal(spec.sg,uvm[1],uvm[2])
            # in-plane radial toward the region's own interior sample point,
            # projected perpendicular to the edge tangent
            Pin=_brep_surface_eval(spec.sg,spec.ruv[1],spec.ruv[2])
            r=_b3(Pin,Pe);r=_b3(r,_b3mul(t,_b3dot(r,t)))
            rl=_b3norm(r);rl<tol && continue
            r=(r[1]/rl,r[2]/rl,r[3]/rl)
            push!(rows,(fi=fi,θ=atan(_b3dot(r,u2),_b3dot(r,u1)),n=n))
        end
        sort!(rows,by=x->x.θ)
        mlen=length(rows)
        for k in 1:mlen
            ra=rows[k];rb2=rows[mod1(k+1,mlen)]
            bisθ=k==mlen ? (ra.θ+rb2.θ+2π)/2 : (ra.θ+rb2.θ)/2
            b=_b3mul(u1,cos(bisθ)) .+ _b3mul(u2,sin(bisθ))
            sa=_b3dot(b,ra.n)>0 ? -1 : 1
            sb=_b3dot(b,rb2.n)>0 ? -1 : 1
            x=useid(ra.fi,sa);y=useid(rb2.fi,sb)
            x!=0 && y!=0 && uni(x,y)
        end
    end
    classes=Dict{Int,Vector{Int}}()
    for k in eachindex(uflat)
        push!(get!(Vector{Int},classes,findx(k)),k)
    end
    cells=NamedTuple[]
    for members in values(classes)
        mask=UInt64(0)
        for k in members
            fi,ui=uflat[k];mask|=kept[fi].uses[ui][1]
        end
        order=minimum(uflat[k][1] for k in members)
        push!(cells,(mask=mask,members=members,order=order))
    end
    sort!(cells,by=c->c.order)
    return cells,uflat
end

# Whether a produced face's geometry lies on sheet operand `i`'s original
# domain — its uv footprint on the partner face record must coincide and the
# sample must land inside the original wire. Credit for shared/coincident
# pieces resolves this way rather than by owner, matching OCC's image maps
# listing the shared face under every operand covering it.
function _occ_spec_on_sheet(repi,spec,P,tol)
    for of in repi.faces
        of.sg.occ===spec.sg.occ || continue
        if of.sg.occ===:plane
            _brep_planes_coincident(of.sg,spec.sg) || continue
        else
            # curved carriers coincide only for verbatim same records here
            spec.stag==of.stag || continue
        end
        uv=_brep_surface_uv(of.sg,P)
        _brep_uv_in_face(repi.edges,of,uv,10tol) && return true
    end
    return false
end

# Sheet-operand credit mask of an emitted face: every sheet operand whose
# original domain covers any of the face's sampled points.
function _occ_face_credits(reps,isvol,spec,tol)
    credit=UInt64(0)
    pts=NTuple{3,Float64}[]
    hasproperty(spec,:ruv) &&
        push!(pts,_brep_surface_eval(spec.sg,spec.ruv[1],spec.ruv[2]))
    for wire in spec.wires, ar in wire
        push!(pts,_brep_curve_eval(ar.crv,(ar.vt0+ar.vt1)/2))
    end
    for i in eachindex(reps)
        isvol[i] && continue
        for P in pts
            _occ_spec_on_sheet(reps[i],spec,P,tol) || continue
            credit|=_brep_bit(i);break
        end
    end
    return credit
end

# A spec's uv interior sample — merged specs (union/intersection unify)
# carry no `ruv`, so rebuild one over the spec's own wires.
function _occ_spec_ruv(spec,tol,caller)
    r=get(spec,:ruv,nothing)
    r===nothing || return r
    arcs=NamedTuple[];cycles=Vector{Tuple{Int,Bool}}[]
    for wire in spec.wires
        cy=Tuple{Int,Bool}[]
        for ar in wire
            push!(arcs,(pc=ar.pc,t0=ar.t0,t1=ar.t1))
            push!(cy,(length(arcs),ar.fwd))
        end
        push!(cycles,cy)
    end
    return _brep_region_sample(arcs,cycles,collect(eachindex(cycles)),
                               tol,caller)
end

# Materialize the mixed Boolean into the scratch model. Returns
# (pieces, vol_image, faces, face_image):
#   pieces     — ordered (vtag, srcs, pseudo, cell) — vtag is the sm volume
#                tag (0 for a pseudo operand kept in place)
#   vol_image  — piece indices crediting each volume operand
#   faces      — ordered free-face entries (stag, credit) — `stag` is the sm
#                surface tag (existing tag when a verbatim operand face is
#                republished unchanged)
#   face_image — ordered sm surface tags crediting each sheet operand,
#                free faces first then shell-absorbed pieces
function _occ_mixed_materialize!(m::GeoModel,op::Symbol,tags::Vector{Int},
        isvol::Vector{Bool},reps,kept,sheet_out,cells,uflat,touched,
        nops::Int,nobj::Int,caller)
    tol=_BREP_TOL
    scale=1.0
    for p in values(m.points)
        scale=max(scale,abs(p[1]),abs(p[2]),abs(p[3]))
    end
    tol3=1e-9*scale
    volmask=UInt64(0)
    for i in 1:nops
        isvol[i] && (volmask|=_brep_bit(i))
    end
    # emitted cells: volume membership that passes the op's keep rule
    emitted=NamedTuple[]
    cellof=Dict{Tuple{Int,Int},Int}()
    for c in cells
        (c.mask&volmask)!=0 || continue
        _brep_label_keep(op,c.mask,nobj) || continue
        push!(emitted,c)
        ci=length(emitted)
        for k in c.members
            cellof[uflat[k]]=ci
        end
    end
    # shell membership per emitted cell — a spec whose every use lands in
    # one single emitted cell is interior free material, not boundary
    member=Dict{Int,Vector{Tuple{Int,Int}}}()   # cell -> [(fi, sign)]
    interior=falses(length(kept))
    for (fi,spec) in enumerate(kept)
        isempty(spec.uses) && continue
        ec=Int[]
        for (ui,u) in enumerate(spec.uses)
            ci=get(cellof,(fi,ui),0);ci!=0 && push!(ec,ci)
        end
        if length(ec)>=2 && length(ec)==length(spec.uses) &&
                all(==(ec[1]),ec)
            interior[fi]=true;continue
        end
        for (ui,u) in enumerate(spec.uses)
            ci=get(cellof,(fi,ui),0);ci==0 && continue
            push!(get!(Vector{Tuple{Int,Int}},member,ci),(fi,u[2]))
        end
    end
    # per cell: componentize members on shared result edges, then order the
    # components by the containment forest into volume groups (outer first)
    pieces=NamedTuple[]   # (srcs,pseudo,comps,order,cell)
    for (ci,c) in enumerate(emitted)
        mem=get(member,ci,Tuple{Int,Int}[])
        isempty(mem) && continue
        mfis=sort!(unique(fi for (fi,_) in mem))
        pos=Dict{Int,Int}(fi=>k for (k,fi) in enumerate(mfis))
        uf=collect(1:length(mfis))
        function findroot(x)
            while uf[x]!=x;uf[x]=uf[uf[x]];x=uf[x] end
            x
        end
        edge_faces=Dict{Any,Int}()
        for fi in mfis, wire in kept[fi].wires, ar in wire
            j=get(edge_faces,ar.key,0)
            if j==0
                edge_faces[ar.key]=fi
            else
                ri,rj=findroot(pos[fi]),findroot(pos[j])
                ri!=rj && (uf[ri]=rj)
            end
        end
        compsd=Dict{Int,Vector{Int}}()
        for fi in mfis
            push!(get!(Vector{Int},compsd,findroot(pos[fi])),fi)
        end
        comp_faces=sort!(collect(values(compsd));by=minimum)
        ncomp=length(comp_faces)
        inside=falses(ncomp,ncomp);inside_count=zeros(Int,ncomp)
        for a in 1:ncomp,b in 1:ncomp
            a==b && continue
            if _brep_point_in_solid(
                    _brep_spec_solid_view(kept,comp_faces[b]),
                    _brep_spec_rep(kept,comp_faces[a]),tol,caller)
                inside[b,a]=true;inside_count[a]+=1
            end
        end
        for vg in _brep_volume_groups(inside,inside_count,caller)
            comps=[comp_faces[x] for x in vg]
            srcs=c.mask&volmask
            push!(pieces,(srcs=srcs,pseudo=0,comps=comps,order=c.order,
                          cell=ci))
        end
    end
    vol_image=[Int[] for _ in 1:nops]
    for (pi,p) in enumerate(pieces)
        for i in 1:nops
            isvol[i] && (p.srcs&_brep_bit(i))!=0 && push!(vol_image[i],pi)
        end
    end
    for (pi,p) in enumerate(pieces)
        for i in 1:nops
            (p.srcs&_brep_bit(i))!=0 || continue
            if !touched[i] && vol_image[i]==[pi]
                pieces[pi]=merge(p,(pseudo=i,));break
            end
        end
    end
    sort!(pieces,by=p->(trailing_zeros(p.srcs),p.srcs,p.order))
    vol_image=[Int[] for _ in 1:nops]
    for (pi,p) in enumerate(pieces)
        for i in 1:nops
            isvol[i] && (p.srcs&_brep_bit(i))!=0 && push!(vol_image[i],pi)
        end
    end
    # faces bound inside emitted shells + free faces (interior, dropped-cell
    # and standalone specs)
    emit=falses(length(kept))
    for p in pieces
        p.pseudo==0 || continue
        for cf in p.comps, fi in cf
            emit[fi]=true
        end
    end
    # a sheet spec emitted into a shell is a modified operand upstream —
    # its result face is a fresh topological entity (Extent ≥ 1), so it
    # must not take the verbatim/preserved path even when its domain
    # survived the arrangement unchanged
    for (fi,spec) in enumerate(kept)
        emit[fi] || continue
        oi=get(spec,:oi,0)
        oi>0 && !isvol[oi] && (touched[oi]=true)
    end
    pseudo_surf=Dict{Int,Set{Int}}()
    pseudo_ctag=Dict{Int,Vector{Int}}()
    for p in pieces
        p.pseudo==0 && continue
        ss=Set{Int}()
        for sh in m.volumes[tags[p.pseudo]], sc in m.surface_loops[sh]
            push!(ss,abs(sc))
        end
        pseudo_surf[p.pseudo]=ss
        pseudo_ctag[p.pseudo]=reps[p.pseudo].ctags
    end
    reuse=falses(length(kept))
    reuse_ctag=Dict{Any,Int}()
    for (fi,spec) in enumerate(kept)
        emit[fi] || continue
        spec.stag>0 || continue
        i=0
        for u in spec.uses
            for (pi,p) in enumerate(pieces)
                if p.pseudo!=0 && (u[1]&_brep_bit(p.pseudo))!=0
                    i=p.pseudo;break
                end
            end
            i==0 || break
        end
        i==0 && continue
        spec.stag in pseudo_surf[i] || continue
        ok=true
        pset=Set(pseudo_ctag[i])
        for w in spec.wires, ar in w
            (ar.key isa Tuple && length(ar.key)>=2 &&
             ar.key[1]===:e && ar.key[2] in pset) || (ok=false;break)
        end
        ok || continue
        reuse[fi]=true;emit[fi]=false
        for w in spec.wires, ar in w
            reuse_ctag[ar.key]=ar.key[2]
        end
    end
    # free faces = interior specs + specs whose surviving uses all landed
    # in dropped cells; an unmodified operand's spec republishes verbatim
    free=Int[]
    for (fi,spec) in enumerate(kept)
        emit[fi] && continue
        reuse[fi] && continue
        any_emitted=false
        for (ui,u) in enumerate(spec.uses)
            get(cellof,(fi,ui),0)!=0 && (any_emitted=true;break)
        end
        isempty(spec.uses) || interior[fi] || !any_emitted || continue
        push!(free,fi)
    end
    # a verbatim spec re-publishes its operand's own sm surface
    verbatim(fi_or_spec)=begin
        spec=fi_or_spec isa Int ? kept[fi_or_spec] : fi_or_spec
        oi=get(spec,:oi,0)
        oi>0 && !isvol[oi] && !touched[oi] && spec.stag==tags[oi]
    end
    # result edge table + 3-D vertex dedupe over everything emitted
    res_edges=NamedTuple[]
    res_verts=NTuple{3,Float64}[]
    emap=Dict{Any,Int}()
    function vertex_of(p)
        for (i,q) in enumerate(res_verts)
            _b3dist(p,q)<tol3 && return i
        end
        push!(res_verts,p);return length(res_verts)
    end
    function scan_arcs(spec)
        for wire in spec.wires, ar in wire
            ei=get!(emap,ar.key) do
                push!(res_edges,(curve=ar.crv,t0=ar.vt0,t1=ar.vt1))
                length(res_edges)
            end
            e=res_edges[ei]
            lo=min(e.t0,ar.vt0,ar.vt1);hi=max(e.t1,ar.vt0,ar.vt1)
            (lo==e.t0 && hi==e.t1) ||
                (res_edges[ei]=(curve=e.curve,t0=lo,t1=hi))
        end
        return nothing
    end
    for (fi,spec) in enumerate(kept)
        # shell members always emit with result-canonical edges — reusing
        # the operand's stored surface would carry its private rim curves
        # where the section table expects the shared `(:s,se)` curve
        (emit[fi] || (fi in free && !verbatim(fi))) && scan_arcs(spec)
    end
    for spec in sheet_out
        verbatim(spec) || scan_arcs(spec)
    end
    edge_verts=Vector{Tuple{Int,Int}}(undef,length(res_edges))
    for (ei,e) in enumerate(res_edges)
        va=vertex_of(_brep_vertex_eval(e.curve,e.t0))
        vb=vertex_of(_brep_vertex_eval(e.curve,e.t1))
        edge_verts[ei]=(va,vb)
    end
    points=Int[];curves=Int[];loops=Int[];surfaces=Int[];shells=Int[]
    pieces_out=Vector{NamedTuple}(undef,length(pieces))
    faces_out=NamedTuple[]    # (stag,credit)
    face_image=[Int[] for _ in 1:nops]
    shell_img=[Int[] for _ in 1:nops]
    try
        ptag=Vector{Int}(undef,length(res_verts))
        for (i,p) in enumerate(res_verts)
            ptag[i]=_add_occ_point!(m,p)
            push!(points,ptag[i])
        end
        for p in points
            delete!(m.point_size,p)
            delete!(m.point_w,p)
        end
        ctag=Vector{Int}(undef,length(res_edges))
        for (i,e) in enumerate(res_edges)
            va,vb=edge_verts[i]
            crv=e.curve
            local c
            if crv.kind===:degenerate
                c=_add_occ_degenerate!(m,ptag[va])
            elseif crv.kind===:line
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_geometry[c]=(occ=:line,origin=crv.o,dir=crv.d,
                                     t0=e.t0,t1=e.t1)
            elseif crv.kind===:circle
                c=_add_occ_circle!(m,ptag[va],ptag[vb],crv.c,crv.n,
                                   crv.x,crv.y,crv.r,e.t0,e.t1)
            elseif crv.kind===:ellipse
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_types[c]=:ellipse
                m.curve_geometry[c]=(occ=:ellipse,center=crv.c,n=crv.n,
                                     X=crv.x,Y=crv.y,rx=crv.rx,ry=crv.ry,
                                     t0=e.t0,t1=e.t1)
            else
                # :bspline — re-emit the evaluator record restamped to the
                # result span; `crv.ctype` preserves the source curve kind
                crv.kind===:bspline || throw(ArgumentError(
                    "$caller: unsupported result curve kind $(crv.kind)"))
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_types[c]=crv.ctype
                m.curve_geometry[c]=merge(crv.g,(t0=e.t0,t1=e.t1))
            end
            ctag[i]=c;push!(curves,c)
        end
        stag=zeros(Int,length(kept))
        function emit_spec(spec)
            pcdict=Dict{Int,NamedTuple}()
            loop_tags=Int[]
            for wire in spec.wires
                ids=Int[]
                for ar in wire
                    c=haskey(reuse_ctag,ar.key) ? reuse_ctag[ar.key] :
                        ctag[emap[ar.key]]
                    efwd=ar.fwd ⊻ ar.neg
                    push!(ids,efwd ? c : -c)
                    pcst=ar.neg ? _brep_flip_pc(ar.pc,caller) : ar.pc
                    ent=get(pcdict,c,nothing)
                    if efwd
                        if !(ent===nothing || ent.fwd===nothing)
                            throw(ErrorException(
                                "$caller: duplicate forward pcurve " *
                                "on Curve[$c]"))
                        end
                        pcdict[c]=(fwd=pcst,
                                   rev=ent===nothing ? nothing : ent.rev)
                    else
                        ent===nothing || ent.rev===nothing ||
                            throw(ErrorException(
                                "$caller: duplicate reverse pcurve " *
                                "on Curve[$c]"))
                        pcdict[c]=(fwd=ent===nothing ? nothing : ent.fwd,
                                   rev=pcst)
                    end
                end
                lt=add_curve_loop!(m,ids)
                push!(loops,lt);push!(loop_tags,lt)
            end
            sg=merge(spec.sg,(pcurves=pcdict,))
            st=_add_occ_surface!(m,spec.sg.occ,loop_tags,sg)
            push!(surfaces,st)
            face=m.surface_geometry[st]
            for wire in spec.wires, ar in wire
                c=haskey(reuse_ctag,ar.key) ? reuse_ctag[ar.key] :
                    ctag[emap[ar.key]]
                ar.crv.kind===:degenerate &&
                    _occ_bind_degenerate!(m,c,face)
            end
            return st
        end
        for (fi,spec) in enumerate(kept)
            if reuse[fi]
                stag[fi]=spec.stag
                continue
            end
            if emit[fi]
                stag[fi]=emit_spec(spec)
                # shell faces absorbed into solids still credit their sheet
                # operands — appended after the free-face images below
                cr=_occ_face_credits(reps,isvol,spec,tol)
                for i in 1:nops
                    (cr&_brep_bit(i))!=0 && push!(shell_img[i],stag[fi])
                end
            end
        end
        # shells per piece component — the shell entry sign is the use's
        # sign recorded for this piece's cell
        for (pi,p) in enumerate(pieces)
            if p.pseudo!=0
                pieces_out[pi]=(vtag=0,srcs=p.srcs,pseudo=p.pseudo,cell=p.cell)
                continue
            end
            sgnof=Dict{Int,Int}()
            for (fi,sgn) in member[p.cell]
                sgnof[fi]=sgn
            end
            pshells=Int[]
            for cf in p.comps
                stags=Int[]
                for fi in cf
                    push!(stags,sgnof[fi]*stag[fi])
                end
                sh=add_surface_loop!(m,stags)
                push!(shells,sh);push!(pshells,sh)
            end
            v=_alloc_tag!(m,3,0,caller)
            m.volumes[v]=pshells
            pieces_out[pi]=(vtag=v,srcs=p.srcs,pseudo=p.pseudo,cell=p.cell)
        end
        # free faces — interior + dropped-cell specs first (kept order),
        # then standalone specs
        for fi in free
            spec=kept[fi]
            stag[fi]=verbatim(fi) ? spec.stag : emit_spec(spec)
            cr=_occ_face_credits(reps,isvol,spec,tol)
            push!(faces_out,(stag=stag[fi],credit=cr))
            for i in 1:nops
                (cr&_brep_bit(i))!=0 && push!(face_image[i],stag[fi])
            end
        end
        for spec in sheet_out
            st=verbatim(spec) ? spec.stag : emit_spec(spec)
            cr=_occ_face_credits(reps,isvol,spec,tol)
            push!(faces_out,(stag=st,credit=cr))
            for i in 1:nops
                (cr&_brep_bit(i))!=0 && push!(face_image[i],st)
            end
        end
        # faces kept bound on a pseudo (unmodified) operand still credit
        # sheet operands covering them — upstream's Generated images list
        # the shared face the same way
        for (fi,spec) in enumerate(kept)
            emit[fi] && continue
            fi in free && continue
            spec.stag>0 || continue
            cr=_occ_face_credits(reps,isvol,spec,tol)
            for i in 1:nops
                (cr&_brep_bit(i))!=0 && push!(shell_img[i],spec.stag)
            end
        end
        for i in 1:nops
            append!(face_image[i],shell_img[i])
        end
    catch
        for sh in shells
            delete!(m.surface_loops,sh)
        end
        _occ_materialize_rollback!(m,points,curves,loops,surfaces,0)
        rethrow()
    end
    return pieces_out,vol_image,faces_out,face_image,Set{Int}(surfaces)
end

# The mixed-dimensional path: at least one operand is a surface. Volume
# results materialize through the same analytic cell decomposition as the
# volume-only engine; surface results come back as free faces. Operation
# rules mirror the oracle: fuse never mixes solids with sheets, and
# difference rejects a solid object minus a surface tool — OCCT's
# BOPAlgo_AlertBOPNotAllowed for both.
function _model_occ_boolean_mixed(m::GeoModel,op::Symbol,odt,tdt,requested,
                                  removeObject,removeTool,caller)
    kernel=m.occ_kernel
    nobj=length(odt)
    in_dt=vcat(odt,tdt)
    N=length(in_dt)
    anyvol=any(d->d==3,first.(in_dt))
    if op===:union && anyvol
        throw(ErrorException("Union failed - BOPAlgo_AlertBOPNotAllowed"))
    end
    if op===:difference && anyvol &&
            any(d->d==2,first.(tdt))
        throw(ErrorException(
            "Difference failed - BOPAlgo_AlertBOPNotAllowed"))
    end
    sm=GeoModel()
    seen=Set{UInt64}()
    matmap=Dict{Tuple{Int,Int},UInt64}()
    snap=NTuple{3,Float64}[]
    for (d,t) in in_dt
        id=kernel.bindings[(d,t)]
        _model_occ_dist_materialize!(sm,kernel,d,id,seen,caller;snap=snap,
            untyped_defaults=true)
        _model_occ_boolean_matmap!(matmap,kernel,id,caller)
    end
    sm.next_tag[1]=max(sm.next_tag[1],maximum(keys(sm.points);init=0))
    sm.next_tag[2]=max(sm.next_tag[2],maximum(keys(sm.curves);init=0))
    sm.next_tag[3]=max(sm.next_tag[3],maximum(keys(sm.surfaces);init=0))
    sm.next_tag[4]=max(sm.next_tag[4],maximum(keys(sm.volumes);init=0))
    tags=Int[t for (_,t) in in_dt]
    isvol=Bool[d==3 for (d,_) in in_dt]
    reps=NamedTuple[]
    for i in 1:N
        push!(reps,isvol[i] ? _brep_extract(sm,tags[i],caller) :
                              _occ_sheet_rep(sm,tags[i],caller))
    end
    kept,sheet_out,touched=_occ_mixed_collect(op,reps,isvol,nobj,caller)
    cells,uflat=_occ_mixed_cells(kept,op,nobj,caller)
    pieces,vol_image,faces,face_image,fresh_faces=_occ_mixed_materialize!(
        sm,op,tags,isvol,reps,kept,sheet_out,cells,uflat,touched,N,nobj,
        caller)
    # ── disposition (upstream preserve-numbering) ──
    # preserved: untouched operand whose result is itself — stays bound
    # deleted:   image empty — unbound when its group's remove flag is set
    # modified:  Extent>=1 — unbind+rebind at tag for a 1:1 image, else the
    #            operand unbinds and images bind fresh (without remove, the
    #            operand stays bound and images bind fresh beside it)
    preserved=falses(N)
    for i in 1:N
        if isvol[i]
            imgs=vol_image[i]
            preserved[i]=!touched[i] && length(imgs)==1 &&
                pieces[imgs[1]].pseudo==i
        else
            preserved[i]=!touched[i]
        end
    end
    # under difference the tool is consumed, never imaged — upstream's
    # mapModified lists it empty and `_multiBind` unbinds it on
    # removeTool; border-touching credits must not leak object faces
    # into the tool's map
    if op===:difference
        for i in nobj+1:N
            preserved[i]=false
            empty!(vol_image[i])
            empty!(face_image[i])
        end
    end
    rem(i)=i<=nobj ? removeObject : removeTool
    # explicit tag — upstream `_multiBind(result, tag, …, highestDimOnly,
    # recursive)`: only the result's single top-dim entity binds, at `tag`
    if requested>0
        top_dim=isempty(pieces) ? 2 : 3
        members=top_dim==3 ? pieces : faces
        length(members)==1 || throw(ErrorException(
            "Cannot bind multiple " *
            (top_dim==3 ? "volumes" : "surfaces") *
            " to single tag $requested"))
        member=only(members)
        pub=_ModelOCCBooleanPub(caller,kernel,sm,matmap,
            collect(kernel.maxima),kernel.next_shape,
            Dict{UInt64,_ModelOCCShape}(),Tuple{Int,Int,UInt64}[],
            Dict{Tuple{Int,Int},UInt64}(),Dict{Tuple{Int,Int},Int}(),
            Dict{Tuple{Int,Int},Int}())
        if top_dim==3
            _model_occ_boolean_solid!(pub,member.vtag,requested,caller)
        else
            _model_occ_boolean_face!(pub,member.stag,caller,requested)
        end
        for (i,(d,t)) in enumerate(in_dt)
            rem(i) || continue
            id=get(kernel.bindings,(d,t),nothing)
            id===nothing || _model_occ_unbind!(kernel,id,true)
        end
        for (dim,ktag,id) in pub.planned
            _model_occ_bind_shape!(kernel,dim,ktag,id,pub.staged[id])
        end
        isempty(pub.planned) ||
            (kernel.next_shape=max(kernel.next_shape,pub.next_id))
        _model_occ_collect_shapes!(kernel)
        out=Tuple{Int,Int}[(top_dim,requested)]
        out_map=Vector{Tuple{Int,Int}}[]
        for (i,(d,t)) in enumerate(in_dt)
            lst=Tuple{Int,Int}[]
            if preserved[i]
                rem(i) || push!(lst,(d,t))
            elseif isvol[i]
                top_dim==3 && !isempty(vol_image[i]) &&
                    push!(lst,(3,requested))
            else
                top_dim==2 && !isempty(face_image[i]) &&
                    push!(lst,(2,requested))
            end
            push!(out_map,lst)
        end
        return out,out_map
    end
    # 1:1 rebind choice: an operand whose whole image is one fresh entity
    # rebinds at its own tag (upstream mapModified with a single image).
    # Operand order decides contention — the first sole claimant wins.
    claim=Dict{Tuple{Int,Int},Int}()    # (dim, sm tag|piece) -> operand
    for i in 1:N
        preserved[i] && continue
        if isvol[i]
            # the Extent==1 rebind counts only fresh (modified) images —
            # a piece reusing an existing operand is not a modification
            imgs=Int[pi for pi in vol_image[i] if pieces[pi].pseudo==0]
            length(imgs)==1 &&
                (claim[(3,pieces[imgs[1]].vtag)]=get(claim,
                    (3,pieces[imgs[1]].vtag),i))
        else
            imgs=Int[st for st in face_image[i] if st in fresh_faces]
            length(imgs)==1 &&
                (claim[(2,imgs[1])]=get(claim,(2,imgs[1]),i))
        end
    end
    # fresh-tag cursors: surviving maxima like upstream's post-removal
    # `_recomputeMaxTag` — bound entities not among the operands plus the
    # operand tags that stay bound (preserved or reclaimed 1:1)
    surviving=[0,0,0,0]
    in_set=Set(in_dt)
    for (d,t) in keys(kernel.bindings)
        d<0 && continue
        (d,t) in in_set && continue
        surviving[d+1]=max(surviving[d+1],t)
    end
    for (i,(d,t)) in enumerate(in_dt)
        preserved[i] && (surviving[d+1]=max(surviving[d+1],t))
        preserved[i] && continue
        rem(i) || continue
        if isvol[i]
            imgs=Int[pi for pi in vol_image[i] if pieces[pi].pseudo==0]
            if length(imgs)==1 &&
                    get(claim,(3,pieces[imgs[1]].vtag),0)==i
                surviving[4]=max(surviving[4],t)
            end
        else
            imgs=Int[st for st in face_image[i] if st in fresh_faces]
            if length(imgs)==1 && get(claim,(2,imgs[1]),0)==i
                surviving[3]=max(surviving[3],t)
            end
        end
    end
    cursor=collect(kernel.maxima)
    # top-level dims reseed to post-removal maxima; sub-entity (dim 0/1)
    # cursors stay on the pre-op maximum — shared sub-shapes can survive an
    # operand unbind through ownership
    cursor[5]=surviving[3];cursor[6]=surviving[4]
    pub=_ModelOCCBooleanPub(caller,kernel,sm,matmap,
        cursor,kernel.next_shape,
        Dict{UInt64,_ModelOCCShape}(),Tuple{Int,Int,UInt64}[],
        Dict{Tuple{Int,Int},UInt64}(),Dict{Tuple{Int,Int},Int}(),
        Dict{Tuple{Int,Int},Int}())
    # stage result solids in `out` order — fresh tags land in traversal
    # order like upstream's SOLID loop — then free faces. Claimed faces
    # pick up their operand tag inside shell staging through `pub.desired`.
    for ((dim,st),i) in claim
        rem(i) && (pub.desired[(dim,st)]=tags[i])
    end
    ktag3=Dict{Int,Int}()   # sm vtag -> kernel tag
    for p in pieces
        p.pseudo==0 || continue
        i=get(claim,(3,p.vtag),0)
        ktag=(i!=0 && rem(i)) ? tags[i] : _model_occ_boolean_tag!(pub,3)
        ktag3[p.vtag]=ktag
        _model_occ_boolean_solid!(pub,p.vtag,ktag,caller)
    end
    ktag2=Dict{Int,Int}()   # sm stag -> kernel tag for free faces
    for fe in faces
        st=fe.stag
        _model_occ_boolean_face!(pub,st,caller)
        ktag2[st]=pub.ktags[(2,st)]
    end
    # kernel mutation: unbind removed operands recursively, then bind the
    # staged shapes — upstream's remove-then-bind order. With
    # remove=false modified operands stay bound beside their images.
    for (i,(d,t)) in enumerate(in_dt)
        preserved[i] && continue
        rem(i) || continue
        id=get(kernel.bindings,(d,t),nothing)
        id===nothing || _model_occ_unbind!(kernel,id,true)
    end
    for (dim,ktag,id) in pub.planned
        _model_occ_bind_shape!(kernel,dim,ktag,id,pub.staged[id])
    end
    isempty(pub.planned) ||
        (kernel.next_shape=max(kernel.next_shape,pub.next_id))
    _model_occ_collect_shapes!(kernel)
    # `out` — upstream `_multiBind` FACE traversal: each result face lists
    # once per compound visit while `!exists || _toPreserve`. A face bound
    # at an operand tag (preserved or 1:1-claimed) pushes per shell
    # occurrence plus its own top-level binding; a fresh face pushes only
    # at its standalone visit — shell-internal fresh faces never list.
    resolve2=st->begin
        st===nothing && return 0
        # verbatim faces of preserved operands keep their original binding
        get(ktag2,st,get(pub.ktags,(2,st),
                         get(matmap,(2,st),0)!=0 ? st : 0))
    end
    out=Tuple{Int,Int}[]
    for p in pieces
        if p.pseudo!=0
            push!(out,(3,tags[p.pseudo]))
        else
            push!(out,(3,ktag3[p.vtag]))
        end
    end
    operand_tags=Set{Int}(tags)
    topreserve(st)=resolve2(st) in operand_tags
    # shell visits in piece order — preserved faces push per occurrence
    for p in pieces
        p.pseudo==0 || continue
        for sh in sm.volumes[p.vtag], ss in sm.surface_loops[sh]
            st=abs(ss)
            topreserve(st) || continue
            push!(out,(2,resolve2(st)))
        end
    end
    # standalone visits: free faces once each; preserved faces bound only
    # inside shells push their own top-level binding once more
    free_stags=Set{Int}(fe.stag for fe in faces)
    for fe in faces
        push!(out,(2,ktag2[fe.stag]))
    end
    pushed_binding=Set{Int}()
    for p in pieces
        p.pseudo==0 || continue
        for sh in sm.volumes[p.vtag], ss in sm.surface_loops[sh]
            st=abs(ss)
            (topreserve(st) && !(st in free_stags) &&
             !(st in pushed_binding)) || continue
            push!(pushed_binding,st);push!(out,(2,resolve2(st)))
        end
    end
    for (i,(d,t)) in enumerate(in_dt)
        preserved[i] || continue
        dt=(d,t)
        dt in out || push!(out,dt)
    end
    # maps: volumes list volume pieces; surfaces list credited faces — free
    # faces first then shell-absorbed pieces
    out_map=Vector{Tuple{Int,Int}}[]
    for (i,(d,t)) in enumerate(in_dt)
        lst=Tuple{Int,Int}[]
        if preserved[i]
            push!(lst,(d,t))
        elseif isvol[i]
            for pi in vol_image[i]
                p=pieces[pi]
                push!(lst,(3,p.pseudo!=0 ? tags[p.pseudo] :
                           ktag3[p.vtag]))
            end
        else
            for st in face_image[i]
                kt=resolve2(st)
                kt!=0 && push!(lst,(2,kt))
            end
        end
        push!(out_map,lst)
    end
    return out,out_map
end
