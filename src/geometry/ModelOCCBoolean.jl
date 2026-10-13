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
    pub.staged[id]=_ModelOCCShape(0,(point=pub.sm.points[s],
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
    g=pub.sm.curve_geometry[s]
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

function _model_occ_boolean_face!(pub::_ModelOCCBooleanPub,s::Int,caller)
    skey=(2,s)
    haskey(pub.ids,skey) && return pub.ids[skey]
    kid=_model_occ_boolean_reuse(pub,2,s)
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
    g=pub.sm.surface_geometry[s]
    if haskey(g,:pcurves)
        pcs=Dict{Int,NamedTuple}()
        for (ct,pc) in g.pcurves
            _model_occ_boolean_edge!(pub,ct,caller)
            pcs[pub.ktags[(1,ct)]]=_model_occ_boolean_pc_entry(pc)
        end
        g=merge(g,(pcurves=pcs,))
    end
    id=_model_occ_boolean_id!(pub)
    tag=_model_occ_boolean_tag!(pub,2)
    pub.staged[id]=_ModelOCCShape(2,(kind=get(pub.sm.surface_types,s,:plane),
        geometry=g),children)
    push!(pub.planned,(2,tag,id))
    pub.ids[skey]=id;pub.ktags[skey]=tag
    return id
end

function _model_occ_boolean_shell!(pub::_ModelOCCBooleanPub,sl::Int,caller)
    uses=Tuple{UInt64,Bool}[]
    for sf in pub.sm.surface_loops[sl]
        fid=_model_occ_boolean_face!(pub,abs(sf),caller)
        push!(uses,(fid,sf<0))
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
    for sl in pub.sm.volumes[v]
        push!(children,_model_occ_boolean_shell!(pub,sl,caller))
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
    explicit!=0 && haskey(kernel.bindings,(3,explicit)) && throw(
        ArgumentError("$caller: OpenCASCADE entity with tag $explicit " *
                      "already exists"))
    for (d,t) in Iterators.flatten((odt,tdt))
        d==3 || throw(ArgumentError(
            "$caller: OCC Boolean operands of dimension $d are not yet " *
            "supported (volumes only)"))
        haskey(kernel.bindings,(3,t)) || throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension 3 with " *
            "tag $t"))
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
    for t in in_tags
        id=kernel.bindings[(3,t)]
        _model_occ_dist_materialize!(sm,kernel,3,id,seen,caller)
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
        Dict{Tuple{Int,Int},UInt64}(),Dict{Tuple{Int,Int},Int}())
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
