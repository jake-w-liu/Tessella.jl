function _model_occ_dimension(value,caller)
    value isa Integer && !(value isa Bool) && -2<=value<=3 ||
        throw(ArgumentError("$caller: OCC dimension must be in -2:3"))
    return Int(value)
end

function _model_occ_requested_tag(value,caller)
    value isa Integer && !(value isa Bool) || throw(ArgumentError(
        "$caller: OCC tag must be an integer"))
    typemin(Int32)<=value<=typemax(Int32) || throw(ArgumentError(
        "$caller: OCC tag exceeds Int32"))
    return Int(value)
end

function _model_occ_next_tag(kernel,dimension,requested,caller)
    requested>=0 && return requested
    maximum=kernel.maxima[dimension+3]
    maximum<typemax(Int32) || throw(ArgumentError(
        "$caller: OCC entity tag allocation exceeds Int32"))
    return maximum+1
end

function _model_occ_reserve_shapes(kernel,count,caller)
    UInt64(count)<=typemax(UInt64)-kernel.next_shape || throw(ArgumentError(
        "$caller: OCC shape identity allocation overflow"))
    return UInt64[kernel.next_shape+UInt64(i) for i in 0:count-1]
end

function _model_occ_bind_shape!(kernel,dimension,tag,id,shape)
    _model_occ_entity_map_bind!(kernel.entity_maps[dimension+3],tag)
    kernel.bindings[dimension,tag]=id
    kernel.shape_tags[id]=(dimension,tag)
    kernel.shapes[id]=shape
    kernel.maxima[dimension+3]=max(kernel.maxima[dimension+3],tag)
    kernel.changed=true
    return nothing
end

function _model_occ_add_point!(m::GeoModel,x,y,z;tag=-1,mesh_size=0.0)
    caller="OCC.add_point"
    point=(_finite_scalar(x,caller,"x"),_finite_scalar(y,caller,"y"),
           _finite_scalar(z,caller,"z"))
    size=_finite_scalar(mesh_size,caller,"mesh size")
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,0,_model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(0,allocated)) && throw(ArgumentError(
        "$caller: OCC Point[$allocated] already exists"))
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    shape=_ModelOCCShape(0,(point=point,mesh_size=size),UInt64[])
    _model_occ_bind_shape!(kernel,0,allocated,id,shape)
    kernel.next_shape=id+UInt64(1)
    return allocated
end

function _model_occ_tensor_surface!(m::GeoModel,points,num_points_u,family;
    tag=-1,degree_u=3,degree_v=3,weights=Float64[],knots_u=Float64[],
    knots_v=Float64[],multiplicities_u=Int[],multiplicities_v=Int[],
    wire_tags=Int[],wire_3d=false)
    caller=family===:bezier ? "OCC.add_bezier_surface" : "OCC.add_bspline_surface"
    kernel=m.occ_kernel
    requested=_model_occ_requested_tag(tag,caller)
    allocated=_model_occ_next_tag(kernel,2,requested,caller)
    haskey(kernel.bindings,(2,allocated)) && throw(ArgumentError(
        "$caller: OCC Surface[$allocated] already exists"))
    # Resolve trim wires BEFORE the scratch build so every validator still
    # precedes publication. `wireTags[1]` is the outer boundary, the rest
    # are holes — same ordering convention as imported bounded faces.
    trim_wids=UInt64[]
    for raw in wire_tags
        wt=abs(_model_occ_requested_tag(raw,caller))
        wid=get(kernel.bindings,(-1,wt),nothing)
        wid===nothing && throw(ArgumentError("$caller: unknown OCC Wire[$wt]"))
        wid in trim_wids && throw(ArgumentError(
            "$caller: OCC Wire[$wt] listed twice"))
        isempty(kernel.shapes[wid].geometry.uses) && throw(ArgumentError(
            "$caller: OCC Wire[$wt] has no edges"))
        push!(trim_wids,wid)
    end
    control_tags=Int[_model_occ_requested_tag(t,caller) for t in points]
    scratch=GeoModel()
    for t in control_tags
        id=get(kernel.bindings,(0,t),nothing)
        id===nothing && throw(ArgumentError("$caller: unknown OCC Point[$t]"))
        # Read the shape's copied pole, never the public GModel Point value.
        scratch.points[t]=kernel.shapes[id].geometry.point
    end
    scratch.next_tag[1]=kernel.maxima[3]
    scratch.next_tag[2]=kernel.maxima[4]
    scratch.next_tag[3]=kernel.maxima[5]
    # The existing pure constructor validates knots/weights/control shape and
    # builds exact copied geometry in a small detached model. Zero is a literal
    # OCC tag; its scratch carrier can use one without changing the shape tag.
    build_tag=allocated==0 ? 1 : allocated
    surface=_model_tensor_surface!(scratch,control_tags,num_points_u,family,
        build_tag,degree_u,degree_v,weights,knots_u,knots_v,multiplicities_u,
        multiplicities_v,Int[],wire_3d,caller)
    geometry=scratch.surface_geometry[surface]
    if !isempty(trim_wids)
        # Curved trim: project every wire edge onto the canonical patch
        # domain — the same BRepBuilderAPI_MakeFace path as
        # `addTrimmedSurface`, including seam-period unwrapping for
        # geometrically closed patches (144/114-style bounded faces).
        tag_of=id->_model_occ_shape_tag(kernel,id,caller)
        pc=Dict{Int,NamedTuple}(); vuv=Dict{Int,NTuple{2,Float64}}()
        for wid in trim_wids, (eid,_) in kernel.shapes[wid].geometry.uses
            etag=tag_of(eid)
            haskey(pc,etag) || (pc[etag]=
                (fwd=_trim_project_pcurve(m,kernel,geometry,eid,caller),
                 rev=nothing))
        end
        for wid in trim_wids, (eid,_) in kernel.shapes[wid].geometry.uses,
                vid in kernel.shapes[eid].children
            vtag=tag_of(vid)
            haskey(vuv,vtag) && continue
            pr=_occ_surface_closest(geometry,kernel.shapes[vid].geometry.point,1.0)
            pr===nothing && throw(ArgumentError(
                "$caller: trim wire vertex does not project onto the surface"))
            vuv[vtag]=(pr[1],pr[2])
        end
        geometry=merge(Base.structdiff(geometry,(boundary_curves=nothing,)),
                       (pcurves=pc,vertex_uv=vuv,trimmed=true))
        face=_ModelOCCShape(2,(kind=scratch.surface_types[surface],
                               geometry=geometry),trim_wids)
        id=only(_model_occ_reserve_shapes(kernel,1,caller))
        _model_occ_bind_shape!(kernel,2,allocated,id,face)
        kernel.next_shape=id+UInt64(1)
        return allocated
    end
    curve_tags=geometry.boundary_curves
    boundary_tags=Int[]
    for c in curve_tags,p in scratch.curves[c]
        p in boundary_tags || push!(boundary_tags,p)
    end
    # Preserve the constructor's explicit corner birth order, not coordinate
    # sorting or equality. The two equal-valued corner poles remain distinct.
    length(boundary_tags)==4 || throw(ArgumentError(
        "$caller: tensor Surface must produce four separate boundary vertices"))
    wire_tag=_model_occ_next_tag(kernel,-1,-1,caller)
    ids=_model_occ_reserve_shapes(kernel,10,caller)
    point_ids=Dict(boundary_tags[i]=>ids[i] for i in 1:4)
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    for (i,t) in enumerate(boundary_tags)
        push!(pending,(0,t,ids[i],_ModelOCCShape(0,
            (point=scratch.points[t],mesh_size=0.0),UInt64[])))
    end
    for (i,t) in enumerate(curve_tags)
        a,b=scratch.curves[t]
        push!(pending,(1,t,ids[4+i],_ModelOCCShape(1,
            (kind=scratch.curve_types[t],geometry=scratch.curve_geometry[t]),
            UInt64[point_ids[a],point_ids[b]])))
    end
    signed=scratch.loops[only(scratch.surfaces[surface])]
    uses=Tuple{UInt64,Bool}[(ids[4+findfirst(==(abs(t)),curve_tags)],t<0) for t in signed]
    wire=_ModelOCCShape(-1,(uses=uses,),UInt64[first(use) for use in uses])
    face=_ModelOCCShape(2,(kind=scratch.surface_types[surface],geometry=geometry),UInt64[ids[9]])
    push!(pending,(-1,wire_tag,ids[9],wire),(2,allocated,ids[10],face))
    # Every validator and birth allocation precedes this publication loop.
    for (dimension,t,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,t,id,shape)
    end
    kernel.next_shape=ids[end]+UInt64(1)
    return allocated
end

_model_occ_add_bezier_surface!(m,points,nu;kwargs...)=
    _model_occ_tensor_surface!(m,points,nu,:bezier;kwargs...)
_model_occ_add_bspline_surface!(m,points,nu;kwargs...)=
    _model_occ_tensor_surface!(m,points,nu,:bspline;kwargs...)

function _model_occ_entities(m::GeoModel,dimension=-1)
    d=_query_dimension(dimension,"OCC.get_entities")
    result=Tuple{Int,Int}[]
    for dim in 0:3
        (d<0 || dim==d) || continue
        append!(result,_model_occ_entity_map_entities(m.occ_kernel.entity_maps[dim+3],dim))
    end
    return result
end

function _model_occ_max_tag(m::GeoModel,dimension)
    d=_model_occ_dimension(dimension,"OCC.get_max_tag")
    return m.occ_kernel.maxima[d+3]
end

function _model_occ_set_max_tag!(m::GeoModel,dimension,tag)
    d=_model_occ_dimension(dimension,"OCC.set_max_tag")
    t=_model_occ_requested_tag(tag,"OCC.set_max_tag")
    m.occ_kernel.maxima[d+3]=max(m.occ_kernel.maxima[d+3],t)
    return nothing
end

function _model_occ_shape_tag(kernel,id,caller)
    entry=get(kernel.shape_tags,id,nothing)
    entry===nothing && throw(ArgumentError("$caller: referenced OCC shape is unbound"))
    return last(entry)
end

function _model_occ_synchronize!(m::GeoModel)
    caller="OCC.synchronize"
    kernel=m.occ_kernel
    removed=sort!([(d,t) for (d,t) in kernel.removed if d>=0];by=e->(-e[1],e[2]))
    isempty(removed) || remove_entities!(m,removed,false)
    empty!(kernel.removed)
    imported=false
    # OCC synchronizes the bound shape graph, independently of the current
    # public entity values. Existing public tags are reused without overwrite.
    for ((dimension,tag),id) in sort!(collect(kernel.bindings);by=first)
        dimension==0 || continue
        g=kernel.shapes[id].geometry
        if !haskey(m.points,tag) && !haskey(m.discrete,(0,tag))
            m.points[tag]=g.point;m.point_size[tag]=0.0
            m.next_tag[1]=max(m.next_tag[1],tag)
            imported=true
        elseif haskey(m.points,tag)
            # A transform may replace a bound vertex's geometry in place —
            # the kernel record is authoritative for bound entities.
            m.points[tag]=g.point
        end
        g.mesh_size>0 && (m.point_size[tag]=g.mesh_size)
    end
    for ((dimension,tag),id) in sort!(collect(kernel.bindings);by=first)
        dimension==1 || continue
        haskey(m.curves,tag) || haskey(m.discrete,(1,tag)) || begin
            shape=kernel.shapes[id]
            a,b=shape.children
            m.curves[tag]=(_model_occ_shape_tag(kernel,a,caller),
                           _model_occ_shape_tag(kernel,b,caller))
            m.curve_types[tag]=shape.geometry.kind
            m.curve_geometry[tag]=shape.geometry.geometry
            m.next_tag[2]=max(m.next_tag[2],tag)
            imported=true
        end
    end
    for ((dimension,tag),id) in sort!(collect(kernel.bindings);by=first)
        dimension==2 || continue
        haskey(m.surfaces,tag) || haskey(m.discrete,(2,tag)) || begin
            shape=kernel.shapes[id]
            loops=Int[]
            for (index,wire_id) in enumerate(shape.children)
                wire=kernel.shapes[wire_id]
                signed=Int[(reverse ? -1 : 1)*_model_occ_shape_tag(kernel,edge,caller)
                           for (edge,reverse) in wire.geometry.uses]
                # BRepBuilderAPI_MakeFace stores hole wires REVERSED: the loop
                # traverses them in reversed order with each edge flipped.
                index>1 && (signed=[-s for s in Iterators.reverse(signed)])
                ordered=try
                    _sort_curve_loop(m,signed,caller)
                catch error
                    error isa ArgumentError || rethrow()
                    # OCCFace::_setup attempts to recompute the live GEdgeLoop,
                    # then reports this failure after lower entities were imported.
                    throw(ArgumentError("$caller: Could not fix wire in surface $tag"))
                end
                loop=_alloc_curve_loop_tag(m,0,caller)
                # This loop is the typed GModel carrier of an explicit OCC wire
                # shape. Its internal slot is not the OCC wire namespace tag.
                m.loops[loop]=ordered
                push!(loops,loop)
            end
            m.surfaces[tag]=loops
            m.surface_types[tag]=shape.geometry.kind
            m.surface_geometry[tag]=shape.geometry.geometry
            m.next_tag[3]=max(m.next_tag[3],tag)
            imported=true
        end
    end
    # OCC shells are dim -2 in the shape namespace; each becomes a GModel
    # surface loop (signed member faces) consumed by the solid arm. The
    # kernel-side map keeps the binding stable across repeated syncs. A
    # shell unbound since the last import leaves with its GModel loop —
    # unless a surviving kernel shell still claims that loop tag: solid
    # rebinds (fillet, chamfer, defeature) orphan the old shell id while
    # the replacement claims the same tag until the shape GC sweeps it.
    shell_loops=kernel.shell_loops
    for (id,loop) in collect(shell_loops)
        haskey(kernel.shapes,id) && continue
        delete!(shell_loops,id)
        loop in values(shell_loops) || delete!(m.surface_loops,loop)
    end
    for ((dimension,tag),id) in sort!(collect(kernel.bindings);by=first)
        dimension==-2 || continue
        haskey(shell_loops,id) && continue
        shape=kernel.shapes[id]
        signed=Int[(reverse ? -1 : 1)*_model_occ_shape_tag(kernel,fid,caller)
                   for (fid,reverse) in shape.geometry.uses]
        loop=_alloc_surface_loop_tag(m,0,caller)
        m.surface_loops[loop]=signed
        shell_loops[id]=loop
        imported=true
    end
    for ((dimension,tag),id) in sort!(collect(kernel.bindings);by=first)
        dimension==3 || continue
        haskey(m.volumes,tag) || haskey(m.discrete,(3,tag)) || begin
            shape=kernel.shapes[id]
            m.volumes[tag]=Int[shell_loops[sid] for sid in shape.children]
            # A primitive solid carries its compact analytic encoding, like
            # the geo-pipeline `add_*!` producers: the boundary stays
            # materialized for queries while volume meshing flows through the
            # native analytic tessellation (`_implicit_volume_surface`).
            primitive=get(shape.geometry,:primitive,nothing)
            if primitive!==nothing
                primitive.occ===:cylinder &&
                    (m.cylinders[tag]=primitive.record)
                primitive.occ===:sphere &&
                    (m.spheres[tag]=primitive.record)
                primitive.occ===:cone &&
                    (m.cones[tag]=primitive.record)
            end
            # GEO-bound encodings rebound through a transform restore the
            # same stores `add_*!`/native Booleans maintain.
            geo_enc=get(shape.geometry,:geo_encodings,nothing)
            if geo_enc!==nothing
                haskey(geo_enc,:cylinder) && (m.cylinders[tag]=geo_enc.cylinder)
                haskey(geo_enc,:sphere) && (m.spheres[tag]=geo_enc.sphere)
                haskey(geo_enc,:cone) && (m.cones[tag]=geo_enc.cone)
                haskey(geo_enc,:booleans) && (m.booleans[tag]=geo_enc.booleans)
                haskey(geo_enc,:boolean_operands) &&
                    (m.boolean_operands[tag]=geo_enc.boolean_operands)
                haskey(geo_enc,:boolean_components) &&
                    (m.boolean_components[tag]=geo_enc.boolean_components)
                if haskey(geo_enc,:box)
                    owned=sort!(collect(_model_volume_owned_points(m,tag)))
                    new=_volume_box_extents(m,owned,geo_enc.box)
                    new===nothing || (m.box_extents[tag]=new)
                end
            end
            m.next_tag[4]=max(m.next_tag[4],tag)
            imported=true
        end
    end
    # `box_extents` tracks materialized shells, not kernel bindings — a
    # rebound or sub-entity move that leaves the corners off the per-axis
    # extrema drops the encoding (`transform_entities!` contract).
    for (tag,ext) in collect(m.box_extents)
        (haskey(m.volumes,tag) && !isempty(m.volumes[tag])) ||
            (delete!(m.box_extents,tag);continue)
        owned=sort!(collect(_model_volume_owned_points(m,tag)))
        new=_volume_box_extents(m,owned,ext)
        new===nothing ? delete!(m.box_extents,tag) : (m.box_extents[tag]=new)
    end
    kernel.changed=false
    return imported || !isempty(removed)
end

function _model_occ_contains_shape(kernel,parent,id)
    for child in kernel.shapes[parent].children
        child==id && return true
        _model_occ_contains_shape(kernel,child,id) && return true
    end
    return false
end

function _model_occ_shape_is_owned(kernel,id)
    dimension=kernel.shapes[id].dimension
    owner_dimension=dimension==0 ? 1 : dimension in (-1,1) ? 2 :
                    dimension==2 ? -2 : 3
    for ((d,_),parent) in kernel.bindings
        d==owner_dimension || continue
        _model_occ_contains_shape(kernel,parent,id) && return true
    end
    return false
end

function _model_occ_unbind!(kernel,id,recursive)
    entry=get(kernel.shape_tags,id,nothing)
    entry===nothing && return nothing
    _model_occ_shape_is_owned(kernel,id) && return nothing
    dimension,tag=entry
    _model_occ_entity_map_unbind!(kernel.entity_maps[first(entry)+3],last(entry))
    delete!(kernel.bindings,entry);delete!(kernel.shape_tags,id)
    push!(kernel.removed,entry)
    kernel.maxima[dimension+3]=maximum((t for (d,t) in keys(kernel.bindings)
                                      if d==dimension);init=0)
    kernel.changed=true
    if recursive
        for child in kernel.shapes[id].children
            _model_occ_unbind!(kernel,child,true)
        end
    end
    return nothing
end

function _model_occ_remove!(m::GeoModel,dim_tags,recursive=false)
    recursive isa Bool || throw(ArgumentError("OCC.remove: recursive must be Bool"))
    entities=_model_removal_dim_tags(dim_tags,"OCC.remove")
    kernel=m.occ_kernel
    for (dimension,tag) in entities
        haskey(kernel.bindings,(dimension,tag)) || throw(ArgumentError(
            "OCC.remove: unknown OCC entity ($dimension,$tag)"))
    end
    for entry in entities
        id=get(kernel.bindings,entry,nothing)
        id===nothing || _model_occ_unbind!(kernel,id,recursive)
    end
    _model_occ_collect_shapes!(kernel)
    return nothing
end

function _model_occ_collect_shapes!(kernel)
    reachable=Set{UInt64}()
    stack=collect(values(kernel.bindings))
    while !isempty(stack)
        id=pop!(stack)
        id in reachable && continue
        push!(reachable,id)
        append!(stack,kernel.shapes[id].children)
    end
    filter!(pair->first(pair) in reachable,kernel.shapes)
    return nothing
end
function _model_occ_geometry_equal(a,b)
    a===b && return true
    typeof(a)===typeof(b) || return false
    if a isa NamedTuple || a isa _ModelNURBSSurface || a isa NURBS.NURBSSurface ||
       a isa _ModelOCCNURBSCurve || a isa _ModelOCCReparamCurve ||
       a isa _ModelOCCConicReparam || a isa NURBS.NURBSCurve
        return all(name->_model_occ_geometry_equal(getfield(a,name),getfield(b,name)),
                   fieldnames(typeof(a)))
    elseif a isa Tuple
        length(a)==length(b) || return false
        return all(_model_occ_geometry_equal(x,y) for (x,y) in zip(a,b))
    elseif a isa AbstractArray
        size(a)==size(b) || return false
        return all(_model_occ_geometry_equal(x,y) for (x,y) in zip(a,b))
    elseif a isa AbstractDict
        keys(a)==keys(b) || return false
        return all(_model_occ_geometry_equal(value,b[key]) for (key,value) in a)
    end
    return a==b
end
