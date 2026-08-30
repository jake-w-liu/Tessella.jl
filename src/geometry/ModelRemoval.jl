function _model_removal_dim_tags(dim_tags,caller::AbstractString)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    result=Tuple{Int,Int}[]
    sizehint!(result,length(dim_tags))
    for entry in dim_tags
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif entry isa Tuple && length(entry)==2
            entry
        else
            throw(ArgumentError(
                "$caller: each dim_tags entry must be a (dimension, tag) pair"))
        end
        dimension=_dimension(pair[1],caller)
        entity_tag=_tag(pair[2],caller,dimension)
        entity_tag>0 || throw(ArgumentError(
            "$caller: entity tags must be positive"))
        push!(result,(dimension,entity_tag))
    end
    return result
end

@inline function _model_removal_exists(
    m::GeoModel,removed::Set{Tuple{Int,Int}},dimension::Int,tag::Int)
    return haskey(_model_entity_dictionary(m,dimension),tag) &&
           !((dimension,tag) in removed)
end

function _model_removal_surface_curves(
    m::GeoModel,surface::Int,caller::AbstractString)
    loops=get(m.surfaces,surface,nothing)
    loops===nothing && throw(ErrorException(
        "$caller: Surface[$surface] disappeared during removal planning"))
    curves=Int[]
    for loop in loops
        signed_curves=get(m.loops,loop,nothing)
        signed_curves===nothing && throw(ErrorException(
            "$caller: Surface[$surface] references missing Loop[$loop]; " *
            "rebuild the model"))
        for signed_curve in signed_curves
            curve=abs(signed_curve)
            haskey(m.curves,curve) || throw(ErrorException(
                "$caller: Loop[$loop] references missing Curve[$curve]; " *
                "rebuild the model"))
            push!(curves,curve)
        end
    end
    return curves
end

function _model_removal_volume_surfaces(
    m::GeoModel,volume::Int,caller::AbstractString)
    shells=get(m.volumes,volume,nothing)
    shells===nothing && throw(ErrorException(
        "$caller: Volume[$volume] disappeared during removal planning"))
    surfaces=Int[]
    for shell in shells
        signed_surfaces=get(m.surface_loops,shell,nothing)
        signed_surfaces===nothing && throw(ErrorException(
            "$caller: Volume[$volume] references missing Surface Loop[$shell]; " *
            "rebuild the model"))
        for signed_surface in signed_surfaces
            surface=abs(signed_surface)
            haskey(m.surfaces,surface) || throw(ErrorException(
                "$caller: Surface Loop[$shell] references missing " *
                "Surface[$surface]; rebuild the model"))
            push!(surfaces,surface)
        end
    end
    return surfaces
end

function _model_removal_has_boundary_owner(
    m::GeoModel,removed::Set{Tuple{Int,Int}},dimension::Int,tag::Int,
    caller::AbstractString)
    if dimension==0
        return any(tag in endpoints
                   for (curve,endpoints) in m.curves
                   if _model_removal_exists(m,removed,1,curve))
    elseif dimension==1
        for surface in keys(m.surfaces)
            _model_removal_exists(m,removed,2,surface) || continue
            tag in _model_removal_surface_curves(m,surface,caller) && return true
        end
    elseif dimension==2
        for volume in keys(m.volumes)
            _model_removal_exists(m,removed,3,volume) || continue
            tag in _model_removal_volume_surfaces(m,volume,caller) && return true
        end
    end
    return false
end

function _model_removal_has_embedding_owner(
    m::GeoModel,removed::Set{Tuple{Int,Int}},dimension::Int,tag::Int,
    caller::AbstractString)
    source=(dimension,tag)
    for (target,sources) in m.embeds
        target in removed && continue
        _model_removal_exists(m,removed,target...) || throw(ErrorException(
            "$caller: embedding target $target does not exist; rebuild the model"))
        source in sources && return true
    end
    return false
end

function _model_removal_boundary(
    m::GeoModel,dimension::Int,tag::Int,caller::AbstractString)
    if dimension==0
        return Tuple{Int,Int}[]
    elseif dimension==1
        endpoints=get(m.curves,tag,nothing)
        endpoints===nothing && throw(ErrorException(
            "$caller: Curve[$tag] disappeared during removal planning"))
        for point in endpoints
            haskey(m.points,point) || throw(ErrorException(
                "$caller: Curve[$tag] references missing Point[$point]; " *
                "rebuild the model"))
        end
        return Tuple{Int,Int}[(0,endpoints[1]),(0,endpoints[2])]
    elseif dimension==2
        return Tuple{Int,Int}[
            (1,curve) for curve in
            _model_removal_surface_curves(m,tag,caller)]
    end
    return Tuple{Int,Int}[
        (2,surface) for surface in
        _model_removal_volume_surfaces(m,tag,caller)]
end

function _model_removal_plan(
    m::GeoModel,requested::Vector{Tuple{Int,Int}},recursive::Bool,
    caller::AbstractString)
    removed=Set{Tuple{Int,Int}}()
    function attempt!(dimension::Int,tag::Int,recurse::Bool)
        _model_removal_exists(m,removed,dimension,tag) || return nothing
        _model_removal_has_boundary_owner(
            m,removed,dimension,tag,caller) && return nothing
        _model_removal_has_embedding_owner(
            m,removed,dimension,tag,caller) && return nothing
        boundary=recurse ?
            _model_removal_boundary(m,dimension,tag,caller) :
            Tuple{Int,Int}[]
        push!(removed,(dimension,tag))
        for (boundary_dimension,boundary_tag) in boundary
            attempt!(boundary_dimension,boundary_tag,true)
        end
        return nothing
    end
    for (dimension,tag) in requested
        attempt!(dimension,tag,recursive)
    end
    return removed
end

function _model_removal_state(
    m::GeoModel,removed::Set{Tuple{Int,Int}})
    removed_tags=ntuple(dimension->Set(
        tag for (entity_dimension,tag) in removed
        if entity_dimension==dimension-1),4)

    points=copy(m.points)
    point_size=copy(m.point_size)
    curves=copy(m.curves)
    surfaces=copy(m.surfaces)
    volumes=copy(m.volumes)
    for tag in removed_tags[1]
        delete!(points,tag);delete!(point_size,tag)
    end
    for tag in removed_tags[2]
        delete!(curves,tag)
    end
    for tag in removed_tags[3]
        delete!(surfaces,tag)
    end
    for tag in removed_tags[4]
        delete!(volumes,tag)
    end

    loops=copy(m.loops)
    for (loop,signed_curves) in m.loops
        any(curve->abs(curve) in removed_tags[2],signed_curves) &&
            delete!(loops,loop)
    end
    surface_loops=copy(m.surface_loops)
    for (loop,signed_surfaces) in m.surface_loops
        any(surface->abs(surface) in removed_tags[3],signed_surfaces) &&
            delete!(surface_loops,loop)
    end
    for (surface,surface_loops_used) in surfaces
        all(loop->haskey(loops,loop),surface_loops_used) || throw(ErrorException(
            "remove_entities!: surviving Surface[$surface] lost a Curve Loop; " *
            "rebuild the model"))
    end
    for (volume,shells) in volumes
        all(shell->haskey(surface_loops,shell),shells) || throw(ErrorException(
            "remove_entities!: surviving Volume[$volume] lost a Surface Loop; " *
            "rebuild the model"))
    end

    entity_names=copy(m.entity_names)
    entity_visibility=copy(m.entity_visibility)
    entity_colors=copy(m.entity_colors)
    for entity in removed
        delete!(entity_names,entity)
        delete!(entity_visibility,entity)
        delete!(entity_colors,entity)
    end

    physical=Dict{Tuple{Int,Int},Vector{Int}}()
    physical_names=copy(m.physical_names)
    for (key,members) in m.physical
        retained=Int[member for member in members
                     if !((key[1],member) in removed)]
        if isempty(retained)
            delete!(physical_names,key)
        else
            physical[key]=retained
        end
    end

    embeds=Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}()
    for (target,sources) in m.embeds
        target in removed && continue
        retained=NTuple{2,Int}[source for source in sources
                               if !(source in removed)]
        isempty(retained) || (embeds[target]=retained)
    end

    periodic=Dict{Tuple{Int,Int},ModelPeriodicConstraint}()
    for (key,constraint) in m.periodic
        slave=(constraint.dim,Int(constraint.slave_entity))
        master=(constraint.dim,Int(constraint.master_entity))
        (slave in removed || master in removed) && continue
        periodic[key]=constraint
    end

    box_extents=copy(m.box_extents)
    cylinders=copy(m.cylinders)
    spheres=copy(m.spheres)
    cones=copy(m.cones)
    booleans=copy(m.booleans)
    boolean_operands=copy(m.boolean_operands)
    for tag in removed_tags[4]
        for encoding in (box_extents,cylinders,spheres,cones,
                         booleans,boolean_operands)
            delete!(encoding,tag)
        end
    end

    return (;points,point_size,curves,loops,surfaces,surface_loops,volumes,
            entity_names,entity_visibility,entity_colors,physical,physical_names,
            box_extents,cylinders,spheres,cones,booleans,boolean_operands,
            periodic,embeds)
end

"""
    remove_entities!(model, dim_tags, recursive=false) -> Int

Atomically process ordered positive `(dimension, entity_tag)` removals and return
the number removed. Missing entities and entities still used as an explicit boundary
or embedding source are skipped. With `recursive=true`, each removed entity's
explicit boundary is processed down to Points; embedded entities are not boundaries.

Entity names, visibility, colors, Physical memberships, embedding targets, periodic
relations, native solid encodings, and Boolean-result snapshots owned by removed
entities are cleaned up. Empty Physical groups and their names are removed,
construction loops that would dangle are discarded, and automatic tag counters
remain monotonic. Primitive and Boolean Volumes have no visible boundary entities,
so recursive removal stops at the Volume.
"""
function remove_entities!(m::GeoModel,dim_tags,recursive=false)
    caller="remove_entities!"
    recursive isa Bool || throw(ArgumentError(
        "$caller: recursive must be Bool"))
    requested=_model_removal_dim_tags(dim_tags,caller)
    removed=_model_removal_plan(m,requested,recursive,caller)
    isempty(removed) && return 0
    state=_model_removal_state(m,removed)

    m.points=state.points
    m.point_size=state.point_size
    m.curves=state.curves
    m.loops=state.loops
    m.surfaces=state.surfaces
    m.surface_loops=state.surface_loops
    m.volumes=state.volumes
    m.entity_names=state.entity_names
    m.entity_visibility=state.entity_visibility
    m.entity_colors=state.entity_colors
    m.physical=state.physical
    m.physical_names=state.physical_names
    m.box_extents=state.box_extents
    m.cylinders=state.cylinders
    m.spheres=state.spheres
    m.cones=state.cones
    m.booleans=state.booleans
    m.boolean_operands=state.boolean_operands
    m.periodic=state.periodic
    m.embeds=state.embeds
    return length(removed)
end

function _remove_volume_entity!(m::GeoModel,tag::Int)
    return remove_entities!(m,[(3,tag)],false)>0
end
