@inline function _model_identity_tag(value,dimension::Int,caller::AbstractString,
                                     what::AbstractString)
    tag=_tag(value,caller,dimension)
    tag>0 || throw(ArgumentError("$caller: $what must be positive"))
    return tag
end

"""
    set_entity_name!(model, dim, tag, name) -> String

Set or replace the name of an existing positive-tag entity. Names need not be
unique. An empty name removes the current name, and a missing entity is a no-op.
Return the resulting name.
"""
function set_entity_name!(m::GeoModel,dim,tag,name)
    caller="set_entity_name!"
    dimension=_dimension(dim,caller)
    entity_tag=_model_identity_tag(tag,dimension,caller,"entity tag")
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    haskey(_model_entity_dictionary(m,dimension),entity_tag) || return ""
    entity_name=String(name)
    key=(dimension,entity_tag)
    if isempty(entity_name)
        delete!(m.entity_names,key)
    else
        m.entity_names[key]=entity_name
    end
    return entity_name
end

"""
    model_entity_name(model, dim, tag) -> String

Return the name of an existing positive-tag entity, or an empty string when the
entity is missing or unnamed.
"""
function model_entity_name(m::GeoModel,dim,tag)
    caller="model_entity_name"
    dimension=_dimension(dim,caller)
    entity_tag=_model_identity_tag(tag,dimension,caller,"entity tag")
    haskey(_model_entity_dictionary(m,dimension),entity_tag) || return ""
    return get(m.entity_names,(dimension,entity_tag),"")
end

"""
    remove_entity_name!(model, name) -> Int

Remove `name` from every model entity carrying it and return the number removed.
A missing or empty name is a no-op.
"""
function remove_entity_name!(m::GeoModel,name)
    caller="remove_entity_name!"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    entity_name=String(name)
    isempty(entity_name) && return 0
    keys_to_remove=Tuple{Int,Int}[
        key for (key,existing) in m.entity_names if existing==entity_name]
    for key in keys_to_remove
        delete!(m.entity_names,key)
    end
    return length(keys_to_remove)
end

function _model_identity_rekey(
    dictionary::Dict{Int,T},old_tag::Int,new_tag::Int,
    caller::AbstractString,what::AbstractString;required::Bool=false) where T
    haskey(dictionary,new_tag) && throw(ArgumentError(
        "$caller: $what already contains target tag $new_tag without a model entity"))
    if !haskey(dictionary,old_tag)
        required && throw(ErrorException(
            "$caller: $what is missing source tag $old_tag; rebuild the model"))
        return copy(dictionary)
    end
    result=copy(dictionary)
    value=pop!(result,old_tag)
    result[new_tag]=value
    return result
end

@inline function _model_identity_signed_tag(value::Int,old_tag::Int,new_tag::Int)
    abs(value)==old_tag || return value
    return value<0 ? -new_tag : new_tag
end

function _model_identity_names(
    m::GeoModel,dimension::Int,old_tag::Int,new_tag::Int,
    caller::AbstractString)
    result=copy(m.entity_names)
    old_key=(dimension,old_tag);new_key=(dimension,new_tag)
    haskey(result,new_key) && throw(ArgumentError(
        "$caller: the target entity name slot $(new_key) is already occupied"))
    if haskey(result,old_key)
        result[new_key]=pop!(result,old_key)
    end
    return result
end

function _model_identity_entity_state(
    dictionary::Dict{Tuple{Int,Int},T},dimension::Int,
    old_tag::Int,new_tag::Int,caller::AbstractString,what::AbstractString) where T
    result=copy(dictionary)
    old_key=(dimension,old_tag);new_key=(dimension,new_tag)
    haskey(result,new_key) && throw(ArgumentError(
        "$caller: the target $what slot $(new_key) is already occupied"))
    if haskey(result,old_key)
        result[new_key]=pop!(result,old_key)
    end
    return result
end

function _model_identity_physical(
    m::GeoModel,dimension::Int,old_tag::Int,new_tag::Int,
    caller::AbstractString)
    result=copy(m.physical)
    for (key,members) in m.physical
        key[1]==dimension || continue
        new_tag in members && throw(ArgumentError(
            "$caller: Physical$(key) already references target entity $new_tag"))
        old_tag in members || continue
        result[key]=Int[member==old_tag ? new_tag : member for member in members]
    end
    return result
end

function _model_identity_embeds(
    m::GeoModel,dimension::Int,old_tag::Int,new_tag::Int,
    caller::AbstractString)
    result=Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}()
    for (key,entities) in m.embeds
        key==(dimension,new_tag) && throw(ArgumentError(
            "$caller: an embedding already targets missing entity $(key)"))
        (dimension,new_tag) in entities && throw(ArgumentError(
            "$caller: an embedding already references missing entity " *
            "($dimension,$new_tag)"))
        new_key=key==(dimension,old_tag) ? (dimension,new_tag) : key
        new_entities=NTuple{2,Int}[
            entity==(dimension,old_tag) ? (dimension,new_tag) : entity
            for entity in entities]
        haskey(result,new_key) && throw(ArgumentError(
            "$caller: retagging would merge embedding targets at $(new_key)"))
        result[new_key]=new_entities
    end
    return result
end

function _model_identity_periodic(
    m::GeoModel,dimension::Int,old_tag::Int,new_tag::Int,
    caller::AbstractString)
    result=Dict{Tuple{Int,Int},ModelPeriodicConstraint}()
    for (key,constraint) in m.periodic
        slave=Int(constraint.slave_entity)
        master=Int(constraint.master_entity)
        key==(constraint.dim,slave) || throw(ErrorException(
            "$caller: periodic key $key disagrees with its slave entity; " *
            "rebuild the model"))
        if constraint.dim==dimension
            (slave==new_tag || master==new_tag) && throw(ArgumentError(
                "$caller: a periodic relation already references missing target " *
                "$(_model_periodic_entity_label(dimension))[$new_tag]"))
            slave=slave==old_tag ? new_tag : slave
            master=master==old_tag ? new_tag : master
        end
        new_key=(constraint.dim,slave)
        haskey(result,new_key) && throw(ArgumentError(
            "$caller: retagging would merge periodic slave relations at $new_key"))
        result[new_key]=ModelPeriodicConstraint(
            constraint.dim,Int32(slave),Int32(master),constraint.affine,
            constraint.reversed,constraint.atol)
    end
    return result
end

function _model_identity_point_state(
    m::GeoModel,old_tag::Int,new_tag::Int,caller::AbstractString)
    points=_model_identity_rekey(
        m.points,old_tag,new_tag,caller,"Point";required=true)
    point_size=_model_identity_rekey(
        m.point_size,old_tag,new_tag,caller,"Point size";required=true)
    curves=copy(m.curves)
    for (curve,(first_point,last_point)) in m.curves
        (first_point==new_tag || last_point==new_tag) && throw(ArgumentError(
            "$caller: Curve[$curve] already references missing target Point[$new_tag]"))
        if first_point==old_tag || last_point==old_tag
            curves[curve]=(first_point==old_tag ? new_tag : first_point,
                           last_point==old_tag ? new_tag : last_point)
        end
    end
    return (points=points,point_size=point_size,curves=curves)
end

function _model_identity_curve_state(
    m::GeoModel,old_tag::Int,new_tag::Int,caller::AbstractString)
    curves=_model_identity_rekey(
        m.curves,old_tag,new_tag,caller,"Curve";required=true)
    loops=copy(m.loops)
    for (loop,signed_curves) in m.loops
        any(value->abs(value)==new_tag,signed_curves) && throw(ArgumentError(
            "$caller: Loop[$loop] already references missing target Curve[$new_tag]"))
        any(value->abs(value)==old_tag,signed_curves) || continue
        loops[loop]=Int[
            _model_identity_signed_tag(value,old_tag,new_tag)
            for value in signed_curves]
    end
    return (curves=curves,loops=loops)
end

function _model_identity_surface_state(
    m::GeoModel,old_tag::Int,new_tag::Int,caller::AbstractString)
    surfaces=_model_identity_rekey(
        m.surfaces,old_tag,new_tag,caller,"Surface";required=true)
    surface_loops=copy(m.surface_loops)
    for (loop,signed_surfaces) in m.surface_loops
        any(value->abs(value)==new_tag,signed_surfaces) && throw(ArgumentError(
            "$caller: Surface Loop[$loop] already references missing target " *
            "Surface[$new_tag]"))
        any(value->abs(value)==old_tag,signed_surfaces) || continue
        surface_loops[loop]=Int[
            _model_identity_signed_tag(value,old_tag,new_tag)
            for value in signed_surfaces]
    end
    return (surfaces=surfaces,surface_loops=surface_loops)
end

function _model_identity_volume_state(
    m::GeoModel,old_tag::Int,new_tag::Int,caller::AbstractString)
    haskey(m.booleans,old_tag)==haskey(m.boolean_operands,old_tag) ||
        throw(ErrorException(
            "$caller: Boolean Volume[$old_tag] has inconsistent snapshot ownership; " *
            "rebuild the model"))
    volumes=_model_identity_rekey(
        m.volumes,old_tag,new_tag,caller,"Volume";required=true)
    box_extents=_model_identity_rekey(
        m.box_extents,old_tag,new_tag,caller,"Box encoding")
    cylinders=_model_identity_rekey(
        m.cylinders,old_tag,new_tag,caller,"Cylinder encoding")
    spheres=_model_identity_rekey(
        m.spheres,old_tag,new_tag,caller,"Sphere encoding")
    cones=_model_identity_rekey(
        m.cones,old_tag,new_tag,caller,"Cone encoding")
    # Boolean operand tags are historical operation-time provenance, not live refs.
    booleans=_model_identity_rekey(
        m.booleans,old_tag,new_tag,caller,"Boolean encoding")
    boolean_operands=_model_identity_rekey(
        m.boolean_operands,old_tag,new_tag,caller,"Boolean operand snapshot")
    return (volumes=volumes,box_extents=box_extents,cylinders=cylinders,
            spheres=spheres,cones=cones,booleans=booleans,
            boolean_operands=boolean_operands)
end

"""
    model_set_tag!(model, dim, tag, new_tag) -> Int

Atomically move an existing positive entity tag to an unused positive tag in the
same dimension. Topology, Physical memberships, embeddings, periodic relations,
entity names, visibility, colors, native solid encodings, and owned Boolean-result
snapshots follow the entity. Automatic tag allocation remains monotonic. Boolean
operand provenance retains the tags recorded when the operation was created.
"""
function model_set_tag!(m::GeoModel,dim,tag,new_tag)
    caller="model_set_tag!"
    dimension=_dimension(dim,caller)
    old_entity_tag=_model_identity_tag(tag,dimension,caller,"source entity tag")
    new_entity_tag=_model_identity_tag(
        new_tag,dimension,caller,"target entity tag")
    entities=_model_entity_dictionary(m,dimension)
    haskey(entities,old_entity_tag) || throw(ArgumentError(
        "$caller: unknown entity ($dimension,$old_entity_tag)"))
    haskey(entities,new_entity_tag) && throw(ArgumentError(
        "$caller: entity ($dimension,$new_entity_tag) already exists"))

    dimension_state=if dimension==0
        _model_identity_point_state(
            m,old_entity_tag,new_entity_tag,caller)
    elseif dimension==1
        _model_identity_curve_state(
            m,old_entity_tag,new_entity_tag,caller)
    elseif dimension==2
        _model_identity_surface_state(
            m,old_entity_tag,new_entity_tag,caller)
    else
        _model_identity_volume_state(
            m,old_entity_tag,new_entity_tag,caller)
    end
    entity_names=_model_identity_names(
        m,dimension,old_entity_tag,new_entity_tag,caller)
    entity_visibility=_model_identity_entity_state(
        m.entity_visibility,dimension,old_entity_tag,new_entity_tag,
        caller,"visibility")
    entity_colors=_model_identity_entity_state(
        m.entity_colors,dimension,old_entity_tag,new_entity_tag,
        caller,"color")
    physical=_model_identity_physical(
        m,dimension,old_entity_tag,new_entity_tag,caller)
    embeds=_model_identity_embeds(
        m,dimension,old_entity_tag,new_entity_tag,caller)
    periodic=_model_identity_periodic(
        m,dimension,old_entity_tag,new_entity_tag,caller)
    next_tag=copy(m.next_tag)
    next_tag[dimension+1]=max(next_tag[dimension+1],new_entity_tag)

    if dimension==0
        m.points=dimension_state.points
        m.point_size=dimension_state.point_size
        m.curves=dimension_state.curves
    elseif dimension==1
        m.curves=dimension_state.curves
        m.loops=dimension_state.loops
    elseif dimension==2
        m.surfaces=dimension_state.surfaces
        m.surface_loops=dimension_state.surface_loops
    else
        m.volumes=dimension_state.volumes
        m.box_extents=dimension_state.box_extents
        m.cylinders=dimension_state.cylinders
        m.spheres=dimension_state.spheres
        m.cones=dimension_state.cones
        m.booleans=dimension_state.booleans
        m.boolean_operands=dimension_state.boolean_operands
    end
    m.entity_names=entity_names
    m.entity_visibility=entity_visibility
    m.entity_colors=entity_colors
    m.physical=physical
    m.embeds=embeds
    m.periodic=periodic
    m.next_tag=next_tag
    return new_entity_tag
end
