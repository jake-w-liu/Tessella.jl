const _MODEL_ENTITY_TYPES=("Point","Line","Plane","Volume")

function _model_metadata_entity(
    m::GeoModel,dim,tag,caller::AbstractString)
    dimension=_dimension(dim,caller)
    entity_tag=_tag(tag,caller,dimension)
    entity_tag>0 || throw(ArgumentError(
        "$caller: entity tag must be positive"))
    haskey(_model_entity_dictionary(m,dimension),entity_tag) ||
        throw(ArgumentError(
            "$caller: unknown entity ($dimension,$entity_tag)"))
    return dimension,entity_tag
end

"""
    model_entity_type(model, dim, tag) -> String

Return the native geometry type of an existing positive-tag entity. Tessella's
explicit entities are `"Point"`, `"Line"`, `"Plane"`, or `"Volume"`; native
primitive and Boolean solids expose only their `"Volume"` entity because their
boundary topology is implicit.
"""
function model_entity_type(m::GeoModel,dim,tag)
    dimension,_=_model_metadata_entity(
        m,dim,tag,"model_entity_type")
    return _MODEL_ENTITY_TYPES[dimension+1]
end

"""
    model_parent(model, dim, tag) -> Tuple{Int,Int}

Return the partition parent of an existing entity. Native `GeoModel` entities are
not partition entities, so the result is always `(-1,-1)`.
"""
function model_parent(m::GeoModel,dim,tag)
    _model_metadata_entity(m,dim,tag,"model_parent")
    return (-1,-1)
end

"""
    model_number_of_partitions(model) -> Int

Return the number of mesh partitions represented by the native geometry model.
`GeoModel` does not own partition entities, so this is always zero.
"""
model_number_of_partitions(::GeoModel)=0

"""
    model_partitions(model, dim, tag) -> Vector{Int}

Return a detached list of partitions containing an existing entity. Native
`GeoModel` entities do not carry partition membership, so the result is empty.
"""
function model_partitions(m::GeoModel,dim,tag)
    _model_metadata_entity(m,dim,tag,"model_partitions")
    return Int[]
end
