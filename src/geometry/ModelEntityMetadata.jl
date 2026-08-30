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

function _model_plane_geometry(
    m::GeoModel,tag::Int,caller::AbstractString)
    loops=m.surfaces[tag]
    isempty(loops) && throw(ErrorException(
        "$caller: Plane[$tag] has no boundary loops; rebuild the model"))
    all_points=_model_points_of(m,[(2,tag)],caller)
    outer_points=unique!(_loop_points(m,first(loops)))
    point_tags=copy(outer_points)
    seen=Set(point_tags)
    for point in all_points
        point in seen && continue
        push!(point_tags,point)
        push!(seen,point)
    end
    coordinates=NTuple{3,Float64}[m.points[point] for point in point_tags]
    anchor,second,third,projection=_model_surface_projection(
        coordinates,point_tags,tag,caller)

    R=Rational{BigInt}
    first_offset=ntuple(index->R(second[index])-R(anchor[index]),3)
    second_offset=ntuple(index->R(third[index])-R(anchor[index]),3)
    normal=(
        first_offset[2]*second_offset[3]-
            first_offset[3]*second_offset[2],
        first_offset[3]*second_offset[1]-
            first_offset[1]*second_offset[3],
        first_offset[1]*second_offset[2]-
            first_offset[2]*second_offset[1],
    )
    squared=sum(component->component^2,normal)
    squared>0 || throw(ErrorException(
        "$caller: Plane[$tag] has a degenerate boundary; rebuild the model"))
    rhs_exact=sum(index->normal[index]*R(anchor[index]),1:3)
    properties=setprecision(BigFloat,256) do
        magnitude=sqrt(BigFloat(squared))
        coefficients=ntuple(
            index->Float64(BigFloat(normal[index])/magnitude),3)
        rhs=Float64(BigFloat(rhs_exact)/magnitude)
        (coefficients...,rhs)
    end
    all(isfinite,properties) || throw(ErrorException(
        "$caller: Plane[$tag] equation is not Float64-representable; " *
        "rebuild the model"))
    return (
        point_tags=point_tags,
        coordinates=coordinates,
        anchor=anchor,
        second=second,
        third=third,
        projection=projection,
        properties=properties,
    )
end

_model_plane_properties(m::GeoModel,tag::Int,caller::AbstractString)=
    _model_plane_geometry(m,tag,caller).properties

"""
    model_entity_properties(model, dim, tag) -> (Vector{Int}, Vector{Float64})

Return detached native geometry properties for an existing entity. For a `Plane`,
the real vector is `[a,b,c,d]` with a unit normal and equation
`a*x+b*y+c*z=d`, oriented by the exterior boundary loop. Visible `Point`, `Line`,
and `Volume` entities have empty property vectors.
"""
function model_entity_properties(m::GeoModel,dim,tag)
    dimension,entity_tag=_model_metadata_entity(
        m,dim,tag,"model_entity_properties")
    dimension==2 || return Int[],Float64[]
    properties=_model_plane_properties(
        m,entity_tag,"model_entity_properties")
    return Int[],collect(properties)
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
