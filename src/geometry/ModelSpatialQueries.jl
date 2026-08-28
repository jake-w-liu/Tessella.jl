const _MODEL_EMPTY_BOUNDS =
    (Inf,Inf,Inf,-Inf,-Inf,-Inf)

@inline function _model_bounds_include(
    bounds::NTuple{6,Float64},point::NTuple{3,Float64})
    return (min(bounds[1],point[1]),min(bounds[2],point[2]),
            min(bounds[3],point[3]),max(bounds[4],point[1]),
            max(bounds[5],point[2]),max(bounds[6],point[3]))
end

@inline function _model_bounds_union(
    first::NTuple{6,Float64},second::NTuple{6,Float64})
    return (min(first[1],second[1]),min(first[2],second[2]),
            min(first[3],second[3]),max(first[4],second[4]),
            max(first[5],second[5]),max(first[6],second[6]))
end

function _model_bounds_checked(
    values,caller::AbstractString,entity::AbstractString)
    bounds=NTuple{6,Float64}(values)
    all(isfinite,bounds) || throw(ArgumentError(
        "$caller: $entity has a non-finite bounding box"))
    (bounds[1]<=bounds[4] && bounds[2]<=bounds[5] &&
     bounds[3]<=bounds[6]) || throw(ErrorException(
        "$caller: $entity has an invalid bounding box; rebuild the model"))
    return bounds
end

function _model_bounds_from_points(
    m::GeoModel,points,caller::AbstractString,entity::AbstractString)
    isempty(points) && throw(ErrorException(
        "$caller: $entity has no boundary Points; rebuild the model"))
    bounds=_MODEL_EMPTY_BOUNDS
    for point in points
        coordinate=get(m.points,point,nothing)
        coordinate===nothing && throw(ErrorException(
            "$caller: $entity references missing Point[$point]; rebuild the model"))
        bounds=_model_bounds_include(bounds,coordinate)
    end
    return _model_bounds_checked(bounds,caller,entity)
end

function _model_bounds_from_surface_mesh(
    mesh::Mesh,caller::AbstractString,entity::AbstractString)
    ntris(mesh)>0 || throw(ErrorException(
        "$caller: $entity has no boundary triangles"))
    bounds=_MODEL_EMPTY_BOUNDS
    node_count=size(mesh.coords,2)
    @inbounds for triangle in axes(mesh.tris,2), local_node in axes(mesh.tris,1)
        node=Int(mesh.tris[local_node,triangle])
        1<=node<=node_count || throw(ErrorException(
            "$caller: $entity has an invalid boundary node; rebuild the model"))
        point=(mesh.coords[1,node],mesh.coords[2,node],mesh.coords[3,node])
        bounds=_model_bounds_include(bounds,point)
    end
    return _model_bounds_checked(bounds,caller,entity)
end

function _model_bounds_unit_axis(axis,caller::AbstractString,entity::AbstractString)
    scale=max(abs(axis[1]),abs(axis[2]),abs(axis[3]))
    (isfinite(scale) && scale>0) || throw(ErrorException(
        "$caller: $entity has an invalid axis; rebuild the model"))
    scaled=(axis[1]/scale,axis[2]/scale,axis[3]/scale)
    magnitude=hypot(scaled...)
    (isfinite(magnitude) && magnitude>0) || throw(ErrorException(
        "$caller: $entity has an invalid axis; rebuild the model"))
    return (scaled[1]/magnitude,scaled[2]/magnitude,scaled[3]/magnitude)
end

function _model_cylinder_bounds(spec,caller::AbstractString,entity::AbstractString)
    center=spec.center
    height=spec.height
    radius=spec.radius
    (isfinite(height) && height>0 && isfinite(radius) && radius>0) ||
        throw(ErrorException(
            "$caller: $entity has invalid cylinder parameters; rebuild the model"))
    axis=_model_bounds_unit_axis(spec.axis,caller,entity)
    endpoint=ntuple(index->center[index]+height*axis[index],3)
    radial=ntuple(index->
        radius*sqrt(max(0.0,muladd(-axis[index],axis[index],1.0))),3)
    return _model_bounds_checked((
        ntuple(index->min(center[index],endpoint[index])-radial[index],3)...,
        ntuple(index->max(center[index],endpoint[index])+radial[index],3)...),
        caller,entity)
end

function _model_cone_bounds(spec,caller::AbstractString,entity::AbstractString)
    center=spec.center
    height=spec.height
    r1=spec.r1
    r2=spec.r2
    (isfinite(height) && height>0 && isfinite(r1) && isfinite(r2) &&
     r1>=0 && r2>=0 && (r1>0 || r2>0)) || throw(ErrorException(
        "$caller: $entity has invalid cone parameters; rebuild the model"))
    axis=_model_bounds_unit_axis(spec.axis,caller,entity)
    endpoint=ntuple(index->center[index]+height*axis[index],3)
    orthogonal=ntuple(index->
        sqrt(max(0.0,muladd(-axis[index],axis[index],1.0))),3)
    return _model_bounds_checked((
        ntuple(index->min(center[index]-r1*orthogonal[index],
                          endpoint[index]-r2*orthogonal[index]),3)...,
        ntuple(index->max(center[index]+r1*orthogonal[index],
                          endpoint[index]+r2*orthogonal[index]),3)...),
        caller,entity)
end

function _model_volume_bounds(
    m::GeoModel,tag::Int,caller::AbstractString)
    entity="Volume[$tag]"
    explicit=!isempty(m.volumes[tag])
    encodings=(haskey(m.box_extents,tag),haskey(m.cylinders,tag),
               haskey(m.spheres,tag),haskey(m.cones,tag),
               haskey(m.booleans,tag))
    count(identity,(explicit,encodings...))<=1 || throw(ErrorException(
        "$caller: $entity has multiple native encodings; rebuild the model"))
    haskey(m.booleans,tag)==haskey(m.boolean_operands,tag) ||
        throw(ErrorException(
            "$caller: $entity has inconsistent Boolean snapshot ownership; " *
            "rebuild the model"))

    if encodings[1]
        x,y,z,dx,dy,dz=m.box_extents[tag]
        (isfinite(dx) && isfinite(dy) && isfinite(dz) &&
         dx>0 && dy>0 && dz>0) || throw(ErrorException(
            "$caller: $entity has invalid box parameters; rebuild the model"))
        return _model_bounds_checked(
            (x,y,z,x+dx,y+dy,z+dz),caller,entity)
    elseif encodings[2]
        return _model_cylinder_bounds(m.cylinders[tag],caller,entity)
    elseif encodings[3]
        sphere=m.spheres[tag]
        center=sphere.center
        radius=sphere.radius
        (isfinite(radius) && radius>0) || throw(ErrorException(
            "$caller: $entity has invalid sphere parameters; rebuild the model"))
        return _model_bounds_checked((
            center[1]-radius,center[2]-radius,center[3]-radius,
            center[1]+radius,center[2]+radius,center[3]+radius),caller,entity)
    elseif encodings[4]
        return _model_cone_bounds(m.cones[tag],caller,entity)
    elseif encodings[5]
        spec=m.booleans[tag]
        first_operand,second_operand=m.boolean_operands[tag]
        surface=mesh_boolean(first_operand,second_operand,spec.op)
        return _model_bounds_from_surface_mesh(surface,caller,entity)
    end

    explicit || throw(ArgumentError(
        "$caller: $entity has no native solid encoding"))
    points=_model_points_of(m,[(3,tag)],caller)
    return _model_bounds_from_points(m,points,caller,entity)
end

function _model_entity_bounding_box(
    m::GeoModel,dimension::Int,tag::Int,caller::AbstractString)
    entities=_model_entity_dictionary(m,dimension)
    haskey(entities,tag) || throw(ArgumentError(
        "$caller: unknown entity ($dimension,$tag)"))
    if dimension==0
        point=m.points[tag]
        return _model_bounds_checked(
            (point[1],point[2],point[3],point[1],point[2],point[3]),
            caller,"Point[$tag]")
    elseif dimension==1
        return _model_bounds_from_points(
            m,m.curves[tag],caller,"Curve[$tag]")
    elseif dimension==2
        points=_model_points_of(m,[(2,tag)],caller)
        return _model_bounds_from_points(m,points,caller,"Surface[$tag]")
    end
    return _model_volume_bounds(m,tag,caller)
end

"""
    model_bounding_box(model, dim, tag) -> NTuple{6,Float64}

Return `(xmin, ymin, zmin, xmax, ymax, zmax)` for one existing positive-tag
entity. Passing `dim=-1, tag=-1` returns the union over the whole nonempty model.
Explicit straight-edge topology and native primitive bounds are analytical;
Boolean bounds belong to the operation-time result snapshot.
"""
function model_bounding_box(m::GeoModel,dim,tag)
    caller="model_bounding_box"
    if dim isa Integer && !(dim isa Bool) && tag isa Integer &&
            !(tag isa Bool) && dim==-1 && tag==-1
        entities=model_entities(m)
        isempty(entities) && throw(ArgumentError(
            "$caller: the model has no entities"))
        bounds=_MODEL_EMPTY_BOUNDS
        for (dimension,entity_tag) in entities
            bounds=_model_bounds_union(bounds,_model_entity_bounding_box(
                m,dimension,entity_tag,caller))
        end
        return _model_bounds_checked(bounds,caller,"model")
    end
    dimension=_dimension(dim,caller)
    entity_tag=_tag(tag,caller,dimension)
    entity_tag>0 || throw(ArgumentError(
        "$caller: entity tag must be positive"))
    return _model_entity_bounding_box(m,dimension,entity_tag,caller)
end

function _model_spatial_coordinate(value,caller::AbstractString,name::AbstractString)
    value isa Real || throw(ArgumentError(
        "$caller: $name must be a real number"))
    value isa Bool && throw(ArgumentError(
        "$caller: $name must not be Bool"))
    coordinate=try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: $name must be Float64-representable"))
    end
    isfinite(coordinate) || throw(ArgumentError(
        "$caller: $name must be finite"))
    return coordinate
end

"""
    model_entities_in_bounding_box(model, xmin, ymin, zmin,
                                   xmax, ymax, zmax, dim=-1)

Return detached, sorted entities whose complete analytical bounding box is inside
the supplied finite box. `dim=-1` selects all dimensions; reversed bounds simply
select no entities. Embeddings do not enlarge their target entity.
"""
function model_entities_in_bounding_box(
    m::GeoModel,xmin,ymin,zmin,xmax,ymax,zmax,dim=-1)
    caller="model_entities_in_bounding_box"
    dimension=_query_dimension(dim,caller)
    values=(xmin,ymin,zmin,xmax,ymax,zmax)
    names=("xmin","ymin","zmin","xmax","ymax","zmax")
    box=ntuple(index->_model_spatial_coordinate(
        values[index],caller,names[index]),6)
    (box[1]<=box[4] && box[2]<=box[5] && box[3]<=box[6]) ||
        return Tuple{Int,Int}[]
    entities=model_entities(m,dimension)
    result=Tuple{Int,Int}[]
    sizehint!(result,length(entities))
    for (entity_dimension,tag) in entities
        bounds=_model_entity_bounding_box(
            m,entity_dimension,tag,caller)
        bounds[1]>=box[1] && bounds[2]>=box[2] && bounds[3]>=box[3] &&
        bounds[4]<=box[4] && bounds[5]<=box[5] && bounds[6]<=box[6] &&
            push!(result,(entity_dimension,tag))
    end
    return result
end
