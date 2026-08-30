function _model_evaluation_values(
    values,stride::Int,caller::AbstractString,what::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: $what must be a vector or tuple"))
    length(values)%stride==0 || throw(ArgumentError(
        "$caller: $what length must be a multiple of $stride"))
    output=Vector{Float64}(undef,length(values))
    for (index,value) in enumerate(values)
        output[index]=_finite_scalar(value,caller,what)
    end
    return output
end

function _model_evaluation_entity(
    m::GeoModel,dim,tag,allowed,caller::AbstractString)
    dimension,entity_tag=_model_metadata_entity(m,dim,tag,caller)
    dimension in allowed || throw(ArgumentError(
        "$caller: dimension $dimension is unsupported; expected " *
        join(allowed," or ")))
    return dimension,entity_tag
end

@inline function _model_rational_float(
    value::Rational{BigInt},caller::AbstractString,what::AbstractString)
    converted=Float64(value)
    isfinite(converted) || throw(ArgumentError(
        "$caller: $what is not Float64-representable"))
    return converted
end

function _model_line_geometry(
    m::GeoModel,tag::Int,caller::AbstractString)
    first_point,last_point=m.curves[tag]
    haskey(m.points,first_point) || throw(ArgumentError(
        "$caller: Line[$tag] references unknown Point[$first_point]"))
    haskey(m.points,last_point) || throw(ArgumentError(
        "$caller: Line[$tag] references unknown Point[$last_point]"))
    first_coordinate=m.points[first_point]
    last_coordinate=m.points[last_point]
    R=Rational{BigInt}
    delta=ntuple(
        axis->R(last_coordinate[axis])-R(first_coordinate[axis]),3)
    squared=sum(component->component^2,delta)
    squared>0 || throw(ArgumentError(
        "$caller: Line[$tag] has coincident endpoints"))
    return (
        first=first_coordinate,
        last=last_coordinate,
        delta=delta,
        squared=squared,
    )
end

function _model_line_derivative(line,tag::Int,caller::AbstractString)
    return ntuple(axis->_model_rational_float(
        line.delta[axis],caller,"Line[$tag] derivative"),3)
end

function _model_line_point(
    line,parameter,caller::AbstractString,point_index::Int)
    R=Rational{BigInt}
    t=R(parameter)
    return ntuple(3) do axis
        _model_rational_float(
            R(line.first[axis])+t*line.delta[axis],caller,
            "Line value for point $point_index")
    end
end

function _model_line_parameter_exact(line,coordinate::NTuple{3,Float64})
    R=Rational{BigInt}
    offset=ntuple(axis->R(coordinate[axis])-R(line.first[axis]),3)
    numerator=sum(axis->offset[axis]*line.delta[axis],1:3)
    return numerator/line.squared,offset
end

function _model_line_contains(line,coordinate::NTuple{3,Float64})
    parameter,offset=_model_line_parameter_exact(line,coordinate)
    cross=(
        line.delta[2]*offset[3]-line.delta[3]*offset[2],
        line.delta[3]*offset[1]-line.delta[1]*offset[3],
        line.delta[1]*offset[2]-line.delta[2]*offset[1],
    )
    return all(iszero,cross) && 0<=parameter<=1
end

@inline function _model_cross3(first,second)
    return (
        first[2]*second[3]-first[3]*second[2],
        first[3]*second[1]-first[1]*second[3],
        first[1]*second[2]-first[2]*second[1],
    )
end

function _model_normalize3(vector,caller::AbstractString,what::AbstractString)
    magnitude=hypot(vector...)
    (isfinite(magnitude) && magnitude>0) || throw(ArgumentError(
        "$caller: $what is degenerate"))
    normalized=ntuple(axis->vector[axis]/magnitude,3)
    all(isfinite,normalized) || throw(ArgumentError(
        "$caller: $what is not Float64-representable"))
    return normalized
end

function _model_plane_frame(
    m::GeoModel,tag::Int,caller::AbstractString)
    geometry=_model_plane_geometry(m,tag,caller)
    normal=(geometry.properties[1],geometry.properties[2],
            geometry.properties[3])
    rhs=geometry.properties[4]
    dominant=if abs(normal[1])>=abs(normal[2]) &&
                abs(normal[1])>=abs(normal[3])
        1
    elseif abs(normal[2])>=abs(normal[1]) &&
            abs(normal[2])>=abs(normal[3])
        2
    else
        3
    end
    R=Rational{BigInt}
    intercept=_model_rational_float(
        R(rhs)/R(normal[dominant]),caller,"Plane[$tag] frame origin")
    origin=ntuple(axis->axis==dominant ? intercept : 0.0,3)
    reference=normal[1]==0 ? (1.0,0.0,0.0) :
              normal[2]==0 ? (0.0,1.0,0.0) : (0.0,0.0,1.0)
    first_direction=_model_normalize3(
        _model_cross3(normal,reference),caller,"Plane[$tag] first direction")
    second_direction=_model_normalize3(
        _model_cross3(first_direction,normal),caller,
        "Plane[$tag] second direction")
    return merge(geometry,(
        normal=normal,
        origin=origin,
        first_direction=first_direction,
        second_direction=second_direction,
    ))
end

function _model_plane_point(
    plane,first_parameter,second_parameter,caller::AbstractString,
    point_index::Int)
    R=Rational{BigInt}
    u=R(first_parameter);v=R(second_parameter)
    return ntuple(3) do axis
        _model_rational_float(
            R(plane.origin[axis])+R(plane.first_direction[axis])*u+
            R(plane.second_direction[axis])*v,caller,
            "Plane value for point $point_index")
    end
end

function _model_plane_parameters_exact(
    plane,coordinate::NTuple{3,Float64})
    R=Rational{BigInt}
    offset=ntuple(axis->R(coordinate[axis])-R(plane.origin[axis]),3)
    first=ntuple(axis->R(plane.first_direction[axis]),3)
    second=ntuple(axis->R(plane.second_direction[axis]),3)
    g11=sum(component->component^2,first)
    g12=sum(axis->first[axis]*second[axis],1:3)
    g22=sum(component->component^2,second)
    b1=sum(axis->offset[axis]*first[axis],1:3)
    b2=sum(axis->offset[axis]*second[axis],1:3)
    determinant=g11*g22-g12^2
    determinant>0 || throw(ErrorException(
        "internal Plane frame is singular"))
    return ((b1*g22-b2*g12)/determinant,
            (b2*g11-b1*g12)/determinant)
end

function _model_plane_parameters(
    plane,coordinate::NTuple{3,Float64},caller::AbstractString,
    point_index::Int)
    exact=_model_plane_parameters_exact(plane,coordinate)
    return ntuple(axis->_model_rational_float(
        exact[axis],caller,"Plane parameter for point $point_index"),2)
end

function _model_plane_polygons(m::GeoModel,tag::Int,plane)
    first_axis,second_axis=plane.projection
    return [NTuple{2,Float64}[
        (m.points[point][first_axis],m.points[point][second_axis])
        for point in _loop_points(m,loop)] for loop in m.surfaces[tag]]
end

function _model_plane_contains(
    polygons,plane,coordinate::NTuple{3,Float64})
    if orient3(plane.anchor,plane.second,plane.third,coordinate)!=0
        return false
    end
    first_axis,second_axis=plane.projection
    projected=(coordinate[first_axis],coordinate[second_axis])
    _model_loop_position2(projected,first(polygons))==1 || return false
    return all(polygon->_model_loop_position2(projected,polygon)==0,
               Iterators.drop(polygons,1))
end

function _model_plane_parameter_bounds(plane,caller::AbstractString)
    parameters=NTuple{2,Float64}[
        _model_plane_parameters(plane,coordinate,caller,index)
        for (index,coordinate) in pairs(plane.coordinates)]
    return (
        minimum(parameter->parameter[1],parameters),
        minimum(parameter->parameter[2],parameters),
    ),(
        maximum(parameter->parameter[1],parameters),
        maximum(parameter->parameter[2],parameters),
    )
end

function _model_append_point!(output::Vector{Float64},point)
    append!(output,point)
    return output
end

"""
    model_value(model, dim, tag, parametric_coordinates) -> Vector{Float64}

Evaluate an explicit Point, straight Line, or Plane parametrization. Line
parameters use `[0,1]`; Plane parameters use the deterministic native orthonormal
frame returned by [`model_parametrization_bounds`](@ref).
"""
function model_value(m::GeoModel,dim,tag,parametric_coordinates)
    caller="model_value"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(0,1,2),caller)
    if dimension==0
        values=_model_evaluation_values(
            parametric_coordinates,1,caller,"parametric coordinates")
        isempty(values) || throw(ArgumentError(
            "$caller: Point parametric coordinates must be empty"))
        return collect(m.points[entity_tag])
    end
    stride=dimension
    values=_model_evaluation_values(
        parametric_coordinates,stride,caller,"parametric coordinates")
    output=Float64[]
    sizehint!(output,3*(length(values)÷stride))
    if dimension==1
        line=_model_line_geometry(m,entity_tag,caller)
        for (index,parameter) in pairs(values)
            _model_append_point!(output,
                _model_line_point(line,parameter,caller,index))
        end
    else
        plane=_model_plane_frame(m,entity_tag,caller)
        for index in 1:2:length(values)
            _model_append_point!(output,_model_plane_point(
                plane,values[index],values[index+1],caller,(index+1)÷2))
        end
    end
    return output
end

"""Evaluate first derivatives for an explicit straight Line or Plane."""
function model_derivative(m::GeoModel,dim,tag,parametric_coordinates)
    caller="model_derivative"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(1,2),caller)
    values=_model_evaluation_values(
        parametric_coordinates,dimension,caller,"parametric coordinates")
    output=Float64[]
    if dimension==1
        line=_model_line_geometry(m,entity_tag,caller)
        derivative=_model_line_derivative(line,entity_tag,caller)
        sizehint!(output,3length(values))
        for _ in values
            append!(output,derivative)
        end
    else
        plane=_model_plane_frame(m,entity_tag,caller)
        sizehint!(output,3length(values))
        for _ in 1:2:length(values)
            append!(output,plane.first_direction)
            append!(output,plane.second_direction)
        end
    end
    return output
end

"""Evaluate second derivatives for an explicit straight Line or Plane."""
function model_second_derivative(m::GeoModel,dim,tag,parametric_coordinates)
    caller="model_second_derivative"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(1,2),caller)
    values=_model_evaluation_values(
        parametric_coordinates,dimension,caller,"parametric coordinates")
    dimension==1 ? _model_line_geometry(m,entity_tag,caller) :
                   _model_plane_frame(m,entity_tag,caller)
    multiplier=dimension==1 ? 3 : 9
    return zeros(Float64,multiplier*(length(values)÷dimension))
end

"""Return zero curvature for an explicit straight Line or Plane."""
function model_curvature(m::GeoModel,dim,tag,parametric_coordinates)
    caller="model_curvature"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(1,2),caller)
    values=_model_evaluation_values(
        parametric_coordinates,dimension,caller,"parametric coordinates")
    dimension==1 ? _model_line_geometry(m,entity_tag,caller) :
                   _model_plane_frame(m,entity_tag,caller)
    return zeros(Float64,length(values)÷dimension)
end

"""
Return the two zero principal curvatures and their native orthonormal directions
for an explicit Plane.
"""
function model_principal_curvatures(m::GeoModel,tag,parametric_coordinates)
    caller="model_principal_curvatures"
    _,entity_tag=_model_metadata_entity(m,2,tag,caller)
    values=_model_evaluation_values(
        parametric_coordinates,2,caller,"parametric coordinates")
    plane=_model_plane_frame(m,entity_tag,caller)
    count=length(values)÷2
    first_directions=Float64[];second_directions=Float64[]
    sizehint!(first_directions,3count);sizehint!(second_directions,3count)
    for _ in 1:count
        append!(first_directions,plane.first_direction)
        append!(second_directions,plane.second_direction)
    end
    return zeros(Float64,count),zeros(Float64,count),
           first_directions,second_directions
end

"""Return the exterior-loop-oriented unit normal of an explicit Plane."""
function model_normal(m::GeoModel,tag,parametric_coordinates)
    caller="model_normal"
    _,entity_tag=_model_metadata_entity(m,2,tag,caller)
    values=_model_evaluation_values(
        parametric_coordinates,2,caller,"parametric coordinates")
    normal=_model_plane_frame(m,entity_tag,caller).normal
    output=Float64[]
    sizehint!(output,3*(length(values)÷2))
    for _ in 1:2:length(values)
        append!(output,normal)
    end
    return output
end

"""
    model_parametrization(model, dim, tag, coordinates) -> Vector{Float64}

Return orthogonal Line or Plane parameters for concatenated 3-D coordinates.
Coordinates need not lie inside the trimmed entity.
"""
function model_parametrization(m::GeoModel,dim,tag,coordinates)
    caller="model_parametrization"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(0,1,2),caller)
    values=_model_evaluation_values(coordinates,3,caller,"coordinates")
    if dimension==0
        length(values)==3 || throw(ArgumentError(
            "$caller: Point parametrization requires exactly one coordinate"))
        return Float64[]
    end
    output=Float64[]
    sizehint!(output,dimension*(length(values)÷3))
    line=dimension==1 ? _model_line_geometry(m,entity_tag,caller) : nothing
    plane=dimension==2 ? _model_plane_frame(m,entity_tag,caller) : nothing
    for index in 1:3:length(values)
        coordinate=(values[index],values[index+1],values[index+2])
        if dimension==1
            parameter,_=_model_line_parameter_exact(line,coordinate)
            push!(output,_model_rational_float(
                parameter,caller,"Line parameter for point $((index+2)÷3)"))
        else
            append!(output,_model_plane_parameters(
                plane,coordinate,caller,(index+2)÷3))
        end
    end
    return output
end

"""Return detached parametric lower and upper bounds for a Point, Line, or Plane."""
function model_parametrization_bounds(m::GeoModel,dim,tag)
    caller="model_parametrization_bounds"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(0,1,2),caller)
    dimension==0 && return Float64[],Float64[]
    dimension==1 && begin
        _model_line_geometry(m,entity_tag,caller)
        return [0.0],[1.0]
    end
    plane=_model_plane_frame(m,entity_tag,caller)
    lower,upper=_model_plane_parameter_bounds(plane,caller)
    return collect(lower),collect(upper)
end

"""
Count concatenated physical or parametric points inside an explicit native entity.
Line endpoints count as inside. Physical Plane coordinates use the trimmed interior,
excluding boundary loops; Plane parameters use the rectangular parameter bounds.
"""
function model_is_inside(m::GeoModel,dim,tag,coordinates,parametric=false)
    caller="model_is_inside"
    parametric isa Bool || throw(ArgumentError(
        "$caller: parametric must be Bool"))
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(0,1,2),caller)
    if dimension==0
        if parametric
            values=_model_evaluation_values(
                coordinates,1,caller,"parametric coordinates")
            isempty(values) || throw(ArgumentError(
                "$caller: Point parametric coordinates must be empty"))
            return 0
        end
        values=_model_evaluation_values(coordinates,3,caller,"coordinates")
        return 0
    end
    stride=parametric ? dimension : 3
    values=_model_evaluation_values(coordinates,stride,caller,
                                    parametric ? "parametric coordinates" :
                                                 "coordinates")
    if dimension==1
        if parametric
            return count(parameter->0<=parameter<=1,values)
        end
        line=_model_line_geometry(m,entity_tag,caller)
        return count(index->_model_line_contains(
            line,(values[index],values[index+1],values[index+2])),
            1:3:length(values))
    end
    plane=_model_plane_frame(m,entity_tag,caller)
    if parametric
        lower,upper=_model_plane_parameter_bounds(plane,caller)
        return count(index->
            lower[1]<=values[index]<=upper[1] &&
            lower[2]<=values[index+1]<=upper[2],1:2:length(values))
    end
    polygons=_model_plane_polygons(m,entity_tag,plane)
    return count(index->_model_plane_contains(
        polygons,plane,
        (values[index],values[index+1],values[index+2])),
        1:3:length(values))
end

"""
    model_closest_point(model, dim, tag, coordinates)

Project concatenated 3-D coordinates onto an explicit Line or Plane. Line
parameters are clamped to the trimmed segment; Plane projections can lie outside
its boundary loops.
"""
function model_closest_point(m::GeoModel,dim,tag,coordinates)
    caller="model_closest_point"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(1,2),caller)
    values=_model_evaluation_values(coordinates,3,caller,"coordinates")
    closest=Float64[];parameters=Float64[]
    sizehint!(closest,length(values))
    sizehint!(parameters,dimension*(length(values)÷3))
    line=dimension==1 ? _model_line_geometry(m,entity_tag,caller) : nothing
    plane=dimension==2 ? _model_plane_frame(m,entity_tag,caller) : nothing
    for index in 1:3:length(values)
        point_index=(index+2)÷3
        coordinate=(values[index],values[index+1],values[index+2])
        if dimension==1
            exact,_=_model_line_parameter_exact(line,coordinate)
            exact=clamp(exact,zero(exact),one(exact))
            parameter=_model_rational_float(
                exact,caller,"Line parameter for point $point_index")
            push!(parameters,parameter)
            _model_append_point!(closest,
                _model_line_point(line,exact,caller,point_index))
        else
            exact=_model_plane_parameters_exact(plane,coordinate)
            parameter=ntuple(axis->_model_rational_float(
                exact[axis],caller,"Plane parameter for point $point_index"),2)
            append!(parameters,parameter)
            _model_append_point!(closest,_model_plane_point(
                plane,exact[1],exact[2],caller,point_index))
        end
    end
    return closest,parameters
end
