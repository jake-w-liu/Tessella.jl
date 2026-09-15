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
    _model_require_line_curve(m,tag,caller,"curve evaluation")
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

# A contains/projection query flattens the boundary to a chord polygon; curved
# boundary curves would silently lose their bulge, so they fail explicitly.
function _model_plane_boundary_polygons(m::GeoModel,tag::Int,plane,
                                        caller::AbstractString)
    for curve in _model_surface_curves(m,tag)
        _model_require_line_curve(m,curve,caller,"plane boundary queries")
    end
    return _model_plane_polygons(m,tag,plane)
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

Evaluate an explicit Point, Line, arc, or Plane parametrization. Curve
parameters use `[0,1]`; Plane parameters use the deterministic native orthonormal
frame returned by [`model_parametrization_bounds`](@ref). `Circle`/`Ellipse`
arcs evaluate with Gmsh's angle parametrization.
"""
function model_value(m::GeoModel,dim,tag,parametric_coordinates)
    caller="model_value"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(0,1,2),caller)
    if dimension==0
        if haskey(m.discrete,(0,entity_tag))
            values=_model_evaluation_values(
                parametric_coordinates,1,caller,"parametric coordinates")
            isempty(values) || throw(ArgumentError(
                "$caller: Point parametric coordinates must be empty"))
            record=m.discrete[(0,entity_tag)]
            isempty(record.node_coords) && throw(ArgumentError(
                "$caller: discrete point has no node"))
            return collect(record.node_coords[:,1])
        end
        values=_model_evaluation_values(
            parametric_coordinates,1,caller,"parametric coordinates")
        isempty(values) || throw(ArgumentError(
            "$caller: Point parametric coordinates must be empty"))
        return collect(m.points[entity_tag])
    end
    if haskey(m.discrete,(dimension,entity_tag))
        values=_model_evaluation_values(
            parametric_coordinates,dimension,caller,"parametric coordinates")
        _,frames=_discrete_eval_setup(m,dimension,entity_tag,caller)
        output=Float64[]
        sizehint!(output,3*(length(values)÷dimension))
        if dimension==1
            for parameter in values
                _model_append_point!(output,
                    _discrete_curve_point(frames,parameter,caller))
            end
        else
            for index in 1:2:length(values)
                weights,points=_discrete_surface_locate(
                    frames,values[index],values[index+1],caller)
                _model_append_point!(output,ntuple(axis->sum(
                    k->weights[k]*points[k][axis],1:3),3))
            end
        end
        return output
    end
    stride=dimension
    values=_model_evaluation_values(
        parametric_coordinates,stride,caller,"parametric coordinates")
    output=Float64[]
    sizehint!(output,3*(length(values)÷stride))
    if dimension==1
        occ=_occ_geometry_checked(m,entity_tag,caller)
        if occ!==nothing && occ.occ===:circle
            for parameter in values
                _model_append_point!(output,
                    _occ_circle_point(occ,parameter))
            end
        elseif occ!==nothing && occ.occ===:degenerate
            point=m.points[m.curves[entity_tag][1]]
            for _ in values
                _model_append_point!(output,point)
            end
        elseif occ!==nothing && occ.occ===:line
            line=_model_line_geometry(m,entity_tag,caller)
            span=occ.t1-occ.t0
            for (index,parameter) in pairs(values)
                u=occ.t0==occ.t1 ? 0.0 : (parameter-occ.t0)/span
                _model_append_point!(output,
                    _model_line_point(line,u,caller,index))
            end
        elseif _curve_type(m,entity_tag)==:line
            line=_model_line_geometry(m,entity_tag,caller)
            for (index,parameter) in pairs(values)
                _model_append_point!(output,
                    _model_line_point(line,parameter,caller,index))
            end
        else
            arc=_arc_geometry(m,entity_tag,caller)
            for parameter in values
                _model_append_point!(output,_arc_point(arc,parameter))
            end
        end
    else
        surface_geometry=get(m.surface_geometry,entity_tag,nothing)
        if surface_geometry!==nothing && hasproperty(surface_geometry,:occ)
            for index in 1:2:length(values)
                _model_append_point!(output,_occ_surface_point(
                    surface_geometry,values[index],values[index+1]))
            end
            return output
        end
        plane=_model_plane_frame(m,entity_tag,caller)
        for index in 1:2:length(values)
            _model_append_point!(output,_model_plane_point(
                plane,values[index],values[index+1],caller,(index+1)÷2))
        end
    end
    return output
end

"""
Evaluate first derivatives for an explicit Line, arc, or Plane. Arc derivatives
are analytic; Gmsh's `getDerivative` uses a finite difference of the same
evaluation, so values agree within that error.
"""
function model_derivative(m::GeoModel,dim,tag,parametric_coordinates)
    caller="model_derivative"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(1,2),caller)
    values=_model_evaluation_values(
        parametric_coordinates,dimension,caller,"parametric coordinates")
    output=Float64[]
    if haskey(m.discrete,(dimension,entity_tag))
        _,frames=_discrete_eval_setup(m,dimension,entity_tag,caller)
        if dimension==1
            for parameter in values
                derivative=_discrete_curve_derivative(
                    frames,parameter,caller)
                append!(output,derivative)
            end
        else
            for index in 1:2:length(values)
                du,dv=_discrete_surface_derivative(
                    frames,values[index],values[index+1],caller)
                append!(output,du);append!(output,dv)
            end
        end
        return output
    end
    if dimension==1
        occ=_occ_geometry_checked(m,entity_tag,caller)
        if occ!==nothing && occ.occ===:circle
            sizehint!(output,3length(values))
            for parameter in values
                append!(output,_occ_circle_derivative(occ,parameter))
            end
        elseif occ!==nothing && occ.occ===:degenerate
            sizehint!(output,3length(values))
            for _ in values
                append!(output,(0.0,0.0,0.0))
            end
        elseif occ!==nothing && occ.occ===:line
            line=_model_line_geometry(m,entity_tag,caller)
            derivative=_model_line_derivative(line,entity_tag,caller)
            span=occ.t1-occ.t0
            span==0.0 || (derivative=(derivative[1]/span,
                                      derivative[2]/span,
                                      derivative[3]/span))
            sizehint!(output,3length(values))
            for _ in values
                append!(output,derivative)
            end
        elseif _curve_type(m,entity_tag)==:line
            line=_model_line_geometry(m,entity_tag,caller)
            derivative=_model_line_derivative(line,entity_tag,caller)
            sizehint!(output,3length(values))
            for _ in values
                append!(output,derivative)
            end
        else
            arc=_arc_geometry(m,entity_tag,caller)
            sizehint!(output,3length(values))
            for parameter in values
                append!(output,_arc_first_derivative(arc,parameter))
            end
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

"""
Evaluate second derivatives for an explicit Line, arc, or Plane. Line and Plane
derivatives are zero; arc derivatives are analytic.
"""
function model_second_derivative(m::GeoModel,dim,tag,parametric_coordinates)
    caller="model_second_derivative"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(1,2),caller)
    values=_model_evaluation_values(
        parametric_coordinates,dimension,caller,"parametric coordinates")
    if haskey(m.discrete,(dimension,entity_tag))
        multiplier=dimension==1 ? 3 : 9
        return zeros(Float64,multiplier*(length(values)÷dimension))
    end
    if dimension==1
        occ=_occ_geometry_checked(m,entity_tag,caller)
        if occ!==nothing && occ.occ===:circle
            output=Float64[]
            sizehint!(output,3*length(values))
            for parameter in values
                append!(output,_occ_circle_second_derivative(occ,parameter))
            end
            return output
        elseif occ!==nothing && occ.occ in (:line,:degenerate)
            return zeros(Float64,3*length(values))
        end
        _curve_type(m,entity_tag)!=:line && begin
            arc=_arc_geometry(m,entity_tag,caller)
            output=Float64[]
            sizehint!(output,3*length(values))
            for parameter in values
                append!(output,_arc_second_derivative(arc,parameter))
            end
            return output
        end
    end
    dimension==1 ? _model_line_geometry(m,entity_tag,caller) :
                   _model_plane_frame(m,entity_tag,caller)
    multiplier=dimension==1 ? 3 : 9
    return zeros(Float64,multiplier*(length(values)÷dimension))
end

"""
Return the curvature of an explicit Line (zero) or arc (`|P' × P''| / |P'|³`,
matching Gmsh's `getCurvature` formula).
"""
function model_curvature(m::GeoModel,dim,tag,parametric_coordinates)
    caller="model_curvature"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(1,2),caller)
    values=_model_evaluation_values(
        parametric_coordinates,dimension,caller,"parametric coordinates")
    if haskey(m.discrete,(dimension,entity_tag))
        return zeros(Float64,length(values)÷dimension)
    end
    if dimension==1
        occ=_occ_geometry_checked(m,entity_tag,caller)
        if occ!==nothing && occ.occ===:circle
            return fill(1.0/occ.r,length(values))
        elseif occ!==nothing && occ.occ in (:line,:degenerate)
            return zeros(Float64,length(values))
        end
        _curve_type(m,entity_tag)!=:line &&
            return Float64[
                _arc_curvature(_arc_geometry(m,entity_tag,caller),parameter)
                for parameter in values]
    end
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
    if haskey(m.discrete,(2,entity_tag))
        _,frames=_discrete_eval_setup(m,2,entity_tag,caller)
        output=Float64[]
        sizehint!(output,3*(length(values)÷2))
        for index in 1:2:length(values)
            _,points=_discrete_surface_locate(
                frames,values[index],values[index+1],caller)
            a,b,c=points
            normal=_model_normalize3(_model_cross3(
                (b[1]-a[1],b[2]-a[2],b[3]-a[3]),
                (c[1]-a[1],c[2]-a[2],c[3]-a[3])),caller,
                "discrete surface normal")
            append!(output,normal)
        end
        return output
    end
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
    if haskey(m.discrete,(dimension,entity_tag))
        _,frames=_discrete_eval_setup(m,dimension,entity_tag,caller)
        output=Float64[]
        sizehint!(output,dimension*(length(values)÷3))
        for index in 1:3:length(values)
            point=(values[index],values[index+1],values[index+2])
            if dimension==1
                push!(output,_discrete_curve_param(frames,point))
            else
                push!(output,_discrete_surface_param(frames,point)...)
            end
        end
        return output
    end
    output=Float64[]
    sizehint!(output,dimension*(length(values)÷3))
    occ=dimension==1 ? _occ_geometry_checked(m,entity_tag,caller) : nothing
    line=dimension==1 && occ===nothing ?
        _model_line_geometry(m,entity_tag,caller) : nothing
    occ!==nothing && occ.occ===:line &&
        (line=_model_line_geometry(m,entity_tag,caller))
    plane=dimension==2 ? _model_plane_frame(m,entity_tag,caller) : nothing
    for index in 1:3:length(values)
        coordinate=(values[index],values[index+1],values[index+2])
        if dimension==1
            if occ===nothing
                parameter,_=_model_line_parameter_exact(line,coordinate)
                push!(output,_model_rational_float(
                    parameter,caller,"Line parameter for point $((index+2)÷3)"))
            elseif occ.occ===:circle
                push!(output,_occ_circle_parameter(occ,coordinate))
            elseif occ.occ===:line
                parameter,_=_model_line_parameter_exact(line,coordinate)
                push!(output,occ.t0+_model_rational_float(
                    parameter,caller,"Line parameter for point $((index+2)÷3)")*
                    (occ.t1-occ.t0))
            else
                push!(output,occ.t0)
            end
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
    if haskey(m.discrete,(dimension,entity_tag))
        record=m.discrete[(dimension,entity_tag)]
        coords=_model_record_node_coords(m)
        frames=_discrete_element_frames(
            m,record,dimension,coords,caller)
        isempty(frames) && throw(ArgumentError(
            "$caller: discrete entity ($dimension,$entity_tag) has no elements"))
        lower=fill(Inf,dimension);upper=fill(-Inf,dimension)
        for (params,_) in frames
            for parameter in params
                for axis in 1:dimension
                    lower[axis]=min(lower[axis],parameter[axis])
                    upper[axis]=max(upper[axis],parameter[axis])
                end
            end
        end
        return lower,upper
    end
    dimension==1 && begin
        # Gmsh reports the stored parameter interval: [0,1] for built-in lines
        # and arcs, the OCC range (arc-length, angle, or degenerate) for
        # materialized primitive edges.
        occ=_occ_geometry_checked(m,entity_tag,caller)
        occ===nothing || return [occ.t0],[occ.t1]
        _curve_type(m,entity_tag)==:line &&
            _model_line_geometry(m,entity_tag,caller)
        return [0.0],[1.0]
    end
    surface_geometry=get(m.surface_geometry,entity_tag,nothing)
    if surface_geometry!==nothing && hasproperty(surface_geometry,:occ)
        # OCC faces report their analytic parameter intervals — azimuth and
        # profile (axis/slant length, latitude, or tube angle).
        return _occ_surface_bounds(surface_geometry)
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
    polygons=_model_plane_boundary_polygons(m,entity_tag,plane,caller)
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
    if haskey(m.discrete,(dimension,entity_tag))
        _,frames=_discrete_eval_setup(m,dimension,entity_tag,caller)
        for index in 1:3:length(values)
            point=(values[index],values[index+1],values[index+2])
            if dimension==1
                parameter=_discrete_curve_param(frames,point)
                push!(parameters,parameter)
                _model_append_point!(closest,
                    _discrete_curve_point(frames,parameter,caller))
            else
                uv=_discrete_surface_param(frames,point)
                append!(parameters,uv)
                weights,points=_discrete_surface_locate(
                    frames,uv[1],uv[2],caller)
                _model_append_point!(closest,ntuple(axis->sum(
                    k->weights[k]*points[k][axis],1:3),3))
            end
        end
        return closest,parameters
    end
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

function _model_reparametrization_selector(which,caller::AbstractString)
    which isa Integer || throw(ArgumentError(
        "$caller: which must be an integer"))
    which isa Bool && throw(ArgumentError(
        "$caller: which must not be Bool"))
    typemin(Int32)<=which<=typemax(Int32) || throw(ArgumentError(
        "$caller: which exceeds the Int32 range"))
    return Int(which)
end

"""
    model_reparametrize_on_surface(
        model, dim, tag, parametric_coordinates, surface_tag, which=0)

Map an explicit Point or straight-Line parametrization into the deterministic
parameters of an explicit Plane. The source need not belong to the Plane and is
orthogonally projected when it is not coplanar. `which` is validated but has no
effect because native Planes are not periodic.
"""
function model_reparametrize_on_surface(
    m::GeoModel,dim,tag,parametric_coordinates,surface_tag,which=0)
    caller="model_reparametrize_on_surface"
    dimension,entity_tag=_model_evaluation_entity(
        m,dim,tag,(0,1),caller)
    _,plane_tag=_model_metadata_entity(m,2,surface_tag,caller)
    _model_reparametrization_selector(which,caller)
    values=_model_evaluation_values(
        parametric_coordinates,1,caller,"parametric coordinates")
    if dimension==0
        isempty(values) || throw(ArgumentError(
            "$caller: Point parametric coordinates must be empty"))
        plane=_model_plane_frame(m,plane_tag,caller)
        return collect(_model_plane_parameters(
            plane,m.points[entity_tag],caller,1))
    end

    line=_model_line_geometry(m,entity_tag,caller)
    plane=_model_plane_frame(m,plane_tag,caller)
    output=Float64[]
    sizehint!(output,2length(values))
    for (index,parameter) in pairs(values)
        coordinate=_model_line_point(line,parameter,caller,index)
        append!(output,_model_plane_parameters(
            plane,coordinate,caller,index))
    end
    return output
end

# ---- discrete entity evaluation ------------------------------------------
#
# Parametrized discrete entities (see `create_geometry!`) evaluate through
# piecewise-linear interpolation: curves interpolate along their segments'
# stored parameters, surfaces interpolate inside the element whose parametric
# footprint contains the query point. Every helper resolves a node's
# parameters from the entity record's `node_params` (owned nodes) or
# `aux_params` (nodes owned by boundary entities).

# Parametric coordinates of `node` with respect to `record`'s entity.
function _record_param_of(record::DiscreteEntity,node::Int32,dim::Int)
    column=findfirst(==(node),record.node_tags)
    if column!==nothing
        size(record.node_params,1)>=dim || return nothing
        return Vector{Float64}(record.node_params[1:dim,column])
    end
    value=get(record.aux_params,node,nothing)
    value===nothing && return nothing
    length(value)>=dim || return nothing
    return value[1:dim]
end

# (corner tags, corner params, corner coords) for one element of `record`.
function _discrete_element_frame(m::GeoModel,record::DiscreteEntity,
                                 msh_type::Int,nodes::Vector{Int32},
                                 coords::Dict{Int32,NTuple{3,Float64}},
                                 caller::AbstractString)
    count=msh_family(msh_type)===:qua ? 4 :
          msh_family(msh_type)===:lin ? 2 : 3
    corners=nodes[1:min(count,length(nodes))]
    params=Vector{Vector{Float64}}(undef,length(corners))
    points=Vector{NTuple{3,Float64}}(undef,length(corners))
    dim=msh_dimension(msh_type)
    for (i,node) in enumerate(corners)
        value=_record_param_of(record,node,dim)
        value===nothing && throw(ArgumentError(
            "$caller: discrete entity has no parametrization; " *
            "run mesh.create_geometry first"))
        params[i]=value
        points[i]=get(coords,node,nothing)===nothing ? throw(ArgumentError(
            "$caller: node $node has no coordinates")) : coords[node]
    end
    return corners,params,points
end

# All (params, coords) corner frames of `record`'s dim-`dim` elements. Quads
# split into two corner triangles so parametric location is uniform.
function _discrete_element_frames(m::GeoModel,record::DiscreteEntity,dim::Int,
                                  coords::Dict{Int32,NTuple{3,Float64}},
                                  caller::AbstractString)
    frames=Tuple{Vector{Vector{Float64}},Vector{NTuple{3,Float64}}}[]
    for (index,msh_type) in enumerate(record.element_types)
        msh_dimension(msh_type)!=dim && continue
        _,params,points=_discrete_element_frame(
            m,record,Int(msh_type),record.element_nodes[index],coords,caller)
        if msh_family(msh_type)===:qua
            push!(frames,([params[1],params[2],params[3]],
                          [points[1],points[2],points[3]]))
            push!(frames,([params[1],params[3],params[4]],
                          [points[1],points[3],points[4]]))
        else
            push!(frames,(params,points))
        end
    end
    return frames
end

# Barycentric coordinates of `p` in triangle (a,b,c) — degenerate triangles
# return nothing.
function _barycentric3(p,a,b,c)
    v0=(b[1]-a[1],b[2]-a[2]);v1=(c[1]-a[1],c[2]-a[2])
    v2=(p[1]-a[1],p[2]-a[2])
    d00=v0[1]*v0[1]+v0[2]*v0[2];d01=v0[1]*v1[1]+v0[2]*v1[2]
    d11=v1[1]*v1[1]+v1[2]*v1[2]
    d20=v2[1]*v0[1]+v2[2]*v0[2];d21=v2[1]*v1[1]+v2[2]*v1[2]
    denom=d00*d11-d01*d01
    denom==0 && return nothing
    v=(d11*d20-d01*d21)/denom
    w=(d00*d21-d01*d20)/denom
    return (1-v-w,v,w)
end

# Closest point of `p` on segment (a,b) → (clamped t, squared distance).
function _segment_closest(p,a,b)
    d=(b[1]-a[1],b[2]-a[2],b[3]-a[3])
    l2=d[1]^2+d[2]^2+d[3]^2
    if l2==0
        return 0.0,(p[1]-a[1])^2+(p[2]-a[2])^2+(p[3]-a[3])^2
    end
    t=((p[1]-a[1])*d[1]+(p[2]-a[2])*d[2]+(p[3]-a[3])*d[3])/l2
    t=clamp(t,0.0,1.0)
    q=(a[1]+t*d[1],a[2]+t*d[2],a[3]+t*d[3])
    return t,(p[1]-q[1])^2+(p[2]-q[2])^2+(p[3]-q[3])^2
end

# Closest point of `p` on triangle (a,b,c) → (barycentric, squared distance).
function _triangle_closest(p,a,b,c)
    ab=(b[1]-a[1],b[2]-a[2],b[3]-a[3])
    ac=(c[1]-a[1],c[2]-a[2],c[3]-a[3])
    ap=(p[1]-a[1],p[2]-a[2],p[3]-a[3])
    d1=ab[1]*ap[1]+ab[2]*ap[2]+ab[3]*ap[3]
    d2=ac[1]*ap[1]+ac[2]*ap[2]+ac[3]*ap[3]
    if d1<=0 && d2<=0
        return (1.0,0.0,0.0),ap[1]^2+ap[2]^2+ap[3]^2
    end
    bp=(p[1]-b[1],p[2]-b[2],p[3]-b[3])
    d3=ab[1]*bp[1]+ab[2]*bp[2]+ab[3]*bp[3]
    d4=ac[1]*bp[1]+ac[2]*bp[2]+ac[3]*bp[3]
    if d3>=0 && d4<=d3
        return (0.0,1.0,0.0),bp[1]^2+bp[2]^2+bp[3]^2
    end
    vc=d1*d4-d3*d2
    if vc<=0 && d1>=0 && d3<=0
        v=d1/(d1-d3)
        q=(a[1]+v*ab[1],a[2]+v*ab[2],a[3]+v*ab[3])
        return (1-v,v,0.0),(p[1]-q[1])^2+(p[2]-q[2])^2+(p[3]-q[3])^2
    end
    cp=(p[1]-c[1],p[2]-c[2],p[3]-c[3])
    d5=ab[1]*cp[1]+ab[2]*cp[2]+ab[3]*cp[3]
    d6=ac[1]*cp[1]+ac[2]*cp[2]+ac[3]*cp[3]
    if d6>=0 && d5<=d6
        return (0.0,0.0,1.0),cp[1]^2+cp[2]^2+cp[3]^2
    end
    vb=d5*d2-d1*d6
    if vb<=0 && d2>=0 && d6<=0
        w=d2/(d2-d6)
        q=(a[1]+w*ac[1],a[2]+w*ac[2],a[3]+w*ac[3])
        return (1-w,0.0,w),(p[1]-q[1])^2+(p[2]-q[2])^2+(p[3]-q[3])^2
    end
    va=d3*d6-d5*d4
    if va<=0 && (d4-d3)>=0 && (d5-d6)>=0
        w=(d4-d3)/((d4-d3)+(d5-d6))
        q=(b[1]+w*(c[1]-b[1]),b[2]+w*(c[2]-b[2]),b[3]+w*(c[3]-b[3]))
        return (0.0,1-w,w),(p[1]-q[1])^2+(p[2]-q[2])^2+(p[3]-q[3])^2
    end
    denom=1/(va+vb+vc)
    v=vb*denom;w=vc*denom
    q=(a[1]+ab[1]*v+ac[1]*w,a[2]+ab[2]*v+ac[2]*w,a[3]+ab[3]*v+ac[3]*w)
    return (1-v-w,v,w),(p[1]-q[1])^2+(p[2]-q[2])^2+(p[3]-q[3])^2
end

# Locate the param-space footprint containing (u,v); falls back to the
# barycentric-clamped nearest triangle when the query lies outside.
function _discrete_surface_locate(frames,u::Float64,v::Float64,
                                  caller::AbstractString)
    isempty(frames) && throw(ArgumentError(
        "$caller: discrete surface has no elements"))
    best=nothing;best_distance=Inf
    for (params,points) in frames
        weights=_barycentric3((u,v),params[1],params[2],params[3])
        weights===nothing && continue
        if all(w->w>=-1e-12,weights)
            return weights,points
        end
        positive=(max(weights[1],0.0),max(weights[2],0.0),max(weights[3],0.0))
        total=positive[1]+positive[2]+positive[3]
        total>0 || continue
        # single assignment (a reassigned `clamped` captured below would be boxed)
        clamped=(positive[1]/total,positive[2]/total,positive[3]/total)
        pu=clamped[1]*params[1][1]+clamped[2]*params[2][1]+clamped[3]*params[3][1]
        pv=clamped[1]*params[1][2]+clamped[2]*params[2][2]+clamped[3]*params[3][2]
        distance=(pu-u)^2+(pv-v)^2
        if distance<best_distance
            best_distance=distance
            best=(clamped,points)
        end
    end
    best===nothing && throw(ArgumentError(
        "$caller: could not locate ($u, $v) on the discrete surface"))
    return best
end

# Shared frame lookup for a parametrized discrete entity.
function _discrete_eval_setup(m::GeoModel,dimension::Int,tag::Int,
                              caller::AbstractString)
    record=get(m.discrete,(dimension,tag),nothing)
    record===nothing && throw(ArgumentError(
        "$caller: entity ($dimension,$tag) is not discrete"))
    coords=_model_record_node_coords(m)
    frames=_discrete_element_frames(m,record,dimension,coords,caller)
    isempty(frames) && throw(ArgumentError(
        "$caller: discrete entity ($dimension,$tag) has no elements"))
    return record,frames
end

# Evaluate a discrete curve at normalized parameter `t` in [0,1].
function _discrete_curve_point(frames,t::Float64,caller::AbstractString)
    best=nothing;best_distance=Inf
    for (params,points) in frames
        pa,pb=params[1][1],params[2][1]
        lo,hi=minmax(pa,pb)
        if lo-1e-12<=t<=hi+1e-12 && hi>lo
            local_t=clamp((t-pa)/(pb-pa),0.0,1.0)
            a,b=points[1],points[2]
            return (a[1]+local_t*(b[1]-a[1]),
                    a[2]+local_t*(b[2]-a[2]),
                    a[3]+local_t*(b[3]-a[3]))
        end
        distance=min(abs(t-pa),abs(t-pb))
        if distance<best_distance
            best_distance=distance
            best=(pa<pb ? (points[1],points[2],pa,pb) :
                         (points[2],points[1],pb,pa))
        end
    end
    best===nothing && throw(ArgumentError(
        "$caller: parameter $t is outside the curve's range"))
    a,b,lo,hi=best
    local_t=clamp((t-lo)/(hi-lo),0.0,1.0)
    return (a[1]+local_t*(b[1]-a[1]),a[2]+local_t*(b[2]-a[2]),
            a[3]+local_t*(b[3]-a[3]))
end

# Inverse map: closest point on the discrete curve → its stored parameter.
function _discrete_curve_param(frames,point)
    best=0.0;best_distance=Inf
    for (params,points) in frames
        t,dist=_segment_closest(point,points[1],points[2])
        if dist<best_distance
            best_distance=dist
            best=params[1][1]+t*(params[2][1]-params[1][1])
        end
    end
    return best
end

# Segment direction d(xyz)/dt of the discrete-curve frame containing `t`.
function _discrete_curve_derivative(frames,t::Float64,caller::AbstractString)
    best=nothing;best_distance=Inf
    for (params,points) in frames
        pa,pb=params[1][1],params[2][1]
        lo,hi=minmax(pa,pb)
        if lo-1e-12<=t<=hi+1e-12 && hi>lo
            a,b=points[1],points[2]
            scale=(pb-pa)
            return ((b[1]-a[1])/scale,(b[2]-a[2])/scale,(b[3]-a[3])/scale)
        end
        distance=min(abs(t-pa),abs(t-pb))
        if distance<best_distance
            best_distance=distance
            best=(params,points)
        end
    end
    best===nothing && throw(ArgumentError(
        "$caller: parameter $t is outside the curve's range"))
    params,points=best
    pa,pb=params[1][1],params[2][1]
    pb==pa && throw(ArgumentError("$caller: degenerate curve segment"))
    a,b=points[1],points[2]
    scale=pb-pa
    return ((b[1]-a[1])/scale,(b[2]-a[2])/scale,(b[3]-a[3])/scale)
end

# d(xyz)/d(u,v) of the discrete-surface element containing (u,v).
function _discrete_surface_derivative(frames,u::Float64,v::Float64,
                                      caller::AbstractString)
    isempty(frames) && throw(ArgumentError(
        "$caller: discrete surface has no elements"))
    best=nothing;best_distance=Inf
    for (params,points) in frames
        weights=_barycentric3((u,v),params[1],params[2],params[3])
        weights===nothing && continue
        distance=weights[1]<0||weights[2]<0||weights[3]<0 ?
            min(abs.(weights)...) : 0.0
        pu=params[1];pv=params[2];pw=params[3]
        du=(pv[1]-pu[1],pw[1]-pu[1]);dv=(pv[2]-pu[2],pw[2]-pu[2])
        determinant=du[1]*dv[2]-du[2]*dv[1]
        determinant==0 && continue
        pa,pb,pc=points
        ex=(pb[1]-pa[1],pb[2]-pa[2],pb[3]-pa[3])
        ey=(pc[1]-pa[1],pc[2]-pa[2],pc[3]-pa[3])
        # Solve [du; dv] * J = [ex; ey] for J rows (dxyz/du, dxyz/dv).
        det=determinant
        d_du=ntuple(axis->( dv[2]*ex[axis]-dv[1]*ey[axis])/det,3)
        d_dv=ntuple(axis->(-du[2]*ex[axis]+du[1]*ey[axis])/det,3)
        if distance<=best_distance
            best_distance=distance
            best=(d_du,d_dv)
            distance==0 && break
        end
    end
    best===nothing && throw(ArgumentError(
        "$caller: could not differentiate at ($u, $v)"))
    return best
end

# Inverse map on a discrete surface: closest xyz point → stored parameters.
function _discrete_surface_param(frames,point)
    best=(0.0,0.0);best_distance=Inf
    for (params,points) in frames
        weights,dist=_triangle_closest(point,points[1],points[2],points[3])
        if dist<best_distance
            best_distance=dist
            best=(sum(k->weights[k]*params[k][1],1:3),
                  sum(k->weights[k]*params[k][2],1:3))
        end
    end
    return best
end
