"""
    Mesh1D

One-dimensional curve meshing under a size field. The default integration,
grading, point-count and close-point filtering policy follows Gmsh 4.15.2
`meshGEdge`: adaptive trapezoidal integration, scalar size smoothing for the
standard 2-D algorithm, and equal increments of the resulting primitive. An
explicitly supplied `nsample` selects a uniform integration grid while keeping
the same grading and count policy.

The curve is a callable `γ(t) -> (x,y,z)`. A callable `derivative(t)` can be
supplied for exact tangent evaluation; otherwise a bounded second-order finite
difference is used. `mesh_segment` always uses its exact constant tangent.
"""
module Mesh1D

using ..SizeField: AbstractSizeField, Metric3, size_at, metric_at,
                   _metric_displacement

export mesh_curve, mesh_segment, curve_length, metric_length

const _GMSH_INTEGRATION_PRECISION = 1.0e-9
const _GMSH_MIN_INTEGRATION_DEPTH = 7
const _GMSH_MAX_INTEGRATION_DEPTH = 26
const _GMSH_SMOOTH_RATIO = 1.8
const _GMSH_SMOOTH_ITERATIONS = 2001
const _DEFAULT_MAX_INTEGRATION_POINTS = 1_000_000
const _DEFAULT_MAX_EDGES = 10_000_000

@inline function _pt3(p, caller::AbstractString="Mesh1D")
    n = try
        length(p)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: a point must be an indexable coordinate collection"))
    end
    n >= 2 || throw(ArgumentError("$caller: a point needs at least two coordinates"))
    q = try
        (Float64(p[1]), Float64(p[2]), n >= 3 ? Float64(p[3]) : 0.0)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: coordinates must be real and Float64-convertible: $(sprint(showerror, err))"))
    end
    (isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3])) ||
        throw(ArgumentError("$caller: coordinates are not finite: $q"))
    return q
end

@inline _dist3(a, b) = hypot(a[1] - b[1], a[2] - b[2], a[3] - b[3])

function _optional_positive_int(value, caller::AbstractString, name::AbstractString)
    value === nothing && return nothing
    value isa Integer || throw(ArgumentError(
        "$caller: $name must be nothing or a positive integer"))
    value > 0 || throw(ArgumentError("$caller: $name must be positive (got $value)"))
    value <= typemax(Int) || throw(ArgumentError(
        "$caller: $name exceeds the platform Int limit (got $value)"))
    return Int(value)
end

function _positive_int(value, caller::AbstractString, name::AbstractString)
    result = _optional_positive_int(value, caller, name)
    result === nothing && throw(ArgumentError("$caller: $name must be a positive integer"))
    return result
end

function _finite_positive(value::Real, caller::AbstractString, name::AbstractString)
    result = try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be representable as Float64"))
    end
    (isfinite(result) && result > 0) || throw(ArgumentError(
        "$caller: $name must be finite and positive (got $value)"))
    return result
end

function _curve_args(t0::Real, t1::Real, nsample, caller::AbstractString)
    a = try
        Float64(t0)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: t0 must be representable as Float64"))
    end
    b = try
        Float64(t1)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: t1 must be representable as Float64"))
    end
    (isfinite(a) && isfinite(b)) || throw(ArgumentError(
        "$caller: t0 and t1 must be finite (got $t0, $t1)"))
    a < b || throw(ArgumentError("$caller: require t0 < t1 (got $a, $b)"))
    isfinite(b - a) || throw(ArgumentError(
        "$caller: parameter span t1-t0 overflowed Float64"))
    return a, b, _optional_positive_int(nsample, caller, "nsample")
end

function _integration_args(integration_precision::Real, max_integration_points,
                           min_integration_depth, max_integration_depth,
                           caller::AbstractString)
    precision = _finite_positive(integration_precision, caller,
                                 "integration_precision")
    maxpoints = _positive_int(max_integration_points, caller,
                              "max_integration_points")
    mindepth = _positive_int(min_integration_depth, caller,
                             "min_integration_depth")
    maxdepth = _positive_int(max_integration_depth, caller,
                             "max_integration_depth")
    mindepth <= maxdepth || throw(ArgumentError(
        "$caller: min_integration_depth must not exceed max_integration_depth"))
    # At greater depths adjacent Float64 parameters can no longer be guaranteed
    # distinct even on an ordinary O(1) interval.
    maxdepth <= 52 || throw(ArgumentError(
        "$caller: max_integration_depth must be at most 52"))
    return precision, maxpoints, mindepth, maxdepth
end

function _endpoint_vertex_entity(value, caller::AbstractString, which::AbstractString)
    value === nothing && return nothing
    (value isa Tuple && length(value) == 2 &&
     value[1] isa Integer && value[2] isa Integer) || throw(ArgumentError(
        "$caller: $which endpoint entity must be nothing or a (0, tag) integer tuple"))
    dim = try
        Int(value[1])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: $which endpoint entity dimension is outside the platform Int range"))
    end
    tag = try
        Int(value[2])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: $which endpoint entity tag is outside the platform Int range"))
    end
    dim == 0 || throw(ArgumentError(
        "$caller: $which endpoint entity must have dimension 0 (got $dim)"))
    tag > 0 || throw(ArgumentError(
        "$caller: $which endpoint entity tag must be positive (got $tag)"))
    return (dim, tag)
end

function _endpoint_vertex_entities(value, caller::AbstractString)
    value === nothing && return (nothing, nothing)
    (value isa Tuple && length(value) == 2) || throw(ArgumentError(
        "$caller: endpoint_entities must be nothing or a (begin, end) tuple"))
    return (_endpoint_vertex_entity(value[1], caller, "begin"),
            _endpoint_vertex_entity(value[2], caller, "end"))
end

@inline function _checked_distance(a, b, caller::AbstractString)
    d = _dist3(a, b)
    isfinite(d) || throw(ArgumentError("$caller: sampled segment length is not finite"))
    return d
end

@inline function _half_product(a::Float64, b::Float64, caller::AbstractString)
    (a >= 0 && b >= 0 && isfinite(a) && isfinite(b)) || throw(ArgumentError(
        "$caller: trapezoidal integrand is not finite and nonnegative"))
    value = a >= b ? (a / 2) * b : a * (b / 2)
    isfinite(value) || throw(ArgumentError("$caller: trapezoidal integral overflowed Float64"))
    return value
end

@inline function _trapezoidal(a, b, caller::AbstractString)
    dt = b.t - a.t
    dt > 0 || throw(ArgumentError("$caller: integration parameters are not increasing"))
    value = _half_product(dt, a.lc, caller) + _half_product(dt, b.lc, caller)
    isfinite(value) || throw(ArgumentError("$caller: trapezoidal integral overflowed Float64"))
    return value
end

function _numeric_tangent(γ, point, t::Float64, t0::Float64, t1::Float64,
                          caller::AbstractString)
    span = t1 - t0
    h = cbrt(eps(Float64)) * span
    if !(h > 0) || t + h == t || t - h == t
        left = _pt3(γ(t0), caller)
        right = _pt3(γ(t1), caller)
        return ((right[1] - left[1]) / span,
                (right[2] - left[2]) / span,
                (right[3] - left[3]) / span)
    end
    if t - h >= t0 && t + h <= t1
        left = _pt3(γ(t - h), caller)
        right = _pt3(γ(t + h), caller)
        width = 2h
        tangent = ((right[1] - left[1]) / width,
                   (right[2] - left[2]) / width,
                   (right[3] - left[3]) / width)
    elseif t + 2h <= t1
        p1 = _pt3(γ(t + h), caller)
        p2 = _pt3(γ(t + 2h), caller)
        width = 2h
        tangent = ((4 * (p1[1] - point[1]) - (p2[1] - point[1])) / width,
                   (4 * (p1[2] - point[2]) - (p2[2] - point[2])) / width,
                   (4 * (p1[3] - point[3]) - (p2[3] - point[3])) / width)
    elseif t - 2h >= t0
        p1 = _pt3(γ(t - h), caller)
        p2 = _pt3(γ(t - 2h), caller)
        width = 2h
        tangent = ((4 * (point[1] - p1[1]) - (point[1] - p2[1])) / width,
                   (4 * (point[2] - p1[2]) - (point[2] - p2[2])) / width,
                   (4 * (point[3] - p1[3]) - (point[3] - p2[3])) / width)
    else
        left_t = max(t0, t - h)
        right_t = min(t1, t + h)
        left_t < right_t || throw(ArgumentError(
            "$caller: cannot form distinct parameters for numerical differentiation"))
        left = _pt3(γ(left_t), caller)
        right = _pt3(γ(right_t), caller)
        width = right_t - left_t
        tangent = ((right[1] - left[1]) / width,
                   (right[2] - left[2]) / width,
                   (right[3] - left[3]) / width)
    end
    (isfinite(tangent[1]) && isfinite(tangent[2]) && isfinite(tangent[3])) ||
        throw(ArgumentError(
        "$caller: numerical curve derivative is not finite at parameter $t"))
    if tangent[1] == 0 && tangent[2] == 0 && tangent[3] == 0
        # Local coordinate differences can underflow even when the derivative
        # is representable (e.g. γ(t)=1e-320*t). Use the full secant only when
        # it contains information that the local stencil lost.
        first_point = _pt3(γ(t0), caller)
        last_point = _pt3(γ(t1), caller)
        if first_point != last_point
            tangent = ((last_point[1] - first_point[1]) / span,
                       (last_point[2] - first_point[2]) / span,
                       (last_point[3] - first_point[3]) / span)
            (isfinite(tangent[1]) && isfinite(tangent[2]) &&
             isfinite(tangent[3])) || throw(ArgumentError(
                "$caller: fallback curve derivative is not finite at parameter $t"))
        end
    end
    return tangent
end

function _curve_tangent(γ, derivative, point, t::Float64, t0::Float64,
                        t1::Float64, caller::AbstractString)
    tangent = derivative === nothing ?
        _numeric_tangent(γ, point, t, t0, t1, caller) :
        _pt3(derivative(t), caller)
    (isfinite(tangent[1]) && isfinite(tangent[2]) && isfinite(tangent[3])) ||
        throw(ArgumentError(
        "$caller: curve derivative is not finite at parameter $t"))
    return tangent
end

struct _IntegrationPoint
    t::Float64
    lc::Float64
    p::Float64
    xp::Float64
    h::Float64
end

@inline _with_primitive(point::_IntegrationPoint, primitive::Float64) =
    _IntegrationPoint(point.t, point.lc, primitive, point.xp, point.h)

@inline _with_size(point::_IntegrationPoint, size::Float64, lc::Float64) =
    _IntegrationPoint(point.t, lc, point.p, point.xp, size)

function _length_point(γ, derivative, t::Float64, t0::Float64, t1::Float64,
                       caller::AbstractString)
    point = _pt3(γ(t), caller)
    tangent = _curve_tangent(γ, derivative, point, t, t0, t1, caller)
    speed = hypot(tangent[1], tangent[2], tangent[3])
    isfinite(speed) || throw(ArgumentError("$caller: curve speed overflowed Float64"))
    return _IntegrationPoint(t, speed, 0.0, speed, 1.0)
end

function _metric_point(γ, derivative, field::AbstractSizeField, t::Float64,
                       t0::Float64, t1::Float64, curve_entity, begin_entity,
                       end_entity, at_begin::Bool, at_end::Bool,
                       anisotropic_metric::Bool, caller::AbstractString)
    point = _pt3(γ(t), caller)
    tangent = _curve_tangent(γ, derivative, point, t, t0, t1, caller)
    speed = hypot(tangent[1], tangent[2], tangent[3])
    isfinite(speed) || throw(ArgumentError("$caller: curve speed overflowed Float64"))

    if anisotropic_metric
        context = at_begin && begin_entity !== nothing ? begin_entity :
                  at_end && end_entity !== nothing ? end_entity : curve_entity
        metric = metric_at(field, point[1], point[2], point[3], context)
        metric isa Metric3 || throw(ArgumentError(
            "$caller: metric_at must return Metric3"))
        lc = _metric_displacement(metric, tangent[1], tangent[2], tangent[3], caller)
        return _IntegrationPoint(t, lc, 0.0, speed, 1.0)
    end

    h = size_at(field, point[1], point[2], point[3], curve_entity)
    at_begin && begin_entity !== nothing &&
        (h = min(h, size_at(field, point[1], point[2], point[3], begin_entity)))
    at_end && end_entity !== nothing &&
        (h = min(h, size_at(field, point[1], point[2], point[3], end_entity)))
    lc = speed == 0 ? 0.0 : speed / h
    isfinite(lc) || throw(ArgumentError(
        "$caller: scalar metric integrand overflowed Float64 at parameter $t"))
    return _IntegrationPoint(t, lc, 0.0, speed, h)
end

function _push_point!(points::Vector{_IntegrationPoint}, point::_IntegrationPoint,
                      maxpoints::Int, caller::AbstractString)
    length(points) < maxpoints || throw(ArgumentError(
        "$caller: adaptive integration exceeded max_integration_points=$maxpoints"))
    push!(points, point)
    return nothing
end

function _recursive_integration!(points::Vector{_IntegrationPoint}, from,
                                 to, evaluate, precision::Float64,
                                 mindepth::Int, maxdepth::Int, maxpoints::Int,
                                 depth::Int, caller::AbstractString)
    midpoint_t = from.t / 2 + to.t / 2
    (from.t < midpoint_t < to.t) || throw(ArgumentError(
        "$caller: adaptive integration cannot refine the Float64 parameter interval further"))
    midpoint = evaluate(midpoint_t, false, false)
    coarse = _trapezoidal(from, to, caller)
    left = _trapezoidal(from, midpoint, caller)
    right = _trapezoidal(midpoint, to, caller)
    fine = left + right
    isfinite(fine) || throw(ArgumentError("$caller: adaptive integral overflowed Float64"))
    error_estimate = abs(coarse - fine)

    if (error_estimate < precision && depth >= mindepth) || depth >= maxdepth
        previous = points[end].p
        midpoint_primitive = previous + left
        isfinite(midpoint_primitive) || throw(ArgumentError(
            "$caller: accumulated integral overflowed Float64"))
        midpoint = _with_primitive(midpoint, midpoint_primitive)
        _push_point!(points, midpoint, maxpoints, caller)
        to_primitive = midpoint.p + right
        isfinite(to_primitive) || throw(ArgumentError(
            "$caller: accumulated integral overflowed Float64"))
        to = _with_primitive(to, to_primitive)
        _push_point!(points, to, maxpoints, caller)
        return nothing
    end

    _recursive_integration!(points, from, midpoint, evaluate, precision,
                            mindepth, maxdepth, maxpoints, depth + 1, caller)
    _recursive_integration!(points, midpoint, to, evaluate, precision,
                            mindepth, maxdepth, maxpoints, depth + 1, caller)
    return nothing
end

function _adaptive_points(evaluate, t0::Float64, t1::Float64,
                          precision::Float64, maxpoints::Int,
                          mindepth::Int, maxdepth::Int,
                          caller::AbstractString)
    maxpoints >= 3 || throw(ArgumentError(
        "$caller: max_integration_points must be at least 3"))
    from = evaluate(t0, true, false)
    to = evaluate(t1, false, true)
    points = _IntegrationPoint[]
    sizehint!(points, min(maxpoints, 257))
    push!(points, from)
    _recursive_integration!(points, from, to, evaluate, precision, mindepth,
                            maxdepth, maxpoints, 1, caller)
    return points
end

function _uniform_points(evaluate, t0::Float64, t1::Float64, nsample::Int,
                         caller::AbstractString)
    nsample <= typemax(Int) - 1 || throw(ArgumentError(
        "$caller: nsample exceeds the allocatable vector length"))
    points = Vector{_IntegrationPoint}(undef, nsample + 1)
    points[1] = evaluate(t0, true, false)
    @inbounds for i in 1:nsample
        t = t0 + (t1 - t0) * i / nsample
        current = evaluate(t, false, i == nsample)
        increment = _trapezoidal(points[i], current, caller)
        primitive = points[i].p + increment
        isfinite(primitive) || throw(ArgumentError(
            "$caller: accumulated integral overflowed Float64"))
        points[i + 1] = _with_primitive(current, primitive)
    end
    return points
end

function _integration_points(evaluate, t0::Float64, t1::Float64, nsample,
                             precision::Float64, maxpoints::Int,
                             mindepth::Int, maxdepth::Int,
                             caller::AbstractString)
    nsample === nothing && return _adaptive_points(
        evaluate, t0, t1, precision, maxpoints, mindepth, maxdepth, caller)
    nsample <= typemax(Int) - 1 || throw(ArgumentError(
        "$caller: nsample exceeds the allocatable vector length"))
    required = nsample + 1
    required <= maxpoints || throw(ArgumentError(
        "$caller: nsample=$nsample requires $required integration points, exceeding max_integration_points=$maxpoints"))
    return _uniform_points(evaluate, t0, t1, nsample, caller)
end

@inline function _point_size(point::_IntegrationPoint)
    point.xp == 0 && return point.h
    value = point.xp / point.lc
    (isfinite(value) && value > 0) || return point.h
    return value
end

function _recompute_primitive!(points::Vector{_IntegrationPoint},
                               caller::AbstractString)
    points[1] = _with_primitive(points[1], 0.0)
    @inbounds for i in 2:length(points)
        primitive = points[i - 1].p +
                    _trapezoidal(points[i - 1], points[i], caller)
        isfinite(primitive) || throw(ArgumentError(
            "$caller: accumulated smoothed integral overflowed Float64"))
        points[i] = _with_primitive(points[i], primitive)
    end
    return points[end].p
end

function _smooth_primitive!(points::Vector{_IntegrationPoint}, smooth_ratio::Real,
                            max_iterations, caller::AbstractString)
    ratio = _finite_positive(smooth_ratio, caller, "smooth_ratio")
    ratio >= 1 || throw(ArgumentError(
        "$caller: smooth_ratio must be at least 1 (got $smooth_ratio)"))
    iterations = _positive_int(max_iterations, caller, "max_smoothing_iterations")
    alpha_minus_one = sqrt(ratio) - 1

    for _ in 1:iterations
        changed = false
        @inbounds for i in 2:length(points)
            previous = points[i - 1]
            current = points[i]
            hprevious = _point_size(previous)
            hcurrent = _point_size(current)
            dt = current.t - previous.t
            limit = dt * alpha_minus_one * current.xp
            if hcurrent - hprevious > limit * 1.01
                hnew = hprevious + limit
                (isfinite(hnew) && hnew > 0) || throw(ArgumentError(
                    "$caller: scalar size smoothing produced an invalid size"))
                lcnew = current.xp == 0 ? 0.0 : current.xp / hnew
                isfinite(lcnew) || throw(ArgumentError(
                    "$caller: smoothed metric integrand overflowed Float64"))
                points[i] = _with_size(current, hnew, lcnew)
                changed = true
            end
        end

        @inbounds for i in length(points):-1:2
            previous = points[i - 1]
            current = points[i]
            hprevious = _point_size(previous)
            hcurrent = _point_size(current)
            dt = current.t - previous.t
            limit = dt * alpha_minus_one * previous.xp
            if hprevious - hcurrent > limit * 1.01
                hnew = hcurrent + limit
                (isfinite(hnew) && hnew > 0) || throw(ArgumentError(
                    "$caller: scalar size smoothing produced an invalid size"))
                lcnew = previous.xp == 0 ? 0.0 : previous.xp / hnew
                isfinite(lcnew) || throw(ArgumentError(
                    "$caller: smoothed metric integrand overflowed Float64"))
                points[i - 1] = _with_size(previous, hnew, lcnew)
                changed = true
            end
        end
        changed || break
    end
    return _recompute_primitive!(points, caller)
end

function _metric_points(γ, field::AbstractSizeField, derivative,
                        t0::Float64, t1::Float64, nsample, curve_entity,
                        begin_entity, end_entity, anisotropic_metric::Bool,
                        integration_precision::Real, max_integration_points,
                        min_integration_depth, max_integration_depth,
                        smooth_ratio, max_smoothing_iterations,
                        caller::AbstractString)
    precision, maxpoints, mindepth, maxdepth = _integration_args(
        integration_precision, max_integration_points, min_integration_depth,
        max_integration_depth, caller)
    evaluate(t, at_begin, at_end) = _metric_point(
        γ, derivative, field, t, t0, t1, curve_entity, begin_entity,
        end_entity, at_begin, at_end, anisotropic_metric, caller)
    points = _integration_points(evaluate, t0, t1, nsample, precision,
                                 maxpoints, mindepth, maxdepth, caller)
    if !anisotropic_metric && smooth_ratio !== nothing
        _smooth_primitive!(points, smooth_ratio, max_smoothing_iterations, caller)
    end
    return points
end

"""
    curve_length(γ; t0=0, t1=1, nsample=nothing, derivative=nothing,
                 integration_precision=1e-9,
                 max_integration_points=1_000_000) -> Float64

Euclidean arc length of `γ`. By default the speed is integrated with Gmsh's
adaptive trapezoidal rule (minimum depth 7, maximum depth 26). Passing an integer
`nsample` uses that many uniform trapezoids, preserving the explicit fixed-grid
mode of earlier Tessella releases. Supply `derivative(t)` when exact tangents are
available; otherwise a bounded second-order finite difference is used.
"""
function curve_length(γ; t0::Real=0.0, t1::Real=1.0, nsample=nothing,
                      derivative=nothing,
                      integration_precision::Real=_GMSH_INTEGRATION_PRECISION,
                      max_integration_points=_DEFAULT_MAX_INTEGRATION_POINTS,
                      min_integration_depth=_GMSH_MIN_INTEGRATION_DEPTH,
                      max_integration_depth=_GMSH_MAX_INTEGRATION_DEPTH)
    t0, t1, nsample = _curve_args(t0, t1, nsample, "curve_length")
    precision, maxpoints, mindepth, maxdepth = _integration_args(
        integration_precision, max_integration_points, min_integration_depth,
        max_integration_depth, "curve_length")
    if nsample !== nothing
        nsample <= typemax(Int) - 1 || throw(ArgumentError(
            "curve_length: nsample exceeds the platform Int limit"))
        nsample + 1 <= maxpoints || throw(ArgumentError(
            "curve_length: nsample=$nsample exceeds the max_integration_points=$maxpoints work bound"))
        # Retain the exact polyline meaning of the public fixed-grid mode.
        previous = _pt3(γ(t0), "curve_length curve")
        total = 0.0
        @inbounds for i in 1:nsample
            t = t0 + (t1 - t0) * i / nsample
            point = _pt3(γ(t), "curve_length curve")
            total += _checked_distance(previous, point, "curve_length")
            isfinite(total) || throw(ArgumentError(
                "curve_length: accumulated length overflowed Float64"))
            previous = point
        end
        return total
    end
    evaluate(t, _, _) = _length_point(
        γ, derivative, t, t0, t1, "curve_length")
    points = _adaptive_points(evaluate, t0, t1, precision, maxpoints,
                              mindepth, maxdepth, "curve_length")
    return points[end].p
end

"""
    metric_length(γ, field; t0=0, t1=1, nsample=nothing, entity=nothing,
                  endpoint_entities=nothing, anisotropic_metric=false,
                  derivative=nothing, smooth_ratio=1.8) -> Float64

Return Gmsh's integrated edge-count primitive. The scalar policy computes
`|γ′|/h` and applies `smoothPrimitive` with `smooth_ratio=1.8`. Set
`smooth_ratio=nothing` to inspect the unsmoothed integral. With
`anisotropic_metric=true`, the integrand is `√(γ′ᵀMγ′)` and smoothing is not
applied, matching Gmsh's BAMG (`Mesh.Algorithm=7`) edge policy.

`entity` classifies the curve. `endpoint_entities=(begin,end)` accepts point
contexts `(0,tag)` or `nothing`. Scalar endpoint sizes are intersected with the
curve size; anisotropic endpoints replace the curve metric, as in Gmsh 4.15.2.
"""
function metric_length(γ, field::AbstractSizeField; t0::Real=0.0,
                       t1::Real=1.0, nsample=nothing, entity=nothing,
                       endpoint_entities=nothing,
                       anisotropic_metric::Bool=false, derivative=nothing,
                       integration_precision::Real=_GMSH_INTEGRATION_PRECISION,
                       max_integration_points=_DEFAULT_MAX_INTEGRATION_POINTS,
                       min_integration_depth=_GMSH_MIN_INTEGRATION_DEPTH,
                       max_integration_depth=_GMSH_MAX_INTEGRATION_DEPTH,
                       smooth_ratio=_GMSH_SMOOTH_RATIO,
                       max_smoothing_iterations=_GMSH_SMOOTH_ITERATIONS)
    t0, t1, nsample = _curve_args(t0, t1, nsample, "metric_length")
    begin_entity, end_entity =
        _endpoint_vertex_entities(endpoint_entities, "metric_length")
    points = _metric_points(
        γ, field, derivative, t0, t1, nsample, entity, begin_entity,
        end_entity, anisotropic_metric, integration_precision,
        max_integration_points, min_integration_depth,
        max_integration_depth, smooth_ratio, max_smoothing_iterations,
        "metric_length")
    return points[end].p
end

function _gmsh_edge_count(metric_total::Float64, minimum_segments,
                          max_edges, closed::Bool, caller::AbstractString)
    minimum = minimum_segments === nothing ? (closed ? 3 : 1) :
              _positive_int(minimum_segments, caller, "minimum_segments")
    maximum = _positive_int(max_edges, caller, "max_edges")
    minimum <= maximum || throw(ArgumentError(
        "$caller: minimum_segments=$minimum exceeds max_edges=$maximum"))
    metric_total >= 0 || throw(ArgumentError(
        "$caller: integrated metric length is negative"))
    count_value = metric_total + 1.99
    (isfinite(count_value) &&
     count_value <= prevfloat(Float64(typemax(Int)))) ||
        throw(ArgumentError(
            "$caller: requested edge count exceeds the platform Int limit; use a larger size field"))
    # meshGEdge uses N = int(a + 1.99), with N counting both endpoints.
    candidate = max(0, trunc(Int, count_value) - 1)
    nedge = max(minimum, candidate)
    nedge <= maximum || throw(ArgumentError(
        "$caller: requested $nedge edges exceeds max_edges=$maximum"))
    return nedge, minimum
end

@inline function _invert_primitive(points::Vector{_IntegrationPoint},
                                   target::Float64)
    target <= points[1].p && return points[1].t
    target >= points[end].p && return points[end].t
    lo = 1
    hi = length(points)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        points[mid].p <= target ? (lo = mid) : (hi = mid)
    end
    p0 = points[lo].p
    p1 = points[hi].p
    weight = p1 == p0 ? 0.0 : (target - p0) / (p1 - p0)
    return points[lo].t + weight * (points[hi].t - points[lo].t)
end

function _check_closed_curve(γ, t0::Float64, t1::Float64)
    first_point = _pt3(γ(t0), "mesh_curve curve")
    last_point = _pt3(γ(t1), "mesh_curve curve")
    closure = _checked_distance(first_point, last_point, "mesh_curve")
    scale = max(maximum(abs, first_point), maximum(abs, last_point), 1.0)
    closure <= 64eps(Float64) * scale || throw(ArgumentError(
        "mesh_curve: closed=true requires γ(t0) == γ(t1) within floating-point tolerance (gap $closure)"))
    return nothing
end

# Port of meshGEdge.cpp filterPoints for the standard (non-BAMG) algorithm.
# Candidate interior vertices closer than 0.3 times the unsmoothed curve size
# are removed as a group, provided the geometry-specific minimum is preserved.
function _filter_close_points!(points::Vector{NTuple{3,Float64}},
                               parameters::Vector{Float64}, γ,
                               field::AbstractSizeField, curve_entity,
                               minimum_segments::Int, closed::Bool,
                               caller::AbstractString)
    length(points) == length(parameters) || throw(ArgumentError(
        "$caller: internal point/parameter length mismatch"))
    last_interior = closed ? length(points) : length(points) - 1
    last_interior >= 2 || return nothing
    interior_count = last_interior - 1
    candidates = Int[]
    sizehint!(candidates, min(interior_count, 256))
    previous_retained = 1
    @inbounds for index in 2:last_interior
        distance = _checked_distance(points[index], points[previous_retained], caller)
        parameter = parameters[index]
        if index != 2
            parameter = parameter / 2 + parameters[previous_retained] / 2
        end
        sample = _pt3(γ(parameter), caller)
        size = size_at(field, sample[1], sample[2], sample[3], curve_entity)
        if distance < size * 0.3
            push!(candidates, index)
        else
            previous_retained = index
        end
    end

    minimum_interior = minimum_segments - 1
    interior_count - length(candidates) >= minimum_interior || return nothing
    isempty(candidates) && return nothing

    candidate_position = 1
    write_position = 1
    @inbounds for read_position in 2:length(points)
        if candidate_position <= length(candidates) &&
           candidates[candidate_position] == read_position
            candidate_position += 1
            continue
        end
        write_position += 1
        points[write_position] = points[read_position]
        parameters[write_position] = parameters[read_position]
    end
    resize!(points, write_position)
    resize!(parameters, write_position)
    return nothing
end

"""
    mesh_curve(γ, field; t0=0, t1=1, closed=false, nsample=nothing,
               entity=nothing, endpoint_entities=nothing,
               anisotropic_metric=false, derivative=nothing,
               smooth_ratio=1.8, minimum_segments=nothing,
               max_edges=10_000_000) -> (points, parameters)

Mesh `γ` at equal increments of the Gmsh-style edge-count primitive. The initial
number of edges follows Gmsh 4.15.2's `int(a + 1.99) - 1` rule, subject to the
requested minimum. Its standard scalar policy then applies `filterPoints`, which
can remove interior points closer than `0.3h` while preserving that minimum; the
BAMG policy does not filter. The default minimum is one open edge or three closed
edges; set `minimum_segments` to reproduce a geometry kernel's larger
curve-specific minimum. Closed curves omit the duplicate endpoint.

Integration, smoothing, derivative and entity-context keywords have the same
meaning as in [`metric_length`](@ref). `max_integration_points` and `max_edges`
bound the two potentially large allocations and fail explicitly when exceeded.
"""
function mesh_curve(γ, field::AbstractSizeField; t0::Real=0.0,
                    t1::Real=1.0, closed::Bool=false, nsample=nothing,
                    entity=nothing, endpoint_entities=nothing,
                    anisotropic_metric::Bool=false, derivative=nothing,
                    integration_precision::Real=_GMSH_INTEGRATION_PRECISION,
                    max_integration_points=_DEFAULT_MAX_INTEGRATION_POINTS,
                    min_integration_depth=_GMSH_MIN_INTEGRATION_DEPTH,
                    max_integration_depth=_GMSH_MAX_INTEGRATION_DEPTH,
                    smooth_ratio=_GMSH_SMOOTH_RATIO,
                    max_smoothing_iterations=_GMSH_SMOOTH_ITERATIONS,
                    minimum_segments=nothing,
                    max_edges=_DEFAULT_MAX_EDGES)
    t0, t1, nsample = _curve_args(t0, t1, nsample, "mesh_curve")
    closed && _check_closed_curve(γ, t0, t1)
    begin_entity, end_entity =
        _endpoint_vertex_entities(endpoint_entities, "mesh_curve")
    integration = _metric_points(
        γ, field, derivative, t0, t1, nsample, entity, begin_entity,
        end_entity, anisotropic_metric, integration_precision,
        max_integration_points, min_integration_depth,
        max_integration_depth, smooth_ratio, max_smoothing_iterations,
        "mesh_curve")
    total = integration[end].p
    total > 0 || throw(ArgumentError(
        "mesh_curve: curve has zero sampled metric length"))
    nedge, effective_minimum = _gmsh_edge_count(
        total, minimum_segments, max_edges, closed, "mesh_curve")
    !closed && nedge == typemax(Int) && throw(ArgumentError(
        "mesh_curve: open-curve node count exceeds the platform Int limit"))
    nnode = closed ? nedge : nedge + 1
    points = Vector{NTuple{3,Float64}}(undef, nnode)
    parameters = Vector{Float64}(undef, nnode)
    @inbounds for k in 0:nnode - 1
        target = total * k / nedge
        parameter = _invert_primitive(integration, target)
        parameters[k + 1] = parameter
        points[k + 1] = _pt3(γ(parameter), "mesh_curve curve")
    end
    if !anisotropic_metric
        _filter_close_points!(points, parameters, γ, field, entity,
                              effective_minimum, closed, "mesh_curve")
    end
    return points, parameters
end

"""
    mesh_segment(a, b, field; kwargs...) -> (points, parameters)

Mesh the straight segment `a → b`. Its exact constant derivative is used, and
all remaining keywords are forwarded to [`mesh_curve`](@ref).
"""
function mesh_segment(a, b, field::AbstractSizeField; kwargs...)
    first_point = _pt3(a, "mesh_segment begin point")
    last_point = _pt3(b, "mesh_segment end point")
    tangent = (last_point[1] - first_point[1],
               last_point[2] - first_point[2],
               last_point[3] - first_point[3])
    all(isfinite, tangent) || throw(ArgumentError(
        "mesh_segment: segment displacement is not finite"))
    γ(t) = (first_point[1] + t * tangent[1],
            first_point[2] + t * tangent[2],
            first_point[3] + t * tangent[3])
    derivative(_) = tangent
    return mesh_curve(γ, field; t0=0.0, t1=1.0, derivative=derivative,
                      kwargs...)
end

end # module Mesh1D
