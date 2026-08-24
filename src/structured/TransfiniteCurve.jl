"""
    TransfiniteCurve

Validated normalized straight-curve node parameters for Gmsh 4.15.2's
`Progression`, `Bump`, and `Beta` transfinite laws and their wall-height
(`HWall`) variants.

This module implements the closed-form distributions on an affine line. It does
not claim bit-for-bit reproduction of Gmsh's adaptive trapezoidal integration,
nor parity for curved CAD parameterizations, `FlexibleTransfinite`, size-map,
periodic, or boundary-layer curves.
"""
module TransfiniteCurve

export transfinite_curve_parameters, transfinite_curve_hwall

const _DEFAULT_MAX_NODES = 10_000_000
const _INT32_MAX = Int(typemax(Int32))
const _LOG_TWO = log(2.0)

function _limit(value, name::AbstractString,
                caller::AbstractString="transfinite_curve_parameters")
    value isa Integer || throw(ArgumentError(
        "$caller: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: $name must not be Bool"))
    value >= 0 || throw(ArgumentError(
        "$caller: $name must be non-negative"))
    value <= _INT32_MAX || throw(ArgumentError(
        "$caller: $name exceeds the Int32 node limit"))
    return Int(value)
end

function _node_count(value, maximum::Int,
                     caller::AbstractString="transfinite_curve_parameters")
    value isa Integer || throw(ArgumentError(
        "$caller: num_nodes must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: num_nodes must not be Bool"))
    value >= 2 || throw(ArgumentError(
        "$caller: num_nodes must be at least 2"))
    value <= _INT32_MAX || throw(ArgumentError(
        "$caller: num_nodes exceeds the Int32 node limit"))
    value <= maximum || throw(ArgumentError(
        "$caller: requested $value nodes exceeds max_nodes=$maximum"))
    return Int(value)
end

function _law(value)
    value isa Symbol || throw(ArgumentError(
        "transfinite_curve_parameters: mesh_type must be a Symbol"))
    value in (:progression, :power, :bump, :beta) || throw(ArgumentError(
        "transfinite_curve_parameters: mesh_type must be :progression, :power, " *
        ":bump, or :beta"))
    return value === :power ? :progression : value
end

function _coefficient(value)
    value isa Real || throw(ArgumentError(
        "transfinite_curve_parameters: coefficient must be real"))
    value isa Bool && throw(ArgumentError(
        "transfinite_curve_parameters: coefficient must not be Bool"))
    coefficient = try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "transfinite_curve_parameters: coefficient must be Float64-representable: " *
            sprint(showerror, err)))
    end
    isfinite(coefficient) || throw(ArgumentError(
        "transfinite_curve_parameters: coefficient must be finite"))
    coefficient != 0.0 || throw(ArgumentError(
        "transfinite_curve_parameters: coefficient must be nonzero"))
    return coefficient
end

@inline function _progression_parameter(k::Int, segments::Int,
                                         log_ratio::Float64)
    if log_ratio > 0.0
        # expm1(k*l)/expm1(n*l), scaled before either exponential can overflow.
        return exp((k - segments) * log_ratio) *
               (-expm1(-k * log_ratio)) / (-expm1(-segments * log_ratio))
    end
    return expm1(k * log_ratio) / expm1(segments * log_ratio)
end

function _fill_progression!(parameters::Vector{Float64}, coefficient::Float64)
    magnitude = abs(coefficient)
    ratio = coefficient > 0.0 ? magnitude : inv(magnitude)
    return _fill_progression_log!(parameters, log(ratio))
end

function _fill_progression_log!(parameters::Vector{Float64},
                                log_ratio::Float64)
    segments = length(parameters) - 1
    if log_ratio == 0.0
        @inbounds for k in 1:segments-1
            parameters[k + 1] = k / segments
        end
        return parameters
    end
    @inbounds for k in 1:segments-1
        parameters[k + 1] = _progression_parameter(k, segments, log_ratio)
    end
    return parameters
end

@inline function _log_sinh_positive(value::Float64)
    value < 20.0 && return log(sinh(value))
    return value - _LOG_TWO + log1p(-exp(-2.0value))
end

@inline function _log_cosh(value::Float64)
    magnitude = abs(value)
    return magnitude - _LOG_TWO + log1p(exp(-2.0magnitude))
end

@inline function _bump_lower_parameter(k::Int, segments::Int,
                                        coefficient::Float64)
    u = k / segments
    if coefficient > 1.0
        q = sqrt(coefficient - 1.0)
        amplitude = atan(q)
        return 0.5 + tan((2.0u - 1.0) * amplitude) / (2.0q)
    end
    q = sqrt(1.0 - coefficient)
    # atanh(q) written this way remains finite when 1 - coefficient rounds to 1.
    amplitude = log1p(q) - 0.5log(coefficient)
    log_parameter = _log_sinh_positive(2.0u * amplitude) - _LOG_TWO -
                    _log_sinh_positive(amplitude) -
                    _log_cosh((1.0 - 2.0u) * amplitude)
    return exp(log_parameter)
end

function _fill_bump!(parameters::Vector{Float64}, coefficient::Float64)
    segments = length(parameters) - 1
    pairs = fld(segments, 2)
    @inbounds for k in 1:pairs
        if 2k == segments
            parameters[k + 1] = 0.5
        else
            value = _bump_lower_parameter(k, segments, coefficient)
            parameters[k + 1] = value
            parameters[segments - k + 1] = 1.0 - value
        end
    end
    return parameters
end

@inline function _beta_parameter(k::Int, segments::Int, coefficient::Float64,
                                 amplitude::Float64,
                                 cosh_amplitude::Float64)
    u = k / segments
    # beta * (tanh(A) - tanh((1-u)A)), with the subtraction rewritten
    # to retain the small parameters produced when beta is close to one.
    return coefficient * sinh(u * amplitude) /
           (cosh_amplitude * cosh((1.0 - u) * amplitude))
end

function _fill_beta!(parameters::Vector{Float64}, coefficient::Float64)
    segments = length(parameters) - 1
    magnitude = abs(coefficient)
    amplitude = 0.5log1p(2.0 / (magnitude - 1.0))
    cosh_amplitude = cosh(amplitude)
    @inbounds for k in 1:segments-1
        parameters[k + 1] = _beta_parameter(
            k, segments, magnitude, amplitude, cosh_amplitude)
    end
    coefficient < 0.0 && _reverse_orientation!(parameters)
    return parameters
end

function _reverse_orientation!(parameters::Vector{Float64})
    count = length(parameters)
    @inbounds for left in 1:fld(count, 2)
        right = count + 1 - left
        left_value = parameters[left]
        right_value = parameters[right]
        parameters[left] = 1.0 - right_value
        parameters[right] = 1.0 - left_value
    end
    if isodd(count)
        middle = fld(count, 2) + 1
        @inbounds parameters[middle] = 1.0 - parameters[middle]
    end
    return parameters
end

function _postcondition(parameters::Vector{Float64}, mesh_type::Symbol,
                        coefficient::Float64,
                        caller::AbstractString="transfinite_curve_parameters")
    previous = parameters[1]
    previous == 0.0 || throw(ErrorException(
        "$caller: internal first-endpoint invariant failed"))
    @inbounds for index in 2:length(parameters)
        current = parameters[index]
        if !(isfinite(current) && previous < current <= 1.0)
            throw(ArgumentError(
                "$caller: $mesh_type coefficient $coefficient " *
                "with $(length(parameters)) nodes is not strictly representable in Float64 " *
                "(node $index is $current after $previous)"))
        end
        previous = current
    end
    previous == 1.0 || throw(ErrorException(
        "$caller: internal last-endpoint invariant failed"))
    return parameters
end

"""
    transfinite_curve_parameters(num_nodes;
        mesh_type=:progression, coefficient=1.0, max_nodes=10_000_000)
        -> Vector{Float64}

Return the `num_nodes` monotonically increasing parameters in `[0, 1]` for a
straight Gmsh-style transfinite curve, including both endpoints.

Supported `mesh_type` values are `:progression` (`:power` is its Gmsh alias),
`:bump`, and `:beta`. The segment ratio for `:progression` is
`abs(coefficient)`; a negative coefficient reverses that distribution. `:bump`
uses the coefficient magnitude and is symmetric, so its sign has no effect.
For `:beta`, magnitudes at most one produce the same uniform fallback as Gmsh;
magnitudes greater than one use the beta law, and a negative coefficient reverses
it. All laws are uniform when `abs(coefficient) == 1`.

`num_nodes` must be at least two and no larger than both `max_nodes` and the
Int32 topology limit. Coefficients must be finite and nonzero. A mathematically
valid distribution that collapses adjacent parameters at Float64 precision is
rejected instead of returning coincident nodes.

This function returns normalized affine-line parameters. Applying the result to
an arbitrary CAD curve does not reproduce Gmsh's derivative-weighted adaptive
integration and is outside this API's contract.
"""
function transfinite_curve_parameters(num_nodes;
        mesh_type=:progression, coefficient=1.0,
        max_nodes=_DEFAULT_MAX_NODES)
    maximum = _limit(max_nodes, "max_nodes")
    count = _node_count(num_nodes, maximum)
    law = _law(mesh_type)
    coef = _coefficient(coefficient)

    parameters = Vector{Float64}(undef, count)
    parameters[1] = 0.0
    parameters[end] = 1.0
    count == 2 && return parameters

    magnitude = abs(coef)
    if magnitude == 1.0 || (law === :beta && magnitude < 1.0)
        segments = count - 1
        @inbounds for k in 1:segments-1
            parameters[k + 1] = k / segments
        end
    elseif law === :progression
        _fill_progression!(parameters, coef)
    elseif law === :bump
        _fill_bump!(parameters, magnitude)
    else
        _fill_beta!(parameters, coef)
    end
    return _postcondition(parameters, law, coef)
end

# log(Σ_{k=0}^{m-1} exp(k*s)), rewritten with negative exponentials so
# neither the geometric sum nor its largest term needs to be represented.
@inline function _log_geometric_sum(segments::Int, s::Float64)
    s == 0.0 && return log(Float64(segments))
    return (segments - 1) * s +
           log(-expm1(-segments * s)) - log(-expm1(-s))
end

# Solve log(Σ_{k=0}^{m-1} exp(k*s)) = log_target for s > 0. The upper
# bound follows from Σ exp(k*s) ≥ exp((m-1)*s).
function _hwall_log_ratio(segments::Int, log_target::Float64)
    log_target > log(Float64(segments)) || return 0.0
    lo = 0.0
    hi = log_target / (segments - 1)
    isfinite(hi) && hi > 0.0 || throw(ArgumentError(
        "transfinite_curve_hwall: wall ratio is not representable in Float64"))
    while _log_geometric_sum(segments, hi) < log_target
        next = 2.0hi
        isfinite(next) || throw(ArgumentError(
            "transfinite_curve_hwall: wall ratio is not representable in Float64"))
        hi = next
    end
    for _ in 1:256
        mid = 0.5 * (lo + hi)
        (mid == lo || mid == hi) && break
        if _log_geometric_sum(segments, mid) < log_target
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

@inline function _bump_hwall_primitive(log_inverse_coefficient::Float64,
                                       wall_fraction::Float64,
                                       count::Int)
    log_inverse_coefficient == 0.0 && return count * wall_fraction
    coefficient = exp(-log_inverse_coefficient)
    q = sqrt(-expm1(-log_inverse_coefficient))
    amplitude = log1p(q) + 0.5log_inverse_coefficient
    primitive = if q < 0.5
        atanh(2q * wall_fraction /
              (coefficient + 2q^2 * wall_fraction))
    else
        log_plus = log(coefficient + 2wall_fraction * q * (q + 1.0))
        log_minus = -log_inverse_coefficient +
                    log1p(-2wall_fraction * q / (1.0 + q))
        0.5 * (log_plus - log_minus)
    end
    return count * primitive / (2amplitude)
end

function _bump_hwall_log_inverse_coefficient(wall_fraction::Float64,
                                             count::Int)
    lo = 0.0
    hi = 1.0
    while _bump_hwall_primitive(hi, wall_fraction, count) < 1.0
        next = 2.0hi
        isfinite(next) || throw(ArgumentError(
            "transfinite_curve_hwall: Bump_HWall coefficient is not " *
            "representable in Float64"))
        hi = next
    end
    for _ in 1:256
        mid = 0.5 * (lo + hi)
        (mid == lo || mid == hi) && break
        if _bump_hwall_primitive(mid, wall_fraction, count) < 1.0
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

function _fill_bump_log_inverse!(parameters::Vector{Float64},
                                 log_inverse_coefficient::Float64)
    log_inverse_coefficient == 0.0 &&
        return _fill_progression_log!(parameters, 0.0)
    segments = length(parameters) - 1
    q = sqrt(-expm1(-log_inverse_coefficient))
    amplitude = log1p(q) + 0.5log_inverse_coefficient
    pairs = fld(segments, 2)
    @inbounds for k in 1:pairs
        if 2k == segments
            parameters[k + 1] = 0.5
        else
            u = k / segments
            log_parameter = _log_sinh_positive(2.0u * amplitude) - _LOG_TWO -
                            _log_sinh_positive(amplitude) -
                            _log_cosh((1.0 - 2.0u) * amplitude)
            value = exp(log_parameter)
            parameters[k + 1] = value
            parameters[segments - k + 1] = 1.0 - value
        end
    end
    return parameters
end

@inline function _log_tanh_positive(value::Float64)
    value < 20.0 && return log(tanh(value))
    exponential = exp(-2.0value)
    return log1p(-exponential) - log1p(exponential)
end

@inline function _beta_hwall_equation(amplitude::Float64,
                                      wall_fraction::Float64,
                                      count::Int)
    ratio = (count - 1) / count
    amplitude == 0.0 && return log((1.0 - wall_fraction) / ratio)
    return log1p(-wall_fraction) + _log_tanh_positive(amplitude) -
           _log_tanh_positive(ratio * amplitude)
end

function _beta_hwall_amplitude(wall_fraction::Float64, count::Int)
    lo = 0.0
    hi = 1.0
    while _beta_hwall_equation(hi, wall_fraction, count) > 0.0
        next = 2.0hi
        isfinite(next) || throw(ArgumentError(
            "transfinite_curve_hwall: Beta_HWall coefficient is not " *
            "representable in Float64"))
        hi = next
    end
    for _ in 1:256
        mid = 0.5 * (lo + hi)
        (mid == lo || mid == hi) && break
        if _beta_hwall_equation(mid, wall_fraction, count) > 0.0
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

function _fill_beta_amplitude!(parameters::Vector{Float64},
                               amplitude::Float64)
    amplitude == 0.0 && return _fill_progression_log!(parameters, 0.0)
    segments = length(parameters) - 1
    log_sinh_amplitude = _log_sinh_positive(amplitude)
    @inbounds for k in 1:segments-1
        u = k / segments
        log_parameter = _log_sinh_positive(u * amplitude) -
                        log_sinh_amplitude -
                        _log_cosh((1.0 - u) * amplitude)
        parameters[k + 1] = exp(log_parameter)
    end
    return parameters
end

"""
    transfinite_curve_hwall(num_nodes; mesh_type=:progression,
                            wall_height, curve_length, orientation=:start,
                            max_nodes=10_000_000)
        -> Vector{Float64}

Return the `num_nodes` monotonically increasing parameters in `[0, 1]` of a
straight `Progression_HWall`, `Bump_HWall`, or `Beta_HWall` transfinite curve.
Use `mesh_type=:progression` (`:power` is an alias), `:bump`, or `:beta`.

For `:progression`, the wall segment measures `wall_height/curve_length` of the
affine curve. The progression ratio `r ≥ 1` is solved from
`wall_height·Σ_{k=0}^{m-1} r^k == curve_length` with `m = num_nodes - 1`
segments, then applied through the same closed-form `Progression` path as
[`transfinite_curve_parameters`](@ref).

For `:bump` and `:beta`, Tessella reproduces Gmsh 4.15.2's nominal HWall
definition: `wall_height/curve_length` locates the first unit of Gmsh's integrated
transfinite primitive, whose inversion uses `num_nodes` (including endpoints).
Consequently the emitted first segment is generally not equal to the nominal wall
height. These laws require `num_nodes ≥ 3` and
`wall_height/curve_length < 1/num_nodes`. `:bump` is symmetric, so orientation
does not change it; `:progression` and `:beta` are mirrored by `orientation=:end`.

All inverse identities are checked to Float64 roundoff, and distributions that
collapse adjacent parameters are rejected. Tessella solves the full representable
root instead of reproducing Gmsh's bounded-search uniform fallbacks.
"""
function transfinite_curve_hwall(num_nodes;
                                 mesh_type=:progression,
                                 wall_height=nothing, curve_length=nothing,
                                 orientation=:start,
                                 max_nodes=_DEFAULT_MAX_NODES)
    caller = "transfinite_curve_hwall"
    mesh_type isa Symbol || throw(ArgumentError(
        "$caller: mesh_type must be a Symbol"))
    mesh_type in (:progression, :power, :bump, :beta) || throw(ArgumentError(
        "$caller: mesh_type must be :progression, :power, :bump, or :beta"))
    law = mesh_type === :power ? :progression : mesh_type
    orientation isa Symbol || throw(ArgumentError(
        "$caller: orientation must be :start or :end"))
    orientation in (:start, :end) || throw(ArgumentError(
        "$caller: orientation must be :start or :end"))
    maximum = _limit(max_nodes, "max_nodes", caller)
    count = _node_count(num_nodes, maximum, caller)

    function positive_float(value, name)
        value === nothing && throw(ArgumentError("$caller: $name is required"))
        value isa Real || throw(ArgumentError("$caller: $name must be real"))
        value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
        converted = try
            Float64(value)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: $name must be Float64-representable: " *
                sprint(showerror, err)))
        end
        isfinite(converted) || throw(ArgumentError(
            "$caller: $name must be finite"))
        converted > 0.0 || throw(ArgumentError(
            "$caller: $name must be positive"))
        return converted
    end
    hw = positive_float(wall_height, "wall_height")
    len = positive_float(curve_length, "curve_length")

    segments = count - 1
    expected_wall = hw / len
    expected_wall > 0.0 || throw(ArgumentError(
        "$caller: wall_height/curve_length is below Float64 resolution"))
    solver_parameter = if law === :progression
        if segments == 1 && hw != len
            throw(ArgumentError(
                "$caller: wall_height must equal curve_length when num_nodes is 2"))
        end
        uniform_wall = len / segments
        uniform_wall > 0.0 || throw(ArgumentError(
            "$caller: curve_length $len divided over $segments segments is below " *
            "Float64 resolution"))
        hw <= uniform_wall || throw(ArgumentError(
            "$caller: wall_height $hw is infeasible for curve_length $len with " *
            "$segments segments (needs wall_height ≤ $uniform_wall)"))
        hw == uniform_wall ?
            0.0 : _hwall_log_ratio(segments, -log(expected_wall))
    else
        count >= 3 || throw(ArgumentError(
            "$caller: $law HWall requires num_nodes to be at least 3"))
        nominal_limit = inv(Float64(count))
        expected_wall < nominal_limit || throw(ArgumentError(
            "$caller: $law HWall requires wall_height/curve_length " *
            "< 1/num_nodes ($nominal_limit)"))
        if law === :bump
            value = _bump_hwall_log_inverse_coefficient(expected_wall, count)
            residual = _bump_hwall_primitive(value, expected_wall, count) - 1.0
            abs(residual) <= 4096eps(Float64) || throw(ErrorException(
                "$caller: Bump_HWall inverse residual is $residual"))
            value
        else
            value = _beta_hwall_amplitude(expected_wall, count)
            residual = _beta_hwall_equation(value, expected_wall, count)
            abs(residual) <= 4096eps(Float64) || throw(ErrorException(
                "$caller: Beta_HWall inverse residual is $residual"))
            value
        end
    end

    parameters = Vector{Float64}(undef, count)
    parameters[1] = 0.0
    parameters[end] = 1.0
    if law === :progression
        _fill_progression_log!(parameters, solver_parameter)
    elseif law === :bump
        _fill_bump_log_inverse!(parameters, solver_parameter)
    else
        _fill_beta_amplitude!(parameters, solver_parameter)
    end
    orientation === :end && law !== :bump && _reverse_orientation!(parameters)
    hwall_law = law === :progression ? :progression_hwall :
                law === :bump ? :bump_hwall : :beta_hwall
    _postcondition(parameters, hwall_law, expected_wall, caller)

    if law === :progression
        wall_parameter = orientation === :start ?
            parameters[2] - parameters[1] : parameters[end] - parameters[end-1]
        tolerance = max(256eps(Float64) * expected_wall, 8eps(expected_wall))
        abs(wall_parameter - expected_wall) <= tolerance ||
            throw(ErrorException(
                "$caller: solved distribution's wall segment measures " *
                "$wall_parameter but the requested wall height implies $expected_wall"))
    end

    return parameters
end

end # module TransfiniteCurve
