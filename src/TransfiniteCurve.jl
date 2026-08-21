"""
    TransfiniteCurve

Validated normalized straight-curve node parameters for Gmsh 4.15.2's
`Progression`, `Bump`, and `Beta` transfinite laws.

This module implements the closed-form distributions on an affine line. It does
not claim bit-for-bit reproduction of Gmsh's adaptive trapezoidal integration,
nor parity for curved CAD parameterizations, `FlexibleTransfinite`, `*HWall`,
size-map, periodic, or boundary-layer curves.
"""
module TransfiniteCurve

export transfinite_curve_parameters

const _DEFAULT_MAX_NODES = 10_000_000
const _INT32_MAX = Int(typemax(Int32))
const _LOG_TWO = log(2.0)

function _limit(value, name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "transfinite_curve_parameters: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "transfinite_curve_parameters: $name must not be Bool"))
    value >= 0 || throw(ArgumentError(
        "transfinite_curve_parameters: $name must be non-negative"))
    value <= _INT32_MAX || throw(ArgumentError(
        "transfinite_curve_parameters: $name exceeds the Int32 node limit"))
    return Int(value)
end

function _node_count(value, maximum::Int)
    value isa Integer || throw(ArgumentError(
        "transfinite_curve_parameters: num_nodes must be an integer"))
    value isa Bool && throw(ArgumentError(
        "transfinite_curve_parameters: num_nodes must not be Bool"))
    value >= 2 || throw(ArgumentError(
        "transfinite_curve_parameters: num_nodes must be at least 2"))
    value <= _INT32_MAX || throw(ArgumentError(
        "transfinite_curve_parameters: num_nodes exceeds the Int32 node limit"))
    value <= maximum || throw(ArgumentError(
        "transfinite_curve_parameters: requested $value nodes exceeds max_nodes=$maximum"))
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
    segments = length(parameters) - 1
    magnitude = abs(coefficient)
    ratio = coefficient > 0.0 ? magnitude : inv(magnitude)
    log_ratio = log(ratio)
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
                        coefficient::Float64)
    previous = parameters[1]
    previous == 0.0 || throw(ErrorException(
        "transfinite_curve_parameters: internal first-endpoint invariant failed"))
    @inbounds for index in 2:length(parameters)
        current = parameters[index]
        if !(isfinite(current) && previous < current <= 1.0)
            throw(ArgumentError(
                "transfinite_curve_parameters: $mesh_type coefficient $coefficient " *
                "with $(length(parameters)) nodes is not strictly representable in Float64 " *
                "(node $index is $current after $previous)"))
        end
        previous = current
    end
    previous == 1.0 || throw(ErrorException(
        "transfinite_curve_parameters: internal last-endpoint invariant failed"))
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

end # module TransfiniteCurve
