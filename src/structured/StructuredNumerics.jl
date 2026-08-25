"""
    StructuredNumerics

Shared numerical certificates for affine structured-volume generators. The
module is internal to Tessella; public entry points live in the individual
transfinite modules.
"""
module StructuredNumerics

const _Dyadic = Rational{BigInt}
# The volume audit tolerates 2^16 output ULPs. Switch to exact interpolation
# while a logical cell can be within 2^40 coordinate ULPs so that three edge
# roundoffs retain headroom below that relative-volume budget.
const _EXACT_CELL_ULPS = 2.0^40
const _VOLUME_ULPS = 65_536
# A determinant whose cancellation consumes at least 12 binary digits is
# recomputed from the exact dyadic values represented by its Float64 inputs.
const _VOLUME_FAST_CONDITION = 2.0^-12

struct _AffineBasis3
    origin::NTuple{3,Float64}
    axis1::NTuple{3,Float64}
    axis2::NTuple{3,Float64}
    axis3::NTuple{3,Float64}
    exact_origin::Union{Nothing,NTuple{3,_Dyadic}}
    exact_axis1::Union{Nothing,NTuple{3,_Dyadic}}
    exact_axis2::Union{Nothing,NTuple{3,_Dyadic}}
    exact_axis3::Union{Nothing,NTuple{3,_Dyadic}}
end

@inline _uses_exact_affine(basis::_AffineBasis3) = basis.exact_origin !== nothing

function _needs_exact_affine(origin::NTuple{3,Float64},
                             axis1::NTuple{3,Float64},
                             axis2::NTuple{3,Float64},
                             axis3::NTuple{3,Float64}, counts::NTuple{3,Int})
    maximum_count = max(counts...)
    @inbounds for dimension in 1:3
        first = axis1[dimension] - origin[dimension]
        second = axis2[dimension] - origin[dimension]
        third = axis3[dimension] - origin[dimension]
        span = abs(first) + abs(second) + abs(third)
        spacing = max(eps(origin[dimension]), eps(axis1[dimension]),
                      eps(axis2[dimension]), eps(axis3[dimension]))
        if !isfinite(span) || (span > 0 &&
           span <= _EXACT_CELL_ULPS * maximum_count * spacing)
            return true
        end
    end
    return false
end

function _affine_basis3(origin::NTuple{3,Float64},
                        axis1::NTuple{3,Float64},
                        axis2::NTuple{3,Float64},
                        axis3::NTuple{3,Float64}, exact::Bool)
    if exact
        ro = ntuple(d -> _Dyadic(origin[d]), 3)
        r1 = ntuple(d -> _Dyadic(axis1[d]) - ro[d], 3)
        r2 = ntuple(d -> _Dyadic(axis2[d]) - ro[d], 3)
        r3 = ntuple(d -> _Dyadic(axis3[d]) - ro[d], 3)
        return _AffineBasis3(origin, axis1, axis2, axis3, ro, r1, r2, r3)
    end
    return _AffineBasis3(origin, axis1, axis2, axis3, nothing, nothing,
                         nothing, nothing)
end

@inline function _affine_point3(basis::_AffineBasis3, first::Float64,
                                second::Float64, third::Float64,
                                caller::AbstractString, logical)
    if !_uses_exact_affine(basis)
        return ntuple(3) do dimension
            d1 = basis.axis1[dimension] - basis.origin[dimension]
            d2 = basis.axis2[dimension] - basis.origin[dimension]
            d3 = basis.axis3[dimension] - basis.origin[dimension]
            muladd(third, d3,
                   muladd(second, d2,
                          muladd(first, d1, basis.origin[dimension])))
        end
    end
    ro = basis.exact_origin::NTuple{3,_Dyadic}
    r1 = basis.exact_axis1::NTuple{3,_Dyadic}
    r2 = basis.exact_axis2::NTuple{3,_Dyadic}
    r3 = basis.exact_axis3::NTuple{3,_Dyadic}
    a = _Dyadic(first); b = _Dyadic(second); c = _Dyadic(third)
    point = ntuple(3) do dimension
        value = Float64(ro[dimension] + a * r1[dimension] +
                        b * r2[dimension] + c * r3[dimension])
        isfinite(value) || throw(ArgumentError(
            "$caller: affine interpolation at logical node $logical is not " *
            "Float64-representable"))
        value
    end
    return point
end

@inline _sub3(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline _maxabs3(a) = max(abs(a[1]), abs(a[2]), abs(a[3]))

function _exact_determinant_measure(a, b, c, d, caller::AbstractString)
    points = ntuple(4) do index
        point = (a, b, c, d)[index]
        ntuple(dimension -> _Dyadic(point[dimension]), 3)
    end
    edge1 = ntuple(dimension -> points[2][dimension] - points[1][dimension], 3)
    edge2 = ntuple(dimension -> points[3][dimension] - points[1][dimension], 3)
    edge3 = ntuple(dimension -> points[4][dimension] - points[1][dimension], 3)
    determinant = abs(
        edge1[1] * (edge2[2] * edge3[3] - edge2[3] * edge3[2]) -
        edge1[2] * (edge2[1] * edge3[3] - edge2[3] * edge3[1]) +
        edge1[3] * (edge2[1] * edge3[2] - edge2[2] * edge3[1]))
    iszero(determinant) && throw(ArgumentError(
        "$caller: exact determinant vanished during volume certification"))
    numerator_value = numerator(determinant)
    denominator_value = denominator(determinant)
    denominator_exponent = trailing_zeros(denominator_value)
    denominator_value == (big(1) << denominator_exponent) ||
        throw(ErrorException(
            "$caller: internal dyadic determinant invariant failed"))
    numerator_bits = ndigits(numerator_value; base=2)
    # At most 53 bits: every retained integer is exactly representable as
    # Float64 and cannot round 2^53-1 upward into the next binade.
    retained_bits = min(numerator_bits, 53)
    shift = numerator_bits - retained_bits
    leading = numerator_value >> shift
    mantissa = ldexp(Float64(leading), -retained_bits)
    (isfinite(mantissa) && 0.5 <= mantissa < 1.0) ||
        throw(ErrorException(
            "$caller: internal exact determinant scaling invariant failed"))
    return mantissa, numerator_bits - denominator_exponent
end

function _determinant_measure(a, b, c, d, caller::AbstractString)
    edge1 = _sub3(b, a)
    edge2 = _sub3(c, a)
    edge3 = _sub3(d, a)
    scale1 = _maxabs3(edge1)
    scale2 = _maxabs3(edge2)
    scale3 = _maxabs3(edge3)
    if !(isfinite(scale1) && isfinite(scale2) && isfinite(scale3) &&
         scale1 > 0 && scale2 > 0 && scale3 > 0)
        return _exact_determinant_measure(a, b, c, d, caller)
    end
    first = (edge1[1] / scale1, edge1[2] / scale1, edge1[3] / scale1)
    second = (edge2[1] / scale2, edge2[2] / scale2, edge2[3] / scale2)
    third = (edge3[1] / scale3, edge3[2] / scale3, edge3[3] / scale3)
    determinant =
        first[1] * (second[2] * third[3] - second[3] * third[2]) -
        first[2] * (second[1] * third[3] - second[3] * third[1]) +
        first[3] * (second[1] * third[2] - second[2] * third[1])
    permanent =
        abs(first[1]) * (abs(second[2] * third[3]) +
                         abs(second[3] * third[2])) +
        abs(first[2]) * (abs(second[1] * third[3]) +
                         abs(second[3] * third[1])) +
        abs(first[3]) * (abs(second[1] * third[2]) +
                         abs(second[2] * third[1]))
    # Independent column scaling avoids exponent under/overflow. The stronger
    # angular-conditioning floor is intentional: this value feeds a relative
    # conservation audit, so a correct sign alone is insufficient after enough
    # cancellation to fabricate or erase material volume.
    if !(isfinite(determinant) && isfinite(permanent) &&
         abs(determinant) > _VOLUME_FAST_CONDITION * permanent)
        return _exact_determinant_measure(a, b, c, d, caller)
    end
    mantissa, exponent = frexp(abs(determinant))
    for scale in (scale1, scale2, scale3)
        scale_mantissa, scale_exponent = frexp(scale)
        mantissa *= scale_mantissa
        exponent += scale_exponent
    end
    adjustment_mantissa, adjustment_exponent = frexp(mantissa)
    return adjustment_mantissa, exponent + adjustment_exponent
end

@inline function _node3(coords, node_id::Int32)
    index = Int(node_id)
    return (coords[1, index], coords[2, index], coords[3, index])
end

function _certify_tet_volume(coords, tets, expected_points::NTuple{4},
                             expected_multiplier::Int,
                             caller::AbstractString,
                             geometry::AbstractString)
    expected_multiplier > 0 || throw(ErrorException(
        "$caller: internal expected-volume multiplier is not positive"))
    expected_mantissa, expected_exponent = _determinant_measure(
        expected_points..., caller)
    expected = expected_multiplier * expected_mantissa
    total = 0.0
    compensation = 0.0
    @inbounds for tet in axes(tets, 2)
        mantissa, exponent = _determinant_measure(
            _node3(coords, tets[1, tet]), _node3(coords, tets[2, tet]),
            _node3(coords, tets[3, tet]), _node3(coords, tets[4, tet]), caller)
        exponent_delta = exponent - expected_exponent
        exponent_delta <= 1023 || throw(ArgumentError(
            "$caller: represented tetrahedron $tet is too large relative " *
            "to the $geometry volume"))
        relative = exponent_delta < -1074 ? 0.0 :
                   ldexp(mantissa, exponent_delta)
        isfinite(relative) || throw(ArgumentError(
            "$caller: represented tetrahedron $tet has a non-finite relative volume"))
        corrected = relative - compensation
        next_total = total + corrected
        compensation = (next_total - total) - corrected
        total = next_total
    end
    isfinite(total) || throw(ArgumentError(
        "$caller: represented tetrahedron volume sum is not finite"))
    tolerance = _VOLUME_ULPS * eps(max(expected, total))
    error = abs(total - expected)
    (isfinite(error) && error <= tolerance) || throw(ArgumentError(
        "$caller: represented tetrahedra do not conserve the $geometry volume " *
        "(sum $total, expected $expected, error $error exceeds tolerance $tolerance)"))
    return nothing
end

end # module StructuredNumerics
