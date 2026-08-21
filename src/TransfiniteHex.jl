"""
    TransfiniteHex

Bounded Gmsh-4.15.2-compatible recombined transfinite meshing for an affine,
six-faced, eight-corner block. The result is an [`Elements.MixedMesh`](@ref)
containing first-order type-3 boundary quadrangles and type-5 hexahedra.

This module deliberately does not claim five-face degeneracies, curved or
warped faces, independently discretized faces, nonuniform curve laws, partial
recombination into prisms or pyramids, QuadTri, holes, multiple blocks,
periodic seams, or high-order elements.
"""
module TransfiniteHex

using ..MeshTypes: tet_volume, triangle_area
using ..Predicates: orient3
import ..Elements
using ..Elements: ElementBlock, MixedMesh

export mesh_transfinite_hex

const _CALLER = "mesh_transfinite_hex"
const _DEFAULT_MAX_NODES = 10_000_000
const _DEFAULT_MAX_HEXAHEDRA = 10_000_000
const _DEFAULT_MAX_BOUNDARY_QUADRANGLES = 20_000_000
const _AFFINE_TOLERANCE = 4096eps(Float64)
const _VOLUME_ULPS = 65_536
const _ARRANGEMENTS = (:left, :right, :alternate_left, :alternate_right)

@inline function _checked_add(a::Int, b::Int, what::AbstractString)
    try
        return Base.checked_add(a, b)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError("$_CALLER: $what count overflows Int"))
    end
end

function _checked_mul(what::AbstractString, values::Int...)
    result = 1
    for value in values
        result = try
            Base.checked_mul(result, value)
        catch err
            err isa InterruptException && rethrow()
            err isa OverflowError || rethrow()
            throw(ArgumentError("$_CALLER: $what count overflows Int"))
        end
    end
    return result
end

function _limit(value, name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$_CALLER: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$_CALLER: $name must not be Bool"))
    value >= 0 || throw(ArgumentError(
        "$_CALLER: $name must be non-negative"))
    value <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $name exceeds the Int32 topology limit"))
    return Int(value)
end

function _count(value, axis::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$_CALLER: $axis cell count must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$_CALLER: $axis cell count must not be Bool"))
    value > 0 || throw(ArgumentError(
        "$_CALLER: $axis cell count must be positive"))
    value <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $axis cell count exceeds the Int32 topology limit"))
    return Int(value)
end

function _three_counts(raw)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: cells must be an indexable collection of three counts"))
    end
    count == 3 || throw(ArgumentError(
        "$_CALLER: cells must contain exactly three counts"))
    values = Vector{Int}(undef, 3)
    cursor = 1
    try
        for value in raw
            cursor <= 3 || throw(ArgumentError(
                "$_CALLER: cells iteration produced more than three counts"))
            values[cursor] = _count(value, ("u", "v", "w")[cursor])
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read cells: $(sprint(showerror, err))"))
    end
    cursor == 4 || throw(ArgumentError(
        "$_CALLER: cells iteration ended before three counts"))
    return values[1], values[2], values[3]
end

function _arrangement(value)
    value isa Symbol || throw(ArgumentError(
        "$_CALLER: arrangement must be a Symbol"))
    value in _ARRANGEMENTS || throw(ArgumentError(
        "$_CALLER: arrangement must be :left, :right, :alternate_left, " *
        "or :alternate_right"))
    return value
end

function _tag(value, name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$_CALLER: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$_CALLER: $name must not be Bool"))
    0 <= value <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $name must lie in 0:$(typemax(Int32))"))
    return Int32(value)
end

function _face_tags(raw)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: face_tags must be an indexable collection"))
    end
    count == 6 || throw(ArgumentError(
        "$_CALLER: face_tags must contain exactly six tags"))
    tags = Vector{Int32}(undef, 6)
    cursor = 1
    try
        for value in raw
            cursor <= 6 || throw(ArgumentError(
                "$_CALLER: face_tags iteration produced more than six tags"))
            tags[cursor] = _tag(value, "face tag $cursor")
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read face_tags: $(sprint(showerror, err))"))
    end
    cursor == 7 || throw(ArgumentError(
        "$_CALLER: face_tags iteration ended before six tags"))
    return tags
end

function _point3(raw, index::Int)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$_CALLER: corner $index is not indexable"))
    end
    count == 3 || throw(ArgumentError(
        "$_CALLER: corner $index must have exactly three coordinates"))
    point = try
        (Float64(raw[1]), Float64(raw[2]), Float64(raw[3]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: corner $index coordinates must be Float64-representable: " *
            sprint(showerror, err)))
    end
    all(isfinite, point) || throw(ArgumentError(
        "$_CALLER: corner $index has a non-finite coordinate"))
    return point
end

function _eight_corners(raw)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: corners must be an indexable collection"))
    end
    count == 8 || throw(ArgumentError(
        "$_CALLER: corners must contain exactly eight points"))
    corners = Vector{NTuple{3,Float64}}(undef, 8)
    cursor = 1
    try
        for value in raw
            cursor <= 8 || throw(ArgumentError(
                "$_CALLER: corners iteration produced more than eight points"))
            corners[cursor] = _point3(value, cursor)
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read corners: $(sprint(showerror, err))"))
    end
    cursor == 9 || throw(ArgumentError(
        "$_CALLER: corners iteration ended before eight points"))
    return corners
end

@inline _sub3(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline _add3(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline _maxabs3(a) = max(abs(a[1]), abs(a[2]), abs(a[3]))
@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                         a[3] * b[1] - a[1] * b[3],
                         a[1] * b[2] - a[2] * b[1])

function _exact_affine_corners(corners)
    @inbounds for dimension in 1:3
        values = ntuple(index -> Rational{BigInt}(corners[index][dimension]), 8)
        origin = values[1]
        values[3] + origin == values[2] + values[4] || return false
        values[6] + origin == values[2] + values[5] || return false
        values[8] + origin == values[4] + values[5] || return false
        values[7] + 2origin == values[2] + values[4] + values[5] || return false
    end
    return true
end

function _certify_affine(corners)
    origin = corners[1]
    scale = 0.0
    deltas = Vector{NTuple{3,Float64}}(undef, 8)
    @inbounds for index in 1:8
        delta = _sub3(corners[index], origin)
        all(isfinite, delta) || throw(ArgumentError(
            "$_CALLER: corner span overflows Float64 at corner $index"))
        deltas[index] = delta
        scale = max(scale, _maxabs3(delta))
    end
    (isfinite(scale) && scale > 0) || throw(ArgumentError(
        "$_CALLER: corners are geometrically degenerate"))
    normalized = NTuple{3,Float64}[
        (delta[1] / scale, delta[2] / scale, delta[3] / scale)
        for delta in deltas]
    all(point -> all(isfinite, point), normalized) || throw(ArgumentError(
        "$_CALLER: normalized corner coordinates are not finite"))

    orientation = orient3(corners[1], corners[2], corners[4], corners[5])
    orientation < 0 || throw(ArgumentError(
        orientation == 0 ?
        "$_CALLER: canonical u/v/w corner directions are coplanar" :
        "$_CALLER: corners must use the positive canonical Gmsh order " *
        "(s0,s1,s2,s3,s4,s5,s6,s7)"))

    u = normalized[2]
    v = normalized[4]
    w = normalized[5]
    determinant = abs(_dot3(u, _cross3(v, w)))
    frobenius_squared = _dot3(u, u) + _dot3(v, v) + _dot3(w, w)
    conditioning_bound = determinant / frobenius_squared
    if !(isfinite(conditioning_bound) && conditioning_bound > 0)
        _exact_affine_corners(corners) && return nothing
        throw(ArgumentError(
            "$_CALLER: affine corner conditioning is not representable in Float64"))
    end
    tolerance = _AFFINE_TOLERANCE * min(1.0, conditioning_bound)
    expected = (_add3(u, v), _add3(u, w),
                _add3(_add3(u, v), w), _add3(v, w))
    indices = (3, 6, 7, 8)
    maximum_error = 0.0
    @inbounds for position in 1:4
        error = _maxabs3(_sub3(normalized[indices[position]], expected[position]))
        isfinite(error) || throw(ArgumentError(
            "$_CALLER: affine residual is not finite"))
        maximum_error = max(maximum_error, error)
    end
    if maximum_error > tolerance && !_exact_affine_corners(corners)
        throw(ArgumentError(
            "$_CALLER: corners do not form an affine parallelepiped " *
            "(normalized residual $maximum_error exceeds the " *
            "conditioning-scaled tolerance $tolerance)"))
    end
    return nothing
end

@inline _lerp(a::Float64, b::Float64, t::Float64) = (1 - t) * a + t * b

@inline function _trilinear(corners, u::Float64, v::Float64, w::Float64)
    ntuple(3) do dimension
        lower0 = _lerp(corners[1][dimension], corners[2][dimension], u)
        lower1 = _lerp(corners[4][dimension], corners[3][dimension], u)
        upper0 = _lerp(corners[5][dimension], corners[6][dimension], u)
        upper1 = _lerp(corners[8][dimension], corners[7][dimension], u)
        _lerp(_lerp(lower0, lower1, v), _lerp(upper0, upper1, v), w)
    end
end

@inline function _node(coords, node_id::Int32)
    index = Int(node_id)
    return (coords[1, index], coords[2, index], coords[3, index])
end

@inline _node_id(i::Int, j::Int, k::Int, npu::Int, npv::Int) =
    Int32(i + 1 + npu * (j + npv * k))

@inline _hex_index(i::Int, j::Int, k::Int, nv::Int, nw::Int) =
    1 + k + nw * (j + nv * i)

function _fill_hexahedra!(hexahedra, nu::Int, nv::Int, nw::Int,
                          npu::Int, npv::Int)
    cursor = 0
    # Gmsh 4.15.2 meshGRegionTransfinite.cpp iterates i, then j, then k.
    @inbounds for i in 0:nu-1, j in 0:nv-1, k in 0:nw-1
        cursor += 1
        hexahedra[1, cursor] = _node_id(i, j, k, npu, npv)
        hexahedra[2, cursor] = _node_id(i + 1, j, k, npu, npv)
        hexahedra[3, cursor] = _node_id(i + 1, j + 1, k, npu, npv)
        hexahedra[4, cursor] = _node_id(i, j + 1, k, npu, npv)
        hexahedra[5, cursor] = _node_id(i, j, k + 1, npu, npv)
        hexahedra[6, cursor] = _node_id(i + 1, j, k + 1, npu, npv)
        hexahedra[7, cursor] = _node_id(i + 1, j + 1, k + 1, npu, npv)
        hexahedra[8, cursor] = _node_id(i, j + 1, k + 1, npu, npv)
    end
    cursor == size(hexahedra, 2) || throw(ErrorException(
        "$_CALLER: internal hexahedron count invariant failed"))
    return nothing
end

@inline function _write_quad!(quadrangles, tags, position::Int,
                              a::Int32, b::Int32, c::Int32, d::Int32,
                              tag::Int32)
    @inbounds begin
        quadrangles[1, position] = a
        quadrangles[2, position] = b
        quadrangles[3, position] = c
        quadrangles[4, position] = d
        tags[position] = tag
    end
    return position + 1
end

function _fill_boundary!(quadrangles, tags, face_tags,
                         nu::Int, nv::Int, nw::Int, npu::Int, npv::Int)
    position = 1
    # Canonical face order: vmin, umax, vmax, umin, wmin, wmax. These are the
    # native surface-element orientations emitted by Gmsh, not an artificial
    # all-outward rewrite. Faces 3, 4 and 5 are consequently reversed relative
    # to the corresponding outward MHexahedron face.
    @inbounds for i in 0:nu-1, k in 0:nw-1
        position = _write_quad!(quadrangles, tags, position,
            _node_id(i, 0, k, npu, npv),
            _node_id(i + 1, 0, k, npu, npv),
            _node_id(i + 1, 0, k + 1, npu, npv),
            _node_id(i, 0, k + 1, npu, npv), face_tags[1])
    end
    @inbounds for j in 0:nv-1, k in 0:nw-1
        position = _write_quad!(quadrangles, tags, position,
            _node_id(nu, j, k, npu, npv),
            _node_id(nu, j + 1, k, npu, npv),
            _node_id(nu, j + 1, k + 1, npu, npv),
            _node_id(nu, j, k + 1, npu, npv), face_tags[2])
    end
    @inbounds for i in 0:nu-1, k in 0:nw-1
        position = _write_quad!(quadrangles, tags, position,
            _node_id(i, nv, k, npu, npv),
            _node_id(i + 1, nv, k, npu, npv),
            _node_id(i + 1, nv, k + 1, npu, npv),
            _node_id(i, nv, k + 1, npu, npv), face_tags[3])
    end
    @inbounds for j in 0:nv-1, k in 0:nw-1
        position = _write_quad!(quadrangles, tags, position,
            _node_id(0, j, k, npu, npv),
            _node_id(0, j + 1, k, npu, npv),
            _node_id(0, j + 1, k + 1, npu, npv),
            _node_id(0, j, k + 1, npu, npv), face_tags[4])
    end
    @inbounds for i in 0:nu-1, j in 0:nv-1
        position = _write_quad!(quadrangles, tags, position,
            _node_id(i, j, 0, npu, npv),
            _node_id(i + 1, j, 0, npu, npv),
            _node_id(i + 1, j + 1, 0, npu, npv),
            _node_id(i, j + 1, 0, npu, npv), face_tags[5])
    end
    @inbounds for i in 0:nu-1, j in 0:nv-1
        position = _write_quad!(quadrangles, tags, position,
            _node_id(i, j, nw, npu, npv),
            _node_id(i + 1, j, nw, npu, npv),
            _node_id(i + 1, j + 1, nw, npu, npv),
            _node_id(i, j + 1, nw, npu, npv), face_tags[6])
    end
    position == size(quadrangles, 2) + 1 || throw(ErrorException(
        "$_CALLER: internal boundary-quadrangle count invariant failed"))
    return nothing
end

# A small outward-rounded interval is used for the common Jacobian path. If
# interval cancellation prevents a proof, the represented Float64 nodes are
# rechecked exactly as dyadic rationals. This keeps the normal path linear and
# allocation-light without turning a floating estimate into a correctness gate.
struct _Interval
    lo::Float64
    hi::Float64
end

struct _ScaledInterval
    bounds::_Interval
    exponent::Int
end

@inline _whole_interval() = _Interval(-Inf, Inf)
@inline _down(value::Float64) = isnan(value) ? -Inf :
    value == -Inf ? -Inf : prevfloat(value)
@inline _up(value::Float64) = isnan(value) ? Inf :
    value == Inf ? Inf : nextfloat(value)

@inline function _iadd(a::_Interval, b::_Interval)
    lo = a.lo + b.lo
    hi = a.hi + b.hi
    (isnan(lo) || isnan(hi)) && return _whole_interval()
    return _Interval(_down(lo), _up(hi))
end

@inline function _isub(a::_Interval, b::_Interval)
    lo = a.lo - b.hi
    hi = a.hi - b.lo
    (isnan(lo) || isnan(hi)) && return _whole_interval()
    return _Interval(_down(lo), _up(hi))
end

@inline function _imul(a::_Interval, b::_Interval)
    p1 = a.lo * b.lo
    p2 = a.lo * b.hi
    p3 = a.hi * b.lo
    p4 = a.hi * b.hi
    any(isnan, (p1, p2, p3, p4)) && return _whole_interval()
    return _Interval(_down(min(p1, p2, p3, p4)),
                     _up(max(p1, p2, p3, p4)))
end

@inline function _idiv_positive(a::_Interval, denominator::Float64)
    lo = a.lo / denominator
    hi = a.hi / denominator
    (isnan(lo) || isnan(hi)) && return _whole_interval()
    return _Interval(_down(lo), _up(hi))
end

@inline function _iadd_positive(a::_Interval, b::_Interval)
    result = _iadd(a, b)
    return _Interval(max(0.0, result.lo), result.hi)
end

@inline function _ildexp_positive(a::_Interval, exponent::Int)
    lo = ldexp(a.lo, exponent)
    hi = ldexp(a.hi, exponent)
    (isnan(lo) || isnan(hi)) && return _Interval(0.0, Inf)
    return _Interval(lo == 0.0 ? 0.0 : max(0.0, _down(lo)), _up(hi))
end

function _scaled_exact_positive(value::Rational{BigInt})
    value > 0 || throw(ErrorException(
        "$_CALLER: internal positive-measure invariant failed"))
    numerator_value = numerator(value)
    denominator_value = denominator(value)
    numerator_bits = ndigits(numerator_value; base=2)
    denominator_bits = ndigits(denominator_value; base=2)
    candidate = numerator_bits - denominator_bits
    at_least_candidate = candidate >= 0 ?
        numerator_value >= (denominator_value << candidate) :
        (numerator_value << -candidate) >= denominator_value
    exponent = candidate + (at_least_candidate ? 1 : 0)

    shift = 53 - exponent
    quotient, remainder = shift >= 0 ?
        divrem(numerator_value << shift, denominator_value) :
        divrem(numerator_value, denominator_value << -shift)
    (big(1) << 52) <= quotient < (big(1) << 53) || throw(ErrorException(
        "$_CALLER: internal exact-measure scaling invariant failed"))
    lo = ldexp(Float64(quotient), -53)
    hi = iszero(remainder) ? lo : ldexp(Float64(quotient + 1), -53)
    return _ScaledInterval(_Interval(lo, hi), exponent)
end

function _normalize_scaled_positive(bounds::_Interval, exponent::Int)
    (bounds.lo > 0 && isfinite(bounds.hi)) || throw(ErrorException(
        "$_CALLER: internal interval-measure invariant failed"))
    _, adjustment = frexp(bounds.hi)
    normalized = _ildexp_positive(bounds, -adjustment)
    (normalized.lo > 0 && isfinite(normalized.hi)) || throw(ErrorException(
        "$_CALLER: internal normalized-measure invariant failed"))
    return _ScaledInterval(normalized, exponent + adjustment)
end

@inline _ivsub(a, b) = (_isub(a[1], b[1]), _isub(a[2], b[2]),
                         _isub(a[3], b[3]))

@inline function _idet3(a, b, c)
    first = _imul(a[1], _isub(_imul(b[2], c[3]), _imul(b[3], c[2])))
    second = _imul(a[2], _isub(_imul(b[1], c[3]), _imul(b[3], c[1])))
    third = _imul(a[3], _isub(_imul(b[1], c[2]), _imul(b[2], c[1])))
    return _iadd(_isub(first, second), third)
end

@inline function _power_to_bernstein_axis!(coefficients, axis::Int)
    half = _Interval(0.5, 0.5)
    if axis == 1
        @inbounds for j in 1:3, k in 1:3
            a0, a1, a2 = coefficients[1, j, k], coefficients[2, j, k],
                         coefficients[3, j, k]
            coefficients[1, j, k] = a0
            coefficients[2, j, k] = _iadd(a0, _imul(half, a1))
            coefficients[3, j, k] = _iadd(_iadd(a0, a1), a2)
        end
    elseif axis == 2
        @inbounds for i in 1:3, k in 1:3
            a0, a1, a2 = coefficients[i, 1, k], coefficients[i, 2, k],
                         coefficients[i, 3, k]
            coefficients[i, 1, k] = a0
            coefficients[i, 2, k] = _iadd(a0, _imul(half, a1))
            coefficients[i, 3, k] = _iadd(_iadd(a0, a1), a2)
        end
    else
        @inbounds for i in 1:3, j in 1:3
            a0, a1, a2 = coefficients[i, j, 1], coefficients[i, j, 2],
                         coefficients[i, j, 3]
            coefficients[i, j, 1] = a0
            coefficients[i, j, 2] = _iadd(a0, _imul(half, a1))
            coefficients[i, j, 3] = _iadd(_iadd(a0, a1), a2)
        end
    end
    return nothing
end

@inline _det3(a, b, c) =
    a[1] * (b[2] * c[3] - b[3] * c[2]) -
    a[2] * (b[1] * c[3] - b[3] * c[1]) +
    a[3] * (b[1] * c[2] - b[2] * c[1])

function _exact_positive_bernstein_jacobian(coords, nodes, cell::Int)
    rational = Rational{BigInt}
    points = ntuple(8) do slot
        index = Int(nodes[slot])
        ntuple(dimension -> rational(coords[dimension, index]), 3)
    end
    origin = points[1]
    deltas = ntuple(slot -> _sub3(points[slot], origin), 8)
    u = deltas[2]
    v = deltas[4]
    w = deltas[5]
    uv = _sub3(_sub3(deltas[3], u), v)
    uw = _sub3(_sub3(deltas[6], u), w)
    vw = _sub3(_sub3(deltas[8], v), w)
    uvw = _sub3(_sub3(_sub3(_sub3(deltas[7], u), v), w),
                 _add3(_add3(uv, uw), vw))
    du = (((0, 0, 0), u), ((0, 1, 0), uv),
          ((0, 0, 1), uw), ((0, 1, 1), uvw))
    dv = (((0, 0, 0), v), ((1, 0, 0), uv),
          ((0, 0, 1), vw), ((1, 0, 1), uvw))
    dw = (((0, 0, 0), w), ((1, 0, 0), uw),
          ((0, 1, 0), vw), ((1, 1, 0), uvw))
    zero_rational = zero(rational)
    coefficients = fill(zero_rational, 3, 3, 3)
    @inbounds for (eu, xu) in du, (ev, xv) in dv, (ew, xw) in dw
        i = eu[1] + ev[1] + ew[1] + 1
        j = eu[2] + ev[2] + ew[2] + 1
        k = eu[3] + ev[3] + ew[3] + 1
        coefficients[i, j, k] += _det3(xu, xv, xw)
    end
    @inbounds for j in 1:3, k in 1:3
        a0, a1, a2 = coefficients[1, j, k], coefficients[2, j, k],
                     coefficients[3, j, k]
        coefficients[2, j, k] = a0 + a1 // 2
        coefficients[3, j, k] = a0 + a1 + a2
    end
    @inbounds for i in 1:3, k in 1:3
        a0, a1, a2 = coefficients[i, 1, k], coefficients[i, 2, k],
                     coefficients[i, 3, k]
        coefficients[i, 2, k] = a0 + a1 // 2
        coefficients[i, 3, k] = a0 + a1 + a2
    end
    @inbounds for i in 1:3, j in 1:3
        a0, a1, a2 = coefficients[i, j, 1], coefficients[i, j, 2],
                     coefficients[i, j, 3]
        coefficients[i, j, 2] = a0 + a1 // 2
        coefficients[i, j, 3] = a0 + a1 + a2
    end
    @inbounds for k in 1:3, j in 1:3, i in 1:3
        coefficients[i, j, k] > 0 || throw(ArgumentError(
            "$_CALLER: hexahedron $cell lacks a global positive trilinear " *
            "Jacobian certificate (Bernstein coefficient ($i,$j,$k) is " *
            "non-positive)"))
    end
    integral = zero_rational
    @inbounds for coefficient in coefficients
        integral += coefficient
    end
    return _scaled_exact_positive(integral / 27)
end

function _certify_positive_jacobian!(coefficients, coords, nodes, cell::Int)
    origin = _node(coords, nodes[1])
    scale = 0.0
    @inbounds for slot in 2:8
        point = _node(coords, nodes[slot])
        delta = _sub3(point, origin)
        all(isfinite, delta) || throw(ArgumentError(
            "$_CALLER: hexahedron $cell coordinate span is not finite"))
        scale = max(scale, _maxabs3(delta))
    end
    (isfinite(scale) && scale > 0) || throw(ArgumentError(
        "$_CALLER: hexahedron $cell has no representable extent"))

    interval_points = ntuple(8) do slot
        point = _node(coords, nodes[slot])
        ntuple(3) do dimension
            numerator = _isub(_Interval(point[dimension], point[dimension]),
                              _Interval(origin[dimension], origin[dimension]))
            _idiv_positive(numerator, scale)
        end
    end
    u = interval_points[2]
    v = interval_points[4]
    w = interval_points[5]
    uv = _ivsub(_ivsub(interval_points[3], u), v)
    uw = _ivsub(_ivsub(interval_points[6], u), w)
    vw = _ivsub(_ivsub(interval_points[8], v), w)
    uvw = _ivsub(
        _ivsub(_ivsub(_ivsub(interval_points[7], u), v), w),
        (_iadd(_iadd(uv[1], uw[1]), vw[1]),
         _iadd(_iadd(uv[2], uw[2]), vw[2]),
         _iadd(_iadd(uv[3], uw[3]), vw[3])))
    du = (((0, 0, 0), u), ((0, 1, 0), uv),
          ((0, 0, 1), uw), ((0, 1, 1), uvw))
    dv = (((0, 0, 0), v), ((1, 0, 0), uv),
          ((0, 0, 1), vw), ((1, 0, 1), uvw))
    dw = (((0, 0, 0), w), ((1, 0, 0), uw),
          ((0, 1, 0), vw), ((1, 1, 0), uvw))
    fill!(coefficients, _Interval(0.0, 0.0))
    @inbounds for (eu, xu) in du, (ev, xv) in dv, (ew, xw) in dw
        i = eu[1] + ev[1] + ew[1] + 1
        j = eu[2] + ev[2] + ew[2] + 1
        k = eu[3] + ev[3] + ew[3] + 1
        coefficients[i, j, k] = _iadd(
            coefficients[i, j, k], _idet3(xu, xv, xw))
    end
    _power_to_bernstein_axis!(coefficients, 1)
    _power_to_bernstein_axis!(coefficients, 2)
    _power_to_bernstein_axis!(coefficients, 3)
    if all(coefficient -> coefficient.lo > 0, coefficients)
        integral = _Interval(0.0, 0.0)
        @inbounds for coefficient in coefficients
            integral = _iadd_positive(integral, coefficient)
        end
        integral = _idiv_positive(integral, 27.0)
        scale_mantissa, scale_exponent = frexp(scale)
        scale_interval = _Interval(scale_mantissa, scale_mantissa)
        physical = _imul(
            _imul(_imul(integral, scale_interval), scale_interval),
            scale_interval)
        if physical.lo > 0 && isfinite(physical.hi)
            return _normalize_scaled_positive(physical, 3scale_exponent)
        end
    end
    return _exact_positive_bernstein_jacobian(coords, nodes, cell)
end

const _CORNER_JACOBIANS = (
    (1, 2, 4, 5, -1), (2, 1, 3, 6, 1),
    (3, 4, 2, 7, -1), (4, 3, 1, 8, 1),
    (5, 6, 8, 1, 1), (6, 5, 7, 2, -1),
    (7, 8, 6, 3, 1), (8, 7, 5, 4, -1))

const _SIX_TET_CERTIFICATE = (
    (1, 2, 4, 5), (2, 4, 5, 6), (5, 6, 4, 8),
    (2, 4, 6, 3), (4, 8, 6, 3), (6, 8, 7, 3))

function _exact_affine_volume(corners)
    rational = Rational{BigInt}
    points = ntuple(4) do position
        corner = corners[(1, 2, 4, 5)[position]]
        ntuple(dimension -> rational(corner[dimension]), 3)
    end
    u = _sub3(points[2], points[1])
    v = _sub3(points[3], points[1])
    w = _sub3(points[4], points[1])
    determinant = _det3(u, v, w)
    determinant > 0 || throw(ErrorException(
        "$_CALLER: internal affine-volume orientation invariant failed"))
    return _scaled_exact_positive(determinant)
end

function _accumulate_scaled!(bins, occupied, measure::_ScaledInterval,
                             reference_exponent::Int, cell::Int)
    relative = _ildexp_positive(
        measure.bounds, measure.exponent - reference_exponent)
    isfinite(relative.hi) || throw(ArgumentError(
        "$_CALLER: represented hexahedron $cell is too large relative to " *
        "the affine volume"))
    level = 1
    while occupied[level]
        relative = _iadd_positive(bins[level], relative)
        occupied[level] = false
        level += 1
        level <= length(bins) || throw(ErrorException(
            "$_CALLER: internal volume-accumulator capacity invariant failed"))
    end
    bins[level] = relative
    occupied[level] = true
    return nothing
end

function _certify_integrated_volume(bins, occupied, expected::_ScaledInterval)
    total = _Interval(0.0, 0.0)
    @inbounds for level in eachindex(bins)
        occupied[level] || continue
        total = _iadd_positive(total, bins[level])
    end
    isfinite(total.hi) || throw(ArgumentError(
        "$_CALLER: represented integrated trilinear-Jacobian sum is not finite"))
    scale = max(expected.bounds.hi, total.hi)
    tolerance = _up(_VOLUME_ULPS * eps(scale))
    error = max(0.0,
                _up(expected.bounds.hi - total.lo),
                _up(total.hi - expected.bounds.lo))
    (isfinite(error) && error <= tolerance) || throw(ArgumentError(
        "$_CALLER: represented hexahedra do not conserve the affine volume " *
        "under exact trilinear-Jacobian integration (sum interval " *
        "[$(total.lo), $(total.hi)], expected interval " *
        "[$(expected.bounds.lo), $(expected.bounds.hi)], error bound $error " *
        "exceeds tolerance $tolerance)"))
    return nothing
end

function _certify_hexahedron_geometry(coords, hexahedra, corners)
    coefficients = Array{_Interval}(undef, 3, 3, 3)
    expected_volume = _exact_affine_volume(corners)
    # Int32 limits cap the carry depth at 31. Extra slots make the invariant
    # explicit and keep accumulation storage independent of the mesh size.
    volume_bins = fill(_Interval(0.0, 0.0), 64)
    occupied_bins = falses(64)
    @inbounds for cell in axes(hexahedra, 2)
        nodes = ntuple(slot -> hexahedra[slot, cell], 8)
        measure = _certify_positive_jacobian!(coefficients, coords, nodes, cell)
        _accumulate_scaled!(volume_bins, occupied_bins, measure,
                            expected_volume.exponent, cell)
        for (corner, u_neighbor, v_neighbor, w_neighbor, expected_sign) in
                _CORNER_JACOBIANS
            a = _node(coords, nodes[corner])
            b = _node(coords, nodes[u_neighbor])
            c = _node(coords, nodes[v_neighbor])
            d = _node(coords, nodes[w_neighbor])
            sign = orient3(a, b, c, d)
            sign == expected_sign || throw(ArgumentError(
                "$_CALLER: hexahedron $cell has a zero or reversed exact " *
                "corner Jacobian at local corner $corner"))
            measure = tet_volume(a, b, c, d)
            (isfinite(measure) && measure > 0) || throw(ArgumentError(
                "$_CALLER: hexahedron $cell corner Jacobian magnitude is " *
                "not a finite positive Float64 at local corner $corner"))
        end
        for slots in _SIX_TET_CERTIFICATE
            a = _node(coords, nodes[slots[1]])
            b = _node(coords, nodes[slots[2]])
            c = _node(coords, nodes[slots[3]])
            d = _node(coords, nodes[slots[4]])
            orient3(a, b, c, d) == -1 || throw(ArgumentError(
                "$_CALLER: hexahedron $cell has a zero or reversed canonical " *
                "six-tetrahedron subcell"))
            measure = tet_volume(a, b, c, d)
            (isfinite(measure) && measure > 0) || throw(ArgumentError(
                "$_CALLER: hexahedron $cell has a canonical subcell whose " *
            "volume is not a finite positive Float64"))
        end
    end
    _certify_integrated_volume(volume_bins, occupied_bins, expected_volume)
    return nothing
end

@inline function _sort4(a::Int32, b::Int32, c::Int32, d::Int32)
    a > b && ((a, b) = (b, a))
    c > d && ((c, d) = (d, c))
    a > c && ((a, c) = (c, a))
    b > d && ((b, d) = (d, b))
    b > c && ((b, c) = (c, b))
    return (a, b, c, d)
end

@inline _quad_key(quadrangles, position::Int) = _sort4(
    quadrangles[1, position], quadrangles[2, position],
    quadrangles[3, position], quadrangles[4, position])

@inline _hex_face_key(hexahedra, cell::Int, slots) = _sort4(
    hexahedra[slots[1], cell], hexahedra[slots[2], cell],
    hexahedra[slots[3], cell], hexahedra[slots[4], cell])

function _certify_quad_geometry(coords, quadrangles, position::Int,
                                interior::Int32, expected_sign::Int)
    a = _node(coords, quadrangles[1, position])
    b = _node(coords, quadrangles[2, position])
    c = _node(coords, quadrangles[3, position])
    d = _node(coords, quadrangles[4, position])
    inside = _node(coords, interior)
    first_sign = orient3(a, b, c, inside)
    second_sign = orient3(a, c, d, inside)
    (first_sign == expected_sign && second_sign == expected_sign) ||
        throw(ArgumentError(
            "$_CALLER: boundary quadrangle $position has a zero or reversed " *
            "native Gmsh orientation"))
    first_area = triangle_area(a, b, c)
    second_area = triangle_area(a, c, d)
    (isfinite(first_area) && first_area > 0 &&
     isfinite(second_area) && second_area > 0) || throw(ArgumentError(
        "$_CALLER: boundary quadrangle $position has a non-finite or zero " *
        "represented area"))
    return nothing
end

function _certify_topology_and_boundary(coords, quadrangles, quad_tags,
                                         hexahedra, face_tags,
                                         nu::Int, nv::Int, nw::Int,
                                         npu::Int, npv::Int)
    # MHexahedron::faces_hexa local slots, converted from zero- to one-based.
    vmin = (1, 2, 6, 5)
    umax = (2, 3, 7, 6)
    vmax = (3, 4, 8, 7)
    umin = (1, 5, 8, 4)
    wmin = (1, 4, 3, 2)
    wmax = (5, 6, 7, 8)

    # Every internal logical face must be represented identically by its two
    # neighboring type-5 cells. This is a constant-memory manifold audit.
    if nu > 1
        @inbounds for i in 0:nu-2, j in 0:nv-1, k in 0:nw-1
            left = _hex_index(i, j, k, nv, nw)
            right = _hex_index(i + 1, j, k, nv, nw)
            _hex_face_key(hexahedra, left, umax) ==
                _hex_face_key(hexahedra, right, umin) || throw(ErrorException(
                    "$_CALLER: internal u-face topology invariant failed"))
        end
    end
    if nv > 1
        @inbounds for i in 0:nu-1, j in 0:nv-2, k in 0:nw-1
            lower = _hex_index(i, j, k, nv, nw)
            upper = _hex_index(i, j + 1, k, nv, nw)
            _hex_face_key(hexahedra, lower, vmax) ==
                _hex_face_key(hexahedra, upper, vmin) || throw(ErrorException(
                    "$_CALLER: internal v-face topology invariant failed"))
        end
    end
    if nw > 1
        @inbounds for i in 0:nu-1, j in 0:nv-1, k in 0:nw-2
            lower = _hex_index(i, j, k, nv, nw)
            upper = _hex_index(i, j, k + 1, nv, nw)
            _hex_face_key(hexahedra, lower, wmax) ==
                _hex_face_key(hexahedra, upper, wmin) || throw(ErrorException(
                    "$_CALLER: internal w-face topology invariant failed"))
        end
    end

    position = 0
    @inbounds for i in 0:nu-1, k in 0:nw-1
        position += 1
        cell = _hex_index(i, 0, k, nv, nw)
        _quad_key(quadrangles, position) ==
            _hex_face_key(hexahedra, cell, vmin) || throw(ErrorException(
                "$_CALLER: vmin boundary topology invariant failed"))
        quad_tags[position] == face_tags[1] || throw(ErrorException(
            "$_CALLER: vmin boundary tag invariant failed"))
        _certify_quad_geometry(coords, quadrangles, position,
            _node_id(i, 1, k, npu, npv), 1)
    end
    @inbounds for j in 0:nv-1, k in 0:nw-1
        position += 1
        cell = _hex_index(nu - 1, j, k, nv, nw)
        _quad_key(quadrangles, position) ==
            _hex_face_key(hexahedra, cell, umax) || throw(ErrorException(
                "$_CALLER: umax boundary topology invariant failed"))
        quad_tags[position] == face_tags[2] || throw(ErrorException(
            "$_CALLER: umax boundary tag invariant failed"))
        _certify_quad_geometry(coords, quadrangles, position,
            _node_id(nu - 1, j, k, npu, npv), 1)
    end
    @inbounds for i in 0:nu-1, k in 0:nw-1
        position += 1
        cell = _hex_index(i, nv - 1, k, nv, nw)
        _quad_key(quadrangles, position) ==
            _hex_face_key(hexahedra, cell, vmax) || throw(ErrorException(
                "$_CALLER: vmax boundary topology invariant failed"))
        quad_tags[position] == face_tags[3] || throw(ErrorException(
            "$_CALLER: vmax boundary tag invariant failed"))
        _certify_quad_geometry(coords, quadrangles, position,
            _node_id(i, nv - 1, k, npu, npv), -1)
    end
    @inbounds for j in 0:nv-1, k in 0:nw-1
        position += 1
        cell = _hex_index(0, j, k, nv, nw)
        _quad_key(quadrangles, position) ==
            _hex_face_key(hexahedra, cell, umin) || throw(ErrorException(
                "$_CALLER: umin boundary topology invariant failed"))
        quad_tags[position] == face_tags[4] || throw(ErrorException(
            "$_CALLER: umin boundary tag invariant failed"))
        _certify_quad_geometry(coords, quadrangles, position,
            _node_id(1, j, k, npu, npv), -1)
    end
    @inbounds for i in 0:nu-1, j in 0:nv-1
        position += 1
        cell = _hex_index(i, j, 0, nv, nw)
        _quad_key(quadrangles, position) ==
            _hex_face_key(hexahedra, cell, wmin) || throw(ErrorException(
                "$_CALLER: wmin boundary topology invariant failed"))
        quad_tags[position] == face_tags[5] || throw(ErrorException(
            "$_CALLER: wmin boundary tag invariant failed"))
        _certify_quad_geometry(coords, quadrangles, position,
            _node_id(i, j, 1, npu, npv), -1)
    end
    @inbounds for i in 0:nu-1, j in 0:nv-1
        position += 1
        cell = _hex_index(i, j, nw - 1, nv, nw)
        _quad_key(quadrangles, position) ==
            _hex_face_key(hexahedra, cell, wmax) || throw(ErrorException(
                "$_CALLER: wmax boundary topology invariant failed"))
        quad_tags[position] == face_tags[6] || throw(ErrorException(
            "$_CALLER: wmax boundary tag invariant failed"))
        _certify_quad_geometry(coords, quadrangles, position,
            _node_id(i, j, nw - 1, npu, npv), 1)
    end
    position == size(quadrangles, 2) || throw(ErrorException(
        "$_CALLER: certified boundary count invariant failed"))
    return nothing
end

"""
    mesh_transfinite_hex(corners, cells=(1,1,1);
        arrangement=:left, volume_tag=0,
        face_tags=(0,0,0,0,0,0), max_nodes=10_000_000,
        max_hexahedra=10_000_000,
        max_boundary_quadrangles=20_000_000) -> MixedMesh

Mesh an affine six-face block with Gmsh 4.15.2's fully recombined transfinite
volume path. `corners` contains eight finite 3-D points in canonical order
`(s0,s1,s2,s3,s4,s5,s6,s7)`: the first four wind around the `w=0` face and
the last four are their `w=1` counterparts. Derived corners must agree with an
affine parallelepiped within a conditioning-scaled `4096eps(Float64)` input
tolerance. `cells=(nu,nv,nw)` contains positive, uniform logical-cell counts.

The output has one Gmsh type-3 boundary block followed by one type-5 volume
block. Hexahedra use Gmsh's exact `CREATE_HEX` local order and its `i,j,k`
element iteration order (`k` varies fastest). Boundary quadrangles use the
native orders of the six transfinite surface entities. Consequently faces 3,
4 and 5 are reversed relative to outward MHexahedron faces; this is deliberate
Gmsh compatibility, not an orientation defect. `face_tags` follow canonical
order `(vmin,umax,vmax,umin,wmin,wmax)`, and `volume_tag` labels every hex.

`:left`, `:right`, `:alternate_left`, and `:alternate_right` are accepted for
`arrangement` and intentionally produce identical output: once all six faces
are recombined, Gmsh emits `(v1,v2,v3,v4)` before consulting the triangle
arrangement. Other arrangements are rejected.

Every represented type-5 cell receives an outward-rounded interval certificate
that all 27 degree-(2,2,2) Bernstein coefficients of its trilinear Jacobian are
strictly positive; uncertain interval cases fall back to exact rational
arithmetic over the represented dyadic nodes. The exact integral of that
polynomial is enclosed for each cell and accumulated in scaled, balanced
intervals. The resulting error enclosure must lie within `65_536` Float64 ulps
of the normalized exact affine-corner volume, even when the physical
determinant is below Float64's range. Exact corner signs, the canonical
six-tetrahedron subcells, finite positive subcell measures and boundary areas,
internal face pairing, and exact boundary equality are checked as independent
postconditions. Counts, Int32 bounds, and caller limits are checked before
corner conversion or output allocation.

Unsupported here: five-face or degenerate volumes, curved/warped or
independently discretized faces, nonuniform curve laws, partial recombination
into prisms or pyramids, QuadTri, holes, multiple blocks, periodic seams,
embedded entities, high-order elements, and coordinate scales whose represented
boundary areas or cell measures are not finite positive Float64 values.
"""
function mesh_transfinite_hex(corners, cells=(1, 1, 1);
                              arrangement=:left,
                              volume_tag=0,
                              face_tags=(0, 0, 0, 0, 0, 0),
                              max_nodes=_DEFAULT_MAX_NODES,
                              max_hexahedra=_DEFAULT_MAX_HEXAHEDRA,
                              max_boundary_quadrangles=
                                  _DEFAULT_MAX_BOUNDARY_QUADRANGLES)::MixedMesh
    _arrangement(arrangement)
    node_limit = _limit(max_nodes, "max_nodes")
    hexahedron_limit = _limit(max_hexahedra, "max_hexahedra")
    quadrangle_limit = _limit(
        max_boundary_quadrangles, "max_boundary_quadrangles")
    nu, nv, nw = _three_counts(cells)

    npu = _checked_add(nu, 1, "node")
    npv = _checked_add(nv, 1, "node")
    npw = _checked_add(nw, 1, "node")
    node_count = _checked_mul("node", npu, npv, npw)
    hexahedron_count = _checked_mul("hexahedron", nu, nv, nw)
    uv = _checked_mul("boundary quadrangle", nu, nv)
    uw = _checked_mul("boundary quadrangle", nu, nw)
    vw = _checked_mul("boundary quadrangle", nv, nw)
    face_cells = _checked_add(
        _checked_add(uv, uw, "boundary quadrangle"),
        vw, "boundary quadrangle")
    quadrangle_count = _checked_mul("boundary quadrangle", 2, face_cells)
    _checked_mul("hexahedron connectivity", 8, hexahedron_count)
    _checked_mul("boundary connectivity", 4, quadrangle_count)

    node_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $node_count nodes exceed Int32 indexing"))
    hexahedron_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $hexahedron_count hexahedra exceed the Int32 topology limit"))
    quadrangle_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $quadrangle_count boundary quadrangles exceed the Int32 " *
        "topology limit"))
    node_count <= node_limit || throw(ArgumentError(
        "$_CALLER: $node_count nodes exceed max_nodes=$node_limit"))
    hexahedron_count <= hexahedron_limit || throw(ArgumentError(
        "$_CALLER: $hexahedron_count hexahedra exceed " *
        "max_hexahedra=$hexahedron_limit"))
    quadrangle_count <= quadrangle_limit || throw(ArgumentError(
        "$_CALLER: $quadrangle_count boundary quadrangles exceed " *
        "max_boundary_quadrangles=$quadrangle_limit"))

    converted_corners = _eight_corners(corners)
    _certify_affine(converted_corners)
    converted_face_tags = _face_tags(face_tags)
    converted_volume_tag = _tag(volume_tag, "volume_tag")

    coords = Matrix{Float64}(undef, 3, node_count)
    @inbounds for k in 0:nw, j in 0:nv, i in 0:nu
        point = _trilinear(
            converted_corners, i / nu, j / nv, k / nw)
        all(isfinite, point) || throw(ArgumentError(
            "$_CALLER: interpolation produced a non-finite coordinate at " *
            "logical node ($i,$j,$k)"))
        index = Int(_node_id(i, j, k, npu, npv))
        coords[1, index] = point[1]
        coords[2, index] = point[2]
        coords[3, index] = point[3]
    end

    hexahedra = Matrix{Int32}(undef, 8, hexahedron_count)
    _fill_hexahedra!(hexahedra, nu, nv, nw, npu, npv)
    quadrangles = Matrix{Int32}(undef, 4, quadrangle_count)
    quadrangle_tags = Vector{Int32}(undef, quadrangle_count)
    _fill_boundary!(quadrangles, quadrangle_tags, converted_face_tags,
                    nu, nv, nw, npu, npv)

    _certify_hexahedron_geometry(coords, hexahedra, converted_corners)
    _certify_topology_and_boundary(
        coords, quadrangles, quadrangle_tags, hexahedra,
        converted_face_tags, nu, nv, nw, npu, npv)

    quadrangle_block = ElementBlock(3, quadrangles, quadrangle_tags)
    hexahedron_block = ElementBlock(
        5, hexahedra, fill(converted_volume_tag, hexahedron_count))
    result = MixedMesh(coords, [quadrangle_block, hexahedron_block])
    diagnostic = Elements.validate(result)
    diagnostic.ok || throw(ErrorException(
        "$_CALLER: internal MixedMesh validation failed — " *
        join(diagnostic.messages, "; ")))
    (size(result.coords, 2) == node_count && length(result.blocks) == 2 &&
     size(result.blocks[1].nodes) == (4, quadrangle_count) &&
     size(result.blocks[2].nodes) == (8, hexahedron_count)) ||
        throw(ErrorException("$_CALLER: output count postcondition failed"))
    return result
end

end # module TransfiniteHex
