"""
    TransfinitePrism

Bounded Gmsh-4.15.2-compatible five-face transfinite meshing for a triangular
prism. The result is a first-order simplex `Mesh` with the exact unrecombined
tetrahedron topology used by Gmsh's volume path — the legacy collapsed-grid
subdivision by default, and with `compact=true` the compact
transfinite-triangle subdivision (`Mesh.TransfiniteTri = 1`, Gmsh's
`transfinite3` branch) whose full-square tab grid leaves the
upper-triangular interior nodes as unreferenced orphans. Corners alone
describe an affine prism; the optional `faces=` input supplies meshed
boundary face grids on the degenerate-hexahedron slot map (`s3≡s0`, `s7≡s4`),
lifting the affine restriction so curved and warped boundaries interpolate
through Gmsh's `transfiniteHex` formula exactly like the six-face path.

This module deliberately does not claim support for independently discretized
mismatched faces, nonuniform curve laws without face grids, recombined prisms
or hexahedra, QuadTri, holes, multiple blocks, periodic seams, or high-order
elements. As required by the finalized `Mesh` contract, represented boundary
areas and tetrahedron volumes must also remain finite Float64 values; finite
input coordinates alone do not imply that their derived measures are finite.
"""
module TransfinitePrism

using ..MeshTypes: Mesh, boundary_faces, nnodes, ntris, ntets, validate,
                   _throw_simplex_validation
using ..Predicates: orient3
using ..StructuredNumerics: _certify_tet_volume
using ..TransfiniteVolume: _WarpedFace, _chord_ratios, _face_node,
                           _transfinite_hex

export mesh_transfinite_prism

const _CALLER = "mesh_transfinite_prism"
const _DEFAULT_MAX_NODES = 10_000_000
const _DEFAULT_MAX_TETS = 60_000_000
const _DEFAULT_MAX_BOUNDARY_TRIANGLES = 20_000_000
const _AFFINE_TOLERANCE = 4096eps(Float64)

function _checked_add(a::Int, b::Int, what::AbstractString)
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
            values[cursor] = _count(value, ("radial", "opposite", "axial")[cursor])
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
    count == 5 || throw(ArgumentError(
        "$_CALLER: face_tags must contain exactly five tags"))
    tags = Vector{Int32}(undef, 5)
    cursor = 1
    try
        for value in raw
            cursor <= 5 || throw(ArgumentError(
                "$_CALLER: face_tags iteration produced more than five tags"))
            tags[cursor] = _tag(value, "face tag $cursor")
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read face_tags: $(sprint(showerror, err))"))
    end
    cursor == 6 || throw(ArgumentError(
        "$_CALLER: face_tags iteration ended before five tags"))
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
    values = try
        (raw[1], raw[2], raw[3])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: corner $index coordinates must be Float64-representable: " *
            sprint(showerror, err)))
    end
    any(value -> value isa Bool, values) && throw(ArgumentError(
        "$_CALLER: corner $index coordinates must not be Bool"))
    point = try
        (Float64(values[1]), Float64(values[2]), Float64(values[3]))
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

function _six_corners(raw)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: corners must be an indexable collection"))
    end
    count == 6 || throw(ArgumentError(
        "$_CALLER: corners must contain exactly six points"))
    corners = Vector{NTuple{3,Float64}}(undef, 6)
    cursor = 1
    try
        for value in raw
            cursor <= 6 || throw(ArgumentError(
                "$_CALLER: corners iteration produced more than six points"))
            corners[cursor] = _point3(value, cursor)
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read corners: $(sprint(showerror, err))"))
    end
    cursor == 7 || throw(ArgumentError(
        "$_CALLER: corners iteration ended before six points"))
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
        values = ntuple(index -> Rational{BigInt}(corners[index][dimension]), 6)
        values[5] + values[1] == values[2] + values[4] || return false
        values[6] + values[1] == values[3] + values[4] || return false
    end
    return true
end

function _certify_affine(corners)
    origin = corners[1]
    scale = 0.0
    deltas = Vector{NTuple{3,Float64}}(undef, 6)
    @inbounds for index in 1:6
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

    orientation = orient3(corners[1], corners[2], corners[3], corners[4])
    orientation < 0 || throw(ArgumentError(
        orientation == 0 ?
        "$_CALLER: canonical base and axial corner directions are coplanar" :
        "$_CALLER: corners must use the positive canonical Gmsh order " *
        "(s0,s1,s2,s4,s5,s6)"))

    u = normalized[2]
    v = normalized[3]
    w = normalized[4]
    determinant = abs(_dot3(u, _cross3(v, w)))
    frobenius_squared = _dot3(u, u) + _dot3(v, v) + _dot3(w, w)
    conditioning_bound = determinant / frobenius_squared
    if !(isfinite(conditioning_bound) && conditioning_bound > 0)
        _exact_affine_corners(corners) && return nothing
        throw(ArgumentError(
            "$_CALLER: affine corner conditioning is not representable in Float64"))
    end
    tolerance = _AFFINE_TOLERANCE * min(1.0, conditioning_bound)
    expected5 = _add3(u, w)
    expected6 = _add3(v, w)
    maximum_error = max(_maxabs3(_sub3(normalized[5], expected5)),
                        _maxabs3(_sub3(normalized[6], expected6)))
    isfinite(maximum_error) || throw(ArgumentError(
        "$_CALLER: affine residual is not finite"))
    if maximum_error > tolerance && !_exact_affine_corners(corners)
        throw(ArgumentError(
            "$_CALLER: corners do not form an affine triangular prism " *
            "(normalized residual $maximum_error exceeds the " *
            "conditioning-scaled tolerance $tolerance)"))
    end
    return nothing
end

@inline _lerp(a::Float64, b::Float64, t::Float64) = (1 - t) * a + t * b

@inline function _prism_point(corners, radial::Float64,
                              opposite::Float64, axial::Float64)
    return ntuple(3) do dimension
        lower_edge = _lerp(corners[2][dimension], corners[3][dimension], opposite)
        upper_edge = _lerp(corners[5][dimension], corners[6][dimension], opposite)
        lower = _lerp(corners[1][dimension], lower_edge, radial)
        upper = _lerp(corners[4][dimension], upper_edge, radial)
        _lerp(lower, upper, axial)
    end
end

@inline function _node(coords, node_id::Int32)
    index = Int(node_id)
    return (coords[1, index], coords[2, index], coords[3, index])
end

function _emit_canonical_tet!(tets, position::Int, coords,
                              a::Int32, b::Int32, c::Int32, d::Int32)
    sign = orient3(_node(coords, a), _node(coords, b),
                   _node(coords, c), _node(coords, d))
    sign == 0 && throw(ArgumentError(
        "$_CALLER: interpolation produced a zero-volume tetrahedron " *
        "at output position $position"))
    sign < 0 || throw(ArgumentError(
        "$_CALLER: interpolation reversed canonical tetrahedron $position; " *
        "the represented grid is folded"))
    @inbounds begin
        tets[1, position] = a
        tets[2, position] = b
        tets[3, position] = c
        tets[4, position] = d
    end
    return nothing
end

function _write_outward_triangle!(tris, tags, position::Int, coords,
                                  a::Int32, b::Int32, c::Int32,
                                  opposite, tag::Int32)
    sign = orient3(_node(coords, a), _node(coords, b), _node(coords, c), opposite)
    sign == 0 && throw(ArgumentError(
        "$_CALLER: interpolation produced a degenerate boundary triangle " *
        "at output position $position"))
    @inbounds begin
        tris[1, position] = a
        if sign > 0
            tris[2, position] = b
            tris[3, position] = c
        else
            tris[2, position] = c
            tris[3, position] = b
        end
        tags[position] = tag
    end
    return position + 1
end

function _canonical_triangles(tris)
    result = Vector{NTuple{3,Int32}}(undef, size(tris, 2))
    @inbounds for triangle in axes(tris, 2)
        a = tris[1, triangle]
        b = tris[2, triangle]
        c = tris[3, triangle]
        result[triangle] = a <= b ?
            (a <= c ? (b <= c ? (a, b, c) : (a, c, b)) : (c, a, b)) :
            (b <= c ? (a <= c ? (b, a, c) : (b, c, a)) : (c, b, a))
    end
    sort!(result)
    return result
end

function _certify_volume(coords, tets, corners)
    return _certify_tet_volume(
        coords, tets, (corners[1], corners[2], corners[3], corners[4]),
        3, _CALLER, "affine prism")
end

# ---------------------------------------------------------------------------
# Boundary face-grid input (Gmsh's oriented transfinite faces).
#
# A five-face prism maps onto the degenerate hexahedron of
# meshGRegionTransfinite.cpp with s3≡s0 and s7≡s4. Face records arrive in
# canonical slot order (f0,f1,f2,f4,f5) — f0=(s0,s1,s5,s4), f1=(s1,s2,s6,s5),
# f2=(s0,s2,s6,s4) are the axial quadrilaterals and f4/f5 the lower/upper
# collapsed triangular fills. Each `points` matrix is the face mesh reindexed
# onto the slot's canonical grid: (nr+1)×(nw+1) for f0/f2, (ns+1)×(nw+1) for
# f1, and the (nr+1)×(ns+1) expanded collapsed grid for f4/f5 whose i=0
# column repeats the apex vertex bitwise.

@inline function _prism_slot_dims(slot::Int, nr::Int, ns::Int, nw::Int)
    slot == 1 && return (nr, nw)
    slot == 2 && return (ns, nw)
    slot == 3 && return (nr, nw)
    slot == 4 && return (nr, ns)
    slot == 5 && return (nr, ns)
    throw(ArgumentError("$_CALLER: face slot must lie in 1:5; got $slot"))
end

function _prism_face(raw, slot::Int, nr::Int, ns::Int, nw::Int,
                     compact::Bool=false)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$_CALLER: faces[$slot] is not indexable"))
    end
    count == 3 || throw(ArgumentError(
        "$_CALLER: faces[$slot] must be a (points, triangles, tags) triple"))
    n1, n2 = _prism_slot_dims(slot, nr, ns, nw)
    expected = _checked_mul("faces[$slot] node", n1 + 1, n2 + 1)

    points = try
        raw[1]
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read faces[$slot] points"))
    end
    points isa AbstractMatrix || throw(ArgumentError(
        "$_CALLER: faces[$slot] points must be a matrix"))
    size(points, 1) == 3 || throw(ArgumentError(
        "$_CALLER: faces[$slot] points must have three rows"))
    size(points, 2) == expected || throw(ArgumentError(
        "$_CALLER: faces[$slot] carries $(size(points, 2)) nodes but the " *
        "slot requires $expected"))
    converted = Matrix{Float64}(undef, 3, expected)
    @inbounds for column in 1:expected, row in 1:3
        value = Float64(points[row, column])
        isfinite(value) || throw(ArgumentError(
            "$_CALLER: faces[$slot] node $column is not finite"))
        converted[row, column] = value
    end
    # The i=0 column of a collapsed triangular grid is the apex vertex
    # repeated — the degenerate hexahedral edge it represents is a single tab
    # node per axial layer.
    if slot >= 4
        apex = _face_node(converted, 0, 0, n1 + 1)
        @inbounds for j in 1:n2
            _face_node(converted, 0, j, n1 + 1) == apex || throw(
                ArgumentError(
                    "$_CALLER: faces[$slot] column i=0 does not repeat the " *
                    "collapsed apex vertex at j=$j"))
        end
    end

    tris = try
        raw[2]
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read faces[$slot] triangles"))
    end
    tris isa AbstractMatrix || throw(ArgumentError(
        "$_CALLER: faces[$slot] triangles must be a matrix"))
    size(tris, 1) == 3 || throw(ArgumentError(
        "$_CALLER: faces[$slot] triangles must have three rows"))
    ntri = size(tris, 2)
    expected_tris = if slot >= 4
        compact ? _checked_mul("faces[$slot] triangle", n1, n2) :
            _checked_mul("faces[$slot] triangle", n2,
                         _checked_add(2 * n1, -1, "triangle column"))
    else
        _checked_mul("faces[$slot] triangle", 2, n1, n2)
    end
    ntri == expected_tris || throw(ArgumentError(
        "$_CALLER: faces[$slot] carries $ntri triangles but the slot " *
        "requires $expected_tris"))
    converted_tris = Matrix{Int32}(undef, 3, ntri)
    @inbounds for column in 1:ntri, row in 1:3
        index = tris[row, column]
        index isa Integer || throw(ArgumentError(
            "$_CALLER: faces[$slot] triangle $column vertex $row is not " *
            "an integer"))
        1 <= index <= expected || throw(ArgumentError(
            "$_CALLER: faces[$slot] triangle $column vertex $row indexes " *
            "node $index outside 1:$expected"))
        converted_tris[row, column] = Int32(index)
    end

    tags = try
        raw[3]
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read faces[$slot] tags"))
    end
    length(tags) == ntri || throw(ArgumentError(
        "$_CALLER: faces[$slot] carries $(length(tags)) tags for $ntri " *
        "triangles"))
    converted_tags = Vector{Int32}(undef, ntri)
    @inbounds for index in 1:ntri
        converted_tags[index] = _tag(tags[index], "faces[$slot] tag $index")
    end
    return _WarpedFace(converted, converted_tris, converted_tags)
end

function _prism_faces(raw, nr::Int, ns::Int, nw::Int, compact::Bool=false)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$_CALLER: faces is not indexable"))
    end
    count == 5 || throw(ArgumentError(
        "$_CALLER: faces must contain exactly five records"))
    faces = Vector{_WarpedFace}(undef, 5)
    cursor = 1
    try
        for entry in raw
            cursor <= 5 || throw(ArgumentError(
                "$_CALLER: faces produced more than five records"))
            faces[cursor] = _prism_face(entry, cursor, nr, ns, nw, compact)
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "$_CALLER: could not read faces: $(sprint(showerror, err))"))
    end
    cursor == 6 || throw(ArgumentError(
        "$_CALLER: faces ended before five records"))
    return faces
end

# The nine real boundary edges shared between adjacent prism face slots — the
# degenerate hexahedron's twelfth edge is the s0s4 collapsed axial edge shared
# by f0 and f2. Adjacent slot grids must agree bitwise so the welded boundary
# nodes are unambiguous.
function _certify_shared_prism_edges(faces, nr::Int, ns::Int, nw::Int)
    f0, f1, f2, f4, f5 = faces
    npu, npv, npw = nr + 1, ns + 1, nw + 1
    function check(name, a, b)
        a == b || throw(ArgumentError(
            "$_CALLER: boundary faces disagree on shared edge $name " *
            "($(a) != $(b))"))
        return nothing
    end
    @inbounds begin
        for i in 0:nr
            check("s0s1", _face_node(f4.points, i, 0, npu),
                  _face_node(f0.points, i, 0, npu))
            check("s0s2", _face_node(f4.points, i, ns, npu),
                  _face_node(f2.points, i, 0, npu))
            check("s4s5", _face_node(f5.points, i, 0, npu),
                  _face_node(f0.points, i, nw, npu))
            check("s4s6", _face_node(f5.points, i, ns, npu),
                  _face_node(f2.points, i, nw, npu))
        end
        for j in 0:ns
            check("s1s2", _face_node(f4.points, nr, j, npu),
                  _face_node(f1.points, j, 0, npv))
            check("s5s6", _face_node(f5.points, nr, j, npu),
                  _face_node(f1.points, j, nw, npv))
        end
        for k in 0:nw
            check("s0s4", _face_node(f0.points, 0, k, npu),
                  _face_node(f2.points, 0, k, npu))
            check("s1s5", _face_node(f0.points, nr, k, npu),
                  _face_node(f1.points, 0, k, npv))
            check("s2s6", _face_node(f2.points, nr, k, npu),
                  _face_node(f1.points, ns, k, npv))
        end
    end
    return nothing
end

# The six canonical corners must coincide bitwise with the collapsed
# triangular grids' corners and keep the positive canonical orientation the
# tetrahedron templates assume.
function _certify_prism_corners(corners, faces, nr::Int, ns::Int)
    f4, f5 = faces[4], faces[5]
    npu = nr + 1
    _face_node(f4.points, 0, 0, npu) == corners[1] &&
        _face_node(f4.points, nr, 0, npu) == corners[2] &&
        _face_node(f4.points, nr, ns, npu) == corners[3] &&
        _face_node(f5.points, 0, 0, npu) == corners[4] &&
        _face_node(f5.points, nr, 0, npu) == corners[5] &&
        _face_node(f5.points, nr, ns, npu) == corners[6] ||
        throw(ArgumentError(
            "$_CALLER: corners do not match the face grids"))
    orientation = orient3(corners[1], corners[2], corners[3], corners[4])
    orientation < 0 || throw(ArgumentError(
        orientation == 0 ?
        "$_CALLER: canonical base and axial corner directions are coplanar" :
        "$_CALLER: corners must use the positive canonical Gmsh order " *
        "(s0,s1,s2,s4,s5,s6)"))
    return nothing
end

# Tab-grid coordinate fill for the face-grid path: boundary nodes reuse the
# face grids bitwise (shared edges were certified identical) and interior
# nodes evaluate Gmsh's `transfiniteHex` on the degenerate hexahedral slot
# layout — the collapsed u=0 face reads the s0s4 axial edge node (c8).
function _fill_warped_prism!(coords, corners, faces, nr::Int, ns::Int,
                             nw::Int, node_id)
    f0, f1, f2, f4, f5 = faces
    npu, npv = nr + 1, ns + 1
    us = _chord_ratios(
        [_face_node(f4.points, i, 0, npu) for i in 0:nr], "s0s1")
    vs = _chord_ratios(
        [_face_node(f1.points, j, 0, npv) for j in 0:ns], "s1s2")
    ws = _chord_ratios(
        [_face_node(f1.points, 0, k, npv) for k in 0:nw], "s1s5")
    s = (corners[1], corners[2], corners[3], corners[1],
         corners[4], corners[5], corners[6], corners[4])
    @inbounds for k in 0:nw
        point = _face_node(f0.points, 0, k, npu)
        index = Int(node_id(0, 0, k))
        coords[1, index] = point[1]
        coords[2, index] = point[2]
        coords[3, index] = point[3]
    end
    @inbounds for i in 1:nr, j in 0:ns, k in 0:nw
        boundary = i == nr || j == 0 || j == ns || k == 0 || k == nw
        point = if boundary
            i == nr ? _face_node(f1.points, j, k, npv) :
            j == 0 ? _face_node(f0.points, i, k, npu) :
            j == ns ? _face_node(f2.points, i, k, npu) :
            k == 0 ? _face_node(f4.points, i, j, npu) :
                     _face_node(f5.points, i, j, npu)
        else
            u = us[i + 1]; v = vs[j + 1]; w = ws[k + 1]
            evaluated = _transfinite_hex(
                (_face_node(f0.points, 0, k, npu),
                 _face_node(f1.points, j, k, npv),
                 _face_node(f0.points, i, k, npu),
                 _face_node(f2.points, i, k, npu),
                 _face_node(f4.points, i, j, npu),
                 _face_node(f5.points, i, j, npu)),
                (_face_node(f4.points, i, 0, npu),
                 _face_node(f5.points, i, 0, npu),
                 _face_node(f4.points, i, ns, npu),
                 _face_node(f5.points, i, ns, npu)),
                (_face_node(f4.points, 0, j, npu),
                 _face_node(f4.points, nr, j, npu),
                 _face_node(f5.points, 0, j, npu),
                 _face_node(f5.points, nr, j, npu)),
                (_face_node(f0.points, 0, k, npu),
                 _face_node(f2.points, 0, k, npu),
                 _face_node(f0.points, nr, k, npu),
                 _face_node(f2.points, nr, k, npu)),
                s, u, v, w)
            all(isfinite, evaluated) || throw(ArgumentError(
                "$_CALLER: interpolation produced a non-finite coordinate " *
                "at logical node ($i,$j,$k)"))
            evaluated
        end
        index = Int(node_id(i, j, k))
        coords[1, index] = point[1]
        coords[2, index] = point[2]
        coords[3, index] = point[3]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Compact `Mesh.TransfiniteTri = 1` subdivision — Gmsh's `transfinite3`
# branch of `MeshTransfiniteVolume`. The two triangular faces keep their
# compact equal-sided kernels expanded onto a *full* (nr+1)×(ns+1) grid —
# slots past the diagonal alias the diagonal vertex — so the tab stays
# square and the tetrahedron loops emit `SIM_10`–`SIM_12` on diagonal cells
# and `SIM_7`–`SIM_12` on strictly lower-triangular cells. Interior slots
# behind the diagonal (`j > i`) become unreferenced orphan nodes exactly like
# upstream.

# Welds the two triangular boundary planes by bitwise coordinate identity —
# the expanded compact grids repeat the aliased diagonal vertices — then
# hands interior slots (i >= 1) one fresh node each. Returns the
# `node_id(i,j,k)` closure plus the exact node count (two compact lattices
# plus one apex and nr*(ns+1) fresh nodes per interior layer).
function _compact_prism_ids(f4::Matrix{Float64}, f5::Matrix{Float64},
                            nr::Int, ns::Int, nw::Int)
    n2 = ns + 1
    plane = _checked_mul("node", nr + 1, n2)
    lattice = div(_checked_mul("node", nr + 1, nr + 2), 2)
    block = _checked_add(_checked_mul("node", nr, n2), 1, "node")
    interior = _checked_mul(
        "node", block, _checked_add(nw, -1, "node"))
    base5 = _checked_add(lattice, interior, "node")
    node_count = _checked_add(base5, lattice, "node")
    # The welded ids index into the coords columns — guard before any Int32
    # conversion below.
    node_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $node_count nodes exceed Int32 indexing"))
    function weld(points)
        ids = Vector{Int32}(undef, plane)
        seen = Dict{NTuple{3,Float64},Int32}()
        next = Ref(0)
        @inbounds for lin in 1:plane
            key = (points[1, lin], points[2, lin], points[3, lin])
            ids[lin] = get!(seen, key) do
                Int32(next[] += 1)
            end
        end
        next[] == lattice || throw(ArgumentError(
            "$_CALLER: compact triangular face expansion carries " *
            "$(next[]) distinct nodes but the prism slot requires $lattice"))
        return ids
    end
    id4 = weld(f4)
    id5 = weld(f5)
    node_id(i::Int, j::Int, k::Int) =
        k == 0 ? id4[i + j * n2 + 1] :
        k == nw ? Int32(base5) + id5[i + j * n2 + 1] :
        i == 0 ? Int32(lattice + (k - 1) * block + 1) :
            Int32(lattice + (k - 1) * block + 2 + (i - 1) * n2 + j)
    return node_id, node_count
end

# Full-square expanded face grids for an affine compact prism: canonical
# slot (i,j) on a triangular face holds the compact lattice vertex
# (i, min(j,i)) — the diagonal fold the kernel welds back together.
function _synthetic_compact_prism_faces(corners, n::Int, nw::Int)
    np = n + 1
    npw = nw + 1
    tri = Matrix{Float64}(undef, 3, _checked_mul("node", np, np))
    f4 = tri
    f5 = Matrix{Float64}(undef, 3, size(tri, 2))
    @inbounds for i in 0:n, j in 0:n
        opposite = i == 0 ? 0.0 : min(j, i) / i
        lower = _prism_point(corners, i / n, opposite, 0.0)
        upper = _prism_point(corners, i / n, opposite, 1.0)
        index = i + j * np + 1
        f4[1, index] = lower[1]
        f4[2, index] = lower[2]
        f4[3, index] = lower[3]
        f5[1, index] = upper[1]
        f5[2, index] = upper[2]
        f5[3, index] = upper[3]
    end
    quads = ntuple(_ -> Matrix{Float64}(undef, 3,
                                        _checked_mul("node", np, npw)), 3)
    f0, f1, f2 = quads
    @inbounds for i in 0:n, k in 0:nw
        index = i + k * np + 1
        a = _prism_point(corners, i / n, 0.0, k / nw)
        f0[1, index] = a[1]; f0[2, index] = a[2]; f0[3, index] = a[3]
        b = _prism_point(corners, i / n, 1.0, k / nw)
        f2[1, index] = b[1]; f2[2, index] = b[2]; f2[3, index] = b[3]
        c = _prism_point(corners, 1.0, i / n, k / nw)
        f1[1, index] = c[1]; f1[2, index] = c[2]; f1[3, index] = c[3]
    end
    quadtris = Matrix{Int32}(undef, 3, _checked_mul(
        "boundary triangle", 2, n, nw))
    cursor = 0
    @inbounds for i in 0:n-1, k in 0:nw-1
        a = i + k * np + 1
        b = a + 1
        d = a + np
        e = b + np
        quadtris[1, cursor += 1] = a
        quadtris[2, cursor] = b
        quadtris[3, cursor] = d
        quadtris[1, cursor += 1] = b
        quadtris[2, cursor] = e
        quadtris[3, cursor] = d
    end
    # Compact triangulation: {(i,j),(i+1,j+1),(i,j+1)} plus
    # {(i,j),(i+1,j),(i+1,j+1)} per strictly-lower cell and the single
    # diagonal triangle {(i,i),(i+1,i),(i+1,i+1)}.
    tritris2 = Matrix{Int32}(undef, 3, _checked_mul(
        "boundary triangle", n, n))
    cursor = 0
    @inbounds for i in 0:n-1, j in 0:i-1
        a = i + j * np + 1
        b = i + 1 + j * np + 1
        cp = i + (j + 1) * np + 1
        g = i + 1 + (j + 1) * np + 1
        tritris2[1, cursor += 1] = a
        tritris2[2, cursor] = g
        tritris2[3, cursor] = cp
        tritris2[1, cursor += 1] = a
        tritris2[2, cursor] = b
        tritris2[3, cursor] = g
    end
    @inbounds for i in 0:n-1
        a = i + i * np + 1
        b = i + 1 + i * np + 1
        c = i + 1 + (i + 1) * np + 1
        tritris2[1, cursor += 1] = a
        tritris2[2, cursor] = b
        tritris2[3, cursor] = c
    end
    tags0 = zeros(Int32, size(quadtris, 2))
    tags4 = zeros(Int32, size(tritris2, 2))
    return (_WarpedFace(f0, quadtris, tags0),
            _WarpedFace(f1, quadtris, tags0),
            _WarpedFace(f2, quadtris, tags0),
            _WarpedFace(f4, tritris2, tags4),
            _WarpedFace(f5, tritris2, tags4))
end

"""
    mesh_transfinite_prism(corners, cells=(1,1,1);
                          volume_tag=0,
                          face_tags=(0,0,0,0,0),
                          faces=nothing,
                          compact=false,
                          max_nodes=10_000_000,
                          max_tets=60_000_000,
                          max_boundary_triangles=20_000_000) -> Mesh

Mesh a triangular prism using Gmsh 4.15.2's five-face transfinite algorithm
with all surfaces unrecombined and using the `Left` triangle arrangement.
`corners` must contain six finite points in canonical order
`(s0,s1,s2,s4,s5,s6)`: the first three define the lower triangular face and
the final three are their axial counterparts. The canonical orientation must
be positive.

With `compact=false` (the default, `Mesh.TransfiniteTri = 0`) the prism uses
the legacy collapsed-grid subdivision: Gmsh's exact three-tetrahedron
collapsed-wedge pattern and six-tetrahedron interior-block pattern without
per-cell vertex reordering. With `compact=true` (`Mesh.TransfiniteTri = 1`,
Gmsh's `transfinite3` branch selected whenever either triangular face carries
the compact transfinite surface) the two triangular faces keep their
equal-sided compact grids expanded onto a full `nr×ns` cell square — slots
past the diagonal alias the diagonal vertex — diagonal cells emit the
`SIM_10`/`SIM_11`/`SIM_12` tetrahedron templates and strictly lower cells
emit `SIM_7` through `SIM_12`; interior grid slots behind the diagonal
(`j > i`) are filled with `transfiniteHex` coordinates but left unreferenced
as orphan nodes exactly like upstream. Compact mode requires `nr == ns`.

`cells=(nr,ns,nw)` gives positive logical-cell counts. The two sides incident
to collapsed corner `s0` each have `nr` cells, side `s1-s2` has `ns`, and all
three axial sides have `nw`; curve spacing is uniform (`Progression 1`). The
returned mesh uses Gmsh's exact three-tetrahedron collapsed-wedge pattern and
six-tetrahedron interior-block pattern without per-cell vertex reordering.
Exact predicates reject any zero or reversed canonical tetrahedron, the emitted
boundary is certified against the tetrahedron boundary, and a compensated
exponent-scaled determinant audit with an exact dyadic fallback rejects
material loss or overlap caused by unrepresentable intermediate coordinates.
`face_tags` follow canonical order
`(f0,f1,f2,f4,f5)`, where `f0=(s0,s1,s5,s4)`,
`f1=(s1,s2,s6,s5)`, `f2=(s0,s2,s6,s4)`, and `f4`/`f5` are the
lower/upper triangular faces.

With `faces=nothing` the prism must be affine: corners `s5`/`s6` are certified
against `s1+(s4-s0)`/`s2+(s4-s0)` and all nodes are trilinear interpolants.
With `faces`, it is a five-tuple of `(points, triangles, tags)` records in the
canonical slot order `(f0,f1,f2,f4,f5)` carrying the meshed boundary grids
reindexed onto the degenerate-hexahedron slots (`s3≡s0`, `s7≡s4`): the axial
quadrilaterals are `(nr+1)×(nw+1)`/`(ns+1)×(nw+1)` grids and each triangular
face is its `(nr+1)×(ns+1)` expanded collapsed grid whose `i=0` column
repeats the apex vertex bitwise. Boundary tab nodes reuse the face grids
bitwise, interior nodes evaluate Gmsh's `transfiniteHex` on the degenerate
slot map with chord-length parameters measured on the `s0s1`/`s1s2`/`s1s5`
reference edges, and each boundary triangle is oriented against its incident
tetrahedron's inward apex — so curved, warped, and ruled boundary faces are
supported. The `corners` need not be affine in this path, but must coincide
bitwise with the face grids' six corner vertices.

Unsupported here: recombination, QuadTri, holes, multiple blocks, periodic
seams, high-order elements, and coordinate scales whose derived triangle
areas or tetrahedron volumes are not finite Float64 values.
"""
function mesh_transfinite_prism(corners, cells=(1, 1, 1);
                                volume_tag=0,
                                face_tags=(0, 0, 0, 0, 0),
                                faces=nothing,
                                compact=false,
                                max_nodes=_DEFAULT_MAX_NODES,
                                max_tets=_DEFAULT_MAX_TETS,
                                max_boundary_triangles=
                                    _DEFAULT_MAX_BOUNDARY_TRIANGLES)
    node_limit = _limit(max_nodes, "max_nodes")
    tet_limit = _limit(max_tets, "max_tets")
    triangle_limit = _limit(max_boundary_triangles, "max_boundary_triangles")
    nr, ns, nw = _three_counts(cells)

    compact isa Bool || throw(ArgumentError(
        "$_CALLER: compact must be a Bool, got $(typeof(compact))"))
    if compact
        return _mesh_transfinite_prism_compact(
            corners, nr, ns, nw, faces;
            volume_tag=volume_tag, face_tags=face_tags,
            node_limit=node_limit, tet_limit=tet_limit,
            triangle_limit=triangle_limit)
    end

    opposite_nodes = _checked_add(ns, 1, "node")
    axial_nodes = _checked_add(nw, 1, "node")
    layer_nodes = _checked_add(
        _checked_mul("node", nr, opposite_nodes), 1, "node")
    node_count = _checked_mul("node", layer_nodes, axial_nodes)
    two_radial_minus_one = _checked_add(
        _checked_mul("tetrahedron", 2, nr), -1, "tetrahedron")
    tet_count = _checked_mul(
        "tetrahedron", 3, ns, nw, two_radial_minus_one)
    base_triangles = _checked_mul(
        "boundary triangle", ns, two_radial_minus_one)
    radial_side_triangles = _checked_mul(
        "boundary triangle", 4, nr, nw)
    opposite_side_triangles = _checked_mul(
        "boundary triangle", 2, ns, nw)
    triangle_count = _checked_add(
        _checked_add(_checked_mul("boundary triangle", 2, base_triangles),
                     radial_side_triangles, "boundary triangle"),
        opposite_side_triangles, "boundary triangle")

    node_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $node_count nodes exceed Int32 indexing"))
    tet_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $tet_count tetrahedra exceed the Int32 topology limit"))
    triangle_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $triangle_count boundary triangles exceed the Int32 topology limit"))
    node_count <= node_limit || throw(ArgumentError(
        "$_CALLER: $node_count nodes exceed max_nodes=$node_limit"))
    tet_count <= tet_limit || throw(ArgumentError(
        "$_CALLER: $tet_count tetrahedra exceed max_tets=$tet_limit"))
    triangle_count <= triangle_limit || throw(ArgumentError(
        "$_CALLER: $triangle_count boundary triangles exceed " *
        "max_boundary_triangles=$triangle_limit"))

    converted_corners = _six_corners(corners)
    converted_face_tags = _face_tags(face_tags)
    converted_volume_tag = _tag(volume_tag, "volume_tag")

    node_id(i::Int, j::Int, k::Int) = Int32(
        k * layer_nodes + (i == 0 ? 1 : 2 + (i - 1) * opposite_nodes + j))

    coords = Matrix{Float64}(undef, 3, node_count)
    if faces === nothing
        _certify_affine(converted_corners)
        @inbounds for k in 0:nw
            axial = k / nw
            collapsed = _prism_point(converted_corners, 0.0, 0.0, axial)
            all(isfinite, collapsed) || throw(ArgumentError(
                "$_CALLER: interpolation produced a non-finite coordinate " *
                "at logical node (0,0,$k)"))
            collapsed_index = Int(node_id(0, 0, k))
            coords[1, collapsed_index] = collapsed[1]
            coords[2, collapsed_index] = collapsed[2]
            coords[3, collapsed_index] = collapsed[3]
            for i in 1:nr, j in 0:ns
                point = _prism_point(
                    converted_corners, i / nr, j / ns, axial)
                all(isfinite, point) || throw(ArgumentError(
                    "$_CALLER: interpolation produced a non-finite " *
                    "coordinate at logical node ($i,$j,$k)"))
                index = Int(node_id(i, j, k))
                coords[1, index] = point[1]
                coords[2, index] = point[2]
                coords[3, index] = point[3]
            end
        end
    else
        wfaces = _prism_faces(faces, nr, ns, nw)
        _certify_shared_prism_edges(wfaces, nr, ns, nw)
        _certify_prism_corners(converted_corners, wfaces, nr, ns)
        _fill_warped_prism!(coords, converted_corners, wfaces,
                            nr, ns, nw, node_id)
    end

    # Boundary-triangle inward witnesses. The affine path keeps the fixed
    # opposite-face corner reference it always used; the face-grid path must
    # follow warped boundary sheets, so it orients each emitted triangle
    # against an inward tab vertex of the incident tetrahedron — strictly
    # off-plane whenever the emitted tet is nondegenerate.
    witness = if faces === nothing
        fixed = (converted_corners[3], converted_corners[1],
                 converted_corners[2], converted_corners[4],
                 converted_corners[1])
        (slot::Int, ::Int, ::Int, ::Int, ::Int) -> fixed[slot]
    else
        # Inward vertex of the emitted triangle's incident tetrahedron —
        # strictly off the triangle's supporting plane on the interior side
        # whenever that tet is nondegenerate, so the sign test orients the
        # boundary correctly on warped sheets. The collapsed wedge and the
        # interior block split the boundary cells differently, so slot 2/3
        # witnesses are half- and radius-dependent.
        function (slot::Int, i::Int, j::Int, k::Int, half::Int)
            slot == 1 && return _node(coords, node_id(max(i, 1), 1, k))
            if slot == 2
                layer = half == 1 && nr > 1 ? k : k + 1
                return _node(coords, node_id(nr - 1, j + 1, layer))
            end
            if slot == 3
                layer = i == 0 && half == 1 ? k : k + 1
                return _node(coords, node_id(i + 1, ns - 1, layer))
            end
            if slot == 4
                i == 0 && return _node(coords, node_id(0, 0, 1))
                return _node(coords, node_id(i + half - 1, j, 1))
            end
            i == 0 && return _node(coords, node_id(1, j + 1, nw - 1))
            return _node(coords, node_id(i + half - 1, j + 1, nw - 1))
        end
    end

    tets = Matrix{Int32}(undef, 4, tet_count)
    tet_position = 0
    @inbounds for j in 0:ns-1, k in 0:nw-1
        a = node_id(0, 0, k)
        b = node_id(1, j, k)
        c = node_id(1, j + 1, k)
        d = node_id(0, 0, k + 1)
        e = node_id(1, j, k + 1)
        f = node_id(1, j + 1, k + 1)
        # Gmsh's collapsed-wedge path uses (d,f,e,c) for its third
        # tetrahedron. This differs from the interior CREATE_SIM_3 template
        # below even though the two tuples are even permutations.
        for vertices in ((a, b, c, d), (b, c, d, e), (d, f, e, c))
            tet_position += 1
            _emit_canonical_tet!(tets, tet_position, coords, vertices...)
        end
    end
    @inbounds for i in 1:nr-1, j in 0:ns-1, k in 0:nw-1
        a = node_id(i, j, k)
        b = node_id(i + 1, j, k)
        c = node_id(i, j + 1, k)
        d = node_id(i, j, k + 1)
        e = node_id(i + 1, j, k + 1)
        f = node_id(i, j + 1, k + 1)
        g = node_id(i + 1, j + 1, k)
        h = node_id(i + 1, j + 1, k + 1)
        for vertices in ((a, b, c, d), (b, c, d, e), (d, e, c, f),
                         (b, c, e, g), (c, f, e, g), (e, f, h, g))
            tet_position += 1
            _emit_canonical_tet!(tets, tet_position, coords, vertices...)
        end
    end
    tet_position == tet_count || throw(ErrorException(
        "$_CALLER: internal tetrahedron count invariant failed"))

    tris = Matrix{Int32}(undef, 3, triangle_count)
    tri_tags = Vector{Int32}(undef, triangle_count)
    position = 1
    # f0: s0-s1 axial side (opposite coordinate j=0).
    @inbounds for i in 0:nr-1, k in 0:nw-1
        a=node_id(i,0,k); b=node_id(i+1,0,k)
        d=node_id(i,0,k+1); e=node_id(i+1,0,k+1)
        opposite=witness(1,i,0,k,1)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,d,opposite,
            converted_face_tags[1])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,b,e,d,opposite,
            converted_face_tags[1])
    end
    # f1: s1-s2 axial side (outer radial row).
    @inbounds for j in 0:ns-1, k in 0:nw-1
        a=node_id(nr,j,k); b=node_id(nr,j+1,k)
        d=node_id(nr,j,k+1); e=node_id(nr,j+1,k+1)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,d,witness(2,nr,j,k,1),
            converted_face_tags[2])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,b,e,d,witness(2,nr,j,k,2),
            converted_face_tags[2])
    end
    # f2: s0-s2 axial side (opposite coordinate j=ns).
    @inbounds for i in 0:nr-1, k in 0:nw-1
        a=node_id(i,ns,k); b=node_id(i+1,ns,k)
        d=node_id(i,ns,k+1); e=node_id(i+1,ns,k+1)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,d,witness(3,i,ns,k,1),
            converted_face_tags[3])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,b,e,d,witness(3,i,ns,k,2),
            converted_face_tags[3])
    end
    # f4: lower triangular face.
    @inbounds for j in 0:ns-1
        a=node_id(0,0,0); b=node_id(1,j,0); c=node_id(1,j+1,0)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,c,witness(4,0,j,0,1),
            converted_face_tags[4])
    end
    @inbounds for i in 1:nr-1, j in 0:ns-1
        a=node_id(i,j,0); b=node_id(i+1,j,0)
        c=node_id(i,j+1,0); g=node_id(i+1,j+1,0)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,c,witness(4,i,j,0,1),
            converted_face_tags[4])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,c,b,g,witness(4,i,j,0,2),
            converted_face_tags[4])
    end
    # f5: upper triangular face.
    @inbounds for j in 0:ns-1
        a=node_id(0,0,nw); b=node_id(1,j,nw); c=node_id(1,j+1,nw)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,c,witness(5,0,j,nw,1),
            converted_face_tags[5])
    end
    @inbounds for i in 1:nr-1, j in 0:ns-1
        a=node_id(i,j,nw); b=node_id(i+1,j,nw)
        c=node_id(i,j+1,nw); g=node_id(i+1,j+1,nw)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,c,witness(5,i,j,nw,1),
            converted_face_tags[5])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,c,b,g,witness(5,i,j,nw,2),
            converted_face_tags[5])
    end
    position == triangle_count + 1 || throw(ErrorException(
        "$_CALLER: internal boundary triangle count invariant failed"))

    extracted_boundary, maximum_incidence = boundary_faces(tets)
    maximum_incidence == 2 || throw(ErrorException(
        "$_CALLER: constructed tet mesh has face incidence $maximum_incidence"))
    sort!(extracted_boundary)
    extracted_boundary == _canonical_triangles(tris) || throw(ErrorException(
        "$_CALLER: emitted boundary triangles do not match the tet boundary"))
    # The enclosed-volume audit compares against the affine corner prism; the
    # warped face-grid path encloses a genuinely different volume and relies
    # on the per-tet orientation certification instead.
    faces === nothing &&
        _certify_volume(coords, tets, converted_corners)

    mesh = Mesh(coords; tris=tris, tets=tets, tri_tag=tri_tags,
                tet_tag=fill(converted_volume_tag, tet_count))
    diagnostic = validate(mesh)
    diagnostic.ok || _throw_simplex_validation(_CALLER, diagnostic.messages)
    (nnodes(mesh), ntris(mesh), ntets(mesh)) ==
        (node_count, triangle_count, tet_count) || throw(ErrorException(
        "$_CALLER: finalized mesh count invariant failed"))
    return mesh
end

# Compact (`Mesh.TransfiniteTri = 1`) pipeline — the full-square tab with
# coordinate-welded triangular planes and the `SIM_7`–`SIM_12` template set.
# Grid slots behind the diagonal (`j > i` at interior layers) stay
# unreferenced orphan nodes, mirroring Gmsh exactly.
function _mesh_transfinite_prism_compact(corners, nr::Int, ns::Int, nw::Int,
                                         faces; volume_tag, face_tags,
                                         node_limit, tet_limit,
                                         triangle_limit)
    ns == nr || throw(ArgumentError(
        "$_CALLER: the compact transfinite-triangle subdivision requires " *
        "equal radial and opposite cell counts (nr=$nr != ns=$ns)"))
    converted_corners = _six_corners(corners)
    converted_face_tags = _face_tags(face_tags)
    converted_volume_tag = _tag(volume_tag, "volume_tag")
    affine = faces === nothing
    affine && _certify_affine(converted_corners)
    wfaces = affine ?
        _synthetic_compact_prism_faces(converted_corners, nr, nw) :
        _prism_faces(faces, nr, ns, nw, true)
    _certify_shared_prism_edges(wfaces, nr, ns, nw)
    _certify_prism_corners(converted_corners, wfaces, nr, ns)

    node_id, node_count = _compact_prism_ids(
        wfaces[4].points, wfaces[5].points, nr, ns, nw)
    n2 = ns + 1
    tet_count = _checked_mul("tetrahedron", 3, nr, nr, nw)
    side_triangles = _checked_mul(
        "boundary triangle", 2, nw,
        _checked_add(_checked_mul("boundary triangle", 2, nr), ns,
                     "boundary triangle"))
    triangle_count = _checked_add(
        side_triangles,
        _checked_mul("boundary triangle", 2, nr, nr),
        "boundary triangle")

    node_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $node_count nodes exceed Int32 indexing"))
    tet_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $tet_count tetrahedra exceed the Int32 topology limit"))
    triangle_count <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $triangle_count boundary triangles exceed the Int32 topology limit"))
    node_count <= node_limit || throw(ArgumentError(
        "$_CALLER: $node_count nodes exceed max_nodes=$node_limit"))
    tet_count <= tet_limit || throw(ArgumentError(
        "$_CALLER: $tet_count tetrahedra exceed max_tets=$tet_limit"))
    triangle_count <= triangle_limit || throw(ArgumentError(
        "$_CALLER: $triangle_count boundary triangles exceed " *
        "max_boundary_triangles=$triangle_limit"))

    coords = Matrix{Float64}(undef, 3, node_count)
    _fill_warped_prism!(coords, converted_corners, wfaces,
                        nr, ns, nw, node_id)

    tets = Matrix{Int32}(undef, 4, tet_count)
    tet_position = 0
    # Diagonal cells (j == i): upstream emits SIM_10/SIM_11/SIM_12.
    @inbounds for i in 0:nr-1, k in 0:nw-1
        a = node_id(i, i, k)
        b = node_id(i + 1, i, k)
        d = node_id(i, i, k + 1)
        e = node_id(i + 1, i, k + 1)
        g = node_id(i + 1, i + 1, k)
        h = node_id(i + 1, i + 1, k + 1)
        for vertices in ((a, b, g, d), (b, g, d, e), (d, e, g, h))
            tet_position += 1
            _emit_canonical_tet!(tets, tet_position, coords, vertices...)
        end
    end
    # Strictly-lower cells (j < i): SIM_7 through SIM_12.
    @inbounds for i in 1:nr-1, j in 0:i-1, k in 0:nw-1
        a = node_id(i, j, k)
        b = node_id(i + 1, j, k)
        c = node_id(i, j + 1, k)
        d = node_id(i, j, k + 1)
        e = node_id(i + 1, j, k + 1)
        f = node_id(i, j + 1, k + 1)
        g = node_id(i + 1, j + 1, k)
        h = node_id(i + 1, j + 1, k + 1)
        for vertices in ((a, c, d, g), (c, f, d, g), (d, f, h, g),
                         (a, b, g, d), (b, g, d, e), (d, e, g, h))
            tet_position += 1
            _emit_canonical_tet!(tets, tet_position, coords, vertices...)
        end
    end
    tet_position == tet_count || throw(ErrorException(
        "$_CALLER: internal tetrahedron count invariant failed"))

    tris = Matrix{Int32}(undef, 3, triangle_count)
    tri_tags = Vector{Int32}(undef, triangle_count)
    position = 1
    # f0 (j = 0): SIM_10/SIM_11 induced faces — diagonal cells included.
    @inbounds for i in 0:nr-1, k in 0:nw-1
        a = node_id(i, 0, k); b = node_id(i + 1, 0, k)
        d = node_id(i, 0, k + 1); e = node_id(i + 1, 0, k + 1)
        inward = _node(coords, node_id(i + 1, 1, k))
        position = _write_outward_triangle!(
            tris, tri_tags, position, coords, a, b, d, inward,
            converted_face_tags[1])
        position = _write_outward_triangle!(
            tris, tri_tags, position, coords, b, d, e, inward,
            converted_face_tags[1])
    end
    # f1 (i = nr): SIM_11's (e,b,g) and SIM_12's (e,g,h) faces; the j = nr-1
    # cell arrives from the diagonal cell of row nr-1.
    @inbounds for j in 0:nr-1, k in 0:nw-1
        b = node_id(nr, j, k); g = node_id(nr, j + 1, k)
        e = node_id(nr, j, k + 1); h = node_id(nr, j + 1, k + 1)
        inward = _node(coords, node_id(nr - 1, j, k + 1))
        position = _write_outward_triangle!(
            tris, tri_tags, position, coords, b, g, e, inward,
            converted_face_tags[2])
        position = _write_outward_triangle!(
            tris, tri_tags, position, coords, e, g, h, inward,
            converted_face_tags[2])
    end
    # f2 (diagonal plane): SIM_10's (a,g,d) and SIM_12's (d,g,h) faces.
    @inbounds for i in 0:nr-1, k in 0:nw-1
        a = node_id(i, i, k); g = node_id(i + 1, i + 1, k)
        d = node_id(i, i, k + 1); h = node_id(i + 1, i + 1, k + 1)
        inward0 = _node(coords, node_id(i + 1, i, k))
        inward1 = _node(coords, node_id(i + 1, i, k + 1))
        position = _write_outward_triangle!(
            tris, tri_tags, position, coords, a, g, d, inward0,
            converted_face_tags[3])
        position = _write_outward_triangle!(
            tris, tri_tags, position, coords, d, g, h, inward1,
            converted_face_tags[3])
    end
    # f4 (k = 0): SIM_10's (a,b,g) on every cell plus SIM_7's (a,c,g) on
    # strictly-lower cells — the compact face triangulation.
    @inbounds for i in 0:nr-1
        a = node_id(i, i, 0); b = node_id(i + 1, i, 0)
        g = node_id(i + 1, i + 1, 0)
        inward = _node(coords, node_id(i, i, 1))
        position = _write_outward_triangle!(
            tris, tri_tags, position, coords, a, b, g, inward,
            converted_face_tags[4])
        for j in 0:i-1
            a = node_id(i, j, 0); b = node_id(i + 1, j, 0)
            c = node_id(i, j + 1, 0); g = node_id(i + 1, j + 1, 0)
            inward = _node(coords, node_id(i, j, 1))
            position = _write_outward_triangle!(
                tris, tri_tags, position, coords, a, c, g, inward,
                converted_face_tags[4])
            position = _write_outward_triangle!(
                tris, tri_tags, position, coords, a, b, g, inward,
                converted_face_tags[4])
        end
    end
    # f5 (k = nw): SIM_12's (d,e,h) on every cell plus SIM_9's (d,f,h) on
    # strictly-lower cells.
    @inbounds for i in 0:nr-1
        d = node_id(i, i, nw); e = node_id(i + 1, i, nw)
        h = node_id(i + 1, i + 1, nw)
        inward = _node(coords, node_id(i + 1, i + 1, nw - 1))
        position = _write_outward_triangle!(
            tris, tri_tags, position, coords, d, e, h, inward,
            converted_face_tags[5])
        for j in 0:i-1
            d = node_id(i, j, nw); e = node_id(i + 1, j, nw)
            f = node_id(i, j + 1, nw); h = node_id(i + 1, j + 1, nw)
            inward = _node(coords, node_id(i + 1, j + 1, nw - 1))
            position = _write_outward_triangle!(
                tris, tri_tags, position, coords, d, f, h, inward,
                converted_face_tags[5])
            position = _write_outward_triangle!(
                tris, tri_tags, position, coords, d, e, h, inward,
                converted_face_tags[5])
        end
    end
    position == triangle_count + 1 || throw(ErrorException(
        "$_CALLER: internal boundary triangle count invariant failed"))

    extracted_boundary, maximum_incidence = boundary_faces(tets)
    maximum_incidence == 2 || throw(ErrorException(
        "$_CALLER: constructed tet mesh has face incidence $maximum_incidence"))
    sort!(extracted_boundary)
    extracted_boundary == _canonical_triangles(tris) || throw(ErrorException(
        "$_CALLER: emitted boundary triangles do not match the tet boundary"))
    affine && _certify_volume(coords, tets, converted_corners)

    mesh = Mesh(coords; tris=tris, tets=tets, tri_tag=tri_tags,
                tet_tag=fill(converted_volume_tag, tet_count))
    diagnostic = validate(mesh)
    diagnostic.ok || _throw_simplex_validation(_CALLER, diagnostic.messages)
    (nnodes(mesh), ntris(mesh), ntets(mesh)) ==
        (node_count, triangle_count, tet_count) || throw(ErrorException(
        "$_CALLER: finalized mesh count invariant failed"))
    return mesh
end

end
