"""
    TransfiniteVolume

Bounded Gmsh-4.15.2-compatible transfinite volume meshing for a six-faced,
eight-corner block.  The result is a first-order simplex `Mesh`: each logical
hexahedron is split into the six tetrahedra used by Gmsh when all six
transfinite boundary faces are unrecombined.

Two input forms exist. With `faces=nothing` the block must be an affine
parallelepiped certified from the eight corners, and the interior follows the
exact dyadic/trilinear affine map. With `faces` the caller supplies the six
canonical-oriented boundary face grids (Gmsh's `(vmin, umax, vmax, umin,
wmin, wmax)` slot order); the interior then follows Gmsh's `transfiniteHex`
Coons interpolation over those grids, so warped and otherwise non-affine
blocks are admitted whenever their faces mesh.

This module deliberately does not claim support for five-faced/prismatic
volumes, recombined hexahedra/prisms, QuadTri, holes, multiple blocks,
periodic seams, or high-order elements. The `faces=nothing` path uses uniform
`Progression 1` spacing; nonuniform laws arrive through the supplied face
grids on the `faces` path. As required by the finalized `Mesh` contract,
represented boundary areas and tetrahedron volumes must remain finite Float64
values; finite input coordinates alone do not imply finite derived measures.
"""
module TransfiniteVolume

using ..MeshTypes: Mesh, boundary_faces, nnodes, ntris, ntets, validate,
                   _throw_simplex_validation
using ..Predicates: orient3
using ..StructuredNumerics: _needs_exact_affine, _affine_basis3,
                            _affine_point3, _uses_exact_affine,
                            _certify_tet_volume
using ..StructuredRecombine: _canonical_tets
import ..StructuredRecombine

export mesh_transfinite_volume

const _DEFAULT_MAX_NODES = 10_000_000
const _DEFAULT_MAX_TETS = 60_000_000
const _DEFAULT_MAX_BOUNDARY_TRIANGLES = 20_000_000
const _AFFINE_TOLERANCE = 4096eps(Float64)

@inline function _checked_add(a::Int, b::Int, what::AbstractString)
    try
        return Base.checked_add(a, b)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError("mesh_transfinite_volume: $what count overflows Int"))
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
            throw(ArgumentError("mesh_transfinite_volume: $what count overflows Int"))
        end
    end
    return result
end

function _limit(value, name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "mesh_transfinite_volume: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "mesh_transfinite_volume: $name must not be Bool"))
    value >= 0 || throw(ArgumentError(
        "mesh_transfinite_volume: $name must be non-negative"))
    value <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $name exceeds the Int32 topology limit"))
    return Int(value)
end

function _count(value, axis::AbstractString)
    value isa Integer || throw(ArgumentError(
        "mesh_transfinite_volume: $axis cell count must be an integer"))
    value isa Bool && throw(ArgumentError(
        "mesh_transfinite_volume: $axis cell count must not be Bool"))
    value > 0 || throw(ArgumentError(
        "mesh_transfinite_volume: $axis cell count must be positive"))
    value <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $axis cell count exceeds the Int32 topology limit"))
    return Int(value)
end

function _three_counts(raw)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: cells must be an indexable collection of three counts"))
    end
    count == 3 || throw(ArgumentError(
        "mesh_transfinite_volume: cells must contain exactly three counts"))
    values = Vector{Int}(undef, 3)
    cursor = 1
    try
        for value in raw
            cursor <= 3 || throw(ArgumentError(
                "mesh_transfinite_volume: cells iteration produced more than three counts"))
            values[cursor] = _count(value, ("u", "v", "w")[cursor])
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: could not read cells: $(sprint(showerror, err))"))
    end
    cursor == 4 || throw(ArgumentError(
        "mesh_transfinite_volume: cells iteration ended before three counts"))
    return values[1], values[2], values[3]
end

function _tag(value, name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "mesh_transfinite_volume: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "mesh_transfinite_volume: $name must not be Bool"))
    0 <= value <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $name must lie in 0:$(typemax(Int32))"))
    return Int32(value)
end

function _face_tags(raw)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: face_tags must be an indexable collection"))
    end
    count == 6 || throw(ArgumentError(
        "mesh_transfinite_volume: face_tags must contain exactly six tags"))
    tags = Vector{Int32}(undef, 6)
    cursor = 1
    try
        for value in raw
            cursor <= 6 || throw(ArgumentError(
                "mesh_transfinite_volume: face_tags iteration produced more than six tags"))
            tags[cursor] = _tag(value, "face tag $cursor")
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: could not read face_tags: $(sprint(showerror, err))"))
    end
    cursor == 7 || throw(ArgumentError(
        "mesh_transfinite_volume: face_tags iteration ended before six tags"))
    return tags
end

function _point3(raw, index::Int)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: corner $index is not indexable"))
    end
    count == 3 || throw(ArgumentError(
        "mesh_transfinite_volume: corner $index must have exactly three coordinates"))
    values = try
        (raw[1], raw[2], raw[3])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: corner $index coordinates must be " *
            "Float64-representable: $(sprint(showerror, err))"))
    end
    any(value -> value isa Bool, values) && throw(ArgumentError(
        "mesh_transfinite_volume: corner $index coordinates must not be Bool"))
    point = try
        (Float64(values[1]), Float64(values[2]), Float64(values[3]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: corner $index coordinates must be " *
            "Float64-representable: $(sprint(showerror, err))"))
    end
    all(isfinite, point) || throw(ArgumentError(
        "mesh_transfinite_volume: corner $index has a non-finite coordinate"))
    return point
end

function _eight_corners(raw)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: corners must be an indexable collection"))
    end
    count == 8 || throw(ArgumentError(
        "mesh_transfinite_volume: corners must contain exactly eight points"))
    corners = Vector{NTuple{3,Float64}}(undef, 8)
    cursor = 1
    try
        for value in raw
            cursor <= 8 || throw(ArgumentError(
                "mesh_transfinite_volume: corners iteration produced more than eight points"))
            corners[cursor] = _point3(value, cursor)
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: could not read corners: $(sprint(showerror, err))"))
    end
    cursor == 9 || throw(ArgumentError(
        "mesh_transfinite_volume: corners iteration ended before eight points"))
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
    # This rare certificate operates on the exact dyadic values represented by
    # the input Float64 coordinates. It avoids rejecting an exactly affine
    # block when normalization and independent rounded sums disagree by an ULP.
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
    @inbounds for i in 1:8
        delta = _sub3(corners[i], origin)
        all(isfinite, delta) || throw(ArgumentError(
            "mesh_transfinite_volume: corner span overflows Float64 at corner $i"))
        deltas[i] = delta
        scale = max(scale, _maxabs3(delta))
    end
    (isfinite(scale) && scale > 0) || throw(ArgumentError(
        "mesh_transfinite_volume: corners are geometrically degenerate"))
    normalized = NTuple{3,Float64}[
        (delta[1] / scale, delta[2] / scale, delta[3] / scale) for delta in deltas]
    all(point -> all(isfinite, point), normalized) || throw(ArgumentError(
        "mesh_transfinite_volume: normalized corner coordinates are not finite"))

    u = normalized[2]
    v = normalized[4]
    w = normalized[5]
    # Certify the represented input corners directly. Normalization can round
    # a very thin but finite affine basis onto a singular floating determinant.
    orientation = orient3(corners[1], corners[2], corners[4], corners[5])
    orientation < 0 || throw(ArgumentError(
        orientation == 0 ?
        "mesh_transfinite_volume: canonical u/v/w corner directions are coplanar" :
        "mesh_transfinite_volume: corners must use the positive canonical " *
        "Gmsh order (s0,s1,s2,s3,s4,s5,s6,s7)"))

    # For the normalized 3x3 basis A, |det(A)| / ||A||_F^2 is a
    # conservative lower bound for its smallest singular value. Scaling the
    # affine residual tolerance by that bound prevents a fixed absolute
    # tolerance from admitting a fold in an arbitrarily thin block.
    determinant = abs(_dot3(u, _cross3(v, w)))
    frobenius_squared = _dot3(u, u) + _dot3(v, v) + _dot3(w, w)
    conditioning_bound = determinant / frobenius_squared
    if !(isfinite(conditioning_bound) && conditioning_bound > 0)
        _exact_affine_corners(corners) && return nothing
        throw(ArgumentError(
            "mesh_transfinite_volume: affine corner conditioning is not " *
            "representable in Float64"))
    end
    affine_tolerance = _AFFINE_TOLERANCE * min(1.0, conditioning_bound)

    expected = (_add3(u, v), _add3(u, w), _add3(_add3(u, v), w), _add3(v, w))
    indices = (3, 6, 7, 8)
    maximum_error = 0.0
    @inbounds for position in 1:4
        actual = normalized[indices[position]]
        target = expected[position]
        error = _maxabs3(_sub3(actual, target))
        isfinite(error) || throw(ArgumentError(
            "mesh_transfinite_volume: affine residual is not finite"))
        maximum_error = max(maximum_error, error)
    end
    if maximum_error > affine_tolerance && !_exact_affine_corners(corners)
        throw(ArgumentError(
            "mesh_transfinite_volume: corners do not form an affine parallelepiped " *
            "(normalized residual $maximum_error exceeds the conditioning-scaled " *
            "tolerance $affine_tolerance)"))
    end
    return nothing
end

@inline _lerp(a::Float64, b::Float64, t::Float64) = (1 - t) * a + t * b

# --- warped (face-input) transfinite volume path ------------------------------
#
# Gmsh's `MeshTransfiniteVolume` interpolates the volume from the six boundary
# *face meshes*, not only from the corners: each interior tab node follows the
# scalar `transfiniteHex` blend of the six face interpolants minus the twelve
# edge interpolants plus the trilinear corner term, and the u/v/w parameters
# are cumulative chord-length ratios measured along the s0s1, s1s2, and s1s5
# edge chains. Boundary tab nodes reuse the face-mesh vertices bitwise, so a
# genuinely non-affine block (warped ruled faces, curved edges) keeps its
# surface discretization instead of collapsing onto the trilinear corner map.
# The emitted boundary is the canonical conforming split induced by the six-
# tet cell subdivision (as in the affine path); the face meshes' own
# triangulations are validated but act as the transfinite-mesh certificate,
# not the emitted boundary. The
# callers hand the kernel six canonical-oriented grids in the face-tag order
# `(vmin, umax, vmax, umin, wmin, wmax)`: every grid is a `3×(N1+1)·(N2+1)`
# matrix whose column `p1 + p2*(N1+1)` is the canonical parametric node — the
# u/v/w axes of each face slot are `(u,w)`, `(v,w)`, `(u,w)`, `(v,w)`,
# `(u,v)`, `(u,v)` respectively.

@inline function _warped_slot_dims(slot::Int, nu::Int, nv::Int, nw::Int)
    slot == 1 && return (nu, nw)
    slot == 2 && return (nv, nw)
    slot == 3 && return (nu, nw)
    slot == 4 && return (nv, nw)
    slot == 5 && return (nu, nv)
    slot == 6 && return (nu, nv)
    throw(ArgumentError(
        "mesh_transfinite_volume: face slot must lie in 1:6; got $slot"))
end

@inline _face_node(points, p1::Int, p2::Int, n1::Int) =
    (points[1, p1 + 1 + p2 * n1], points[2, p1 + 1 + p2 * n1],
     points[3, p1 + 1 + p2 * n1])

# Six converted face records in canonical slot order. Each record carries the
# point matrix, the face mesh's own triangles (indices into the same canonical
# layout), and the boundary-element tags.
struct _WarpedFace
    points::Matrix{Float64}
    tris::Matrix{Int32}
    tags::Vector{Int32}
end

function _warped_face(raw, slot::Int, nu::Int, nv::Int, nw::Int)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: faces[$slot] is not indexable"))
    end
    count == 3 || throw(ArgumentError(
        "mesh_transfinite_volume: faces[$slot] must be a " *
        "(points, triangles, tags) triple"))
    n1, n2 = _warped_slot_dims(slot, nu, nv, nw)
    expected = _checked_mul("faces[$slot] node", n1 + 1, n2 + 1)

    points = try
        raw[1]
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: could not read faces[$slot] points"))
    end
    points isa AbstractMatrix || throw(ArgumentError(
        "mesh_transfinite_volume: faces[$slot] points must be a matrix"))
    size(points, 1) == 3 || throw(ArgumentError(
        "mesh_transfinite_volume: faces[$slot] points must have three rows"))
    size(points, 2) == expected || throw(ArgumentError(
        "mesh_transfinite_volume: faces[$slot] carries $(size(points, 2)) " *
        "nodes but the slot requires $expected"))
    converted = Matrix{Float64}(undef, 3, expected)
    @inbounds for column in 1:expected, row in 1:3
        value = Float64(points[row, column])
        isfinite(value) || throw(ArgumentError(
            "mesh_transfinite_volume: faces[$slot] node $column is not finite"))
        converted[row, column] = value
    end

    tris = try
        raw[2]
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: could not read faces[$slot] triangles"))
    end
    tris isa AbstractMatrix || throw(ArgumentError(
        "mesh_transfinite_volume: faces[$slot] triangles must be a matrix"))
    size(tris, 1) == 3 || throw(ArgumentError(
        "mesh_transfinite_volume: faces[$slot] triangles must have three rows"))
    ntri = size(tris, 2)
    expected_tris = _checked_mul("faces[$slot] triangle", 2, n1, n2)
    ntri == expected_tris || throw(ArgumentError(
        "mesh_transfinite_volume: faces[$slot] carries $ntri triangles but " *
        "the slot requires $expected_tris"))
    converted_tris = Matrix{Int32}(undef, 3, ntri)
    @inbounds for column in 1:ntri, row in 1:3
        index = tris[row, column]
        index isa Integer || throw(ArgumentError(
            "mesh_transfinite_volume: faces[$slot] triangle $column vertex " *
            "$row is not an integer"))
        1 <= index <= expected || throw(ArgumentError(
            "mesh_transfinite_volume: faces[$slot] triangle $column vertex " *
            "$row indexes node $index outside 1:$expected"))
        converted_tris[row, column] = Int32(index)
    end

    tags = try
        raw[3]
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: could not read faces[$slot] tags"))
    end
    length(tags) == ntri || throw(ArgumentError(
        "mesh_transfinite_volume: faces[$slot] carries $(length(tags)) tags " *
        "for $ntri triangles"))
    converted_tags = Vector{Int32}(undef, ntri)
    @inbounds for index in 1:ntri
        converted_tags[index] = _tag(tags[index], "faces[$slot] tag $index")
    end
    return _WarpedFace(converted, converted_tris, converted_tags)
end

function _warped_faces(raw, nu::Int, nv::Int, nw::Int)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: faces is not indexable"))
    end
    count == 6 || throw(ArgumentError(
        "mesh_transfinite_volume: faces must contain exactly six records"))
    faces = Vector{_WarpedFace}(undef, 6)
    cursor = 1
    try
        for entry in raw
            cursor <= 6 || throw(ArgumentError(
                "mesh_transfinite_volume: faces produced more than six records"))
            faces[cursor] = _warped_face(entry, cursor, nu, nv, nw)
            cursor += 1
        end
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_volume: could not read faces: " *
            "$(sprint(showerror, err))"))
    end
    cursor == 7 || throw(ArgumentError(
        "mesh_transfinite_volume: faces ended before six records"))
    return faces
end

# The twelve boundary edges shared between adjacent face slots. Each entry is
# (face_a, p1a, p2a_start, step_a) style data is too rigid for the mixed row/
# column accesses, so the audit lists them explicitly.
function _certify_shared_edges(faces, nu::Int, nv::Int, nw::Int)
    f1, f2, f3, f4, f5, f6 = faces
    npu, npv, npw = nu + 1, nv + 1, nw + 1
    function check(name, a, b)
        a == b || throw(ArgumentError(
            "mesh_transfinite_volume: boundary faces disagree on shared " *
            "edge $name ($(a) != $(b))"))
        return nothing
    end
    @inbounds begin
        for i in 0:nu
            check("s0s1", _face_node(f5.points, i, 0, npu),
                  _face_node(f1.points, i, 0, npu))
            check("s3s2", _face_node(f5.points, i, nv, npu),
                  _face_node(f3.points, i, 0, npu))
            check("s4s5", _face_node(f6.points, i, 0, npu),
                  _face_node(f1.points, i, nw, npu))
            check("s7s6", _face_node(f6.points, i, nv, npu),
                  _face_node(f3.points, i, nw, npu))
        end
        for j in 0:nv
            check("s0s3", _face_node(f5.points, 0, j, npu),
                  _face_node(f4.points, j, 0, npv))
            check("s1s2", _face_node(f5.points, nu, j, npu),
                  _face_node(f2.points, j, 0, npv))
            check("s4s7", _face_node(f6.points, 0, j, npu),
                  _face_node(f4.points, j, nw, npv))
            check("s5s6", _face_node(f6.points, nu, j, npu),
                  _face_node(f2.points, j, nw, npv))
        end
        for k in 0:nw
            check("s0s4", _face_node(f1.points, 0, k, npu),
                  _face_node(f4.points, 0, k, npv))
            check("s1s5", _face_node(f1.points, nu, k, npu),
                  _face_node(f2.points, 0, k, npv))
            check("s3s7", _face_node(f3.points, 0, k, npu),
                  _face_node(f4.points, nv, k, npv))
            check("s2s6", _face_node(f3.points, nu, k, npu),
                  _face_node(f2.points, nv, k, npv))
        end
    end
    return nothing
end

# Cumulative chord-length ratios along one edge chain: gmsh parameterizes the
# interpolation by the physical edge mesh, not by uniform logical index.
function _chord_ratios(chain, what::AbstractString)
    rows = length(chain)
    ratios = Vector{Float64}(undef, rows)
    total = 0.0
    ratios[1] = 0.0
    @inbounds for i in 2:rows
        delta = _sub3(chain[i], chain[i - 1])
        total += sqrt(_dot3(delta, delta))
        ratios[i] = total
    end
    (isfinite(total) && total > 0) || throw(ArgumentError(
        "mesh_transfinite_volume: reference edge $what is geometrically " *
        "degenerate"))
    @inbounds for i in 1:rows
        ratios[i] /= total
    end
    return ratios
end

# The scalar `transfiniteHex` interpolation from Gmsh's
# meshGRegionTransfinite.cpp, applied coordinate-wise with the same term
# order so the result tracks upstream bitwise where inputs are identical.
# `faces_uvw` = (umin, umax, vmin, vmax, wmin, wmax) interpolants at the
# canonical face params. `edges_w` lists the w-parallel edges in Gmsh's
# literal order (u0v0=s0s4, u0v1=s3s7, u1v0=s1s5, u1v1=s2s6), `edges_u` the
# u-parallel edges (v0w0=s0s1, v0w1=s4s5, v1w0=s3s2, v1w1=s7s6), `edges_v` the
# v-parallel edges (u0w0=s0s3, u1w0=s1s2, u0w1=s4s7, u1w1=s5s6), and `s` the
# eight corner vertices.
@inline function _transfinite_hex(faces_uvw, edges_u, edges_v, edges_w, s,
                                  u::Float64, v::Float64, w::Float64)
    fu0, fu1 = faces_uvw[1], faces_uvw[2]
    fv0, fv1 = faces_uvw[3], faces_uvw[4]
    fw0, fw1 = faces_uvw[5], faces_uvw[6]
    return ntuple(3) do d
        (1 - u) * fu0[d] + u * fu1[d] + (1 - v) * fv0[d] + v * fv1[d] +
        (1 - w) * fw0[d] + w * fw1[d] -
        ((1 - u) * (1 - v) * edges_w[1][d] + (1 - u) * v * edges_w[2][d] +
         u * (1 - v) * edges_w[3][d] + u * v * edges_w[4][d]) -
        ((1 - v) * (1 - w) * edges_u[1][d] + (1 - v) * w * edges_u[2][d] +
         v * (1 - w) * edges_u[3][d] + v * w * edges_u[4][d]) -
        ((1 - u) * (1 - w) * edges_v[1][d] + (1 - w) * u * edges_v[2][d] +
         w * (1 - u) * edges_v[3][d] + u * w * edges_v[4][d]) +
        (1 - u) * (1 - v) * (1 - w) * s[1][d] + u * (1 - v) * (1 - w) * s[2][d] +
        u * v * (1 - w) * s[3][d] + (1 - u) * v * (1 - w) * s[4][d] +
        (1 - u) * (1 - v) * w * s[5][d] + u * (1 - v) * w * s[6][d] +
        u * v * w * s[7][d] + (1 - u) * v * w * s[8][d]
    end
end

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

# The shared recombined-emission machinery lives in `StructuredRecombine`
# (loaded before this module); these delegations keep the diagnostics
# prefixed with this kernel's public entry point.
_emit_canonical_tet!(tets, position::Int, coords,
                     a::Int32, b::Int32, c::Int32, d::Int32) =
    StructuredRecombine._emit_canonical_tet!(
        "mesh_transfinite_volume", tets, position, coords, a, b, c, d)

_write_outward_triangle!(tris, tags, position::Int, coords,
                         a::Int32, b::Int32, c::Int32, opposite, tag::Int32) =
    StructuredRecombine._write_outward_triangle!(
        "mesh_transfinite_volume", tris, tags, position, coords, a, b, c,
        opposite, tag)

_write_outward_quad!(quads, tags, position::Int, coords,
                     a::Int32, b::Int32, c::Int32, d::Int32, opposite,
                     tag::Int32) =
    StructuredRecombine._write_outward_quad!(
        quads, tags, position, coords, a, b, c, d, opposite, tag)

_emit_recombined_cell!(cells::AbstractMatrix{Int32}, column::Int,
                       shadow::AbstractMatrix{Int32}, position::Int,
                       coords, vertices, decomp) =
    StructuredRecombine._emit_recombined_cell!(
        "mesh_transfinite_volume", cells, column, shadow, position, coords,
        vertices, decomp)

_certify_recombined_boundary(shadow, tris, quads) =
    StructuredRecombine._certify_recombined_boundary(
        "mesh_transfinite_volume", shadow, tris, quads)

_recombined_mixed_mesh(coords, tris, tri_tags, quads, quad_tags, hexes,
                       prisms, volume_tag::Int32) =
    StructuredRecombine._recombined_mixed_mesh(
        "mesh_transfinite_volume", coords, tris, tri_tags, quads, quad_tags,
        hexes, prisms, volume_tag)

function _emit_positive_tet!(tets, position::Int, coords,
                             a::Int32, b::Int32, c::Int32, d::Int32)
    sign = orient3(_node(coords, a), _node(coords, b),
                   _node(coords, c), _node(coords, d))
    sign == 0 && throw(ArgumentError(
        "mesh_transfinite_volume: interpolation produced a zero-volume " *
        "tetrahedron at output position $position"))
    @inbounds begin
        tets[1, position] = a
        tets[2, position] = b
        if sign < 0
            tets[3, position] = c
            tets[4, position] = d
        else
            tets[3, position] = d
            tets[4, position] = c
        end
    end
    return nothing
end

# Gmsh's unrecombined six-tetrahedron cell subdivision
# (`CREATE_SIM_1` through `CREATE_SIM_6`); also the reference partition the
# recombined shadow decompositions are certified against.
function _emit_six_simplex_tets!(tets, coords, node_id,
                                 nu::Int, nv::Int, nw::Int)
    tet_position = 0
    @inbounds for k in 0:nw-1, j in 0:nv-1, i in 0:nu-1
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
            _emit_positive_tet!(tets, tet_position, coords, vertices...)
        end
    end
    return tet_position
end

@inline function _write_triangle!(tris, tags, position::Int,
                                  a::Int32, b::Int32, c::Int32, tag::Int32)
    @inbounds begin
        tris[1, position] = a
        tris[2, position] = b
        tris[3, position] = c
        tags[position] = tag
    end
    return position + 1
end

_canonical_triangles(tris) = StructuredRecombine._canonical_triangles(tris)

# The canonical face-slot parametric directions into the logical tab grid:
# slot s maps its (p1,p2) node to tab coordinates (i,j,k) via these tables.
# (param axis 1 index, param axis 2 index, fixed axis index, fixed value)
# expressed per slot in the order (vmin, umax, vmax, umin, wmin, wmax).
const _WARPED_SLOT_PLANE =
    ((1, 3, 2, 0),   # vmin:  (p1,p2) = (u,w),   fixed v = 0
     (2, 3, 1, :hi), # umax:  (p1,p2) = (v,w),   fixed u = nu
     (1, 3, 2, :hi), # vmax:  (p1,p2) = (u,w),   fixed v = nv
     (2, 3, 1, 0),   # umin:  (p1,p2) = (v,w),   fixed u = 0
     (1, 2, 3, 0),   # wmin:  (p1,p2) = (u,v),   fixed w = 0
     (1, 2, 3, :hi)) # wmax:  (p1,p2) = (u,v),   fixed w = nw

# The four inward tab nodes adjacent to a boundary cell: the cell at canonical
# (I,J) has its neighbors offset one layer along the slot's fixed axis.
# Order: (I,J), (I+1,J), (I,J+1), (I+1,J+1) in the slot's parametric axes.
@inline function _inward_candidates(node_id, slot::Int, I::Int, J::Int,
                                    nu::Int, nv::Int, nw::Int)
    a1, a2, fixed_axis, fixed = _WARPED_SLOT_PLANE[slot]
    inner_coord = fixed == 0 ? 1 :
        fixed_axis == 1 ? nu - 1 : fixed_axis == 2 ? nv - 1 : nw - 1
    return ntuple(4) do corner
        offset = corner - 1
        coord = (0, 0, 0)
        coord = Base.setindex(coord, I + (offset & 1), a1)
        coord = Base.setindex(coord, J + (offset >> 1), a2)
        coord = Base.setindex(coord, inner_coord, fixed_axis)
        node_id(coord[1], coord[2], coord[3])
    end
end

# Cell-letter space for the boundary fold audit: the eight corners of a cell
# are named a..h for (i,j,k),(i+1,j,k),(i,j+1,k),(i,j,k+1),(i+1,j,k+1),
# (i,j+1,k+1),(i+1,j+1,k),(i+1,j+1,k+1) — the vertex order of the fixed
# six-tet subdivision (a,b,c,d),(b,c,d,e),(d,e,c,f),(b,c,e,g),(c,f,e,g),
# (e,f,h,g). Per slot, `_WARPED_SLOT_LETTERS` maps each letter to a positive
# index into the boundary quad tabs (a,b,c,d) or a negative index into the
# `_inward_candidates` tuple; `_WARPED_SLOT_TRIS` gives the two emitted
# boundary tris as letter triples — (a,b,d),(b,e,d) on vmin say — and
# `_WARPED_SLOT_APEX_LETTER` the interior apex of the incident tet for each.
# The incident apex is the canonical inward witness: strictly on the inward
# side certifies the emitted orientation; exactly on the triangle's plane
# means the incident tet is degenerate and cannot certify.
const _WARPED_SLOT_LETTERS =
    (( 1, 2,-1, 3, 4,-3,-2,-4),   # vmin:  boundary a,b,d,e
     (-1, 1,-2,-3, 3,-4, 2, 4),   # umax:  boundary b,e,g,h
     (-1,-2, 1,-3,-4, 3, 2, 4),   # vmax:  boundary c,f,g,h
     ( 1,-1, 2, 3,-3, 4,-2,-4),   # umin:  boundary a,c,d,f
     ( 1, 2, 3,-1,-2,-3, 4,-4),   # wmin:  boundary a,b,c,g
     (-1,-2,-3, 1, 2, 3,-4, 4))   # wmax:  boundary d,e,f,h
const _WARPED_SLOT_TRIS =
    (((1, 2, 4), (2, 5, 4)),      # vmin
     ((2, 7, 5), (7, 8, 5)),      # umax
     ((3, 7, 6), (7, 8, 6)),      # vmax
     ((1, 3, 4), (3, 6, 4)),      # umin
     ((1, 2, 3), (2, 7, 3)),      # wmin
     ((4, 5, 6), (5, 8, 6)))      # wmax
const _WARPED_SLOT_APEX_LETTER =
    ((3, 3), (3, 6), (5, 5), (2, 5), (4, 5), (3, 7))

# Every tet edge of the six-tet cell subdivision, in letter indices: a real
# fold is some edge piercing an emitted boundary triangle's interior — the
# discriminating test when the corner signs straddle the supporting plane.
const _WARPED_CELL_EDGES =
    ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4),(2,5),(3,5),(4,5),(4,6),(5,6),
     (3,6),(2,7),(3,7),(5,7),(6,7),(5,8),(6,8),(7,8))

@inline function _warped_cell_corner_ids(map, tabs, cands)
    return ntuple(8) do letter
        source = map[letter]
        source > 0 ? tabs[source] : cands[-source]
    end
end

# True when the open segment (p,q) pierces the interior or boundary of
# triangle (a,b,c): the endpoints already straddle the supporting plane
# (nonzero, opposite), so the three edge-plane orientations decide whether
# the crossing point lands in the closed triangle. A grazing contact — the
# pierce point on a triangle edge — counts as a crossing too: a cell edge
# touching the emitted boundary sheet is not a certifiable configuration.
@inline function _segment_crosses_triangle(p, q, a, b, c)
    o1 = orient3(p, q, a, b)
    o2 = orient3(p, q, b, c)
    o3 = orient3(p, q, c, a)
    return (o1 >= 0 && o2 >= 0 && o3 >= 0) ||
           (o1 <= 0 && o2 <= 0 && o3 <= 0)
end

function _transfinite_volume_warped_mesh(converted_corners, raw_faces,
                                         raw_face_tags,
                                         nu::Int, nv::Int, nw::Int,
                                         npu::Int, npv::Int, npw::Int,
                                         node_count::Int, tet_count::Int,
                                         triangle_count::Int,
                                         converted_volume_tag::Int32)
    faces = _warped_faces(raw_faces, nu, nv, nw)
    _certify_shared_edges(faces, nu, nv, nw)
    converted_face_tags = _face_tags(raw_face_tags)
    node_id(i::Int, j::Int, k::Int) = Int32(i + 1 + npu * (j + npv * k))
    coords = Matrix{Float64}(undef, 3, node_count)
    tets = Matrix{Int32}(undef, 4, tet_count)
    tri_tags = Vector{Int32}(undef, triangle_count)
    tris = _mesh_transfinite_volume_warped(converted_corners, nu, nv, nw,
                                           faces, node_id, coords, tets,
                                           tri_tags, converted_face_tags)
    extracted_boundary, maximum_incidence = boundary_faces(tets)
    maximum_incidence == 2 || throw(ErrorException(
        "mesh_transfinite_volume: constructed tet mesh has face incidence " *
        "$maximum_incidence"))
    sort!(extracted_boundary)
    extracted_boundary == _canonical_triangles(tris) || throw(ErrorException(
        "mesh_transfinite_volume: emitted boundary triangles do not match " *
        "the tet boundary"))
    mesh = Mesh(coords; tris=tris, tets=tets, tri_tag=tri_tags,
                tet_tag=fill(converted_volume_tag, tet_count))
    diagnostic = validate(mesh)
    diagnostic.ok || _throw_simplex_validation(
        "mesh_transfinite_volume", diagnostic.messages)
    (nnodes(mesh), ntris(mesh), ntets(mesh)) ==
        (node_count, triangle_count, tet_count) || throw(ErrorException(
        "mesh_transfinite_volume: finalized mesh count invariant failed"))
    return mesh
end

function _mesh_transfinite_volume_warped(corners, nu::Int, nv::Int, nw::Int,
                                         faces, node_id, coords, tets,
                                         tri_tags_source, face_tags)
    _fill_warped_volume_coords!(coords, corners, faces, nu, nv, nw, node_id)
    _emit_six_simplex_tets!(tets, coords, node_id, nu, nv, nw) ==
        size(tets, 2) || throw(ErrorException(
        "mesh_transfinite_volume: internal tetrahedron count invariant failed"))
    tris = Matrix{Int32}(undef, 3, length(tri_tags_source))
    _emit_warped_volume_boundary!(tris, tri_tags_source,
                                  Matrix{Int32}(undef, 4, 0), Int32[],
                                  coords, node_id, nu, nv, nw, face_tags,
                                  nothing)
    return tris
end

# Warped-path grid fill shared by the simplex and recombined paths: the
# caller-provided corners must be the eight face-grid corners bitwise, the
# u/v/w parameters follow cumulative chord-length ratios along Gmsh's three
# reference edges, and boundary tab nodes reuse the face grids bitwise.
function _fill_warped_volume_coords!(coords, corners, faces,
                                     nu::Int, nv::Int, nw::Int, node_id)
    f1, f2, f3, f4, f5, f6 = faces
    npu, npv = nu + 1, nv + 1

    # The caller-provided corners must be the eight face-grid corners bitwise.
    _face_node(f5.points, 0, 0, npu) == corners[1] &&
        _face_node(f5.points, nu, 0, npu) == corners[2] &&
        _face_node(f5.points, nu, nv, npu) == corners[3] &&
        _face_node(f5.points, 0, nv, npu) == corners[4] &&
        _face_node(f6.points, 0, 0, npu) == corners[5] &&
        _face_node(f6.points, nu, 0, npu) == corners[6] &&
        _face_node(f6.points, nu, nv, npu) == corners[7] &&
        _face_node(f6.points, 0, nv, npu) == corners[8] ||
        throw(ArgumentError(
            "mesh_transfinite_volume: corners do not match the face grids"))
    orientation = orient3(corners[1], corners[2], corners[4], corners[5])
    orientation < 0 || throw(ArgumentError(
        orientation == 0 ?
        "mesh_transfinite_volume: canonical u/v/w corner directions are " *
        "coplanar" :
        "mesh_transfinite_volume: corners must use the positive canonical " *
        "Gmsh order (s0,s1,s2,s3,s4,s5,s6,s7)"))

    # Chord-length parameterization along Gmsh's three reference edges.
    us = _chord_ratios(
        [_face_node(f5.points, i, 0, npu) for i in 0:nu], "s0s1")
    vs = _chord_ratios(
        [_face_node(f2.points, j, 0, npv) for j in 0:nv], "s1s2")
    ws = _chord_ratios(
        [_face_node(f2.points, 0, k, npv) for k in 0:nw], "s1s5")

    # Boundary tab nodes reuse the face grids bitwise (shared edges were
    # certified identical, so any owning slot writes the same value).
    @inbounds for k in 0:nw, j in 0:nv, i in 0:nu
        boundary = i == 0 || i == nu || j == 0 || j == nv || k == 0 || k == nw
        if boundary
            point = i == 0 ? _face_node(f4.points, j, k, npv) :
                    i == nu ? _face_node(f2.points, j, k, npv) :
                    j == 0 ? _face_node(f1.points, i, k, npu) :
                    j == nv ? _face_node(f3.points, i, k, npu) :
                    k == 0 ? _face_node(f5.points, i, j, npu) :
                    _face_node(f6.points, i, j, npu)
        else
            u = us[i + 1]; v = vs[j + 1]; w = ws[k + 1]
            point = _transfinite_hex(
                (_face_node(f4.points, j, k, npv),
                 _face_node(f2.points, j, k, npv),
                 _face_node(f1.points, i, k, npu),
                 _face_node(f3.points, i, k, npu),
                 _face_node(f5.points, i, j, npu),
                 _face_node(f6.points, i, j, npu)),
                (_face_node(f5.points, i, 0, npu),
                 _face_node(f6.points, i, 0, npu),
                 _face_node(f5.points, i, nv, npu),
                 _face_node(f6.points, i, nv, npu)),
                (_face_node(f5.points, 0, j, npu),
                 _face_node(f5.points, nu, j, npu),
                 _face_node(f6.points, 0, j, npu),
                 _face_node(f6.points, nu, j, npu)),
                (_face_node(f1.points, 0, k, npu),
                 _face_node(f3.points, 0, k, npu),
                 _face_node(f1.points, nu, k, npu),
                 _face_node(f3.points, nu, k, npu)),
                corners, u, v, w)
            all(isfinite, point) || throw(ArgumentError(
                "mesh_transfinite_volume: interpolation produced a " *
                "non-finite coordinate at logical node ($i,$j,$k)"))
        end
        index = Int(node_id(i, j, k))
        coords[1, index] = point[1]
        coords[2, index] = point[2]
        coords[3, index] = point[3]
    end
    return nothing
end

# Audits one canonical boundary half-triangle of cell (p1,p2) in face `slot`:
# signs of the cell's eight corners against the triangle's supporting plane,
# then a pierce test for every straddling cell edge. Returns the orientation
# sign of the incident tet's interior apex (nonzero on success).
function _warped_boundary_half_sign(coords, ids, ltri, apex_letter::Int,
                                    p1::Int, p2::Int, slot::Int)
    pa = _node(coords, ids[ltri[1]])
    pb = _node(coords, ids[ltri[2]])
    pc = _node(coords, ids[ltri[3]])
    # Signs of the cell's eight corners against the triangle's supporting
    # plane; a straddling cell edge is then tested against the triangle's
    # interior — only a true pierce is a fold, so a tilted boundary band that
    # leaves a non-incident corner on the outward side stays legal (Gmsh
    # parity).
    signs = ntuple(8) do letter
        letter == ltri[1] || letter == ltri[2] || letter == ltri[3] ? 0 :
            orient3(pa, pb, pc, _node(coords, ids[letter]))
    end
    for (x, y) in _WARPED_CELL_EDGES
        (signs[x] == 0 || signs[y] == 0 || signs[x] == signs[y]) && continue
        _segment_crosses_triangle(_node(coords, ids[x]),
                                  _node(coords, ids[y]),
                                  pa, pb, pc) && throw(ArgumentError(
            "mesh_transfinite_volume: boundary cell ($p1,$p2) of face " *
            "slot $slot folds over its inward neighbors"))
    end
    sign = signs[apex_letter]
    sign == 0 && throw(ArgumentError(
        "mesh_transfinite_volume: boundary cell ($p1,$p2) of face slot " *
        "$slot cannot certify its orientation (inward reference degenerates)"))
    return sign
end

# Boundary emission for the warped path. With `mask === nothing` every face
# cell emits the canonical tet-boundary split (A,B,C),(B,D,C) — one fixed
# diagonal per face cell in every slot's canonical params — into `tris`.
# With a six-slot `mask`, recombined face cells emit one quadrangle wound
# like the audited first canonical half (its reverse when that half needed
# flipping) into `quads` while unrecombined face cells keep the triangle
# split. Tessella's volume meshes carry the conforming boundary, which
# `_consume_volume_attributes` also re-derives downstream; each surface's own
# arrangement is a property of its surface record, not of this mesh. Every
# emitted cell is certified strictly outward against the inward-adjacent tab
# nodes of its cell.
function _emit_warped_volume_boundary!(tris, tri_tags_source,
                                       quads, quad_tags,
                                       coords, node_id,
                                       nu::Int, nv::Int, nw::Int,
                                       face_tags, mask)
    tri_position = 1
    quad_position = 1
    @inbounds for slot in 1:6
        a1, a2, fixed_axis, fixed = _WARPED_SLOT_PLANE[slot]
        n1, n2 = _warped_slot_dims(slot, nu, nv, nw)
        fixed_coord = fixed == 0 ? 0 :
            fixed_axis == 1 ? nu : fixed_axis == 2 ? nv : nw
        tag = face_tags[slot]
        function tab_node(p1::Int, p2::Int)
            coord = (0, 0, 0)
            coord = Base.setindex(coord, p1, a1)
            coord = Base.setindex(coord, p2, a2)
            coord = Base.setindex(coord, fixed_coord, fixed_axis)
            return node_id(coord[1], coord[2], coord[3])
        end
        letter_map = _WARPED_SLOT_LETTERS[slot]
        tri_letters = _WARPED_SLOT_TRIS[slot]
        apex_letters = _WARPED_SLOT_APEX_LETTER[slot]
        recombined = mask !== nothing && mask[slot]
        for p2 in 0:n2-1, p1 in 0:n1-1
            a = tab_node(p1, p2)
            b = tab_node(p1 + 1, p2)
            c = tab_node(p1, p2 + 1)
            d = tab_node(p1 + 1, p2 + 1)
            ids = _warped_cell_corner_ids(
                letter_map, (a, b, c, d),
                _inward_candidates(node_id, slot, p1, p2, nu, nv, nw))
            if recombined
                # The canonical halves share the cyclic quad (a,b,d,c);
                # half 1 (a,b,c) governs the emitted winding.
                sign = _warped_boundary_half_sign(
                    coords, ids, tri_letters[1], apex_letters[1], p1, p2, slot)
                _warped_boundary_half_sign(
                    coords, ids, tri_letters[2], apex_letters[2], p1, p2, slot)
                @inbounds begin
                    quads[1, quad_position] = a
                    if sign > 0
                        quads[2, quad_position] = b
                        quads[3, quad_position] = d
                        quads[4, quad_position] = c
                    else
                        quads[2, quad_position] = c
                        quads[3, quad_position] = d
                        quads[4, quad_position] = b
                    end
                    quad_tags[quad_position] = tag
                end
                quad_position += 1
            else
                for half in 1:2
                    ltri = tri_letters[half]
                    sign = _warped_boundary_half_sign(
                        coords, ids, ltri, apex_letters[half], p1, p2, slot)
                    if sign < 0
                        tris[1, tri_position] = ids[ltri[2]]
                        tris[2, tri_position] = ids[ltri[1]]
                        tris[3, tri_position] = ids[ltri[3]]
                    else
                        tris[1, tri_position] = ids[ltri[1]]
                        tris[2, tri_position] = ids[ltri[2]]
                        tris[3, tri_position] = ids[ltri[3]]
                    end
                    tri_tags_source[tri_position] = tag
                    tri_position += 1
                end
            end
        end
    end
    tri_position == size(tris, 2) + 1 || throw(ErrorException(
        "mesh_transfinite_volume: internal boundary triangle count " *
        "invariant failed"))
    quad_position == size(quads, 2) + 1 || throw(ErrorException(
        "mesh_transfinite_volume: internal boundary quadrangle count " *
        "invariant failed"))
    return nothing
end

"""
    mesh_transfinite_volume(corners, cells=(1,1,1);
                            volume_tag=0,
                            face_tags=(0,0,0,0,0,0),
                            faces=nothing,
                            recombine=nothing,
                            max_nodes=10_000_000,
                            max_tets=60_000_000,
                            max_boundary_triangles=20_000_000) -> Mesh | MixedMesh

Mesh a six-face block with the Gmsh 4.15.2 transfinite volume subdivision.
`corners` must contain eight finite 3-D points in Gmsh's canonical order
`(s0,s1,s2,s3,s4,s5,s6,s7)`: the first four wind around the `w=0` face and
the final four are their `w=1` counterparts. `cells=(nu,nv,nw)` gives
positive logical-cell counts; curve laws are uniformly spaced
(`Progression 1`).

By default (`faces=nothing`) the block must be affine: the four derived
corners must agree with an affine parallelepiped to a normalized,
conditioning-scaled `4096eps(Float64)` tolerance, and every interior node
follows the exact dyadic or trilinear affine map. Passing `faces` lifts that
restriction: each entry is a `(points, triangles, tags)` triple for one
boundary face slot in the `face_tags` order `(vmin, umax, vmax, umin, wmin,
wmax)`, where `points` is a `3×(N1+1)(N2+1)` matrix in the slot's canonical
parametric order — `(u,w)` for the v-faces, `(v,w)` for the u-faces, `(u,v)`
for the w-faces — and `triangles`/`tags` carry that face mesh's own
triangulation and element tags as the transfinite-mesh certificate (the
triangulation is validated but the emitted boundary is always canonical,
below). The interior then follows Gmsh's `transfiniteHex` Coons interpolation
over the six face interpolants minus the twelve edge interpolants plus the
trilinear corner term, parameterized by chord-length ratios along the `s0s1`,
`s1s2`, and `s1s5` edge chains. Boundary tab nodes reuse the face vertices
bitwise — shared edges must agree exactly or the call throws — and `corners`
must be the eight face-grid corners bitwise in positive canonical order.

With `recombine=nothing` (or `false`, or an all-false mask) the returned
simplex `Mesh` uses the exact six-tetrahedron connectivity pattern from
Gmsh's `CREATE_SIM_1` through `CREATE_SIM_6`, and the boundary triangles are
the conforming split that subdivision induces: one fixed diagonal per face
cell in every slot, certified strictly outward against the inward-adjacent
interior nodes. In the affine path, exact dyadic affine interpolation
protects remote, narrow blocks whose nested Float64 interpolation would lose
material, and a compensated exponent-scaled determinant audit certifies
conservation of the corner-defined volume. `face_tags` follow Gmsh's
canonical face order `(vmin, umax, vmax, umin, wmin, wmax)`; `volume_tag`
labels every tet.

`recombine` selects Gmsh's six-face recombination decision tree and returns
a `MixedMesh`: `true` (or all-true mask) emits one hexahedron per cell with
quadrangles on every face; a mask with exactly one opposite face pair
unrecombined emits the corresponding prism pair per cell —
`(false,true,false,true,true,true)` and `(true,false,true,false,true,true)`
leave the v- or u-face pair triangular with prisms spanning that axis, and
`(true,true,true,true,false,false)` emits Gmsh's `CREATE_PRISM_1`/
`CREATE_PRISM_2` pair with quadrangles on the four side faces. Every
emitted volume cell carries a positional tetrahedral shadow decomposition
certified strictly positively oriented; when the shadow reproduces the
unrecombined partition it is checked tet-for-tet against the reference
simplex subdivision, and the emitted boundary is audited against the
shadow's exterior faces for exact coverage. Any other partial mask throws —
matching Gmsh's "wrong surface recombination in transfinite volume"
rejection. `max_tets` bounds the decomposed shadow-tetrahedron count and
`max_boundary_triangles` bounds the emitted boundary cells (triangles plus
quadrangles).

Unsupported here: five-face degeneracies (`mesh_transfinite_prism`),
nonuniform curve laws on the `faces=nothing` path (nonuniform spacing
arrives through `faces` grids), QuadTri, holes, multiple blocks, periodic
seams, high-order elements, and coordinate scales whose derived boundary
areas or tetrahedron volumes are not finite Float64 values.
"""
function mesh_transfinite_volume(corners, cells=(1, 1, 1);
                                 volume_tag=0,
                                 face_tags=(0, 0, 0, 0, 0, 0),
                                 faces=nothing,
                                 recombine=nothing,
                                 max_nodes=_DEFAULT_MAX_NODES,
                                 max_tets=_DEFAULT_MAX_TETS,
                                 max_boundary_triangles=
                                     _DEFAULT_MAX_BOUNDARY_TRIANGLES)
    node_limit = _limit(max_nodes, "max_nodes")
    tet_limit = _limit(max_tets, "max_tets")
    triangle_limit = _limit(max_boundary_triangles, "max_boundary_triangles")
    nu, nv, nw = _three_counts(cells)

    npu = _checked_add(nu, 1, "node")
    npv = _checked_add(nv, 1, "node")
    npw = _checked_add(nw, 1, "node")
    node_count = _checked_mul("node", npu, npv, npw)
    tet_count = _checked_mul("tetrahedron", 6, nu, nv, nw)
    uv = _checked_mul("boundary triangle", nu, nv)
    uw = _checked_mul("boundary triangle", nu, nw)
    vw = _checked_mul("boundary triangle", nv, nw)

    node_count <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $node_count nodes exceed Int32 indexing"))
    tet_count <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $tet_count tetrahedra exceed the Int32 topology limit"))
    node_count <= node_limit || throw(ArgumentError(
        "mesh_transfinite_volume: $node_count nodes exceed max_nodes=$node_limit"))
    tet_count <= tet_limit || throw(ArgumentError(
        "mesh_transfinite_volume: $tet_count tetrahedra exceed max_tets=$tet_limit"))

    converted_corners = _eight_corners(corners)
    mask = _recombine_volume_mask(recombine)
    mask === nothing || return _transfinite_volume_recombined(
        converted_corners, faces, face_tags, mask, nu, nv, nw,
        npu, npv, npw, node_count, tet_count, uv, uw, vw,
        node_limit, tet_limit, triangle_limit,
        _tag(volume_tag, "volume_tag"))

    face_cells = _checked_add(_checked_add(uv, uw, "boundary triangle"),
                              vw, "boundary triangle")
    triangle_count = _checked_mul("boundary triangle", 4, face_cells)
    triangle_count <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $triangle_count boundary triangles exceed " *
        "the Int32 topology limit"))
    triangle_count <= triangle_limit || throw(ArgumentError(
        "mesh_transfinite_volume: $triangle_count boundary triangles exceed " *
        "max_boundary_triangles=$triangle_limit"))

    faces === nothing || return _transfinite_volume_warped_mesh(
        converted_corners, faces, face_tags, nu, nv, nw,
        npu, npv, npw, node_count, tet_count, triangle_count,
        _tag(volume_tag, "volume_tag"))
    _certify_affine(converted_corners)
    converted_face_tags = _face_tags(face_tags)
    converted_volume_tag = _tag(volume_tag, "volume_tag")

    node_id(i::Int, j::Int, k::Int) = Int32(i + 1 + npu * (j + npv * k))
    coords = Matrix{Float64}(undef, 3, node_count)
    _fill_affine_volume_coords!(coords, converted_corners, nu, nv, nw,
                                node_id)

    tets = Matrix{Int32}(undef, 4, tet_count)
    _emit_six_simplex_tets!(tets, coords, node_id, nu, nv, nw) ==
        tet_count || throw(ErrorException(
        "mesh_transfinite_volume: internal tetrahedron count invariant failed"))

    tris = Matrix{Int32}(undef, 3, triangle_count)
    tri_tags = Vector{Int32}(undef, triangle_count)
    position = 1
    # Face 0: vmin.
    @inbounds for k in 0:nw-1, i in 0:nu-1
        a=node_id(i,0,k); b=node_id(i+1,0,k)
        d=node_id(i,0,k+1); e=node_id(i+1,0,k+1)
        position=_write_triangle!(tris,tri_tags,position,a,b,d,converted_face_tags[1])
        position=_write_triangle!(tris,tri_tags,position,b,e,d,converted_face_tags[1])
    end
    # Face 1: umax.
    @inbounds for k in 0:nw-1, j in 0:nv-1
        b=node_id(nu,j,k); g=node_id(nu,j+1,k)
        e=node_id(nu,j,k+1); h=node_id(nu,j+1,k+1)
        position=_write_triangle!(tris,tri_tags,position,b,g,e,converted_face_tags[2])
        position=_write_triangle!(tris,tri_tags,position,e,g,h,converted_face_tags[2])
    end
    # Face 2: vmax.
    @inbounds for k in 0:nw-1, i in 0:nu-1
        c=node_id(i,nv,k); g=node_id(i+1,nv,k)
        f=node_id(i,nv,k+1); h=node_id(i+1,nv,k+1)
        position=_write_triangle!(tris,tri_tags,position,c,f,g,converted_face_tags[3])
        position=_write_triangle!(tris,tri_tags,position,f,h,g,converted_face_tags[3])
    end
    # Face 3: umin.
    @inbounds for k in 0:nw-1, j in 0:nv-1
        a=node_id(0,j,k); c=node_id(0,j+1,k)
        d=node_id(0,j,k+1); f=node_id(0,j+1,k+1)
        position=_write_triangle!(tris,tri_tags,position,a,d,c,converted_face_tags[4])
        position=_write_triangle!(tris,tri_tags,position,c,d,f,converted_face_tags[4])
    end
    # Face 4: wmin.
    @inbounds for j in 0:nv-1, i in 0:nu-1
        a=node_id(i,j,0); b=node_id(i+1,j,0)
        c=node_id(i,j+1,0); g=node_id(i+1,j+1,0)
        position=_write_triangle!(tris,tri_tags,position,a,c,b,converted_face_tags[5])
        position=_write_triangle!(tris,tri_tags,position,b,c,g,converted_face_tags[5])
    end
    # Face 5: wmax.
    @inbounds for j in 0:nv-1, i in 0:nu-1
        d=node_id(i,j,nw); e=node_id(i+1,j,nw)
        f=node_id(i,j+1,nw); h=node_id(i+1,j+1,nw)
        position=_write_triangle!(tris,tri_tags,position,d,e,f,converted_face_tags[6])
        position=_write_triangle!(tris,tri_tags,position,e,h,f,converted_face_tags[6])
    end
    position == triangle_count + 1 || throw(ErrorException(
        "mesh_transfinite_volume: internal boundary triangle count invariant failed"))

    # A rounded trilinear center can land exactly on a face of a thin valid
    # block. Certify each face against a represented corner on its interior
    # side instead; exact orient3 then remains decisive at every finite scale.
    face_triangle_counts = (2uw, 2vw, 2uw, 2vw, 2uv, 2uv)
    opposite_corners = (converted_corners[4], converted_corners[1],
                        converted_corners[1], converted_corners[2],
                        converted_corners[5], converted_corners[1])
    triangle = 0
    @inbounds for face in 1:6
        opposite = opposite_corners[face]
        for _ in 1:face_triangle_counts[face]
            triangle += 1
            sign = orient3(_node(coords, tris[1, triangle]),
                           _node(coords, tris[2, triangle]),
                           _node(coords, tris[3, triangle]), opposite)
            sign > 0 || throw(ArgumentError(
                "mesh_transfinite_volume: boundary triangle $triangle is not " *
                "strictly outward-oriented"))
        end
    end
    triangle == triangle_count || throw(ErrorException(
        "mesh_transfinite_volume: boundary orientation count invariant failed"))

    extracted_boundary, maximum_incidence = boundary_faces(tets)
    maximum_incidence == 2 || throw(ErrorException(
        "mesh_transfinite_volume: constructed tet mesh has face incidence $maximum_incidence"))
    sort!(extracted_boundary)
    extracted_boundary == _canonical_triangles(tris) || throw(ErrorException(
        "mesh_transfinite_volume: emitted boundary triangles do not match the tet boundary"))
    _certify_tet_volume(
        coords, tets,
        (converted_corners[1], converted_corners[2],
         converted_corners[4], converted_corners[5]),
        6, "mesh_transfinite_volume", "affine block")

    mesh = Mesh(coords; tris=tris, tets=tets, tri_tag=tri_tags,
                tet_tag=fill(converted_volume_tag, tet_count))
    diagnostic = validate(mesh)
    diagnostic.ok || _throw_simplex_validation(
        "mesh_transfinite_volume", diagnostic.messages)
    (nnodes(mesh), ntris(mesh), ntets(mesh)) ==
        (node_count, triangle_count, tet_count) || throw(ErrorException(
        "mesh_transfinite_volume: finalized mesh count invariant failed"))
    return mesh
end

# ============================ Recombination ============================
#
# Gmsh's six-face transfinite recombination decision tree
# (meshGRegionTransfinite.cpp). The mask follows the canonical face order
# (vmin, umax, vmax, umin, wmin, wmax) — the same order as `face_tags`:
#   (T,T,T,T,T,T)  → one `CREATE_HEX` per cell;
#   (F,T,F,T,T,T)  → a prism pair spanning the v direction;
#   (T,F,T,F,T,T)  → a prism pair spanning the u direction;
#   (T,T,T,T,F,F)  → `CREATE_PRISM_1`/`CREATE_PRISM_2` spanning w;
#   all-false      → the unrecombined six-tet subdivision (simplex `Mesh`).
# Every other partial mask is rejected — Gmsh logs "Wrong surface
# recombination in transfinite volume" and fails the volume mesh; this
# kernel throws the equivalent ArgumentError.
function _recombine_volume_mask(recombine)
    recombine === nothing && return nothing
    mask = if recombine isa Bool
        ntuple(_ -> recombine, 6)
    elseif recombine isa Union{Tuple, AbstractVector}
        length(recombine) == 6 || throw(ArgumentError(
            "mesh_transfinite_volume: recombine mask must hold six entries " *
            "in canonical face order (vmin,umax,vmax,umin,wmin,wmax), got " *
            "$(length(recombine))"))
        ntuple(6) do slot
            flag = recombine[slot]
            flag isa Bool || throw(ArgumentError(
                "mesh_transfinite_volume: recombine mask entry $slot must " *
                "be Bool, got $(typeof(flag))"))
            flag
        end
    else
        throw(ArgumentError(
            "mesh_transfinite_volume: recombine must be a Bool or a " *
            "six-element Bool mask in canonical face order " *
            "(vmin,umax,vmax,umin,wmin,wmax), got $(typeof(recombine))"))
    end
    any(mask) || return nothing
    all(mask) ||
        mask == (false, true, false, true, true, true) ||
        mask == (true, false, true, false, true, true) ||
        mask == (true, true, true, true, false, false) ||
        throw(ArgumentError(
            "mesh_transfinite_volume: wrong surface recombination in " *
            "transfinite volume — Gmsh's six-face pattern admits only all " *
            "six faces recombined or a single opposite face pair left " *
            "unrecombined"))
    return mask
end

# Positional simplex decompositions of the emitted volume cells, indexed in
# each cell tuple's own vertex order. `_SHADOW_HEX` and the PRISM_1/PRISM_2
# pair reproduce the unrecombined six-tet partition exactly; the R variants
# decompose the v-free/u-free prism pair, indexed in the emitted vertex
# order — which is Gmsh's literal macro order with positions 1<->2 and
# 4<->5 swapped, matching the orientation fixup Gmsh applies to the
# otherwise negative-oriented partial-mask prisms.
const _SHADOW_HEX = ((1, 2, 4, 5), (2, 4, 5, 6), (5, 6, 4, 8),
                     (2, 4, 6, 3), (4, 8, 6, 3), (6, 8, 7, 3))
const _SHADOW_PRISM_1 = ((1, 2, 3, 4), (2, 3, 4, 5), (4, 5, 3, 6))
const _SHADOW_PRISM_2 = ((3, 2, 6, 1), (2, 5, 6, 1), (6, 5, 4, 1))
const _SHADOW_PRISM_1R = ((2, 1, 5, 3), (1, 3, 4, 5), (5, 4, 6, 3))
const _SHADOW_PRISM_2R = ((3, 1, 2, 6), (1, 4, 2, 6), (6, 4, 2, 5))

# Affine-path grid fill shared by the simplex and recombined paths.
function _fill_affine_volume_coords!(coords, converted_corners,
                                     nu::Int, nv::Int, nw::Int, node_id)
    affine_points = (converted_corners[1], converted_corners[2],
                     converted_corners[4], converted_corners[5])
    exact_interpolation = _needs_exact_affine(
        affine_points..., (nu, nv, nw)) &&
        _exact_affine_corners(converted_corners)
    affine_basis = _affine_basis3(affine_points..., exact_interpolation)
    @inbounds for k in 0:nw, j in 0:nv, i in 0:nu
        u = i / nu; v = j / nv; w = k / nw
        point = _uses_exact_affine(affine_basis) ?
            _affine_point3(affine_basis, u, v, w,
                           "mesh_transfinite_volume", (i, j, k)) :
            _trilinear(converted_corners, u, v, w)
        all(isfinite, point) || throw(ArgumentError(
            "mesh_transfinite_volume: interpolation produced a non-finite " *
            "coordinate at logical node ($i,$j,$k)"))
        index = Int(node_id(i, j, k))
        coords[1, index] = point[1]
        coords[2, index] = point[2]
        coords[3, index] = point[3]
    end
    return nothing
end

# Volume-cell emission in Gmsh's i-outermost loop order. `mask` is already
# validated to one of the four admissible patterns by `_recombine_volume_mask`.
function _emit_six_recombined_cells!(hexes, prisms, shadow, coords,
                                     node_id, nu::Int, nv::Int, nw::Int, mask)
    hex_position = 0
    prism_position = 0
    shadow_position = 0
    all_recombined = all(mask)
    v_free = !mask[1] && !mask[3]
    u_free = !mask[2] && !mask[4]
    @inbounds for i in 0:nu-1, j in 0:nv-1, k in 0:nw-1
        a = node_id(i, j, k); b = node_id(i + 1, j, k)
        c = node_id(i, j + 1, k); d = node_id(i, j, k + 1)
        e = node_id(i + 1, j, k + 1); f = node_id(i, j + 1, k + 1)
        g = node_id(i + 1, j + 1, k); h = node_id(i + 1, j + 1, k + 1)
        if all_recombined
            hex_position += 1
            shadow_position = _emit_recombined_cell!(
                hexes, hex_position, shadow, shadow_position, coords,
                (a, b, g, c, d, e, h, f), _SHADOW_HEX)
        elseif v_free
            # Gmsh's literal macro order is negative-oriented; the emitted
            # order applies the MPrism orientation fixup (1<->2, 4<->5).
            prism_position += 1
            shadow_position = _emit_recombined_cell!(
                prisms, prism_position, shadow, shadow_position, coords,
                (b, a, d, g, c, f), _SHADOW_PRISM_1R)
            prism_position += 1
            shadow_position = _emit_recombined_cell!(
                prisms, prism_position, shadow, shadow_position, coords,
                (d, e, b, f, h, g), _SHADOW_PRISM_2R)
        elseif u_free
            prism_position += 1
            shadow_position = _emit_recombined_cell!(
                prisms, prism_position, shadow, shadow_position, coords,
                (g, b, e, c, a, d), _SHADOW_PRISM_1R)
            prism_position += 1
            shadow_position = _emit_recombined_cell!(
                prisms, prism_position, shadow, shadow_position, coords,
                (e, h, g, d, f, c), _SHADOW_PRISM_2R)
        else
            prism_position += 1
            shadow_position = _emit_recombined_cell!(
                prisms, prism_position, shadow, shadow_position, coords,
                (a, b, c, d, e, f), _SHADOW_PRISM_1)
            prism_position += 1
            shadow_position = _emit_recombined_cell!(
                prisms, prism_position, shadow, shadow_position, coords,
                (g, c, b, h, f, e), _SHADOW_PRISM_2)
        end
    end
    hex_position == size(hexes, 2) &&
        prism_position == size(prisms, 2) &&
        shadow_position == size(shadow, 2) || throw(ErrorException(
        "mesh_transfinite_volume: internal recombined-cell count invariant " *
        "failed"))
    return nothing
end

# Affine-path boundary emission: unrecombined slots emit the canonical
# tet-boundary triangle split, recombined slots emit one quadrangle per face
# cell in the parametric cyclic order (A,B,D,C). Every cell is wound strictly
# outward against the face's opposite corner.
function _emit_affine_recombined_boundary!(tris, tri_tags, quads, quad_tags,
        coords, node_id, opposite_corners, nu::Int, nv::Int, nw::Int,
        face_tags, mask)
    tri_position = 1
    quad_position = 1
    @inbounds for slot in 1:6
        a1, a2, fixed_axis, fixed = _WARPED_SLOT_PLANE[slot]
        n1, n2 = _warped_slot_dims(slot, nu, nv, nw)
        fixed_coord = fixed == 0 ? 0 :
            fixed_axis == 1 ? nu : fixed_axis == 2 ? nv : nw
        tag = face_tags[slot]
        opposite = opposite_corners[slot]
        function tab_node(p1::Int, p2::Int)
            coord = (0, 0, 0)
            coord = Base.setindex(coord, p1, a1)
            coord = Base.setindex(coord, p2, a2)
            coord = Base.setindex(coord, fixed_coord, fixed_axis)
            return node_id(coord[1], coord[2], coord[3])
        end
        if mask[slot]
            for p2 in 0:n2-1, p1 in 0:n1-1
                quad_position = _write_outward_quad!(
                    quads, quad_tags, quad_position, coords,
                    tab_node(p1, p2), tab_node(p1 + 1, p2),
                    tab_node(p1 + 1, p2 + 1), tab_node(p1, p2 + 1),
                    opposite, tag)
            end
        else
            for p2 in 0:n2-1, p1 in 0:n1-1
                a = tab_node(p1, p2); b = tab_node(p1 + 1, p2)
                c = tab_node(p1, p2 + 1); d = tab_node(p1 + 1, p2 + 1)
                tri_position = _write_outward_triangle!(
                    tris, tri_tags, tri_position, coords, a, b, c,
                    opposite, tag)
                tri_position = _write_outward_triangle!(
                    tris, tri_tags, tri_position, coords, b, d, c,
                    opposite, tag)
            end
        end
    end
    tri_position == size(tris, 2) + 1 &&
        quad_position == size(quads, 2) + 1 || throw(ErrorException(
        "mesh_transfinite_volume: internal boundary cell count invariant " *
        "failed"))
    return nothing
end

function _transfinite_volume_recombined(converted_corners, raw_faces,
        raw_face_tags, mask::NTuple{6,Bool}, nu::Int, nv::Int, nw::Int,
        npu::Int, npv::Int, npw::Int, node_count::Int, tet_count::Int,
        uv::Int, uw::Int, vw::Int,
        node_limit::Int, tet_limit::Int, boundary_limit::Int,
        converted_volume_tag::Int32)
    all_recombined = all(mask)
    cell_count = _checked_mul("cell", nu, nv, nw)
    hex_count = all_recombined ? cell_count : 0
    prism_count = all_recombined ? 0 : _checked_mul("prism", 2, cell_count)
    slot_cells = (uw, vw, uw, vw, uv, uv)
    tri_count = 0
    quad_count = 0
    @inbounds for slot in 1:6
        if mask[slot]
            quad_count = _checked_add(
                quad_count, slot_cells[slot], "boundary cell")
        else
            tri_count = _checked_add(
                tri_count, _checked_mul("boundary triangle", 2,
                                        slot_cells[slot]),
                "boundary triangle")
        end
    end
    boundary_cells = _checked_add(tri_count, quad_count, "boundary cell")
    boundary_cells <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $boundary_cells boundary cells exceed " *
        "the Int32 topology limit"))
    boundary_cells <= boundary_limit || throw(ArgumentError(
        "mesh_transfinite_volume: $boundary_cells boundary cells exceed " *
        "max_boundary_triangles=$boundary_limit"))

    node_id(i::Int, j::Int, k::Int) = Int32(i + 1 + npu * (j + npv * k))
    coords = Matrix{Float64}(undef, 3, node_count)
    faces = nothing
    if raw_faces === nothing
        _certify_affine(converted_corners)
        _fill_affine_volume_coords!(coords, converted_corners, nu, nv, nw,
                                    node_id)
    else
        faces = _warped_faces(raw_faces, nu, nv, nw)
        _certify_shared_edges(faces, nu, nv, nw)
        _fill_warped_volume_coords!(coords, converted_corners, faces,
                                    nu, nv, nw, node_id)
    end

    hexes = Matrix{Int32}(undef, 8, hex_count)
    prisms = Matrix{Int32}(undef, 6, prism_count)
    shadow = Matrix{Int32}(undef, 4, tet_count)
    _emit_six_recombined_cells!(hexes, prisms, shadow, coords, node_id,
                                nu, nv, nw, mask)

    # When the shadow reproduces the unrecombined partition (all-hex and the
    # w-free prism pair), certify it against the reference simplex
    # subdivision; the partial v/u-free tilings are instead certified by
    # orientation, conformity, and boundary coverage below.
    if all_recombined || (mask[1] && mask[2] && mask[3] && mask[4])
        reference = Matrix{Int32}(undef, 4, tet_count)
        _emit_six_simplex_tets!(reference, coords, node_id, nu, nv, nw) ==
            tet_count || throw(ErrorException(
            "mesh_transfinite_volume: internal reference partition count " *
            "invariant failed"))
        _canonical_tets(shadow) == _canonical_tets(reference) ||
            throw(ErrorException(
            "mesh_transfinite_volume: recombined cells do not partition " *
            "the reference simplex mesh"))
    end

    tris = Matrix{Int32}(undef, 3, tri_count)
    quads = Matrix{Int32}(undef, 4, quad_count)
    tri_tags = Vector{Int32}(undef, tri_count)
    quad_tags = Vector{Int32}(undef, quad_count)
    converted_face_tags = _face_tags(raw_face_tags)
    if faces === nothing
        opposite_corners = (converted_corners[4], converted_corners[1],
                            converted_corners[1], converted_corners[2],
                            converted_corners[5], converted_corners[1])
        _emit_affine_recombined_boundary!(tris, tri_tags, quads, quad_tags,
            coords, node_id, opposite_corners, nu, nv, nw,
            converted_face_tags, mask)
        # Every emitted affine boundary cell must be strictly outward against
        # the interior opposite corner; witness emission guarantees the
        # winding but this independent pass keeps the certification absolute.
        triangle = 0
        @inbounds for slot in 1:6
            mask[slot] && continue
            opposite = opposite_corners[slot]
            for _ in 1:2*slot_cells[slot]
                triangle += 1
                orient3(_node(coords, tris[1, triangle]),
                        _node(coords, tris[2, triangle]),
                        _node(coords, tris[3, triangle]), opposite) > 0 ||
                    throw(ArgumentError(
                        "mesh_transfinite_volume: boundary triangle " *
                        "$triangle is not strictly outward-oriented"))
            end
        end
        quad = 0
        @inbounds for slot in 1:6
            mask[slot] || continue
            opposite = opposite_corners[slot]
            for _ in 1:slot_cells[slot]
                quad += 1
                a = _node(coords, quads[1, quad])
                b = _node(coords, quads[2, quad])
                c = _node(coords, quads[3, quad])
                d = _node(coords, quads[4, quad])
                orient3(a, b, c, opposite) > 0 &&
                    orient3(a, c, d, opposite) > 0 || throw(ArgumentError(
                        "mesh_transfinite_volume: boundary quadrangle " *
                        "$quad is not strictly outward-oriented"))
            end
        end
    else
        _emit_warped_volume_boundary!(tris, tri_tags, quads, quad_tags,
                                      coords, node_id, nu, nv, nw,
                                      converted_face_tags, mask)
    end

    if faces === nothing
        _certify_tet_volume(
            coords, shadow,
            (converted_corners[1], converted_corners[2],
             converted_corners[4], converted_corners[5]),
            6, "mesh_transfinite_volume", "affine block")
    end
    _certify_recombined_boundary(shadow, tris, quads)
    mesh = _recombined_mixed_mesh(coords, tris, tri_tags, quads, quad_tags,
                                  hexes, prisms, converted_volume_tag)
    size(mesh.coords, 2) == node_count || throw(ErrorException(
        "mesh_transfinite_volume: finalized mesh count invariant failed"))
    return mesh
end

end
