"""
    TransfinitePrism

Bounded Gmsh-4.15.2-compatible five-face transfinite meshing for an affine
triangular prism. The result is a first-order simplex `Mesh` with the exact
unrecombined tetrahedron topology used by Gmsh's legacy collapsed-grid
(`Mesh.TransfiniteTri = 0`) volume path.

This module deliberately does not claim support for the specific compact
triangular algorithm (`Mesh.TransfiniteTri = 1`), curved or independently
discretized faces, non-affine prisms, nonuniform curve laws, recombined prisms
or hexahedra, QuadTri, holes, multiple blocks, periodic seams, or high-order
elements. As required by the finalized `Mesh` contract, represented boundary
areas and tetrahedron volumes must also remain finite Float64 values; finite
input coordinates alone do not imply that their derived measures are finite.
"""
module TransfinitePrism

using ..MeshTypes: Mesh, boundary_faces, nnodes, ntris, ntets, validate
using ..Predicates: orient3
using ..StructuredNumerics: _certify_tet_volume, _throw_simplex_validation

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

"""
    mesh_transfinite_prism(corners, cells=(1,1,1);
                          volume_tag=0,
                          face_tags=(0,0,0,0,0),
                          max_nodes=10_000_000,
                          max_tets=60_000_000,
                          max_boundary_triangles=20_000_000) -> Mesh

Mesh an affine triangular prism using Gmsh 4.15.2's five-face legacy
collapsed-grid transfinite algorithm (`Mesh.TransfiniteTri = 0`) with all
surfaces unrecombined and using the `Left` triangle arrangement. `corners` must
contain six finite points in canonical order `(s0,s1,s2,s4,s5,s6)`: the first
three define the lower triangular face and the final three are their common
affine translation. The canonical orientation must be positive.

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

Unsupported here: `Mesh.TransfiniteTri = 1`, curved/warped or independently
discretized faces, nonuniform curve laws, recombination, QuadTri, holes,
multiple blocks, periodic seams, high-order elements, and coordinate scales
whose derived triangle areas or tetrahedron volumes are not finite Float64
values.
"""
function mesh_transfinite_prism(corners, cells=(1, 1, 1);
                                volume_tag=0,
                                face_tags=(0, 0, 0, 0, 0),
                                max_nodes=_DEFAULT_MAX_NODES,
                                max_tets=_DEFAULT_MAX_TETS,
                                max_boundary_triangles=
                                    _DEFAULT_MAX_BOUNDARY_TRIANGLES)
    node_limit = _limit(max_nodes, "max_nodes")
    tet_limit = _limit(max_tets, "max_tets")
    triangle_limit = _limit(max_boundary_triangles, "max_boundary_triangles")
    nr, ns, nw = _three_counts(cells)

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
    _certify_affine(converted_corners)
    converted_face_tags = _face_tags(face_tags)
    converted_volume_tag = _tag(volume_tag, "volume_tag")

    node_id(i::Int, j::Int, k::Int) = Int32(
        k * layer_nodes + (i == 0 ? 1 : 2 + (i - 1) * opposite_nodes + j))

    coords = Matrix{Float64}(undef, 3, node_count)
    @inbounds for k in 0:nw
        axial = k / nw
        collapsed = _prism_point(converted_corners, 0.0, 0.0, axial)
        all(isfinite, collapsed) || throw(ArgumentError(
            "$_CALLER: interpolation produced a non-finite coordinate at " *
            "logical node (0,0,$k)"))
        collapsed_index = Int(node_id(0, 0, k))
        coords[1, collapsed_index] = collapsed[1]
        coords[2, collapsed_index] = collapsed[2]
        coords[3, collapsed_index] = collapsed[3]
        for i in 1:nr, j in 0:ns
            point = _prism_point(
                converted_corners, i / nr, j / ns, axial)
            all(isfinite, point) || throw(ArgumentError(
                "$_CALLER: interpolation produced a non-finite coordinate at " *
                "logical node ($i,$j,$k)"))
            index = Int(node_id(i, j, k))
            coords[1, index] = point[1]
            coords[2, index] = point[2]
            coords[3, index] = point[3]
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
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,d,converted_corners[3],
            converted_face_tags[1])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,b,e,d,converted_corners[3],
            converted_face_tags[1])
    end
    # f1: s1-s2 axial side (outer radial row).
    @inbounds for j in 0:ns-1, k in 0:nw-1
        a=node_id(nr,j,k); b=node_id(nr,j+1,k)
        d=node_id(nr,j,k+1); e=node_id(nr,j+1,k+1)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,d,converted_corners[1],
            converted_face_tags[2])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,b,e,d,converted_corners[1],
            converted_face_tags[2])
    end
    # f2: s0-s2 axial side (opposite coordinate j=ns).
    @inbounds for i in 0:nr-1, k in 0:nw-1
        a=node_id(i,ns,k); b=node_id(i+1,ns,k)
        d=node_id(i,ns,k+1); e=node_id(i+1,ns,k+1)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,d,converted_corners[2],
            converted_face_tags[3])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,b,e,d,converted_corners[2],
            converted_face_tags[3])
    end
    # f4: lower triangular face.
    @inbounds for j in 0:ns-1
        a=node_id(0,0,0); b=node_id(1,j,0); c=node_id(1,j+1,0)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,c,converted_corners[4],
            converted_face_tags[4])
    end
    @inbounds for i in 1:nr-1, j in 0:ns-1
        a=node_id(i,j,0); b=node_id(i+1,j,0)
        c=node_id(i,j+1,0); g=node_id(i+1,j+1,0)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,c,converted_corners[4],
            converted_face_tags[4])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,c,b,g,converted_corners[4],
            converted_face_tags[4])
    end
    # f5: upper triangular face.
    @inbounds for j in 0:ns-1
        a=node_id(0,0,nw); b=node_id(1,j,nw); c=node_id(1,j+1,nw)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,c,converted_corners[1],
            converted_face_tags[5])
    end
    @inbounds for i in 1:nr-1, j in 0:ns-1
        a=node_id(i,j,nw); b=node_id(i+1,j,nw)
        c=node_id(i,j+1,nw); g=node_id(i+1,j+1,nw)
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,a,b,c,converted_corners[1],
            converted_face_tags[5])
        position=_write_outward_triangle!(
            tris,tri_tags,position,coords,c,b,g,converted_corners[1],
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

end
