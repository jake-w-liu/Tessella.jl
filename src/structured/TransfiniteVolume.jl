"""
    TransfiniteVolume

Bounded Gmsh-4.15.2-compatible transfinite volume meshing for an affine,
six-faced, eight-corner block.  The result is a first-order simplex `Mesh`:
each logical hexahedron is split into the six tetrahedra used by Gmsh when all
six transfinite boundary faces are unrecombined.

This module deliberately does not claim support for five-faced/prismatic
volumes, curved or independently discretized faces, non-affine hexahedra,
recombined hexahedra/prisms, QuadTri, holes, multiple blocks, periodic seams,
or high-order elements. As required by the finalized `Mesh` contract,
represented boundary areas and tetrahedron volumes must remain finite Float64
values; finite input coordinates alone do not imply finite derived measures.
"""
module TransfiniteVolume

using ..MeshTypes: Mesh, boundary_faces, nnodes, ntris, ntets, validate
using ..Predicates: orient3
using ..StructuredNumerics: _needs_exact_affine, _affine_basis3,
                            _affine_point3, _uses_exact_affine,
                            _certify_tet_volume, _throw_simplex_validation

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

"""
    mesh_transfinite_volume(corners, cells=(1,1,1);
                            volume_tag=0,
                            face_tags=(0,0,0,0,0,0),
                            max_nodes=10_000_000,
                            max_tets=60_000_000,
                            max_boundary_triangles=20_000_000) -> Mesh

Mesh an affine six-face block with the Gmsh 4.15.2 unrecombined transfinite
volume subdivision. `corners` must contain eight finite 3-D points in Gmsh's
canonical order `(s0,s1,s2,s3,s4,s5,s6,s7)`: the first four wind around the
`w=0` face and the final four are their `w=1` counterparts. The four derived
corners must agree with an affine parallelepiped to a normalized,
conditioning-scaled `4096eps(Float64)` tolerance. `cells=(nu,nv,nw)` gives
positive logical-cell counts; curve laws are uniformly spaced (`Progression 1`).

The returned simplex mesh uses the exact six-tetrahedron connectivity pattern
from Gmsh's `CREATE_SIM_1` through `CREATE_SIM_6`. Boundary triangles are the
induced conforming face arrangements. Exact dyadic affine interpolation protects
remote, narrow blocks whose nested Float64 interpolation would lose material, and
a compensated exponent-scaled determinant audit certifies conservation of the
corner-defined volume. `face_tags` follow Gmsh's canonical face order
`(vmin, umax, vmax, umin, wmin, wmax)`; `volume_tag` labels every tet.

Unsupported here: five-face degeneracies, curved/warped or independently
discretized faces, nonuniform curve laws, recombination into hexahedra/prisms,
QuadTri, holes, multiple blocks, periodic seams, high-order elements, and
coordinate scales whose derived boundary areas or tetrahedron volumes are not
finite Float64 values.
"""
function mesh_transfinite_volume(corners, cells=(1, 1, 1);
                                 volume_tag=0,
                                 face_tags=(0, 0, 0, 0, 0, 0),
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
    face_cells = _checked_add(_checked_add(uv, uw, "boundary triangle"),
                              vw, "boundary triangle")
    triangle_count = _checked_mul("boundary triangle", 4, face_cells)

    node_count <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $node_count nodes exceed Int32 indexing"))
    tet_count <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $tet_count tetrahedra exceed the Int32 topology limit"))
    triangle_count <= typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_volume: $triangle_count boundary triangles exceed " *
        "the Int32 topology limit"))
    node_count <= node_limit || throw(ArgumentError(
        "mesh_transfinite_volume: $node_count nodes exceed max_nodes=$node_limit"))
    tet_count <= tet_limit || throw(ArgumentError(
        "mesh_transfinite_volume: $tet_count tetrahedra exceed max_tets=$tet_limit"))
    triangle_count <= triangle_limit || throw(ArgumentError(
        "mesh_transfinite_volume: $triangle_count boundary triangles exceed " *
        "max_boundary_triangles=$triangle_limit"))

    converted_corners = _eight_corners(corners)
    _certify_affine(converted_corners)
    converted_face_tags = _face_tags(face_tags)
    converted_volume_tag = _tag(volume_tag, "volume_tag")
    affine_points = (converted_corners[1], converted_corners[2],
                     converted_corners[4], converted_corners[5])
    exact_interpolation = _needs_exact_affine(
        affine_points..., (nu, nv, nw)) &&
        _exact_affine_corners(converted_corners)
    affine_basis = _affine_basis3(affine_points..., exact_interpolation)

    node_id(i::Int, j::Int, k::Int) = Int32(i + 1 + npu * (j + npv * k))
    coords = Matrix{Float64}(undef, 3, node_count)
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

    tets = Matrix{Int32}(undef, 4, tet_count)
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
    tet_position == tet_count || throw(ErrorException(
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

end
