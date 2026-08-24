"""
    TransfiniteQuad

Validated recombined four-sided transfinite quadrangle patches for an affine
planar surface. Node placement and boundary certification are delegated to the
public four-sided [`Transfinite.mesh_transfinite_patch`](@ref) contract; the
certified logical grid is emitted as linear Gmsh type-1 boundary lines and
type-3 quadrangles in an [`Elements.MixedMesh`](@ref).

This module does not discretize curves, apply curve laws, map general CAD
parameterizations, repair unequal opposite-side counts, smooth, handle holes or
periodic seams, generate high-order/recombined three-sided meshes, or construct
transfinite volumes.
"""
module TransfiniteQuad

using ..MeshTypes: Mesh, nnodes, nsegs
using ..Transfinite: mesh_transfinite_patch
using ..Predicates: orient2
import ..Elements
using ..Elements: ElementBlock, MixedMesh

export mesh_transfinite_quad_patch

const _CALLER = "mesh_transfinite_quad_patch"
const _DEFAULT_MAX_NODES = 10_000_000
const _DEFAULT_MAX_QUADRANGLES = 10_000_000
const _INT32_MAX = Int(typemax(Int32))

function _checked_add(a::Int, b::Int, what::AbstractString)
    try
        return Base.checked_add(a, b)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError("$_CALLER: $what count overflows Int"))
    end
end

function _checked_mul(a::Int, b::Int, what::AbstractString)
    try
        return Base.checked_mul(a, b)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError("$_CALLER: $what count overflows Int"))
    end
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

function _tag(value, name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$_CALLER: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$_CALLER: $name must not be Bool"))
    0 <= value <= typemax(Int32) || throw(ArgumentError(
        "$_CALLER: $name must lie in 0:$(typemax(Int32))"))
    return Int32(value)
end

function _arrangement(value)
    value isa Symbol || throw(ArgumentError(
        "$_CALLER: arrangement must be a Symbol"))
    value in (:left, :right, :alternate_left, :alternate_right) ||
        throw(ArgumentError(
            "$_CALLER: arrangement must be :left, :right, :alternate_left, " *
            "or :alternate_right"))
    return value
end

@inline _node(i::Int, j::Int, width::Int) = Int32(i + 1 + j * width)

@inline function _exact_delta(coords, index::Int, origin)
    return ntuple(3) do dimension
        Rational{BigInt}(coords[dimension, index]) - origin[dimension]
    end
end

@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                         a[3] * b[1] - a[1] * b[3],
                         a[1] * b[2] - a[2] * b[1])

function _exact_dropped_axis(mesh::Mesh, origin)
    exact_origin = ntuple(d -> Rational{BigInt}(origin[d]), 3)
    zero_exact = Rational{BigInt}(0)
    axis = (zero_exact, zero_exact, zero_exact)
    axis_norm = zero_exact
    @inbounds for segment in axes(mesh.segs, 2), endpoint in 1:2
        delta = _exact_delta(
            mesh.coords, Int(mesh.segs[endpoint, segment]), exact_origin)
        delta_norm = _dot3(delta, delta)
        if delta_norm > axis_norm
            axis = delta
            axis_norm = delta_norm
        end
    end
    axis_norm > 0 || throw(ErrorException(
        "$_CALLER: certified grid has no exact boundary extent"))

    best = (zero_exact, zero_exact, zero_exact)
    best_norm = zero_exact
    @inbounds for segment in axes(mesh.segs, 2), endpoint in 1:2
        delta = _exact_delta(
            mesh.coords, Int(mesh.segs[endpoint, segment]), exact_origin)
        candidate = _cross3(axis, delta)
        candidate_norm = _dot3(candidate, candidate)
        if candidate_norm > best_norm
            best = candidate
            best_norm = candidate_norm
        end
    end
    best_norm > 0 || throw(ErrorException(
        "$_CALLER: certified boundary has no exact plane normal"))
    magnitude = max(abs(best[1]), abs(best[2]), abs(best[3]))
    return abs(best[1]) == magnitude ? 1 : abs(best[2]) == magnitude ? 2 : 3
end

function _fill_quadrangles!(quadrangles, width::Int, horizontal::Int,
                            vertical::Int)
    cursor = 0
    @inbounds for i in 0:horizontal-1
        for j in 0:vertical-1
            cursor += 1
            quadrangles[1, cursor] = _node(i, j, width)
            quadrangles[2, cursor] = _node(i + 1, j, width)
            quadrangles[3, cursor] = _node(i + 1, j + 1, width)
            quadrangles[4, cursor] = _node(i, j + 1, width)
        end
    end
    cursor == size(quadrangles, 2) || throw(ErrorException(
        "$_CALLER: internal quadrangle count invariant failed"))
    return nothing
end

function _projection(mesh::Mesh)
    first_node = Int(mesh.segs[1, 1])
    origin = (mesh.coords[1, first_node], mesh.coords[2, first_node],
              mesh.coords[3, first_node])
    scale = 0.0
    @inbounds for index in axes(mesh.coords, 2)
        dx = mesh.coords[1, index] - origin[1]
        dy = mesh.coords[2, index] - origin[2]
        dz = mesh.coords[3, index] - origin[3]
        (isfinite(dx) && isfinite(dy) && isfinite(dz)) || throw(ArgumentError(
            "$_CALLER: certified coordinate span is not finite"))
        scale = max(scale, abs(dx), abs(dy), abs(dz))
    end
    (isfinite(scale) && scale > 0) || throw(ErrorException(
        "$_CALLER: certified grid has no representable extent"))

    nx = 0.0
    ny = 0.0
    nz = 0.0
    accumulated_magnitude = 0.0
    @inbounds for segment in axes(mesh.segs, 2)
        first = Int(mesh.segs[1, segment])
        second = Int(mesh.segs[2, segment])
        px = (mesh.coords[1, first] - origin[1]) / scale
        py = (mesh.coords[2, first] - origin[2]) / scale
        pz = (mesh.coords[3, first] - origin[3]) / scale
        qx = (mesh.coords[1, second] - origin[1]) / scale
        qy = (mesh.coords[2, second] - origin[2]) / scale
        qz = (mesh.coords[3, second] - origin[3]) / scale
        tx = (py - qy) * (pz + qz)
        ty = (pz - qz) * (px + qx)
        tz = (px - qx) * (py + qy)
        nx += tx
        ny += ty
        nz += tz
        accumulated_magnitude += abs(tx) + abs(ty) + abs(tz)
    end
    all(isfinite, (nx, ny, nz, accumulated_magnitude)) ||
        throw(ErrorException(
            "$_CALLER: certified boundary normal is not finite"))
    magnitude = max(abs(nx), abs(ny), abs(nz))
    error_bound = 128eps(Float64) * size(mesh.segs, 2) *
                  accumulated_magnitude
    dropped = if magnitude > error_bound
        abs(nx) == magnitude ? 1 : abs(ny) == magnitude ? 2 : 3
    else
        _exact_dropped_axis(mesh, origin)
    end
    return dropped
end

@inline function _project(coords, node::Int32, dropped::Int)
    index = Int(node)
    x = coords[1, index]
    y = coords[2, index]
    z = coords[3, index]
    all(isfinite, (x, y, z)) || throw(ArgumentError(
        "$_CALLER: quadrangle node cannot be projected"))
    return dropped == 1 ? (y, z) : dropped == 2 ? (z, x) : (x, y)
end

function _validate_quadrangle_geometry(mesh::Mesh, quadrangles)
    dropped = _projection(mesh)
    reference = 0
    @inbounds for cell in axes(quadrangles, 2)
        p1 = _project(mesh.coords, quadrangles[1, cell], dropped)
        p2 = _project(mesh.coords, quadrangles[2, cell], dropped)
        p3 = _project(mesh.coords, quadrangles[3, cell], dropped)
        p4 = _project(mesh.coords, quadrangles[4, cell], dropped)
        # A planar bilinear quadrangle's Jacobian determinant is affine in its
        # reference coordinates. These four turns are its four corner values,
        # so one nonzero common exact sign certifies the whole reference cell.
        turns = (orient2(p1, p2, p3), orient2(p2, p3, p4),
                 orient2(p3, p4, p1), orient2(p4, p1, p2))
        all(!=(0), turns) || throw(ArgumentError(
            "$_CALLER: quadrangle $cell has a zero corner Jacobian"))
        all(==(turns[1]), turns) || throw(ArgumentError(
            "$_CALLER: quadrangle $cell is concave or folded"))
        if reference == 0
            reference = turns[1]
        elseif turns[1] != reference
            throw(ArgumentError(
                "$_CALLER: quadrangle $cell reverses patch orientation"))
        end
    end
    return reference
end

"""
    mesh_transfinite_quad_patch(side1, side2, side3, side4;
        arrangement=:left, face_tag=0, side_tags=(0,0,0,0),
        max_nodes=10_000_000, max_quadrangles=10_000_000) -> MixedMesh

Construct a recombined four-sided planar transfinite patch. The four boundary
chains are already-discretized finite 3-D points oriented cyclically as
`c1→c2`, `c2→c3`, `c3→c4`, and `c4→c1`; opposite chains must have
matching node counts. Gmsh type-1 line elements preserve `side_tags`, while
Gmsh type-3 quadrangles preserve `face_tag`.

The accepted arrangements are `:left`, `:right`, `:alternate_left`, and
`:alternate_right`. They intentionally produce identical output: Gmsh 4.15.2's
four-sided recombined path emits `(v1,v2,v3,v4)` before consulting triangle
arrangement. The operation certifies the Coons grid through the public
four-sided simplex implementation, then additionally checks all four exact
projected corner Jacobians of every linear quadrangle. Concave, folded,
degenerate, inverted, non-manifold, or invalid boundaries are rejected without
fallback.

Counts, Int32 bounds, and caller limits are checked before output allocation.
The internal certification uses two triangles per logical quadrangle, so this
bounded API additionally requires `2 * nquadrangles <= typemax(Int32)`.

This function does not discretize curves, apply curve laws or fields, repair
unequal counts, smooth, handle holes/periodic seams/embedded entities, map a
general CAD/ruled/spherical parameterization, generate high-order elements,
recombine three-sided patches, create lossless CAD entity metadata, or construct
transfinite volumes.
"""
function mesh_transfinite_quad_patch(side1::AbstractVector,
                                     side2::AbstractVector,
                                     side3::AbstractVector,
                                     side4::AbstractVector;
                                     arrangement=:left,
                                     face_tag=0,
                                     side_tags=(0, 0, 0, 0),
                                     max_nodes=_DEFAULT_MAX_NODES,
                                     max_quadrangles=_DEFAULT_MAX_QUADRANGLES)::MixedMesh
    _arrangement(arrangement)
    node_limit = _limit(max_nodes, "max_nodes")
    quadrangle_limit = _limit(max_quadrangles, "max_quadrangles")
    side_tags isa Tuple && length(side_tags) == 4 || throw(ArgumentError(
        "$_CALLER: side_tags must be a four-integer tuple"))
    physical_side_tags = ntuple(
        index -> _tag(side_tags[index], "side_tags[$index]"), 4)
    physical_face_tag = _tag(face_tag, "face_tag")

    lengths = (length(side1), length(side2), length(side3), length(side4))
    @inbounds for side in 1:4
        lengths[side] >= 2 || throw(ArgumentError(
            "$_CALLER: side $side needs at least two points"))
    end
    lengths[1] == lengths[3] || throw(ArgumentError(
        "$_CALLER: opposite sides 1 and 3 have non-matching node counts " *
        "$(lengths[1]) and $(lengths[3])"))
    lengths[2] == lengths[4] || throw(ArgumentError(
        "$_CALLER: opposite sides 2 and 4 have non-matching node counts " *
        "$(lengths[2]) and $(lengths[4])"))

    horizontal = lengths[1] - 1
    vertical = lengths[2] - 1
    nodes = _checked_mul(lengths[1], lengths[2], "node")
    quadrangles = _checked_mul(horizontal, vertical, "quadrangle")
    certification_triangles = _checked_mul(2, quadrangles,
                                            "certification triangle")
    segments = _checked_mul(
        2, _checked_add(horizontal, vertical, "segment"), "segment")
    nodes <= _INT32_MAX || throw(ArgumentError(
        "$_CALLER: $nodes nodes exceed the Int32 indexing limit"))
    quadrangles <= _INT32_MAX || throw(ArgumentError(
        "$_CALLER: $quadrangles quadrangles exceed the Int32 topology limit"))
    certification_triangles <= _INT32_MAX || throw(ArgumentError(
        "$_CALLER: $quadrangles quadrangles exceed the bounded Int32 " *
        "certification-triangle limit $(_INT32_MAX ÷ 2)"))
    segments <= _INT32_MAX || throw(ArgumentError(
        "$_CALLER: $segments segments exceed the Int32 topology limit"))
    nodes <= node_limit || throw(ArgumentError(
        "$_CALLER: $nodes nodes exceed max_nodes=$node_limit"))
    quadrangles <= quadrangle_limit || throw(ArgumentError(
        "$_CALLER: $quadrangles quadrangles exceed " *
        "max_quadrangles=$quadrangle_limit"))

    # A fixed :left split is only a transient validity certificate. Gmsh's
    # recombined four-sided path ignores arrangement, as does the returned mesh.
    certified = mesh_transfinite_patch(
        side1, side2, side3, side4;
        arrangement=:left,
        face_tag=physical_face_tag,
        side_tags=physical_side_tags,
        max_nodes=node_limit,
        max_triangles=certification_triangles)
    (nnodes(certified) == nodes && nsegs(certified) == segments) ||
        throw(ErrorException(
            "$_CALLER: certified grid count postcondition failed"))

    quadrangle_topology = Matrix{Int32}(undef, 4, quadrangles)
    _fill_quadrangles!(quadrangle_topology, lengths[1], horizontal, vertical)
    _validate_quadrangle_geometry(certified, quadrangle_topology)
    line_block = ElementBlock(1, certified.segs, certified.seg_tag)
    quadrangle_block = ElementBlock(
        3, quadrangle_topology, fill(physical_face_tag, quadrangles))
    result = MixedMesh(certified.coords, [line_block, quadrangle_block])
    diagnostic = Elements.validate(result)
    diagnostic.ok || throw(ErrorException(
        "$_CALLER: internal MixedMesh validation failed — " *
        join(diagnostic.messages, "; ")))
    (size(result.coords, 2) == nodes && length(result.blocks) == 2 &&
     size(result.blocks[1].nodes, 2) == segments &&
     size(result.blocks[2].nodes, 2) == quadrangles) || throw(ErrorException(
        "$_CALLER: output count postcondition failed"))
    return result
end

end # module TransfiniteQuad
