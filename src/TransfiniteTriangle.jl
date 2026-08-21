"""
    TransfiniteTriangle

Validated planar, three-sided transfinite triangle patches matching the
unrecombined `Mesh.TransfiniteTri = 1` path in Gmsh 4.15.2. Boundary chains are
supplied already discretized, cyclically oriented, and must have the same node
count. Gmsh's triangular interpolation, chord parameters, compact triangular
connectivity, and arrangement-independent triangle path are reproduced for an
affine planar surface.

This module does not discretize curves, implement Gmsh's legacy
`Mesh.TransfiniteTri = 0` collapsed-grid algorithm, quasi-transfinite repair,
CAD/ruled/spherical parameterizations, recombination or quadrangles, smoothing,
holes, periodic seams, embedded entities, size/quality fields, or transfinite
volumes.
"""
module TransfiniteTriangle

using ..MeshTypes: Mesh, boundary_edges, nnodes, nsegs, ntris, validate
using ..Predicates: orient2

export mesh_transfinite_triangle

const _CALLER = "mesh_transfinite_triangle"
const _DEFAULT_MAX_NODES = 10_000_000
const _DEFAULT_MAX_TRIANGLES = 20_000_000
const _INT32_MAX = Int(typemax(Int32))
const _BOUNDARY_AUDIT_MULTIPLIER = 64
const _BOUNDARY_AUDIT_FLOOR = 4096

struct _PlaneFrame
    u::NTuple{3,Float64}
    v::NTuple{3,Float64}
    n::NTuple{3,Float64}
end

struct _SegmentBox
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    segment::Int
end

struct _BoundaryNode
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    segment::Int
    left::Int
    right::Int
    count::Int
end

@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                         a[3] * b[1] - a[1] * b[3],
                         a[1] * b[2] - a[2] * b[1])
@inline _norm3(a) = hypot(a[1], a[2], a[3])
@inline _edge_key(a::Int32, b::Int32) = a < b ? (a, b) : (b, a)

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

function _point3(raw, side::Int, index::Int)
    count = try
        length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: side $side point $index is not indexable"))
    end
    count == 3 || throw(ArgumentError(
        "$_CALLER: side $side point $index must have exactly three coordinates"))
    point = try
        (Float64(raw[1]), Float64(raw[2]), Float64(raw[3]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$_CALLER: side $side point $index coordinates must be " *
            "Float64-representable: $(sprint(showerror, err))"))
    end
    all(isfinite, point) || throw(ArgumentError(
        "$_CALLER: side $side point $index has a non-finite coordinate"))
    return point
end

function _convert_side(side::AbstractVector, number::Int)
    result = Vector{NTuple{3,Float64}}(undef, length(side))
    destination = 1
    for raw in side
        destination <= length(result) || throw(ArgumentError(
            "$_CALLER: side $number iteration produced more than its declared length"))
        result[destination] = _point3(raw, number, destination)
        destination += 1
    end
    destination == length(result) + 1 || throw(ArgumentError(
        "$_CALLER: side $number iteration length changed during conversion"))
    return result
end

function _distance(a, b, description::AbstractString)
    dx = a[1] - b[1]
    dy = a[2] - b[2]
    dz = a[3] - b[3]
    (isfinite(dx) && isfinite(dy) && isfinite(dz)) || throw(ArgumentError(
        "$_CALLER: $description coordinate span overflows Float64"))
    distance = hypot(dx, dy, dz)
    (isfinite(distance) && distance > 0) || throw(ArgumentError(
        "$_CALLER: $description has zero or non-finite length"))
    return distance
end

function _validate_side_edges(side, number::Int)
    @inbounds for index in 1:length(side)-1
        _distance(side[index + 1], side[index], "side $number segment $index")
    end
    return nothing
end

function _boundary_ring(sides)
    count = sum(length(side) - 1 for side in sides; init=0)
    ring = Vector{NTuple{3,Float64}}(undef, count)
    cursor = 0
    @inbounds for side in sides
        for index in 1:length(side)-1
            cursor += 1
            ring[cursor] = side[index]
        end
    end
    cursor == count || throw(ErrorException(
        "$_CALLER: internal boundary count invariant failed"))
    return ring
end

function _normalization(ring)
    origin = ring[1]
    scale = 0.0
    @inbounds for (index, point) in pairs(ring)
        dx = point[1] - origin[1]
        dy = point[2] - origin[2]
        dz = point[3] - origin[3]
        (isfinite(dx) && isfinite(dy) && isfinite(dz)) || throw(ArgumentError(
            "$_CALLER: boundary coordinate span overflows Float64 at node $index"))
        scale = max(scale, abs(dx), abs(dy), abs(dz))
    end
    (isfinite(scale) && scale > 0) || throw(ArgumentError(
        "$_CALLER: boundary is geometrically degenerate"))
    return origin, scale
end

@inline _normalize(point, origin, scale) =
    ((point[1] - origin[1]) / scale,
     (point[2] - origin[2]) / scale,
     (point[3] - origin[3]) / scale)

function _normalized_side(side, origin, scale)
    result = Vector{NTuple{3,Float64}}(undef, length(side))
    @inbounds for index in eachindex(side)
        point = _normalize(side[index], origin, scale)
        all(isfinite, point) || throw(ArgumentError(
            "$_CALLER: normalized boundary coordinate is not finite"))
        result[index] = point
    end
    return result
end

function _plane_frame(ring, origin, scale)
    nx = 0.0
    ny = 0.0
    nz = 0.0
    @inbounds for index in eachindex(ring)
        point = _normalize(ring[index], origin, scale)
        next = _normalize(ring[mod1(index + 1, length(ring))], origin, scale)
        nx += (point[2] - next[2]) * (point[3] + next[3])
        ny += (point[3] - next[3]) * (point[1] + next[1])
        nz += (point[1] - next[1]) * (point[2] + next[2])
    end
    normal_length = hypot(nx, ny, nz)
    (isfinite(normal_length) && normal_length > 0) || throw(ArgumentError(
        "$_CALLER: boundary has no representable plane normal"))
    normal = (nx / normal_length, ny / normal_length, nz / normal_length)
    reference = abs(normal[1]) <= abs(normal[2]) ?
        (abs(normal[1]) <= abs(normal[3]) ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)) :
        (abs(normal[2]) <= abs(normal[3]) ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0))
    axis_u = _cross3(reference, normal)
    axis_length = _norm3(axis_u)
    (isfinite(axis_length) && axis_length > 0) || throw(ArgumentError(
        "$_CALLER: could not construct an in-plane frame"))
    axis_u = (axis_u[1] / axis_length,
              axis_u[2] / axis_length,
              axis_u[3] / axis_length)
    axis_v = _cross3(normal, axis_u)
    frame = _PlaneFrame(axis_u, axis_v, normal)
    tolerance = 256eps(Float64)
    @inbounds for (index, point) in pairs(ring)
        normalized = _normalize(point, origin, scale)
        distance = abs(_dot3(normalized, normal))
        (isfinite(distance) && distance <= tolerance) || throw(ArgumentError(
            "$_CALLER: boundary node $index is not coplanar " *
            "(normalized distance $distance exceeds $tolerance)"))
    end
    return frame
end

@inline _project(frame::_PlaneFrame, point) =
    (_dot3(point, frame.u), _dot3(point, frame.v))

@inline function _on_segment(a, b, point)
    orient2(a, b, point) == 0 || return false
    return min(a[1], b[1]) <= point[1] <= max(a[1], b[1]) &&
           min(a[2], b[2]) <= point[2] <= max(a[2], b[2])
end

@inline function _segments_intersect(a, b, c, d)
    o1 = orient2(a, b, c)
    o2 = orient2(a, b, d)
    o3 = orient2(c, d, a)
    o4 = orient2(c, d, b)
    return (o1 == 0 && _on_segment(a, b, c)) ||
           (o2 == 0 && _on_segment(a, b, d)) ||
           (o3 == 0 && _on_segment(c, d, a)) ||
           (o4 == 0 && _on_segment(c, d, b)) ||
           (o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && o1 != o2 && o3 != o4)
end

@inline _boundary_adjacent(i::Int, j::Int, count::Int) =
    j == i + 1 || (i == 1 && j == count)

function _adjacent_overlap(a, b, c, d)
    shared = a == c || a == d ? a : b == c || b == d ? b : nothing
    shared === nothing && return true
    for point in (a, b)
        point == shared || !_on_segment(c, d, point) || return true
    end
    for point in (c, d)
        point == shared || !_on_segment(a, b, point) || return true
    end
    return false
end

@inline _box_overlap(a, b) =
    a.xmin <= b.xmax && b.xmin <= a.xmax &&
    a.ymin <= b.ymax && b.ymin <= a.ymax

function _build_boundary_tree!(nodes, order, boxes, lo::Int, hi::Int)
    xmin = Inf
    xmax = -Inf
    ymin = Inf
    ymax = -Inf
    @inbounds for position in lo:hi
        box = boxes[order[position]]
        xmin = min(xmin, box.xmin)
        xmax = max(xmax, box.xmax)
        ymin = min(ymin, box.ymin)
        ymax = max(ymax, box.ymax)
    end
    count = hi - lo + 1
    node_index = length(nodes) + 1
    push!(nodes, _BoundaryNode(xmin, xmax, ymin, ymax, 0, 0, 0, count))
    if lo == hi
        @inbounds segment = boxes[order[lo]].segment
        nodes[node_index] = _BoundaryNode(
            xmin, xmax, ymin, ymax, segment, 0, 0, 1)
        return node_index
    end
    xspan = xmax - xmin
    yspan = ymax - ymin
    if xspan >= yspan
        sort!(@view(order[lo:hi]);
              by=index -> begin
                  box = boxes[index]
                  box.xmin / 2 + box.xmax / 2
              end,
              alg=QuickSort)
    else
        sort!(@view(order[lo:hi]);
              by=index -> begin
                  box = boxes[index]
                  box.ymin / 2 + box.ymax / 2
              end,
              alg=QuickSort)
    end
    middle = lo + (hi - lo) ÷ 2
    left = _build_boundary_tree!(nodes, order, boxes, lo, middle)
    right = _build_boundary_tree!(nodes, order, boxes, middle + 1, hi)
    nodes[node_index] = _BoundaryNode(
        xmin, xmax, ymin, ymax, 0, left, right, count)
    return node_index
end

function _audit_boundary_pair(points, i::Int, j::Int, count::Int)
    a = points[i]
    b = points[mod1(i + 1, count)]
    c = points[j]
    d = points[mod1(j + 1, count)]
    _segments_intersect(a, b, c, d) || return nothing
    adjacent = _boundary_adjacent(min(i, j), max(i, j), count)
    if !adjacent || _adjacent_overlap(a, b, c, d)
        throw(ArgumentError(
            "$_CALLER: boundary segments $i and $j intersect"))
    end
    return nothing
end

function _validate_simple_boundary(points)
    count = length(points)
    boxes = Vector{_SegmentBox}(undef, count)
    @inbounds for index in 1:count
        point = points[index]
        next = points[mod1(index + 1, count)]
        (isfinite(point[1]) && isfinite(point[2])) || throw(ArgumentError(
            "$_CALLER: projected boundary node $index is not finite"))
        boxes[index] = _SegmentBox(
            min(point[1], next[1]), max(point[1], next[1]),
            min(point[2], next[2]), max(point[2], next[2]), index)
    end
    order = collect(1:count)
    nodes = _BoundaryNode[]
    node_capacity = _checked_add(
        _checked_mul(2, count, "boundary audit node"), -1,
        "boundary audit node")
    sizehint!(nodes, node_capacity)
    root = _build_boundary_tree!(nodes, order, boxes, 1, count)
    limit = max(_BOUNDARY_AUDIT_FLOOR,
                _checked_mul(_BOUNDARY_AUDIT_MULTIPLIER, count,
                             "boundary-intersection audit"))
    stack = Tuple{Int,Int}[(root, root)]
    visits = 0
    candidates = 0
    while !isempty(stack)
        first, second = pop!(stack)
        visits += 1
        visits <= limit || throw(ArgumentError(
            "$_CALLER: boundary intersection audit exceeded its bounded " *
            "traversal limit $limit"))
        node1 = nodes[first]
        node2 = nodes[second]
        _box_overlap(node1, node2) || continue
        if first == second
            node1.segment != 0 && continue
            push!(stack, (node1.left, node1.left),
                         (node1.left, node1.right),
                         (node1.right, node1.right))
        elseif node1.segment != 0 && node2.segment != 0
            candidates += 1
            candidates <= limit || throw(ArgumentError(
                "$_CALLER: boundary intersection audit exceeded its bounded " *
                "candidate limit $limit"))
            _audit_boundary_pair(points, node1.segment, node2.segment, count)
        elseif node2.segment != 0 ||
               (node1.segment == 0 && node1.count >= node2.count)
            push!(stack, (node1.left, second), (node1.right, second))
        else
            push!(stack, (first, node2.left), (first, node2.right))
        end
    end
    return nothing
end

function _averaged_chords(first, opposite)
    length(first) == length(opposite) || throw(ErrorException(
        "$_CALLER: internal opposing-side count invariant failed"))
    cumulative = zeros(Float64, length(first))
    total = 0.0
    @inbounds for index in 1:length(first)-1
        first_length = _distance(
            first[index + 1], first[index], "first-side interval $index")
        opposite_length = _distance(
            opposite[index + 1], opposite[index],
            "opposing-side interval $index")
        increment = 0.5first_length + 0.5opposite_length
        (isfinite(increment) && increment > 0) || throw(ArgumentError(
            "$_CALLER: averaged chord $index is not finite and positive"))
        next = total + increment
        (isfinite(next) && next > total) || throw(ArgumentError(
            "$_CALLER: averaged chord accumulation is not strictly representable"))
        total = next
        cumulative[index + 1] = total
    end
    _validate_chord_ratios(cumulative, total, "averaged")
    return cumulative, total
end

function _side_chords(side)
    cumulative = zeros(Float64, length(side))
    total = 0.0
    @inbounds for index in 1:length(side)-1
        increment = _distance(
            side[index + 1], side[index], "second-side interval $index")
        next = total + increment
        (isfinite(next) && next > total) || throw(ArgumentError(
            "$_CALLER: second-side chord accumulation is not strictly representable"))
        total = next
        cumulative[index + 1] = total
    end
    _validate_chord_ratios(cumulative, total, "second-side")
    return cumulative, total
end

function _validate_chord_ratios(cumulative, total, description::AbstractString)
    (isfinite(total) && total > 0) || throw(ArgumentError(
        "$_CALLER: $description chord total is not finite and positive"))
    previous = 0.0
    @inbounds for index in 2:length(cumulative)-1
        current = cumulative[index] / total
        (isfinite(current) && previous < current < 1.0) ||
            throw(ArgumentError(
                "$_CALLER: $description chord parameters are not strictly representable"))
        previous = current
    end
    return nothing
end

function _chord_sample(side, cumulative, total::Float64, parameter::Float64)
    0.0 < parameter < 1.0 || throw(ErrorException(
        "$_CALLER: internal side interpolation parameter invariant failed"))
    lo = 2
    hi = length(cumulative)
    while lo < hi
        middle = lo + (hi - lo) ÷ 2
        if parameter <= cumulative[middle] / total
            hi = middle
        else
            lo = middle + 1
        end
    end
    right = lo
    left = right - 1
    fraction = (parameter * total - cumulative[left]) /
               (cumulative[right] - cumulative[left])
    isfinite(fraction) || throw(ArgumentError(
        "$_CALLER: side interpolation fraction is not finite"))
    first = side[left]
    second = side[right]
    return ntuple(coordinate ->
        first[coordinate] + fraction * (second[coordinate] - first[coordinate]), 3)
end

# Gmsh 4.15.2 meshGFaceTransfinite.cpp lines 592-611 use j/i as the
# triangular ray coordinate, interpolate side 2 by chord length at that value,
# and apply TRAN_TRI. Corner 1 is the normalization origin, so its zero term is
# omitted from this algebraically identical form.
@inline function _transfinite_point(first, second, third, corner2, corner3,
                                    u::Float64, v::Float64)
    one_v = 1.0 - v
    return ntuple(3) do coordinate
        u * second[coordinate] + one_v * first[coordinate] +
        v * third[coordinate] -
        (u * one_v * corner2[coordinate] + u * v * corner3[coordinate])
    end
end

@inline function _physical_point(normalized, origin, scale)
    point = (origin[1] + scale * normalized[1],
             origin[2] + scale * normalized[2],
             origin[3] + scale * normalized[3])
    all(isfinite, point) || throw(ArgumentError(
        "$_CALLER: generated coordinate is not finite"))
    return point
end

@inline _node(i::Int, j::Int) = Int32((i * (i + 1)) ÷ 2 + j + 1)

function _fill_segments!(segments, tags, divisions::Int, side_tags)
    cursor = 0
    @inbounds for i in 0:divisions-1
        cursor += 1
        segments[1, cursor] = _node(i, 0)
        segments[2, cursor] = _node(i + 1, 0)
        tags[cursor] = side_tags[1]
    end
    @inbounds for j in 0:divisions-1
        cursor += 1
        segments[1, cursor] = _node(divisions, j)
        segments[2, cursor] = _node(divisions, j + 1)
        tags[cursor] = side_tags[2]
    end
    @inbounds for i in divisions:-1:1
        cursor += 1
        segments[1, cursor] = _node(i, i)
        segments[2, cursor] = _node(i - 1, i - 1)
        tags[cursor] = side_tags[3]
    end
    cursor == size(segments, 2) || throw(ErrorException(
        "$_CALLER: internal segment count invariant failed"))
    return nothing
end

function _fill_triangles!(triangles, divisions::Int)
    cursor = 0
    @inbounds for i in 0:divisions-1
        for j in 0:i
            v1 = _node(i, j)
            v2 = _node(i + 1, j)
            v3 = _node(i + 1, j + 1)
            if i > 0 && j < i
                v4 = _node(i, j + 1)
                cursor += 1
                triangles[1, cursor] = v1
                triangles[2, cursor] = v3
                triangles[3, cursor] = v4
            end
            cursor += 1
            triangles[1, cursor] = v1
            triangles[2, cursor] = v2
            triangles[3, cursor] = v3
        end
    end
    cursor == size(triangles, 2) || throw(ErrorException(
        "$_CALLER: internal triangle count invariant failed"))
    return nothing
end

@inline function _project_output(coords, node_id::Int32, origin, scale, frame)
    index = Int(node_id)
    normalized = ((coords[1, index] - origin[1]) / scale,
                  (coords[2, index] - origin[2]) / scale,
                  (coords[3, index] - origin[3]) / scale)
    all(isfinite, normalized) || throw(ArgumentError(
        "$_CALLER: generated coordinate cannot be projected"))
    return _project(frame, normalized)
end

function _validate_triangle_orientation(coords, triangles, origin, scale, frame)
    reference = 0
    @inbounds for triangle in axes(triangles, 2)
        a = _project_output(coords, triangles[1, triangle], origin, scale, frame)
        b = _project_output(coords, triangles[2, triangle], origin, scale, frame)
        c = _project_output(coords, triangles[3, triangle], origin, scale, frame)
        orientation = orient2(a, b, c)
        orientation != 0 || throw(ArgumentError(
            "$_CALLER: cell triangle $triangle is folded or degenerate"))
        if reference == 0
            reference = orientation
        elseif orientation != reference
            throw(ArgumentError(
                "$_CALLER: cell triangle $triangle reverses patch orientation"))
        end
    end
    return reference
end

function _validate_boundary_postcondition(mesh::Mesh)
    actual, max_incidence = boundary_edges(mesh.tris)
    expected_incidence = ntris(mesh) == 1 ? 1 : 2
    max_incidence == expected_incidence || throw(ErrorException(
        "$_CALLER: internal triangle incidence postcondition failed"))
    expected = Vector{NTuple{2,Int32}}(undef, nsegs(mesh))
    @inbounds for segment in 1:nsegs(mesh)
        expected[segment] = _edge_key(
            mesh.segs[1, segment], mesh.segs[2, segment])
    end
    sort!(actual)
    sort!(expected)
    actual == expected || throw(ErrorException(
        "$_CALLER: triangle boundary does not equal emitted segments"))
    return nothing
end

"""
    mesh_transfinite_triangle(side1, side2, side3;
        arrangement=:left, face_tag=0, side_tags=(0,0,0),
        max_nodes=10_000_000, max_triangles=20_000_000) -> Mesh

Construct a planar three-sided structured triangle patch using Gmsh 4.15.2's
specific `Mesh.TransfiniteTri = 1` algorithm. Each side is an
already-discretized vector of finite 3-D points, cyclically oriented as
`c1→c2`, `c2→c3`, and `c3→c1`. Adjacent endpoints must match exactly
after conversion to `Float64`, and all three sides must have the same node
count.

The accepted arrangements are `:left`, `:right`, `:alternate_left`, and
`:alternate_right`. They intentionally produce identical triangle meshes:
Gmsh's unrecombined specific three-sided path does not consult the arrangement
when creating triangles. Arrangement-dependent triangular/quadrilateral mixes
from Gmsh's recombined path are outside this function's simplex-only contract.

The returned mesh contains all three boundary segment chains. `face_tag` is
copied to every triangle and each entry of `side_tags` to the corresponding
chain. Counts, Int32 topology bounds, and caller limits are checked before
output allocation. The completed output is checked for finite coordinates,
consistent nonzero orientation, manifold topology, and exact boundary
conservation; invalid or folded inputs are rejected without an unstructured
fallback.

This operation does not discretize curves, apply curve laws, reproduce Gmsh's
legacy `Mesh.TransfiniteTri = 0` collapsed-grid algorithm, repair unequal side
counts, map a general CAD/ruled/spherical parameterization, recombine into
quadrangles, smooth, handle holes/periodic seams/embedded entities, apply size
or quality fields, or construct volumes. A boundary whose spatial intersection
audit exceeds its linear candidate budget is rejected.
"""
function mesh_transfinite_triangle(side1::AbstractVector,
                                   side2::AbstractVector,
                                   side3::AbstractVector;
                                   arrangement=:left,
                                   face_tag=0,
                                   side_tags=(0, 0, 0),
                                   max_nodes=_DEFAULT_MAX_NODES,
                                   max_triangles=_DEFAULT_MAX_TRIANGLES)::Mesh
    # Gmsh 4.15.2 ignores this value for unrecombined transfinite3 triangles;
    # validating it still prevents silent misspellings and preserves API parity.
    _arrangement(arrangement)
    node_limit = _limit(max_nodes, "max_nodes")
    triangle_limit = _limit(max_triangles, "max_triangles")
    side_tags isa Tuple && length(side_tags) == 3 || throw(ArgumentError(
        "$_CALLER: side_tags must be a three-integer tuple"))
    physical_side_tags = ntuple(
        index -> _tag(side_tags[index], "side_tags[$index]"), 3)
    physical_face_tag = _tag(face_tag, "face_tag")

    lengths = (length(side1), length(side2), length(side3))
    @inbounds for side in 1:3
        lengths[side] >= 2 || throw(ArgumentError(
            "$_CALLER: side $side needs at least two points"))
    end
    lengths[1] == lengths[2] == lengths[3] || throw(ArgumentError(
        "$_CALLER: all three sides need matching node counts; got " *
        "$(lengths[1]), $(lengths[2]), and $(lengths[3])"))

    divisions = lengths[1] - 1
    triangles = _checked_mul(divisions, divisions, "triangle")
    twice_nodes = _checked_mul(
        _checked_add(divisions, 1, "node factor"),
        _checked_add(divisions, 2, "node factor"), "twice-node")
    nodes = twice_nodes ÷ 2
    segments = _checked_mul(3, divisions, "segment")
    nodes <= _INT32_MAX || throw(ArgumentError(
        "$_CALLER: $nodes nodes exceed the Int32 indexing limit"))
    triangles <= _INT32_MAX || throw(ArgumentError(
        "$_CALLER: $triangles triangles exceed the Int32 topology limit"))
    segments <= _INT32_MAX || throw(ArgumentError(
        "$_CALLER: $segments segments exceed the Int32 topology limit"))
    nodes <= node_limit || throw(ArgumentError(
        "$_CALLER: $nodes nodes exceed max_nodes=$node_limit"))
    triangles <= triangle_limit || throw(ArgumentError(
        "$_CALLER: $triangles triangles exceed max_triangles=$triangle_limit"))

    sides = (_convert_side(side1, 1),
             _convert_side(side2, 2),
             _convert_side(side3, 3))
    @inbounds for side in 1:3
        _validate_side_edges(sides[side], side)
        next = mod1(side + 1, 3)
        sides[side][end] == sides[next][1] || throw(ArgumentError(
            "$_CALLER: side $side endpoint does not exactly match side $next start point"))
    end
    corners = (sides[1][1], sides[2][1], sides[3][1])
    length(Set(corners)) == 3 || throw(ArgumentError(
        "$_CALLER: the three corners must be distinct"))

    ring = _boundary_ring(sides)
    origin, scale = _normalization(ring)
    frame = _plane_frame(ring, origin, scale)
    projected_ring = NTuple{2,Float64}[
        _project(frame, _normalize(point, origin, scale)) for point in ring]
    _validate_simple_boundary(projected_ring)

    first = _normalized_side(sides[1], origin, scale)
    second = _normalized_side(sides[2], origin, scale)
    third = reverse(_normalized_side(sides[3], origin, scale))
    radial_chords, radial_total = _averaged_chords(first, third)
    second_chords, second_total = _side_chords(second)
    corner2 = first[end]
    corner3 = second[end]

    coordinates = Matrix{Float64}(undef, 3, nodes)
    @inbounds for i in 0:divisions
        u = radial_chords[i + 1] / radial_total
        for j in 0:i
            point = if j == 0
                sides[1][i + 1]
            elseif i == divisions
                sides[2][j + 1]
            elseif j == i
                sides[3][divisions - i + 1]
            else
                v = j / i
                side_point = _chord_sample(
                    second, second_chords, second_total, v)
                normalized = _transfinite_point(
                    first[i + 1], side_point, third[i + 1],
                    corner2, corner3, u, v)
                all(isfinite, normalized) || throw(ArgumentError(
                    "$_CALLER: transfinite interpolation generated a non-finite coordinate"))
                _physical_point(normalized, origin, scale)
            end
            node = Int(_node(i, j))
            coordinates[1, node] = point[1]
            coordinates[2, node] = point[2]
            coordinates[3, node] = point[3]
        end
    end

    segment_topology = Matrix{Int32}(undef, 2, segments)
    segment_tags = Vector{Int32}(undef, segments)
    _fill_segments!(
        segment_topology, segment_tags, divisions, physical_side_tags)
    triangle_topology = Matrix{Int32}(undef, 3, triangles)
    _fill_triangles!(triangle_topology, divisions)
    _validate_triangle_orientation(
        coordinates, triangle_topology, origin, scale, frame)
    triangle_tags = fill(physical_face_tag, triangles)

    mesh = Mesh(coordinates;
                segs=segment_topology,
                tris=triangle_topology,
                seg_tag=segment_tags,
                tri_tag=triangle_tags)
    diagnostic = validate(mesh)
    diagnostic.ok || throw(ErrorException(
        "$_CALLER: internal output validation failed — " *
        join(diagnostic.messages, "; ")))
    (nnodes(mesh) == nodes && nsegs(mesh) == segments &&
     ntris(mesh) == triangles) || throw(ErrorException(
        "$_CALLER: internal output count postcondition failed"))
    _validate_boundary_postcondition(mesh)
    return mesh
end

end # module TransfiniteTriangle
