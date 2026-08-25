"""
    Refine

Deterministic one-level uniform refinement of finalized linear simplex meshes.
Every referenced edge receives one shared midpoint; segments, triangles and
tetrahedra are replaced by 2, 4 and 8 children, respectively.
"""
module Refine

using ..Predicates: orient3
import ..MeshTypes
using ..MeshTypes: Mesh

export refine_uniform

const _INT32_MAX = Int(typemax(Int32))
const _MIDPOINT_SMALL = 2floatmin(Float64)
const _MIDPOINT_LARGE = floatmax(Float64) / 2

@inline function _limit(value, name::AbstractString)
    (value isa Integer && !(value isa Bool)) || throw(ArgumentError(
        "refine_uniform: $name must be an integer other than Bool"))
    0 <= value <= typemax(Int32) || throw(ArgumentError(
        "refine_uniform: $name must lie in 0:$(typemax(Int32))"))
    return Int(value)
end

@inline function _checked_mul(x::Int, y::Int, what::AbstractString)
    try
        return Base.checked_mul(x, y)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "refine_uniform: $what count overflows the platform Int limit"))
    end
end

@inline function _checked_add(x::Int, y::Int, what::AbstractString)
    try
        return Base.checked_add(x, y)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "refine_uniform: $what count overflows the platform Int limit"))
    end
end

@inline _edge(a::Int32, b::Int32) = a < b ? (a, b) : (b, a)

@inline function _midpoint_coordinate(a::Float64, b::Float64)
    # Arrange the calculation so that no branch overflows and at most one
    # operation rounds. This avoids the double rounding in `a + (b-a)/2`, while
    # retaining contributions from subnormal endpoints next to huge endpoints.
    aa = abs(a)
    ab = abs(b)
    if aa <= _MIDPOINT_LARGE && ab <= _MIDPOINT_LARGE
        return (a + b) / 2
    elseif aa < _MIDPOINT_SMALL
        return a + b / 2
    elseif ab < _MIDPOINT_SMALL
        return a / 2 + b
    end
    return a / 2 + b / 2
end

@inline function _point(coords::Matrix{Float64}, node::Int32)
    i = Int(node)
    @inbounds return (coords[1, i], coords[2, i], coords[3, i])
end

function _write_positive_tet!(tetrahedra::Matrix{Int32}, column::Int,
                              vertices::NTuple{4,Int32},
                              coords::Matrix{Float64}, parent::Int,
                              child::Int)
    a, b, c, d = vertices
    orientation = orient3(_point(coords, a), _point(coords, b),
                          _point(coords, c), _point(coords, d))
    orientation == 0 && throw(ArgumentError(
        "refine_uniform: tetrahedron $parent child $child has zero exact volume"))
    # `orient3` follows Shewchuk's convention, opposite to MeshTypes' positive
    # signed-volume convention. Swap once if the raw template is negative.
    orientation > 0 && ((a, b) = (b, a))
    @inbounds begin
        tetrahedra[1, column] = a
        tetrahedra[2, column] = b
        tetrahedra[3, column] = c
        tetrahedra[4, column] = d
    end
    return nothing
end

function _mark_referenced!(used::BitVector, mesh::Mesh)
    @inbounds for cell in axes(mesh.segs, 2), local_node in 1:2
        used[Int(mesh.segs[local_node, cell])] = true
    end
    @inbounds for cell in axes(mesh.tris, 2), local_node in 1:3
        used[Int(mesh.tris[local_node, cell])] = true
    end
    @inbounds for cell in axes(mesh.tets, 2), local_node in 1:4
        used[Int(mesh.tets[local_node, cell])] = true
    end
    return nothing
end

function _compact_map(used::BitVector)
    remap = zeros(Int32, length(used))
    count = 0
    @inbounds for old in eachindex(used)
        used[old] || continue
        count += 1
        remap[old] = Int32(count)
    end
    return remap, count
end

@inline function _remapped(remap::Vector{Int32}, node::Int32)
    @inbounds return remap[Int(node)]
end

function _edge_records(mesh::Mesh, remap::Vector{Int32}, count::Int)
    records = Vector{NTuple{2,Int32}}(undef, count)
    record = 0
    @inbounds for cell in axes(mesh.segs, 2)
        a = _remapped(remap, mesh.segs[1, cell])
        b = _remapped(remap, mesh.segs[2, cell])
        record += 1
        records[record] = _edge(a, b)
    end
    @inbounds for cell in axes(mesh.tris, 2)
        a = _remapped(remap, mesh.tris[1, cell])
        b = _remapped(remap, mesh.tris[2, cell])
        c = _remapped(remap, mesh.tris[3, cell])
        for edge in (_edge(a, b), _edge(b, c), _edge(c, a))
            record += 1
            records[record] = edge
        end
    end
    @inbounds for cell in axes(mesh.tets, 2)
        a = _remapped(remap, mesh.tets[1, cell])
        b = _remapped(remap, mesh.tets[2, cell])
        c = _remapped(remap, mesh.tets[3, cell])
        d = _remapped(remap, mesh.tets[4, cell])
        for edge in (_edge(a, b), _edge(a, c), _edge(a, d),
                     _edge(b, c), _edge(b, d), _edge(c, d))
            record += 1
            records[record] = edge
        end
    end
    record == count || error("refine_uniform: internal edge-record count mismatch")
    sort!(records; alg=QuickSort)
    isempty(records) && return records

    unique_count = 1
    @inbounds for read in 2:length(records)
        records[read] == records[unique_count] && continue
        unique_count += 1
        records[unique_count] = records[read]
    end
    resize!(records, unique_count)
    return records
end

function _coordinates_and_midpoints(mesh::Mesh, used::BitVector,
                                    remap::Vector{Int32}, compact_nodes::Int,
                                    edges::Vector{NTuple{2,Int32}}, final_nodes::Int)
    coords = Matrix{Float64}(undef, 3, final_nodes)
    @inbounds for old in eachindex(used)
        used[old] || continue
        new = Int(remap[old])
        coords[1, new] = mesh.coords[1, old]
        coords[2, new] = mesh.coords[2, old]
        coords[3, new] = mesh.coords[3, old]
    end

    midpoint_ids = Dict{NTuple{2,Int32},Int32}()
    sizehint!(midpoint_ids, length(edges))
    @inbounds for (index, edge) in pairs(edges)
        a, b = edge
        pa = _point(coords, a)
        pb = _point(coords, b)
        midpoint = (_midpoint_coordinate(pa[1], pb[1]),
                    _midpoint_coordinate(pa[2], pb[2]),
                    _midpoint_coordinate(pa[3], pb[3]))
        all(isfinite, midpoint) || throw(ArgumentError(
            "refine_uniform: edge $edge has a non-finite midpoint"))
        (midpoint != pa && midpoint != pb) || throw(ArgumentError(
            "refine_uniform: edge $edge midpoint is below Float64 coordinate resolution"))
        node = compact_nodes + index
        coords[1, node] = midpoint[1]
        coords[2, node] = midpoint[2]
        coords[3, node] = midpoint[3]
        midpoint_ids[edge] = Int32(node)
    end
    return coords, midpoint_ids
end

@inline function _midpoint(midpoint_ids::Dict{NTuple{2,Int32},Int32},
                           a::Int32, b::Int32)
    return midpoint_ids[_edge(a, b)]
end

function _refined_segments(mesh::Mesh, remap::Vector{Int32},
                           midpoint_ids::Dict{NTuple{2,Int32},Int32}, count::Int)
    segments = Matrix{Int32}(undef, 2, count)
    tags = Vector{Int32}(undef, count)
    @inbounds for parent in axes(mesh.segs, 2)
        a = _remapped(remap, mesh.segs[1, parent])
        b = _remapped(remap, mesh.segs[2, parent])
        midpoint = _midpoint(midpoint_ids, a, b)
        first = 2parent - 1
        segments[1, first] = a
        segments[2, first] = midpoint
        segments[1, first + 1] = midpoint
        segments[2, first + 1] = b
        tags[first] = mesh.seg_tag[parent]
        tags[first + 1] = mesh.seg_tag[parent]
    end
    return segments, tags
end

function _refined_triangles(mesh::Mesh, remap::Vector{Int32},
                            midpoint_ids::Dict{NTuple{2,Int32},Int32}, count::Int)
    triangles = Matrix{Int32}(undef, 3, count)
    tags = Vector{Int32}(undef, count)
    @inbounds for parent in axes(mesh.tris, 2)
        a = _remapped(remap, mesh.tris[1, parent])
        b = _remapped(remap, mesh.tris[2, parent])
        c = _remapped(remap, mesh.tris[3, parent])
        mab = _midpoint(midpoint_ids, a, b)
        mbc = _midpoint(midpoint_ids, b, c)
        mca = _midpoint(midpoint_ids, c, a)
        children = ((a, mab, mca),
                    (mab, mbc, mca),
                    (mab, b, mbc),
                    (mca, mbc, c))
        base = 4(parent - 1)
        for child in 1:4
            column = base + child
            triangle = children[child]
            triangles[1, column] = triangle[1]
            triangles[2, column] = triangle[2]
            triangles[3, column] = triangle[3]
            tags[column] = mesh.tri_tag[parent]
        end
    end
    return triangles, tags
end

function _refined_tetrahedra(mesh::Mesh, remap::Vector{Int32},
                             midpoint_ids::Dict{NTuple{2,Int32},Int32},
                             coords::Matrix{Float64}, count::Int)
    tetrahedra = Matrix{Int32}(undef, 4, count)
    tags = Vector{Int32}(undef, count)
    @inbounds for parent in axes(mesh.tets, 2)
        a = _remapped(remap, mesh.tets[1, parent])
        b = _remapped(remap, mesh.tets[2, parent])
        c = _remapped(remap, mesh.tets[3, parent])
        d = _remapped(remap, mesh.tets[4, parent])
        mab = _midpoint(midpoint_ids, a, b)
        mbc = _midpoint(midpoint_ids, b, c)
        mca = _midpoint(midpoint_ids, c, a)
        mad = _midpoint(midpoint_ids, a, d)
        mbd = _midpoint(midpoint_ids, b, d)
        mcd = _midpoint(midpoint_ids, c, d)
        children = ((a, mab, mca, mad),
                    (mab, b, mbc, mbd),
                    (mca, mbc, c, mcd),
                    (mad, mbd, mcd, d),
                    (mab, mca, mad, mbd),
                    (mab, mbd, mbc, mca),
                    (mca, mad, mbd, mcd),
                    (mca, mcd, mbd, mbc))
        base = 8(parent - 1)
        for child in 1:8
            column = base + child
            _write_positive_tet!(tetrahedra, column, children[child], coords,
                                 parent, child)
            tags[column] = mesh.tet_tag[parent]
        end
    end
    return tetrahedra, tags
end

"""
    refine_uniform(mesh::Mesh;
                   max_nodes::Integer=typemax(Int32),
                   max_cells::Integer=typemax(Int32)) -> Mesh

Uniformly refine every linear simplex in `mesh` once. Unreferenced input nodes
are removed first. Referenced original nodes retain their relative order, and
one shared midpoint is appended for every undirected edge in lexicographic edge
order. Parent segments, triangles and tetrahedra produce 2, 4 and 8 children,
with the parent's tag copied to every child.

The input and result must satisfy [`MeshTypes.validate`](@ref). Output resource
counts are checked before output allocation: `max_nodes` bounds the final node
count and `max_cells` bounds the combined segment, triangle and tetrahedron
count. Both limits must lie in `0:typemax(Int32)` and must not be `Bool`.
"""
function refine_uniform(mesh::Mesh;
                        max_nodes=typemax(Int32),
                        max_cells=typemax(Int32))::Mesh
    node_limit = _limit(max_nodes, "max_nodes")
    cell_limit = _limit(max_cells, "max_cells")

    diagnostic = MeshTypes.validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "refine_uniform: input mesh is invalid — " * join(diagnostic.messages, "; ")))

    input_segments = size(mesh.segs, 2)
    input_triangles = size(mesh.tris, 2)
    input_tetrahedra = size(mesh.tets, 2)
    output_segments = _checked_mul(2, input_segments, "output segment")
    output_triangles = _checked_mul(4, input_triangles, "output triangle")
    output_tetrahedra = _checked_mul(8, input_tetrahedra, "output tetrahedron")
    output_segments <= _INT32_MAX || throw(ArgumentError(
        "refine_uniform: output segment count exceeds the Int32 topology limit"))
    output_triangles <= _INT32_MAX || throw(ArgumentError(
        "refine_uniform: output triangle count exceeds the Int32 topology limit"))
    output_tetrahedra <= _INT32_MAX || throw(ArgumentError(
        "refine_uniform: output tetrahedron count exceeds the Int32 topology limit"))
    output_cells = _checked_add(
        output_segments,
        _checked_add(output_triangles, output_tetrahedra, "output cell"),
        "output cell")
    output_cells <= _INT32_MAX || throw(ArgumentError(
        "refine_uniform: combined output cell count exceeds the Int32 limit"))
    output_cells <= cell_limit || throw(ArgumentError(
        "refine_uniform: output requires $output_cells cells, exceeding max_cells=$cell_limit"))

    edge_records = _checked_add(
        input_segments,
        _checked_add(_checked_mul(3, input_triangles, "edge-record"),
                     _checked_mul(6, input_tetrahedra, "edge-record"),
                     "edge-record"),
        "edge-record")
    # Preflight every dense output-array extent before any dense output array is
    # allocated. Allocation failures themselves are deliberately not caught.
    _checked_mul(edge_records, sizeof(NTuple{2,Int32}), "edge-record byte")
    segment_entries = _checked_mul(2, output_segments, "segment connectivity entry")
    triangle_entries = _checked_mul(3, output_triangles, "triangle connectivity entry")
    tetrahedron_entries = _checked_mul(4, output_tetrahedra,
                                       "tetrahedron connectivity entry")
    _checked_mul(segment_entries, sizeof(Int32), "segment connectivity byte")
    _checked_mul(triangle_entries, sizeof(Int32), "triangle connectivity byte")
    _checked_mul(tetrahedron_entries, sizeof(Int32), "tetrahedron connectivity byte")
    _checked_mul(output_segments, sizeof(Int32), "segment-tag byte")
    _checked_mul(output_triangles, sizeof(Int32), "triangle-tag byte")
    _checked_mul(output_tetrahedra, sizeof(Int32), "tetrahedron-tag byte")
    _checked_mul(size(mesh.coords, 2), sizeof(Int32), "node-remap byte")

    used = falses(size(mesh.coords, 2))
    _mark_referenced!(used, mesh)
    remap, compact_nodes = _compact_map(used)
    compact_nodes <= node_limit || throw(ArgumentError(
        "refine_uniform: $compact_nodes referenced nodes exceed max_nodes=$node_limit"))
    minimum_midpoints = edge_records == 0 ? 0 : 1
    minimum_nodes = _checked_add(compact_nodes, minimum_midpoints, "output node")
    minimum_nodes <= _INT32_MAX || throw(ArgumentError(
        "refine_uniform: output node count exceeds Int32 indexing"))
    minimum_nodes <= node_limit || throw(ArgumentError(
        "refine_uniform: output requires more than max_nodes=$node_limit"))

    edges = _edge_records(mesh, remap, edge_records)
    final_nodes = _checked_add(compact_nodes, length(edges), "output node")
    final_nodes <= _INT32_MAX || throw(ArgumentError(
        "refine_uniform: output node count exceeds Int32 indexing"))
    final_nodes <= node_limit || throw(ArgumentError(
        "refine_uniform: output requires $final_nodes nodes, exceeding max_nodes=$node_limit"))
    coordinate_entries = _checked_mul(3, final_nodes, "coordinate entry")
    _checked_mul(coordinate_entries, sizeof(Float64), "coordinate byte")

    coords, midpoint_ids = _coordinates_and_midpoints(
        mesh, used, remap, compact_nodes, edges, final_nodes)
    segments, segment_tags = _refined_segments(
        mesh, remap, midpoint_ids, output_segments)
    triangles, triangle_tags = _refined_triangles(
        mesh, remap, midpoint_ids, output_triangles)
    tetrahedra, tetrahedron_tags = _refined_tetrahedra(
        mesh, remap, midpoint_ids, coords, output_tetrahedra)

    result = Mesh(coords; segs=segments, tris=triangles, tets=tetrahedra,
                  seg_tag=segment_tags, tri_tag=triangle_tags,
                  tet_tag=tetrahedron_tags)
    output_diagnostic = MeshTypes.validate(result)
    output_diagnostic.ok || throw(ArgumentError(
        "refine_uniform: refinement produced an invalid mesh — " *
        join(output_diagnostic.messages, "; ")))
    return result
end

end # module Refine
