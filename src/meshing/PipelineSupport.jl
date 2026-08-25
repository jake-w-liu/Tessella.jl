"""
    PipelineSupport

Internal input conversion, resource accounting, and postcondition certificates
shared by Tessella's top-level planar and volume-meshing pipelines.
"""
module PipelineSupport

using ..MeshTypes: Mesh
using ..Mesh1D: _convex_coordinate

const PIPELINE_DEFAULT_MAX_NODES = 10_000_000
const PIPELINE_DEFAULT_MAX_TETS = 60_000_000
const _PipelineDyadic = Rational{BigInt}

function pipeline_float(value, caller::AbstractString, name::AbstractString;
                        allow_positive_inf::Bool=false)
    value isa Real || throw(ArgumentError("$caller: $name must be real"))
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    converted = try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: $name must be Float64-representable: $(sprint(showerror, err))"))
    end
    if allow_positive_inf
        (isfinite(converted) || converted == Inf) || throw(ArgumentError(
            "$caller: $name must be finite or positive Inf"))
    else
        isfinite(converted) || throw(ArgumentError("$caller: $name must be finite"))
    end
    return converted
end

function pipeline_limit(value, caller::AbstractString, name::AbstractString,
                        maximum::Integer=typemax(Int))
    value isa Integer || throw(ArgumentError("$caller: $name must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    0 <= value <= maximum || throw(ArgumentError(
        "$caller: $name must lie in 0:$maximum"))
    return Int(value)
end

function pipeline_seed(value, caller::AbstractString, maximum::Integer)
    value isa Integer || throw(ArgumentError("$caller: rng_seed must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: rng_seed must not be Bool"))
    0 <= value <= maximum || throw(ArgumentError(
        "$caller: rng_seed must lie in 0:$maximum"))
    return value
end

function pipeline_checked_mul(caller::AbstractString, what::AbstractString,
                              values::Int...)
    result = 1
    for value in values
        result = try
            Base.checked_mul(result, value)
        catch err
            err isa InterruptException && rethrow()
            err isa OverflowError || rethrow()
            throw(ArgumentError("$caller: $what count overflows Int"))
        end
    end
    return result
end

function pipeline_prism_step(hmax::Float64)
    step = hmax / sqrt(2.0)
    while isfinite(step) && step > 0 && hypot(step, step) > hmax
        step = prevfloat(step)
    end
    (isfinite(step) && step > 0) || throw(ArgumentError(
        "mesh_sized_extrude: hmax is below Float64 spacing resolution"))
    return step
end

function pipeline_layer_count(z0::Float64, z1::Float64, step::Float64)
    ratio = (_PipelineDyadic(z1) - _PipelineDyadic(z0)) /
            _PipelineDyadic(step)
    count = max(big(1), cld(numerator(ratio), denominator(ratio)))
    count <= typemax(Int) || throw(ArgumentError(
        "mesh_sized_extrude: axial layer count exceeds the platform Int limit"))
    return Int(count)
end

function pipeline_z_levels(z0::Float64, z1::Float64, count::Int,
                           edge_limit::Float64)
    levels = Vector{Float64}(undef, count + 1)
    levels[1] = z0
    @inbounds for k in 1:count-1
        value = _convex_coordinate(z0, z1, k / count)
        (isfinite(value) && levels[k] < value < z1) || throw(ArgumentError(
            "mesh_sized_extrude: axial level $k of $count is not a " *
            "distinct interior Float64 coordinate"))
        gap = value - levels[k]
        (isfinite(gap) && gap <= edge_limit) || throw(ArgumentError(
            "mesh_sized_extrude: represented axial layer $k exceeds hmax " *
            "at Float64 resolution"))
        levels[k + 1] = value
    end
    final_gap = z1 - levels[count]
    (isfinite(final_gap) && final_gap > 0 && final_gap <= edge_limit) ||
        throw(ArgumentError(
            "mesh_sized_extrude: represented final axial layer exceeds " *
            "hmax at Float64 resolution"))
    levels[end] = z1
    return levels
end

function certify_pipeline_planar_edges(mesh::Mesh, limit::Float64)
    @inbounds for triangle in axes(mesh.tris, 2)
        vertices = (mesh.tris[1, triangle], mesh.tris[2, triangle],
                    mesh.tris[3, triangle])
        for (first, second) in ((1, 2), (2, 3), (3, 1))
            a = Int(vertices[first])
            b = Int(vertices[second])
            edge_length = hypot(mesh.coords[1, b] - mesh.coords[1, a],
                                mesh.coords[2, b] - mesh.coords[2, a])
            (isfinite(edge_length) && edge_length <= limit) ||
                throw(ErrorException(
                    "mesh_sized_extrude: transverse refinement postcondition " *
                    "failed on triangle $triangle edge ($a,$b)"))
        end
    end
    return nothing
end

function certify_pipeline_tet_edges(coords::Matrix{Float64},
                                    tets::Matrix{Int32}, limit::Float64)
    @inbounds for tet in axes(tets, 2), first in 1:3, second in first+1:4
        a = Int(tets[first, tet])
        b = Int(tets[second, tet])
        edge_length = hypot(coords[1, b] - coords[1, a],
                            coords[2, b] - coords[2, a],
                            coords[3, b] - coords[3, a])
        (isfinite(edge_length) && edge_length <= limit) ||
            throw(ErrorException(
                "mesh_sized_extrude: maximum-edge postcondition failed on " *
                "tetrahedron $tet edge ($a,$b)"))
    end
    return nothing
end

function pipeline_points2(xs::AbstractVector, ys::AbstractVector,
                          caller::AbstractString)
    count = length(xs)
    length(ys) == count || throw(ArgumentError(
        "$caller: xs and ys must have the same length " *
        "(got $count and $(length(ys)))"))
    count >= 3 || throw(ArgumentError("$caller: at least three points are required"))
    count <= typemax(Int32) || throw(ArgumentError(
        "$caller: $count input points exceed Int32 indexing"))
    converted_x = Vector{Float64}(undef, count)
    converted_y = Vector{Float64}(undef, count)
    cursor = 0
    for (raw_x, raw_y) in zip(xs, ys)
        cursor += 1
        cursor <= count || throw(ArgumentError(
            "$caller: coordinate iteration produced more values than length reported"))
        converted_x[cursor] = pipeline_float(
            raw_x, caller, "point $cursor x coordinate")
        converted_y[cursor] = pipeline_float(
            raw_y, caller, "point $cursor y coordinate")
    end
    cursor == count || throw(ArgumentError(
        "$caller: coordinate iteration ended after $cursor of $count points"))
    return converted_x, converted_y
end

function pipeline_segments(records::AbstractVector, point_count::Int,
                           caller::AbstractString)
    length(records) <= typemax(Int32) || throw(ArgumentError(
        "$caller: segment count exceeds the Int32 topology limit"))
    segments = Vector{Tuple{Int,Int}}(undef, length(records))
    cursor = 0
    for record in records
        cursor += 1
        cursor <= length(segments) || throw(ArgumentError(
            "$caller: segment iteration produced more records than length reported"))
        (record isa Tuple && length(record) == 2) || throw(ArgumentError(
            "$caller: segment $cursor must be a two-integer tuple"))
        first, second = record
        (first isa Integer && second isa Integer) || throw(ArgumentError(
            "$caller: segment $cursor endpoints must be integers"))
        (first isa Bool || second isa Bool) && throw(ArgumentError(
            "$caller: segment $cursor endpoints must not be Bool"))
        1 <= first <= point_count || throw(ArgumentError(
            "$caller: segment $cursor endpoint $first is outside 1:$point_count"))
        1 <= second <= point_count || throw(ArgumentError(
            "$caller: segment $cursor endpoint $second is outside 1:$point_count"))
        first != second || throw(ArgumentError(
            "$caller: segment $cursor has identical endpoints $first"))
        segments[cursor] = (Int(first), Int(second))
    end
    cursor == length(segments) || throw(ArgumentError(
        "$caller: segment iteration ended after $cursor of $(length(segments)) records"))
    return segments
end

end # module PipelineSupport
