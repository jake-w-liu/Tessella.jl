"""
    StructuredRecombine

Shared emission and certification machinery for the recombined transfinite
volume kernels (six-face `TransfiniteVolume`, five-face `TransfinitePrism`).
Every routine takes the calling kernel's name as `caller` so diagnostics keep
their public entry-point prefix. The helpers mirror the verification model
upstream Gmsh documents implicitly: emitted hex/prism cells carry a positional
tetrahedral *shadow* decomposition — each shadow tetrahedron certified
strictly positively oriented — and the emitted boundary triangles and
quadrangles are audited against the shadow's exterior faces for exact
coverage. Gmsh is used only as a differential oracle; this module is
production code.
"""
module StructuredRecombine

using ..MeshTypes: boundary_faces, validate, _throw_simplex_validation
using ..Predicates: orient3
using ..Elements: ElementBlock, MixedMesh

@inline function _node(coords, node_id::Int32)
    index = Int(node_id)
    return (coords[1, index], coords[2, index], coords[3, index])
end

@inline _canon3(a::Int32, b::Int32, c::Int32) =
    a <= b ? (a <= c ? (b <= c ? (a, b, c) : (a, c, b)) : (c, a, b)) :
             (b <= c ? (a <= c ? (b, a, c) : (b, c, a)) : (c, b, a))

@inline function _canon4(a::Int32, b::Int32, c::Int32, d::Int32)
    a, b = minmax(a, b); c, d = minmax(c, d)
    a, c = minmax(a, c); b, d = minmax(b, d)
    b, c = minmax(b, c)
    return (a, b, c, d)
end

function _canonical_tets(tets::AbstractMatrix{Int32})
    result = Vector{NTuple{4,Int32}}(undef, size(tets, 2))
    @inbounds for tet in axes(tets, 2)
        result[tet] = _canon4(
            tets[1, tet], tets[2, tet], tets[3, tet], tets[4, tet])
    end
    return sort!(result)
end

function _canonical_triangles(tris::AbstractMatrix{Int32})
    result = Vector{NTuple{3,Int32}}(undef, size(tris, 2))
    @inbounds for triangle in axes(tris, 2)
        a = tris[1, triangle]
        b = tris[2, triangle]
        c = tris[3, triangle]
        result[triangle] = _canon3(a, b, c)
    end
    sort!(result)
    return result
end

# Strict analogue of the simplex path's auto-correcting emitter: the fixed
# shadow templates must reproduce the canonical orientation, otherwise the
# represented grid is folded and the recombined cell cannot be certified.
function _emit_canonical_tet!(caller::AbstractString, tets, position::Int,
                              coords, a::Int32, b::Int32, c::Int32, d::Int32)
    sign = orient3(_node(coords, a), _node(coords, b),
                   _node(coords, c), _node(coords, d))
    sign == 0 && throw(ArgumentError(
        "$caller: interpolation produced a zero-volume tetrahedron " *
        "at output position $position"))
    sign < 0 || throw(ArgumentError(
        "$caller: interpolation reversed canonical tetrahedron $position; " *
        "the represented grid is folded"))
    @inbounds begin
        tets[1, position] = a
        tets[2, position] = b
        tets[3, position] = c
        tets[4, position] = d
    end
    return nothing
end

function _write_outward_triangle!(caller::AbstractString, tris, tags,
                                  position::Int, coords,
                                  a::Int32, b::Int32, c::Int32,
                                  opposite, tag::Int32)
    sign = orient3(_node(coords, a), _node(coords, b), _node(coords, c), opposite)
    sign == 0 && throw(ArgumentError(
        "$caller: interpolation produced a degenerate boundary triangle " *
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

# Quadrilateral analogue of `_write_outward_triangle!`: emits the boundary
# quadrangle wound so its first corner triangle faces away from `opposite`.
# A degenerate first triangle (non-planar warped quad) keeps the incoming
# winding — the coverage audit still certifies the cell's placement.
function _write_outward_quad!(quads, tags, position::Int, coords,
                              a::Int32, b::Int32, c::Int32, d::Int32,
                              opposite, tag::Int32)
    sign = orient3(_node(coords, a), _node(coords, b),
                   _node(coords, c), opposite)
    @inbounds begin
        quads[1, position] = a
        if sign >= 0
            quads[2, position] = b
            quads[3, position] = c
            quads[4, position] = d
        else
            quads[2, position] = d
            quads[3, position] = c
            quads[4, position] = b
        end
        tags[position] = tag
    end
    return position + 1
end

# Writes the cell tuple into `cells` at `column` and emits its tetrahedron
# decomposition into `shadow`, certifying every constituent tetrahedron
# through the canonical orientation check. Returns the updated shadow
# position.
@inline function _emit_recombined_cell!(caller::AbstractString,
        cells::AbstractMatrix{Int32}, column::Int,
        shadow::AbstractMatrix{Int32}, position::Int,
        coords, vertices, decomp)
    @inbounds for i in eachindex(vertices)
        cells[i, column] = vertices[i]
    end
    @inbounds for tet in decomp
        position += 1
        _emit_canonical_tet!(caller, shadow, position, coords,
                             vertices[tet[1]], vertices[tet[2]],
                             vertices[tet[3]], vertices[tet[4]])
    end
    return position
end

# Audits the emitted boundary cells against the shadow partition's boundary:
# every emitted triangle must be a shadow boundary face and every emitted
# quadrangle must contain exactly the two boundary triangles it covers —
# once each, with nothing left over. This certifies exact coverage even for
# shifted quadrangle splits, whose covered triangle pairs are always subsets
# of the quad's vertex set.
function _certify_recombined_boundary(caller::AbstractString, shadow, tris,
                                      quads)
    extracted_boundary, maximum_incidence = boundary_faces(shadow)
    maximum_incidence == 2 || throw(ErrorException(
        "$caller: recombined volume decomposition produced face incidence " *
        "$maximum_incidence"))
    remaining = Set{NTuple{3,Int32}}()
    sizehint!(remaining, length(extracted_boundary))
    for face in extracted_boundary
        push!(remaining, face)
    end
    @inbounds for t in axes(tris, 2)
        key = _canon3(tris[1, t], tris[2, t], tris[3, t])
        key in remaining || throw(ErrorException(
            "$caller: emitted boundary triangle $t is not on the volume " *
            "boundary"))
        delete!(remaining, key)
    end
    @inbounds for q in axes(quads, 2)
        a = quads[1, q]; b = quads[2, q]; c = quads[3, q]; d = quads[4, q]
        covered = 0
        for face in (_canon3(a, b, c), _canon3(a, c, d),
                     _canon3(a, b, d), _canon3(b, c, d))
            if face in remaining
                delete!(remaining, face)
                covered += 1
            end
        end
        covered == 2 || throw(ErrorException(
            "$caller: emitted boundary quadrangle $q does not cover " *
            "exactly two boundary faces"))
    end
    isempty(remaining) || throw(ErrorException(
        "$caller: emitted boundary leaves $(length(remaining)) volume " *
        "boundary faces uncovered"))
    return nothing
end

# Assembles the recombined mixed mesh: only nonempty element blocks are
# materialized, volume cells carry `volume_tag`, and the result is validated
# before returning.
function _recombined_mixed_mesh(caller::AbstractString, coords,
                                tris, tri_tags, quads, quad_tags,
                                hexes, prisms, volume_tag::Int32)
    blocks = ElementBlock[]
    size(tris, 2) > 0 &&
        push!(blocks, ElementBlock(2, tris, tri_tags))
    size(quads, 2) > 0 &&
        push!(blocks, ElementBlock(3, quads, quad_tags))
    size(hexes, 2) > 0 &&
        push!(blocks, ElementBlock(5, hexes,
                                   fill(volume_tag, size(hexes, 2))))
    size(prisms, 2) > 0 &&
        push!(blocks, ElementBlock(6, prisms,
                                   fill(volume_tag, size(prisms, 2))))
    mesh = MixedMesh(coords, blocks)
    diagnostic = validate(mesh)
    diagnostic.ok || _throw_simplex_validation(caller, diagnostic.messages)
    return mesh
end

end
