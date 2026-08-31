"""
    MeshEntityTopology

Deterministic global edge and triangular-face identifiers for finalized
linear-simplex [`Mesh`](@ref) values. Numbering follows first encounter in
segment, triangle, then tetrahedron order, using Gmsh's local simplex patterns.
Lookups return detached Gmsh-shaped tag and orientation arrays.
"""
module MeshEntityTopology

using ..MeshTypes: Mesh, nnodes

const _SEGMENT_EDGES=((1,2),)
const _TRIANGLE_EDGES=((1,2),(2,3),(3,1))
const _TETRAHEDRON_EDGES=((1,2),(2,3),(3,1),(4,1),(4,3),(4,2))
const _TRIANGLE_FACES=((1,2,3),)
const _TETRAHEDRON_FACES=((1,3,2),(1,2,4),(1,4,3),(4,2,3))

@inline function _simplex_edge_patterns(element_type::Int)
    element_type==1 && return _SEGMENT_EDGES
    element_type==2 && return _TRIANGLE_EDGES
    element_type==4 && return _TETRAHEDRON_EDGES
    return ()
end

@inline function _simplex_face_patterns(element_type::Int,face_type::Int)
    face_type==3 || return ()
    element_type==2 && return _TRIANGLE_FACES
    element_type==4 && return _TETRAHEDRON_FACES
    return ()
end

struct MeshEdgeTopology
    node_count::Int
    nodes::Vector{NTuple{2,Int32}}
    tags::Dict{NTuple{2,Int32},UInt64}
end

struct MeshFaceTopology
    node_count::Int
    triangle_nodes::Vector{NTuple{3,Int32}}
    triangle_tags::Dict{NTuple{3,Int32},UInt64}
end

@inline _edge_key(a::Int32,b::Int32)=a<b ? (a,b) : (b,a)

@inline function _face_key(a::Int32,b::Int32,c::Int32)
    a>b && ((a,b)=(b,a))
    b>c && ((b,c)=(c,b))
    a>b && ((a,b)=(b,a))
    return (a,b,c)
end

function _add_edge!(nodes,tags,a::Int32,b::Int32,
                    element_type::Int,cell::Int)
    a!=b || throw(ArgumentError(
        "mesh_edge_topology: type-$element_type cell $cell has a repeated " *
        "edge node $a"))
    key=_edge_key(a,b)
    haskey(tags,key) && return nothing
    push!(nodes,(a,b))
    tags[key]=UInt64(length(nodes))
    return nothing
end

function _add_face!(nodes,tags,a::Int32,b::Int32,c::Int32,
                    element_type::Int,cell::Int)
    (a!=b && a!=c && b!=c) || throw(ArgumentError(
        "mesh_face_topology: type-$element_type cell $cell has repeated " *
        "face nodes ($a, $b, $c)"))
    key=_face_key(a,b,c)
    haskey(tags,key) && return nothing
    push!(nodes,(a,b,c))
    tags[key]=UInt64(length(nodes))
    return nothing
end

function _checked_edge_candidate_count(mesh::Mesh)
    try
        triangle_count=Base.checked_mul(3,size(mesh.tris,2))
        tetrahedron_count=Base.checked_mul(6,size(mesh.tets,2))
        return Base.checked_add(
            size(mesh.segs,2),Base.checked_add(
                triangle_count,tetrahedron_count))
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "mesh_edge_topology: candidate count exceeds the platform Int range"))
    end
end

function _checked_face_candidate_count(mesh::Mesh)
    try
        tetrahedron_count=Base.checked_mul(4,size(mesh.tets,2))
        return Base.checked_add(size(mesh.tris,2),tetrahedron_count)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "mesh_face_topology: candidate count exceeds the platform Int range"))
    end
end

function _mesh_edge_topology(mesh::Mesh)
    candidate_count=_checked_edge_candidate_count(mesh)
    nodes=NTuple{2,Int32}[]
    tags=Dict{NTuple{2,Int32},UInt64}()
    sizehint!(nodes,candidate_count)
    sizehint!(tags,candidate_count)
    for (element_type,cells) in
        ((1,mesh.segs),(2,mesh.tris),(4,mesh.tets))
        patterns=_simplex_edge_patterns(element_type)
        @inbounds for cell in axes(cells,2),pattern in patterns
            _add_edge!(nodes,tags,cells[pattern[1],cell],
                       cells[pattern[2],cell],element_type,cell)
        end
    end
    return MeshEdgeTopology(nnodes(mesh),nodes,tags)
end

function _mesh_face_topology(mesh::Mesh)
    candidate_count=_checked_face_candidate_count(mesh)
    nodes=NTuple{3,Int32}[]
    tags=Dict{NTuple{3,Int32},UInt64}()
    sizehint!(nodes,candidate_count)
    sizehint!(tags,candidate_count)
    for (element_type,cells) in ((2,mesh.tris),(4,mesh.tets))
        patterns=_simplex_face_patterns(element_type,3)
        @inbounds for cell in axes(cells,2),pattern in patterns
            _add_face!(nodes,tags,cells[pattern[1],cell],
                       cells[pattern[2],cell],cells[pattern[3],cell],
                       element_type,cell)
        end
    end
    return MeshFaceTopology(nnodes(mesh),nodes,tags)
end

function _checked_node_sequence(values,stride::Int,node_count::Int,
                                caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: node_tags must be a vector or tuple of integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    count=length(values)
    count%stride==0 || throw(ArgumentError(
        "$caller: node_tags length must be divisible by $stride; got $count"))
    result=Vector{Int32}(undef,count)
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: node_tags[$index] must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: node_tags[$index] must not be Bool"))
        tag=try
            Int(value)
        catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: node_tags[$index] exceeds the platform Int range"))
        end
        1<=tag<=node_count || throw(ArgumentError(
            "$caller: unknown node tag $value; expected a dense tag in " *
            "1:$node_count"))
        result[index]=Int32(tag)
    end
    return result
end

function _checked_face_type(value,caller::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$caller: face_type must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: face_type must not be Bool"))
    face_type=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError ||
         err isa MethodError) || rethrow()
        throw(ArgumentError(
            "$caller: face_type exceeds the platform Int range"))
    end
    face_type in (3,4) || throw(ArgumentError(
        "$caller: face_type must be 3 (triangle) or 4 (quadrangle)"))
    return face_type
end

function _mesh_edges(topology::Union{Nothing,MeshEdgeTopology},
                     mesh::Mesh,node_tags,
                     caller::AbstractString="mesh_edges")
    nodes=_checked_node_sequence(node_tags,2,nnodes(mesh),caller)
    isempty(nodes) && return UInt64[],Int32[]
    topology===nothing && throw(ArgumentError(
        "$caller: global edges have not been created; call create_edges first"))
    topology.node_count==nnodes(mesh) || error(
        "$caller: internal edge topology does not match the cached mesh")
    count=length(nodes)÷2
    tags=Vector{UInt64}(undef,count)
    orientations=Vector{Int32}(undef,count)
    @inbounds for edge in 1:count
        first_node=nodes[2edge-1]
        second_node=nodes[2edge]
        first_node!=second_node || throw(ArgumentError(
            "$caller: edge $edge repeats node $first_node"))
        key=_edge_key(first_node,second_node)
        tag=get(topology.tags,key,UInt64(0))
        tag!=0 || throw(ArgumentError(
            "$caller: unknown mesh edge ($first_node, $second_node)"))
        tags[edge]=tag
        orientations[edge]=first_node<second_node ? Int32(1) : Int32(-1)
    end
    return tags,orientations
end

function _mesh_faces(topology::Union{Nothing,MeshFaceTopology},
                     mesh::Mesh,face_type_value,node_tags,
                     caller::AbstractString="mesh_faces")
    face_type=_checked_face_type(face_type_value,caller)
    nodes=_checked_node_sequence(node_tags,face_type,nnodes(mesh),caller)
    isempty(nodes) && return UInt64[],Int32[]
    topology===nothing && throw(ArgumentError(
        "$caller: global faces have not been created; call create_faces first"))
    topology.node_count==nnodes(mesh) || error(
        "$caller: internal face topology does not match the cached mesh")
    count=length(nodes)÷face_type
    if face_type==4
        first_node=nodes[1]
        second_node=nodes[2]
        third_node=nodes[3]
        fourth_node=nodes[4]
        (first_node!=second_node && first_node!=third_node &&
         first_node!=fourth_node && second_node!=third_node &&
         second_node!=fourth_node && third_node!=fourth_node) ||
            throw(ArgumentError("$caller: face 1 contains repeated node tags"))
        throw(ArgumentError(
            "$caller: unknown quadrangular mesh face " *
            "($first_node, $second_node, $third_node, $fourth_node)"))
    end
    tags=Vector{UInt64}(undef,count)
    orientations=zeros(Int32,count)
    @inbounds for face in 1:count
        offset=face_type*(face-1)
        first_node=nodes[offset+1]
        second_node=nodes[offset+2]
        third_node=nodes[offset+3]
        (first_node!=second_node && first_node!=third_node &&
         second_node!=third_node) || throw(ArgumentError(
            "$caller: face $face contains repeated node tags"))
        key=_face_key(first_node,second_node,third_node)
        tag=get(topology.triangle_tags,key,UInt64(0))
        tag!=0 || throw(ArgumentError(
            "$caller: unknown triangular mesh face " *
            "($first_node, $second_node, $third_node)"))
        tags[face]=tag
        # Gmsh 4.15.2 leaves face orientations as zero in its public API.
        orientations[face]=Int32(0)
    end
    return tags,orientations
end

function _mesh_all_edges(topology::Union{Nothing,MeshEdgeTopology})
    topology===nothing && return UInt64[],UInt64[]
    count=length(topology.nodes)
    tags=UInt64.(1:count)
    nodes=Vector{UInt64}(undef,Base.checked_mul(2,count))
    @inbounds for edge in 1:count
        nodes[2edge-1]=UInt64(topology.nodes[edge][1])
        nodes[2edge]=UInt64(topology.nodes[edge][2])
    end
    return tags,nodes
end

function _mesh_all_faces(topology::Union{Nothing,MeshFaceTopology},
                         face_type_value,
                         caller::AbstractString="mesh_all_faces")
    face_type=_checked_face_type(face_type_value,caller)
    (topology===nothing || face_type==4) && return UInt64[],UInt64[]
    count=length(topology.triangle_nodes)
    tags=UInt64.(1:count)
    nodes=Vector{UInt64}(undef,Base.checked_mul(3,count))
    @inbounds for face in 1:count
        nodes[3face-2]=UInt64(topology.triangle_nodes[face][1])
        nodes[3face-1]=UInt64(topology.triangle_nodes[face][2])
        nodes[3face]=UInt64(topology.triangle_nodes[face][3])
    end
    return tags,nodes
end

end # module MeshEntityTopology
