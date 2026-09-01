"""
    MeshEntityTopology

Deterministic global edge and face identifiers for finalized linear-simplex
[`Mesh`](@ref) values. Automatic numbering follows first encounter in segment,
triangle, then tetrahedron order, using Gmsh's local simplex patterns. Explicit
positive identifiers can also be attached atomically to node pairs, triangles,
or quadrangles. Lookups return detached Gmsh-shaped tag and orientation arrays.
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
    identifiers::Vector{UInt64}
    tags::Dict{NTuple{2,Int32},UInt64}
    used_tags::Set{UInt64}
end

struct MeshFaceTopology
    node_count::Int
    triangle_nodes::Vector{NTuple{3,Int32}}
    triangle_identifiers::Vector{UInt64}
    triangle_tags::Dict{NTuple{3,Int32},UInt64}
    quadrangle_nodes::Vector{NTuple{4,Int32}}
    quadrangle_identifiers::Vector{UInt64}
    quadrangle_tags::Dict{NTuple{4,Int32},UInt64}
    used_tags::Set{UInt64}
end

@inline _edge_key(a::Int32,b::Int32)=a<b ? (a,b) : (b,a)

@inline function _triangle_key(a::Int32,b::Int32,c::Int32)
    a>b && ((a,b)=(b,a))
    b>c && ((b,c)=(c,b))
    a>b && ((a,b)=(b,a))
    return (a,b,c)
end

@inline function _quadrangle_key(
    a::Int32,b::Int32,c::Int32,d::Int32)
    a>b && ((a,b)=(b,a))
    c>d && ((c,d)=(d,c))
    a>c && ((a,c)=(c,a))
    b>d && ((b,d)=(d,b))
    b>c && ((b,c)=(c,b))
    return (a,b,c,d)
end

@inline _distinct_triangle(a,b,c)=a!=b && a!=c && b!=c
@inline _distinct_quadrangle(a,b,c,d)=
    a!=b && a!=c && a!=d && b!=c && b!=d && c!=d

function _next_topology_tag(used_tags::Set{UInt64},count::Int,
                            caller::AbstractString)
    candidate=UInt64(count)+UInt64(1)
    while candidate in used_tags
        candidate==typemax(UInt64) && throw(ArgumentError(
            "$caller: no positive UInt64 identifier remains"))
        candidate+=UInt64(1)
    end
    return candidate
end

@inline function _record_edge!(nodes,identifiers,tags,used_tags,
                               a::Int32,b::Int32,tag::UInt64)
    push!(nodes,(a,b))
    push!(identifiers,tag)
    tags[_edge_key(a,b)]=tag
    push!(used_tags,tag)
    return nothing
end

function _add_generated_edge!(nodes,identifiers,tags,used_tags,
                              a::Int32,b::Int32,
                              element_type::Int,cell::Int)
    a!=b || throw(ArgumentError(
        "mesh_edge_topology: type-$element_type cell $cell has a repeated " *
        "edge node $a"))
    key=_edge_key(a,b)
    haskey(tags,key) && return nothing
    tag=_next_topology_tag(
        used_tags,length(identifiers),"mesh_edge_topology")
    _record_edge!(nodes,identifiers,tags,used_tags,a,b,tag)
    return nothing
end

@inline function _record_triangle!(nodes,identifiers,tags,used_tags,
                                   a::Int32,b::Int32,c::Int32,tag::UInt64)
    push!(nodes,(a,b,c))
    push!(identifiers,tag)
    tags[_triangle_key(a,b,c)]=tag
    push!(used_tags,tag)
    return nothing
end

@inline function _record_quadrangle!(nodes,identifiers,tags,used_tags,
                                     a::Int32,b::Int32,c::Int32,d::Int32,
                                     tag::UInt64)
    push!(nodes,(a,b,c,d))
    push!(identifiers,tag)
    tags[_quadrangle_key(a,b,c,d)]=tag
    push!(used_tags,tag)
    return nothing
end

function _add_generated_face!(nodes,identifiers,tags,used_tags,
                              a::Int32,b::Int32,c::Int32,
                              element_type::Int,cell::Int)
    _distinct_triangle(a,b,c) || throw(ArgumentError(
        "mesh_face_topology: type-$element_type cell $cell has repeated " *
        "face nodes ($a, $b, $c)"))
    key=_triangle_key(a,b,c)
    haskey(tags,key) && return nothing
    tag=_next_topology_tag(
        used_tags,length(used_tags),"mesh_face_topology")
    _record_triangle!(nodes,identifiers,tags,used_tags,a,b,c,tag)
    return nothing
end

function _edge_topology_copy(
    topology::Union{Nothing,MeshEdgeTopology},node_count::Int)
    topology===nothing && return MeshEdgeTopology(
        node_count,NTuple{2,Int32}[],UInt64[],
        Dict{NTuple{2,Int32},UInt64}(),Set{UInt64}())
    topology.node_count==node_count || error(
        "mesh_edge_topology: internal topology does not match the mesh")
    return MeshEdgeTopology(
        node_count,copy(topology.nodes),copy(topology.identifiers),
        copy(topology.tags),copy(topology.used_tags))
end

function _face_topology_copy(
    topology::Union{Nothing,MeshFaceTopology},node_count::Int)
    topology===nothing && return MeshFaceTopology(
        node_count,NTuple{3,Int32}[],UInt64[],
        Dict{NTuple{3,Int32},UInt64}(),
        NTuple{4,Int32}[],UInt64[],
        Dict{NTuple{4,Int32},UInt64}(),Set{UInt64}())
    topology.node_count==node_count || error(
        "mesh_face_topology: internal topology does not match the mesh")
    return MeshFaceTopology(
        node_count,copy(topology.triangle_nodes),
        copy(topology.triangle_identifiers),copy(topology.triangle_tags),
        copy(topology.quadrangle_nodes),copy(topology.quadrangle_identifiers),
        copy(topology.quadrangle_tags),copy(topology.used_tags))
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

function _checked_topology_capacity(base::Int,additional::Int,
                                    caller::AbstractString)
    try
        return Base.checked_add(base,additional)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "$caller: topology count exceeds the platform Int range"))
    end
end

function _mesh_edge_topology(
    mesh::Mesh,topology::Union{Nothing,MeshEdgeTopology}=nothing)
    candidate_count=_checked_edge_candidate_count(mesh)
    replacement=_edge_topology_copy(topology,nnodes(mesh))
    capacity=_checked_topology_capacity(
        length(replacement.nodes),candidate_count,"mesh_edge_topology")
    sizehint!(replacement.nodes,capacity)
    sizehint!(replacement.identifiers,capacity)
    sizehint!(replacement.tags,capacity)
    sizehint!(replacement.used_tags,capacity)
    for (element_type,cells) in
        ((1,mesh.segs),(2,mesh.tris),(4,mesh.tets))
        patterns=_simplex_edge_patterns(element_type)
        @inbounds for cell in axes(cells,2),pattern in patterns
            _add_generated_edge!(
                replacement.nodes,replacement.identifiers,replacement.tags,
                replacement.used_tags,cells[pattern[1],cell],
                cells[pattern[2],cell],element_type,cell)
        end
    end
    return replacement
end

function _mesh_face_topology(
    mesh::Mesh,topology::Union{Nothing,MeshFaceTopology}=nothing)
    candidate_count=_checked_face_candidate_count(mesh)
    replacement=_face_topology_copy(topology,nnodes(mesh))
    capacity=_checked_topology_capacity(
        length(replacement.triangle_nodes),candidate_count,
        "mesh_face_topology")
    sizehint!(replacement.triangle_nodes,capacity)
    sizehint!(replacement.triangle_identifiers,capacity)
    sizehint!(replacement.triangle_tags,capacity)
    total_capacity=_checked_topology_capacity(
        length(replacement.used_tags),candidate_count,"mesh_face_topology")
    sizehint!(replacement.used_tags,total_capacity)
    for (element_type,cells) in ((2,mesh.tris),(4,mesh.tets))
        patterns=_simplex_face_patterns(element_type,3)
        @inbounds for cell in axes(cells,2),pattern in patterns
            _add_generated_face!(
                replacement.triangle_nodes,
                replacement.triangle_identifiers,
                replacement.triangle_tags,replacement.used_tags,
                cells[pattern[1],cell],cells[pattern[2],cell],
                cells[pattern[3],cell],element_type,cell)
        end
    end
    return replacement
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

function _checked_identifier_sequence(values,name::AbstractString,
                                      caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: $name must be a vector or tuple of positive integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    result=Vector{UInt64}(undef,length(values))
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: $name[$index] must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: $name[$index] must not be Bool"))
        value>0 || throw(ArgumentError(
            "$caller: $name[$index] must be positive; got $value"))
        result[index]=try
            UInt64(value)
        catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: $name[$index] exceeds the UInt64 range"))
        end
    end
    return result
end

function _checked_group_nodes(values,stride::Int,group_count::Int,
                              node_count::Int,name::AbstractString,
                              caller::AbstractString)
    nodes=_checked_node_sequence(values,stride,node_count,caller)
    expected=try
        Base.checked_mul(stride,group_count)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "$caller: $name count exceeds the platform Int range"))
    end
    length(nodes)==expected || throw(ArgumentError(
        "$caller: expected $expected node_tags for $group_count $name, " *
        "got $(length(nodes))"))
    return nodes
end

function _add_manual_edge!(nodes,identifiers,tags,used_tags,
                           a::Int32,b::Int32,tag::UInt64,index::Int,
                           caller::AbstractString)
    a!=b || throw(ArgumentError(
        "$caller: edge $index repeats node $a"))
    key=_edge_key(a,b)
    existing=get(tags,key,UInt64(0))
    if existing!=0
        existing==tag || throw(ArgumentError(
            "$caller: edge $index ($a, $b) already has tag $existing; " *
            "cannot assign tag $tag"))
        return nothing
    end
    tag in used_tags && throw(ArgumentError(
        "$caller: edge tag $tag is already assigned to another edge"))
    _record_edge!(nodes,identifiers,tags,used_tags,a,b,tag)
    return nothing
end

function _add_manual_triangle!(nodes,identifiers,tags,used_tags,
                               a::Int32,b::Int32,c::Int32,tag::UInt64,
                               index::Int,caller::AbstractString)
    _distinct_triangle(a,b,c) || throw(ArgumentError(
        "$caller: face $index contains repeated node tags"))
    key=_triangle_key(a,b,c)
    existing=get(tags,key,UInt64(0))
    if existing!=0
        existing==tag || throw(ArgumentError(
            "$caller: triangular face $index ($a, $b, $c) already has " *
            "tag $existing; cannot assign tag $tag"))
        return nothing
    end
    tag in used_tags && throw(ArgumentError(
        "$caller: face tag $tag is already assigned to another face"))
    _record_triangle!(nodes,identifiers,tags,used_tags,a,b,c,tag)
    return nothing
end

function _add_manual_quadrangle!(nodes,identifiers,tags,used_tags,
                                 a::Int32,b::Int32,c::Int32,d::Int32,
                                 tag::UInt64,index::Int,
                                 caller::AbstractString)
    _distinct_quadrangle(a,b,c,d) || throw(ArgumentError(
        "$caller: face $index contains repeated node tags"))
    key=_quadrangle_key(a,b,c,d)
    existing=get(tags,key,UInt64(0))
    if existing!=0
        existing==tag || throw(ArgumentError(
            "$caller: quadrangular face $index ($a, $b, $c, $d) already " *
            "has tag $existing; cannot assign tag $tag"))
        return nothing
    end
    tag in used_tags && throw(ArgumentError(
        "$caller: face tag $tag is already assigned to another face"))
    _record_quadrangle!(nodes,identifiers,tags,used_tags,a,b,c,d,tag)
    return nothing
end

function _mesh_add_edges(
    topology::Union{Nothing,MeshEdgeTopology},mesh::Mesh,
    edge_tag_values,edge_node_values,
    caller::AbstractString="mesh_add_edges")
    identifiers=_checked_identifier_sequence(
        edge_tag_values,"edge_tags",caller)
    nodes=_checked_group_nodes(
        edge_node_values,2,length(identifiers),nnodes(mesh),"edges",caller)
    replacement=_edge_topology_copy(topology,nnodes(mesh))
    capacity=_checked_topology_capacity(
        length(replacement.nodes),length(identifiers),caller)
    sizehint!(replacement.nodes,capacity)
    sizehint!(replacement.identifiers,capacity)
    sizehint!(replacement.tags,capacity)
    sizehint!(replacement.used_tags,capacity)
    @inbounds for index in eachindex(identifiers)
        _add_manual_edge!(
            replacement.nodes,replacement.identifiers,replacement.tags,
            replacement.used_tags,nodes[2index-1],nodes[2index],
            identifiers[index],index,caller)
    end
    return replacement
end

function _mesh_add_faces(
    topology::Union{Nothing,MeshFaceTopology},mesh::Mesh,face_type_value,
    face_tag_values,face_node_values,
    caller::AbstractString="mesh_add_faces")
    face_type=_checked_face_type(face_type_value,caller)
    identifiers=_checked_identifier_sequence(
        face_tag_values,"face_tags",caller)
    nodes=_checked_group_nodes(
        face_node_values,face_type,length(identifiers),nnodes(mesh),
        "faces",caller)
    replacement=_face_topology_copy(topology,nnodes(mesh))
    total_capacity=_checked_topology_capacity(
        length(replacement.used_tags),length(identifiers),caller)
    sizehint!(replacement.used_tags,total_capacity)
    if face_type==3
        capacity=_checked_topology_capacity(
            length(replacement.triangle_nodes),length(identifiers),caller)
        sizehint!(replacement.triangle_nodes,capacity)
        sizehint!(replacement.triangle_identifiers,capacity)
        sizehint!(replacement.triangle_tags,capacity)
        @inbounds for index in eachindex(identifiers)
            offset=3index-3
            _add_manual_triangle!(
                replacement.triangle_nodes,
                replacement.triangle_identifiers,
                replacement.triangle_tags,replacement.used_tags,
                nodes[offset+1],nodes[offset+2],nodes[offset+3],
                identifiers[index],index,caller)
        end
    else
        capacity=_checked_topology_capacity(
            length(replacement.quadrangle_nodes),length(identifiers),caller)
        sizehint!(replacement.quadrangle_nodes,capacity)
        sizehint!(replacement.quadrangle_identifiers,capacity)
        sizehint!(replacement.quadrangle_tags,capacity)
        @inbounds for index in eachindex(identifiers)
            offset=4index-4
            _add_manual_quadrangle!(
                replacement.quadrangle_nodes,
                replacement.quadrangle_identifiers,
                replacement.quadrangle_tags,replacement.used_tags,
                nodes[offset+1],nodes[offset+2],nodes[offset+3],
                nodes[offset+4],identifiers[index],index,caller)
        end
    end
    return replacement
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
        "$caller: global edges have not been populated; call add_edges or " *
        "create_edges first"))
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
        "$caller: global faces have not been populated; call add_faces or " *
        "create_faces first"))
    topology.node_count==nnodes(mesh) || error(
        "$caller: internal face topology does not match the cached mesh")
    count=length(nodes)÷face_type
    tags=Vector{UInt64}(undef,count)
    orientations=zeros(Int32,count)
    @inbounds for face in 1:count
        offset=face_type*(face-1)
        first_node=nodes[offset+1]
        second_node=nodes[offset+2]
        third_node=nodes[offset+3]
        if face_type==3
            _distinct_triangle(first_node,second_node,third_node) ||
                throw(ArgumentError(
                    "$caller: face $face contains repeated node tags"))
            key=_triangle_key(first_node,second_node,third_node)
            tag=get(topology.triangle_tags,key,UInt64(0))
            tag!=0 || throw(ArgumentError(
                "$caller: unknown triangular mesh face " *
                "($first_node, $second_node, $third_node)"))
            tags[face]=tag
        else
            fourth_node=nodes[offset+4]
            _distinct_quadrangle(
                first_node,second_node,third_node,fourth_node) ||
                throw(ArgumentError(
                    "$caller: face $face contains repeated node tags"))
            key=_quadrangle_key(
                first_node,second_node,third_node,fourth_node)
            tag=get(topology.quadrangle_tags,key,UInt64(0))
            tag!=0 || throw(ArgumentError(
                "$caller: unknown quadrangular mesh face " *
                "($first_node, $second_node, $third_node, $fourth_node)"))
            tags[face]=tag
        end
        # Gmsh 4.15.2 leaves face orientations as zero in its public API.
        orientations[face]=Int32(0)
    end
    return tags,orientations
end

function _mesh_all_entities(
    identifiers::Vector{UInt64},nodes::Vector{NTuple{N,Int32}}) where {N}
    count=length(identifiers)
    length(nodes)==count || error(
        "mesh entity topology: internal tag/node lengths differ")
    order=sortperm(identifiers)
    result_tags=Vector{UInt64}(undef,count)
    node_value_count=try
        Base.checked_mul(N,count)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "mesh entity topology: flattened node count exceeds the " *
            "platform Int range"))
    end
    result_nodes=Vector{UInt64}(undef,node_value_count)
    @inbounds for output_index in 1:count
        source_index=order[output_index]
        result_tags[output_index]=identifiers[source_index]
        source_nodes=nodes[source_index]
        offset=N*(output_index-1)
        for local_node in 1:N
            result_nodes[offset+local_node]=UInt64(source_nodes[local_node])
        end
    end
    return result_tags,result_nodes
end

function _mesh_all_edges(topology::Union{Nothing,MeshEdgeTopology})
    topology===nothing && return UInt64[],UInt64[]
    return _mesh_all_entities(topology.identifiers,topology.nodes)
end

function _mesh_all_faces(topology::Union{Nothing,MeshFaceTopology},
                         face_type_value,
                         caller::AbstractString="mesh_all_faces")
    face_type=_checked_face_type(face_type_value,caller)
    topology===nothing && return UInt64[],UInt64[]
    if face_type==3
        return _mesh_all_entities(
            topology.triangle_identifiers,topology.triangle_nodes)
    end
    return _mesh_all_entities(
        topology.quadrangle_identifiers,topology.quadrangle_nodes)
end

end # module MeshEntityTopology
