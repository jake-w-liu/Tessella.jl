# SPDX-FileCopyrightText: 1997-2025 C. Geuzaine, J.-F. Remacle
# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: MIT AND GPL-2.0-or-later
# Contains Gmsh 4.15.2 adaptations: geo/GModelCreateTopologyFromMesh.cpp.
# Julia adaptation/changes: Tessella contributors, 2026.
# Licenses apply to their respective portions; see LICENSES/README.md and NOTICES.md#gmsh.

# `GModel::createTopologyFromMesh` on the discrete records: the boundary of
# each region, face and curve is first looked up among the elements the model
# already has (a surface whose triangles bound a region becomes one of its
# faces), and only the rest becomes new entities, one per group of adjacent
# higher-dimensional entities. Gmsh orders some sets by pointer; pointers are
# taken to follow creation, so entity and element tags are Gmsh's, while
# boundary lists come in tag order.

# Gmsh's face and edge tables (1-based): MTetrahedron::faces_tetra,
# MHexahedron::faces_hexa, MPrism::faces_prism, MPyramid::faces_pyramid,
# MTriangle/MQuadrangle edges.
const _TOPO_FACES=Dict{Int32,Vector{Vector{Int}}}(
    Int32(4)=>[[1,3,2],[1,2,4],[1,4,3],[4,2,3]],
    Int32(5)=>[[1,4,3,2],[1,2,6,5],[1,5,8,4],[2,3,7,6],[3,4,8,7],[5,6,7,8]],
    Int32(6)=>[[1,3,2],[4,5,6],[1,2,5,4],[1,4,6,3],[2,3,6,5]],
    Int32(7)=>[[1,2,5],[4,1,5],[2,3,5],[3,4,5],[1,4,3,2]])
const _TOPO_EDGES=Dict{Int32,Vector{NTuple{2,Int}}}(
    Int32(2)=>[(1,2),(2,3),(3,1)],
    Int32(3)=>[(1,2),(2,3),(3,4),(4,1)])
# GRegion/GFace::getMeshElement order of the families.
const _TOPO_RANK=Dict{Int32,Int}(Int32(4)=>1,Int32(5)=>2,Int32(6)=>3,Int32(7)=>4,
                                 Int32(2)=>1,Int32(3)=>2,Int32(1)=>1,Int32(15)=>1)

# The element that first produced a face or an edge: record, element, side.
struct _TopoParent
    owner::Int
    element::Int
    side::Int
end

mutable struct _TopoFacePair
    first::Int          # 0 when no second region has the face
    second::Int
    parent::_TopoParent
end

mutable struct _TopoEdgeFaces
    faces::Vector{Int}  # sorted
    parent::_TopoParent
end

function _topo_supported(m::GeoModel)
    for (_,record) in m.discrete
        for msh in record.element_types
            haskey(_TOPO_RANK,msh) || return false
        end
    end
    return true
end

function _topo_exists(m::GeoModel)
    for (key,record) in m.discrete
        key[1]>=1 && isempty(record.boundary) && return false
    end
    return true
end

function _topo_tags(m::GeoModel,dim::Int)
    tags=Int[]
    for key in keys(m.discrete)
        key[1]==dim && push!(tags,key[2])
    end
    return sort!(tags)
end

function _topo_max_tag(m::GeoModel,dim::Int)
    result=0
    for key in keys(m.discrete)
        key[1]==dim && (result=max(result,key[2]))
    end
    for tag in keys(_model_entity_dictionary(m,dim))
        result=max(result,Int(tag))
    end
    return result
end

# The record's elements of the given dimension in getMeshElement order.
function _topo_elements(record::DiscreteEntity,dim::Int)
    order=Int[]
    for (index,msh) in enumerate(record.element_types)
        msh_dimension(Int(msh))==dim && push!(order,index)
    end
    sort!(order;by=index->_TOPO_RANK[record.element_types[index]],alg=Base.Sort.MergeSort)
    return order
end

function _topo_face_key(nodes::Vector{UInt64},face::Vector{Int})
    sorted=sort!(UInt64[nodes[slot] for slot in face])
    return ntuple(i->i<=length(sorted) ? sorted[i] : UInt64(0),4)
end

function _topo_new_entity!(m::GeoModel,dim::Int)
    tag=_topo_max_tag(m,dim)+1
    add_discrete_entity!(m,dim,tag)
    return tag
end

function _topo_3d!(m::GeoModel,next_element::Base.RefValue{UInt64})
    face_owner=Dict{NTuple{4,UInt64},Int}()
    for surface in _topo_tags(m,2)
        record=m.discrete[(2,surface)]
        for index in _topo_elements(record,2)
            get!(face_owner,_topo_face_key(record.element_nodes[index],
                collect(1:length(record.element_nodes[index]))),surface)
        end
    end
    pairs=Dict{NTuple{4,UInt64},_TopoFacePair}()
    region_faces=Dict{Int,Set{Int}}()
    for region in _topo_tags(m,3)
        record=m.discrete[(3,region)]
        for index in _topo_elements(record,3)
            nodes=record.element_nodes[index]
            faces=_TOPO_FACES[record.element_types[index]]
            for side in eachindex(faces)
                key=_topo_face_key(nodes,faces[side])
                surface=get(face_owner,key,0)
                if surface!=0
                    push!(get!(Set{Int},region_faces,region),surface)
                    continue
                end
                pair=get(pairs,key,nothing)
                if pair===nothing
                    pairs[key]=_TopoFacePair(0,region,_TopoParent(region,index,side))
                else
                    pair.first=region
                end
            end
        end
    end
    keys_sorted=sort!(collect(keys(pairs)))
    pair_face=Dict{NTuple{2,Int},Int}()
    for key in keys_sorted
        pair=pairs[key]
        pair.first==pair.second && continue
        normalized=pair.first==0 ? (0,pair.second) : minmax(pair.first,pair.second)
        haskey(pair_face,normalized) && continue
        surface=_topo_new_entity!(m,2)
        pair.first==0 || push!(get!(Set{Int},region_faces,pair.first),surface)
        push!(get!(Set{Int},region_faces,pair.second),surface)
        pair_face[normalized]=surface
    end
    # Gmsh looks the pair up as stored, (later region, earlier region): only
    # a boundary pair (null, region) finds its face and gets its elements.
    for key in keys_sorted
        pair=pairs[key]
        surface=get(pair_face,(pair.first,pair.second),0)
        surface==0 && continue
        parent=m.discrete[(3,pair.parent.owner)]
        msh=parent.element_types[pair.parent.element]
        face=_TOPO_FACES[msh][pair.parent.side]
        nodes=parent.element_nodes[pair.parent.element]
        tag=next_element[];next_element[]+=1
        msh==7 && continue   # a pyramid's face element is made, not kept
        _record_append_element!(m.discrete[(2,surface)],length(face)==3 ? 2 : 3,tag,
            UInt64[nodes[slot] for slot in face])
    end
    for (region,faces) in region_faces
        m.discrete[(3,region)].boundary=NTuple{2,Int}[(2,f) for f in sort!(collect(faces))]
    end
    return nothing
end

# ensureSimplyConnectedEdge: the curve keeps its first connected part (lines
# in the order the depth-first walk meets them); every further part becomes a
# new curve. Returns the curve and its new parts.
function _topo_simply_connected!(m::GeoModel,curve::Int)
    record=m.discrete[(1,curve)]
    lines=Int[]
    for (index,msh) in enumerate(record.element_types)
        msh==1 && push!(lines,index)
    end
    result=Int[curve]
    isempty(lines) && return result
    connection=Dict{UInt64,NTuple{2,Int}}()
    for line in lines, j in 1:2
        node=record.element_nodes[line][j]
        existing=get(connection,node,nothing)
        connection[node]=existing===nothing ? (line,0) : (existing[1],line)
    end
    remaining=Set{Int}(lines)
    by_tag=sort(lines;by=line->record.element_tags[line])
    parts=Vector{Vector{Int}}()
    cursor=1
    while !isempty(remaining)
        while !(by_tag[cursor] in remaining)
            cursor+=1
        end
        stack=Int[by_tag[cursor]]
        part=Int[]
        while !isempty(stack)
            line=pop!(stack)
            delete!(remaining,line)
            (isempty(part) || part[end]!=line) && push!(part,line)
            for j in 1:2
                first,second=connection[record.element_nodes[line][j]]
                if first==line && second!=0 && second in remaining
                    push!(stack,second)
                elseif second==line && first in remaining
                    push!(stack,first)
                end
            end
        end
        push!(parts,part)
    end
    others=Int[]
    for (index,msh) in enumerate(record.element_types)
        msh==1 || push!(others,index)
    end
    types=record.element_types;tags=record.element_tags;nodes=record.element_nodes
    for (number,part) in enumerate(parts)
        target=if number==1
            record
        else
            tag=_topo_new_entity!(m,1)
            push!(result,tag)
            m.discrete[(1,tag)]
        end
        keep=number==1 ? vcat(part,others) : part
        new_types=Int32[];new_tags=UInt64[];new_nodes=Vector{UInt64}[]
        for index in keep
            push!(new_types,types[index]);push!(new_tags,tags[index]);push!(new_nodes,nodes[index])
        end
        target.element_types=new_types;target.element_tags=new_tags;target.element_nodes=new_nodes
    end
    return result
end

function _topo_2d!(m::GeoModel,next_element::Base.RefValue{UInt64})
    edge_owner=Dict{NTuple{2,UInt64},Int}()
    for curve in _topo_tags(m,1)
        record=m.discrete[(1,curve)]
        for index in _topo_elements(record,1)
            nodes=record.element_nodes[index]
            edge_owner[minmax(nodes[1],nodes[2])]=curve
        end
    end
    edge_faces=Dict{NTuple{2,UInt64},_TopoEdgeFaces}()
    face_edges=Dict{Int,Set{Int}}()
    for surface in _topo_tags(m,2)
        record=m.discrete[(2,surface)]
        for index in _topo_elements(record,2)
            nodes=record.element_nodes[index]
            edges=_TOPO_EDGES[record.element_types[index]]
            for side in eachindex(edges)
                a,b=edges[side]
                key=minmax(nodes[a],nodes[b])
                curve=get(edge_owner,key,0)
                if curve!=0
                    push!(get!(Set{Int},face_edges,surface),curve)
                    continue
                end
                entry=get(edge_faces,key,nothing)
                if entry!==nothing && insorted(surface,entry.faces)
                    delete!(edge_faces,key)        # in bulk
                elseif entry===nothing
                    edge_faces[key]=_TopoEdgeFaces(Int[surface],_TopoParent(surface,index,side))
                else
                    insert!(entry.faces,searchsortedfirst(entry.faces,surface),surface)
                end
            end
        end
    end
    keys_sorted=sort!(collect(keys(edge_faces)))
    faces_edge=Dict{Vector{Int},Int}()
    for key in keys_sorted
        faces=edge_faces[key].faces
        haskey(faces_edge,faces) && continue
        curve=_topo_new_entity!(m,1)
        for surface in faces
            push!(get!(Set{Int},face_edges,surface),curve)
        end
        faces_edge[faces]=curve
    end
    for key in keys_sorted
        entry=edge_faces[key]
        curve=faces_edge[entry.faces]
        parent=m.discrete[(2,entry.parent.owner)]
        a,b=_TOPO_EDGES[parent.element_types[entry.parent.element]][entry.parent.side]
        nodes=parent.element_nodes[entry.parent.element]
        tag=next_element[];next_element[]+=1
        _record_append_element!(m.discrete[(1,curve)],1,tag,UInt64[nodes[a],nodes[b]])
    end
    splits=Dict{Int,Vector{Int}}()
    visited=Set{Int}()
    while true
        pending=Int[]
        for curve in _topo_tags(m,1)
            curve in visited || push!(pending,curve)
        end
        isempty(pending) && break
        for curve in pending
            push!(visited,curve)
            parts=_topo_simply_connected!(m,curve)
            length(parts)==1 || (splits[curve]=parts)
            for part in parts
                push!(visited,part)
            end
        end
    end
    for (_,edges) in face_edges
        added=Int[]
        for curve in edges
            parts=get(splits,curve,nothing)
            parts===nothing || append!(added,parts)
        end
        union!(edges,added)
    end
    for (surface,edges) in face_edges
        m.discrete[(2,surface)].boundary=NTuple{2,Int}[(1,c) for c in sort!(collect(edges))]
    end
    return nothing
end

function _topo_point_node(record::DiscreteEntity)
    for (index,msh) in enumerate(record.element_types)
        msh==15 && return record.element_nodes[index][1]
    end
    isempty(record.node_tags) && return nothing
    return record.node_tags[1]
end

function _topo_1d!(m::GeoModel,next_element::Base.RefValue{UInt64})
    point_of_node=Dict{UInt64,Int}()
    for point in _topo_tags(m,0)
        node=_topo_point_node(m.discrete[(0,point)])
        node===nothing || (point_of_node[node]=point)
    end
    node_edges=Dict{UInt64,Vector{Int}}()
    edge_points=Dict{Int,Vector{Int}}()
    for curve in _topo_tags(m,1)
        record=m.discrete[(1,curve)]
        for (index,msh) in enumerate(record.element_types)
            msh==1 || continue
            for j in 1:2
                node=record.element_nodes[index][j]
                point=get(point_of_node,node,0)
                if point!=0
                    points=get!(Vector{Int},edge_points,curve)
                    insorted(point,points) || insert!(points,searchsortedfirst(points,point),point)
                    continue
                end
                curves=get(node_edges,node,nothing)
                if curves!==nothing && insorted(curve,curves)
                    delete!(node_edges,node)       # in bulk
                elseif curves===nothing
                    node_edges[node]=Int[curve]
                else
                    insert!(curves,searchsortedfirst(curves,curve),curve)
                end
            end
        end
    end
    for node in sort!(collect(keys(node_edges)))
        point=_topo_new_entity!(m,0)
        record=m.discrete[(0,point)]
        tag=next_element[];next_element[]+=1
        _record_append_element!(record,15,tag,UInt64[node])
        point_of_node[node]=point
        for curve in node_edges[node]
            points=get!(Vector{Int},edge_points,curve)
            insorted(point,points) || insert!(points,searchsortedfirst(points,point),point)
        end
    end
    for curve in sort!(collect(keys(edge_points)))
        points=edge_points[curve]
        record=m.discrete[(1,curve)]
        if length(points)==2
            record.boundary=NTuple{2,Int}[(0,points[1]),(0,points[2])]
        elseif length(_topo_simply_connected!(m,curve))==1
            record.boundary=NTuple{2,Int}[(0,points[1]),(0,points[1])]
        end
        # Otherwise Gmsh reports a single- or multiply-ended curve and leaves
        # it without vertices; the next loop gives it one.
    end
    for curve in _topo_tags(m,1)
        record=m.discrete[(1,curve)]
        count(entry->entry[1]==0,record.boundary)>=1 && continue
        line=findfirst(==(Int32(1)),record.element_types)
        line===nothing && continue
        node=record.element_nodes[line][1]
        point=_topo_new_entity!(m,0)
        tag=next_element[];next_element[]+=1
        _record_append_element!(m.discrete[(0,point)],15,tag,UInt64[node])
        record.boundary=NTuple{2,Int}[(0,point),(0,point)]
    end
    return nothing
end

function _topo_node_coordinates(m::GeoModel)
    coordinates=Dict{UInt64,NTuple{3,Float64}}()
    for (_,record) in m.discrete
        for (index,tag) in enumerate(record.node_tags)
            coordinates[tag]=(record.node_coords[1,index],record.node_coords[2,index],
                              record.node_coords[3,index])
        end
    end
    return coordinates
end

# pruneMeshVertexAssociations: every node used by an element goes to the
# lowest-dimensional entity using it (the lowest tag among equals), unused
# nodes are dropped; a node that stays in its record keeps its parameters.
function _topo_prune_nodes!(m::GeoModel)
    owner=Dict{UInt64,NTuple{2,Int}}()
    for dim in 3:-1:0, tag in _topo_tags(m,dim)
        record=m.discrete[(dim,tag)]
        for nodes in record.element_nodes, node in nodes
            current=get(owner,node,nothing)
            (current===nothing || current[1]>dim) && (owner[node]=(dim,tag))
        end
    end
    source=Dict{UInt64,Tuple{NTuple{2,Int},Int}}()
    for (key,record) in m.discrete, (index,tag) in enumerate(record.node_tags)
        source[tag]=(key,index)
    end
    owned=Dict{NTuple{2,Int},Vector{UInt64}}()
    for (node,key) in owner
        haskey(source,node) && push!(get!(Vector{UInt64},owned,key),node)
    end
    rebuilt=Dict{NTuple{2,Int},DiscreteEntity}()
    for (key,record) in m.discrete
        nodes=sort!(get(owned,key,UInt64[]))
        coords=Matrix{Float64}(undef,3,length(nodes))
        width=size(record.node_params,1)
        params=width>0 ? zeros(width,length(nodes)) : Matrix{Float64}(undef,0,0)
        native=Dict{UInt64,Vector{Float64}}()
        for (column,node) in enumerate(nodes)
            from,index=source[node]
            other=m.discrete[from]
            coords[:,column]=@view other.node_coords[:,index]
            if from==key
                width>0 && size(other.node_params,2)==length(other.node_tags) &&
                    (params[:,column]=@view other.node_params[:,index])
                stored=get(other.node_native_params,node,nothing)
                stored===nothing || (native[node]=stored)
            end
        end
        copyrecord=DiscreteEntity(record.boundary,nodes,coords,params,record.element_types,
            record.element_tags,record.element_nodes,record.aux_params,native,
            record.geometry_parametrized)
        rebuilt[key]=copyrecord
    end
    for (key,record) in rebuilt
        m.discrete[key]=record
    end
    return nothing
end

function _create_topology_from_mesh!(m::GeoModel)
    _topo_exists(m) && return nothing
    dim=0
    for key in keys(m.discrete)
        dim=max(dim,key[1])
    end
    for d in 1:3
        isempty(_model_entity_dictionary(m,d)) || (dim=max(dim,d))
    end
    next_element=Ref(_model_fresh_element_tag(m))
    dim>=3 && _topo_3d!(m,next_element)
    dim>=2 && _topo_2d!(m,next_element)
    dim>=1 && _topo_1d!(m,next_element)
    _topo_prune_nodes!(m)
    return nothing
end
