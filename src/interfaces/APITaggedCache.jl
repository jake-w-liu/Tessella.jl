# Public mesh tags are labels, while native mesh connectivity and kernels use
# dense indices. The element array follows stored block order, not MSH type
# order or the order of public element tags.
struct _CachePublicTags
    mesh::Union{Mesh,MixedMesh}
    node_tags::Vector{UInt64}
    element_tags::Vector{UInt64}
    node_indices::Dict{UInt64,Int}
    element_indices::Dict{UInt64,Int}
    authority::Set{Tuple{Int,Int32}}
end

function _cache_tag_inverse(tags,caller,kind)
    indices=Dict{UInt64,Int}()
    sizehint!(indices,length(tags))
    for (index,tag) in enumerate(tags)
        tag>0 || throw(ArgumentError("$caller: $kind tags must be positive"))
        haskey(indices,tag) && throw(ArgumentError(
            "$caller: duplicate public $kind tag $tag"))
        indices[tag]=index
    end
    return indices
end

function _cache_public_tags(mesh::Union{Mesh,MixedMesh};authority=(),
                            node_tags=nothing,element_tags=nothing,
                            node_count::Integer=nnodes(mesh),
                            caller::AbstractString="API mesh tag table")
    if node_tags!==nothing || element_tags!==nothing
        node_tags===nothing && (node_tags=collect(UInt64,1:node_count))
        if element_tags===nothing
            _,_,count=_mesh_element_offsets(mesh)
            element_tags=collect(UInt64,1:count)
        end
        node_tags=UInt64[_cache_checked_public_tag(tag,caller,"node") for tag in node_tags]
        element_tags=UInt64[_cache_checked_public_tag(tag,caller,"element") for tag in element_tags]
        length(node_tags)==node_count || throw(ArgumentError(
            "$caller: public node tags do not align with actual mesh coordinates"))
        _,_,count=_mesh_element_offsets(mesh)
        length(element_tags)==count || throw(ArgumentError(
            "$caller: public element tags do not align with mesh cells"))
    elseif mesh isa MixedMesh && mesh.entity_data!==nothing
        data=mesh.entity_data
        length(data.external_node_tags)==nnodes(mesh) || throw(ArgumentError(
            "$caller: public node tags do not align with mesh coordinates"))
        length(data.external_element_tags)==length(mesh.blocks) || throw(ArgumentError(
            "$caller: public element tags do not align with mesh blocks"))
        node_tags=copy(data.external_node_tags)
        element_tags=UInt64[]
        for (block,tags) in zip(mesh.blocks,data.external_element_tags)
            length(tags)==length(block.tags) || throw(ArgumentError(
                "$caller: public element tags do not align with block cells"))
            append!(element_tags,tags)
        end
    else
        node_tags=collect(UInt64,1:node_count)
        count=mesh isa Mesh ? size(mesh.segs,2)+size(mesh.tris,2)+size(mesh.tets,2) :
                            sum(length(block.tags) for block in mesh.blocks;init=0)
        element_tags=collect(UInt64,1:count)
    end
    covered=Set{Tuple{Int,Int32}}()
    for entity in authority
        entity isa Tuple && length(entity)==2 &&
            entity[1] isa Integer && !(entity[1] isa Bool) &&
            entity[2] isa Integer && !(entity[2] isa Bool) &&
            0<=entity[1]<=3 && 0<entity[2]<=typemax(Int32) || throw(ArgumentError(
                "$caller: authority must contain valid (dimension, entity tag) pairs"))
        push!(covered,(Int(entity[1]),Int32(entity[2])))
    end
    return _CachePublicTags(mesh,node_tags,element_tags,
        _cache_tag_inverse(node_tags,caller,"node"),
        _cache_tag_inverse(element_tags,caller,"element"),covered)
end

@inline _cache_node_tag(tags::_CachePublicTags,index::Integer)=tags.node_tags[index]
@inline _cache_node_tag(::Nothing,index::Integer)=UInt64(index)
@inline _cache_element_tag(tags::_CachePublicTags,index::Integer)=tags.element_tags[index]
@inline _cache_element_tag(::Nothing,index::Integer)=UInt64(index)
@inline _cache_node_tag(mesh::Mesh,index::Integer)=_cache_node_tag(_mesh_public_tags(mesh),index)
@inline function _cache_node_tag(mesh::MixedMesh,index::Integer)
    public=_mesh_public_tags(mesh)
    public===nothing || return public.node_tags[index]
    data=mesh.entity_data
    return data===nothing ? UInt64(index) : data.external_node_tags[index]
end
function _cache_element_tag(mesh::MixedMesh,index::Integer)
    public=_mesh_public_tags(mesh)
    public===nothing || return public.element_tags[index]
    data=mesh.entity_data
    data===nothing && return UInt64(index)
    offset=0
    for tags in data.external_element_tags
        index<=offset+length(tags) && return tags[index-offset]
        offset+=length(tags)
    end
    throw(BoundsError(data.external_element_tags,index))
end
_cache_node_tags(tags,indices)=UInt64[_cache_node_tag(tags,index) for index in indices]
_cache_element_tags(tags,indices)=UInt64[_cache_element_tag(tags,index) for index in indices]

function _cache_checked_public_tag(value,caller,kind)
    value isa Integer && !(value isa Bool) && value>0 || throw(ArgumentError(
        "$caller: $kind tag must be a positive integer"))
    return try
        UInt64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $kind tag exceeds UInt64"))
    end
end

function _cache_node_index(tags::_CachePublicTags,value,count::Integer,caller)
    tag=_cache_checked_public_tag(value,caller,"node")
    index=get(tags.node_indices,tag,0)
    index!=0 || throw(ArgumentError("$caller: unknown node tag $tag"))
    return index
end

function _cache_element_index(tags::_CachePublicTags,value,count::Integer,caller)
    tag=_cache_checked_public_tag(value,caller,"element")
    index=get(tags.element_indices,tag,0)
    index!=0 || throw(ArgumentError("$caller: unknown element tag $tag"))
    return index
end

function _cache_dense_index(value,count,caller,kind)
    tag=_cache_checked_public_tag(value,caller,kind)
    tag<=count || throw(ArgumentError("$caller: unknown $kind tag $tag"))
    return Int(tag)
end
_cache_node_index(::Nothing,value,count::Integer,caller)=
    _cache_dense_index(value,count,caller,"node")
_cache_element_index(::Nothing,value,count::Integer,caller)=
    _cache_dense_index(value,count,caller,"element")
_cache_node_index(tags::_CachePublicTags,value,caller)=
    _cache_node_index(tags,value,length(tags.node_tags),caller)
_cache_element_index(tags::_CachePublicTags,value,caller)=
    _cache_element_index(tags,value,length(tags.element_tags),caller)

@inline _cache_covers_record(::Nothing,dim::Integer,tag::Integer)=false
@inline _cache_covers_record(tags::_CachePublicTags,dim::Integer,tag::Integer)=
    (Int(dim),Int32(tag)) in tags.authority

# Only query output is sorted by numeric MSH family. Internal slots retain the
# stored block order so locator, quality, and Jacobian indices stay valid.
_cache_numeric_family_order(types)=sortperm(types;by=Int,alg=Base.Sort.MergeSort)

function _classification_with_public_tags(class,tags::_CachePublicTags)
    tags.mesh===class.mesh || throw(ArgumentError(
        "API mesh classification: public tags describe a different mesh cache"))
    return _MeshClassification(class.mesh,class.entity,class.entities,
        class.node_entities,class.boundaries,class.seg_entities,class.tri_entities,
        class.tet_entities,class.cell_entities,tags,class.edge_entities,class.face_entities,
        class.quad_entities)
end

function _classification_with_public_tags(class;authority=())
    return _classification_with_public_tags(class,
        _cache_public_tags(class.mesh;authority=authority))
end

# Explicit stored identities accompany selection and node compaction. Primary
# edges identify overlay interpolation rows even when coordinates coincide.
function _native_rebind_public_tags(class,mesh::Mesh;primary_map=nothing,selections=nothing)
    public=class.public_tags
    public===nothing && return nothing
    source=class.mesh
    node_tags=UInt64[]
    if primary_map===nothing
        nnodes(source)==nnodes(mesh) || throw(ArgumentError(
            "API mesh tag rebind: changed primary rows require an explicit map"))
        append!(node_tags,@view public.node_tags[1:nnodes(source)])
    else
        resize!(node_tags,nnodes(mesh));fill!(node_tags,0)
        for old in 1:nnodes(source)
            target=primary_map[old]
            target==0 && continue
            node_tags[target]==0 && (node_tags[target]=public.node_tags[old])
        end
        all(!iszero,node_tags) || throw(ArgumentError(
            "API mesh tag rebind: new primary rows lack a stored identity"))
    end
    overlay=_high_order_overlay(source)
    if overlay!==nothing
        old_mids=Dict{Tuple{Int32,Int32},UInt64}()
        for column in axes(_p2_cells(overlay),2),(slot,i,j) in _p2_edge_slots(overlay)
            a,b=_p2_cells(overlay)[i,column],_p2_cells(overlay)[j,column]
            primary_map===nothing || ((a,b)=(primary_map[a],primary_map[b]))
            (a==0 || b==0) && continue
            get!(old_mids,minmax(a,b),public.node_tags[_p2_cells(overlay)[slot,column]])
        end
        skeleton=ntets(mesh)>0 ? mesh.tets : mesh.tris
        slots=ntets(mesh)>0 ? HighOrder._P2_EDGE_SLOTS : HighOrder._P2TRI_EDGE_SLOTS
        seen=Set{Tuple{Int32,Int32}}()
        for column in axes(skeleton,2),(_,i,j) in slots
            key=minmax(skeleton[i,column],skeleton[j,column])
            key in seen && continue
            push!(seen,key)
            push!(node_tags,old_mids[key])
        end
    end
    element_tags=UInt64[]
    offset=0
    for (block,(_,cells,_)) in enumerate(_cache_native_blocks(source))
        positions=selections===nothing ? collect(axes(cells,2)) : selections[block]
        append!(element_tags,public.element_tags[offset .+ positions])
        offset+=size(cells,2)
    end
    return _cache_public_tags(mesh;authority=public.authority,node_tags,element_tags,
        node_count=length(node_tags))
end

# Record connectivity can refer to overlay rows as well as primary rows. Edge
# identities carry those references through primary-node compaction, including
# two formerly separate edges that now share one interpolation row.
function _native_public_node_remap(class,new_class,geometry;primary_map,source_mesh=nothing)
    public=class===nothing ? nothing : class.public_tags
    rebound=new_class===nothing ? nothing : new_class.public_tags
    remap=Dict{Int32,Int32}()
    for old in eachindex(primary_map)
        target=primary_map[old]
        target==0 && continue
        remap[Int32(_cache_node_tag(public,old))]=Int32(_cache_node_tag(rebound,target))
    end
    source=_high_order_overlay(class===nothing ? source_mesh : class.mesh)
    source===nothing && return remap
    support=Dict{Tuple{Int32,Int32},Int32}()
    for column in axes(_p2_cells(geometry),2),(slot,i,j) in _p2_edge_slots(geometry)
        cells=_p2_cells(geometry)
        support[minmax(cells[i,column],cells[j,column])]=
            Int32(_cache_node_tag(rebound,cells[slot,column]))
    end
    for column in axes(_p2_cells(source),2),(slot,i,j) in _p2_edge_slots(source)
        cells=_p2_cells(source)
        a,b=primary_map[cells[i,column]],primary_map[cells[j,column]]
        (a==0 || b==0) && continue
        remap[Int32(_cache_node_tag(public,cells[slot,column]))]=support[minmax(a,b)]
    end
    return remap
end

function _native_dedup_node_identity_check(overlay,primary_map,class,selected,scan_cache,caller)
    (overlay===nothing || !scan_cache) && return nothing
    count=length(primary_map)
    identities=fill((0,Int32(0),Int32(0)),nnodes(overlay))
    for node in 1:count
        identities[node]=(0,primary_map[node],Int32(0))
    end
    cells=_p2_cells(overlay)
    for column in axes(cells,2),(slot,i,j) in _p2_edge_slots(overlay)
        a,b=minmax(primary_map[cells[i,column]],primary_map[cells[j,column]])
        identities[cells[slot,column]]=(1,a,b)
    end
    seen=Dict{NTuple{3,Float64},Tuple{Int,Int32,Int32}}()
    for node in axes(overlay.coords,2)
        owner=node<=count ? (class===nothing ? nothing : class.node_entities[node]) :
            LAST_MESH_HIGH_ORDER_MIDS[][node-count]
        (selected===nothing || owner in selected) || continue
        key=Tuple(overlay.coords[:,node]);identity=identities[node]
        previous=get!(seen,key,identity)
        previous==identity || throw(ArgumentError(
            "$caller: coincident actual quadratic nodes on distinct primary or " *
            "support identities cannot be compacted in the native overlay"))
    end
    return nothing
end
