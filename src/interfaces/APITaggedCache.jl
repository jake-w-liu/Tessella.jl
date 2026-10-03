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
                            caller::AbstractString="API mesh tag table")
    if mesh isa MixedMesh && mesh.entity_data!==nothing
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
        node_tags=collect(UInt64,1:nnodes(mesh))
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
@inline function _cache_node_tag(mesh::MixedMesh,index::Integer)
    data=mesh.entity_data
    return data===nothing ? UInt64(index) : data.external_node_tags[index]
end
function _cache_element_tag(mesh::MixedMesh,index::Integer)
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
        class.tet_entities,class.cell_entities,tags)
end

function _classification_with_public_tags(class;authority=())
    return _classification_with_public_tags(class,
        _cache_public_tags(class.mesh;authority=authority))
end
