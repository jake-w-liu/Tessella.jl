# OCCT 7.8.1 NCollection_DataMap<Standard_Integer,TopoDS_Shape> iteration.
# Each OCC namespace owns its bucket state through bind/unbind operations.
# A removed key does not reset capacity; rebinding joins the current bucket head.
mutable struct _ModelOCCEntityMap
    heads::Vector{Int}
    links::Dict{Int,Int}
end
_ModelOCCEntityMap()=_ModelOCCEntityMap(Int[],Dict{Int,Int}())
Base.:(==)(first::_ModelOCCEntityMap,last::_ModelOCCEntityMap)=
    first.heads==last.heads && first.links==last.links

const _MODEL_OCC_MAP_PRIMES=(101,1009,2003,5003,10007,20011,37003,57037,
    65003,100019,209953,472393,995329,2359297,4478977,9437185,17915905,
    35831809,71663617,150994945,301989889,573308929,1019215873,2038431745)

function _model_occ_entity_map_bind!(map::_ModelOCCEntityMap,tag::Int)
    if isempty(map.links) || length(map.links)>length(map.heads)
        index=findfirst(prime->prime>length(map.links),_MODEL_OCC_MAP_PRIMES)
        index===nothing && throw(ArgumentError("OCC namespace map exceeds native capacity"))
        buckets=_MODEL_OCC_MAP_PRIMES[index]
        if isempty(map.heads) || buckets>length(map.heads)
            heads=zeros(Int,buckets)
            for encoded in map.heads
                while encoded!=0
                    current=encoded-1
                    following=map.links[current]
                    bucket=mod(current,buckets)+1
                    map.links[current]=heads[bucket]
                    heads[bucket]=encoded
                    encoded=following
                end
            end
            map.heads=heads
        end
    end
    haskey(map.links,tag) && return false
    bucket=mod(tag,length(map.heads))+1
    map.links[tag]=map.heads[bucket]
    map.heads[bucket]=tag+1
    return true
end

function _model_occ_entity_map_unbind!(map::_ModelOCCEntityMap,tag::Int)
    haskey(map.links,tag) || return false
    bucket=mod(tag,length(map.heads))+1
    encoded=map.heads[bucket]
    previous=0
    while encoded!=0
        current=encoded-1
        following=map.links[current]
        if current==tag
            if previous==0
                map.heads[bucket]=following
            else
                map.links[previous-1]=following
            end
            delete!(map.links,current)
            return true
        end
        previous=encoded
        encoded=following
    end
    throw(ArgumentError("OCC namespace map loses an existing shape tag"))
end

function _model_occ_entity_map_entities(map::_ModelOCCEntityMap,dimension::Int)
    result=Tuple{Int,Int}[]
    sizehint!(result,length(map.links))
    for encoded in map.heads
        while encoded!=0
            tag=encoded-1
            push!(result,(dimension,tag))
            encoded=map.links[tag]
        end
    end
    return result
end
