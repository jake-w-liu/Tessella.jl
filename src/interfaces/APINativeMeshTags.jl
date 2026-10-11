# Mesh public labels use size_t; CAD entity tags and dense indices retain their own parsers.
@inline _mesh_query_uid(value,caller,label)=_cache_checked_public_tag(value,caller,label)
@inline _native_element_constructor_uid(uid::UInt64)=
    reinterpret(UInt64,Int64(reinterpret(Int32,UInt32(uid & UInt64(0xffffffff)))))
@inline function _native_mesh_next_uid(maximum_tag::UInt64,caller)
    return maximum_tag+UInt64(1)
end

# Input rejection is a native Msg event before any mesh constructor. Its
# process-global count controls later optimization until generation resets it.
function _native_mesh_generation_begin_locked()
    NATIVE_LOOKUP_QUERY_ERRORS[]=0
    NATIVE_MSG_LAST_ERROR[]=""
    return nothing
end

function _native_geo_abort_option_locked(value)
    OPTIONS["General.AbortOnError"]=value
    return nothing
end

function _native_geo_numeric_option_locked(family,member,value)
    name=string(family,'.',member)
    key=name=="Mesh.MinimumCirclePoints" ? "Mesh.MinCircleNodes" :
        name=="Mesh.MinimumCurvePoints" ? "Mesh.MinCurveNodes" :
        get(_OPTION_ALIASES,name,name)
    haskey(OPTIONS,key) || return nothing
    # GEO's option callback has already applied the native coercion. Publish
    # it at the assignment, including assignments preceding an aborted open.
    OPTIONS[key]=value
    if key in ("Mesh.MeshSizeFactor","Mesh.CharacteristicLengthFactor")
        OPTIONS[key=="Mesh.MeshSizeFactor" ?
            "Mesh.CharacteristicLengthFactor" : "Mesh.MeshSizeFactor"]=value
    end
    return nothing
end

function _native_mesh_input_error_locked(message::String,
                                         abort=OPTIONS["General.AbortOnError"])
    NATIVE_LOOKUP_QUERY_ERRORS[]+=1
    NATIVE_MSG_LAST_ERROR[]=message
    _native_diagnostic_emit_locked(:error,message)
    (abort==2 || abort==3) && throw(ErrorException(message))
    abort==4 && exit(1)
    return nothing
end

function _native_mesh_last_error()
    return lock(STATE_LOCK) do
        NATIVE_MSG_LAST_ERROR[]
    end
end

@inline function _native_mesh_element_bucket(kind)
    family=msh_spec(kind).family
    return family===:pnt ? 0 : family===:lin ? 1 : family===:tri ? 2 :
        family===:qua ? 3 : family===:tet ? 4 : family===:hex ? 5 :
        family===:pri ? 6 : family===:pyr ? 7 : 8
end

function _native_record_renumber_plan(records,old_tags,new_tags,node::Bool,caller;
        model=nothing)
    old_values,new_values=_tagged_renumber_lists(old_tags,new_tags,caller,node ? "node" : "element")
    positions=collect(eachindex(records))
    if !node && isempty(old_values) && model!==nothing
        physical=_native_empty_element_physical_owners(model)
        if !isempty(physical)
            owners=IdDict{DiscreteEntity,Tuple{Int,Int32}}()
            for (owner,record) in Model._discrete_mesh_records_model(model)
                owners[record]=(owner[1],Int32(owner[2]))
            end
            # The input records already follow CAD entity order; a stable
            # physical-first partition keeps each entity's stored objects.
            sort!(positions;by=i->!(get(owners,records[i],(-1,Int32(-1))) in physical),
                alg=Base.Sort.MergeSort)
        end
    end
    order=UInt64[]
    for position in positions
        record=records[position]
        if node
            append!(order,record.node_tags)
        else
            positions=sortperm(record.element_types;by=_native_mesh_element_bucket,
                alg=Base.Sort.MergeSort)
            append!(order,record.element_tags[positions])
        end
    end
    label=node ? "node" : "element"
    mapping=_tagged_renumber_mapping(order,old_values,new_values,caller,label)
    events=_tagged_renumber_mapping_events(order,old_values,new_values,label)
    return (;mapping,events)
end


_native_record_renumber_mapping(records,old_tags,new_tags,node::Bool,caller;model=nothing)=
    _native_record_renumber_plan(records,old_tags,new_tags,node,caller;model).mapping

# Validate view rows and prepare every detached replacement before sparse
# records change. The unchanged record mutation owns its parameter/connectivity
# remap; session views consume exactly the same complete native label map.
function _native_record_renumber_labels_locked!(model,records,old_tags,new_tags,node::Bool,caller)
    numbering=_native_record_renumber_plan(records,old_tags,new_tags,node,caller;model)
    mapping=numbering.mapping
    plan=(node_tag_map=node ? mapping : nothing,element_tag_map=node ? nothing : mapping)
    views=_native_recombine_views_plan(plan;renumber=true,kind=node ? :node : :element)
    _record_renumber!(records,mapping,node ? "node" : "element")
    node && Model._model_periodic_uid_remap!(model,mapping)
    maximum_tag=maximum(values(mapping);init=UInt64(0))
    counter=node ? NODE_TAG_MAX : ELEMENT_TAG_MAX
    counter[]=max(counter[],maximum_tag)
    for (kind,message) in numbering.events
        _native_diagnostic_emit_locked(kind,message)
    end
    _native_recombine_views_publish!(views)
    return nothing
end
