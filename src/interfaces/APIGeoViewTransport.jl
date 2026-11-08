# GEO executes against owned view data. Model contexts retain their actual
# source slots, while aliases retain one shared staged owner. Publication is
# prepared completely before any held session object is changed.
function _geo_import_copy_view(source,memo)
    haskey(memo,source) && return memo[source]
    target=_SessionView(source.name);memo[source]=target
    for field in fieldnames(typeof(source))
        value=getfield(source,field)
        copied=if field===:alias_of
            value===nothing ? nothing : _geo_import_copy_view(value,memo)
        elseif field===:model_step_ranges
            copy(value)
        elseif field===:model_contexts
            contexts=typeof(value)()
            for (key,context) in value
                entries=map(keys(context)) do name
                    entry=getfield(context,name)
                    name===:model_slot ? entry : deepcopy(entry)
                end
                contexts[key]=NamedTuple{keys(context)}(entries)
            end
            contexts
        else
            deepcopy(value)
        end
        setfield!(target,field,copied)
    end
    return target
end

function _geo_import_view_seed()
    memo=IdDict{Any,Any}();originals=Dict{Int,Any}()
    views=Any[]
    for tag in sort!(collect(keys(SESSION_VIEWS)))
        source=SESSION_VIEWS[tag];originals[tag]=source
        push!(views,GeoExec._GeoSessionView(tag,_geo_import_copy_view(source,memo)))
    end
    return (;views,originals,memo,maximum=VIEW_TAG_MAX[],next=Ref(VIEW_TAG_MAX[]))
end

_geo_session_view(source,caller)=source isa GeoExec._GeoSessionView ?
    source.data : _session_view_from_pos(source,caller)

function _geo_import_view_plan(sources,seed,caller)
    original_by_copy=IdDict{Any,Any}(copy=>source for (source,copy) in seed.memo)
    result=Tuple{Int,Any}[]
    for source in sources
        data=_geo_session_view(source,caller)
        push!(result,(Int(source.tag),data))
    end
    tags=first.(result)
    allunique(tags) || throw(ArgumentError("$caller: GEO view tags are not unique"))
    all(tag->0<=tag<=typemax(Int32),tags) || throw(ArgumentError("$caller: GEO view tag is out of range"))
    # Rebind staged aliases to the original held owner, including owners
    # whose public tag was deleted by the script. Their data remain alive.
    for (_,view) in result
        view.alias_of===nothing || (view.alias_of=get(original_by_copy,view.alias_of,view.alias_of))
    end
    # Unregistered owners can still be held by aliases. Publish their staged
    # data as well, even when the script removed the owner's public tag.
    updates=Tuple{Any,Any}[(original,staged) for (original,staged) in seed.memo]
    added=Tuple{Int,Any}[]
    for (tag,view) in result
        if !haskey(seed.originals,tag)
            push!(added,(tag,view))
        end
    end
    removed=Int[tag for tag in keys(seed.originals) if !(tag in tags)]
    return (;updates,added,removed,maximum=max(seed.next[],maximum(tags;init=seed.maximum)))
end

function _geo_import_publish_views!(plan)
    for tag in plan.removed;delete!(SESSION_VIEWS,tag);end
    for (original,staged) in plan.updates,field in fieldnames(typeof(original))
        setfield!(original,field,getfield(staged,field))
    end
    for (tag,view) in plan.added;SESSION_VIEWS[tag]=view;end
    return nothing
end

# View-only operations use the existing API implementations against an owned
# table, without changing the active model or any source-model identity.
function _geo_view_operation!(operation,store,maximum_tag,caller)
    saved=copy(SESSION_VIEWS);saved_max=VIEW_TAG_MAX[]
    saved_errors=NATIVE_LOOKUP_QUERY_ERRORS[];saved_message=NATIVE_MSG_LAST_ERROR[]
    try
        empty!(SESSION_VIEWS)
        for source in store;SESSION_VIEWS[Int(source.tag)]=_geo_session_view(source,caller);end
        VIEW_TAG_MAX[]=maximum_tag
        operation()
        result=Any[GeoExec._GeoSessionView(tag,SESSION_VIEWS[tag])
            for tag in sort!(collect(keys(SESSION_VIEWS)))]
        errors=NATIVE_LOOKUP_QUERY_ERRORS[]>saved_errors ? [NATIVE_MSG_LAST_ERROR[]] : String[]
        return result,VIEW_TAG_MAX[],errors
    finally
        empty!(SESSION_VIEWS);merge!(SESSION_VIEWS,saved);VIEW_TAG_MAX[]=saved_max
        NATIVE_LOOKUP_QUERY_ERRORS[]=saved_errors;NATIVE_MSG_LAST_ERROR[]=saved_message
    end
end
