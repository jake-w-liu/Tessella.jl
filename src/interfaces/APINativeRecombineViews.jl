# Allocated native steps retain the model that created them, including empty
# gaps. Keep ranges rather than storage proportional to the last step number.
struct _NativeViewStepRange
    first::Int
    last::Int
    model_slot::Any
    components::Int
end

function _native_view_step_range(ranges,step::Int)
    for range in ranges
        range.first<=step<=range.last && return range
    end
    return nothing
end

function _native_view_model_steps(view,slot)
    count=0
    for range in view.model_step_ranges
        range.model_slot===slot || continue
        count=Base.checked_add(count,Base.checked_add(range.last-range.first,1))
    end
    return count
end

# A type change starts a new sequence with the new model and default gap times.
# Same-kind updates retain the actual allocated ranges and populated tables.
function _native_view_storage(view,kind::String,slot,step::Int,components::Int)
    step>=0 || throw(ArgumentError("model view: step must be non-negative"))
    step<typemax(Int) || throw(ArgumentError("model view: step count exceeds the platform Int range"))
    reset=view.model_kind!==nothing && view.model_kind!=kind
    ranges=copy(view.model_step_ranges)
    data=reset ? empty(view.model_data) : copy(view.model_data)
    contexts=reset ? empty(view.model_contexts) : copy(view.model_contexts)
    if reset && !isempty(ranges)
        empty!(ranges)
    end
    extent=isempty(ranges) ? 0 : Base.checked_add(last(ranges).last,1)
    if step>=extent
        if !isempty(ranges) && last(ranges).model_slot===slot && last(ranges).components==components
            previous=pop!(ranges)
            push!(ranges,_NativeViewStepRange(previous.first,step,slot,components))
        else
            push!(ranges,_NativeViewStepRange(extent,step,slot,components))
        end
    end
    return (;ranges,data,contexts,reset)
end

# Native default scalar nodal data stays scalar; NodeData and ElementData
# infer their component count from the first caller row. Explicit counts
# keep the existing separate admission path.
@inline function _native_view_default_components(kind,components::Int,row_width::Integer)
    components==-1 || return components
    return kind=="ElementNodeData" ? 1 : Int(row_width)
end

# A view's sparse tag storage has a logical extent, including unused slots.
# Keep that extent as a scalar instead of allocating storage proportional to
# its greatest UID. New entries beyond the extent grow it by 100 slots.
function _native_model_view_extent(tags,initial::Integer;previous::UInt64=UInt64(0))
    initial>=0 || throw(ArgumentError("model view: negative storage extent"))
    extent=max(previous,UInt64(initial))
    for tag in tags
        tag isa Integer && !(tag isa Bool) && tag>=0 || throw(ArgumentError(
            "model view: row UIDs must be nonnegative integers"))
        value=UInt64(tag)
        value>=extent && (extent=Base.checked_add(value,UInt64(100)))
    end
    return extent
end

function _native_view_model_slot(model_name)
    _slot_capture_locked!()
    if isempty(model_name)
        index=_current_slot_index_locked()
        return index==0 ? nothing : MODEL_SLOTS[index]
    end
    index=findfirst(slot->slot.name==model_name,MODEL_SLOTS)
    return index===nothing ? nothing : MODEL_SLOTS[index]
end

function _native_view_model_count(slot,kind)
    slot===nothing && return 0
    tags=Set{UInt64}()
    if slot.mesh!==nothing
        public=slot.class===nothing ? _mesh_public_tags(slot.mesh) : slot.class.public_tags
        if public!==nothing
            if kind=="NodeData"
                union!(tags,public.node_tags)
            else
                union!(tags,public.element_tags)
            end
        end
    end
    for (_,record) in Model._discrete_mesh_records_model(slot.model)
        if kind=="NodeData"
            union!(tags,record.node_tags)
        else
            union!(tags,record.element_tags)
        end
    end
    return length(tags)
end

function _native_view_element_node_count(slot,tag::UInt64)
    slot===nothing && return nothing
    for (_,record) in Model._discrete_mesh_records_model(slot.model)
        index=findfirst(==(tag),record.element_tags)
        index===nothing || return Elements._msh_primary_nodes(msh_spec(record.element_types[index]).family)
    end
    if slot.mesh!==nothing
        public=slot.class===nothing ? _mesh_public_tags(slot.mesh) : slot.class.public_tags
        if public!==nothing
            index=findfirst(==(tag),public.element_tags)
            if index!==nothing
                offset=0
                for (msh,cells,_) in _cache_native_blocks(slot.mesh)
                    count=size(cells,2)
                    index<=offset+count && return Elements._msh_primary_nodes(msh_spec(msh).family)
                    offset=Base.checked_add(offset,count)
                end
            end
        end
    end
    return nothing
end

# A newly allocated nodal row is checked against its live element. Repeated
# writes reuse that row's storage; their retained tails are checked by the
# merge operation rather than producing another allocation warning.
function _native_view_nodal_row_warning(slot,old_tags,tags,counts,components::Int)
    components>0 || return nothing
    seen=Set{UInt64}(old_tags)
    for (tag,count) in zip(tags,counts)
        tag in seen && continue
        push!(seen,tag)
        width=_native_view_element_node_count(slot,tag)
        if width!==nothing && count<Base.checked_mul(width,components)
            return "Some elements in ElementNodeData have less values than number of nodes"
        end
    end
    return nothing
end

# Native repeated writes update complete rows and preserve other UIDs in the
# same step. A repeated UID in one write uses that write's final row.
function _native_view_merge_rows(old_tags,old_counts,old_values,tags,counts,values;
                                 one_node_width=nothing)
    selected=Dict{UInt64,Vector{Float64}}()
    for (row_tags,row_counts,row_values) in ((old_tags,old_counts,old_values),(tags,counts,values))
        length(row_tags)==length(row_counts) || throw(ArgumentError("model view: incomplete row counts"))
        offset=1
        for row in eachindex(row_tags)
            count=Int(row_counts[row]);count>0 || throw(ArgumentError("model view: empty row"))
            last=Base.checked_add(offset,count)-1
            last<=length(row_values) || throw(ArgumentError("model view: incomplete row values"))
            incoming=copy(@view row_values[offset:last]);tag=UInt64(row_tags[row])
            previous=get(selected,tag,nothing)
            if one_node_width!==nothing && count==one_node_width && previous!==nothing && length(previous)>count
                retained=copy(previous);retained[1:count]=incoming;selected[tag]=retained
            else
                selected[tag]=incoming
            end
            offset=last+1
        end
        offset-1==length(row_values) || throw(ArgumentError("model view: incomplete row values"))
    end
    result_tags=sort!(collect(keys(selected)))
    result_counts=Int32[length(selected[tag]) for tag in result_tags]
    result_values=Float64[value for tag in result_tags for value in selected[tag]]
    return (tags=result_tags,counts=result_counts,values=result_values)
end

function _native_recombine_views_plan(plan;renumber::Bool,kind=nothing)
    kind in (nothing,:node,:element) || throw(ArgumentError("model view: invalid numbering kind"))
    renumber && plan!==nothing && get(plan,:publication_required,true) &&
        get(plan,:native_recombine_finished,true) || return (updates=NamedTuple[],events=Tuple{Symbol,String}[])
    index=_current_slot_index_locked()
    index==0 && return (updates=NamedTuple[],events=Tuple{Symbol,String}[])
    slot=MODEL_SLOTS[index]
    updates=NamedTuple[];events=Tuple{Symbol,String}[]
    # Aliased views reach the same shared data through different tags;
    # renumber each data owner once.
    seen=Set{_SessionView}()
    for view_tag in sort!(collect(keys(SESSION_VIEWS)))
        view=_view_data_owner(SESSION_VIEWS[view_tag])
        view in seen && continue
        push!(seen,view)
        data_kind=view.model_kind
        data_kind in ("NodeData","ElementData","ElementNodeData") || continue
        kind===nothing || (data_kind=="NodeData")===(kind===:node) || continue
        steps=_native_view_model_steps(view,slot)
        steps>0 || continue
        description=data_kind=="NodeData" ? "nodal" : "element"
        suffix=steps==1 ? "step" : "steps"
        push!(events,(:info,"Renumbering $description model data ($steps $suffix)"))
        mapping=data_kind=="NodeData" ? get(plan,:native_view_node_tag_map,plan.node_tag_map) :
            get(plan,:native_view_element_tag_map,plan.element_tag_map)
        for key in sort!(collect(keys(view.model_contexts)))
            context=view.model_contexts[key]
            context.model_slot===slot || continue
            tags,values,components=view.model_data[key]
            result=_native_model_view_rows(tags,context.counts,values,mapping;
                extent=context.extent,renumber=true)
            updated_context=merge(context,(counts=result.counts,extent=result.extent))
            push!(updates,(view=view,key=key,data=(result.tags,result.values,components),
                context=updated_context))
            result.warning===nothing || push!(events,(:warning,result.warning))
        end
    end
    return (;updates,events)
end

function _native_recombine_views_publish!(plan)
    for update in plan.updates
        update.view.model_data[update.key]=update.data
        update.view.model_contexts[update.key]=update.context
    end
    for (kind,message) in plan.events
        _native_diagnostic_emit_locked(kind,message)
    end
    return nothing
end

# Renumbering is a detached operation over complete view rows. A source UID
# outside the logical storage extent leaves the original table untouched;
# otherwise the new table contains only rows referenced by the live map.
# This also preserves rows referring to consumed cells when numbering is off.
function _native_model_view_rows(tags,counts,values,mapping;
        extent::UInt64,renumber::Bool)
    length(tags)==length(counts) || throw(ArgumentError(
        "model view: row counts must align with row UIDs"))
    allunique(tags) || throw(ArgumentError("model view: repeated row UID"))
    all(count->count isa Integer && !(count isa Bool) && count>0,counts) ||
        throw(ArgumentError("model view: row value counts must be positive integers"))
    offsets=Vector{Int}(undef,length(counts)+1);offsets[1]=1
    for row in eachindex(counts)
        offsets[row+1]=Base.checked_add(offsets[row],Int(counts[row]))
    end
    offsets[end]-1==length(values) || throw(ArgumentError(
        "model view: values must align with complete rows"))
    original()=(tags=UInt64.(tags),counts=Int32.(counts),values=Float64.(values),
        extent=extent,warning=nothing)
    renumber || return original()
    outside=UInt64[old for old in keys(mapping) if UInt64(old)>=extent]
    if !isempty(outside)
        first_uid=minimum(outside)
        return merge(original(),(warning="Wrong source index $first_uid in step data renumbering",))
    end
    selected=Int[row for row in eachindex(tags) if haskey(mapping,tags[row])]
    sort!(selected;by=row->mapping[tags[row]],alg=MergeSort)
    result_tags=UInt64[mapping[tags[row]] for row in selected]
    allunique(result_tags) || throw(ArgumentError("model view: numbering map repeats a target UID"))
    result_counts=Int32[counts[row] for row in selected]
    result_values=Float64[]
    sizehint!(result_values,sum(Int,result_counts;init=0))
    for row in selected
        append!(result_values,@view(values[offsets[row]:offsets[row+1]-1]))
    end
    result_extent=isempty(mapping) ? UInt64(1) :
        Base.checked_add(maximum(UInt64(value) for value in Base.values(mapping)),UInt64(1))
    return (tags=result_tags,counts=result_counts,values=result_values,
        extent=result_extent,warning=nothing)
end
