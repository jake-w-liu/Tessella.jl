# View selection is display state. Actual list/model data determine the legal
# step count; updates do not silently change an existing displayed selection.
function _view_plain_list_steps(types,counts,data)
    steps=nothing
    for type in types
        count=counts[type];values=data[type]
        if count==0
            isempty(values) || return nothing
            continue
        end
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && return nothing
        _,nodes,components=layout
        length(values)%count==0 || return nothing
        width=length(values)÷count-3nodes
        width>0 && width%(nodes*components)==0 || return nothing
        current=width÷(nodes*components)
        steps=steps===nothing ? current : min(steps,current)
    end
    return steps
end

function _view_list_time_plan(view,type,count,values;reset_model=false)
    view=_view_data_owner(view)
    reset_model && !haskey(Post._POS_RECORD_TYPES,type) && throw(ArgumentError(
        "API.view.add_list_data: model-to-list conversion supports plain linear list records"))
    types=reset_model ? String[] : copy(view.list_types)
    type in types || push!(types,type)
    counts=reset_model ? empty(view.list_counts) : copy(view.list_counts);counts[type]=count
    data=reset_model ? empty(view.list_data) : copy(view.list_data);data[type]=values
    steps=_view_plain_list_steps(types,counts,data)
    if reset_model && steps===nothing
        count==0 && isempty(values) || throw(ArgumentError(
            "API.view.add_list_data: model-to-list conversion requires complete linear list payloads"))
        steps=0
    end
    time=reset_model ? Float64[] : view.list_time
    if steps!==nothing && !isempty(time) && steps>length(time)
        # Retain declared physical-time history; only new slots get index time.
        time=vcat(time,Float64[Float64(index) for index in length(time):steps-1])
    end
    return (;types,counts,data,steps,time)
end

function _view_data_steps(view)
    view=_view_data_owner(view)
    if !isempty(view.model_step_ranges)
        return Base.checked_add(last(view.model_step_ranges).last,1)
    end
    steps=_view_plain_list_steps(view.list_types,view.list_counts,view.list_data)
    steps===nothing && return Int(get(view.numbers,"NbTimeStep",0.0))
    view.numbers["NbTimeStep"]=Float64(steps)
    return steps
end

function _view_time_at(view,step)
    view=_view_data_owner(view)
    step<0 && return 0.0
    if !isempty(view.model_step_ranges)
        context=get(view.model_contexts,(step,view.model_kind),nothing)
        return context===nothing ? 0.0 : context.time
    end
    step<length(view.list_time) && return view.list_time[step+1]
    return Float64(step)
end

function _view_select_step!(view,requested)
    steps=_view_data_steps(view)
    step=requested<0 ? steps-1 : requested>=steps ? 0 : requested
    view.numbers["TimeStep"]=Float64(step)
    view.numbers["Time"]=steps==0 ? 0.0 : _view_time_at(view,step)
    return nothing
end

function _view_closest_step(view,time,steps)
    view=_view_data_owner(view)
    steps==0 && return 0
    candidates=Int[]
    if !isempty(view.model_step_ranges)
        for (step,kind) in keys(view.model_contexts)
            kind==view.model_kind && 0<=step<steps && push!(candidates,step)
        end
        sort!(unique!(candidates))
        # Every unpopulated allocated gap has time zero. Its earliest index
        # represents all such gaps without scanning/allocating the full extent.
        gap=0
        for step in candidates
            if step==gap
                gap+=1
            elseif step>gap
                break
            end
        end
        gap<steps && push!(candidates,gap)
        sort!(candidates)
    else
        # Actual list values already own every step; no sparse model expansion.
        candidates=0:steps-1
    end
    selected=first(candidates);distance=abs(_view_time_at(view,selected)-time)
    for step in Iterators.drop(candidates,1)
        current=abs(_view_time_at(view,step)-time)
        if current<distance;selected=step;distance=current;end
    end
    return selected
end

function _view_set_time_option!(view,name,value,caller)
    name=="NbTimeStep" && return nothing
    if name=="TimeStep"
        step=_geo_signed_gmsh_int_value(value,"$caller: TimeStep")
        return _view_select_step!(view,step)
    end
    value<0 && return nothing
    step=_view_closest_step(view,value,_view_data_steps(view))
    return _view_select_step!(view,step)
end
