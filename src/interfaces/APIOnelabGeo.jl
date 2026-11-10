mutable struct _OnelabGeoRunner
    store::_OnelabServerStore
    numbers::Dict{String,Float64}
    strings::Dict{String,String}
    seed_numbers::Dict{String,Float64}
    seed_strings::Dict{String,String}
    epoch::Union{Nothing,UInt64}
    remote::Union{Nothing,_OnelabRemoteSession}
end

const _ONELAB_LIVE_GEO=Ref{Any}(nothing)

function _onelab_geo_namespace_stage(context,allocator)
    numbers=Dict{String,Vector{Float64}}(name=>[value] for (name,value) in context.values)
    for name in context.list_variables
        haskey(context.lists,name) && (numbers[name]=copy(context.lists[name]))
    end
    return (;context,allocator,numbers,strings=deepcopy(context.strings))
end

function _onelab_geo_adopt_namespace!(stage)
    _ONELAB_LIVE_GEO[]=(;context=stage.context,allocator=stage.allocator)
    empty!(_PARSER_NUMBERS);merge!(_PARSER_NUMBERS,stage.numbers)
    empty!(_PARSER_STRINGS);merge!(_PARSER_STRINGS,stage.strings)
    return nothing
end

function _onelab_geo_runner(;session=false)
    store=if session
        source=_ONELAB_SESSION_STORE
        # Existing records are immutable within this operation: scalar writes
        # and received wire records replace a copied record before mutation.
        _OnelabServerStore(copy(source.parameters),copy(source.changed),
            copy(source.typed_changed),deepcopy(source.wire),copy(source.client_changed),
            copy(source.commands),ReentrantLock())
    else
        _OnelabServerStore(Dict{_OnelabParameterKey,_OnelabParameter}(),Set{String}(),
            Set{_OnelabParameterKey}(),Dict{_OnelabParameterKey,_OnelabWireParameter}(),
            Dict{String,Bool}(),Dict{String,String}(),ReentrantLock())
    end
    remote=_onelab_remote_current()
    remote===nothing || _onelab_remote_refresh!(store,remote)
    runner=_OnelabGeoRunner(store,Dict{String,Float64}(),Dict{String,String}(),
        Dict{String,Float64}(),Dict{String,String}(),session ? _ONELAB_SESSION_EPOCH[] : nothing,remote)
    _onelab_geo_refresh!(runner)
    return runner
end

function _onelab_geo_refresh!(runner)
    empty!(runner.numbers);empty!(runner.strings)
    for ((_,name),parameter) in runner.store.parameters
        if parameter.kind===:number
            runner.numbers[name]=isempty(parameter.values_f) ? 0. : first(parameter.values_f)
        else
            runner.strings[name]=isempty(parameter.values_s) ? "" : first(parameter.values_s)
        end
    end
    empty!(runner.seed_numbers);merge!(runner.seed_numbers,runner.numbers)
    empty!(runner.seed_strings);merge!(runner.seed_strings,runner.strings)
    return nothing
end

function _onelab_geo_sync!(runner)
    numbers=runner.numbers;strings=runner.strings;store=runner.store
    for (kind,values,seeds) in ((:number,numbers,runner.seed_numbers),
                               (:string,strings,runner.seed_strings))
        for (name,value) in values
            key=(kind,name);old=get(store.parameters,key,nothing)
            changed=!haskey(seeds,name) || seeds[name]!=value
            if old===nothing || changed
                parameter=old===nothing ? _OnelabParameter(name,kind) : deepcopy(old)
                if kind===:number
                    parameter.values_f=Float64[value]
                else
                    parameter.values_s=String[value]
                end
                store.parameters[key]=parameter
                _onelab_store_mark_changed!(store,key)
                runner.remote===nothing || _onelab_remote_publish!(store,runner.remote,key)
            end
        end
        for name in keys(seeds)
            if !haskey(values,name)
                _onelab_store_delete!(store,(kind,name))
                runner.remote===nothing || _onelab_remote_send(runner.remote,_ONELAB_MSG.PARAMETER_CLEAR,name)
            end
        end
    end
    # Empty and multi-value records retain their complete values when their
    # scalar projection is unchanged.
    empty!(runner.seed_numbers);merge!(runner.seed_numbers,numbers)
    empty!(runner.seed_strings);merge!(runner.seed_strings,strings)
    return nothing
end

function _onelab_geo_adopt_runner!(runner)
    runner.epoch==_ONELAB_SESSION_EPOCH[] || throw(ArgumentError("ONELAB: session was reset"))
    target=_ONELAB_SESSION_STORE;source=runner.store
    for field in (:parameters,:changed,:typed_changed,:wire,:client_changed,:commands)
        destination=getfield(target,field)
        empty!(destination);merge_or_union=destination isa AbstractSet ? union! : merge!
        merge_or_union(destination,getfield(source,field))
    end
    return nothing
end

function _onelab_geo_source!(model,context,allocator,source;filename=nothing)
    statements,locations=filename===nothing ?
        GeoExec._geo_exec_statements(IOBuffer(source)) : GeoExec._geo_exec_statements(filename)
    context.exec_hook=src->GeoExec._geo_exec_value_term(model,src,context,allocator)
    context.entity_name_lookup=(dim,tag,kind)->kind===:physical ?
        get(model.physical_names,(dim,tag),"") : get(model.entity_names,(dim,tag),"")
    context.soft_unknown_reads=true
    result=GeoExec._exec_geo_statements!(model,statements,locations,firstindex(statements),
        lastindex(statements),context,allocator,Ref(0))
    _onelab_geo_status!(context.native_onelab_runner,context.stop,context.exit_code)
    context.stop===:toomany && throw(ArgumentError("ONELAB: too many GEO parse errors"))
    context.stop===:exit && context.exit_code!=0 && throw(ArgumentError("ONELAB: GEO Exit $(context.exit_code)"))
    GeoExec._geo_sync_physical_view!(model,context)
    GeoExec._geo_sync_mesh_algorithms!(model,context)
    GeoExec._geo_sync_native_surface_controls!(model,context)
    return result
end

function _onelab_geo_command!(runner,model,context,allocator,kind,payload,client)
    _onelab_geo_refresh!(runner)
    if kind==_ONELAB_MSG.PARSE_STRING
        Base.invokelatest(_onelab_geo_source!,model,context,allocator,payload)
    elseif kind==_ONELAB_MSG.MERGE_FILE || kind==_ONELAB_MSG.OPEN_PROJECT
        kind==_ONELAB_MSG.OPEN_PROJECT && GeoExec._exec_line!(model,"NewModel;",context,allocator)
        Base.invokelatest(GeoExec._geo_exec_include!,model,payload,context,
            "ONELAB client $client",allocator)
    elseif kind==_ONELAB_MSG.OLPARSE
        return Base.invokelatest(_onelab_preprocess!,runner.store,payload,runner.epoch;
            directory=dirname(abspath(context.file_name)),
            diagnostics=message->GeoExec._geo_msg_error!(context,message))
    end
    _onelab_geo_sync!(runner)
    return true
end

function _onelab_geo_context_runner!(context)
    runner=context.native_onelab_runner
    if runner===nothing
        runner=_onelab_geo_runner()
        runner.numbers=context.onelab_numbers;runner.strings=context.onelab_strings
        context.native_onelab_runner=runner
    end
    return runner
end

function _onelab_geo_run!(model,context,name,command,allocator)
    runner=_onelab_geo_context_runner!(context)
    runner.remote===nothing || return _onelab_remote_run!(runner,name,command)
    _onelab_geo_sync!(runner)
    logger=(kind,client,message)->begin
        kind==_ONELAB_MSG.ERROR ? GeoExec._geo_msg_error!(context,"$client: $message") :
        kind==_ONELAB_MSG.WARNING ? GeoExec._geo_yywarn!(context,"$client: $message") :
            GeoExec._geo_yyinfo!(context,"$client: $message")
        nothing
    end
    try
        _onelab_run_store(runner.store,name,command;epoch=runner.epoch,
            directory=dirname(abspath(context.file_name)),quiet=!context.stderr_diagnostics,
            command_hook=(kind,payload,client)->_onelab_geo_command!(runner,model,context,allocator,kind,payload,client),
            log_hook=logger)
    finally
        _onelab_geo_refresh!(runner)
        context.onelab_numbers=runner.numbers;context.onelab_strings=runner.strings
    end
    return nothing
end

function _onelab_geo_seed_namespace!(context)
    allocator=GeoExec._GeoTagAllocatorState()
    previous=_ONELAB_LIVE_GEO[]
    if previous!==nothing
        for field in (:values,:lists,:list_variables,:unavailable,:unavailable_lists,
                :stored_list_items,:functions,:structs,:struct_maxtag,:strings)
            setfield!(context,field,deepcopy(getfield(previous.context,field)))
        end
        allocator=deepcopy(previous.allocator)
    end
    for field in (:values,:lists,:unavailable,:unavailable_lists)
        table=getfield(context,field)
        for name in collect(keys(table))
            haskey(_PARSER_NUMBERS,name) || delete!(table,name)
        end
    end
    filter!(name->haskey(_PARSER_NUMBERS,name),context.list_variables)
    for (name,values) in _PARSER_NUMBERS
        unchanged=name in context.list_variables ? get(context.lists,name,nothing)==values :
            haskey(context.values,name) && length(values)==1 && context.values[name]==first(values)
        unchanged && continue
        context.lists[name]=copy(values);push!(context.list_variables,name)
        isempty(values) ? delete!(context.values,name) : (context.values[name]=first(values))
        delete!(context.unavailable,name);delete!(context.unavailable_lists,name)
    end
    context.stored_list_items=sum(length,Base.values(context.lists);init=0)
    for name in collect(keys(context.strings))
        haskey(_PARSER_STRINGS,name) || delete!(context.strings,name)
    end
    merge!(context.strings,deepcopy(_PARSER_STRINGS))
    return allocator
end

function _onelab_live_parse_string!(source;filename=nothing)
    caller=filename===nothing ? "ONELAB ParseString" : "API.parser.parse"
    current=_model_locked();index=_current_slot_index_locked()
    original=MODEL_SLOTS[index];model=Model._model_mesh_transaction_copy(current)
    complete=_model_merge_complete(current,LAST_MESH[],LAST_MESH_CLASS[],caller)
    model.meshing.mesh_authority=_model_merge_authority(model,complete,caller)
    seed=_model_merge_seed(model,original,caller)
    context=GeoExec._GeoNumericContext()
    allocator=_onelab_geo_seed_namespace!(context)
    GeoExec._geo_seed_mesh_context!(context,seed)
    GeoExec._geo_allocator_resync_model!(allocator,model)
    options=copy(OPTIONS);diagnostics=_model_import_diagnostics(options)
    for (name,value) in options
        family,member=split(name,'.';limit=2)
        context.option_numbers[(family,0,member)]=value
    end
    context.native_numeric_option_callback=(family,member,value)->_model_import_option!(options,family,member,value)
    context.native_context_lc=NATIVE_CONTEXT_LC[];context.native_context_forced_bbox=NATIVE_CONTEXT_FORCED_BBOX[]
    context.native_mesh_error_callback=diagnostics.input_error
    context.native_mesh_begin_callback=diagnostics.begin_mesh
    context.native_recombine_diagnostics=diagnostics.hooks
    fields=_model_merge_fields();context.fields=fields.fields
    context.background_field=fields.background_field;context.boundary_layer_fields=fields.boundary_layer_fields
    allocator.live_fields=Set(Int(tag) for tag in keys(context.fields))
    allocator.physical_group_max=model.physical_tag_max
    for ((dimension,entity),physicals) in model.entity_physicals,physical in physicals
        push!(get!(() -> Int[],context.raw_physicals,(dimension,physical)),entity)
    end
    context.file_name=filename===nothing ?
        (isempty(MODEL_FILE_NAME[]) ? joinpath(pwd(),"onelab.geo") : MODEL_FILE_NAME[]) : filename
    runner=_onelab_geo_runner(;session=true)
    runner.remote===nothing || (context.native_onelab_action=runner.remote.action)
    context.onelab_numbers=runner.numbers;context.onelab_strings=runner.strings
    context.native_onelab_runner=runner
    views=_geo_import_view_seed()
    GeoExec._GEO_EXEC_VIEWS[context]=copy(views.views)
    GeoExec._GEO_EXEC_VIEW_MAX[context]=views.maximum
    GeoExec._GEO_EXEC_VIEW_CALLBACK[context]=tag->(views.next[]=max(views.next[],tag);nothing)
    transfinite=Base.invokelatest(_onelab_geo_source!,model,context,allocator,source;filename)
    _onelab_geo_sync!(runner)
    namespace=_onelab_geo_namespace_stage(context,allocator)
    params=(;fields=context.fields,background_field=context.background_field,
        boundary_layer_fields=context.boundary_layer_fields)
    result=(;model,mesh=context.mesh,params,transfinite_tri=transfinite,
        mesh_node_max=context.mesh_node_max,mesh_element_max=context.mesh_element_max,
        native_context_lc=context.native_context_lc,native_context_forced_bbox=context.native_context_forced_bbox,
        native_element_lookup_receipt=context.native_element_lookup_receipt,
        views=copy(get(GeoExec._GEO_EXEC_VIEWS,context,Any[])),onelab_stage=runner)
    final=model.meshing.mesh_authority
    incoming=_model_merge_incoming_sections(complete,final===nothing ? nothing : final.mesh,caller)
    plan=_model_import_prepare_result(result,context.file_name,options,diagnostics;
        caller,detached_result=true,view_sections=incoming,geo_views=views)
    plan.slot.file_name=MODEL_FILE_NAME[];plan.slot.name=model.name
    _model_import_publish!(plan;legacy=true)
    _onelab_geo_adopt_namespace!(namespace)
    return true
end
