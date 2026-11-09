mutable struct _OnelabServerStore
    parameters::Dict{String,_OnelabParameter}
    changed::Set{String}
    wire::Dict{String,_OnelabWireParameter}
    client_changed::Dict{String,Bool}
    commands::Dict{String,String}
    lock::Union{Nothing,ReentrantLock}
end

const _ONELAB_SESSION_STORE=_OnelabServerStore(SESSION_ONELAB,
    SESSION_ONELAB_CHANGED,Dict{String,_OnelabWireParameter}(),
    Dict{String,Bool}(),Dict{String,String}(),nothing)

function _onelab_reset_server_state!()
    Threads.atomic_add!(_ONELAB_SESSION_EPOCH,UInt64(1))
    empty!(_ONELAB_SESSION_STORE.wire)
    empty!(_ONELAB_SESSION_STORE.client_changed)
    empty!(_ONELAB_SESSION_STORE.commands)
    _ONELAB_LIVE_GEO[]=nothing
    return nothing
end

function _onelab_store_locked(action,store,epoch=nothing)
    return lock(store.lock===nothing ? STATE_LOCK : store.lock) do
        epoch===nothing || epoch==_ONELAB_SESSION_EPOCH[] || throw(ArgumentError(
            "API.onelab.run: session was reset"))
        action()
    end
end

function _onelab_store_mark_changed!(store,name)
    push!(store.changed,name)
    metadata=get(store.wire,name,nothing)
    if metadata!==nothing
        for client in keys(metadata.clients)
            metadata.clients[client]=1
            store.client_changed[client]=true
        end
    end
    return nothing
end

_onelab_mark_changed!(name)=_onelab_store_mark_changed!(_ONELAB_SESSION_STORE,name)

function _onelab_store_delete!(store,name)
    _onelab_store_mark_changed!(store,name)
    delete!(store.parameters,name);delete!(store.changed,name);delete!(store.wire,name)
    return nothing
end

function _onelab_store_register_client!(store,parameter,client)
    metadata=get(store.wire,parameter.name,nothing)
    if metadata===nothing || metadata.kind!==parameter.kind
        metadata=_OnelabWireParameter(parameter.kind,parameter.name)
        store.wire[parameter.name]=metadata
    end
    if !isempty(client) && !haskey(metadata.clients,client)
        _onelab_dict_put!(metadata.clients,metadata.client_order,client,0)
        store.client_changed[client]=get(store.client_changed,client,false) ||
            parameter.name in store.changed
    end
    return metadata
end

function _onelab_store_wire(store,parameter,client)
    result=deepcopy(_onelab_store_register_client!(store,parameter,client))
    result.name=parameter.name;result.kind=parameter.kind
    result.label=something(parameter.label,"")
    result.help=something(parameter.short_help,"")
    result.visible=parameter.visible;result.read_only=parameter.read_only
    result.attributes=copy(parameter.attributes)
    result.values_f=copy(parameter.values_f);result.values_s=copy(parameter.values_s)
    result.choices_f=copy(parameter.choices_f);result.choices_s=copy(parameter.choices_s)
    result.min=something(parameter.min,-1e200)
    result.max=something(parameter.max,1e200)
    result.step=something(parameter.step,0.)
    result.kind_name=parameter.kind_name=="parameter" ? "generic" : parameter.kind_name
    result.value_labels=Dict{Float64,String}()
    for (key,label) in parameter.value_labels
        value=tryparse(Float64,key)
        value===nothing || !isfinite(value) || (result.value_labels[value]=label)
    end
    return result
end

function _onelab_store_receive!(store,wire,client)
    isempty(wire.name) && throw(ArgumentError("ONELAB server: parameter name must be nonempty"))
    all(isfinite,wire.values_f) && all(isfinite,wire.choices_f) &&
        all(isfinite,(wire.min,wire.max,wire.step)) &&
        all(isfinite,keys(wire.value_labels)) || throw(ArgumentError(
            "ONELAB server: parameter numbers must be finite"))
    old=get(store.parameters,wire.name,nothing)
    parameter=old===nothing || old.kind!==wire.kind ?
        _OnelabParameter(wire.name,wire.kind) : deepcopy(old)
    changed=old===nothing || old.kind!==wire.kind ||
        old.values_f!=wire.values_f || old.values_s!=wire.values_s
    parameter.values_f=wire.values_f;parameter.values_s=wire.values_s
    parameter.label=isempty(wire.label) ? nothing : wire.label
    parameter.short_help=isempty(wire.help) ? nothing : wire.help
    parameter.visible=wire.visible;parameter.read_only=wire.read_only
    parameter.attributes=copy(wire.attributes)
    parameter.min=wire.min;parameter.max=wire.max;parameter.step=wire.step
    parameter.choices_f=wire.choices_f;parameter.choices_s=wire.choices_s
    parameter.value_labels=Dict(_onelab_wire_f64(value)=>label for (value,label) in wire.value_labels)
    parameter.kind_name=wire.kind_name
    previous=get(store.wire,wire.name,nothing)
    if previous!==nothing
        for (name,value) in previous.clients
            haskey(wire.clients,name) || _onelab_dict_put!(wire.clients,wire.client_order,name,value)
        end
    end
    haskey(wire.clients,client) || _onelab_dict_put!(wire.clients,wire.client_order,client,0)
    # The metadata keeps ownership/ordering without a duplicate values buffer.
    wire.values_f=Float64[];wire.values_s=String[]
    wire.choices_f=Float64[];wire.choices_s=String[]
    store.parameters[parameter.name]=parameter
    store.wire[parameter.name]=wire
    changed && _onelab_store_mark_changed!(store,parameter.name)
    return nothing
end

function _onelab_store_set_client_changed!(store,name,value)
    store.client_changed[name]=value
    for metadata in values(store.wire)
        haskey(metadata.clients,name) && (metadata.clients[name]=Int(value))
    end
    return nothing
end

function _onelab_server_parameter!(store,client,socket,kind,payload,epoch=nothing)
    return _onelab_store_locked(store,epoch) do
        if kind==_ONELAB_MSG.PARAMETER_QUERY
            query=_onelab_wire_decode(payload,"ONELAB server")
            parameter=get(store.parameters,query.name,nothing)
            if parameter===nothing || parameter.kind!==query.kind
                _onelab_server_send(socket,_ONELAB_MSG.PARAMETER_NOT_FOUND,query.name)
            else
                result=_onelab_store_wire(store,parameter,client)
                _onelab_server_send(socket,_ONELAB_MSG.PARAMETER,_onelab_wire_encode(result))
            end
        elseif kind==_ONELAB_MSG.PARAMETER || kind==_ONELAB_MSG.PARAMETER_UPDATE
            _onelab_store_receive!(store,_onelab_wire_decode(payload,"ONELAB server"),client)
        elseif kind==_ONELAB_MSG.PARAMETER_QUERY_ALL
            fields=split(payload,'\0';keepempty=true)
            wanted=if length(fields)==2
                fields[1]==_ONELAB_WIRE_VERSION && fields[2] in ("number","string") ||
                    throw(ArgumentError("ONELAB server: malformed parameter-list query"))
                Symbol(fields[2])
            else
                # Native clients may send the complete dummy parameter record;
                # the two-field version/kind request is also interoperable.
                _onelab_wire_decode(payload,"ONELAB server").kind
            end
            for name in sort!(collect(keys(store.parameters)))
                parameter=store.parameters[name]
                parameter.kind===wanted || continue
                result=_onelab_store_wire(store,parameter,client)
                _onelab_server_send(socket,_ONELAB_MSG.PARAMETER_QUERY_ALL,_onelab_wire_encode(result))
            end
            _onelab_server_send(socket,_ONELAB_MSG.PARAMETER_QUERY_END)
        elseif kind==_ONELAB_MSG.PARAMETER_CLEAR
            if isempty(payload)
                for name in collect(keys(store.parameters));_onelab_store_delete!(store,name);end
            else
                _onelab_store_delete!(store,payload)
            end
        elseif kind==_ONELAB_MSG.CLIENT_CHANGED
            fields=split(payload,'\0';keepempty=true)
            if length(fields)==2 && fields[1]=="get"
                value=get(store.client_changed,fields[2],fields[2] in store.changed)
                _onelab_server_send(socket,kind,value ? "true" : "false")
            elseif length(fields)==3 && fields[1]=="set" && fields[3] in ("true","false")
                _onelab_store_set_client_changed!(store,String(fields[2]),fields[3]=="true")
            else
                throw(ArgumentError("ONELAB server: malformed changed-flag request"))
            end
        else
            throw(ArgumentError("ONELAB server: unsupported message type $kind"))
        end
    end
end

function _onelab_run_store(store,name,command;command_hook,log_hook,epoch=nothing,
        directory=pwd(),quiet=false,kwargs...)
    caller="API.onelab.run"
    name=_onelab_command_text(name,caller,"client name")
    command=_onelab_command_text(command,caller,"command")
    selected=_onelab_store_locked(store,epoch) do
        if isempty(name) && isempty(command)
            return sort!(collect(store.commands);by=first)
        end
        actual=isempty(name) ? "Tessella client" : name
        executable=isempty(strip(command)) ? get(store.commands,actual,"") : command
        isempty(executable) && throw(ArgumentError("$caller: no command registered for client $(repr(actual))"))
        return [actual=>executable]
    end
    for (client,executable) in selected
        # Registration precedes execution so a client can be retried by name;
        # active process/socket state is always removed by the engine's finally.
        _onelab_store_locked(store,epoch) do
            store.commands[client]=executable
        end
        _onelab_run_engine(store,client,executable;command_hook,log_hook,
            epoch,directory,quiet,kwargs...)
    end
    return nothing
end

function _onelab_live_log(kind,name,message,epoch)
    lock(STATE_LOCK) do
        epoch==_ONELAB_SESSION_EPOCH[] || throw(ArgumentError("API.onelab.run: session was reset"))
        level=kind==_ONELAB_MSG.ERROR ? :error : kind==_ONELAB_MSG.WARNING ? :warning : :info
        _native_diagnostic_emit_locked(level,"$name: $message")
    end
    return nothing
end

function _onelab_live_command(kind,payload,name,epoch)
    return lock(STATE_LOCK) do
        epoch==_ONELAB_SESSION_EPOCH[] || throw(ArgumentError("API.onelab.run: session was reset"))
        if kind==_ONELAB_MSG.MERGE_FILE
            Base.invokelatest(merge,payload)
        elseif kind==_ONELAB_MSG.OPEN_PROJECT
            Base.invokelatest(open,payload)
        elseif kind==_ONELAB_MSG.PARSE_STRING
            Base.invokelatest(_onelab_live_parse_string!,payload)
        elseif kind==_ONELAB_MSG.OLPARSE
            return Base.invokelatest(_onelab_preprocess!,_ONELAB_SESSION_STORE,payload,epoch;
                diagnostics=message->_native_diagnostic_emit_locked(:error,message))
        end
        return true
    end
end

function _onelab_run_session(name,command;kwargs...)
    epoch,quiet=lock(STATE_LOCK) do
        _model_locked()
        (_ONELAB_SESSION_EPOCH[],iszero(get(OPTIONS,"General.Terminal",1.)))
    end
    return _onelab_run_store(_ONELAB_SESSION_STORE,name,command;
        epoch,quiet,command_hook=(kind,payload,client)->_onelab_live_command(kind,payload,client,epoch),
        log_hook=(kind,client,message)->_onelab_live_log(kind,client,message,epoch),kwargs...)
end
