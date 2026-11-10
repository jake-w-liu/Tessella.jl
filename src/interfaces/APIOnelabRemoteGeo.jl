function _onelab_remote_message_sink(session)
    return function (prefix,message)
        session.receiving_log && return nothing
        kind=startswith(prefix,"Error") ? _ONELAB_MSG.ERROR :
            startswith(prefix,"Warning") ? _ONELAB_MSG.WARNING : _ONELAB_MSG.INFO
        _onelab_remote_send(session,kind,String(message))
        return nothing
    end
end

function _onelab_remote_scope(action,name,address,epoch;timeout=60.)
    session=_onelab_remote_begin(name,address;timeout,epoch)
    storage=task_local_storage()
    previous=get(storage,_ONELAB_REMOTE_KEY,nothing)
    previous_sink=get(storage,IO._MSG_TASK_SINK_KEY,nothing)
    storage[_ONELAB_REMOTE_KEY]=session
    storage[IO._MSG_TASK_SINK_KEY]=_onelab_remote_message_sink(session)
    try
        return action(session)
    finally
        previous===nothing ? delete!(storage,_ONELAB_REMOTE_KEY) : (storage[_ONELAB_REMOTE_KEY]=previous)
        previous_sink===nothing ? delete!(storage,IO._MSG_TASK_SINK_KEY) : (storage[IO._MSG_TASK_SINK_KEY]=previous_sink)
        _onelab_remote_end!(session)
    end
end

function _onelab_remote_model_info!(session;path=nothing,source=nothing)
    _onelab_remote_valid(session)
    file=path===nothing ? session.model_file : String(path)
    session.model_file=file
    for (name,value,kind) in (("Gmsh/Model name",file,"file"),
            ("Gmsh/Model absolute path",isempty(file) ? "" : dirname(abspath(file))*"/","generic"))
        parameter=_OnelabWireParameter(:string,name)
        parameter.values_s=[value];parameter.kind_name=kind
        parameter.visible=false;parameter.read_only=true
        _onelab_remote_send(session,_ONELAB_MSG.PARAMETER,_onelab_wire_encode(parameter))
    end
    dimension=source===nothing ? model.get_dimension() : Model.model_dimension(source)
    physicals=source===nothing ? length(model.get_physical_groups()) : length(Model.model_physical_groups(source))
    for (name,value) in (("Gmsh/Model dimension",dimension),("Gmsh/Number of physical groups",physicals))
        parameter=_OnelabWireParameter(:number,name)
        parameter.values_f=Float64[value];parameter.visible=false;parameter.read_only=true
        _onelab_remote_send(session,_ONELAB_MSG.PARAMETER,_onelab_wire_encode(parameter))
    end
    session.model_info_sent=true
    return nothing
end

function _onelab_remote_default_file(session)
    parameter=_onelab_remote_fetch(session,:string,"Gmsh/Model name")
    (parameter===nothing || isempty(parameter.values_s) || isempty(first(parameter.values_s))) &&
        return nothing
    file=first(parameter.values_s)
    if !isabspath(file)
        directory=_onelab_remote_fetch(session,:string,"Gmsh/Model absolute path")
        if directory!==nothing && !isempty(directory.values_s)
            file=joinpath(first(directory.values_s),file)
        end
    end
    return file
end

function _onelab_remote_command!(session,kind,payload)
    try
        lock(STATE_LOCK) do
            _onelab_remote_valid(session)
            _onelab_live_command(kind,payload,session.name,session.epoch)
            kind==_ONELAB_MSG.OPEN_PROJECT && (session.model_file=model.get_file_name())
            _onelab_remote_model_info!(session)
        end
    catch error
        error isa InterruptException && rethrow()
        _onelab_remote_valid(session)
        _onelab_remote_send(session,_ONELAB_MSG.ERROR,sprint(showerror,error))
        session.exit_requested && rethrow()
    end
    return nothing
end

function _onelab_remote_service!(session)
    while !session.exit_requested
        frame=if isempty(session.commands)
            _onelab_remote_frame(session;idle=true)
        else
            command=popfirst!(session.commands)
            session.queued_bytes-=ncodeunits(command[2])
            command
        end
        frame===nothing && return 0
        kind,payload=frame
        if kind in _ONELAB_REMOTE_COMMANDS
            _onelab_remote_command!(session,kind,payload)
        elseif kind in _ONELAB_ASYNC
            _onelab_remote_log(session,kind,payload)
        elseif kind==_ONELAB_MSG.STOP
            return 0
        else
            _onelab_remote_send(session,_ONELAB_MSG.ERROR,"Unsupported remote message type $kind")
        end
    end
    return session.exit_code
end
