mutable struct _OnelabRemoteSession
    name::String
    address::String
    socket::Union{Nothing,Base.LibuvStream}
    owner::Task
    epoch::UInt64
    timeout::Float64
    cancelled::Threads.Atomic{Bool}
    finished::Threads.Atomic{Bool}
    timed_out::Threads.Atomic{Bool}
    deadline::Threads.Atomic{Float64}
    commands::Vector{Tuple{Int32,String}}
    queued_bytes::Int
    watchdog::Union{Nothing,Task}
    exit_requested::Bool
    exit_code::Int
    action::String
    receiving_log::Bool
    model_info_sent::Bool
    model_file::String
end

const _ONELAB_REMOTE_KEY=gensym(:tessella_onelab_remote)
const _ONELAB_ACTIVE_REMOTE=IdDict{_OnelabRemoteSession,Nothing}()
const _ONELAB_REMOTE_CAPACITY=64*1024*1024
const _ONELAB_REMOTE_COMMANDS=(_ONELAB_MSG.MERGE_FILE,_ONELAB_MSG.PARSE_STRING,_ONELAB_MSG.OPEN_PROJECT)

function _onelab_remote_current()
    session=get(task_local_storage(),_ONELAB_REMOTE_KEY,nothing)
    return session!==nothing && session.owner===current_task() ? session : nothing
end

function _onelab_remote_valid(session)
    if session.timed_out[]
        throw(ArgumentError("ONELAB remote client: socket operation timed out"))
    end
    !session.cancelled[] && !session.finished[] && session.epoch==_ONELAB_SESSION_EPOCH[] ||
        throw(ArgumentError("ONELAB remote client: session was reset or closed"))
    return nothing
end

function _onelab_remote_cancel!(session)
    session.cancelled[]=true
    socket=session.socket
    socket===nothing || _onelab_abort_socket!(socket)
    return nothing
end

function _onelab_cancel_remote_sessions!()
    sessions=lock(_ONELAB_RUN_LOCK) do
        collect(keys(_ONELAB_ACTIVE_REMOTE))
    end
    foreach(_onelab_remote_cancel!,sessions)
    return nothing
end

function _onelab_remote_watchdog(session)
    while !session.finished[]
        if session.cancelled[] || session.epoch!=_ONELAB_SESSION_EPOCH[]
            _onelab_remote_cancel!(session)
            return nothing
        elseif _onelab_monotonic_time()>=session.deadline[]
            session.timed_out[]=true
            _onelab_remote_cancel!(session)
            return nothing
        end
        sleep(.01)
    end
    return nothing
end

function _onelab_remote_io(action,session;idle=false)
    _onelab_remote_valid(session)
    previous=session.deadline[]
    session.deadline[]=idle ? previous : min(previous,_onelab_monotonic_time()+session.timeout)
    try
        result=action()
        _onelab_remote_valid(session)
        return result
    catch error
        error isa InterruptException && rethrow()
        _onelab_remote_valid(session)
        rethrow()
    finally
        session.deadline[]=previous
    end
end

function _onelab_remote_send(session,kind,payload="")
    ncodeunits(payload)<=_ONELAB_REMOTE_CAPACITY || throw(ArgumentError(
        "ONELAB remote client: message exceeds byte capacity"))
    _onelab_valid_text(payload) || throw(ArgumentError("ONELAB remote client: message is not valid UTF-8"))
    return _onelab_remote_io(session) do
        _onelab_server_send(session.socket,Int32(kind),payload)
    end
end

function _onelab_remote_frame(session;idle=false)
    first_byte=try
        _onelab_remote_io(session;idle) do
            read(session.socket,UInt8)
        end
    catch error
        error isa EOFError && return nothing
        rethrow()
    end
    return _onelab_remote_io(session) do
        header=UInt8[first_byte]
        try
            append!(header,_onelab_client_recv_exact(session.socket,7))
            kind=ltoh(reinterpret(Int32,header[1:4])[1])
            count=Int(ltoh(reinterpret(Int32,header[5:8])[1]))
            0<=count<=_ONELAB_REMOTE_CAPACITY || throw(ArgumentError(
                "ONELAB remote client: invalid message byte count"))
            payload=count==0 ? "" : String(_onelab_client_recv_exact(session.socket,count))
            _onelab_valid_text(payload) || throw(ArgumentError(
                "ONELAB remote client: message is not valid UTF-8"))
            return kind,payload
        catch error
            error isa EOFError && throw(ArgumentError("ONELAB remote client: truncated message"))
            rethrow()
        end
    end
end

function _onelab_remote_queue!(session,kind,payload)
    length(session.commands)<1024 && ncodeunits(payload)<=_ONELAB_REMOTE_CAPACITY-session.queued_bytes ||
        throw(ArgumentError("ONELAB remote client: command queue exceeds capacity"))
    push!(session.commands,(kind,payload))
    session.queued_bytes+=ncodeunits(payload)
    return nothing
end

function _onelab_remote_log(session,kind,message)
    previous=session.receiving_log;session.receiving_log=true
    try
        _onelab_live_log(kind,session.name,message,session.epoch)
    finally
        session.receiving_log=previous
    end
    return nothing
end

function _onelab_remote_reply(session;idle=false)
    while true
        frame=_onelab_remote_frame(session;idle)
        frame===nothing && throw(ArgumentError("ONELAB remote client: server closed before replying"))
        kind,payload=frame
        if kind in _ONELAB_REMOTE_COMMANDS
            _onelab_remote_queue!(session,kind,payload)
        elseif kind in _ONELAB_ASYNC
            _onelab_remote_log(session,kind,payload)
        else
            return kind,payload
        end
    end
end

function _onelab_remote_fetch(session,kind,name)
    return _onelab_remote_io(session) do
    request=_OnelabWireParameter(kind,name)
    _onelab_remote_send(session,_ONELAB_MSG.PARAMETER_QUERY,
        _onelab_wire_encode(request;max_bytes=_ONELAB_REMOTE_CAPACITY))
    message,payload=_onelab_remote_reply(session)
    message==_ONELAB_MSG.PARAMETER_NOT_FOUND && return nothing
    message==_ONELAB_MSG.PARAMETER || throw(ArgumentError("ONELAB remote client: unexpected parameter reply"))
    return _onelab_client_match_reply(_onelab_wire_decode(payload,"ONELAB remote client"),request,
        "ONELAB remote client")
    end
end

function _onelab_remote_all(session)
    return _onelab_remote_io(session) do
    result=_OnelabWireParameter[];bytes=0
    for kind in (:number,:string)
        query=_onelab_wire_encode(_OnelabWireParameter(kind,"dummy"))
        _onelab_remote_send(session,_ONELAB_MSG.PARAMETER_QUERY_ALL,query)
        while true
            message,payload=_onelab_remote_reply(session)
            message==_ONELAB_MSG.PARAMETER_QUERY_END && break
            message==_ONELAB_MSG.PARAMETER_QUERY_ALL || throw(ArgumentError(
                "ONELAB remote client: unexpected parameter list reply"))
            ncodeunits(payload)<=_ONELAB_REMOTE_CAPACITY-bytes || throw(ArgumentError(
                "ONELAB remote client: parameter list exceeds byte capacity"))
            bytes+=ncodeunits(payload)
            parameter=_onelab_wire_decode(payload,"ONELAB remote client")
            parameter.kind===kind || throw(ArgumentError("ONELAB remote client: parameter list has wrong type"))
            push!(result,parameter)
        end
    end
    return result
    end
end

function _onelab_remote_refresh!(store,session)
    parameters=_onelab_remote_all(session)
    stage=_OnelabServerStore(Dict{_OnelabParameterKey,_OnelabParameter}(),Set{String}(),
        Set{_OnelabParameterKey}(),Dict{_OnelabParameterKey,_OnelabWireParameter}(),
        Dict{String,Bool}(),copy(store.commands),ReentrantLock())
    for parameter in parameters
        _onelab_store_receive!(stage,parameter,session.name;mark_changed=false)
    end
    for metadata in values(stage.wire),(name,value) in metadata.clients
        stage.client_changed[name]=get(stage.client_changed,name,false) || value!=0
    end
    for field in (:parameters,:changed,:typed_changed,:wire,:client_changed)
        target=getfield(store,field);empty!(target)
        target isa AbstractSet ? union!(target,getfield(stage,field)) : merge!(target,getfield(stage,field))
    end
    return nothing
end

function _onelab_remote_publish!(store,session,key)
    record=get(store.parameters,key,nothing)
    record===nothing && return nothing
    wire=_onelab_store_wire(store,record,session.name)
    _onelab_remote_send(session,_ONELAB_MSG.PARAMETER,_onelab_wire_encode(wire;max_bytes=_ONELAB_REMOTE_CAPACITY))
    return nothing
end

function _onelab_remote_query!(runner,kind,name)
    session=runner.remote
    session===nothing && return nothing
    _onelab_geo_sync!(runner)
    parameter=_onelab_remote_fetch(session,kind,name)
    key=(kind,name)
    if parameter===nothing
        haskey(runner.store.parameters,key) && _onelab_store_delete!(runner.store,key)
    else
        _onelab_store_receive!(runner.store,parameter,session.name;mark_changed=false)
    end
    _onelab_geo_refresh!(runner)
    return nothing
end

function _onelab_remote_subclient!(runner,name,command)
    session=runner.remote
    _onelab_geo_sync!(runner)
    _onelab_remote_send(session,_ONELAB_MSG.CONNECT,name*"\0"*command)
    while true
        message,payload=_onelab_remote_reply(session;idle=true)
        message==_ONELAB_MSG.STOP && break
        throw(ArgumentError("ONELAB remote client: unexpected subclient completion"))
    end
    _onelab_remote_refresh!(runner.store,session)
    _onelab_geo_refresh!(runner)
    return nothing
end

function _onelab_remote_run!(runner,name,command)
    caller="ONELAB remote client"
    name=_onelab_command_text(name,caller,"client name")
    command=_onelab_command_text(command,caller,"command")
    selected=if isempty(name) && isempty(command)
        sort!(collect(runner.store.commands);by=first)
    else
        actual=isempty(name) ? "Tessella client" : name
        executable=isempty(strip(command)) ? get(runner.store.commands,actual,"") : command
        isempty(executable) || (runner.store.commands[actual]=executable)
        [actual=>executable]
    end
    for (name,command) in selected
        _onelab_remote_subclient!(runner,name,command)
    end
    return nothing
end

function _onelab_remote_attach(name,address,socket,timeout,epoch)
    session=_OnelabRemoteSession(String(name),String(address),socket,current_task(),epoch,timeout,
        Threads.Atomic{Bool}(false),Threads.Atomic{Bool}(false),Threads.Atomic{Bool}(false),
        Threads.Atomic{Float64}(Inf),Tuple{Int32,String}[],0,nothing,false,0,"",false,false,"")
    try
        lock(_ONELAB_RUN_LOCK) do
            _ONELAB_ACTIVE_REMOTE[session]=nothing
        end
        session.watchdog=@async _onelab_remote_watchdog(session)
        return session
    catch
        _onelab_remote_end!(session;notify=false)
        rethrow()
    end
end

function _onelab_remote_begin(name,address;timeout=60.,epoch=_ONELAB_SESSION_EPOCH[])
    caller="ONELAB remote client"
    timeout=_onelab_run_limit(timeout,caller,"socket timeout")
    tcp=occursin(r":\d+$",String(address))
    host,port,canonical=tcp ? _onelab_client_parse_addr(address,nothing,caller) :
        (String(address),nothing,String(address))
    socket=tcp ? TCPSocket() : Base.PipeEndpoint()
    session=_onelab_remote_attach(name,canonical,socket,timeout,epoch)
    try
        _onelab_remote_io(session) do
            port===nothing ? Sockets.connect(socket,host) : Sockets.connect(socket,host,port)
        end
        _onelab_remote_send(session,_ONELAB_MSG.START,string(getpid()))
        action=_onelab_remote_fetch(session,:string,session.name*"/Action")
        session.action=action===nothing || isempty(action.values_s) ? "" : first(action.values_s)
        return session
    catch
        _onelab_remote_end!(session;notify=false)
        rethrow()
    end
end

function _onelab_geo_status!(runner,stop,code)
    runner===nothing && return nothing
    session=runner.remote===nothing ? _onelab_listener_current() : runner.remote
    session===nothing && return nothing
    if stop===:exit
        session.exit_requested=true
        session.exit_code=Int(code)
    end
    return nothing
end

function _onelab_remote_end!(session;notify=true)
    try
        if notify && !session.cancelled[] && !session.finished[]
            try _onelab_remote_send(session,_ONELAB_MSG.STOP,"Goodbye!") catch end
        end
    finally
        session.finished[]=true
        _onelab_remote_cancel!(session)
        socket=session.socket
        socket===nothing || try close(socket) catch end
        session.socket=nothing
        watchdog=session.watchdog
        watchdog===nothing || try wait(watchdog) catch end
        lock(_ONELAB_RUN_LOCK) do
            delete!(_ONELAB_ACTIVE_REMOTE,session)
        end
        empty!(session.commands);session.queued_bytes=0
    end
    return nothing
end
