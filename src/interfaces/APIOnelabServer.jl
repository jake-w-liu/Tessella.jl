include("APIOnelabProcess.jl")

mutable struct _OnelabRunJob
    name::String
    command::String
    listener::Sockets.TCPServer
    socket::Union{Nothing,TCPSocket}
    process::Union{Nothing,_OnelabOwnedProcess}
    reader::Union{Nothing,Task}
    parent_socket::Union{Nothing,TCPSocket}
    connected::Bool
    stopped::Bool
    eof_seen::Bool
    reader_done::Bool
    notified::Bool
    started::Float64
end

mutable struct _OnelabRunEngine
    jobs::Vector{_OnelabRunJob}
    events::Channel{Any}
    cancelled::Threads.Atomic{Bool}
    timed_out::Threads.Atomic{Bool}
    connect_expired::Threads.Atomic{Bool}
    connect_name::Base.RefValue{String}
    pending_bytes::Threads.Atomic{Int}
    maximum_message::Int
    maximum_clients::Int
    connect_timeout::Float64
    deadline::Float64
    directory::String
    quiet::Bool
    store::Any
    command_hook::Any
    log_hook::Any
    epoch::Union{Nothing,UInt64}
end

const _ONELAB_RUN_LOCK=ReentrantLock()
const _ONELAB_ACTIVE_ENGINES=IdDict{_OnelabRunEngine,Nothing}()
const _ONELAB_SESSION_EPOCH=Threads.Atomic{UInt64}(0)
_onelab_monotonic_time()=Float64(time_ns())*1.e-9

function _onelab_run_limit(value,caller,what;infinite=false)
    value isa Real && !(value isa Bool) || throw(ArgumentError("$caller: $what must be a number"))
    result=Float64(value)
    ((isfinite(result) && result>0) || infinite && result==Inf) ||
        throw(ArgumentError("$caller: $what must be positive$(infinite ? " or Inf" : " and finite")"))
    return result
end

function _onelab_server_read_exact(socket::TCPSocket,count::Int)
    # Client and server share the incremental exact reader: a truncated frame
    # allocates according to bytes actually received, not its declared length.
    return _onelab_client_recv_exact(socket,count)
end

function _onelab_server_send(socket,type::Int32,payload::AbstractString="")
    bytes=codeunits(payload)
    length(bytes)<=typemax(Int32) || throw(ArgumentError("ONELAB server: response exceeds wire length"))
    write(socket,htol(type));write(socket,htol(Int32(length(bytes))))
    write(socket,bytes);flush(socket)
    return nothing
end

function _onelab_frame_reserve!(engine,count)
    while !engine.cancelled[]
        used=engine.pending_bytes[]
        if used<=engine.maximum_message-count &&
                Threads.atomic_cas!(engine.pending_bytes,used,used+count)==used
            return true
        end
        sleep(.001)
    end
    return false
end

function _onelab_read_job!(engine,index)
    job=engine.jobs[index]
    try
        while !engine.cancelled[]
            socket=Sockets.accept(job.listener)
            if engine.cancelled[]
                close(socket);break
            end
            job.socket=socket
            if engine.cancelled[]
                close(socket);break
            end
            handed_stop=false
            try
                while !engine.cancelled[]
                    header=_onelab_server_read_exact(socket,8)
                    kind=ltoh(reinterpret(Int32,header[1:4])[1])
                    count=Int(ltoh(reinterpret(Int32,header[5:8])[1]))
                    0<=count<=engine.maximum_message || throw(ArgumentError(
                        "ONELAB server: invalid frame length $count"))
                    _onelab_frame_reserve!(engine,count) || break
                    queued=false
                    try
                        payload=String(_onelab_server_read_exact(socket,count))
                        _onelab_valid_text(payload) || throw(ArgumentError("ONELAB server: frame is not valid UTF-8"))
                        put!(engine.events,(:frame,index,socket,kind,payload,count))
                        queued=true
                    finally
                        queued || Threads.atomic_sub!(engine.pending_bytes,count)
                    end
                    if kind==_ONELAB_MSG.STOP
                        # The caller must finish all earlier queued replies
                        # before closing this socket. A client can send STOP
                        # as soon as it receives the final reply bytes, while
                        # that write's flush has not yet returned in the pump.
                        handed_stop=true
                        return nothing
                    end
                end
            catch exception
                if exception isa EOFError || exception isa Base.IOError &&
                        exception.code in (Base.UV_ECONNRESET,Base.UV_ECONNABORTED,Base.UV_EPIPE)
                    put!(engine.events,(:eof,index,socket))
                else
                    rethrow()
                end
            finally
                handed_stop || (isopen(socket) && close(socket))
            end
            # The native client can reconnect after a failed write. Keep the
            # listener until STOP, process exit, cancellation or a deadline.
        end
    catch exception
        engine.cancelled[] || try put!(engine.events,(:error,index,exception)) catch end
    finally
        job.reader_done=true
    end
    return nothing
end

function _onelab_spawn_job!(engine,name::String,command::String,parent_socket=nothing)
    return lock(_ONELAB_RUN_LOCK) do
    engine.cancelled[] && throw(ArgumentError("API.onelab.run: client was cancelled"))
    any(job->job.name==name && !job.notified,engine.jobs) && throw(ArgumentError(
        "API.onelab.run: client $(repr(name)) is already running"))
    for other in keys(_ONELAB_ACTIVE_ENGINES)
        other===engine && continue
        other.store===engine.store && _onelab_engine_valid(other) &&
            any(job->job.name==name && !job.notified,other.jobs) && throw(ArgumentError(
                "API.onelab.run: client $(repr(name)) is already running"))
    end
    reusable=findfirst(job->job.notified && job.reader_done,engine.jobs)
    reusable!==nothing || length(engine.jobs)<engine.maximum_clients || throw(ArgumentError(
        "API.onelab.run: active client count exceeds $(engine.maximum_clients)"))
    listener=Sockets.listen(Sockets.localhost,0)
    _,port=Sockets.getsockname(listener)
    job=_OnelabRunJob(name,command,listener,nothing,nothing,nothing,
        parent_socket,false,false,false,false,false,_onelab_monotonic_time())
    index=reusable===nothing ? length(engine.jobs)+1 : reusable
    reusable===nothing ? push!(engine.jobs,job) : (engine.jobs[index]=job)
    try
        engine.cancelled[] && throw(InterruptException())
        job.process=_onelab_process_start(command,
            String["-onelab",name,"127.0.0.1:$port"];
            directory=engine.directory,quiet=engine.quiet)
        job.reader=@async _onelab_read_job!(engine,index)
    catch
        isopen(listener) && close(listener)
        job.process===nothing || _onelab_process_close!(job.process)
        reusable===nothing ? pop!(engine.jobs) : (job.notified=true;job.reader_done=true)
        rethrow()
    end
    return index
    end
end

function _onelab_engine_cancel!(engine)
    jobs=lock(_ONELAB_RUN_LOCK) do
        engine.cancelled[]=true
        copy(engine.jobs)
    end
    # Kill the owned process trees first. A pump blocked writing a large
    # response to a client that never reads is released by the peer closing;
    # closing our own socket first can wait for that blocked write's stream
    # lock indefinitely (independently reproduced on Julia 1.13.3/Windows).
    for job in jobs
        job.process===nothing || try _onelab_process_close!(job.process) catch end
    end
    for job in jobs
        socket=job.socket
        socket===nothing || try isopen(socket) && close(socket) catch end
        try isopen(job.listener) && close(job.listener) catch end
    end
    return nothing
end

# Enforces the absolute run deadline and per-client connect deadline even while
# the pump is blocked inside a response write. Cancellation is idempotent; the
# pump converts the resulting IO error into the matching ArgumentError.
function _onelab_watchdog!(engine)
    while !engine.cancelled[]
        now=_onelab_monotonic_time()
        if now>=engine.deadline
            engine.timed_out[]=true
            _onelab_engine_cancel!(engine)
            return nothing
        end
        jobs=lock(_ONELAB_RUN_LOCK) do
            copy(engine.jobs)
        end
        for job in jobs
            if !job.notified && !job.connected &&
                    now-job.started>engine.connect_timeout
                engine.connect_name[]=job.name
                engine.connect_expired[]=true
                _onelab_engine_cancel!(engine)
                return nothing
            end
        end
        wait_until=min(engine.deadline,now+0.05)
        for job in jobs
            job.notified && continue
            job.connected || (wait_until=min(wait_until,
                job.started+engine.connect_timeout))
        end
        sleep(max(0.001,wait_until-now))
    end
    return nothing
end

function _onelab_cancel_session_runs!()
    engines=lock(_ONELAB_RUN_LOCK) do
        [engine for engine in keys(_ONELAB_ACTIVE_ENGINES) if engine.epoch!==nothing]
    end
    foreach(_onelab_engine_cancel!,engines)
    return nothing
end

function _onelab_engine_valid(engine)
    !engine.cancelled[] && (engine.epoch===nothing || engine.epoch==_ONELAB_SESSION_EPOCH[])
end

function _onelab_dispatch!(engine,index,socket,kind,payload)
    job=engine.jobs[index]
    if kind==_ONELAB_MSG.START
        !job.connected || throw(ArgumentError("ONELAB server: duplicate client START"))
        pid=tryparse(Int,payload)
        pid!==nothing && pid>0 || throw(ArgumentError("ONELAB server: invalid client START"))
        job.connected=true
        job.eof_seen=false
    elseif kind==_ONELAB_MSG.STOP
        job.connected || throw(ArgumentError("ONELAB server: STOP before START"))
        any(child->child.parent_socket===socket && !child.notified,engine.jobs) &&
            throw(ArgumentError("ONELAB server: client stopped with unfinished sub-clients"))
        job.stopped=true
        isopen(socket) && close(socket)
    else
        job.connected || throw(ArgumentError("ONELAB server: message before START"))
        job.eof_seen=false
        if kind==_ONELAB_MSG.CONNECT
            fields=split(payload,'\0';keepempty=true)
            length(fields)==2 || throw(ArgumentError("ONELAB server: malformed sub-client request"))
            name=_onelab_command_text(fields[1],"API.onelab.run","client name")
            command=_onelab_command_text(fields[2],"API.onelab.run","command")
            isempty(name) && throw(ArgumentError("ONELAB server: empty sub-client name"))
            _onelab_spawn_job!(engine,name,command,socket)
        elseif kind in _ONELAB_ASYNC
            Base.invokelatest(engine.log_hook,kind,job.name,payload)
        elseif kind in (_ONELAB_MSG.MERGE_FILE,_ONELAB_MSG.PARSE_STRING,
                        _ONELAB_MSG.OPEN_PROJECT,_ONELAB_MSG.OLPARSE)
            result=Base.invokelatest(engine.command_hook,kind,payload,job.name)
            kind==_ONELAB_MSG.OLPARSE && _onelab_server_send(socket,kind,result ? "true" : "false")
        else
            Base.invokelatest(_onelab_server_parameter!,engine.store,job.name,socket,kind,payload,engine.epoch)
        end
    end
    return nothing
end

include("APIOnelabStore.jl")
include("APIOnelabPreprocess.jl")
include("APIOnelabGeo.jl")

function _onelab_engine_failure(engine,exception)
    engine.timed_out[] && return ArgumentError(
        "API.onelab.run: client run timed out")
    if engine.connect_expired[]
        name=engine.connect_name[]
        return ArgumentError("API.onelab.run: client $(repr(name)) did not connect before its deadline")
    end
    engine.cancelled[] && return ArgumentError(
        "API.onelab.run: session was reset or client was cancelled")
    return exception
end

function _onelab_pump!(engine)
    # Give an already finished, very short client its pending accept/read turn
    # before deciding that the process exited without its START frame.
    sleep(.001)
    while true
        engine.timed_out[] && throw(ArgumentError("API.onelab.run: client run timed out"))
        if engine.connect_expired[]
            name=engine.connect_name[]
            throw(ArgumentError("API.onelab.run: client $(repr(name)) did not connect before its deadline"))
        end
        _onelab_engine_valid(engine) || throw(ArgumentError("API.onelab.run: session was reset or client was cancelled"))
        _onelab_monotonic_time()<=engine.deadline || throw(ArgumentError("API.onelab.run: client run timed out"))
        while isready(engine.events)
            event=take!(engine.events)
            if event[1]===:frame
                try
                    _onelab_dispatch!(engine,event[2],event[3],event[4],event[5])
                catch exception
                    throw(_onelab_engine_failure(engine,exception))
                finally
                    Threads.atomic_sub!(engine.pending_bytes,event[6])
                end
            elseif event[1]===:error
                throw(event[3])
            elseif event[1]===:eof
                # A reconnect is still possible while the owned process lives.
                engine.jobs[event[2]].eof_seen=true
            end
        end
        complete=true
        for job in engine.jobs
            job.notified && continue
            process=job.process
            process===nothing && throw(ArgumentError("API.onelab.run: missing client process"))
            if _onelab_process_exited(process)
                # Let already queued STOP/parameter frames drain before judging
                # an exit. The reader reports EOF or STOP without waiting on us.
                if isready(engine.events)
                    complete=false;continue
                end
                process.exit_code==0 || throw(ArgumentError(
                    "API.onelab.run: client $(repr(job.name)) exited with status $(process.exit_code)"))
                job.connected || throw(ArgumentError(
                    "API.onelab.run: client $(repr(job.name)) exited without connecting"))
                if !job.stopped && !job.eof_seen && !job.reader_done
                    complete=false;continue
                end
                job.stopped || throw(ArgumentError(
                    "API.onelab.run: client $(repr(job.name)) disconnected without STOP"))
                _onelab_process_close!(process)
                isopen(job.listener) && close(job.listener)
                job.notified=true
                job.parent_socket===nothing || _onelab_server_send(
                    job.parent_socket,_ONELAB_MSG.STOP,job.name)
            else
                complete=false
                !job.connected && _onelab_monotonic_time()-job.started>engine.connect_timeout && throw(ArgumentError(
                    "API.onelab.run: client $(repr(job.name)) did not connect before its deadline"))
            end
        end
        complete && return nothing
        sleep(.001)
    end
end

function _onelab_run_engine(store,name::String,command::String;
        command_hook,log_hook,epoch=nothing,directory=pwd(),quiet=false,
        timeout=Inf,connect_timeout=60.,max_message_bytes=64*1024*1024,max_clients=32)
    caller="API.onelab.run"
    duration=_onelab_run_limit(timeout,caller,"timeout";infinite=true)
    startup=_onelab_run_limit(connect_timeout,caller,"connect_timeout")
    max_message_bytes isa Integer && !(max_message_bytes isa Bool) &&
        1<=max_message_bytes<=typemax(Int32) || throw(ArgumentError("$caller: invalid max_message_bytes"))
    max_clients isa Integer && !(max_clients isa Bool) &&
        1<=max_clients<=1024 || throw(ArgumentError("$caller: invalid max_clients"))
    engine=_OnelabRunEngine(_OnelabRunJob[],Channel{Any}(64),Threads.Atomic{Bool}(false),
        Threads.Atomic{Bool}(false),Threads.Atomic{Bool}(false),Ref(""),
        Threads.Atomic{Int}(0),Int(max_message_bytes),Int(max_clients),startup,
        isfinite(duration) ? _onelab_monotonic_time()+duration : Inf,String(directory),quiet,store,
        command_hook,log_hook,epoch)
    lock(_ONELAB_RUN_LOCK) do
        _ONELAB_ACTIVE_ENGINES[engine]=nothing
    end
    watchdog=@async _onelab_watchdog!(engine)
    try
        _onelab_engine_valid(engine) || throw(ArgumentError("$caller: session was reset"))
        _onelab_spawn_job!(engine,name,command)
        _onelab_pump!(engine)
    finally
        _onelab_engine_cancel!(engine)
        close(engine.events)
        for job in engine.jobs
            job.reader===nothing || try wait(job.reader) catch end
        end
        try wait(watchdog) catch end
        lock(_ONELAB_RUN_LOCK) do
            delete!(_ONELAB_ACTIVE_ENGINES,engine)
        end
    end
    return nothing
end
