# External listener peers own sockets. Their START process identifier does
# not grant process ownership; only CONNECT-launched solver jobs are owned.
mutable struct _OnelabListener
    listener::Any
    address::String
    path::String
    path_identity::Any
    owner::Task
    epoch::UInt64
    cancelled::Threads.Atomic{Bool}
    finished::Threads.Atomic{Bool}
    timed_out::Threads.Atomic{Bool}
    stopped::Threads.Atomic{Bool}
    option_controlled::Bool
    deadline::Float64
    input::Union{Nothing,_OnelabRemoteSession}
    output::Union{Nothing,_OnelabRemoteSession}
    engine::Union{Nothing,_OnelabRunEngine}
    watchdog::Union{Nothing,Task}
    exit_requested::Bool
    exit_code::Int
end

struct _OnelabListenerStopped <: Exception end

const _ONELAB_LISTENER_KEY=gensym(:tessella_onelab_listener)
const _ONELAB_ACTIVE_LISTENERS=IdDict{_OnelabListener,Nothing}()

mutable struct _OnelabFrontendStatus
    owner::Task
    exit_requested::Bool
    exit_code::Int
end

function _onelab_listener_current()
    listener=get(task_local_storage(),_ONELAB_LISTENER_KEY,nothing)
    return listener!==nothing && listener.owner===current_task() ? listener : nothing
end

function _onelab_listener_cancel!(listener)
    engine,input,output=lock(_ONELAB_RUN_LOCK) do
        listener.cancelled[]=true
        (listener.engine,listener.input,listener.output)
    end
    engine===nothing || _onelab_engine_cancel!(engine)
    for transport in (input,output)
        transport===nothing || _onelab_remote_cancel!(transport)
    end
    try close(listener.listener) catch end
    return nothing
end

function _onelab_cancel_listeners!()
    listeners=lock(_ONELAB_RUN_LOCK) do
        collect(keys(_ONELAB_ACTIVE_LISTENERS))
    end
    foreach(_onelab_listener_cancel!,listeners)
    return nothing
end

function _onelab_listener_valid(listener)
    listener.timed_out[] && throw(ArgumentError("ONELAB listener: service timed out"))
    listener.epoch==_ONELAB_SESSION_EPOCH[] ||
        throw(ArgumentError("ONELAB listener: session was reset or closed"))
    listener.stopped[] && throw(_OnelabListenerStopped())
    listener.cancelled[] && throw(ArgumentError("ONELAB listener: session was reset or closed"))
    return nothing
end

function _onelab_listener_disabled(listener;nonblocking=false)
    listener.option_controlled || return false
    if nonblocking
        trylock(STATE_LOCK) || return false
    else
        lock(STATE_LOCK)
    end
    try
        return iszero(get(OPTIONS,"Solver.AlwaysListen",0.))
    finally
        unlock(STATE_LOCK)
    end
end

function _onelab_listener_watchdog(listener)
    while !listener.finished[]
        if listener.epoch!=_ONELAB_SESSION_EPOCH[] || listener.cancelled[]
            _onelab_listener_cancel!(listener)
            return nothing
        elseif _onelab_monotonic_time()>=listener.deadline
            listener.timed_out[]=true
            _onelab_listener_cancel!(listener)
            return nothing
        elseif _onelab_listener_disabled(listener;nonblocking=true)
            listener.stopped[]=true
            _onelab_listener_cancel!(listener)
            return nothing
        end
        sleep(.01)
    end
    return nothing
end

_onelab_server_send(writer::_OnelabBoundWriter,kind::Int32,payload::AbstractString="")=
    _onelab_remote_send(writer.transport,kind,payload)

function _onelab_listener_reader(input,engine)
    first_frame=true
    try
        while _onelab_engine_valid(engine)
            frame=_onelab_remote_frame(input;idle=!first_frame)
            if frame===nothing
                put!(engine.events,(:listener_eof,))
                return nothing
            end
            first_frame=false
            kind,payload=frame;count=ncodeunits(payload)
            _onelab_frame_reserve!(engine,count) || return nothing
            queued=false
            try
                put!(engine.events,(:listener_frame,kind,payload,count))
                queued=true
            finally
                queued || Threads.atomic_sub!(engine.pending_bytes,count)
            end
            kind==_ONELAB_MSG.STOP && return nothing
        end
    catch error
        engine.cancelled[] || try put!(engine.events,(:listener_error,error)) catch end
    end
    return nothing
end

function _onelab_listener_model_info!(listener)
    lock(STATE_LOCK) do
        _onelab_listener_valid(listener)
        file=model.get_file_name()
        for (name,value,kind) in (("Gmsh/Model name",file,"file"),
                ("Gmsh/Model absolute path",isempty(file) ? "" : dirname(abspath(file))*"/","generic"))
            parameter=_OnelabWireParameter(:string,name)
            parameter.values_s=[value];parameter.kind_name=kind
            parameter.visible=false;parameter.read_only=true
            _onelab_store_receive!(_ONELAB_SESSION_STORE,parameter,"Gmsh";mark_changed=false)
        end
        for (name,value) in (("Gmsh/Model dimension",model.get_dimension()),
                ("Gmsh/Number of physical groups",length(model.get_physical_groups())))
            parameter=_OnelabWireParameter(:number,name)
            parameter.values_f=Float64[value];parameter.visible=false;parameter.read_only=true
            _onelab_store_receive!(_ONELAB_SESSION_STORE,parameter,"Gmsh";mark_changed=false)
        end
    end
    return nothing
end

function _onelab_listener_command!(listener,kind,payload)
    try
        _onelab_live_command(kind,payload,"Listen",listener.epoch)
        _onelab_listener_model_info!(listener)
    catch error
        error isa InterruptException && rethrow()
        _onelab_listener_valid(listener)
        _onelab_remote_send(listener.output,_ONELAB_MSG.ERROR,sprint(showerror,error))
        listener.exit_requested && rethrow()
    end
    return nothing
end

mutable struct _OnelabPeerResources
    input::Union{Nothing,_OnelabRemoteSession}
    output::Union{Nothing,_OnelabRemoteSession}
    engine::Union{Nothing,_OnelabRunEngine}
    watchdog::Union{Nothing,Task}
    reader::Union{Nothing,Task}
    replies::IdDict{TCPSocket,_OnelabBoundWriter}
end

function _onelab_listener_peer!(listener,socket,io_timeout,connect_timeout,max_clients)
    resources=_OnelabPeerResources(nothing,nothing,nothing,nothing,nothing,
        IdDict{TCPSocket,_OnelabBoundWriter}())
    storage=task_local_storage();previous=get(storage,_ONELAB_LISTENER_KEY,nothing)
    previous_sink=get(storage,IO._MSG_TASK_SINK_KEY,nothing)
    try
        directory=pwd()
        quiet=lock(STATE_LOCK) do
            _onelab_listener_valid(listener)
            iszero(get(OPTIONS,"General.Terminal",1.))
        end
        # Cleanup is armed before either ownership transfer. Publication is
        # short and locked; factory calls and watchdog waits never hold it.
        resources.input=_onelab_remote_attach("Listen",listener.address,socket,io_timeout,listener.epoch)
        lock(_ONELAB_RUN_LOCK) do
            listener.input=resources.input
        end
        _onelab_listener_valid(listener)
        resources.output=_onelab_remote_attach("Listen",listener.address,socket,io_timeout,listener.epoch)
        lock(_ONELAB_RUN_LOCK) do
            listener.output=resources.output
        end
        _onelab_listener_valid(listener)
        reply_hook=raw->get!(resources.replies,raw) do
            transport=_onelab_remote_attach("Listen solver",listener.address,raw,io_timeout,listener.epoch)
            _OnelabBoundWriter(transport)
        end
        resources.engine,resources.watchdog=_onelab_engine_start(_ONELAB_SESSION_STORE;
            command_hook=(kind,payload,name)->_onelab_listener_command!(listener,kind,payload),
            log_hook=(kind,name,message)->_onelab_live_log(kind,name,message,listener.epoch),
            epoch=listener.epoch,directory,quiet,connect_timeout,max_clients,reply_hook)
        lock(_ONELAB_RUN_LOCK) do
            listener.engine=resources.engine
        end
        _onelab_listener_valid(listener)
        input=resources.input;output=resources.output;engine=resources.engine
        writer=_OnelabBoundWriter(output);replies=resources.replies
        resources.reader=@async _onelab_listener_reader(input,engine)
        storage[_ONELAB_LISTENER_KEY]=listener
        storage[IO._MSG_TASK_SINK_KEY]=_onelab_remote_message_sink(output)
        while !listener.exit_requested
            _onelab_listener_valid(listener)
            _onelab_engine_valid(engine) || throw(_onelab_engine_failure(engine,
                ArgumentError("ONELAB listener: solver engine stopped")))
            while isready(engine.events)
                event=take!(engine.events)
                if event[1]===:listener_frame
                    _,kind,payload,count=event
                    try
                        if kind==_ONELAB_MSG.START
                            pid=tryparse(Int,payload)
                            pid!==nothing && pid>0 || throw(ArgumentError("ONELAB listener: invalid START"))
                        elseif kind==_ONELAB_MSG.STOP
                            any(job->job.parent_socket===writer && !job.notified,engine.jobs) &&
                                throw(ArgumentError("ONELAB listener: peer stopped with unfinished solvers"))
                            return nothing
                        elseif kind==_ONELAB_MSG.CONNECT
                            fields=split(payload,'\0';keepempty=true)
                            length(fields)==2 || throw(ArgumentError("ONELAB listener: malformed solver request"))
                            name=_onelab_command_text(fields[1],"ONELAB listener","client name")
                            command=_onelab_command_text(fields[2],"ONELAB listener","command")
                            isempty(name) && throw(ArgumentError("ONELAB listener: empty solver name"))
                            _onelab_spawn_job!(engine,name,command,writer)
                        elseif kind in _ONELAB_ASYNC
                            _onelab_remote_log(output,kind,payload)
                        elseif kind in _ONELAB_REMOTE_COMMANDS
                            _onelab_listener_command!(listener,kind,payload)
                        elseif kind==_ONELAB_MSG.OLPARSE
                            result=_onelab_live_command(kind,payload,"Listen",listener.epoch)
                            _onelab_server_send(writer,kind,result ? "true" : "false")
                        else
                            _onelab_server_parameter!(_ONELAB_SESSION_STORE,"Listen",writer,kind,payload,listener.epoch)
                        end
                    finally
                        Threads.atomic_sub!(engine.pending_bytes,count)
                    end
                elseif event[1]===:listener_eof
                    return nothing
                elseif event[1]===:listener_error
                    throw(event[2])
                elseif event[1]===:frame
                    try
                        _onelab_dispatch!(engine,event[2],event[3],event[4],event[5])
                    finally
                        Threads.atomic_sub!(engine.pending_bytes,event[6])
                    end
                elseif event[1]===:error
                    throw(event[3])
                elseif event[1]===:eof
                    engine.jobs[event[2]].eof_seen=true
                end
                if _onelab_listener_disabled(listener)
                    listener.stopped[]=true
                    throw(_OnelabListenerStopped())
                end
                listener.exit_requested && break
            end
            listener.exit_requested && break
            _onelab_jobs_complete!(engine)
            for socket in collect(keys(replies))
                isopen(socket) && continue
                _onelab_remote_end!(pop!(replies,socket).transport;notify=false)
            end
            sleep(.001)
        end
    finally
        previous===nothing ? delete!(storage,_ONELAB_LISTENER_KEY) : (storage[_ONELAB_LISTENER_KEY]=previous)
        previous_sink===nothing ? delete!(storage,IO._MSG_TASK_SINK_KEY) : (storage[IO._MSG_TASK_SINK_KEY]=previous_sink)
        resources.engine===nothing || _onelab_engine_finish!(resources.engine,resources.watchdog)
        for reply in values(resources.replies)
            _onelab_remote_end!(reply.transport;notify=false)
        end
        resources.output===nothing || _onelab_remote_end!(resources.output;notify=false)
        resources.input===nothing || _onelab_remote_end!(resources.input;notify=false)
        resources.reader===nothing || try wait(resources.reader) catch end
        lock(_ONELAB_RUN_LOCK) do
            listener.input=nothing;listener.output=nothing;listener.engine=nothing
        end
    end
    return nothing
end

_onelab_socket_identity(record)=(record.device,record.inode,record.mode & 0xf000)

function _onelab_listener_remove_path(path,identity)
    if !isempty(path) && identity!==nothing && (ispath(path) || islink(path))
        record=try
            lstat(path)
        catch error
            (ispath(path) || islink(path)) && rethrow()
            return nothing
        end
        _onelab_socket_identity(record)==identity && rm(path;force=true)
    end
    return nothing
end

function _onelab_unix_listener(path)
    # uv_pipe_bind records the pathname and unlinks it on close, including
    # an unrelated replacement. Bind the native descriptor ourselves, then
    # transfer only descriptor ownership through Julia's PipeServer(fd).
    # Path cleanup consequently remains under the inode/type certificate.
    bsd=Sys.KERNEL in (:Darwin,:FreeBSD,:OpenBSD,:NetBSD,:DragonFly)
    bytes=codeunits(path);capacity=bsd ? 104 : 108
    length(bytes)<capacity || throw(ArgumentError("ONELAB listener: Unix socket path is too long"))
    address=zeros(UInt8,2+length(bytes)+1)
    if bsd
        address[1]=UInt8(length(address));address[2]=0x01
    else
        copyto!(address,1,reinterpret(UInt8,Cushort[1]),1,2)
    end
    copyto!(address,3,bytes,1,length(bytes))
    descriptor_ref=Ref{Cint}(-1)
    server=Ref{Union{Nothing,Sockets.PipeServer}}(nothing)
    adopted=Ref(false);identity=Ref{Any}(nothing)
    try
        # Serialize creation/FD_CLOEXEC with Julia's process spawn IO lock;
        # Linux additionally sets SOCK_CLOEXEC in the socket syscall itself.
        Base.disable_sigint() do
            Base.iolock_begin()
            try
                kind=Sys.islinux() ? Cint(1|0x80000) : Cint(1)
                descriptor_ref[]=ccall(:socket,Cint,(Cint,Cint,Cint),1,kind,0)
                Base.systemerror("ONELAB listener: socket",descriptor_ref[]<0)
                Base.systemerror("ONELAB listener: close-on-exec",
                    ccall(:fcntl,Cint,(Cint,Cint,Cint),descriptor_ref[],2,1)<0)
            finally
                Base.iolock_end()
            end
        end
        descriptor=descriptor_ref[]
        Base.disable_sigint() do
            bind_result=GC.@preserve address ccall(:bind,Cint,(Cint,Ptr{UInt8},Cuint),
                descriptor,address,length(address))
            Base.systemerror("ONELAB listener: bind",bind_result<0)
            path_record=lstat(path)
            path_record.mode & 0xf000 == 0xc000 || throw(ArgumentError(
                "ONELAB listener: socket pathname changed while binding"))
            identity[]=_onelab_socket_identity(path_record)
        end
        server[]=Sockets.PipeServer()
        Base.disable_sigint() do
            # PipeServer(fd) was added after Julia 1.12. Use its public
            # libuv adoption operation directly, without a Base extension.
            Base.iolock_begin()
            try
                open_result=GC.@preserve server ccall(:uv_pipe_open,Cint,
                    (Ptr{Cvoid},Cint),server[].handle,descriptor)
                Base.uv_error("ONELAB listener: pipe_open",open_result)
                adopted[]=true
                server[].status=Base.StatusOpen
            finally
                Base.iolock_end()
            end
        end
        Sockets.listen(server[])
        return server[],path,path,identity[]
    catch
        server[]===nothing || close(server[])
        if !adopted[] && descriptor_ref[]>=0
            ccall(:close,Cint,(Cint,),descriptor_ref[])
        end
        _onelab_listener_remove_path(path,identity[])
        rethrow()
    end
end

function _onelab_listener_bind(address)
    caller="ONELAB listener"
    address=_onelab_command_text(address,caller,"address")
    isempty(address) && throw(ArgumentError("$caller: empty address"))
    if occursin(r":\d+$",address)
        host,port,_=_onelab_client_parse_addr(address,nothing,caller)
        listener=Sockets.listen(Sockets.getaddrinfo(host),port)
        try
            ip,actual=Sockets.getsockname(listener)
            return listener,string(ip)*":"*string(actual),"",nothing
        catch
            close(listener);rethrow()
        end
    end
    path=Sys.iswindows() ? address : abspath(address)
    !Sys.iswindows() && (ispath(path) || islink(path)) &&
        throw(ArgumentError("$caller: socket path already exists"))
    Sys.iswindows() || return _onelab_unix_listener(path)
    listener=Sockets.listen(path)
    try
        return listener,path,path,nothing
    catch
        close(listener);rethrow()
    end
end

function _onelab_listen(address;timeout=Inf,io_timeout=60.,connect_timeout=60.,
        max_clients=32,ready=nothing,option_controlled=false)
    caller="ONELAB listener"
    duration=_onelab_run_limit(timeout,caller,"timeout";infinite=true)
    io_limit=_onelab_run_limit(io_timeout,caller,"IO timeout")
    startup=_onelab_run_limit(connect_timeout,caller,"solver connect timeout")
    max_clients isa Integer && !(max_clients isa Bool) && 1<=max_clients<=1024 ||
        throw(ArgumentError("$caller: invalid max_clients"))
    option_controlled isa Bool || throw(ArgumentError("$caller: invalid option_controlled"))
    epoch=lock(STATE_LOCK) do
        _model_locked();_ONELAB_SESSION_EPOCH[]
    end
    server,actual,path,identity=_onelab_listener_bind(address)
    listener=_OnelabListener(server,actual,path,identity,current_task(),epoch,
        Threads.Atomic{Bool}(false),Threads.Atomic{Bool}(false),Threads.Atomic{Bool}(false),
        Threads.Atomic{Bool}(false),option_controlled,
        isfinite(duration) ? _onelab_monotonic_time()+duration : Inf,
        nothing,nothing,nothing,nothing,false,0)
    try
        lock(_ONELAB_RUN_LOCK) do
            _ONELAB_ACTIVE_LISTENERS[listener]=nothing
        end
        listener.watchdog=@async _onelab_listener_watchdog(listener)
        _onelab_listener_valid(listener)
        _onelab_listener_model_info!(listener)
        lock(STATE_LOCK) do
            _native_diagnostic_emit_locked(:info,"ONELAB listening at $actual")
        end
        ready===nothing || ready(actual)
        while !listener.exit_requested
            _onelab_listener_valid(listener)
            socket=try
                Sockets.accept(server)
            catch error
                _onelab_listener_valid(listener);rethrow()
            end
            try
                _onelab_listener_peer!(listener,socket,io_limit,startup,Int(max_clients))
            catch error
                error isa InterruptException && rethrow()
                _onelab_listener_valid(listener)
                listener.exit_requested && rethrow()
                lock(STATE_LOCK) do
                    _native_diagnostic_emit_locked(:warning,"ONELAB peer: "*sprint(showerror,error))
                end
            finally
                _onelab_abort_socket!(socket);try close(socket) catch end
            end
        end
        return listener.exit_code
    catch error
        error isa _OnelabListenerStopped && return 0
        rethrow()
    finally
        listener.finished[]=true;_onelab_listener_cancel!(listener)
        listener.watchdog===nothing || try wait(listener.watchdog) catch end
        lock(_ONELAB_RUN_LOCK) do
            delete!(_ONELAB_ACTIVE_LISTENERS,listener)
        end
        _onelab_listener_remove_path(path,identity)
    end
end
