module OnelabRemoteTransportChecks
using Test,Tessella,Sockets
const A=Tessella.API
const SocketHandle=Sys.iswindows() ? Base.Libc.WindowsRawSocket : Base.RawFD

function send_frame(peer,kind,payload="")
    write(peer,htol(Int32(kind)));write(peer,htol(Int32(ncodeunits(payload))))
    write(peer,payload);flush(peer)
end
function receive_frame(peer)
    kind=Int(ltoh(read(peer,Int32)));count=Int(ltoh(read(peer,Int32)))
    0<=count<=64*1024*1024 || error("fixture byte count")
    return kind,String(read(peer,count))
end

function with_remote(action,handler;timeout=1.,small_buffer=false)
    listener=listen(Sockets.localhost,0);peer_ref=Ref{Any}(nothing)
    if small_buffer
        buffer=Ref{Cint}(4096)
        code=Sys.iswindows() ? ccall((:setsockopt,"ws2_32"),Cint,
            (SocketHandle,Cint,Cint,Ref{Cint},Cint),fd(listener),0xffff,0x1002,buffer,4) :
            ccall(:setsockopt,Cint,(SocketHandle,Cint,Cint,Ref{Cint},Cuint),fd(listener),
                Sys.islinux() ? 1 : 0xffff,Sys.islinux() ? 8 : 0x1002,buffer,4)
        @test code==0
    end
    address="127.0.0.1:"*string(getsockname(listener)[2])
    server=@async try
        peer=accept(listener);peer_ref[]=peer
        @test first(receive_frame(peer))==1
        kind,payload=receive_frame(peer)
        @test kind==24
        @test A._onelab_wire_decode(payload,"fixture").name=="fixture/Action"
        send_frame(peer,29,"fixture/Action")
        handler(peer)
    catch error
        error isa EOFError || error isa Base.IOError || rethrow()
    end
    emergency=Timer(10.) do timer
        peer_ref[]===nothing || close(peer_ref[])
        close(listener)
    end
    session=nothing
    try
        session=A._onelab_remote_begin("fixture",address;timeout)
        return action(session)
    finally
        session===nothing || A._onelab_remote_end!(session;notify=false)
        peer_ref[]===nothing || close(peer_ref[])
        close(listener);close(emergency);wait(server)
    end
end

@testset "Remote ONELAB frames and bounded command queue" begin
    A.initialize();A.option("General.Terminal",0)
    try
        for bytes in (UInt8[0x15],reinterpret(UInt8,Int32[21,100]) |> collect,
                reinterpret(UInt8,Int32[21,-1]) |> collect,
                reinterpret(UInt8,Int32[21,typemax(Int32)]) |> collect)
            with_remote(peer->begin
                write(peer,bytes);flush(peer);close(peer)
            end) do session
                @test_throws ArgumentError A._onelab_remote_frame(session)
            end
        end
        with_remote(peer->begin
            send_frame(peer,10,"progress")
            send_frame(peer,21,"x=1;")
            p=A._OnelabWireParameter(:number,"answer");p.values_f=[7.]
            send_frame(peer,23,A._onelab_wire_encode(p));read(peer)
        end) do session
            p=A._onelab_remote_fetch(session,:number,"answer")
            @test p.values_f==[7.]
            @test session.commands==[(Int32(21),"x=1;")]
            @test session.queued_bytes==4
            empty!(session.commands);session.queued_bytes=0
            for _ in 1:1024;A._onelab_remote_queue!(session,Int32(21),"");end
            @test_throws ArgumentError A._onelab_remote_queue!(session,Int32(21),"")
            @test length(session.commands)==1024 && session.queued_bytes==0
        end
        for kind in (:wrong_name,:wrong_type)
            with_remote(peer->begin
                receive_frame(peer)
                p=A._OnelabWireParameter(kind===:wrong_type ? :string : :number,
                    kind===:wrong_name ? "other" : "answer")
                send_frame(peer,23,A._onelab_wire_encode(p));read(peer)
            end) do session
                @test_throws ArgumentError A._onelab_remote_fetch(session,:number,"answer")
            end
        end
        commands=String[]
        with_remote(peer->begin
            for _ in 1:2
                kind,payload=receive_frame(peer)
                @test kind==27
                push!(commands,payload)
                send_frame(peer,2,"nested")
                for _ in 1:2
                    @test first(receive_frame(peer))==25
                    send_frame(peer,26)
                end
            end
            read(peer)
        end) do session
            runner=A._onelab_geo_runner();runner.remote=session
            A._onelab_remote_run!(runner,"nested","actual command")
            A._onelab_remote_run!(runner,"nested","")
        end
        @test commands==["nested\0actual command","nested\0actual command"]
        messages=Tuple{Int,String}[]
        with_remote(peer->begin
            while true
                kind,payload=receive_frame(peer)
                kind==24 && (send_frame(peer,29,"barrier");break)
                push!(messages,(kind,payload))
            end
            read(peer)
        end) do session
            storage=task_local_storage()
            storage[A._ONELAB_REMOTE_KEY]=session
            storage[Tessella.IO._MSG_TASK_SINK_KEY]=A._onelab_remote_message_sink(session)
            try
                wait(@async begin
                    @test A._onelab_remote_current()===nothing
                    Tessella.IO._msg_log_write("Info    : ","unbound task")
                end)
                Tessella.IO._msg_log_write("Info    : ","bound task")
                @test A._onelab_remote_fetch(session,:number,"barrier")===nothing
            finally
                delete!(storage,A._ONELAB_REMOTE_KEY)
                delete!(storage,Tessella.IO._MSG_TASK_SINK_KEY)
            end
        end
        @test messages==[(10,"bound task")]
        @test isempty(A._ONELAB_ACTIVE_REMOTE)
    finally
        A.finalize()
    end
end

function pending_write(session)
    socket=session.socket;socket===nothing && return false
    Base.iolock_begin()
    try
        socket.handle==C_NULL && return false
        return GC.@preserve socket ccall(:uv_stream_get_write_queue_size,
            Csize_t,(Ptr{Cvoid},),socket.handle)>0
    finally
        Base.iolock_end()
    end
end

@testset "Remote ONELAB deadlines and surviving nonreading peers" begin
    A.initialize();A.option("General.Terminal",0)
    try
        for partial in (false,true)
            with_remote(peer->begin
                partial && (write(peer,UInt8(21));flush(peer))
                read(peer)
            end;timeout=.15) do session
                start=time()
                error=try A._onelab_remote_frame(session);nothing catch e;e end
                @test error isa ArgumentError
                @test occursin("timed out",sprint(showerror,error))
                @test time()-start<3.
            end
        end
        for reset in (false,true)
            ready=Channel{Any}(1);reading_stopped=Channel{Nothing}(1)
            job=@async with_remote(peer->begin
                # Repeated one-MiB frames fill the connection buffers while
                # this peer stays alive without reading. Reuse one payload.
                buffer=Ref{Cint}(4096)
                code=Sys.iswindows() ? ccall((:setsockopt,"ws2_32"),Cint,
                    (SocketHandle,Cint,Cint,Ref{Cint},Cint),fd(peer),
                    0xffff,0x1002,buffer,4) : ccall(:setsockopt,Cint,
                    (SocketHandle,Cint,Cint,Ref{Cint},Cuint),fd(peer),
                    Sys.islinux() ? 1 : 0xffff,Sys.islinux() ? 8 : 0x1002,buffer,4)
                @test code==0
                Base.stop_reading(peer)
                put!(reading_stopped,nothing)
                while isopen(peer);sleep(.01);end
            end;timeout=reset ? 30. : 1.,small_buffer=true) do session
                take!(reading_stopped)
                buffer=Ref{Cint}(4096)
                code=Sys.iswindows() ? ccall((:setsockopt,"ws2_32"),Cint,
                    (SocketHandle,Cint,Cint,Ref{Cint},Cint),fd(session.socket),
                    0xffff,0x1001,buffer,4) : ccall(:setsockopt,Cint,
                    (SocketHandle,Cint,Cint,Ref{Cint},Cuint),fd(session.socket),
                    Sys.islinux() ? 1 : 0xffff,Sys.islinux() ? 7 : 0x1001,buffer,4)
                @test code==0
                put!(ready,session)
                payload=repeat("x",1024*1024)
                try
                    A._onelab_remote_io(session) do
                        for _ in 1:256;A._onelab_remote_send(session,23,payload);end
                    end
                    nothing
                catch e
                    e
                end
            end
            session=take!(ready)
            @test timedwait(()->pending_write(session),2.;pollint=.01)==:ok
            start=time();reset && A.initialize()
            @test timedwait(()->istaskdone(job),4.;pollint=.01)==:ok
            error=fetch(job)
            @test error isa ArgumentError
            @test occursin(reset ? "reset" : "timed out",sprint(showerror,error))
            @test time()-start<3.
            @test session.finished[] && session.socket===nothing
            @test isempty(A._ONELAB_ACTIVE_REMOTE)
        end
    finally
        A.finalize()
    end
end
end
