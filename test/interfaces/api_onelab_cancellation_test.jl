module OnelabCancellationChecks
using Test,Tessella,Sockets
const A=Tessella.API
const SocketHandle=Sys.iswindows() ? Base.Libc.WindowsRawSocket : Base.RawFD

function relay_peer(directory)
    output=joinpath(directory,"separate client.txt")
    fixture=joinpath(@__DIR__,"..","fixtures","onelab_native_client.jl")
    quote_arg=Sys.iswindows() ? A._onelab_windows_quote : Base.shell_escape_posixly
    command=join(quote_arg.(String[Base.julia_cmd().exec[1],"--startup-file=no",
        "--compile=min","-O0",fixture,"relay",output])," ")
    return output,command
end

function connect_peer(output)
    @test timedwait(()->isfile(output),15.;pollint=.01)==:ok
    host,port=split(read(output,String),':')
    peer=connect(host,parse(Int,port))
    buffer=Ref{Cint}(4096)
    result=if Sys.iswindows()
        ccall((:setsockopt,"ws2_32"),Cint,
            (SocketHandle,Cint,Cint,Ref{Cint},Cint),fd(peer),0xffff,0x1002,buffer,4)
    else
        ccall(:setsockopt,Cint,(Base.RawFD,Cint,Cint,Ref{Cint},Cuint),fd(peer),
            Sys.islinux() ? 1 : 0xffff,Sys.islinux() ? 8 : 0x1002,buffer,4)
    end
    @test result==0
    for (kind,payload) in ((1,string(getpid())),(25,"1.3\0string"))
        write(peer,htol(Int32(kind)));write(peer,htol(Int32(ncodeunits(payload))))
        write(peer,payload);flush(peer)
    end
    return peer
end

function pending_write(engine)
    Base.iolock_begin()
    try
        for job in engine.jobs
            socket=job.socket
            socket===nothing && continue
            socket.handle==C_NULL && continue
            queued=GC.@preserve socket ccall(:uv_stream_get_write_queue_size,
                Csize_t,(Ptr{Cvoid},),socket.handle)
            queued>0 && return true
        end
        return false
    finally
        Base.iolock_end()
    end
end

@testset "ONELAB cancellation with surviving socket peer" begin
    A.initialize();A.option("General.Terminal",0)
    try
        mktempdir() do directory
            for reset in (false,true)
                A.onelab.set_string("huge",[repeat("x",4*1024*1024)])
                output,command=relay_peer(directory)
                rm(output;force=true)
                runner=@async try
                    A.onelab.run("separate client",command;
                        timeout=reset ? 30. : 3.,connect_timeout=10.)
                    nothing
                catch exception
                    exception
                end
                peer=connect_peer(output)
                engine=lock(A._ONELAB_RUN_LOCK) do
                    only(keys(A._ONELAB_ACTIVE_ENGINES))
                end
                emergency_fired=Ref(false)
                emergency=Timer(6.) do timer
                    emergency_fired[]=true
                    close(peer)
                end
                try
                    @test timedwait(()->pending_write(engine),2.;pollint=.01)==:ok
                    start=time()
                    reset && A.initialize()
                    @test timedwait(()->istaskdone(runner),5.;pollint=.01)==:ok
                    result=fetch(runner)
                    @test result isa ArgumentError
                    @test occursin(reset ? "cancelled" : "timed out",sprint(showerror,result))
                    @test time()-start<4.
                    @test !emergency_fired[]
                    @test isempty(A._ONELAB_ACTIVE_ENGINES)
                    @test all(job->job.process.closed && istaskdone(job.reader),engine.jobs)
                    @test all(job->!isopen(job.socket) && !isopen(job.listener),engine.jobs)
                finally
                    close(emergency);isopen(peer) && close(peer)
                end
            end
        end
    finally
        A.finalize()
    end
end
end
