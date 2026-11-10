module CliOnelabListenerChecks
using Test,Tessella,Sockets,Random
const A=Tessella.API
const C=Tessella.CLI

function frame(peer,kind,payload="")
    write(peer,htol(Int32(kind)));write(peer,htol(Int32(ncodeunits(payload))))
    write(peer,payload);flush(peer)
end
function reply(peer,messages=Tuple{Int,String}[])
    while true
        kind=Int(ltoh(read(peer,Int32)));count=Int(ltoh(read(peer,Int32)))
        0<=count<=64*1024*1024 || error("fixture byte count")
        payload=String(read(peer,count))
        if kind in (10,11,12,13)
            push!(messages,(kind,payload))
        else
            return kind,payload
        end
    end
end
function parameter(peer,kind,name,messages=Tuple{Int,String}[])
    request=A._OnelabWireParameter(kind,name)
    frame(peer,24,A._onelab_wire_encode(request))
    message,payload=reply(peer,messages)
    message==29 && return nothing
    message==23 || error("fixture parameter reply")
    return A._onelab_wire_decode(payload,"listener fixture")
end
function active_listener()
    values=lock(A._ONELAB_RUN_LOCK) do
        collect(keys(A._ONELAB_ACTIVE_LISTENERS))
    end
    return isempty(values) ? nothing : only(values)
end
function stop_server(server)
    istaskdone(server) || A.initialize()
    wait(server);A.finalize()
end

function production_case(transport)
    mktempdir() do directory
        input=joinpath(directory,"initial.geo");output=joinpath(directory,"listener.msh")
        write(input,"n=7;SetNumber(\"Listen/flag\",Solver.AlwaysListen);Point(1)={0,0,0,1};\nMacro later\n"*
            "Point(newp)={n,1,0,1};\nReturn\n")
        requested=transport===:tcp ? "127.0.0.1:0" : Sys.iswindows() ?
            "\\\\.\\pipe\\tessella-listen-"*randstring(12) : joinpath(directory,"listen.sock")
        args=[input,"-setnumber","unused","2","-listen",requested]
        server=@async try C.main(args) catch error;error end
        @test timedwait(()->active_listener()!==nothing,30.;pollint=.01)==:ok
        state=active_listener();peer=nothing
        emergency=Timer(60.) do timer;A.initialize();end
        try
            peer=transport===:tcp ? connect(Sockets.localhost,
                parse(Int,last(split(state.address,':')))) : connect(requested)
            # START is metadata, not process ownership. This is the listener
            # process's own PID and must never cause it to be terminated.
            frame(peer,1,string(getpid()))
            @test parameter(peer,:string,"Gmsh/Model name").values_s==[input]
            @test parameter(peer,:number,"Listen/flag").values_f==[1.]
            p=A._OnelabWireParameter(:number,"shared");p.values_f=[2.,3.]
            frame(peer,23,A._onelab_wire_encode(p))
            q=A._OnelabWireParameter(:string,"shared");q.values_s=["λ","tail"]
            frame(peer,23,A._onelab_wire_encode(q))
            @test parameter(peer,:number,"shared").values_f==[2.,3.]
            @test parameter(peer,:string,"shared").values_s==["λ","tail"]
            frame(peer,21,"n+=1;Call later;SetNumber(\"Listen/result\",n);"*
                "SetNumber(\"Listen/rpc flag\",Solver.AlwaysListen);"*
                "Line(1)={1,2};Transfinite Curve{1}=3;Mesh 1;"*
                "Save \""*replace(output,'\\'=>'/')*"\";Printf(\"LISTEN_RPC\");")
            messages=Tuple{Int,String}[]
            @test parameter(peer,:number,"Listen/result",messages).values_f==[8.]
            @test parameter(peer,:number,"Listen/rpc flag").values_f==[1.]
            @test any(message->message[1]==10 && message[2]=="LISTEN_RPC",messages)
            mesh=Tessella.Elements.read_mixed_msh(output)
            @test Tessella.validate(mesh).ok
            @test any(row->mesh.coords[:,row]==[8.,1.,0.],axes(mesh.coords,2))
            @test sum(size(block.nodes,2) for block in mesh.blocks if block.msh==1)==2
            frame(peer,2,"Goodbye!");close(peer);peer=nothing
            @test !istaskdone(server)
            peer=transport===:tcp ? connect(Sockets.localhost,
                parse(Int,last(split(state.address,':')))) : connect(requested)
            # Native listeners also accept a reconnect's parameter request
            # without another START; the native client retry uses this path.
            @test parameter(peer,:number,"Listen/result").values_f==[8.]
            quote_arg=Sys.iswindows() ? A._onelab_windows_quote : Base.shell_escape_posixly
            solver=joinpath(@__DIR__,"..","fixtures","onelab_native_client.jl")
            result=joinpath(directory,"direct solver.txt")
            command=join(quote_arg.(String[Base.julia_cmd().exec[1],"--startup-file=no",
                "--compile=min","-O0",solver,"child",result])," ")
            frame(peer,27,"direct solver\0"*command)
            @test reply(peer)==(2,"direct solver")
            @test readlines(result)==["direct solver"]
            @test parameter(peer,:number,"Output/child").values_f==[9.]
            frame(peer,21,"Exit;")
            @test timedwait(()->istaskdone(server),5.;pollint=.01)==:ok
            @test fetch(server)==input
            @test state.finished[] && !isopen(state.listener)
            @test isempty(A._ONELAB_ACTIVE_LISTENERS) && isempty(A._ONELAB_ACTIVE_REMOTE) &&
                isempty(A._ONELAB_ACTIVE_ENGINES)
            !Sys.iswindows() && transport===:pipe && @test !ispath(requested)
        finally
            close(emergency);peer===nothing || close(peer);stop_server(server)
        end
    end
end

function pending_write(state)
    output=state.output;output===nothing && return false
    socket=output.socket;socket===nothing && return false
    Base.iolock_begin()
    try
        socket.handle==C_NULL && return false
        return GC.@preserve socket ccall(:uv_stream_get_write_queue_size,
            Csize_t,(Ptr{Cvoid},),socket.handle)>0
    finally
        Base.iolock_end()
    end
end

@testset "Listener deadlines and reset release a surviving nonreading controller" begin
    for reset in (false,true)
        A.initialize();A.option("General.Terminal",0)
        payload=repeat("x",1024*1024)
        for index in 1:16;A.onelab.set_string("large/$index",[payload]);end
        server=@async try A._onelab_listen("127.0.0.1:0";
            timeout=15.,io_timeout=reset ? 30. : .75) catch e;e end
        @test timedwait(()->active_listener()!==nothing,5.;pollint=.01)==:ok
        state=active_listener();port=parse(Int,last(split(state.address,':')))
        peer=connect(Sockets.localhost,port)
        @test timedwait(()->state.output!==nothing,5.;pollint=.01)==:ok
        emergency_fired=Ref(false)
        emergency=Timer(5.) do timer;emergency_fired[]=true;close(peer);end
        try
            query=A._OnelabWireParameter(:string,"dummy")
            frame(peer,25,A._onelab_wire_encode(query))
            @test timedwait(()->pending_write(state),2.;pollint=.01)==:ok
            start=time()
            if reset
                A.initialize();A.onelab.set_number("replacement",[42.])
                @test timedwait(()->istaskdone(server),3.;pollint=.01)==:ok
                @test fetch(server) isa ArgumentError
                @test A.onelab.get_number("replacement")==[42.]
            else
                @test timedwait(()->state.output===nothing,3.;pollint=.01)==:ok
                @test !istaskdone(server)
                second=connect(Sockets.localhost,port)
                try
                    @test parameter(second,:number,"Gmsh/Model dimension")!==nothing
                    frame(second,21,"Exit;")
                finally
                    close(second)
                end
                @test timedwait(()->istaskdone(server),3.;pollint=.01)==:ok
                @test fetch(server)==0
            end
            @test time()-start<3.
            @test !emergency_fired[]
            @test state.finished[] && isempty(A._ONELAB_ACTIVE_LISTENERS) &&
                isempty(A._ONELAB_ACTIVE_REMOTE) && isempty(A._ONELAB_ACTIVE_ENGINES)
        finally
            close(emergency);close(peer);stop_server(server)
        end
    end
    A.initialize();A.option("General.Terminal",0)
    state=Ref{Any}(nothing)
    try
        @test_throws ArgumentError A._onelab_listen("127.0.0.1:0";timeout=.15,
            ready=address->(state[]=active_listener()))
        @test state[].finished[] && !isopen(state[].listener)
        @test isempty(A._ONELAB_ACTIVE_LISTENERS)
    finally
        A.finalize()
    end
end

@testset "CLI listener preloads geometry, meshes RPCs and retains reconnect state" begin
    production_case(:tcp);production_case(:pipe)
end

@testset "Listener default address and reset preserve replacement state" begin
    A.initialize()
    @test A.option("Solver.AlwaysListen")==0.
    @test A.option("Solver.AlwaysListen",2.9)==2.
    @test_throws ArgumentError A.option("Solver.AlwaysListen",true)
    A.finalize()
    selected,ordinary=C._onelab_cli_flags(["-listen"])
    @test selected.listening && selected.address===nothing && isempty(ordinary)
    @test_throws ArgumentError C._onelab_cli_flags(["-listen","-socket","x"])
    server=@async try C.main(["-listen","-setstring","Solver.SocketName","127.0.0.1:0"]) catch e;e end
    @test timedwait(()->active_listener()!==nothing,30.;pollint=.01)==:ok
    state=active_listener()
    try
        @test startswith(state.address,"127.0.0.1:") && !endswith(state.address,":0")
        A.initialize();A.onelab.set_number("replacement",[42.])
        @test timedwait(()->istaskdone(server),5.;pollint=.01)==:ok
        @test fetch(server) isa ArgumentError
        @test A.onelab.get_number("replacement")==[42.]
        @test state.finished[] && isempty(A._ONELAB_ACTIVE_LISTENERS)
    finally
        stop_server(server)
    end
end

@testset "Listener startup Exit and explicit disable never leave a service" begin
    mktempdir() do directory
        input=joinpath(directory,"exit.geo")
        write(input,"Exit;\n")
        @test C.main([input,"-listen","127.0.0.1:0"])==input
        @test isempty(A._ONELAB_ACTIVE_LISTENERS)
        write(input,"Exit 7;\n")
        @test_throws ArgumentError C.main([input,"-listen","127.0.0.1:0"])
        @test isempty(A._ONELAB_ACTIVE_LISTENERS)
        write(input,"Solver.AlwaysListen=0;\n")
        @test C.main([input,"-listen","127.0.0.1:0"])==input
        @test isempty(A._ONELAB_ACTIVE_LISTENERS)
    end
    for connected in (false,true)
        server=@async try C.main(["-listen","127.0.0.1:0"]) catch e;e end
        @test timedwait(()->active_listener()!==nothing,30.;pollint=.01)==:ok
        state=active_listener();peer=nothing
        try
            if connected
                peer=connect(Sockets.localhost,parse(Int,last(split(state.address,':'))))
                frame(peer,21,"Solver.AlwaysListen=0;")
            else
                A.option("Solver.AlwaysListen",0.)
            end
            @test timedwait(()->istaskdone(server),5.;pollint=.01)==:ok
            @test fetch(server)==""
            @test state.stopped[] && state.finished[]
            @test isempty(A._ONELAB_ACTIVE_LISTENERS) && isempty(A._ONELAB_ACTIVE_REMOTE)
        finally
            peer===nothing || close(peer);stop_server(server)
        end
    end
end

@testset "Listener partial frames, malformed counts and socket path rollback" begin
    A.initialize();A.option("General.Terminal",0)
    server=@async try A._onelab_listen("127.0.0.1:0";timeout=15.,io_timeout=.5) catch e;e end
    @test timedwait(()->active_listener()!==nothing,5.;pollint=.01)==:ok
    state=active_listener();port=parse(Int,last(split(state.address,':')));peer=nothing
    try
        peer=connect(Sockets.localhost,port)
        write(peer,UInt8(23));flush(peer)
        emergency_fired=Ref(false)
        emergency=Timer(3.) do timer;emergency_fired[]=true;close(peer);end
        try
            @test_throws EOFError read(peer,UInt8)
            @test !emergency_fired[]
        finally
            close(emergency);close(peer);peer=nothing
        end
        peer=connect(Sockets.localhost,port)
        write(peer,htol(Int32(23)));write(peer,htol(typemax(Int32)));flush(peer)
        @test_throws EOFError read(peer,UInt8)
        close(peer);peer=connect(Sockets.localhost,port)
        @test parameter(peer,:number,"Gmsh/Model dimension")!==nothing
        frame(peer,21,"Exit;")
        @test timedwait(()->istaskdone(server),5.;pollint=.01)==:ok
        @test fetch(server)==0
    finally
        peer===nothing || close(peer);stop_server(server)
    end
    if !Sys.iswindows()
        mktempdir() do directory
            path=joinpath(directory,"occupied.sock");write(path,"keep")
            @test_throws ArgumentError A._onelab_listener_bind(path)
            @test read(path,String)=="keep"
        end
    end
end

if Sys.isunix()
    @testset "Listener setup failure releases ownership and accepts the next peer" begin
        A.initialize();A.option("General.Terminal",0)
        server=@async try A._onelab_listen("127.0.0.1:0";timeout=15.) catch e;e end
        @test timedwait(()->active_listener()!==nothing,5.;pollint=.01)==:ok
        state=active_listener();port=parse(Int,last(split(state.address,':')))
        original=pwd()
        try
            mktempdir() do directory
                target=joinpath(directory,"removed-working-directory");mkdir(target)
                @test dirname(realpath(target))==realpath(directory)
                cd(target);rm(target) # Empty directory created by this fixture.
                peer=connect(Sockets.localhost,port)
                try
                    @test_throws EOFError read(peer,UInt8)
                finally
                    close(peer);cd(original)
                end
                @test isempty(A._ONELAB_ACTIVE_REMOTE) && isempty(A._ONELAB_ACTIVE_ENGINES)
                peer=connect(Sockets.localhost,port)
                try
                    @test parameter(peer,:number,"Gmsh/Model dimension")!==nothing
                    frame(peer,21,"Exit;")
                finally
                    close(peer)
                end
                @test timedwait(()->istaskdone(server),5.;pollint=.01)==:ok
                @test fetch(server)==0
            end
        finally
            cd(original);stop_server(server)
        end
    end
    @testset "Listener closes only its owned Unix socket pathname" begin
        A.initialize();A.option("General.Terminal",0)
        mktempdir() do directory
            path=joinpath(directory,"owned.sock");old=joinpath(directory,"moved.sock")
            server=@async try A._onelab_listen(path;timeout=15.) catch e;e end
            @test timedwait(()->active_listener()!==nothing,5.;pollint=.01)==:ok
            peer=connect(path)
            try
                mv(path,old);write(path,"replacement file")
                frame(peer,21,"Exit;")
                @test timedwait(()->istaskdone(server),5.;pollint=.01)==:ok
                @test fetch(server)==0
                @test isfile(path) && read(path,String)=="replacement file"
            finally
                close(peer);stop_server(server)
            end
        end
    end
end

@testset "Listener launches an actual remote mesher and delegated solver" begin
    A.initialize();A.option("General.Terminal",0)
    server=@async try A._onelab_listen("127.0.0.1:0";timeout=180.,connect_timeout=120.) catch e;e end
    @test timedwait(()->active_listener()!==nothing,5.;pollint=.01)==:ok
    state=active_listener();client=nothing
    try
        mktempdir() do directory
            quote_arg=Sys.iswindows() ? A._onelab_windows_quote : Base.shell_escape_posixly
            geo_quote(value)=Sys.iswindows() ? A._onelab_windows_quote(String(value)) :
                "\""*replace(String(value),'\\'=>"\\\\",'"'=>"\\\"")*"\""
            result=joinpath(directory,"solver.txt")
            solver=joinpath(@__DIR__,"..","fixtures","onelab_native_client.jl")
            command=join(geo_quote.([String(Base.julia_cmd().exec[1]),"--startup-file=no",
                "--compile=min","-O0",solver,"child",result])," ")
            source=joinpath(directory,"mesher.geo")
            write(source,"OnelabRun(\"delegated solver\", '"*command*"');"*
                "SetNumber(\"Listener/delegated\",GetNumber(\"Output/child\"));\n")
            mesher=joinpath(@__DIR__,"..","fixtures","onelab_remote_cli_client.jl")
            executable=join(quote_arg.(String[Base.julia_cmd().exec[1],"--startup-file=no",
                "--compiled-modules=no","--compile=min","-O0","--project="*dirname(Base.active_project()),
                mesher,source,"-parse_and_exit"])," ")
            client=A.onelab.client_connect(state.address;name="controller")
            @test A.onelab.client_run(client,"remote mesher",executable)===nothing
            @test readlines(result)==["delegated solver"]
            @test A.onelab.client_get_number(client,"Listener/delegated")==[9.]
            @test A.onelab.client_get_number(client,"Output/child")==[9.]
            A.onelab.client_send_command(client,"Exit;")
            @test timedwait(()->istaskdone(server),5.;pollint=.01)==:ok
            @test fetch(server)==0
            @test isempty(A._ONELAB_ACTIVE_ENGINES) && isempty(A._ONELAB_ACTIVE_REMOTE)
        end
    finally
        client===nothing || try A.onelab.client_disconnect(client) catch end
        stop_server(server)
    end
end
end
