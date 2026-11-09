using Test
using Tessella
import Sockets

function _onelab_truncated_allocation(reader,count)
    listener=Sockets.listen(Sockets.localhost,0)
    client=Sockets.connect(Sockets.localhost,Sockets.getsockname(listener)[2])
    peer=Sockets.accept(listener);close(listener)
    write(peer,"short");close(peer)
    try
        return @allocated try reader(client,count) catch exception
            exception isa EOFError || rethrow()
        end
    finally
        close(client)
    end
end

@testset "ONELAB bounded wire decoder" begin
    A=Tessella.API
    @test A._onelab_exit_status(0,0)==0
    @test A._onelab_exit_status(7,0)==7
    @test A._onelab_exit_status(0,15)==143
    negative=join(["1.3","string","bad","","","31","1","0","-1","-1","-1","generic","-1"],'\0')
    @test_throws ArgumentError A._onelab_wire_decode(negative,"test")
    parameter=A._OnelabWireParameter(:string,"empty")
    parameter.values_s=[""]
    wire=A._onelab_wire_encode(parameter)
    @test A._onelab_wire_decode(wire,"test").values_s==[""]
    @test_throws ArgumentError A._onelab_wire_decode(wire*"\0extra","test")
    parameter.values_s=["embedded\0nul"]
    @test_throws ArgumentError A._onelab_wire_encode(parameter)
    # Strict UTF-8 classification, including overlongs, surrogates and
    # out-of-range code points, on both the ASCII fast path and the fallback.
    for bytes in (UInt8[], UInt8[0x68,0x69], UInt8[0xc3,0xa9],
            UInt8[0xf0,0x90,0x80,0x80], UInt8[0xff], UInt8[0xc3],
            UInt8[0xc0,0x80], UInt8[0xe0,0x80,0x80],
            UInt8[0xed,0xa0,0x80], UInt8[0xf4,0x90,0x80,0x80],
            UInt8[0xf5,0x80,0x80,0x80])
        text=String(bytes)
        @test A._onelab_valid_text(text)==isvalid(text)
    end
    A._onelab_dict_put!(parameter.attributes,parameter.attr_order,"removed","blue")
    empty!(parameter.attributes);parameter.values_s=[""]
    @test isempty(A._onelab_wire_decode(A._onelab_wire_encode(parameter),"test").attributes)
    for reader in (A._onelab_client_recv_exact,A._onelab_server_read_exact)
        _onelab_truncated_allocation(reader,8)
        small=_onelab_truncated_allocation(reader,1024)
        huge=_onelab_truncated_allocation(reader,64*1024*1024)
        @test huge<small+32768
        @test huge<200000
    end
end

@testset "native ONELAB process server" begin
    A=Tessella.API
    quote_arg=Sys.iswindows() ? A._onelab_windows_quote : Base.shell_escape_posixly
    fixture=joinpath(@__DIR__,"..","fixtures","onelab_native_client.jl")
    command(args...)=join(quote_arg.(String[Base.julia_cmd().exec[1],
        "--startup-file=no","--compile=min","-O0",fixture,args...])," ")
    A.initialize();A.option("General.Terminal",0)
    try
        mktempdir() do directory
            output=joinpath(directory,"client result λ.txt")
            A.onelab.set_number("Input/seed",[2.,3.])
            @test A.onelab.run("solver λ\"tail\\",command("basic",output);timeout=30)===nothing
            @test read(output,String)=="solver λ\"tail\\\n2"
            @test A.onelab.get_number("Output/λ")==[4.,5.]
            @test A.onelab.get_string("Output/text")==["hello λ","second"]
            @test occursin("result λ",A.onelab.get("Output/λ"))
            @test occursin("Blue",A.onelab.get("Output/λ"))
            @test !A.onelab.get_changed("solver λ\"tail\\")
            A.onelab.set_number("Input/seed",6.)
            @test A.onelab.get_changed("solver λ\"tail\\")
            A.onelab.clear_changed()
            @test !A.onelab.get_changed("solver λ\"tail\\")
            @test A.onelab.run("solver λ\"tail\\";timeout=30)===nothing
            @test isempty(A._ONELAB_ACTIVE_ENGINES)

            # Use the product's native client against the product's process
            # server. The tiny owned child keeps the process lifetime open.
            relay=joinpath(directory,"native relay.txt")
            runner=@async A.onelab.run("native relay",command("relay",relay);timeout=30)
            @test timedwait(()->isfile(relay),15.;pollint=.01)==:ok
            client=A.onelab.client_connect(read(relay,String);name="native relay")
            try
                @test A.onelab.client_get_number(client,"Input/seed")==[6.]
                A.onelab.set_number("Output/λ",[4.,5.])
                replacement=A.onelab.client_get_parameter(client,:number,"Output/λ")
                @test isempty(replacement.attributes) && isempty(replacement.value_labels)
                @test A.onelab.client_define_number(client,"Input/new";values=[1.,2.],min=0.,max=8.)==[1.,2.]
                A.onelab.client_set_number(client,"Input/new",[3.,4.])
                @test A.onelab.client_get_number(client,"Input/new")==[3.,4.]
                A.onelab.client_set_string(client,"Input/string",["one","two"])
                @test A.onelab.client_get_string(client,"Input/string")==["one","two"]
                @test "Input/new" in A.onelab.client_get_names(client)
                A.onelab.client_set_changed(client,"native relay",false)
                @test !A.onelab.client_get_changed(client,"native relay")
                A.onelab.client_clear(client,"Input/string")
                @test !("Input/string" in A.onelab.client_get_names(client))
            finally
                A.onelab.client_disconnect(client)
                write(relay*".done","")
            end
            @test fetch(runner)===nothing
            @test A.onelab.get_number("Input/new")==[3.,4.]

            child=command("child",joinpath(directory,"children.txt"))
            @test A.onelab.run("parent",command("nested",output,child);
                timeout=30,max_clients=2)===nothing
            @test read(output,String)=="nested complete"
            @test readlines(joinpath(directory,"children.txt"))==["nested child","nested child"]
            @test A.onelab.get_number("Output/child")==[9.]

            # ParseString uses the full native GEO evaluator and publishes its
            # actual model/parameter changes before the next wire query.
            source="n=3; For i In {1:n} Point(i)={i,0,0,1}; EndFor SetNumber(\"Output/parsed\",n);"
            @test A.onelab.run("parse",command("parse",output,source);timeout=30)===nothing
            @test read(output,String)=="parsed"
            @test A.model.get_entities(0)==[(0,1),(0,2),(0,3)]
            @test A.onelab.get_number("Output/parsed")==[3.]
            @test A.onelab.run("parse",command("parse",output,
                "Point(4)={n+1,0,0,1}; SetNumber(\"Output/parsed\",n+1);");timeout=30)===nothing
            @test A.onelab.get_number("Output/parsed")==[4.]
            @test A.onelab.run("parse",command("parse",output,"Delete n;");timeout=30)===nothing
            @test isempty(A.parser.get_number("n"))
            @test A.onelab.run("parse",command("parse",output,
                "SetNumber(\"Output/deleted\",Exists(n));");timeout=30)===nothing
            @test A.onelab.get_number("Output/deleted")==[0.]
            A.parser.set_number("external",8.)
            @test A.onelab.run("parse",command("parse",output,
                "SetNumber(\"Output/external\",external);");timeout=30)===nothing
            @test A.onelab.get_number("Output/external")==[8.]
            A.parser.clear("external")
            @test A.onelab.run("parse",command("parse",output,
                "SetNumber(\"Output/external\",Exists(external));");timeout=30)===nothing
            @test A.onelab.get_number("Output/external")==[0.]

            file=joinpath(directory,"incoming.geo")
            write(file,"Point(20)={20,0,0,1}; SetNumber(\"Output/file\",12);\n")
            @test A.onelab.run("merge file",command("file",output,"20",file);timeout=30)===nothing
            @test A.model.get_entities(0)==[(0,1),(0,2),(0,3),(0,4),(0,20)]
            @test A.onelab.get_number("Output/file")==[12.]
            write(file,"Point(30)={30,0,0,1};\n")
            @test A.onelab.run("open file",command("file",output,"33",file);timeout=30)===nothing
            @test A.model.get_entities(0)==[(0,30)]

            # File execution stages the full store. Merely reading the first
            # value preserves vectors, and a later failure adopts no writes.
            A.onelab.set_number("Input/seed",[2.,3.])
            geo=joinpath(directory,"runner.geo")
            run_command=command("basic",output)
            write(geo,"OnelabRun(\"geo client\",'$run_command');\n"*
                "SetNumber(\"Output/geo\",GetNumber(\"Output/λ\"));\n")
            @test A.merge(geo)===nothing
            @test startswith(read(output,String),"geo client\n")
            @test A.onelab.get_number("Input/seed")==[2.,3.]
            @test A.onelab.get_number("Output/λ")==[4.,5.]
            @test A.onelab.get_number("Output/geo")==[4.]
            write(geo,"OnelabRun \"geo client\";\n")
            @test A.merge(geo)===nothing
            A.onelab.set_number("Output/λ",[14.,15.])
            before=A.onelab.get()
            write(geo,"SetNumber(\"Input/seed\",99); OnelabRun(\"geo client\"); Exit 7;\n")
            @test_throws ArgumentError A.merge(geo)
            @test A.onelab.get()==before

            A.onelab.set_number("Input/empty",Float64[])
            write(geo,"value=GetNumber(\"Input/seed\");\n")
            result=Base.invokelatest(A.open_geo!,geo)
            @test result isa Tessella.GeoExec.GeoExecution
            @test A.onelab.get_number("Input/seed")==[2.,3.]
            @test A.onelab.get_number("Input/empty")==Float64[]
            A.onelab.set_string("Input/empty string",String[])
            write(geo,"x=DefineNumber(9,Name \"Input/empty\");\n"*
                "y=GetNumber(\"Input/empty\",5);\n"*
                "s=DefineString(\"fallback\",Name \"Input/empty string\");\n")
            result=Base.invokelatest(A.open_geo!,geo)
            @test result.values["x"]==9. && result.values["y"]==0.
            @test result.strings["s"]==[""]
            @test isempty(A.onelab.get_number("Input/empty"))
            @test isempty(A.onelab.get_string("Input/empty string"))

            # The native preprocessor consumes plain solver text, not GEO.
            # Its reply reports whether it resolved a parameter.
            template=joinpath(directory,"solver input.txt.ol")
            write(template,"value=OL.get(Input/seed) text=OL.get(Output/text)\n")
            @test A.onelab.run("preprocess",command("preprocess",output,template);timeout=30)===nothing
            @test read(output,String)=="(28, \"true\")"
            @test read(chop(template;tail=3),String)=="value=2 text=hello λ\n\n"
            write(template,"plain solver text\n")
            @test A.onelab.run("preprocess";timeout=30)===nothing
            @test read(output,String)=="(28, \"false\")"
            @test read(chop(template;tail=3),String)=="plain solver text\n\n"
            A.onelab.set_number("Input/zero",0.)
            write(joinpath(directory,"included.ol"),"included=OL.get(Input/seed)\n")
            write(template,"OL.iftrue(Input/zero)\nwrong\nOL.else\nOL.include(included.ol)\nOL.endif\n"*
                "OL.ifntrue(Input/zero)\nright\nOL.endif\n")
            @test A.onelab.run("preprocess";timeout=30)===nothing
            @test read(chop(template;tail=3),String)=="included=2\n\nright\n\n"
            retained=read(chop(template;tail=3),String)
            write(template,"OL.include(solver input.txt.ol)\n")
            @test_throws ArgumentError A.onelab.run("preprocess";timeout=30)
            @test read(chop(template;tail=3),String)==retained
            @test isempty(A._ONELAB_ACTIVE_ENGINES)

            for mode in ("fast","fast","fast")
                @test A.onelab.run("quick",command(mode,output);timeout=30)===nothing
            end
            for (mode,pattern) in (("fail","status 7"),("nostart","without connecting"),
                    ("nostop","without STOP"),("badlength","frame length"),
                    ("partial","without STOP"))
                exception=try A.onelab.run(mode,command(mode,output);timeout=30,
                    max_message_bytes=2*1024*1024);nothing catch caught;caught end
                @test exception isa ArgumentError
                @test occursin(pattern,sprint(showerror,exception))
                @test isempty(A._ONELAB_ACTIVE_ENGINES)
            end
            if Sys.isunix()
                @test_throws ArgumentError A.onelab.run("signal",command("signal",output);timeout=30)
            end
            @test_throws ArgumentError A.onelab.run("missing","tessella-no-such-program-20261109";timeout=3)
            @test isempty(A._ONELAB_ACTIVE_ENGINES)
            @test_throws ArgumentError A.onelab.run("hung",command("hang",output);timeout=1.)
            @test isempty(A._ONELAB_ACTIVE_ENGINES)
            @test_throws ArgumentError A.onelab.run("limits",command("fast",output);max_clients=0)

            # Reset cancels only this runner's exact process tree, even while
            # its event pump is waiting. Retain handles before the reset.
            rm(output;force=true)
            task=@async try
                A.onelab.run("tree",command("tree",output);timeout=30)
                nothing
            catch exception
                exception
            end
            @test timedwait(()->isfile(output*".descendant"),15.;pollint=.01)==:ok
            pids=parse.(Int,[read(output,String),read(output*".descendant",String)])
            handles=Sys.iswindows() ? [ccall((:OpenProcess,"kernel32"),Ptr{Cvoid},
                (UInt32,Int32,UInt32),0x00100000,0,pid) for pid in pids] : Ptr{Cvoid}[]
            A.initialize()
            @test timedwait(()->istaskdone(task),10.;pollint=.01)==:ok
            @test fetch(task) isa ArgumentError
            @test isempty(A._ONELAB_ACTIVE_ENGINES)
            for handle in handles
                try
                    @test handle!=C_NULL
                    @test ccall((:WaitForSingleObject,"kernel32"),UInt32,
                        (Ptr{Cvoid},UInt32),handle,5000)==0
                finally
                    handle==C_NULL || ccall((:CloseHandle,"kernel32"),Int32,(Ptr{Cvoid},),handle)
                end
            end
            @test A.onelab.run()===nothing

            # Backpressure: a client that requests a large value and then
            # never reads must not extend the run deadline. The watchdog kills
            # the owned process before closing its socket, which releases the
            # pump's blocked response write.
            A.onelab.set_string("Input/huge",[repeat("x",4*1024*1024)])
            before=time()
            exception=try
                A.onelab.run("no-read",command("noread",output);
                    timeout=1.,connect_timeout=3.)
                nothing
            catch caught
                caught
            end
            @test exception isa ArgumentError
            @test occursin("timed out",sprint(showerror,exception))
            @test time()-before<4.
            @test isempty(A._ONELAB_ACTIVE_ENGINES)

            # Reset while the pump is blocked in that response write must also
            # finish promptly and clean up the owned client.
            rm(output;force=true)
            task=@async try
                A.onelab.run("no-read reset",command("noread",output);timeout=30)
                nothing
            catch caught
                caught
            end
            @test timedwait(()->!isempty(A._ONELAB_ACTIVE_ENGINES),15.;
                pollint=.01)==:ok
            sleep(.5)
            start=time()
            A.initialize()
            @test timedwait(()->istaskdone(task),10.;pollint=.01)==:ok
            @test fetch(task) isa ArgumentError
            @test time()-start<5.
            @test isempty(A._ONELAB_ACTIVE_ENGINES)

            # The standalone frontend owns its own parameter store and does
            # not require or initialize the process-global API session.
            A.finalize()
            write(geo,"SetNumber(\"Input/seed\",2); OnelabRun(\"standalone\",'$run_command');\n"*
                "result=GetNumber(\"Output/λ\");\n")
            result=Base.invokelatest(Tessella.GeoExec.execute_geo,geo;mesh_dim=0)
            @test result.values["result"]==4.
            @test startswith(read(output,String),"standalone\n")
            @test A.CURRENT[]===nothing
        end
    finally
        A.finalize()
    end
end
