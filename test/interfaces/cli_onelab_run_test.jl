module CliOnelabRunChecks
using Test,Tessella
const A=Tessella.API
const C=Tessella.CLI

@testset "CLI ONELAB metamodel action lifecycle" begin
    @test C.main(["-run"])==""
    @test_throws ArgumentError C.main(["-run","-run"])
    @test_throws ArgumentError C.main(["-run","-2"])
    @test_throws ArgumentError C.main(["-run","-parse_and_exit"])
    A.initialize();A.option("General.Terminal",0)
    try
        A.onelab.set_number("session parameter",[13.,14.])
        before=A.onelab.get()
        mktempdir() do directory
            path=joinpath(directory,"metamodel.geo");output=joinpath(directory,"phases.txt")
            quote_arg(value)=Sys.iswindows() ? A._onelab_windows_quote(String(value)) :
                "\""*replace(String(value),'\\'=>"\\\\",'\"'=>"\\\"")*"\""
            fixture=joinpath(@__DIR__,"..","fixtures","onelab_native_client.jl")
            command(args...)=join(quote_arg.([String(Base.julia_cmd().exec[1]),
                "--startup-file=no","--compile=min","-O0",fixture,args...])," ")
            phase_command=command("phase",output)
            write(path,"Point(1)={0,0,0,1}; SetString(\"Phase\",OnelabAction);"*
                "n=DefineNumber(2,Name \"Controls/N\",Min 1,Max 10,Step 0.5);"*
                "OnelabRun(\"phase\",'"*phase_command*"');\n")
            bytes=read(path);files_before=sort(readdir(directory))
            @test C.main([path,"-run"])==path
            @test readlines(output)==[":2","check:7"]
            @test read(path)==bytes
            @test sort(readdir(directory))==sort([files_before;"phases.txt";"metamodel.db"])
            @test isempty(A._ONELAB_ACTIVE_ENGINES)
            @test A.onelab.get()==before

            # Option preludes and phase-sensitive field scans use the same
            # action as the geometry execution, including the check pass.
            write(path,"If(StrCmp(OnelabAction,\"check\")==0)"*
                "Mesh.MeshSizeMax=4; Else Mesh.MeshSizeMax=8; EndIf\n")
            normal=Tessella.IO.read_geo_params(path)
            check=Tessella.IO.read_geo_params(path;_onelab_action="check")
            @test normal.mesh_size_max==8. && check.mesh_size_max==4.
            result=Tessella.GeoExec.execute_geo(path;_initial_onelab_action="check")
            @test result.params.mesh_size_max==4.
            @test C.main([path,"-run","-setnumber","General.Terminal","0"])==path
            @test sort(readdir(directory))==sort([files_before;"phases.txt";"metamodel.db"])

            # Exit prevents the later action; a failed owned child leaves
            # neither an active engine nor any session parameter mutation.
            write(path,"Exit 0;\n")
            @test C.main([path,"-run"])==path
            write(path,"OnelabRun(\"bad\",'"*command("fail",output)*"');\n")
            @test_throws ArgumentError C.main([path,"-run"])
            @test isempty(A._ONELAB_ACTIVE_ENGINES)
            @test A.onelab.get()==before
            @test_throws ArgumentError C.main([path,"-run","-o",joinpath(directory,"out.msh")])
            mesh=joinpath(directory,"mesh.msh");write(mesh,"")
            @test_throws ArgumentError C.main([mesh,"-run"])
            @test_throws ArgumentError C.main([path,"-run","-merge",path])
        end
    finally
        A.finalize()
    end
end
end
