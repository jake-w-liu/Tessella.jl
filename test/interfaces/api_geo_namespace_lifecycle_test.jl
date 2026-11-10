module GeoNamespaceLifecycleChecks
using Test,Tessella
const A=Tessella.API
@testset "GEO merge preserves its language and allocator; open resets them atomically" begin
    A.initialize();A.option("General.Terminal",0)
    try
        mktempdir() do directory
            first=joinpath(directory,"initial.geo")
            write(first,"n=7;xs[]={1,2};msg=\"λ\";Point(1)={0,0,0,1};\n"*
                "Macro later\nPoint(newp)={n,0,0,1};\nReturn\n")
            second=joinpath(directory,"merge.geo")
            write(second,"n+=1;v=external[1];xs[1]=9;msg=StrCat(msg,\"tail\");Call later;\n")
            @test A.open(first)===nothing
            A.parser.set_number("external",[11.,12.]);A.onelab.set_number("keep",[4.])
            @test A.merge(second)===nothing
            @test A.parser.get_number("n")==[8.]
            @test A.parser.get_number("v")==[12.]
            @test A.parser.get_number("xs")==[1.,9.]
            @test A.parser.get_string("msg")==["λtail"]
            @test A.model.get_entities()==[(0,1),(0,2)]
            @test A.model.get_value(0,2,Float64[])==[8.,0.,0.]
            before=A.model.get_entities()
            write(second,"n=99;SetNumber(\"keep\",99);Point(99)={99,0,0,1};Exit 7;\n")
            @test_throws ArgumentError A.merge(second)
            @test A.parser.get_number("n")==[8.] && A.model.get_entities()==before
            @test A.onelab.get_number("keep")==[4.]
            mesh=joinpath(directory,"point.msh");write(mesh,"invalid")
            @test_throws ArgumentError A.open(mesh)
            @test A.parser.get_number("n")==[8.] && A.model.get_entities()==before
            write(mesh,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n1\n"*
                "1 0 0 0\n\$EndNodes\n\$Elements\n1\n1 15 2 0 1 1\n\$EndElements\n")
            @test A.open(mesh)===nothing
            @test isempty(A.parser.get_number("n")) && isempty(A.parser.get_string("msg"))
            @test isnothing(A._ONELAB_LIVE_GEO[])
            @test A.onelab.get_number("keep")==[4.]
            write(second,"Point(newp)={2,0,0,1};\n")
            @test A.parser.parse(second)===nothing
            @test A.model.get_entities()==[(0,1),(0,2)]
            @test A.model.get_value(0,2,Float64[])==[2.,0.,0.]
        end
    finally
        A.finalize()
    end
end

@testset "Matched MSH open resets GEO namespace while merge preserves it" begin
    for operation in (A.open,A.merge)
        A.initialize();A.option("General.Terminal",0)
        try
            mktempdir() do directory
                geo=joinpath(directory,"initial.geo")
                write(geo,"n=7;msg=\"λ\";Point(1)={0,0,0,1};\n"*
                    "Macro later\nPoint(newp)={n,0,0,1};\nReturn\n")
                mesh=joinpath(directory,"point.msh")
                A.open(geo);A.option("Geometry.MatchGeomAndMesh",1)
                A.parser.set_number("external",[11.,12.])
                A.onelab.set_number("keep",[4.])
                previous_namespace=A._ONELAB_LIVE_GEO[]
                write(mesh,"invalid")
                @test_throws ArgumentError operation(mesh)
                @test A.parser.get_number("n")==[7.]
                @test A.parser.get_string("msg")==["λ"]
                @test A._ONELAB_LIVE_GEO[]===previous_namespace
                write(mesh,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n1\n"*
                    "1 0 0 0\n\$EndNodes\n\$Elements\n1\n1 15 2 0 1 1\n\$EndElements\n")
                before=A.model.list()
                @test operation(mesh)===nothing
                @test A.model.list()==before
                @test A.model.get_entities()==[(0,1)]
                @test first(A.mesh.get_nodes())==[1]
                @test A.onelab.get_number("keep")==[4.]
                second=joinpath(directory,"next.geo")
                if operation===A.open
                    @test isempty(A.parser.get_number("n"))
                    @test isempty(A.parser.get_number("external"))
                    @test isempty(A.parser.get_string("msg"))
                    @test isnothing(A._ONELAB_LIVE_GEO[])
                    write(second,"Point(newp)={2,0,0,1};\n")
                else
                    @test A.parser.get_number("n")==[7.]
                    @test A.parser.get_number("external")==[11.,12.]
                    @test A.parser.get_string("msg")==["λ"]
                    @test A._ONELAB_LIVE_GEO[]===previous_namespace
                    write(second,"n+=1;Call later;\n")
                end
                @test A.parser.parse(second)===nothing
                @test A.model.get_entities()==[(0,1),(0,2)]
                @test A.model.get_value(0,2,Float64[])==
                    [operation===A.open ? 2. : 8.,0.,0.]
            end
        finally
            A.finalize()
        end
    end
end
end
