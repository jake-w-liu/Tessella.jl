module ParserLiveChecks
using Test,Tessella
const A=Tessella.API

@testset "parser file execution shares the current model and namespace" begin
    A.initialize();A.option("General.Terminal",0);A.parser.clear()
    try
        mktempdir() do directory
            path=joinpath(directory,"parse.geo")
            included=joinpath(directory,"included.geo")
            A.model.geo.add_point(2.,0.,0.;meshSize=1.,tag=2)
            A.model.geo.synchronize();A.parser.set_number("external",7.)
            write(path,"y=external+1; Point(1)={1,0,0,1};\n")
            @test A.parser.parse(path)===nothing
            @test A.model.get_entities()==[(0,1),(0,2)]
            @test A.parser.get_number("external")==[7.]
            @test A.parser.get_number("y")==[8.]
            @test A.model.get_file_name()==""

            write(included,"y+=2; label=CurrentFileName;\n"*
                "Macro create_point\nPoint(9)={y,0,0,1};\nReturn\n")
            write(path,"Include \"included.geo\";\n")
            @test A.parser.parse(path)===nothing
            @test A.parser.get_number("y")==[10.]
            @test A.parser.get_string("label")==["included.geo"]
            write(path,"Call create_point;\nDelete external;\n")
            @test A.parser.parse(path)===nothing
            @test (0,9) in A.model.get_entities()
            @test A.model.get_value(0,9,Float64[])==[10.,0,0]
            @test isempty(A.parser.get_number("external"))

            write(path,"Line(3)={1,2}; Transfinite Curve{3}=3; Mesh 1;\n")
            @test A.parser.parse(path)===nothing
            nodes=A.mesh.get_nodes();elements=A.mesh.get_elements()
            line_elements=A.mesh.get_elements(1,3)
            @test line_elements[1]==[1] && length(only(line_elements[2]))==2
            line_nodes=A.mesh.get_nodes(1,3)
            @test length(line_nodes[1])==1 && line_nodes[2]==[1.5,0.,0.]
            write(path,"y+=1; SetNumber(\"Parsed/value\",y);\n")
            @test A.parser.parse(path)===nothing
            @test A.parser.get_number("y")==[11.]
            @test A.onelab.get_number("Parsed/value")==[11.]
            @test A.mesh.get_nodes()==nodes
            @test A.mesh.get_elements()==elements

            before=A.model.get_entities();numbers=A.parser.get_number("y")
            parameters=A.onelab.get()
            write(path,"Point(90)={0,0,0,1}; y=99; SetNumber(\"Parsed/value\",99); Exit 7;\n")
            @test_throws ArgumentError A.parser.parse(path)
            @test A.model.get_entities()==before
            @test A.parser.get_number("y")==numbers
            @test A.onelab.get()==parameters
            @test A.mesh.get_nodes()==nodes && A.mesh.get_elements()==elements
            @test_throws ArgumentError A.parser.parse(joinpath(directory,"missing.geo"))
        end
    finally
        A.finalize()
    end
end

@testset "parser variables belong to the API session" begin
    A.initialize()
    try
        for finalize_first in (true,false)
            A.parser.set_number("session number",7.)
            A.parser.set_string("session string","old")
            finalize_first && A.finalize()
            A.initialize()
            @test isempty(A.parser.get_number("session number"))
            @test isempty(A.parser.get_string("session string"))
        end
    finally
        A.finalize()
    end
end
end
