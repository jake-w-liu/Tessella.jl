using Test
using Tessella
using Tessella.CLI: main
using Tessella.IO: read_msh
using Tessella.Elements: read_mixed_msh, mixed_crc
using Tessella.MeshTypes: mesh_crc, ntris, validate

const _SQUARE_GEO="""
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
"""

const _PERIODIC_CLI_GEO=replace(
    _SQUARE_GEO,"0, 1};"=>"0, 0.5};") *
    "Periodic Curve {2} = {4} Translate {1, 0, 0};\n"
const _PERIODIC_VOLUME_CLI_GEO=
    _PERIODIC_CLI_GEO * "Box(1) = {0, 0, 0, 1, 1, 1};\n"
const _UNUSED_PERIODIC_CLI_GEO=_SQUARE_GEO * """
Point(5) = {0, 2, 0, 0.5};
Point(6) = {1, 2, 0, 0.5};
Point(7) = {0, 3, 0, 0.5};
Point(8) = {1, 3, 0, 0.5};
Line(5) = {5, 6};
Line(6) = {7, 8};
Periodic Curve {5} = {6} Translate {0, -1, 0};
"""
const _EMBEDDED_CLI_GEO=_SQUARE_GEO * """
Point(5) = {0.25, 0.5, 0, 0.5};
Point(6) = {0.75, 0.5, 0, 0.5};
Point(7) = {0.5, 0.25, 0, 0.5};
Line(5) = {5, 6};
Point{7} In Surface{1};
Line{5} In Surface{1};
Physical Point("embedded points", 31) = {5, 6, 7};
Physical Curve("embedded line", 32) = {5};
Physical Surface("domain", 33) = {1};
"""
const _EMBEDDED_VOLUME_CLI_GEO="""
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
Point(101) = {0.2, 0.2, 0.5, 0.5};
Point(102) = {0.8, 0.2, 0.5, 0.5};
Point(103) = {0.5, 0.8, 0.5, 0.5};
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 101};
Line Loop(101) = {101, 102, 103};
Plane Surface(101) = {101};
Surface{101} In Volume{1};
Physical Point("sheet points", 51) = {101, 102, 103};
Physical Curve("sheet boundary", 52) = {101, 102, 103};
Physical Surface("sheet", 53) = {101};
Physical Volume("domain", 54) = {1};
"""

@testset "bounded non-destructive CLI" begin
    mktempdir() do directory
        input=joinpath(directory,"square")
        write(input,_SQUARE_GEO)

        @test main([input])==input
        @test read(input,String)==_SQUARE_GEO

        output=main([input,"-2"])
        @test output==input*".msh"
        @test isfile(output)
        @test read(input,String)==_SQUARE_GEO
        mesh=read_msh(output).mesh
        @test validate(mesh).ok && ntris(mesh)>0
        @test mesh_crc(mesh).sha==
              "92e578bac6d8feb3f0f845f100665dcc145edf924965ad72be716da933f34461"

        explicit=joinpath(directory,"explicit.msh")
        wrapped="x"*input*"x"
        input_view=SubString(wrapped,nextind(wrapped,firstindex(wrapped)),
                            prevind(wrapped,lastindex(wrapped)))
        @test main([input_view])==input
        @test main([input,"-o",explicit,"-2"])==explicit
        @test mesh_crc(read_msh(explicit).mesh)==mesh_crc(mesh)

        @test_throws ArgumentError main([input,"-2","-o",input])
        @test read(input,String)==_SQUARE_GEO
        link=joinpath(directory,"input-link")
        hardlink(input,link)
        @test_throws ArgumentError main([input,"-2","-o",link])
        @test read(input,String)==_SQUARE_GEO

        @test_throws ArgumentError main(String[])
        @test_throws ArgumentError main([input,input])
        @test_throws ArgumentError main([input,"--bad"])
        @test_throws ArgumentError main([input,"-2","-2"])
        @test_throws ArgumentError main([input,"-2","-3"])
        @test_throws ArgumentError main([input,"-2","-o","a","-o","b"])
        @test_throws ArgumentError main([input,"-2","-o"])
        @test_throws ArgumentError main([input,"-2","-o","-3"])
        @test_throws ArgumentError main([input,"-o","unused.msh"])
        @test_throws ArgumentError main([""])
        @test_throws ArgumentError main([input*"\0bad"])
        @test_throws ArgumentError main(fill("x",10_001))
        @test_throws ArgumentError main([repeat("x",1_000_001)])

        periodic_input=joinpath(directory,"periodic.geo")
        periodic_output=joinpath(directory,"periodic.msh")
        write(periodic_input,_PERIODIC_CLI_GEO)
        @test main([periodic_input,"-2","-o",periodic_output])==
              periodic_output
        @test read(periodic_input,String)==_PERIODIC_CLI_GEO
        periodic=read_mixed_msh(periodic_output)
        @test validate(periodic).ok
        @test periodic.entity_data!==nothing
        @test periodic.elementary_entities===nothing
        @test mixed_crc(periodic).sha==
              "cf03be1a36427f1ef0fbc4e852996bd65d2630b5ac384fa0267dd14e46ea6280"
        @test length(periodic.periodic_links)==3
        @test sort([(Int(link.slave_entity),Int(link.master_entity))
                    for link in periodic.periodic_links if link.dim==0])==
              [(2,1),(3,4)]
        periodic_curve=only(filter(
            link->link.dim==1,periodic.periodic_links))
        @test periodic_curve.slave_entity==2
        @test periodic_curve.master_entity==4
        @test length(periodic_curve.slave_nodes)==5

        embedded_input=joinpath(directory,"embedded.geo")
        embedded_output=joinpath(directory,"embedded.msh")
        write(embedded_input,_EMBEDDED_CLI_GEO)
        @test main([embedded_input,"-2","-o",embedded_output])==embedded_output
        @test read(embedded_input,String)==_EMBEDDED_CLI_GEO
        embedded=read_mixed_msh(embedded_output)
        @test validate(embedded).ok
        @test embedded.entity_data!==nothing
        @test embedded.entity_data.entities[(2,1)].embedded_curves==Int32[5]
        point_block=only(findall(block->block.msh==15,embedded.blocks))
        line_block=only(findall(block->block.msh==1,embedded.blocks))
        @test Set(embedded.entity_data.block_entities[point_block])==
              Set(Int32.(1:7))
        @test Set(embedded.entity_data.block_entities[line_block])==
              Set(Int32.(1:5))
        @test embedded.physical_names==Dict(
            (0,31)=>"embedded points",(1,32)=>"embedded line",(2,33)=>"domain")
        @test mixed_crc(embedded).sha==
              "6025846e0f58581418401081f092630d2e49999a26e608bf377d1bae4c51dc4b"

        embedded_volume_input=joinpath(directory,"embedded-volume.geo")
        embedded_volume_output=joinpath(directory,"embedded-volume.msh")
        write(embedded_volume_input,_EMBEDDED_VOLUME_CLI_GEO)
        @test main([embedded_volume_input,"-3","-o",embedded_volume_output])==
              embedded_volume_output
        @test read(embedded_volume_input,String)==_EMBEDDED_VOLUME_CLI_GEO
        embedded_volume=read_mixed_msh(embedded_volume_output)
        @test validate(embedded_volume).ok
        @test embedded_volume.entity_data!==nothing
        @test haskey(embedded_volume.entity_data.entities,(2,101))
        @test haskey(embedded_volume.entity_data.entities,(3,1))
        @test isempty(embedded_volume.entity_data.entities[(3,1)].boundaries)
        @test Set(block.msh for block in embedded_volume.blocks)==Set([15,1,2,4])
        @test embedded_volume.physical_names==Dict(
            (0,51)=>"sheet points",(1,52)=>"sheet boundary",
            (2,53)=>"sheet",(3,54)=>"domain")
        @test mixed_crc(embedded_volume).sha==
              "2fc633ca160054b1c8f86b9981febc79c83ff248afe409aad7fbb8e1e3f27ef4"

        periodic_volume_input=joinpath(directory,"periodic-volume.geo")
        periodic_volume_output=joinpath(directory,"periodic-volume.msh")
        write(periodic_volume_input,_PERIODIC_VOLUME_CLI_GEO)
        write(periodic_volume_output,"unchanged")
        projection_error=try
            main([periodic_volume_input,"-3","-o",periodic_volume_output])
            nothing
        catch err
            err
        end
        @test projection_error isa ArgumentError
        @test occursin(
            "periodic metadata projection is limited to surface meshes",
            sprint(showerror,projection_error))
        @test read(periodic_volume_input,String)==_PERIODIC_VOLUME_CLI_GEO
        @test read(periodic_volume_output,String)=="unchanged"

        unused_input=joinpath(directory,"unused-periodic.geo")
        unused_output=joinpath(directory,"unused-periodic.msh")
        write(unused_input,_UNUSED_PERIODIC_CLI_GEO)
        write(unused_output,"unchanged")
        unused_error=try
            main([unused_input,"-2","-o",unused_output])
            nothing
        catch err
            err
        end
        @test unused_error isa ArgumentError
        @test occursin(
            "selected Surface[1] does not contain periodic slave Curve tags [5]",
            sprint(showerror,unused_error))
        @test read(unused_input,String)==_UNUSED_PERIODIC_CLI_GEO
        @test read(unused_output,String)=="unchanged"
    end

    @test isempty(Docs.undocumented_names(Tessella.CLI;private=false))
end
