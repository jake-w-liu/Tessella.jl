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
const _PERIODIC_SURFACE_VOLUME_CLI_GEO=read(normpath(joinpath(
    @__DIR__,"..","fixtures","periodic_surface_volume.geo")),String)
const _GEOMETRY_EXPRESSION_CLI_GEO=read(normpath(joinpath(
    @__DIR__,"..","fixtures","geo_geometry_expressions.geo")),String)
const _LIST_VARIABLE_CLI_GEO=read(normpath(joinpath(
    @__DIR__,"..","fixtures","geo_list_variables.geo")),String)
const _DYNAMIC_TAG_CLI_GEO=read(normpath(joinpath(
    @__DIR__,"..","fixtures","geo_dynamic_tags.geo")),String)
const _SET_MAX_TAG_CLI_GEO=read(normpath(joinpath(
    @__DIR__,"..","fixtures","geo_set_max_tags.geo")),String)
const _POINT_MESH_SIZE_CLI_GEO=read(normpath(joinpath(
    @__DIR__,"..","fixtures","geo_point_mesh_sizes.geo")),String)
const _POINT_MESH_SIZE_POINTS_OF_CLI_GEO=read(normpath(joinpath(
    @__DIR__,"..","fixtures","geo_point_mesh_size_points_of.geo")),String)
const _PERIODIC_EMBEDDED_CLI_GEO=_SQUARE_GEO * """
Point(5) = {0.25, 0.25, 0, 1};
Point(6) = {0.75, 0.25, 0, 1};
Point(7) = {0.25, 0.75, 0, 1};
Point(8) = {0.75, 0.75, 0, 1};
Point(9) = {0.5, 0.25, 0, 1};
Line(5) = {5, 6};
Line(6) = {7, 8};
Point{9} In Surface{1};
Line{5, 6} In Surface{1};
Periodic Curve {6} = {5} Translate {0, 0.5, 0};
Physical Curve("periodic traces", 21) = {5, 6};
Physical Surface("domain", 22) = {1};
"""
const _PERIODIC_GRAPH_CLI_GEO="""
Point(1) = {0, 0, 0, 0.4};
Point(2) = {1, 0, 0, 0.4};
Point(3) = {1, 1, 0, 0.4};
Point(4) = {0, 1, 0, 0.4};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Point(101) = {0.2, 0.2, 0, 0.4};
Point(102) = {0.8, 0.2, 0, 0.4};
Point(103) = {0.2, 0.5, 0, 0.4};
Point(104) = {0.8, 0.5, 0, 0.4};
Point(105) = {0.2, 0.8, 0, 0.4};
Point(106) = {0.8, 0.8, 0, 0.4};
Point(107) = {0.425, 0.8, 0, 0.4};
Line(30) = {101, 102};
Line(20) = {103, 104};
Line(10) = {106, 105};
Point{107} In Surface{1};
Line{30, 20, 10} In Surface{1};
graphSlaveBegin = Sqrt(100);
graphSlaveEnd = 4 * 5;
graphMasterBegin = graphSlaveEnd;
graphMasterEnd = 3 * 10;
graphCurveStep = 20 / 2;
graphShift = 3 / 10;
graphZero = Atan2(0, 1);
Periodic Curve {graphSlaveBegin:graphSlaveEnd:graphCurveStep} =
  {graphMasterBegin:graphMasterEnd:graphCurveStep}
  Translate {graphZero, graphShift, Sin(graphZero)};
Physical Curve("periodic traces", 41) = {30, 20, 10};
Physical Surface("domain", 42) = {1};
"""
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
Point(104) = {0.3, 0.35, 0.5, 0.5};
Point(105) = {0.7, 0.35, 0.5, 0.5};
Point(106) = {0.5, 0.35, 0.5, 0.5};
Line(104) = {104, 105};
Point{106} In Surface{101};
Line{104} In Surface{101};
Surface{101} In Volume{1};
Physical Point("sheet points", 51) = {101, 102, 103, 104, 105, 106};
Physical Curve("sheet curves", 52) = {101, 102, 103, 104};
Physical Surface("sheet", 53) = {101};
Physical Volume("domain", 54) = {1};
"""
const _EXPLICIT_SHELL_CLI_GEO="""
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {0, 1, 0, 1};
Point(4) = {0, 0, 1, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {2, 4};
Line(6) = {3, 4};
Curve Loop(1) = {1, 2, 3};
Curve Loop(2) = {1, 5, -4};
Curve Loop(3) = {2, 6, -5};
Curve Loop(4) = {3, 4, -6};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Surface Loop(1) = {1, 2, 3, 4};
Volume(1) = {1};
Physical Surface("boundary", 71) = {1, 2, 3, 4};
Physical Volume("domain", 72) = {1};
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

        periodic_embedded_input=joinpath(
            directory,"periodic-embedded.geo")
        periodic_embedded_output=joinpath(
            directory,"periodic-embedded.msh")
        write(periodic_embedded_input,_PERIODIC_EMBEDDED_CLI_GEO)
        @test main([periodic_embedded_input,"-2","-o",
                    periodic_embedded_output])==periodic_embedded_output
        @test read(periodic_embedded_input,String)==
              _PERIODIC_EMBEDDED_CLI_GEO
        periodic_embedded=read_mixed_msh(periodic_embedded_output)
        @test validate(periodic_embedded).ok
        @test periodic_embedded.entity_data.entities[(2,1)].embedded_curves==
              Int32[5,6]
        @test periodic_embedded.physical_names==Dict(
            (1,21)=>"periodic traces",(2,22)=>"domain")
        @test sort([(Int(link.slave_entity),Int(link.master_entity))
                    for link in periodic_embedded.periodic_links
                    if link.dim==0])==[(7,5),(8,6)]
        embedded_periodic_curve=only(filter(
            link->link.dim==1,periodic_embedded.periodic_links))
        @test embedded_periodic_curve.slave_entity==6
        @test embedded_periodic_curve.master_entity==5
        @test length(embedded_periodic_curve.slave_nodes)==3
        @test mixed_crc(periodic_embedded).sha==
              "b02e4da6aa4910bfb488f6b3eebc4ae0bf22f91867506a973390cc65f450d9ec"

        periodic_graph_input=joinpath(directory,"periodic-graph.geo")
        periodic_graph_output=joinpath(directory,"periodic-graph.msh")
        write(periodic_graph_input,_PERIODIC_GRAPH_CLI_GEO)
        @test main([periodic_graph_input,"-2","-o",
                    periodic_graph_output])==periodic_graph_output
        @test read(periodic_graph_input,String)==_PERIODIC_GRAPH_CLI_GEO
        periodic_graph=read_mixed_msh(periodic_graph_output)
        @test validate(periodic_graph).ok
        @test periodic_graph.entity_data.entities[(2,1)].embedded_curves==
              Int32[10,20,30]
        @test periodic_graph.physical_names==Dict(
            (1,41)=>"periodic traces",(2,42)=>"domain")
        @test sort([(Int(link.slave_entity),Int(link.master_entity))
                    for link in periodic_graph.periodic_links
                    if link.dim==0])==
              [(103,101),(104,102),(105,103),(106,104)]
        graph_curve_links=Dict(
            Int(link.slave_entity)=>link for link in
            periodic_graph.periodic_links if link.dim==1)
        @test Dict(slave=>Int(link.master_entity)
                   for (slave,link) in graph_curve_links)==
              Dict(10=>20,20=>30)
        @test all(link->length(link.slave_nodes)==9,
                  values(graph_curve_links))
        @test mixed_crc(periodic_graph).sha==
              "3f98267cc70f9326ebe490c854cb59a9987c638e6aaabcba326d086bfb887ab1"

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
              "ac8238a9530b41f1cffbb596f304f77e07a047f2f7db0695bc25770af117448e"

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
        @test embedded_volume.entity_data.entities[(2,101)].embedded_curves==
              Int32[104]
        @test isempty(embedded_volume.entity_data.entities[(3,1)].boundaries)
        @test Set(block.msh for block in embedded_volume.blocks)==Set([15,1,2,4])
        @test embedded_volume.physical_names==Dict(
            (0,51)=>"sheet points",(1,52)=>"sheet curves",
            (2,53)=>"sheet",(3,54)=>"domain")
        @test mixed_crc(embedded_volume).sha==
              "745bc23ab2aa7c0824006a94ef279514a1c6fa97d3cd95960943797b85c6336c"

        explicit_shell_input=joinpath(directory,"explicit-shell.geo")
        explicit_shell_output=joinpath(directory,"explicit-shell.msh")
        write(explicit_shell_input,_EXPLICIT_SHELL_CLI_GEO)
        @test main([explicit_shell_input,"-3","-o",explicit_shell_output])==
              explicit_shell_output
        @test read(explicit_shell_input,String)==_EXPLICIT_SHELL_CLI_GEO
        explicit_shell=read_mixed_msh(explicit_shell_output)
        @test validate(explicit_shell).ok
        @test explicit_shell.entity_data!==nothing
        @test explicit_shell.entity_data.entities[(3,1)].boundaries==
              Int32[1,2,3,4]
        @test Set(block.msh for block in explicit_shell.blocks)==Set([15,1,2,4])
        @test explicit_shell.physical_names==Dict(
            (2,71)=>"boundary",(3,72)=>"domain")
        @test mixed_crc(explicit_shell).sha==
              "fe4bc18f3b9156c654c5ee1433b43c87cb9b8d7f372d9fbb72689f500584623a"

        geometry_expression_input=joinpath(
            directory,"geometry-expressions.geo")
        geometry_expression_output=joinpath(
            directory,"geometry-expressions.msh")
        write(geometry_expression_input,_GEOMETRY_EXPRESSION_CLI_GEO)
        @test main([geometry_expression_input,"-3","-o",
                    geometry_expression_output])==geometry_expression_output
        @test read(geometry_expression_input,String)==
              _GEOMETRY_EXPRESSION_CLI_GEO
        geometry_expression=read_mixed_msh(geometry_expression_output)
        @test validate(geometry_expression).ok
        @test geometry_expression.physical_names==Dict(
            (0,70)=>"corners",(0,71)=>"probe",(1,72)=>"edges",
            (2,73)=>"boundary",(3,74)=>"domain")
        @test mixed_crc(geometry_expression).sha==
              "89ee7d39873b202e264917e98fce2756b038d3e6f2f06d1bdbce9c46f7e628cd"

        list_variable_input=joinpath(directory,"list-variables.geo")
        list_variable_output=joinpath(directory,"list-variables.msh")
        write(list_variable_input,_LIST_VARIABLE_CLI_GEO)
        @test main([list_variable_input,"-3","-o",list_variable_output])==
              list_variable_output
        @test read(list_variable_input,String)==_LIST_VARIABLE_CLI_GEO
        list_variable_mesh=read_mixed_msh(list_variable_output)
        @test validate(list_variable_mesh).ok
        @test mixed_crc(list_variable_mesh).sha==
              "27417f652cf93e0d6aad41c2f1b6c65af3751dfb3cb3166432d2e798f25a6493"
        @test list_variable_mesh.physical_names==Dict(
            (0,61)=>"corners",(0,65)=>"face probes",(1,62)=>"edges",
            (2,63)=>"boundary",(3,64)=>"domain")

        dynamic_tag_input=joinpath(directory,"dynamic-tags.geo")
        dynamic_tag_output=joinpath(directory,"dynamic-tags.msh")
        write(dynamic_tag_input,_DYNAMIC_TAG_CLI_GEO)
        @test main([dynamic_tag_input,"-3","-o",dynamic_tag_output])==
              dynamic_tag_output
        @test read(dynamic_tag_input,String)==_DYNAMIC_TAG_CLI_GEO
        dynamic_tag_mesh=read_mixed_msh(dynamic_tag_output)
        @test validate(dynamic_tag_mesh).ok
        @test mixed_crc(dynamic_tag_mesh).sha==
              "99aeefc2e269090b518f5896f2388bd419c8fb44d06e4743579d50d61aaedf81"
        @test dynamic_tag_mesh.physical_names==Dict(
            (0,61)=>"corners",(0,65)=>"face probes",(1,62)=>"edges",
            (2,63)=>"boundary",(3,64)=>"domain")
        @test sort([(Int(link.slave_entity),Int(link.master_entity),
                     length(link.slave_nodes))
                    for link in dynamic_tag_mesh.periodic_links
                    if link.dim==2])==[(22,24,5),(23,21,4)]

        set_max_tag_input=joinpath(directory,"set-max-tags.geo")
        set_max_tag_output=joinpath(directory,"set-max-tags.msh")
        write(set_max_tag_input,_SET_MAX_TAG_CLI_GEO)
        @test main([set_max_tag_input,"-3","-o",set_max_tag_output])==
              set_max_tag_output
        @test read(set_max_tag_input,String)==_SET_MAX_TAG_CLI_GEO
        set_max_tag_mesh=read_mixed_msh(set_max_tag_output)
        @test validate(set_max_tag_mesh).ok
        @test mixed_crc(set_max_tag_mesh).sha==
              "aa3127e5c1a302ebdb98890a4d7a62d5cb385e3cf73aa024372f809a7542a45d"
        @test set_max_tag_mesh.physical_names==Dict(
            (0,603)=>"corners",(1,604)=>"edges",
            (2,605)=>"boundary",(3,606)=>"domain")

        point_mesh_size_input=joinpath(directory,"point-mesh-sizes.geo")
        point_mesh_size_output=joinpath(directory,"point-mesh-sizes.msh")
        write(point_mesh_size_input,_POINT_MESH_SIZE_CLI_GEO)
        @test main([point_mesh_size_input,"-2","-o",point_mesh_size_output])==
              point_mesh_size_output
        @test read(point_mesh_size_input,String)==_POINT_MESH_SIZE_CLI_GEO
        point_mesh_size=read_msh(point_mesh_size_output).mesh
        @test validate(point_mesh_size).ok
        @test mesh_crc(point_mesh_size).sha==
              "b3f1bf410e917d050eacceab998b0fdf7b4cd61d1d9f263805b5120c06f1f4df"

        points_of_input=joinpath(directory,"point-mesh-size-points-of.geo")
        points_of_output=joinpath(directory,"point-mesh-size-points-of.msh")
        write(points_of_input,_POINT_MESH_SIZE_POINTS_OF_CLI_GEO)
        @test main([points_of_input,"-3","-o",points_of_output])==
              points_of_output
        @test read(points_of_input,String)==_POINT_MESH_SIZE_POINTS_OF_CLI_GEO
        points_of_mesh=read_mixed_msh(points_of_output)
        @test validate(points_of_mesh).ok
        @test mixed_crc(points_of_mesh).sha==
              "a03e62cdb5b5049f0c3ac79f829707cab1ce7575ade8a2b246808a7e222e50ca"
        @test points_of_mesh.physical_names==
              Dict((0,11)=>"vertices",(3,12)=>"domain")

        periodic_surface_volume_input=joinpath(
            directory,"periodic-surface-volume.geo")
        periodic_surface_volume_output=joinpath(
            directory,"periodic-surface-volume.msh")
        write(periodic_surface_volume_input,_PERIODIC_SURFACE_VOLUME_CLI_GEO)
        @test main([periodic_surface_volume_input,"-3","-o",
                    periodic_surface_volume_output])==
              periodic_surface_volume_output
        @test read(periodic_surface_volume_input,String)==
              _PERIODIC_SURFACE_VOLUME_CLI_GEO
        periodic_surface_volume=read_mixed_msh(
            periodic_surface_volume_output)
        @test validate(periodic_surface_volume).ok
        @test mixed_crc(periodic_surface_volume).sha==
              "27417f652cf93e0d6aad41c2f1b6c65af3751dfb3cb3166432d2e798f25a6493"
        @test length(periodic_surface_volume.periodic_links)==15
        @test sort([(Int(link.slave_entity),Int(link.master_entity),
                     length(link.slave_nodes))
                    for link in periodic_surface_volume.periodic_links
                    if link.dim==2])==[(4,6,5),(5,3,4)]
        @test periodic_surface_volume.physical_names==Dict(
            (0,61)=>"corners",(0,65)=>"face probes",(1,62)=>"edges",
            (2,63)=>"boundary",(3,64)=>"domain")

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
            "selected volume projection omits periodic relations [(1, 2, 4)]",
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
            "selected surface projection omits periodic relations [(1, 5, 6)]",
            sprint(showerror,unused_error))
        @test read(unused_input,String)==_UNUSED_PERIODIC_CLI_GEO
        @test read(unused_output,String)=="unchanged"
    end

    @test isempty(Docs.undocumented_names(Tessella.CLI;private=false))
end
