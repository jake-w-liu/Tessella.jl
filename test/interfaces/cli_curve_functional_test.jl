module CLICurveFunctionalChecks
using Test, Tessella
const E=Tessella.Elements

const CURVES="""
    Point(1) = {0, 0, 0, 1};
    Point(2) = {1, 0, 0, 1};
    Point(3) = {0, 1, 0, 1};
    Point(4) = {1, 1, 0, 1};
    Point(5) = {3, 2, 0, 1};
    Line(1) = {1, 2};
    Line(2) = {3, 4};
    Transfinite Curve {1, 2} = 3;
    Periodic Curve {2} = {1} Translate {0, 1, 0};
    Physical Curve("pair", 11) = {1, 2};
    Physical Point("orphan", 12) = {5};
    Mesh.SaveAll = 1;
    """

@testset "classified one-dimensional CLI" begin
    mktempdir() do directory
        source=joinpath(directory,"curves.geo")
        output=joinpath(directory,"curves.msh")
        write(source,CURVES)
        execution=Tessella.execute_geo(source;mesh_dim=1)
        @test execution.mesh!==nothing
        @test Tessella.MeshTypes.nsegs(execution.mesh)==4
        @test all(part->part[1] in (0,1),execution.actual_mesh_parts)
        @test Tessella.CLI.main([source,"-1","-o",output])==output
        mesh=E.read_mixed_msh(output)
        @test Tessella.validate(mesh).ok
        @test Set(Int(block.msh) for block in mesh.blocks)==Set([15,1])
        @test sum(size(block.nodes,2) for block in mesh.blocks if block.msh==1)==4
        @test mesh.physical_names[(1,11)]=="pair" && mesh.physical_names[(0,12)]=="orphan"
        @test Set(mesh.entity_data.node_entities)==Set([(0,Int32(i)) for i in 1:5]) ∪
            Set([(1,Int32(1)),(1,Int32(2))])
        link=only(filter(link->link.dim==1,mesh.periodic_links))
        @test (link.slave_entity,link.master_entity)==(2,1)
        @test length(link.slave_nodes)==3
        @test all(mesh.coords[:,s]≈mesh.coords[:,m]+[0.0,1.0,0.0]
                  for (s,m) in zip(link.slave_nodes,link.master_nodes))
        @test read(source,String)==CURVES

        # The same producer supports configured high-order Curve elements.
        quadratic=joinpath(directory,"quadratic.geo")
        write(quadratic,CURVES * "Mesh.ElementOrder = 2;\n")
        high=Tessella.CLI.main([quadratic,"-1"])
        elevated=E.read_mixed_msh(high)
        @test Tessella.validate(elevated).ok
        @test Set(Int(block.msh) for block in elevated.blocks)==Set([15,8])
        @test sum(size(block.nodes,2) for block in elevated.blocks if block.msh==8)==4
        highlink=only(filter(link->link.dim==1,elevated.periodic_links))
        @test length(highlink.slave_nodes)==5

        for format in ("msh2","msh4")
            binary_path=joinpath(directory,format*".msh")
            @test Tessella.CLI.main([source,"-1","-format",format,"-bin","-o",binary_path])==binary_path
            recovered=E.read_mixed_msh(binary_path)
            @test Tessella.validate(recovered).ok
            @test sum(size(block.nodes,2) for block in recovered.blocks if block.msh==1)==4
            @test recovered.physical_names==mesh.physical_names
            @test length(only(filter(link->link.dim==1,recovered.periodic_links)).slave_nodes)==3
            open(binary_path) do io
                @test readline(io)=="\$MeshFormat"
                @test split(readline(io))[1:2]==[format=="msh2" ? "2.2" : "4.1","1"]
            end
        end

        triangle_source=joinpath(directory,"triangle.geo")
        write(triangle_source,"""
            Point(1) = {0, 0, 0, 1}; Point(2) = {1, 0, 0, 1};
            Point(3) = {0, 1, 0, 1};
            Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 1};
            Curve Loop(1) = {1, 2, 3}; Plane Surface(1) = {1};
            Transfinite Curve {1, 2, 3} = 2; Transfinite Surface {1};
            """)
        @test Tessella.execute_geo(triangle_source;mesh_dim=2).mesh isa Tessella.Mesh
        triangle_output=joinpath(directory,"triangle.msh")
        Tessella.CLI.main([triangle_source,"-2","-bin","-o",triangle_output])
        triangle=E.read_mixed_msh(triangle_output)
        @test Tessella.validate(triangle).ok
        @test sum(size(block.nodes,2) for block in triangle.blocks if block.msh==2)==1

        before=read(output)
        @test_throws ArgumentError Tessella.CLI.main([source,"-1","-2","-o",output])
        @test_throws ArgumentError Tessella.CLI.main([source,"-1","-o",source])
        @test_throws ArgumentError Tessella.CLI.main([source,"-1","-format","bad","-o",output])
        @test_throws ArgumentError Tessella.CLI.main([source,"-1","-bin","-bin","-o",output])
        @test_throws ArgumentError Tessella.CLI.main([source,"-format","msh2"])
        @test read(output)==before && read(source,String)==CURVES
        empty_source=joinpath(directory,"empty.geo")
        write(empty_source,"// no geometry\n")
        @test_throws ArgumentError Tessella.CLI.main([empty_source,"-1","-o",output])
        @test read(output)==before
    end
end
end
