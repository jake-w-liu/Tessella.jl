module CLINativeImportDispatchChecks
using Test,Tessella
const A=Tessella.API
const C=Tessella.CLI
const E=Tessella.Elements

@testset "CLI native CAD and mesh inputs use API import" begin
    mktempdir() do directory
        stl=joinpath(directory,"triangle.stl")
        write(stl,"solid triangle\nfacet normal 0 0 1\nouter loop\n"*
            "vertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\n"*
            "endloop\nendfacet\nendsolid triangle\n")
        vtk=joinpath(directory,"triangle.vtk")
        write(vtk,"# vtk DataFile Version 2.0\ntriangle\nASCII\n"*
            "DATASET UNSTRUCTURED_GRID\nPOINTS 3 double\n0 0 0\n1 0 0\n0 1 0\n"*
            "CELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n5\n")
        for path in (stl,vtk)
            @test C.main([path,"-parse_and_exit"])==path
            output=joinpath(directory,basename(path)*".msh")
            @test C.main([path,"-save","-o",output])==output
            mesh=E.read_mixed_msh(output)
            @test Tessella.validate(mesh).ok
            @test size(mesh.coords,2)==3
            @test sum(size(block.nodes,2) for block in mesh.blocks if
                E.msh_dimension(Int(block.msh))==2)==1
            @test_throws ArgumentError C.main([path,"-string","n=7;"])
        end
        for fixture in ("brep_native_plane.brep","brep_contact_vertex.step",
                        "brep_iges_msbo_tet.igs")
            path=normpath(joinpath(@__DIR__,"..","fixtures",fixture))
            @test C.main([path,"-parse_and_exit"])==path
        end
        # Extensionless GEO sources retain their standalone execution path.
        path=joinpath(directory,"script")
        write(path,"Point(1)={0,0,0,1};\n")
        @test C.main([path,"-parse_and_exit"])==path
    end
end
end
