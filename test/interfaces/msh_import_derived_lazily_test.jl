module MshImportDerivedLazilyTest

using Test, Tessella

# Two surfaces (a triangle and a quadrangle) and a curve, one physical group.
const MSH = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$PhysicalNames
1
2 5 "plate"
\$EndPhysicalNames
\$Entities
0 1 2 0
1 0 0 0 1 0 0 0 0
1 0 0 0 2 1 0 1 5 0
2 0 0 0 2 1 0 1 5 0
\$EndEntities
\$Nodes
2 5 1 5
1 1 0 2
1
2
0 0 0
1 0 0
2 1 0 3
3
4
5
1 1 0
0 1 0
2 1 0
\$EndNodes
\$Elements
3 3 1 3
1 1 1 1
1 1 2
2 1 2 1
2 1 2 5
2 2 3 1
3 2 5 3 4
\$EndElements
"""

@testset "Import derives its records and parts from the arrays" begin
    mktempdir() do dir
        path = joinpath(dir, "m.msh")
        write(path, MSH)
        GE = Tessella.GeoExec
        function discrete(fast)
            GE._IMPORT_FAST_MERGE[] = fast
            try
                model = GE._geo_stage_mesh_file(path).model
                return (model.discrete, model.physical, model.entity_physicals, model.physical_names)
            finally
                GE._IMPORT_FAST_MERGE[] = true
            end
        end
        # The counting-sort merge builds the same records as the general function.
        @test discrete(true) == discrete(false)
        @test length(discrete(true)[1]) == 3

        API = Tessella.API
        API.initialize()
        try
            API.open(path)
            @test sum(length, API.mesh.get_elements(2)[2]) == 2
            @test API.model.get_physical_name(2, 5) == "plate"
            out = joinpath(dir, "o.msh")
            API.option("Mesh.MshFileVersion", 4.1)
            API.write(out)
            API.initialize()
            API.open(out)
            @test sum(length, API.mesh.get_elements(2)[2]) == 2
        finally
            API.finalize()
        end
    end
end

@testset "Writing packs binary runs and validates fixed blocks without a key per cell" begin
    E = Tessella.Elements
    coords = [0.0 1 0 0 1; 0 0 1 0 1; 0 0 0 1 0]
    nodes = Int32[1 1 2; 2 2 3; 3 3 4; 4 4 5]
    # The second tetrahedron repeats the first.
    mesh = E.MixedMesh(coords, [E.ElementBlock(4, nodes, Int32[1, 1, 1])])
    messages = E.validate(mesh).messages
    @test any(occursin("duplicate type 4 cell at block 1 cell 2", m) for m in messages)
    bad = E.MixedMesh(coords, [E.ElementBlock(4, Int32[1 2; 1 3; 2 4; 3 5], Int32[0, 0])])
    @test any(occursin("repeats a node index", m) for m in E.validate(bad).messages)
    mktempdir() do dir
        a = joinpath(dir, "a.msh"); b = joinpath(dir, "b.msh")
        ok = E.MixedMesh(coords, [E.ElementBlock(4, Int32[1 1; 2 2; 3 3; 4 5], Int32[1, 2])])
        E.write_mixed_msh(a, ok; version = 4.1, binary = true)
        E.write_mixed_msh(b, E.read_mixed_msh(a); version = 4.1, binary = true)
        @test read(a) == read(b)
    end
end

end # module
