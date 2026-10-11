@testset "MSH and STL correctness and failure atomicity" begin
    mktempdir() do directory
        source = unit_tetrahedron()
        for version in (2.2, 4.1)
            path = joinpath(directory, "tet-$version.msh")
            Tessella.IO.write_msh(path, source; version,
                                 physical_names=Dict((3, 7) => "air  region"))
            result = Tessella.IO.read_msh(path)
            @test result.mesh.coords == source.coords
            @test result.mesh.tets == source.tets && result.mesh.tet_tag == source.tet_tag
            @test result.physical_names[(3, 7)] == "air  region"
            @test validate(result.mesh).ok
        end

        path = joinpath(directory, "quad.msh")
        mixed = recombine_triangles(unit_square())
        write_mixed_msh(path, mixed)
        recovered = read_mixed_msh(path)
        @test validate(recovered).ok && recovered.coords == mixed.coords
        @test [(block.msh, block.nodes) for block in recovered.blocks] ==
              [(block.msh, block.nodes) for block in mixed.blocks]

        write(path, "previous output")
        @test_throws ArgumentError Tessella.IO.write_msh(path, source; version=3.0)
        @test read(path, String) == "previous output"
        write(path, "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n" *
                    "\$Nodes\n1\n1 0 0 0\n\$EndNodes\n" *
                    "\$Elements\n1\n1 1 0 1 99\n\$EndElements\n")
        @test_throws ArgumentError Tessella.IO.read_msh(path)
        # A failed read releases its handle, including on Windows.
        rm(path)
        @test !isfile(path)

        path = joinpath(directory, "square.stl")
        write(path, """
            solid square
            facet normal 0 0 1
              outer loop
                vertex 0 0 0
                vertex 1 0 0
                vertex 1 1 0
              endloop
            endfacet
            facet normal 0 0 1
              outer loop
                vertex 0 0 0
                vertex 1 1 0
                vertex 0 1 0
              endloop
            endfacet
            endsolid square
            """)
        surface = Tessella.IO.read_stl(path)
        @test MT.nnodes(surface) == 4 && MT.ntris(surface) == 2
        @test validate(surface).ok && surface_area(surface) ≈ 1.0
        @test_throws ArgumentError Tessella.IO.read_stl(path; max_facets=1)
    end
end
