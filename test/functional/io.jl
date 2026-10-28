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

@testset "MSH 4.1 external tags beyond the Int32 range" begin
    # MSH4 node and element tags are `size_t`: the external tag space is the
    # full unsigned 64-bit range, while entity/physical tags stay `int`
    # (upstream `addDiscreteEntity` auto-assigns tags beyond Int32). Verified
    # against live Gmsh 4.15.2: node tag 3000000005 round-trips verbatim;
    # element tags are truncated to 32 bits by upstream's own writer.
    mktempdir() do directory
        node_tags = UInt64[2_147_483_748, 2_147_483_749, 3_000_000_005]
        element_tags = UInt64[4_294_967_297]  # beyond typemax(UInt32)
        mesh = MixedMesh(
            Float64[0 1 0; 0 0 1; 0 0 0],
            Any[ElementBlock(2, reshape(Int32[1, 2, 3], 3, 1), Int32[3])];
            entity_data=MixedEntityData(
                Dict{Tuple{Int,Int},MixedEntity}((2, 7) => MixedEntity(
                    2, 7, (0.0, 0.0, 0.0, 1.0, 1.0, 0.0);
                    physical_tags=Int32[3]));
                node_entities=[(2, Int32(7)), (2, Int32(7)), (2, Int32(7))],
                node_parametric=[nothing, nothing, nothing],
                external_node_tags=node_tags,
                block_entities=[Int32[7]],
                external_element_tags=[element_tags]))
        @test validate(mesh).ok
        crc = mixed_crc(mesh).sha
        for binary in (false, true)
            path = joinpath(directory, "wide-tags-$binary.msh")
            write_mixed_msh(path, mesh; version=4.1, binary)
            back = read_mixed_msh(path)
            @test validate(back).ok
            @test back.entity_data.external_node_tags == node_tags
            @test back.entity_data.external_element_tags == [element_tags]
            @test mixed_crc(back).sha == crc
            binary ||
                @test occursin("3000000005", read(path, String))
        end
    end
end
