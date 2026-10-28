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

# ── MSH 2.2 binary edge cases ────────────────────────────────────────────────
# Layout verified against live Gmsh 4.15.2: an equivalent hand-built binary
# file (same framing) was loaded by the oracle without errors. The binary
# `$Elements` section is a sequence of (type, count, numTags) block headers
# followed by full element records — types may transition freely and the same
# type may repeat in later blocks. Binary view-data sections keep an ASCII
# header and carry binary rows: Int32 tag + ncomp Float64 for NodeData/
# ElementData, plus an Int32 per-row node count for ElementNodeData.
function _io_msh22bin(path; sections=true)
    open(path, "w") do io
        print(io, "\$MeshFormat\n2.2 1 8\n")
        write(io, Int32(1))
        print(io, "\n\$EndMeshFormat\n")
        print(io, "\$PhysicalNames\n2\n1 102 \"edge\"\n2 103 \"surf\"\n",
            "\$EndPhysicalNames\n")
        print(io, "\$Nodes\n5\n")
        for (i, xyz) in enumerate(((0.0,0,0),(1,0,0),(1,1,0),(0,1,0),(0.5,0.5,0)))
            write(io, Int32(i)); for x in xyz; write(io, Float64(x)); end
        end
        print(io, "\n\$EndNodes\n\$Elements\n6\n")
        # tri block(1) → line block(1) → tri block(2) → line block(2): type
        # transitions and a repeated type header inside one binary section.
        write(io, Int32[2,1,2]); write(io, Int32[1,103,1,1,2,5])
        write(io, Int32[1,1,2]); write(io, Int32[2,102,1,1,2])
        write(io, Int32[2,2,2]); write(io, Int32[3,103,1,2,3,5])
        write(io, Int32[4,103,1,3,4,5])
        write(io, Int32[1,2,2]); write(io, Int32[5,102,2,2,3])
        write(io, Int32[6,102,2,3,4])
        print(io, "\n\$EndElements\n")
        sections || return path
        print(io, "\$NodeData\n1\n\"nodal\"\n1\n0.25\n3\n0\n1\n5\n")
        for (i, v) in enumerate((0.5,1.5,2.5,3.5,4.5))
            write(io, Int32(i)); write(io, Float64(v))
        end
        print(io, "\n\$EndNodeData\n",
            "\$ElementData\n1\n\"elscal\"\n1\n0.5\n3\n0\n1\n6\n")
        for (i, v) in enumerate((10.,20.,30.,40.,50.,60.))
            write(io, Int32(i)); write(io, Float64(v))
        end
        print(io, "\n\$EndElementData\n",
            "\$ElementNodeData\n1\n\"elnode\"\n1\n0.75\n3\n0\n3\n2\n")
        write(io, Int32(1)); write(io, Int32(3))
        for v in 1.0:9.0; write(io, v); end
        write(io, Int32(3)); write(io, Int32(3))
        for v in 10.0:18.0; write(io, v); end
        print(io, "\n\$EndElementNodeData\n")
    end
    return path
end

@testset "MSH 2.2 binary edges" begin
    mktempdir() do directory
        path = _io_msh22bin(joinpath(directory, "mixed.msh"))
        m = read_mixed_msh(path)
        @test size(m.coords, 2) == 5
        @test m.physical_names == Dict((1,102)=>"edge", (2,103)=>"surf")
        blocks = Dict(Int(b.msh) => b for b in m.blocks)
        @test sort([Int(b.msh) for b in m.blocks]) == [1, 2]
        # Both line blocks merged; connectivity keeps file order.
        @test blocks[1].nodes == Int32[1 2 3; 2 3 4]
        @test blocks[2].nodes == Int32[1 2 3; 2 3 4; 5 5 5]
        by_name = Dict(s.name => s for s in m.data_sections)
        nd = by_name["NodeData"]
        @test nd.header == [0, 1, 5] && nd.nodes == Int32[1:5;]
        @test nd.values == collect(0.5:1.0:4.5)
        ed = by_name["ElementData"]
        # `elements` holds flat internal cell indices in block order: the two
        # line blocks group first (tags 2,5,6 → cells 1-3), then the tris
        # (tags 1,3,4 → cells 4-6).
        @test ed.header == [0, 1, 6] && ed.elements == Int32[4,1,5,6,2,3]
        @test ed.values == collect(10.0:10.0:60.0)
        en = by_name["ElementNodeData"]
        @test en.header == [0, 3, 2] && en.elements == Int32[4, 5]
        @test en.row_nodes == Int32[3, 3]
        @test en.values == collect(1.0:18.0)

        # A block header promising two records but carrying one rejects
        # loudly instead of fabricating a second row.
        truncated = joinpath(directory, "trunc.msh")
        open(truncated, "w") do io
            print(io, "\$MeshFormat\n2.2 1 8\n")
            write(io, Int32(1))
            print(io, "\n\$EndMeshFormat\n\$Nodes\n2\n")
            for (i, xyz) in enumerate(((0.0,0,0),(1.0,0,0)))
                write(io, Int32(i)); for x in xyz; write(io, Float64(x)); end
            end
            print(io, "\n\$EndNodes\n\$Elements\n2\n")
            write(io, Int32[1,2,2])                 # declares 2 line elements
            write(io, Int32[1,0,0,1,2])             # only one record follows
            print(io, "\n\$EndElements\n")
        end
        @test_throws ArgumentError read_mixed_msh(truncated)
    end
end
