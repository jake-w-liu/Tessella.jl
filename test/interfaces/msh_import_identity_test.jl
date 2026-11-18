using Test
using Tessella

# Regression tests for MSH import scaling and binary-read robustness (no Gmsh needed).
const IOA_API = Tessella.API
const IOA_E = Tessella.Elements

# A MSH 4.1 ASCII file whose quadrangle surface (entity 1) has a lower tag than
# its triangle surface (entity 2). Gmsh numbers MSH 2.2 elements type by type
# (triangles before quadrangles), whatever the entity tags.
const IOA_QUAD_TRI_MSH4 = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
0 0 2 0
1 0 0 0 1 1 0 0 0
2 1 0 0 2 1 0 0 0
\$EndEntities
\$Nodes
2 5 1 5
2 1 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
2 2 0 1
5
2 0 0
\$EndNodes
\$Elements
2 2 1 2
2 1 3 1
1 1 2 3 4
2 2 2 1
2 2 3 5
\$EndElements
"""

@testset "MSH import identity and binary reads" begin
    mktempdir() do dir
        src = joinpath(dir, "qt.msh")
        write(src, IOA_QUAD_TRI_MSH4)

        @testset "identity receipts are vectors, not giant tuples" begin
            _, identities = IOA_E._read_mixed_msh_with_identities(src)
            # `Tuple(vector)` made the receipt's type grow with the element
            # count and overflowed the stack for meshes of ~5e5 cells.
            @test identities.cell_events isa Vector
            @test length(identities.cell_events) == 2
        end

        @testset "truncated binary scalars still raise a clear error" begin
            bin = joinpath(dir, "bin.msh")
            m = IOA_E.read_mixed_msh(src)
            IOA_E.write_mixed_msh(bin, m; version = 4.1, binary = true)
            bytes = read(bin)
            for cut in (length(bytes) - 40, length(bytes) - 7, length(bytes) ÷ 2)
                cutfile = joinpath(dir, "cut_$cut.msh")
                write(cutfile, bytes[1:cut])
                err = try
                    IOA_E.read_mixed_msh(cutfile); nothing
                catch e
                    e
                end
                @test err isa ArgumentError
            end
        end

        @testset "sparse and huge node tags keep reading correctly" begin
            # A tag far above 4x the node count bypasses the dense lookup vector.
            huge = "1000000000000"
            sparse = replace(IOA_QUAD_TRI_MSH4,
                "2 5 1 5
" => "2 5 1 $huge
",
                "5
2 0 0" => "$huge
2 0 0",
                "2 2 3 5" => "2 2 3 $huge")
            file = joinpath(dir, "sparse.msh")
            write(file, sparse)
            m = IOA_E.read_mixed_msh(file)
            tri = only(b for b in m.blocks if b.msh == 2)
            quad = only(b for b in m.blocks if b.msh == 3)
            @test vec(tri.nodes) == Int32[2, 3, 5]
            @test vec(quad.nodes) == Int32[1, 2, 3, 4]
            @test m.entity_data.external_node_tags[5] == UInt64(10)^12
        end
    end
end
