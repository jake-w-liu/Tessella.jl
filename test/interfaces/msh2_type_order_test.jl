using Test
using Tessella

# Regression test for Gmsh-compatible MSH 2.2 element numbering (no Gmsh needed).
const IOB_API = Tessella.API
const IOB_E = Tessella.Elements

# A MSH 4.1 ASCII file whose quadrangle surface (entity 1) has a lower tag than
# its triangle surface (entity 2). Gmsh numbers MSH 2.2 elements type by type
# (triangles before quadrangles), whatever the entity tags.
const IOB_QUAD_TRI_MSH4 = """
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

@testset "MSH 2.2 element order" begin
    mktempdir() do dir
        src = joinpath(dir, "qt.msh")
        write(src, IOB_QUAD_TRI_MSH4)

        @testset "MSH 2.2 export orders elements type by type like Gmsh" begin
            IOB_API.initialize()
            try
                IOB_API.open(src)
                IOB_API.option("Mesh.MshFileVersion", 2.2)
                IOB_API.option("Mesh.Binary", 0)
                out = joinpath(dir, "qt2.msh")
                IOB_API.write(out)
                text = read(out, String)
                lines = split(replace(text, "\r\n" => "\n"), '\n')
                first_element = findfirst(==("\$Elements"), lines)
                records = [split(l) for l in lines[first_element+2:first_element+3]]
                @test parse(Int, records[1][2]) == 2   # triangle first (tag 1)
                @test parse(Int, records[2][2]) == 3   # then the quadrangle (tag 2)
                @test records[1][1] == "1" && records[2][1] == "2"
            finally
                IOB_API.finalize()
            end
        end
    end
end
