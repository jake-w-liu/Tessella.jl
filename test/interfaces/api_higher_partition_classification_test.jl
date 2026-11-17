module HigherPartitionClassificationTests
using Tessella,Test,SHA,Serialization
const E=Tessella.Elements
const A=Tessella.API

const PARTITION_SOURCE=raw"""
$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
2 1 0 0
1 0 0 0 0
2 1 0 0 0
11 0 0 0 1 0 0 1 5 2 1 -2
$EndEntities
$PartitionedEntities
2
1
21 2
2 1 0 0
21 0 1 1 1 0 0 0 0
22 0 2 1 2 1 0 0 0
21 1 11 2 1 2 0 0 0 1 0 0 1 5 2 21 -22
$EndPartitionedEntities
$Nodes
3 3 10 30
0 21 0 1
10
0 0 0
0 22 0 1
20
1 0 0
1 21 0 1
30
0.5 0 0
$EndNodes
$Elements
1 2 101 102
1 21 1 2
101 10 30
102 30 20
$EndElements
$GhostElements
1
102 2 1 1
$EndGhostElements
"""

function digest(mesh)
    output=IOBuffer();serialize(output,mesh)
    return bytes2hex(sha256(take!(output)))
end

@testset "Imported partition entities retain their complete node closure" begin
    mktempdir() do directory
        source=joinpath(directory,"partition-ascii.msh")
        write(source,PARTITION_SOURCE)
        original=E.read_mixed_msh(source)
        for binary in (false,true)
            path=joinpath(directory,"partition-$binary.msh")
            E.write_mixed_msh(path,original;version=4.1,binary)
            mesh=E.read_mixed_msh(path)
            before=digest(mesh)
            class=A._higher_classification_from_data(mesh,"partition declaration closure")
            @test class.entities==Tuple{Int,Int32}[
                (0,1),(0,2),(0,21),(0,22),(1,11),(1,21)]
            @test class.boundaries[(1,Int32(11))]==Int32[1,2]
            @test class.boundaries[(1,Int32(21))]==Int32[21,22]
            @test A._mesh_entity_nodes_with_boundary(class,1,Int32(21))==[3,1,2]
            @test class.public_tags.node_tags==UInt64[10,20,30]
            @test class.public_tags.element_tags==UInt64[101,102]
            @test class.node_entities==[(0,Int32(21)),(0,Int32(22)),(1,Int32(21))]
            @test mesh.partition_data.entities[(1,21)].boundaries==Int32[21,-22]
            @test digest(mesh)==before
        end
    end
end
end
