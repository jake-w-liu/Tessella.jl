module HigherPartitionClassificationTests
using Tessella,Test,SHA,Serialization
const E=Tessella.Elements
const A=Tessella.API

function digest(mesh)
    output=IOBuffer();serialize(output,mesh)
    return bytes2hex(sha256(take!(output)))
end

@testset "Imported partition entities retain their complete node closure" begin
    # These original native MSH packets also back the full import regression.
    data=joinpath(@__DIR__,"..","data","native_msh_import_uint64_v1")
    for name in ("partition_ascii4.msh","partition_binary4.msh")
        mesh=E.read_mixed_msh(joinpath(data,name))
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
