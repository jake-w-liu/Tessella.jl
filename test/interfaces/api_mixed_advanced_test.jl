using Test
using Tessella
using Tessella.Elements: ElementBlock, MixedMesh, lagrange_nodes

function _mixed_advanced_fixture(owners=Int32[1,1])
    api=Tessella.API
    api.finalize();api.initialize()
    for owner in unique(owners)
        api.model.add_discrete_entity(3,owner)
    end
    coordinates=hcat(lagrange_nodes(5),lagrange_nodes(5))
    cells=hcat(Int32.(1:8),Int32.(9:16))
    mesh=MixedMesh(coordinates,[ElementBlock(5,cells,Int32[7,8])])
    entities=[(3,owner) for owner in unique(owners)]
    nodes=vcat(fill((3,owners[1]),8),fill((3,owners[2]),8))
    class=api._mixed_classification(mesh,first(entities),entities,nodes,
        Dict{Tuple{Int,Int32},Vector{Int32}}(),Dict(5=>copy(owners)))
    lock(api.STATE_LOCK) do
        api._replace_mesh_cache_locked!(mesh,class)
    end
    return api
end

function _mixed_isolated_partition_fixture(count)
    coords=zeros(3,2count)
    coords[1,:]=collect(1:2count)
    return MixedMesh(coords,[ElementBlock(1,reshape(Int32.(1:2count),2,count))])
end

@testset "native mixed duplicate removal and partitioning" begin
    api=Tessella.API
    try
        _mixed_advanced_fixture()
        api.mesh.remove_duplicate_nodes()
        mesh=api.mesh.get()
        @test mesh isa MixedMesh
        @test size(mesh.coords,2)==8
        @test only(mesh.blocks).msh==5
        @test only(mesh.blocks).nodes[:,1]==only(mesh.blocks).nodes[:,2]
        @test api.mesh.get_element(2)[3:4]==(3,1)
        api.mesh.remove_duplicate_elements()
        @test size(only(api.mesh.get().blocks).nodes,2)==1
        @test only(api.mesh.get().blocks).tags==Int32[7]
        @test api.mesh.get_element(1)[1]==5

        _mixed_advanced_fixture(Int32[1,2])
        api.mesh.remove_duplicate_nodes([(3,1)])
        @test size(api.mesh.get().coords,2)==16
        api.mesh.remove_duplicate_nodes()
        @test size(api.mesh.get().coords,2)==8
        api.mesh.remove_duplicate_elements()
        @test size(only(api.mesh.get().blocks).nodes,2)==2
        @test api.mesh.get_element(1)[4]==1
        @test api.mesh.get_element(2)[4]==2
        before=api.mesh.get()
        @test_throws ArgumentError api.mesh.remove_duplicate_nodes([(3,99)])
        @test api.mesh.get().coords==before.coords
        @test only(api.mesh.get().blocks).nodes==only(before.blocks).nodes

        _mixed_advanced_fixture()
        api.mesh.partition(2)
        record=api.LAST_MESH_PARTITION[]
        @test record.element_partitions==Int32[1,2]
        @test record.mesh===api.LAST_MESH[]
        api.mesh.unpartition()
        @test api.LAST_MESH_PARTITION[]===nothing
        api.mesh.partition(0,[1,2],[2,1])
        @test api.LAST_MESH_PARTITION[].num_partitions==2
        @test api.LAST_MESH_PARTITION[].element_partitions==Int32[2,1]
        @test_throws ArgumentError api.mesh.partition(2,[1],[1])
        @test api.LAST_MESH_PARTITION[].element_partitions==Int32[2,1]
        api.mesh.remove_duplicate_nodes()
        @test api.LAST_MESH_PARTITION[]===nothing
        small=_mixed_isolated_partition_fixture(2000)
        large=_mixed_isolated_partition_fixture(4000)
        @test api._partition_elements(small,1)==ones(Int32,2000)
        @test api._partition_elements(large,3)==vcat(
            fill(Int32(1),1334),fill(Int32(2),1334),fill(Int32(3),1332))
        api._partition_elements(large,1)
        a=minimum(@allocated(api._partition_elements(small,1)) for _ in 1:3)
        b=minimum(@allocated(api._partition_elements(large,1)) for _ in 1:3)
        @test b<=2a+262144
    finally
        api.finalize()
    end
end
