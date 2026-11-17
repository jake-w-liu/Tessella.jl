module ElementFamilyGroupingTests
using Tessella,Test
const A=Tessella.API
const E=Tessella.Elements
@testset "Parent-family node queries retain interleaved polynomial order" begin
    A.initialize()
    try
        A.model.add_discrete_entity(1,7)
        coords=reshape(Float64[0,-0.,0,1,0,0,4,0,0,5,0,0,4.5,0,0,2,-0.,0,3,0,0],3,:)
        blocks=[E.ElementBlock(1,reshape(Int32[1,2],2,1)),
            E.ElementBlock(8,reshape(Int32[3,4,5],3,1)),
            E.ElementBlock(1,reshape(Int32[6,7],2,1))]
        data=E.MixedEntityData(Dict((1,7)=>E.MixedEntity(1,7,[0.,0.,0.,5.,0.,0.]));
            node_entities=fill((1,Int32(7)),7),node_parametric=fill(nothing,7),
            external_node_tags=UInt64[17,29,53,61,71,41,43],
            block_entities=[Int32[7],Int32[7],Int32[7]],
            external_element_tags=[UInt64[701],UInt64[901],UInt64[1201]])
        cache=E.MixedMesh(coords,blocks;entity_data=data)
        A._replace_mesh_cache_locked!(cache,A._higher_classification_from_data(cache,"mixed-order family fixture"))
            # Isolated Gmsh4.15.2 returns the same actual parent-family incidences
            # for requested MSH1/MSH8 and the aggregate/entity7 selectors.
            expected_nodes=UInt64[17,29,53,61,71,41,43]
            expected_xyz=Float64[0,-0.,0,1,0,0,4,0,0,5,0,0,4.5,0,0,2,-0.,0,3,0,0]
            cache=A.LAST_MESH[];before=deepcopy(cache)
            for msh in (1,8),entity in (-1,7),parametric in (false,true)
                nodes,xyz,uv=A.mesh.get_nodes_by_element_type(msh,entity,parametric)
                @test nodes==expected_nodes
                @test isequal(xyz,expected_xyz)
                @test isempty(uv)
            end
            @test A.mesh.get_elements_by_type(1,7)==(UInt64[701,1201],UInt64[17,29,41,43])
            @test A.mesh.get_elements_by_type(8,7)==(UInt64[901],UInt64[53,61,71])
            @test A.LAST_MESH[]===cache
            @test cache.entity_data.external_element_tags==before.entity_data.external_element_tags
            @test [block.nodes for block in cache.blocks]==[block.nodes for block in before.blocks]
    finally A.finalize() end
end

@testset "Family grouping preserves sparse labels, dense slots and task slices" begin
    A.initialize()
    try
        A.model.add_discrete_entity(1,7)
        coords=Float64[0 1 2 3;-0. 0. -0. 0.;0 0 0 0]
        blocks=[E.ElementBlock(1,reshape(Int32[1,2],2,1)),
            E.ElementBlock(1,reshape(Int32[3,4],2,1))]
        metadata=E.MixedEntityData(Dict((1,7)=>E.MixedEntity(1,7,[0.,0.,0.,3.,0.,0.]));
            node_entities=fill((1,Int32(7)),4),node_parametric=fill(nothing,4),
            external_node_tags=UInt64[17,29,41,43],
            block_entities=[Int32[7],Int32[7]],
            external_element_tags=[UInt64[701],UInt64[1201]])
        cache=E.MixedMesh(coords,blocks;entity_data=metadata)
        class=A._higher_classification_from_data(cache,"family grouping fixture")
        A._replace_mesh_cache_locked!(cache,class)
        before=(copy(cache.coords),[copy(block.nodes) for block in cache.blocks],
            deepcopy(cache.entity_data.external_element_tags))
        expected=(Int32[1],[UInt64[701,1201]],[UInt64[17,29,41,43]])
        @test A.mesh.get_elements(1,7)==expected
        @test A.mesh.get_elements(1)==expected
        @test A.mesh.get_elements_by_type(1,7)==(expected[2][1],expected[3][1])
        @test A.mesh.get_elements_by_type(1)==(expected[2][1],expected[3][1])
        for tasks in (2,3,4),task in 0:tasks-1
            chosen=A._mesh_task_range(2,task,tasks)
            @test A.mesh.get_elements_by_type(1,7,task,tasks)==
                (expected[2][1][chosen],vec(reshape(expected[3][1],2,2)[:,chosen]))
        end
        @test A.mesh.get_element(701)==(1,UInt64[17,29],1,7)
        @test A.mesh.get_element(1201)==(1,UInt64[41,43],1,7)
        @test A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===class
        @test isequal(cache.coords,before[1])
        @test [block.nodes for block in cache.blocks]==before[2]
        @test cache.entity_data.external_element_tags==before[3]
        result=A.mesh.get_elements(1,7)
        result[2][1][1]=999;result[3][1][1]=999
        @test A.mesh.get_elements(1,7)==expected

        # Another family occupies dense slot2 between the two matching blocks;
        # its public label is independent of the family-relative positions.
        A.model.add_discrete_entity(2,9)
        interleaved=[blocks[1],E.ElementBlock(2,reshape(Int32[1,2,3],3,1)),blocks[2]]
        metadata2=E.MixedEntityData(merge(metadata.entities,
            Dict((2,9)=>E.MixedEntity(2,9,[0.,0.,0.,3.,0.,0.])));
            node_entities=metadata.node_entities,node_parametric=metadata.node_parametric,
            external_node_tags=metadata.external_node_tags,
            block_entities=[Int32[7],Int32[9],Int32[7]],
            external_element_tags=[UInt64[701],UInt64[901],UInt64[1201]])
        cache2=E.MixedMesh(coords,interleaved;entity_data=metadata2)
        A._replace_mesh_cache_locked!(cache2,A._higher_classification_from_data(cache2,"interleaved fixture"))
        A._native_lookup_invalidate_locked!()
        @test A.mesh.get_elements_by_type(1,7)==(expected[2][1],expected[3][1])
        @test A.mesh.get_elements_by_type(1)==(expected[2][1],expected[3][1])
        @test A.mesh.get_elements_by_type(2,9)==(UInt64[901],UInt64[17,29,41])
        @test A.mesh.get_elements_by_type(1,7,1,2)==(UInt64[1201],UInt64[41,43])
        @test A.mesh.get_element(901)==(2,UInt64[17,29,41],2,9)
        A.model.add_discrete_entity(1,8)
        @test A.mesh.get_elements_by_type(1,8)==(UInt64[],UInt64[])
    finally
        A.finalize()
    end
end

if Base.JLOptions().compile_enabled!=3
@testset "Repeated family query allocation growth stays linear" begin
    A.initialize()
    allocations=Int[]
    try
        A.model.add_discrete_entity(1,7)
        for count in (1000,2000,4000)
            block=E.ElementBlock(1,repeat(reshape(Int32[1,2],2,1),1,count))
            labels=[collect(UInt64,1:count),collect(UInt64,count+1:2count)]
            data=E._native_mixed_entity_data(Dict((1,7)=>E.MixedEntity(1,7,[0.,0.,0.,1.,0.,0.]));
                node_entities=fill((1,Int32(7)),2),node_parametric=fill(nothing,2),
                external_node_tags=UInt64[17,29],block_entities=[fill(Int32(7),count),fill(Int32(7),count)],
                external_element_tags=labels)
            cache=E.MixedMesh(Float64[0 1;0 0;0 0],[block,block];entity_data=data)
            A._replace_mesh_cache_locked!(cache,A._higher_classification_from_data(cache,"repeated family resource fixture"))
            getter=()->A.mesh.get_elements_by_type(1,7)
            getter()
            push!(allocations,minimum(@allocated(getter()) for _ in 1:3))
            tags,nodes=getter()
            @test tags==vcat(labels...)
            @test nodes==repeat(UInt64[17,29],2count)
            @test A.mesh.get_elements_by_type(1,7,2count-1,2count)==(UInt64[2count],UInt64[17,29])
            @test A.LAST_MESH[]===cache
        end
        @test allocations[2]<=2.15*allocations[1]+65536
        @test allocations[3]<=2.15*allocations[2]+65536
        println("Repeated family public query allocation sweep: ",allocations)
    finally A.finalize() end
end
end
end
