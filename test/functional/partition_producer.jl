@testset "owned Tet partition interface and ghost aliases" begin
    api=Tessella.API
    api.initialize()
    try
        api.option("General.Terminal",0);api.option("Mesh.Renumber",0)
        api.model.add_discrete_entity(3,11);api.model.add_discrete_entity(3,22)
        api.mesh.add_nodes(3,11,[11,12,13,14,15],
            [0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,-1.])
        api.mesh.add_elements_by_type(11,4,[501],[11,12,13,14])
        api.mesh.add_elements_by_type(22,4,[502],[11,13,12,15])
        api.mesh.generate(3)
        api.option("Mesh.PartitionCreateTopology",1);api.option("Mesh.PartitionCreateGhostCells",1)
        api.mesh.partition(2,[501,502],[1,2])
        @test api.model.get_boundary([(3,23)],false,true,false)==[(2,1)] &&
              api.model.get_boundary([(3,24)],false,true,false)==[(2,-1)]
        owned=api.mesh.get()
        @test sum(size(block.nodes,2) for block in owned.blocks if block.msh==4)==2 &&
              api.mesh.get_ghost_elements(3,25)==(UInt64[502],Int32[2]) &&
              api.mesh.get_ghost_elements(3,26)==(UInt64[501],Int32[1])
    finally
        api.finalize()
    end
end
