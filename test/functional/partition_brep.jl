@testset "Closed lower BRep partition" begin
    A=Tessella.API
    try
        A.finalize();A.initialize();A.model.add("closed partition")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        A.model.add_discrete_entity(2,30);A.model.add_discrete_entity(3,11,[30])
        coordinates=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
        A.mesh.add_nodes(3,11,11:15,vec(coordinates))
        A.mesh.add_elements_by_type(30,2,601:606,
            [11,12,14,12,13,14,13,11,14,11,13,15,13,12,15,12,11,15])
        A.mesh.add_elements_by_type(11,4,[501,502],[11,12,13,14,11,13,12,15])
        A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
        A.mesh.partition(2,[601,602,603,604,605,606,501,502],[1,1,1,2,2,2,1,2])
        @test Set(A.model.get_boundary([(3,12)],false,true,false))==Set([(2,31),(2,33)]) &&
              Set(A.model.get_boundary([(3,13)],false,true,false))==Set([(2,32),(2,-33)])
        types,ids,_=A.mesh.get_elements(1,1)
        @test types==[1] && length(only(ids))==3 &&
              isempty(A.model.get_boundary([(1,1)],false,true,false)) &&
              A.model.get_parent(1,1)==(2,30) && A.model.get_partitions(1,1)==[1,2]
    finally
        A.finalize()
    end
end
