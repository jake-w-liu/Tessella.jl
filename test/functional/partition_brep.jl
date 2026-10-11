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

@testset "Declared lower BRep partition" begin
    A=Tessella.API
    try
        A.finalize();A.initialize();A.model.add("declared partition")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        A.model.add_discrete_entity(0,1);A.model.add_discrete_entity(0,3)
        A.model.add_discrete_entity(1,41,[1,-3]);A.model.add_discrete_entity(1,42,[3,-1])
        A.model.add_discrete_entity(2,30,[41,42]);A.model.add_discrete_entity(2,31,[-41,-42])
        A.model.add_discrete_entity(3,11,[30,31])
        coordinates=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
        A.mesh.add_nodes(0,1,[11],coordinates[:,1])
        A.mesh.add_nodes(0,3,[13],coordinates[:,3])
        A.mesh.add_nodes(1,41,[12],coordinates[:,2])
        A.mesh.add_nodes(3,11,14:15,vec(coordinates[:,4:5]))
        A.mesh.add_elements_by_type(30,2,601:603,[11,12,14,12,13,14,13,11,14])
        A.mesh.add_elements_by_type(31,2,604:606,[11,13,15,13,12,15,12,11,15])
        A.mesh.add_elements_by_type(41,1,[701,702],[11,12,12,13])
        A.mesh.add_elements_by_type(42,1,[703],[13,11])
        A.mesh.add_elements_by_type(1,15,[801],[11])
        A.mesh.add_elements_by_type(3,15,[802],[13])
        A.mesh.add_elements_by_type(11,4,[501,502],[11,12,13,14,11,13,12,15])
        A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
        A.mesh.partition(2,[601:606...,701:703...,801,802,501,502],
            [1,1,1,2,2,2,1,1,1,1,1,1,2])
        children=Dict(dim=>Dict(parent=>only(tag for (_,tag) in A.model.get_entities(dim)
            if A.model.get_parent(dim,tag)==(dim,parent)) for parent in parents)
            for (dim,parents) in ((0,[1,3]),(1,[41,42]),(2,[30,31])))
        @test A.model.get_partitions(0,children[0][1])==[1,2] &&
              A.model.get_partitions(0,children[0][3])==[1,2] &&
              A.model.get_partitions(1,children[1][41])==[1,2] &&
              A.model.get_partitions(1,children[1][42])==[1,2]
        @test Set(A.model.get_boundary([(1,children[1][41])],false,true,false))==
            Set([(0,children[0][1]),(0,-children[0][3])])
        @test Set(A.model.get_boundary([(2,children[2][30])],false,true,false))==
            Set([(1,children[1][41]),(1,children[1][42])])
        types,ids,nodes=A.mesh.get_elements(1,children[1][41])
        @test only(ids)==[701,702] && vec(only(nodes))==[11,12,12,13]
        types,ids,nodes=A.mesh.get_elements(0,children[0][1])
        @test only(ids)==[801] && vec(only(nodes))==[11]
        interface=only(tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)==(3,11))
        bounds=A.model.get_boundary([(2,interface)],false,true,false)
        @test Set(abs.(last.(bounds)))==Set([children[1][41],children[1][42]]) &&
            allequal(sign.(last.(bounds)))
        A.mesh.unpartition()
        @test A.model.get_entities()==[(0,1),(0,3),(1,41),(1,42),(2,30),(2,31),(3,11)]
        # MSH4 transport keeps signed non-partition bounds; the imported
        # declared topology partitions through the unchanged certificate.
        A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.SaveAll",1)
        path=tempname()*".msh";A.write(path)
        A.model.add("imported declared partition");A.open(path);rm(path)
        @test A.model.get_boundary([(1,41)],false,true,false)==[(0,1),(0,-3)] &&
            A.model.get_boundary([(2,30)],false,true,false)==[(1,41),(1,42)]
        A.option("Mesh.PartitionCreateTopology",1)
        A.mesh.partition(2,[601:606...,701:703...,801,802,501,502],
            [1,1,1,2,2,2,1,1,1,1,1,1,2])
        child41=only(tag for (_,tag) in A.model.get_entities(1)
            if A.model.get_parent(1,tag)==(1,41))
        @test A.model.get_number_of_partitions()==2 &&
            A.model.get_partitions(1,child41)==[1,2]
        A.mesh.unpartition()
        @test A.model.get_boundary([(1,41)],false,true,false)==[(0,1),(0,-3)]
    finally
        A.finalize()
    end
end
