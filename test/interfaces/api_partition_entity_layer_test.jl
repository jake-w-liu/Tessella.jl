module PartitionEntityLayerChecks
# `Mesh.PartitionCreateTopology=0` + `Mesh.PartitionCreateGhostCells=0` still
# materialize upstream's `createNewEntities` layer: one partition child per
# (dimension, parent entity, partition) owner, with the parent's physical
# memberships, the child as its cells' owner, and no interface topology or
# ghost metadata. Verified against the gmsh 4.15.2 oracle: a two-region Tet4
# fixture lists `(3,23)` parented to `(3,11)` with partition [1] and
# `(3,24)` parented to `(3,22)` with partition [2] after
# `partitionMesh(2)` under both options off.
using Test, Tessella
const A=Tessella.API
const E=Tessella.Elements

function box_fixture!()
    A.finalize();A.initialize();A.model.add("two tets")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    A.model.add_discrete_entity(3,11);A.model.add_discrete_entity(3,22)
    A.model.add_physical_group(3,[11];tag=101,name="left")
    A.model.add_physical_group(3,[22];tag=202,name="right")
    coordinates=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
    A.mesh.add_nodes(3,11,[11,12,13,14,15],vec(coordinates))
    A.mesh.add_elements_by_type(11,4,[501],[11,12,13,14])
    A.mesh.add_elements_by_type(22,4,[502],[11,13,12,15])
    A.mesh.generate(3)
    return nothing
end

function square_fixture!()
    A.finalize();A.initialize();A.model.add("square")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    for p in 1:4
        A.model.add_discrete_entity(0,p)
    end
    for (c,(a,b)) in enumerate([(1,2),(2,3),(3,4),(4,1)])
        A.model.add_discrete_entity(1,10+c,[a,b])
    end
    A.model.add_discrete_entity(2,20,[11,12,13,14])
    A.model.add_physical_group(2,[20];tag=100,name="plate")
    coordinates=[0. 1. 1. 0.;0. 0. 1. 1.;0. 0. 0. 0.]
    for p in 1:4
        A.mesh.add_nodes(0,p,[p],coordinates[:,p])
    end
    A.mesh.add_nodes(2,20,[5],[0.5,0.5,0.0])
    for (c,(a,b)) in enumerate([(1,2),(2,3),(3,4),(4,1)])
        A.mesh.add_elements_by_type(10+c,1,[100+c],[a,b])
    end
    for p in 1:4
        A.mesh.add_elements_by_type(p,15,[200+p],[p])
    end
    A.mesh.add_elements_by_type(20,2,[301,302],[1,2,3, 1,4,5])
    A.mesh.generate(2)
    return nothing
end

@testset "topology=0 ghosts=0 still materializes the 3-D entity layer" begin
    try
        box_fixture!()
        @test A.option("Mesh.PartitionCreateTopology")==0.0 &&
              A.option("Mesh.PartitionCreateGhostCells")==0.0
        A.mesh.partition(2,[501,502],[1,2])
        mesh=A.mesh.get()
        pdata=mesh.partition_data
        # Real partition entities — not a bare label record.
        @test pdata!==nothing && Int(pdata.num_partitions)==2
        @test length(pdata.entities)==2
        @test isempty(pdata.ghost_entities) && isempty(pdata.ghost_elements)
        @test pdata.entities[(3,23)].parent==(3,11) &&
              pdata.entities[(3,23)].partitions==Int32[1]
        @test pdata.entities[(3,24)].parent==(3,22) &&
              pdata.entities[(3,24)].partitions==Int32[2]
        # Bare ownership children: no boundary/interface wiring, physical
        # memberships transported (verified layout of gmsh 4.15.2).
        @test isempty(pdata.entities[(3,23)].boundaries) &&
              isempty(pdata.entities[(3,24)].boundaries)
        @test pdata.entities[(3,23)].physical_tags==Int32[101] &&
              pdata.entities[(3,24)].physical_tags==Int32[202]
        @test A.model.get_entities()==[(3,11),(3,22),(3,23),(3,24)]
        @test A.model.get_number_of_partitions()==2
        @test A.model.get_parent(3,23)==(3,11) && A.model.get_partitions(3,23)==[1]
        @test A.model.get_parent(3,24)==(3,22) && A.model.get_partitions(3,24)==[2]
        @test A.model.get_partitions(3,11)==[1] && A.model.get_partitions(3,22)==[2]
        @test Set(A.model.get_entities_for_physical_group(3,101))==Set([11,23])
        @test Set(A.model.get_entities_for_physical_group(3,202))==Set([22,24])
        # Cells moved to their partition children; parents own none.
        families,uids,cells=A.mesh.get_elements(3,23)
        i=findfirst(==(4),families)
        @test uids[i]==UInt64[501]
        families,uids,cells=A.mesh.get_elements(3,24)
        i=findfirst(==(4),families)
        @test uids[i]==UInt64[502]
        families,uids,cells=A.mesh.get_elements(3,11)
        @test isempty(uids) || all(isempty,uids)
        # The mesh itself round-trips the partition layer through MSH4.
        mktempdir() do directory
            path=joinpath(directory,"part.msh")
            E.write_mixed_msh(path,mesh;version=4.1)
            text=read(path,String)
            @test occursin("\$PartitionedEntities",text)
            reread=E.read_mixed_msh(path)
            @test reread.partition_data!==nothing &&
                  length(reread.partition_data.entities)==2
        end
        A.mesh.unpartition()
        @test A.model.get_entities()==[(3,11),(3,22)] &&
              A.model.get_number_of_partitions()==0
    finally
        A.finalize()
    end
end

@testset "topology=0 ghosts=0 still materializes the 2-D entity layer" begin
    try
        square_fixture!()
        @test A.option("Mesh.PartitionCreateTopology")==0.0 &&
              A.option("Mesh.PartitionCreateGhostCells")==0.0
        A.mesh.partition(2)
        mesh=A.mesh.get()
        pdata=mesh.partition_data
        @test pdata!==nothing && Int(pdata.num_partitions)==2
        @test !isempty(pdata.entities)
        @test isempty(pdata.ghost_entities) && isempty(pdata.ghost_elements)
        # One child per (dim,parent,partition) owner; every child is bare.
        for (key,entity) in pdata.entities
            @test entity.parent[1]==key[1]
            @test isempty(entity.boundaries)
            @test !isempty(entity.partitions)
            @test A.model.get_parent(key[1],key[2])==entity.parent
            @test sort(Int.(entity.partitions))==
                  sort(A.model.get_partitions(key[1],key[2]))
        end
        # Surface children appear in the public entity listing, like
        # `PartitionMesh`/upstream `createNewEntities` produce.
        surfaces=[key for key in keys(pdata.entities) if key[1]==2]
        @test !isempty(surfaces)
        for (dim,tag) in surfaces
            @test (dim,tag) in A.model.get_entities(2)
        end
        @test A.model.get_number_of_partitions()==2
        # Physical membership transported to the surface children.
        members=A.model.get_entities_for_physical_group(2,100)
        @test 20 in members &&
              any(key->key[1]==2 && key[2] in members,
                  [(k[1],k[2]) for k in surfaces])
        A.mesh.unpartition()
        @test A.model.get_entities(2)==[(2,20)] &&
              A.model.get_number_of_partitions()==0
    finally
        A.finalize()
    end
end
end
