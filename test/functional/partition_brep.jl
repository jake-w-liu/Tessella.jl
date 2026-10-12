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

@testset "Mixed linear family BRep partition" begin
    A=Tessella.API
    try
        # Tet4 + Prism6 share a Tri3 face inside one declared Region; the
        # partition interface is a Tri3 cell between Quad4 lateral surfaces.
        A.finalize();A.initialize();A.model.add("mixed partition")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        for tag in 31:37;A.model.add_discrete_entity(2,tag);end
        A.model.add_discrete_entity(3,11,[31,32,33,34,35,36,37])
        coordinates=[0. 1. 0. 0. 2. 1. 1.;0. 0. 1. 0. 0. 0. 1.;0. 0. 0. 1. 0. 1. 0.]
        A.mesh.add_nodes(3,11,1:7,vec(coordinates))
        A.mesh.add_elements_by_type(31,2,[601],[1,3,2])
        A.mesh.add_elements_by_type(32,2,[602],[1,2,4])
        A.mesh.add_elements_by_type(33,2,[603],[1,4,3])
        A.mesh.add_elements_by_type(34,2,[604],[5,7,6])
        A.mesh.add_elements_by_type(35,3,[605],[5,6,4,2])
        A.mesh.add_elements_by_type(36,3,[606],[5,2,3,7])
        A.mesh.add_elements_by_type(37,3,[607],[6,7,3,4])
        A.mesh.add_elements_by_type(11,4,[501],[1,2,3,4])
        A.mesh.add_elements_by_type(11,6,[502],[5,6,7,2,4,3])
        A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
        A.mesh.partition(2,[601:607...,501,502],[1,1,1,2,2,2,2,1,2])
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        types,ids,nodes=A.mesh.get_elements(2,interface)
        @test A.model.get_number_of_partitions()==2 &&
            A.model.get_partitions(2,interface)==[1,2]
        @test types==Int32[2] && only(ids)==[608] && vec(only(nodes))==[4,2,3]
        @test length(A.model.get_entities(1))==3 &&
            all(A.model.get_partitions(1,seam)==[1,2] for seam in 1:3)
        A.mesh.unpartition()
        @test A.model.get_entities()==[(2,31),(2,32),(2,33),(2,34),(2,35),(2,36),(2,37),(3,11)]
    finally
        A.finalize()
    end
end

@testset "Second-order lower BRep partition" begin
    A=Tessella.API
    try
        # Two Tet10 share the {1,2,3} Tri6; edge-node rows are shared between
        # volume cells and their surface carriers. Tags 6:14 are edge nodes.
        A.finalize();A.initialize();A.model.add("p2 partition")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        A.model.add_discrete_entity(2,30);A.model.add_discrete_entity(3,11,[30])
        coordinates=[0. 1. 0. 0. 0. .5 .5 0. 0. 0. .5 .5 0. 0.;
                     0. 0. 1. 0. 0. 0. .5 .5 0. .5 0. 0. .5 0.;
                     0. 0. 0. 1. -1. 0. 0. 0. .5 .5 .5 -.5 -.5 -.5]
        A.mesh.add_nodes(3,11,1:14,vec(coordinates))
        A.mesh.add_elements_by_type(30,9,601:606,
            [1,2,4,6,11,9, 2,3,4,7,10,11, 3,1,4,8,9,10,
             1,3,5,8,13,14, 3,2,5,7,12,13, 2,1,5,6,14,12])
        A.mesh.add_elements_by_type(11,11,[501,502],
            [1,2,3,4,6,7,8,9,10,11, 1,3,2,5,8,7,6,14,12,13])
        A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
        A.mesh.partition(2,[601:606...,501,502],[1,1,1,2,2,2,1,2])
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        types,ids,nodes=A.mesh.get_elements(2,interface)
        @test types==[9] && length(only(ids))==1 &&
            sort(vec(only(nodes)))==[1,2,3,6,7,8]
        seam=only(tag for (dim,tag) in A.model.get_entities()
            if dim==1 && A.model.get_parent(1,tag)==(2,30))
        types,ids,nodes=A.mesh.get_elements(1,seam)
        @test types==[8] && length(only(ids))==3 &&
            A.model.get_partitions(1,seam)==[1,2]
        A.mesh.unpartition()
        @test A.model.get_entities()==[(2,30),(3,11)]
    finally
        A.finalize()
    end
end

@testset "Junction BRep partition" begin
    A=Tessella.API
    try
        # Three volumes meet at rim edge {1,2}: two pairwise interfaces plus
        # exterior carriers close the junction fan; the shared seam and the
        # dim-0 junction entities carry all three partitions.
        A.finalize();A.initialize();A.model.add("junction partition")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        A.model.add_discrete_entity(2,30);A.model.add_discrete_entity(3,11,[30])
        coordinates=[0. 0. 1. 0. -1. -0.5;0. 0. 0. 1. -1. -1.;0. 1. 0. 0. 0.2 0.2]
        A.mesh.add_nodes(3,11,1:6,vec(coordinates))
        A.mesh.add_elements_by_type(30,2,601:608,
            [1,4,3, 4,2,3, 1,5,2, 1,3,5, 5,3,2, 1,6,4, 1,2,6, 6,2,4])
        A.mesh.add_elements_by_type(11,4,[501,502,503],
            [1,2,3,4, 1,3,2,5, 1,2,4,6])
        A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
        A.mesh.partition(3,[601:608...,501,502,503],[1,1,2,2,2,3,3,3,1,2,3])
        jpoints=[tag for (_,tag) in A.model.get_entities(0)]
        jseam=only(tag for (_,tag) in A.model.get_entities(1)
            if A.model.get_partitions(1,tag)==[1,2,3])
        @test length(jpoints)==2 && A.model.get_parent(1,jseam)==(2,30) &&
            all(A.model.get_partitions(0,tag)==[1,2,3] for tag in jpoints)
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test A.model.get_entities()==[(2,30),(3,11)]
    finally
        A.finalize()
    end
end

@testset "Third-order lower BRep partition" begin
    A=Tessella.API
    try
        # Two Tet20 share the {1,2,3} Tri10; edge nodes sit at i/3 fractions
        # and face nodes at centroids, so volume cells and surface carriers
        # share literal node rows. Tags 6:17 are t1 edge nodes, 18:21 its face
        # nodes; 22:27 are t2-only edge nodes and 28:30 its face nodes.
        A.finalize();A.initialize();A.model.add("p3 partition")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        A.model.add_discrete_entity(2,30);A.model.add_discrete_entity(3,11,[30])
        t=1/3
        coordinates=[0. 1 0 0 0  t 2t 2t t 0 0  0 0 0 0 t 2t  t t 0 t  0 0 t 2t 0 0  t 0 t;
                     0. 0 1 0 0  0 0 t 2t 2t t  0 0 t 2t 0 0  t 0 t t  0 0 0 0 t 2t  0 t t;
                     0. 0 0 1 -1 0 0 0 0 0 0  2t t 2t t 2t t  0 t t t  -2t -t -2t -t -2t -t  -t -t -t]
        A.mesh.add_nodes(3,11,1:30,vec(coordinates))
        A.mesh.add_elements_by_type(30,21,601:606,
            [1,2,4,6,7,17,16,12,13,19, 1,4,3,13,12,14,15,10,11,20,
             4,2,3,16,17,8,9,15,14,21, 1,3,5,11,10,27,26,22,23,29,
             1,5,2,23,22,24,25,7,6,28, 5,3,2,26,27,9,8,25,24,30])
        A.mesh.add_elements_by_type(11,29,[501,502],
            [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
             1,3,2,5,11,10,9,8,7,6,22,23,24,25,26,27,18,29,28,30])
        A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
        A.mesh.partition(2,[601:606...,501,502],[1,1,1,2,2,2,1,2])
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        types,ids,nodes=A.mesh.get_elements(2,interface)
        @test types==[21] && length(only(ids))==1 &&
            sort(vec(only(nodes)))==[1,2,3,6,7,8,9,10,11,18]
        seam=only(tag for (dim,tag) in A.model.get_entities()
            if dim==1 && A.model.get_parent(1,tag)==(2,30))
        types,ids,nodes=A.mesh.get_elements(1,seam)
        @test types==[26] && length(only(ids))==3 &&
            A.model.get_partitions(1,seam)==[1,2]
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test A.model.get_entities()==[(2,30),(3,11)]
    finally
        A.finalize()
    end
end

@testset "Coincident Surface partition" begin
    A=Tessella.API
    try
        # Two Regions share facet {1,2,4} carried by declared Surface 30: the
        # partition child keeps the single cell, merges [1,2] membership and
        # binds both Region children — no duplicate interface entity or cell.
        A.finalize();A.initialize();A.model.add("coincident partition")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        for tag in 30:36;A.model.add_discrete_entity(2,tag);end
        A.model.add_discrete_entity(3,11,[30,31,32,33])
        A.model.add_discrete_entity(3,12,[-30,34,35,36])
        coordinates=[0. 1. 0. 0. .5;0. 0. 1. 0. -.5;0. 0. 0. 1. .2]
        A.mesh.add_nodes(3,11,1:5,vec(coordinates))
        A.mesh.add_elements_by_type(30,2,[601],[1,2,4])
        A.mesh.add_elements_by_type(31,2,[602],[1,3,2])
        A.mesh.add_elements_by_type(32,2,[603],[1,4,3])
        A.mesh.add_elements_by_type(33,2,[604],[4,2,3])
        A.mesh.add_elements_by_type(34,2,[605],[1,2,5])
        A.mesh.add_elements_by_type(35,2,[606],[1,5,4])
        A.mesh.add_elements_by_type(36,2,[607],[5,2,4])
        A.mesh.add_elements_by_type(11,4,[501],[1,2,3,4])
        A.mesh.add_elements_by_type(12,4,[502],[1,2,4,5])
        A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
        A.mesh.partition(2,[601:607...,501,502],[1,1,1,1,2,2,2,1,2])
        shared=only(tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)==(2,30))
        region1=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11))
        region2=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,12))
        bounds1=Set(A.model.get_boundary([(3,region1)],false,true,false))
        bounds2=Set(A.model.get_boundary([(3,region2)],false,true,false))
        @test A.model.get_partitions(2,shared)==[1,2] &&
            (2,shared) in bounds1 && (2,-shared) in bounds2 &&
            isempty(A.model.get_entities(1))
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test A.model.get_entities()==[(2,s) for s in 30:36] ∪ [(3,11),(3,12)]
    finally
        A.finalize()
    end
end

@testset "Periodic and parametric partition transport" begin
    A=Tessella.API
    try
        # Slave Surface 32 is the translated master Surface 31; every boundary
        # vertex carries a Surface UV parameter. Partition children keep the
        # link bound to the original entities (Gmsh convention) and parametric
        # rows reclassify to same-parent children.
        A.finalize();A.initialize();A.model.add("periodic prism")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        for tag in 31:35;A.model.add_discrete_entity(2,tag);end
        A.model.add_discrete_entity(3,11,[31,32,33,34,35])
        A.mesh.add_nodes(2,31,[11,12,13],[0.,0.,0., 1.,0.,0., 0.,1.,0.],
            [0.,0., 1.,0., 0.,1.])
        A.mesh.add_nodes(2,32,[14,15,16],[0.,0.,1., 1.,0.,1., 0.,1.,1.],
            [0.,0., 1.,0., 0.,1.])
        A.mesh.add_elements_by_type(31,2,[631],[11,13,12])
        A.mesh.add_elements_by_type(32,2,[632],[14,15,16])
        A.mesh.add_elements_by_type(33,2,[633,634],[11,12,15, 11,15,14])
        A.mesh.add_elements_by_type(34,2,[635,636],[13,11,16, 11,14,16])
        A.mesh.add_elements_by_type(35,2,[637,638],[12,13,16, 12,16,15])
        A.mesh.add_elements_by_type(11,4,[501,502,503],
            [11,12,13,16, 11,12,16,15, 11,14,15,16])
        A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
        path=joinpath(mktempdir(),"prism.msh");A.write(path)
        lines=readlines(path)
        append!(lines,["\$Periodic","1","2 32 31","16",
            "1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 1.0",
            "3","14 11","15 12","16 13","\$EndPeriodic"])
        write(path,join(lines,"\n")*"\n")
        A.model.add("imported periodic prism");A.open(path);rm(path)
        A.mesh.partition(2,[631:638...,501,502,503],
            Int32[1,1,2,1,1,1,1,2,1,2,1])
        link=only(A._cached_mesh_locked("test").periodic_links)
        @test (link.dim,link.slave_entity,link.master_entity)==(2,32,31) &&
            length(link.slave_nodes)==3
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test length(A._cached_mesh_locked("test").periodic_links)==1 &&
            count(!isnothing,A._cached_mesh_locked("test").entity_data.node_parametric)==6
    finally
        A.finalize()
    end
end
