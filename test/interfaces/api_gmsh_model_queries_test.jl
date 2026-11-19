# Model queries on meshes loaded from MSH files answer what Gmsh 4.15.2 answers.
# The test/fixtures/gmsh_4152_api_*.msh files were written by Gmsh 4.15.2 and
# every expected value below is what `gmsh` returned after reloading the file.
module GmshModelQueryParityChecks
using Test, Tessella
const A=Tessella.API
const DATA=joinpath(@__DIR__,"..","fixtures")

function load(name)
    A.finalize();A.initialize();A.option("General.Terminal",0)
    A.open(joinpath(DATA,"gmsh_4152_api_"*name*".msh"))
    return nothing
end

@testset "get_boundary follows GModel::getBoundaryTags" begin
    try
        load("square_o1")
        # A Curve lists begin then end, never signed; a negative tag swaps them.
        @test A.model.get_boundary([(1,1)],false,true,false)==[(0,1),(0,2)]
        @test A.model.get_boundary([(1,-1)],false,true,false)==[(0,2),(0,1)]
        # A Surface tag is taken by magnitude.
        @test A.model.get_boundary([(2,1)],false,true,false)==[(1,1),(1,2),(1,3),(1,4)]
        @test A.model.get_boundary([(2,-1)],false,true,false)==[(1,1),(1,2),(1,3),(1,4)]
        # Recursive: tag order for a Surface, begin/end for a Curve, and an
        # unknown Point is not an error.
        @test A.model.get_boundary([(2,1)],false,false,true)==[(0,1),(0,2),(0,3),(0,4)]
        @test A.model.get_boundary([(1,3)],false,false,true)==[(0,3),(0,4)]
        @test isempty(A.model.get_boundary([(0,9)],false,false,true))
        @test A.model.get_adjacencies(1,1)==([1],[1,2])
        load("cube_o1")
        @test A.model.get_boundary([(3,1)],false,false,true)==[(0,tag) for tag in 1:8]
        # A closed Curve lists its Point twice; the pair cancels when combined.
        load("disk_o2")
        @test A.model.get_boundary([(1,1)],false,true,false)==[(0,1),(0,1)]
        @test A.model.get_boundary([(1,1)],false,false,true)==[(0,1),(0,1)]
        @test isempty(A.model.get_boundary([(1,1)],true,false,false))
        @test A.model.get_adjacencies(1,1)==([1],[1,1])
    finally
        A.finalize()
    end
end

@testset "add_discrete_entity bounds a Curve like gmsh::model::addDiscreteEntity" begin
    try
        A.finalize();A.initialize();A.option("General.Terminal",0);A.model.add("declared")
        foreach(tag->A.model.add_discrete_entity(0,tag),1:3)
        # Gmsh looks the two tags up as Points: -3 names none and is dropped.
        A.model.add_discrete_entity(1,41,[1,-3]);A.model.add_discrete_entity(1,42,[1,3])
        A.model.add_discrete_entity(1,43,[3,1]);A.model.add_discrete_entity(1,44,[2,2])
        A.model.add_discrete_entity(2,30,[42,-43])
        @test A.model.get_boundary([(1,41)],false,true,false)==[(0,1)]
        @test A.model.get_boundary([(1,42)],false,true,false)==[(0,1),(0,3)]
        @test A.model.get_boundary([(1,43)],false,true,false)==[(0,3),(0,1)]
        @test A.model.get_boundary([(1,44)],false,true,false)==[(0,2),(0,2)]
        @test A.model.get_boundary([(1,-42)],false,true,false)==[(0,3),(0,1)]
        @test A.model.get_boundary([(2,30)],false,true,false)==[(1,42),(1,-43)]
        @test A.model.get_boundary([(2,-30)],false,true,false)==[(1,42),(1,-43)]
        @test A.model.get_boundary([(2,30)],false,false,true)==[(0,1),(0,3)]
        @test A.model.get_adjacencies(1,44)==(Int[],[2,2])
        A.mesh.add_nodes(0,1,[11],[0.,0.,0.]);A.mesh.add_nodes(0,3,[13],[1.,0.,0.])
        A.mesh.add_nodes(1,42,[12],[0.5,0.,0.])
        # Boundary nodes of a discrete Curve sit at its parametric bounds, and
        # `discreteEdge::parBounds` of a Curve without geometry ends at
        # double(size_t(0)-1); those of a discrete Surface get (-1, -1).
        tags,_,parameters=A.mesh.get_nodes(1,42,true,true)
        @test tags==[12,11,13] && parameters==[0.0,1.8446744073709552e19]
        tags,_,parameters=A.mesh.get_nodes(2,30,true,true)
        @test tags==[12,11,13] && parameters==fill(-1.0,6)
    finally
        A.finalize()
    end
end

@testset "get_nodes appends boundary nodes in Gmsh's order" begin
    try
        load("square_o1")
        tags,coordinates,parameters=A.mesh.get_nodes(2,1,true,true)
        @test tags==[1,2,3,4] && parameters==fill(-1.0,8)
        @test coordinates==[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.]
        tags,_,parameters=A.mesh.get_nodes(1,2,true,true)
        @test tags==[2,3] && parameters==[0.0,1.8446744073709552e19]
        load("cube_o1")
        # Points of a Volume or Surface come last, in tag order.
        @test A.mesh.get_nodes(3,1,true,false)[1]==1:8
        @test A.mesh.get_nodes(2,1,true,false)[1]==1:4
        load("disk_o2")
        # The single Point of a closed Curve is its begin and its end.
        tags,_,parameters=A.mesh.get_nodes(1,1,true,true)
        @test tags==[2,3,4,5,6,7,8,1,1] && parameters==[0.0,0.0]
        @test A.mesh.get_nodes(2,1,true,false)[1]==[9,2,3,4,5,6,7,8,1]
        @test A.mesh.get_nodes(2,1,true,true)[3]==fill(-1.0,16)
    finally
        A.finalize()
    end
end

@testset "periodic queries on a loaded mesh" begin
    try
        load("periodic")
        affine=[1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.]
        # Relations of the $Periodic section, for Points and Curves alike.
        link=A.mesh.get_periodic_nodes(1,2)
        @test link.master_entity==4 && link.affine==affine
        @test sort!(collect(zip(Int.(link.slave_nodes),Int.(link.master_nodes))))==[(2,1),(3,4),(6,8)]
        link=A.mesh.get_periodic_nodes(0,3)
        @test link.master_entity==4 && link.slave_nodes==[3] && link.master_nodes==[4]
        # An entity that is its own master answers (tag, [], [], []).
        for (dim,tag) in ((0,1),(1,1),(1,3),(2,1))
            own=A.mesh.get_periodic_nodes(dim,tag)
            @test own.master_entity==tag && isempty(own.slave_nodes) &&
                  isempty(own.master_nodes) && isempty(own.affine)
        end
        @test A.mesh.get_periodic(1,[1,2,3,4])==[1,4,3,4]
        @test A.mesh.get_periodic(0,[1,2,3,4])==[1,1,4,4]
        # Keys: Gmsh writes the master coordinates into the first array.
        master,types,master_types,keys,master_keys,coordinates,own=
            A.mesh.get_periodic_keys(1,"Lagrange",2)
        @test master==4 && types==zeros(4) && master_types==zeros(4)
        @test keys==[2,6,6,3] && master_keys==[1,8,8,4]
        @test coordinates==[0.,0.,0.,0.,0.5,0.,0.,0.5,0.,0.,1.,0.]
        @test own==[1.,0.,0.,1.,0.5,0.,1.,0.5,0.,1.,1.,0.]
        # A non-periodic entity answers before the space is even looked at.
        for space in ("H1Legendre2","Nonsense")
            answer=A.mesh.get_periodic_keys(1,space,1)
            @test answer[1]==1 && all(isempty,answer[2:end])
        end
        @test_throws ArgumentError A.mesh.get_periodic_keys(1,"H1Legendre2",2)
        load("pyramid")
        own=A.mesh.get_periodic_nodes(3,1)
        @test own.master_entity==1 && isempty(own.slave_nodes) && isempty(own.affine)
    finally
        A.finalize()
    end
end

@testset "changing a loaded mesh that has physical groups" begin
    try
        load("square_phys")
        # The file carries the groups (0,30), (1,20) and (2,10). An automatic
        # tag follows the largest tag of its own dimension until groups have
        # been added through the API, whose counter spans the dimensions; a
        # negative tag is automatic too.
        @test A.model.add_physical_group(2,[1])==11
        @test A.model.add_physical_group(1,[1,2])==21
        @test A.model.add_physical_group(0,[2])==31
        @test A.model.add_physical_group(2,[1])==32
        @test A.model.add_physical_group(1,[3];tag=-7)==33
        @test A.model.get_physical_groups()==
            [(0,30),(0,31),(1,20),(1,21),(1,33),(2,10),(2,11),(2,32)]
        # renumberMeshVertices parks every node on count+1 when physical
        # groups may restrict what is saved, and getMaxNodeTag keeps it.
        load("square_phys")
        @test A.mesh.get_max_node_tag()==4
        A.mesh.renumber_nodes()
        @test A.mesh.get_max_node_tag()==5 && A.mesh.get_nodes()[1]==1:4
        A.mesh.renumber_nodes()
        @test A.mesh.get_max_node_tag()==5
        load("square_phys")
        A.mesh.renumber_nodes([1,2],[2,1])
        @test A.mesh.get_max_node_tag()==5
        @test A.mesh.get_nodes(0,1)[1]==[2] && A.mesh.get_nodes(0,2)[1]==[1]
        load("square_phys");A.option("Mesh.SaveAll",1)
        A.mesh.renumber_nodes()
        @test A.mesh.get_max_node_tag()==4
        load("square_phys");A.option("Mesh.SaveAll",0)
        A.model.remove_physical_groups()
        A.mesh.renumber_nodes()
        @test A.mesh.get_max_node_tag()==4
        # Discrete entities embed and are embedded; Points are listed first,
        # by tag, then the Curves in the order they were embedded. A Point or
        # a Curve embeds nothing, whether it exists or not.
        load("square_phys")
        @test isempty(A.mesh.get_embedded(2,1)) && isempty(A.mesh.get_embedded(0,1)) &&
              isempty(A.mesh.get_embedded(1,99))
        A.mesh.embed(1,[2],2,1);A.mesh.embed(0,[2],2,1);A.mesh.embed(0,[1],2,1)
        @test A.mesh.get_embedded(2,1)==[(0,1),(0,2),(1,2)]
        A.mesh.embed(1,[1],2,1)
        @test A.mesh.get_embedded(2,1)==[(0,1),(0,2),(1,2),(1,1)]
        A.model.add_discrete_entity(2,7)
        @test isempty(A.mesh.get_embedded(2,7))
        @test_throws ArgumentError A.mesh.get_embedded(2,99)
        # A discrete Point keeps a prescribed size.
        @test A.mesh.get_sizes([(0,1),(0,2)])==[0.0,0.0]
        A.mesh.set_size([(0,1)],0.1)
        @test A.mesh.get_sizes([(0,1),(0,2)])==[0.1,0.0]
    finally
        A.finalize()
    end
end

end # module
