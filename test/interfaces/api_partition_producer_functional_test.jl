module PartitionProducerFunctionalChecks
using Test, LinearAlgebra, Tessella
const A=Tessella.API
const E=Tessella.Elements

function fixture(;lower=false)
    A.finalize();A.initialize();A.model.add("two tets")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    A.model.add_discrete_entity(3,11);A.model.add_discrete_entity(3,22)
    A.model.add_physical_group(3,[11];tag=101,name="left")
    A.model.add_physical_group(3,[22];tag=202,name="right")
    coordinates=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
    A.mesh.add_nodes(3,11,[11,12,13,14,15],vec(coordinates))
    A.mesh.add_elements_by_type(11,4,[501],[11,12,13,14])
    A.mesh.add_elements_by_type(22,4,[502],[11,13,12,15])
    if lower
        A.model.add_discrete_entity(2,30)
        A.mesh.add_elements_by_type(30,2,[503],[11,13,12])
    end
    A.mesh.generate(3)
    return nothing
end

snapshot()=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
    A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag(),A.model.get_number_of_partitions())
function family_cells(dim,owner,kind)
    families,uids,vertices=A.mesh.get_elements(dim,owner)
    index=findfirst(==(kind),families)
    index===nothing && return (UInt64[],zeros(UInt64,E.msh_num_nodes(kind),0))
    return uids[index],reshape(vertices[index],E.msh_num_nodes(kind),length(uids[index]))
end
function records(mesh)
    data=mesh.entity_data
    Dict(uid=>(Int(block.msh),Int(owner),Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (block,tags,owners) in zip(mesh.blocks,data.external_element_tags,data.block_entities)
        for (column,(uid,owner)) in enumerate(zip(tags,owners)))
end
function partition_records(mesh)
    part=mesh.partition_data
    Dict(key=>(entity.parent,entity.partitions,entity.boundaries,entity.physical_tags)
        for (key,entity) in part.entities)
end
function ghost_records(mesh)
    tags=reduce(vcat,mesh.entity_data.external_element_tags;init=UInt64[])
    [(tags[record.element],record.partition,record.ghost_partitions)
        for record in mesh.partition_data.ghost_elements]
end

@testset "actual two-Tet partition interface, owners and ghost aliases" begin
    try
        fixture()
        @test A.option("Mesh.PartitionCreateTopology")==0.0 && A.option("Mesh.PartitionCreateGhostCells")==0.0
        @test A.mesh.get_element(501)[3:4]==(3,11)
        A.option("Mesh.PartitionCreateTopology",1);A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[501,502],[1,2])
        @test A.model.get_number_of_partitions()==2
        @test A.model.get_entities(2)==[(2,1)]
        @test A.model.get_parent(2,1)==(3,11) && A.model.get_partitions(2,1)==[1,2]
        @test A.model.get_parent(3,23)==(3,11) && A.model.get_partitions(3,23)==[1]
        @test A.model.get_parent(3,24)==(3,22) && A.model.get_partitions(3,24)==[2]
        @test A.model.get_boundary([(3,23)],false,true,false)==[(2,1)]
        @test A.model.get_boundary([(3,24)],false,true,false)==[(2,-1)]
        @test isempty(A.model.get_boundary([(2,1)],false,true,false))
        @test family_cells(2,1,2)==(UInt64[503],reshape(UInt64[11,13,12],3,1))
        @test family_cells(3,23,4)==(UInt64[501],reshape(UInt64[11,12,13,14],4,1))
        @test family_cells(3,24,4)==(UInt64[502],reshape(UInt64[11,13,12,15],4,1))
        @test isempty(first(family_cells(3,11,4))) && isempty(first(family_cells(3,22,4)))
        @test Set(first(A.mesh.get_nodes(2,1)))==Set(UInt64[11,12,13])
        @test first(A.mesh.get_nodes(3,23))==UInt64[14] && first(A.mesh.get_nodes(3,24))==UInt64[15]
        @test Set(A.model.get_entities_for_physical_group(3,101))==Set([11,23])
        @test Set(A.model.get_entities_for_physical_group(3,202))==Set([22,24])
        @test A.mesh.get_ghost_elements(3,25)==(UInt64[502],Int32[2])
        @test A.mesh.get_ghost_elements(3,26)==(UInt64[501],Int32[1])
        @test family_cells(3,25,4)==family_cells(3,24,4)
        @test family_cells(3,26,4)==family_cells(3,23,4)
        @test first(family_cells(3,-1,4))==UInt64[501,502,502,501]
        # Lookup owner cache has its own public lifecycle; partition publication
        # does not silently change a previously cached element binding.
        @test A.mesh.get_element(501)[3:4]==(3,11)
        A.mesh.rebuild_element_cache(false)
        @test A.mesh.get_element(501)[3:4]==(3,26) && A.mesh.get_element(502)[3:4]==(3,25)
        owned=A.mesh.get();material=[block for block in owned.blocks if block.msh==4]
        @test sum(size(block.nodes,2) for block in material)==2
        volumes=[det(hcat(owned.coords[:,cell[2]]-owned.coords[:,cell[1]],
            owned.coords[:,cell[3]]-owned.coords[:,cell[1]],owned.coords[:,cell[4]]-owned.coords[:,cell[1]]))/6
            for block in material for cell in eachcol(block.nodes)]
        @test all(>(0),volumes) && sum(volumes)==1/3
        before=snapshot();owned.coords[1,1]=99.
        @test snapshot()==before
        owned=A.mesh.get()
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"owned-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                reread=E.read_mixed_msh(path)
                @test records(reread)==records(owned) && partition_records(reread)==partition_records(owned)
                @test ghost_records(reread)==ghost_records(owned)
            end
            ids,_,_=A.mesh.get_nodes();view=A.view.add("partition values")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,[[Float64(uid)] for uid in ids])
            for binary in (0,1,0)
                A.option("Mesh.Binary",binary);A.option("PostProcessing.Binary",binary)
                path=joinpath(directory,"view-$(binary).msh");A.view.write(view,path)
                format=open(path) do io;readline(io);split(readline(io))[2];end
                @test format==string(binary)
                reread=E.read_mixed_msh(path)
                @test records(reread)==records(owned) && partition_records(reread)==partition_records(owned)
                @test ghost_records(reread)==ghost_records(owned)
                @test length(reread.data_sections)==1 && Set(only(reread.data_sections).nodes)==Set(1:5)
            end
            source=joinpath(directory,"roundtrip.geo")
            write(source,"Merge \""*replace(joinpath(directory,"view-1.msh"),'\\'=>'/')*"\";\n")
            A.model.add("roundtrip");A.open_geo!(source)
            @test A.model.get_parent(3,23)==(3,11) && A.model.get_partitions(2,1)==[1,2]
            @test A.mesh.get_ghost_elements(3,25)==(UInt64[502],Int32[2]) &&
                  A.mesh.get_ghost_elements(3,26)==(UInt64[501],Int32[1])
            @test first(family_cells(3,-1,4))==UInt64[501,502,502,501]
            A.model.set_current("two tets")
            A.option("Mesh.MshFileVersion",2.2)
            protected=joinpath(directory,"unsupported.msh");write(protected,"retained")
            @test_throws ArgumentError A.view.write(view,protected)
            @test read(protected,String)=="retained"
        end
        A.mesh.unpartition()
        @test A.model.get_number_of_partitions()==0 && A.model.get_entities()==[(3,11),(3,22)]
        @test family_cells(3,11,4)==(UInt64[501],reshape(UInt64[11,12,13,14],4,1))
        @test family_cells(3,22,4)==(UInt64[502],reshape(UInt64[11,13,12,15],4,1))
        @test first(family_cells(3,-1,4))==UInt64[501,502]
    finally
        A.finalize()
    end
end

@testset "ghost-only production and existing assignment default" begin
    try
        fixture();A.mesh.partition(2,[501,502],[1,2])
        # The entity layer materializes one child per (dim,parent,partition)
        # even with topology/ghosts off (`9417a32` / upstream
        # `createNewEntities`) — the children list beside the parents.
        @test A.model.get_entities()==[(3,11),(3,22),(3,23),(3,24)] &&
              A.model.get_number_of_partitions()==2
        @test A.model.get_partitions(3,11)==[1] && A.model.get_partitions(3,22)==[2]
        A.mesh.unpartition();A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[501,502],[1,2])
        @test isempty(A.model.get_entities(2)) && A.model.get_parent(3,23)==(3,11)
        @test A.mesh.get_ghost_elements(3,25)==(UInt64[502],Int32[2])
        @test first(family_cells(3,-1,4))==UInt64[501,502,502,501]
    finally
        A.finalize()
    end
end

@testset "unsupported topology requests reject before publication" begin
    try
        fixture();A.option("Mesh.PartitionCreateTopology",1)
        before=snapshot()
        @test_throws ArgumentError A.mesh.partition(3,[501,502],[1,2])
        @test snapshot()==before
        fixture(;lower=true);A.option("Mesh.PartitionCreateTopology",1)
        before=snapshot()
        @test_throws ArgumentError A.mesh.partition(2,[501,502,503],[1,2,1])
        @test snapshot()==before
        # The existing assignment default still handles the same mixed cache.
        A.option("Mesh.PartitionCreateTopology",0)
        @test A.mesh.partition(2,[501,502,503],[1,2,1])===nothing
        @test A.model.get_number_of_partitions()==2
    finally
        A.finalize()
    end
end
end
