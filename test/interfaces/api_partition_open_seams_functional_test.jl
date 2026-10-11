module PartitionOpenSeamFunctionalChecks
using Test, LinearAlgebra, Tessella
const A=Tessella.API
const E=Tessella.Elements
const FACES=((1,3,2),(1,2,4),(2,3,4),(3,1,4))
function fixture()
    A.finalize();A.initialize();A.model.add("open seams")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    for tag in 31:36;A.model.add_discrete_entity(2,tag);end
    A.model.add_discrete_entity(3,11,collect(31:36))
    A.model.add_physical_group(2,collect(31:36);tag=310,name="skin")
    A.model.add_physical_group(3,[11];tag=111,name="material")
    xyz=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
    A.mesh.add_nodes(3,11,11:15,vec(xyz))
    triangles=UInt64[11 12 13 11 13 12;12 13 11 13 12 11;14 14 14 15 15 15]
    for column in 1:6;A.mesh.add_elements_by_type(30+column,2,[600+column],triangles[:,column]);end
    A.mesh.add_elements_by_type(11,4,[501,502],[11,12,13,14,11,13,12,15])
    A.mesh.generate(3)
end
function family(dim,owner,kind)
    types,ids,nodes=A.mesh.get_elements(dim,owner);index=findfirst(==(kind),types)
    index===nothing && return UInt64[],zeros(UInt64,E.msh_num_nodes(kind),0)
    ids[index],reshape(nodes[index],E.msh_num_nodes(kind),length(ids[index]))
end
support(nodes)=Tuple(sort(collect(nodes)))
cycle_equal(a,b)=any(Tuple(a)==ntuple(i->b[mod1(i+shift,3)],3) for shift in 0:2)
function closed_region(tag)
    _,tets=family(3,tag,4);hits=Dict{Tuple,Vector{Tuple}}()
    for cell in eachcol(tets),face in FACES
        uids=ntuple(i->cell[face[i]],3);push!(get!(()->Tuple[],hits,support(uids)),uids)
    end
    exterior=Dict(key=>only(faces) for (key,faces) in hits if length(faces)==1);actual=Dict{Tuple,Tuple}()
    for (_,signed) in A.model.get_boundary([(3,tag)],false,true,false),cell in eachcol(last(family(2,abs(signed),2)))
        uids=signed>0 ? Tuple(cell) : (cell[1],cell[3],cell[2]);key=support(uids)
        haskey(actual,key) && return false
        actual[key]=uids
    end
    keys(actual)==keys(exterior) && all(cycle_equal(actual[key],exterior[key]) for key in keys(actual))
end
function records(mesh)
    data=mesh.entity_data
    Dict(uid=>(Int(block.msh),Int(owner),Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (block,ids,owners) in zip(mesh.blocks,data.external_element_tags,data.block_entities)
        for (column,(uid,owner)) in enumerate(zip(ids,owners)))
end
part_records(mesh)=Dict(key=>(entity.parent,entity.partitions,entity.boundaries,entity.physical_tags)
    for (key,entity) in mesh.partition_data.entities)

@testset "actual open partition seams retain original Surface borders" begin
    try
        fixture();original=A.mesh.get_elements()
        A.option("Mesh.PartitionCreateTopology",1);A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[601,602,603,604,605,606,501,502],[1,1,1,2,2,2,1,2])
        @test A.model.get_number_of_partitions()==2 && isempty(A.model.get_entities(0))
        @test A.model.get_entities(1)==[(1,1),(1,2),(1,3)]
        @test [A.model.get_parent(1,tag) for tag in 1:3]==[(2,33),(2,32),(2,31)]
        @test all(A.model.get_partitions(1,tag)==[1,2] for tag in 1:3)
        @test all(isempty(A.model.get_boundary([(1,tag)],false,true,false)) for tag in 1:3)
        actual_edges=[Tuple(vec(last(family(1,tag,1)))) for tag in 1:3]
        @test actual_edges==[(13,11),(12,13),(11,12)]
        @test [only(first(family(1,tag,1))) for tag in 1:3]==UInt64[608,609,610]
        @test Set(A.model.get_boundary([(2,43)],false,true,false))==Set([(1,-1),(1,-2),(1,-3)])
        @test Set(support(edge) for edge in actual_edges)==Set([(11,12),(11,13),(12,13)])
        @test [A.model.get_parent(2,tag) for tag in 37:42]==[(2,tag) for tag in 31:36]
        @test [only(first(family(2,tag,2))) for tag in 37:42]==UInt64[601,602,603,604,605,606]
        @test all(A.model.get_partitions(2,tag)==[tag<=39 ? 1 : 2] for tag in 37:42)
        @test Set(A.model.get_entities_for_physical_group(2,310))==Set(31:42)
        @test A.model.get_entities_for_physical_group(3,111)==[11,12,13]
        @test closed_region(12) && closed_region(13)
        @test [A.model.get_boundary([(2,tag)],false,true,false) for tag in 37:42]==
            [[(1,3)],[(1,2)],[(1,1)],[(1,-1)],[(1,-2)],[(1,-3)]]
        @test A.mesh.get_ghost_elements(3,14)==(UInt64[502],Int32[2]) &&
              A.mesh.get_ghost_elements(3,15)==(UInt64[501],Int32[1])
        owned=A.mesh.get();material=[block for block in owned.blocks if block.msh==4]
        volumes=[det(hcat(owned.coords[:,cell[2]]-owned.coords[:,cell[1]],
            owned.coords[:,cell[3]]-owned.coords[:,cell[1]],owned.coords[:,cell[4]]-owned.coords[:,cell[1]]))/6
            for block in material for cell in eachcol(block.nodes)]
        @test volumes==[1/6,1/6] && sum(volumes)==1/3 && Tessella.validate(owned).ok
        before=A.mesh.get_nodes();owned.coords[1,1]=99.;@test A.mesh.get_nodes()==before
        owned=A.mesh.get()
        mktempdir() do directory
            path=joinpath(directory,"open.msh");E.write_mixed_msh(path,owned;version=4.1)
            imported=E.read_mixed_msh(path)
            @test records(imported)==records(owned) && part_records(imported)==part_records(owned)
            merge=joinpath(directory,"open.geo");write(merge,"Merge \""*replace(path,'\\'=>'/')*"\";\n")
            A.model.add("open roundtrip");A.open_geo!(merge)
            @test closed_region(12) && closed_region(13)
            @test [A.model.get_parent(1,tag) for tag in 1:3]==[(2,33),(2,32),(2,31)] &&
                  [Tuple(vec(last(family(1,tag,1)))) for tag in 1:3]==actual_edges
            @test A.mesh.get_ghost_elements(3,14)==(UInt64[502],Int32[2])
            A.model.set_current("open seams")
        end
        A.mesh.unpartition()
        @test A.model.get_entities()==vcat([(2,tag) for tag in 31:36],[(3,11)])
        @test A.mesh.get_elements()==original && A.model.get_number_of_partitions()==0
        @test A.model.get_boundary([(3,11)],false,true,false)==[(2,tag) for tag in 31:36]
    finally
        A.finalize()
    end
end
end
