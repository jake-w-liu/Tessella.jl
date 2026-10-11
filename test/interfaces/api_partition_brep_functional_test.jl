module PartitionBRepFunctionalChecks
using Test, LinearAlgebra, Tessella
const A=Tessella.API
const E=Tessella.Elements
const CORNERS=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
const EXTERIOR=UInt64[11 12 13 11 13 12;12 13 11 13 12 11;14 14 14 15 15 15]
const TETS=UInt64[11 11;12 13;13 12;14 15]
const FACES=((1,3,2),(1,2,4),(2,3,4),(3,1,4))

function fixture(;shared_original=false,retained_data=false)
    A.finalize();A.initialize();A.model.add("lower BRep")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    if shared_original
        for tag in 30:36;A.model.add_discrete_entity(2,tag);end
        A.model.add_discrete_entity(3,11,[-30,31,32,33])
        A.model.add_discrete_entity(3,22,[30,34,35,36])
        A.mesh.add_nodes(3,11,11:15,vec(CORNERS))
        A.mesh.add_elements_by_type(30,2,[601],[11,12,13])
        for column in 1:6
            A.mesh.add_elements_by_type(30+column,2,[601+column],EXTERIOR[:,column])
        end
        A.mesh.add_elements_by_type(11,4,[501],TETS[:,1])
        A.mesh.add_elements_by_type(22,4,[502],TETS[:,2])
    else
        A.model.add_discrete_entity(2,30)
        A.model.add_discrete_entity(3,11,[30])
        A.model.add_physical_group(2,[30];tag=310,name="skin")
        A.model.add_physical_group(3,[11];tag=111,name="material")
        A.mesh.add_nodes(3,11,11:15,vec(CORNERS))
        A.mesh.add_elements_by_type(30,2,601:606,vec(EXTERIOR))
        A.mesh.add_elements_by_type(11,4,[501,502],vec(TETS))
    end
    A.mesh.generate(3)
    if retained_data
        mktempdir() do directory
            view=A.view.add("retained cell data")
            A.view.add_model_data(view,0,A.model.get_current(),"ElementData",[501,502],[[10.],[20.]])
            path=joinpath(directory,"input.msh");A.view.write(view,path)
            merge=joinpath(directory,"input.geo");write(merge,"Merge \""*replace(path,'\\'=>'/')*"\";\n")
            A.model.add("lower BRep data");A.open_geo!(merge)
        end
    end
end
partition_fixture()=A.mesh.partition(2,[601,602,603,604,605,606,501,502],[1,1,1,2,2,2,1,2])
snapshot()=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
    A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag(),A.model.get_number_of_partitions())
function family(dim,owner,kind)
    types,ids,nodes=A.mesh.get_elements(dim,owner);index=findfirst(==(kind),types)
    index===nothing && return UInt64[],zeros(UInt64,E.msh_num_nodes(kind),0)
    return ids[index],reshape(nodes[index],E.msh_num_nodes(kind),length(ids[index]))
end
support(nodes)=Tuple(sort(collect(nodes)))
cycle_equal(a,b)=any(Tuple(a)==ntuple(i->b[mod1(i+shift,3)],3) for shift in 0:2)
function closed_region(tag)
    _,tets=family(3,tag,4)
    hits=Dict{Tuple,Vector{Tuple}}()
    for cell in eachcol(tets),face in FACES
        uids=ntuple(i->cell[face[i]],3);push!(get!(()->Tuple[],hits,support(uids)),uids)
    end
    exterior=Dict(key=>only(faces) for (key,faces) in hits if length(faces)==1)
    actual=Dict{Tuple,Tuple}()
    for (_,signed) in A.model.get_boundary([(3,tag)],false,true,false)
        _,tris=family(2,abs(signed),2)
        for cell in eachcol(tris)
            uids=signed>0 ? Tuple(cell) : (cell[1],cell[3],cell[2]);key=support(uids)
            haskey(actual,key) && return false
            actual[key]=uids
        end
    end
    keys(actual)==keys(exterior) && all(cycle_equal(actual[key],exterior[key]) for key in keys(actual))
end
function records(mesh)
    data=mesh.entity_data
    Dict(uid=>(Int(block.msh),Int(owner),Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (block,ids,owners) in zip(mesh.blocks,data.external_element_tags,data.block_entities)
        for (column,(uid,owner)) in enumerate(zip(ids,owners)))
end
partition_records(mesh)=Dict(key=>(entity.parent,entity.partitions,entity.boundaries,entity.physical_tags)
    for (key,entity) in mesh.partition_data.entities)
function ghost_records(mesh)
    ids=reduce(vcat,mesh.entity_data.external_element_tags;init=UInt64[])
    [(ids[row.element],row.partition,row.ghost_partitions) for row in mesh.partition_data.ghost_elements]
end
function cell_values(mesh)
    ids=reduce(vcat,mesh.entity_data.external_element_tags;init=UInt64[])
    section=only(filter(section->section.name=="ElementData",mesh.data_sections))
    Dict(ids[row]=>value for (row,value) in zip(section.elements,section.values))
end

@testset "lower Tri3 BRep partition has actual closed children and seam" begin
    try
        fixture(;retained_data=true);original=A.mesh.get_elements();source_model=A.model.get_current()
        A.option("Mesh.PartitionCreateTopology",1);A.option("Mesh.PartitionCreateGhostCells",1)
        partition_fixture()
        @test A.model.get_number_of_partitions()==2
        @test A.model.get_parent(2,31)==(2,30) && A.model.get_partitions(2,31)==[1]
        @test A.model.get_parent(2,32)==(2,30) && A.model.get_partitions(2,32)==[2]
        @test A.model.get_parent(3,12)==(3,11) && A.model.get_partitions(3,12)==[1]
        @test A.model.get_parent(3,13)==(3,11) && A.model.get_partitions(3,13)==[2]
        @test A.model.get_parent(2,33)==(3,11) && A.model.get_partitions(2,33)==[1,2]
        @test A.model.get_parent(1,1)==(2,30) && A.model.get_partitions(1,1)==[1,2]
        @test family(2,31,2)==(UInt64[601,602,603],EXTERIOR[:,1:3])
        @test family(2,32,2)==(UInt64[604,605,606],EXTERIOR[:,4:6])
        @test family(3,12,4)==(UInt64[501],TETS[:,1:1]) && family(3,13,4)==(UInt64[502],TETS[:,2:2])
        @test only(first(family(2,33,2)))==607 && support(vec(last(family(2,33,2))))==(11,12,13)
        @test Set(A.model.get_boundary([(3,12)],false,true,false))==Set([(2,31),(2,33)])
        @test Set(A.model.get_boundary([(3,13)],false,true,false))==Set([(2,32),(2,-33)])
        @test closed_region(12) && closed_region(13)
        line_ids,lines=family(1,1,1)
        @test length(line_ids)==3 && Set(line_ids)==Set(UInt64[608,609,610])
        @test Set(support(cell) for cell in eachcol(lines))==Set([(11,12),(12,13),(11,13)])
        @test all(count(==(uid),lines[1,:])==1 && count(==(uid),lines[2,:])==1 for uid in 11:13)
        @test isempty(A.model.get_boundary([(1,1)],false,true,false))
        @test A.model.get_boundary([(2,31)],false,true,false)==[(1,1)]
        @test A.model.get_boundary([(2,32)],false,true,false)==[(1,-1)]
        @test A.model.get_boundary([(2,33)],false,true,false)==[(1,-1)]
        @test Set(first(A.mesh.get_nodes(1,1)))==Set(UInt64[11,12,13])
        @test first(A.mesh.get_nodes(2,31))==UInt64[14] && first(A.mesh.get_nodes(2,32))==UInt64[15]
        @test Set(A.model.get_entities_for_physical_group(2,310))==Set([30,31,32])
        @test Set(A.model.get_entities_for_physical_group(3,111))==Set([11,12,13])
        @test A.mesh.get_ghost_elements(3,14)==(UInt64[502],Int32[2]) && A.mesh.get_ghost_elements(3,15)==(UInt64[501],Int32[1])
        @test family(3,14,4)==family(3,13,4) && family(3,15,4)==family(3,12,4)
        owned=A.mesh.get();material=[block for block in owned.blocks if block.msh==4]
        volumes=[det(hcat(owned.coords[:,cell[2]]-owned.coords[:,cell[1]],
            owned.coords[:,cell[3]]-owned.coords[:,cell[1]],owned.coords[:,cell[4]]-owned.coords[:,cell[1]]))/6
            for block in material for cell in eachcol(block.nodes)]
        @test volumes==[1/6,1/6] && sum(volumes)==1/3
        @test Tessella.validate(owned).ok && sum(size(block.nodes,2) for block in owned.blocks)==12
        before=snapshot();owned.coords[1,1]=99.;@test snapshot()==before
        owned=A.mesh.get()
        @test ghost_records(owned)==[(UInt64(501),Int32(1),Int32[2]),(UInt64(502),Int32(2),Int32[1])]
        @test cell_values(owned)==Dict(UInt64(501)=>10.,UInt64(502)=>20.)
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"lower-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                imported=E.read_mixed_msh(path)
                @test records(imported)==records(owned) && partition_records(imported)==partition_records(owned)
                @test ghost_records(imported)==ghost_records(owned)
                @test cell_values(imported)==Dict(UInt64(501)=>10.,UInt64(502)=>20.)
            end
            ids,_,_=A.mesh.get_nodes();view=A.view.add("lower values")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,[[Float64(uid)] for uid in ids])
            A.option("Mesh.MshFileVersion",4.1);A.option("PostProcessing.Binary",1)
            path=joinpath(directory,"view.msh");A.view.write(view,path)
            imported=E.read_mixed_msh(path)
            @test records(imported)==records(owned) && partition_records(imported)==partition_records(owned)
            @test ghost_records(imported)==ghost_records(owned)
            merge=joinpath(directory,"roundtrip.geo");write(merge,"Merge \""*replace(path,'\\'=>'/')*"\";\n")
            A.model.add("lower roundtrip");A.open_geo!(merge)
            @test A.model.get_parent(1,1)==(2,30) && A.model.get_partitions(2,33)==[1,2]
            @test closed_region(12) && closed_region(13)
            @test family(1,1,1)==(line_ids,lines) && A.mesh.get_ghost_elements(3,14)==(UInt64[502],Int32[2])
            A.model.set_current(source_model)
        end
        A.mesh.unpartition()
        @test A.model.get_entities()==[(2,30),(3,11)] && A.model.get_number_of_partitions()==0
        @test A.mesh.get_elements()==original
        @test A.model.get_boundary([(3,11)],false,true,false)==[(2,30)]
        @test isempty(A.model.get_boundary([(2,30)],false,true,false)) && closed_region(11)
        @test A.model.get_entities_for_physical_group(2,310)==[30] && A.model.get_entities_for_physical_group(3,111)==[11]
    finally
        A.finalize()
    end
end

@testset "lower BRep ghost-only and duplicate-role rejection are atomic" begin
    try
        fixture();A.option("Mesh.PartitionCreateGhostCells",1);partition_fixture()
        @test isempty(A.model.get_entities(1)) && length(A.model.get_entities(2))==3
        @test A.model.get_boundary([(3,12)],false,true,false)==[(2,31)]
        @test A.mesh.get_ghost_elements(3,14)==(UInt64[502],Int32[2])
        fixture(;shared_original=true);A.option("Mesh.PartitionCreateTopology",1)
        @test A.model.get_boundary([(3,11)],false,true,false)==[(2,-30),(2,31),(2,32),(2,33)]
        before=snapshot()
        caught=try
            A.mesh.partition(2,[601,602,603,604,605,606,607,501,502],[1,1,1,1,2,2,2,1,2]);nothing
        catch error
            error
        end
        @test caught isa ArgumentError && occursin("distinct geometric roles",sprint(showerror,caught))
        @test snapshot()==before
    finally
        A.finalize()
    end
end
end
