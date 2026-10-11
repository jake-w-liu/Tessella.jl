module ModelExportFunctionalChecks
using Test,Tessella
const A=Tessella.API
const E=Tessella.Elements

state()=(A.model.get_current(),A.model.get_file_name(),A.mesh.get_nodes(),
    A.mesh.get_elements(),A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag(),A.view.get_tags())
function lines_model()
    A.model.add("current-export")
    for (entity,x) in ((10,0.),(20,3.))
        A.model.add_discrete_entity(1,entity)
        tags=UInt64[entity*10+3,entity*10+1,entity*10+8]
        A.mesh.add_nodes(1,entity,tags,[x,0,0,x+1,0,0,x+2,0,0])
        A.mesh.add_elements_by_type(entity,1,[entity*100],tags[1:2])
    end
    A.model.add_physical_group(1,[10];tag=7);A.model.set_physical_name(1,7,"owned wire")
end

function retained_fixture()
    coordinates=[0. 1 0 0 0 0 1 3 4;0 0 1 0 0 3 3 3 3;0 0 0 1 -1 0 0 0 0]
    bounds=(0.,0.,-1.,4.,3.,1.)
    entities=Dict((3,1)=>E.MixedEntity(3,1,bounds),
        (3,99)=>E.MixedEntity(3,99,bounds),
        (1,21)=>E.MixedEntity(1,21,(0.,3.,0.,1.,3.,0.)),
        (1,22)=>E.MixedEntity(1,22,(3.,3.,0.,4.,3.,0.)))
    parts=Dict((3,11)=>E.MixedEntity(3,11,bounds;partitions=[1],parent=(3,1)),
        (3,12)=>E.MixedEntity(3,12,bounds;partitions=[2],parent=(3,1)))
    blocks=[E.ElementBlock(4,reshape(Int32[1,2,3,4],4,1)),
        E.ElementBlock(1,Int32[6 8;7 9]),E.ElementBlock(4,reshape(Int32[1,3,2,5],4,1))]
    owners=[Int32[11],Int32[21,22],Int32[12]]
    metadata=E.MixedEntityData(entities;node_entities=vcat(fill((3,1),5),[(1,21),(1,21),(1,22),(1,22)]),node_parametric=fill(nothing,9),
        external_node_tags=UInt64.(1001:1009),block_entities=owners,
        external_element_tags=[UInt64[1017],UInt64[3001,3002],UInt64[1029]])
    partition=E.MixedPartitionData(2;entities=parts,ghost_entities=[(99,1)],
        ghost_elements=[E.MixedGhostElement(4,2,[1])])
    affine=[1. 0 0 3;0 1 0 0;0 0 1 0;0 0 0 1]
    periodic=[E.MixedPeriodicLink(1,22,21,Int32[8,9],Int32[6,7];affine)]
    sections=[E.MshDataSection("NodeData",["\"owned node\""],[.5],Int64[0,1,1],Int32[],Int32[6],Int32[],false,[3.25],4),
        E.MshDataSection("ElementData",["\"owned cell\""],[.5],Int64[0,1,1],Int32[4],Int32[],Int32[],false,[42.],4)]
    ancillary=[E.MshAncillarySection("Comments",false,collect(codeunits("retained export receipt\n")),4)]
    E.MixedMesh(coordinates,blocks;entity_data=metadata,elementary_entities=owners,
        partition_data=partition,periodic_links=periodic,data_sections=sections,ancillary_sections=ancillary)
end

function exercise(directory)
    @testset "current model configuration and detachment" begin
        A.initialize()
        try
            A.option("General.Terminal",0);lines_model()
            for version in (4.1,2.2),binary in (0,1),save_all in (0,1)
                A.option("Mesh.MshFileVersion",version);A.option("Mesh.Binary",binary);A.option("Mesh.SaveAll",save_all)
                before=state();cache=A.LAST_MESH[];classification=A.LAST_MESH_CLASS[]
                path=joinpath(directory,"current-$version-$binary-$save_all.msh")
                @test A.write(path)===nothing
                mesh=E.read_mixed_msh(path);data=mesh.entity_data
                expected_nodes=version==4.1 ? (save_all==1 ? 6 : 3) : (save_all==1 ? 4 : 2)
                @test size(mesh.coords,2)==expected_nodes && sum(size(b.nodes,2) for b in mesh.blocks)==(save_all==1 ? 2 : 1)
                @test mesh.physical_names[(1,7)]=="owned wire"
                @test version==4.1 ? Set(data.external_node_tags)==Set(UInt64.(save_all==1 ? [101,103,108,201,203,208] : [101,103,108])) &&
                    Set(reduce(vcat,data.external_element_tags;init=UInt64[]))==Set(UInt64.(save_all==1 ? [1000,2000] : [1000])) :
                    all(==(save_all==1 ? 0 : 7),vcat([b.tags for b in mesh.blocks]...))
                @test state()==before && A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===classification
            end
            # The API-local filename method must preserve existing stream IO.
            A.option("Mesh.MshFileVersion",4.1);A.option("PostProcessing.Binary",1)
            view=A.view.add("binary continuity");A.view.add_model_data(view,0,"current-export","NodeData",[103],[[2.5]],0.,1)
            path=joinpath(directory,"view-continuity.msh");A.view.write(view,path)
            mesh=E.read_mixed_msh(path)
            @test only(mesh.data_sections).values==[2.5]
            @test mesh.entity_data.external_node_tags[only(only(mesh.data_sections).nodes)]==103
            before=state();path=joinpath(directory,"atomic.msh");Base.write(path,"sentinel")
            A.option("Mesh.MshFileVersion",3.0)
            @test_throws ArgumentError A.write(path)
            @test read(path,String)=="sentinel" && state()==before
            @test_throws ArgumentError A.write(joinpath(directory,"unsupported.vtu"))
            @test_throws ArgumentError A.write("bad\0.msh")
        finally;A.finalize();end
    end
    @testset "large literal node UIDs" begin
        A.initialize()
        try
            A.option("General.Terminal",0);A.model.add("uid-export");A.model.add_discrete_entity(1,7)
            tags=(UInt64(1)<<53).+UInt64[19,3]
            A.mesh.add_nodes(1,7,tags,[0.,0,0,1,0,0]);A.mesh.add_elements_by_type(7,1,[1002],tags)
            for binary in (0,1)
                A.option("Mesh.Binary",binary);path=joinpath(directory,"uids-$binary.msh");before=state()
                A.write(path);mesh=E.read_mixed_msh(path)
                @test Set(mesh.entity_data.external_node_tags)==Set(tags) &&
                    mesh.entity_data.external_node_tags[only(mesh.blocks).nodes[:,1]]==tags
                @test state()==before
            end
        finally;A.finalize();end
    end
    @testset "retained periodic ghost data and literal UID transport" begin
        source=retained_fixture();path=joinpath(directory,"retained.msh")
        E.write_mixed_msh(path,source)
        script=joinpath(directory,"retained.geo");Base.write(script,"Merge \""*replace(path,'\\'=>'/')*"\";\n")
        A.initialize()
        try
            A.option("General.Terminal",0);A.open_geo!(script);A.option("Mesh.SaveAll",1)
            original_ancillary=deepcopy(A.LAST_MESH[].ancillary_sections)
            for binary in (0,1)
                A.option("Mesh.Binary",binary);before=state();path=joinpath(directory,"metadata-$binary.msh")
                A.write(path);mesh=E.read_mixed_msh(path);data=mesh.entity_data
                flat=reduce(vcat,data.external_element_tags;init=UInt64[])
                @test Set(data.external_node_tags)==Set(UInt64.(1001:1009)) && Set(flat)==Set(UInt64[1017,1029,3001,3002])
                link=only(mesh.periodic_links)
                @test data.external_node_tags[link.slave_nodes]==UInt64[1008,1009] && data.external_node_tags[link.master_nodes]==UInt64[1006,1007]
                @test mesh.partition_data.entities[(3,12)].parent==(3,1) && mesh.partition_data.entities[(3,12)].partitions==[2]
                ghost=only(mesh.partition_data.ghost_elements)
                @test flat[ghost.element]==1029 && ghost.partition==2 && ghost.ghost_partitions==[1]
                node_section=only(filter(s->s.name=="NodeData",mesh.data_sections))
                cell_section=only(filter(s->s.name=="ElementData",mesh.data_sections))
                @test data.external_node_tags[only(node_section.nodes)]==1006 && node_section.values==[3.25]
                @test flat[only(cell_section.elements)]==1029 && cell_section.values==[42.]
                # Binary readers distinguish the separator newline from the
                # captured payload. Compare actual published bytes with the
                # literal imported cache payload instead of a guessed suffix.
                original=only(original_ancillary)
                @test length(mesh.ancillary_sections)==1 && occursin(
                    "\$"*original.name*"\n"*String(original.payload),read(path,String))
                @test state()==before
            end
            A.option("Mesh.MshFileVersion",2.2);before=state();path=joinpath(directory,"partition-atomic.msh");Base.write(path,"sentinel")
            @test_throws ArgumentError A.write(path)
            @test read(path,String)=="sentinel" && state()==before
        finally;A.finalize();end
    end
end

@testset "public current-model MSH export" begin
    retained=get(ENV,"TESSELLA_MODEL_EXPORT_RETAIN_DIR",nothing)
    if retained===nothing;mktempdir(exercise);else;mkpath(retained);exercise(retained);end
end
end
