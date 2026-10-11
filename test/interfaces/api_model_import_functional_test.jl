module ModelImportFunctionalChecks
using Test,Tessella
const A=Tessella.API
const E=Tessella.Elements

state()=(models=A.model.list(),current=A.model.get_current(),file=A.model.get_file_name(),
    entities=A.model.get_entities(),nodes=A.mesh.get_nodes(),cells=A.mesh.get_elements(),
    views=A.view.get_tags(),factor=A.option("Mesh.MeshSizeFactor"),
    alias=A.option("Mesh.CharacteristicLengthFactor"),abort=A.option("General.AbortOnError"),
    context=A.option("General.BoundingBoxSize"),logger=A.logger.get())

function wire_fixture(;large=true,orphan=true)
    all_uids=(large ? UInt64(1)<<53 : UInt64(100)).+UInt64[17,3,19,7,29]
    count=orphan ? 5 : 4;uids=all_uids[1:count];data_node=orphan ? 5 : 1
    entities=Dict((1,7)=>E.MixedEntity(1,7,(0.,0,0,2,0,0);physical_tags=[11]),
        (1,8)=>E.MixedEntity(1,8,(3.,0,0,4,0,0)))
    data=E.MixedEntityData(entities;node_entities=[(1,7),(1,7),(1,8),(1,8),(1,7)][1:count],
        node_parametric=fill(nothing,count),external_node_tags=uids,
        block_entities=[Int32[7,8]],external_element_tags=[UInt64[1001,1002]])
    links=[E.MixedPeriodicLink(1,8,7,[3,4],[1,2];
        affine=[1. 0 0 3;0 1 0 0;0 0 1 0;0 0 0 1])]
    sections=[
        E.MshDataSection("NodeData",["\"temperature\""],[0.],[0,1,1],Int32[],Int32[data_node],Int32[],false,[3.25],4),
        E.MshDataSection("NodeData",["\"temperature\""],[.5],[1,1,1],Int32[],Int32[data_node],Int32[],false,[4.25],4),
        E.MshDataSection("ElementData",["\"cell data\""],[0.],[0,1,1],Int32[2],Int32[],Int32[],false,[42.],4),
        E.MshDataSection("ElementNodeData",["\"linear nodal\""],[0.],[0,1,1],Int32[1],Int32[],Int32[2],true,[1.,2.],4)]
    coords=[0. 1 3 4 2;0 0 0 0 0;0 0 0 0 0]
    return E.MixedMesh(coords[:,1:count],
        [E.ElementBlock(1,Int32[1 3;2 4],Int32[11,0])];entity_data=data,
        physical_names=Dict((1,11)=>"owned wire"),periodic_links=links,data_sections=sections)
end

function partition_fixture()
    coords=[0. 1 0 0 2 3 2 2;0 0 1 0 0 0 1 0;0 0 0 1 0 0 0 1]
    entities=Dict((3,1)=>E.MixedEntity(3,1,(0.,0,0,3,1,1)),
        (3,99)=>E.MixedEntity(3,99,(0.,0,0,3,1,1)))
    children=Dict((3,11)=>E.MixedEntity(3,11,(0.,0,0,1,1,1);partitions=[1],parent=(3,1)),
        (3,12)=>E.MixedEntity(3,12,(2.,0,0,3,1,1);partitions=[2],parent=(3,1)))
    data=E.MixedEntityData(entities;node_entities=vcat(fill((3,11),4),fill((3,12),4)),
        node_parametric=fill(nothing,8),external_node_tags=UInt64.(1001:1008),
        block_entities=[Int32[11,12]],external_element_tags=[UInt64[1017,1029]])
    partition=E.MixedPartitionData(2;entities=children,ghost_entities=[(99,1)],
        ghost_elements=[E.MixedGhostElement(2,2,[1])])
    return E.MixedMesh(coords,[E.ElementBlock(4,Int32[1 5;2 6;3 7;4 8],Int32[0,0])];
        entity_data=data,partition_data=partition)
end

function exercise(directory)
    @testset "actual file imports and held model views" begin
        A.initialize()
        try
            A.option("General.Terminal",0);A.option("Mesh.SaveAll",1)
            A.model.add("kept model");A.model.set_file_name("kept-associated.msh")
            A.model.add_discrete_entity(1,99);A.mesh.add_nodes(1,99,[501,502],[0.,0,0,1,0,0])
            A.mesh.add_elements_by_type(99,1,[901],[501,502])
            held=A.view.add("held old values");A.view.add_model_data(held,0,"kept model","NodeData",[501],[[7.25]])
            before=state()
            for (version,binary) in ((4.1,false),(4.1,true),(2.2,false),(2.2,true))
                # Binary data-section tags have an actual signed 32-bit format
                # width; the ASCII case separately retains literal large UIDs.
                source=wire_fixture(;large=!binary,orphan=version==4.1)
                path=joinpath(directory,"wire-$version-$binary.msh")
                E.write_mixed_msh(path,source;version,binary)
                literal=E.read_mixed_msh(path);data=literal.entity_data
                node_uids=data===nothing ? UInt64.(1:E.nnodes(literal)) : data.external_node_tags
                flat=data===nothing ? UInt64.(1:sum(size(block.nodes,2) for block in literal.blocks)) :
                    reduce(vcat,data.external_element_tags)
                nodal=first(filter(section->section.name=="NodeData",literal.data_sections))
                element=only(filter(section->section.name=="ElementData",literal.data_sections))
                element_nodal=only(filter(section->section.name=="ElementNodeData",literal.data_sections))
                source_link=only(literal.periodic_links)
                previous_models=A.model.list();previous_views=A.view.get_tags()
                @test A.open(path)===nothing
                @test A.model.list()==vcat(previous_models,[splitext(basename(path))[1]]) && A.model.get_file_name()==path
                @test Set(A.mesh.get_nodes()[1])==Set(node_uids)
                @test A.mesh.get_elements(1,7)[3][1]==node_uids[literal.blocks[1].nodes[:,1]]
                imported=setdiff(A.view.get_tags(),previous_views)
                @test length(imported)==3 && A.view.get_model_data(imported[1],0)==(["NodeData"],node_uids[nodal.nodes],[3.25])
                @test A.view.get_model_data(imported[1],1)==(["NodeData"],node_uids[nodal.nodes],[4.25])
                @test A.view.get_model_data(imported[2],0)==(["ElementData"],flat[element.elements],[42.])
                @test A.view.get_model_data(imported[3],0)==(["ElementNodeData"],flat[element_nodal.elements],[1.,2.])
                A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.Binary",binary ? 1 : 0)
                output=joinpath(directory,"roundtrip-$version-$binary.msh");A.write(output)
                restored=E.read_mixed_msh(output);link=only(restored.periodic_links)
                @test restored.entity_data.external_node_tags[link.slave_nodes]==node_uids[source_link.slave_nodes] &&
                    restored.entity_data.external_node_tags[link.master_nodes]==node_uids[source_link.master_nodes]
                @test length(restored.data_sections)==4 && restored.physical_names==literal.physical_names
                view_file=joinpath(directory,"view-$version-$binary.msh");A.view.write(imported[1],view_file)
                saved=E.read_mixed_msh(view_file)
                @test sort([section.reals[1] for section in saved.data_sections])==[0.,.5]
                A.model.set_current("kept model")
                @test A.mesh.get_nodes()==before.nodes && A.mesh.get_elements()==before.cells &&
                    A.model.get_file_name()==before.file && A.view.get_model_data(held,0)==(["NodeData"],UInt64[501],[7.25])
            end
            # A held old model view still writes its own model after new models open.
            path=joinpath(directory,"held-old-view.msh");A.option("PostProcessing.Binary",1);A.view.write(held,path)
            @test Set(E.read_mixed_msh(path).entity_data.external_node_tags)==Set(UInt64[501,502])
            partition_path=joinpath(directory,"partitioned.msh");E.write_mixed_msh(partition_path,partition_fixture();binary=true)
            @test A.open(partition_path)===nothing
            @test A.model.get_parent(3,11)==(3,1) && A.model.get_partitions(3,12)==[2]
            output=joinpath(directory,"partition-roundtrip.msh");A.write(output);restored=E.read_mixed_msh(output)
            ghost=only(restored.partition_data.ghost_elements);uids=reduce(vcat,restored.entity_data.external_element_tags)
            @test uids[ghost.element]==1029 && ghost.partition==2 && ghost.ghost_partitions==[1]
        finally;A.finalize();end
    end
    @testset "GEO publication and atomic preparation" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            held=A.view.add("list lifetime");A.view.add_list_data(held,"SP",1,[0.,0,0,9.])
            geo=joinpath(directory,"a name.geo");Base.write(geo,"Point(61)={2,3,4,1};\n")
            opened=A.open(geo)
            println("GEO initial open: ",(models=A.model.list(),current=A.model.get_current(),views=A.view.get_tags(),held=held))
            @test opened===nothing && A.model.list()==["a name"] && A.view.get_tags()==[held]
            @test A.model.get_entities()==[(0,61)] && A.model.get_file_name()==geo
            @test A.view.get_list_data(held)==(["SP"],[[0.,0,0,9.]])
            @test A.open(geo)===nothing && A.model.list()==["a name","a name"]
            named=joinpath(directory,"named.geo");Base.write(named,"SetName \"declared name\";Point(71)={1,0,0,1};\n")
            @test A.open(named)===nothing && A.model.get_current()=="declared name"
            empty_name=joinpath(directory,"explicit-empty.geo");Base.write(empty_name,"SetName \"\";Point(72)={1,1,0,1};\n")
            @test A.open(empty_name)===nothing && A.model.get_current()==""
            A.model.add("legacy current");legacy=joinpath(directory,"legacy.geo")
            Base.write(legacy,"Mesh.MeshSizeFactor=2.5;Point(81)={0,0,0,1};\n")
            count=length(A.model.list());result=A.open_geo!(legacy)
            @test length(A.model.list())==count && A.model.get_current()=="legacy current"
            result.model.points[81]=(9.,9.,9.)
            @test A.model.get_value(0,81,Float64[])==[0.,0.,0.]
            @test A.option("Mesh.MeshSizeFactor")==2.5 && A.option("Mesh.CharacteristicLengthFactor")==2.5
            A.logger.start();before=state()
            failure=joinpath(directory,"failure.geo")
            Base.write(failure,"Mesh.MeshSizeFactor=8;General.AbortOnError=0;RecombineMesh;Point(91)={0,0,0,1};Point(92)={4,0,0,1};Line(2)={99,100};\n")
            # Empty RecombineMesh produces diagnostic events without generating a mesh.
            for endpoint in (A.open,A.open_geo!)
                @test_throws ArgumentError endpoint(failure)
                @test state()==before
            end
            malformed=joinpath(directory,"malformed.msh");Base.write(malformed,"\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Nodes\nnot a header\n")
            for path in (malformed,joinpath(directory,"missing.msh"),joinpath(directory,"unsupported.step"),"bad\0.geo")
                @test_throws ArgumentError A.open(path)
                @test state()==before
            end
            # A named empty current model is reused, leaving older models alive.
            A.logger.stop();A.model.add("named empty");count=length(A.model.list())
            @test A.open(geo)===nothing && length(A.model.list())==count && A.model.get_current()=="a name"
        finally;A.finalize();end
    end
end

@testset "public atomic model file import" begin
    retained=get(ENV,"TESSELLA_MODEL_IMPORT_RETAIN_DIR",nothing)
    if retained===nothing;mktempdir(exercise);else;mkpath(retained);exercise(retained);end
end
end
