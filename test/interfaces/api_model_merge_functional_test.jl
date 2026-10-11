module ModelMergeFunctionalChecks
using Test,Tessella
const A=Tessella.API
const E=Tessella.Elements

function line_mesh(owner,nodes,cell;value=nothing)
    entities=Dict((1,owner)=>E.MixedEntity(1,owner,(5.,0,0,6,0,0)))
    data=E.MixedEntityData(entities;node_entities=fill((1,owner),2),
        node_parametric=fill(nothing,2),external_node_tags=UInt64.(nodes),
        block_entities=[Int32[owner]],external_element_tags=[UInt64[cell]])
    sections=value===nothing ? E.MshDataSection[] : [E.MshDataSection(
        "NodeData",["\"old cache values\""],[0.],[0,1,1],Int32[],Int32[1],Int32[],false,[value],4)]
    return E.MixedMesh([5. 6;0 0;0 0],[E.ElementBlock(1,reshape(Int32[1,2],2,1),Int32[0])];
        entity_data=data,data_sections=sections)
end

function incoming_wire(;orphan=true)
    count=orphan ? 5 : 4
    uids=UInt64[117,103,119,107,129][1:count]
    entities=Dict((1,7)=>E.MixedEntity(1,7,(0.,0,0,2,0,0);physical_tags=[11]),
        (1,8)=>E.MixedEntity(1,8,(3.,0,0,4,0,0)))
    data=E.MixedEntityData(entities;node_entities=[(1,7),(1,7),(1,8),(1,8),(1,7)][1:count],
        node_parametric=fill(nothing,count),external_node_tags=uids,
        block_entities=[Int32[7,8]],external_element_tags=[UInt64[1001,1002]])
    node=orphan ? 5 : 1
    sections=[E.MshDataSection("NodeData",["\"incoming values\""],[0.],[0,1,1],
            Int32[],Int32[node],Int32[],false,[3.25],4),
        E.MshDataSection("NodeData",["\"incoming values\""],[.5],[1,1,1],
            Int32[],Int32[node],Int32[],false,[4.25],4),
        E.MshDataSection("ElementData",["\"incoming cells\""],[0.],[0,1,1],
            Int32[2],Int32[],Int32[],false,[42.],4),
        E.MshDataSection("ElementNodeData",["\"incoming nodal\""],[0.],[0,1,1],
            Int32[1],Int32[],Int32[2],true,[1.,2.],4)]
    links=[E.MixedPeriodicLink(1,8,7,[3,4],[1,2];affine=[1. 0 0 3;0 1 0 0;0 0 1 0;0 0 0 1])]
    return E.MixedMesh([0. 1 3 4 2;0 0 0 0 0;0 0 0 0 0][:,1:count],
        [E.ElementBlock(1,Int32[1 3;2 4],Int32[11,0])];entity_data=data,
        physical_names=Dict((1,11)=>"incoming physical"),periodic_links=links,data_sections=sections)
end

function partition_mesh(;count=2,offset=0)
    entities=Dict((3,1)=>E.MixedEntity(3,1,(0.,0,0,3,1,1)),
        (3,99)=>E.MixedEntity(3,99,(0.,0,0,3,1,1)))
    children=Dict((3,11)=>E.MixedEntity(3,11,(0.,0,0,1,1,1);partitions=[1],parent=(3,1)),
        (3,12)=>E.MixedEntity(3,12,(2.,0,0,3,1,1);partitions=[2],parent=(3,1)))
    data=E.MixedEntityData(entities;node_entities=vcat(fill((3,11),4),fill((3,12),4)),
        node_parametric=fill(nothing,8),external_node_tags=UInt64.(2001:2008).+UInt64(offset),
        block_entities=[Int32[11,12]],external_element_tags=[UInt64[2017,2029].+UInt64(offset)])
    parts=E.MixedPartitionData(count;entities=children,ghost_entities=[(99,1)],
        ghost_elements=[E.MixedGhostElement(2,2,[1])])
    return E.MixedMesh([0. 1 0 0 2 3 2 2;0 0 1 0 0 0 1 0;0 0 0 1 0 0 0 1],
        [E.ElementBlock(4,Int32[1 5;2 6;3 7;4 8],Int32[0,0])];entity_data=data,partition_data=parts)
end

snapshot()=(models=A.model.list(),name=A.model.get_current(),file=A.model.get_file_name(),
    entities=A.model.get_entities(),nodes=A.mesh.get_nodes(),cells=A.mesh.get_elements(),
    views=A.view.get_tags(),factor=A.option("Mesh.MeshSizeFactor"),
    alias=A.option("Mesh.CharacteristicLengthFactor"),bounds=A.option("General.BoundingBoxSize"),
    log=A.logger.get())

function exercise(directory)
    @testset "literal mesh merge with cache, records and held views" begin
        for (version,binary) in ((4.1,false),(4.1,true),(2.2,false),(2.2,true))
            A.initialize()
            try
                A.option("General.Terminal",0);A.option("Mesh.SaveAll",1)
                A.model.add("older");A.model.add_discrete_entity(1,88)
                A.mesh.add_nodes(1,88,[751,752],[8.,0,0,9,0,0])
                A.mesh.add_elements_by_type(88,1,[701],[751,752])
                older=A.view.add("older view");A.view.add_model_data(older,0,"older","NodeData",[751],[[9.25]])
                cache=joinpath(directory,"current.msh")
                E.write_mixed_msh(cache,line_mesh(98,[401,402],801;value=4.75);version=4.1)
                A.open(cache);A.model.set_file_name("kept-associated.msh")
                cached=last(A.view.get_tags())
                A.model.add_discrete_entity(1,99)
                A.mesh.add_nodes(1,99,[501,502,503],[0.,0,0,1,0,0,2,0,0])
                A.mesh.add_elements_by_type(99,1,[901],[501,502])
                held=A.view.add("held current");A.view.add_model_data(held,0,"current","NodeData",[501],[[7.25]])
                before=snapshot();old_views=A.view.get_tags()
                path=joinpath(directory,"incoming-$version-$binary.msh")
                E.write_mixed_msh(path,incoming_wire(;orphan=version==4.1);version,binary)
                literal=E.read_mixed_msh(path)
                node_uids=literal.entity_data===nothing ? UInt64.(1:E.nnodes(literal)) : literal.entity_data.external_node_tags
                cell_uids=literal.entity_data===nothing ? UInt64[1,2] : reduce(vcat,literal.entity_data.external_element_tags)
                @test A.merge(path)===nothing && A.model.list()==before.models &&
                    A.model.get_current()==before.name && A.model.get_file_name()==before.file
                @test Set(A.mesh.get_nodes()[1])==union(Set(UInt64[401,402,501,502,503]),Set(node_uids)) &&
                    A.mesh.get_elements(1,99)[3]==[UInt64[501,502]] && A.mesh.get_nodes(1,99)[1]==[501,502,503]
                added=setdiff(A.view.get_tags(),old_views)
                @test length(added)==3 && A.view.get_model_data(cached,0)==(["NodeData"],UInt64[401],[4.75]) &&
                    A.view.get_model_data(held,0)==(["NodeData"],UInt64[501],[7.25]) &&
                    A.view.get_model_data(older,0)==(["NodeData"],UInt64[751],[9.25])
                section=first(filter(s->s.name=="NodeData",literal.data_sections))
                @test A.view.get_model_data(added[1],0)==(["NodeData"],node_uids[section.nodes],[3.25]) &&
                    A.view.get_model_data(added[1],1)==(["NodeData"],node_uids[section.nodes],[4.25]) &&
                    A.view.get_model_data(added[2],0)==(["ElementData"],cell_uids[2:2],[42.]) &&
                    A.view.get_model_data(added[3],0)==(["ElementNodeData"],cell_uids[1:1],[1.,2.])
                A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.Binary",binary ? 1 : 0)
                output=joinpath(directory,"merged-$version-$binary.msh");A.write(output)
                restored=E.read_mixed_msh(output);link=only(restored.periodic_links);original=only(literal.periodic_links)
                @test length(restored.data_sections)==5 && restored.entity_data.external_node_tags[link.slave_nodes]==node_uids[original.slave_nodes] &&
                    restored.entity_data.external_node_tags[link.master_nodes]==node_uids[original.master_nodes]
                # Native read-time caches can omit caller nodes added after
                # import; the public SDK requires the same explicit rebuild.
                A.mesh.rebuild_node_cache(false)
                view_file=joinpath(directory,"held-$version-$binary.msh");A.view.write(held,view_file)
                view_mesh=E.read_mixed_msh(view_file)
                @test only(view_mesh.data_sections).values==[7.25] &&
                    Set(view_mesh.entity_data.external_node_tags)==Set(A.mesh.get_nodes()[1])
            finally;A.finalize();end
        end
    end

    @testset "seeded existing geometry, fields and genuine tiny Mesh1" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            initial=joinpath(directory,"initial.geo")
            Base.write(initial,"Point(20)={0,0,0,1};Point(30)={1,0,0,1};Line(40)={20,30};Transfinite Curve{40}=3;Mesh 1;\n")
            A.open(initial);A.model.set_file_name("original-associated.geo")
            A.mesh.field.add("Box",5);A.mesh.field.set_number(5,"VIn",0.25)
            A.mesh.field.set_number(5,"VOut",0.25);A.mesh.field.set_as_background_mesh(5)
            before=snapshot()
            added=joinpath(directory,"existing.geo")
            Base.write(added,"p=newp;Point(p)={2,0,0,1};Line(41)={30,p};Field[5].VIn=0.125;Field[6]=Box;Field[6].VIn=0.5;Field[6].VOut=0.5;Background Field=6;SetName \"renamed merge\";\n")
            @test A.merge(added)===nothing && length(A.model.list())==length(before.models) &&
                A.model.get_current()=="renamed merge" && A.model.get_file_name()==before.file
            @test Set(A.model.get_entities())==Set([(0,20),(0,30),(0,31),(1,40),(1,41)]) &&
                A.model.get_value(0,31,Float64[])==[2.,0,0] && A.mesh.get_nodes()==before.nodes
            @test A.mesh.field.list()==[5,6] && A.mesh.field.get_number(5,"VIn")==0.125 &&
                A.mesh.field.get_number(5,"VOut")==0.25 && A.mesh.field.get_number(6,"VIn")==0.5
            generated=joinpath(directory,"mesh-existing.geo")
            Base.write(generated,"Transfinite Curve{41}=3;Mesh 1;\n")
            @test A.merge(generated)===nothing && length(A.mesh.get_nodes()[1])==5 &&
                length(A.mesh.get_nodes(1,40)[1])==1 && length(A.mesh.get_nodes(1,41)[1])==1
            @test Set(A.mesh.get_elements()[1])==Set([1,15]) &&
                length(reduce(vcat,A.mesh.get_elements()[2]))==7 && A.model.get_file_name()==before.file
            @test A.merge((left=1,),(right=2,))==(left=1,right=2) &&
                A.merge(Dict(:a=>1),Dict(:b=>2),Dict(:c=>3))==Dict(:a=>1,:b=>2,:c=>3)
        finally;A.finalize();end
    end

    @testset "unmeshed periodic geometry and explicit mesh/model clearing" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            geometry=joinpath(directory,"periodic-geometry.geo")
            Base.write(geometry,"Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={2,0,0,1};Point(4)={3,0,0,1};Line(11)={1,2};Line(12)={3,4};Periodic Curve{12}={11} Translate{2,0,0};\n")
            A.open(geometry);name=A.model.get_current();file=A.model.get_file_name()
            added=joinpath(directory,"periodic-added.geo");Base.write(added,"Point(50)={4,0,0,1};\n")
            @test A.merge(added)===nothing && A.model.get_current()==name &&
                A.model.get_file_name()==file && (0,50) in A.model.get_entities()
            @test isempty(A.mesh.get_nodes()[1]) && A.mesh.get_periodic(1,[12])==[11]
            cache=joinpath(directory,"clear-current.msh")
            E.write_mixed_msh(cache,line_mesh(98,[401,402],801;value=4.75);version=4.1)
            A.open(cache);view=last(A.view.get_tags());count=length(A.model.list());file=A.model.get_file_name()
            clear=joinpath(directory,"clear-mesh.geo");Base.write(clear,"Delete Meshes;\n")
            @test A.merge(clear)===nothing && isempty(A.mesh.get_nodes()[1]) &&
                isempty(reduce(vcat,A.mesh.get_elements()[2];init=UInt64[])) && A.model.get_entities()==[(1,98)]
            remove=joinpath(directory,"clear-model.geo");Base.write(remove,"Delete Model;\n")
            @test A.merge(remove)===nothing && isempty(A.model.get_entities()) &&
                length(A.model.list())==count && A.model.get_file_name()==file &&
                A.view.get_model_data(view,0)==(["NodeData"],UInt64[401],[4.75])
        finally;A.finalize();end
    end

    @testset "actual retained transfinite Quad constructor" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            geometry=joinpath(directory,"quad.geo")
            Base.write(geometry,"Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={1,1,0,1};Point(4)={0,1,0,1};Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};Transfinite Curve{:}=2;Transfinite Surface{1};Recombine Surface{1};Mesh 2;\n")
            A.open(geometry);nodes=A.mesh.get_nodes();cells=A.mesh.get_elements();surface=A.mesh.get_elements(2,1)
            added=joinpath(directory,"quad-extra.geo");Base.write(added,"Point(20)={2,0,0,1};\n")
            @test A.merge(added)===nothing && surface[1]==[3] && A.mesh.get_elements(2,1)==surface
            @test A.mesh.get_nodes()==nodes && A.mesh.get_elements()==cells && (0,20) in A.model.get_entities()
        finally;A.finalize();end
    end

    @testset "partition metadata and atomic file failures" begin
        A.initialize()
        try
            A.option("General.Terminal",0);A.option("Mesh.SaveAll",1)
            path=joinpath(directory,"partition.msh");E.write_mixed_msh(path,partition_mesh();version=4.1)
            @test A.merge(path)===nothing && A.model.get_parent(3,11)==(3,1) && A.mesh.get_ghost_elements(3,99)==(UInt64[2029],Int32[2])
            A.model.add_discrete_entity(1,77);A.mesh.add_nodes(1,77,[301,302],[0.,3,0,1,3,0])
            A.mesh.add_elements_by_type(77,1,[401],[301,302])
            line=joinpath(directory,"extra.msh");E.write_mixed_msh(line,line_mesh(78,[311,312],411);version=4.1)
            @test A.merge(line)===nothing && A.model.get_parent(3,12)==(3,1) && Set(A.mesh.get_nodes(1,77)[1])==Set(UInt64[301,302])
            A.option("Mesh.MshFileVersion",4.1);output=joinpath(directory,"partition-merged.msh");A.write(output)
            saved=E.read_mixed_msh(output);ghost=only(saved.partition_data.ghost_elements)
            @test reduce(vcat,saved.entity_data.external_element_tags)[ghost.element]==2029 && ghost.partition==2 && ghost.ghost_partitions==[1]
            A.logger.start();before=snapshot()
            @test_throws ArgumentError A.merge(line)
            @test snapshot()==before
            namespace=joinpath(directory,"namespace.msh")
            E.write_mixed_msh(namespace,partition_mesh(;count=3,offset=10000);version=4.1)
            @test_throws ArgumentError A.merge(namespace)
            missing=joinpath(directory,"missing.msh");@test_throws ArgumentError A.merge(missing)
            bad=joinpath(directory,"bad.msh");Base.write(bad,"\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Nodes\nbad\n")
            @test_throws ArgumentError A.merge(bad)
            failed=joinpath(directory,"failed.geo")
            Base.write(failed,"Mesh.MeshSizeFactor=9;RecombineMesh;Line(900)={9998,9999};\n")
            @test_throws ArgumentError A.merge(failed)
            @test snapshot()==before && A.mesh.get_ghost_elements(3,99)==(UInt64[2029],Int32[2])
        finally;A.finalize();end
    end
end

retained=get(ENV,"TESSELLA_MODEL_MERGE_RETAIN","")
isempty(retained) ? mktempdir(exercise) : (mkpath(retained);exercise(retained))
end
