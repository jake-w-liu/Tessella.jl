module ReclassificationLifecycleTests
using Test,Tessella
using Tessella.Elements: MixedMesh,ElementBlock,MixedEntity,MixedEntityData,MshDataSection
const API=Tessella.API

failures=String[]
try
@testset "Malformed reclassification leaves published state intact" begin
    API.initialize()
    try
        API.model.add_discrete_entity(1,1)
        API.mesh.add_nodes(1,1,[11,12],[0.,0.,0.,1.,0.,0.],[.1,.2])
        API.mesh.add_elements_by_type(1,1,[42],[11,12])
        record=API.CURRENT[].discrete[(1,1)]
        record.element_nodes[1][2]=999
        model=API.CURRENT[];parameters=deepcopy(record.node_native_params)
        maxima=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])
        @test_throws r"unknown source node 999" API.mesh.reclassify_nodes()
        @test API.CURRENT[]===model
        @test API.LAST_MESH[]===nothing
        @test API.CURRENT[].discrete[(1,1)]===record
        @test record.node_native_params==parameters
        @test record.element_nodes==[Int32[11,999]]
        @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==maxima
        record.element_nodes[1][2]=12
        API.mesh.reclassify_nodes()
        API.mesh.set_visibility([42],0)
        API.mesh.reclassify_nodes()
        @test API.mesh.get_visibility([42])==Int32[0]
        @test API.mesh.get_elements(1,1)==([1],[UInt64[42]],[UInt64[11,12]])
        @test API.mesh.get_node(12)==([1.,0.,0.],[.2],1,1)
    finally
        API.finalize()
    end
end

catch exception
    push!(failures,sprint(showerror,exception))
    showerror(stderr,exception,catch_backtrace());println(stderr)
end
try
@testset "Reclassification remaps view rows by actual node and element identity" begin
    API.initialize()
    try
        API.model.add_discrete_entity(1,1)
        API.model.add_discrete_entity(2,1)
        coords=Float64[.5 0 1 9;1 0 0 9;0 0 0 9]
        blocks=[ElementBlock(2,reshape(Int32[2,3,1],3,1)),
            ElementBlock(1,reshape(Int32[2,3],2,1))]
        metadata=Dict((1,1)=>MixedEntity(1,1,(0.,0.,0.,1.,0.,0.)),(2,1)=>MixedEntity(2,1,(0.,0.,0.,9.,9.,0.)))
        data=MixedEntityData(metadata;node_entities=fill((2,Int32(1)),4),
            node_parametric=[[.5,1.],[0.,0.],[1.,0.],[9.,9.]],
            external_node_tags=[13,11,12,14],block_entities=[Int32[1],Int32[1]],
            external_element_tags=[UInt64[301],UInt64[201]])
        sections=[MshDataSection("NodeData",["nodes"],[0.],[0,1,2],Int32[],
            Int32[1,3],Int32[],false,[130.,120.],0),
            MshDataSection("ElementData",["elements"],[0.],[0,1,2],Int32[1,2],
                Int32[],Int32[],false,[3010.,2010.],0),
            MshDataSection("ElementNodeData",["supports"],[0.],[0,1,1],Int32[1],
                Int32[2,3,1],Int32[3],false,[110.,120.,130.],0),
            MshDataSection("NodeData",["orphan and vector"],[2.],[7,2,2,8],Int32[],
                Int32[4,1],Int32[],false,[90.,91.,130.,131.],3),
            MshDataSection("ElementNodeData",["implicit vector"],[3.],[8,2,1],Int32[1],
                Int32[],Int32[3],true,[110.,111.,120.,121.,130.,131.],4)]
        partition=Tessella.Elements.MixedPartitionData(2;ghost_entities=[(8,2)],
            ghost_elements=[Tessella.Elements.MixedGhostElement(1,1,[2]),
                Tessella.Elements.MixedGhostElement(2,2,[1])])
        ancillary=[Tessella.Elements.MshAncillarySection("Comments",false,
            collect(codeunits("retained metadata\n")),2)]
        mesh=MixedMesh(coords,blocks;entity_data=data,data_sections=sections,
            elementary_entities=data.block_entities,partition_data=partition,
            ancillary_sections=ancillary)
        class=API._higher_classification_from_data(mesh,"reclassification lifecycle fixture")
        API._replace_mesh_cache_locked!(mesh,class)
        API.mesh.reclassify_nodes()
        result=API.LAST_MESH[]
        @test result.entity_data.external_node_tags==UInt64[11,12,13]
        @test [block.msh for block in result.blocks]==[1,2]
        @test result.data_sections[1].nodes==Int32[3,2]
        @test result.data_sections[1].values==[130.,120.]
        @test result.data_sections[2].elements==Int32[2,1]
        @test result.data_sections[2].values==[3010.,2010.]
        @test result.data_sections[3].elements==Int32[2]
        @test result.data_sections[3].nodes==Int32[1,2,3]
        @test result.data_sections[3].values==[110.,120.,130.]
        @test result.data_sections[4].nodes==Int32[3]
        @test result.data_sections[4].values==[130.,131.]
        @test result.data_sections[4].header==[7,2,1,8]
        @test result.data_sections[4].strings==["orphan and vector"]
        @test result.data_sections[4].reals==[2.]
        @test result.data_sections[4].anchor==3
        @test result.data_sections[5].elements==Int32[2]
        @test result.data_sections[5].implicit_nodes
        @test isempty(result.data_sections[5].nodes)
        @test result.data_sections[5].row_nodes==Int32[3]
        @test result.data_sections[5].values==[110.,111.,120.,121.,130.,131.]
        @test result.partition_data.num_partitions==2
        @test result.partition_data.ghost_entities==[(Int32(8),Int32(2))]
        @test [(g.element,g.partition,g.ghost_partitions) for g in
            result.partition_data.ghost_elements]==[(Int32(2),Int32(1),Int32[2]),
                (Int32(1),Int32(2),Int32[1])]
        @test [(s.name,s.binary,s.payload,s.anchor) for s in result.ancillary_sections]==
            [(s.name,s.binary,s.payload,s.anchor) for s in ancillary]
        @test Tessella.MeshTypes.validate(result).ok
        API.mesh.reclassify_nodes()
        repeated=API.LAST_MESH[]
        @test [(s.name,s.nodes,s.elements,s.row_nodes,s.implicit_nodes,s.values,s.header)
            for s in repeated.data_sections]==
            [(s.name,s.nodes,s.elements,s.row_nodes,s.implicit_nodes,s.values,s.header)
                for s in result.data_sections]
        @test [(g.element,g.partition,g.ghost_partitions) for g in
            repeated.partition_data.ghost_elements]==[(g.element,g.partition,g.ghost_partitions)
                for g in result.partition_data.ghost_elements]
    finally
        API.finalize()
    end
end

catch exception
    push!(failures,sprint(showerror,exception))
    showerror(stderr,exception,catch_backtrace());println(stderr)
end
try
@testset "Reclassification retains a native carrier through actual UID bindings" begin
    API.initialize()
    try
        execution=mktemp() do path,io
            write(io,raw"""
Mesh.Renumber=0;
Point(1)={0,0,0,.2};Point(2)={1,0,0,.2};
Point(3)={1,1,0,.2};Point(4)={0,1,0,.2};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
Transfinite Curve{1:4}=3;Transfinite Surface{1};Mesh 2;
""")
            close(io)
            execute_geo(path)
        end
        # Exercise the shared model/authority protocol using a carrier and
        # explicit UID row bindings produced by the native GEO frontend.
        API.initialize()
        slot=API.MODEL_SLOTS[API._current_slot_index_locked()]
        slot.model=execution.model
        API.CURRENT[]=execution.model
        initial=execution.model.meshing.mesh_authority
        class=API._higher_classification_from_data(initial.mesh,"native carrier lifecycle")
        API._replace_mesh_cache_locked!(initial.mesh,class)
        previous=API.CURRENT[].meshing.mesh_authority
        @test previous.carrier isa Tessella.MeshTypes.Mesh
        carrier=previous.carrier
        labels=previous.mesh.entity_data.external_node_tags[previous.carrier_rows]
        public_before=API.mesh.get_nodes()
        node_values=Dict(tag=>API.mesh.get_node(tag) for tag in public_before[1])
        elements_before=API.mesh.get_elements()
        API.mesh.reclassify_nodes()
        result=API.CURRENT[].meshing.mesh_authority
        @test result.carrier isa Tessella.MeshTypes.Mesh
        @test result.carrier.coords==carrier.coords
        @test result.carrier.segs==carrier.segs
        @test result.carrier.tris==carrier.tris
        @test result.mesh.entity_data.external_node_tags[result.carrier_rows]==labels
        @test result.mesh.coords[:,result.carrier_rows]==result.carrier.coords
        @test all(tag->API.mesh.get_node(tag)==node_values[tag],public_before[1])
        @test API.mesh.get_elements()==elements_before
    finally
        API.finalize()
    end
end

catch exception
    push!(failures,sprint(showerror,exception))
    showerror(stderr,exception,catch_backtrace());println(stderr)
end
isempty(failures) || error("Reclassification lifecycle failures: ",failures)
end
