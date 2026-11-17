module APINativeRecombineLifecycleTests
using Test,Tessella,Serialization
const A=Tessella.API
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
fingerprint(x)=begin io=IOBuffer();serialize(io,x);take!(io);end
function start()
    call(A.initialize)
    # Native Msg errors persist across initialize/finalize. An empty generation
    # resets that process state before this independent fixture is constructed.
    call(A.mesh.generate,0)
    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
    @test call(A.logger.get_last_error)==""
end
function plane_fixture(case)
    corners=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.))
    for (tag,point) in enumerate(corners);call(A.model.add_point,point...;tag,meshSize=.5);end
    for tag in 1:4;call(A.model.add_line,tag,mod1(tag+1,4);tag);end
    if case["surface"]
        call(A.model.add_curve_loop,[1,2,3,4];tag=1)
        call(A.model.add_plane_surface,[1];tag=1)
    end
    for (tag,point) in enumerate(corners);call(A.mesh.add_nodes,0,tag,[100+tag],collect(point));end
    for tag in 1:4
        call(A.mesh.add_elements_by_type,tag,1,[2000+tag],[100+tag,100+mod1(tag+1,4)])
    end
    case["surface"] && !case["empty_surface"] && call(A.mesh.add_elements_by_type,1,2,
        [3001,3002],[101,102,103,101,103,104])
end
@testset "detached numbering preserves constructor receipts and raw empty matrices" begin
    start()
    try
        case=Dict("surface"=>true,"empty_surface"=>false)
        plane_fixture(case)
        @test call(A.option,"Mesh.RecombineNodeRepositioning")==1
        @test call(A.option,"Mesh.RecombineOptimizeTopology")==5
        @test call(A.option,"Mesh.Smoothing")==1
        @test call(A.option,"Mesh.RecombineMinimumQuality")==.01
        call(A.option,"Mesh.Renumber",0)
        call(A.mesh.recombine)
        mesh=A.LAST_MESH[]
        class=A.LAST_MESH_CLASS[]
        model=A.CURRENT[]
        options=copy(A.OPTIONS)
        plan=call(A._api_native_recombine_plan,model,mesh,class,"constructor receipt";
            options,initial_max_node_tag=A.NODE_TAG_MAX[],initial_max_element_tag=A.ELEMENT_TAG_MAX[])
        authority=plan.model.meshing.mesh_authority
        receipt=call(Tessella.Model._model_mesh_product,authority,2,1)
        receipt.blocks[1].nodes[:,:]=reverse(receipt.blocks[1].nodes;dims=1)
        plan.model.meshing.mesh_authority=call(Tessella.Model._model_mesh_with_surface_constructor_inputs,
            authority,Dict((2,Int32(1))=>receipt))
        for shape in (0,1)
            input=deepcopy(plan)
            for (_,record) in Tessella.Model._discrete_mesh_records_model(input.model)
                record.node_params=zeros(0,shape==0 ? 0 : length(record.node_tags))
            end
            before=fingerprint(input)
            output=call(A._api_native_recombine_renumber_plan,input,"constructor receipt";save_all=false)
            @test fingerprint(input)==before
            @test output.model.meshing.mesh_authority.mesh===output.mesh
            @test output.max_node_tag==104 && output.max_element_tag==3004
            @test output.model.meshing.mesh_authority.private_rank==input.model.meshing.mesh_authority.private_rank
            @test isequal(output.native_runtime_parameters,input.native_runtime_parameters)
            updated=only(values(output.model.meshing.mesh_authority.surface_constructor_inputs))
            @test updated.blocks[1].nodes==receipt.blocks[1].nodes
            @test updated.coords==receipt.coords && updated.coords!==receipt.coords
            @test updated.entity_data.external_node_tags==UInt64[output.node_tag_map[tag]
                for tag in receipt.entity_data.external_node_tags]
            @test updated.entity_data.external_element_tags==[UInt64[output.element_tag_map[tag]
                for tag in labels] for labels in receipt.entity_data.external_element_tags]
            for (_,record) in Tessella.Model._discrete_mesh_records_model(output.model)
                @test size(record.node_params)==(0,shape==0 ? 0 : length(record.node_tags))
            end
            output.mesh.coords[1,1]=9.
            @test fingerprint(input)==before
        end
    finally
        call(A.finalize)
    end
end
@testset "invalid private selector context rejects before publication" begin
    start()
    try
        case=Dict("surface"=>true,"empty_surface"=>false)
        plane_fixture(case)
        call(A.option,"Mesh.Renumber",1)
        options=copy(A.OPTIONS)
        options["Mesh.RecombinationAlgorithm"]=5.0
        previous=fingerprint((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],A.NODE_TAG_MAX[],
            A.ELEMENT_TAG_MAX[],A.NATIVE_LOOKUP_QUERY_ERRORS[],A.NATIVE_MSG_LAST_ERROR[]))
        @test_throws ArgumentError call(A._api_native_recombine_plan,A.CURRENT[],A.LAST_MESH[],
            A.LAST_MESH[]===nothing ? nothing : A._cached_classification_locked(A.LAST_MESH[]),"invalid private selector";
            options,initial_max_node_tag=A.NODE_TAG_MAX[],initial_max_element_tag=A.ELEMENT_TAG_MAX[])
        @test previous==fingerprint((A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[],A.NODE_TAG_MAX[],
            A.ELEMENT_TAG_MAX[],A.NATIVE_LOOKUP_QUERY_ERRORS[],A.NATIVE_MSG_LAST_ERROR[]))
    finally
        call(A.finalize)
    end
end

end # module
