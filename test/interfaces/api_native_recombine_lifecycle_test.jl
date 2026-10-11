module APINativeRecombineLifecycleTests
using Test,Tessella,SHA,TOML,Serialization
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
const G=Tessella.GmshLibm
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
bits(x)=reinterpret(UInt64,Float64.(x))
fingerprint(x)=begin io=IOBuffer();serialize(io,x);take!(io);end
const ALGORITHM_TOML=joinpath(@__DIR__,"../artifacts/native_recombine_algorithm_api_primary.toml")
const ALGORITHM_JSON=joinpath(@__DIR__,"../artifacts/native_recombine_algorithm_api_primary.json")
const LIFECYCLE_TOML=joinpath(@__DIR__,"../artifacts/native_recombine_lifecycle_primary.toml")
const LIFECYCLE_JSON=joinpath(@__DIR__,"../artifacts/native_recombine_lifecycle_primary.json")
const LOOKUP_TOML=joinpath(@__DIR__,"../artifacts/native_recombine_lookup_primary.toml")
const LOOKUP_JSON=joinpath(@__DIR__,"../artifacts/native_recombine_lookup_primary.json")
function start()
    call(A.initialize)
    # Native Msg errors persist across initialize/finalize. An empty generation
    # resets that process state before this independent fixture is constructed.
    call(A.mesh.generate,0)
    @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
    @test call(A.logger.get_last_error)==""
end
function check_state(expected)
    types,tags,nodes=call(A.mesh.get_elements)
    @test types==expected["types"]
    @test tags==[UInt64.(row) for row in expected["tags"]]
    @test nodes==[UInt64.(row) for row in expected["nodes"]]
    ids=first(call(A.mesh.get_nodes))
    @test Set(ids)==Set(UInt64.(expected["node_tags"]))
    for (index,uid) in enumerate(expected["node_tags"])
        xyz,parameters,dimension,owner=call(A.mesh.get_node,uid)
        @test bits(xyz)==bits(expected["xyz"][index])
        @test bits(parameters)==bits(expected["parameters"][index])
        @test (dimension,owner)==Tuple(expected["owners"][index])
    end
    @test call(A.mesh.get_max_node_tag)==expected["max_node"]
    @test call(A.mesh.get_max_element_tag)==expected["max_element"]
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
@testset "pinned raw public lifecycle projections" begin
    @test bytes2hex(sha256(read(ALGORITHM_JSON)))=="20d5296f0ace96e3e1fd8de95c61383116a8d826850835dae59e0b44df881671"
    @test bytes2hex(sha256(read(ALGORITHM_TOML)))=="41a351f9f372e52d86b4af3d678c256572a6514505d5a2066e5434f0ee238f98"
    @test bytes2hex(sha256(read(LIFECYCLE_JSON)))=="b453193fe0a75ff9c3fb1519961bf1afc9127d116f1e9fc8046b51fe6a66b69d"
    @test bytes2hex(sha256(read(LIFECYCLE_TOML)))=="e5ffa12005eebaa0e62cf06041a47db3af607a64ba13f6b6421984de3fe56c5f"
    @test bytes2hex(sha256(read(LOOKUP_JSON)))=="12a1afae5689c27c7a6d3842a0d0d4b961f16dc7ba363cadd8a286686b46eb11"
    @test bytes2hex(sha256(read(LOOKUP_TOML)))=="c88f302c9c02680186476ebb7d9d08191859baa49d5a467924ec3668cc1112ff"
end
@testset "all pinned public selectors preserve native births and directed products" begin
    corpus=TOML.parsefile(ALGORITHM_TOML)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(ALGORITHM_JSON)))
    @test length(corpus["cases"])==20
    for case in corpus["cases"]
        start()
        try
            call(A.model.add_discrete_entity,2,1)
            before=case["before"]
            coordinates=vec(reduce(hcat,Float64.(x) for x in before["xyz"]))
            call(A.mesh.add_nodes,2,1,before["node_tags"],coordinates)
            for (msh,tags,nodes) in zip(before["types"],before["tags"],before["nodes"])
                call(A.mesh.add_elements_by_type,1,msh,tags,nodes)
            end
            call(A.mesh.set_recombine,2,1,0.)
            for (name,value) in (("Mesh.Renumber",0),("Mesh.RecombinationAlgorithm",case["algorithm"]),
                ("Mesh.RecombineNodeRepositioning",0),("Mesh.RecombineOptimizeTopology",0),("Mesh.RecombineMinimumQuality",0))
                call(A.option,name,value)
            end
            check_state(before)
            word=G._win_x87_control_word()
            call(A.mesh.recombine)
            check_state(case["after"])
            @test G._win_x87_control_word()==word
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
        finally
            call(A.finalize)
        end
    end
end
@testset "prior Msg skip, dense numbering, no-Surface and empty-Surface lifecycle" begin
    corpus=TOML.parsefile(LIFECYCLE_TOML)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(LIFECYCLE_JSON)))
    @test length(corpus["cases"])==18
    for case in corpus["cases"]
        start()
        try
            plane_fixture(case)
            call(A.option,"Mesh.Renumber",case["renumber"])
            call(A.option,"Mesh.RecombinationAlgorithm",1)
            call(A.option,"General.AbortOnError",case["abort"])
            @test call(A.option,"Mesh.RecombineNodeRepositioning")==1
            @test call(A.option,"Mesh.RecombineOptimizeTopology")==5
            @test call(A.option,"Mesh.Smoothing")==1
            @test call(A.option,"Mesh.RecombineMinimumQuality")==.01
            check_state(case["before"])
            if case["error"]
                if case["abort"]==0
                    call(A.mesh.get_node,999999)
                else
                    @test_throws ErrorException call(A.mesh.get_node,999999)
                end
            end
            @test call(A.logger.get_last_error)==case["last_before"]
            errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
            @test errors==Int(case["error"])
            word=G._win_x87_control_word()
            call(A.mesh.recombine)
            check_state(case["after"])
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==errors
            @test call(A.logger.get_last_error)==case["last_after"]
            call(A.mesh.recombine)
            check_state(case["repeat"])
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==errors
            @test call(A.logger.get_last_error)==case["last_after"]
            @test G._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
    end
end
@testset "native primed lookup policy follows automatic numbering" begin
    corpus=TOML.parsefile(LOOKUP_TOML)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(LOOKUP_JSON)))
    @test length(corpus["cases"])==2
    for case in corpus["cases"]
        start()
        try
            call(A.model.add_discrete_entity,2,1)
            call(A.mesh.add_nodes,2,1,[101,102,103,104],[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.])
            call(A.mesh.add_elements_by_type,1,2,[3001,3002],[101,102,103,101,103,104])
            for (name,value) in (("Mesh.Renumber",case["renumber"]),("Mesh.RecombineNodeRepositioning",0),
                ("Mesh.RecombineOptimizeTopology",0),("Mesh.RecombineMinimumQuality",0))
                call(A.option,name,value)
            end
            @test call(A.mesh.get_node,101)[3:4]==(2,1)
            @test first(call(A.mesh.get_element,3001))==2
            call(A.mesh.recombine)
            @test call(A.mesh.get_elements)==(case["types"],[UInt64.(x) for x in case["labels"]],
                [UInt64.(x) for x in case["cells"]])
            xyz,parameters,dimension,owner=call(A.mesh.get_node,case["node_tag"])
            @test bits(xyz)==bits(case["node"][1])
            @test bits(parameters)==bits(case["node"][2])
            @test (dimension,owner)==(case["node"][3],case["node"][4])
            if case["element_accepted"]
                msh,nodes,dimension,owner=call(A.mesh.get_element,case["element_tag"])
                @test (msh,nodes,dimension,owner)==(case["element"][1],UInt64.(case["element"][2]),
                    case["element"][3],case["element"][4])
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
            else
                @test_throws ErrorException call(A.mesh.get_element,case["element_tag"])
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==1
            end
            @test call(A.logger.get_last_error)==case["last_error"]
        finally
            call(A.finalize)
        end
    end
end
@testset "detached numbering preserves constructor receipts and raw empty matrices" begin
    start()
    try
        case=first(TOML.parsefile(LIFECYCLE_TOML)["cases"])
        plane_fixture(case)
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
        case=first(TOML.parsefile(LIFECYCLE_TOML)["cases"])
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

@testset "native Msg errors cross initialize/finalize until empty generation" begin
    raw=joinpath(@__DIR__,"../artifacts/native_recombine_legacy_order_primary.json")
    projection=joinpath(@__DIR__,"../artifacts/native_recombine_legacy_order_primary.toml")
    @test bytes2hex(sha256(read(raw)))=="a19555221379428a65e42c403a1ea8d9e631ce5a75fb570d7b683c1bcce7cd13"
    @test bytes2hex(sha256(read(projection)))=="4b86227af0f9999b107a2a2540c2986ea086af357c4255a5552fde5482f066c1"
    control=TOML.parsefile(projection)
    @test control["raw_sha256"]==bytes2hex(sha256(read(raw)))
    @test control["fixture_cases"]==24
    fixture=joinpath(@__DIR__,"../artifacts/native_recombine_degree_two_primary.tsv")
    @test bytes2hex(sha256(read(fixture)))=="46dc14373b2d83dc691513130bf3e051d7d7848ab05963a608b0ab839d2a552d"
    source=readlines(fixture)[2:end]
    function check_projection(expected)
        @test call(A.mesh.get_elements,2,1)==(expected["types"],
            [UInt64.(row) for row in expected["tags"]],
            [UInt64.(row) for row in expected["cells"]])
        @test first(call(A.mesh.get_nodes))==UInt64.(expected["nodes"])
        @test call(A.mesh.get_max_node_tag)==expected["max_node"]
        @test call(A.mesh.get_max_element_tag)==expected["max_element"]
    end
    for (index,expected) in enumerate(control["cases"])
        index==1 ? start() : call(A.initialize)
        try
            @test call(A.logger.get_last_error)==expected["last_on_initialize"]
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==index-1
            fields=split(source[expected["index"]+1],'\t')
            labels=parse.(UInt64,split(fields[4],','))
            coordinates=reshape(parse.(Float64,split(fields[2],',')),3,:)
            owners=[parse.(Int,split(value,',')) for value in split(fields[3],';')]
            cells=reduce(vcat,(parse.(Int,split(value,',')) for value in split(fields[5],';')))
            call(A.option,"Mesh.Renumber",0)
            call(A.model.add_discrete_entity,2,1)
            for (row,(dimension,owner)) in enumerate(owners)
                dimension==0 && call(A.model.add_discrete_entity,0,owner)
                call(A.mesh.add_nodes,dimension,owner,[labels[row]],coordinates[:,row])
            end
            call(A.mesh.add_elements_by_type,1,3,parse.(UInt64,split(fields[6],',')),labels[cells])
            check_projection(expected["before"])
            call(A.mesh.recombine)
            check_projection(expected["after"])
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==index-1
            if expected["deleted_accepted"]
                @test call(A.mesh.get_node,expected["deleted_node"])[3:4]==(2,1)
            else
                @test_throws ErrorException call(A.mesh.get_node,expected["deleted_node"])
            end
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==1
            @test call(A.logger.get_last_error)==expected["last_after"]
        finally
            call(A.finalize)
        end
    end
    call(A.initialize)
    try
        @test call(A.logger.get_last_error)==control["last_before_empty_generation_reset"]
        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==1
        call(A.mesh.generate,0)
        @test call(A.logger.get_last_error)==control["last_after_empty_generation_reset"]
        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
        @test call(A.option,"Mesh.Renumber")==control["default_renumber"]
        call(A.model.add_discrete_entity,2,1)
        call(A.mesh.add_nodes,2,1,[101,102,103,104],[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.])
        call(A.mesh.add_elements_by_type,1,3,[301],[101,102,103,104])
        check_projection(control["sparse_quad_before"])
        call(A.mesh.recombine)
        check_projection(control["sparse_quad_after"])
        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
    finally
        call(A.finalize)
    end
end
end
end
