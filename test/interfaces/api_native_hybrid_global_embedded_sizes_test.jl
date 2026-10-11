module NativeHybridGlobalEmbeddedSizesChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,validate
using Tessella.Elements: ElementBlock,MixedMesh,MixedEntityData,MixedEntity
const M=Tessella.Model
const artifacts=joinpath(@__DIR__,"..","artifacts")
const native=TOML.parsefile(joinpath(artifacts,"native_hybrid_foreign_embedded_sizes_4152_v1.toml"))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
function helper(expr,name)
    expr isa Expr || return nothing
    if expr.head===:function && expr.args[1] isa Expr && expr.args[1].head===:call && expr.args[1].args[1]===name
        return expr
    end
    for child in expr.args
        found=helper(child,name);found===nothing || return found
    end
    return nothing
end
path=joinpath(@__DIR__,"api_native_hybrid_actual_model_test.jl")
for name in (:box_model,:actual_lower)
    found=helper(Meta.parseall(read(path,String);filename=path),name)
    @assert found!==nothing;Core.eval(@__MODULE__,found)
end
@testset "native global embedded Point and Curve sizing preserves selected Source" begin
    observed=Dict{String,Vector{UInt64}}()
    for case in native["case"]
        m=box_model();m.meshing.lc_extend_from_boundary=2
        # The SDK receipt declares an empty hidden discrete GRegion[2]. It
        # supplies sizing seeds and is outside the selected constructor group.
        m.volumes[2]=Int[]
        if startswith(case["name"],"point")
            m.point_size[1]=.1
        else
            M.add_point!(m,-.05,0,0;tag=99,mesh_size=.4)
            M.add_line!(m,1,99;tag=99)
        end
        embedded=NTuple{2,Int}[Tuple(Int.(entity)) for entity in case["embedded"]]
        isempty(embedded) || (m.embeds[(3,2)]=embedded)
        input=actual_lower(m,case);saved=deepcopy(input.mesh)
        receipt=M._model_native_hybrid_global_size_seeds(m,input,"native global seeds")
        @test receipt.input===input
        @test length(receipt.records)==(isempty(embedded) ? 0 : startswith(case["name"],"point") ? 1 : 2)
        @test M._model_native_hybrid_seed_records(receipt)===receipt.records
        options=M._model_default_mesh1d_options(m,"native global seeds")
        options=M._ModelMesh1DOptions((name===:optimize ? 0 : getfield(options,name)
            for name in fieldnames(M._ModelMesh1DOptions))...)
        options.random_state[]=UInt32(1)
        component=M._model_native_hybrid_component_receipt(m,[1],input,"native global seeds")
        births=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),UInt64(case["current_lower_element_high_water"]))
        batch=M._model_native_hybrid_connected_production(m,component,options,births,"native global seeds")
        production=only(values(batch.products));expected=only(case["expected"])
        @test validate(production.mesh).ok
        @test batch.node_high_water==UInt64(case["node_high_water"])
        @test batch.element_high_water==UInt64(case["element_high_water"])
        ids=UInt64[uid for cell in expected["cell"] for uid in cell["labels"]]
        corners=Tuple.(UInt64.(nodes) for cell in expected["cell"] for nodes in cell["nodes"])
        @test production.element_tags==ids
        @test [Tuple(production.boundary_node_tags[rows]) for block in production.mesh.blocks for rows in eachcol(block.nodes)]==corners
        positions=Dict(uid=>row for (row,uid) in enumerate(production.boundary_node_tags))
        @test isequal(production.mesh.coords[:,[positions[UInt64(uid)] for uid in expected["node_tags"]]],
            reduce(hcat,hexpoint.(expected["xyz_hex"])))
        @test M._model_native_hybrid_volume_certificate(production.mesh,production.boundary_node_tags,production.constructor,"native global seeds")===nothing
        @test isequal(input.mesh.coords,saved.coords)
        @test isequal(input.mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
        @test input.mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
        observed[case["name"]]=production.element_tags
        println("HYBRID_GLOBAL_SEEDS ",case["name"]," records=",length(receipt.records),
            " clocks=",(batch.node_high_water,batch.element_high_water));flush(stdout)
        if !isempty(embedded)
            changed=deepcopy(input)
            row=findfirst(==(first(receipt.records).uid),changed.mesh.entity_data.external_node_tags)
            changed.mesh.coords[1,row]=nextfloat(changed.mesh.coords[1,row])
            bad=M._ModelNativeHybridSizeSeeds(m,changed,receipt.records)
            @test_throws ArgumentError M._model_native_hybrid_seed_records(bad)
            constructor=component.constructor
            source=Mesh(constructor.coordinates;tets=zeros(Int32,4,0))
            @test_throws ArgumentError M.ModelNativeHybridInsertion.native_sizes(source,constructor.node_uids,[],"changed seed";minimize=true,seeds=bad)
            foreign=deepcopy(m);foreign.embeds[(3,99)]=[(0,1)]
            before=options.random_state[]
            @test_throws ArgumentError M._model_native_hybrid_connected_production(foreign,component,options,births,"foreign target")
            @test options.random_state[]==before
            @test births.node_high_water==UInt64(case["node_clock"])
            @test births.element_high_water==UInt64(case["current_lower_element_high_water"])
            if startswith(case["name"],"curve")
                @test !(UInt64(9) in production.boundary_node_tags)
                @test only(filter(record->record.uid==UInt64(1),receipt.records)).size==.05
            end
        end
    end
    @test observed["point_control"]!=observed["point_foreign_seed"]
    @test observed["curve_control"]!=observed["curve_foreign_seed"]
end
end
