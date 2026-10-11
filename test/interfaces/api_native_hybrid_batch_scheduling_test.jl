module NativeHybridBatchSchedulingChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,validate
using Tessella.Elements: ElementBlock,MixedMesh,MixedEntityData,MixedEntity
const M=Tessella.Model
const artifacts=joinpath(@__DIR__,"..","artifacts")
const native=TOML.parsefile(joinpath(artifacts,"native_hybrid_batch_scheduling_4152_v1.toml"))
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
body=replace(read(path,String),"for tag in 1:6"=>"for tag in sort!(collect(keys(m.surfaces)))")
for name in (:box_model,:actual_lower)
    found=helper(Meta.parseall(body;filename=path),name)
    @assert found!==nothing;Core.eval(@__MODULE__,found)
end
path=joinpath(@__DIR__,"api_native_hybrid_connected_components_test.jl")
Core.eval(@__MODULE__,helper(Meta.parseall(read(path,String);filename=path),:adjacent_model))
function disjoint_model()
    m=box_model()
    points=((3,0,0),(4,0,0),(4,1,0),(3,1,0),(3,0,1),(4,0,1),(4,1,1),(3,1,1))
    for (tag,p) in zip(101:108,points);M.add_point!(m,p...;tag,mesh_size=.4);end
    curves=((101,102),(102,103),(103,104),(104,101),(105,106),(106,107),(107,108),(108,105),(101,105),(102,106),(103,107),(104,108))
    for (tag,ends) in zip(101:112,curves);M.add_line!(m,ends...;tag);end
    loops=((101,102,103,104),(105,106,107,108),(101,110,-105,-109),(102,111,-106,-110),(103,112,-107,-111),(104,109,-108,-112))
    for (tag,curves) in zip(101:106,loops)
        M.add_curve_loop!(m,collect(curves);tag);M.add_plane_surface!(m,[tag];tag)
    end
    M.add_surface_loop!(m,collect(101:106);tag=101);M.add_volume!(m,[101];tag=2)
    return m
end
@testset "native connected and disjoint hybrid birth scheduling before default optimization" begin
    results=Dict{String,Any}()
    for case in native["case"]
        println("HYBRID_BATCH_CASE ",case["name"]);flush(stdout)
        m=startswith(case["name"],"connected") ? adjacent_model() : disjoint_model()
        input=actual_lower(m,case);saved=deepcopy(input.mesh)
        options=M._model_default_mesh1d_options(m,"native global batch")
        options=M._ModelMesh1DOptions((name===:optimize ? Int(case["optimize"]) : getfield(options,name)
            for name in fieldnames(M._ModelMesh1DOptions))...)
        # Native Generator resets its RNG at this call. The producer receives
        # and carries the same Ref through every component without reseeding.
        options.random_state[]=UInt32(1)
        births=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),UInt64(case["current_lower_element_high_water"]))
        before=options.random_state[]
        @test_throws ArgumentError M._model_native_hybrid_batch_production(m,[1,1],input,options,births,"duplicate Region")
        @test_throws ArgumentError M._model_native_hybrid_batch_production(m,[1,99],input,options,births,"foreign Region")
        missing=deepcopy(input);delete!(missing.products,(2,Int32(1)))
        @test_throws ArgumentError M._model_native_hybrid_batch_production(m,[1,2],missing,options,births,"missing lower Face")
        changed=deepcopy(input);changed.mesh.coords[1,1]=nextfloat(changed.mesh.coords[1,1])
        @test_throws ArgumentError M._model_native_hybrid_batch_production(m,[1,2],changed,options,births,"one ULP Source")
        @test options.random_state[]==before
        @test births.node_high_water==UInt64(case["node_clock"])
        @test births.element_high_water==UInt64(case["current_lower_element_high_water"])
        batch=M._model_native_hybrid_batch_production(m,[2,1],input,options,births,"native global batch")
        @test batch isa M._ModelNativeHybridBatchProduction
        @test batch.input===input
        @test batch.regions==[1,2]
        @test batch.groups==case["groups"]
        @test batch.random_state==options.random_state[]
        @test batch.random_state!=before
        @test batch.node_high_water==UInt64(case["node_high_water"])
        @test batch.element_high_water==UInt64(case["element_high_water"])
        @test Set(keys(batch.products))==Set([1,2])
        native_owners=Dict{UInt64,Int}()
        for expected in case["expected"],uid in expected["node_tags"]
            @assert !haskey(native_owners,UInt64(uid));native_owners[UInt64(uid)]=Int(expected["owner"][2])
        end
        owned=Dict{UInt64,Int}()
        for expected in case["expected"]
            tag=Int(expected["owner"][2]);production=batch.products[tag]
            @test production.input===input
            @test validate(production.mesh).ok
            @test production.node_high_water==batch.node_high_water
            @test production.element_high_water==batch.element_high_water
            ids=UInt64[uid for cell in expected["cell"] for uid in cell["labels"]]
            corners=Tuple.(UInt64.(nodes) for cell in expected["cell"] for nodes in cell["nodes"])
            @test production.element_tags==ids
            @test [Tuple(production.boundary_node_tags[rows]) for block in production.mesh.blocks for rows in eachcol(block.nodes)]==corners
            @test allunique(production.boundary_node_tags)
            positions=Dict(uid=>row for (row,uid) in enumerate(production.boundary_node_tags))
            observed=production.mesh.coords[:,[positions[UInt64(uid)] for uid in expected["node_tags"]]]
            wanted=reduce(hcat,hexpoint.(expected["xyz_hex"]))
            differences=[(axis,expected["node_tags"][row]) for row in axes(observed,2),axis in 1:3 if !isequal(observed[axis,row],wanted[axis,row])]
            println("HYBRID_BATCH_FINAL ",case["name"]," Region=",tag," families=",[(block.msh,size(block.nodes,2)) for block in production.mesh.blocks],
                " clocks=",(production.node_high_water,production.element_high_water)," coordinate_bit_differences=",length(differences));flush(stdout)
            @test isequal(observed,wanted)
            @test M._model_native_hybrid_volume_certificate(production.mesh,production.boundary_node_tags,production.constructor,"native global batch")===nothing
            source=Set(input.mesh.entity_data.external_node_tags)
            registered=production.boundary_node_tags[production.registration_rows]
            @test UInt64[uid for uid in registered if !(uid in source)]==UInt64.(expected["node_tags"])
            newids=UInt64[uid for uid in production.boundary_node_tags if !(uid in source)]
            @test Set(newids)==Set(UInt64.(expected["node_tags"]))
            for uid in newids
                @assert !haskey(owned,uid);owned[uid]=tag
            end
        end
        @test owned==native_owners
        @test allunique(vcat((product.element_tags for product in values(batch.products))...))
        @test isequal(input.mesh.coords,saved.coords)
        @test isequal(input.mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
        @test input.mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
        @test input.mesh.entity_data.external_element_tags==saved.entity_data.external_element_tags
        results[case["name"]]=batch
    end
    # The default pass changes surviving Tet births after both components;
    # Pyramid factory identities remain those of the unoptimized operation.
    plain=results["disjoint_original_opt0"];optimized=results["disjoint_original_default"]
    @test optimized.node_high_water==plain.node_high_water
    @test optimized.element_high_water>plain.element_high_water
    for tag in (1,2)
        pyramid_ids(product)=product.element_tags[sum(size(block.nodes,2) for block in product.mesh.blocks if block.msh==4)+1:end]
        @test pyramid_ids(plain.products[tag])==pyramid_ids(optimized.products[tag])
        @test plain.products[tag].element_tags!=optimized.products[tag].element_tags
    end
end
end
