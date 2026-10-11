module NativeHybridQuadBirthChecks
using Test,Tessella,TOML,SHA
const M=Tessella.Model
const TRANSPORT=joinpath(@__DIR__,"..","artifacts","native_hybrid_quad_birth_4152_v2.toml")
const PRIMARY=joinpath(@__DIR__,"..","artifacts","native_hybrid_quad_birth_4152_v2.json")
data=TOML.parsefile(TRANSPORT)
@assert data["primary_sha256"]==bytes2hex(sha256(read(PRIMARY)))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
@testset "native hybrid Quad birth and ordered Pyramid factory" begin
    for case in data["case"]
        @testset "$(case["name"])" begin
            ids=UInt64.(case["source_uids"])
            coordinates=reduce(hcat,hexpoint.(case["source_xyz_hex"]))
            quads=[(owner=(2,Int32(owner)),uid=UInt64(uid),nodes=Tuple(UInt64.(nodes))) for
                (owner,uid,nodes) in zip(case["owners"],case["source_quad_uids"],case["source_quad_nodes"])]
            snapshot=(copy(ids),copy(coordinates),deepcopy(quads))
            plan=M._model_native_hybrid_quad_plan(ids,coordinates,quads,UInt64(case["node_clock"]),"native hybrid")
            @test [birth.apex_uid for birth in plan.births]==UInt64.(case["apex_uids"])
            @test [birth.element_uid for birth in plan.births]==UInt64.(case["source_quad_uids"])
            @test plan.node_high_water==UInt64(case["node_clock"]+length(quads))
            @test plan.node_uids[1:length(ids)]==ids
            @test plan.coordinates[:,1:length(ids)]==coordinates
            @test size(plan.fan_facets)==(3,4length(quads))
            for (ordinal,birth) in enumerate(plan.births)
                @test plan.node_uids[birth.apex_row]==birth.apex_uid
                @test plan.source_uids[collect(birth.base_rows)]==collect(birth.base_uids)
                @test plan.fan_facets[:,4(ordinal-1)+1:4ordinal]==
                    Int32[birth.base_rows[1] birth.base_rows[2] birth.base_rows[3] birth.base_rows[4];
                          birth.base_rows[2] birth.base_rows[3] birth.base_rows[4] birth.base_rows[1];
                          birth.apex_row birth.apex_row birth.apex_row birth.apex_row]
            end
            final=copy(plan.coordinates)
            final[:,length(ids)+1:end].=reduce(hcat,hexpoint.(case["apex_xyz_hex"]))
            nodes=M._model_native_hybrid_pyramid_nodes(plan,final,"native hybrid")
            expected=permutedims(reduce(hcat,UInt64.(cell) for cell in case["pyramid_nodes"]))
            @test permutedims(plan.node_uids[nodes])==expected
            clock=minimum(UInt64.(case["pyramid_uids"]))-UInt64(1)
            production=M._model_native_hybrid_pyramid_production(plan,final,clock,"native Pyramid births")
            @test production.nodes==nodes
            @test production.element_tags==UInt64.(case["pyramid_uids"])
            @test production.element_high_water==maximum(UInt64.(case["pyramid_uids"]))
            @test production.quad_uids==UInt64.(case["source_quad_uids"])
            @test_throws OverflowError M._model_native_hybrid_pyramid_production(plan,final,typemax(UInt64),"overflow")
            order=reverse(collect(eachindex(plan.births)))
            @test M._model_native_hybrid_pyramid_nodes(plan,final,"restored face vector";birth_order=order)==nodes[:,order]
            @test_throws ArgumentError M._model_native_hybrid_pyramid_nodes(plan,final,"missing Quad";birth_order=order[2:end])
            duplicate_order=copy(order);duplicate_order[end]=duplicate_order[1]
            @test_throws ArgumentError M._model_native_hybrid_pyramid_nodes(plan,final,"duplicate Quad";birth_order=duplicate_order)
            @test snapshot==(ids,coordinates,quads)
            @test_throws ArgumentError M._model_native_hybrid_pyramid_nodes(plan,plan.coordinates,"flat apex")
            changed=copy(final);changed[1,1]=nextfloat(changed[1,1])
            @test_throws ArgumentError M._model_native_hybrid_pyramid_nodes(plan,changed,"Source ULP")
            duplicate=copy(quads);push!(duplicate,first(quads))
            @test_throws ArgumentError M._model_native_hybrid_quad_plan(ids,coordinates,duplicate,UInt64(case["node_clock"]),"duplicate")
            foreign=copy(quads);foreign[1]=merge(first(quads),(owner=(3,Int32(1)),))
            @test_throws ArgumentError M._model_native_hybrid_quad_plan(ids,coordinates,foreign,UInt64(case["node_clock"]),"foreign")
            missing=copy(quads);missing[1]=merge(first(quads),(nodes=(UInt64(999999),first(quads).nodes[2:end]...),))
            @test_throws ArgumentError M._model_native_hybrid_quad_plan(ids,coordinates,missing,UInt64(case["node_clock"]),"missing")
            @test snapshot==(ids,coordinates,quads)
        end
    end
end
end
