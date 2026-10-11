module NativeHybridInteriorChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,ntets
const M=Tessella.Model
const artifacts=joinpath(@__DIR__,"..","artifacts")
data=TOML.parsefile(joinpath(artifacts,"native_hybrid_quad_birth_4152_v2.toml"))
const primary="c4d9dbba6c819623781f53be4d532f1b4167847dba40101d7e76da1f522e222e"
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]

@testset "native hybrid recovery preserves complete ordered Tet birth vectors" begin
    @test bytes2hex(sha256(read(joinpath(artifacts,"native_hybrid_quad_birth_4152_v2.json"))))==primary
    for name in ("coarse","original")
        native=TOML.parsefile(joinpath(artifacts,"native_hybrid_initial_"*name*"_4152_v2.toml"))
        case=only(filter(c->c["name"]==name*"_no_tet_optimizer",data["case"]))
        ids=UInt64.(case["source_uids"])
        coords=reduce(hcat,hexpoint.(case["source_xyz_hex"]))
        quads=[(owner=(2,Int32(owner)),uid=UInt64(uid),nodes=Tuple(UInt64.(nodes))) for
            (owner,uid,nodes) in zip(case["owners"],case["source_quad_uids"],case["source_quad_nodes"])]
        plan=M._model_native_hybrid_quad_plan(ids,coords,quads,UInt64(case["node_clock"]),"native hybrid")
        surface=Mesh(plan.coordinates;tris=plan.fan_facets)
        frozen=copy(surface.coords)
        order=sortperm(plan.node_uids);inverse=invperm(order)
        births=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),
            UInt64(native["current_lower_element_high_water"]))
        births.source=Mesh(coords);births.node_tags=copy(plan.node_uids)
        kernel=Ref{M.ModelNativeInitialTet.Kernel}()
        options=(random_state=Ref(UInt32(1)),random_factor3d=1e-12)
        M._model_native_initial_delaunay_fill(surface,order,options,"native hybrid";births,kernel_output=kernel)
        facets=[Tuple(Int32.(inverse[plan.fan_facets[:,column]])) for column in axes(plan.fan_facets,2)]
        # splitQuadRecovery::getTri traverses its actual MFaceLessThan map.
        sort!(facets;by=face->Tuple(sort!(UInt64[plan.node_uids[order[row]] for row in face])))
        before=[(cell.vertices,cell.neighbors) for cell in kernel[].cells]
        product=M.ModelNativeHybridInterior.recover(kernel[],facets)
        mapped=reduce(hcat,Int32[order[row] for row in cell] for cell in product.cells)
        mesh=Mesh(plan.coordinates;tets=mapped)
        M._model_native_volume_cell_births!(births,mesh)
        expected=Dict(UInt64(uid)=>Tuple(UInt64.(nodes)) for
            (uid,nodes) in zip(native["tet_uids"],native["tet_nodes"]))
        observed=Dict(uid=>Tuple(plan.node_uids[mesh.tets[:,column]]) for
            (column,uid) in enumerate(births.element_tags))
        @test native["primary_sha256"]==primary
        @test Set(keys(observed))==Set(keys(expected))
        for uid in sort!(collect(keys(expected)))
            @test observed[uid]==expected[uid]
        end
        @test births.temporary_elements==(name=="coarse" ? 88 : 1392)
        @test births.node_high_water==(name=="coarse" ? 22 : 270)
        @test births.element_high_water==maximum(native["tet_uids"])
        @test surface.coords==frozen
        @test births.node_tags==plan.node_uids
        @test options.random_state[]==kernel[].random_state
        @test [(cell.vertices,cell.neighbors) for cell in kernel[].cells]==before
        @test issorted(product.pool_slots)
        @test M.ModelNativeHybridUnflip.metric(product.arena,first(product.cells))==0.0
        @test_throws M.ModelNativeTetgenArena.Unsupported M.ModelNativeHybridInterior.recover(kernel[],facets[1:end-1])
        @test_throws M.ModelNativeTetgenArena.Unsupported M.ModelNativeHybridInterior.recover(kernel[],vcat(facets,facets[1:1]))
        invalid=copy(facets);invalid[1]=(invalid[1][1],invalid[1][2],Int32(kernel[].count+1))
        @test_throws M.ModelNativeTetgenArena.Unsupported M.ModelNativeHybridInterior.recover(kernel[],invalid)
    end
end
end
