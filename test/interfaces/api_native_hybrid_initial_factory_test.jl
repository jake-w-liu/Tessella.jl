module NativeHybridInitialFactoryChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,ntets
const M=Tessella.Model
data=TOML.parsefile(joinpath(@__DIR__,"..","artifacts/native_hybrid_quad_birth_4152_v2.toml"))
native=TOML.parsefile(joinpath(@__DIR__,"..","artifacts/native_hybrid_initial_coarse_4152_v2.toml"))
case=only(filter(c->c["name"]=="coarse_no_tet_optimizer",data["case"]))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
ids=UInt64.(case["source_uids"]);coords=reduce(hcat,hexpoint.(case["source_xyz_hex"]))
quads=[(owner=(2,Int32(owner)),uid=UInt64(uid),nodes=Tuple(UInt64.(nodes))) for
    (owner,uid,nodes) in zip(case["owners"],case["source_quad_uids"],case["source_quad_nodes"])]
plan=M._model_native_hybrid_quad_plan(ids,coords,quads,UInt64(case["node_clock"]),"native initial")
surface=Mesh(plan.coordinates;tris=plan.fan_facets)
births=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),UInt64(native["current_lower_element_high_water"]))
births.source=Mesh(coords);births.node_tags=copy(plan.node_uids)
kernel=Ref{M.ModelNativeInitialTet.Kernel}()
mesh=M._model_native_initial_delaunay_fill(surface,sortperm(plan.node_uids),nothing,"native initial";
    births,kernel_output=kernel)
M._model_native_volume_cell_births!(births,mesh)
wanted=Dict(UInt64(uid)=>Tuple(UInt64.(nodes)) for (uid,nodes) in zip(native["tet_uids"],native["tet_nodes"]))
observed=Dict(uid=>Tuple(plan.node_uids[mesh.tets[:,column]]) for (column,uid) in enumerate(births.element_tags))
@testset "native coarse hybrid initial Tet factories retain exact births" begin
    @test ntets(mesh)==24
    @test births.temporary_elements==88
    @test births.node_high_water==22
    @test births.element_high_water==156
    @test observed==wanted
    @test mesh.coords==plan.coordinates
    @test births.node_tags==plan.node_uids
    @test !M.ModelNativeTetgenArena.has_bad_quality(kernel[],M.ModelNativeTetgenArena.QualityWorkspace())
end
end
