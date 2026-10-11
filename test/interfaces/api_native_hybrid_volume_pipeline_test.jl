module NativeHybridVolumePipelineChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,ntets,validate
const M=Tessella.Model
const artifacts=joinpath(@__DIR__,"..","artifacts")
data=TOML.parsefile(joinpath(artifacts,"native_hybrid_quad_birth_4152_v2.toml"))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
@testset "actual Full Model hybrid insertion presentation and relocated births" begin
    for name in ("coarse","original")
        println("HYBRID_FULL_MODEL_CASE ",name);flush(stdout)
        native=TOML.parsefile(joinpath(artifacts,"native_hybrid_initial_"*name*"_4152_v2.toml"))
        case=only(filter(c->c["name"]==name*"_no_tet_optimizer",data["case"]))
        ids=UInt64.(case["source_uids"]);coords=reduce(hcat,hexpoint.(case["source_xyz_hex"]))
        quads=[(owner=(2,Int32(owner)),uid=UInt64(uid),nodes=Tuple(UInt64.(nodes))) for
            (owner,uid,nodes) in zip(case["owners"],case["source_quad_uids"],case["source_quad_nodes"])]
        plan=M._model_native_hybrid_quad_plan(ids,coords,quads,UInt64(case["node_clock"]),"native hybrid")
        surface=Mesh(plan.coordinates;tris=plan.fan_facets)
        order=sortperm(plan.node_uids);inverse=invperm(order)
        births=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),UInt64(native["current_lower_element_high_water"]))
        births.source=Mesh(coords);births.node_tags=copy(plan.node_uids)
        kernel=Ref{M.ModelNativeInitialTet.Kernel}()
        options=(random_state=Ref(UInt32(1)),random_factor3d=1e-12)
        M._model_native_initial_delaunay_fill(surface,order,options,"native hybrid";births,kernel_output=kernel)
        facets=[Tuple(Int32.(inverse[plan.fan_facets[:,column]])) for column in axes(plan.fan_facets,2)]
        sort!(facets;by=face->Tuple(sort!(UInt64[plan.node_uids[order[row]] for row in face])))
        interior=M.ModelNativeHybridInterior.recover(kernel[],facets)
        tets=reduce(hcat,Int32[order[row] for row in cell] for cell in interior.cells)
        mesh=Mesh(plan.coordinates;tets=tets)
        @test validate(mesh).ok
        M._model_native_volume_cell_births!(births,mesh)
        println("HYBRID_FULL_MODEL_INTERIOR ",name," count ",ntets(mesh));flush(stdout)
        refined=M.ModelNativeHybridInsertion.quiescent(mesh,births.node_tags,births.element_tags,NamedTuple[],"native hybrid")
        @test refined.node_numbers==plan.node_uids
        @test maximum(refined.radii)<1.0
        @test maximum(refined.element_numbers)==births.element_high_water
        expected=Dict(UInt64(uid)=>Tuple(UInt64.(nodes)) for (uid,nodes) in zip(native["tet_uids"],native["tet_nodes"]))
        observed=Dict(uid=>Tuple(plan.node_uids[refined.mesh.tets[:,column]]) for (column,uid) in enumerate(refined.element_numbers))
        @test observed==expected
        @test refined.element_numbers==UInt64.(native["tet_uids"])
        println("HYBRID_FULL_MODEL_INSERTION ",name," ordered_uid_vector_equal ",refined.element_numbers==UInt64.(native["tet_uids"]));flush(stdout)
        relocated=M.ModelNativeHybridRelocate.relocate(plan,refined.mesh,"native hybrid")
        wanted=reduce(hcat,hexpoint.(case["apex_xyz_hex"]))
        actual=relocated.coords[:,length(ids)+1:end]
        differences=[(axis,row,actual[axis,row],wanted[axis,row]) for row in axes(actual,2),axis in 1:3 if !isequal(actual[axis,row],wanted[axis,row])]
        println("HYBRID_FULL_MODEL_RELOCATION ",name," coordinate_bit_differences ",length(differences)," first ",isempty(differences) ? nothing : first(differences));flush(stdout)
        @test isequal(actual,wanted)
        @test isequal(relocated.coords[:,1:length(ids)],coords)
        @test relocated.tets==refined.mesh.tets
        pyramids=M._model_native_hybrid_pyramid_production(plan,relocated.coords,births.element_high_water,"native hybrid")
        @test pyramids.element_tags==UInt64.(case["pyramid_uids"])
        @test [Tuple(plan.node_uids[cell]) for cell in eachcol(pyramids.nodes)]==Tuple.(UInt64.(nodes) for nodes in case["pyramid_nodes"])
        final=M._gmsh_optimize_tetrahedra(relocated;movable=falses(nnodes(relocated)),
            element_numbers=refined.element_numbers,element_high_water=pyramids.element_high_water,threshold=0.3)
        default=TOML.parsefile(joinpath(artifacts,"native_hybrid_default_"*name*"_4152_v2.toml"))
        wanted_default=Dict(UInt64(uid)=>Tuple(UInt64.(nodes)) for (uid,nodes) in zip(default["tet_uids"],default["tet_nodes"]))
        observed_default=Dict(uid=>Tuple(plan.node_uids[final.mesh.tets[:,column]]) for (column,uid) in enumerate(final.element_numbers))
        differences_default=[uid for uid in union(keys(observed_default),keys(wanted_default)) if get(observed_default,uid,nothing)!=get(wanted_default,uid,nothing)]
        println("HYBRID_FULL_MODEL_DEFAULT ",name," count ",ntets(final.mesh)," HWM ",final.element_high_water," ordered_uid_differences ",length(differences_default));flush(stdout)
        @test observed_default==wanted_default
        @test final.element_numbers==UInt64.(default["tet_uids"])
        @test final.element_high_water==UInt64(default["element_high_water"])
        @test final.mesh.coords==relocated.coords

    end
end
end
