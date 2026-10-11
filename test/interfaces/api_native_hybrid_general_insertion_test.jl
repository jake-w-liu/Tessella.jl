module NativeHybridGeneralInsertionChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,ntets,validate
using Tessella.Elements: ElementBlock,MixedMesh,MixedEntityData,MixedEntity
const M=Tessella.Model
const artifacts=joinpath(@__DIR__,"..","artifacts")
const native=TOML.parsefile(joinpath(artifacts,"native_hybrid_general_insertion_4152_v2.toml"))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]

# Reuse only the actual owned-lower fixture constructors. The production
# testsets in that file are separate controls and are not evaluated here.
function helper(expr,name)
    expr isa Expr || return nothing
    if expr.head===:function && expr.args[1] isa Expr && expr.args[1].head===:call &&
            expr.args[1].args[1]===name
        return expr
    end
    for child in expr.args
        found=helper(child,name);found===nothing || return found
    end
    return nothing
end
fixture_helper=joinpath(@__DIR__,"api_native_hybrid_actual_model_test.jl")
for name in (:box_model,:actual_lower)
    found=helper(Meta.parseall(read(fixture_helper,String);filename=fixture_helper),name)
    @assert found!==nothing
    Core.eval(@__MODULE__,found)
end

function expected_cells(expected,family)
    cell=only(filter(c->c["msh"]==family,expected["cell"]))
    return UInt64.(cell["labels"]),Tuple.(UInt64.(ids) for ids in cell["nodes"])
end

@testset "actual native hybrid nonquiescent insertion and owned Pyramid relocation" begin
    for case in native["case"]
        println("HYBRID_GENERAL_CASE ",case["name"]);flush(stdout)
        m=box_model();input=actual_lower(m,case);saved=deepcopy(input.mesh)
        constructor=M._model_native_hybrid_volume_input(m,1,input,"native general hybrid")
        @test constructor!==nothing
        plan=M._model_native_hybrid_quad_plan(constructor.node_uids,constructor.coordinates,
            constructor.quads,UInt64(case["node_clock"]),"native general hybrid")
        initial_positions=Dict(UInt64(uid)=>row for (row,uid) in enumerate(case["initial_node_tags"]))
        actual_initial=reduce(hcat,hexpoint.(case["initial_xyz_hex"]))
        # Algorithm3 retains the flat Quad apex only through Tet incidence;
        # native getNodes/getNode does not expose those unclassified objects.
        @test Set(plan.source_uids)==Set(keys(initial_positions))
        @test isequal(plan.source_coordinates,actual_initial[:,[initial_positions[uid] for uid in plan.source_uids]])
        @test Set(UInt64(uid) for cell in case["initial_tet_nodes"] for uid in cell)==Set(plan.node_uids)
        positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(plan.node_uids))
        cells=reduce(hcat,Int32[positions[UInt64(uid)] for uid in cell] for cell in case["initial_tet_nodes"])
        initial=Mesh(copy(plan.coordinates);tets=cells)
        @test validate(initial).ok
        options=M._model_default_mesh1d_options(m,"native general hybrid")
        # Qualify fresh construction independently before consuming the SDK
        # initial cells in the isolated insertion control below.
        options.random_state[]=UInt32(1)
        fan=copy(plan.fan_facets)
        for triangle in constructor.triangles;fan=hcat(fan,collect(triangle.rows));end
        order=sortperm(plan.node_uids)
        births=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),UInt64(case["current_lower_element_high_water"]))
        births.source=Mesh(plan.source_coordinates);births.node_tags=copy(plan.node_uids)
        kernel=Ref{M.ModelNativeInitialTet.Kernel}()
        M._model_native_initial_delaunay_fill(Mesh(plan.coordinates;tris=fan),order,options,"native general hybrid";
            births,kernel_output=kernel)
        facets=M._model_native_hybrid_recovery_facets(plan,constructor,order,"native general hybrid")
        interior=M.ModelNativeHybridInterior.recover(kernel[],facets)
        fresh=Mesh(plan.coordinates;tets=reduce(hcat,Int32[order[row] for row in cell] for cell in interior.cells))
        M._model_native_volume_cell_births!(births,fresh)
        @test births.node_high_water==UInt64(case["initial_node_high_water"])
        @test births.element_high_water==UInt64(case["initial_element_high_water"])
        @test births.element_tags==UInt64.(case["initial_tet_uids"])
        @test [Tuple(plan.node_uids[cell]) for cell in eachcol(fresh.tets)]==Tuple.(UInt64.(nodes) for nodes in case["initial_tet_nodes"])
        println("HYBRID_GENERAL_FRESH_INITIAL ",case["name"]," native_ordered_cells=",ntets(fresh));flush(stdout)
        background=M._ModelAlgorithmBackgroundSize(m,3,1,options,nothing,"native general hybrid")
        refined=M.ModelNativeHybridInsertion.refine(initial,plan.node_uids,
            UInt64.(case["initial_tet_uids"]),constructor.triangles,background,"native general hybrid";
            node_high_water=UInt64(case["initial_node_high_water"]),
            element_high_water=UInt64(case["initial_element_high_water"]))
        println("HYBRID_GENERAL_INSERTION ",case["name"]," rows=",nnodes(refined.mesh),
            " cells=",ntets(refined.mesh)," clocks=",(refined.node_high_water,refined.element_high_water));flush(stdout)
        wanted_tet_uids,wanted_tets=expected_cells(case["expected"],4)
        @test refined.element_numbers==wanted_tet_uids
        @test [Tuple(refined.node_numbers[cell]) for cell in eachcol(refined.mesh.tets)]==wanted_tets
        @test refined.node_high_water==UInt64(case["node_high_water"])
        @test refined.element_high_water+UInt64(length(plan.births))==UInt64(case["element_high_water"])
        @test isequal(refined.mesh.coords[:,1:length(plan.node_uids)],plan.coordinates)
        @test refined.node_numbers[1:length(plan.node_uids)]==plan.node_uids
        for (key,passes) in (("no_smoothing",0),("expected",m.meshing.native_smoothing))
            smoothed=M.ModelNativeHybridInsertion.smooth(refined,length(plan.node_uids),passes,"native general hybrid")
            relocated=M.ModelNativeHybridRelocate.relocate(plan,smoothed.mesh,"native general hybrid";
                birth_order=constructor.pyramid_order,node_uids=refined.node_numbers)
            pyramids=M._model_native_hybrid_pyramid_production(plan,relocated.coords,
                refined.element_high_water,"native general hybrid";birth_order=constructor.pyramid_order)
            expected=case[key];wanted_pyramid_uids,wanted_pyramids=expected_cells(expected,7)
            @test pyramids.element_tags==wanted_pyramid_uids
            @test [Tuple(refined.node_numbers[cell]) for cell in eachcol(pyramids.nodes)]==wanted_pyramids
            wanted_positions=Dict(UInt64(uid)=>row for (row,uid) in enumerate(expected["node_tags"]))
            wanted_coords=reduce(hcat,hexpoint.(expected["xyz_hex"]))
            new_ids=refined.node_numbers[length(plan.source_uids)+1:end]
            @test Set(new_ids)==Set(keys(wanted_positions))
            actual=relocated.coords[:,length(plan.source_uids)+1:end]
            wanted=wanted_coords[:,[wanted_positions[uid] for uid in new_ids]]
            differences=[(axis,new_ids[row]) for row in axes(actual,2),axis in 1:3 if !isequal(actual[axis,row],wanted[axis,row])]
            println("HYBRID_GENERAL_RELOCATION ",case["name"]," ",key,
                " coordinate_bit_differences=",length(differences)," first=",isempty(differences) ? nothing : first(differences));flush(stdout)
            @test isequal(actual,wanted)
            @test isequal(relocated.coords[:,1:length(plan.source_uids)],plan.source_coordinates)
            mixed=MixedMesh(relocated.coords,[ElementBlock(4,relocated.tets),ElementBlock(7,pyramids.nodes)])
            @test validate(mixed).ok
            @test M._model_native_hybrid_volume_certificate(mixed,refined.node_numbers,constructor,"native general hybrid")===nothing
        end
        @test isequal(input.mesh.coords,saved.coords)
        @test input.mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
        @test isequal(input.mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
    end
end
end
