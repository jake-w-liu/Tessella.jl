module NativeHybridConnectedComponentChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,ntets,validate
using Tessella.Elements: ElementBlock,MixedMesh,MixedEntityData,MixedEntity
const M=Tessella.Model
const artifacts=joinpath(@__DIR__,"..","artifacts")
const native=TOML.parsefile(joinpath(artifacts,"native_hybrid_connected_components_4152_v2.toml"))
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
function adjacent_model()
    m=box_model()
    for (tag,p) in zip(101:104,((2,0,0),(2,1,0),(2,0,1),(2,1,1)))
        M.add_point!(m,p...;tag,mesh_size=.4)
    end
    for (tag,ends) in zip(101:108,((2,101),(101,102),(102,3),(6,103),(103,104),(104,7),(101,103),(102,104)))
        M.add_line!(m,ends...;tag)
    end
    for (tag,curves) in zip(101:105,((101,102,103,-2),(104,105,106,-6),(101,107,-104,-10),(102,108,-105,-107),(103,11,-106,-108)))
        M.add_curve_loop!(m,collect(curves);tag);M.add_plane_surface!(m,[tag];tag)
    end
    M.add_surface_loop!(m,[101,102,103,104,105,-4];tag=101)
    M.add_volume!(m,[101];tag=2)
    return m
end
@testset "connected native hybrid initial recovery and exact declared Region partition" begin
    for case in native["case"]
        println("HYBRID_CONNECTED_CASE ",case["name"]);flush(stdout)
        m=adjacent_model();input=actual_lower(m,case);saved=deepcopy(input.mesh)
        component=M._model_native_hybrid_component_receipt(m,[2,1],input,"native connected")
        @test component.regions==[1,2]
        @test component.face_regions[4]==[1,2]
        @test all(length(owners)==1 for (face,owners) in component.face_regions if face!=4)
        constructor=component.constructor
        plan=M._model_native_hybrid_quad_plan(constructor.node_uids,constructor.coordinates,
            constructor.quads,UInt64(case["node_clock"]),"native connected")
        work=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),UInt64(case["current_lower_element_high_water"]))
        work.source=Mesh(plan.source_coordinates);work.node_tags=copy(plan.node_uids)
        work.node_high_water=plan.node_high_water
        options=M._model_default_mesh1d_options(m,"native connected")
        options.random_state[]=UInt32(1)
        order=sortperm(plan.node_uids);kernel=Ref{M.ModelNativeInitialTet.Kernel}()
        surface=Mesh(plan.coordinates;tris=plan.fan_facets)
        M._model_native_initial_delaunay_fill(surface,order,options,"native connected";births=work,kernel_output=kernel)
        facets=M._model_native_hybrid_recovery_facets(plan,constructor,order,"native connected")
        certificate=M._model_native_hybrid_component_recovery_receipt(component,plan,order,kernel[],"native connected")
        @test certificate.kernel===kernel[]
        @test count(==(2),values(certificate.multiplicities))==length(constructor.triangles)
        @test_throws M.ModelNativeTetgenArena.Unsupported M.ModelNativeHybridInterior.recover(kernel[],facets)
        @test_throws ArgumentError M.ModelNativeHybridInterior.recover(deepcopy(kernel[]),facets;component=certificate)
        changed=copy(facets);reverse!(changed)
        @test_throws ArgumentError M.ModelNativeHybridInterior.recover(kernel[],changed;component=certificate)
        interior=M.ModelNativeHybridInterior.recover(kernel[],facets;component=certificate)
        mesh=Mesh(plan.coordinates;tets=reduce(hcat,Int32[order[row] for row in cell] for cell in interior.cells))
        M._model_native_volume_cell_births!(work,mesh)
        @test validate(mesh).ok
        @test work.node_high_water==UInt64(case["initial_node_high_water"])
        @test work.element_high_water==UInt64(case["initial_element_high_water"])
        expected=only(filter(product->product["owner"]==[3,1],case["initial"]))
        expected_cells=only(expected["cell"])
        @test work.element_tags==UInt64.(expected_cells["labels"])
        observed=[Tuple(plan.node_uids[rows]) for rows in eachcol(mesh.tets)]
        @test observed==Tuple.(UInt64.(ids) for ids in expected_cells["nodes"])
        owners=M._model_native_hybrid_component_partition(mesh,plan.node_uids,work.element_tags,component,plan,"native connected")
        @test [count(==(tag),owners) for tag in component.regions]==case["expected_region_tets"]
        # Coordinates verify the already selected exact Face-incidence owner.
        # They never choose a UID, panel, or Region for the production.
        @test all((sum(mesh.coords[1,mesh.tets[:,column]])/4<1)==(owners[column]==1) for column in axes(mesh.tets,2))
        changed_mesh=Mesh(copy(mesh.coords);tets=copy(mesh.tets))
        changed_mesh.coords[1,1]=nextfloat(changed_mesh.coords[1,1])
        @test_throws ArgumentError M._model_native_hybrid_component_partition(changed_mesh,plan.node_uids,work.element_tags,component,plan,"changed Source")
        @test_throws ArgumentError M._model_native_hybrid_component_partition(mesh,plan.node_uids[1:end-1],work.element_tags,component,plan,"missing UID")
        @test_throws ArgumentError M._model_native_hybrid_component_partition(mesh,plan.node_uids,fill(first(work.element_tags),ntets(mesh)),component,plan,"duplicate cell UID")
        repeated=Mesh(mesh.coords;tets=hcat(mesh.tets,mesh.tets[:,1]))
        @test_throws ArgumentError M._model_native_hybrid_component_partition(repeated,plan.node_uids,vcat(work.element_tags,last(work.element_tags)+1),component,plan,"duplicate native cell")
        @test_throws ArgumentError M._model_native_hybrid_component_receipt(m,[1,1],input,"duplicate Region")
        @test_throws ArgumentError M._model_native_hybrid_component_receipt(m,[1,99],input,"foreign Region")
        ambiguous=deepcopy(m);ambiguous.volumes[2]=copy(ambiguous.volumes[1])
        @test_throws ArgumentError M._model_native_hybrid_component_receipt(ambiguous,[1,2],input,"ambiguous Region")
        changed_input=deepcopy(input)
        first_row=findfirst(==(first(constructor.node_uids)),changed_input.mesh.entity_data.external_node_tags)
        changed_input.mesh.coords[1,first_row]=nextfloat(changed_input.mesh.coords[1,first_row])
        changed_component=M._ModelNativeHybridComponentReceipt(changed_input,component.regions,
            component.boundaries,component.face_regions,merge(constructor,(input=changed_input,)))
        @test_throws ArgumentError M._model_native_hybrid_component_recovery_receipt(changed_component,plan,order,kernel[],"one ULP Source")
        @test isequal(input.mesh.coords,saved.coords)
        @test isequal(input.mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
        @test input.mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
        println("HYBRID_CONNECTED_INITIAL ",case["name"]," cells=",ntets(mesh)," Regions=",[count(==(tag),owners) for tag in component.regions]);flush(stdout)
        options=M._ModelMesh1DOptions((name===:optimize ? 0 : getfield(options,name)
            for name in fieldnames(M._ModelMesh1DOptions))...)
        options.random_state[]=UInt32(1)
        births=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),UInt64(case["current_lower_element_high_water"]))
        before=options.random_state[]
        @test_throws ArgumentError M._model_native_hybrid_connected_production(m,changed_component,options,births,"changed batch Source")
        @test options.random_state[]==before
        @test births.node_high_water==UInt64(case["node_clock"])
        @test births.element_high_water==UInt64(case["current_lower_element_high_water"])
        batch=M._model_native_hybrid_connected_production(m,component,options,births,"native connected final")
        @test batch.input===input
        @test Set(keys(batch.products))==Set([1,2])
        @test batch.node_high_water==UInt64(case["final_node_high_water"])
        @test batch.element_high_water==UInt64(case["final_element_high_water"])
        for expected in case["expected"]
            tag=Int(expected["owner"][2]);production=batch.products[tag]
            @test production.input===input
            @test validate(production.mesh).ok
            ids=UInt64[uid for cell in expected["cell"] for uid in cell["labels"]]
            cells=Tuple.(UInt64.(ids) for cell in expected["cell"] for ids in cell["nodes"])
            @test production.element_tags==ids
            @test [Tuple(production.boundary_node_tags[rows]) for block in production.mesh.blocks for rows in eachcol(block.nodes)]==cells
            positions=Dict(uid=>row for (row,uid) in enumerate(production.boundary_node_tags))
            observed=production.mesh.coords[:,[positions[UInt64(uid)] for uid in expected["node_tags"]]]
            wanted=reduce(hcat,hexpoint.(expected["xyz_hex"]))
            differences=[(axis,expected["node_tags"][row]) for row in axes(observed,2),axis in 1:3 if !isequal(observed[axis,row],wanted[axis,row])]
            println("HYBRID_CONNECTED_FINAL ",case["name"]," Region=",tag," families=",[(block.msh,size(block.nodes,2)) for block in production.mesh.blocks],
                " clocks=",(production.node_high_water,production.element_high_water)," coordinate_bit_differences=",length(differences));flush(stdout)
            @test isequal(observed,wanted)
            source_uids=Set(input.mesh.entity_data.external_node_tags)
            registered=production.boundary_node_tags[production.registration_rows]
            @test UInt64[uid for uid in registered if !(uid in source_uids)]==UInt64.(expected["node_tags"])
            @test M._model_native_hybrid_volume_certificate(production.mesh,production.boundary_node_tags,production.constructor,"native connected final")===nothing
        end
        @test isequal(input.mesh.coords,saved.coords)
        @test isequal(input.mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
        @test input.mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
    end
end
end
