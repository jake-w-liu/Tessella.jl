module NativeHybridGeneralModelChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,validate
using Tessella.Elements: ElementBlock,MixedMesh,MixedEntityData,MixedEntity
const M=Tessella.Model
const A=Tessella.API
const artifacts=joinpath(@__DIR__,"..","artifacts")
const native=TOML.parsefile(joinpath(artifacts,"native_hybrid_general_model_4152_v2.toml"))
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
helper_path=joinpath(@__DIR__,"api_native_hybrid_actual_model_test.jl")
for name in (:box_model,:actual_lower)
    found=helper(Meta.parseall(read(helper_path,String);filename=helper_path),name)
    @assert found!==nothing;Core.eval(@__MODULE__,found)
end
function foreign_face!(m)
    M.add_point!(m,-.05,0,0;tag=99,mesh_size=.4)
    M.add_point!(m,0,-.05,0;tag=100,mesh_size=.4)
    M.add_line!(m,1,99;tag=99);M.add_line!(m,99,100;tag=100);M.add_line!(m,100,1;tag=101)
    M.add_curve_loop!(m,[99,100,101];tag=99);M.add_plane_surface!(m,[99];tag=99)
    return m
end
@testset "actual Model nonquiescent hybrid and global GFace sizing" begin
    for case in native["case"]
        println("HYBRID_GENERAL_MODEL_CASE ",case["name"]);flush(stdout)
        m=box_model();case["name"]=="foreign_shared_point" && foreign_face!(m)
        input=actual_lower(m,case);saved=deepcopy(input.mesh)
        options=M._model_default_mesh1d_options(m,"native general Model")
        options.random_state[]=UInt32(case["native_phase_random_state"])
        options=M._ModelMesh1DOptions((name===:optimize ? 0 : getfield(options,name) for name in fieldnames(M._ModelMesh1DOptions))...)
        production=M.mesh_model_volume(m,1;_mesh_input=input,_mesh1d_options=options,
            _return_volume_production=true,_mesh_node_high_water=UInt64(case["node_clock"]),
            _mesh_element_high_water=UInt64(case["current_lower_element_high_water"]))
        @test production isa M._ModelNativeHybridVolumeProduction
        @test production.input===input
        @test validate(production.mesh).ok
        @test production.node_high_water==UInt64(case["node_high_water"])
        @test production.element_high_water==UInt64(case["element_high_water"])
        expected=case["expected"]
        wanted_ids=UInt64[uid for cell in expected["cell"] for uid in cell["labels"]]
        wanted_nodes=[Tuple(UInt64.(nodes)) for cell in expected["cell"] for nodes in cell["nodes"]]
        @test production.element_tags==wanted_ids
        @test [Tuple(production.boundary_node_tags[nodes]) for block in production.mesh.blocks for nodes in eachcol(block.nodes)]==wanted_nodes
        rows=Dict(uid=>row for (row,uid) in enumerate(production.boundary_node_tags))
        wanted_coords=reduce(hcat,hexpoint.(expected["xyz_hex"]))
        observed=production.mesh.coords[:,[rows[UInt64(uid)] for uid in expected["node_tags"]]]
        differences=[(axis,expected["node_tags"][row]) for row in axes(observed,2),axis in 1:3 if !isequal(observed[axis,row],wanted_coords[axis,row])]
        println("HYBRID_GENERAL_MODEL_OUTPUT ",case["name"]," families=",[(b.msh,size(b.nodes,2)) for b in production.mesh.blocks],
            " clocks=",(production.node_high_water,production.element_high_water)," coordinate_bit_differences=",length(differences));flush(stdout)
        @test isequal(observed,wanted_coords)
        source_uids=Set(input.mesh.entity_data.external_node_tags)
        registered=production.boundary_node_tags[production.registration_rows]
        @test UInt64[uid for uid in registered if !(uid in source_uids)]==UInt64.(expected["node_tags"])
        @test M._model_native_hybrid_volume_certificate(production.mesh,production.boundary_node_tags,production.constructor,"native general Model")===nothing
        @test isequal(input.mesh.coords,saved.coords)
        @test input.mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
        @test isequal(input.mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
        # Native lower Tri vectors are reordered after recovery. Exercise the
        # actual detached helper, preserving constructor/event identities.
        data=input.mesh.entity_data
        nodes=A._Dim01Node[A._Dim01Node(uid,ntuple(axis->input.mesh.coords[axis,row],3),
            data.node_entities[row],data.node_parametric[row]) for (row,uid) in enumerate(data.external_node_tags)]
        cells=A._Dim01Cell[]
        for (block,owners,labels) in zip(input.mesh.blocks,data.block_entities,data.external_element_tags),column in axes(block.nodes,2)
            owner=(Tessella.Elements.msh_dimension(block.msh),owners[column])
            push!(cells,A._Dim01Cell(Int(block.msh),owner,collect(block.nodes[:,column]),labels[column],labels[column],Int32(0)))
        end
        detached=(nodes=nodes,cells=cells,caller="native hybrid lower order")
        A._higher_native_region_surface_order!(detached,m,1)
        for product in case["post_product"],block in product["cell"]
            selected_cells=filter(cell->cell.owner==(Int(product["owner"][1]),Int32(product["owner"][2])) && cell.msh==block["msh"],cells)
            @test [cell.tag for cell in selected_cells]==UInt64.(block["labels"])
            @test [Tuple(nodes[row].tag for row in cell.nodes) for cell in selected_cells]==Tuple.(UInt64.(ids) for ids in block["nodes"])
        end
        @test isequal(input.mesh.coords,saved.coords)
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
        triangles=M._model_native_hybrid_global_triangles(m,input,"native general Model")
        @test !isempty(triangles)
        selected=M._model_native_hybrid_volume_input(m,1,input,"native general Model")
        first_triangle=first(filter(t->any(uid->uid in selected.node_uids,t.nodes),triangles))
        coordinates=collect(first_triangle.coordinates)
        slot=findfirst(uid->uid in selected.node_uids,first_triangle.nodes)
        coordinates[slot]=(nextfloat(coordinates[slot][1]),coordinates[slot][2],coordinates[slot][3])
        bad=merge(first_triangle,(coordinates=Tuple(coordinates),))
        @test_throws ArgumentError M.ModelNativeHybridInsertion.native_sizes(Mesh(selected.coordinates;tets=zeros(Int32,4,0)),
            selected.node_uids,[bad],"changed selected UID")
        if case["name"]=="foreign_shared_point"
            @test any(t->t.owner==(2,Int32(99)),triangles)
            @test !(UInt64(9) in production.boundary_node_tags) && !(UInt64(10) in production.boundary_node_tags)
            undeclared=deepcopy(m);delete!(undeclared.surfaces,99)
            state=options.random_state[]
            @test_throws ArgumentError M.mesh_model_volume(undeclared,1;_mesh_input=input,_mesh1d_options=options,
                _return_volume_production=true,_mesh_node_high_water=UInt64(case["node_clock"]),
                _mesh_element_high_water=UInt64(case["current_lower_element_high_water"]))
            @test options.random_state[]==state
            @test isequal(input.mesh.coords,saved.coords)
        end
    end
end
end
