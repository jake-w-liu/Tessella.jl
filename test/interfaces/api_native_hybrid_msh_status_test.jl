module NativeHybridMshStatusChecks
using Test,TOML,SHA,Tessella
const E=Tessella.Elements
const M=Tessella.Model
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const native=TOML.parsefile(joinpath(ROOT,"test/artifacts/native_hybrid_msh_status_4152_v6.toml"))
hexpoint(values)=Float64[Core.eval(@__MODULE__,Meta.parse(value)) for value in values]
function box_model()
    m=M.GeoModel()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(points);M.add_point!(m,p...;tag,mesh_size=.4);end
    curves=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(curves);M.add_line!(m,ends...;tag);end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,curves) in enumerate(loops)
        M.add_curve_loop!(m,collect(curves);tag);M.add_plane_surface!(m,[tag];tag)
    end
    M.add_surface_loop!(m,collect(1:6);tag=1);M.add_volume!(m,[1];tag=1)
    return m
end
@testset "owned MSH factory separates lower completion from imported presence" begin
    @test bytes2hex(sha256(read(joinpath(ROOT,native["primary"]))))==native["primary_sha256"]
    @test length(native["case"])==16
    @test M._model_mesh_cleared_authority(nothing)===nothing
    for case in native["case"]
        m=box_model();history=case["history"]
        # These are explicit preceding native stage states, rather than states
        # inferred from the imported cells: Generate1 resets Faces to pending.
        done=Set{Tuple{Int,Int32}}()
        if history!="cold"
            union!(done,((1,Int32(tag)) for tag in 1:12))
            history=="generate2_clear" && union!(done,((2,Int32(tag)) for tag in 1:6))
        end
        previous=M._ModelMeshAuthority(nothing,Dict{Tuple{Int,Int32},M._ModelMeshProduct}(),done,UInt64[])
        cleared=M._model_mesh_cleared_authority(previous)
        @test cleared!==previous && cleared.done==done && cleared.done!==done
        @test cleared.mesh===nothing && isempty(cleared.products) && isempty(cleared.native_inputs)
        @test isempty(cleared.surface_constructor_inputs) && isempty(cleared.producer_rows)
        @test cleared.emitted_geometry===nothing && M._model_mesh_status(m,cleared)==-1
        input=E._mixed_msh_factory_context(case["initial_node_max"],case["initial_element_max"])
        path=joinpath(ROOT,"test/data/native_hybrid_import_events_4152_v1",case["variant"]*".msh")
        if case["collision"]
            @test_throws ArgumentError E._read_mixed_msh_factory_with_identities(path,input)
            @test (input.node_max,input.element_max)==(UInt64(case["initial_node_max"]),UInt64(case["initial_element_max"]))
            @test cleared.done==done && cleared.mesh===nothing
            continue
        end
        mesh,receipt=E._read_mixed_msh_factory_with_identities(path,input)
        plan=M._model_mesh_factory_authority(m,mesh,receipt,cleared)
        authority=plan.authority;data=authority.mesh.entity_data
        @test plan.factory===receipt.factory_events && authority.mesh===mesh
        @test_throws ArgumentError M._model_mesh_factory_authority(m,mesh,receipt,authority)
        @test authority.done==done && authority.done!==cleared.done && cleared.done==done
        @test M._model_mesh_status(m,authority)==case["status"]
        phases=M._model_mesh_generation_phases(3,case["status"])
        @test phases==(one=case["status"]<1,two=case["status"]<2,three=true)
        retained=!phases.one && (!phases.two || case["only_empty"]==1)
        @test retained==case["retained"]
        for product in case["product"]
            owner=Tuple(product["owner"]);rows=findall(==(owner),data.node_entities)
            @test data.external_node_tags[rows]==UInt64.(product["node_tags"])
            @test isequal(vec(mesh.coords[:,rows]),Float64[value for point in product["xyz_hex"] for value in hexpoint(point)])
            @test isequal(reduce(vcat,(data.node_parametric[row]===nothing ? Float64[] : data.node_parametric[row] for row in rows);init=Float64[]),hexpoint(product["uv_hex"]))
            observed=Any[]
            for (block,labels,owners) in zip(mesh.blocks,data.external_element_tags,data.block_entities)
                E.msh_dimension(block.msh)==owner[1] || continue
                columns=findall(==(Int32(owner[2])),owners);isempty(columns) && continue
                push!(observed,(Int(block.msh),labels[columns],[data.external_node_tags[block.nodes[:,column]] for column in columns]))
            end
            @test observed==[(Int(cell["msh"]),UInt64.(cell["labels"]),[UInt64.(nodes) for nodes in cell["nodes"]]) for cell in get(product,"cell",[])]
            owner[1]==2 || continue
            constructor=authority.surface_constructor_inputs[(2,Int32(owner[2]))]
            actual=M._model_mesh_product(authority,2,Int(owner[2]))
            @test M._model_surface_constructor_binding(actual,constructor,"native import")===nothing
            @test constructor.entity_data.external_element_tags==actual.entity_data.external_element_tags
        end
        @test_throws ArgumentError M._model_mesh_factory_authority(m,deepcopy(mesh),receipt,cleared)
        for changed_factory in (
                merge(plan.factory,(final_node_max=plan.factory.final_node_max+UInt64(1),)),
                merge(plan.factory,(final_element_max=plan.factory.final_element_max+UInt64(1),)),
                merge(plan.factory,(node_events=Base.front(plan.factory.node_events),)),
                merge(plan.factory,(cell_events=Base.front(plan.factory.cell_events),)),
                merge(plan.factory,(node_events=reverse(plan.factory.node_events),)),
                merge(plan.factory,(cell_events=reverse(plan.factory.cell_events),)))
            @test_throws ArgumentError M._model_mesh_factory_authority(m,mesh,
                merge(receipt,(factory_events=changed_factory,)),cleared)
        end
        saved=mesh.coords[1,1];mesh.coords[1,1]=nextfloat(saved)
        @test_throws ArgumentError M._model_mesh_factory_authority(m,mesh,receipt,cleared)
        mesh.coords[1,1]=saved
        changed=deepcopy(m);delete!(changed.surfaces,1)
        @test_throws ArgumentError M._model_mesh_factory_authority(changed,mesh,receipt,cleared)
        @test cleared.done==done && cleared.mesh===nothing
        after_clear=M._model_mesh_cleared_authority(authority)
        @test after_clear.done==done && after_clear.mesh===nothing
        @test isempty(after_clear.surface_constructor_inputs) && isempty(after_clear.products)
    end
end
end
