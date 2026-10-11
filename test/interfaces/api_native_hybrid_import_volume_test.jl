module NativeHybridImportVolumeChecks
using Test,TOML,SHA,Tessella
const M=Tessella.Model
const E=Tessella.Elements
const A=Tessella.API
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const transition=TOML.parsefile(joinpath(ROOT,"test/artifacts/native_hybrid_import_first_order_4152_v8.toml"))
const native=TOML.parsefile(joinpath(ROOT,"test/artifacts/native_hybrid_import_volume_4152_v7.toml"))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
function box_model()
    m=M.GeoModel()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(points);M.add_point!(m,p...;tag,mesh_size=.4);end
    curves=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(curves)
        M.add_line!(m,ends...;tag);M.set_transfinite_curve!(m,tag,3)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,curves) in enumerate(loops)
        M.add_curve_loop!(m,collect(curves);tag);M.add_plane_surface!(m,[tag];tag)
        M.set_transfinite_surface!(m,tag);tag==1 || M.set_recombine!(m,2,tag)
    end
    M.add_surface_loop!(m,collect(1:6);tag=1);M.add_volume!(m,[1];tag=1)
    return m
end
function first_order_input(m,plan)
    source=plan.authority.mesh;data=source.entity_data
    nodes=A._Dim01Node[A._Dim01Node(uid,Tuple(source.coords[:,row]),data.node_entities[row],
        data.node_parametric[row]===nothing ? nothing : copy(data.node_parametric[row]))
        for (row,uid) in enumerate(data.external_node_tags)]
    cells=A._Dim01Cell[A._Dim01Cell(Int(block.msh),(E.msh_dimension(block.msh),owners[column]),copy(block.nodes[:,column]),uid,uid,block.tags[column])
        for (block,owners,labels) in zip(source.blocks,data.block_entities,data.external_element_tags)
        for (column,uid) in enumerate(labels)]
    sort!(cells;by=A._dim01_cell_order)
    assembly=A._Dim01Assembly(nodes,cells,Dict(n.tag=>Int32(row) for (row,n) in enumerate(nodes)),
        plan.node_high_water,plan.element_high_water,2_000_000,4_000_000,"native import linearization",false,false)
    @test !any(A._dim01_first_order!(assembly,m))
    mapping=Dict(cell.source_tag=>cell.tag for cell in assembly.cells)
    @test length(mapping)==length(assembly.cells) && allunique(values(mapping))
    actual=deepcopy(source)
    for labels in actual.entity_data.external_element_tags
        labels.=[mapping[uid] for uid in labels]
    end
    authority=M._model_mesh_rebind_authority(plan.authority,actual)
    expected=Set(uid for constructor in values(plan.authority.surface_constructor_inputs)
        for labels in constructor.entity_data.external_element_tags for uid in labels)
    surface_mapping=Dict(uid=>mapping[uid] for uid in expected)
    inputs=M._model_mesh_first_order_surface_constructor_inputs(plan.authority,authority,surface_mapping)
    @test length(inputs)==length(plan.authority.surface_constructor_inputs)
    for (owner,constructor) in inputs
        @test [uid for labels in constructor.entity_data.external_element_tags for uid in labels]==
            [mapping[uid] for labels in plan.authority.surface_constructor_inputs[owner].entity_data.external_element_tags for uid in labels]
    end
    before=deepcopy(plan.authority.mesh);after=deepcopy(authority.mesh)
    for bad in (Dict{UInt64,UInt64}(),merge(surface_mapping,Dict(typemax(UInt64)=>UInt64(1))),
            Dict(uid=>UInt64(1) for uid in keys(surface_mapping)))
        @test_throws ArgumentError M._model_mesh_first_order_surface_constructor_inputs(plan.authority,authority,bad)
        @test isequal(before.coords,plan.authority.mesh.coords) && isequal(after.coords,authority.mesh.coords)
        @test before.entity_data.external_element_tags==plan.authority.mesh.entity_data.external_element_tags &&
            after.entity_data.external_element_tags==authority.mesh.entity_data.external_element_tags
    end
    original=authority.mesh.coords[1,1];authority.mesh.coords[1,1]=nextfloat(original)
    @test_throws ArgumentError M._model_mesh_first_order_surface_constructor_inputs(plan.authority,authority,surface_mapping)
    authority.mesh.coords[1,1]=original
    authority=M._model_mesh_with_surface_constructor_inputs(authority,inputs)
    return (authority=authority,node_high_water=assembly.node_max,element_high_water=assembly.cell_max,mapping=mapping)
end
function check_first_order(actual,case)
    data=actual.entity_data
    for product in case["product"]
        owner=Tuple((Int(product["owner"][1]),Int32(product["owner"][2])))
        rows=findall(==(owner),data.node_entities)
        @test data.external_node_tags[rows]==UInt64.(product["node_tags"])
        @test isequal(actual.coords[:,rows],isempty(rows) ? zeros(3,0) : reduce(hcat,hexpoint.(product["xyz_hex"])))
        @test isequal(reduce(vcat,(data.node_parametric[row]===nothing ? Float64[] : data.node_parametric[row] for row in rows);init=Float64[]),hexpoint(product["uv_hex"]))
        observed=NamedTuple[]
        for (block,owners,labels) in zip(actual.blocks,data.block_entities,data.external_element_tags)
            E.msh_dimension(block.msh)==owner[1] || continue
            columns=findall(==(owner[2]),owners);isempty(columns) && continue
            push!(observed,(msh=Int(block.msh),labels=labels[columns],nodes=[Tuple(data.external_node_tags[block.nodes[:,column]]) for column in columns]))
        end
        @test observed==[(msh=Int(cell["msh"]),labels=UInt64.(cell["labels"]),nodes=Tuple.(Vector{UInt64}.(cell["nodes"]))) for cell in get(product,"cell",[])]
    end
end
@testset "actual Model consumes ingestion-owned mixed Face constructor vectors" begin
    @test bytes2hex(sha256(read(joinpath(ROOT,native["primary"]))))==native["primary_sha256"]
    @test length(native["case"])==6
    for case in native["case"]
        println("IMPORTED_HYBRID_VOLUME_CASE ",case["name"]);flush(stdout)
        m=box_model();done=Set{Tuple{Int,Int32}}((dim,Int32(tag)) for (dim,tag) in M.model_entities(m) if dim in (1,2))
        previous=M._ModelMeshAuthority(nothing,Dict{Tuple{Int,Int32},M._ModelMeshProduct}(),done,UInt64[])
        path=joinpath(ROOT,"test/data/native_hybrid_import_events_4152_v1",case["name"]*".msh")
        mesh,identity=E._read_mixed_msh_factory_with_identities(path,
            E._mixed_msh_factory_context(case["initial_node_max"],case["initial_element_max"]))
        plan=M._model_mesh_factory_authority(m,mesh,identity,previous)
        @test M._model_mesh_status(m,plan.authority)==2
        @test (plan.node_high_water,plan.element_high_water)==
            (UInt64(case["import_node_max"]),UInt64(case["import_element_max"]))
        saved=deepcopy(mesh)
        source_plan=first_order_input(m,plan)
        first_case=only(filter(row->row["name"]==case["name"],transition["case"]))
        @test (source_plan.node_high_water,source_plan.element_high_water)==(UInt64(first_case["node_max"]),UInt64(first_case["element_max"]))
        check_first_order(source_plan.authority.mesh,first_case)
        options=M._model_default_mesh1d_options(m,"native imported hybrid")
        options.random_state[]=UInt32(1)
        options=M._ModelMesh1DOptions((name===:optimize ? 0 : getfield(options,name) for name in fieldnames(M._ModelMesh1DOptions))...)
        production=M.mesh_model_volume(m,1;_mesh_input=source_plan.authority,_mesh1d_options=options,
            _return_volume_production=true,_mesh_node_high_water=source_plan.node_high_water,
            _mesh_element_high_water=source_plan.element_high_water)
        @test production isa M._ModelNativeHybridVolumeProduction
        @test production.input===source_plan.authority
        @test production.node_high_water==UInt64(case["final_node_max"])
        @test production.element_high_water==UInt64(case["final_element_max"])
        @test production.element_tags==UInt64[uid for cell in case["cell"] for uid in cell["labels"]]
        @test [Tuple(production.boundary_node_tags[nodes]) for block in production.mesh.blocks for nodes in eachcol(block.nodes)]==
            [Tuple(UInt64.(nodes)) for cell in case["cell"] for nodes in cell["nodes"]]
        positions=Dict(uid=>row for (row,uid) in enumerate(production.boundary_node_tags))
        @test isequal(production.mesh.coords[:,[positions[UInt64(uid)] for uid in case["node_tags"]]],reduce(hcat,hexpoint.(case["xyz_hex"])))
        original=Set(mesh.entity_data.external_node_tags)
        @test UInt64[uid for uid in production.boundary_node_tags[production.registration_rows] if !(uid in original)]==UInt64.(case["node_tags"])
        @test M._model_native_hybrid_volume_certificate(production.mesh,production.boundary_node_tags,production.constructor,"native imported hybrid")===nothing
        @test isequal(mesh.coords,saved.coords)
        @test isequal(mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
        @test [b.nodes for b in mesh.blocks]==[b.nodes for b in saved.blocks]
        @test mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
    end
end
end
