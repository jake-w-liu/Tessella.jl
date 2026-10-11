module NativeHybridMshFactoryChecks
using Test,TOML,SHA,Tessella
const E=Tessella.Elements
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const native=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","native_hybrid_msh_factory_4152_v4.toml"))
hexpoint(values)=Float64[Core.eval(@__MODULE__,Meta.parse(value)) for value in values]
@testset "actual MSH4 factory replays literal native file events into final identity" begin
    @test bytes2hex(sha256(read(joinpath(ROOT,native["primary"]))))==native["primary_sha256"]
    @test length(native["case"])==11
    snapshots=Dict{String,Any}()
    for case in native["case"]
        path=joinpath(ROOT,case["path"])
        @test bytes2hex(sha256(read(path)))==case["sha256"]
        input=E._mixed_msh_factory_context(case["initial_node_max"],case["initial_element_max"])
        initial=(input.node_max,input.element_max)
        mesh,receipt=E._read_mixed_msh_factory_with_identities(path,input)
        factory=receipt.factory_events
        @test (input.node_max,input.element_max)==initial
        @test factory.source_sha256==case["sha256"]
        @test (factory.initial_node_max,factory.initial_element_max)==UInt64.(initial)
        @test (factory.final_node_max,factory.final_element_max)==
            (UInt64(case["final_node_max"]),UInt64(case["final_element_max"]))
        @test factory.node_events==Tuple((Int(row[1]),Int32(row[2]),UInt64(row[3]),UInt64(row[4])) for row in case["node_events"])
        @test factory.cell_events==Tuple((Int(row[1]),Int32(row[2]),UInt64(row[3]),UInt64(row[4])) for row in case["cell_events"])
        @test receipt.cell_events==Tuple((Int(row[1]),Int32(row[2]),UInt64(row[3])) for row in case["cell_events"])
        @test receipt.node_tags==UInt64.(case["final_node_tags"])
        @test mesh.entity_data.external_node_tags==UInt64.(case["final_node_tags"])
        @test Int[block.msh for block in mesh.blocks]==Int.(case["bucket_types"])
        @test receipt.element_tags==[UInt64.(tags) for tags in case["bucket_tags"]]
        @test mesh.entity_data.external_element_tags==[UInt64.(tags) for tags in case["bucket_tags"]]
        @test E._mixed_msh_factory_binding(mesh,factory)
        @test Tessella.MeshTypes.validate(mesh).ok
        for product in case["product"]
            owner=Tuple(product["owner"]);rows=findall(==(owner),mesh.entity_data.node_entities)
            @test mesh.entity_data.external_node_tags[rows]==UInt64.(product["node_tags"])
            @test isequal(vec(mesh.coords[:,rows]),Float64[value for point in product["xyz_hex"] for value in hexpoint(point)])
            @test isequal(reduce(vcat,(mesh.entity_data.node_parametric[row]===nothing ? Float64[] : mesh.entity_data.node_parametric[row] for row in rows);init=Float64[]),hexpoint(product["uv_hex"]))
            observed=Any[]
            for (block,labels,owners) in zip(mesh.blocks,mesh.entity_data.external_element_tags,mesh.entity_data.block_entities)
                E.msh_dimension(block.msh)==owner[1] || continue
                columns=findall(==(Int32(owner[2])),owners);isempty(columns) && continue
                push!(observed,(Int(block.msh),labels[columns],[mesh.entity_data.external_node_tags[block.nodes[:,column]] for column in columns]))
            end
            @test observed==[(Int(cell["msh"]),UInt64.(cell["labels"]),[UInt64.(nodes) for nodes in cell["nodes"]]) for cell in get(product,"cell",[])]
        end
        snapshots[case["name"]]=(mesh=mesh,receipt=receipt)
        if any(iszero,case["raw_node_tags"]) || any(row->iszero(row[3]),case["cell_events"])
            @test_throws ArgumentError E._read_mixed_msh_with_identities(path)
        else
            generic,identities=E._read_mixed_msh_with_identities(path)
            @test identities.factory_events===nothing
            @test isequal(generic.coords,mesh.coords)
            @test [block.nodes for block in generic.blocks]==[block.nodes for block in mesh.blocks]
            @test generic.entity_data.external_element_tags==receipt.element_tags
        end
        @test_throws ArgumentError E._read_mixed_msh_factory_with_identities(path,input;max_elements=1)
        @test (input.node_max,input.element_max)==initial
    end
    ordinary=snapshots["ordinary"];mesh=ordinary.mesh;factory=ordinary.receipt.factory_events
    @test !E._mixed_msh_factory_binding(deepcopy(mesh),factory)
    saved=mesh.coords[1,1];mesh.coords[1,1]=nextfloat(saved)
    @test !E._mixed_msh_factory_binding(mesh,factory)
    mesh.coords[1,1]=saved
    saved=mesh.entity_data.external_node_tags[1];mesh.entity_data.external_node_tags[1]=UInt64(900)
    @test !E._mixed_msh_factory_binding(mesh,factory)
    mesh.entity_data.external_node_tags[1]=saved
    saved=mesh.entity_data.node_entities[1];mesh.entity_data.node_entities[1]=(0,Int32(2))
    @test !E._mixed_msh_factory_binding(mesh,factory)
    mesh.entity_data.node_entities[1]=saved
    param=findfirst(row->row!==nothing && !isempty(row),mesh.entity_data.node_parametric)
    saved=mesh.entity_data.node_parametric[param][1];mesh.entity_data.node_parametric[param][1]=nextfloat(saved)
    @test !E._mixed_msh_factory_binding(mesh,factory)
    mesh.entity_data.node_parametric[param][1]=saved
    quad=findfirst(block->block.msh==3,mesh.blocks)
    saved=copy(mesh.blocks[quad].nodes[:,1]);mesh.blocks[quad].nodes[1,1],mesh.blocks[quad].nodes[2,1]=mesh.blocks[quad].nodes[2,1],mesh.blocks[quad].nodes[1,1]
    @test !E._mixed_msh_factory_binding(mesh,factory)
    mesh.blocks[quad].nodes[:,1]=saved
    saved=mesh.entity_data.external_element_tags[quad][1];mesh.entity_data.external_element_tags[quad][1]=UInt64(900)
    @test !E._mixed_msh_factory_binding(mesh,factory)
    mesh.entity_data.external_element_tags[quad][1]=saved
    saved=mesh.entity_data.block_entities[quad][1];mesh.entity_data.block_entities[quad][1]=Int32(2)
    @test !E._mixed_msh_factory_binding(mesh,factory)
    mesh.entity_data.block_entities[quad][1]=saved
    saved=mesh.blocks[quad].tags[1];mesh.blocks[quad].tags[1]=Int32(900)
    @test !E._mixed_msh_factory_binding(mesh,factory)
    mesh.blocks[quad].tags[1]=saved
    @test E._mixed_msh_factory_binding(mesh,factory)
    @test ordinary.receipt.factory_events.cell_events!=snapshots["family_block_reverse"].receipt.factory_events.cell_events
    @test snapshots["zero"].receipt.factory_events.final_element_max==UInt64(66)
    @test snapshots["zero_family_block_reverse"].receipt.factory_events.final_element_max==UInt64(65)
    @test snapshots["zero"].receipt.node_tags[1]==UInt64(27)
    @test snapshots["cold_zero"].receipt.node_tags[1]==UInt64(1)
    @test_throws ArgumentError E._mixed_msh_factory_context(true,0)
    @test_throws ArgumentError E._mixed_msh_factory_context(0,false)
    @test_throws ArgumentError E._mixed_msh_factory_context(-1,0)
    @test_throws ArgumentError E._mixed_msh_factory_context(0,big(typemax(UInt64))+1)
    zero_case=only(filter(case->case["name"]=="zero",native["case"]))
    path=joinpath(ROOT,zero_case["path"])
    @test_throws ArgumentError E._read_mixed_msh_factory_with_identities(path,E._mixed_msh_factory_context(typemax(UInt64),60))
    @test_throws ArgumentError E._read_mixed_msh_factory_with_identities(path,E._mixed_msh_factory_context(1,60))
    collision=E._mixed_msh_factory_context(20,32)
    @test_throws ArgumentError E._read_mixed_msh_factory_with_identities(path,collision)
    @test (collision.node_max,collision.element_max)==(UInt64(20),UInt64(32))
    @test_throws ArgumentError E._read_mixed_msh_factory_with_identities(path,E._mixed_msh_factory_context(26,typemax(UInt64)))
end
end
