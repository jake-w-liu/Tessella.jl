module NativeHybridMshParserEventChecks
using Test,TOML,SHA,Tessella
const E=Tessella.Elements
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const native=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","native_hybrid_msh_parser_events_4152_v1.toml"))
hexpoint(values)=Float64[Core.eval(@__MODULE__,Meta.parse(value)) for value in values]
@testset "actual MSH parser records constructor encounter order before family buckets" begin
    @test bytes2hex(sha256(read(joinpath(ROOT,native["primary"]))))==native["primary_sha256"]
    @test length(native["case"])==11
    receipts=Dict{String,Any}()
    for case in native["case"]
        path=joinpath(ROOT,case["path"])
        @test bytes2hex(sha256(read(path)))==case["sha256"]
        # File zero is a constructor event needing a counter context, whereas
        # this default reader produces literal positive identities. Keep that
        # distinction strict until the ingestion operator replays the events.
        if case["zero_constructor_input"]
            reference=Ref{Any}(:untouched)
            @test_throws ArgumentError E.read_mixed_msh(path;_source_identities=reference)
            @test reference[]===:untouched
            continue
        end
        mesh,receipt=E._read_mixed_msh_with_identities(path)
        @test receipt.cell_events isa Tuple
        @test receipt.cell_events==Tuple((Int(row[1]),Int32(row[2]),UInt64(row[3])) for row in case["events"])
        @test length(receipt.cell_events)==sum(size(block.nodes,2) for block in mesh.blocks)
        @test receipt.node_tags==UInt64.(case["source_node_tags"])
        @test mesh.entity_data.external_node_tags==UInt64.(case["source_node_tags"])
        @test Int[block.msh for block in mesh.blocks]==Int.(case["source_bucket_types"])
        @test receipt.element_tags==[UInt64.(tags) for tags in case["source_bucket_tags"]]
        @test mesh.entity_data.external_element_tags==[UInt64.(tags) for tags in case["source_bucket_tags"]]
        @test receipt.is_v4
        @test Tessella.MeshTypes.validate(mesh).ok
        receipts[case["name"]]=receipt
        data=mesh.entity_data
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
        end
        generic=E.read_mixed_msh(path)
        @test isequal(generic.coords,mesh.coords)
        @test [block.nodes for block in generic.blocks]==[block.nodes for block in mesh.blocks]
        @test generic.entity_data.external_element_tags==receipt.element_tags
        reference=Ref{Any}(:untouched)
        @test_throws ArgumentError E.read_mixed_msh(path;_source_identities=reference,max_elements=1)
        @test reference[]===:untouched
    end
    @test receipts["ordinary"].cell_events!=receipts["family_block_reverse"].cell_events
    @test receipts["ordinary"].cell_events!=receipts["within_family_reverse"].cell_events
    face_quads(receipt)=UInt64[uid for (kind,owner,uid) in receipt.cell_events if kind==3 && owner==1]
    @test face_quads(receipts["ordinary"])==UInt64[61,62,63]
    @test face_quads(receipts["family_block_reverse"])==UInt64[61,62,63]
    @test face_quads(receipts["within_family_reverse"])==UInt64[63,62,61]
end
end
