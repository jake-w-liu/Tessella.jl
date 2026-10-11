module GeometryExpressionNativeControls
using Tessella,Test,SHA
const E=Tessella.Elements
const MT=Tessella.MeshTypes
const NATIVE_API_LITERAL=joinpath(@__DIR__,"..","..","validation","native_initial_tet","quality32_native_v1.jl")
include(NATIVE_API_LITERAL)
const NATIVE_MSH=joinpath(@__DIR__,"..","artifacts","geo_geometry_expression_native_parametric_v1.msh")
const NATIVE_MSH_SHA="d3f3e1887b6919353c91539582f8cf62085d13513e8212c439016bc5940161ff"

function expected_api_projection()
    state=INITIAL_TET_QUALITY32_NATIVE_V1["original_default"]
    @test bytes2hex(sha256(read(NATIVE_API_LITERAL)))==
        "b51d83698a2ef85f43d4755b7592618ea4f6de0f7ead795fdf11dd356b6f512e"
    @test bytes2hex(sha256(read(joinpath(@__DIR__,"..","fixtures","geo_geometry_expressions.geo"))))==
        INITIAL_TET_QUALITY32_NATIVE_V1["scene_sha256"]
    @test bytes2hex(sha256(read(NATIVE_MSH)))==NATIVE_MSH_SHA
    serialized=read_mixed_msh(NATIVE_MSH)
    @test E.mixed_crc(serialized).sha==
        "e0c23dc0310041f95b2b22a26bab6319fac2f00efcce1025a2dfc0d50beba2b0"
    uids=UInt64.(state["nodes"]["node_tags"])
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(uids))
    coords=reshape([reinterpret(Float64,parse(UInt64,value;base=16))
        for value in state["nodes"]["xyz_bits"]],3,:)
    owners=Vector{Tuple{Int,Int32}}(undef,length(uids))
    parameters=Vector{Union{Nothing,Vector{Float64}}}(undef,length(uids))
    element_owners=Dict{UInt64,Int32}();physical=Dict{UInt64,Int32}()
    for entity in state["entities"]
        dim,tag=entity["dimension"],entity["tag"]
        for (index,uid) in enumerate(entity["nodes"]["node_tags"])
            row=positions[uid];owners[row]=(dim,Int32(tag))
            uv=entity["nodes"]["uv_bits"]
            # Native API reports empty Point vectors. The MSH parser uses
            # nothing for its absent parametric marker; these are distinct
            # valid storage routes with no coordinates at dimension zero.
            parameters[row]=dim==0 ? Float64[] :
                [reinterpret(Float64,parse(UInt64,value;base=16))
                 for value in uv[dim*(index-1)+1:dim*index]]
        end
        group=serialized.entity_data.entities[(dim,tag)].physical_tags
        @test length(group)==1
        for tags in entity["cells"]["element_tags"],uid in tags
            element_owners[uid]=Int32(tag);physical[uid]=only(group)
        end
    end
    blocks=E.ElementBlock[];block_owners=Vector{Int32}[];element_tags=Vector{UInt64}[]
    for (kind,tags,corners) in zip(state["cells"]["types"],
            state["cells"]["element_tags"],state["cells"]["connectivity"])
        width=E.msh_num_nodes(kind)
        cells=reshape(Int32[positions[uid] for uid in corners],width,:)
        push!(blocks,E.ElementBlock(kind,cells,Int32[physical[uid] for uid in tags]))
        push!(block_owners,Int32[element_owners[uid] for uid in tags])
        push!(element_tags,UInt64.(tags))
    end
    data=E.MixedEntityData(serialized.entity_data.entities;node_entities=owners,
        node_parametric=parameters,external_node_tags=uids,block_entities=block_owners,
        external_element_tags=element_tags)
    expected=E.MixedMesh(coords,blocks;physical_names=serialized.physical_names,entity_data=data)
    @test E.mixed_crc(expected).sha==
        "7e1b27f1cd12b46c765b72e222f33f30d60a62759077746f5a747f296ed5fb90"
    return expected
end

function check(execution,projected)
    expected=expected_api_projection()
    @test E.mixed_crc(projected)==E.mixed_crc(expected)
    authority=execution.model.meshing.mesh_authority
    uids=authority.mesh.entity_data.external_node_tags[authority.carrier_rows]
    @test allunique(uids) && Set(uids)==Set(expected.entity_data.external_node_tags)
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(uids))
    cells=Dict{Int,Matrix{Int32}}()
    for block in expected.blocks
        block.msh in (1,2,4) || continue
        cells[block.msh]=map(row->positions[expected.entity_data.external_node_tags[row]],block.nodes)
    end
    # The renderer retains its own row/coordinate storage. Only explicit
    # authority UID bindings transport native cells into that storage.
    native_positions=Dict(uid=>row for (row,uid) in enumerate(expected.entity_data.external_node_tags))
    carrier_coords=expected.coords[:,[native_positions[uid] for uid in uids]]
    carrier=MT.Mesh(carrier_coords;segs=cells[1],tris=cells[2],tets=cells[4])
    @test MT.mesh_crc(carrier).sha==MT.mesh_crc(execution.mesh).sha
    @test MT.mesh_crc(carrier).sha==
        "a648e22f1ae8d58bce2638ed1533b5826755fdd8c67f53c9fb4b1ed243cbccd9"
    return nothing
end
end
