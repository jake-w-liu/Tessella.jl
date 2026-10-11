module DefaultBoxOptimizerFullNativeV1
using Test,Tessella
include(joinpath(@__DIR__,"..","artifacts","default_box_optimizer_full_native_v1.jl"))

bits(values)=[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
row(unused_points)=only(filter(value->value["unused_points"]==unused_points,
    DEFAULT_BOX_OPTIMIZER_FULL_NATIVE_V1["rows"]))
function node_bits(native,tags)
    nodes=Dict(node["uid"]=>node["XYZ_bits"] for node in native["state"]["node_records"])
    return reduce(vcat,(nodes[Int(tag)] for tag in tags);init=String[])
end

function check(api,carrier;unused_points=false)
    native=row(unused_points);state=native["state"]
    global_nodes=state["global_nodes"];global_elements=state["global_elements"]
    tags,xyz,uv=api.mesh.get_nodes()
    @test tags==UInt64.(global_nodes["uids"])
    @test bits(xyz)==global_nodes["XYZ_bits"]
    @test bits(uv)==global_nodes["UV_bits"]
    types,elements,cells=api.mesh.get_elements()
    @test types==Int32.(global_elements["types"])
    @test elements==[UInt64.(block) for block in global_elements["uids"]]
    @test cells==[UInt64.(block) for block in global_elements["connectivity"]]
    @test api.mesh.get_max_node_tag()==UInt64(state["maxima"][1])
    @test api.mesh.get_max_element_tag()==UInt64(state["maxima"][2])
    for entity in state["per_entity"]
        dim,tag=entity["owner"]
        tags,xyz,uv=api.mesh.get_nodes(dim,tag,false,true)
        @test tags==UInt64.(entity["node_uids"])
        @test bits(xyz)==entity["XYZ_bits"]
        @test bits(uv)==entity["UV_bits"]
        types,elements,cells=api.mesh.get_elements(dim,tag)
        @test types==Int32.(entity["types"])
        @test elements==[UInt64.(block) for block in entity["element_uids"]]
        @test cells==[UInt64.(block) for block in entity["connectivity"]]
    end
    for expected in state["node_records"]
        xyz,uv,dim,tag=api.mesh.get_node(expected["uid"])
        @test bits(xyz)==expected["XYZ_bits"]
        @test bits(uv)==expected["UV_bits"]
        @test [Int(dim),Int(tag)]==expected["owner"]
    end
    for expected in state["element_records"]
        kind,vertices,dim,tag=api.mesh.get_element(expected["uid"])
        @test Int(kind)==expected["type"]
        @test vertices==UInt64.(expected["connectivity"])
        @test [Int(dim),Int(tag)]==expected["owner"]
    end
    # Returned geometry has its own stable storage. The independently proved
    # incidence map binds those rows to the complete actual public products.
    expected=native["carrier"]
    @test carrier isa Tessella.MeshTypes.Mesh
    @test bits(vec(carrier.coords))==expected["coordinate_bits"]
    for (kind,cells) in ((1,carrier.segs),(2,carrier.tris),(4,carrier.tets))
        block=only(filter(value->value["msh"]==kind,expected["blocks"]))
        @test vec(cells)==block["connectivity"]
    end
    @test Tessella.MeshTypes.mesh_crc(carrier).sha==expected["crc_sha"]
    authority=api.CURRENT[].meshing.mesh_authority
    @test authority.mesh===api.LAST_MESH[]
    actual_tags=authority.mesh.entity_data.external_node_tags
    @test actual_tags[authority.carrier_rows]==UInt64.(native["carrier_rows_to_native_UIDs"])
    @test bits(vec(authority.mesh.coords[:,authority.carrier_rows]))==expected["coordinate_bits"]
    return native
end
end
