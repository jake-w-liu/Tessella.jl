module PointClosedContextTestsV1
using Tessella,Test
const A=Tessella.API
include("../../validation/api_point_position/point_closed_context_primary_v1.jl")
bits(values)=string.(reinterpret(UInt64,Float64[values...]);base=16,pad=16)
function open_text(text)
    mktempdir() do directory
        path=joinpath(directory,"closed.geo")
        write(path,text);A.open_geo!(path;mesh_dim=0)
    end
end
function check(expected)
    tags=A.mesh.get_nodes()[1]
    @test Set(tags)==Set(UInt64(node["uid"]) for node in expected["nodes"])
    for node in expected["nodes"]
        xyz,param,dimension,owner=A.mesh.get_node(node["uid"])
        @test (dimension,owner)==Tuple(Int.(node["owner"]))
        @test bits(xyz)==node["xyz_bits"]
        @test bits(param)==node["parameter_bits"]
    end
    @test A.mesh.get_nodes(0,1)[1]==UInt64.(expected["point_nodes"])
    types,ids,nodes=A.mesh.get_elements()
    @test types==Int32[cell["type"] for cell in expected["cells"]]
    @test ids==[UInt64.(cell["tags"]) for cell in expected["cells"]]
    @test nodes==[UInt64.(cell["nodes"]) for cell in expected["cells"]]
    @test A.mesh.get_max_node_tag()==expected["max_node"]
    @test A.mesh.get_max_element_tag()==expected["max_element"]
    @test A.model.get_value(0,1,Float64[])==Float64.(expected["point"])
    @test only(bits([A.option("General.BoundingBoxSize")]))==expected["lc_bits"]
end
@testset "Closed Curve endpoint receipts follow legitimate Point mutation" begin
    for row in POINT_CLOSED_CONTEXT_ROWS_V1
        @testset "$(row["name"])" begin
            try
                A.initialize()
                A.option("Mesh.Renumber",row["renumber"])
                A.option("Mesh.MeshOnlyEmpty",row["only_empty"])
                open_text(row["geo"])
                A.mesh.generate(1)
                check(row["before"])
                current=A.CURRENT[]
                @test haskey(current.meshing.native_curve_endpoints,1)
                previous=current.meshing.native_curve_endpoints[1]
                @test previous.curve==11 && previous.node==row["first_uid"]
                cache=A.LAST_MESH[];class=A.LAST_MESH_CLASS[]
                A.model.set_coordinates(1,.125,-.25,.375)
                check(row["after"])
                @test A.CURRENT[]===current && A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===class
                @test !haskey(current.meshing.native_curve_endpoints,1)
                A.mesh.generate(1)
                check(row["regenerated"])
                if row["only_empty"]==0
                    receipt=A.CURRENT[].meshing.native_curve_endpoints[1]
                    @test receipt.curve==11 && receipt.node==first(A.mesh.get_nodes(0,1)[1])
                    @test collect(receipt.point)==A.mesh.get_node(receipt.node)[1]
                end
                A.model.synchronize()
                check(row["synchronized"])
            finally
                A.finalize()
            end
        end
    end
end
end
