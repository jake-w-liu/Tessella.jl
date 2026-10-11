module NativePointPositionTestsV1
using Tessella,Test,TOML
const A=Tessella.API
include("../../validation/api_point_position/native_point_position_primary_v1.jl")

function setup(spec)
    A.initialize()
    A.option("Mesh.Renumber",spec["renumber"])
    A.option("Mesh.ElementOrder",1)
    for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        A.model.add_point(x,y,0.;tag)
    end
    for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        A.model.add_line(first,last;tag)
    end
    A.model.add_curve_loop([1,2,3,4];tag=1)
    A.model.add_plane_surface([1];tag=1)
    for tag in 1:4
        A.mesh.set_transfinite_curve(tag,3)
    end
    A.mesh.set_transfinite_surface(1)
    raw=spec["raw"]
    if raw=="single"
        A.mesh.add_nodes(0,1,[111],[0.,0.,0.])
    elseif raw=="multiple"
        A.mesh.add_nodes(0,1,[900,7,400],[0.,0.,0.,.01,.02,.03,.04,.05,.06])
        A.mesh.add_elements_by_type(1,15,[777,778],[900,7])
    elseif raw=="incidence"
        for (owner,uid,point) in ((1,111,(0.,0.,0.)),(2,222,(1.,0.,0.)),(3,333,(1.,1.,0.)))
            A.mesh.add_nodes(0,owner,[uid],collect(point))
        end
        A.mesh.add_elements_by_type(1,15,[777],[111])
        A.mesh.add_elements_by_type(1,1,[888],[111,222])
        A.mesh.add_elements_by_type(1,2,[999],[111,222,333])
    elseif raw in ("cellfree","cells")
        A.mesh.add_nodes(0,1,[101,102,103],zeros(9))
        raw=="cells" && A.mesh.add_elements_by_type(1,15,[777,778],[101,102])
    end
    if spec["dimension"]>0
        A.mesh.generate(spec["dimension"])
        spec["order"]==2 && A.mesh.set_order(2)
    end
    if get(spec,"sparse",false)
        tags=A.mesh.get_nodes()[1]
        A.mesh.renumber_nodes(tags,100000 .+ 17 .* tags)
    end
end

function snapshot()
    tags=A.mesh.get_nodes()[1]
    nodes=Dict{UInt64,Any}(tag=>A.mesh.get_node(tag) for tag in tags)
    return (nodes=nodes,cells=A.mesh.get_elements(),max_node=A.mesh.get_max_node_tag(),
        max_element=A.mesh.get_max_element_tag(),point=A.model.get_value(0,1,Float64[]))
end

function native_nodes(row)
    Dict(UInt64(node["uid"])=>node for node in row["nodes_before"])
end

const observations=Dict{String,Any}[]
@testset "Point coordinate setters update the current first attached vertex" begin
    for row in NATIVE_POINT_POSITION_ROWS_V1
        spec=row["spec"]
        @testset "$(spec["name"])" begin
            try
                setup(spec)
                before=snapshot()
                cache=A.LAST_MESH[];classification=A.LAST_MESH_CLASS[]
                authority=A.CURRENT[].meshing.mesh_authority
                point_nodes=A.mesh.get_nodes(0,1)[1]
                @test point_nodes==UInt64.(row["point_nodes_before"])
                native=native_nodes(row)
                @test Set(keys(before.nodes))==Set(keys(native))
                for (uid,node) in before.nodes
                    expected=native[uid]
                    @test node[3:4]==Tuple(Int.(expected["owner"]))
                    @test isapprox(node[1],Float64.(expected["xyz"]);atol=2e-11,rtol=0)
                    @test isapprox(node[2],Float64.(expected["parameters"]);atol=2e-11,rtol=0)
                end
                @test before.max_node==row["original_max_node"]
                @test before.max_element==row["original_max_element"]
                A.model.set_coordinates(1,.125,-.25,.375)
                after=snapshot()
                @test after.point==[.125,-.25,.375]
                @test Set(keys(after.nodes))==Set(keys(before.nodes))
                @test after.cells==before.cells
                @test (after.max_node,after.max_element)==(before.max_node,before.max_element)
                @test A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===classification
                @test A.CURRENT[].meshing.mesh_authority===authority
                changed=UInt64[]
                for (uid,previous) in before.nodes
                    current=after.nodes[uid]
                    current[1]==previous[1] || push!(changed,uid)
                    @test current[2:4]==previous[2:4]
                    expected=uid==row["first_uid"] ? [.125,-.25,.375] : previous[1]
                    @test current[1]==expected
                end
                @test changed==UInt64.(row["changed_node_uids"])
                A.model.set_coordinates(1,.125,-.25,.375)
                @test snapshot()==after
                @test_throws ArgumentError A.model.set_coordinates(999,.5,.5,.5)
                @test snapshot()==after && A.LAST_MESH[]===cache
                @test_throws ArgumentError A.model.set_coordinates(1,Inf,.5,.5)
                @test snapshot()==after && A.LAST_MESH[]===cache
                push!(observations,Dict("name"=>spec["name"],"first_uid"=>row["first_uid"],
                    "node_count"=>length(after.nodes),"max_node"=>after.max_node,
                    "max_element"=>after.max_element))
            finally
                A.finalize()
            end
        end
    end
end
open(joinpath(@__DIR__,"../tmp/point_position_observed_v1_$(VERSION).toml"),"w") do stream
    TOML.print(stream,Dict("observations"=>observations))
end
end
