module NativePointPositionTestsV1
using Tessella,Test
const A=Tessella.API

function setup()
    A.initialize()
    A.option("Mesh.Renumber",0)
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
    A.mesh.add_nodes(0,1,[900,7,400],[0.,0.,0.,.01,.02,.03,.04,.05,.06])
    A.mesh.add_elements_by_type(1,15,[777,778],[900,7])
end

function snapshot()
    tags=A.mesh.get_nodes()[1]
    nodes=Dict{UInt64,Any}(tag=>A.mesh.get_node(tag) for tag in tags)
    return (nodes=nodes,cells=A.mesh.get_elements(),max_node=A.mesh.get_max_node_tag(),
        max_element=A.mesh.get_max_element_tag(),point=A.model.get_value(0,1,Float64[]))
end

@testset "Point setters retain mesh identities and reject invalid input" begin
    try
        setup()
        before=snapshot()
        cache=A.LAST_MESH[];classification=A.LAST_MESH_CLASS[]
        authority=A.CURRENT[].meshing.mesh_authority
        A.model.set_coordinates(1,.125,-.25,.375)
        after=snapshot()
        @test after.point==[.125,-.25,.375]
        @test Set(keys(after.nodes))==Set(keys(before.nodes))
        @test after.cells==before.cells
        @test (after.max_node,after.max_element)==(before.max_node,before.max_element)
        @test A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===classification
        @test A.CURRENT[].meshing.mesh_authority===authority
        A.model.set_coordinates(1,.125,-.25,.375)
        @test snapshot()==after
        @test_throws ArgumentError A.model.set_coordinates(999,.5,.5,.5)
        @test snapshot()==after && A.LAST_MESH[]===cache
        @test_throws ArgumentError A.model.set_coordinates(1,Inf,.5,.5)
        @test snapshot()==after && A.LAST_MESH[]===cache
    finally
        A.finalize()
    end
end
end
