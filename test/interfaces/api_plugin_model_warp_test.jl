module PluginModelWarpTests
using Test,Tessella
const A=Tessella.API

function fixture()
    A.model.add("warp")
    A.model.add_discrete_entity(2,1)
    A.mesh.add_nodes(2,1,[1,2,3,4],
        [0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.])
    A.mesh.add_elements_by_type(1,2,[11,12],[1,2,3,2,4,3])
    scalar=A.view.add("scalar")
    A.view.add_model_data(scalar,0,"warp","NodeData",[1,2,3,4],
        [[1.],[2.],[3.],[4.]],0.,1)
    vector=A.view.add("vector")
    A.view.add_model_data(vector,0,"warp","NodeData",[1,2,3,4],
        [1.,0.,0., 0.,1.,0., 0.,0.,1., 1.,1.,0.],0.,3)
    return scalar,vector
end

function coords()
    tags,flat,_=A.mesh.get_nodes()
    return Dict(tag=>[flat[3(i-1)+a] for a in 1:3] for (i,tag) in enumerate(tags))
end

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        f()
    finally
        A.finalize()
    end
end

@testset "Model Warp displaces actual mesh nodes from a vector OtherView" begin
    session() do
        scalar,vector=fixture()
        nodes=A.mesh.get_nodes()
        elements=A.mesh.get_elements()
        model=A.model.get_current()
        A.plugin.set_number("Warp","View",scalar)
        A.plugin.set_number("Warp","OtherView",vector)
        A.plugin.set_number("Warp","Factor",2.)
        A.plugin.run("Warp")
        actual=coords()
        @test actual[UInt64(1)]==[2.,0.,0.]
        @test actual[UInt64(2)]==[1.,2.,0.]
        @test actual[UInt64(3)]==[0.,1.,2.]
        @test actual[UInt64(4)]==[3.,3.,0.]
        @test A.model.get_current()==model
        @test A.mesh.get_nodes()[1]==nodes[1]
        @test A.mesh.get_elements()==elements
        # View payloads are unchanged by Warp.
        @test A.view.get_model_data(scalar,0)[3]==[1.,2.,3.,4.]
        @test A.view.get_model_data(vector,0)[3]==
            [1.,0.,0.,0.,1.,0.,0.,0.,1.,1.,1.,0.]
    end
end

@testset "Model Warp honours aliases, Factor zero and missing OtherView" begin
    session() do
        scalar,vector=fixture()
        alias=A.view.add_alias(vector)
        A.plugin.set_number("Warp","View",scalar)
        A.plugin.set_number("Warp","OtherView",alias)
        A.plugin.set_number("Warp","Factor",1.)
        A.plugin.run("Warp")
        @test coords()[UInt64(1)]==[1.,0.,0.]
        A.plugin.set_number("Warp","Factor",0.)
        A.plugin.run("Warp")
        @test coords()[UInt64(1)]==[1.,0.,0.]
    end
    session() do
        scalar,vector=fixture()
        before=coords()
        # A model view without OtherView is upstream's successful no-op.
        A.plugin.set_number("Warp","View",scalar)
        A.plugin.set_number("Warp","OtherView",-1.)
        A.plugin.run("Warp")
        @test coords()==before
    end
end

const AFF_X=Float64[1,0,0,1, 0,1,0,0, 0,0,1,0, 0,0,0,1]

function periodic_box()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,point) in enumerate(points)
        A.model.add_point(point...;tag,meshSize=0.5)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(edges)
        A.model.add_line(ends...;tag)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        A.model.add_curve_loop(collect(lines);tag)
        A.model.add_plane_surface([tag];tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
end

@testset "Model Warp preserves periodic links" begin
    session() do
        A.option("Mesh.Renumber",0)
        periodic_box()
        A.mesh.set_periodic(2,[4],[6],AFF_X)
        A.mesh.generate(3)
        relation=A.mesh.get_periodic_nodes(2,4)
        @test !isempty(relation.slave_nodes)
        tags=A.mesh.get_nodes()[1]
        scalar=A.view.add("scalar")
        A.view.add_model_data(scalar,0,"", "NodeData",tags,
            [[0.] for _ in tags],0.,1)
        vector=A.view.add("vector")
        A.view.add_model_data(vector,0,"", "NodeData",tags,
            repeat([1.,0.,0.],length(tags)),0.,3)
        A.plugin.set_number("Warp","View",scalar)
        A.plugin.set_number("Warp","OtherView",vector)
        A.plugin.run("Warp")
        updated=A.mesh.get_periodic_nodes(2,4)
        @test updated.slave_nodes==relation.slave_nodes
        @test updated.master_nodes==relation.master_nodes
        @test updated.affine==relation.affine
    end
end

@testset "Malformed Model Warp is atomic" begin
    session() do
        scalar,vector=fixture()
        A.plugin.set_number("Warp","View",scalar)
        A.plugin.set_number("Warp","OtherView",vector)
        before=coords()
        A.plugin.set_number("Warp","TimeStep",5.)
        @test_throws ArgumentError A.plugin.run("Warp")
        @test coords()==before
        A.plugin.set_number("Warp","TimeStep",1.5)
        @test_throws ArgumentError A.plugin.run("Warp")
        @test coords()==before
        A.plugin.set_number("Warp","TimeStep",0.)
        bad=A.view.add("bad")
        A.view.add_model_data(bad,0,"warp","NodeData",[999],
            [1.,2.,3.],0.,3)
        A.plugin.set_number("Warp","OtherView",bad)
        @test_throws ArgumentError A.plugin.run("Warp")
        @test coords()==before
        # Component count other than three is rejected before any move.
        short=A.view.add("short")
        A.view.add_model_data(short,0,"warp","NodeData",[1,2,3,4],
            [1.,2.,3.,4.],0.,1)
        A.plugin.set_number("Warp","OtherView",short)
        @test_throws ArgumentError A.plugin.run("Warp")
        @test coords()==before
    end
end
end
