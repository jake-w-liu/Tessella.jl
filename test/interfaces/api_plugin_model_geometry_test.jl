module PluginModelGeometryTests
using Test,Tessella
const A=Tessella.API

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        f()
    finally
        A.finalize()
    end
end

function triangle(;cached=false,highorder=false)
    A.model.add("modelgeometry")
    A.model.add_discrete_entity(2,7)
    tags=highorder ? [11,21,31,41,51,61] : [11,21,31]
    coords=highorder ? [1.,0.,0.,2.,0.,0.,1.,1.,0.,1.5,0.,0.,1.5,.5,0.,1.,.5,0.] :
        [1.,0.,0.,2.,0.,0.,1.,1.,0.]
    A.mesh.add_nodes(2,7,tags,coords)
    A.mesh.add_elements_by_type(7,highorder ? 9 : 2,[301],tags)
    A.model.add_physical_group(2,[7];tag=23)
    if cached
        mktempdir() do dir
            path=joinpath(dir,"cached.msh")
            A.write(path)
            A.clear()
            A.open(path)
        end
    end
    return A.model.get_current(),tags,coords
end

@testset "Model geometry plugins preserve source model, UIDs, values and ownership" begin
    for plugin in ("Transform","SphericalRaise"), cached in (false,true), highorder in (false,true)
        session() do
            name,tags,coordinates=triangle(;cached,highorder)
            view=A.view.add("source")
            A.view.add_model_data(view,0,name,"NodeData",tags,[[Float64(i)] for i in eachindex(tags)])
            values=A.view.get_model_data(view,0)
            elements=A.mesh.get_elements()
            A.model.add("foreign")
            A.plugin.set_number(plugin,"View",view)
            if plugin=="Transform"
                A.plugin.set_number(plugin,"A11",2)
                A.plugin.set_number(plugin,"Tx",.5)
            end
            A.plugin.run(plugin)
            @test A.model.get_current()=="foreign"
            @test A.view.get_model_data(view,0)==values
            A.model.set_current(name)
            got_tags,got,_=A.mesh.get_nodes()
            expected=copy(coordinates)
            expected[1:9]=plugin=="Transform" ? [2.5,0.,0.,4.5,0.,0.,2.5,1.,0.] :
                [2.,0.,0.,4.,0.,0.,3.1213203435596424,3.1213203435596424,0.]
            bytag=Dict(Int(tag)=>got[3i-2:3i] for (i,tag) in enumerate(got_tags))
            @test sort!(Int.(got_tags))==tags
            @test all(bytag[tag]≈expected[3i-2:3i] for (i,tag) in enumerate(tags))
            @test A.mesh.get_elements()==elements
            @test A.model.get_physical_groups()==[(2,23)]
            @test A.model.get_entities_for_physical_group(2,23)==[7]
        end
    end
end

function shared_model(kind,steps)
    A.model.add("sharedgeometry");A.model.add_discrete_entity(2,1)
    A.mesh.add_nodes(2,1,[1,2,3,4],[1.,0.,0.,2.,0.,0.,1.,1.,0.,2.,1.,0.])
    A.mesh.add_elements_by_type(1,2,[11,12],[1,2,3,2,4,3])
    view=A.view.add("shared")
    for step in 0:steps-1
        tags=kind=="NodeData" ? [1,2,3,4] : [11,12]
        data=kind=="NodeData" ? [[1.],[2.],[3.],[4.]] :
            kind=="ElementData" ? [[1.],[2.]] : [[1.,2.,3.],[4.,5.,6.]]
        A.view.add_model_data(view,step,"sharedgeometry",kind,tags,data,Float64(step),1)
    end
    return view
end

@testset "Model Transform follows nodal and incident-element visit semantics" begin
    for kind in ("NodeData","ElementData","ElementNodeData"), steps in (1,2)
        session() do
            view=shared_model(kind,steps)
            before=[A.view.get_model_data(view,step) for step in 0:steps-1]
            A.plugin.set_number("Transform","View",view)
            A.plugin.set_number("Transform","Tx",.5)
            A.plugin.run("Transform")
            repeat=kind=="NodeData" ? 1 : steps
            expected=[1.0 + .5repeat,0.,0.,2.0 + repeat,0.,0.,1.0 + repeat,1.,0.,2.0 + .5repeat,1.,0.]
            kind=="NodeData" && (expected=[1.5,0.,0.,2.5,0.,0.,1.5,1.,0.,2.5,1.,0.])
            @test A.mesh.get_nodes()[2]≈expected
            @test [A.view.get_model_data(view,step) for step in 0:steps-1]==before
            @test A.mesh.get_elements()[3]==[UInt64[1,2,3,2,4,3]]
        end
    end
end

@testset "Model Transform orientation is applied once to nodal elements" begin
    for kind in ("NodeData","ElementData","ElementNodeData"), steps in (1,2)
        session() do
            view=shared_model(kind,steps)
            A.plugin.set_number("Transform","View",view)
            A.plugin.set_number("Transform","SwapOrientation",1)
            A.plugin.run("Transform")
            expected=kind=="NodeData" ? UInt64[1,3,2,2,3,4] : UInt64[1,2,3,2,4,3]
            @test A.mesh.get_elements()[3]==[expected]
        end
    end
end

@testset "Incomplete model rows are skipped and invalid coordinates are atomic" begin
    session() do
        name,tags,coords=triangle()
        view=A.view.add("partial")
        A.view.add_model_data(view,0,name,"NodeData",tags[1:1],[[1.]])
        before=A.mesh.get_nodes()
        for plugin in ("Transform","SphericalRaise")
            A.plugin.set_number(plugin,"View",view)
            A.plugin.run(plugin)
            @test A.mesh.get_nodes()==before
        end
        full=A.view.add("complete")
        A.view.add_model_data(full,0,name,"NodeData",tags,[[1.],[2.],[3.]])
        A.plugin.set_number("Transform","View",full)
        A.plugin.set_number("Transform","A11",1e308)
        @test_throws ArgumentError A.plugin.run("Transform")
        @test A.mesh.get_nodes()==before
        @test A.model.get_current()==name
    end
end

@testset "Model radial raise follows all incident step visits" begin
 for kind in ("NodeData","ElementData","ElementNodeData"),selected in (0,1)
  A.initialize()
  try
   A.option("General.Terminal",0);A.model.add("a");A.model.add_discrete_entity(2,1)
   A.mesh.add_nodes(2,1,[1,2,3],[1.,0.,0.,2.,0.,0.,1.,1.,0.]);A.mesh.add_elements_by_type(1,2,[11],[1,2,3])
   view=A.view.add("v")
   for step in (0,1)
    tags=kind=="NodeData" ? [1,2,3] : [11]
    data=kind=="NodeData" ? [[1.0*10^step],[2.0*10^step],[3.0*10^step]] : kind=="ElementData" ? [[1.0*10^step]] : [[1.0*10^step,2.0*10^step,3.0*10^step]]
    A.view.add_model_data(view,step,"a",kind,tags,data,Float64(step),1)
   end
   A.plugin.set_number("SphericalRaise","View",view);A.plugin.set_number("SphericalRaise","TimeStep",selected);A.plugin.run("SphericalRaise")
   expected=kind=="NodeData" ? [2.,0.,0.,4.,0.,0.,3.1213203435596424,3.1213203435596424,0.] : kind=="ElementData" ? [12.,0.,0.,13.,0.,0.,8.778174593052023,8.778174593052023,0.] : [12.,0.,0.,24.,0.,0.,24.33452377915607,24.33452377915607,0.]
   @test A.mesh.get_nodes()[2]≈expected
  finally;A.finalize();end
 end
end


function two_models(scales=(1.0,2.0))
    for (name,scale) in zip(("a","b"),scales)
        A.model.add(name);A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[1,2,3],[scale,0.,0.,2scale,0.,0.,scale,scale,0.])
        A.mesh.add_elements_by_type(1,2,[11],[1,2,3])
    end
    view=A.view.add("two meshes")
    for (step,name) in enumerate(("a","b"))
        A.view.add_model_data(view,step-1,name,"NodeData",[1,2,3],
            [[1.0*10^(step-1)],[2.0*10^(step-1)],[3.0*10^(step-1)]],Float64(step-1),1)
    end
    return view
end

@testset "Model geometry visits distinct source meshes and stages atomically" begin
    for (plugin,selected) in (("Transform",0),("SphericalRaise",0),("SphericalRaise",1))
        session() do
            view=two_models()
            before_data=[A.view.get_model_data(view,step) for step in (0,1)]
            A.model.add("foreign")
            A.plugin.set_number(plugin,"View",view)
            if plugin=="Transform"
                A.plugin.set_number(plugin,"Tx",.5)
            else
                A.plugin.set_number(plugin,"TimeStep",selected)
            end
            A.plugin.run(plugin)
            @test A.model.get_current()=="foreign"
            @test [A.view.get_model_data(view,step) for step in (0,1)]==before_data
            for name in ("a","b")
                A.model.set_current(name)
                expected=plugin=="Transform" ? (name=="a" ?
                    [1.5,0.,0.,2.5,0.,0.,1.5,1.,0.] : [2.5,0.,0.,4.5,0.,0.,2.5,2.,0.]) :
                    (name=="a" ? [2.,0.,0.,4.,0.,0.,3.1213203435596424,3.1213203435596424,0.] :
                    [12.,0.,0.,24.,0.,0.,23.213203435596423,23.213203435596423,0.])
                @test A.mesh.get_nodes()[2]≈expected
                @test A.mesh.get_nodes()[1]==UInt64[1,2,3]
            end
        end
    end
    session() do
        view=two_models((1e-308,1.0))
        before=Dict{String,Any}()
        for name in ("a","b")
            A.model.set_current(name);before[name]=A.mesh.get_nodes()
        end
        A.plugin.set_number("Transform","View",view)
        A.plugin.set_number("Transform","A11",1e308)
        @test_throws ArgumentError A.plugin.run("Transform")
        @test A.model.get_current()=="b"
        for name in ("a","b")
            A.model.set_current(name);@test A.mesh.get_nodes()==before[name]
        end
        A.option("General.AbortOnError",0)
        A.plugin.set_number("SphericalRaise","View",view)
        A.plugin.set_number("SphericalRaise","TimeStep",4)
        A.plugin.run("SphericalRaise")
        for name in ("a","b")
            A.model.set_current(name);@test A.mesh.get_nodes()==before[name]
        end
    end
end

end
