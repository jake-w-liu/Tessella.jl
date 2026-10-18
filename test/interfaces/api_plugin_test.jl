# `plugin` module UI01: option registry validation plus real plugin
# computations on small session-view fixtures (list data and model data).
module PluginChecks
using Test, Tessella
const A=Tessella.API

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "plugin registry validation" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # Unknown plugin names and options raise ArgumentError.
        @test_throws ArgumentError A.plugin.run("NoSuchPlugin")
        @test_throws ArgumentError A.plugin.set_number("MinMax","Bogus",1.0)
        @test_throws ArgumentError A.plugin.set_string("NewView","Bogus","x")
        @test occursin("unknown plugin",_errtext(()->A.plugin.run("NoSuchPlugin")))
        @test occursin("unknown numeric option",
            _errtext(()->A.plugin.set_number("MinMax","Bogus",1.0)))
        @test occursin("unknown string option",
            _errtext(()->A.plugin.set_string("NewView","Bogus","x")))
        # Numeric options reject set_string and vice versa.
        @test_throws ArgumentError A.plugin.set_string("MinMax","View","0")
        @test_throws ArgumentError A.plugin.set_number("NewView","Type",1.0)
        # Value type/range validation follows the API's numeric conventions.
        @test_throws ArgumentError A.plugin.set_number("MinMax","View","x")
        @test_throws ArgumentError A.plugin.set_number("MinMax","View",true)
        @test_throws ArgumentError A.plugin.set_number("MinMax","View",Inf)
        @test_throws ArgumentError A.plugin.set_string("NewView","Type",3)
        @test_throws ArgumentError A.plugin.run(3)
        @test_throws ArgumentError A.plugin.set_number("MinMax",3,1.0)
        @test_throws ArgumentError A.plugin.set_string("NewView","Type",1)
        # Running on an empty session reports the missing view.
        @test occursin("no views exist",_errtext(()->A.plugin.run("MinMax")))
        # Unknown view tags are rejected.
        A.plugin.set_number("MinMax","View",77)
        @test occursin("unknown view tag 77",_errtext(()->A.plugin.run("MinMax")))
    finally
        A.finalize()
    end
end

@testset "list-data plugins compute real output" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # One scalar triangle: nodes (0,0,0),(1,0,0),(0,1,0), values 1,3,5.
        field=A.view.add("field")
        A.view.add_list_data(field,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        latest()=last(A.view.get_tags())

        # MinMax produces *_Min and *_Max scalar-point views with the extrema.
        A.plugin.set_number("MinMax","View",field)
        A.plugin.run("MinMax")
        tags=A.view.get_tags()
        min_tag=tags[end-1];max_tag=tags[end]
        types,data=A.view.get_list_data(min_tag)
        @test types==["SP"]
        @test data[1][1:3]≈[0.5,0.5,0.0]
        @test data[1][4]≈1.0
        types,data=A.view.get_list_data(max_tag)
        @test types==["SP"]
        @test data[1][4]≈5.0
        # Argument=1 anchors each record at the extremum location.
        A.plugin.set_number("MinMax","Argument",1)
        A.plugin.run("MinMax")
        tags=A.view.get_tags()
        arg_min=tags[end-1];arg_max=tags[end]
        _,data=A.view.get_list_data(arg_min)
        @test data[1]≈[0.,0.,0.,1.]
        _,data=A.view.get_list_data(arg_max)
        @test data[1]≈[0.,1.,0.,5.]
        A.plugin.set_number("MinMax","Argument",0)

        # Probe interpolates at (0.25,0.25,0): N=(0.5,0.25,0.25).
        A.plugin.set_number("Probe","X",0.25)
        A.plugin.set_number("Probe","Y",0.25)
        A.plugin.set_number("Probe","Z",0.0)
        A.plugin.set_number("Probe","View",field)
        A.plugin.run("Probe")
        probe_tag=latest()
        types,data=A.view.get_list_data(probe_tag)
        @test types==["SP"]
        @test data[1]≈[0.25,0.25,0.0,2.5]

        # Integrate: mean(1,3,5)*area(0.5)=1.5.
        A.plugin.set_number("Integrate","View",field)
        A.plugin.run("Integrate")
        integ_tag=latest()
        types,data=A.view.get_list_data(integ_tag)
        @test types==["SP"]
        @test data[1][1:3]≈[0.5,0.5,0.0]
        @test data[1][4]≈1.5

        # Skin of a single triangle is its 3 edges (scalar line records).
        A.plugin.set_number("Skin","View",field)
        A.plugin.run("Skin")
        skin_tag=latest()
        types,data=A.view.get_list_data(skin_tag)
        @test types==["SL"]
        # Each SL record: x1 x2 y1 y2 z1 z2 v1 v2.
        @test length(data[1])==3*8
        @test sort([data[1][7],data[1][8],data[1][15],data[1][16],
                    data[1][23],data[1][24]])≈[1.,1.,3.,3.,5.,5.]

        # ExtractElements keeps the triangle when mean(1,3,5)=3 is inside
        # [0,4) and drops it for [4,10).
        A.plugin.set_number("ExtractElements","View",field)
        A.plugin.set_number("ExtractElements","MinVal",0.0)
        A.plugin.set_number("ExtractElements","MaxVal",4.0)
        A.plugin.run("ExtractElements")
        keep_tag=latest()
        types,data=A.view.get_list_data(keep_tag)
        @test types==["ST"]
        @test data[1]≈[0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,3.,5.]
        A.plugin.set_number("ExtractElements","MinVal",4.0)
        A.plugin.set_number("ExtractElements","MaxVal",10.0)
        A.plugin.run("ExtractElements")
        drop_tag=latest()
        types,_=A.view.get_list_data(drop_tag)
        @test isempty(types)

        # Transform translates the element coordinates in place.
        movable=A.view.add("movable")
        A.view.add_list_data(movable,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        A.plugin.set_number("Transform","Tx",10.0)
        A.plugin.set_number("Transform","View",movable)
        A.plugin.run("Transform")
        _,data=A.view.get_list_data(movable)
        @test data[1][1:3]≈[10.,11.,10.]
        @test data[1][10:12]≈[1.,3.,5.]

        # Remove drops the selected (shape, component) classes in place.
        mixed=A.view.add("mixed")
        A.view.add_list_data(mixed,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        A.view.add_list_data(mixed,"SL",1,
            [0.,1., 0.,0., 0.,0., 2.,4.])
        A.plugin.set_number("Remove","Triangles",1.0)
        A.plugin.set_number("Remove","View",mixed)
        A.plugin.run("Remove")
        types,data=A.view.get_list_data(mixed)
        @test types==["SL"]
        @test data[1]≈[0.,1.,0.,0.,0.,0.,2.,4.]

        # Smooth averages coincident node values across elements.
        shared=A.view.add("shared")
        A.view.add_list_data(shared,"ST",2,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,2.,4.,
             1.,1.,0., 0.,1.,1., 0.,0.,0., 10.,20.,30.])
        A.plugin.set_number("Smooth","View",shared)
        A.plugin.run("Smooth")
        _,data=A.view.get_list_data(shared)
        # Shared (1,0,0) node: (2+10)/2=6; shared (0,1,0): (4+30)/2=17.
        @test data[1][10:12]≈[0.,6.,17.]
        @test data[1][22:24]≈[6.,20.,17.]

        # The "Run" option executes the plugin immediately.
        A.plugin.set_number("MinMax","View",field)
        A.plugin.set_number("MinMax","Run",1)
        @test length(A.view.get_tags())>0
    finally
        A.finalize()
    end
end

@testset "model-data plugins and NewView" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # A coarse unit-box surface mesh supplies real element geometry.
        A.model.add("box")
        A.model.add_box(0,0,0,1,1,1;tag=1)
        A.option("Mesh.MeshSizeMin",1.0)
        A.option("Mesh.MeshSizeMax",1.0)
        A.mesh.generate(2)
        node_ids,node_xyz,_=A.mesh.get_nodes()
        types,labels,_=A.mesh.get_elements(2,-1)
        @test !isempty(node_ids) && !isempty(labels)

        # NodeData view: MinMax + Probe vs the session probe.
        nodal=A.view.add("nodal")
        nodal_rows=[[Float64(i)] for (i,_) in enumerate(node_ids)]
        A.view.add_model_data(nodal,0,"box","NodeData",
            node_ids,nodal_rows)
        A.plugin.set_number("MinMax","View",nodal)
        A.plugin.run("MinMax")
        tags=A.view.get_tags()
        _,data=A.view.get_list_data(tags[end-1])
        @test data[1][4]≈1.0
        _,data=A.view.get_list_data(tags[end])
        @test data[1][4]≈Float64(length(node_ids))

        # Probe at a node position returns that node's value exactly.
        px,py,pz=node_xyz[1],node_xyz[2],node_xyz[3]
        expected=A.view.probe("box",nodal,px,py,pz)
        A.plugin.set_number("Probe","X",px)
        A.plugin.set_number("Probe","Y",py)
        A.plugin.set_number("Probe","Z",pz)
        A.plugin.set_number("Probe","View",nodal)
        A.plugin.run("Probe")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SP"]
        @test data[1][4:3+length(expected)]≈expected
        @test data[1][4]≈1.0

        # ElementData view: constant per-element value 2 integrates to
        # 2 * total surface area = 2 * 6.0 = 12.0.
        etags=UInt64[tag for block in labels for tag in block]
        elem=A.view.add("elemental")
        A.view.add_model_data(elem,0,"box","ElementData",
            etags,[Float64[2.0] for _ in etags])
        A.plugin.set_number("Integrate","View",elem)
        A.plugin.run("Integrate")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SP"]
        @test data[1][4]≈12.0

        # NewView fills a constant NodeData view over the mesh nodes.
        A.plugin.set_number("NewView","Value",2.5)
        A.plugin.run("NewView")
        new_tag=last(A.view.get_tags())
        dtype,ntags,ndata=A.view.get_model_data(new_tag,0)
        @test dtype==["NodeData"]
        @test sort(ntags)==sort(node_ids)
        @test Set(ndata)==Set([2.5])
        # The Type string option validates like upstream.
        A.plugin.set_string("NewView","Type","ElementData")
        @test occursin("unknown data type",
            _errtext(()->A.plugin.run("NewView")))
        A.plugin.set_string("NewView","Type","NodeData")

        # Session reset restores plugin defaults: after clear() the view
        # option no longer points at the removed views.
        A.plugin.set_number("MinMax","View",nodal)
        A.clear()
        @test occursin("no views exist",_errtext(()->A.plugin.run("MinMax")))
    finally
        A.finalize()
    end
end

end # module PluginChecks
