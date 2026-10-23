# `plugin` module UI01: CutSphere sphere-levelset sections, Warp coordinate
# displacement, MathEval per-node expression evaluation and Scal2Vec scalar
# packing on small list/model view fixtures.
module PluginComputeChecks
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

@testset "CutSphere plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # One scalar tetrahedron (0,0,0),(1,0,0),(0,1,0),(0,0,1), values
        # 0,10,20,30. The r=0.5 sphere at the origin cuts the three edges
        # leaving node 1 at t = 0.25 -> points (0.25,0,0),(0,0.25,0),
        # (0,0,0.25) and one ST record; values interpolate to
        # (2.5, 5, 7.5) in angle order.
        tet=A.view.add("tet")
        A.view.add_list_data(tet,"SS",1,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1., 0.,10.,20.,30.])
        for (opt,v) in (("Xc",0.0),("Yc",0.0),("Zc",0.0),("R",0.5),
                        ("View",Float64(tet)))
            A.plugin.set_number("CutSphere",opt,v)
        end
        A.plugin.run("CutSphere")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        rec=data[1]
        @test length(rec)==3*3+3
        pts=[(rec[node+1],rec[3+node+1],rec[6+node+1]) for node in 0:2]
        @test sort!([collect(p) for p in pts])≈
            sort!([[0.25,0.,0.],[0.,0.25,0.],[0.,0.,0.25]])
        @test sort(rec[10:12])≈[2.5,5.0,7.5]
        # Upstream faces the section's normal along the outward radial
        # direction (away from the sphere centre at the origin here), which
        # the (p2-p1)x(p3-p1) winding measures directly.
        e1=(pts[2][1]-pts[1][1],pts[2][2]-pts[1][2],pts[2][3]-pts[1][3])
        e2=(pts[3][1]-pts[1][1],pts[3][2]-pts[1][2],pts[3][3]-pts[1][3])
        n=(e1[2]*e2[3]-e1[3]*e2[2],e1[3]*e2[1]-e1[1]*e2[3],
           e1[1]*e2[2]-e1[2]*e2[1])
        @test n[1]+n[2]+n[3]>0

        # A line through a unit sphere emits one scalar point at the linear
        # interpolation position t = 0.25 -> (0.5,0,0), value 1.5.
        lin=A.view.add("lin")
        A.view.add_list_data(lin,"SL",1,[0.,2., 0.,0., 0.,0., 1.,3.])
        A.plugin.set_number("CutSphere","R",1.0)
        A.plugin.set_number("CutSphere","View",lin)
        A.plugin.run("CutSphere")
        types,data=A.view.get_list_data(latest())
        @test types==["SP"]
        @test data[1]≈[0.5,0.0,0.0,1.5]

        # A triangle entirely inside the sphere cuts nothing: the emitted
        # view stays empty.
        tri=A.view.add("tri")
        A.view.add_list_data(tri,"ST",1,
            [0.,0.1,0., 0.,0.,0.1, 0.,0.,0., 1.,2.,3.])
        A.plugin.set_number("CutSphere","View",tri)
        A.plugin.run("CutSphere")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)

        # ExtractVolume=-1 keeps elements fully inside; +1 keeps outside.
        A.plugin.set_number("CutSphere","View",tri)
        A.plugin.set_number("CutSphere","ExtractVolume",-1)
        A.plugin.run("CutSphere")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        @test data[1][10:12]≈[1.,2.,3.]
        A.plugin.set_number("CutSphere","ExtractVolume",1)
        A.plugin.run("CutSphere")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)
        A.plugin.set_number("CutSphere","ExtractVolume",0)

        # A tet with two steps: section values interpolate per step (back to
        # the r = 0.5 sphere).
        tet2=A.view.add("tet2")
        A.view.add_list_data(tet2,"SS",1,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
             0.,10.,20.,30., 0.,20.,40.,60.])
        A.plugin.set_number("CutSphere","R",0.5)
        A.plugin.set_number("CutSphere","View",tet2)
        A.plugin.run("CutSphere")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        rec=data[1]
        @test length(rec)==9+6
        @test sort(rec[10:12])≈[2.5,5.0,7.5]
        @test sort(rec[13:15])≈[5.0,10.0,15.0]
    finally
        A.finalize()
    end
end

@testset "Warp plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)

        # Scalar triangle with both View and OtherView < 0: the smoothed
        # surface normal (0,0,1) times each node's value displaces z by the
        # node value.
        tri=A.view.add("warped")
        A.view.add_list_data(tri,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        A.plugin.run("Warp")
        _,data=A.view.get_list_data(tri)
        @test data[1]≈[0.,1.,0., 0.,0.,1., 1.,3.,5., 1.,3.,5.]

        # A vector view warped by itself (View >= 0, OtherView < 0): each
        # node translates by its own vector at TimeStep, scaled by Factor.
        vec=A.view.add("vec")
        A.view.add_list_data(vec,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             1.,0.,0., 0.,1.,0., 0.,0.,1.])
        A.plugin.set_number("Warp","Factor",2.0)
        A.plugin.set_number("Warp","View",vec)
        A.plugin.run("Warp")
        _,data=A.view.get_list_data(vec)
        @test data[1]≈[2.,1.,0., 0.,2.,1., 0.,0.,2.,
                       1.,0.,0., 0.,1.,0., 0.,0.,1.]
        A.plugin.set_number("Warp","Factor",1.0)

        # A scalar view warped by a separate same-shape vector view.
        target=A.view.add("target")
        A.view.add_list_data(target,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        field=A.view.add("field")
        A.view.add_list_data(field,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             1.,0.,0., 0.,1.,0., 0.,0.,1.])
        A.plugin.set_number("Warp","View",target)
        A.plugin.set_number("Warp","OtherView",field)
        A.plugin.run("Warp")
        _,data=A.view.get_list_data(target)
        @test data[1]≈[1.,1.,0., 0.,1.,1., 0.,0.,1., 1.,3.,5.]

        # TimeStep selects which step of OtherView supplies the field;
        # element coordinates still move once per target step.
        field2=A.view.add("field2")
        A.view.add_list_data(field2,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             1.,0.,0., 0.,0.,0., 0.,0.,0., 0.,1.,0., 0.,0.,0., 0.,0.,0.])
        target2=A.view.add("target2")
        A.view.add_list_data(target2,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        A.plugin.set_number("Warp","View",target2)
        A.plugin.set_number("Warp","OtherView",field2)
        A.plugin.set_number("Warp","TimeStep",1)
        A.plugin.run("Warp")
        _,data=A.view.get_list_data(target2)
        @test data[1]≈[0.,1.,0., 1.,0.,1., 0.,0.,0., 1.,3.,5.]
        A.plugin.set_number("Warp","TimeStep",0)

        # A mismatched element count logs "Incompatible views" and leaves
        # the target unchanged.
        other=A.view.add("other")
        A.view.add_list_data(other,"VT",2,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,1.,0., 0.,0.,1., 0.,0.,0.,
             0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,1.,0., 0.,0.,1., 0.,0.,0.])
        target3=A.view.add("target3")
        A.view.add_list_data(target3,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        A.plugin.set_number("Warp","View",target3)
        A.plugin.set_number("Warp","OtherView",other)
        A.plugin.run("Warp")
        _,data=A.view.get_list_data(target3)
        @test data[1]≈[0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.]

        # Invalid TimeStep (past OtherView's step count) also aborts
        # without touching the view.
        A.plugin.set_number("Warp","OtherView",field)
        A.plugin.set_number("Warp","TimeStep",4)
        A.plugin.run("Warp")
        _,data=A.view.get_list_data(target3)
        @test data[1]≈[0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.]
        A.plugin.set_number("Warp","TimeStep",0)

        # Non-integer and unknown OtherView tags throw.
        A.plugin.set_number("Warp","OtherView",2.5)
        @test occursin("OtherView tag must be an integer",
            _errtext(()->A.plugin.run("Warp")))
        A.plugin.set_number("Warp","OtherView",99)
        @test occursin("unknown view tag 99",
            _errtext(()->A.plugin.run("Warp")))
    finally
        A.finalize()
    end
end

@testset "MathEval plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # Scalar triangle: v0*2 doubles each node value.
        tri=A.view.add("tri")
        A.view.add_list_data(tri,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        A.plugin.set_string("MathEval","Expression0","v0*2")
        A.plugin.set_number("MathEval","View",tri)
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        @test data[1][10:12]≈[2.,6.,10.]

        # Coordinate variables evaluate per node.
        A.plugin.set_string("MathEval","Expression0","x+2*y")
        A.plugin.run("MathEval")
        _,data=A.view.get_list_data(latest())
        @test data[1][10:12]≈[0.,1.,2.]

        # Default expression on vector data: the Euclidean norm of each
        # node's (v0,v1,v2) -> scalar record.
        vec=A.view.add("vec")
        A.view.add_list_data(vec,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             3.,0.,0., 0.,4.,0., 0.,0.,12.])
        for i in 0:8
            A.plugin.set_string("MathEval","Expression$i","")
        end
        A.plugin.set_string("MathEval","Expression0",
                            "Sqrt(v0^2+v1^2+v2^2)")
        A.plugin.set_number("MathEval","View",vec)
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        @test data[1][10:12]≈[3.,4.,12.]

        # Expression1/2 -> a 3-component record; missing expressions are 0.
        A.plugin.set_string("MathEval","Expression0","v0")
        A.plugin.set_string("MathEval","Expression1","x")
        A.plugin.set_string("MathEval","Expression2","v0+v1")
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(latest())
        @test types==["VT"]
        @test data[1][10:18]≈[3.,0.,3., 0.,1.,4., 0.,0.,0.]

        # Expression3..8 -> a 9-component tensor record padded with zeros.
        for i in 1:3
            A.plugin.set_string("MathEval","Expression$i","")
        end
        A.plugin.set_string("MathEval","Expression0","v0")
        A.plugin.set_string("MathEval","Expression3","v1")
        A.plugin.set_string("MathEval","Expression8","v2")
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(latest())
        @test types==["TT"]
        @test data[1][10:36]≈[3.,0.,0.,0.,0.,0.,0.,0.,0.,
                              0.,0.,0.,4.,0.,0.,0.,0.,0.,
                              0.,0.,0.,0.,0.,0.,0.,0.,12.]

        # w-variables read OtherView at the same (element,node,step) slot
        # when grids match.
        other=A.view.add("other")
        A.view.add_list_data(other,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 10.,20.,30.])
        for i in 0:8
            A.plugin.set_string("MathEval","Expression$i","")
        end
        A.plugin.set_string("MathEval","Expression0","v0+w0")
        A.plugin.set_number("MathEval","View",tri)
        A.plugin.set_number("MathEval","OtherView",other)
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        @test data[1][10:12]≈[11.,23.,35.]

        # ForceInterpolation evaluates w through the containing element
        # instead — node coordinates coincide here, so the result matches.
        A.plugin.set_number("MathEval","ForceInterpolation",1)
        A.plugin.run("MathEval")
        _,data=A.view.get_list_data(latest())
        @test data[1][10:12]≈[11.,23.,35.]

        # Dimension > 0 restricts the interpolation search to that element
        # dimension: no line records exist in `other`, so w stays 0.
        A.plugin.set_string("MathEval","Expression0","w0")
        A.plugin.set_number("MathEval","Dimension",1)
        A.plugin.run("MathEval")
        _,data=A.view.get_list_data(latest())
        @test data[1][10:12]≈[0.,0.,0.]
        A.plugin.set_number("MathEval","Dimension",-1)
        A.plugin.set_number("MathEval","ForceInterpolation",0)
        A.plugin.set_string("MathEval","Expression0","v0+w0")

        # Two steps: TimeStep=-1 emits both, TimeStep=1 emits only step 1.
        tri2=A.view.add("tri2")
        A.view.add_list_data(tri2,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5., 2.,4.,6.])
        A.plugin.set_string("MathEval","Expression0","v0*10")
        A.plugin.set_number("MathEval","View",tri2)
        A.plugin.set_number("MathEval","OtherView",-1)
        A.plugin.run("MathEval")
        _,data=A.view.get_list_data(latest())
        @test data[1][10:15]≈[10.,30.,50.,20.,40.,60.]
        A.plugin.set_number("MathEval","TimeStep",1)
        A.plugin.run("MathEval")
        _,data=A.view.get_list_data(latest())
        @test data[1][10:12]≈[20.,40.,60.]
        A.plugin.set_number("MathEval","TimeStep",-1)

        # PhysicalGroup > 0 on a list view (no entity tags) drops all
        # elements like the upstream entity-less `ge` path.
        A.plugin.set_number("MathEval","PhysicalGroup",2)
        A.plugin.run("MathEval")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)
        A.plugin.set_number("MathEval","PhysicalGroup",-1)

        # Unknown identifiers and malformed expressions raise ArgumentError.
        A.plugin.set_string("MathEval","Expression0","v9")
        @test occursin("unknown identifier",_errtext(()->A.plugin.run("MathEval")))
        A.plugin.set_string("MathEval","Expression0","v0+*")
        @test_throws ArgumentError A.plugin.run("MathEval")
        A.plugin.set_string("MathEval","Expression0","F2")
        @test occursin("unknown identifier",_errtext(()->A.plugin.run("MathEval")))
    finally
        A.finalize()
    end
end

@testset "Scal2Vec plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        coords=[0.,1.,0., 0.,0.,1., 0.,0.,0.]
        vx=A.view.add("vx")
        A.view.add_list_data(vx,"ST",1,[coords;[1.,3.,5.]])
        vy=A.view.add("vy")
        A.view.add_list_data(vy,"ST",1,[coords;[10.,30.,50.]])
        vz=A.view.add("vz")
        A.view.add_list_data(vz,"ST",1,[coords;[100.,300.,500.]])

        # Three scalar views pack into one vector record.
        A.plugin.set_number("Scal2Vec","ViewX",vx)
        A.plugin.set_number("Scal2Vec","ViewY",vy)
        A.plugin.set_number("Scal2Vec","ViewZ",vz)
        A.plugin.set_string("Scal2Vec","NameNewView","packed")
        A.plugin.run("Scal2Vec")
        tag=latest()
        types,data=A.view.get_list_data(tag)
        @test types==["VT"]
        @test data[1]≈[coords;[1.,10.,100.,3.,30.,300.,5.,50.,500.]]
        @test A.view.get_string(tag,"Name")=="packed"

        # Unset components emit zeros; the first set view is the reference.
        A.plugin.set_number("Scal2Vec","ViewY",-1)
        A.plugin.set_number("Scal2Vec","ViewZ",-1)
        A.plugin.run("Scal2Vec")
        types,data=A.view.get_list_data(latest())
        @test types==["VT"]
        @test data[1]≈[coords;[1.,0.,0.,3.,0.,0.,5.,0.,0.]]

        # Two steps in the reference pack both steps.
        vx2=A.view.add("vx2")
        A.view.add_list_data(vx2,"ST",1,[coords;[1.,3.,5.,2.,4.,6.]])
        A.plugin.set_number("Scal2Vec","ViewX",vx2)
        A.plugin.run("Scal2Vec")
        _,data=A.view.get_list_data(latest())
        @test data[1][10:18]≈[1.,0.,0.,3.,0.,0.,5.,0.,0.]
        @test data[1][19:27]≈[2.,0.,0.,4.,0.,0.,6.,0.,0.]

        # No view selected -> upstream "could not find any view" failure.
        A.plugin.set_number("Scal2Vec","ViewX",-1)
        @test occursin("could not find any view",
            _errtext(()->A.plugin.run("Scal2Vec")))
        A.plugin.set_number("Scal2Vec","ViewX",88)
        @test occursin("unknown view tag 88",
            _errtext(()->A.plugin.run("Scal2Vec")))

        # A component view with fewer elements than the reference is a
        # structural mismatch (upstream would read out-of-range), not a
        # silent zero-fill.
        small=A.view.add("small")
        A.view.add_list_data(small,"ST",1,[coords;[2.,4.,6.]])
        big=A.view.add("big")
        A.view.add_list_data(big,"ST",2,
            [coords;[1.,3.,5.];coords;[2.,4.,6.]])
        A.plugin.set_number("Scal2Vec","ViewX",big)
        A.plugin.set_number("Scal2Vec","ViewY",small)
        @test occursin("incompatible element structure",
            _errtext(()->A.plugin.run("Scal2Vec")))
    finally
        A.finalize()
    end
end

@testset "MathEval and Scal2Vec on model data" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("box")
        A.model.add_box(0,0,0,1,1,1;tag=1)
        A.option("Mesh.MeshSizeMin",1.0)
        A.option("Mesh.MeshSizeMax",1.0)
        A.mesh.generate(2)
        node_ids,node_xyz,_=A.mesh.get_nodes()
        @test !isempty(node_ids)
        latest()=last(A.view.get_tags())

        # NodeData u = x+y+z evaluated through MathEval into a list view:
        # per element node, result = v0 (same data) so records echo x+y+z.
        field=A.view.add("modeldata")
        rows=[[node_xyz[3(i-1)+1]+node_xyz[3(i-1)+2]+node_xyz[3(i-1)+3]]
              for (i,_) in enumerate(node_ids)]
        A.view.add_model_data(field,0,"box","NodeData",node_ids,rows)
        for i in 0:8
            A.plugin.set_string("MathEval","Expression$i","")
        end
        A.plugin.set_string("MathEval","Expression0","v0")
        A.plugin.set_number("MathEval","View",field)
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(latest())
        # The generated surface mesh may embed corner point elements.
        @test sort(types)==["SL","SP","ST"]
        block=Dict(type=>data[i] for (i,type) in enumerate(types))
        for (type,arity) in (("SL",2),("SP",1),("ST",3))
            vec=block[type];stride=3arity+arity
            for element in 0:(length(vec)÷stride)-1
                base=element*stride
                for node in 0:arity-1
                    x,y,z=vec[base+node+1],vec[base+arity+node+1],
                          vec[base+2arity+node+1]
                    @test vec[base+3arity+node+1]≈x+y+z atol=1e-10
                end
            end
        end

        # w-variables on a model view interpolate through the mesh locator
        # when OtherView is list-based (different grid -> probe path).
        probe_view=A.view.add("probe_src")
        A.view.add_list_data(probe_view,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 100.,100.,100.])
        A.plugin.set_string("MathEval","Expression0","w0")
        A.plugin.set_number("MathEval","OtherView",probe_view)
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(latest())
        @test !isempty(types)
        A.plugin.set_number("MathEval","OtherView",-1)
        A.plugin.set_string("MathEval","Expression0","v0")

        # Scal2Vec packs two NodeData model views into a vector list view.
        ex=A.view.add("mx")
        A.view.add_model_data(ex,0,"box","NodeData",node_ids,
                              [[node_xyz[3(i-1)+1]] for i in eachindex(node_ids)])
        ey=A.view.add("my")
        A.view.add_model_data(ey,0,"box","NodeData",node_ids,
                              [[node_xyz[3(i-1)+2]] for i in eachindex(node_ids)])
        A.plugin.set_number("Scal2Vec","ViewX",ex)
        A.plugin.set_number("Scal2Vec","ViewY",ey)
        A.plugin.set_number("Scal2Vec","ViewZ",-1)
        A.plugin.run("Scal2Vec")
        types,data=A.view.get_list_data(latest())
        @test sort(types)==["VL","VP","VT"]
        block=Dict(type=>data[i] for (i,type) in enumerate(types))
        for (type,arity) in (("VL",2),("VP",1),("VT",3))
            vec=block[type];stride=3arity+3arity
            for element in 0:(length(vec)÷stride)-1
                base=element*stride
                for node in 0:arity-1
                    x,y=vec[base+node+1],vec[base+arity+node+1]
                    g=vec[base+3arity+3node+1:base+3arity+3node+3]
                    @test g≈[x,y,0.0] atol=1e-10
                end
            end
        end

        # PhysicalGroup filters by the element's owning entity: assign a
        # physical group to one surface and check only its triangles pass.
        # (Adding the group invalidates the current model slot, so this runs
        # last on the view's captured slot.)
        surfaces=A.model.get_entities(2)
        A.model.add_physical_group(2,[surfaces[1][2]];tag=7)
        A.plugin.set_number("MathEval","PhysicalGroup",7)
        A.plugin.set_number("MathEval","View",field)
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(latest())
        @test !isempty(types)
        A.plugin.set_number("MathEval","PhysicalGroup",-1)
    finally
        A.finalize()
    end
end

end # module PluginComputeChecks
