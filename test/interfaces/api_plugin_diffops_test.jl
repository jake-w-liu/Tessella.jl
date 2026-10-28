# `plugin` module UI01: Divergence/Curl corner-gradient differential
# operators, Lambda2 eigenvalues of S*S+Omega*Omega, SphericalRaise in-place
# radial displacement and Annotate text records on small list fixtures.
module PluginDiffOpsChecks
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

@testset "Divergence plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # Linear vector field v = (x, y, 0) on a unit right triangle:
        # div = 1 + 1 + 0 = 2 at every corner node.
        tri=A.view.add("tri")
        A.view.add_list_data(tri,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             0.,0.,0., 1.,0.,0., 0.,1.,0.])
        A.plugin.set_number("Divergence","View",tri)
        A.plugin.run("Divergence")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        @test data[1]≈[0.,1.,0., 0.,0.,1., 0.,0.,0., 2.,2.,2.]
        @test A.view.get_string(latest(),"Name")=="tri_Divergence"

        # v = (x, 2y, 3z) on a tetrahedron: div = 6.
        tet=A.view.add("tet")
        A.view.add_list_data(tet,"VS",1,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
             0.,0.,0., 1.,0.,0., 0.,2.,0., 0.,0.,3.])
        A.plugin.set_number("Divergence","View",tet)
        A.plugin.run("Divergence")
        types,data=A.view.get_list_data(latest())
        @test types==["SS"]
        @test data[1][13:16]≈[6.,6.,6.,6.]

        # Two steps differentiate independently: v=(x,y,0) then v=(2x,2y,0)
        # give div 2 then 4.
        tri2=A.view.add("tri2")
        A.view.add_list_data(tri2,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             0.,0.,0., 1.,0.,0., 0.,1.,0.,
             0.,0.,0., 2.,0.,0., 0.,2.,0.])
        A.plugin.set_number("Divergence","View",tri2)
        A.plugin.run("Divergence")
        _,data=A.view.get_list_data(latest())
        @test data[1][10:15]≈[2.,2.,2.,4.,4.,4.]

        # Tangential divergence on a line: v = (x,0,0) along [0,0,0]-[2,0,0]
        # gives t.dv/ds = 1.
        lin=A.view.add("lin")
        A.view.add_list_data(lin,"VL",1,
            [0.,2., 0.,0., 0.,0., 0.,0.,0., 2.,0.,0.])
        A.plugin.set_number("Divergence","View",lin)
        A.plugin.run("Divergence")
        types,data=A.view.get_list_data(latest())
        @test types==["SL"]
        @test data[1][7:8]≈[1.,1.]

        # A vector point record emits the upstream zero gradient.
        pnt=A.view.add("pnt")
        A.view.add_list_data(pnt,"VP",1,[1., 2., 3., 4.,5.,6.])
        A.plugin.set_number("Divergence","View",pnt)
        A.plugin.run("Divergence")
        types,data=A.view.get_list_data(latest())
        @test types==["SP"]
        @test data[1]≈[1.,2.,3.,0.]

        # Scalar views are skipped wholesale like upstream's numComp gate.
        scalar=A.view.add("scalar")
        A.view.add_list_data(scalar,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,2.,3.])
        A.plugin.set_number("Divergence","View",scalar)
        A.plugin.run("Divergence")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)
    finally
        A.finalize()
    end
end

@testset "Curl plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # Rigid rotation v = (-y, x, 0): curl = (0,0,2) at every node.
        tri=A.view.add("tri")
        A.view.add_list_data(tri,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             0.,0.,0., 0.,1.,0., -1.,0.,0.])
        A.plugin.set_number("Curl","View",tri)
        A.plugin.run("Curl")
        types,data=A.view.get_list_data(latest())
        @test types==["VT"]
        @test data[1][10:18]≈[0.,0.,2., 0.,0.,2., 0.,0.,2.]
        @test A.view.get_string(latest(),"Name")=="tri_Curl"

        # v = (0, 0, y) on a tetrahedron: curl = (1,0,0).
        tet=A.view.add("tet")
        A.view.add_list_data(tet,"VS",1,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
             0.,0.,0., 0.,0.,0., 0.,0.,1., 0.,0.,0.])
        A.plugin.set_number("Curl","View",tet)
        A.plugin.run("Curl")
        types,data=A.view.get_list_data(latest())
        @test types==["VS"]
        for node in 0:3
            @test data[1][13+3node:15+3node]≈[1.,0.,0.]
        end

        # Scalar views produce no output.
        scalar=A.view.add("scalar")
        A.view.add_list_data(scalar,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,2.,3.])
        A.plugin.set_number("Curl","View",scalar)
        A.plugin.run("Curl")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)
    finally
        A.finalize()
    end
end

@testset "Lambda2 plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # Pure rotation v = (-y, x, 0): the velocity gradient is
        # antisymmetric, S = 0 and A = Omega^2 = diag(-1,-1,0), so the
        # descending eigenvalues are (0,-1,-1): lambda2 = -1 (a vortex).
        tri=A.view.add("tri")
        A.view.add_list_data(tri,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             0.,0.,0., 0.,1.,0., -1.,0.,0.])
        A.plugin.set_number("Lambda2","View",tri)
        A.plugin.run("Lambda2")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        @test data[1][10:12]≈[-1.,-1.,-1.]
        @test A.view.get_string(latest(),"Name")=="tri_Lambda2"

        # Pure strain v = (x, -y, 0): S = diag(1,-1,0), A = diag(1,1,0),
        # eigenvalues (1,1,0) -> lambda2 = 1.
        tri2=A.view.add("tri2")
        A.view.add_list_data(tri2,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             0.,0.,0., 1.,0.,0., 0.,-1.,0.])
        A.plugin.set_number("Lambda2","View",tri2)
        A.plugin.run("Lambda2")
        types,data=A.view.get_list_data(latest())
        @test data[1][10:12]≈[1.,1.,1.]

        # 9-component records supply the tensor directly: g = diag(1,4,9)
        # gives A = diag(1,16,81), eigenvalues (81,16,1).
        tp=A.view.add("tp")
        A.view.add_list_data(tp,"TP",1,
            [1., 2., 3., 1.,0.,0., 0.,4.,0., 0.,0.,9.])
        A.plugin.set_number("Lambda2","View",tp)
        A.plugin.run("Lambda2")
        types,data=A.view.get_list_data(latest())
        @test types==["SP"]
        @test data[1]≈[1.,2.,3.,16.]
        A.plugin.set_number("Lambda2","Eigenvalue",1)
        A.plugin.run("Lambda2")
        _,data=A.view.get_list_data(latest())
        @test data[1]≈[1.,2.,3.,81.]
        A.plugin.set_number("Lambda2","Eigenvalue",3)
        A.plugin.run("Lambda2")
        _,data=A.view.get_list_data(latest())
        @test data[1]≈[1.,2.,3.,1.]

        # Out-of-range Eigenvalue logs the upstream error and emits
        # nothing; `Msg::Error` only continues when `General.AbortOnError`
        # is 0.
        A.option("General.AbortOnError",0)
        before=length(A.view.get_tags())
        A.plugin.set_number("Lambda2","Eigenvalue",4)
        A.plugin.run("Lambda2")
        @test length(A.view.get_tags())==before
        A.plugin.set_number("Lambda2","Eigenvalue",2)
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

@testset "SphericalRaise plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)

        # Nodes on the axes at distances 1,2,3 from the origin move along
        # their radial unit vectors by (Offset + Raise*value): values 1,2,3
        # extend the radii to 1+1, 2+2, 3+3.
        tri=A.view.add("tri")
        A.view.add_list_data(tri,"ST",1,
            [1.,0.,0., 0.,2.,0., 0.,0.,3., 1.,2.,3.])
        A.plugin.set_number("SphericalRaise","View",tri)
        A.plugin.run("SphericalRaise")
        _,data=A.view.get_list_data(tri)
        @test data[1]≈[2.,0.,0., 0.,4.,0., 0.,0.,6., 1.,2.,3.]

        # coef = Offset + Raise*value = 1 - value: the first node stays, the
        # others retract along their radial unit vectors.
        tri2=A.view.add("tri2")
        A.view.add_list_data(tri2,"ST",1,
            [1.,0.,0., 0.,2.,0., 0.,0.,3., 1.,2.,3.])
        A.plugin.set_number("SphericalRaise","View",tri2)
        A.plugin.set_number("SphericalRaise","Raise",-1.0)
        A.plugin.set_number("SphericalRaise","Offset",1.0)
        A.plugin.run("SphericalRaise")
        _,data=A.view.get_list_data(tri2)
        @test data[1]≈[1.,0.,0., 0.,1.,0., 0.,0.,1., 1.,2.,3.]
        A.plugin.set_number("SphericalRaise","Raise",1.0)
        A.plugin.set_number("SphericalRaise","Offset",0.0)

        # TimeStep picks the step whose values displace the coordinates:
        # step-1 values (10,20,30) move the unit-radius axis nodes to
        # radii 11, 21, 31.
        tri3=A.view.add("tri3")
        A.view.add_list_data(tri3,"ST",1,
            [1.,0.,0., 0.,1.,0., 0.,0.,1.,
             1.,2.,3., 10.,20.,30.])
        A.plugin.set_number("SphericalRaise","View",tri3)
        A.plugin.set_number("SphericalRaise","TimeStep",1)
        A.plugin.run("SphericalRaise")
        _,data=A.view.get_list_data(tri3)
        @test data[1]≈[11.,0.,0., 0.,21.,0., 0.,0.,31., 1.,2.,3.,10.,20.,30.]

        # An invalid TimeStep logs the upstream error and leaves the
        # view; `Msg::Error` only continues when `General.AbortOnError` is 0.
        A.option("General.AbortOnError",0)
        tri4=A.view.add("tri4")
        A.view.add_list_data(tri4,"ST",1,
            [1.,0.,0., 0.,2.,0., 0.,0.,3., 1.,2.,3.])
        A.plugin.set_number("SphericalRaise","View",tri4)
        A.plugin.set_number("SphericalRaise","TimeStep",5)
        A.plugin.run("SphericalRaise")
        _,data=A.view.get_list_data(tri4)
        @test data[1]≈[1.,0.,0., 0.,2.,0., 0.,0.,3., 1.,2.,3.]
        A.plugin.set_number("SphericalRaise","TimeStep",0)

        # Model views are rejected like the other in-place writers.
        model_view=A.view.add("model")
        A.model.add("box")
        A.model.add_box(0,0,0,1,1,1;tag=1)
        A.option("Mesh.MeshSizeMin",2.0)
        A.option("Mesh.MeshSizeMax",2.0)
        A.mesh.generate(2)
        node_ids,_,_=A.mesh.get_nodes()
        A.view.add_model_data(model_view,0,"box","NodeData",node_ids,
                              [[1.0] for _ in node_ids])
        A.plugin.set_number("SphericalRaise","View",model_view)
        @test occursin("rewrites list data only",
            _errtext(()->A.plugin.run("SphericalRaise")))
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

@testset "Annotate plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # ThreeD=1 anchors the string at the model position (X,Y,Z).
        tri=A.view.add("tri")
        A.view.add_list_data(tri,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,2.,3.])
        A.plugin.set_string("Annotate","Text","hello")
        A.plugin.set_number("Annotate","ThreeD",1)
        A.plugin.set_number("Annotate","X",1.0)
        A.plugin.set_number("Annotate","Y",2.0)
        A.plugin.set_number("Annotate","Z",3.0)
        A.plugin.set_number("Annotate","View",tri)
        A.plugin.run("Annotate")
        types,strings,coords=A.view.get_list_data_strings(tri,1,-1,true)
        @test types==["T3"]
        @test strings==["hello"]
        @test coords≈[1.,2.,3.,0.]

        # A second annotation appends to the same record block.
        A.plugin.set_string("Annotate","Text","again")
        A.plugin.run("Annotate")
        types,strings,coords=A.view.get_list_data_strings(tri,2,-1,true)
        @test types==["T3"]
        @test strings==["hello","again"]
        @test coords≈[1.,2.,3.,0., 1.,2.,3.,0.]

        # ThreeD=0 stores screen coordinates with z = 0 under "T2".
        A.plugin.set_number("Annotate","ThreeD",0)
        A.plugin.run("Annotate")
        types,strings,coords=A.view.get_list_data_strings(tri,1,-1,true)
        @test types==["T2"]
        @test strings==["again"]
        @test coords≈[1.,2.,0.,0.]

        # Model input creates a new `*_Annotate` view carrying the string
        # and the source time table.
        A.model.add("box")
        A.model.add_box(0,0,0,1,1,1;tag=1)
        A.option("Mesh.MeshSizeMin",2.0)
        A.option("Mesh.MeshSizeMax",2.0)
        A.mesh.generate(2)
        node_ids,_,_=A.mesh.get_nodes()
        md=A.view.add("md")
        A.view.add_model_data(md,0,"box","NodeData",node_ids,
                              [[1.0] for _ in node_ids])
        A.plugin.set_number("Annotate","View",md)
        A.plugin.run("Annotate")
        tag=latest()
        @test tag!=md
        @test A.view.get_string(tag,"Name")=="md_Annotate"
        types,strings,_=A.view.get_list_data_strings(tag,1,-1,true)
        @test types==["T2"]
        @test strings==["again"]
    finally
        A.finalize()
    end
end

end # module
