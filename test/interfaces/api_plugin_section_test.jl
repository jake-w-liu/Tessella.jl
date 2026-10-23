# `plugin` module UI01: Gradient (per-node isoparametric gradient of scalar
# data) and CutPlane (`A x + B y + C z + D = 0` levelset sections plus the
# ExtractVolume whole-element filter) on small list and model view fixtures.
module PluginSectionChecks
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

@testset "Gradient plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # Scalar triangle: u = x + 2y + 5 on (0,0),(1,0),(0,1) → ∇u = (1,2,0)
        # replicated at every node (linear data: constant gradient).
        tri=A.view.add("gtri")
        A.view.add_list_data(tri,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 5.,6.,7.])
        A.plugin.set_number("Gradient","View",tri)
        A.plugin.run("Gradient")
        types,data=A.view.get_list_data(latest())
        @test types==["VT"]
        @test data[1]≈[0.,1.,0., 0.,0.,1., 0.,0.,0.,
                       1.,2.,0., 1.,2.,0., 1.,2.,0.]

        # Line along z with u = 3z: tangential gradient (0,0,3).
        lin=A.view.add("glin")
        A.view.add_list_data(lin,"SL",1,
            [0.,0., 0.,0., 0.,2., 0.,6.])
        A.plugin.set_number("Gradient","View",lin)
        A.plugin.run("Gradient")
        types,data=A.view.get_list_data(latest())
        @test types==["VL"]
        @test data[1][7:12]≈[0.,0.,3., 0.,0.,3.]

        # Tetrahedron: u = x + 2y + 3z → (1,2,3) at all four corners.
        tet=A.view.add("gtet")
        A.view.add_list_data(tet,"SS",1,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1., 0.,1.,2.,3.])
        A.plugin.set_number("Gradient","View",tet)
        A.plugin.run("Gradient")
        types,data=A.view.get_list_data(latest())
        @test types==["VS"]
        for node in 0:3
            @test data[1][13+3node:15+3node]≈[1.,2.,3.]
        end

        # Bilinear quad on the unit square, u = x*y: the per-node gradient
        # varies — (0,0), (0,1), (1,1), (1,0) in the corners' xy-plane.
        qua=A.view.add("gqua")
        A.view.add_list_data(qua,"SQ",1,
            [0.,1.,1.,0., 0.,0.,1.,1., 0.,0.,0.,0., 0.,0.,1.,0.])
        A.plugin.set_number("Gradient","View",qua)
        A.plugin.run("Gradient")
        types,data=A.view.get_list_data(latest())
        @test types==["VQ"]
        @test data[1][13:24]≈[0.,0.,0., 0.,1.,0., 1.,1.,0., 1.,0.,0.]

        # Two steps differentiate independently: step1 u = x, step2 u = 2y.
        two=A.view.add("gtwo")
        A.view.add_list_data(two,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,1.,0., 0.,0.,2.])
        A.plugin.set_number("Gradient","View",two)
        A.plugin.run("Gradient")
        types,data=A.view.get_list_data(latest())
        @test types==["VT"]
        @test data[1][10:18]≈[1.,0.,0., 1.,0.,0., 1.,0.,0.]
        @test data[1][19:27]≈[0.,2.,0., 0.,2.,0., 0.,2.,0.]

        # Non-scalar records and degenerate/point elements are skipped.
        vec=A.view.add("gvec")
        A.view.add_list_data(vec,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,0.,0., 0.,1.,0., 0.,0.,1.])
        A.plugin.set_number("Gradient","View",vec)
        A.plugin.run("Gradient")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)
        flat=A.view.add("gflat")
        A.view.add_list_data(flat,"ST",1,
            [0.,2.,1., 0.,0.,0., 0.,0.,0., 1.,2.,3.])
        A.plugin.set_number("Gradient","View",flat)
        A.plugin.run("Gradient")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)
        pnt=A.view.add("gpnt")
        A.view.add_list_data(pnt,"SP",1,[0.5,0.5,0.5, 7.0])
        A.plugin.set_number("Gradient","View",pnt)
        A.plugin.run("Gradient")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)

        # Registry plumbing: unknown options and missing views are errors.
        @test_throws ArgumentError A.plugin.set_number("Gradient","Bogus",1.0)
        @test_throws ArgumentError A.plugin.set_string("Gradient","View","x")
        @test_throws ArgumentError A.plugin.run("Gradient2")
    finally
        A.finalize()
    end
end

@testset "CutPlane plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # Triangle cut by x = 0.5: crossings on edges (1,2) at (0.5,0) and
        # (2,3) at (0.5,0.5) emit a line with interpolated values 2 and 4.
        tri=A.view.add("ctri")
        A.view.add_list_data(tri,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        for (opt,v) in (("A",1.0),("B",0.0),("C",0.0),("D",-0.5),
                        ("View",Float64(tri)))
            A.plugin.set_number("CutPlane",opt,v)
        end
        A.plugin.run("CutPlane")
        types,data=A.view.get_list_data(latest())
        @test types==["SL"]
        @test data[1]≈[0.5,0.5, 0.,0.5, 0.,0., 2.,4.]

        # Tetrahedron cut by z = 0.5: three midpoints on the apex edges form
        # a section triangle with interpolated corner values.
        tet=A.view.add("ctet")
        A.view.add_list_data(tet,"SS",1,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1., 1.,2.,3.,4.])
        for (opt,v) in (("A",0.0),("B",0.0),("C",1.0),("D",-0.5),
                        ("View",Float64(tet)))
            A.plugin.set_number("CutPlane",opt,v)
        end
        A.plugin.run("CutPlane")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        got=Set{NTuple{4,Float64}}()
        rec=data[1]
        for node in 0:2
            push!(got,(rec[1+node],rec[4+node],rec[7+node],rec[10+node]))
        end
        @test got==Set([(0.0,0.0,0.5,2.5),(0.5,0.0,0.5,3.0),(0.0,0.5,0.5,3.5)])

        # Unit hexahedron cut by z = 0.5: upstream decomposes the hex into
        # six tetrahedra whose sections emit 4 ST + 2 SQ records tiling the
        # unit square exactly once (verified byte-identical to gmsh 4.15.2).
        hex=A.view.add("chex")
        A.view.add_list_data(hex,"SH",1,
            [0.,1.,1.,0., 0.,1.,1.,0.,   # x1..x8
             0.,0.,1.,1., 0.,0.,1.,1.,   # y
             0.,0.,0.,0., 1.,1.,1.,1.,   # z
             zeros(4)..., ones(4)...])
        A.plugin.set_number("CutPlane","View",hex)
        A.plugin.run("CutPlane")
        types,data=A.view.get_list_data(latest())
        @test types==["ST","SQ"]
        tris,quads=data
        @test length(tris)==4*12 && length(quads)==2*16
        area=0.0
        for element in 0:3
            rec=tris[element*12+1:element*12+12]
            pts=[(rec[i],rec[3+i],rec[6+i]) for i in 1:3]
            for p in pts
                @test p[3]≈0.5
            end
            e1=(pts[2].-pts[1]);e2=(pts[3].-pts[1])
            area+=0.5hypot(e1[2]*e2[3]-e1[3]*e2[2],
                           e1[3]*e2[1]-e1[1]*e2[3],
                           e1[1]*e2[2]-e1[2]*e2[1])
            @test rec[10:12]≈fill(0.5,3)
        end
        for element in 0:1
            rec=quads[element*16+1:element*16+16]
            pts=[(rec[i],rec[4+i],rec[8+i]) for i in 1:4]
            for p in pts
                @test p[3]≈0.5
            end
            e1=(pts[3].-pts[1]);e2=(pts[4].-pts[2])
            area+=0.5abs(e1[1]*e2[2]-e1[2]*e2[1])
            @test rec[13:16]≈fill(0.5,4)
        end
        @test area≈1.0

        # A plane touching a single vertex (x + y = 0 hits only node 1, whose
        # incident edges both report it) emits nothing upstream — sections
        # below the element's codimension are discarded.
        apex=A.view.add("capex")
        A.view.add_list_data(apex,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 9.,3.,5.])
        for (opt,v) in (("A",1.0),("B",1.0),("C",0.0),("D",0.0),
                        ("View",Float64(apex)))
            A.plugin.set_number("CutPlane",opt,v)
        end
        A.plugin.run("CutPlane")
        types,data=A.view.get_list_data(latest())
        @test isempty(types)

        # A triangle lying in the cut plane emits its first edge as an `SL`
        # record (upstream's all-on-levelset degenerate case).
        flat=A.view.add("cflat")
        A.view.add_list_data(flat,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 4.,5.,6.])
        for (opt,v) in (("A",0.0),("B",0.0),("C",1.0),("D",0.0),
                        ("View",Float64(flat)))
            A.plugin.set_number("CutPlane",opt,v)
        end
        A.plugin.run("CutPlane")
        types,data=A.view.get_list_data(latest())
        @test types==["SL"]
        @test data[1]≈[0.,1., 0.,0., 0.,0., 4.,5.]

        # ExtractVolume keeps whole elements on the selected side: the
        # straddling tetrahedron goes to neither side, the low one to -1.
        pair=A.view.add("cpair")
        A.view.add_list_data(pair,"SS",2,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1., 1.,2.,3.,4.,
             0.,1.,0.,0., 0.,0.,1.,0., -2.,-2.,-2.,-1., 1.,2.,3.,4.])
        for (opt,v) in (("A",0.0),("B",0.0),("C",1.0),("D",-0.5),
                        ("ExtractVolume",-1.0),("View",Float64(pair)))
            A.plugin.set_number("CutPlane",opt,v)
        end
        A.plugin.run("CutPlane")
        types,data=A.view.get_list_data(latest())
        @test types==["SS"]
        @test data[1]≈[0.,1.,0.,0., 0.,0.,1.,0., -2.,-2.,-2.,-1., 1.,2.,3.,4.]
        A.plugin.set_number("CutPlane","ExtractVolume",1.0)
        A.plugin.run("CutPlane")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)
        A.plugin.set_number("CutPlane","ExtractVolume",0.0)

        # Two-step data interpolates each step independently.
        two=A.view.add("ctwo")
        A.view.add_list_data(two,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 2.,6.,4., 20.,60.,40.])
        for (opt,v) in (("A",1.0),("B",0.0),("C",0.0),("D",-0.5),
                        ("View",Float64(two)))
            A.plugin.set_number("CutPlane",opt,v)
        end
        A.plugin.run("CutPlane")
        types,data=A.view.get_list_data(latest())
        @test types==["SL"]
        @test data[1][7:8]≈[4.,5.]
        @test data[1][9:10]≈[40.,50.]

        # A degenerate plane normal rejects loudly.
        for (opt,v) in (("A",0.0),("B",0.0),("C",0.0))
            A.plugin.set_number("CutPlane",opt,v)
        end
        @test occursin("plane normal",_errtext(()->A.plugin.run("CutPlane")))
        @test_throws ArgumentError A.plugin.set_number("CutPlane","Bogus",1.0)
        @test occursin("unknown view tag",
            _errtext(()->(A.plugin.set_number("CutPlane","View",77);
                          A.plugin.set_number("CutPlane","A",1.0);
                          A.plugin.run("CutPlane"))))
    finally
        A.finalize()
    end
end

@testset "Gradient and CutPlane on model data" begin
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

        # NodeData u = x: the gradient is the orthogonal projection of x̂
        # onto each element's tangent space — (1,0,0) on y/z-normal faces and
        # edges along x, zero on x-normal faces and y/z edges.
        field=A.view.add("xcoord")
        rows=[[node_xyz[3(i-1)+1]] for (i,_) in enumerate(node_ids)]
        A.view.add_model_data(field,0,"box","NodeData",node_ids,rows)
        A.plugin.set_number("Gradient","View",field)
        A.plugin.run("Gradient")
        types,data=A.view.get_list_data(latest())
        @test sort(types)==["VL","VT"]
        block=Dict(type=>data[i] for (i,type) in enumerate(types))
        let vec=block["VL"],arity=2
            for element in 0:(length(vec)÷(3arity+3arity))-1
                base=element*(3arity+3arity)
                e=(vec[base+2]-vec[base+1],
                   vec[base+arity+2]-vec[base+arity+1],
                   vec[base+2arity+2]-vec[base+2arity+1])
                s=hypot(e...)
                s>0 || continue
                ex=e[1]/s
                for node in 0:arity-1
                    g=vec[base+3arity+3node+1:base+3arity+3node+3]
                    @test g≈[ex*ex,ex*e[2]/s,ex*e[3]/s] atol=1e-10
                end
            end
        end
        let vec=block["VT"],arity=3
            for element in 0:(length(vec)÷(3arity+3arity))-1
                base=element*(3arity+3arity)
                p1=(vec[base+1],vec[base+arity+1],vec[base+2arity+1])
                e1=(vec[base+2]-p1[1],vec[base+arity+2]-p1[2],
                    vec[base+2arity+2]-p1[3])
                e2=(vec[base+3]-p1[1],vec[base+arity+3]-p1[2],
                    vec[base+2arity+3]-p1[3])
                n=(e1[2]*e2[3]-e1[3]*e2[2],e1[3]*e2[1]-e1[1]*e2[3],
                   e1[1]*e2[2]-e1[2]*e2[1])
                s=hypot(n...)
                s>0 || continue
                nx,ny,nz=n./s
                expected=(1-nx*nx,-nx*ny,-nx*nz)
                for node in 0:arity-1
                    g=vec[base+3arity+3node+1:base+3arity+3node+3]
                    @test g≈collect(expected) atol=1e-10
                end
            end
        end

        # CutPlane x = 0.5: every section point and interpolated value sits
        # on the plane.
        for (opt,v) in (("A",1.0),("B",0.0),("C",0.0),("D",-0.5),
                        ("View",Float64(field)))
            A.plugin.set_number("CutPlane",opt,v)
        end
        A.plugin.run("CutPlane")
        types,data=A.view.get_list_data(latest())
        @test !isempty(types)
        block=Dict(type=>data[i] for (i,type) in enumerate(types))
        for (type,arity) in (("SP",1),("SL",2),("ST",3))
            haskey(block,type) || continue
            vec=block[type];stride=3arity+arity
            for element in 0:(length(vec)÷stride)-1
                base=element*stride
                for node in 0:arity-1
                    @test vec[base+node+1]≈0.5 atol=1e-10
                    @test vec[base+3arity+node+1]≈0.5 atol=1e-10
                end
            end
        end
        @test haskey(block,"SL") || haskey(block,"ST")
    finally
        A.finalize()
    end
end

end # module PluginSectionChecks
