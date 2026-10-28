# `plugin` module UI01: MakeSimplex in-place list decomposition,
# ModifyComponents expression-driven in-place component rewrites,
# FieldFromAmplitudePhase two-step NodeData synthesis on an imported mesh,
# ExtractEdges sharp/boundary edge records on the current model mesh, and the
# Vorticity Tessella extension (Gmsh 4.15.2 registers Curl, not Vorticity).
module PluginBatch2Checks
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

@testset "MakeSimplex in-place decomposition" begin
    A.initialize()
    try
        A.option("General.Terminal",0)

        # Quadrangle → 2 triangles through the upstream maps (1,2,3)/(1,3,4),
        # with coordinates and per-node values reordered accordingly.
        quad=A.view.add("quad")
        A.view.add_list_data(quad,"SQ",1,
            [0.,1.,1.,0., 0.,0.,1.,1., 0.,0.,0.,0., 1.,3.,5.,7.])
        A.plugin.set_number("MakeSimplex","View",quad)
        A.plugin.run("MakeSimplex")
        types,data=A.view.get_list_data(quad)
        @test types==["ST"]
        @test data[1]≈[0.,1.,1., 0.,0.,1., 0.,0.,0., 1.,3.,5.,
                        0.,1.,0., 0.,1.,1., 0.,0.,0., 1.,5.,7.]

        # Hexahedron → 6 tetrahedra (first map (1,2,4,8), second (1,5,2,8));
        # a pre-existing ST record keeps its payload ahead of the appended
        # decomposed records.
        vol=A.view.add("vol")
        A.view.add_list_data(vol,"ST",1,
            [9.,8.,7., 0.,0.,0., 0.,0.,0., 10.,20.,30.])
        A.view.add_list_data(vol,"SH",1,
            [0.,1.,1.,0.,0.,1.,1.,0., 0.,0.,1.,1.,0.,0.,1.,1.,
             0.,0.,0.,0.,1.,1.,1.,1., Float64.(10:10:80)...])
        A.plugin.set_number("MakeSimplex","View",vol)
        A.plugin.run("MakeSimplex")
        types,data=A.view.get_list_data(vol)
        @test types==["ST","SS"]
        st=data[1]
        @test st[1:12]≈[9.,8.,7.,0.,0.,0.,0.,0.,0.,10.,20.,30.]
        ss=data[2]
        # 6 tetrahedra × (3*4 coords + 4 values) = 96 entries.
        @test length(ss)==6*16
        # First decomposed tetrahedron uses map (1,2,4,8): node values
        # 10, 20, 40, 80.
        @test ss[13:16]≈[10.,20.,40.,80.]
        # Second map (1,5,2,8): values 10, 50, 20, 80.
        @test ss[29:32]≈[10.,50.,20.,80.]

        # Prism → 3 tetrahedra and pyramid → 2 tetrahedra, scalar payloads.
        mix=A.view.add("mix")
        A.view.add_list_data(mix,"SI",1,
            [0.,1.,0.,0.,1.,0., 0.,0.,1.,0.,0.,1., 0.,0.,0.,1.,1.,1.,
             Float64.(1:6)...])
        A.view.add_list_data(mix,"SY",1,
            [0.,1.,1.,0.,0., 0.,0.,1.,1.,0., 0.,0.,0.,0.,1.,
             Float64.(1:5)...])
        A.plugin.set_number("MakeSimplex","View",mix)
        A.plugin.run("MakeSimplex")
        types,data=A.view.get_list_data(mix)
        @test types==["SS"]
        @test length(data[1])==(3+2)*16
        # Prism first map (1,2,3,5): values 1,2,3,5.
        @test data[1][13:16]≈[1.,2.,3.,5.]
        # Pyramid first map (1,2,4,5): values 1,2,4,5 → record 4 of 5.
        @test data[1][3*16+13:3*16+16]≈[1.,2.,4.,5.]

        # Vector records decompose per component (VQ → VT).
        vec=A.view.add("vec")
        A.view.add_list_data(vec,"VQ",1,
            [0.,1.,1.,0., 0.,0.,1.,1., 0.,0.,0.,0.,
             Float64.(1:12)...])
        A.plugin.set_number("MakeSimplex","View",vec)
        A.plugin.run("MakeSimplex")
        types,data=A.view.get_list_data(vec)
        @test types==["VT"]
        # tri1 (nodes 1,2,3) keeps the first three value triples.
        @test data[1][10:18]≈Float64.(1:9)
        # tri2 (nodes 1,3,4): values 1..3, 7..9, 10..12.
        @test data[1][28:36]≈[1.,2.,3.,7.,8.,9.,10.,11.,12.]

        # Lines and points pass through unchanged.
        keep=A.view.add("keep")
        A.view.add_list_data(keep,"SL",1,
            [0.,1., 0.,0., 0.,0., 2.,4.])
        A.plugin.set_number("MakeSimplex","View",keep)
        A.plugin.run("MakeSimplex")
        types,data=A.view.get_list_data(keep)
        @test types==["SL"] && data[1]≈[0.,1.,0.,0.,0.,0.,2.,4.]

        # Model views are rejected like upstream's list-only contract.
        nodal=A.view.add("nodal")
        A.model.add("m")
        A.model.add_box(0,0,0,1,1,1)
        A.option("Mesh.MeshSizeMin",2.0)
        A.option("Mesh.MeshSizeMax",2.0)
        A.mesh.generate(2)
        node_ids,_,_=A.mesh.get_nodes()
        A.view.add_model_data(nodal,0,"m","NodeData",
            node_ids,[[1.0] for _ in node_ids])
        A.plugin.set_number("MakeSimplex","View",nodal)
        @test occursin("list data only",
            _errtext(()->A.plugin.run("MakeSimplex")))
    finally
        A.finalize()
    end
end

@testset "ModifyComponents in-place expressions" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # Out-of-range step requests emit upstream `Msg::Error` ("Number of
        # time steps don't match: using step 0") and continue only when
        # `General.AbortOnError` is 0.
        A.option("General.AbortOnError",0)

        # Default Expression0 = "v0 * Sin(x)" rewrites the scalar component
        # in place; the view tag and record layout are unchanged.
        v=A.view.add("v")
        A.view.add_list_data(v,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        A.plugin.set_number("ModifyComponents","View",v)
        A.plugin.run("ModifyComponents")
        types,data=A.view.get_list_data(v)
        @test types==["ST"]
        @test data[1][10:12]≈[0.0,3sin(1.0),0.0]
        @test A.view.get_number(v,"NbTimeStep")==1

        # A custom expression plus an empty slot on a vector view: only the
        # expressed component changes.
        w=A.view.add("w")
        A.view.add_list_data(w,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,2.,3., 4.,5.,6., 7.,8.,9.])
        A.plugin.set_string("ModifyComponents","Expression0","w0 * v1")
        A.plugin.set_string("ModifyComponents","Expression1","")
        A.plugin.set_string("ModifyComponents","Expression2","")
        A.plugin.set_number("ModifyComponents","View",w)
        A.plugin.set_number("ModifyComponents","OtherView",w)
        A.plugin.run("ModifyComponents")
        types,data=A.view.get_list_data(w)
        @test types==["VT"]
        # comp0 = w0*v1 = v0*v1 (OtherView=self, same step): node-wise
        @test data[1][10:18]≈[2., 2., 3., 20., 5., 6., 56., 8., 9.]

        # OtherTimeStep: w reads the other view's step 0 while step 1 runs.
        two=A.view.add("two")
        A.view.add_list_data(two,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5., 2.,4.,6.])
        A.plugin.set_string("ModifyComponents","Expression0","w0")
        A.plugin.set_number("ModifyComponents","View",two)
        A.plugin.set_number("ModifyComponents","OtherView",two)
        A.plugin.set_number("ModifyComponents","OtherTimeStep",0)
        A.plugin.set_number("ModifyComponents","TimeStep",1)
        A.plugin.run("ModifyComponents")
        _,data=A.view.get_list_data(two)
        # Only step 1 (0-based index 1) changed; it now holds step 0 values.
        @test data[1][10:12]≈[1.,3.,5.]
        @test data[1][13:15]≈[1.,3.,5.]

        # TimeStep variable: on step 1 the values become the step index.
        A.plugin.set_string("ModifyComponents","Expression0","TimeStep")
        A.plugin.set_number("ModifyComponents","TimeStep",1)
        A.plugin.set_number("ModifyComponents","OtherTimeStep",-1)
        A.plugin.run("ModifyComponents")
        _,data=A.view.get_list_data(two)
        @test data[1][13:15]≈[1.,1.,1.]

        # Empty expressions preserve components; out-of-range TimeStep emits
        # the upstream diagnostic and rewrites step 0 instead.
        three=A.view.add("three")
        A.view.add_list_data(three,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,3.,5.])
        A.plugin.set_string("ModifyComponents","Expression0","7.")
        A.plugin.set_number("ModifyComponents","View",three)
        A.plugin.set_number("ModifyComponents","TimeStep",-1)
        A.plugin.run("ModifyComponents")
        _,data=A.view.get_list_data(three)
        @test data[1][10:12]≈[7.,7.,7.]

        # Invalid TimeStep (> last) falls back to step 0 without raising.
        A.plugin.set_string("ModifyComponents","Expression0","9.")
        A.plugin.set_number("ModifyComponents","TimeStep",4)
        A.plugin.run("ModifyComponents")
        _,data=A.view.get_list_data(three)
        @test data[1][10:12]≈[9.,9.,9.]

        # A vector OtherView on a different grid interpolates like the
        # upstream octree path: v0 gets the interpolated w1 of a constant
        # vector field defined on a coarse line spanning the triangle.
        src=A.view.add("src")
        A.view.add_list_data(src,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 2.,4.,6., 2.,4.,6., 2.,4.,6.])
        dst=A.view.add("dst")
        A.view.add_list_data(dst,"SL",2,
            [0.,1., 0.,0., 0.,0., 8.,8., 9.,9.,
             0.,0., 0.,1., 0.,0., 8.,8., 9.,9.])
        A.plugin.set_string("ModifyComponents","Expression0","w1")
        A.plugin.set_number("ModifyComponents","View",dst)
        A.plugin.set_number("ModifyComponents","OtherView",src)
        A.plugin.set_number("ModifyComponents","OtherTimeStep",-1)
        A.plugin.set_number("ModifyComponents","TimeStep",-1)
        A.plugin.run("ModifyComponents")
        _,data=A.view.get_list_data(dst)
        @test data[1][7:10]≈[4.,4.,4.,4.]
        @test data[1][17:20]≈[4.,4.,4.,4.]
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

@testset "ModifyComponents on model NodeData" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.MeshSizeMin",2.0)
        A.option("Mesh.MeshSizeMax",2.0)
        A.model.add("m")
        A.model.add_box(0,0,0,1,1,1)
        A.mesh.generate(2)
        node_ids,_,_=A.mesh.get_nodes()

        # Every node written once (the upstream tagNode dedup): v0*2 doubles
        # exactly once even though surface elements share nodes.
        md=A.view.add("md")
        A.view.add_model_data(md,0,"m","NodeData",
            node_ids,[[1.5] for _ in node_ids])
        A.plugin.set_string("ModifyComponents","Expression0","v0 * 2")
        A.plugin.set_number("ModifyComponents","View",md)
        A.plugin.run("ModifyComponents")
        _,_,vals=A.view.get_model_data(md,0)
        @test Set(vals)==Set([3.0])

        # w0 from a same-grid NodeData other view at the same step.
        other=A.view.add("other")
        A.view.add_model_data(other,0,"m","NodeData",
            node_ids,[[9.0] for _ in node_ids])
        A.plugin.set_string("ModifyComponents","Expression0","w0 + v0")
        A.plugin.set_number("ModifyComponents","OtherView",other)
        A.plugin.set_number("ModifyComponents","View",md)
        A.plugin.run("ModifyComponents")
        _,_,vals=A.view.get_model_data(md,0)
        @test Set(vals)==Set([12.0])
    finally
        A.finalize()
    end
end

@testset "FieldFromAmplitudePhase complex rotation" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.MeshSizeMin",2.0)
        A.option("Mesh.MeshSizeMax",2.0)
        A.model.add("m")
        A.model.add_box(0,0,0,1,1,1)
        A.mesh.generate(2)
        node_ids,_,_=A.mesh.get_nodes()

        mktempdir() do dir
            meshfile=joinpath(dir,"fine.msh")
            A.write(meshfile)

            # Amplitude view: two NodeData steps (real=1, imag=2); phase view:
            # one scalar step phi=pi/6 — u = a * exp(j*k*phi) with k=4.
            amp=A.view.add("amp")
            A.view.add_model_data(amp,0,"m","NodeData",
                node_ids,[[1.0] for _ in node_ids])
            A.view.add_model_data(amp,1,"m","NodeData",
                node_ids,[[2.0] for _ in node_ids])
            ph=A.view.add("phase")
            phi=pi/6
            A.view.add_model_data(ph,0,"m","NodeData",
                node_ids,[[phi] for _ in node_ids])
            A.plugin.set_number("FieldFromAmplitudePhase","Wavenumber",4.0)
            A.plugin.set_number("FieldFromAmplitudePhase","AmplitudeView",amp)
            A.plugin.set_number("FieldFromAmplitudePhase","PhaseView",ph)
            A.plugin.set_string("FieldFromAmplitudePhase","MeshFile",meshfile)
            A.plugin.run("FieldFromAmplitudePhase")
            tag=last(A.view.get_tags())
            @test A.view.get_string(tag,"Name")=="FieldFromAPhi"
            dtype,tags,vals=A.view.get_model_data(tag,0)
            @test dtype==["NodeData"]
            @test sort(tags)==sort(node_ids)
            k=4.0
            @test Set(vals)==Set([cos(k*phi)-2sin(k*phi)])
            _,tags1,vals1=A.view.get_model_data(tag,1)
            @test tags1==tags
            @test Set(vals1)==Set([sin(k*phi)+2cos(k*phi)])
            # The imported-mesh model joins the list and the current model
            # stays put (upstream keeps the current-model index).
            @test "fine" in A.model.list()
            @test A.model.get_current()=="m"
        end

        # A wrong amplitude step count reports the upstream error and emits
        # no extra view; `Msg::Error` only continues when
        # `General.AbortOnError` is 0.
        A.option("General.AbortOnError",0)
        bad=A.view.add("bad")
        A.view.add_model_data(bad,0,"m","NodeData",
            node_ids,[[0.0] for _ in node_ids])
        ntags=length(A.view.get_tags())
        A.plugin.set_number("FieldFromAmplitudePhase","AmplitudeView",bad)
        A.plugin.run("FieldFromAmplitudePhase")
        @test length(A.view.get_tags())==ntags
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

@testset "ExtractEdges sharp and boundary edges" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # Two unit square panels sharing the edge x=1, z=0 (a 90° fold):
        # the shared mesh edges are sharp, the outline is boundary.  The mesh
        # arrives via a file import so the fixture needs no meshing.
        msh=joinpath(mktempdir(),"fold.msh")
        open(msh,"w") do io
            write(io,"""
\$MeshFormat
2.2 0 8
\$EndMeshFormat
\$Nodes
6
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
5 1 0 1
6 1 1 1
\$EndNodes
\$Elements
4
1 2 2 1 1 1 2 3
2 2 2 1 1 1 3 4
3 2 2 1 1 2 5 6
4 2 2 1 1 2 6 3
\$EndElements
""")
        end
        A.open(msh)

        A.plugin.set_number("ExtractEdges","IncludeBoundary",0)
        A.plugin.run("ExtractEdges")
        tag=last(A.view.get_tags())
        types,data=A.view.get_list_data(tag)
        @test types==["SL"]
        nedge=length(data[1])÷8
        @test nedge>=1
        @test all(≈(1.0),data[1][7:8:end])
        # At least one emitted edge lies along the shared fold (x≈1, z≈0).
        fold=0
        for e in 0:nedge-1
            x1,y1,z1,x2,y2,z2=data[1][8*e+1:8*e+6]
            x1≈1.0 && x2≈1.0 && z1≈0.0 && z2≈0.0 && (fold+=1)
        end
        @test fold>=1

        # Boundary edges join the output when IncludeBoundary=1.
        A.plugin.set_number("ExtractEdges","IncludeBoundary",1)
        A.plugin.run("ExtractEdges")
        tag2=last(A.view.get_tags())
        _,data2=A.view.get_list_data(tag2)
        @test length(data2[1])>length(data[1])

        # A current model with no mesh triangles → the upstream error text.
        A.model.add("empty")
        @test occursin("No triangles",_errtext(()->A.plugin.run("ExtractEdges")))
    finally
        A.finalize()
    end
end

@testset "Vorticity extension reuses Curl" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # Linear field v = (-y, x, 0) on a unit triangle: curl = (0,0,2).
        v=A.view.add("v")
        A.view.add_list_data(v,"VT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             0.,0.,0., 0.,1.,0., -1.,0.,0.])
        A.plugin.set_number("Vorticity","View",v)
        A.plugin.run("Vorticity")
        tag=last(A.view.get_tags())
        @test A.view.get_string(tag,"Name")=="v_Vorticity"
        types,data=A.view.get_list_data(tag)
        @test types==["VT"]
        @test data[1][10:18]≈[0.,0.,2., 0.,0.,2., 0.,0.,2.]

        # Scalar views contribute nothing, matching the Curl numComp gate.
        s=A.view.add("s")
        A.view.add_list_data(s,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,2.,3.])
        A.plugin.set_number("Vorticity","View",s)
        A.plugin.run("Vorticity")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test isempty(types)
    finally
        A.finalize()
    end
end

end # module PluginBatch2Checks
