# `plugin` module UI01: Eigenvalues tensor eigentriple views, ModulusPhase
# in-place complex-to-modulus/phase rewrite, HarmonicToTime harmonic
# resampling and Density per-node element-incidence counts on small
# list/model view fixtures.
module PluginFieldOpsChecks
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

@testset "Eigenvalues plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # One 9-component tensor point record holding diag(1,4,9): the
        # descending eigenvalues are (9,4,1), so the three emitted views
        # carry 1, 4 and 9.
        tp=A.view.add("tp")
        A.view.add_list_data(tp,"TP",1,
            [1., 2., 3., 1.,0.,0., 0.,4.,0., 0.,0.,9.])
        A.plugin.set_number("Eigenvalues","View",tp)
        A.plugin.run("Eigenvalues")
        tags=A.view.get_tags()
        @test length(tags)==4
        for (tag,expected,suffix) in
                ((tags[2],1.0,"_MinEigenvalues"),
                 (tags[3],4.0,"_MidEigenvalues"),
                 (tags[4],9.0,"_MaxEigenvalues"))
            types,data=A.view.get_list_data(tag)
            @test types==["SP"]
            @test data[1]≈[1.,2.,3.,expected]
            @test A.view.get_string(tag,"Name")=="tp"*suffix
        end

        # A 2-step tensor triangle: step 0 = diag(1,4,9) on every corner,
        # step 1 = [[2,1,0],[1,2,0],[0,0,5]] (eigenvalues 5,3,1).
        tt=A.view.add("tt")
        A.view.add_list_data(tt,"TT",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0.,
             1.,0.,0.,0.,4.,0.,0.,0.,9.,
             1.,0.,0.,0.,4.,0.,0.,0.,9.,
             1.,0.,0.,0.,4.,0.,0.,0.,9.,
             2.,1.,0.,1.,2.,0.,0.,0.,5.,
             2.,1.,0.,1.,2.,0.,0.,0.,5.,
             2.,1.,0.,1.,2.,0.,0.,0.,5.])
        A.plugin.set_number("Eigenvalues","View",tt)
        A.plugin.run("Eigenvalues")
        tags=A.view.get_tags()
        _,data=A.view.get_list_data(tags[end-2])
        @test data[1][10:15]≈[1.,1.,1., 1.,1.,1.]
        _,data=A.view.get_list_data(tags[end-1])
        @test data[1][10:15]≈[4.,4.,4., 3.,3.,3.]
        _,data=A.view.get_list_data(tags[end])
        @test data[1][10:15]≈[9.,9.,9., 5.,5.,5.]

        # Records without 9 components are skipped; the three views are
        # still created (empty) like upstream.
        tri=A.view.add("tri")
        A.view.add_list_data(tri,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,2.,3.])
        before=length(A.view.get_tags())
        A.plugin.set_number("Eigenvalues","View",tri)
        A.plugin.run("Eigenvalues")
        tags=A.view.get_tags()
        @test length(tags)==before+3
        for tag in tags[end-2:end]
            types,_=A.view.get_list_data(tag)
            @test isempty(types)
        end
    finally
        A.finalize()
    end
end

@testset "ModulusPhase plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # A 2-step scalar triangle: step 0 = real parts (3,-3,0), step 1 =
        # imaginary parts (4,0,-1). In place, step 0 becomes the moduli
        # (5,3,1) and step 1 the phases atan2(im,re).
        mp=A.view.add("mp")
        A.view.add_list_data(mp,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 3.,-3.,0., 4.,0.,-1.])
        A.plugin.set_number("ModulusPhase","View",mp)
        A.plugin.run("ModulusPhase")
        types,data=A.view.get_list_data(mp)
        @test types==["ST"]
        @test data[1][10:12]≈[5.,3.,1.]
        @test data[1][13:15]≈[atan(4,3),atan(0,-3),atan(-1,0)]
        @test A.view.get_string(mp,"Name")=="mp_ModulusPhase"

        # Vector records transform per component.
        vp=A.view.add("vp")
        A.view.add_list_data(vp,"VP",1,
            [5., 6., 7., 3.,0.,4., 4.,0.,-3.])
        A.plugin.set_number("ModulusPhase","View",vp)
        A.plugin.run("ModulusPhase")
        _,data=A.view.get_list_data(vp)
        @test data[1][4:6]≈[5.,0.,5.]
        @test data[1][7:9]≈[atan(4,3),0.0,atan(-3,4)]

        # Out-of-range indices emit the upstream diagnostic and leave
        # data; `Msg::Error` only continues when `General.AbortOnError` is 0.
        A.option("General.AbortOnError",0)
        bad=A.view.add("bad")
        A.view.add_list_data(bad,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 3.,-3.,0., 4.,0.,-1.])
        A.plugin.set_number("ModulusPhase","ImaginaryPart",2)
        A.plugin.set_number("ModulusPhase","View",bad)
        A.plugin.run("ModulusPhase")
        _,data=A.view.get_list_data(bad)
        @test data[1][10:15]≈[3.,-3.,0.,4.,0.,-1.]
        A.plugin.set_number("ModulusPhase","ImaginaryPart",1)

        # Model NodeData: stored rows are per tag, so two steps' rows are
        # rewritten in place — (3,4) -> modulus 5, phase atan(4,3).
        A.model.add("box")
        md=A.view.add("md")
        A.view.add_model_data(md,0,"box","NodeData",
            UInt64[1,2,3],[[3.],[0.],[6.]])
        A.view.add_model_data(md,1,"box","NodeData",
            UInt64[1,2,3],[[4.],[2.],[-8.]])
        A.plugin.set_number("ModulusPhase","View",md)
        A.plugin.run("ModulusPhase")
        dtype,tags0,vals0=A.view.get_model_data(md,0)
        dtype1,tags1,vals1=A.view.get_model_data(md,1)
        @test dtype==["NodeData"] && dtype1==["NodeData"]
        @test sort(tags0)==[1,2,3] && sort(tags1)==[1,2,3]
        @test sort(vals0)≈[2.,5.,10.]
        @test sort(vals1)≈sort([atan(4,3),atan(2,0),atan(-8,6)])
        @test A.view.get_string(md,"Name")=="md_ModulusPhase"
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

@testset "HarmonicToTime plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # re=1, im=0 at every corner; NumSteps=4 samples cos(2 pi k/4) =
        # 1, 0, -1, 0 -> a 4-step scalar view named *_HarmonicToTime.
        h=A.view.add("h")
        A.view.add_list_data(h,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 1.,1.,1., 0.,0.,0.])
        A.plugin.set_number("HarmonicToTime","NumSteps",4)
        A.plugin.set_number("HarmonicToTime","View",h)
        A.plugin.run("HarmonicToTime")
        tag=last(A.view.get_tags())
        types,data=A.view.get_list_data(tag)
        @test types==["ST"]
        @test data[1][10:21]≈[1.,1.,1., 0.,0.,0., -1.,-1.,-1., 0.,0.,0.]
        @test A.view.get_number(tag,"NbTimeStep")==4
        @test A.view.get_string(tag,"Name")=="h_HarmonicToTime"

        # re=0, im=1 with the default TimeSign=-1 gives -sin(p): 0,-1,0,1.
        h2=A.view.add("h2")
        A.view.add_list_data(h2,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,0.,0., 1.,1.,1.])
        A.plugin.set_number("HarmonicToTime","View",h2)
        A.plugin.run("HarmonicToTime")
        _,data=A.view.get_list_data(last(A.view.get_tags()))
        @test data[1][10:21]≈[0.,0.,0., -1.,-1.,-1., 0.,0.,0., 1.,1.,1.]

        # TimeSign=+1 flips the sine term: 0,1,0,-1.
        h3=A.view.add("h3")
        A.view.add_list_data(h3,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,0.,0., 1.,1.,1.])
        A.plugin.set_number("HarmonicToTime","TimeSign",1)
        A.plugin.set_number("HarmonicToTime","View",h3)
        A.plugin.run("HarmonicToTime")
        _,data=A.view.get_list_data(last(A.view.get_tags()))
        @test data[1][10:21]≈[0.,0.,0., 1.,1.,1., 0.,0.,0., -1.,-1.,-1.]
        A.plugin.set_number("HarmonicToTime","TimeSign",-1)

        # Frequency=0 collapses the phase to p=0: every sample is the real
        # part and all sample times are 0.
        h4=A.view.add("h4")
        A.view.add_list_data(h4,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 2.,2.,2., 1.,1.,1.])
        A.plugin.set_number("HarmonicToTime","Frequency",0)
        A.plugin.set_number("HarmonicToTime","NumSteps",3)
        A.plugin.set_number("HarmonicToTime","View",h4)
        A.plugin.run("HarmonicToTime")
        _,data=A.view.get_list_data(last(A.view.get_tags()))
        @test data[1][10:18]≈[2.,2.,2., 2.,2.,2., 2.,2.,2.]
        A.plugin.set_number("HarmonicToTime","Frequency",1)
        A.plugin.set_number("HarmonicToTime","NumSteps",20)

        # NumSteps <= 0 and out-of-range step indices log the upstream
        # errors and create nothing; `Msg::Error` only continues when
        # `General.AbortOnError` is 0.
        A.option("General.AbortOnError",0)
        before=length(A.view.get_tags())
        A.plugin.set_number("HarmonicToTime","NumSteps",0)
        A.plugin.set_number("HarmonicToTime","View",h4)
        A.plugin.run("HarmonicToTime")
        @test length(A.view.get_tags())==before
        A.plugin.set_number("HarmonicToTime","NumSteps",4)
        A.plugin.set_number("HarmonicToTime","ImaginaryPart",5)
        A.plugin.run("HarmonicToTime")
        @test length(A.view.get_tags())==before
        A.plugin.set_number("HarmonicToTime","ImaginaryPart",1)
        A.plugin.set_number("HarmonicToTime","NumSteps",20)
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

@testset "Density plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # Two triangles sharing the (1,0)-(0,1) edge: shared corner
        # positions count 2 incident elements, the free corners count 1.
        # Element blocks are [x1..x3, y1..y3, z1..z3, values].
        d=A.view.add("d")
        A.view.add_list_data(d,"ST",2,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 9.,9.,9.,
             1.,1.,0., 0.,1.,1., 0.,0.,0., 8.,8.,8.])
        A.plugin.set_number("Density","View",d)
        A.plugin.run("Density")
        tag=last(A.view.get_tags())
        types,data=A.view.get_list_data(tag)
        @test types==["ST"]
        # Single-step scalar records: coords then one value per node.
        @test data[1][1:9]≈[0.,1.,0., 0.,0.,1., 0.,0.,0.]
        @test data[1][10:12]≈[1.,2.,2.]
        @test data[1][13:21]≈[1.,1.,0., 0.,1.,1., 0.,0.,0.]
        @test data[1][22:24]≈[2.,1.,2.]
        @test A.view.get_number(tag,"NbTimeStep")==1
        @test A.view.get_string(tag,"Name")=="d_Density"

        # Lines count their endpoints: [0,0,0]-[1,0,0] and [1,0,0]-[2,0,0]
        # share the midpoint at incidence 2.
        l=A.view.add("l")
        A.view.add_list_data(l,"SL",2,
            [0.,1., 0.,0., 0.,0., 1.,2.,
             1.,2., 0.,0., 0.,0., 3.,4.])
        A.plugin.set_number("Density","View",l)
        A.plugin.run("Density")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SL"]
        @test data[1][7:8]≈[1.,2.]
        @test data[1][15:16]≈[2.,1.]
    finally
        A.finalize()
    end
end

end # module PluginFieldOpsChecks
