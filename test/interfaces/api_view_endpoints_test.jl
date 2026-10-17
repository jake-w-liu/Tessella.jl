# `model.view` VIEW02 endpoints — per-element string list data,
# option copying, per-window visibility and homogeneous model data.
module ViewEndpointsChecks
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

@testset "view string lists, options, visibility, homogeneous data" begin
    A.initialize()
    try
        A.option("General.Terminal",0)

        view=A.view.add("labels")
        labels=["elem A","elem B","elem C"]
        anchors=[0.0,0.0,0.0, 0.5,0.5,0.0, 1.0,1.0,0.0]
        A.view.add_list_data_string(view,"ST",labels,anchors)
        types,data=A.view.get_list_data_strings(view,3)
        @test types==["ST"]
        @test data==labels
        types2,data2,coords=A.view.get_list_data_strings(
            view,3,-1,true)
        @test types2==["ST"] && data2==labels && coords==anchors

        # Second block with a different element count selects by numEle.
        A.view.add_list_data_string(view,"SQ",["q1"],
                                    [0.0,0.0,0.0])
        t1,d1=A.view.get_list_data_strings(view,1)
        @test t1==["SQ"] && d1==["q1"]

        # Validation.
        @test occursin("only strings",
            _errtext(()->A.view.add_list_data_string(
                view,"ST",[1,2],[0.0,0.0,0.0])))
        @test occursin("3 coordinates",
            _errtext(()->A.view.add_list_data_string(
                view,"ST",["x"],[0.0,0.0])))
        @test occursin("numEle",
            _errtext(()->A.view.add_list_data_string(
                view,"ST",["x","y"],zeros(6),1)))
        @test occursin("no string list data",
            _errtext(()->A.view.get_list_data_strings(view,7)))
        @test occursin("unknown view tag",
            _errtext(()->A.view.get_list_data_strings(99999,1)))

        # copyOptions moves numbers and strings (name stays per-view).
        src=A.view.add("source")
        A.view.set_number(src,"Scale",2.5)
        A.view.set_string(src,"Format","binary")
        dst=A.view.add("target")
        A.view.copy_options(src,dst)
        @test A.view.get_number(dst,"Scale")==2.5
        @test A.view.get_string(dst,"Format")=="binary"
        @test A.view.get_string(dst,"Name")=="target"

        # Per-window visibility persists as a numbered option row.
        A.view.set_visibility_per_window(dst,0.0,2)
        @test A.view.get_number(dst,"VisibilityPerWindow:2")==0.0
        @test A.view.get_number(dst,"VisibilityPerWindow:0")==0.0

        # Homogeneous model data round-trips through model storage.
        box=A.model.occ.add_box(0,0,0,1,1,1)
        A.model.occ.synchronize()
        mname=A.model.get_current()
        md=A.view.add("model scalars")
        A.view.add_homogeneous_model_data(
            md,0,mname,"ElementData",[1,2],[1.5,2.5])
        ht,htags,hdata=A.view.get_homogeneous_model_data(md,0,mname)
        @test ht==["ElementData"]
        @test htags==[1,2]
        @test hdata==[1.5,2.5]
        @test occursin("Unknown model",
            _errtext(()->A.view.get_homogeneous_model_data(
                md,0,"no such model")))
    finally
        A.finalize()
    end
end
end
