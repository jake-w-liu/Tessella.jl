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

        # setColor/getColor pack upstream's `A<<24|B<<16|G<<8|R` word into the
        # same numeric option table, so copyOptions carries colors too.
        A.view.set_color(src,"Color.Ambient",10,20,30,40)
        @test A.view.get_color(src,"Color.Ambient")==(10,20,30,40)
        A.view.set_color(src,"Color.Points",1,2,3)   # default alpha 255
        @test A.view.get_color(src,"Color.Points")==(1,2,3,255)
        A.view.copy_options(src,dst)
        @test A.view.get_color(dst,"Color.Ambient")==(10,20,30,40)
        @test A.view.get_color(dst,"Color.Missing")==(0,0,0,0)
        @test_throws ArgumentError A.view.set_color(src,"C",-1,0,0)
        @test_throws ArgumentError A.view.set_color(src,"C",0,0,256)
        @test occursin("unknown view tag",
            _errtext(()->A.view.set_color(99999,"C",0,0,0)))

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

@testset "view.probe interpolates nodal data inside elements" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.MeshSizeMax",0.5)
        sq=A.model.occ.add_rectangle(0,0,0,1,1)
        A.model.occ.synchronize()
        A.mesh.generate(2)
        mname=A.model.get_current()
        ntags,ncoords,_=A.mesh.get_nodes(-1,-1,true,false)
        coord_of=Dict(ntags[i]=>ncoords[3i-2:3i] for i in eachindex(ntags))

        # ElementNodeData: f = x+2y at every corner — linear field is
        # interpolated exactly by the barycentric basis.
        etypes,etags,enodes=A.mesh.get_elements(2,-1)
        @test etypes==Int32[2]
        npe=3
        rows=Vector{Float64}[]
        flat_etags=reduce(vcat,etags)
        for (btags,bnodes) in zip(etags,enodes)
            for e in eachindex(btags)
                push!(rows,[let xy=coord_of[bnodes[(e-1)*npe+j]]
                                xy[1]+2.0*xy[2]
                            end for j in 1:npe])
            end
        end
        v=A.view.add("field")
        A.view.add_model_data(v,0,mname,"ElementNodeData",
                              flat_etags,rows,0.0,1,0)
        for (x,y) in ((0.3,0.4),(0.5,0.5),(0.1,0.8))
            got=A.view.probe(mname,v,x,y,0.0)
            @test length(got)==1
            @test got[1]≈x+2.0*y atol=1e-10
        end
        # NodeData path and ElementData path.
        v2=A.view.add("nodal")
        nvals=[coord_of[t][1]+2.0*coord_of[t][2] for t in ntags]
        A.view.add_model_data(v2,0,mname,"NodeData",ntags,nvals)
        @test A.view.probe(mname,v2,0.3,0.4,0.0)[1]≈0.3+0.8 atol=1e-10
        v3=A.view.add("per-element")
        A.view.add_model_data(v3,0,mname,"ElementData",
                              flat_etags,fill(7.5,length(flat_etags)))
        @test A.view.probe(mname,v3,0.3,0.4,0.0)==[7.5]
        # Outside the mesh and unknown model reject.
        @test occursin("no element contains",
            _errtext(()->A.view.probe(mname,v,5.0,5.0,0.0)))
        @test occursin("Unknown model",
            _errtext(()->A.view.probe("nope",v,0.3,0.4,0.0)))
        @test occursin("unknown view tag",
            _errtext(()->A.view.probe(mname,99999,0.3,0.4,0.0)))
    finally
        A.finalize()
    end
end
end
