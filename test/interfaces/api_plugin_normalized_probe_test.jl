module PluginNormalizedProbeChecks
using Test,Tessella
const A=Tessella.API
const E=Tessella.Elements

function geometry(kind)
    if kind===:quadrangle
        return [0. 1. 1. 0.;0. 0. 1. 1.;0. 0. 0. 0.],"SQ"
    end
    msh=kind===:hexahedron ? 5 : kind===:prism ? 6 : 7
    coords=E.lagrange_nodes(msh)
    if kind!==:prism;coords[1:2,:].=(coords[1:2,:].+1)./2;end
    kind!==:pyramid && (coords[3,:].=(coords[3,:].+1)./2)
    return coords,kind===:hexahedron ? "SH" : kind===:prism ? "SI" : "SY"
end

function source_view(kind,scale,offset)
    coords,record=geometry(kind)
    values=coords[1,:].+2coords[2,:].+3coords[3,:]
    coords.*=scale;coords[1:2,:].+=offset
    payload=vcat(vec(coords[1,:]),vec(coords[2,:]),vec(coords[3,:]),values)
    view=A.view.add("source")
    A.view.add_list_data(view,record,1,payload)
    return view,coords,payload
end

function cutgrid(view,point)
    for (key,value) in (("View",view),("NumPointsU",1),("NumPointsV",1),
                        ("ConnectPoints",0),("X0",point[1]),("Y0",point[2]),("Z0",point[3]))
        A.plugin.set_number("CutGrid",key,value)
    end
    before=A.view.get_tags();A.plugin.run("CutGrid")
    output=only(setdiff(A.view.get_tags(),before))
    return A.view.get_list_data(output)
end

@testset "Plugin resampling uses normalized certified element maps" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        for kind in (:quadrangle,:hexahedron,:prism,:pyramid),
                (scale,offset) in ((1.,0.),(1.,1e12),(1e-150,0.),(1e150,0.))
            source,coords,payload=source_view(kind,scale,offset)
            saved=copy(payload)
            z=kind===:quadrangle ? 0. : .25scale
            point=(offset+.25scale,offset+.25scale,z)
            types,data=cutgrid(source,point)
            @test types==["SP"]
            @test last(only(data))≈(kind===:quadrangle ? .75 : 1.5) atol=2e-13
            @test only(data)[1:3]==collect(point)
            @test payload==saved
            _,outside=cutgrid(source,(offset+3scale,offset+.25scale,z))
            @test last(only(outside))==0.
        end
        warped=[0. 1. 1. 0.;0. 0. 1. 1.;0. 0. .2 0.]
        view=A.view.add("warped quad")
        values=warped[1,:].+2warped[2,:].+3warped[3,:]
        A.view.add_list_data(view,"SQ",1,vcat(vec(warped[1,:]),vec(warped[2,:]),vec(warped[3,:]),values))
        _,data=cutgrid(view,(.25,.25,.0125))
        @test last(only(data))≈.7875 atol=2e-13

        # MathEval's OtherView spatial resampling reaches the same inverse.
        other,_,_=source_view(:quadrangle,1.,1e12)
        target=A.view.add("sampling point")
        A.view.add_list_data(target,"SP",1,[1e12+.25,1e12+.25,0.,0.])
        A.plugin.set_number("MathEval","View",target)
        A.plugin.set_number("MathEval","OtherView",other)
        A.plugin.set_number("MathEval","ForceInterpolation",1)
        A.plugin.set_string("MathEval","Expression0","w0")
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SP"]
        @test last(only(data))≈.75 atol=2e-13

        # RK4 integration of dx/dt=1+x has a known one-step polynomial.
        offset=1e9;coords,_=geometry(:quadrangle);local_x=copy(coords[1,:])
        coords[1:2,:].+=offset
        vectors=zeros(3,4);vectors[1,:].=1 .+local_x
        vector_view=A.view.add("translated linear flow")
        A.view.add_list_data(vector_view,"VQ",1,vcat(vec(coords[1,:]),vec(coords[2,:]),vec(coords[3,:]),vec(vectors)))
        for (key,value) in (("View",vector_view),("NumPointsU",1),("NumPointsV",1),
                            ("X0",offset+.25),("Y0",offset+.25),("Z0",0.),
                            ("DT",.1),("MaxIter",3))
            A.plugin.set_number("StreamLines",key,value)
        end
        A.plugin.run("StreamLines")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["VP"]
        factor=1+.1+.1^2/2+.1^3/6+.1^4/24
        @test only(data)[4:3:10]≈[1.25*(factor^step-1) for step in 1:3] atol=2e-6
    finally A.finalize() end
end
end
