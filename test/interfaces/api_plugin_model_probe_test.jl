module PluginModelProbeTests
using Test,Tessella
const A=Tessella.API

function fixture(dir,name,msh,coords)
    path=joinpath(dir,name*".msh")
    ids=UInt64[10i+1 for i in eachindex(coords)]
    open(path,"w") do io
        println(io,raw"$MeshFormat");println(io,"2.2 0 8");println(io,raw"$EndMeshFormat")
        println(io,raw"$Nodes");println(io,length(coords))
        for (uid,p) in zip(ids,coords);println(io,uid," ",join(p," "));end
        println(io,raw"$EndNodes");println(io,raw"$Elements");println(io,1)
        println(io,"77 ",msh," 2 0 1 ",join(ids," "));println(io,raw"$EndElements")
    end
    A.open(path)
    return A.model.get_current(),ids
end
function probe(view,point)
    A.plugin.set_number("Probe","View",view)
    for (name,value) in zip(("X","Y","Z"),point)
        A.plugin.set_number("Probe",name,value)
    end
    A.plugin.run("Probe")
    return A.view.get_list_data(last(A.view.get_tags()))
end
function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir(f)
    finally
        A.finalize()
    end
end

@testset "Model Probe uses actual mixed-element shape functions and UIDs" begin
    fixtures=(
        ("line",1,[(0.,0.,0.),(1.,0.,0.)],(.25,0.,0.),p->1+2p[1]),
        ("quad",3,[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)],
            (.25,.5,0.),p->1+2p[1]+3p[2]+4p[1]*p[2]),
        ("hex",5,[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
                  (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.)],
            (.25,.5,.25),p->1+2p[1]+3p[2]+4p[3]+5p[1]*p[2]+6p[1]*p[3]+7p[2]*p[3]+8prod(p)),
        ("prism",6,[(0.,0.,0.),(1.,0.,0.),(0.,1.,0.),
                    (0.,0.,1.),(1.,0.,1.),(0.,1.,1.)],
            (.2,.3,.4),p->1+2p[1]+3p[2]+4p[3]+5p[1]*p[3]+6p[2]*p[3]),
        ("pyramid",7,[(-1.,-1.,0.),(1.,-1.,0.),(1.,1.,0.),(-1.,1.,0.),(0.,0.,1.)],
            (.1,.2,.25),p->1+2p[1]+3p[2]+4p[3]),
        # Curved Line3: x=(1+u)/2, y=(1-u^2)/2. The query is on the
        # parabola, away from its straight corner hull, with f=1+x^2.
        ("curved_line3",8,[(0.,0.,0.),(1.,0.,0.),(.5,.5,0.)],
            (.25,.375,0.),p->1+p[1]^2))
    for (name,msh,coords,point,value) in fixtures
        session() do dir
            model,ids=fixture(dir,name,msh,coords)
            rows=[[value(p)] for p in coords]
            for kind in ("NodeData","ElementNodeData")
                view=A.view.add(name*kind)
                A.view.add_model_data(view,0,model,kind,
                    kind=="NodeData" ? ids : [77],
                    kind=="NodeData" ? rows : [reduce(vcat,rows)])
                types,data=probe(view,point)
                @test types==["SP"]
                @test only(data)≈[point...,value(point)]
            end
            # A dense index does not identify a sparse-UID element.
            unknown=A.view.add("unknown UID")
            A.view.add_model_data(unknown,0,model,"ElementData",[1],[[42.]])
            types,data=probe(unknown,point)
            @test types==["SP"] && only(data)≈[point...,0.]
        end
    end
end

@testset "Model Probe follows each step's owned model without current-model dependence" begin
    session() do dir
        model,ids=fixture(dir,"first",3,[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)])
        view=A.view.add("two model steps")
        A.view.add_model_data(view,0,model,"NodeData",ids,[[1.],[3.],[10.],[4.]],2.0)
        A.model.add("second")
        model2,ids2=fixture(dir,"second",3,[(0.,0.,0.),(2.,0.,0.),(2.,2.,0.),(0.,2.,0.)])
        A.view.add_model_data(view,1,model2,"NodeData",ids2,[[10.],[14.],[20.],[16.]],3.0)
        # Same UIDs in different models denote different geometry and rows.
        types,data=probe(view,(.25,.5,0.))
        @test types==["SP"]
        @test only(data)≈[.25,.5,0.,3.5,14.0]
        alias=A.view.add_alias(view)
        A.view.remove(view)
        @test probe(alias,(.25,.5,0.))==(types,data)
        types,data=probe(alias,(4.,4.,0.))
        @test isempty(types) && isempty(data)
    end
end
@testset "Model probe snapshots are reused and released with their run" begin
    session() do dir
        model,ids=fixture(dir,"owned",3,[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)])
        view=A.view.add("owned field")
        A.view.add_model_data(view,0,model,"NodeData",ids,[[1.],[3.],[10.],[4.]])
        owner=A.SESSION_VIEWS[view]
        contexts=IdDict{Any,Any}()
        key=A._PLUGIN_MODEL_PROBE_CONTEXT_KEY
        lock(A.STATE_LOCK) do
            task_local_storage(key,contexts) do
                first_values,_=A._plugin_probe_model(owner,(.25,.5,0.),"probe cache check")
                snapshot=only(values(contexts))
                second_values,_=A._plugin_probe_model(owner,(.5,.25,0.),"probe cache check")
                @test only(first_values)≈[3.5]
                @test only(second_values)≈[3.25]
                @test length(contexts)==1 && only(values(contexts))===snapshot
                slot=owner.model_contexts[(0,"NodeData")].model_slot
                @test snapshot.mesh.coords!==slot.mesh.coords
                @test snapshot.public.element_tags==[77]
                # A nested public plugin run owns a distinct cache and restores
                # the caller's context, including the failure cleanup path.
                @test only(last(probe(view,(.25,.5,0.))))≈[.25,.5,0.,3.5]
                @test current_task().storage[key]===contexts
                bad=A.view.add("unsupported components")
                A.view.add_model_data(bad,0,model,"NodeData",ids,fill([1.,2.],4))
                @test_throws ArgumentError probe(bad,(.25,.5,0.))
                @test current_task().storage[key]===contexts
            end
        end
        storage=current_task().storage
        @test storage===nothing || !haskey(storage,key)
    end
end
@testset "Model interpolation honors dimension and inactive UID authorities" begin
    session() do dir
        path=joinpath(dir,"mixed.msh")
        write(path,raw"""$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
11 0 0 0
21 1 0 0
31 0 1 0
41 0 0 1
$EndNodes
$Elements
2
66 2 2 0 1 11 21 31
77 4 2 0 2 11 21 31 41
$EndElements
""")
        A.open(path);model=A.model.get_current()
        other=A.view.add("dimension values")
        A.view.add_model_data(other,0,model,"ElementData",[66,77],[[2.],[10.]])
        source=A.view.add("surface sampling points")
        A.view.add_list_data(source,"ST",1,[.2,.4,.2, .2,.2,.4, 0.,0.,0., 1.,1.,1.])
        A.plugin.set_number("MathEval","View",source)
        A.plugin.set_number("MathEval","OtherView",other)
        A.plugin.set_number("MathEval","Dimension",2)
        A.plugin.set_string("MathEval","Expression0","w0")
        A.plugin.run("MathEval")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["ST"]
        @test only(data)[10:12]≈[2.,2.,2.]

        unknown=A.view.add("unknown element")
        A.view.add_model_data(unknown,0,model,"ElementData",[1],[[42.]])
        A.plugin.set_number("Integrate","View",unknown)
        A.plugin.run("Integrate")
        _,data=A.view.get_list_data(last(A.view.get_tags()))
        @test last(only(data))==0.0

        nodal=A.view.add("inactive node UIDs")
        A.view.add_model_data(nodal,0,model,"NodeData",[11,21,31,41],fill([1.],4))
        A.model.add("unrelated current model")
        A.plugin.set_number("Integrate","Dimension",2)
        for (view,integral) in ((other,1.0),(nodal,0.5))
            A.plugin.set_number("Integrate","View",view)
            A.plugin.run("Integrate")
            _,data=A.view.get_list_data(last(A.view.get_tags()))
            @test last(only(data))≈integral
        end
    end
end
end
