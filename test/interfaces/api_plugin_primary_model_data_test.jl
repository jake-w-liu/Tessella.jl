module PluginPrimaryModelDataTests
using Test,Tessella
const A=Tessella.API

const POINTS=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(.5,0.,.2),(.5,.5,.1),(0.,.5,.1))
const IDS=UInt64.(11:16)
const VALUES=[0.,1.,1.,4.,4.,4.]
const CORNERS=[0.,1.,0.,0.,0.,1.,0.,0.,0.]

function fixture(dir,materialized)
    if materialized
        path=joinpath(dir,"triangle.msh")
        open(path,"w") do io
            println(io,raw"$MeshFormat");println(io,"2.2 0 8");println(io,raw"$EndMeshFormat")
            println(io,raw"$Nodes");println(io,length(POINTS))
            for (id,p) in zip(IDS,POINTS);println(io,id," ",join(p," "));end
            println(io,raw"$EndNodes");println(io,raw"$Elements");println(io,1)
            println(io,"301 9 2 0 7 ",join(IDS," "));println(io,raw"$EndElements")
        end
        A.open(path)
    else
        A.model.add("builder")
        A.model.add_discrete_entity(2,7)
        A.mesh.add_nodes(2,7,IDS,collect(Iterators.flatten(POINTS)))
        A.mesh.add_elements_by_type(7,9,[301],IDS)
    end
    return A.model.get_current()
end

function run_plugin(name,source)
    A.plugin.set_number(name,"View",source)
    before=A.view.get_tags()
    A.plugin.run(name)
    return [A.view.get_list_data(tag) for tag in setdiff(A.view.get_tags(),before)]
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

@testset "Higher-order model primitive plugins retain primary data and full bounds" begin
    # SDK 4.15.2 controls: Tri6 midside values do not enter these primitive
    # operators, but curved midside coordinates still determine the SP anchor.
    for materialized in (false,true),kind in ("NodeData","ElementNodeData")
        @testset "cached=$materialized $kind" begin
            session() do dir
                model=fixture(dir,materialized)
                source=A.view.add("curved source")
                tags=kind=="NodeData" ? IDS : UInt64[301]
                rows=kind=="NodeData" ? [[v] for v in VALUES] : [VALUES]
                A.view.add_model_data(source,0,model,kind,tags,rows)
                before=(A.view.get_model_data(source,0),A.mesh.get_nodes(),A.mesh.get_elements())
                @test only(A.view.probe(source,.2,.2,.16)[1])≈4.24
                A.plugin.set_string("MathEval","Expression0","v0")
                types,data=only(run_plugin("MathEval",source))
                @test types==["ST"]
                @test only(data)==[CORNERS;VALUES[1:3]]
                types,data=only(run_plugin("Gradient",source))
                @test types==["VT"]
                @test only(data)≈[CORNERS;repeat([1.,1.,0.],3)]
                extrema=run_plugin("MinMax",source)
                @test length(extrema)==2
                @test all(first(record)==["SP"] for record in extrema)
                @test sort([only(last(record))[end] for record in extrema])==[0.,1.]
                @test all(only(last(record))[1:3]≈[.5,.5,.1] for record in extrema)
                types,data=only(run_plugin("Integrate",source))
                @test types==["SP"]
                @test only(data)≈[.5,.5,.1,1/3]
                @test (A.view.get_model_data(source,0),A.mesh.get_nodes(),A.mesh.get_elements())==before
                @test only(A.view.probe(source,.2,.2,.16)[1])≈4.24
            end
        end
    end
end

@testset "Primitive projection preserves component layout and populated steps" begin
    for kind in ("NodeData","ElementNodeData"),nc in (3,9)
        session() do dir
            model=fixture(dir,false)
            source=A.view.add("component source")
            rows=[[v+10c for c in 1:nc] for v in VALUES]
            tags=kind=="NodeData" ? IDS : UInt64[301]
            payload=kind=="NodeData" ? rows : [reduce(vcat,rows)]
            A.view.add_model_data(source,0,model,kind,tags,payload,2.,nc)
            A.view.add_model_data(source,2,model,kind,tags,
                kind=="NodeData" ? [2r for r in rows] : [2only(payload)],5.,nc)
            for c in 0:8
                A.plugin.set_string("MathEval","Expression$c",c<nc ? "v$c" : "")
            end
            types,data=only(run_plugin("MathEval",source))
            @test types==[nc==3 ? "VT" : "TT"]
            expected=reduce(vcat,rows[1:3])
            @test only(data)≈[CORNERS;expected;2expected]
            @test A.view.get_model_data(source,2)[3]==reduce(vcat,
                kind=="NodeData" ? [2r for r in rows] : [2only(payload)])
        end
    end
end

@testset "Builder model snapshots are reused by primitive and interpolation paths" begin
    session() do dir
        model=fixture(dir,false)
        source=A.view.add("owned builder source")
        A.view.add_model_data(source,0,model,"NodeData",IDS,[[v] for v in VALUES])
        owner=A.SESSION_VIEWS[source]
        slot=owner.model_contexts[(0,"NodeData")].model_slot
        contexts=IdDict{Any,Any}()
        key=A._PLUGIN_MODEL_PROBE_CONTEXT_KEY
        lock(A.STATE_LOCK) do
            task_local_storage(key,contexts) do
                first_elements,_=A._plugin_model_elements(owner,(0,"NodeData"),"snapshot check")
                snapshot=only(values(contexts))
                second_elements,_=A._plugin_model_elements(owner,(0,"NodeData"),"snapshot check")
                @test first_elements[1].vals==second_elements[1].vals==VALUES[1:3]
                @test only(values(contexts))===snapshot
                probe_context=A._plugin_model_probe_context(slot,contexts,"snapshot check")
                @test probe_context.mesh===snapshot.mesh
                @test probe_context.public===snapshot.public
                @test A._plugin_model_probe_context(slot,contexts,"snapshot check")===probe_context
                @test A._plugin_model_mesh(slot,"snapshot check")===snapshot.mesh
                @test snapshot.public.element_tags==[301]
                @test snapshot.public.node_tags==IDS
            end
        end
        storage=current_task().storage
        @test storage===nothing || !haskey(storage,key)
    end
end
end
