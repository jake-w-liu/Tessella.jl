module APIGeoViewTransportChecks
using Test,Tessella
const A=Tessella.API
const P=Tessella.Post

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir(f)
    finally
        A.finalize()
    end
end

@testset "GEO imports publish owned list schemes, time and shared aliases" begin
    session() do dir
        pos=joinpath(dir,"views.pos")
        write(pos,"""
View "time" { TIME{2,5};SL(0,0,0,1,0,0){.25,.25,.5,.5}; };
View "quadratic" { SL2(0,0,0,1,0,0,.5,0,0){1,2,1.25}; };
""")
        geo=joinpath(dir,"views.geo")
        write(geo,"Merge \"views.pos\";Alias View[0];")
        A.open(geo)
        tags=A.view.get_tags()
        @test length(tags)==3 && allunique(tags)
        @test [A.view.get_string(tag,"Name") for tag in tags]==["time","quadratic","time"]
        @test first(A.view.probe(tags[1],.25,0.,0.))==[.25,.5]
        @test first(A.view.probe(tags[2],.25,0.,0.))≈[1.0625]
        exported=joinpath(dir,"time.pos");A.view.write(tags[1],exported)
        @test only(P.read_pos(exported)).time==[2.,5.]
        A.view.add_list_data(tags[1],"SL",1,[0.,1.,0.,0.,0.,0.,.6,.6,.7,.7])
        @test A.view.get_list_data(tags[3])==A.view.get_list_data(tags[1])
        A.view.remove(tags[1])
        @test first(A.view.probe(tags[3],.25,0.,0.))==[.6,.7]
        @test read(pos,String)=="""
View "time" { TIME{2,5};SL(0,0,0,1,0,0){.25,.25,.5,.5}; };
View "quadratic" { SL2(0,0,0,1,0,0,.5,0,0){1,2,1.25}; };
"""
    end
end

@testset "GEO view selection uses actual session order, tags and model sources" begin
    session() do dir
        A.model.add("source")
        p=A.model.add_point(0.,0.,0.;meshSize=.5)
        q=A.model.add_point(1.,0.,0.;meshSize=.5)
        A.model.add_line(p,q);A.mesh.generate(1)
        ids,_,_=A.mesh.get_nodes()
        held=A.view.add("held model",100)
        A.view.add_model_data(held,0,"source","NodeData",ids,fill([.5],length(ids)))
        alias=A.view.add_alias(held,false,200)
        held_object=A.SESSION_VIEWS[held]
        source_slot=only(values(held_object.model_contexts)).model_slot
        pos=joinpath(dir,"incoming.pos")
        write(pos,"View \"incoming\" { SL(0,0,0,1,0,0){.25,.25}; };")
        geometry="Point(11)={0,0,0,1};Point(21)={1,0,0,1};Line(31)={11,21};"
        options="Mesh.MeshSizeFromPoints=0;Mesh.MeshSizeExtendFromBoundary=0;"
        geo=joinpath(dir,"target.geo")
        write(geo,"Merge \"incoming.pos\";"*geometry*options*
            "Field[1]=PostView;Field[1].ViewIndex=0;Background Field=1;")
        A.open(geo)
        @test A.view.get_tags()==[100,200,201]
        @test A.SESSION_VIEWS[held]===held_object
        @test only(values(held_object.model_contexts)).model_slot===source_slot
        A.mesh.generate(1)
        @test length(first(A.mesh.get_nodes()))==3
        @test A.view.probe("source",held,.25,0.,0.)==[.5]
        # The second input sees all existing views. Literal ViewTag=202
        # selects its newly allocated source, independent of list position.
        write(geo,"Merge \"incoming.pos\";"*geometry*options*
            "Field[1]=PostView;Field[1].ViewTag=202;Background Field=1;")
        A.open(geo);A.mesh.generate(1)
        @test A.view.get_tags()==[100,200,201,202]
        @test length(first(A.mesh.get_nodes()))==5
        @test A.SESSION_VIEWS[alias].alias_of===held_object
        A.view.add_model_data(held,0,"source","NodeData",ids,fill([.75],length(ids)))
        @test first(A.view.probe(alias,.25,0.,0.))==[.75]
        # A fatally failed script has already staged a deletion and a new
        # view. Its failure must leave held identities, data and allocation
        # alone. `Error(...)` is a recoverable diagnostic, so the rollback
        # path is reached through a nonzero `Exit`.
        tags=A.view.get_tags();data=deepcopy(A.view.get_model_data(held,0))
        current=A.model.get_current();entities=A.model.get_entities()
        failure=joinpath(dir,"failure.geo")
        write(failure,"Delete View[0];Merge \"incoming.pos\";Exit 2;")
        @test_throws ArgumentError A.merge(failure)
        @test A.view.get_tags()==tags
        @test A.view.get_model_data(held,0)==data
        @test A.model.get_current()==current && A.model.get_entities()==entities
        @test A.SESSION_VIEWS[held]===held_object
        @test A.view.add("after failure")==203
    end
end

@testset "GEO Delete Empty Views tests owned session records" begin
    session() do dir
        empty=A.view.add("session empty")
        used=A.view.add("session used")
        A.view.add_list_data(used,"SL",1,[0.,0.,0.,1.,0.,0.,.4,.6])
        geo=joinpath(dir,"empty.geo")
        write(geo,"Delete Empty Views;")
        A.merge(geo)
        @test A.view.get_tags()==[used]
        @test A.view.get_list_data(used)==(["SL"],
            [[0.,0.,0.,1.,0.,0.,.4,.6]])
        @test empty!=used
    end
end

@testset "GEO plugin output retains held and newly created list aliases" begin
    session() do dir
        source=A.view.add("source",100)
        A.view.add_list_data(source,"SP",1,[0.,0.,0.,1.])
        alias=A.view.add_alias(source,false,200)
        owner=A.SESSION_VIEWS[source];alias_object=A.SESSION_VIEWS[alias]
        geo=joinpath(dir,"plugin.geo")
        write(geo,"Plugin(Probe).View=100;Plugin(Probe).X=0;" *
            "Plugin(Probe).Y=0;Plugin(Probe).Z=0;Plugin(Probe).Run;Alias View[0];")
        A.merge(geo)
        @test A.view.get_tags()==[100,200,201,202]
        @test A.SESSION_VIEWS[source]===owner && A.SESSION_VIEWS[alias]===alias_object
        @test A.SESSION_VIEWS[alias].alias_of===owner
        @test A.SESSION_VIEWS[202].alias_of===owner
        @test first(A.view.probe(201,0.,0.,0.))==[1.]
        A.view.add_list_data(source,"SP",1,[0.,0.,0.,2.])
        @test A.view.get_list_data(alias)==A.view.get_list_data(source)
        @test A.view.get_list_data(202)==A.view.get_list_data(source)
        # A fatal script can mutate staged owners and create plugin output;
        # neither operation may detach a held alias or publish partial data.
        before=A.view.get_tags();values=A.view.get_list_data(source)
        write(geo,"Plugin(Probe).View=100;Plugin(Probe).Run;Delete View[0];Exit 2;")
        @test_throws ArgumentError A.merge(geo)
        @test A.view.get_tags()==before && A.view.get_list_data(source)==values
        @test A.SESSION_VIEWS[alias].alias_of===owner
        A.view.remove(source)
        @test first(A.view.probe(alias,0.,0.,0.))==[2.]
        @test first(A.view.probe(202,0.,0.,0.))==[2.]
    end
end

@testset "GEO plugin model views publish replacements and live source slots" begin
    line=raw"""
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
2
11 0 0 0
21 1 0 0
$EndNodes
$Elements
1
31 1 2 0 1 11 21
$EndElements
"""
    for import_file in (A.open,A.merge)
        session() do dir
            retained=nothing;retained_owner=nothing;retained_slot=nothing
            if import_file===A.merge
                # Merge into a populated model whose already held model view
                # must keep its source slot while the new model view is bound.
                existing=replace(line,"11 0 0 0"=>"101 0 1 0",
                    "21 1 0 0"=>"121 1 1 0","31 1 2 0 1 11 21"=>"131 1 2 0 1 101 121")
                existing_path=joinpath(dir,"existing.msh");write(existing_path,existing)
                A.open(existing_path)
                retained=A.view.add("retained model",250)
                A.view.add_model_data(retained,0,A.model.get_current(),"NodeData",
                    [101,121],[[3.],[3.]])
                retained_owner=A.SESSION_VIEWS[retained]
                retained_slot=only(values(retained_owner.model_contexts)).model_slot
            end
            source=A.view.add("old source",100)
            A.view.add_list_data(source,"SP",1,[0.,0.,0.,1.])
            alias=A.view.add_alias(source,false,200);old=A.SESSION_VIEWS[source]
            write(joinpath(dir,"line.msh"),line)
            geo=joinpath(dir,"model-view.geo")
            create="Merge \"line.msh\";Plugin(NewView).ViewTag=100;" *
                "Plugin(NewView).Value=7;Plugin(NewView).Run;"
            # Replacing a tag publishes a new view at that tag; a held alias
            # keeps the old owner's data, as direct Plugin(NewView) does.
            write(geo,create);import_file(geo)
            @test A.SESSION_VIEWS[source]!==old
            @test A.SESSION_VIEWS[source].model_kind=="NodeData"
            @test Set(A.view.get_model_data(source,0)[2])==
                Set(import_file===A.merge ? UInt64[11,21,101,121] : UInt64[11,21])
            @test A.SESSION_VIEWS[alias].alias_of===old
            @test first(A.view.probe(alias,0.,0.,0.))==[1.]
            context=only(values(A.SESSION_VIEWS[source].model_contexts))
            @test context.model_slot===A.MODEL_SLOTS[A._current_slot_index_locked()]
            @test first(A.view.probe(source,.25,0.,0.))==[7.]
            if retained!==nothing
                @test A.SESSION_VIEWS[retained]===retained_owner
                @test only(values(retained_owner.model_contexts)).model_slot===retained_slot
                @test first(A.view.probe(retained,.25,1.,0.))==[3.]
            end
            held=A.SESSION_VIEWS[source];data=A.view.get_model_data(source,0)
            write(geo,"Plugin(NewView).ViewTag=100;Plugin(NewView).Value=9;" *
                "Plugin(NewView).Run;Exit 2;")
            @test_throws ArgumentError A.merge(geo)
            @test A.SESSION_VIEWS[source]===held && A.view.get_model_data(source,0)==data
        end
        session() do dir
            write(joinpath(dir,"line.msh"),line)
            geo=joinpath(dir,"model-alias.geo")
            write(geo,"Merge \"line.msh\";Plugin(NewView).Value=5;" *
                "Plugin(NewView).Run;Alias View[0];Delete View[0];")
            import_file(geo)
            alias=only(A.view.get_tags());owner=A._view_data_owner(A.SESSION_VIEWS[alias])
            @test A.SESSION_VIEWS[alias].alias_of===owner
            @test only(values(owner.model_contexts)).model_slot===
                A.MODEL_SLOTS[A._current_slot_index_locked()]
            @test first(A.view.probe(alias,.25,0.,0.))==[5.]
        end
    end
end
end
