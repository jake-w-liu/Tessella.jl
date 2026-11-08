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
end
