module ViewRepresentationFunctionalChecks
using Test,Tessella
const A=Tessella.API
const P=Tessella.Post
metadata(tag)=Tuple(A.view.get_number(tag,key) for key in ("NbTimeStep","TimeStep","Time"))
geometry()=(A.model.list(),A.model.get_current(),A.model.get_entities(),
    A.mesh.get_nodes(),A.mesh.get_elements())

function point_model()
    A.model.add("point");A.model.add_discrete_entity(0,17)
    A.mesh.add_nodes(0,17,[501],[2.,3,4])
    A.mesh.add_elements_by_type(17,15,[701],[501])
    A.mesh.rebuild_node_cache(false)
end

function exercise(directory)
    @testset "list and model representation replacement retains held slots" begin
        A.initialize()
        try
            A.option("General.Terminal",0);point_model()
            held=A.view.add("held")
            A.view.add_model_data(held,0,"point","NodeData",[501],[[99.]],.25)
            before=geometry();tag=A.view.add("transition")
            A.view.set_string(tag,"Name","actual name")
            A.view.set_number(tag,"Visible",0)
            payload=[0.,0,0,10,20];A.view.add_list_data(tag,"SP",1,payload)
            payload[end]=999.
            @test A.view.get_list_data(tag)==(["SP"],[[0.,0,0,10,20]])
            A.view.set_number(tag,"TimeStep",1)
            @test metadata(tag)==(2,1,1)
            A.view.add_model_data(tag,3,"point","NodeData",[501],[[30.]],7.5)
            @test metadata(tag)==(4,1,1) &&
                A.view.get_model_data(tag,3)==(["NodeData"],UInt64[501],[30.])
            @test_throws ArgumentError A.view.get_list_data(tag)
            @test A.view.get_model_data(tag,1)==(["NodeData"],UInt64[],Float64[])
            for binary in (0,1)
                A.option("PostProcessing.Binary",binary)
                path=joinpath(directory,"model-$binary.msh");A.view.write(tag,path)
                mesh=Tessella.read_mixed_msh(path);section=only(mesh.data_sections)
                @test mesh.coords==reshape([2.,3,4],3,1) &&
                    mesh.entity_data.node_entities==[(0,17)] &&
                    mesh.entity_data.external_node_tags==UInt64[501] &&
                    section.name=="NodeData" && section.reals==[7.5] &&
                    section.header[1]==3 && section.values==[30.]
            end
            A.view.set_number(tag,"TimeStep",3)
            A.view.add_list_data(tag,"SP",1,[0.,0,0,40,50])
            @test metadata(tag)==(2,3,7.5) &&
                A.view.get_list_data(tag)==(["SP"],[[0.,0,0,40,50]]) &&
                A.view.get_string(tag,"Name")=="actual name" &&
                A.view.get_number(tag,"Visible")==0
            @test_throws ArgumentError A.view.get_model_data(tag,3)
            A.option("PostProcessing.Binary",0)
            path=joinpath(directory,"replaced.pos");A.view.write(tag,path)
            parsed=only(P.read_pos(path))
            @test parsed.name=="actual name" && parsed.time==[0.,1] &&
                only(parsed.elements).values==reshape([40.,50.],1,1,2)
            path=joinpath(directory,"replaced.csv");A.view.write(tag,path)
            rows=[parse.(Float64,split(line,',')) for line in readlines(path)]
            @test length(rows)==2 && first.(rows)==[0.,1] && last.(rows)==[40.,50.]
            returned=A.view.get_list_data(tag);returned[2][1][end]=999.
            A.view.set_number(tag,"TimeStep",1)
            @test metadata(tag)==(2,1,1) && last(only(A.view.get_list_data(tag)[2]))==50
            @test geometry()==before && A.view.get_tags()==[held,tag] &&
                A.view.get_model_data(held,0)==(["NodeData"],UInt64[501],[99.])
        finally;A.finalize();end
    end
    @testset "physical list history resets only on representation replacement" begin
        A.initialize()
        try
            A.option("General.Terminal",0);point_model()
            path=joinpath(directory,"physical.pos")
            Base.write(path,"View \"physical\" {TIME{0.25,1.5};SP(0,0,0){10,20};};\n")
            A.merge(path);tag=only(A.view.get_tags())
            A.view.set_number(tag,"TimeStep",1)
            A.view.add_model_data(tag,0,"point","ElementData",[701],[[42.]],4.)
            @test metadata(tag)==(1,1,1.5) &&
                A.view.get_model_data(tag,0)==(["ElementData"],UInt64[701],[42.])
            A.view.add_list_data(tag,"SP",1,[0.,0,0,40,50])
            path=joinpath(directory,"reset-history.pos");A.view.write(tag,path)
            @test only(P.read_pos(path)).time==[0.,1]
            A.view.add_model_data(tag,3,"point","NodeData",[501],[[30.]],7.5)
            A.view.set_number(tag,"TimeStep",3)
            A.view.add_list_data(tag,"SP",0,Float64[])
            @test metadata(tag)==(0,3,7.5) && all(isempty,A.view.get_list_data(tag)[2])
            @test_throws ArgumentError A.view.get_model_data(tag,3)
        finally;A.finalize();end
    end
    @testset "failed representation admission preserves actual view data" begin
        A.initialize()
        try
            A.option("General.Terminal",0);point_model();tag=A.view.add("atomic")
            A.view.add_list_data(tag,"SP",1,[0.,0,0,10,20])
            A.view.set_number(tag,"TimeStep",1)
            before=(metadata(tag),A.view.get_list_data(tag),geometry())
            @test_throws ArgumentError A.view.add_model_data(tag,0,"missing","NodeData",[501],[[1.]],0.)
            @test (metadata(tag),A.view.get_list_data(tag),geometry())==before
            A.view.add_model_data(tag,3,"point","NodeData",[501],[[30.]],7.5)
            A.view.set_number(tag,"TimeStep",3)
            before=(metadata(tag),A.view.get_model_data(tag,3),geometry(),A.view.get_tags())
            @test_throws ArgumentError A.view.add_list_data(tag,"XX",1,[0.,0,0,1])
            @test_throws ArgumentError A.view.add_list_data(tag,"ST2",1,zeros(12))
            @test_throws ArgumentError A.view.add_list_data(tag,"SP",1,[0.,0,0])
            @test_throws ArgumentError A.view.add_list_data(tag,"SP",1,[0.,0,0,NaN])
            @test (metadata(tag),A.view.get_model_data(tag,3),geometry(),A.view.get_tags())==before
        finally;A.finalize();end
    end
end

directory=get(ENV,"TESSELLA_VIEW_REPRESENTATION_DIR","")
if isempty(directory);mktempdir(exercise);else;mkpath(directory);exercise(directory);end
end
