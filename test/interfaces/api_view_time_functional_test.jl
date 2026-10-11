module ViewTimeFunctionalChecks
using Test,Tessella
const A=Tessella.API
const P=Tessella.Post
metadata(tag)=Tuple(A.view.get_number(tag,key) for key in ("NbTimeStep","TimeStep","Time"))
list_state(tag)=(metadata(tag),A.view.get_list_data(tag),A.model.list(),
    A.model.get_current(),A.model.get_file_name(),A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_tags())

function exercise(directory)
    @testset "derived list metadata and display selection" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            path=joinpath(directory,"physical.pos")
            Base.write(path,"View \"physical\" {TIME{0.25,1.5,4};SP(0,0,0){10,20,30};};\n")
            A.open(path);tag=only(A.view.get_tags());original=A.view.get_list_data(tag)
            @test metadata(tag)==(3,0,0)
            A.view.set_number(tag,"TimeStep",1.7);@test metadata(tag)==(3,1,1.5)
            A.view.set_number(tag,"TimeStep",-2);@test metadata(tag)==(3,2,4)
            A.view.set_number(tag,"TimeStep",9);@test metadata(tag)==(3,0,.25)
            A.view.set_number(tag,"TimeStep",-.5);@test metadata(tag)==(3,0,.25)
            A.view.set_number(tag,"Time",.875);@test metadata(tag)==(3,0,.25)
            A.view.set_number(tag,"Time",2.75);@test metadata(tag)==(3,1,1.5)
            A.view.set_number(tag,"Time",4);A.view.set_number(tag,"Time",-1)
            @test metadata(tag)==(3,2,4)
            A.view.set_number(tag,"NbTimeStep",99)
            @test metadata(tag)==(3,2,4) && A.view.get_list_data(tag)==original
            A.view.add_list_data(tag,"SP",1,[0.,0,0,100,200])
            @test metadata(tag)==(2,2,4)
            A.view.set_number(tag,"Time",4);@test metadata(tag)==(2,1,1.5)
            A.view.add_list_data(tag,"SP",1,[0.,0,0,10,20,30,40,50])
            output=joinpath(directory,"expanded.pos");A.view.write(tag,output)
            @test only(P.read_pos(output)).time==[.25,1.5,4,3,4]
            A.view.set_number(tag,"TimeStep",4);@test metadata(tag)==(5,4,4)
            A.view.add_list_data(tag,"SL",1,[0.,1,0,0,0,0,1,2,3,4,5,6,7,8])
            @test metadata(tag)==(4,4,4)
            A.view.add_list_data(tag,"SP",0,Float64[]);A.view.add_list_data(tag,"SL",0,Float64[])
            A.view.set_number(tag,"TimeStep",-1)
            @test metadata(tag)==(4,3,3)
        finally;A.finalize();end
    end
    @testset "empty views and physical-time ties" begin
        A.initialize()
        try
            A.option("General.Terminal",0);empty=A.view.add("empty")
            @test metadata(empty)==(0,0,0)
            A.view.set_number(empty,"TimeStep",-2);@test metadata(empty)==(0,-1,0)
            path=joinpath(directory,"tie.pos")
            Base.write(path,"View \"unordered\" {TIME{4,0.25,1.5};SP(0,0,0){10,20,30};};\n"*
                "View \"duplicate\" {TIME{0.25,1.5,1.5};SP(0,0,0){10,20,30};};\n"*
                "View \"negative\" {TIME{-4,-1,0.5};SP(0,0,0){10,20,30};};\n")
            A.merge(path);unordered,duplicate,negative=setdiff(A.view.get_tags(),[empty])
            A.view.set_number(unordered,"Time",.875);A.view.set_number(duplicate,"Time",1.5)
            @test metadata(unordered)==(3,1,.25) && metadata(duplicate)==(3,1,1.5)
            A.view.set_number(negative,"TimeStep",1);A.view.set_number(negative,"Time",-1)
            @test metadata(negative)==(3,1,-1)
            A.view.set_number(negative,"Time",0);@test metadata(negative)==(3,2,.5)
        finally;A.finalize();end
    end
    @testset "actual sparse model steps and data-kind replacement" begin
        A.initialize()
        try
            A.option("General.Terminal",0);A.model.add("point")
            A.model.add_discrete_entity(0,17);A.mesh.add_nodes(0,17,[501],[2.,3,4])
            A.mesh.add_elements_by_type(17,15,[701],[501]);tag=A.view.add("sparse")
            geometry=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements())
            A.view.add_model_data(tag,3,"point","NodeData",[501],[[30.]],7.5)
            @test metadata(tag)==(4,0,0)
            gap=A.view.get_model_data(tag,1)
            @test gap==(["NodeData"],UInt64[],Float64[])
            A.view.set_number(tag,"TimeStep",3)
            A.view.add_model_data(tag,3,"point","NodeData",[501],[[31.]],9.25)
            @test metadata(tag)==(4,3,7.5) && A.view.get_model_data(tag,3)[3]==[31.]
            A.view.add_model_data(tag,0,"point","NodeData",[501],[[10.]],.25)
            A.view.set_number(tag,"Time",0);@test metadata(tag)==(4,1,0)
            A.view.set_number(tag,"TimeStep",-1);@test metadata(tag)==(4,3,9.25)
            A.view.add_model_data(tag,2,"point","ElementData",[701],[[42.]],4.)
            @test metadata(tag)==(3,3,9.25)
            @test_throws ArgumentError A.view.get_model_data(tag,3)
            actual=A.view.get_model_data(tag,2);gap=A.view.get_model_data(tag,1)
            A.view.set_number(tag,"TimeStep",1)
            @test actual==(["ElementData"],UInt64[701],[42.]) &&
                gap==(["ElementData"],UInt64[],Float64[]) && metadata(tag)==(3,1,0) &&
                (A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements())==geometry
            # Owned range storage and closest-zero-gap selection stay bounded.
            A.view.add_model_data(tag,1_000_000,"point","NodeData",[501],[[99.]],9.)
            A.view.set_number(tag,"Time",0)
            @test metadata(tag)==(1_000_001,0,0)
            A.view.set_number(tag,"TimeStep",-1)
            @test metadata(tag)==(1_000_001,1_000_000,9)
            before=(metadata(tag),A.view.get_model_data(tag,1_000_000),geometry)
            @test_throws ArgumentError A.view.add_model_data(tag,typemax(Int),"point","NodeData",[501],[[1.]],1.)
            @test (metadata(tag),A.view.get_model_data(tag,1_000_000),
                (A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements()))==before
        finally;A.finalize();end
    end
    @testset "finite and integer-range errors are atomic" begin
        A.initialize()
        try
            A.option("General.Terminal",0);tag=A.view.add("finite")
            A.view.add_list_data(tag,"SP",1,[0.,0,0,10,20]);A.view.set_number(tag,"TimeStep",1)
            before=list_state(tag)
            @test_throws ArgumentError A.view.set_number(tag,"Time",NaN)
            @test_throws ArgumentError A.view.set_number(tag,"TimeStep",Float64(typemax(Int32))+1)
            @test list_state(tag)==before
        finally;A.finalize();end
    end
end

directory=get(ENV,"TESSELLA_VIEW_TIME_DIR","")
if isempty(directory);mktempdir(exercise);else;mkpath(directory);exercise(directory);end
end
