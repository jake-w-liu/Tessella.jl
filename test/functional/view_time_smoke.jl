@testset "view step metadata" begin
    A=Tessella.API;A.initialize()
    try
        tag=A.view.add("steps");A.view.add_list_data(tag,"SP",1,[0.,0,0,3,7])
        @test A.view.get_number(tag,"NbTimeStep")==2
        A.view.set_number(tag,"TimeStep",1.7)
        @test A.view.get_number(tag,"TimeStep")==1 && A.view.get_number(tag,"Time")==1
    finally;A.finalize();end
end
