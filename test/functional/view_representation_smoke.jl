@testset "view representation replacement" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0);A.model.add("point")
        A.model.add_discrete_entity(0,17);A.mesh.add_nodes(0,17,[501],[2.,3,4])
        tag=A.view.add("replacement");A.view.add_list_data(tag,"SP",1,[0.,0,0,3,7])
        A.view.add_model_data(tag,0,"point","NodeData",[501],[[11.]],.25)
        @test A.view.get_number(tag,"NbTimeStep")==1 &&
            A.view.get_model_data(tag,0)==(["NodeData"],UInt64[501],[11.])
        A.view.add_list_data(tag,"SP",0,Float64[])
        @test A.view.get_number(tag,"NbTimeStep")==0 && all(isempty,A.view.get_list_data(tag)[2])
    finally;A.finalize();end
end
