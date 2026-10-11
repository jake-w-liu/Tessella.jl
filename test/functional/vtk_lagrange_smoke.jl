@testset "public cubic Triangle VTK export" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0);A.model.add("cubic-vtk")
        _,_,_,count,local_coordinates,_=A.mesh.get_element_properties(21)
        coordinates=vcat(reshape(local_coordinates,2,count),zeros(1,count))
        A.model.add_discrete_entity(2,7);tags=UInt64.(1:count)
        A.mesh.add_nodes(2,7,tags,vec(coordinates))
        A.mesh.add_elements_by_type(7,21,[5000],tags)
        view=A.view.add("values")
        A.view.add_model_data(view,0,"cubic-vtk","NodeData",tags,[[Float64(tag)] for tag in tags],0.,1)
        mktempdir() do directory
            path=joinpath(directory,"tri.vtk");A.view.write(view,path)
            text=replace(read(path,String),"\r\n"=>"\n")
            @test occursin("CELL_TYPES 1\n69\n",text)
            @test occursin("POINTS 10 double\n",text) && occursin("step_0 1 10 double\n",text)
        end
    finally
        A.finalize()
    end
end
