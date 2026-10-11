using Test, Tessella

@testset "VTK retains every complete quadratic hex node and value" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("hex27");A.model.add_discrete_entity(3,7)
        _,_,_,count,reference,_=A.mesh.get_element_properties(12)
        tags=collect(101:100+count)
        A.mesh.add_nodes(3,7,tags,reference);A.mesh.add_elements_by_type(7,12,[301],tags)
        view=A.view.add("hex node values")
        A.view.add_model_data(view,0,"hex27","NodeData",tags,[[Float64(index)] for index in 1:count])
        mktempdir() do directory
            path=joinpath(directory,"hex.vtk");A.view.write(view,path);lines=readlines(path)
            cell=split(lines[findfirst(startswith("CELLS "),lines)+1]);slots=parse.(Int,cell[2:end]).+1
            field(name,T)=let row=findfirst(startswith(name*" "),lines)
                parse.(T,split(join(lines[row+1:row+count],' ')))
            end
            order=[1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16,23,24,22,25,21,26,27]
            @test lines[findfirst(startswith("CELL_TYPES "),lines)+1]=="29" && parse(Int,first(cell))==27
            @test field("node_uid",UInt64)[slots]==UInt64.(tags[order])
            @test field("step_0",Float64)[slots]==Float64.(order)
        end
    finally
        A.finalize()
    end
end
