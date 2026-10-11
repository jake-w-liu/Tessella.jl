using Test, Tessella

@testset "VTK preserves every Pyramid13 node" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("pyramid13");A.model.add_discrete_entity(3,7)
        _,_,_,count,coordinates,_=A.mesh.get_element_properties(19)
        tags=collect(101:100+count)
        A.mesh.add_nodes(3,7,tags,coordinates);A.mesh.add_elements_by_type(7,19,[301],tags)
        view=A.view.add("pyramid values")
        A.view.add_model_data(view,0,"pyramid13","NodeData",tags,[[Float64(index)] for index in 1:count])
        mktempdir() do directory
            path=joinpath(directory,"pyramid.vtk");A.view.write(view,path);lines=readlines(path)
            cell=split(lines[findfirst(startswith("CELLS "),lines)+1]);slots=parse.(Int,cell[2:end]).+1
            field(name,T)=let row=findfirst(startswith(name*" "),lines)
                parse.(T,split(join(lines[row+1:row+count],' ')))
            end
            expected_order=[1,2,3,4,5,6,9,11,7,8,10,12,13]
            @test lines[findfirst(startswith("CELL_TYPES "),lines)+1]=="27" && parse(Int,first(cell))==13
            @test field("node_uid",UInt64)[slots]==UInt64.(tags[expected_order]) && field("step_0",Float64)[slots]==Float64.(expected_order)
        end
    finally
        A.finalize()
    end
end
