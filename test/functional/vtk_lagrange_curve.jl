using Test, Tessella

@testset "VTK keeps all cubic Curve nodes and values" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("cubic curve");A.model.add_discrete_entity(1,7)
        coordinates=[0.,0.,0.,1.,1.,0.,1/3,1/9,0.,2/3,4/9,0.]
        tags=[101,102,103,104]
        A.mesh.add_nodes(1,7,tags,coordinates);A.mesh.add_elements_by_type(7,26,[301],tags)
        view=A.view.add("curve nodal values")
        A.view.add_model_data(view,0,"cubic curve","NodeData",tags,[[value] for value in (0.,1.,1/3,2/3)])
        mktempdir() do directory
            path=joinpath(directory,"curve.vtk");A.view.write(view,path);lines=readlines(path)
            cell=split(lines[findfirst(startswith("CELLS "),lines)+1]);slots=parse.(Int,cell[2:end]).+1
            field(name,T)=let row=findfirst(startswith(name*" "),lines)
                parse.(T,split(join(lines[row+1:row+4],' ')))
            end
            @test lines[findfirst(startswith("CELL_TYPES "),lines)+1]=="68" && parse(Int,first(cell))==4
            @test field("node_uid",UInt64)[slots]==UInt64.(tags) && field("step_0",Float64)[slots]==[0.,1.,1/3,2/3]
        end
    finally
        A.finalize()
    end
end
