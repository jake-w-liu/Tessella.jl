module ModelViewExportSmokeChecks
using Test, Tessella
const A=Tessella.API
@testset "model CSV and VTK keep sparse quadratic data" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("model smoke");A.model.add_discrete_entity(1,7)
        A.mesh.add_nodes(1,7,[101,102,103],[0.,0.,0.,2.,0.,0.,1.,0.,0.])
        A.mesh.add_elements_by_type(7,8,[301],[101,102,103])
        view=A.view.add("partial line values")
        A.view.add_model_data(view,0,"model smoke","NodeData",[101],[[1.25]])
        before=(A.mesh.get_nodes(),A.view.get_model_data(view,0))
        mktempdir() do directory
            path=joinpath(directory,"line.vtk");A.view.write(view,path)
            text=replace(read(path,String),"\r\n"=>"\n")
            @test occursin("CELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n21\n",text)
            @test occursin("step_0 1 3 double\n1.25\nNaN\nNaN\n",text) &&
                occursin("present_0 1 3 int\n1\n0\n0\n",text)
            path=joinpath(directory,"line.csv");A.view.write(view,path)
            lines=readlines(path);header=split(first(lines),',')
            rows=[Dict(zip(header,split(line,','))) for line in lines[2:end]]
            @test count(row->row["record"]=="node",rows)==3 &&
                only(row["connectivity"] for row in rows if row["record"]=="cell")=="101;102;103" &&
                only(row["value_1"] for row in rows if row["record"]=="data")=="1.25"
        end
        @test (A.mesh.get_nodes(),A.view.get_model_data(view,0))==before
    finally
        A.finalize()
    end
end
end
