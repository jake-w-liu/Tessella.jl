module ListViewExportSmokeChecks
using Test, Tessella
const A=Tessella.API
@testset "list export formats and validation" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        view=A.view.add("two-step point")
        A.view.add_list_data(view,"SP",1,[.125,.5,.75,1.25,2.5])
        bad=A.view.add("incomplete triangle")
        A.view.add_list_data(bad,"ST",1,zeros(10))
        mktempdir() do directory
            txt=joinpath(directory,"point.txt");csv=joinpath(directory,"point.csv")
            vtk=joinpath(directory,"point.vtk")
            for path in (txt,csv,vtk);A.view.write(view,path);end
            expected=[[0.,0.,0.,0.,.125,.5,.75,1.25],[1.,1.,0.,0.,.125,.5,.75,2.5]]
            txt_rows=[parse.(Float64,split(line)) for line in readlines(txt)]
            csv_rows=[parse.(Float64,split(line,',')) for line in readlines(csv)]
            @test txt_rows==expected && csv_rows==expected
            text=replace(read(vtk,String),"\r\n"=>"\n")
            @test startswith(text,"# vtk DataFile Version 2.0\n") &&
                occursin("CELLS 1 2\n1 0\nCELL_TYPES 1\n1\n",text) &&
                occursin("SCALARS step_1 double 1\nLOOKUP_TABLE default\n2.5\n",text)
            write(csv,"sentinel")
            @test_throws ArgumentError A.view.write(bad,csv)
            @test read(csv,String)=="sentinel"
        end
    finally
        A.finalize()
    end
end
end
