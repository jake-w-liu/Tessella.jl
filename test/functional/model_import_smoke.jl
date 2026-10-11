@testset "current model file open" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            path=joinpath(directory,"opened.geo");Base.write(path,"Point(17)={1,2,3,1};\n")
            @test A.open(path)===nothing && A.model.list()==["opened"] && A.model.get_entities()==[(0,17)]
            before=(A.model.list(),A.model.get_current(),A.model.get_entities(),A.option("Mesh.MeshSizeFactor"))
            @test_throws ArgumentError A.open(joinpath(directory,"missing.msh"))
            @test (A.model.list(),A.model.get_current(),A.model.get_entities(),A.option("Mesh.MeshSizeFactor"))==before
        end
    finally;A.finalize();end
end
