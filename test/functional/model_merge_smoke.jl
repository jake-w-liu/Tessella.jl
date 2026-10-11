@testset "current model file merge" begin
    Tessella.API.initialize()
    try
        Tessella.API.option("General.Terminal",0)
        Tessella.API.model.add("merge smoke")
        mktempdir() do directory
            file=joinpath(directory,"point.geo");Base.write(file,"Point(17)={1,2,3,1};\n")
            @test Tessella.API.merge(file)===nothing && Tessella.API.model.get_current()=="merge smoke"
            @test Tessella.API.model.get_entities()==[(0,17)] && Tessella.API.model.get_file_name()==""
        end
    finally;Tessella.API.finalize();end
end
