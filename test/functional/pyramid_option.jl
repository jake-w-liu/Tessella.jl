let api=Tessella.API
    api.initialize()
    try
        @test api.option("Mesh.OptimizePyramids",-0.1)==-0.1
        @test api.option("Mesh.OptimizePyramids",0.7)==0.7
    finally
        api.finalize()
    end
end
