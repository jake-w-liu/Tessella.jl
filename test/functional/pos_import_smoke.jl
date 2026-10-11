@testset "session POS file import" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            path=joinpath(directory,"list.pos");Base.write(path,"View \"named\" {SP(1,2,3){4};};\n")
            A.open(path);tag=only(A.view.get_tags())
            @test A.model.list()==["list"] && A.view.get_string(tag,"Name")=="named" &&
                A.view.get_list_data(tag)==(["SP"],[[1.,2,3,4]])
            empty=joinpath(directory,"empty.pos");Base.write(empty,"// no views\n")
            before=(A.model.list(),A.model.get_current(),A.model.get_file_name(),A.view.get_tags())
            A.merge(empty)
            @test (A.model.list(),A.model.get_current(),A.model.get_file_name(),A.view.get_tags())==before
        end
    finally;A.finalize();end
end
