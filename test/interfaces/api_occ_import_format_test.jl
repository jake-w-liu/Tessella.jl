module OCCImportFormatChecks
using Test, Tessella
const A=Tessella.API
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

@testset "OCC explicit format selects the native reader" begin
    mktempdir() do directory
        for (source,format,extension) in (("iges_occ_curves.igs","iges",".stp"),
                                          ("brep_occ_box.step","step",".igs"))
            for suffix in (".dat",extension)
                path=joinpath(directory,"payload"*suffix)
                cp(joinpath(FIXTURES,source),path)
                A.initialize()
                try
                    A.option("General.Terminal",0)
                    out=A.model.occ.import_shapes(path;format)
                    @test !isempty(out)
                    @test all(first(dt)==(format=="iges" ? 1 : 3) for dt in out)
                finally
                    A.finalize()
                end
                rm(path)
            end
        end
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.model.occ.add_point(2.,3.,4.)
            before=A.model.occ.get_entities()
            for invalid in ("unknown","igs","stp","IGES",".step",7)
                @test_throws ArgumentError A.model.occ.import_shapes(
                    joinpath(FIXTURES,"brep_occ_box.step");format=invalid)
                @test A.model.occ.get_entities()==before
            end
            @test_throws ArgumentError A.model.occ.import_shapes(
                joinpath(FIXTURES,"iges_occ_curves.igs");format="step")
            @test A.model.occ.get_entities()==before
            @test_throws ArgumentError A.model.occ.import_shapes(
                joinpath(FIXTURES,"brep_occ_box.step");format="iges")
            @test A.model.occ.get_entities()==before
        finally
            A.finalize()
        end
    end
end
end
