# `ArcCos`/`ArcSin`/`ArcTan`/`ArcTan2` are lexer aliases of `Acos`/`Asin`/
# `Atan`/`Atan2` in Gmsh (used by examples/simple_geo/antenna.geo and
# machine.i1). Values printed by Gmsh 4.15.2 with `%.15g`.
using Test
using Tessella

@testset "geo Arc* function aliases" begin
    dir=mktempdir()
    path=joinpath(dir,"alias.geo")
    write(path,"a = ArcTan(1); b = ArcCos(0.5); c = ArcSin(1) + ArcTan2(1, -1);\n" *
               "t = ArcTan(8.e-3/2);\n")
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.values["a"]≈0.785398163397448
    @test ex.values["b"]≈1.0471975511966
    @test ex.values["c"]≈3.92699081698724
    @test ex.values["t"]≈0.00399997866687146
end
