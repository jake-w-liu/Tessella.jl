# `Field[i].F = Sprintf(...)` (tutorials/t10.geo) is a StringExpr: Gmsh stores
# the formatted text. Tessella stored the raw `Sprintf(...)` source, which the
# MathEval parser then rejected ("unexpected character '\"'").
using Test
using Tessella

@testset "geo Field string options from StringExpr" begin
    dir=mktempdir()
    path=joinpath(dir,"matheval.geo")
    write(path,"""
        lc = 0.15;
        s = "F2 + x";
        Field[4] = Distance;
        Field[5] = MathEval;
        Field[5].F = Sprintf("F4^3 + %g", lc / 100);
        Field[6] = MathEval;
        Field[6].F = "cos(x)";
        Field[7] = MathEval;
        Field[7].F = StrCat(s, " + 1");
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    f=ex.params.fields
    @test f[5].options["F"]=="\"F4^3 + 0.0015\""
    @test f[6].options["F"]=="\"cos(x)\""
    @test f[7].options["F"]=="\"F2 + x + 1\""
end
