# `Include "file"` parses the file as `.geo` whatever its extension (Gmsh's
# `tInclude` always calls `ParseFile`), as in examples/simple_geo/tower.geo
# (`Include "tower.i1"`).
using Test
using Tessella

@testset "geo Include of non-.geo extensions" begin
    dir=mktempdir()
    write(joinpath(dir,"part.i1"),"a = 3;\nPoint(1) = {a,0,0,0.1};\n")
    path=joinpath(dir,"main.geo")
    write(path,"Include \"part.i1\";\nb = a + 1;\n")
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.values["b"]==4
    @test ex.model.points[1]==(3.0,0.0,0.0)
end
