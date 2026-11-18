# `DefineConstant` option tails: `Min -5` is an option item (not `Min - 5`
# glued onto the value), and a value read like `Geometry.Curves` is the
# value, not an option word (examples/boolean/coherence.geo, compsolid*.geo,
# transform.geo, examples/simple_geo/indheat.geo). Values printed by Gmsh 4.15.2.
using Test
using Tessella

@testset "geo DefineConstant option tails" begin
    dir=mktempdir()
    path=joinpath(dir,"defconst.geo")
    write(path,"""
        DefineConstant[
          x = {0, Min -5, Max 5, Step 0.1, Name "Bloc 1/0x"}
          dx = {2, Min 0.1, Max 5, Step 0.1, Name "Bloc 1/dx"}
          x2 = {x+dx, Min -5, Max 5, Step 0.1, Name "Bloc 2/0x"}
        ];
        DefineConstant[
          showLines = {Geometry.Curves, Choices {0,1}, Name "Options/Show curves?",
                       GmshOption "Geometry.Curves", AutoCheck 0}
          macro = {"title.script", Name "Macros/Add title",
                   Macro "Gmsh", AutoCheck 0 /*, Highlight "Turquoise"*/ }
          a4 = {dx+1, Name "A/a4"}
        ];
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test (ex.values["x"],ex.values["dx"],ex.values["x2"])==(0.0,2.0,2.0)
    @test ex.values["showLines"]==1.0
    @test ex.strings["macro"]==["title.script"]
    @test ex.values["a4"]==3.0
end
