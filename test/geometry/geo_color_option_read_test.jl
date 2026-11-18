# A color option may be assigned from another color option
# (`Geometry.Color.Surfaces = Geometry.Color.Points;` in tutorials/t3.geo);
# Gmsh 4.15.2 accepts it silently, Tessella reported "Unknown color".
using Test
using Tessella

@testset "geo color option read in ColorExpr" begin
    dir=mktempdir()
    path=joinpath(dir,"colors.geo")
    write(path,"""
        Geometry.Color.Points = {10,20,30};
        Geometry.Color.Surfaces = Geometry.Color.Points;
        Geometry.Color.Curves = Geometry.Color.Volumes;
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
end
