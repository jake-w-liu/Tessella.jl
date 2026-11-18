# `Field[i] = Extend` accepts Gmsh 4.15.2's ExtendField options (CurvesList,
# SurfacesList, DistMax, SizeMax, Power); DistMax/Power were rejected as
# unknown options and silently dropped (examples/boolean/full_quad.geo,
# extend_field.geo), while BoundaryLayer options were accepted.
using Test
using Tessella

@testset "geo Extend field options" begin
    dir=mktempdir()
    path=joinpath(dir,"extend.geo")
    write(path,"""
        Point(1) = {0,0,0,0.1}; Point(2) = {1,0,0,0.1};
        Line(1) = {1,2};
        Field[1] = Extend;
        Field[1].CurvesList = {1};
        Field[1].DistMax = 0.5;
        Field[1].SizeMax = 0.2;
        Field[1].Power = 2;
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    opts=ex.params.fields[1].options
    @test parse(Float64,opts["DistMax"])==0.5
    @test parse(Float64,opts["Power"])==2.0
    bad=joinpath(dir,"extend_bad.geo")
    write(bad,"Field[1] = Extend;\nField[1].Thickness = 0.5;\n")
    ex=execute_geo(bad)
    @test any(e->occursin("Unknown option 'Thickness' in field 1 of type 'Extend'",e),
              ex.errors)
end
