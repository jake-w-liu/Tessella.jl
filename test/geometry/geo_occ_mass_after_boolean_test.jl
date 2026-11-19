module GeoOCCMassAfterBooleanChecks
using Test,Tessella

# `Mass Volume{…}` after a `.geo` Boolean threw
# `MethodError: no method matching haskey(::Nothing, ::Symbol)`
# (examples/boolean/number_of_tets.geo): the result faces keep their pcurves
# per use, so an edge used only reversed has an empty forward slot, and the
# faces were bound into the OCC registry without the normalization the API
# path applies. A planar Boolean shows it; the value is exact.
@testset "geo Mass after an OCC Boolean" begin
    dir=mktempdir()
    path=joinpath(dir,"mass_after_boolean.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,2,2,2};
        Box(2) = {1,1,1,2,2,2};
        BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
        v = Mass Volume{3};
        c() = CenterOfMass Volume{3};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.values["v"]≈7.0 rtol=1e-12
    @test ex.lists["c"]≈fill((8-1.5)/7,3) rtol=1e-12
end
end
