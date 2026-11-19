module GeoOCCSurfaceBooleanTagChecks
using Test,Tessella

# `BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };` —
# a surface Boolean with an explicit result tag — threw
# `MethodError: no method matching _model_occ_boolean_face!(…, ::String, ::Int64)`:
# the face publisher was called with the tag as a fourth argument it never
# had. The plate with a square hole has area 9 − 1.
@testset "geo surface Boolean with an explicit tag" begin
    dir=mktempdir()
    path=joinpath(dir,"plate.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Rectangle(1) = {0,0,0,3,3};
        Rectangle(2) = {1,1,0,1,1};
        BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
        a = Mass Surface{3};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.values["a"]≈8.0 rtol=1e-12
end
end
