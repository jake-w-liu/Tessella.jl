# Field ids/options computed from model-dependent lists (`all() = Volume{:};
# Field[#all() + 1] = Min;` in examples/boolean/mesh_size_per_volume.geo) used
# to abort `execute_geo` inside its static parameter pre-scan. Execution builds
# the field map itself; `read_geo_params` alone still reports the list as
# unavailable.
using Test
using Tessella

@testset "geo Field ids from model lists" begin
    dir=mktempdir()
    path=joinpath(dir,"fields.geo")
    write(path,"""
        Point(1) = {0,0,0,0.1}; Point(2) = {1,0,0,0.1}; Point(3) = {1,1,0,0.1};
        all() = Point{:};
        For i In {0 : #all() - 1}
          Field[i + 1] = Constant;
          Field[i + 1].VIn = 0.1*(i+1);
        EndFor
        Field[#all() + 1] = Min;
        Field[#all() + 1].FieldsList = {1:#all()};
        Background Field = #all() + 1;
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test sort(collect(keys(ex.params.fields)))==[1,2,3,4]
    @test ex.params.fields[4].kind=="Min"
    @test_throws ArgumentError Tessella.IO.read_geo_params(path)
end
