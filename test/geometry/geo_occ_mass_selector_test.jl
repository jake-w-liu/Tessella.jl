# `Mass`/`CenterOfMass`/`MatrixOfInertia X{t}` under the OpenCASCADE factory
# read the OCC mass properties (examples/boolean/mesh_size_per_volume.geo,
# number_of_tets.geo) instead of aborting. Values printed by Gmsh 4.15.2.
using Test
using Tessella

@testset "geo OCC Mass/CenterOfMass/MatrixOfInertia" begin
    dir=mktempdir()
    path=joinpath(dir,"mass.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,1,2,3};
        v = Mass Volume{1};
        s = Mass Surface{1};
        c() = CenterOfMass Volume{1};
        i() = MatrixOfInertia Volume{1};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.values["v"]≈6.0
    @test ex.values["s"]≈6.0
    @test ex.lists["c"]≈[0.5,1.0,1.5]
    i=ex.lists["i"]
    @test length(i)==9
    @test [i[1],i[5],i[9]]≈[6.5,5.0,2.5]
    @test abs(i[2])<1e-12
end
