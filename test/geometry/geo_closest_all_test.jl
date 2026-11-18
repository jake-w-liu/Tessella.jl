# `Closest {x,y,z}{ shapes }` returns every given entity sorted by distance
# (`getClosestEntities(..., n = dimTags.size())`), not just the nearest one —
# tutorials/t16.geo reads `closest(1)`. Expected values from Gmsh 4.15.2.
using Test
using Tessella

@testset "geo Closest returns all entities sorted" begin
    dir=mktempdir()
    path=joinpath(dir,"closest.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,1,2,3};
        c() = Closest {0.5,1,10} { Surface{1:6}; };
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    c=ex.lists["c"]
    @test length(c)==6
    @test sort(c)==collect(1.0:6.0)
    @test (c[1],c[end])==(6.0,5.0)
end
