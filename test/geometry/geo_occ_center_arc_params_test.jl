# OCC `Circle(t) = {x,y,z,r,angle2}` (5 parameters) and
# `Ellipse(t) = {x,y,z,rx,ry,angle2}` (6 parameters) take angle1 = 0, as in
# Gmsh 4.15.2's grammar (examples/boolean/primitives.geo). Values printed by
# Gmsh on the same script.
using Test
using Tessella

@testset "geo OCC center-form arcs with one angle" begin
    dir=mktempdir()
    path=joinpath(dir,"arcs.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Circle(1) = {0,0,0, 0.3, Pi/3};
        Ellipse(2) = {2,0,0, 0.4,0.1, Pi/3};
        l1 = Mass Curve{1};
        l2 = Mass Curve{2};
        b1() = Boundary{ Curve{1}; };
        p() = Point{b1(1)};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.values["l1"]≈0.31415927 atol=1e-7
    @test ex.values["l2"]≈0.22830658 atol=1e-7
    @test length(ex.lists["b1"])==2
    @test ex.lists["p"][1:2]≈[0.15,0.25980762] atol=1e-7
end
