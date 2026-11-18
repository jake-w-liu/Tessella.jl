# Under SetFactory("OpenCASCADE"):
# - `Curve Loop` ignores curve signs and orients each connected edge itself
#   (`OCC_Internals::addWire` / `BRepBuilderAPI_MakeWire`), so `{1,2,3}`
#   with Line(2) = {3,2} is a valid loop;
# - the point form `Ellipse(t) = {start, center, major, end}` builds an OCC
#   ellipse arc record (examples/boolean/twist.geo).
# Values printed by Gmsh 4.15.2 on the same script.
using Test
using Tessella

@testset "geo OCC wire orientation and point-form ellipse" begin
    dir=mktempdir()
    path=joinpath(dir,"occwire.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Point(1) = {0,0,0}; Point(2) = {1,0,0}; Point(3) = {0,1,0};
        Line(1) = {1,2}; Line(2) = {3,2}; Line(3) = {3,1};
        Curve Loop(1) = {1,2,3};
        Plane Surface(1) = {1};
        a1 = Mass Surface{1};
        Point(5) = {0, 0, 2};
        Point(6) = {1, 0, 2};
        Point(7) = {0, 0.5, 2};
        Ellipse(5) = {6, 5, 6, 7};
        l5 = Mass Curve{5};
        b() = Boundary{ Curve{5}; };
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.model.loops[1]==[1,-2,3]
    @test ex.values["a1"]≈0.5
    @test ex.values["l5"]≈1.211056 atol=1e-6
    @test ex.lists["b"]==[6.0,7.0]
end
