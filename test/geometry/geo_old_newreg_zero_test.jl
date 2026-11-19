# `Geometry.OldNewReg = 0`: Extrude/Duplicata copies draw each dimension's
# own `getMaxTag(dim) + 1` and coherence merges no longer reset the counters.
# Expected tags printed by Gmsh 4.15.2 for the same script (with the default
# OldNewReg = 1 the same script gives top surface 72 and newp 15).
module GeoOldNewRegZeroTest
using Test, Tessella

@testset "geo Geometry.OldNewReg = 0 tag allocation" begin
    dir=mktempdir()
    path=joinpath(dir,"oldnewreg0.geo")
    write(path,"""
        Geometry.OldNewReg = 0;
        Point(1) = {0,0,0,0.1}; Point(2) = {1,0,0,0.1};
        Point(3) = {1,1,0,0.1}; Point(4) = {0,1,0,0.1};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Curve Loop(10) = {1,2,3,4};
        Plane Surface(1) = {10};
        Physical Surface(50) = {1};
        a[] = Extrude {0,0,1} { Surface{1}; };
        n1[] = {newp, newl, news, newv, newreg};
        b[] = Translate {3,0,0} { Duplicata { Surface{1}; } };
        c[] = Translate {0,3,0} { Duplicata { Volume{a[1]}; } };
        n2[] = {newp, newl, news, newv, newreg};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.lists["a"]==[7.0,1.0,3.0,4.0,5.0,6.0]
    @test ex.lists["n1"]==[45.0,21.0,8.0,2.0,51.0]
    @test ex.lists["b"]==[8.0]
    @test ex.lists["c"]==[2.0]
    @test ex.lists["n2"]==[157.0,49.0,15.0,3.0,51.0]
end
end
