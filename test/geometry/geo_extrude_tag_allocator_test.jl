# `newp`-family reads after built-in `Extrude`/`Duplicata`/transform
# statements, and the automatic tags those statements allocate. Expected
# values were printed by Gmsh 4.15.2 on the same scripts.
using Test
using Tessella

function _geo_alloc_run(source::AbstractString)
    dir=mktempdir()
    path=joinpath(dir,"script.geo")
    write(path,source)
    return execute_geo(path)
end

@testset "geo extrude/duplicata tag allocators" begin
    # The statement observer cannot replay these statements; the counters
    # come from the model afterwards instead of becoming unavailable.
    ex=_geo_alloc_run("""
        Point(1) = {0,0,0,0.1}; Point(2) = {1,0,0,0.1};
        Line(1) = {1,2};
        e[] = Extrude {0,1,0} { Curve{1}; };
        p1 = newp; l1 = newl; s1 = news; r1 = newreg;
        q[] = Translate {5,0,0} { Duplicata { Surface{e[1]}; } };
        p2 = newp; l2 = newl; s2 = news;
        Rotate {{0,0,1},{0,0,0},Pi/2} { Duplicata { Point{1}; } }
        p3 = newp;
        """)
    @test isempty(ex.errors)
    @test [ex.values[k] for k in ("p1","l1","s1","r1")]==[5,6,6,6]
    @test ex.lists["q"]==[6.0]
    @test [ex.values[k] for k in ("p2","l2","s2","p3")]==[15,11,11,15]

    # Extruded/duplicated surfaces with holes: Gmsh's derived surfaces carry
    # no curve loops, so their hole wires must not bump `NEWREG()`, and an
    # explicit `Curve Loop(t)` may reuse a tag the model used internally.
    ex=_geo_alloc_run("""
        Point(1) = {0,0,0,0.1}; Point(2) = {3,0,0,0.1}; Point(3) = {3,3,0,0.1};
        Point(4) = {0,3,0,0.1}; Point(5) = {1,1,0,0.1}; Point(6) = {2,1,0,0.1};
        Point(7) = {2,2,0,0.1}; Point(8) = {1,2,0,0.1};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
        Curve Loop(1) = {1,2,3,4};
        Curve Loop(2) = {5,6,7,8};
        Plane Surface(1) = {1,2};
        e[] = Extrude {0,0,1} { Surface{1}; };
        r1 = newreg; ll1 = newll;
        d[] = Translate {10,0,0} { Duplicata { Surface{1}; } };
        r2 = newreg;
        Curve Loop(e[2]) = {5,6,7,8};
        Plane Surface(newreg) = {e[2]};
        s = newreg - 1;
        b[] = Boundary{ Surface{e[2]}; };
        """)
    @test isempty(ex.errors)
    @test ex.lists["e"][1:2]==[50.0,1.0]
    @test ex.lists["e"][3]==21.0
    @test (ex.values["r1"],ex.values["ll1"])==(51,51)
    @test ex.lists["d"]==[51.0]
    @test (ex.values["r2"],ex.values["s"])==(60,60)
    @test length(ex.lists["b"])==4
    @test ex.model.loops[21]==[5,6,7,8]
    @test ex.model.surfaces[60]==[21]
end
