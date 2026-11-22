# examples/boolean/primitives.geo in small, against Gmsh 4.15.2.
#
# The allocator observer re-evaluated a primitive's parameter list after the
# statement had run, so `Box(newv) = {x++,...}` advanced `x` twice (the file's
# shapes ran to x = 27.8 where Gmsh stops at 14.8); `Spline`/`Bezier` under
# the OpenCASCADE factory were not recorded, so the next `newl` named the
# curve just made ("OpenCASCADE curve with tag 131 already exists"). The
# side planes of a partial cylinder measured 0.375 where Gmsh has 0.25 (their
# edges are stored out of traversal order and the shoelace added a closing
# chord); `occ.get_mass` did not reach shapes a `.geo` statement left in the
# model only, and a `Plane Surface` with a Bezier side had no plane record.
module ApiOccPrimitivesGeoGmshTest
using Test, Tessella

@testset "primitives.geo in small (Gmsh 4.15.2)" begin
    A=Tessella.API; occ=A.model.occ
    path=tempname()*".geo"
    write(path,"""
        SetFactory("OpenCASCADE");
        x = 0; y = 0;
        Box(newv) = {x++,y,0, 0.5,0.5,0.5};
        Cylinder(newv) = {x++,y,0, 0.5,0,0, 0.5, Pi/3};
        p = newp;
          Point(p) = {x++,y,0}; Point(p+1) = {x-0.7,y+0.5,0};
          Point(p+2) = {x-0.3,y+0.5,0};
          Point(p+3) = {x-0.1,y,0}; Point(p+4) = {x-0.9,y-0.2,0};
        l = newl;
          Bezier(l) = {p,p+4,p+3,p+2}; Line(l+1) = {p+2,p+1};
          Line(l+2) = {p+1,p};
        ll = newll;
          Curve Loop(ll) = {l:l+2};
        Plane Surface(news) = {ll};
        p = newp; Point(p) = {x++,y,0}; Point(p+1) = {x-0.5,y+0.3,0};
        Point(p+2) = {x-0.2,y,0};
        Spline(newl) = {p:p+2};
        p = newp; Point(p) = {x++,y,0}; Point(p+1) = {x-0.5,y+0.3,0};
        Point(p+2) = {x-0.2,y,0};
        Bezier(newl) = {p:p+2};
        """)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        @test [length(A.model.get_entities(d)) for d in 0:3]==[25,26,12,2]
        @test A.model.get_bounding_box(-1,-1)[4]≈4.8 atol=1e-6
        @test occ.get_mass(3,1)≈0.125 rtol=1e-12
        @test occ.get_mass(3,13)≈0.0654498469497874 rtol=1e-12
        @test [occ.get_mass(2,t) for t in 7:11]≈[0.261799387799149,
            0.130899693899575,0.130899693899575,0.25,0.25] rtol=1e-12
        @test occ.get_mass(2,26)≈0.2905 rtol=1e-12
        @test sort!([t for (_,t) in A.model.get_entities(1)])[end-1:end]==[27,28]
    finally
        A.finalize()
        rm(path;force=true)
    end
end
end
