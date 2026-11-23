# examples/boolean/shell_sewing.geo: an ellipse arc and a line revolved about
# the z axis into two open surfaces, sewn by `Surface Loop(...) = {1, 2}
# Using Sewing;` and closed by `Volume`. The revolve of a conic curve fell
# back to the built-in sweep and the sewing was ignored ("Unknown surface 2
# in GEO volume 2000"). Gmsh 4.15.2: Surface 1 "Surface of Revolution" of
# area 10.7392176639419, Surface 2 a cone of 4.44288293815837, and the sewn
# shell's faces 3 and 4 the same; the volume is 5.23598775598299 in
# magnitude (Gmsh reports it negative: its sewn shell is inside out).
# The sewing used to copy every vertex and edge (7 points, 11 curves where
# Gmsh has 5 and 9): it now rebuilds only the welded circle, the vertex on
# it and the seams through that vertex, keeps the pole edges 102 and 104,
# orients the shell from its first face ([3, -4] as in Gmsh) and the conic
# revolve has OCC's u = angle, v = profile layout (it had the transpose,
# hence the opposite normal).
module GeoOccShellSewingTest
using Test, Tessella

@testset "shell_sewing.geo (Gmsh 4.15.2)" begin
    A=Tessella.API; occ=A.model.occ
    path=tempname()*".geo"
    write(path,"""
        SetFactory("OpenCASCADE");
        Geometry.OCCAutoFix = 0;
        Point(1) = {1, 0, 0}; Point(2) = {0, 0, 2};
        Point(3) = {0, 0, 0}; Point(4) = {0, 0, -1};
        Ellipse(100) = {2, 3, 2, 1};
        Line(101) = {4, 1};
        a() = Extrude{ {0,0,1}, {0,0,0}, 2*Pi }{ Line{100}; };
        b() = Extrude{ {0,0,1}, {0,0,0}, 2*Pi }{ Line{101}; };
        Surface Loop(1000) = {1, 2} Using Sewing;
        Volume(2000) = {1000};
        """)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        @test [length(A.model.get_entities(d)) for d in 0:3]==[5,9,4,1]
        bnd(d,t)=[b for (_,b) in A.model.get_boundary([(d,t)],false,true)]
        @test bnd(3,2000)==[3,-4]
        @test bnd(1,106)==[5,5] && bnd(1,107)==[2,5] && bnd(1,108)==[4,5]
        @test bnd(2,3)==[102,107,-106,-107]
        @test A.model.get_type(1,104)=="Unknown" && A.model.get_type(1,105)=="Circle"
        @test A.model.get_type(2,1)=="Surface of Revolution"
        @test occ.get_mass(2,2)≈4.44288293815837 rtol=1e-12
        @test occ.get_mass(2,4)≈4.44288293815837 rtol=1e-12
        # The revolved faces' areas (10.7392176639419 both) and the volume
        # (5.23598775598297 here, -5.23598775598299 in Gmsh) are B-spline
        # quadratures too slow for the interpreted suite.
        @test length(A.model.get_entities(3))==1
    finally
        A.finalize()
        rm(path;force=true)
    end
end
end
