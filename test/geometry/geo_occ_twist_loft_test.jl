# examples/boolean/twist.geo up to its Boolean: three sections of four
# ellipse arcs, each rotated by Pi/8, lofted by `ThruSections(1)`. The arcs
# of one section are wound both ways about the cap's normal, and the cap's
# plane fit stopped with "boundary conics have inconsistent orientations";
# a conic wound against the face axis now traces its projection backwards
# (a `neg2d` pcurve). Gmsh 4.15.2: 6 faces, 1 volume of 0.0706312232625922
# (OpenCASCADE's approximated sections; Tessella interpolates them
# exactly: the caps are pi*0.3*0.075, Gmsh's 1.7e-8 below).
module GeoOccTwistLoftTest
using Test, Tessella

@testset "twist.geo loft (Gmsh 4.15.2)" begin
    A=Tessella.API; occ=A.model.occ
    path=tempname()*".geo"
    write(path,"""
        SetFactory("OpenCASCADE");
        a = 0.300; b = 0.075;
        Point(1) = {a*Cos(Pi/6), a*Sin(Pi/6), 0};
        Point(2) = {-b*Sin(Pi/6), b*Cos(Pi/6), 0};
        Point(3) = {-a*Cos(Pi/6), -a*Sin(Pi/6), 0};
        Point(4) = {b*Sin(Pi/6), -b*Cos(Pi/6), 0};
        Point(5) = {0, 0, 0};
        Ellipse(1) = {3, 5, 3, 2}; Ellipse(2) = {3, 5, 3, 4};
        Ellipse(3) = {1, 5, 1, 4}; Ellipse(4) = {1, 5, 1, 2};
        l~{1}() = {1:4};
        Curve Loop(1) = l~{1}();
        N = 3; angle = Pi/4;
        For i In {2:N}
          l~{i}() = Translate{0,0,1/(N-1)}{ Duplicata{ Curve{l~{i-1}()}; } };
          Rotate {{0, 0, 1}, {0, 0, 0}, angle/(N-1)} { Curve{l~{i}()}; }
          Curve Loop(i) = l~{i}();
        EndFor
        ThruSections(1) = {1:N};
        """)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        @test length(A.model.get_entities(2))==6
        @test length(A.model.get_entities(3))==1
        @test occ.get_mass(3,1)≈0.0706312232625922 rtol=1e-7
        @test occ.get_mass(2,5)≈π*0.3*0.075 rtol=1e-12
        @test occ.get_mass(2,1)≈0.323743630051411 rtol=1e-7
    finally
        A.finalize()
        rm(path;force=true)
    end
end
end
