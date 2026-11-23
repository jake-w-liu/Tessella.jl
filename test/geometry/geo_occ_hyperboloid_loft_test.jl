# examples/boolean/hyperboloid.geo and twist.geo pieces, against Gmsh 4.15.2.
#
# - A plane face bounded by lines and a spline, one line on the axis, revolved
#   a full turn: the OCC revolve rejected the spline edge and the statement
#   fell back to the built-in sweep (a Volume of mass 0, the profile kept).
#   Gmsh: 2 disks and a "Surface of Revolution". (hyperboloid.geo's 20-point
#   spline: volume 16.75519028 by OpenCASCADE's GProp, 16.7551738395 by
#   2000x20-point Gauss quadrature of pi*x^2 dy; Tessella gives the latter.)
# - A smooth loft through three sections or more bounds its face with new
#   B-spline edges at the first and last section (the section curves stay
#   free); Tessella reused the section curves.
# - `Translate` rebinds one input at a time and a copy's sub-shapes take their
#   tags in the order `_bind` walks them: twist.geo's second section has
#   points 14-17 (was 6-9), translating two boxes renumbers their points
#   17-32.
# - A .geo OCC curve loop runs from its first curve (twist.geo's {1:4} is
#   1,-4,3,-2), so the loft's faces follow that order.
module GeoOccHyperboloidLoftTest
using Test, Tessella
const A = Tessella.API

function opened(src)
    path = tempname() * ".geo"
    write(path, src)
    A.initialize(); A.option("General.Terminal", 0)
    A.open(path); rm(path; force=true)
    return A.model
end

@testset "hyperboloid.geo / twist.geo pieces (Gmsh 4.15.2)" begin
    try
        m = opened("""
            SetFactory("OpenCASCADE");
            For i In {0:4}
              u = -1 + 2 / 4 * i;
              Point(1+i) = {Sqrt(1 + u^2), 2 * u, 0};
            EndFor
            Point(100) = {0, -2, 0};
            Point(101) = {0, 2, 0};
            Line(1) = {100, 1};
            Spline(2) = {1:5};
            Line(3) = {5, 101};
            Line(4) = {101, 100};
            Curve Loop(1) = {1:4};
            Plane Surface(1) = {1};
            Extrude{{0,1,0}, {0,0,0}, 2*Pi} { Surface{1}; }
            """)
        @test [length(m.get_entities(d)) for d in 0:3] == [7, 6, 4, 1]
        @test m.get_type(1, 2) == "BSpline"
        @test [m.get_type(2, t) for t in 2:4] ==
            ["Plane", "Surface of Revolution", "Plane"]
        # Gauss quadrature of pi*x^2 dy on the spline; Gmsh 16.7332366562
        @test m.occ.get_mass(3, 1) ≈ 16.7332366424 rtol=1e-9
    finally
        A.finalize()
    end
    try
        m = opened("""
            SetFactory("OpenCASCADE");
            For i In {0:2}
              Point(2*i+1) = {0, 0, i};
              Point(2*i+2) = {1 + 0.1 * i, 0, i};
              Line(i+1) = {2*i+1, 2*i+2};
              Wire(i+1) = {i+1};
            EndFor
            ThruSections{1:3}
            """)
        @test [length(m.get_entities(d)) for d in 0:2] == [6, 7, 1]
        @test [sort(abs.([b for (_, b) in m.get_boundary([(1, t)], false, false)]))
               for t in 4:7] == [[1, 2], [2, 6], [5, 6], [1, 5]]
        @test all(m.get_type(1, t) == "BSpline" for t in 4:7)
        @test [abs(b) for (_, b) in m.get_boundary([(2, 1)], false, false)] ==
            [4, 5, 6, 7]
    finally
        A.finalize()
    end
    try
        m = opened("""
            SetFactory("OpenCASCADE");
            a = 0.300; b = 0.075;
            Point(1) = {a*Cos(Pi/6), a*Sin(Pi/6), 0};
            Point(2) = {-b*Sin(Pi/6), b*Cos(Pi/6), 0};
            Point(3) = {-a*Cos(Pi/6), -a*Sin(Pi/6), 0};
            Point(4) = {b*Sin(Pi/6), -b*Cos(Pi/6), 0};
            Point(5) = {0, 0, 0};
            Ellipse(1) = {3, 5, 3, 2};
            Ellipse(2) = {3, 5, 3, 4};
            Ellipse(3) = {1, 5, 1, 4};
            Ellipse(4) = {1, 5, 1, 2};
            Curve Loop(1) = {1:4};
            l2() = Translate{0,0,0.5}{ Duplicata{ Curve{1:4}; } };
            Rotate {{0, 0, 1}, {0, 0, 0}, Pi/8} { Curve{l2()}; }
            Curve Loop(2) = l2();
            l3() = Translate{0,0,0.5}{ Duplicata{ Curve{l2()}; } };
            Rotate {{0, 0, 1}, {0, 0, 0}, Pi/8} { Curve{l3()}; }
            Curve Loop(3) = l3();
            ThruSections(1) = {1:3};
            """)
        @test [t for (_, t) in m.get_entities(0)] ==
            [1, 2, 3, 4, 5, 14, 15, 16, 17, 26, 27, 28, 29]
        bnd(t) = [b for (_, b) in m.get_boundary([(1, t)], false, false)]
        @test bnd(5) == [14, 15] && bnd(8) == [17, 15]
        # faces along curves 1, -4, 3, -2: first-section edges 13, 17, 20, 23
        @test [sort(abs.(bnd(t))) for t in (13, 17, 20, 23)] ==
            [[2, 3], [1, 2], [1, 4], [3, 4]]
        @test m.get_type(2, 6) == "Plane"
    finally
        A.finalize()
    end
    A.initialize()
    try
        occ = A.model.occ
        occ.add_box(0, 0, 0, 1, 1, 1); occ.add_box(2, 0, 0, 1, 1, 1)
        occ.translate([(3, 1), (3, 2)], 0, 0, 1); occ.synchronize()
        @test [t for (_, t) in A.model.get_entities(0)] == collect(17:32)
        @test [t for (_, t) in A.model.get_entities(2)] == collect(13:24)
        # a mirror is a GTransform upstream: curves TrimmedCurve, faces
        # "BSpline surface" (a translation keeps Line and Plane)
        occ.mirror([(3, 1)], 1, 0, 0, 0); occ.synchronize()
        f1 = abs(first(A.model.get_boundary([(3, 1)], false, false))[2])
        f2 = abs(first(A.model.get_boundary([(3, 2)], false, false))[2])
        c1 = abs(first(A.model.get_boundary([(2, f1)], false, false))[2])
        c2 = abs(first(A.model.get_boundary([(2, f2)], false, false))[2])
        @test (f1, c1, f2, c2) == (25, 49, 19, 37)
        @test A.model.get_type(1, c1) == "TrimmedCurve"
        @test A.model.get_type(2, f1) == "BSpline surface"
        @test A.model.get_type(1, c2) == "Line" && A.model.get_type(2, f2) == "Plane"
    finally
        A.finalize()
    end
end

end
