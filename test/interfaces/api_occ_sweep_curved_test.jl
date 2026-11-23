# OCC extrusion and revolution of plane faces bounded by conics, the copy's
# entity order, and OCC faces imported next to derived built-in loops — all
# against Gmsh 4.15.2.
#
# `occ.extrude` covered only all-line faces and full disks along their
# normal, `occ.revolve` only all-line faces: examples/boolean/extrude.geo
# and revolve.geo (fragments of squares and an ellipse) fell back to the
# built-in sweep, whose records do not match their end points ("Curve[19]
# geometry is inconsistent"). A lateral face is now the exact tensor patch
# of the edge swept by the translation or the rotation (Gmsh: "Unknown" and
# "Surface of Revolution"), and an edge two inputs share sweeps one face.
# `occ.copy` numbers the copy's sub-shapes in the order of the originals'
# tags (BRepBuilderAPI_Copy keeps the structure), and the OCC import no
# longer overwrites the curve loop of a built-in extruded surface, which had
# kept four mirrored points of symmetry_test_occ_builtin.geo at their old
# place. Counts, output lists, volumes and point coordinates are Gmsh's.
module ApiOccSweepCurvedTest
using Test, Tessella

@testset "occ.extrude / revolve of conic-bounded faces (Gmsh 4.15.2)" begin
    A=Tessella.API; occ=A.model.occ
    A.initialize()
    try
        p=[occ.add_point(c...) for c in
           [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(2,0,0),(2,1,0)]]
        l=[occ.add_line(p[1],p[2]),occ.add_line(p[2],p[3]),
           occ.add_line(p[3],p[4]),occ.add_line(p[4],p[1]),
           occ.add_line(p[2],p[5]),occ.add_line(p[5],p[6]),
           occ.add_line(p[6],p[3])]
        s1=occ.add_plane_surface([occ.add_curve_loop([l[1],l[2],l[3],l[4]])])
        s2=occ.add_plane_surface([occ.add_curve_loop([l[5],l[6],l[7],-l[2]])])
        e=occ.add_disk(4,0.5,0,0.5,0.3)
        out=occ.extrude([(2,s1),(2,s2),(2,e)],0,0,0.3)
        occ.synchronize()
        @test out==[(2,8),(3,1),(2,4),(2,5),(2,6),(2,7),(2,12),(3,2),(2,9),
                    (2,10),(2,11),(2,5),(2,14),(3,3),(2,13)]
        @test [length(A.model.get_entities(d)) for d in 0:3]==[14,23,14,3]
        @test occ.get_mass(3,1)≈0.3 rtol=1e-12
        @test occ.get_mass(3,2)≈0.3 rtol=1e-12
        @test occ.get_mass(3,3)≈0.141371669411541 rtol=1e-12
        @test A.model.get_type(2,13)=="Unknown"
        @test count(t->A.model.get_type(2,t)=="Plane",1:14)==13
    finally
        A.finalize()
    end
    A.initialize()
    try
        e=occ.add_disk(1,0,0,0.5,0.3;zAxis=[0,1,0],xAxis=[1,0,0])
        out=occ.revolve([(2,e)],0,0,0,0,0,1,π/2)
        occ.synchronize()
        @test out==[(2,3),(3,1),(2,2)]
        @test [length(A.model.get_entities(d)) for d in 0:3]==[2,3,3,1]
        @test occ.get_mass(3,1)≈π*0.5*0.3*π/2 rtol=1e-12
        @test A.model.get_type(2,2)=="Surface of Revolution"
        # (the lateral area, 16 s interpreted, is checked compiled: Tessella
        # 4.009771604872344, (π/2)·perimeter 4.009771604872327, Gmsh
        # 4.00977967464682)
    finally
        A.finalize()
    end
    A.initialize()
    try
        occ.add_box(0.1,2,0,1,1,1); c=occ.copy([(3,1)])
        occ.mirror(c,1,0,0,0); occ.synchronize()
        @test [A.model.get_value(0,t,Float64[]) for t in 9:16]≈
            [[-0.1,2,1],[-0.1,2,0],[-0.1,3,1],[-0.1,3,0],
             [-1.1,2,1],[-1.1,2,0],[-1.1,3,1],[-1.1,3,0]]
    finally
        A.finalize()
    end
    path=tempname()*".geo"
    write(path,"""
        Point(1) = {0.1,0,0}; Point(2) = {1.1,0,0};
        Point(3) = {1.1,1,0}; Point(4) = {0.1,1,0};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Line Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
        e() = Extrude {0,0,1} { Surface{1}; };
        SetFactory("OpenCASCADE");
        Box(100) = {0.1, 2, 0, 1, 1, 1};
        Symmetry {1, 0, 0, 0} { Duplicata { Volume{100}; } }
        """)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        @test [length(A.model.get_entities(d)) for d in 0:3]==[24,36,18,3]
        @test [A.model.get_value(0,t,Float64[]) for t in 23:30]≈
            [[-0.1,2,1],[-0.1,2,0],[-0.1,3,1],[-0.1,3,0],
             [-1.1,2,1],[-1.1,2,0],[-1.1,3,1],[-1.1,3,0]]
    finally
        A.finalize()
        rm(path;force=true)
    end
end
end
