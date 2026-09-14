using Test
using Tessella

function _execute_transform_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _transform_error(source::AbstractString)
    try
        _execute_transform_source(source)
        return nothing
    catch err
        return err
    end
end

const _SQUARE_GEO = """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Point(3) = {1,1,0,1};
    Point(4) = {0,1,0,1};
    Line(1) = {1,2};
    Line(2) = {2,3};
    Line(3) = {3,4};
    Line(4) = {4,1};
    Curve Loop(1) = {1,2,3,4};
    Plane Surface(1) = {1};
    """

@testset ".geo Translate/Dilate/Rotate/Symmetry on explicit entities" begin
    # Points move in place.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Translate {5,0,0} { Point{1}; }
        Symmetry {1,0,0,-5} { Point{2}; }
        Dilate {{2,0,0}, 3} { Point{2}; }
        """)
    @test r.model.points[1]==(5.0,0.0,0.0)
    # (1,0,0) reflects across x=5 to (9,0,0), then dilates about x=2 by 3 → 23.
    @test r.model.points[2]==(23.0,0.0,0.0)

    # Curves transform their endpoints; shared vertices move exactly once.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {2,0,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Translate {10,0,0} { Line{1}; Line{2}; Point{1}; }
        """)
    @test r.model.points[1]==(10.0,0.0,0.0)
    @test r.model.points[2]==(11.0,0.0,0.0)   # shared endpoint moved once, not twice
    @test r.model.points[3]==(12.0,0.0,0.0)

    # Surfaces recurse through boundary curves; points stay shared.
    r=_execute_transform_source(_SQUARE_GEO * """
        Translate {5,0,0} { Surface{1}; }
        """)
    @test [r.model.points[t] for t in 1:4]==
        [(5.0,0.0,0.0),(6.0,0.0,0.0),(6.0,1.0,0.0),(5.0,1.0,0.0)]
    @test r.model.surfaces[1]==[1]
    @test r.model.loops[1]==[1,2,3,4]

    # Rotate by an arbitrary (non-π/2) angle about an arbitrary axis.
    r=_execute_transform_source("""
        Point(1) = {1,0,0,1};
        Rotate {{0,0,1},{0,0,0}, Pi/4} { Point{1}; }
        """)
    @test collect(r.model.points[1])≈[sqrt(2)/2,sqrt(2)/2,0.0]

    # Rotate about a non-origin point and a non-axis-aligned axis.
    r=_execute_transform_source("""
        Point(1) = {2,0,0,1};
        Rotate {{0,0,1},{1,0,0}, Pi/2} { Point{1}; }
        """)
    @test collect(r.model.points[1])≈[1.0,1.0,0.0]

    # Per-axis Dilate.
    r=_execute_transform_source("""
        Point(1) = {2,3,4,1};
        Dilate {{1,1,1},{2,4,8}} { Point{1}; }
        """)
    @test r.model.points[1]==(1+2*1.0,1+4*2.0,1+8*3.0)
end

@testset ".geo transform lists: Physical, Boundary, nested transforms" begin
    # Physical groups expand to their members.
    r=_execute_transform_source(_SQUARE_GEO * """
        Physical Point(7) = {1};
        Physical Curve(8) = {1};
        Translate {5,0,0} { Physical Point{7}; Physical Curve{8}; }
        """)
    @test r.model.points[1]==(5.0,0.0,0.0)
    @test r.model.points[2]==(6.0,0.0,0.0)   # via Physical Curve{8}
    @test r.model.points[3]==(1.0,1.0,0.0)

    # Boundary of a surface selects its curves → their endpoints move.
    r=_execute_transform_source(_SQUARE_GEO * """
        Translate {0,9,0} { Boundary{Surface{1};} }
        """)
    @test r.model.points[3]==(1.0,10.0,0.0)
    @test r.model.curves[1]==(1,2)

    # CombinedBoundary cancels shared interior curves.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1}; Point(2) = {1,0,0,1}; Point(3) = {2,0,0,1};
        Point(4) = {0,1,0,1}; Point(5) = {1,1,0,1}; Point(6) = {2,1,0,1};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {4,5}; Line(4) = {5,6};
        Line(5) = {1,4}; Line(6) = {2,5}; Line(7) = {3,6};
        Curve Loop(1) = {1,6,-3,-5};
        Curve Loop(2) = {2,7,-4,-6};
        Plane Surface(1) = {1}; Plane Surface(2) = {2};
        Translate {0,0,9} { CombinedBoundary{Surface{1}; Surface{2};} }
        """)
    # Shared curve 6 cancels; boundary curves move their endpoints.
    @test r.model.points[1]==(0.0,0.0,9.0)
    @test r.model.points[2]==(1.0,0.0,9.0)   # also endpoint of cancelled 6? no — 6's endpoints are 2,5 which are endpoints of surviving 1,3
    @test r.model.points[3]==(2.0,0.0,9.0)
    @test r.model.points[4]==(0.0,1.0,9.0)
    @test r.model.points[5]==(1.0,1.0,9.0)
    @test r.model.points[6]==(2.0,1.0,9.0)

    # Nested transform executes first; its entities are then transformed again.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Translate {10,0,0} { Translate {0,1,0} { Point{1}; } }
        """)
    @test r.model.points[1]==(10.0,1.0,0.0)

    # PointsOf inside a transform list.
    r=_execute_transform_source(_SQUARE_GEO * """
        Translate {0,0,7} { PointsOf{Curve{1};} }
        """)
    @test r.model.points[1]==(0.0,0.0,7.0)
    @test r.model.points[2]==(1.0,0.0,7.0)
    @test r.model.points[3]==(1.0,1.0,0.0)
end

@testset ".geo Duplicata" begin
    # Bare Duplicata copies without merging (Gmsh copy() has no coherence
    # pass). Gmsh allocates the curve tag from the shared NEWREG counter, then
    # two coincident control-point copies (stored on the copy), then the wired
    # beg/end vertex copies.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Duplicata { Line{1}; }
        """)
    @test sort!(collect(keys(r.model.points)))==[1,2,3,4,5,6]
    @test r.model.curves[2]==(5,6)
    @test r.model.curve_control_points[2]==[3,4]
    @test r.model.points[3]==(0.0,0.0,0.0)
    @test r.model.points[5]==(0.0,0.0,0.0)

    # Translate + Duplicata: the copy's whole vertex set (endpoints plus
    # control points) moves, then the post-transform coherence merge collapses
    # coincident pairs onto the lowest tag — exactly Gmsh's pipeline.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Translate {5,0,0} { Duplicata{Line{1};} }
        """)
    @test sort!(collect(keys(r.model.points)))==[1,2,3,4]
    @test r.model.curves[2]==(3,4)
    @test r.model.points[3]==(5.0,0.0,0.0)
    @test r.model.points[4]==(6.0,0.0,0.0)
    @test r.model.points[1]==(0.0,0.0,0.0)  # originals stay put

    # Surface copy: shared allocator tag, deep-copied boundary.
    r=_execute_transform_source(_SQUARE_GEO * """
        Translate {5,0,0} { Duplicata{Surface{1};} }
        """)
    @test sort!(collect(keys(r.model.surfaces)))==[1,5]
    @test r.model.loops[5]==[6,7,8,9]
    @test [r.model.curves[c] for c in 6:9]==
        [(5,6),(6,10),(10,14),(14,5)]
    @test r.model.points[5]==(5.0,0.0,0.0)
    @test r.model.points[14]==(5.0,1.0,0.0)

    # Point copy inside a transform: only the copy moves.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Translate {5,0,0} { Duplicata{Point{1};} }
        """)
    @test r.model.points[1]==(0.0,0.0,0.0)
    @test r.model.points[3]==(5.0,0.0,0.0)

    # Nested Duplicata: copies the inner copy's result.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Duplicata { Duplicata{Point{1};} }
        """)
    @test sort!(collect(keys(r.model.points)))==[1,2,3]
    @test r.model.points[2]==(0.0,0.0,0.0)
    @test r.model.points[3]==(0.0,0.0,0.0)

    # Gmsh grammar: a transform/action is the whole MultipleShape — it cannot
    # be mixed with other list entries, in either order.
    @test _transform_error("""
        Point(1)={0,0,0,1}; Line(1)={1,1};
        Translate {5,0,0} { Duplicata{Point{1};} Line{1}; }
        """) isa ArgumentError
    @test _transform_error("""
        Point(1)={0,0,0,1}; Line(1)={1,1};
        Translate {5,0,0} { Line{1}; Duplicata{Point{1};} }
        """) isa ArgumentError

    # Inline Shape definitions in a transform list create then transform.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Translate {5,0,0} { Point(7) = {0,0,0,1}; }
        """)
    @test r.model.points[7]==(5.0,0.0,0.0)
    @test r.model.points[1]==(0.0,0.0,0.0)
end

@testset ".geo Coherence statements" begin
    # Coherence merges coincident entities globally.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {0,0,0,1};
        Line(1) = {1,2};
        Coherence;
        """)
    @test sort!(collect(keys(r.model.points)))==[1,2]

    # Coherence Point{...} snaps every listed point to the first tag's
    # position, then merges.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Coherence Point{1,2};
        """)
    @test sort!(collect(keys(r.model.points)))==[1]
    @test r.model.points[1]==(0.0,0.0,0.0)

    # Coherence Geometry == Coherence; Coherence Mesh is a no-op mid-script.
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {0,0,0,1};
        Coherence Geometry;
        Coherence Mesh;
        """)
    @test sort!(collect(keys(r.model.points)))==[1]
    @test _transform_error("Point(1)={0,0,0,1}; Coherence Foo;") isa ArgumentError
end

@testset ".geo transforms on volume encodings" begin
    # Materialized Box: a π/2 rotation keeps the encoding; a π/4 rotation
    # drops it while preserving the explicit shell.
    r=_execute_transform_source("""
        Box(1) = {0,0,0, 1,2,1};
        Translate {1,1,1} { Volume{1}; }
        """)
    @test haskey(r.model.box_extents,1)
    @test r.model.box_extents[1]==(1.0,1.0,1.0,1.0,2.0,1.0)
    # Gmsh addBox corner order: point 1 is the (x0,y0,z1) corner.
    @test r.model.points[1]==(1.0,1.0,2.0)
    @test r.model.points[2]==(1.0,1.0,1.0)

    r=_execute_transform_source("""
        Box(1) = {0,0,0, 1,2,1};
        Rotate {{0,0,1},{0,0,0}, Pi/2} { Volume{1}; }
        """)
    @test haskey(r.model.box_extents,1)
    x0,y0,z0,dx,dy,dz=r.model.box_extents[1]
    @test sort!([dx,dy,dz])≈[1.0,1.0,2.0]

    r=_execute_transform_source("""
        Box(1) = {0,0,0, 1,2,1};
        Rotate {{0,0,1},{0,0,0}, Pi/4} { Volume{1}; }
        """)
    @test !haskey(r.model.box_extents,1)   # not axis-aligned → shell only
    @test length(r.model.volumes[1])==1
    @test collect(r.model.points[1])≈[0.0,0.0,1.0]   # (x0,y0,z1) corner
    @test collect(r.model.points[2])≈[0.0,0.0,0.0]
    @test collect(r.model.points[6])≈[sqrt(2)/2,sqrt(2)/2,0.0]  # (1,0,0) corner

    # Sub-entity transform on a materialized box drops the encoding when the
    # shell no longer matches the encoded corners.
    r=_execute_transform_source("""
        Box(1) = {0,0,0, 1,1,1};
        Translate {0,0,1} { Point{1}; }
        """)
    @test !haskey(r.model.box_extents,1)
    @test r.model.points[1]==(0.0,0.0,2.0)  # (x0,y0,z1) corner

    # Implicit primitives update their encodings.
    r=_execute_transform_source("""
        Cylinder(1) = {0,0,0, 0,0,2, 0.5};
        Sphere(2) = {1,1,1, 0.25};
        Cone(3) = {0,0,0, 1,0,0, 0.4,0.2};
        Translate {1,0,0} { Volume{1}; Volume{2}; Volume{3}; }
        Dilate {{0,0,0}, 2} { Volume{2}; }
        Rotate {{1,0,0},{0,0,0}, Pi/2} { Volume{1}; }
        """)
    @test r.model.cylinders[1].center==(1.0,0.0,0.0)
    @test collect(r.model.cylinders[1].axis)≈[0.0,-2.0,0.0]  # Rx(π/2): (0,0,2)→(0,-2,0)
    @test r.model.cylinders[1].radius≈0.5
    @test r.model.cylinders[1].height≈2.0
    @test r.model.spheres[2].center==(4.0,2.0,2.0)  # translated then ×2 dilated
    @test r.model.spheres[2].radius≈0.5
    @test r.model.cones[3].center==(1.0,0.0,0.0)
    @test r.model.cones[3].r1≈0.4

    # Non-representable transforms fail before mutating.
    err=_transform_error("""
        Sphere(1) = {0,0,0,1};
        Dilate {{0,0,0},{2,1,1}} { Volume{1}; }
        """)
    @test err isa ArgumentError
    @test occursin("not representable",err.msg)

    # Boolean operand snapshots follow the Boolean result, not later operand
    # transforms.
    r=_execute_transform_source("""
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.5,0,0, 1,1,1};
        BooleanUnion(3) = { Volume{1}; }{ Volume{2}; };
        Translate {10,0,0} { Volume{3}; }
        Translate {100,0,0} { Volume{1}; }
        """)
    A,B=r.model.boolean_operands[3]
    @test minimum(A.coords[1,:])≈10.0     # snapshot moved with the result
    @test maximum(B.coords[1,:])≈11.5
    @test r.model.box_extents[1][1]≈100.0  # operand transform independent
end

@testset ".geo transform error paths" begin
    @test _transform_error("Point(1)={0,0,0,1}; Translate {1,0,0} { Point{9}; }") isa ArgumentError
    @test occursin("unknown Point",_transform_error(
        "Point(1)={0,0,0,1}; Translate {1,0,0} { Point{9}; }").msg)
    @test occursin("must end with `;`",_transform_error(
        "Point(1)={0,0,0,1}; Translate {1,0,0} { Point{1} }").msg)
    # An empty shape list is legal and transforms nothing (the global
    # coherence merge still runs).
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {0,0,0,1};
        Translate {1,0,0} { }
        """)
    @test sort!(collect(keys(r.model.points)))==[1]
    # A degenerate zero symmetry plane takes Gmsh's p -> 1e-12 floor and acts
    # as the identity (the global merge still runs).
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {0,0,0,1};
        Symmetry {0,0,0,5} { Point{1}; }
        """)
    @test r.model.points[1]==(0.0,0.0,0.0)
    @test sort!(collect(keys(r.model.points)))==[1]
    @test occursin("nonzero",_transform_error(
        "Point(1)={0,0,0,1}; Rotate {{0,0,0},{0,0,0},Pi/2} { Point{1}; }").msg)
    @test occursin("OpenCASCADE",_transform_error(
        "Point(1)={0,0,0,1}; Affine {1,0,0,0, 0,1,0,0, 0,0,1,0} { Point{1}; }").msg)
    @test _transform_error(
        "Point(1)={0,0,0,1}; Translate {1,0,0} { Frobnicate{1}; }") isa ArgumentError
    @test occursin("unknown action on multiple shapes",_transform_error(
        "Point(1)={0,0,0,1}; Translate {1,0,0} { Frobnicate{Point{1};} }").msg)
    @test occursin("unknown Physical",_transform_error(
        "Point(1)={0,0,0,1}; Translate {1,0,0} { Physical Point{99}; }").msg)
    # Boundary of a Point returns nothing (Gmsh behavior, no error).
    r=_execute_transform_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {0,0,0,1};
        Boundary { Point{1}; }
        """)
    @test sort!(collect(keys(r.model.points)))==[1,2]
    # Malformed parameter counts.
    @test _transform_error("Point(1)={0,0,0,1}; Translate {1,0} { Point{1}; }") isa ArgumentError
    @test _transform_error("Point(1)={0,0,0,1}; Symmetry {1,0,0} { Point{1}; }") isa ArgumentError
end
