using Test
using Tessella

function _execute_extrude_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _extrude_error(source::AbstractString)
    try
        _execute_extrude_source(source)
        return nothing
    catch err
        return err
    end
end

const _EXTRUDE_SQUARE = """
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

@testset ".geo translational Extrude" begin
    # Point extrusion: top point then connecting curve (Gmsh-verified
    # layout and output ordering).
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        out[] = Extrude {0,0,1} { Point{1}; };
        """)
    @test r.model.points[1]==(0.0,0.0,0.0)
    @test r.model.points[2]==(0.0,0.0,1.0)
    @test r.model.curves[1]==(1,2)
    @test r.lists["out"]==[2.0,1.0]

    # A standalone (non-assignment) Extrude statement executes for its
    # side effects only.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        Extrude {0,0,1} { Point{1}; };
        """)
    @test r.model.points[2]==(0.0,0.0,1.0)
    @test r.model.curves[1]==(1,2)

    # Curve extrusion: top curve, lateral surface, then the lateral
    # generatrices (all but source/top). Gmsh: out=[2,5,4,-3].
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude {0,0,1} { Curve{1}; };
        """)
    @test r.lists["out"]==[2.0,5.0,4.0,-3.0]
    @test r.model.curves[2]==(3,4)       # chapeau copy
    @test r.model.curves[3]==(1,3)       # connecting curve at beg
    @test r.model.curves[4]==(2,4)       # connecting curve at end
    @test r.model.loops[only(r.model.surfaces[5])]==[1,4,-2,-3]

    # Negative generatrix: the reversed record extrudes; the signed input
    # tag is kept in the lateral tail (Gmsh: out=[2,5,-1,4,-3]).
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude {0,0,1} { Curve{-1}; };
        """)
    @test r.lists["out"]==[2.0,5.0,-1.0,4.0,-3.0]
    # The chapeau copies the reversed record: point 3 is the copy of
    # point 2, so (3,4) runs (1,0,1)->(0,0,1) — reversed vs Curve{1}.
    @test r.model.curves[2]==(3,4)
    @test r.model.points[3]==(1.0,0.0,1.0)
    @test r.model.points[4]==(0.0,0.0,1.0)
    @test r.model.loops[only(r.model.surfaces[5])]==[-1,4,-2,-3]

    # Surface extrusion materializes the volume with the exact Gmsh tag
    # layout: laterals 13,17,21,25, re-tagged top 26, dedicated volume 1.
    r=_execute_extrude_source(_EXTRUDE_SQUARE * """
        out[] = Extrude {0,0,1} { Surface{1}; };
        """)
    @test r.lists["out"]==[26.0,1.0,13.0,17.0,21.0,25.0]
    @test [r.model.loops[l] for l in r.model.surfaces[13]]==[[1,12,-6,-11]]
    @test [r.model.loops[l] for l in r.model.surfaces[17]]==[[2,16,-7,-12]]
    @test [r.model.loops[l] for l in r.model.surfaces[21]]==[[3,20,-8,-16]]
    @test [r.model.loops[l] for l in r.model.surfaces[25]]==[[4,11,-9,-20]]
    @test [r.model.loops[l] for l in r.model.surfaces[26]]==[[6,7,8,9]]
    @test [r.model.surface_loops[sl] for sl in r.model.volumes[1]]==
        [[-1,26,13,17,21,25]]
    @test sort!(collect(keys(r.model.points)))==[1,2,3,4,5,6,10,14]

    # Surface{-1}: the sign is metadata-only for the built-in geometry; the
    # source tag is kept in the lateral tail (Gmsh: out=[26,1,1,13,...]).
    r=_execute_extrude_source(_EXTRUDE_SQUARE * """
        out[] = Extrude {0,0,1} { Surface{-1}; };
        """)
    @test r.lists["out"]==[26.0,1.0,1.0,13.0,17.0,21.0,25.0]
    @test [r.model.surface_loops[sl] for sl in r.model.volumes[1]]==
        [[-1,26,13,17,21,25]]

    # Multiple entities per list run sequentially: Gmsh's allocator gives
    # out=[2,5,4,-3,8,6] for {Curve{1}; Point{5}}.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Point(5) = {5,0,0,1};
        out[] = Extrude {0,0,1} { Curve{1}; Point{5}; };
        """)
    @test r.lists["out"]==[2.0,5.0,4.0,-3.0,8.0,6.0]

    # `x = Extrude{..}` (scalar assignment form) stores the flat list.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        x = Extrude {0,0,1} { Point{1}; };
        """)
    @test r.lists["x"]==[2.0,1.0]

    # `GeoEntity{dim}{tags}` is the generic selector inside shape lists.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        out[] = Extrude {0,0,1} { GeoEntity{0}{1}; };
        """)
    @test r.lists["out"]==[2.0,1.0]
    @test r.model.curves[1]==(1,2)

    # Geometry.ExtrudeReturnLateralEntities = 0 keeps only the top/body
    # pair per input (Gmsh: surface -> [26,1], curve -> [2,5],
    # point -> [2,1]).
    r=_execute_extrude_source(_EXTRUDE_SQUARE * """
        Geometry.ExtrudeReturnLateralEntities = 0;
        out[] = Extrude {0,0,1} { Surface{1}; };
        """)
    @test r.lists["out"]==[26.0,1.0]
    # The geometry itself is still fully materialized.
    @test haskey(r.model.volumes,1)
    @test haskey(r.model.surfaces,26)

    r=_execute_extrude_source("""
        Geometry.ExtrudeReturnLateralEntities = 0;
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude {0,0,1} { Curve{1}; };
        """)
    @test r.lists["out"]==[2.0,5.0]

    # Merged lateral: a second extrusion of the reversed curve produces
    # the same lateral surface, which coherence merges — the body is then
    # dropped from the output (Gmsh: out2=[-2]).
    r=_execute_extrude_source("""
        Geometry.ExtrudeReturnLateralEntities = 0;
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        out[] = Extrude {0,0,1} { Curve{1}; };
        out2[] = Extrude {0,0,1} { Curve{-1}; };
        """)
    @test r.lists["out"]==[2.0,5.0]
    @test r.lists["out2"]==[-2.0]
    @test sort!(collect(keys(r.model.surfaces)))==[5]
    @test sort!(collect(keys(r.model.curves)))==[1,2,3,4]
end

@testset ".geo Extrude parameters" begin
    # Layers{counts}{heights}... documented Gmsh form is
    # Layers{{counts...},{heights...}}; Recombine, ScaleLast, QuadTri and
    # Using modifiers attach to every created entity.
    r=_execute_extrude_source(_EXTRUDE_SQUARE * """
        out[] = Extrude {0,0,1} { Surface{1};
            Layers{{2,4},{0.2,0.8}};
            ScaleLast;
            QuadTriAddVerts RecombLaterals;
            Using Index[2];
        };
        """)
    p=r.model.meshing.extrude
    @test haskey(p,(2,26)) && haskey(p,(3,1)) && haskey(p,(2,13))
    q=p[(2,26)]
    @test q.layers==[2,4]
    @test q.heights==[0.2,0.8]
    @test q.scale_last==true
    @test q.recombine==false
    @test q.quad_to_tri==:add_verts
    @test q.recomb_laterals==true

    # Layers{n}: one layer of n elements, height 1.0.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Extrude {0,0,1} { Curve{1}; Layers{3}; Recombine; };
        """)
    q=r.model.meshing.extrude[(1,2)]
    @test q.layers==[3]
    @test q.heights==[1.0]
    @test q.recombine==true

    # Layers{counts,heights}: flat one-layer form (Gmsh splits at the
    # first comma of the unbraced lists).
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Extrude {0,0,1} { Curve{1}; Layers{3,5}; };
        """)
    q=r.model.meshing.extrude[(1,2)]
    @test q.layers==[3]
    @test q.heights==[5.0]

    # Using View[i] is an inert translational modifier like Using Index.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        out[] = Extrude {0,0,1} { Point{1}; Using View[3]; };
        """)
    @test r.lists["out"]==[2.0,1.0]

    # QuadTriNoNewVerts selects the other quad-to-tri mode.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Extrude {0,0,1} { Curve{1}; QuadTriNoNewVerts; };
        """)
    @test r.model.meshing.extrude[(2,5)].quad_to_tri==:no_new_verts

    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude {0,0,1} { Point{1}; Layers{2,3,4}; };
        """)
    @test err isa ArgumentError
    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude {0,0,1} { Point{1}; Layers{{2,3},{0.5}}; };
        """)
    @test err isa ArgumentError
    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude {0,0,1} { Point{1}; Hole(2); };
        """)
    @test err isa ArgumentError
end

@testset ".geo Extrude degenerate and error paths" begin
    # Zero delta on a point: the copy coincides — Gmsh returns the source
    # tag, creates no curve, and leaves the unmerged copy behind.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        out[] = Extrude {0,0,0} { Point{1}; };
        """)
    @test r.lists["out"]==[1.0]
    @test isempty(r.model.curves)
    @test sort!(collect(keys(r.model.points)))==[1,2]

    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude {{0,0,1},{0,0,0},Pi/2} { Point{1}; };
        """)
    @test err isa ArgumentError
    @test occursin("rotational",sprint(showerror,err))

    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude {{0,0,1},{0,0,0},{1,0,0},Pi/2} { Point{1}; };
        """)
    @test err isa ArgumentError
    @test occursin("twist",sprint(showerror,err))

    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude { Point{1}; };
        """)
    @test err isa ArgumentError
    @test occursin("boundary-layer",sprint(showerror,err))

    # The pipe form splits at the shape-list brace; the trailing
    # `Using Wire` must reject the whole statement before side effects.
    err=_extrude_error(_EXTRUDE_SQUARE * """
        Extrude { Surface{1}; } Using Wire {3};
        """)
    @test err isa ArgumentError
    @test occursin("pipe",sprint(showerror,err))

    err=_extrude_error(_EXTRUDE_SQUARE * """
        out[] = Extrude {0,0,1} { Surface{1}; };
        out2[] = Extrude {0,0,1} { Volume{1}; };
        """)
    @test err isa ArgumentError
    @test occursin("Volume",sprint(showerror,err))

    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude {0,0,1} { Point{99}; };
        """)
    @test err isa ArgumentError
    @test occursin("Point[99]",sprint(showerror,err))

    # Every entity block in the shape list must end with a semicolon.
    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude {0,0,1} { Point{1} };
        """)
    @test err isa ArgumentError
    @test occursin("semicolon",sprint(showerror,err))

    err=_extrude_error("""
        Point(1) = {0,0,0,1};
        Extrude {0,0,1} { Extrude {0,0,1} { Point{1}; }; };
        """)
    @test err isa ArgumentError
end

@testset ".geo Extrude allocator and variable plumbing" begin
    # The chapeau point uses the dedicated point counter while every
    # curve/surface tag goes through the shared NEWREG maximum — matching
    # Gmsh's mixed allocation in the multi-entity probe above.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(1) = {1,2};
        Point(5) = {5,0,0,1};
        Extrude {0,0,1} { Curve{1}; };
        out[] = Extrude {0,0,1} { Point{5}; };
        """)
    @test r.lists["out"]==[8.0,6.0]
    @test r.model.points[8]==(5.0,0.0,1.0)
    @test r.model.curves[6]==(5,8)

    # Extrude output lists feed later list expressions like any other
    # variable.
    r=_execute_extrude_source("""
        Point(1) = {0,0,0,1};
        out[] = Extrude {0,0,1} { Point{1}; };
        n = #out[];
        """)
    @test r.lists["out"]==[2.0,1.0]

    # read_geo_params degrades an unexecutable Extrude term to an
    # unavailable list instead of throwing.
    gp=mktemp() do path,io
        write(io,"""
            Point(1) = {0,0,0,1};
            out[] = Extrude {0,0,1} { Point{1}; };
            Mesh.MeshSizeMin = 0.25;
            """)
        close(io)
        Tessella.IO.read_geo_params(path)
    end
    @test gp.mesh_size_min==0.25
end
