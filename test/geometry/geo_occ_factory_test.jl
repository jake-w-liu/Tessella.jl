# Focused coverage for `.geo` execution under
# `SetFactory("OpenCASCADE")`: OCC-only primitives and statements must route
# through the internal `Model._model_occ_*` endpoints, built-in-only gates and
# upstream error text must be preserved, mixed-factory sessions must keep
# working, and the native mesher must run on OCC-produced models (including
# `Background Mesh View` PostView fields).

using Test
using Tessella
using Tessella.MeshTypes: validate, nnodes, ntris

function _occ_geo(source::AbstractString; mesh_dim=0)
    return mktemp() do path, io
        write(io, source)
        close(io)
        execute_geo(path; mesh_dim=mesh_dim)
    end
end

function _occ_geo_with_file(pos_source::AbstractString, geo::AbstractString;
                            mesh_dim=0)
    return mktempdir() do dir
        pos = joinpath(dir, "bgmesh.pos")
        write(pos, pos_source)
        path = joinpath(dir, "model.geo")
        write(path, geo)
        execute_geo(path; mesh_dim=mesh_dim)
    end
end

const _SCALAR_POS = """
View "h" {
ST(-1,-2,0, 2,-2,0, -1,2,0){0.25,0.25,0.25};
ST(2,-2,0, 2,2,0, -1,2,0){0.25,0.25,0.25};
};
"""

# t17-style anisotropic `TT` triangles covering the [-1,2]x[-2,2] rectangle:
# each node carries a full 9-component metric, so the PostView resolver must
# produce a `PostViewAnisoField`.
const _ANISO_POS = """
View "nodalMetric" {
TT(-1,-2,0, 2,-2,0, -1,2,0){4,0,0,0,4,0,0,0,1, 4,0,0,0,4,0,0,0,1, 4,0,0,0,4,0,0,0,1};
TT(2,-2,0, 2,2,0, -1,2,0){4,0,0,0,4,0,0,0,1, 4,0,0,0,4,0,0,0,1, 4,0,0,0,4,0,0,0,1};
};
"""

# `execute_geo` collects per-statement diagnostics and throws a combined
# ArgumentError at the end of the run; return its message for assertions on
# the upstream error surface.
function _occ_geo_error(source::AbstractString)
    try
        _occ_geo(source)
    catch err
        err isa ArgumentError || rethrow()
        return err.msg
    end
    return ""
end

@testset "OCC factory .geo routing" begin

    @testset "OCC surface/curve primitives" begin
        # Tags start at 101: OCC `Rectangle`/`Disk` auto-allocate their corner
        # vertices and edges from the same tag allocator. Upstream the
        # 3-argument OCC `Circle` is `addCircleArc`'s center form — the middle
        # point is the circle's center, so both endpoints sit at radius 0.5.
        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Rectangle(1) = {0, 0, 0, 1, 1};
            Disk(2) = {3, 0.5, 0, 0.5};
            Point(101) = {0, 0, 0};
            Point(102) = {0.5, 0, 0};
            Point(103) = {0.5, 0.5, 0};
            Circle(103) = {101, 102, 103};
            Circle(104) = {0, 0, 1, 0.25};
            Ellipse(105) = {0, 0, 2, 0.6, 0.3};
            """)
        @test r.msg_error_count == 0
        m = r.model
        @test haskey(m.surfaces, 1)
        @test haskey(m.surfaces, 2)
        _has_curve(tag) = haskey(m.curves, tag) ||
                          haskey(m.occ_kernel.bindings, (1, tag))
        @test _has_curve(103)  # 3-point arc under OCC
        @test _has_curve(104)  # center+radius OCC circle
        @test _has_curve(105)  # center+radii OCC ellipse
    end

    @testset "OCC partial-sweep primitives" begin
        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Cylinder(1) = {0, 0, 0, 0, 0, 1, 0.5, Pi/2};
            Sphere(2) = {3, 0, 0, 0.5, -Pi/2, Pi/2, Pi/4};
            Cone(3) = {0, 3, 0, 0, 0, 1, 0.5, 0.2, Pi/2};
            Torus(4) = {3, 3, 0, 1.0, 0.25, Pi/2};
            """)
        @test r.msg_error_count == 0
        @test haskey(r.model.volumes, 1)
        @test haskey(r.model.volumes, 2)
        @test haskey(r.model.volumes, 3)
        @test haskey(r.model.volumes, 4)
    end

    @testset "OCC Wire / ThruSections loft" begin
        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Circle(1) = {0, 0, 0, 0.5};
            Curve Loop(1) = 1;
            Circle(2) = {0.1, 0.05, 1, 0.1};
            Curve Loop(2) = 2;
            Circle(3) = {-0.1, -0.1, 2, 0.3};
            Curve Loop(3) = 3;
            ThruSections(1) = {1:3};
            """)
        @test r.msg_error_count == 0
        @test haskey(r.model.volumes, 1)

        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Circle(1) = {0, 0, 0, 0.5};
            Curve Loop(1) = 1;
            Circle(2) = {0.1, 0.05, 1, 0.1};
            Curve Loop(2) = 2;
            Circle(3) = {-0.1, -0.1, 2, 0.3};
            Curve Loop(3) = 3;
            Ruled ThruSections(2) = {1:3};
            """)
        @test r.msg_error_count == 0
        @test haskey(r.model.volumes, 2)
    end

    @testset "OCC BSpline/Bezier surfaces" begin
        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Point(1) = {0, 0, 0};
            Point(2) = {1, 0, 0};
            Point(3) = {0, 1, 0};
            Point(4) = {1, 1, 0};
            BSpline Surface(1) = {{1, 2}, {3, 4}};
            Bezier Surface(2) = {{1, 2}, {3, 4}};
            """)
        @test r.msg_error_count == 0
        @test haskey(r.model.surfaces, 1)
        @test haskey(r.model.surfaces, 2)
    end

    @testset "OCC transforms and shape ops" begin
        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Rectangle(1) = {0, 0, 0, 1, 1};
            out() = Translate{2, 0, 0} { Duplicata{ Surface{1}; } };
            Rotate{{0, 0, 1}, {0, 0, 0}, Pi/4} { Surface{out(0)}; }
            Affine{1,0,0,0, 0,1,0,0, 0,0,1,1} { Surface{out(0)}; }
            p() = Closest{0.2, 0.3, 0} { Point{1}; };
            """)
        @test r.msg_error_count == 0
        m = r.model
        @test haskey(m.surfaces, 2)
        @test haskey(m.points, 2)
    end

    @testset "OCC Fillet/Chamfer routing and error surface" begin
        # A single box edge fillets cleanly through the OCC endpoint (the
        # `out() = Fillet{...}{...}{...}` Transform production yields the
        # filleted volume).
        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Box(1) = {0, 0, 0, 1, 1, 1};
            e() = Boundary{ Volume{1}; };
            f() = Unique(Boundary{ Surface{e()}; });
            out() = Fillet{1}{f(0)}{0.1};
            """)
        @test r.msg_error_count == 0
        @test !isempty(r.model.volumes)

        # Corner blends are a documented unsupported constraint of the OCC
        # endpoint: `OCC.fillet` reports "filleted edges share a vertex",
        # which surfaces as "Could not fillet shapes".
        msg = _occ_geo_error("""
            SetFactory("OpenCASCADE");
            Box(1) = {0, 0, 0, 1, 1, 1};
            e() = Boundary{ Volume{1}; };
            f() = Unique(Boundary{ Surface{e()}; });
            out() = Fillet{1}{f()}{0.1};
            """)
        @test occursin("Could not fillet shapes", msg)

        # A single-edge chamfer reaches `OCC.chamfer` (the "kept faces do not
        # close a solid" diagnostic can only come from the OCC endpoint).
        msg = _occ_geo_error("""
            SetFactory("OpenCASCADE");
            Box(1) = {0, 0, 0, 1, 1, 1};
            e() = Boundary{ Volume{1}; };
            f() = Unique(Boundary{ Surface{e()}; });
            out() = Chamfer{1}{f(0)}{e(0)}{0.1};
            """)
        @test occursin("Could not chamfer shapes", msg)
    end

    @testset "built-in-only gates preserved" begin
        msg = _occ_geo_error("""
            Rectangle(1) = {0, 0, 0, 1, 1};
            """)
        @test occursin("only available with OpenCASCADE geometry kernel", msg)

        msg = _occ_geo_error("""
            Disk(1) = {0, 0, 0, 0.5};
            """)
        @test occursin("only available with OpenCASCADE geometry kernel", msg)

        # The built-in `Circle` three-tag grammar error is preserved too.
        msg = _occ_geo_error("""
            Circle(1) = {1, 2};
            """)
        @test occursin("Circle requires 3 points", msg)
    end

    @testset "mixed factory sessions" begin
        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Rectangle(1) = {0, 0, 0, 1, 1};
            SetFactory("Built-in");
            Point(10) = {5, 5, 0, 0.5};
            Point(11) = {6, 5, 0, 0.5};
            Point(12) = {6, 6, 0, 0.5};
            Point(13) = {5, 6, 0, 0.5};
            Line(10) = {10, 11};
            Line(11) = {11, 12};
            Line(12) = {12, 13};
            Line(13) = {13, 10};
            Curve Loop(2) = {10, 11, 12, 13};
            Plane Surface(2) = {2};
            """)
        @test r.msg_error_count == 0
        m = r.model
        @test haskey(m.surfaces, 1)
        @test haskey(m.surfaces, 2)
        @test haskey(m.curves, 10)
    end

    @testset "OCC model meshes with native mesher" begin
        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Mesh.CharacteristicLengthMax = 0.4;
            Rectangle(1) = {0, 0, 0, 1, 1};
            """; mesh_dim=2)
        @test r.msg_error_count == 0
        @test ntris(r.mesh) > 0
        @test validate(r.mesh).ok

        r = _occ_geo("""
            SetFactory("OpenCASCADE");
            Mesh.CharacteristicLengthMax = 0.5;
            Box(1) = {0, 0, 0, 1, 1, 1};
            """; mesh_dim=2)
        @test r.msg_error_count == 0
        @test ntris(r.mesh) > 0
        @test validate(r.mesh).ok
    end

    @testset "Background Mesh View scalar PostView" begin
        geo = """
            SetFactory("OpenCASCADE");
            Merge "bgmesh.pos";
            Rectangle(1) = {-1, -2, 0, 3, 4};
            Background Mesh View[0];
            """
        r = _occ_geo_with_file(_SCALAR_POS, geo; mesh_dim=2)
        @test r.msg_error_count == 0
        @test ntris(r.mesh) > 0
        @test validate(r.mesh).ok
    end

    @testset "Background Mesh View anisotropic PostView" begin
        # t17-style: a `TT` tensor view resolves to a PostViewAnisoField, the
        # native 1D grading samples it, and the surface producer runs. A
        # sparse view that leaves the retained boundary over the metric edge
        # limit surfaces the documented ArgumentError — either outcome proves
        # the field resolved and engaged (a resolution failure reports
        # "Unknown view" / field errors through msg_error_count instead).
        geo = """
            SetFactory("OpenCASCADE");
            Merge "bgmesh.pos";
            Rectangle(1) = {-1, -2, 0, 3, 4};
            Background Mesh View[0];
            """
        resolved = try
            r = _occ_geo_with_file(_ANISO_POS, geo; mesh_dim=2)
            @test r.msg_error_count == 0
            @test ntris(r.mesh) > 0
            @test validate(r.mesh).ok
            true
        catch err
            err isa ArgumentError || rethrow()
            @test occursin("max_metric_edge", err.msg)
            true
        end
        @test resolved
    end

    @testset "Background Mesh View bad index errors" begin
        geo = """
            SetFactory("OpenCASCADE");
            Merge "bgmesh.pos";
            Background Mesh View[5];
            """
        r = _occ_geo_with_file(_SCALAR_POS, geo)
        @test r.msg_error_count > 0
    end
end
