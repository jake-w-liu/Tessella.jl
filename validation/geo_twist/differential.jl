#!/usr/bin/env julia
# Entity-level differential for .geo twist (`TRANSLATE_ROTATE`) `Extrude`.
# Each script is executed by Tessella's bounded GeoExec interpreter and
# independently opened by Gmsh 4.15.2's built-in parser; the resulting
# entity tags per dimension, every point's coordinates, curve endpoint
# wiring, surface-boundary and volume-shell absolute tags, the `out[]`
# result lists, and spline-curve evaluations at the generatrix knot
# parameters are compared. Gmsh's twist semantics — decoded from
# `GEO_Internals::twist`/`ExtrudePoint` — are `{{delta}, {axis}, {point},
# angle}`: the first vector is the translation, the second the rotation
# axis direction, the third a point on the axis (not the revolve order).
# The generatrix is a `MSH_SEGM_SPLN` through `Geometry.ExtrudeSplinePoints`
# (default 5) vertices, each one `DuplicateVertex`-copied and stepped by
# `angle/d` about the axis and `delta/d` along it.
#
# Point coordinates and curve evaluations are bit-exact for cases whose
# motion is a pure axis-aligned screw. Cases mixing arbitrary-axis
# rotations contract different FMA chains than Gmsh's C++ and are compared
# at a 64-ulp absolute tolerance, the same convention geo_transforms uses.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella

const TARGET_GMSH_VERSION = "4.15.2"

function find_gmsh_executable()
    explicit = get(ENV, "GMSH_EXECUTABLE", "")
    if !isempty(explicit)
        isfile(explicit) || error(
            "GMSH_EXECUTABLE does not name a file: $explicit")
        return realpath(explicit)
    end
    executable = Sys.which("gmsh")
    executable !== nothing && return realpath(executable)
    fallback = "/opt/homebrew/bin/gmsh"
    isfile(fallback) && return realpath(fallback)
    error("Gmsh $TARGET_GMSH_VERSION is required; install it or set GMSH_EXECUTABLE")
end

function find_gmsh_api(executable)
    explicit = get(ENV, "GMSH_JULIA_API", "")
    if !isempty(explicit)
        isfile(explicit) || error(
            "GMSH_JULIA_API does not name a file: $explicit")
        return realpath(explicit)
    end
    prefix = dirname(dirname(executable))
    candidates = (joinpath(prefix, "lib", "gmsh.jl"),
                  joinpath(prefix, "lib64", "gmsh.jl"),
                  joinpath(prefix, "share", "gmsh", "api", "julia", "gmsh.jl"))
    for candidate in candidates
        isfile(candidate) && return realpath(candidate)
    end
    error("could not locate gmsh.jl for $executable; set GMSH_JULIA_API")
end

const SQUARE = """
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

const CASES = (
    # Point twist: the translation is parallel to the axis — a clean helix.
    # {{T}, {A}, {X}} = translate (0,0,1), rotate about z through (0,0,1).
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Point{1}; };
    """,
    # Screw with translation NOT parallel to the axis and a center off the
    # origin: each step's translation is rotated into the current frame.
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {{1,0,0},{5,0,0},{0,2,0},Pi/2} { Point{1}; };
    """,
    # Axis parallel to translation, center off the axis.
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {{0,0,2},{0,0,1},{0,0,0},Pi/3} { Point{1}; };
    """,
    # ExtrudeSplinePoints = 1: a single-step spline (two control points).
    """
    Geometry.ExtrudeSplinePoints = 1;
    Point(1) = {1,0,0,1};
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Point{1}; };
    """,
    # ExtrudeSplinePoints = 3.
    """
    Geometry.ExtrudeSplinePoints = 3;
    Point(1) = {1,0,0,1};
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Point{1}; };
    """,
    # Curve twist: spline generatrices from both endpoints, ruled lateral.
    """
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Line(1) = {1,2};
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Curve{1}; };
    """,
    # Negative generatrix tag.
    """
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Line(1) = {1,2};
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Curve{-1}; };
    """,
    # Surface twist: full-turn-free volume with spline generatrices.
    SQUARE * """
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Surface{1}; };
    """,
    # Signed surface generatrix.
    SQUARE * """
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/3} { Surface{-1}; };
    """,
    # ExtrudeReturnLateralEntities = 0 on a twisted curve.
    """
    Geometry.ExtrudeReturnLateralEntities = 0;
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Line(1) = {1,2};
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Curve{1}; };
    """,
    # A second twist of the top surface closes a toroidal segment.
    SQUARE * """
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Surface{1}; };
    out2[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Surface{-1}; };
    """,
    # Mixed entity list: curve and point through the same twist.
    """
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Line(1) = {1,2};
    Point(5) = {5,0,0,1};
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Curve{1}; Point{5}; };
    """,
    # Non-axis-aligned twist: Gram-Schmidt rotation frame (ULP case).
    """
    Point(1) = {1,2,3,1};
    out[] = Extrude {{0.5,-0.5,1},{1,1,1},{-1,0.5,2},Pi/3} { Point{1}; };
    """,
    # Extrude parameters are mesh metadata — geometry must not change.
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2}
            { Point{1}; Layers{2}; };
    """,
    # Twist after a translate: allocator interleaving across forms.
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {0,0,1} { Point{1}; };
    out2[] = Extrude {{0,0,1},{0,0,1},{0,0,1},Pi/2} { Point{2}; };
    """,
)

# Cases whose Gram-Schmidt rotation mixes axis components: the Julia and
# C++ builds contract different FMAs, so coordinates compare at 64 ulp.
const ULP_CASES = Set((2, 13))
const ULP_TOL = 64 * eps(Float64)

const GMSH_EXECUTABLE = find_gmsh_executable()
const GMSH_CLI_VERSION = strip(read(`$GMSH_EXECUTABLE --version`, String))
(GMSH_CLI_VERSION == TARGET_GMSH_VERSION ||
 startswith(GMSH_CLI_VERSION, TARGET_GMSH_VERSION * "-")) ||
    error("expected Gmsh $TARGET_GMSH_VERSION, got $GMSH_CLI_VERSION")
const GMSH_API_FILE = find_gmsh_api(GMSH_EXECUTABLE)
include(GMSH_API_FILE)
gmsh.GMSH_API_VERSION == TARGET_GMSH_VERSION || error(
    "expected Gmsh API $TARGET_GMSH_VERSION, got $(gmsh.GMSH_API_VERSION)")

const LIST_VARS = ("out", "out2", "x")

_samples_close(a, b) = all(abs.(a .- b) .<= ULP_TOL)

function _points_equal(actual, expected, ulp)
    ulp && return _samples_close(collect(actual), collect(expected))
    return all(reinterpret(UInt64, collect(actual)) .==
               reinterpret(UInt64, collect(expected)))
end

samples = Ref(0)
gmsh.initialize([GMSH_EXECUTABLE, "-nopopup"], false, false)
try
    gmsh.option.setNumber("General.Terminal", 0)
    runtime_version = gmsh.option.getString("General.Version")
    (runtime_version == TARGET_GMSH_VERSION ||
     startswith(runtime_version, TARGET_GMSH_VERSION * "-")) || error(
        "expected Gmsh runtime $TARGET_GMSH_VERSION, got $runtime_version")

    mktempdir() do directory
        for (case_index, source) in enumerate(CASES)
            path = joinpath(directory, "case$case_index.geo")
            write(path, source)

            parsed = Tessella.GeoExec.execute_geo(path)
            model = parsed.model

            gmsh.clear()
            # `gmsh.clear` does not restore Geometry.* option state, so a
            # `Geometry.ExtrudeReturnLateralEntities = 0` or
            # `Geometry.ExtrudeSplinePoints` case would leak into later
            # ones; reset both explicitly before each open.
            gmsh.option.setNumber(
                "Geometry.ExtrudeReturnLateralEntities", 1)
            gmsh.option.setNumber("Geometry.ExtrudeSplinePoints", 5)
            gmsh.open(path)

            ulp = case_index in ULP_CASES

            for dim in 0:3
                expected = sort!(last.(gmsh.model.getEntities(dim)))
                actual = sort!(collect(keys(
                    dim == 0 ? model.points : dim == 1 ? model.curves :
                    dim == 2 ? model.surfaces : model.volumes)))
                actual == expected || error(
                    "case $case_index dim $dim entities differ: " *
                    "Tessella=$actual Gmsh=$expected")
            end
            for tag in sort!(collect(keys(model.points)))
                expected = gmsh.model.getValue(0, tag, Float64[])
                actual = collect(model.points[tag])
                _points_equal(actual, expected, ulp) || error(
                    "case $case_index Point($tag) differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
            for tag in sort!(collect(keys(model.curves)))
                expected = sort!(last.(gmsh.model.getBoundary(
                    [(1, tag)], false, false, false)))
                actual = sort!(collect(model.curves[tag]))
                actual == expected || error(
                    "case $case_index Curve($tag) endpoints differ: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
                # Spline generatrices: evaluate both curves at the knot
                # parameters — Catmull-Rom interpolates the control
                # points, so this covers CP order and count.
                cps = get(model.curve_control_points, tag, nothing)
                cps === nothing && continue
                get(model.curve_types, tag, :line) === :spline || continue
                g = Tessella.Model._spline_geometry(
                    model, tag, "geo_twist differential")
                for i in 0:length(cps)-1
                    u = i / (length(cps) - 1)
                    u = clamp(u, 0.0, 1.0)
                    expected_p = gmsh.model.getValue(1, tag, [u])
                    actual_p = collect(
                        Tessella.Model._spline_point(g, u))
                    _points_equal(actual_p, expected_p, ulp) || error(
                        "case $case_index Curve($tag) at u=$u differs: " *
                        "Tessella=$actual_p Gmsh=$expected_p")
                    samples[] += 1
                end
            end
            for tag in sort!(collect(keys(model.surfaces)))
                expected = sort!(abs.(last.(gmsh.model.getBoundary(
                    [(2, tag)], false, false, false))))
                actual = sort!(unique!(abs.(vcat(
                    (model.loops[l] for l in model.surfaces[tag])...))))
                actual == expected || error(
                    "case $case_index Surface($tag) boundary differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
            for tag in sort!(collect(keys(model.volumes)))
                expected = sort!(abs.(last.(gmsh.model.getBoundary(
                    [(3, tag)], false, false, false))))
                actual = sort!(unique!(abs.(vcat(
                    (model.surface_loops[sl]
                     for sl in model.volumes[tag])...))))
                actual == expected || error(
                    "case $case_index Volume($tag) shell differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
            # Extrusion result lists compare bit-for-bit through the
            # .geo parser variables.
            for name in LIST_VARS
                haskey(parsed.lists, name) || continue
                expected = gmsh.parser.getNumber(name)
                actual = parsed.lists[name]
                actual == expected || println(
                    "case $case_index `$name`: Tessella=$actual Gmsh=$expected")
                all(reinterpret(UInt64, actual) .==
                    reinterpret(UInt64, expected)) || error(
                    "case $case_index list `$name` differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
        end
    end

    println("GEO_TWIST_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(CASES)) samples=$(samples[]) " *
            "bit_exact=$(length(CASES) - length(ULP_CASES))/" *
            "$(length(CASES))")
finally
    gmsh.finalize()
end
