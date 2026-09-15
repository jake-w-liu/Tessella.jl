#!/usr/bin/env julia
# Entity-level differential for .geo translational `Extrude`. Each script is
# executed by Tessella's bounded GeoExec interpreter and independently opened
# by Gmsh 4.15.2's built-in parser; the resulting entity tags per dimension,
# every point's coordinates, curve endpoint wiring, surface-boundary and
# volume-shell absolute tags, and the `out[]`/`out2[]` result lists must
# match bit-for-bit (loop and surface-loop records have no Gmsh entity
# counterpart, so topology is compared through sorted absolute boundary
# tags; `Boundary`-API signs are not part of Gmsh's public surface/volume
# boundary query, so the comparison is orientation-insensitive there).

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
                  "/opt/homebrew/lib/gmsh.jl",
                  "/opt/homebrew/opt/gmsh/lib/gmsh.jl",
                  "/usr/local/opt/gmsh/lib/gmsh.jl")
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
    # Point extrusion: out=[top point, connecting curve].
    """
    Point(1) = {0,0,0,1};
    out[] = Extrude {0,0,1} { Point{1}; };
    """,
    # Curve extrusion: out=[top, lateral surface, laterals...].
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    out[] = Extrude {0,0,1} { Curve{1}; };
    """,
    # Negative generatrix: reversed record, signed lateral tail.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    out[] = Extrude {0,0,1} { Curve{-1}; };
    """,
    # Surface extrusion: Gmsh's lateral/top/volume tag layout.
    SQUARE * """
    out[] = Extrude {0,0,1} { Surface{1}; };
    """,
    # Surface{-1}: sign is metadata-only; source tag kept in the tail.
    SQUARE * """
    out[] = Extrude {0,0,1} { Surface{-1}; };
    """,
    # A second extrusion of the reversed curve: the lateral merges.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    out[] = Extrude {0,0,1} { Curve{1}; };
    out2[] = Extrude {0,0,1} { Curve{-1}; };
    """,
    # Multi-entity list with mixed allocation: out=[2,5,4,-3,8,6].
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    Point(5) = {5,0,0,1};
    out[] = Extrude {0,0,1} { Curve{1}; Point{5}; };
    """,
    # ExtrudeReturnLateralEntities = 0 on every dimension.
    """
    Geometry.ExtrudeReturnLateralEntities = 0;
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    out[] = Extrude {0,0,1} { Curve{1}; };
    out2[] = Extrude {0,0,1} { Curve{-1}; };
    """,
    # ExtrudeReturnLateralEntities = 0 on a surface.
    """
    Geometry.ExtrudeReturnLateralEntities = 0;
    """ * SQUARE * """
    out[] = Extrude {0,0,1} { Surface{1}; };
    """,
    # A second surface extrusion merges the top and laterals into the
    # first volume's boundary; volumes themselves never merge.
    SQUARE * """
    out[] = Extrude {0,0,1} { Surface{1}; };
    out2[] = Extrude {0,0,1} { Surface{-1}; };
    """,
    # Extrude parameters are mesh metadata — geometry must not change.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    out[] = Extrude {0,0,1} { Curve{1}; Layers{3}; Recombine; };
    """,
    # Non-axis-aligned delta and a degenerate (zero) point extrusion.
    """
    Point(1) = {0,0,0,1};
    out[] = Extrude {1.5,-2,0.25} { Point{1}; };
    out2[] = Extrude {0,0,0} { Point{1}; };
    """,
    # Scalar assignment and the GeoEntity selector.
    """
    Point(1) = {0,0,0,1};
    x = Extrude {0,0,1} { GeoEntity{0}{1}; };
    """,
    # Rotational point extrusion: a Circle arc [start, axis-center, end].
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {{0,0,1},{0,0,0},Pi/2} { Point{1}; };
    """,
    # The rotation axis is a line — an off-origin point moves it.
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {{0,0,1},{1,1,0},Pi/2} { Point{1}; };
    """,
    # Curve revolve: shared axis-center point, out=[2,5,4,-3].
    """
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Line(1) = {1,2};
    out[] = Extrude {{0,0,1},{0,0,0},Pi/2} { Curve{1}; };
    """,
    # Negative generatrix and negative angle.
    """
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Line(1) = {1,2};
    out[] = Extrude {{0,0,1},{0,0,0},-Pi/3} { Curve{-1}; };
    """,
    # An endpoint on the axis collapses: TRIC lateral, out=[3,5,4].
    """
    Point(1) = {2,0,0,1};
    Point(2) = {1,0,0,1};
    Point(3) = {0,0,0,1};
    Circle(1) = {1,2,3};
    Line(2) = {3,1};
    Curve Loop(1) = {1,2};
    Plane Surface(1) = {1};
    out[] = Extrude {{0,0,1},{0,0,0},Pi/2} { Curve{2}; };
    """,
    # Arc generatrix: the circle's control points rotate with the copy.
    """
    Point(1) = {2,0,0,1};
    Point(2) = {1,0,0,1};
    Point(3) = {0,0,0,1};
    Circle(1) = {1,2,3};
    out[] = Extrude {{0,0,1},{0,0,0},Pi/2} { Curve{1}; };
    """,
    # Surface revolve: volume shell, laterals 13/17/21/25, top 26.
    """
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Point(3) = {1,1,1,1};
    Point(4) = {1,0,1,1};
    Line(1) = {1,2};
    Line(2) = {2,3};
    Line(3) = {3,4};
    Line(4) = {4,1};
    Curve Loop(1) = {1,2,3,4};
    Plane Surface(1) = {1};
    out[] = Extrude {{0,0,1},{0,0,0},Pi/2} { Surface{1}; };
    """,
    # Signed surface generatrix.
    """
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Point(3) = {1,1,1,1};
    Point(4) = {1,0,1,1};
    Line(1) = {1,2};
    Line(2) = {2,3};
    Line(3) = {3,4};
    Line(4) = {4,1};
    Curve Loop(1) = {1,2,3,4};
    Plane Surface(1) = {1};
    out[] = Extrude {{0,0,1},{0,0,0},Pi/3} { Surface{-1}; };
    """,
    # Collapsed revolves: on-axis point, zero angle, full turn.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    out[] = Extrude {{0,0,1},{0,0,0},Pi/2} { Point{1}; };
    out2[] = Extrude {{0,0,1},{0,0,0},0} { Point{2}; };
    x[] = Extrude {{0,0,1},{0,0,0},2*Pi} { Point{2}; };
    """,
    # An arc spanning more than Pi stays a single Circle entity.
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {{0,0,1},{0,0,0},3*Pi/2} { Point{1}; };
    """,
    # ExtrudeReturnLateralEntities = 0 on a revolved curve.
    """
    Geometry.ExtrudeReturnLateralEntities = 0;
    Point(1) = {1,0,0,1};
    Point(2) = {1,1,0,1};
    Line(1) = {1,2};
    out[] = Extrude {{0,0,1},{0,0,0},Pi/2} { Curve{1}; };
    """,
    # Non-axis-aligned revolve: Gmsh's Gram-Schmidt rotation frame.
    """
    Point(1) = {1,2,3,1};
    Point(2) = {0,1,0,1};
    Line(1) = {1,2};
    out[] = Extrude {{1,1,1},{-1,0.5,2},Pi/3} { Curve{1}; };
    """,
    # A revolve after a translate: allocator interleaving across forms.
    """
    Point(1) = {1,0,0,1};
    out[] = Extrude {0,0,1} { Point{1}; };
    out2[] = Extrude {{0,0,1},{0,0,0},Pi/2} { Point{2}; };
    """,
)

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
            # `Geometry.ExtrudeReturnLateralEntities = 0` case would leak
            # into later ones; reset it explicitly before each open.
            gmsh.option.setNumber(
                "Geometry.ExtrudeReturnLateralEntities", 1)
            gmsh.open(path)

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
                all(reinterpret(UInt64, actual) .==
                    reinterpret(UInt64, expected)) || error(
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

    println("GEO_EXTRUDE_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(CASES)) samples=$(samples[])")
finally
    gmsh.finalize()
end
