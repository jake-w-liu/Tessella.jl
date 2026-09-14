#!/usr/bin/env julia
# Entity-level differential for .geo curved entities. Each script is executed
# by Tessella's bounded GeoExec interpreter and independently opened by Gmsh
# 4.15.2's built-in parser; the resulting entity tags per dimension, every
# point's coordinates, curve endpoint wiring, ordered signed surface
# boundaries (which exercise Tessella's SortEdgesInLoop port), `getType`
# strings, and arc evaluations must match bit-for-bit. Tessella rejects
# `EndCurve`-invalid arcs (zero radius, unsolvable ellipse, non-cocircular
# circle) at parse time where Gmsh only logs an error and keeps the entity,
# so every case here is a valid arc by construction.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.Model: model_value, model_entity_type

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

const QUARTER_CIRCLE = """
    Point(1) = {1,0,0};
    Point(2) = {0,0,0};
    Point(3) = {0,1,0};
    Circle(1) = {1,2,3};
    """

const CASES = (
    # Quarter circle in the z=0 plane.
    QUARTER_CIRCLE,
    # Semicircle: collinear endpoints make the point-derived normal
    # degenerate, so the `Plane{..}` normal picks the -y half arc.
    """
    Point(1) = {1,0,0};
    Point(2) = {0,0,0};
    Point(3) = {-1,0,0};
    Circle(1) = {1,2,3} Plane {0,0,-1};
    """,
    # Circle in a non-xy plane; the normal comes only from the points.
    """
    Point(1) = {0,0,1};
    Point(2) = {0,0,0};
    Point(3) = {0,1,0};
    Circle(1) = {1,2,3};
    """,
    # Four-tag ellipse: start=(2,0), major direction +x, end=(0,1)
    # solves to a=2, b=1.
    """
    Point(1) = {2,0,0};
    Point(2) = {0,0,0};
    Point(5) = {3,0,0};
    Point(3) = {0,1,0};
    Ellipse(1) = {1,2,5,3};
    """,
    # Three-tag ellipse: Gmsh duplicates start as the major point, which
    # solves to the unit circle.
    QUARTER_CIRCLE * "Ellipse(2) = {1,2,3};\n",
    # Ruled Surface over an out-of-order 4-curve loop and Surface over a
    # 3-curve loop; both report `getType` "Surface". The loop members are
    # deliberately scrambled so the boundary comparison exercises
    # SortEdgesInLoop chaining.
    QUARTER_CIRCLE * """
    Point(4) = {-1,0,0};
    Point(5) = {0,-1,0};
    Circle(2) = {3,2,4};
    Circle(3) = {4,2,5};
    Line(4) = {5,1};
    Line(5) = {4,1};
    Curve Loop(1) = {3,1,4,2};
    Curve Loop(2) = {5,1,2};
    Ruled Surface(1) = {1};
    Surface(2) = {2} In Sphere {2};
    """,
    # `Using Point` constraint variant and a bare Surface.
    QUARTER_CIRCLE * """
    Point(4) = {-1,0,0};
    Point(9) = {0,0,0.5};
    Circle(2) = {3,2,4};
    Line(3) = {4,1};
    Curve Loop(1) = {1,2,3};
    Surface(1) = {1} Using Point {9};
    """,
    # Extruding a Circle: the top copy keeps the "Circle" type, the
    # lateral is a "Surface" (MSH_SURF_REGL).
    QUARTER_CIRCLE * "out[] = Extrude {0,0,1} { Curve{1}; };\n",
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

# Parameters sampled along every arc in each case.
const SAMPLE_US = (0.0, 0.25, 0.5, 0.75, 1.0)

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
                expected_type = gmsh.model.getType(1, tag)
                actual_type = model_entity_type(model, 1, tag)
                actual_type == expected_type || error(
                    "case $case_index Curve($tag) type differs: " *
                    "Tessella=$actual_type Gmsh=$expected_type")
                expected = sort!(last.(gmsh.model.getBoundary(
                    [(1, tag)], false, false, false)))
                actual = sort!(collect(model.curves[tag]))
                actual == expected || error(
                    "case $case_index Curve($tag) endpoints differ: " *
                    "Tessella=$actual Gmsh=$expected")
                for u in SAMPLE_US
                    expected = gmsh.model.getValue(1, tag, [u])
                    actual = collect(model_value(model, 1, tag, [u]))
                    all(reinterpret(UInt64, actual) .==
                        reinterpret(UInt64, expected)) || error(
                        "case $case_index Curve($tag) eval u=$u differs: " *
                        "Tessella=$actual Gmsh=$expected")
                    samples[] += 1
                end
            end
            for tag in sort!(collect(keys(model.surfaces)))
                expected_type = gmsh.model.getType(2, tag)
                actual_type = model_entity_type(model, 2, tag)
                actual_type == expected_type || error(
                    "case $case_index Surface($tag) type differs: " *
                    "Tessella=$actual_type Gmsh=$expected_type")
                # Ordered generatrices: Gmsh's Curve Loop records keep the
                # declaration order after SortEdgesInLoop chaining.
                # Extrusion-created laterals report unsigned tags from
                # getBoundary, so the comparison is on absolute tags.
                expected = abs.(last.(gmsh.model.getBoundary(
                    [(2, tag)], false, false, false)))
                actual = abs.(vcat(
                    (model.loops[l] for l in model.surfaces[tag])...))
                actual == expected || error(
                    "case $case_index Surface($tag) boundary differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
            if haskey(parsed.lists, "out")
                expected = gmsh.parser.getNumber("out")
                actual = parsed.lists["out"]
                all(reinterpret(UInt64, actual) .==
                    reinterpret(UInt64, expected)) || error(
                    "case $case_index list `out` differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
        end
    end

    println("GEO_CURVED_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(CASES)) samples=$(samples[])")
finally
    gmsh.finalize()
end
