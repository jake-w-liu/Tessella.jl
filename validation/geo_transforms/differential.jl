#!/usr/bin/env julia
# Entity-level differential for .geo explicit-entity transforms and shape
# actions. Each script is executed by Tessella's bounded GeoExec interpreter
# and independently opened by Gmsh 4.15.2's built-in parser; the resulting
# entity tags per dimension, every point's coordinates, physical-group
# memberships, curve endpoint wiring, and surface-loop topology must match
# bit-for-bit (loop and surface-loop records have no Gmsh entity
# counterpart, so their tags differ deliberately and are compared
# structurally through the surfaces that reference them).
#
# Cases marked `exact = false` contain arbitrary-angle `Rotate`s: Tessella
# reproduces Gmsh's Gram-Schmidt matrix and three-pass pipeline, but the
# compiled C++ contracts some multiply-adds into FMAs in a build-specific
# pattern, so those coordinates are compared at a 64-ulp absolute
# tolerance instead of bit-for-bit.

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

const CASES = (
    # Translate on an explicit point; shared-vertex in-place semantics.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    Translate {5,0,0} { Point{1}; }
    """,
    # A curve transform moves its endpoints; a shared vertex moves once.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Point(3) = {1,1,0,1};
    Line(1) = {1,2};
    Line(2) = {2,3};
    Translate {0,5,0} { Curve{1}; Curve{2}; }
    """,
    # A surface transform recurses through boundary curves.
    """
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
    Dilate {{0,0,0},{2,1,1}} { Surface{1}; }
    """,
    # Arbitrary-axis rotation and plane symmetry.
    """
    Point(1) = {1,0,0,1};
    Point(2) = {0,1,0,1};
    Point(3) = {0,0,1,1};
    Rotate {{1,1,1},{0,0,0}, 2*Pi/3} { Point{1}; Point{2}; Point{3}; }
    Point(4) = {1,0,0,1};
    Symmetry {1,0,0,-0.5} { Point{4}; }
    """,
    # Physical selectors inside a transform list.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    Physical Point(7) = {1};
    Physical Curve(8) = {1};
    Translate {0,9,0} { Physical Point{7}; Physical Curve{8}; }
    """,
    # Boundary action as the whole shape list.
    """
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
    Translate {0,0,7} { Boundary{Surface{1};} }
    """,
    # Nested transforms compose; PointsOf selects boundary points.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    Translate {1,0,0} { Translate{0,2,0} { Point{1}; } }
    Dilate {{0,0,0}, 3} { PointsOf{Curve{1};} }
    """,
    # Duplicata of a point inside a transform: only the copy moves.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    Translate {5,0,0} { Duplicata{Point{1};} }
    """,
    # Duplicata of a curve inside a transform: the copy's whole vertex set
    # (endpoint copies plus control-point copies) moves, then the global
    # post-transform coherence merge collapses coincident pairs onto the
    # lowest tag.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    Translate {5,0,0} { Duplicata{Curve{1};} }
    """,
    # Bare Duplicata: coincident copies survive (no coherence pass).
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    Duplicata { Line{1}; }
    """,
    # Surface Duplicata: shared NEWREG tag, deep-copied boundary, then the
    # translate merge leaves the classic 5/6/10/14 corner wiring.
    """
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
    Translate {5,0,0} { Duplicata{Surface{1};} }
    """,
    # The post-transform coherence merge is global: an unrelated coincident
    # pair collapses even when it is not in the transformed set.
    """
    Point(1) = {0,0,0,1};
    Point(3) = {0,0,0,1};
    Point(4) = {5,0,0,1};
    Translate {5,0,0} { Point{4}; }
    """,
    # Coherence statement: global merge, and Coherence Point snaps listed
    # points onto the first tag before merging.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Point(3) = {0,0,0,1};
    Line(1) = {1,2};
    Coherence;
    Point(4) = {2,0,0,1};
    Point(5) = {9,0,0,1};
    Coherence Point{4,5};
    """,
    # Wildcards and an inline shape definition inside a transform list.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Translate {1,0,0} { Point{:}; }
    Translate {5,0,0} { Point(7) = {0,0,0,1}; }
    """,
    # A transform on nothing is legal; the merge still runs.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {0,0,0,1};
    Translate {1,0,0} { }
    """,
    # Volume transform on a materialized box recurses through the shell.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Point(3) = {1,1,0,1};
    Point(4) = {0,1,0,1};
    Point(5) = {0,0,1,1};
    Point(6) = {1,0,1,1};
    Point(7) = {1,1,1,1};
    Point(8) = {0,1,1,1};
    Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
    Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
    Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};
    Curve Loop(1) = {1,2,3,4};
    Curve Loop(2) = {5,6,7,8};
    Curve Loop(3) = {1,10,-5,-9};
    Curve Loop(4) = {2,11,-6,-10};
    Curve Loop(5) = {3,12,-7,-11};
    Curve Loop(6) = {4,9,-8,-12};
    Plane Surface(1) = {1}; Plane Surface(2) = {2};
    Plane Surface(3) = {3}; Plane Surface(4) = {4};
    Plane Surface(5) = {5}; Plane Surface(6) = {6};
    Surface Loop(1) = {1,2,3,4,5,6};
    Volume(1) = {1};
    Rotate {{0,0,1},{0,0,0}, Pi/2} { Volume{1}; }
    """,
    # OrientedBoundary selects the same absolute entities for transforms.
    """
    Point(1) = {0,0,0,1};
    Point(2) = {1,0,0,1};
    Line(1) = {1,2};
    Translate {0,0,3} { OrientedBoundary{Curve{1};} }
    """,
)

# Indices of cases containing arbitrary-angle Rotate statements; their point
# coordinates compare at COORD_ATOL instead of bit-for-bit (see header).
const ULP_CASES = Set((4, 16))
const COORD_ATOL = 64 * eps(1.0)

const GMSH_EXECUTABLE = find_gmsh_executable()
const GMSH_CLI_VERSION = strip(read(`$GMSH_EXECUTABLE --version`, String))
(GMSH_CLI_VERSION == TARGET_GMSH_VERSION ||
 startswith(GMSH_CLI_VERSION, TARGET_GMSH_VERSION * "-")) ||
    error("expected Gmsh $TARGET_GMSH_VERSION, got $GMSH_CLI_VERSION")
const GMSH_API_FILE = find_gmsh_api(GMSH_EXECUTABLE)
include(GMSH_API_FILE)
gmsh.GMSH_API_VERSION == TARGET_GMSH_VERSION || error(
    "expected Gmsh API $TARGET_GMSH_VERSION, got $(gmsh.GMSH_API_VERSION)")

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
                if case_index in ULP_CASES
                    all(abs.(actual .- expected) .<= COORD_ATOL) || error(
                        "case $case_index Point($tag) differs: " *
                        "Tessella=$actual Gmsh=$expected")
                else
                    all(reinterpret(UInt64, actual) .==
                        reinterpret(UInt64, expected)) || error(
                        "case $case_index Point($tag) differs: " *
                        "Tessella=$actual Gmsh=$expected")
                end
                samples[] += 1
            end
            # Curve endpoint wiring (orientation-insensitive) must match.
            for tag in sort!(collect(keys(model.curves)))
                expected = sort!(last.(gmsh.model.getBoundary(
                    [(1, tag)], false, false, false)))
                actual = sort!(collect(model.curves[tag]))
                actual == expected || error(
                    "case $case_index Curve($tag) endpoints differ: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
            # Surface boundary curves compare as sorted absolute tags; Gmsh
            # has no first-class loop entity.
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
            # Physical-group memberships must match.
            for (dim, tag) in gmsh.model.getPhysicalGroups()
                expected = sort!(gmsh.model.getEntitiesForPhysicalGroup(
                    dim, tag))
                actual = sort!(get(model.physical, (dim, tag), Int[]))
                actual == expected || error(
                    "case $case_index Physical dim $dim tag $tag differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
        end
    end

    println("GEO_TRANSFORMS_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(CASES)) samples=$(samples[]) " *
            "bit_exact=$(length(CASES) - length(ULP_CASES))/" *
            "$(length(CASES))")
finally
    gmsh.finalize()
end
