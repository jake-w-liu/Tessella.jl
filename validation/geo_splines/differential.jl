#!/usr/bin/env julia
# Entity-level differential for .geo spline-family curves. Each script is
# executed by Tessella's bounded GeoExec interpreter and independently opened
# by Gmsh 4.15.2's built-in parser; the resulting entity tags per dimension,
# every point's coordinates, curve endpoint wiring, `getType` strings,
# `getValue`/`getDerivative`/`getSecondDerivative`/`getCurvature` samples,
# parametrization bounds, 10-sample bounding boxes, `getParametrization`
# (parFromPoint), `getClosestPoint`, and `isInside` results must match
# bit-for-bit. Tessella rejects malformed Nurbs knot vectors and inferred
# degrees outside `0 <= degree <= npoints-1` at parse time — outside that
# range upstream `findSpan`/`basisFuns` index the control-point list out of
# bounds where Gmsh silently reads garbage — so every case here is a valid
# record by construction. NaN components compare positionally rather than
# bit-wise: both engines report the same degenerate queries undefined.
#
# The derivative/curvature and parFromPoint/closestPoint battery runs only on
# curves Gmsh reports as "Nurb" (the built-in kernel's type for the whole
# spline family). Gmsh differentiates a built-in Line by a 1e-5 forward
# difference and inverts it by XYZToU Newton iterates, while Tessella returns
# the exact rational delta and projects exactly — a pre-existing deliberate
# divergence on plain lines that this spline-family file does not re-litigate.
# A reversed Nurbs copy in a surface loop is likewise out of scope: Gmsh
# 4.15.2 segfaults inside `findSpan` on the mirrored knot vector while
# synchronizing, so only the record layout is compared (see
# `test/geometry/geo_spline_test.jl`).

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.Model: model_value, model_derivative, model_second_derivative,
                      model_curvature, model_parametrization_bounds,
                      model_parametrization, model_closest_point,
                      model_is_inside, model_entity_type, model_bounding_box

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

const ZIGZAG = """
    Point(1) = {0,0,0};
    Point(2) = {1,1,0};
    Point(3) = {2,0,0};
    Point(4) = {3,1,0};
    Point(5) = {4,0,0};
    """

const CASES = (
    # Catmull-Rom spline through the zig-zag polygon.
    ZIGZAG * "Spline(1) = {1,2,3,4,5};\n",
    # Two-point spline — extrapolated ghost ends collapse to a segment.
    "Point(1) = {0,0,0};\nPoint(2) = {1,1,1};\nSpline(1) = {1,2};\n",
    # Closed (periodic) spline: the first and last control points coincide,
    # taking Catmull-Rom's cyclic ghost branch (`c->beg == c->end`).
    """
    Point(1) = {1,0,0}; Point(2) = {0,1,0}; Point(3) = {-1,0,0};
    Point(4) = {0,-1,0};
    Spline(1) = {1,2,3,4,1};
    """,
    # Uniform B-splines at every matrix branch: N=5 (matext left/right),
    # N=6 (mat6_1/mat6_2), N=7 (mat7_2), N=8 (interior matbs/6 plus the
    # NbControlPoints>6 right-extremity reflection).
    ZIGZAG * "BSpline(1) = {1,2,3,4,5};\n",
    """
    Point(1) = {0,0,0}; Point(2) = {1,1,0}; Point(3) = {2,-1,0};
    Point(4) = {3,1,0}; Point(5) = {4,-1,0}; Point(6) = {5,0,0};
    BSpline(1) = {1,2,3,4,5,6};
    """,
    """
    Point(1) = {0,0,0}; Point(2) = {1,1,0}; Point(3) = {2,-1,0};
    Point(4) = {3,1,0}; Point(5) = {4,-1,0}; Point(6) = {5,0,0};
    Point(7) = {6,1,0};
    BSpline(1) = {1,2,3,4,5,6,7};
    """,
    """
    Point(1) = {0,0,0}; Point(2) = {1,1,0}; Point(3) = {2,-1,0};
    Point(4) = {3,1,0}; Point(5) = {4,-1,0}; Point(6) = {5,0,0};
    Point(7) = {6,1,0}; Point(8) = {7,-1,0};
    BSpline(1) = {1,2,3,4,5,6,7,8};
    """,
    # Degenerate B-spline counts: 2 cps -> linear, 3 -> quadratic Bezier,
    # 4 -> cubic Bezier matrices.
    "Point(1) = {0,0,0};\nPoint(2) = {2,1,0};\nBSpline(1) = {1,2};\n",
    """
    Point(1) = {0,0,0}; Point(2) = {1,2,0}; Point(3) = {2,0,0};
    BSpline(1) = {1,2,3};
    """,
    ZIGZAG * "BSpline(1) = {1,2,3,4};\n",
    # Closed (periodic) B-spline: the first and last control points coincide.
    """
    Point(1) = {1,0,0}; Point(2) = {0,1,0}; Point(3) = {-1,0,0};
    Point(4) = {0,-1,0};
    BSpline(1) = {1,2,3,4,1};
    """,
    # Full-list De Casteljau Bezier — 4 and 6 control points.
    ZIGZAG * "Bezier(1) = {1,2,3,4};\n",
    """
    Point(1) = {0,0,0}; Point(2) = {1,2,0}; Point(3) = {2,-1,0};
    Point(4) = {3,1,0}; Point(5) = {4,-2,0}; Point(6) = {5,0,0};
    Bezier(1) = {1,2,3,4,5,6};
    """,
    # Clamped cubic Nurbs on the raw [0,2] interval — same polygon as the
    # B-spline case, so the record is comparable but parametrized differently.
    ZIGZAG * "Nurbs(1) = {1,2,3,4,5} Knots {0,0,0,0,1,2,2,2,2} Order 3;\n",
    # Degree-2 Nurbs on a non-unit interval, and an Order expression that the
    # built-in kernel must evaluate but ignore.
    """
    Point(1) = {0,0,0}; Point(2) = {1,1,0}; Point(3) = {2,0,0};
    Point(4) = {3,1,0};
    Nurbs(1) = {1,2,3,4} Knots {0.5,0.5,0.5,1.25,2.5,2.5,2.5} Order 1+1;
    """,
    # Degree-1 Nurbs — piecewise linear through the control polygon
    # (5 points, 7 knots -> degree 1).
    ZIGZAG * "Nurbs(1) = {1,2,3,4,5} Knots {0,0,0.25,0.5,0.75,1,1} Order 1;\n",
    # Empty `Knots {}` routes to the plain BSpline record.
    ZIGZAG * "Nurbs(1) = {1,2,3,4,5} Knots {} Order 9;\n",
    # List variables and numeric expressions inside the control list, and a
    # `Nurbs` whose `ListOfDouble` positions are list variables on both sides
    # of `Knots`.
    """
    Point(1) = {0,0,0}; Point(2) = {1,1,0}; Point(3) = {2,0,0};
    Point(4) = {3,1,0}; Point(5) = {4,0,0};
    pts[] = {1,2,3,4,5};
    kk[] = {0,0,0,0,1,1,1,1};
    Spline(7) = pts[];
    Bezier(8) = {1+0,2,3,4};
    Nurbs(9) = pts[] Knots kk[] Order 7;
    """,
    # `ListOfDouble` entity-list forms: a negated `Curve Loop` and a
    # multiplier-scaled `Spline` control list.
    """
    Point(1) = {0,0,0}; Point(2) = {1,1,0}; Point(3) = {2,0,0};
    Point(4) = {2,-1,0};
    Spline(1) = 1*{1,2,3};
    Line(2) = {3,4};
    Line(3) = {4,1};
    Curve Loop(1) = -{3,2,1};
    Plane Surface(1) = {1};
    """,
    # Duplicata of a surface whose loop traverses the spline backwards: the
    # copy carries the reversed record under a positive tag.
    """
    Point(1) = {0,0,0}; Point(2) = {1,1,0}; Point(3) = {2,0,0};
    Point(4) = {2,-1,0};
    Spline(1) = {1,2,3};
    Line(2) = {3,4};
    Line(3) = {4,1};
    Curve Loop(1) = {-3,-2,-1};
    Plane Surface(1) = {1};
    Duplicata { Surface{1}; }
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

# Query points for projection/containment — on-curve, near-curve, and far
# off-curve coordinates exercised on every spline-family curve.
const QUERIES = ([2.0,0.5,0.0],[1.0,0.9,0.0],[9.0,9.0,9.0],[-0.4,0.6,0.0])

# Bit-identical comparison per component, except that NaN matches NaN:
# degenerate queries (a 0/0 inside `basisFuns` at a clamped end, say) produce
# NaN on both sides through different arithmetic paths whose payloads need not
# agree — the parity claim is that both engines report the query undefined at
# the same positions.
function bit_equal(actual, expected)
    a=collect(actual);e=collect(expected)
    length(a)==length(e) || return false
    for i in eachindex(a)
        (isnan(a[i]) && isnan(e[i])) && continue
        reinterpret(UInt64,a[i])==reinterpret(UInt64,e[i]) || return false
    end
    return true
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
                bit_equal(actual, expected) || error(
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
                # Sampled on the entity's own parameter interval so Nurbs
                # curves exercise their raw knot range.
                lo, hi = gmsh.model.getParametrizationBounds(1, tag)
                alo, ahi = model_parametrization_bounds(model, 1, tag)
                bit_equal(alo, lo) && bit_equal(ahi, hi) || error(
                    "case $case_index Curve($tag) bounds differ: " *
                    "Tessella=($alo,$ahi) Gmsh=($lo,$hi)")
                samples[] += 1
                nurb = expected_type == "Nurb"
                us = range(lo[1], hi[1]; length = 7)
                for u in us
                    expected = gmsh.model.getValue(1, tag, [u])
                    actual = model_value(model, 1, tag, [u])
                    bit_equal(actual, expected) || error(
                        "case $case_index Curve($tag) eval u=$u differs: " *
                        "Tessella=$actual Gmsh=$expected")
                    if nurb
                        expected = gmsh.model.getDerivative(1, tag, [u])
                        actual = model_derivative(model, 1, tag, [u])
                        bit_equal(actual, expected) || error(
                            "case $case_index Curve($tag) der u=$u " *
                            "differs: Tessella=$actual Gmsh=$expected")
                        expected = gmsh.model.getSecondDerivative(
                            1, tag, [u])
                        actual = model_second_derivative(model, 1, tag, [u])
                        bit_equal(actual, expected) || error(
                            "case $case_index Curve($tag) der2 u=$u " *
                            "differs: Tessella=$actual Gmsh=$expected")
                        expected = gmsh.model.getCurvature(1, tag, [u])
                        actual = model_curvature(model, 1, tag, [u])
                        bit_equal(actual, expected) || error(
                            "case $case_index Curve($tag) curv u=$u " *
                            "differs: Tessella=$actual Gmsh=$expected")
                        samples[] += 3
                    end
                    samples[] += 1
                end
                eb = gmsh.model.getBoundingBox(1, tag)
                ab = model_bounding_box(model, 1, tag)
                bit_equal(collect(ab), collect(eb)) || error(
                    "case $case_index Curve($tag) bounding box differs: " *
                    "Tessella=$ab Gmsh=$eb")
                samples[] += 1
                if nurb
                    for q in QUERIES
                        expected = gmsh.model.getParametrization(1, tag, q)
                        actual = model_parametrization(model, 1, tag, q)
                        bit_equal(actual, expected) || error(
                            "case $case_index Curve($tag) parFromPoint($q) " *
                            "differs: Tessella=$actual Gmsh=$expected")
                        ecoord, epar = gmsh.model.getClosestPoint(1, tag, q)
                        acoord, apar = model_closest_point(model, 1, tag, q)
                        bit_equal(acoord, ecoord) && bit_equal(apar, epar) ||
                            error(
                                "case $case_index Curve($tag) closest($q) " *
                                "differs: Tessella=($acoord,$apar) " *
                                "Gmsh=($ecoord,$epar)")
                        expected = gmsh.model.isInside(1, tag, q)
                        actual = model_is_inside(model, 1, tag, q)
                        actual == expected || error(
                            "case $case_index Curve($tag) isInside($q) " *
                            "differs: Tessella=$actual Gmsh=$expected")
                        samples[] += 3
                    end
                end
            end
            for tag in sort!(collect(keys(model.surfaces)))
                expected_type = gmsh.model.getType(2, tag)
                actual_type = model_entity_type(model, 2, tag)
                actual_type == expected_type || error(
                    "case $case_index Surface($tag) type differs: " *
                    "Tessella=$actual_type Gmsh=$expected_type")
                # Oriented boundary — reversed generatrices must arrive as
                # positive-tag copies of the reversed records.
                expected = last.(gmsh.model.getBoundary(
                    [(2, tag)], false, true, false))
                actual = vcat((model.loops[l]
                               for l in model.surfaces[tag])...)
                actual == expected || error(
                    "case $case_index Surface($tag) oriented boundary " *
                    "differs: Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
        end
    end

    println("GEO_SPLINES_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(CASES)) samples=$(samples[])")
finally
    gmsh.finalize()
end
