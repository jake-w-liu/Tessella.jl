#!/usr/bin/env julia
# Entity-level differential for materialized Cylinder/Sphere/Cone boundaries.
# Gmsh 4.15.2's built-in .geo kernel rejects solid primitives (they are
# OCC-only), so each Tessella `.geo` execution is compared against the
# equivalent `gmsh.model.occ.addCylinder/addSphere/addCone` model: entity
# tags per dimension, point coordinates, curve/surface type names, signed
# boundary wiring, shell signs, parametrization bounds, curve evaluations,
# and bounding boxes must all match — coordinates/evaluations bitwise where
# the same IEEE operations produce them, within 1e-12 otherwise.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.Model: model_boundary, model_entity_type,
                      model_parametrization_bounds, model_value,
                      model_bounding_box

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

# Each case: the Tessella .geo source and the matching OCC builder calls.
const CASES = (
    (source="Cylinder(1) = {0,0,0,0,0,2,1};",
     occ=()->gmsh.model.occ.addCylinder(0,0,0,0,0,2,1,1)),
    # Tilted axis exercises the OCC reference-direction frame.
    (source="Cylinder(1) = {1,2,3,-2,4,5,3};",
     occ=()->gmsh.model.occ.addCylinder(1,2,3,-2,4,5,3,1)),
    (source="Cylinder(1) = {0,0,0,0,2,0,1};",
     occ=()->gmsh.model.occ.addCylinder(0,0,0,0,2,0,1,1)),
    (source="Sphere(1) = {0,0,0,1};",
     occ=()->gmsh.model.occ.addSphere(0,0,0,1,1)),
    (source="Sphere(1) = {1,-2,0.5,2.5};",
     occ=()->gmsh.model.occ.addSphere(1,-2,0.5,2.5,1)),
    (source="Cone(1) = {0,0,0,0,0,2,2,1};",
     occ=()->gmsh.model.occ.addCone(0,0,0,0,0,2,2,1,1)),
    (source="Cone(1) = {1,2,3,-2,4,5,6,2};",
     occ=()->gmsh.model.occ.addCone(1,2,3,-2,4,5,6,2,1)),
    # r2 == 0: degenerate apex edge on top, no top cap.
    (source="Cone(1) = {0,0,0,0,0,2,2,0};",
     occ=()->gmsh.model.occ.addCone(0,0,0,0,0,2,2,0,1)),
    # r1 == 0: degenerate apex edge at the base, top cap only.
    (source="Cone(1) = {0,0,0,0,0,2,0,2};",
     occ=()->gmsh.model.occ.addCone(0,0,0,0,0,2,0,2,1)),
    (source="Torus(1) = {0,0,0,3,1};",
     occ=()->gmsh.model.occ.addTorus(0,0,0,3,1,1)),
    (source="Torus(1) = {1,-2,0.5,5,1.5};",
     occ=()->gmsh.model.occ.addTorus(1,-2,0.5,5,1.5,1)),
    # Spindle torus (r1 < r2): self-intersecting tube, same OCC layout.
    (source="Torus(1) = {0,0,0,2,3};",
     occ=()->gmsh.model.occ.addTorus(0,0,0,2,3,1)),
    # angle == 2π spelled out still materializes the full-torus layout.
    (source="Torus(1) = {0,0,0,3,1,6.283185307179586};",
     occ=()->gmsh.model.occ.addTorus(0,0,0,3,1,1,2*pi)),
    # Partial torus: two rim vertices, trimmed equator, closed meridians,
    # Plane caps under shell [torus,+start,-end].
    (source="Torus(1) = {0,0,0,3,1,1.5707963267948966};",
     occ=()->gmsh.model.occ.addTorus(0,0,0,3,1,1,pi/2)),
    (source="Torus(1) = {1,2,3,4,0.75,4.2};",
     occ=()->gmsh.model.occ.addTorus(1,2,3,4,0.75,1,4.2)),
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

samples = Ref(0)
gmsh.initialize([GMSH_EXECUTABLE, "-nopopup"], false, false)
try
    gmsh.option.setNumber("General.Terminal", 0)
    runtime_version = gmsh.option.getString("General.Version")
    (runtime_version == TARGET_GMSH_VERSION ||
     startswith(runtime_version, TARGET_GMSH_VERSION * "-")) || error(
        "expected Gmsh runtime $TARGET_GMSH_VERSION, got $runtime_version")

    mktempdir() do directory
        for (case_index, case) in enumerate(CASES)
            path = joinpath(directory, "case$case_index.geo")
            write(path, case.source)

            parsed = Tessella.GeoExec.execute_geo(path)
            model = parsed.model

            gmsh.clear()
            case.occ()
            gmsh.model.occ.synchronize()

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
                all(abs.(actual .- expected) .< 1e-12) || error(
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
                expected = last.(gmsh.model.getBoundary(
                    [(1, tag)], false, true, false))
                actual = collect(model_boundary(
                    model, [(1, tag)], false, true, false) .|> last)
                actual == expected || error(
                    "case $case_index Curve($tag) boundary differs: " *
                    "Tessella=$actual Gmsh=$expected")
                expected_bounds = gmsh.model.getParametrizationBounds(1, tag)
                actual_bounds = model_parametrization_bounds(model, 1, tag)
                (actual_bounds[1] ≈ expected_bounds[1] &&
                 actual_bounds[2] ≈ expected_bounds[2]) || error(
                    "case $case_index Curve($tag) bounds differ: " *
                    "Tessella=$actual_bounds Gmsh=$expected_bounds")
                expected_type == "Unknown" && continue
                t0, t1 = only(actual_bounds[1]), only(actual_bounds[2])
                for t in (t0, t0 + (t1 - t0) / 3, t0 + 2(t1 - t0) / 3, t1)
                    expected = gmsh.model.getValue(1, tag, [t])
                    actual = model_value(model, 1, tag, [t])
                    all(abs.(actual .- expected) .< 1e-12) || error(
                        "case $case_index Curve($tag) eval($t) differs: " *
                        "Tessella=$actual Gmsh=$expected")
                    samples[] += 1
                end
                samples[] += 1
            end
            for tag in sort!(collect(keys(model.surfaces)))
                expected_type = gmsh.model.getType(2, tag)
                actual_type = model_entity_type(model, 2, tag)
                actual_type == expected_type || error(
                    "case $case_index Surface($tag) type differs: " *
                    "Tessella=$actual_type Gmsh=$expected_type")
                expected = last.(gmsh.model.getBoundary(
                    [(2, tag)], false, true, false))
                actual = collect(model_boundary(
                    model, [(2, tag)], false, true, false) .|> last)
                actual == expected || error(
                    "case $case_index Surface($tag) boundary differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
                # Analytic OCC faces share parametrization; Plane fits are
                # basis-dependent and stay unchecked.
                expected_type in ("Cylinder", "Cone", "Sphere", "Torus") ||
                    continue
                expected_bounds = gmsh.model.getParametrizationBounds(2, tag)
                actual_bounds = model_parametrization_bounds(model, 2, tag)
                (all(actual_bounds[1] .≈ expected_bounds[1]) &&
                 all(actual_bounds[2] .≈ expected_bounds[2])) || error(
                    "case $case_index Surface($tag) bounds differ: " *
                    "Tessella=$actual_bounds Gmsh=$expected_bounds")
                u0, u1 = actual_bounds[1][1], actual_bounds[2][1]
                v0, v1 = actual_bounds[1][2], actual_bounds[2][2]
                for (fu, fv) in ((0.0, 0.5), (0.5, 0.0), (1 / 3, 2 / 3),
                                 (1.0, 0.5))
                    uv = [u0 + (u1 - u0) * fu, v0 + (v1 - v0) * fv]
                    expected = gmsh.model.getValue(2, tag, uv)
                    actual = model_value(model, 2, tag, uv)
                    all(abs.(actual .- expected) .< 1e-12) || error(
                        "case $case_index Surface($tag) eval($uv) differs: " *
                        "Tessella=$actual Gmsh=$expected")
                    samples[] += 1
                end
            end
            for tag in sort!(collect(keys(model.volumes)))
                expected_shell = last.(gmsh.model.getBoundary(
                    [(3, tag)], false, true, false))
                actual = last.(model_boundary(
                    model, [(3, tag)], true, true, false))
                actual == expected_shell || error(
                    "case $case_index Volume($tag) shell differs: " *
                    "Tessella=$actual Gmsh=$expected_shell")
                # `getBoundingBox` reports OCC's padded box (1e-7 absolute
                # tolerance); Tessella's analytic box must sit inside it
                # while staying tight itself. Torus faces fall back to a
                # coarse polyhedral bound in OCC, so torus-shelled volumes
                # are checked for containment only.
                expected = collect(gmsh.model.getBoundingBox(3, tag))
                actual = collect(model_bounding_box(model, 3, tag))
                polyhedral = any(gmsh.model.getType(2, s) == "Torus"
                                 for s in abs.(expected_shell))
                (all(actual[1:3] .>= expected[1:3] .- 1e-12) &&
                 all(actual[4:6] .<= expected[4:6] .+ 1e-12) &&
                 (polyhedral ||
                  all(abs.(actual .- expected) .< 1e-6))) || error(
                    "case $case_index Volume($tag) bbox differs: " *
                    "Tessella=$actual Gmsh=$expected")
                samples[] += 1
            end
        end
    end

    println("GEO_PRIMITIVES_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(CASES)) samples=$(samples[])")
finally
    gmsh.finalize()
end
