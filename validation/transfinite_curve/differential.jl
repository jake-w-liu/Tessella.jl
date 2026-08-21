#!/usr/bin/env julia
# In-memory straight-line differential for Gmsh 4.15.2 transfinite curve laws.
# No geometry or mesh file is read or written, and Gmsh is always finalized.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using SHA: sha256

if !isdefined(Tessella, :TransfiniteCurve)
    Base.include(Tessella,
                 joinpath(@__DIR__, "..", "..", "src", "TransfiniteCurve.jl"))
end
using Tessella.TransfiniteCurve: transfinite_curve_parameters

const TARGET_GMSH_VERSION = "4.15.2"
const PARAMETER_CRC =
    "b525628945d55f49bc5151ad313d697f9c32f7b3325015e5d0080a69260ffd0e"

function find_gmsh_api()
    configured = get(ENV, "GMSH_JULIA_API", "")
    if !isempty(configured)
        isfile(configured) || error(
            "GMSH_JULIA_API does not name a file: $configured")
        return realpath(configured)
    end
    candidates = String[]
    executable = Sys.which("gmsh")
    if executable !== nothing
        prefix = dirname(dirname(realpath(executable)))
        append!(candidates, (joinpath(prefix, "lib", "gmsh.jl"),
                             joinpath(prefix, "lib64", "gmsh.jl")))
    end
    append!(candidates, ("/opt/homebrew/lib/gmsh.jl",
                         "/opt/homebrew/opt/gmsh/lib/gmsh.jl",
                         "/usr/local/opt/gmsh/lib/gmsh.jl"))
    for candidate in unique(candidates)
        isfile(candidate) && return realpath(candidate)
    end
    error("Gmsh 4.15.2 Julia API not found; set GMSH_JULIA_API")
end

include(find_gmsh_api())
gmsh.GMSH_API_VERSION == TARGET_GMSH_VERSION || error(
    "expected Gmsh API $TARGET_GMSH_VERSION, got $(gmsh.GMSH_API_VERSION)")

const CASES = (
    (name=:progression_2_n6, law=:progression, gmsh_type="Progression",
     coefficient=2.0, count=6),
    (name=:progression_negative_2_n6, law=:progression,
     gmsh_type="Progression", coefficient=-2.0, count=6),
    (name=:progression_half_n6, law=:progression, gmsh_type="Progression",
     coefficient=0.5, count=6),
    (name=:progression_negative_half_n6, law=:progression,
     gmsh_type="Progression", coefficient=-0.5, count=6),
    (name=:progression_1_5_n9, law=:progression, gmsh_type="Progression",
     coefficient=1.5, count=9),
    (name=:power_alias_n6, law=:power, gmsh_type="Power",
     coefficient=2.0, count=6),
    (name=:bump_2_n6, law=:bump, gmsh_type="Bump",
     coefficient=2.0, count=6),
    (name=:bump_negative_2_n6, law=:bump, gmsh_type="Bump",
     coefficient=-2.0, count=6),
    (name=:bump_0_05_n9, law=:bump, gmsh_type="Bump",
     coefficient=0.05, count=9),
    (name=:bump_0_05_n17, law=:bump, gmsh_type="Bump",
     coefficient=0.05, count=17),
    (name=:bump_half_n17, law=:bump, gmsh_type="Bump",
     coefficient=0.5, count=17),
    (name=:bump_20_n17, law=:bump, gmsh_type="Bump",
     coefficient=20.0, count=17),
    (name=:beta_2_n6, law=:beta, gmsh_type="Beta",
     coefficient=2.0, count=6),
    (name=:beta_negative_2_n6, law=:beta, gmsh_type="Beta",
     coefficient=-2.0, count=6),
    (name=:beta_1_01_n9, law=:beta, gmsh_type="Beta",
     coefficient=1.01, count=9),
    (name=:beta_1_01_n17, law=:beta, gmsh_type="Beta",
     coefficient=1.01, count=17),
    (name=:beta_1_2_n17, law=:beta, gmsh_type="Beta",
     coefficient=1.2, count=17),
    (name=:beta_20_n17, law=:beta, gmsh_type="Beta",
     coefficient=20.0, count=17),
    (name=:beta_half_n17, law=:beta, gmsh_type="Beta",
     coefficient=0.5, count=17),
    (name=:beta_negative_half_n17, law=:beta, gmsh_type="Beta",
     coefficient=-0.5, count=17),
    (name=:uniform_two_nodes, law=:progression, gmsh_type="Progression",
     coefficient=1.0, count=2),
)

function parameter_crc(groups)
    stream = IOBuffer()
    for values in groups
        write(stream, htol(UInt64(length(values))))
        for value in values
            write(stream, htol(reinterpret(UInt64, value)))
        end
    end
    return bytes2hex(sha256(take!(stream)))
end

function gmsh_parameters(case)
    gmsh.clear()
    gmsh.model.add(String(case.name))
    first = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    last = gmsh.model.geo.addPoint(1.0, 0.0, 0.0)
    curve = gmsh.model.geo.addLine(first, last)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setTransfiniteCurve(
        curve, case.count, case.gmsh_type, case.coefficient)
    gmsh.model.mesh.generate(1)

    node_tags, coordinates, parameters =
        gmsh.model.mesh.getNodes(1, curve, true, true)
    length(node_tags) == case.count || error(
        "$(case.name): expected $(case.count) Gmsh nodes, got $(length(node_tags))")
    length(coordinates) == 3case.count || error(
        "$(case.name): malformed Gmsh coordinate block")
    length(parameters) == case.count || error(
        "$(case.name): malformed Gmsh parameter block")
    order = sortperm(parameters)
    result = Vector{Float64}(undef, case.count)
    for (destination, source) in pairs(order)
        x = coordinates[3source - 2]
        y = coordinates[3source - 1]
        z = coordinates[3source]
        isfinite(x) && abs(y) <= 8eps(Float64) && abs(z) <= 8eps(Float64) ||
            error("$(case.name): Gmsh returned an invalid unit-line coordinate")
        abs(x - parameters[source]) <= 64eps(Float64) || error(
            "$(case.name): line coordinate and CAD parameter differ")
        result[destination] = x
    end
    result[1] == 0.0 && result[end] == 1.0 || error(
        "$(case.name): Gmsh did not preserve both endpoints")
    all(>(0.0), diff(result)) || error(
        "$(case.name): Gmsh parameters are not strictly increasing")

    element_types, element_tags, element_nodes =
        gmsh.model.mesh.getElements(1, curve)
    line_block = findfirst(==(Int32(1)), element_types)
    line_block === nothing && error("$(case.name): Gmsh emitted no line elements")
    length(element_tags[line_block]) == case.count - 1 || error(
        "$(case.name): Gmsh line-element count differs")
    length(element_nodes[line_block]) == 2(case.count - 1) || error(
        "$(case.name): malformed Gmsh line connectivity")
    return result
end

function maximum_difference(left, right)
    length(left) == length(right) || error("parameter-count mismatch")
    return maximum((abs(left[index] - right[index]) for index in eachindex(left));
                   init=0.0)
end

gmsh.initialize(String[], false, false)
try
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.option.setNumber("Mesh.FlexibleTransfinite", 0)
    gmsh.option.setNumber("Mesh.LcIntegrationPrecision", 1e-9)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    runtime_version = gmsh.option.getString("General.Version")
    (runtime_version == TARGET_GMSH_VERSION ||
     startswith(runtime_version, TARGET_GMSH_VERSION * "-")) || error(
        "expected Gmsh runtime $TARGET_GMSH_VERSION, got $runtime_version")

    gmsh_results = Dict{Symbol,Vector{Float64}}()
    tessella_results = Dict{Symbol,Vector{Float64}}()
    maximum_error = 0.0
    for case in CASES
        tessella = transfinite_curve_parameters(
            case.count; mesh_type=case.law, coefficient=case.coefficient,
            max_nodes=case.count)
        reference = gmsh_parameters(case)
        error_value = maximum_difference(tessella, reference)
        error_value <= 1e-7 || error(
            "$(case.name): analytic/Gmsh error $error_value exceeds 1e-7")
        maximum_error = max(maximum_error, error_value)
        tessella_results[case.name] = tessella
        gmsh_results[case.name] = reference
    end

    maximum_difference(gmsh_results[:progression_half_n6],
                       gmsh_results[:progression_negative_2_n6]) == 0.0 || error(
        "Gmsh progression coefficient direction/order semantics changed")
    maximum_difference(gmsh_results[:progression_2_n6],
                       gmsh_results[:power_alias_n6]) == 0.0 || error(
        "Gmsh Power alias differs from Progression")
    maximum_difference(gmsh_results[:bump_2_n6],
                       gmsh_results[:bump_negative_2_n6]) == 0.0 || error(
        "Gmsh Bump unexpectedly depends on coefficient sign")
    maximum_difference(
        gmsh_results[:beta_half_n17], collect(range(0.0, 1.0; length=17))) <=
        2e-12 || error("Gmsh Beta coefficient <= 1 is no longer uniform")
    maximum_difference(
        gmsh_results[:beta_negative_half_n17],
        gmsh_results[:beta_half_n17]) == 0.0 || error(
        "Gmsh subunit Beta fallback unexpectedly depends on coefficient sign")

    checksum_groups = (
        tessella_results[:progression_2_n6],
        tessella_results[:progression_negative_2_n6],
        tessella_results[:bump_2_n6],
        tessella_results[:beta_2_n6],
        tessella_results[:beta_negative_2_n6],
        tessella_results[:bump_0_05_n9],
        tessella_results[:beta_1_01_n9])
    checksum = parameter_crc(checksum_groups)
    checksum == PARAMETER_CRC || error(
        "Tessella transfinite-curve checksum changed: $checksum")

    println("TRANSFINITE_CURVE_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "laws=3 cases=$(length(CASES)) coordinate_samples=" *
            "$(sum(case.count for case in CASES; init=0)) " *
            "max_abs_error=$maximum_error sha=$checksum")
finally
    gmsh.isInitialized() != 0 && gmsh.finalize()
end

# Intentional nonclaims: non-affine CAD curves and their derivative-weighted
# adaptive integration, FlexibleTransfinite, *HWall/size-map laws, closed,
# periodic, extruded, degenerate and boundary-layer curves, and bitwise equality
# with Gmsh's adaptive-trapezoid/linear-inversion roundoff.
