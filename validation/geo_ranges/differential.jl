#!/usr/bin/env julia
# Pointwise differential for bounded constant Gmsh list ranges. The same field
# option lists are parsed independently by Tessella and Gmsh 4.15.2 and compared
# bit-for-bit; no geometry or mesh file is generated.

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

function tessella_values(params, tag, option)
    raw = params.fields[tag].options[option]
    entries = Tessella.IO._geo_split_list(
        raw, "geo-range differential Field[$tag].$option")
    return parse.(Float64, entries)
end

function bit_equal(first, second)
    length(first) == length(second) || return false
    return all(reinterpret(UInt64, first[index]) ==
               reinterpret(UInt64, second[index])
               for index in eachindex(first))
end

const FLOAT_CASES = (
    "1:5",
    "5:1",
    "3:3",
    "1:5:2",
    "5:1:-2",
    "5:-2:1",
    "-2:2",
    "0:0.25:0.1",
    "1 / 2:3 / 2:1 / 4",
    "Min(2, 3):Max(4, 5)",
    "1e-300:3e-300:1e-300",
    "(1 + 1):(3 * 2):2",
    "(5.5 % 2.2):(5.5 % 2.2)",
)

const INTEGER_CASES = (
    "1:3:0.5",
    "5:1:-0.75",
    "2:2",
    "4:-1:1",
)

const WRAPPED_CASES = (
    ("BoundaryLayer", "SizesList", "-{-0.1:-0.3:-0.1}"),
    ("Distance", "PointsList", "-{-1:-3}"),
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
        path = joinpath(directory, "ranges.geo")
        open(path, "w") do io
            tag = 0
            for expression in FLOAT_CASES
                tag += 1
                println(io, "Field[$tag] = BoundaryLayer;")
                println(io, "Field[$tag].SizesList = {$expression};")
            end
            for expression in INTEGER_CASES
                tag += 1
                println(io, "Field[$tag] = Distance;")
                println(io, "Field[$tag].PointsList = {$expression};")
            end
            for (kind, option, expression) in WRAPPED_CASES
                tag += 1
                println(io, "Field[$tag] = $kind;")
                println(io, "Field[$tag].$option = $expression;")
            end
        end

        params = Tessella.IO.read_geo_params(path)
        gmsh.open(path)
        tag = 0
        for _ in FLOAT_CASES
            tag += 1
            expected = Float64.(gmsh.model.mesh.field.getNumbers(tag, "SizesList"))
            actual = tessella_values(params, tag, "SizesList")
            bit_equal(actual, expected) || error(
                "Field[$tag].SizesList differs: Tessella=$actual Gmsh=$expected")
            samples[] += length(expected)
        end
        for _ in INTEGER_CASES
            tag += 1
            expected = Float64.(gmsh.model.mesh.field.getNumbers(tag, "PointsList"))
            actual = tessella_values(params, tag, "PointsList")
            bit_equal(actual, expected) || error(
                "Field[$tag].PointsList differs: Tessella=$actual Gmsh=$expected")
            samples[] += length(expected)
        end
        for (_, option, _) in WRAPPED_CASES
            tag += 1
            expected = Float64.(gmsh.model.mesh.field.getNumbers(tag, option))
            actual = tessella_values(params, tag, option)
            bit_equal(actual, expected) || error(
                "Field[$tag].$option differs: Tessella=$actual Gmsh=$expected")
            samples[] += length(expected)
        end
    end

    println("GEO_RANGE_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "float_cases=$(length(FLOAT_CASES)) " *
            "integer_cases=$(length(INTEGER_CASES)) " *
            "wrapped_cases=$(length(WRAPPED_CASES)) " *
            "samples=$(samples[]) bit_exact=1")
finally
    gmsh.finalize()
end
