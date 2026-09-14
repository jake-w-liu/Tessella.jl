#!/usr/bin/env julia
# Entity-level differential for .geo control flow. Each script is executed by
# Tessella's bounded GeoExec interpreter and independently opened by Gmsh
# 4.15.2's built-in parser; the resulting entity tags per dimension and every
# point's coordinates must match bit-for-bit.

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
    # If/ElseIf/Else selection, same-line form, nesting.
    """
    x = 1;
    If (x > 0)
      Point(50) = {5, 0, 0, 1};
    ElseIf (x > 2)
      Point(51) = {6, 0, 0, 1};
    Else
      Point(52) = {7, 0, 0, 1};
    EndIf
    If (x < 0)
      Point(53) = {0, 0, 0, 1};
    ElseIf (x == 1)
      Point(54) = {54, 0, 0, 1};
    Else
      Point(55) = {0, 0, 0, 1};
    EndIf
    If (0) Point(70) = {0,0,0,1}; Else Point(71) = {1,0,0,1}; EndIf
    """,
    # For ranges: implicit +1, explicit step, descending step, empty range,
    # loop-variable persistence at the first out-of-range value.
    """
    i = 7;
    For i In {0:2}
      Point(i+1) = {i, 0, 0, 1};
    EndFor
    Point(99) = {i, 0, 0, 1};
    For k In {0:4:2}
      Point(k+10) = {k, 0, 0, 1};
    EndFor
    For k In {5:0:-1}
      Point(80+k) = {k, 0, 0, 1};
    EndFor
    Point(91) = {k, 0, 0, 1};
    For m In {5:0}
      Point(200+m) = {m, 0, 0, 1};
    EndFor
    Point(92) = {m, 0, 0, 1};
    """,
    # Nested For/If and loop bounds evaluated once at entry.
    """
    For i In {0:1}
      If (i == 0)
        For j In {10:11}
          Point(i*100+j) = {j, 0, 0, 1};
        EndFor
      Else
        Point(9) = {9, 0, 0, 1};
      EndIf
    EndFor
    n = 2;
    For i In {0:n}
      n = 100;
      Point(300+i) = {i, 0, 0, 1};
    EndFor
    """,
    # Conditional operators inside geometry expressions, plus allocator reads.
    """
    x = 1;
    For i In {0:2}
      Point(newp) = {(i > 0) ? i*10 : -1, !(x == 1), (x >= 1 && i < 2) || 0, 1};
    EndFor
    Point(newp) = {x != 2, 1 < 2 == 1, 0, 1};
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
        end
    end

    println("GEO_CONTROL_FLOW_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(CASES)) samples=$(samples[]) bit_exact=1")
finally
    gmsh.finalize()
end
