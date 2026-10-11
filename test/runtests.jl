using Test
using LinearAlgebra
using Tessella
import Tessella.MeshTypes as MT

# Small public operations and analytic invariants; opt-in campaigns: README.md.
include("functional/fixtures.jl")

@testset verbose=true "Tessella functionality" begin
    for feature in ("core", "geometry", "meshing", "io", "api")
        println("Testing ", feature, "...")
        flush(stdout)
        elapsed = @elapsed include(joinpath("functional", feature * ".jl"))
        println("Completed ", feature, " in ", round(elapsed; digits=2), " s")
        flush(stdout)
    end
    include("meshing/ordered_matching_seeded_test.jl")
    include("geometry/model_compact_atomic_coverage_test.jl")
    include("functional/list_view_export.jl")
end
