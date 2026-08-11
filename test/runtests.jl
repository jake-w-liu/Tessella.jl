using Test
using Tessella

@testset "Tessella" begin
    # Stage 0 — Foundations (CRC-gated, DEVELOPMENT.md discipline).
    include("predicates_test.jl")   # exact predicates vs exact-rational oracle
    include("meshtypes_test.jl")    # mesh container, topology, quality, checksum
    include("io_test.jl")           # .msh v2/v4 round-trip, STL, .geo scan

    # Stage 1 — 2-D meshing (CRC-gated).
    include("mesh2d_test.jl")       # Delaunay: exact empty-circumcircle oracle

    # Stage 2 — 1-D edge meshing + surface meshing (CRC-gated).
    include("mesh1d_test.jl")       # size fields + graded edge meshing vs arc length
    include("meshsurface_test.jl")  # planar/cylinder/parametric surface meshing

    # Stage 3 — 3-D meshing (CRC-gated; boundary recovery lands incrementally).
    include("mesh3d_test.jl")       # 3-D Delaunay + volume filling + coax junction

    # Stage 4 — optimization + quality (CRC-gated; sliver removal WIP).
    include("optimize_test.jl")     # tet quality report + Laplacian & ODT smoothing

    # Stage 5 — heal + native primitives + top-level pipeline.
    include("heal_test.jl")         # surface diagnostics (open/non-manifold/degenerate…)
    include("geometry_test.jl")     # box / cylinder / box-tunnel primitives
    include("pipeline_test.jl")     # mesh_volume: validated-or-explicit-blocker

    # Stage 6 — high-order elements.
    include("highorder_test.jl")    # quadratic (P2) tet generation + type-11 I/O

    @testset "stage banner" begin
        @test Tessella.stage() isa Int
        @test Tessella.stage() >= 1     # Stage 0 gate is green
    end
end
