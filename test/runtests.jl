using Test
using Tessella

@testset "Tessella" begin
    # Stage 0 — Foundations (CRC-gated, DEVELOPMENT.md discipline).
    include("predicates_test.jl")   # exact predicates vs exact-rational oracle
    include("meshtypes_test.jl")    # mesh container, topology, quality, checksum
    include("transform_test.jl")    # validated affine transforms + orientation preservation
    include("io_test.jl")           # .msh v2/v4 round-trip, STL, .geo scan
    include("elements_test.jl")     # fixed/special Gmsh 4.15.2 records + mixed entity I/O
    include("recombine_test.jl")    # deterministic triangle-to-quad recombination
    include("refine_test.jl")       # deterministic one-level uniform simplex refinement
    include("transfinite_test.jl")  # validated four-sided planar structured patches
    include("transfinite_curve_test.jl") # straight-curve Progression/Bump/Beta laws

    # Stage 1 — 2-D meshing (CRC-gated).
    include("mesh2d_test.jl")       # Delaunay: exact empty-circumcircle oracle

    # Stage 2 — Gmsh-compatible size fields + 1-D/surface meshing (CRC-gated).
    include("sizefield_test.jl") # Distance/Threshold/Box/Min fields + local 3-D sizing
    include("mesh1d_test.jl")       # size fields + graded edge meshing vs arc length
    include("meshsurface_test.jl")  # planar/cylinder/parametric surface meshing

    # Stage 3 — 3-D meshing and certified boundary recovery (CRC-gated).
    include("mesh3d_test.jl")       # 3-D Delaunay + volume filling + coax junction

    # Stage 4 — optimization, quality, and sliver removal (CRC-gated).
    include("optimize_test.jl")     # tet quality report + Laplacian & ODT smoothing

    # Stage 5 — heal + native primitives + top-level pipeline.
    include("heal_test.jl")         # surface diagnostics (open/non-manifold/degenerate…)
    include("geometry_test.jl")     # box / cylinder / box-tunnel primitives
    include("cad_test.jl")          # native analytical geometry (surfaces + exact imprints), no OCC
    include("pipeline_test.jl")     # mesh_volume: validated-or-explicit-blocker

    # Stage 6 — high-order elements.
    include("highorder_test.jl")    # quadratic (P2) tet generation + type-11 I/O

    # Application: native meshes for all 22 HFSS User Guide case geometries (no gmsh/OCC).
    include("hfss_cases_test.jl")   # STATUS #12 meshing half — valid+watertight+conforming

    @testset "stage banner" begin
        @test Tessella.stage() isa Int
        @test Tessella.stage() == 6     # all package stages through P2/I/O are green
    end
end
