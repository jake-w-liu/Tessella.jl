using Test
using Tessella

@testset "Tessella" begin
    @testset "repository Julia-file layout" begin
        repository_root = normpath(joinpath(@__DIR__, ".."))
        function julia_files(relative_root)
            root = joinpath(repository_root, relative_root)
            paths = String[
                relpath(joinpath(directory, filename), root)
                for (directory, _, filenames) in walkdir(root)
                for filename in filenames if endswith(filename, ".jl")
            ]
            return sort!(paths)
        end

        source_files = julia_files("src")
        @test filter(path -> length(splitpath(path)) == 1, source_files) ==
              ["Tessella.jl"]
        source_domains = Set(("core", "fields", "geometry", "interfaces",
                              "meshing", "structured"))
        source_members = filter(!=("Tessella.jl"), source_files)
        @test all(length(splitpath(path)) == 2 for path in source_members)
        @test Set(first(splitpath(path)) for path in source_members) == source_domains

        test_files = julia_files("test")
        @test filter(path -> length(splitpath(path)) == 1, test_files) ==
              ["runtests.jl"]
        test_domains = Set(("core", "fields", "geometry", "integration",
                            "interfaces", "meshing", "structured"))
        test_members = filter(!=("runtests.jl"), test_files)
        @test all(length(splitpath(path)) == 2 for path in test_members)
        discovered_test_domains = Set(
            first(splitpath(path)) for path in test_members if
            first(splitpath(path)) != "tmp")
        @test discovered_test_domains == test_domains

        validation_files = julia_files("validation")
        @test filter(path -> length(splitpath(path)) == 1, validation_files) ==
              ["run_all.jl"]
        validation_members = filter(!=("run_all.jl"), validation_files)
        @test all(length(splitpath(path)) == 2 for path in validation_members)
        @test "support" in Set(first(splitpath(path)) for path in validation_members)
    end

    # Stage 0 — Foundations (CRC-gated, DEVELOPMENT.md discipline).
    include("core/predicates_test.jl")   # exact predicates vs exact-rational oracle
    include("core/meshtypes_test.jl")    # mesh container, topology, quality, checksum
    include("core/transform_test.jl")    # validated affine transforms + orientation preservation
    include("interfaces/io_test.jl")     # .msh v2/v4 round-trip, STL, .geo scan
    include("core/elements_test.jl")     # fixed/special Gmsh 4.15.2 records + mixed entity I/O
    include("meshing/recombine_test.jl") # deterministic triangle-to-quad recombination
    include("meshing/refine_test.jl")    # deterministic one-level uniform simplex refinement
    include("structured/transfinite_test.jl") # validated four-sided planar structured patches
    include("structured/transfinite_curve_test.jl") # straight Progression/Bump/Beta/HWall laws
    include("structured/transfinite_triangle_test.jl") # specific three-sided triangle/quad patches
    include("structured/transfinite_quad_test.jl") # recombined four-sided quadrangle patches
    include("structured/transfinite_volume_test.jl") # affine six-face structured volumes
    include("structured/transfinite_prism_test.jl") # affine five-face transfinite prisms
    include("structured/transfinite_hex_test.jl") # affine six-face recombined hexahedra
    include("geometry/model_test.jl")    # entity kernel + .geo execution
    include("geometry/geo_periodic_test.jl") # literal periodic-curve transforms
    include("geometry/model_periodic_io_test.jl") # classified periodic MSH projection
    include("geometry/nurbs_test.jl")    # De Boor vs Bernstein/circle oracles
    include("meshing/boundarylayer_test.jl") # prismatic layer extrusion
    include("meshing/periodic_test.jl")  # periodic identification
    include("interfaces/api_test.jl") # synchronized session, detached mesh cache
    include("interfaces/cli_test.jl") # bounded parser + non-destructive output
    include("interfaces/gui_test.jl") # validated headless command/state machine
    include("interfaces/post_test.jl") # owned scalar views + synchronized plugins

    # Stage 1 — 2-D meshing (CRC-gated).
    include("meshing/mesh2d_test.jl") # Delaunay: exact empty-circumcircle oracle

    # Stage 2 — Gmsh-compatible size fields + 1-D/surface meshing (CRC-gated).
    include("fields/sizefield_test.jl") # Distance/Threshold/Box/Min fields + local 3-D sizing
    include("meshing/mesh1d_test.jl") # size fields + graded edge meshing vs arc length
    include("meshing/meshsurface_test.jl") # planar/cylinder/parametric surface meshing

    # Stage 3 — 3-D meshing and certified boundary recovery (CRC-gated).
    include("meshing/mesh3d_test.jl") # 3-D Delaunay + volume filling + coax junction

    # Stage 4 — optimization, quality, and sliver removal (CRC-gated).
    include("meshing/optimize_test.jl") # tet quality report + Laplacian & ODT smoothing

    # Stage 5 — heal + native primitives + top-level pipeline.
    include("geometry/heal_test.jl") # surface diagnostics (open/non-manifold/degenerate…)
    include("geometry/geometry_test.jl") # box / cylinder / box-tunnel primitives
    include("geometry/cad_test.jl") # native analytical geometry (surfaces + exact imprints), no OCC
    include("geometry/brep_test.jl") # ISO-10303-21 STEP / IGES classified-solid import
    include("integration/pipeline_test.jl") # mesh_volume: validated-or-explicit-blocker

    # Stage 6 — high-order elements.
    include("meshing/highorder_test.jl") # quadratic (P2) tet generation + type-11 I/O

    # Application: native meshes for all 22 HFSS User Guide case geometries (no gmsh/OCC).
    include("integration/hfss_cases_test.jl") # STATUS #12 meshing half — valid+watertight+conforming

    @testset "stage banner" begin
        @test Tessella.stage() isa Int
        @test Tessella.stage() == 6     # all package stages through P2/I/O are green
    end
end
