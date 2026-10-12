# Default checks avoid optimizing the large production dispatchers. Compiler
# coverage is opt-in with TESSELLA_TEST_COMPILER=optimized (see README.md).
let mode=get(ENV,"TESSELLA_TEST_COMPILER","minimal")
    mode in ("minimal","optimized") || error(
        "TESSELLA_TEST_COMPILER must be minimal or optimized")
    options=Base.JLOptions()
    coverage_file=options.output_code_coverage==C_NULL ? "" : unsafe_string(options.output_code_coverage)
    # Base.julia_cmd deliberately omits fixed coverage filenames in children.
    fixed_coverage=!isempty(coverage_file) && !occursin("%p",coverage_file)
    if mode=="minimal" && options.compile_enabled!=3 && !fixed_coverage
        project=Base.active_project()
        project===nothing && error("Run functional tests with --project=.")
        command=`$(Base.julia_cmd()) --project=$(dirname(project)) --startup-file=no --check-bounds=yes --compile=min -O0 $(abspath(@__FILE__)) $ARGS`
        println("Running small functional checks with minimal compilation.")
        flush(stdout)
        run(command)
        exit(0)
    end
end

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
    include("functional/brep_assembly_smoke.jl")
    include("functional/brep_mapped_occurrences_smoke.jl")
    include("functional/brep_nested_relationship_smoke.jl")
    include("functional/brep_mixed_relationship_smoke.jl")
    include("functional/brep_unit_conversion_smoke.jl")
    include("functional/brep_sphere_smoke.jl")
    include("functional/brep_cylinder_smoke.jl")
    include("functional/brep_cone_smoke.jl")
    include("functional/brep_torus_segment_smoke.jl")
    include("functional/model_view_export.jl")
    include("functional/model_export_smoke.jl")
    include("functional/model_import_smoke.jl")
    include("functional/model_merge_smoke.jl")
    include("functional/pos_import_smoke.jl")
    include("functional/view_time_smoke.jl")
    include("functional/view_representation_smoke.jl")
    include("functional/cli_information.jl")
    include("functional/cli_check_smoke.jl")
    include("functional/vtk_complete_hex.jl")
    include("functional/vtk_lagrange_curve.jl")
    include("functional/vtk_lagrange_smoke.jl")
    include("functional/vtk_pyramid13.jl")
    include("functional/brep_planar_faceholes_smoke.jl")
    include("functional/constant_metric_box.jl")
    include("functional/constant_metric_regions.jl")
    include("functional/pyramid_option.jl")
    include("functional/brep_bounded143_smoke.jl")
    include("functional/brep_plane190_uv_smoke.jl")
    include("functional/brep_plane190_transform_smoke.jl")
    include("functional/brep_plane190_chain_smoke.jl")
    include("functional/partition_producer.jl")
    include("functional/partition_brep.jl")
    include("functional/rigid_periodic_smoke.jl")
    include("functional/occ_circle_smoke.jl")
    include("functional/occ_plane_smoke.jl")
    include("functional/occ_ellipse_smoke.jl")
    include("functional/occ_volume_smoke.jl")
    include("functional/occ_primitives_smoke.jl")
    include("functional/occ_transform_smoke.jl")
end
