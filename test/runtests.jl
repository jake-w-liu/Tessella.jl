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
    for path in (
        "meshing/ordered_matching_seeded_test.jl",
        "core/transform_axis_scale_test.jl",
        "fields/distance_triangle_translation_test.jl",
        "geometry/model_plane_closed_curve_test.jl",
        "geometry/model_compact_atomic_coverage_test.jl",
        "functional/list_view_export.jl",
        "functional/brep_assembly_smoke.jl",
        "functional/brep_mapped_occurrences_smoke.jl",
        "functional/brep_nested_relationship_smoke.jl",
        "functional/brep_mixed_relationship_smoke.jl",
        "functional/brep_unit_conversion_smoke.jl",
        "functional/brep_sphere_smoke.jl",
        "functional/brep_cylinder_smoke.jl",
        "functional/brep_cone_smoke.jl",
        "functional/brep_torus_segment_smoke.jl",
        "functional/model_view_export.jl",
        "functional/model_export_smoke.jl",
        "functional/model_import_smoke.jl",
        "functional/model_merge_smoke.jl",
        "functional/pos_import_smoke.jl",
        "functional/view_time_smoke.jl",
        "functional/view_representation_smoke.jl",
        "functional/cli_information.jl",
        "functional/cli_check_smoke.jl",
        "functional/vtk_complete_hex.jl",
        "functional/vtk_lagrange_curve.jl",
        "functional/vtk_lagrange_smoke.jl",
        "functional/vtk_pyramid13.jl",
        "functional/brep_planar_faceholes_smoke.jl",
        "functional/brep_curved_trims_smoke.jl",
        "functional/brep_nurbs_trim_smoke.jl",
        "functional/brep_nurbs_mixed_smoke.jl",
        "functional/brep_ellipse_trims_smoke.jl",
        "functional/constant_metric_box.jl",
        "functional/constant_metric_regions.jl",
        "functional/pyramid_option.jl",
        "functional/brep_bounded143_smoke.jl",
        "functional/brep_plane190_uv_smoke.jl",
        "functional/brep_plane190_transform_smoke.jl",
        "functional/brep_plane190_chain_smoke.jl",
        "functional/partition_producer.jl",
        "functional/partition_brep.jl",
        "functional/rigid_periodic_smoke.jl",
        "functional/occ_circle_smoke.jl",
        "functional/occ_plane_smoke.jl",
        "functional/occ_ellipse_smoke.jl",
        "functional/occ_volume_smoke.jl",
        "functional/occ_primitives_smoke.jl",
        "functional/occ_transform_smoke.jl",
        "functional/occ_queries_smoke.jl",
        "functional/occ_measures_smoke.jl",
        "functional/occ_wedge_smoke.jl",
        "functional/occ_inertia_smoke.jl",
        "functional/occ_spline_smoke.jl",
        "functional/occ_closest_smoke.jl",
        "functional/occ_sweep_smoke.jl",
        "functional/occ_revolve_smoke.jl",
        "functional/occ_corner_smoke.jl",
        "functional/occ_offset_smoke.jl",
        "functional/occ_nurbs_smoke.jl",
        "functional/occ_boolean_smoke.jl",
        "functional/api_surface_smoke.jl",
        "interfaces/api_model_postview_field_test.jl",
        "interfaces/api_plugin_model_probe_test.jl",
    )
        println("Testing ", path, "...")
        flush(stdout)
        elapsed = @elapsed include(path)
        println("Completed ", path, " in ", round(elapsed; digits=2), " s")
        flush(stdout)
    end
end
