using Test,Tessella

Base.JLOptions().compile_enabled==3 && error(
    "Coverage requires normal compilation; omit --compile=min.")
Base.JLOptions().code_coverage==0 && error(
    "Run coverage tests with --code-coverage=user.")

# Separate modules keep each standalone test's helpers and imports isolated.
@testset "Compiled native kernel coverage" begin
    for (index,path) in enumerate((
        "core/predicates_test.jl",
        "core/meshtypes_test.jl",
        "core/mesh_entity_topology_test.jl",
        "core/mesh_point_location_test.jl",
        "core/mesh_reference_geometry_test.jl",
        "core/mesh_quadrature_test.jl",
        "core/mesh_element_quality_test.jl",
        "core/transform_test.jl",
        "core/allocation_audit_test.jl",
        "meshing/mesh1d_test.jl",
        "meshing/mesh2d_test.jl",
        "meshing/refine_test.jl",
        "meshing/periodic_test.jl",
        "coverage/native_geometry.jl",
        "structured/transfinite_curve_test.jl",
        "geometry/nurbs_test.jl",
        "geometry/heal_test.jl",
        "fields/scalar_finite_interpolation_test.jl",
        "fields/metric_spectral_range_test.jl",
        "fields/distance_triangle_translation_test.jl",
        "fields/finite_difference_range_test.jl",
        "fields/restrict_anisotropic_test.jl",
        "fields/metric_spectrum_certificate_test.jl",
        "fields/structured_extend_lonlat_range_test.jl",
        "fields/axial_field_geometry_test.jl",
        "fields/postview_translated_line_test.jl",
        "fields/postview_compact_tolerance_test.jl",
        "interfaces/api_hybrid_winslow_allocation_test.jl",
        "interfaces/io_test.jl",
    ))
        println("Coverage tests: ",path);flush(stdout)
        scope=Core.eval(Main,:(module $(Symbol("KernelCoverage",index)) end))
        Base.include_string(scope,"using Test,LinearAlgebra,Tessella; import Tessella.MeshTypes as MT")
        Base.include(scope,joinpath(@__DIR__,"..",path))
    end
    module_io=Core.eval(Main,:(module IOCoverage end))
    Base.include_string(module_io,"using Test,Tessella; import Tessella.MeshTypes as MT")
    Base.include(module_io,joinpath(@__DIR__,"..","functional","fixtures.jl"))
    Base.include(module_io,joinpath(@__DIR__,"..","functional","io.jl"))
end

# Coverage counters are flushed only when this Julia process exits.
report_dir=joinpath(@__DIR__,"..","..","coverage")
mkpath(report_dir)
write(joinpath(report_dir,"runner-pid"),string(getpid()))
