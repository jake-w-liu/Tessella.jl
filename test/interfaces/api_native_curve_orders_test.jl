using Test
@testset "Native Curve conversion, generation, mutation and raw queries" begin
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","public_curve_order_probe_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","public_curve_generate_visible_probe_v3.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","native_closed_curve_endpoint_probe_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","native_closed_curve_mutation_probe_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","public_curve_discrete_mixed_v7_probe_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","public_curve_query_probe_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","native_ordinary_closed_curve_probe_v11.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","native_ordinary_curve_grading_probe_v3.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","native_curve_dense_geometry_probe_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","native_occ_center_circle_broad_probe_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","native_curve_point_storage_probe_v1.jl"))
    mktempdir() do output_dir
        withenv("TESSELLA_GEOMETRY_DERIVATIVES_ACTUAL"=>joinpath(output_dir,"derivatives.toml")) do
            include(joinpath(@__DIR__,"..","..","validation","native_curve_orders","native_curve_geometry_derivatives_probe_v1.jl"))
        end
    end
end
