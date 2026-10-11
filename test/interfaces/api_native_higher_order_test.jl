# Native higher-order coverage: factories/public identities, queries and raw
# identities, cached queries, quality algorithms and retained public state,
# and legacy-composition query authority.
# (Merged api_native_higher_{order,queries,cached_queries,quality,
# legacy_composition}_test.jl; embedded actual_ast/resources probes dropped —
# the global allocation audit covers that surface.)

using Tessella, Test

@testset "Native higher-order factories and public identities" begin
    fixture_dir=joinpath(@__DIR__,"..","..","validation","higher_order_native")
    for name in ("raw_sol1","lifecycle","native_plane","generation")
        include(joinpath(fixture_dir,"higher_order_"*name*"_probe_v1.jl"))
    end
    include(joinpath(fixture_dir,"higher_order_metadata_probe_v1.jl"))
    include(joinpath(fixture_dir,"higher_raw_volume_generation_probe_v1.jl"))
end

@testset "Native higher-order queries and original raw identities" begin
    fixtures=joinpath(@__DIR__,"..","..","validation","higher_order_native")
    for name in ("higher_order_query_probe_v2.jl","higher_order_bare_query_probe_v1.jl","higher_order_tet_measure_probe_v1.jl","higher_legacy_p2_native_payload_probe_v1.jl","higher_legacy_p2_determinant_probe_v1.jl")
        include(joinpath(fixtures,name))
    end
end

@testset "Original N records beside a generated cache" begin
    include(joinpath(@__DIR__,"..","..","validation","higher_order_native",
        "higher_order_cached_query_probe_v2.jl"))
end

@testset "Native higher-order quality algorithms and retained public state" begin
    validation=joinpath(@__DIR__,"..","..","validation","higher_order_native")
    for name in ("higher_native_bounds_probe_v1.jl","higher_native_bounds_components_v4.jl",
                 "higher_native_remaining_quality_probe_v5.jl","higher_native_quality_edge_probe_v5.jl")
        include(joinpath(validation,name))
    end
end

@testset "Legacy cache and independent native records retain complete query authority" begin
    include(joinpath(@__DIR__,"..","..","validation","higher_order_native",
        "higher_legacy_cached_query_probe_v2.jl"))
    include(joinpath(@__DIR__,"..","..","validation","higher_order_native",
        "higher_legacy_composition_supports_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","higher_order_native",
        "higher_legacy_metadata_numeric_witness_v5.jl"))
end
