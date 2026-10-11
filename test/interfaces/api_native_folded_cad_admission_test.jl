using Test, Tessella
@testset "Native folded CAD placement and SOL1 admission" begin
    include(joinpath(@__DIR__, "..", "..", "validation", "native_folded_cad", "native_folded_cad_probe_v1.jl"))
    include(joinpath(@__DIR__, "..", "..", "validation", "native_folded_cad", "native_folded_cad_receipt_audit_v1.jl"))
    include(joinpath(@__DIR__, "..", "..", "validation", "native_folded_cad", "native_folded_curve_sol1_probe_v1.jl"))
    include(joinpath(@__DIR__, "..", "..", "validation", "native_folded_cad", "native_folded_surface_sol1_probe_v1.jl"))
end
