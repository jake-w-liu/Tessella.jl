module APINativeQuadSplitTests
using Tessella,Test
const VALIDATION=joinpath(@__DIR__,"..","..","validation","native_geometry")
include(joinpath(VALIDATION,"native_quad_split_complete_v1.jl"))
include(joinpath(VALIDATION,"native_quad_split_carrier_controls_v1.jl"))
const OUTPUT_PREFIX=haskey(ENV,"TESSELLA_NATIVE_QUAD_OUTPUT_PREFIX") ?
    ENV["TESSELLA_NATIVE_QUAD_OUTPUT_PREFIX"] : joinpath(mktempdir(),"native_quad_split")
@testset "Complete native Quad split scopes" begin
NativeQuadSplitCompleteV1.run(OUTPUT_PREFIX*".functional.toml")
NativeQuadSplitCarrierControlsV1.run(OUTPUT_PREFIX*".carrier.toml")
end
end
