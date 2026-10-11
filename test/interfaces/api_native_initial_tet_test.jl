module NativeInitialTetProductionTests
using Tessella,Test
const M=Tessella.Model
const A=Tessella.API
const K=M.ModelNativeInitialTet
const NATIVE_INITIAL_TET_VALIDATION=joinpath(@__DIR__,"..","..","validation","native_initial_tet")
include(joinpath(NATIVE_INITIAL_TET_VALIDATION,"fixtures_v1.jl"))
include(joinpath(NATIVE_INITIAL_TET_VALIDATION,"extended_fixtures_v2.jl"))
include(joinpath(NATIVE_INITIAL_TET_VALIDATION,"kernel_probe_v3.jl"))
include(joinpath(NATIVE_INITIAL_TET_VALIDATION,"certificates_v2.jl"))
include(joinpath(NATIVE_INITIAL_TET_VALIDATION,"public_production_v7.jl"))
end
module NativeEmbeddedTetProductionTests
include(joinpath(@__DIR__,"..","..","validation","native_initial_tet","embedded_public_v1.jl"))
end
module NativeInitialTetQualityTests
include(joinpath(@__DIR__,"..","..","validation","native_initial_tet","quality32_controls_v1.jl"))
end
