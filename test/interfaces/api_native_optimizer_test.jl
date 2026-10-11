# Native optimizer coverage: ring/cube/generation kernels, model options,
# factory-UID rounds and the legacy P2 optimizer.
# (Merged api_native_factory_optimizer_uid_test.jl and
# api_native_legacy_p2_optimizer_test.jl; the resource-byte probes were
# dropped — global allocation audit covers that surface.)

module NativeOptimizerRingTests
include(joinpath(@__DIR__,"..","..","validation","native_optimizer","ring_kernel_v1.jl"))
end
module NativeOptimizerCubeTests
using Tessella,Test
const K=Tessella.Mesh3D
include(joinpath(@__DIR__,"..","..","validation","native_optimizer","cube_fixtures_v1.jl"))
include(joinpath(@__DIR__,"..","..","validation","native_optimizer","cube_controls_v1.jl"))
end
module NativeOptimizerGenerationTests
include(joinpath(@__DIR__,"..","..","validation","native_optimizer","tetrahedron_generation_v1.jl"))
end

module NativeOptimizerModelOptionsTests
include(joinpath(@__DIR__,"..","..","validation","native_optimizer","model_options_v4.jl"))
end

module NativeFactoryOptimizerUidR1Tests
include(joinpath(@__DIR__,"..","..","validation","native_optimizer","factory_uid_r1_controls_v1.jl"))
end

module NativeFactoryOptimizerUidR0Tests
include(joinpath(@__DIR__,"..","..","validation","native_optimizer","factory_uid_r0_controls_v1.jl"))
end

module APINativeLegacyP2OptimizerTests
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_legacy_p2_optimizer_v1.jl"))
prefix=get(ENV,"TESSELLA_NATIVE_LEGACY_P2_OUTPUT_PREFIX",joinpath(mktempdir(),"native_legacy_p2"))
NativeLegacyP2OptimizerV1.run(prefix)
end
