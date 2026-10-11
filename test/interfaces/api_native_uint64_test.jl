# Native UInt64 tag coverage: quad birth, post-wrap generation, counter wrap
# lifecycle/builder and mesh-tag producers/publication/cached mutations.
# (Merged from api_native_{quad_uint64_birth,uint64_post_wrap_generation,
# uint64_counter_wrap,uint64_mesh_tags}_test.jl; AST-structure and byte-budget
# probes dropped — those are covered by the global allocation audit.)

module APINativeQuadUInt64BirthTestV1
using Tessella,Test
include("../../validation/native_geometry/native_quad_uint64_birth_controls_v1.jl")
prefix=get(ENV,"TESSELLA_NATIVE_QUAD_UINT64_OUTPUT_PREFIX",joinpath(mktempdir(),"quad_uint64"))
NativeQuadUInt64BirthV1.run(prefix*".functional.toml")
end

module APINativeUInt64PostWrapGenerationTests
using Tessella,Test
const VALIDATION=joinpath(@__DIR__,"..","..","validation","native_geometry")
include(joinpath(VALIDATION,"native_uint64_post_wrap_birth_v2.jl"))
include(joinpath(VALIDATION,"native_uint64_post_wrap_generate_v1.jl"))
const PREFIX=get(ENV,"TESSELLA_NATIVE_NODE_ERROR_OUTPUT_PREFIX",joinpath(mktempdir(),"native_post_wrap"))
NativeUInt64PostWrapBirthProbeV2.run(PREFIX*".birth.toml")
NativeUInt64PostWrapGenerateProbeV1.run(PREFIX*".no_birth.toml")
end

module APINativeUInt64CounterWrapTests
using Tessella,Test
const VALIDATION=joinpath(@__DIR__,"..","..","validation","native_geometry")
include(joinpath(VALIDATION,"native_uint64_counter_wrap_v1.jl"))
include(joinpath(VALIDATION,"native_uint64_counter_wrap_lifecycle_v1.jl"))
include(joinpath(VALIDATION,"native_uint64_counter_wrap_builder_v1.jl"))
const PREFIX=get(ENV,"TESSELLA_NATIVE_UINT64_WRAP_OUTPUT_PREFIX",joinpath(mktempdir(),"native_wrap"))
@testset "Defined native UInt64 wrap, model lifetime and detached builder" begin
    NativeUInt64WrapProbeV1.run(PREFIX*".wrap.toml")
    NativeUInt64WrapLifecycleProbeV1.run(PREFIX*".lifecycle.toml")
    NativeUInt64WrapBuilderV1.run(PREFIX*".builder.toml")
end
include("../../validation/native_geometry/native_uint64_node_prevalidation_v1.jl")
NativeUInt64NodePrevalidationProbeV1.run(PREFIX*".prevalidation.toml")
end

module APINativeUInt64MeshTagsTests
using Tessella,Test
const VALIDATION=joinpath(@__DIR__,"..","..","validation","native_geometry")
include(joinpath(VALIDATION,"native_uint64_producer_v1.jl"))
include(joinpath(VALIDATION,"native_uint64_partial_publication_v1.jl"))
include(joinpath(VALIDATION,"native_uint64_cached_mutations_v1.jl"))
include(joinpath(VALIDATION,"native_uint64_carrier_controls_v1.jl"))
include(joinpath(VALIDATION,"native_lookup_large_tag_v1.jl"))
const PREFIX=get(ENV,"TESSELLA_NATIVE_UINT64_OUTPUT_PREFIX",joinpath(mktempdir(),"native_uint64"))
@testset "Native UInt64 mesh producers, zero identity, partial publication and cached mutation" begin
    NativeUInt64ProducerProbeV4.run(Tessella.API,PREFIX*".producer.toml")
    NativeUInt64PartialErrorProbeV4.run(Tessella.API,PREFIX*".partial.toml")
    NativeUInt64CachedMutationProbeV3.run(Tessella.API,PREFIX*".cached.toml")
    NativeUInt64CarrierControlsV1.run(PREFIX*".carrier.toml")
    NativeLookupLargeElementTagsProbeV2.run(PREFIX*".query16.toml")
end
end
