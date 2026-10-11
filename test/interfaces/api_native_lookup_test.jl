# Native lookup coverage: lifetime, checked tags, raw Tri3 geometry and the
# abort-policy controls. (Merged api_native_lookup_abort_policy_test.jl; the
# AST-structure probe was dropped — global allocation audit covers that.)

module NativeLookupPermanentTestV1
using Test,Tessella
const SUPPORT=normpath(joinpath(@__DIR__,"..","..","validation","native_geometry"))
include(joinpath(SUPPORT,"native_lookup_large_tag_v1.jl"))
include(joinpath(SUPPORT,"native_lookup_lifecycle_v1.jl"))
include(joinpath(SUPPORT,"native_tri3_mapped_point_v1.jl"))

function run(prefix)
    @testset "Native lookup lifetime, checked tags and raw Tri3 geometry" begin
        NativeLookupLargeElementTagsProbeV2.run(prefix*".large_tags.toml")
        NativeLookupLifecycleProbeV2.run(prefix*".lookup.toml")
        NativeTri3MappedPointProbeV1.run(prefix*".tri3.toml")
    end
end

prefix=get(ENV,"TESSELLA_NATIVE_LOOKUP_OUTPUT_PREFIX","")
if isempty(prefix)
    mktempdir() do directory
        run(joinpath(directory,"native_lookup"))
    end
else
    run(prefix)
end
end

module APINativeLookupAbortPolicyTestV1
include("../../validation/native_geometry/native_lookup_abort_policy_controls_v1.jl")
include("../../validation/native_geometry/native_lookup_abort_neighbors_controls_v1.jl")
prefix=get(ENV,"TESSELLA_NATIVE_LOOKUP_ABORT_OUTPUT_PREFIX",joinpath(mktempdir(),"lookup_abort"))
NativeLookupAbortPolicyV1.run(prefix)
NativeLookupAbortNeighborsV1.run(prefix)
include("../../validation/native_geometry/native_lookup_abort_fractional_controls_v1.jl")
NativeLookupAbortFractionalV1.run(prefix)
end
