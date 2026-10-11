# Native public Relocate3D controls plus the P2 Msg-lifetime probe.
# (Merged api_native_public_relocate3d_p2_msg_test.jl; the byte-budget and
# AST probes that used to accompany these were dropped — the global
# allocation audit covers that surface.)

module APINativePublicRelocate3DTests
using Tessella,Test
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_public_relocate3d_v1.jl"))
const PREFIX=get(ENV,"TESSELLA_NATIVE_RELOCATE3D_OUTPUT_PREFIX",joinpath(mktempdir(),"native_relocate3d"))
NativePublicRelocate3DV1.run(PREFIX*".relocate3d.toml")
end

include(joinpath(@__DIR__,"..","..","validation","native_geometry",
    "native_public_relocate3d_p2_msg_v1.jl"))
NativePublicRelocate3DP2MsgV1.run(get(ENV,"TESSELLA_NATIVE_RELOCATE3D_OUT",tempname()))
