# Native OCC coverage: literal OCC shapes/copied-pole registry gates and the
# curve namespace identity gates.
# (Merged api_occ_shape_registry_test.jl and api_occ_curve_namespace_test.jl.)

module NativeOCCShapeRegistryPermanentV1
include(joinpath(@__DIR__,"..","..","validation","occ_shape_registry","occ_native_zero_metadata_gate_v14.jl"))
include(joinpath(@__DIR__,"..","..","validation","occ_shape_registry","occ_shape_registry_extended_gate_v7.jl"))
include(joinpath(@__DIR__,"..","..","validation","occ_shape_registry","occ_native_zero_identity_gate_v11.jl"))
end

module OCCCurveNamespacePermanentV1
using Test,Tessella
const asset_dir=normpath(joinpath(@__DIR__,"..","..","validation","occ_curve_identity"))
for path in ("occ_copied_curve_gate_v3.jl","occ_namespace_identity_gate_v6.jl",
             "occ_curve_endpoint_compiled_v8.jl")
    include(joinpath(asset_dir,path))
end
end
