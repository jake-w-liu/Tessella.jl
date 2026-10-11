# Native NURBS coverage: cached surface projection, distant initial-state
# probe, and tensor births / mixed-kernel context corrections.
# (Merged api_native_nurbs_{projection,projection_distant,surface}_test.jl;
# the NURBS ast/resources probes were dropped — global audit covers them.)

module NativeNURBSProjectionPermanentV1
using Test
const _projection_validation=joinpath(@__DIR__,"..","..","validation","nurbs_surface")
@testset "Native cached Surface projection and original root arithmetic" begin
    include(joinpath(_projection_validation,"projection_public_numeric_v1.jl"))
    include(joinpath(_projection_validation,"closest48_v4.jl"))
    include(joinpath(_projection_validation,"projection_cache_v2.jl"))
    include(joinpath(_projection_validation,"projection_root_v2.jl"))
end
end

using Test,Tessella
include(joinpath(@__DIR__,"..","..","validation","native_nurbs_projection","initial_state_probe_v1.jl"))

@testset "Native NURBS tensor births and mixed-kernel context correction" begin
    include(joinpath(@__DIR__,"..","..","validation","nurbs_surface","first_slice_v5.jl"))
    include(joinpath(@__DIR__,"..","..","validation","nurbs_surface","lifecycle_v5.jl"))
    include(joinpath(@__DIR__,"..","..","validation","nurbs_surface","curved_boundary_v12.jl"))
    include(joinpath(@__DIR__,"..","..","validation","nurbs_surface","query_lifecycle_v5.jl"))
    include(joinpath(@__DIR__,"..","..","validation","nurbs_surface","constructor_admission_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","nurbs_surface","mesh_routes_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","nurbs_surface","coincident_identity_v1.jl"))
    include(joinpath(@__DIR__,"..","..","validation","nurbs_surface","tf_context_v8.jl"))
end
