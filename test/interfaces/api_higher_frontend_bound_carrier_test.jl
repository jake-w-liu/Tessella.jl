module HigherFrontendBoundCarrierTestsV1
using Test, Tessella
include(joinpath(@__DIR__,"..","..","validation","native_geometry",
    "higher_frontend_bound_carrier_branches_v1.jl"))
include(joinpath(@__DIR__,"..","..","validation","native_geometry",
    "higher_frontend_bound_carrier_scene_v1.jl"))
mktempdir() do directory
    HigherFrontendBoundCarrierBranchesV1.run(joinpath(directory,"branches.toml"))
    HigherFrontendBoundCarrierSceneV1.run(directory)
end
end
