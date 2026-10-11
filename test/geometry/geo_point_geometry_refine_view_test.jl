using Test,Tessella
include(joinpath(@__DIR__,"..","..","validation","native_geometry","point_geometry_refine_controls_v1.jl"))
mktempdir() do directory
    GeometryPointRefineControlsV1.run(joinpath(directory,"point_refine"))
end
