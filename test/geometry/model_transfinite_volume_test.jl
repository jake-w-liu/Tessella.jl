# Transfinite-volume coverage: actual emission controls, the full native cube
# public controls, and retained boundary proposal bindings.
# (Merged model_transfinite_volume_{emission,full_native_cube,boundary_native}_test.jl.)

using Test,Tessella
include(joinpath(@__DIR__,"..","..","validation","native_geometry","transfinite_volume_actual_emission_controls_v5.jl"))
mktempdir() do directory
    TFVolumeActualEmissionControlsV5.run(joinpath(directory,"tf_emission"))
end

include(joinpath(@__DIR__,"..","..","validation","native_geometry","full_native_cube_controls_v1.jl"))
mktempdir() do directory
    TFExactCubeNativePublicControlsV1.run(joinpath(directory,"full_native_cube"))
end

include(joinpath(@__DIR__,"..","..","validation","native_geometry","transfinite_volume_actual_binding_controls_v1.jl"))
@testset "TF retained boundary proposal actual controls" begin
    mktempdir() do directory
        TFVolumeActualBoundaryControlsV1.run(joinpath(directory,"tf_boundary"))
    end
end
