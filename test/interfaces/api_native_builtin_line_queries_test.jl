using Test
@testset "Native builtin Line inverse and retained query state" begin
    directory=joinpath(@__DIR__,"..","..","validation","native_geometry")
    include(joinpath(directory,"native_builtin_line_queries_probe_v2.jl"))
end
