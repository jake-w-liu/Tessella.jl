using Test
using Tessella

@testset "Tessella" begin
    # Stage 0 — Foundations (CRC-gated, DEVELOPMENT.md discipline).
    include("predicates_test.jl")   # exact predicates vs exact-rational oracle
    include("meshtypes_test.jl")    # mesh container, topology, quality, checksum
    include("io_test.jl")           # .msh v2/v4 round-trip, STL, .geo scan

    @testset "stage banner" begin
        @test Tessella.stage() isa Int
        @test Tessella.stage() >= 1     # Stage 0 gate is green
    end
end
