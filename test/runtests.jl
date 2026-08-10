using Test
using Tessella

@testset "Tessella" begin
    # Stage-0 skeleton smoke test. Real CRC suites arrive with each stage
    # (see DEVELOPMENT.md): predicates-vs-exact-rational, .msh round-trip vs gmsh,
    # Delaunay invariants, boundary-recovery on the enclosure coax junction, etc.
    @test Tessella.stage() == 0
end
