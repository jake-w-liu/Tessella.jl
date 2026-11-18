# Regression: mesh_planar used refine!'s fixed 300k-step budget, so any request
# finer than ~600k triangles failed with "reached maxsteps" (unit square with a
# hole, max_area=2e-6).  The budget now scales with bounding-box area / max_area.
using Test
using Tessella

@testset "mesh_planar refinement step budget scales with max_area" begin
    xs = [0.0, 1.0, 1.0, 0.0]; ys = [0.0, 0.0, 1.0, 1.0]
    @test Tessella._planar_refine_steps(xs, ys, Inf) == 300_000
    @test Tessella._planar_refine_steps(xs, ys, 1e-3) == 300_000
    @test Tessella._planar_refine_steps(xs, ys, 2e-6) == 2_000_000
    @test Tessella._planar_refine_steps(xs .* 1e200, ys .* 1e200, 1e-300) == 300_000  # overflowing estimate
    @test Tessella._planar_refine_steps(xs, ys, 1e-20) == 1_000_000_000_000
    m = mesh_planar(xs, ys, [(1, 2), (2, 3), (3, 4), (4, 1)]; max_area=1e-3)
    @test validate(m).ok
end
