# The suite includes every test file into Main; this file keeps its aliases to itself.
module PcurveSplineLocalTests
using Test
using Random
using Tessella
import Tessella.Model as M

# A spline pcurve is evaluated on the poles of one knot span. The values must
# be those of the full-row evaluation of the explicitly derived splines, bit
# for bit: meshes of Boolean section curves depend on them.

function pcurve_reference(knots, p, cps, t)
    point = M._occ_spl2d_eval_dense(knots, p, cps, t)
    k1, p1, cps1 = M._occ_spl2d_derive(knots, p, cps)
    first = M._occ_spl2d_eval_dense(k1, p1, cps1, t)
    second = (0.0, 0.0)
    if p1 > 0 && length(cps1) > 1
        k2, p2, cps2 = M._occ_spl2d_derive(k1, p1, cps1)
        second = M._occ_spl2d_eval_dense(k2, p2, cps2, t)
    end
    return point, first, second
end

@testset "Spline pcurves evaluate locally with the dense values" begin
    rng = MersenneTwister(20261118)
    checked = 0
    for _ in 1:300
        p = rand(rng, 1:5)
        n = rand(rng, p+1:p+30)
        interior = sort!(rand(rng, n - p - 1))
        for k in 2:length(interior)
            rand(rng) < 0.2 && (interior[k] = interior[k-1])      # repeated knots
        end
        scale = 10.0^rand(rng, -3:3)
        knots = vcat(fill(0.0, p + 1), interior, fill(1.0, p + 1)) .* scale .+ rand(rng, (-5.0, 0.0, 100.0))
        cps = [(3randn(rng), 3randn(rng)) for _ in 1:n]
        pc = (spl2d=(knots, p, cps),)
        span = knots[end] - knots[1]
        samples = vcat(knots[1], knots[end], knots[1] - 0.1span, knots[end] + 0.1span,
                       rand(rng, knots, 4), knots[1] .+ span .* rand(rng, 12))
        for t in samples
            expected = pcurve_reference(knots, p, cps, t)
            @test isequal(M._occ_pcurve_d2(pc, t), expected)
            @test isequal(M._occ_pcurve_d1(pc, t), (expected[1], expected[2]))
            @test isequal(M._occ_spl2d_eval(knots, p, cps, t), expected[1])
            checked += 1
        end
    end
    @test checked == 300 * 20

    # Degrees beyond the stack tables keep the dense path.
    p = 9
    knots = vcat(fill(0.0, p + 1), [0.25, 0.5, 0.75], fill(1.0, p + 1))
    cps = [(randn(rng), randn(rng)) for _ in 1:length(knots)-p-1]
    @test isequal(M._occ_pcurve_d2((spl2d=(knots, p, cps),), 0.4), pcurve_reference(knots, p, cps, 0.4))
end
end
