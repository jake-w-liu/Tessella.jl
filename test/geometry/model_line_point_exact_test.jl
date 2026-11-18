# The suite includes every test file into Main; this file keeps its aliases to itself.
module LinePointExactTests
using Test
using Random
using Tessella
import Tessella.Model as M

# A point on a straight Line is its exact rational value rounded once. The
# Float64 evaluation must return that value bit for bit wherever it answers,
# and decline where it cannot.

line_reference(a, b, t) =
    Float64(Rational{BigInt}(a) + Rational{BigInt}(t) * (Rational{BigInt}(b) - Rational{BigInt}(a)))

function line_coordinate(rng)
    kind = rand(rng, 1:7)
    kind == 1 && return rand(rng)
    kind == 2 && return randn(rng) * 10.0^rand(rng, -12:12)
    kind == 3 && return Float64(rand(rng, -8:8))
    kind == 4 && return rand(rng, (0.0, -0.0, 1.0, 0.5, 0.1, 0.2, 0.3))
    kind == 5 && return ldexp(1.0 + rand(rng, 0:7) / 8, rand(rng, -40:40))
    kind == 6 && return nextfloat(rand(rng, (1.0, 0.5, 2.0, 3.0)), rand(rng, -3:3))
    return rand(rng, (0.1, 0.2, 0.3, 0.7)) * rand(rng, (1, 3, 5, 7, 1e6, 1e-6))
end

function line_parameter(rng)
    kind = rand(rng, 1:4)
    kind == 1 && return rand(rng)
    kind == 2 && return rand(rng, 0:128) / 128          # the integration's dyadic points
    kind == 3 && return rand(rng, (0.0, 1.0, 0.5, 1 / 3, 2 / 3, 0.1))
    return nextfloat(rand(rng, 0:16) / 16, rand(rng, -2:2))
end

@testset "Line points are the rounded exact value" begin
    rng = MersenneTwister(20261118)
    answered = 0
    for _ in 1:20_000
        a = line_coordinate(rng); b = line_coordinate(rng); t = line_parameter(rng)
        value = M._model_exact_lerp(a, b, t)
        value === nothing && continue
        answered += 1
        expected = line_reference(a, b, t)
        value == expected || @test value == expected
    end
    @test answered > 15_000

    # Ties between two doubles are broken by the terms below them.
    @test M._model_exact_sum5(1.0, eps(1.0) / 2, floatmin(Float64), 0.0, 0.0) == nextfloat(1.0)
    @test M._model_exact_sum5(1.0, eps(1.0) / 2, -floatmin(Float64), 0.0, 0.0) == 1.0
    @test M._model_exact_sum5(1.0, eps(1.0) / 2, 0.0, 0.0, 0.0) == 1.0                  # even
    @test M._model_exact_sum5(nextfloat(1.0), eps(1.0) / 2, 0.0, 0.0, 0.0) == nextfloat(1.0, 2)
    @test M._model_exact_sum5(1e308, 1e308, 0.0, 0.0, 0.0) === nothing

    # Outside its range the evaluation declines and the rational path answers.
    @test M._model_exact_lerp(1e300, 0.75, 0.9453125) === nothing        # terms 2^990 apart
    @test M._model_exact_lerp(1e-200, 1e-200, 1e-200) === nothing        # product error underflows
    @test M._model_exact_lerp(0.0, 0.0, 0.5) == 0.0

    m = GeoModel()
    add_point!(m, 0.1, 0.2, 0.3; tag=1)
    add_point!(m, 0.7, -1.3, 2.9; tag=2)
    add_point!(m, 1e300, 0.2, 0.3; tag=3)
    add_line!(m, 1, 2; tag=1)
    add_line!(m, 1, 3; tag=2)
    for t in (0.0, 1.0, 0.5, 0.3, 1 / 3, 0.9453125)
        for curve in (1, 2)
            line = M._model_line_geometry(m, curve, "test")
            @test M._model_curve_point(m, curve, t, "test") ===
                  M._model_line_point(line, t, "test", 1)
        end
    end
    @test M._model_line_point_fast(m, 1, 0.3) !== nothing
    @test M._model_line_point_fast(m, 2, 0.9453125) === nothing
end
end
