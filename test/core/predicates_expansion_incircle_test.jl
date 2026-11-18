# incircle/insphere: inputs the A-stage filter cannot certify are evaluated exactly
# with Float64 expansions (previously always Rational{BigInt}, ~2 µs / ~6-9 µs and
# thousands of bytes per call, which dominated Delaunay time on cospherical input).
# Signs must equal the independent rational oracle; the exact-zero path must not
# allocate.
using Test
using Tessella.Predicates
import Tessella.Predicates as P
include("oracles.jl")
using .Oracles

@noinline _ic_alloc(a, b, c, d) = @allocated incircle(a, b, c, d)
@noinline _is_alloc(a, b, c, d, e) = @allocated insphere(a, b, c, d, e)

@testset "incircle/insphere exact expansion stage" begin
    s = UInt64(0x5DEECE66D)
    rnd() = (s ⊻= s << 13; s ⊻= s >> 7; s ⊻= s << 17; s)
    function coord()
        r = rnd() % 6
        r == 0 && return Float64(Int(rnd() % 9) - 4)
        r == 1 && return (Int(rnd() % 9) - 4) / 10          # inexact differences
        r == 2 && return cospi(Int(rnd() % 24) / 12)
        r == 3 && return nextfloat(1.0, Int(rnd() % 7) - 3)
        r == 4 && return ldexp(Float64(rnd() >> 11) / 2.0^53, Int(rnd() % 120) - 60)
        return Float64(rnd() >> 11) / 2.0^53
    end
    sph = ((1, 2, 2), (2, 1, -2), (-2, 2, 1), (3, 0, 0), (0, -3, 0), (0, 0, 3), (-1, -2, 2))
    mism = 0; zeros = 0; nexp = 0
    for k in 1:250
        if isodd(k)
            c = (coord(), coord(), coord()); h = ldexp(1.0, Int(rnd() % 7) - 3)
            pts = ntuple(i -> ntuple(j -> c[j] + h * sph[mod1(i + k, 7)][j], 3), 5)
        else
            pts = ntuple(_ -> (coord(), coord(), coord()), 5)
        end
        flat = Tuple(x for p in pts for x in p)
        P._insphere_expansion(flat...) === nothing || (nexp += 1)
        w = oracle_insphere(pts...)
        w == 0 && (zeros += 1)
        insphere(pts...) == w || (mism += 1)
        q = ntuple(i -> (pts[i][1], pts[i][2]), 4)
        incircle(q...) == oracle_incircle(q...) || (mism += 1)
        q = ((pts[1][1], 0.0), (0.0, pts[1][1]), (-pts[1][1], 0.0), (pts[2][3], pts[3][3]))
        incircle(q...) == oracle_incircle(q...) || (mism += 1)
    end
    @test mism == 0
    @test zeros > 10
    @test nexp > 200
    # exact zero decided by expansions, without allocation
    a = (0.0, 0.0, 0.0); b = (1.0, 0.0, 0.0); c = (0.0, 1.0, 0.0); d = (1.0, 1.0, 0.0); e = (0.0, 0.0, 1.0)
    @test insphere(a, b, c, d, e) == 0
    @test _is_alloc(a, b, c, d, e) == 0
    @test incircle((1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)) == 0
    @test _ic_alloc((1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)) == 0
end
