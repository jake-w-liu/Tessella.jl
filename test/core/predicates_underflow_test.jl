# Regression: the Shewchuk A-stage filters assume no underflow.  A partial product
# that underflows is later multiplied by a large factor (a z-difference or a lift),
# so its lost bits can dominate a determinant whose computed value and permanent
# are comfortably normal; the old filters then certified the wrong sign.
using Test
using Tessella.Predicates
include("oracles.jl")
using .Oracles

@testset "Predicates — gradual-underflow filter soundness" begin
    # orient3: true det = adz*(bdx*cdy) - cdz*(adx*bdy) ≈ 1e-100 - 1e-210 > 0, but
    # bdx*cdy = 1e-400 underflows to zero and the filter saw only -1e-210.
    a = (1e-50, 0.0, 1e300); b = (1e-200, -1e-60, 0.0)
    c = (0.0, 1e-200, 1e-100); d = (0.0, 0.0, 0.0)
    @test oracle_orient3(a, b, c, d) == 1
    @test orient3(a, b, c, d) == 1
    @test orient3_sos(a, b, c, d, 1, 2, 3, 4) == 1
    # Same configuration found by randomized fuzzing (old result -1).
    a = (7.578713009293459e-183, -8.373544579578191e119, 1.5439778652904985e270)
    b = (-1.650340909568883e-308, 5.882002648422471e-32, -1.3170839363e-313)
    c = (-1.287315065638103e-301, 1.34855715853e-313, 1.819828852943443e-152)
    d = (9.966e-320, -8.7220304e-314, 4.741035008215046e-10)
    @test orient3(a, b, c, d) == oracle_orient3(a, b, c, d) == 1

    # incircle: a lift ~1e300 multiplies an underflowed 2x2 minor (old result -1).
    P = ((1.2377063243218752e150, -7.931815252332647e-272),
         (1.3686587225979084e-301, -1.3691098009441343e-181),
         (2.695891806360233e-181, 2.7014285000421523e-121),
         (0.0, 6.651341422178277e-293))
    @test incircle(P...) == oracle_incircle(P...) == 1
    @test incircle_sos(P..., 1, 2, 3, 4) == 1
    P = ((5.491523538419388e-182, 1.255e-319), (1.896647119891547e-294, -2.7908418106333283e90),
         (1.2119218833290388e-69, -3.46011255169221e-302), (-0.0, -1.1862806016635119e-293))
    @test incircle(P...) == oracle_incircle(P...) == 1
    P = ((-2.874627766125948e144, -4.2886552031479354e144), (-5.708522375146372e-70, 2.3162330977017516e-272),
         (-9.066816376580749e-91, 6.61254305740444e-293), (-3.3425367472600816e-91, 5.343594551819552e-302))
    @test incircle(P...) == oracle_incircle(P...) == -1

    # Fixed-seed mixed-exponent sweep (coordinates spanning 2^-1000..2^1000, with
    # zeros): every sign must match the independent oracle.
    s = UInt64(0x9E3779B97F4A7C15)
    function nextcoord()
        s ⊻= s << 13; s ⊻= s >> 7; s ⊻= s << 17
        r = Int(s % 7)
        r == 0 && return 0.0
        m = (Float64(s >> 11) / 2.0^53) - 0.5
        e = Int((s >> 3) % 61) * 32 - 960
        return ldexp(m, e)
    end
    bad = 0
    for _ in 1:150
        p2 = ntuple(_ -> (nextcoord(), nextcoord()), 4)
        p3 = ntuple(_ -> (nextcoord(), nextcoord(), nextcoord()), 5)
        orient2(p2[1], p2[2], p2[3]) == oracle_orient2(p2[1], p2[2], p2[3]) || (bad += 1)
        incircle(p2...) == oracle_incircle(p2...) || (bad += 1)
        orient3(p3[1], p3[2], p3[3], p3[4]) == oracle_orient3(p3[1], p3[2], p3[3], p3[4]) || (bad += 1)
        insphere(p3...) == oracle_insphere(p3...) || (bad += 1)
    end
    @test bad == 0
end
