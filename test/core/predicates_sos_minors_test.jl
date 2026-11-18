# Regression: the exact +ε SoS evaluator used to enumerate every generalized
# cofactor as a freshly allocated Rational{BigInt} matrix (orient3_sos on four
# coplanar points: ~20 µs and ~1400 allocations; insphere_sos on five cospherical
# points: ~110 µs).  Minors that keep only Float64 coordinate columns are now
# decided by the exact Float64 predicates.  The signs must equal the rational
# predicates (whose minors all take the BigInt path) on every degenerate input.
using Test
using Tessella.Predicates
import Tessella.Predicates as P

@noinline _o3sos_alloc(a, b, c, d) = @allocated orient3_sos(a, b, c, d, 4, 2, 7, 1)
@noinline _o2sos_alloc(a, b, c) = @allocated orient2_sos(a, b, c, 3, 1, 2)

@testset "SoS degenerate minors: Float64 fast minors match exact rational minors" begin
    R = Rational{BigInt}
    q2(p) = (R(p[1]), R(p[2])); q3(p) = (R(p[1]), R(p[2]), R(p[3]))
    s = UInt64(0x2545F4914F6CDD1D)
    nexti() = (s ⊻= s << 13; s ⊻= s >> 7; s ⊻= s << 17; Int(s % 5) - 2)
    mism = 0; ndeg = 0
    for _ in 1:600
        p3 = ntuple(_ -> (Float64(nexti()), Float64(nexti()), 0.5 * nexti()), 5)
        p2 = ntuple(_ -> (Float64(nexti()), 0.5 * nexti()), 4)
        ix = (3, 9, 1, 6, 4)
        orient3(p3[1:4]...) == 0 && (ndeg += 1)
        orient3_sos(p3[1:4]..., ix[1:4]...) == orient3_rat(map(q3, p3[1:4])..., ix[1:4]...) || (mism += 1)
        orient2_sos(p2[1:3]..., ix[1:3]...) == orient2_rat(map(q2, p2[1:3])..., ix[1:3]...) || (mism += 1)
        incircle_sos(p2..., ix[1:4]...) == incircle_rat(map(q2, p2)..., ix[1:4]...) || (mism += 1)
        insphere_sos(p3..., ix...) == insphere_rat(map(q3, p3)..., ix...) || (mism += 1)
    end
    @test ndeg > 50          # the sweep really exercises the degenerate path
    @test mism == 0
    # cube corners: exactly cospherical with coplanar faces
    c = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    @test insphere(c...) == 0
    @test insphere_sos(c..., 1, 2, 3, 4, 5) == insphere_rat(map(q3, c)..., 1, 2, 3, 4, 5)
    @test insphere_sos(c..., 5, 4, 3, 2, 1) == insphere_rat(map(q3, c)..., 5, 4, 3, 2, 1)
    # degenerate orient2_sos / orient3_sos allocate nothing
    a = (0.0, 0.0, 0.0); b = (1.0, 0.0, 0.0); cc = (0.0, 1.0, 0.0); d = (0.25, 0.5, 0.0)
    orient3_sos(a, b, cc, d, 4, 2, 7, 1)
    @test _o3sos_alloc(a, b, cc, d) == 0
    orient2_sos((0.0, 0.0), (1.0, 1.0), (2.0, 2.0), 3, 1, 2)
    @test _o2sos_alloc((0.0, 0.0), (1.0, 1.0), (2.0, 2.0)) == 0
end
