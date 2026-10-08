using Test
using Tessella

@testset "Windows integer powers preserve large exponent direction" begin
    # Powers of two have an independent exact exponent oracle. Large integer
    # exponents belong to the primary log2l/exp2l path; their product must retain
    # its overflow/underflow direction rather than wrap an integer accumulator.
    for power in (-1023,-1000,-538,-537,-512,-58,-20,-2,-1,0,1,2,20,58,
                  512,537,538,1000,1023),
        count in (-9e18,-4e18,-3425408785282518016.0,-1e18,-1075.0,-1074.0,
                  -1024.0,-2.0,-1.0,0.0,1.0,2.0,1023.0,1024.0,1074.0,
                  1075.0,1e18,3425408785282518016.0,4e18,9e18),
        negative in (false,true)
        exponent_exact=big(power)*BigInt(count)
        magnitude=exponent_exact>1023 ? Inf : exponent_exact < -1074 ? 0.0 :
                  ldexp(1.0,Int(exponent_exact))
        expected=negative && isodd(BigInt(count)) ? -magnitude : magnitude
        base=negative ? -ldexp(1.0,power) : ldexp(1.0,power)
        @test Tessella.GmshLibm._win_pow(base,count) === expected
    end
    # A non-power-of-two base also retains its reciprocal subnormal; the
    # independent 256-bit reciprocal and Gmsh 4.15.2 agree on this rounding.
    @test Tessella.GmshLibm._win_pow(1e155,-2.0) ===
          Float64(inv(BigFloat(1e155;precision=256)^2))
    for base in (5.647527718625593e161,-5.647527718625593e161,1e161,1e162)
        @test Tessella.GmshLibm._win_pow(base,-2.0) ===
              Float64(inv(BigFloat(base;precision=256)^2))
    end
    # Large integral exponents use native extended registers without allocation.
    measure(f)=(f(); @allocated f())
    for count in (3425408785282518016.0,-3425408785282518016.0,4500000000000001.0)
        @test measure(()->Tessella.GmshLibm._win_pow(2.0^58,count))==0
    end
end
