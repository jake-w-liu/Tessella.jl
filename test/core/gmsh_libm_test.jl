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

module WindowsAtan2RegressionTests
using Test,Tessella
const G=Tessella.GmshLibm

@static if Sys.iswindows() && Sys.ARCH===:x86_64
    @inline function x87_stack_top()
        Base.llvmcall(raw"""
            %word = alloca i16, align 2
            call void asm sideeffect "fnstsw ($0)", "r,~{memory}"(ptr %word)
            %answer = load i16, ptr %word, align 2
            %top = and i16 %answer, 14336
            ret i16 %top
            """,UInt16,Tuple{})
    end

    # Actual pinned Gmsh 4.15.2 public Atan2 parser outputs under each RC.
    const pairs=(
        (0x3fc0000000000000,0xc004000000000000),
        (0xc00a000000000000,0x3fd8000000000000),
        (0x3ff0000000000000,0x3ff0000000000000),
        (0x3ff0000000000000,0xbff0000000000000),
        (0x0000000000000001,0x7fefffffffffffff),
        (0x7fefffffffffffff,0x0000000000000001),
        (0x0000000000000000,0x8000000000000000),
        (0x8000000000000000,0x8000000000000000),
        (0x7ff0000000000000,0xfff0000000000000))
    const expected=(
        (0x4008bbaabde5e29c,0xbff74b727a4eb80a,0x3fe921fb54442d18,0x4002d97c7f3321d2,0x0000000000000000,0x3ff921fb54442d18,0x400921fb54442d18,0xc00921fb54442d18,0x4002d97c7f3321d2),
        (0x4008bbaabde5e29b,0xbff74b727a4eb80a,0x3fe921fb54442d18,0x4002d97c7f3321d2,0x0000000000000000,0x3ff921fb54442d18,0x400921fb54442d18,0xc00921fb54442d19,0x4002d97c7f3321d2),
        (0x4008bbaabde5e29c,0xbff74b727a4eb809,0x3fe921fb54442d19,0x4002d97c7f3321d3,0x0000000000000001,0x3ff921fb54442d19,0x400921fb54442d19,0xc00921fb54442d18,0x4002d97c7f3321d3),
        (0x4008bbaabde5e29b,0xbff74b727a4eb809,0x3fe921fb54442d18,0x4002d97c7f3321d2,0x0000000000000000,0x3ff921fb54442d18,0x400921fb54442d18,0xc00921fb54442d18,0x4002d97c7f3321d2))

    function mode_values(word)
        values=Vector{UInt64}(undef,length(pairs))
        saved=G._win_x87_control_word();before=x87_stack_top()
        cw=UInt16(0);top=UInt16(0)
        try
            G._win_x87_set_control_word!(word)
            @inbounds for index in eachindex(pairs)
                y,x=pairs[index]
                values[index]=reinterpret(UInt64,G._win_atan2(
                    reinterpret(Float64,y),reinterpret(Float64,x)))
            end
            cw=G._win_x87_control_word();top=x87_stack_top()
        finally
            G._win_x87_set_control_word!(saved)
        end
        return values,cw,top==before
    end

    @noinline function repeated_angles(count)
        total=0.0
        for index in 1:count
            total+=G._win_atan2(1.,Float64(index))
        end
        return total
    end

    @testset "Windows atan2 exact rounding, state and bounded workspace" begin
        saved=G._win_x87_control_word()
        G._win_atan2(1.,1.);x87_stack_top()
        for pc in (UInt16(0),UInt16(0x0200),UInt16(0x0300)),
            (index,rc) in enumerate((UInt16(0),UInt16(0x0400),UInt16(0x0800),UInt16(0x0c00)))
            word=(saved&UInt16(0xf0ff))|pc|rc
            values,cw,top_held=mode_values(word)
            @test values==collect(expected[index])
            @test cw==word
            @test top_held
            @test G._win_x87_control_word()==saved
        end
        @test isequal(G._win_atan2(-0.,1.),-0.)
        @test isnan(G._win_atan2(NaN,1.))
        @test isnan(G._win_atan2(1.,NaN))
        repeated_angles(1000)
        for count in (1000,2000,4000)
            @test minimum(@allocated(repeated_angles(count)) for _ in 1:3)==0
        end
        @test G._win_x87_control_word()==saved
    end
end
end # module WindowsAtan2RegressionTests
