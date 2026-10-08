using Test
using SHA
using Tessella

# FXSAVE is a test-only snapshot: production only reads/writes the two-byte
# control word. The abridged tag identifies leaked x87 stack entries.
module WindowsPowerTest
using Tessella
const Libm=Tessella.GmshLibm
@static if Sys.iswindows() && Sys.ARCH===:x86_64
    @inline function state()
        Base.llvmcall(raw"""
            %state = alloca [512 x i8], align 16
            call void asm sideeffect "fxsave ($0)", "r,~{memory}"(ptr %state)
            %swp = getelementptr [512 x i8], ptr %state, i64 0, i64 2
            %twp = getelementptr [512 x i8], ptr %state, i64 0, i64 4
            %cw = load i16, ptr %state, align 2
            %sw = load i16, ptr %swp, align 2
            %tw = load i8, ptr %twp, align 1
            %c = zext i16 %cw to i64
            %s = zext i16 %sw to i64
            %t = zext i8 %tw to i64
            %ss = shl i64 %s, 16
            %tt = shl i64 %t, 32
            %a = or i64 %c, %ss
            %r = or i64 %a, %tt
            ret i64 %r
            """,UInt64,Tuple{})
    end
    unchanged(before,after)=UInt16(before&0xffff)==UInt16(after&0xffff) &&
        (before>>27)&7==(after>>27)&7 && (before>>32)&0xff==(after>>32)&0xff==0
    function compare(inputs,expected)
        mismatches=Int[];state_errors=Int[]
        for i in eachindex(inputs)
            x,y=inputs[i]
            before=state()
            actual=Libm._gm87_pow(x,y)
            after=state()
            reinterpret(UInt64,actual)==expected[i] || push!(mismatches,i)
            unchanged(before,after) || push!(state_errors,i)
        end
        return mismatches,state_errors
    end
    function repeated()
        before=state()
        total=0.0
        for i in 1:10_000
            total+=Libm._gm87_pow(isodd(i) ? nextfloat(1.0) : prevfloat(1.0),2147483648.0)
        end
        return total,unchanged(before,state())
    end
    function batch!(output,inputs)
        for i in eachindex(output)
            x,y=inputs[mod1(i,length(inputs))]
            output[i]=Libm._gm87_pow(x,y)
        end
        return output
    end
    # Specializing on the callback keeps harness dispatch and return boxing
    # outside the measurement; the output arrays are retained by the caller.
    measure(callback::F) where F=(callback(); @allocated callback())
end
end

@testset "Windows power primary corpus and FPU resource state" begin
    @static if Sys.iswindows() && Sys.ARCH===:x86_64
        libm=Tessella.GmshLibm
        fixture=normpath(joinpath(@__DIR__,"../../validation/windows_power/corpus.tsv"))
        @test bytes2hex(sha256(read(fixture)))==
            "71bc0963d24b8b409b84c126005d66c48ca521b58f7dc79a870743335af46bee"
        rows=[split(row,'\t') for row in readlines(fixture)[2:end]]
        inputs=[(reinterpret(Float64,parse(UInt64,row[2];base=16)),
                 reinterpret(Float64,parse(UInt64,row[3];base=16))) for row in rows]
        @test length(inputs)==2874
        saved=libm._win_x87_control_word()
        try
            for (column,mode) in enumerate((0x0200,0x0600,0x0a00,0x0e00,
                                            0x0300,0x0700,0x0b00,0x0f00,
                                            0x0000,0x0400,0x0800,0x0c00))
                expected=[parse(UInt64,row[column+3];base=16) for row in rows]
                caller=(saved&0xf0ff)|UInt16(mode)
                libm._win_x87_set_control_word!(caller)
                mismatches,state_errors=WindowsPowerTest.compare(inputs,expected)
                @test isempty(mismatches)
                @test isempty(state_errors)
                @test libm._win_x87_control_word()==caller
            end
            # CLI task context uses the PC64 primary while leaving the PC53
            # caller word untouched between operations, in every rounding mode.
            for (index,mode) in enumerate((0x0200,0x0600,0x0a00,0x0e00))
                expected=[parse(UInt64,row[index+7];base=16) for row in rows]
                caller=(saved&0xf0ff)|UInt16(mode)
                libm._win_x87_set_control_word!(caller)
                mismatches,state_errors=libm._with_win_cli_precision() do
                    WindowsPowerTest.compare(inputs,expected)
                end
                @test isempty(mismatches)
                @test isempty(state_errors)
                @test libm._win_x87_control_word()==caller
                @test !libm._win_cli_precision()
            end
            libm._win_x87_set_control_word!((saved&0xf0ff)|UInt16(0x0200))
            @test WindowsPowerTest.repeated()[2]
            @test WindowsPowerTest.measure(WindowsPowerTest.repeated)==0
            @test libm._with_win_cli_precision() do
                WindowsPowerTest.measure(WindowsPowerTest.repeated)==0
            end
            for count in (1000,2000,4000)
                output=Vector{Float64}(undef,count)
                @test WindowsPowerTest.measure(() -> WindowsPowerTest.batch!(output,inputs))==0
                @test length(output)==count
                @test sizeof(output)==8count
            end
        finally
            libm._win_x87_set_control_word!(saved)
        end
        @test libm._win_x87_control_word()==saved

        # FPU control is per OS thread. Each worker restores its own caller;
        # no yielding or Julia global precision mutation occurs in the kernel.
        states=zeros(Int,Threads.nthreads())
        Threads.@threads :static for worker in eachindex(states)
            thread_saved=libm._win_x87_control_word()
            try
                libm._win_x87_set_control_word!((thread_saved&0xf0ff)|UInt16(0x0300))
                states[worker]=WindowsPowerTest.repeated()[2] ? 10_000 : -1
            finally
                libm._win_x87_set_control_word!(thread_saved)
            end
            libm._win_x87_control_word()==thread_saved || (states[worker]=-1)
        end
        @test sum(states)==10_000Threads.nthreads()

        # Audit the actual production definitions, including the CLI callback.
        for function_object in (libm._win_pow,libm._win_pow_native,libm._win_cli_precision,
                                libm._win_powi,libm._win_powi_product,
                                libm._win_x87_powlog,libm._with_win_cli_precision,
                                Tessella.CLI.main,Tessella.CLI._main)
            for method in methods(function_object)
                code=Base.uncompressed_ast(method)
                @test !any(code.code) do statement
                    statement isa Expr && statement.head==:(=) || return false
                    rhs=statement.args[2]
                    return rhs isa Expr && rhs.head==:call && rhs.args[1]==GlobalRef(Core,:Box)
                end
            end
        end
    else
        # Inline x87 is deliberately unavailable off Windows x86_64.
        @test !Tessella.GmshLibm._WIN_X87_NATIVE
        for (x,y) in ((2.0,3.0),(0.5,-4.0),(1.25,0.3))
            @test Tessella.GmshLibm._win_pow(x,y)===Tessella.GmshLibm._gm_pow(x,y)
        end
    end
end
