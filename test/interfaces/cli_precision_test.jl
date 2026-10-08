using Test
using Tessella

@testset "Windows standalone CLI precision boundary" begin
    @static if Sys.iswindows() && Sys.ARCH===:x86_64
        libm=Tessella.GmshLibm
        saved=libm._win_x87_control_word()
        mktempdir() do folder
            source=joinpath(folder,"precision.geo")
            output=joinpath(folder,"precision.msh")
            script="""
                Point(1)={0,0,0};
                a=(1.0000000000000002)^9007199254740992;
                Point(2)={a,0,0}; Point(3)={0,1,0};
                Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};
                Curve Loop(1)={1,2,3}; Plane Surface(1)={1};
                Transfinite Curve{1,2,3}=2;
                Transfinite Surface{1}={1,2,3};
                """
            write(source,script)
            invalid=joinpath(folder,"invalid.geo")
            write(invalid,"Point(1)={unknown_scalar,0,0};")
            try
                caller53=(saved&0xf0ff)|UInt16(0x0200)
                libm._win_x87_set_control_word!(caller53)
                direct=Tessella.execute_geo(source;mesh_dim=2)
                @test maximum(direct.mesh.coords[1,:])===7.389056098930648
                @test libm._win_x87_control_word()==caller53
                @test Tessella.CLI.main([source,"-2","-o",output])==output
                mesh=Tessella.IO.read_msh(output).mesh
                @test maximum(mesh.coords[1,:])===7.389056098930649
                @test size(mesh.coords,2)==3
                @test validate(mesh).ok
                @test read(source,String)==script
                @test libm._win_x87_control_word()==caller53

                for precision in (0x0200,0x0300),mode in (0,0x0400,0x0800,0x0c00)
                    caller=(saved&0xf0ff)|UInt16(precision|mode)
                    libm._win_x87_set_control_word!(caller)
                    for args in (["--unknown"],[joinpath(folder,"missing.geo")],
                                 [invalid],[source,"-2","-o",folder])
                        @test_throws Exception Tessella.CLI.main(args)
                        @test libm._win_x87_control_word()==caller
                    end
                    value=libm._with_win_cli_precision() do
                        scoped=libm._win_x87_control_word()
                        @test scoped==caller
                        @test libm._win_cli_precision()
                        @test libm._with_win_cli_precision(() -> 17)==17
                        @test libm._win_x87_control_word()==scoped
                        31
                    end
                    @test value==31
                    @test libm._win_x87_control_word()==caller
                    @test !libm._win_cli_precision()
                end

                # Julia task switches do not save/restore the x87 control word.
                # A yielding CLI scope must therefore store only task context.
                libm._win_x87_set_control_word!(caller53)
                channel=Channel{Nothing}(1)
                resume=Channel{Nothing}(1)
                observations=Tuple{Bool,UInt16,Float64}[]
                @sync begin
                    @async libm._with_win_cli_precision() do
                        put!(channel,nothing)
                        take!(resume)
                        push!(observations,(libm._win_cli_precision(),
                            libm._win_x87_control_word(),
                            libm._gm87_pow(nextfloat(1.0),9007199254740992.0)))
                    end
                    @async begin
                        take!(channel)
                        push!(observations,(libm._win_cli_precision(),
                            libm._win_x87_control_word(),
                            libm._gm87_pow(nextfloat(1.0),9007199254740992.0)))
                        put!(resume,nothing)
                    end
                end
                @test observations==[(false,caller53,7.389056098930648),
                                     (true,caller53,7.389056098930649)]
                @test !libm._win_cli_precision()
                @test libm._win_x87_control_word()==caller53

                # A spawned scope may resume on another worker. At every resume
                # only its task-local context follows it; that worker's caller
                # word is preserved by the individual power operation.
                results=fetch.([Threads.@spawn begin
                    original=libm._win_x87_control_word()
                    value=libm._with_win_cli_precision() do
                        for _ in 1:100
                            yield()
                            before=libm._win_x87_control_word()
                            libm._win_cli_precision() || return false
                            libm._gm87_pow(nextfloat(1.0),9007199254740992.0)===
                                7.389056098930649 || return false
                            libm._win_x87_control_word()==before || return false
                        end
                        true
                    end
                    value && !libm._win_cli_precision() &&
                        libm._win_x87_control_word()==original
                end for _ in 1:4])
                @test all(results)

                @test_throws ErrorException libm._with_win_cli_precision() do
                    @test libm._win_cli_precision()
                    error("scope failure")
                end
                @test !libm._win_cli_precision()
                @test libm._win_x87_control_word()==caller53
            finally
                libm._win_x87_set_control_word!(saved)
            end
        end
        @test libm._win_x87_control_word()==saved
    else
        @test Tessella.GmshLibm._with_win_cli_precision(() -> 31)==31
        @test_throws ArgumentError Tessella.CLI.main(String[])
    end
end
