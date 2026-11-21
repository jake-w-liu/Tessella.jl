module HighOrderArcMidpointTests
using Test, Tessella
const A = Tessella.API

# The second-order node of a Circle arc is placed by Gmsh's Newton iteration on
# arc lengths. Its Jacobian is a finite difference with a step of 1e-5 of the
# parameter range, so a change in the last bit of a sine or cosine moves the
# converged node by about 1e-10.
#
# Which platforms reproduce Gmsh 4.15.2's digits:
# - x86 Windows, compiled: bit for bit. Gmsh's Windows build evaluates the
#   paired sine and cosine with the x87 FSINCOS instruction, and so does
#   Tessella there (`Model._model_high_order_fsincos_exact()` is true).
# - x86 Windows, interpreted (`--compile=min`, as this suite runs): the
#   instruction is evaluated in software, because `llvmcall` needs the
#   compiler. That evaluation is correctly rounded and the instruction is
#   not: they differ by one ulp on about 3 in 100000 arguments. On these two
#   arcs they agree and the nodes come out bit-identical, but that is not
#   guaranteed for every arc.
# - every other platform: the system math library, whose sine and cosine
#   differ from FSINCOS by an ulp more often and by more next to multiples of
#   π/2, where the FPU's 66-bit π shows.
# The bound for the last two is 5e-10. The largest deviation observed is
# 6.5e-11, for the radius-2 arc when the C runtime's cosine is used: its
# cosine of the Float64 nearest π/2 is 6.123233995736766e-17 against the FPU's
# 6.123031769111886e-17, and the node came out at x = 1.4142135587734221
# (this was the interpreted result before the software evaluation existed).
# 5e-10 leaves a factor of 8 over that.
const EXACT = Tessella.Model._model_high_order_fsincos_exact()
const BOUND = 5e-10
close_to(p, q) = EXACT ? p == q : maximum(abs.(p .- q)) <= BOUND

@testset "Second-order arc nodes are Gmsh's" begin
    A.finalize(); A.initialize()
    try
        A.option("General.Terminal", 0)
        mktempdir() do directory
            path = joinpath(directory, "arc_p2.geo")
            write(path, """
                Point(1)={2,0,0,1};Point(2)={0,2,0,1};Point(5)={0,0,0,1};
                Point(3)={1,0,0,1};Point(4)={0,1,0,1};
                Circle(2)={1,5,2};
                Circle(4)={4,5,3};
                Transfinite Curve{2,4}=2;
                """)
            A.open_geo!(path)
        end
        A.option("Mesh.ElementOrder", 2)
        A.mesh.generate(1)
        tags, xyz, _ = A.mesh.get_nodes()
        nodes = [(xyz[3i-2], xyz[3i-1], xyz[3i]) for i in eachindex(tags)]
        @test length(nodes) == 7
        outer = only(p for p in nodes if p[1] > 1.2 && p[2] > 1.2)
        inner = only(p for p in nodes if 0.5 < p[1] < 0.9 && 0.5 < p[2] < 0.9)
        # Gmsh 4.15.2, Windows build: nodes 6 and 7 of the same model.
        @test close_to(outer, (1.4142135587089986, 1.4142135660371915, 0.0))
        @test close_to(inner, (0.70710678301859575, 0.70710677935449928, 0.0))
        @test abs(hypot(outer[1], outer[2]) - 2) <= 4eps(2.0)
        @test abs(hypot(inner[1], inner[2]) - 1) <= 4eps(1.0)
    finally
        A.finalize()
    end
end

# The software evaluation reduces by quadrants of the FPU's 66-bit π: the
# cosine of the Float64 nearest π/2 is what `fcos` returns, not the exact one.
@testset "Software FSINCOS reproduces the FPU near π/2" begin
    quarter = 1.5707963267948966
    s, c = Tessella.GmshLibm._x87_fsincos_emulated(quarter)
    @test s == 1.0
    @test c == reinterpret(Float64, 0x3c91a60000000000)   # exact cosine: 0x3c91a62633145c07
    @test Tessella.GmshLibm._x87_fsincos_emulated(0.5) == (sin(0.5), cos(0.5))
    # On x86 Windows both ways of evaluating the pair give that value: the
    # instruction when compiled, the software evaluation when interpreted.
    if Sys.iswindows() && Sys.ARCH in (:x86_64, :i686)
        @test Tessella.Model._model_high_order_curve_sincos(quarter) == (s, c)
    end
    @test EXACT == (Sys.iswindows() && Sys.ARCH in (:x86_64, :i686) &&
                    Base.JLOptions().compile_enabled in (1, 2))
end
end
