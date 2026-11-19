# `Printf` arguments may contain comparisons (`StrLen(CurrentFileName) > 0`):
# the statement is split at the argument list's own closing parenthesis, so
# `> 0)` is no longer taken for a `> file` redirection. Quoted strings may
# span lines (tutorials/t4.geo's DoubleClickedCurveCommand). Gmsh 4.15.2
# prints `c 1`, `f 7` and the two-line string for this script.
module GeoPrintfArgsTest
using Test, Tessella

@testset "geo Printf arguments and multi-line strings" begin
    dir=mktempdir()
    path=joinpath(dir,"cf2.geo")
    write(path,"""
        Printf("c %g", StrLen("abc") > 0);
        Printf("f %g", StrLen(CurrentFileName));
        Printf("x %g", 3) > "pr.txt";
        s = "a
        b";
        Geometry.DoubleClickedCurveCommand = "Printf('Curve %g',
          Geometry.DoubleClickedEntityTag);";
        """)
    output=mktemp() do out,io
        ex=redirect_stdout(io) do
            execute_geo(path)
        end
        flush(io)
        @test isempty(ex.errors)
        @test ex.strings["s"]==["a\nb"]
        read(out,String)
    end
    @test occursin("c 1",output)
    @test occursin("f 7",output)
    @test strip(read(joinpath(dir,"pr.txt"),String))=="x 3"
end
end
