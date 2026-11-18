# Control headers may span several lines (`If(a &&\n b)` in tutorials/t18.geo):
# Gmsh's flex scanner has no line structure. Output matches Gmsh 4.15.2.
using Test
using Tessella

@testset "geo multi-line control headers" begin
    dir=mktempdir()
    path=joinpath(dir,"control.geo")
    write(path,"""
        a = 1; b = 2; r = 0;
        If (a == 1 &&
            b == 2) // trailing comment
          r += 1;
        EndIf
        If (a == 2 ||
            // a comment line
            b == 3)
          r += 10;
        ElseIf (a +
                b == 3)
          r += 100;
        EndIf
        n = 0;
        For i In {1:
                  3}
          n += i;
        EndFor
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.values["r"]==101
    @test ex.values["n"]==6
    # The parameter scanner splits statements the same way.
    params=Tessella.IO.read_geo_params(path)
    @test isempty(params.scan_errors)
end
