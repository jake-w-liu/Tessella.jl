using Test
using Tessella

# A small full-parser regression. Its native execution dispatcher has substantial
# cold compilation cost, so it is explicitly opt-in alongside the parser suites.
@testset "bounded geo execution and error handling" begin
    mktempdir() do directory
        path = joinpath(directory, "expressions.geo")
        write(path, """
            n = 2;
            For i In {1:n}
              Point(i) = {i - 1, 0, 0, 1};
            EndFor
            If (n == 2)
              Line(1) = {1, 2};
            EndIf
            """)
        result = execute_geo(path)
        @test result.mesh === nothing && isempty(result.warnings)
        @test result.model.points[2] == (1.0, 0.0, 0.0)
        @test Tessella.Model.model_entities(result.model, 1) == [(1, 1)]
        write(path, "Line(1) = {1, 99};")
        @test_throws ArgumentError execute_geo(path)
    end
end
