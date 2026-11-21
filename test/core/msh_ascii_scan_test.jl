module MshAsciiScanTest

using Test, Tessella

const E = Tessella.Elements

# The ASCII element readers translate well-formed one-record lines in place and
# hand every other line to the general tokenizer. Both routes must agree.

const V2 = """
\$MeshFormat
2.2 0 8
\$EndMeshFormat
\$Nodes
5
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
5 1 1 0
\$EndNodes
\$Elements
3
1 2 2 7 1 1 2 3
2 4 2 8 2 1 2 3 4
3 3 2 9 3 2 3 5 1
\$EndElements
"""

const V4 = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
0 0 2 1
1 0 0 0 1 1 0 1 7 0
2 0 0 0 1 1 1 1 8 0
3 0 0 0 1 1 1 0 0
\$EndEntities
\$Nodes
1 5 1 5
3 3 0 5
1
2
3
4
5
0 0 0
1 0 0
0 1 0
0 0 1
1 1 0
\$EndNodes
\$Elements
2 3 1 3
2 1 2 2
1 1 2 3
2 2 3 5
3 3 4 1
3 1 2 3 4
\$EndElements
"""

# Same content with every irregularity the scanner must not mishandle.
variant(text) = replace(text, "\n" => "\r\n")
spread(text) = replace(text, "1 1 2 3\n" => "1\n1 2\n3\n", "1 2 3 4\n" => "\n 1  2\t3 4  \n")
# MSH 2.2 records are one line each: only the spacing may vary.
spaced(text) = replace(text, "1 2 2 7 1 1 2 3" => "  1  2\t2 7 1  1 2 3 ",
                       "3 3 2 9 3 2 3 5 1" => "3 3 2 9 3 2 3 5 1\t")
signed(text) = replace(text, "1 2 2 7 1 1 2 3" => "1 2 2 +7 1 1 2 3")

function crc(path)
    m = E.read_mixed_msh(path)
    return (size(m.coords), [(b.msh, copy(b.nodes), copy(b.tags)) for b in m.blocks])
end

@testset "ASCII element lines: fast scan and general tokenizer agree" begin
    mktempdir() do dir
        for (name, text, derived) in (("v2", V2, (variant, spaced, signed)),
                                      ("v4", V4, (variant, spread)))
            base = joinpath(dir, name * ".msh")
            write(base, text)
            reference = crc(base)
            for f in derived
                path = joinpath(dir, name * "_" * string(f) * ".msh")
                write(path, f(text))
                @test crc(path) == reference
            end
        end
    end
end

@testset "ASCII element diagnostics are unchanged" begin
    mktempdir() do dir
        function fails(text; message = nothing, kwargs...)
            path = joinpath(dir, "bad.msh")
            write(path, text)
            err = try E.read_mixed_msh(path; kwargs...); nothing catch e; e end
            return err isa ArgumentError && (message === nothing || occursin(message, err.msg))
        end
        @test fails(replace(V2, "3 3 2 9 3 2 3 5 1" => "2 3 2 9 3 2 3 5 1"); message = "duplicate element tag")
        @test fails(replace(V2, "2 4 2 8 2 1 2 3 4" => "2 4 2 8 2 1 2 3 9"); message = "unknown node tag")
        @test fails(replace(V2, "1 2 2 7 1 1 2 3" => "0 2 2 7 1 1 2 3"); message = "positive")
        @test fails(replace(V2, "1 2 2 7 1 1 2 3" => "1 2 2 7 1 1 2"))
        @test fails(replace(V4, "2 2 3 5" => "2 2 3 9"); message = "unknown node tag")
        @test fails(replace(V4, "1 1 2 3\n" => "2 1 2 3\n"); message = "duplicate element tag")
        # A record of 300 000 tokens is refused without reading it as a line.
        long = replace(V2, "3 3 2 9 3 2 3 5 1" => "3 34 0 300000" * " 1"^300_000)
        @test fails(long; max_connectivity = 0)
    end
end

end # module
