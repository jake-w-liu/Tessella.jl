# `parser` module UI01: the `.geo` parser variable store — get/set/clear,
# regex name queries, and `.geo` file adoption through `parser.parse`.
module ParserChecks
using Test, Tessella
const A=Tessella.API

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "parser variable store" begin
    A.initialize()
    try
        A.parser.clear()
        # set/get scalars and vectors.
        A.parser.set_number("x",1.5)
        A.parser.set_number("y",[2.0,3.0,4.0])
        A.parser.set_string("s","hello")
        A.parser.set_string("t",["a","b"])
        @test A.parser.get_number("x")==[1.5]
        @test A.parser.get_number("y")==[2.0,3.0,4.0]
        @test A.parser.get_string("s")==["hello"]
        @test A.parser.get_string("t")==["a","b"]
        # Missing variables return empty vectors (upstream semantics).
        @test A.parser.get_number("missing")==Float64[]
        @test A.parser.get_string("missing")==String[]
        # Update keeps one record per name across types.
        A.parser.set_number("x",9.0)
        @test A.parser.get_number("x")==[9.0]
        # Name queries take a regular expression; empty returns everything.
        @test A.parser.get_names()==["s","t","x","y"]
        @test A.parser.get_names("^x")==["x"]
        @test A.parser.get_names("^(x|y)\$")==["x","y"]
        @test A.parser.get_names("zzz")==String[]
        @test occursin("regular expression",
            _errtext(()->A.parser.get_names("[")))
        # clear(name) removes one record; clear() empties the store.
        A.parser.clear("x")
        @test A.parser.get_number("x")==Float64[]
        @test A.parser.get_names()==["s","t","y"]
        A.parser.clear()
        @test A.parser.get_names()==String[]
        # Input validation mirrors the API's conventions.
        @test_throws ArgumentError A.parser.set_number("",1.0)
        @test_throws ArgumentError A.parser.set_number("n","s")
        @test_throws ArgumentError A.parser.set_number("n",Inf)
        @test_throws ArgumentError A.parser.set_number("n",[1.0,NaN])
        @test_throws ArgumentError A.parser.set_string("s",1.0)
    finally
        A.finalize()
    end
end

@testset "parser.parse adopts .geo variables" begin
    mktempdir() do dir
        geo=joinpath(dir,"vars.geo")
        write(geo,"""
a = 3.25;
b[] = {1.0, 2.0, 5.5};
c = "text";
""")
        A.initialize()
        try
            A.parser.clear()
            A.parser.parse(geo)
            @test A.parser.get_number("a")==[3.25]
            @test A.parser.get_number("b")==[1.0,2.0,5.5]
            @test A.parser.get_string("c")==["text"]
            # Parsed variables merge cumulatively.
            A.parser.set_number("keep",1.0)
            @test A.parser.get_names() |> issorted
            @test "keep" in A.parser.get_names()
            @test "a" in A.parser.get_names()
            # Parsing does not create model geometry upstream-style — the
            # current model carries no entities from the file.
            @test isempty(A.model.get_entities())
            # Missing files fail loudly.
            @test_throws Exception A.parser.parse(joinpath(dir,"none.geo"))
        finally
            A.finalize()
        end
    end
end

@testset "parser lifecycle requirement" begin
    # Like the rest of the session API, parser calls need `API.initialize()`.
    A.initialize(); A.finalize()
    @test_throws Exception A.parser.set_number("x",1.0)
    @test_throws Exception A.parser.get_names()
end

end # module
