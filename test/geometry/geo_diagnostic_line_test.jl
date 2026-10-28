# `.geo` diagnostic file/line prefixes — Gmsh 4.15.2's `yyerror`/`yymsg`
# print `'gmsh_yyname', line gmsh_yylineno - 1: <msg>` where the reported
# number counts the line feeds the lexer has *consumed*. Because `YY_INPUT`
# refills one source line per buffer and bumps `yylineno` only on `'\n'`,
# the same statement reports different numbers on LF and CRLF input: a LF
# line's feed is inside its own refill (tokens report their true line)
# while a CRLF line stops at `'\r'` and defers the feed into the next
# refill (tokens report one less), exactly like an unterminated last line.
# Statement actions (`Could not add line`, affectation errors) report the
# terminator line; mid-expression errors report the lookahead line; a
# `Call`ed function body reports inside the file that *defined* it —
# `FunctionManager::enterFunction` swaps `gmsh_yyname`/`gmsh_yylineno` for
# the duration of the call. Every expectation below was pinned against the
# pinned Gmsh 4.15.2 binary on equivalent scripts.
using Test
using Tessella

function _capture_stderr(f)
    return mktemp() do errpath,errstream
        result=redirect_stderr(errstream) do
            try f();nothing catch error;error end
        end
        flush(errstream)
        (result,read(errpath,String))
    end
end

function _execute_captured(source::AbstractString;
                           files=Dict{String,String}(),mesh_dim=0,
                           filename="script.geo")
    dir=mktempdir()
    path=joinpath(dir,filename)
    write(path,source)
    for (name,content) in files
        write(joinpath(dir,name),content)
    end
    return (path,_capture_stderr(()->execute_geo(path;mesh_dim=mesh_dim)))
end

_err_lines(err::AbstractString)=
    [line for line in split(err,'\n') if startswith(line,"Error   :")]
_warn_lines(err::AbstractString)=
    [line for line in split(err,'\n') if startswith(line,"Warning :")]

@testset "LF lines report their physical line" begin
    path,(result,err)=_execute_captured(
        "Point(1) = {0, 0, 0, 1};\nbad = undef;\nPoint(2) = {1, 0, 0, 1};\n")
    @test result===nothing
    @test _err_lines(err)==
        ["Error   : '$path', line 2: Unknown variable 'undef'"]
    # A second failing statement keeps reporting its own line.
    path,(result,err)=_execute_captured(
        "Point(1) = {0, 0, 0, 1};\nbad_a = undef_a;\nPoint(2) = {1, 0, 0, 1};"*
        "\nbad_b = undef_b;\nPoint(3) = {2, 0, 0, 1};\n")
    @test _err_lines(err)==
        ["Error   : '$path', line 2: Unknown variable 'undef_a'",
         "Error   : '$path', line 4: Unknown variable 'undef_b'"]
end

@testset "CRLF and unterminated lines report one less" begin
    # `YY_INPUT` leaves a CRLF line's `'\n'` in the next refill, so the
    # feed count lags one line behind the LF convention.
    dir=mktempdir();path=joinpath(dir,"script.geo")
    write(path,"Point(1) = {0, 0, 0, 1};\r\nbad_a = undef_a;\r\n"*
               "Point(2) = {1, 0, 0, 1};\r\nbad_b = undef_b;\r\n"*
               "Point(3) = {2, 0, 0, 1};\r\n")
    (result,err)=_capture_stderr(()->execute_geo(path))
    @test _err_lines(err)==
        ["Error   : '$path', line 1: Unknown variable 'undef_a'",
         "Error   : '$path', line 3: Unknown variable 'undef_b'"]
    # An unterminated final line has no feed to consume either.
    write(path,"Point(1) = {0, 0, 0, 1};\nbad = undef;")
    (result,err)=_capture_stderr(()->execute_geo(path))
    @test _err_lines(err)==
        ["Error   : '$path', line 1: Unknown variable 'undef'"]
end

@testset "statement actions report the terminator line" begin
    # `Could not add line` fires at the `tEND`-adjacent reduction upstream —
    # on the `;` line, not on the last expression token's line.
    path,(result,err)=_execute_captured(
        "Point(1) = {0, 0, 0, 1};\nLine(1) = {\n  999, 998\n};\n"*
        "Point(2) = {1, 0, 0, 1};\n")
    @test "Error   : '$path', line 4: Could not add line" in _err_lines(err)
    # Mid-expression errors keep their lookahead position instead: the
    # reduce of `undef2` sees `+` on its own line; `undef_var`'s reduce
    # already sees `+` on the terminator line.
    path,(result,err)=_execute_captured(
        "Point(1) = {0, 0, 0, 1};\nbad = 1 +\n  undef_var\n  + 3;\n"*
        "bad2 = undef2 +\n  4;\nPoint(2) = {1, 0, 0, 1};\n")
    @test _err_lines(err)==
        ["Error   : '$path', line 4: Unknown variable 'undef_var'",
         "Error   : '$path', line 5: Unknown variable 'undef2'"]
end

@testset "Include and Call report their own file and line" begin
    # `ParseFile` swaps `gmsh_yyname`/`gmsh_yylineno` for the child file
    # and restores them when it returns.
    dir=mktempdir()
    write(joinpath(dir,"inner.geo"),
          "// comment\ninner_bad = undef_inner;\nPoint(5) = {2, 0, 0, 1};\n")
    path=joinpath(dir,"outer.geo")
    write(path,"Point(1) = {0, 0, 0, 1};\nInclude \"inner.geo\";\n"*
               "Point(2) = {1, 0, 0, 1};\nouter_bad = undef_outer;\n")
    (result,err)=_capture_stderr(()->execute_geo(path))
    @test _err_lines(err)==
        ["Error   : '$(joinpath(dir,"inner.geo"))', line 2: Unknown variable 'undef_inner'",
         "Error   : '$path', line 4: Unknown variable 'undef_outer'"]
    # A `Call`ed body reports inside the file that defined it — upstream's
    # `enterFunction` swaps `gmsh_yyname` for the body's definition site.
    write(joinpath(dir,"inner.geo"),
          "// comment\nFunction mk\n  inner_bad = undef_inner;\nReturn\n")
    write(path,"Include \"inner.geo\";\nCall mk;\nPoint(2) = {1, 0, 0, 1};\n")
    (result,err)=_capture_stderr(()->execute_geo(path))
    @test _err_lines(err)==
        ["Error   : '$(joinpath(dir,"inner.geo"))', line 3: Unknown variable 'undef_inner'"]
end

@testset "syntax errors, warnings and Msg channels" begin
    # `yyerror` diagnostics carry the same prefix as `yymsg(0)`.
    path,(result,err)=_execute_captured(
        "Point(1) = {0, 0, 0, 1};\na = 1 +;\nPoint(2) = {1, 0, 0, 1};\n")
    @test _err_lines(err)==
        ["Error   : '$path', line 2: syntax error (;)"]
    # `Warning("...")`/`Error("...")` run through `Msg::Warning`/`Msg::Error`
    # — plain channels with no parser prefix.
    path,(result,err)=_execute_captured(
        "Warning(\"hello\");\nError(\"bomb\");\nPoint(1) = {0, 0, 0, 1};\n")
    @test _err_lines(err)==["Error   : bomb"]
    @test _warn_lines(err)==["Warning : hello"]
    # `SetFactory` rejects an unknown factory through `yymsg(1)` — prefixed.
    path,(result,err)=_execute_captured(
        "SetFactory(\"bogus\");\nPoint(1) = {0, 0, 0, 1};\n")
    @test _warn_lines(err)==
        ["Warning : '$path', line 1: Unknown factory \"bogus\" - using " *
         "\"Built-in\" instead"]
end

@testset "API .geo execution keeps the same continuation stream" begin
    # `API.open` routes through the same `execute_geo` engine — identical
    # file/line diagnostics and continuation.
    dir=mktempdir()
    write(joinpath(dir,"inner.geo"),"inner_bad = undef_inner;\n")
    path=joinpath(dir,"outer.geo")
    write(path,"Include \"inner.geo\";\nouter_bad = undef_outer;\n"*
               "Point(1) = {0, 0, 0, 1};\n")
    (result,err)=_capture_stderr() do
        Tessella.API.initialize()
        Tessella.API.open(path)
    end
    @test _err_lines(err)==
        ["Error   : '$(joinpath(dir,"inner.geo"))', line 1: Unknown variable 'undef_inner'",
         "Error   : '$path', line 2: Unknown variable 'undef_outer'"]
end
