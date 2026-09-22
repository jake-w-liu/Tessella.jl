using Test
using Tessella

const _GE=Tessella.GeoExec
const _GI=Tessella.IO

# `execute_geo` throws accumulated `yymsg(0)` diagnostics; the internal driver
# inspects the context state recoverable statements left behind (warnings,
# `Msg::Error` counts, struct tables).
function _exec_ctx(source::AbstractString)
    return mktemp() do path,io
        write(io,source)
        close(io)
        params=_GI.read_geo_params(path)
        model=Tessella.Model.GeoModel()
        context=_GI._GeoNumericContext()
        context.file_name=String(path)
        context.fields=empty!(params.fields)
        context.soft_unknown_reads=true
        context.stderr_diagnostics=false
        context.entity_name_lookup=(dim,tag,kind)->
            kind===:physical ? get(model.physical_names,(dim,tag),"") :
                get(model.entity_names,(dim,tag),"")
        allocator_state=_GE._GeoTagAllocatorState()
        context.exec_hook=src->_GE._geo_exec_value_term(
            model,src,context,allocator_state)
        statements=_GE._geo_exec_statements(path)
        _GE._exec_geo_statements!(
            model,statements,firstindex(statements),lastindex(statements),
            context,allocator_state,Ref(0))
        _GE._geo_sync_physical_view!(model,context)
        return model,context
    end
end

function _struct_error(source::AbstractString)
    try
        mktemp() do path,io
            write(io,source)
            close(io)
            execute_geo(path)
        end
        return nothing
    catch err
        return err
    end
end

@testset "bounded .geo Struct definitions and member reads" begin
    # Gmsh 4.15.2: `Struct name[key value, ...]` — key/value pairs without
    # `=`; a `DefineStruct` expression returns the auto tag.
    _,ctx=_exec_ctx("""
        x = Struct s2 [a 5];
        Struct s [a 1, b 2];
        t = s;
        av = s.a;
        bv = s.b;
        """)
    @test ctx.values["x"]==1
    @test ctx.values["t"]==2   # second struct gets tag 2
    @test ctx.values["av"]==1
    @test ctx.values["bv"]==2

    # Numeric list members index with `()` or `[]`; `#m()` is the member dim.
    _,ctx2=_exec_ctx("""
        Struct s [a 1, b {2,3}];
        b0 = s.b(0);
        b1 = s.b[1];
        d = #s.b();
        """)
    @test ctx2.values["b0"]==2
    @test ctx2.values["b1"]==3
    @test ctx2.values["d"]==2

    # Char members are stored via `key "str"`/`key Str[...]`/`key Str[{...}]`
    # but are invisible to *numeric* member reads: `x = s.name` takes the
    # FExpr production upstream — `Unknown member` + 0, never a string
    # assignment. `StringExprVar` positions (`StrCat[...]`, `Printf` format)
    # reach the member, including `(i)` indexing.
    _,ctx3=_exec_ctx("""
        Struct s [name Str["hello"], lst Str[{"a","b"}]];
        n = s.name;
        l0 = s.lst(0);
        c = StrCat[s.name];
        e = StrCat[s.lst(1)];
        """)
    @test ctx3.values["n"]==0
    @test ctx3.values["l0"]==0
    @test ctx3.strings["c"]==["hello"]
    @test ctx3.strings["e"]==["b"]
    @test length(filter(e->occursin("Unknown member",e),ctx3.exec_errors))==2

    # A quoted string inside a `{...}` list is a syntax error upstream —
    # string lists need `Str[{...}]` (`lst {"a","b"}` → `syntax error (")`).
    err=_struct_error("Struct s [lst {\"a\",\"b\"}];")
    @test err isa ArgumentError
    @test occursin("syntax error (\")",sprint(showerror,err))
end

@testset "bounded .geo Struct redefinition and Append" begin
    # Redefinition keeps the first body and emits `Redefinition of Struct
    # '::s'`; `(Append 0)` is not an append — same diagnostic.
    _,ctx=_exec_ctx("""
        Struct s [a 1];
        Struct s [b 2];
        Struct s (Append 0) [c 3];
        """)
    @test any(e->occursin("Redefinition of Struct '::s'",e),ctx.exec_errors)
    @test length(filter(e->occursin("Redefinition of Struct",e),
                        ctx.exec_errors))==2
    s=ctx.structs[""]["s"]
    @test haskey(s.fopt,"a") && !haskey(s.fopt,"b") && !haskey(s.fopt,"c")

    # `(Append)` merges only new members; existing members are preserved.
    _,ctx2=_exec_ctx("""
        Struct s [a 1];
        Struct s (Append) [a 9, c 3];
        av = s.a;
        cv = s.c;
        """)
    @test ctx2.values["av"]==1
    @test ctx2.values["cv"]==3
    @test isempty(ctx2.exec_errors)

    # An unknown member read reports `Unknown member 'b' of Struct s` and
    # yields the default 0.
    _,ctx3=_exec_ctx("""
        Struct s [a 1];
        b = s.b;
        """)
    @test any(e->occursin("Unknown member 'b' of Struct s",e),ctx3.exec_errors)
    @test ctx3.values["b"]==0

    # Struct members are read-only: `s.a = 5` parses as a number-option
    # write on the unknown category `s` (upstream `Msg::Error`).
    _,ctx4=_exec_ctx("""
        Struct s [a 1];
        s.a = 5;
        s.b += 2;
        """)
    @test ctx4.msg_error_count==2
    @test ctx4.structs[""]["s"].fopt["a"]==[1.0]
end

@testset "bounded .geo namespaces, DimNameSpace and NameStruct" begin
    # `Struct ns::name` stores under `ns`; `ns::name` reads the tag,
    # `ns::name.member` the member; `DimNameSpace(ns)` counts defs.
    _,ctx=_exec_ctx("""
        Struct ns::s [a 1];
        Struct ns::t [a 2];
        tag = ns::s;
        av = ns::s.a;
        d = DimNameSpace(ns);
        d0 = DimNameSpace();
        """)
    @test ctx.values["tag"]==1
    @test ctx.values["av"]==1
    @test ctx.values["d"]==2
    @test ctx.values["d0"]==0

    # `DimNameSpace("ns")` — a quoted argument is a syntax error at `"`.
    err=_struct_error("""Struct ns::s [a 1];\nd = DimNameSpace("ns");""")
    @test err isa ArgumentError
    @test occursin("syntax error (\")",sprint(showerror,err))

    # `NameStruct(ns::#tag)` resolves the name inside `ns`; `NameStruct(#tag)`
    # looks in the root namespace "" only — `Unknown NameSpace ''` warning
    # when no root struct exists; `Unknown Struct of index n` for a missing
    # tag in an existing namespace.
    _,ctx2=_exec_ctx("""
        Struct ns::s [a 1];
        Struct ns::t [a 2];
        n1 = NameStruct(ns::#1);
        n2 = NameStruct(ns::#2);
        """)
    @test ctx2.strings["n1"]==["s"]
    @test ctx2.strings["n2"]==["t"]

    _,ctx3=_exec_ctx("""
        Struct ns::s [a 1];
        n = NameStruct(#1);
        """)
    @test ctx3.strings["n"]==[""]
    @test any(w->occursin("Unknown NameSpace '' of Struct",w),
              ctx3.exec_warnings)

    _,ctx4=_exec_ctx("""
        Struct s [a 1];
        n = NameStruct(#9);
        """)
    @test ctx4.strings["n"]==[""]
    @test any(w->occursin("Unknown Struct of index 9",w),ctx4.exec_warnings)

    # Malformed `NameStruct` args — upstream bison offender tokens: a bare
    # identifier reduces as `String__Index` first, so `foo)`/`ns)`/`ns::`/
    # `ns::#` all report `)`; a non-name head reports itself.
    for (src,tok) in (
            ("Printf(NameStruct(5));","5"),
            ("Printf(NameStruct());",")"),
            ("Printf(NameStruct(foo));",")"),
            ("Printf(NameStruct(ns));",")"),
            ("Printf(NameStruct(ns::));",")"),
            ("Printf(NameStruct(ns::x));","x"),
            ("Printf(NameStruct(ns::#));",")"),
            ("Printf(NameStruct(#));",")"),
            ("Printf(NameStruct(ns::5));","5"),
            ("Printf(NameStruct(#1,2));",","))
        err=_struct_error("Struct s [a 1];\n"*src)
        @test err isa ArgumentError
        @test occursin("syntax error ($tok)",sprint(showerror,err))
    end

    # `NameStruct` is a `StringExpr` production — it cannot appear in a
    # numeric argument position (`Printf("x %s", NameStruct(#1))` →
    # `syntax error (NameStruct)` upstream).
    err=_struct_error("""Struct s [a 1];\nPrintf("x %s", NameStruct(#1));""")
    @test err isa ArgumentError
    @test occursin("syntax error (NameStruct)",sprint(showerror,err))

    # `ns::` reads resolve only struct names — a namespace-qualified
    # assignment target is a syntax error at `::` upstream.
    err2=_struct_error("ns::x = 5;")
    @test err2 isa ArgumentError
    @test occursin("syntax error (::)",sprint(showerror,err2))

    # `x = ns::#1` — the `#` tag form is `NameStruct_Arg` only, never an
    # FExpr (`syntax error (#)` upstream).
    err3=_struct_error("Struct ns::s [a 1];\nx = ns::#1;")
    @test err3 isa ArgumentError
    @test occursin("syntax error (#)",sprint(showerror,err3))

    # `Exists`/`GetForced` accept `ns::` full names.
    _,ctx5=_exec_ctx("""
        Struct ns::s [a 1];
        e1 = Exists(ns::s.a);
        e2 = Exists(ns::s);
        g = GetForced(ns::s.a, 9);
        """)
    @test ctx5.values["e1"]==1
    @test ctx5.values["e2"]==1
    @test ctx5.values["g"]==1

    # `DimNameSpace(ns);` as a bare statement — expression tokens are not
    # statements (`syntax error (DimNameSpace)` upstream).
    err4=_struct_error("Struct ns::s [a 1];\nDimNameSpace(ns);")
    @test err4 isa ArgumentError
    @test occursin("syntax error (DimNameSpace)",sprint(showerror,err4))
end

@testset "bounded .geo EOF semantics" begin
    # Upstream lexer rules verified against Gmsh 4.15.2: an EOF reached while
    # a dead-region `skip`/`skipTest` or a function-body `skip` is mid-flight
    # reports `Unexpected end of file` and fails the parse; open *active*
    # `If`/`For` levels are silent (their last pass already streamed).
    for source in ("If (0)\nPoint(1)={0,0,0};",
                   "If (1)\nElse\nPoint(1)={0,0,0};",
                   "If (0)\nIf (1)\nPoint(1)={0,0,0};",
                   "Function f\nPoint(1)={0,0,0};")
        err=_struct_error(source)
        @test err isa ArgumentError
        @test occursin("Unexpected end of file",sprint(showerror,err))
    end

    # Active unclosed `If`/`For` — silent; the body still ran once.
    m1,_=_exec_ctx("If (1)\nPoint(1)={0,0,0};")
    @test sort!(collect(keys(m1.points)))==[1]
    m2,c2=_exec_ctx("For i In {0:2}\nPoint(1)={0,0,0};")
    @test sort!(collect(keys(m2.points)))==[1]
    @test c2.values["i"]==0

    # A `For` whose range never matches skips the body — the skip reports
    # the EOF like any dead region.
    err=_struct_error("For i In {3:1}\nPoint(1)={0,0,0};")
    @test err isa ArgumentError
    @test occursin("Unexpected end of file",sprint(showerror,err))

    # An unclosed `Function` still registered the name at its header
    # (upstream `createFunction` runs before the body skip) — a following
    # `Call` in the parent stream is not `Unknown function`... but the EOF
    # error itself is what `execute_geo` reports.
    err2=_struct_error("Function f\nPoint(1)={0,0,0};")
    @test err2 isa ArgumentError
    @test occursin("Unexpected end of file",sprint(showerror,err2))
end
