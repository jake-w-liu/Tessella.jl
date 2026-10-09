using Test
using Tessella

function _execute_control_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _control_error(source::AbstractString)
    try
        _execute_control_source(source)
        return nothing
    catch err
        return err
    end
end

@testset "bounded .geo If/ElseIf/Else/EndIf" begin
    parsed=_execute_control_source(raw"""
        x = 1;
        If (x > 0)
          Point(50) = {5, 0, 0, 1};
        ElseIf (x > 2)
          Point(51) = {6, 0, 0, 1};
        Else
          Point(52) = {7, 0, 0, 1};
        EndIf
        If (x < 0)
          Point(53) = {0, 0, 0, 1};
        ElseIf (x == 1)
          Point(54) = {54, 0, 0, 1};
        Else
          Point(55) = {0, 0, 0, 1};
        EndIf
        """)
    # Matches the pinned-Gmsh unrolled file: only the true branches execute.
    @test sort!(collect(keys(parsed.model.points)))==[50,54]
    @test parsed.model.points[50]==(5.0,0.0,0.0)
    @test parsed.model.points[54]==(54.0,0.0,0.0)

    single_line=_execute_control_source(
        "If (0) Point(70) = {0,0,0,1}; Else Point(71) = {1,0,0,1}; EndIf")
    @test sort!(collect(keys(single_line.model.points)))==[71]

    nested=_execute_control_source(raw"""
        If (1)
          If (0)
            Point(1) = {0,0,0,1};
          Else
            Point(2) = {2,0,0,1};
          EndIf
        EndIf
        Point(3) = {3,0,0,1};
        """)
    @test sort!(collect(keys(nested.model.points)))==[2,3]

    # Orphan `Else`/`ElseIf` are `yymsg(0)` errors; an orphan `EndIf` is only
    # a `yymsg(1)` warning upstream (Gmsh.y:4094-4099) and does not fail.
    for source in ("Else",
                   "ElseIf (1)",
                   "If 1)\nPoint(1)={0,0,0,1};\nEndIf",
                   "If (1\nPoint(1)={0,0,0,1};\nEndIf")
        @test _control_error(source) isa ArgumentError
    end
    orphan=_execute_control_source("EndIf")
    @test orphan.warnings==["Orphan EndIf"]

    # Upstream's `Else`/`ElseIf` on an already-taken `If` runs
    # `skip("If","EndIf")`: everything through the matching `EndIf` is
    # invisible text — later `Else`/`ElseIf` headers and their bodies do not
    # parse at all, and the level closes without diagnostics.
    for source in ("If (1)\nElse\nElse\nPoint(1)={0,0,0,1};\nEndIf",
                   "If (1)\nElse\nElseIf (1)\nPoint(1)={0,0,0,1};\nEndIf")
        silent=_execute_control_source(source)
        @test isempty(silent.model.points)
        @test isempty(silent.warnings)
    end
    # An `Else` that follows an untaken chain does not set the taken flag —
    # a later `ElseIf` can still activate (Gmsh.y:4082-4093 leaves status 0).
    else_elseif=_execute_control_source(raw"""
        If (0)
          Point(1) = {0,0,0,1};
        Else
          Point(2) = {2,0,0,1};
        ElseIf (1)
          Point(3) = {3,0,0,1};
        EndIf
        """)
    @test sort!(collect(keys(else_elseif.model.points)))==[2,3]

    # Control families interleave through the flat stream: an `EndIf` inside
    # a `For` body decrements the outer `ImbricatedTest` on every iteration —
    # the second one drives it negative (`Orphan EndIf` warning) while the
    # loop itself completes silently (verified against pinned Gmsh 4.15.2).
    interleaved=_execute_control_source(raw"""
        If (1)
        For i In {0:1}
        EndIf
        EndFor
        Point(1) = {1,0,0,1};
        """)
    @test haskey(interleaved.model.points,1)
    @test interleaved.warnings==["Orphan EndIf"]
    # An `Else` inside the body `skip`s through the `EndIf` that belongs to
    # the outer `If` — the loop frame stays open and finishes normally.
    interleaved2=_execute_control_source(raw"""
        If (1)
        For i In {0:1}
        Else
          Point(9) = {9,0,0,1};
        EndFor
        EndIf
        Printf("done %g", i);
        """)
    @test !haskey(interleaved2.model.points,9)
    @test interleaved2.values["i"]==0.0

    # `EndWhile` with no open While is a silent marker in the bounded
    # extension (upstream lexes `EndWhile` as a plain identifier).
    noerr=_execute_control_source("If (1)\nEndWhile\nEndIf\nPoint(1)={0,0,0,1};")
    @test haskey(noerr.model.points,1)

    # Upstream runs an unclosed taken `If` to end-of-file without complaint
    # (Gmsh 4.15.2): the trailing statements execute normally.
    unclosed=_execute_control_source("If (1)\nPoint(1)={0,0,0,1};")
    @test haskey(unclosed.model.points,1)
    # An unclosed *untaken* `If` dies in `skipTest` — `Unexpected end of
    # file` fails the parse upstream (gmsh exits 1) without executing the
    # body.
    dead=_control_error("If (0)\nPoint(1)={0,0,0,1};")
    @test dead isa ArgumentError
    @test occursin("Unexpected end of file",sprint(showerror,dead))
end

@testset "bounded .geo For/EndFor" begin
    parsed=_execute_control_source(raw"""
        For i In {0:2}
          Point(i+1) = {i, 0, 0, 1};
        EndFor
        For k In {0:4:2}
          Point(k+10) = {k, 0, 0, 1};
        EndFor
        """)
    # Identical to the pinned-Gmsh unrolled file for the same script.
    @test sort!(collect(keys(parsed.model.points)))==[1,2,3,10,12,14]
    @test parsed.model.points[3]==(2.0,0.0,0.0)
    @test parsed.model.points[14]==(4.0,0.0,0.0)

    # The implicit two-term increment is +1: {5:0} iterates zero times, and the
    # loop variable still lands on the (unexecuted) start value.
    descending=_execute_control_source(raw"""
        For k In {5:0}
          Point(80+k) = {k,0,0,1};
        EndFor
        Point(90) = {k,0,0,1};
        """)
    @test sort!(collect(keys(descending.model.points)))==[90]
    @test descending.model.points[90]==(5.0,0.0,0.0)

    # Explicit negative increments iterate down; the variable persists at the
    # first out-of-range value, as in Gmsh's unrolled output.
    explicit=_execute_control_source(raw"""
        i = 7;
        For i In {0:1}
          Point(i+1) = {i,0,0,1};
        EndFor
        Point(99) = {i,0,0,1};
        For k In {5:0:-1}
          Point(80+k) = {k,0,0,1};
        EndFor
        Point(91) = {k,0,0,1};
        """)
    @test sort!(collect(keys(explicit.model.points)))==[1,2,80,81,82,83,84,85,91,99]
    @test explicit.model.points[99]==(2.0,0.0,0.0)
    @test explicit.model.points[91]==(-1.0,0.0,0.0)

    # Range bounds are evaluated once at loop entry.
    once=_execute_control_source(raw"""
        n = 2;
        For i In {0:n}
          n = 100;
          Point(i+1) = {i,0,0,1};
        EndFor
        """)
    @test sort!(collect(keys(once.model.points)))==[1,2,3]

    # Loop variables can drive entity tags and allocator reads per iteration.
    alloc=_execute_control_source(raw"""
        For i In {0:2}
          Point(newp) = {i,0,0,1};
        EndFor
        Point(newp) = {9,0,0,1};
        """)
    @test sort!(collect(keys(alloc.model.points)))==[1,2,3,4]

    nested=_execute_control_source(raw"""
        For i In {0:1}
          If (i == 0)
            For j In {10:11}
              Point(i*100+j) = {j,0,0,1};
            EndFor
          Else
            Point(9) = {9,0,0,1};
          EndIf
        EndFor
        """)
    @test sort!(collect(keys(nested.model.points)))==[9,10,11]

    # A zero step is legal upstream (Gmsh.y:3943-3963): the body runs once
    # and the `EndFor` range check pops the level immediately.
    zerostep=_execute_control_source(raw"""
        For i In {0:5:0}
          Point(i+1) = {i,0,0,1};
        EndFor
        Point(99) = {i,0,0,1};
        """)
    @test sort!(collect(keys(zerostep.model.points)))==[1,99]
    @test zerostep.model.points[99]==(0.0,0.0,0.0)

    # The anonymous `For (a:b[:c])` range form exists upstream (Gmsh.y:3904)
    # and iterates without binding a variable.
    anon=_execute_control_source(raw"""
        For (0:1)
          Point(newp) = {0,0,0,1};
        EndFor
        """)
    @test sort!(collect(keys(anon.model.points)))==[1,2]

    for source in ("EndFor",
                   "For i In {1,2,3}\nEndFor",
                   "For i In {5}\nEndFor",
                   "For (i = 0; i < 3; i++)\nEndFor",
                   "For i {0:2}\nEndFor",
                   "For In {0:2}\nEndFor",
                   "For Pi In {0:1}\nEndFor",
                   "For newp In {0:1}\nEndFor",
                   "For i In {0:1\nEndFor",
                   "For i In {0:2}\nElse\nEndFor")
        @test _control_error(source) isa ArgumentError
    end

    # Upstream runs an unclosed `For` body once, to end-of-file, without
    # complaint (Gmsh 4.15.2): the loop variable is consumed and the trailing
    # statements execute a single time.
    unclosed=_execute_control_source("For i In {0:2}\nPoint(i+1)={i,0,0,1};")
    @test haskey(unclosed.model.points,1)
    @test !haskey(unclosed.model.points,2)

    # `EndFor` advances the loop variable through the live symbol, so a
    # body's own write to `i` propagates into the next-iteration check.
    mutated=_execute_control_source(raw"""
        For i In {0:4}
          i = i + 3;
          Point(i+1) = {i,0,0,1};
        EndFor
        """)
    @test sort!(collect(keys(mutated.model.points)))==[4,8]
end

@testset "bounded .geo While/EndWhile (Tessella extension)" begin
    # Pinned Gmsh 4.15.2 has no While keyword; this is a bounded extension.
    parsed=_execute_control_source(raw"""
        j = 0;
        While (j < 3)
          Point(60+j) = {j*10, 0, 0, 1};
          j = j + 1;
        EndWhile
        """)
    @test sort!(collect(keys(parsed.model.points)))==[60,61,62]
    @test parsed.model.points[62]==(20.0,0.0,0.0)

    false_entry=_execute_control_source(raw"""
        While (0)
          Point(1) = {0,0,0,1};
        EndWhile
        Point(2) = {2,0,0,1};
        """)
    @test sort!(collect(keys(false_entry.model.points)))==[2]

    nested=_execute_control_source(raw"""
        For i In {0:1}
          j = i;
          While (j < i + 2)
            Point(i*10+j+10) = {j,0,0,1};
            j = j + 1;
          EndWhile
        EndFor
        """)
    @test sort!(collect(keys(nested.model.points)))==[10,11,21,22]

    cap=_control_error("While (1)\nEndWhile")
    @test cap isa ArgumentError
    @test occursin("iterations",sprint(showerror,cap))
    # A bare `EndWhile` with no open While is a silent marker (upstream has
    # no While keyword — it lexes the word as a plain identifier).
    @test _control_error("EndWhile")===nothing
    # An unclosed `While` whose header condition held streams its body once
    # and stays silent — the extension mirrors upstream's `For` EOF rule.
    unclosed=_execute_control_source("While (1)\nPoint(1)={0,0,0,1};")
    @test sort!(collect(keys(unclosed.model.points)))==[1]
    for source in ("While (0)\nPoint(1)={0,0,0,1};",
                   "While 1)\nEndWhile",
                   "While (1)\nEndFor\nEndWhile")
        @test _control_error(source) isa ArgumentError
    end
end

@testset "bounded .geo Function/Call/Return" begin
    # Gmsh 4.15.2: `Function name ... Return` registers a zero-argument body
    # executed by `Call name;` in the shared variable scope; each call
    # re-executes the body.
    parsed=_execute_control_source(raw"""
        Function makeLine
          p = newp;
          Point(p) = {0, 0, 0};
          Point(p+1) = {1, 0, 0};
          Line(newl) = {p, p+1};
        Return
        Call makeLine;
        x = 42;
        Call makeLine;
        Point(99) = {x, 0, 0};
        """)
    @test sort!(collect(keys(parsed.model.points)))==[1,2,3,4,99]
    @test sort!(collect(keys(parsed.model.curves)))==[1,2]
    @test parsed.model.points[99]==(42.0,0.0,0.0)

    # `Call` is legal inside If/For blocks and re-runs per iteration.
    loop=_execute_control_source(raw"""
        Function mark
          Point(newp) = {0, 0, 0};
        Return
        If (1)
          Call mark;
        EndIf
        For i In {0:2}
          Call mark;
        EndFor
        """)
    @test sort!(collect(keys(loop.model.points)))==[1,2,3,4]

    # Quoted string-literal names match Gmsh's string-expression header.
    quoted=_execute_control_source(raw"""
        Function "quoted name"
          Point(7) = {0, 0, 0};
        Return
        Call "quoted name";
        """)
    @test sort!(collect(keys(quoted.model.points)))==[7]

    # Functions compose: a body may `Call` an already-registered function.
    # (Gmsh's token-level capture ends a `Function` body at the first `Return`,
    # so a `Function` inside another body always steals the outer body's
    # terminator — a nested definition is not a useful form.)
    nested=_execute_control_source(raw"""
        Function inner
          Point(1) = {0, 0, 0};
        Return
        Function outer
          Call inner;
          Point(2) = {1, 0, 0};
        Return
        Call outer;
        """)
    @test sort!(collect(keys(nested.model.points)))==[1,2]

    # A `Function` inside an If block registers only when the branch runs.
    branch=_execute_control_source(raw"""
        If (1)
          Function g
            Point(5) = {0, 0, 0};
          Return
        EndIf
        Call g;
        """)
    @test sort!(collect(keys(branch.model.points)))==[5]

    cond=_execute_control_source(raw"""
        If (0)
          Function hidden
            Point(1) = {0, 0, 0};
          Return
        EndIf
        """)
    @test isempty(cond.model.points)
    @test _control_error(raw"""
        If (0)
          Function hidden
            Point(1) = {0, 0, 0};
          Return
        EndIf
        Call hidden;
        """) isa ArgumentError

    err=_control_error("Call makePoint;\nFunction makePoint\n  Point(1)={0,0,0};\nReturn\n")
    @test err isa ArgumentError
    @test occursin("Unknown function 'makePoint'",sprint(showerror,err))
    err=_control_error("Function f\nReturn\nFunction f\nReturn\n")
    @test err isa ArgumentError
    @test occursin("Redefinition of function f",sprint(showerror,err))
    err=_control_error("Function r\n  Call r;\nReturn\nCall r;\n")
    @test err isa ArgumentError
    @test occursin("Call depth",sprint(showerror,err))
    err=_control_error("Point(1)={0,0,0};\nReturn\n")
    @test err isa ArgumentError
    @test occursin("Error while exiting function",sprint(showerror,err))
    err=_control_error("Function f\n  Point(1)={0,0,0};\n")
    @test err isa ArgumentError
    @test occursin("Unexpected end of file",sprint(showerror,err))
    for source in ("Function f\nReturn\nCall f(1,2);",
                   "Function\nReturn",
                   "Call;",
                   "Call f;",
                   "Function 7f\nReturn",
                   "Function f\n  If (1)\n    Return;\n  EndIf\nReturn\nCall f;")
        @test _control_error(source) isa ArgumentError
    end

    # A `Function` definition inside a `For` body re-registers every
    # iteration: the second registration reports `Redefinition` and keeps
    # the original body (Gmsh.y:4000-4006), verified against pinned Gmsh.
    loopdef=_control_error(raw"""
        For i In {0:1}
          Function f
            Point(newp) = {0,0,0};
          Return
        EndFor
        """)
    @test loopdef isa ArgumentError
    @test occursin("Redefinition of function f",sprint(showerror,loopdef))
    # `Call` before the definition reports `Unknown function` on every
    # iteration — the body registered later is never seen by earlier calls.
    callfirst=_control_error(raw"""
        For i In {0:1}
          Call g;
        EndFor
        Function g
          Point(1) = {0,0,0};
        Return
        """)
    @test callfirst isa ArgumentError
    @test occursin("Unknown function 'g'",sprint(showerror,callfirst))
end

@testset "bounded .geo error tEND recovery" begin
    # A syntax error in a control header triggers upstream's `error tEND`
    # recovery (Gmsh.y:276): the discard eats following control keywords and
    # the first `;`-terminated statement, then parsing resumes — `c` still
    # runs and the `EndIf` reports `Orphan EndIf` as a *warning*, so the only
    # recorded error is the header's own `syntax error (y)`.
    err=_control_error(raw"""
        If (1 y)
          b = 5;
          c = 6;
        EndIf
        """)
    @test err isa ArgumentError
    @test sprint(showerror,err)=="ArgumentError: execute_geo: syntax error (y)"
    # When the discard swallows the `EndIf` itself, no orphan diagnostic
    # follows — the first `;`-statement ends the recovery.
    err2=_control_error(raw"""
        If (x y)
        EndIf
        Printf("done");
        """)
    @test err2 isa ArgumentError
    @test occursin("syntax error (y)",sprint(showerror,err2))
    @test !occursin("Orphan",sprint(showerror,err2))
    # A malformed `For` range errors at the offending token, the discard
    # eats the body statement, and the stray `EndFor` then reports
    # `Invalid For/EndFor loop`.
    err3=_control_error("For i In {missing}\nPoint(1)={0,0,0};\nEndFor")
    @test err3 isa ArgumentError
    @test occursin("Unknown variable 'missing'",sprint(showerror,err3))
    @test occursin("syntax error (})",sprint(showerror,err3))
    @test occursin("Invalid For/EndFor loop",sprint(showerror,err3))
    @test !occursin("Unexpected end of file",sprint(showerror,err3))
end

@testset "bounded .geo comparison, logical, and ternary operators" begin
    parsed=_execute_control_source(raw"""
        x = 1;
        a = (x > 0) ? 42 : -1;
        b = !(x == 1);
        c = (x >= 1) && (x < 2) || 0;
        d = 1 < 2 == 1;
        e = x != 2;
        Point(1) = {a, b, c, 1};
        Point(2) = {d, e, (0 ? 9 : 3), 1};
        """)
    @test parsed.model.points[1]==(42.0,0.0,1.0)
    @test parsed.model.points[2]==(1.0,1.0,3.0)

    # A ':' answering a pending '?' inside a brace list binds to the ternary,
    # not to a range — matching the pinned-Gmsh unrolled output.
    inline=_execute_control_source(raw"""
        Point(1) = {1 ? 2 : 3, 0, 0, 1};
        Point(2) = {0 ? 9 : 3, 0, 0, 1};
        a = 2;
        Point(3) = {a > 1 ? a : -a, 0, 0, 1};
        """)
    @test inline.model.points[1]==(2.0,0.0,0.0)
    @test inline.model.points[2]==(3.0,0.0,0.0)
    @test inline.model.points[3]==(2.0,0.0,0.0)

    mixed=_execute_control_source(
        "If (2 > 1 && 3 <= 3) Point(1) = {1,0,0,1}; EndIf")
    @test haskey(mixed.model.points,1)

    for source in ("x = 1 ? 2;",
                   "x = 1 ? 2 : ;",
                   "x = 1 = 2;",
                   "x = 1 ? 2 ? 3 : 4;")
        @test _control_error(source) isa ArgumentError
    end
    # `&`/`|` are integer-truncating bitwise operators upstream (Gmsh 4.15.2:
    # `5 & 3.7` = 1, `0.5 | 2` = 2).
    bits=_execute_control_source("x = 6 & 3; y = 4 | 1; z = 5 & 3.7; w = 0.5 | 2;")
    @test bits.values["x"]==2.0 && bits.values["y"]==5.0
    @test bits.values["z"]==1.0 && bits.values["w"]==2.0

    # Bitwise/shift operands and counts take the shipped binary's x86
    # semantics (Gmsh 4.15.2 verified): `cvttsd2si` truncation, the
    # 0x80000000 "integer indefinite" value on out-of-range/non-finite
    # operands, and a 5-bit masked shift count.
    hw=_execute_control_source(
        "a = 6 >> 255; b = 1 << 32; c = 255 & 1e30; d = -3 | 2; " *
        "e = 5 >> Sqrt(-1); f = 1 << -1; g = 6.9 >> 1.9; h = 1e30 >> 1; " *
        "i = 4294967295 & -1; j = -2147483648 >> 31; " *
        "k = (Pi - 12) ^ +0.25 % 100; l = -7 % 3; m = 7 % -3;")
    @test hw.values["a"]==0.0      # 255 & 31 = 31; 6 >> 31 = 0
    @test hw.values["b"]==1.0      # 32 & 31 = 0; 1 << 0 = 1
    @test hw.values["c"]==0.0      # int-indefinite & 0xFF = 0
    @test hw.values["d"]==-1.0     # -3 | 2 = -1
    @test hw.values["e"]==5.0      # NaN count -> 0x80000000 & 31 = 0
    @test hw.values["f"]==-2147483648.0  # -1 & 31 = 31; 1 << 31 wraps
    @test hw.values["g"]==3.0      # trunc(6.9) >> trunc(1.9) = 6 >> 1
    @test hw.values["h"]==-1073741824.0  # int-indefinite >> 1
    @test hw.values["i"]==-2147483648.0  # out-of-range -> int-indefinite
    @test hw.values["j"]==-1.0     # arithmetic right shift
    @test hw.values["k"]==-48.0    # int-indefinite % 100 = -48
    @test hw.values["l"]==-1.0     # C remainder: sign follows dividend
    @test hw.values["m"]==1.0
    # `% 0` and `INT_MIN % -1` trap the x86 `idiv` upstream — the process
    # dies — so an explicit diagnostic is the faithful equivalent here.
    @test _control_error("x = 5 % 0;") isa ArgumentError
    @test _control_error("x = -2147483648 % -1;") isa ArgumentError

    # `^` edge semantics (Gmsh 4.15.2 verified): `x^0 = 1` even for NaN x;
    # `(±1)^±Inf = 1`; integer exponents use binary squaring — parity flips
    # the sign of negative bases and ±Inf bases (odd → -Inf/-0, even →
    # +Inf/+0) — while integer exponents past the Int64 range are all even
    # and square `|x|`; positive overflow of the y·log2(x) chain is `Inf`.
    pe=_execute_control_source(
        "a = Sqrt(-1) ^ 0; b = (-1) ^ Exp(999); c = 1 ^ (-Exp(999)); " *
        "d = (-Exp(999)) ^ 0.5; e = (-Exp(999)) ^ 3; " *
        "f = (-Exp(999)) ^ (-3); g = (-Exp(999)) ^ (-2); " *
        "h = (-1) ^ 1e30; i = (-2) ^ 1e30; j = 2 ^ 1e19; k = 0.5 ^ (-1e30); " *
        "l = (-1) ^ (4.5e15 + 1); m = (-1) ^ 4.5e15; n = 2 ^ (-1074);")
    @test pe.values["a"]==1.0        # NaN^0 = 1
    @test pe.values["b"]==1.0        # (-1)^Inf = 1
    @test pe.values["c"]==1.0        # 1^(-Inf) = 1
    @test pe.values["d"]==Inf        # (-Inf)^non-integer = +Inf
    @test pe.values["e"]==-Inf       # odd integer parity
    @test pe.values["f"]===-0.0      # odd negative parity
    @test pe.values["g"]==0.0 && !signbit(pe.values["g"])  # even: +0
    @test pe.values["h"]==1.0        # |y|≥2^63 integer → even → |x|^y
    @test pe.values["i"]==Inf        # (-2)^1e30 = +Inf
    @test pe.values["j"]==Inf        # positive chain overflow → Inf
    @test pe.values["k"]==Inf        # 0.5^(-1e30) overflows + 
    @test pe.values["l"]==-1.0       # odd squaring parity on -1
    @test pe.values["m"]==1.0        # even squaring parity on -1
    @test pe.values["n"]==5e-324     # tracked-exponent subnormal squaring

    # Finite integer exponents still exceed a machine Int's accumulated
    # binary scaling range. Exact powers of two prove the overflow direction;
    # the pinned Windows Gmsh parser independently returns Inf/0/0/Inf.
    large_powers=_execute_control_source(
        "a = (2 ^ 58) ^ 3425408785282518016; " *
        "b = (2 ^ 58) ^ (-3425408785282518016); " *
        "c = (2 ^ (-58)) ^ 3425408785282518016; " *
        "d = (2 ^ (-58)) ^ (-3425408785282518016); " *
        "e = (2 ^ 512) ^ (-2); f = (2 ^ 537) ^ (-2);")
    @test large_powers.values["a"]===Inf
    @test large_powers.values["b"]===0.0
    @test large_powers.values["c"]===0.0
    @test large_powers.values["d"]===Inf
    @test large_powers.values["e"]===ldexp(1.0,-1024)
    @test large_powers.values["f"]===ldexp(1.0,-1074)

    # Trig past the x87 FPU range (|x| ≥ 2^63) takes mingw's assembly
    # fallback — `fldpi; fadd; fprem1` — reducing modulo twice π_hw rounded
    # at the current hardware precision, unlike in-range `fsin`/`fcos`.
    # All values verified bit-for-bit against the Windows gmsh.exe; other
    # platforms link the system libm upstream, so the pins are Windows-only.
    # These CLI pins use the standalone executable's PC64 precision. Direct
    # DLL/API calls instead retain the caller's precision and rounding.
    if Sys.iswindows()
        x87_source=
            "a = Sin(1e300); b = Cos(1e300); c = Sin(9.5e18); " *
            "d = Cos(9.5e18); e = Sin(-9.3e18); f = Cos(-9.3e18); " *
            "g = Sin(1e19); h = Cos(1e19); " *
            "i = Sin(4.5e15); j = Cos(4.5e15);"
        gmsh_math=Tessella.GmshLibm
        x87_cli_context=gmsh_math._win_cli_precision()
        x87_control_word=Sys.ARCH===:x86_64 ? gmsh_math._win_x87_control_word() : nothing
        x87=gmsh_math._with_win_cli_precision(() -> _execute_control_source(x87_source))
        @test gmsh_math._win_cli_precision()===x87_cli_context
        if Sys.ARCH===:x86_64
            @test gmsh_math._win_x87_control_word()===x87_control_word
        end
        @test x87.values["a"]===0.9790015909522538
        @test x87.values["b"]===0.20385260585274806
        @test x87.values["c"]===-0.9819889907942527
        @test x87.values["d"]===0.18893814320799576
        @test x87.values["e"]===0.7158654763735572
        @test x87.values["f"]===-0.6982382256339594
        @test x87.values["g"]===-0.8556574595436565
        @test x87.values["h"]===-0.517542570159495
        @test x87.values["i"]===0.06881493122288082
        @test x87.values["j"]===0.9976294428497939

        if Sys.ARCH===:x86_64
            # Literal bits independently captured from the pinned Gmsh 4.15.2
            # DLL parser at PC53/nearest for this exact program. The first eight
            # differ from the standalone PC64 executable; i/j are controls.
            api53_bits=(
                0xbfe52ec076f4df9b,0x3fe7fc3f18dde1db,
                0xbfeff980b1087879,0xbfa46346948342ef,
                0xbfe9ffde57a992ae,0xbfe2a7cd22d091a2,
                0xbfefc07dab441d63,0xbfbfd0b0ce36e482,
                0x3fb19ddaf71429c6,0x3fefec9494d2214f,
            )
            api53_word=(x87_control_word & UInt16(0xf0ff)) | UInt16(0x0200)
            try
                gmsh_math._win_x87_set_control_word!(api53_word)
                api53=_execute_control_source(x87_source)
                @test gmsh_math._win_x87_control_word()===api53_word
                @test gmsh_math._win_cli_precision()===x87_cli_context
                for (name,expected) in zip(("a","b","c","d","e","f","g","h","i","j"),api53_bits)
                    @test reinterpret(UInt64,api53.values[name])===expected
                end
            finally
                gmsh_math._win_x87_set_control_word!(x87_control_word)
            end
            @test gmsh_math._win_x87_control_word()===x87_control_word
            @test gmsh_math._win_cli_precision()===x87_cli_context
        end
    end
end
