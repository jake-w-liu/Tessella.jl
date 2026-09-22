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
end
