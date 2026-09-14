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

    for source in ("EndIf",
                   "Else",
                   "ElseIf (1)",
                   "If (1)\nPoint(1)={0,0,0,1};",
                   "If (1)\nElse\nElse\nPoint(1)={0,0,0,1};\nEndIf",
                   "If (1)\nElse\nElseIf (1)\nPoint(1)={0,0,0,1};\nEndIf",
                   "If (1)\nFor i In {0:1}\nEndIf\nEndFor",
                   "If (1)\nEndWhile\nEndIf",
                   "If 1)\nPoint(1)={0,0,0,1};\nEndIf",
                   "If (1\nPoint(1)={0,0,0,1};\nEndIf")
        @test _control_error(source) isa ArgumentError
    end
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

    for source in ("EndFor",
                   "For i In {0:2}\nPoint(1)={0,0,0,1};",
                   "For i In {1,2,3}\nEndFor",
                   "For i In {5}\nEndFor",
                   "For (i = 0; i < 3; i++)\nEndFor",
                   "For i In {0:5:0}\nEndFor",
                   "For i {0:2}\nEndFor",
                   "For In {0:2}\nEndFor",
                   "For Pi In {0:1}\nEndFor",
                   "For newp In {0:1}\nEndFor",
                   "For i In {0:1\nEndFor",
                   "For i In {0:2}\nElse\nEndFor")
        @test _control_error(source) isa ArgumentError
    end
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
    for source in ("EndWhile",
                   "While (1)\nPoint(1)={0,0,0,1};",
                   "While 1)\nEndWhile",
                   "While (1)\nEndFor\nEndWhile")
        @test _control_error(source) isa ArgumentError
    end
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
                   "x = 1 & 2;",
                   "x = 1 | 2;",
                   "x = 1 = 2;",
                   "x = 1 ? 2 ? 3 : 4;")
        @test _control_error(source) isa ArgumentError
    end
end
