# `DefineConstant[...]` entry grammar — Gmsh 4.15.2 accepts comma-separated
# AND whitespace/newline-separated `name = spec` entries (t20/t21 write the
# latter). An entry ends where a new entry can start: an identifier directly
# following a completed atom (name, literal, `)`, `]`, `}` or string close)
# at bracket depth zero.
using Test
using Tessella

function _execute_geo_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _point_coords(result)
    points=result.model.points
    return sort!(collect(values(points)))
end

@testset "DefineConstant separators (Gmsh 4.15.2 grammar)" begin
    # Comma-separated entries — the long-standing form.
    r=_execute_geo_source("""
        DefineConstant[A = 1, B = 2, C = 3];
        Point(1) = {A, B, C, 1};""")
    @test _point_coords(r)==[(1.0,2.0,3.0)]

    # Newline-separated annotated entries — verbatim t21 style. Commas live
    # inside the braces; the entries themselves carry none.
    r=_execute_geo_source("""
        DefineConstant[
          partitioner = {0, Choices{0="Metis", 1="SimplePartition"},
            Name "Parameters/0Mesh partitioner"}
          N = {3, Min 1, Max 256, Step 1,
            Name "Parameters/1Number of partitions"}
          topology = {1, Choices{0, 1},
            Name "Parameters/2Create partition topology (BRep)?"}
          ghosts = {0, Choices{0, 1},
            Name "Parameters/3Create ghost cells?"}
        ];
        Point(1) = {N, topology, ghosts, 1};""")
    @test _point_coords(r)==[(3.0,1.0,0.0)]

    # Comma-less entries on one line — whitespace alone separates.
    r=_execute_geo_source("""
        DefineConstant[A = 1 B = 2];
        Point(1) = {A, B, 0, 1};""")
    @test _point_coords(r)==[(1.0,2.0,0.0)]

    # Bare-name entries and an `=`-expression boundary mid-line:
    # `C` and `E` parse as separate bare entries (upstream reads `D = A + 1`
    # then `E`, not a longer expression).
    r=_execute_geo_source("""
        DefineConstant[A = 1 B = 2, C D = A + 1 E];
        Point(1) = {A, B, C, 1};
        Point(2) = {D, E, 0, 1};""")
    @test _point_coords(r)==[(1.0,2.0,0.0),(2.0,0.0,0.0)]

    # Whitespace-separated bare names.
    r=_execute_geo_source("""
        DefineConstant[X Y];
        Point(1) = {X, Y, 0, 1};""")
    @test _point_coords(r)==[(0.0,0.0,0.0)]

    # A `}`-closed annotation followed directly by the next entry.
    r=_execute_geo_source("""
        DefineConstant[A = {4, Name "a"} B = {5, Name "b"}];
        Point(1) = {A, B, 0, 1};""")
    @test _point_coords(r)==[(4.0,5.0,0.0)]

    # String constants keep the brace options in place.
    r=_execute_geo_source("""
        DefineConstant[txt = {"hello", Name "greeting"} num = 7];
        Point(1) = {num, 0, 0, 1};""")
    @test _point_coords(r)==[(7.0,0.0,0.0)]

    # `UndefineConstant` shares the whitespace-separated name grammar.
    r=_execute_geo_source("""
        DefineConstant[A = 1 B = 2];
        UndefineConstant[A B];
        Point(1) = {A, B, 0, 1};""")
    @test _point_coords(r)==[(1.0,2.0,0.0)]
end
