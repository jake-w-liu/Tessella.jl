# `.geo` statement-level error continuation — Gmsh 4.15.2 reports most
# statement failures as `Error   : ...` diagnostics and keeps executing the
# rest of the script (the bison `error tEND` recovery for syntax errors,
# `Msg::Error`/`yymsg(0)` for runtime failures). `execute_geo` mirrors that:
# recoverable diagnostics land on `GeoExecution.errors`/`msg_error_count`
# and the model carries every side effect upstream would have produced.
# Fatal input — unterminated statements/comments/strings, resource bounds —
# still rejects outright.
using Test
using Tessella

function _execute_geo_source(source::AbstractString;files=Dict{String,String}(),mesh_dim=0)
    dir=mktempdir()
    path=joinpath(dir,"script.geo")
    write(path,source)
    for (name,content) in files
        write(joinpath(dir,name),content)
    end
    return execute_geo(path;mesh_dim=mesh_dim)
end

@testset "unknown commands are recoverable diagnostics" begin
    # `Foobar;` reduces to a bare `tSTRING` — `syntax error (;)` upstream.
    ex=_execute_geo_source("Point(1)={0,0,0};\nFoobar;\nPrintf(\"tail\");")
    @test "syntax error (;)" in ex.errors
    @test haskey(ex.model.points,1)
    # `Foobar 1`/`Foobar(1)` evaluate the argument then report the command.
    ex=_execute_geo_source("Foobar 1;\nPoint(1)={0,0,0};")
    @test ex.errors==["Unknown command 'Foobar'"]
    @test haskey(ex.model.points,1)
    ex=_execute_geo_source("Foobar(1);\nPoint(1)={0,0,0};")
    @test ex.errors==["Unknown command 'Foobar'"]
    @test haskey(ex.model.points,1)
    ex=_execute_geo_source("SendToServer \"x\";\nPoint(1)={0,0,0};")
    @test ex.errors==["Unknown command 'SendToServer'"]
end

@testset "missing Merge/Include files are recoverable" begin
    ex=_execute_geo_source("""
        Point(1)={0,0,0};
        Merge "missing.pos";
        Include "nothere.geo";
        Point(2)={1,0,0};""")
    @test "Unable to open file 'missing.pos'" in ex.errors
    @test "Unable to open file 'nothere.geo'" in ex.errors
    @test haskey(ex.model.points,1) && haskey(ex.model.points,2)
    # An included file's syntax error reports in place and returns to the
    # caller's stream.
    ex=_execute_geo_source(
        "Include \"bad.geo\";\nPoint(2)={1,0,0};";
        files=Dict("bad.geo"=>"Point(1)={0,0,0};\n}};\nPoint(3)={2,0,0};\n"))
    @test haskey(ex.model.points,2)
end

@testset "transform/Duplicata/Boundary unknown entities are partial" begin
    # Upstream applies the transform to resolvable entities and reports
    # `Unknown GEO <kind> with tag <n>` + `Could not translate shapes`.
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Translate {1,0,0} { Point{9}; Curve{1}; }
        Point(3)={0,1,0};""")
    @test "Unknown GEO point with tag 9" in ex.errors
    @test "Could not translate shapes" in ex.errors
    @test haskey(ex.model.points,3)
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Duplicata { Point{9}; Curve{1}; }
        Point(30)={0,1,0};""")
    @test "Unknown GEO point 9" in ex.errors
    @test "Could not apply operation on shapes" in ex.errors
    # The valid entity was still duplicated.
    @test length(ex.model.curves)==2
    @test haskey(ex.model.points,30)
    # `Boundary` skips points silently and reports other unknown entities.
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Boundary { Point{9}; }
        Boundary { Curve{9}; }
        Point(3)={0,1,0};""")
    @test "Unknown model curve with tag 9" in ex.errors
    @test "Could not apply operation on shapes" in ex.errors
    @test !("Unknown model vertex with tag 9" in ex.errors)
    @test haskey(ex.model.points,3)
end

@testset "Extrude silently skips unresolved entities" begin
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={0,1,0};Point(4)={1,1,0};
        Line(1)={1,2};Line(2)={2,4};Line(3)={4,3};Line(4)={3,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Extrude {0,0,1} { Surface{1}; Surface{9}; }
        Point(100)={5,5,5};""")
    @test isempty(ex.errors)
    @test !isempty(ex.model.volumes)
    @test haskey(ex.model.points,100)
end

@testset "geometry construction failures are diagnostics" begin
    ex=_execute_geo_source("""
        Point(1)={0,0,0};
        Line(1)={1,99};
        Point(2)={1,0,0};""")
    @test any(err->occursin("Unknown control point 99",err),ex.errors)
    @test "Could not add line" in ex.errors
    @test haskey(ex.model.points,2)
end

@testset "Delete operations report and continue" begin
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Delete { Point{9}; Curve{1}; }
        Point(3)={0,1,0};""")
    @test "Unknown model vertex with tag 9" in ex.errors
    @test !haskey(ex.model.curves,1)
    @test haskey(ex.model.points,3)
    ex=_execute_geo_source("Delete View[5];\nPoint(1)={0,0,0};")
    @test ex.errors==["Unknown view 5"]
    ex=_execute_geo_source("Delete Point{9};\nPoint(1)={0,0,0};")
    @test ex.errors==["syntax error (Point)"]
    @test haskey(ex.model.points,1)
    ex=_execute_geo_source("Delete Foo Bar;\nPoint(1)={0,0,0};")
    @test ex.errors==["Unknown command 'Delete Foo Bar'"]
    ex=_execute_geo_source("Delete Widget[0];\nPoint(1)={0,0,0};")
    @test ex.errors==["Unknown command 'Delete Widget'"]
end

@testset "embeds report missing parents/members and continue" begin
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={1,1,0};Point(4)={0,1,0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={0.5,0.5,0};
        Point{5} In Surface{9};
        Point{9} In Surface{1};
        Point{5} In Surface{1};
        Line(5)={1,3};
        Point(6)={2,2,2};""")
    @test "Unknown model surface with tag 9" in ex.errors
    @test "Unknown model vertex with tag 9" in ex.errors
    @test get(ex.model.embeds,(2,1),[])==[(0,5)]
    @test haskey(ex.model.points,6)
end

@testset "plugins, fields and views" begin
    ex=_execute_geo_source("""
        Plugin(NoSuchPlugin).Run;
        Plugin(Integrate).NoSuchOption = 1;
        Plugin(Integrate).View;
        Point(1)={0,0,0};""")
    @test "Failed action 'Run' or unknown plugin 'NoSuchPlugin'" in ex.errors
    @test "Unknown option 'NoSuchOption' or plugin 'Integrate'" in ex.errors
    @test "Failed action 'View' or unknown plugin 'Integrate'" in ex.errors
    @test haskey(ex.model.points,1)
    ex=_execute_geo_source("""
        Field[1] = BogusField;
        Point(1)={0,0,0};""")
    @test "Unknown field type \"BogusField\"" in ex.errors
    @test "Cannot create field 1 of type 'BogusField'" in ex.errors
    ex=_execute_geo_source("""
        Save View[0] "x.pos";
        Print View[0] "y.pos";
        Alias View[0];
        Point(1)={0,0,0};""")
    @test "Unknown view 0" in ex.errors
    # `Print` is not a view verb upstream — `Unknown view` wins when the
    # index does not resolve.
    @test count(==("Unknown view 0"),ex.errors)==3
end

@testset "stray closing braces recover like upstream" begin
    # A `}` mid-expression kills the whole statement upstream — the bison
    # reduction needs `tEND`, so the point is never added, and the error
    # token is the stray brace itself.
    ex=_execute_geo_source("""
        Point(1)={0,0,0,1}};
        Point(2)={1,0,0,1};""")
    @test ex.errors==["syntax error (})"]
    @test !haskey(ex.model.points,1)
    @test haskey(ex.model.points,2)
    ex=_execute_geo_source("""
        Point(1)={0,0,0,1};
        } extra tokens;
        Point(2)={1,0,0,1};""")
    @test ex.errors==["syntax error (})"]
    @test haskey(ex.model.points,1) && haskey(ex.model.points,2)
end

@testset "mesh statements on empty inputs continue" begin
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Mesh 3;
        PartitionMesh(0);
        Delete Meshes;
        PartitionMesh(2);
        RefineMesh;
        Point(3)={0,1,0};""")
    @test isempty(ex.errors)
    @test "No mesh nodes were found" in ex.warnings
    @test haskey(ex.model.points,3)
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Mesh 1;
        SetOrder 0;
        MeshSize {9} = 0.1;
        MeshSize {1} = 0.2;
        Point(3)={0,1,0};""")
    @test isempty(ex.errors)
    @test ex.model.point_size[1]≈0.2
    @test !haskey(ex.model.point_size,9)
end

@testset "boolean and periodic statements report and continue" begin
    # Built-in factory: the OCC gate is a recoverable `yymsg(0)`.
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        BooleanUnion { Curve{1}; }{ Curve{1}; }
        Point(3)={0,1,0};""")
    @test "Boolean operators only available with OpenCASCADE " *
          "geometry kernel" in ex.errors
    @test haskey(ex.model.points,3)
    # `Periodic Volume`/`Periodic Point` are `syntax error` diagnostics.
    ex=_execute_geo_source("""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Periodic Volume {1} = {1} Translate {1,0,0};
        Point(3)={0,1,0};""")
    @test any(err->occursin("syntax error (Volume)",err),ex.errors)
    @test haskey(ex.model.points,3)
end

@testset "Save on odd targets continues" begin
    dir=mktempdir()
    geo=joinpath(dir,"t.geo")
    write(geo,"""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Save "x.bogus";
        Save "x.pos";
        Save "x.vtk";
        Point(3)={0,1,0};""")
    ex=execute_geo(geo)
    @test ex.errors==["Unknown output file format"]
    @test isfile(joinpath(dir,"x.pos"))
    @test isfile(joinpath(dir,"x.vtk"))
    @test haskey(ex.model.points,3)
end

@testset "state around failed statements is preserved" begin
    ex=_execute_geo_source("""
        a = 1;
        Point(1)={0,0,0};
        Line(9)={1,99};
        a = 2;
        Point(2)={1,0,0};
        Merge "missing.pos";
        a = 3;""")
    @test ex.values["a"]≈3.0
    @test haskey(ex.model.points,1) && haskey(ex.model.points,2)
    @test !isempty(ex.errors)
    # Errors inside a `For` body do not stop the loop.
    ex=_execute_geo_source("""
        For i In {1:3}
          Line(i)={i,99};
          Point(i+10)={i,0,0};
        EndFor
        Point(20)={0,0,0};""")
    @test haskey(ex.model.points,11) && haskey(ex.model.points,12) &&
          haskey(ex.model.points,13)
    @test haskey(ex.model.points,20)
end

@testset "fatal input still rejects" begin
    @test_throws ArgumentError _execute_geo_source("/* unterminated\n")
    @test_throws ArgumentError _execute_geo_source(
        "Point(1)={0,0,0\n")
    @test_throws ArgumentError _execute_geo_source(
        "Point(1)={0,0,0}; stray text")
    @test_throws ArgumentError _execute_geo_source(
        "Point(1)={0,0,0}; \"unterminated")
end

@testset "Abort stops the stream without an error" begin
    # `Abort` is a voluntary stop (upstream `Msg::Exit`): nothing after it
    # runs, and the execution still returns.
    ex=_execute_geo_source("""
        Point(1)={0,0,0};
        Abort;
        Point(2)={1,0,0};""")
    @test haskey(ex.model.points,1)
    @test !haskey(ex.model.points,2)
end

@testset "Exit semantics are preserved" begin
    ex=_execute_geo_source("Exit 0;\nPoint(1)={0,0,0};")
    @test ex.exit_code==0
    @test !haskey(ex.model.points,1)
    # Errors recorded before `Exit 0` are dropped like upstream's
    # `exit(0)` leaving before the errorstate check.
    ex=_execute_geo_source("Foobar 1;\nExit 0;")
    @test ex.exit_code==0
    @test_throws ArgumentError _execute_geo_source("Exit 2;")
end
