# `.geo` Boolean operand groups — the upstream
# `'{' ListOfShapes BooleanOption '}'` production (Gmsh.y `Boolean`/
# `BooleanShape`). Operand members are `ListOfShapes` items of any
# dimension (`Point`/`Curve`/`Line`/`Surface`/`Volume` blocks, `Line` as
# the `Curve` lexer alias, `GeoEntity{d}` clauses, `Physical`/`Parent`
# selectors, `:` wildcards and inline `Kind(tag) = ...` definitions), then
# the trailing `Delete`/`Recursive Delete` option with an optional `FExpr`
# gate. Operand lists reduce — side effects included — before the
# OpenCASCADE factory gate upstream, so resolution happens even under
# built-in. All-volume operations keep the native arrangement path; any
# sub-volume operand routes through the OCC kernel bridge like the
# `model.occ.*` boolean endpoints. Operand dims the bounded kernel cannot
# map (0/1) surface the endpoint's own message plus the recoverable
# `Could not apply boolean operator`, and the statement stream continues.

using Test
using Tessella

function _execute_geo_source(source::AbstractString;mesh_dim=0)
    dir=mktempdir()
    path=joinpath(dir,"script.geo")
    write(path,source)
    return execute_geo(path;mesh_dim=mesh_dim)
end

# A unit box plus a `1×1` wall exactly covering its mid-plane — the
# `occ.fragment` volume-by-surface fixture from the API suite.
const _BOX_WALL="""
    SetFactory("OpenCASCADE");
    Box(1) = {0,0,0, 1,1,1};
    Rectangle(7) = {0,0,0.5, 1,1};
"""

@testset "volume fragmented by a surface wall" begin
    ex=_execute_geo_source(_BOX_WALL * """
        BooleanFragments{ Volume{1}; Delete; }{ Surface{7}; Delete; }
        """)
    @test isempty(ex.errors)
    # two half-volumes; the wall rebounds at its own tag
    @test sort(collect(keys(ex.model.volumes)))==[1,2]
    @test haskey(ex.model.surfaces,7)
    # the standalone `Boolean` term form takes the same operands — its
    # production carries no `tEND`, so a trailing `;` is the stray-token
    # syntax error upstream reports
    ex=_execute_geo_source(_BOX_WALL * """
        BooleanFragments{ Volume{1}; Delete; }{ Surface{7}; Delete; };
        """)
    @test ex.errors==["syntax error (;)"]
    @test sort(collect(keys(ex.model.volumes)))==[1,2]
    # tagged form binds the single result — an overhanging wall leaves
    # one volume so `BooleanX(t)` has one top-dim entity to bind
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Rectangle(7) = {-0.5,-0.5,0.5, 2,2};
        BooleanFragments(9) = { Volume{1}; Delete; }{ Surface{7}; Delete; };
        """)
    @test isempty(ex.errors)
    @test sort(collect(keys(ex.model.volumes)))==[9]
    # the deleted wall's image pieces bind fresh — upstream `_multiBind`
    # applies the requested tag to the single top-dim result only
    @test length(ex.model.surfaces)==10
end

# A unit box plus a `2×2` wall overhanging its footprint — the sheet's
# rim survives a `Surface - Volume` cut as the outside frame, and its
# interior patch is the `Volume ∩ Surface` section.
const _BOX_OVERHANG="""
    SetFactory("OpenCASCADE");
    Box(1) = {0,0,0, 1,1,1};
    Rectangle(7) = {-0.5,-0.5,0.5, 2,2};
"""

@testset "surface object operands" begin
    # sheet minus volume keeps the outside frame, rebound at its own tag
    ex=_execute_geo_source(_BOX_OVERHANG * """
        BooleanDifference{ Surface{7}; Delete; }{ Volume{1}; Delete; }
        """)
    @test isempty(ex.errors)
    @test isempty(ex.model.volumes)
    @test haskey(ex.model.surfaces,7)
    # volume ∩ sheet yields the section sheet
    ex=_execute_geo_source(_BOX_OVERHANG * """
        BooleanIntersection{ Volume{1}; Delete; }{ Surface{7}; Delete; }
        """)
    @test isempty(ex.errors)
    @test isempty(ex.model.volumes)
    @test haskey(ex.model.surfaces,7)
end

@testset "volume-only operands stay on the native path" begin
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.5,0,0, 1,1,1};
        BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; }
        """)
    @test isempty(ex.errors)
    @test length(ex.model.volumes)==1
    # `Recursive Delete` is accepted and treated like `Delete` upstream
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.5,0,0, 1,1,1};
        BooleanUnion{ Volume{1}; Recursive Delete; }
                    { Volume{2}; Recursive Delete; }
        """)
    @test isempty(ex.errors)
    @test length(ex.model.volumes)==1
    # `Delete FExpr` gates: `Delete 0` keeps the object, `Delete 1`
    # removes the tool — the union rebinds at the post-removal maximum
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.5,0,0, 1,1,1};
        BooleanUnion{ Volume{1}; Delete 0; }{ Volume{2}; Delete 1; }
        """)
    @test isempty(ex.errors)
    @test sort(collect(keys(ex.model.volumes)))==[1,2]
    # volume 1 is the kept operand (original box faces), volume 2 the
    # union piece on fresh surfaces
    @test all(t->abs(t)<=6,ex.model.surface_loops[only(ex.model.volumes[1])])
    @test all(t->abs(t)>12,ex.model.surface_loops[only(ex.model.volumes[2])])
    # disjoint difference honours the same flags
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {5,0,0, 1,1,1};
        BooleanDifference{ Volume{1}; Delete 0; }{ Volume{2}; Delete 1; }
        """)
    @test isempty(ex.errors)
    @test sort(collect(keys(ex.model.volumes)))==[1]
end

@testset "operand member forms" begin
    # `Line` is the upstream `tCurve` lexer alias — it resolves to a
    # dim-1 operand, which the bounded kernel then rejects atomically
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Point(13) = {0,0,0}; Point(14) = {1,0,0};
        Curve(13) = {13,14};
        BooleanUnion{ Line{13}; }{ Line{13}; }
        Point(30) = {2,0,0};
        """)
    @test any(err->occursin("dimension 1 are not yet supported",err),
              ex.errors)
    @test "Could not apply boolean operator" in ex.errors
    @test haskey(ex.model.points,30)
    # `GeoEntity{d}{tags}` clauses — dim-2 GeoEntity tool
    ex=_execute_geo_source(_BOX_WALL * """
        BooleanFragments{ Volume{1}; Delete; }{ GeoEntity{2}{7}; Delete; }
        """)
    @test isempty(ex.errors)
    @test sort(collect(keys(ex.model.volumes)))==[1,2]
    # `:` wildcard and `Physical` selectors resolve to operand entities —
    # the physical group member unions with the explicit tool volume
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.5,0,0, 1,1,1};
        Physical Volume(5) = {1};
        BooleanUnion{ Physical Volume{5}; Delete; }{ Volume{2}; Delete; }
        """)
    @test isempty(ex.errors)
    @test length(ex.model.volumes)==1
    @test all(t->abs(t)>12,
              ex.model.surface_loops[only(only(values(ex.model.volumes)))])
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.5,0,0, 1,1,1};
        BooleanUnion{ Volume{:}; Delete; }{ }
        """)
    @test isempty(ex.errors)
    @test length(ex.model.volumes)==1
    # an inline `Kind(tag) = ...` definition is a `Shape` member — its
    # side effect lands and the entity joins the operand list
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanUnion{ Volume{1}; Point(9) = {9,9,9}; }{ }
        """)
    @test haskey(ex.model.points,9)
    @test any(err->occursin("dimension 0 are not yet supported",err),
              ex.errors)
    @test "Could not apply boolean operator" in ex.errors
end

@testset "operand-group grammar tail" begin
    # a `Delete` clause anywhere but the tail cannot be the BooleanOption
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanUnion{ Volume{1}; Delete; Volume{1}; }{ Volume{1}; }
        Point(30) = {3,0,0};
        """)
    @test any(err->occursin("must close its",err),ex.errors)
    @test haskey(ex.model.points,30)
    # `Delete{..}` is the entity-delete action — not a BooleanOption
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanUnion{ Volume{1}; }{ Volume{1}; Delete{1}; }
        """)
    @test any(err->occursin("malformed operand clause",err),ex.errors)
    # `;;`/leading `;` and a missing member `tEND` are syntax errors
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanUnion{ Volume{1};; }{ Volume{1}; }
        """)
    @test "syntax error (;)" in ex.errors
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanUnion{ Volume{1} }{ Volume{1}; }
        """)
    @test "syntax error (})" in ex.errors
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanUnion{ Volume{1}; Delete }{ Volume{1}; }
        """)
    @test "syntax error (})" in ex.errors
    # a `Name{..}` transform/action is not a `ListOfShapes` member
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanUnion{ Duplicata{ Volume{1}; }; }{ Volume{1}; }
        """)
    @test any(err->occursin("not a Boolean operand member",err),ex.errors)
end

@testset "operand application diagnostics" begin
    # tags that resolve to no live entity report the OCC lookup pair
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanFragments{ Volume{1}; }{ Surface{9}; }
        Point(30) = {3,0,0};
        """)
    @test "Unknown OpenCASCADE entity of dimension 2 with tag 9" in
          ex.errors
    @test "Could not apply boolean operator" in ex.errors
    @test haskey(ex.model.points,30)
    # an empty object group is upstream's early `return true` no-op
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanFragments{ }{ Volume{1}; }
        """)
    @test isempty(ex.errors)
    @test sort(collect(keys(ex.model.volumes)))==[1]
    # a lone operand maps to itself (A ∪ ∅, A \ ∅, fragments of one arg)
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanUnion{ Volume{1}; }{ }
        """)
    @test isempty(ex.errors)
    @test sort(collect(keys(ex.model.volumes)))==[1]
end

@testset "built-in factory gating" begin
    # the bare `Boolean` term reports the OCC-kernel diagnostic…
    ex=_execute_geo_source("""
        Point(1)={0,0,0};
        BooleanUnion{Volume{1};}{Volume{1};};
        Point(2)={1,0,0};
        """)
    @test "Boolean operators only available with OpenCASCADE " *
          "geometry kernel" in ex.errors
    @test haskey(ex.model.points,2)
    # …while the tagged `BooleanShape` statement silently no-ops upstream
    ex=_execute_geo_source("""
        Point(1)={0,0,0};
        BooleanUnion(3) = {Volume{1};}{Volume{1};};
        Point(2)={1,0,0};
        """)
    @test isempty(ex.errors)
    @test isempty(ex.model.volumes)
    # operand `ListOfShapes` side effects still land before the gate —
    # the inline `Shape` definition executes either way upstream
    ex=_execute_geo_source("""
        BooleanUnion{ Point(7)={1,0,0}; }{ };
        """)
    @test haskey(ex.model.points,7)
    @test "Boolean operators only available with OpenCASCADE " *
          "geometry kernel" in ex.errors
end
