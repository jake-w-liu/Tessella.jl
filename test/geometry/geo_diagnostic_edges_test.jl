# `.geo` executor diagnostic/operand-resolution edge cases:
#
# * A degenerate `BooleanFragments` operand (a volume fragmented by its own
#   face) makes the bounded OCC engine raise `KeyError`. Upstream
#   `OCC_Internals::booleanOperator` turns every kernel failure into
#   `Msg::Error` + `false`, and the grammar reports the recoverable
#   `Could not apply boolean operator` — the exception must never escape
#   `execute_geo`, both messages must land in `errors`, and parsing must
#   continue. The engine-side root cause lives in ModelBoolean and is
#   tracked separately; this pins the executor's recoverability contract.
#
# * `Physical Kind(t) = {..}` members resolve against the synchronized
#   model. Upstream runs `OCC_Internals::synchronize` before the GEO
#   internals sync at every physical-resolution point, so a lazily bound
#   OCC kernel entity (an API-built model merged into `.geo`, before any
#   `.geo` OCC op materializes it) counts as an existing member — not a
#   warn-skipped unknown.
#
# * `Physical Kind(t) op {..}` inside a `BooleanX{..}{..}` operand group —
#   upstream's `tPhysical GeoEntity '(' PhysicalId ')' NumericAffectation
#   ListOfDouble tEND` Shape production: the declaration runs as a member
#   side effect and contributes no operand entities. It must not abort the
#   operand list with `syntax error (`.

using Test
using Tessella

function _execute_geo_source(source::AbstractString;_initial_model=nothing)
    dir=mktempdir()
    path=joinpath(dir,"script.geo")
    write(path,source)
    return execute_geo(path;_initial_model=_initial_model)
end

@testset "degenerate self-fragment reports recoverable diagnostics" begin
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        BooleanFragments(7) = { Volume{1}; Delete; }{ Surface{2}; Delete; };
        Point(99) = {9,9,9};
        """)
    # The engine's own message is recorded, followed by the recoverable
    # `Could not apply boolean operator` pair — nothing escapes the call.
    @test any(e->occursin("not found",e),ex.errors)
    @test "Could not apply boolean operator" in ex.errors
    # The statement stream continues after the failed boolean.
    @test haskey(ex.model.points,99)
end

@testset "physical declaration resolves OCC-bound members" begin
    # An OCC-built model merged into `.geo` carries lazily bound kernel
    # entities that no native record has materialized yet — upstream
    # materializes them at the sync that precedes physical resolution.
    m=Tessella.Model.GeoModel()
    Tessella.Model._model_occ_add_box!(m,0,0,0,1,1,1;tag=1)
    Tessella.Model._model_occ_add_box!(m,5,0,0,1,1,1;tag=2)
    @test haskey(m.occ_kernel.bindings,(3,1))
    @test haskey(m.occ_kernel.bindings,(3,2))
    @test isempty(m.volumes)
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Physical Volume(5) = {1,2};
        Physical Volume(9) = {77};
        """;_initial_model=m)
    @test ex.model.physical[(3,5)]==[1,2]
    # The bindings materialized through the sync — genuine entities, no
    # warn-skip; a truly unknown member still reports upstream's warning.
    @test haskey(ex.model.volumes,1) && haskey(ex.model.volumes,2)
    @test ex.warnings==["Skipping unknown volume 77 in physical volume 9"]
    @test isempty(ex.errors)
end

@testset "nested physical declaration inside a boolean operand group" begin
    # The declaration member runs as a side effect and contributes no
    # operand entities; with no `Delete` option the tool volume survives
    # so the declared group resolves at the sync — upstream emits no
    # diagnostic at all here.
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.5,0,0, 1,1,1};
        BooleanUnion(9) = { Volume{1}; }{ Physical Volume(5)={2}; Volume{2}; };
        Point(97) = {7,7,7};
        """)
    @test isempty(ex.errors)
    @test isempty(ex.warnings)
    @test ex.model.physical[(3,5)]==[2]
    @test haskey(ex.model.volumes,9)
    @test haskey(ex.model.points,97)
    # With `Delete` the operand list consumes tool volume 2 before the
    # physical view resolves — upstream reports the member skip the same
    # way (`Skipping unknown volume 2 in physical volume 5`).
    ex=_execute_geo_source("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.5,0,0, 1,1,1};
        BooleanUnion(9) = { Volume{1}; Delete; }{ Physical Volume(5)={2}; Volume{2}; Delete; };
        Point(98) = {8,8,8};
        """)
    @test isempty(ex.errors)
    @test ex.warnings==["Skipping unknown volume 2 in physical volume 5"]
    @test haskey(ex.model.volumes,9)
    @test haskey(ex.model.points,98)
end
