using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, ntets
using Tessella.Elements: mixed_crc

const _GEO_LIST_VARIABLE_FIXTURE=normpath(joinpath(
    @__DIR__,"..","fixtures","geo_list_variables.geo"))

function _execute_list_variable_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _list_variable_error(source::AbstractString)
    try
        _execute_list_variable_source(source)
        return nothing
    catch err
        return err
    end
end

@testset "bounded .geo numeric list variables" begin
    parsed=execute_geo(_GEO_LIST_VARIABLE_FIXTURE)
    model=parsed.model
    @test sort!(collect(keys(model.points)))==vcat(collect(1:8),[101,102])
    @test sort!(collect(keys(model.curves)))==collect(1:12)
    @test sort!(collect(keys(model.loops)))==collect(1:6)
    @test sort!(collect(keys(model.surfaces)))==collect(1:6)
    @test sort!(collect(keys(model.surface_loops)))==[1]
    @test sort!(collect(keys(model.volumes)))==[1]
    @test model.loops[1]==collect(1:4)
    @test model.loops[3]==[1,10,-5,-9]
    @test model.surface_loops[1]==collect(1:6)
    @test model.physical==Dict(
        (0,61)=>collect(1:8),(1,62)=>collect(1:12),
        (2,63)=>collect(1:6),(3,64)=>[1],(0,65)=>[101,102])
    @test model.physical_names==Dict(
        (0,61)=>"corners",(1,62)=>"edges",(2,63)=>"boundary",
        (3,64)=>"domain",(0,65)=>"face probes")
    @test parsed.params.fields[201].options["PointsList"]=="{101, 102}"
    @test get(model.embeds,(2,6),NTuple{2,Int}[])==[(0,101)]
    @test get(model.embeds,(2,4),NTuple{2,Int}[])==[(0,102)]
    @test [(constraint.dim,Int(constraint.slave_entity),
            Int(constraint.master_entity))
           for constraint in model_periodic_constraints(model)]==
          [(2,4,6),(2,5,3)]

    meshed=execute_geo(_GEO_LIST_VARIABLE_FIXTURE;mesh_dim=3)
    @test validate(meshed.mesh).ok
    @test nnodes(meshed.mesh)==125
    @test ntets(meshed.mesh)==384
    @test mesh_crc(meshed.mesh).sha==
          "a58374071a4c485a339e1c5b48b8b0f3e69bf362ff0e41f57ca1a665139e81df"
    @test length(model_periodic_nodes(meshed.model,meshed.mesh,2,4).slave_nodes)==25
    @test length(model_periodic_nodes(meshed.model,meshed.mesh,2,5).slave_nodes)==25
    projected=model_to_mixed(meshed.model,meshed.mesh,3,1)
    @test validate(projected).ok
    @test projected.physical_names==model.physical_names
    @test mixed_crc(projected).sha==
          "9a608febd598ab7090d3fb64cb9ff7c7b90c26dc4cc8442e7f54fa121a8fa9ad"

    # Diagnostics match Gmsh 4.15.2's recoverable `yymsg` text: the entity
    # statements still execute (e.g. `Point(missing[0])` creates `Point(0)`
    # after the unknown-variable read softens to 0) and the accumulated
    # errors surface at end of parse.
    invalid_sources=(
        "a[] = {1}; Point(missing[0]) = {0,0,0,1};"=>"Unknown variable",
        "a[] = {1}; Point(a[-1]) = {0,0,0,1};"=>"Uninitialized variable",
        "a[] = {1,2}; a[{0,1}] = {3};"=>
            "Incompatible array dimensions in affectation",
        "a[] = {1,2}; a[] *= 2;"=>"Operators *= and /= not available for lists",
        "newp[] = {1};"=>"read-only",
        "a[] = {1:65537};"=>"expanded list exceeds 65536 entries",
        "a[] = {}; Point(1)={0,0,0,1}; Line(1)=a[];"=>"Could not add line",
        "Point(1)={0,0,0,1}; Point(2)={1,0,0,1}; Line(1)=1,2;"=>
            "syntax error (,)",
        "coords[] = {0,0,0}; Point(1) = {coords[], 1};"=>
            "syntax error (])",
    )
    for (source,message) in invalid_sources
        err=_list_variable_error(source)
        @test err isa ArgumentError
        @test occursin(message,sprint(showerror,err))
    end

    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end

@testset "Gmsh 4.15.2 ++/-- postfix increments" begin
    # Statement-level scalar postfix mutates in place (Affectation
    # `String__Index NumericIncrement`); expression-level postfix returns the
    # OLD value then increments (FExpr `String__Index NumericIncrement`).
    parsed=_execute_list_variable_source("""
        x = 5;
        y = x++;
        z = x--;
        w = 0;
        w++;
        w++;
    """)
    @test parsed.values["x"]==5.0
    @test parsed.values["y"]==5.0
    @test parsed.values["z"]==6.0
    @test parsed.values["w"]==2.0

    # Indexed statement form auto-resizes with zeros (`incrementVariable`);
    # the expression indexed form returns the old element and mutates it,
    # with `(int)` truncation on the index.
    parsed=_execute_list_variable_source("""
        a[] = {1,2};
        a[5]++;
        b() = {9,9};
        c = b(1)++;
        u() = {3,4};
        u(1.9)++;
    """)
    @test parsed.lists["a"]==[1.0,2.0,0.0,0.0,0.0,1.0]
    @test parsed.values["c"]==9.0
    @test parsed.lists["b"]==[9.0,10.0]
    @test parsed.lists["u"]==[3.0,5.0]

    # The expression indexed form works on scalars too (`s.value[0]` reads
    # element zero of the scalar payload) while the statement form requires a
    # list — Gmsh's `s.list` check inside `incrementVariable`.
    parsed=_execute_list_variable_source("""
        s = 7;
        t = s[0]++;
    """)
    @test parsed.values["t"]==7.0
    @test parsed.values["s"]==8.0

    err=_list_variable_error("d = 1; d[2]++;")
    @test err isa ArgumentError
    @test occursin("Variable 'd' is not a list",sprint(showerror,err))

    err=_list_variable_error("v[] = {1}; v++;")
    @test err isa ArgumentError
    @test occursin("Variable 'v' is a list",sprint(showerror,err))

    err=_list_variable_error("missing[0]++;")
    @test err isa ArgumentError
    @test occursin("Unknown variable 'missing'",sprint(showerror,err))

    err=_list_variable_error("unknownvar++;")
    @test err isa ArgumentError
    @test occursin("Unknown variable 'unknownvar'",sprint(showerror,err))

    # Prefix `++`/`--` is a syntax error in Gmsh's grammar too (`++` is the
    # offending token); a negative index reaches `s.value[-1]` UB upstream,
    # so Tessella fails explicitly.
    err=_list_variable_error("x = 1; y = ++x;")
    @test err isa ArgumentError
    @test occursin("syntax error (++)",sprint(showerror,err))
    err=_list_variable_error("a[] = {1}; a[-1]++;")
    @test err isa ArgumentError
    @test occursin("negative index",sprint(showerror,err))
end
