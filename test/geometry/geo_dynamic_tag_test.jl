using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, ntets
using Tessella.Elements: mixed_crc

const _GEO_DYNAMIC_TAG_FIXTURE=normpath(joinpath(
    @__DIR__,"..","fixtures","geo_dynamic_tags.geo"))

function _execute_dynamic_tag_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _dynamic_tag_error(source::AbstractString)
    try
        _execute_dynamic_tag_source(source)
        return nothing
    catch err
        return err
    end
end

# Gmsh emits the `Msg::Error` detail on stderr and the `yymsg` caller line as
# the thrown diagnostic — capture both channels.
function _dynamic_tag_stderr(source::AbstractString)
    err=nothing
    text=mktemp() do path,io
        redirect_stdio(;stderr=io) do
            try
                _execute_dynamic_tag_source(source)
            catch e
                err=e
            end
        end
        flush(io)
        read(path,String)
    end
    return (err,text)
end

@testset "bounded .geo geometry and Physical tag allocators" begin
    parsed=execute_geo(_GEO_DYNAMIC_TAG_FIXTURE)
    model=parsed.model
    @test sort!(collect(keys(model.points)))==collect(1:10)
    @test sort!(collect(keys(model.curves)))==collect(1:12)
    @test sort!(collect(keys(model.loops)))==collect(13:18)
    @test sort!(collect(keys(model.surfaces)))==collect(19:24)
    @test sort!(collect(keys(model.surface_loops)))==[25]
    @test sort!(collect(keys(model.volumes)))==[26]
    @test model.loops[13]==collect(1:4)
    @test model.loops[15]==[1,10,-5,-9]
    @test model.surface_loops[25]==collect(19:24)
    @test model.physical==Dict(
        (0,61)=>collect(1:8),(1,62)=>collect(1:12),
        (2,63)=>collect(19:24),(3,64)=>[26],(0,65)=>[9,10])
    @test model.physical_names==Dict(
        (0,61)=>"corners",(1,62)=>"edges",(2,63)=>"boundary",
        (3,64)=>"domain",(0,65)=>"face probes")
    @test parsed.params.fields[1].options["PointsList"]=="{9, 10}"
    @test get(model.embeds,(2,24),NTuple{2,Int}[])==[(0,9)]
    @test get(model.embeds,(2,22),NTuple{2,Int}[])==[(0,10)]
    @test [(constraint.dim,Int(constraint.slave_entity),
            Int(constraint.master_entity))
           for constraint in model_periodic_constraints(model)]==
          [(2,22,24),(2,23,21)]

    meshed=execute_geo(_GEO_DYNAMIC_TAG_FIXTURE;mesh_dim=3)
    @test validate(meshed.mesh).ok
    @test nnodes(meshed.mesh)==125
    @test ntets(meshed.mesh)==384
    @test mesh_crc(meshed.mesh).sha==
          "a58374071a4c485a339e1c5b48b8b0f3e69bf362ff0e41f57ca1a665139e81df"
    @test length(model_periodic_nodes(meshed.model,meshed.mesh,2,22).slave_nodes)==25
    @test length(model_periodic_nodes(meshed.model,meshed.mesh,2,23).slave_nodes)==25
    projected=model_to_mixed(meshed.model,meshed.mesh,3,26)
    @test validate(projected).ok
    @test projected.physical_names==model.physical_names
    @test mixed_crc(projected).sha==
          "8a7d8009ce298b69e9f15cad9927d24446854ac5bb527ed81f5ec0125ee2a713"

    primitive_source=raw"""
        SetFactory("OpenCASCADE");
        first = newv;
        firstAgain = newv;
        Box(first) = {0,0,0,1,1,1};
        second = newv;
        Cylinder(second) = {2,0,0,0,0,1,0.5};
        third = newv;
        Sphere(third) = {4,0,0,0.5};
        fourth = newv;
        Cone(fourth) = {6,0,0,0,0,1,0.5,0.25};
        point1 = newp;
        Point(point1) = {8,0,0,1};
        point2 = newp;
        Point(point2) = {9,0,0,1};
        curve = newl;
        Line(curve) = {point1,point2};
        Mesh.MeshSizeMin = newreg / 100;
        Mesh.MeshSizeMax = newp / 10;
        """
    primitives=_execute_dynamic_tag_source(primitive_source)
    @test sort!(collect(keys(primitives.model.volumes)))==[1,13,16,19]
    # Box materializes its eight corner Points and twelve edge Curves, and
    # the curved primitives materialize their OCC layouts the same way —
    # two rim/pole Points and three Curves each, exactly the child counts
    # the allocator reserves.
    @test sort!(collect(keys(primitives.model.points)))==[1:16;]
    @test sort!(collect(keys(primitives.model.curves)))==[1:22;]
    @test primitives.params.mesh_size_min==0.23
    @test primitives.params.mesh_size_max==1.7

    physical_source=raw"""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Physical Point("witness", 65) = {1};
        Line(newreg) = {1,2};
        """
    physical=_execute_dynamic_tag_source(physical_source)
    @test sort!(collect(keys(physical.model.curves)))==[66]

    automatic_physical_source=raw"""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Line(3) = {1,2};
        Physical Point("auto point") = {1};
        Physical Curve("auto curve") = {3};
        Line(newreg) = {1,2};
        Physical Point("explicit", 20) = {2};
        Physical Curve("later") = {3};
        Physical Point("endpoints") = CombinedBoundary{Line{3};};
        Line(newreg) = {1,2};
        """
    automatic_physical=_execute_dynamic_tag_source(automatic_physical_source)
    @test sort!(collect(keys(automatic_physical.model.curves)))==[3,4,23]
    @test automatic_physical.model.physical==Dict(
        (0,1)=>[1],(1,2)=>[3],(0,20)=>[2],(1,21)=>[3],(0,22)=>[1,2])
    @test automatic_physical.model.physical_names==Dict(
        (0,1)=>"auto point",(1,2)=>"auto curve",(0,20)=>"explicit",
        (1,21)=>"later",(0,22)=>"endpoints")
    @test automatic_physical.params.physical_groups==
          automatic_physical.model.physical_names
    @test automatic_physical.model.physical_tag_max==22

    field_source=raw"""
        Field[newf] = Distance;
        Field[newf] = Min;
        Point(newf) = {0,0,0,1};
        """
    fields=_execute_dynamic_tag_source(field_source)
    @test sort!(collect(keys(fields.params.fields)))==[1,2]
    @test sort!(collect(keys(fields.model.points)))==[3]

    cone_tip_source=raw"""
        SetFactory("OpenCASCADE");
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {0,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Line(3) = {3,1};
        Curve Loop(4) = {1,2,3};
        Plane Surface(100) = {4};
        Cone(newv) = {2,0,0,0,0,1,0.5,0};
        Mesh.MeshSizeMin = newreg / 100;
        """
    cone_tip=_execute_dynamic_tag_source(cone_tip_source)
    @test sort!(collect(keys(cone_tip.model.volumes)))==[101]
    @test cone_tip.params.mesh_size_min==1.03

    # A partial torus consumes two rim Points, three Curves, and three
    # Surfaces of hidden topology; the full torus consumes one, two, and one.
    torus_partial_source=raw"""
        SetFactory("OpenCASCADE");
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {0,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Line(3) = {3,1};
        Curve Loop(4) = {1,2,3};
        Plane Surface(100) = {4};
        Torus(newv) = {2,0,0,3,1,1.5707963267948966};
        Mesh.MeshSizeMin = newreg / 100;
        """
    torus_partial=_execute_dynamic_tag_source(torus_partial_source)
    @test sort!(collect(keys(torus_partial.model.volumes)))==[101]
    @test torus_partial.params.mesh_size_min==1.04

    torus_full_source=raw"""
        SetFactory("OpenCASCADE");
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {0,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Line(3) = {3,1};
        Curve Loop(4) = {1,2,3};
        Plane Surface(100) = {4};
        Torus(newv) = {2,0,0,3,1};
        Mesh.MeshSizeMin = newreg / 100;
        """
    torus_full=_execute_dynamic_tag_source(torus_full_source)
    @test sort!(collect(keys(torus_full.model.volumes)))==[101]
    @test torus_full.params.mesh_size_min==1.02

    boolean_source=raw"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,1,1,1};
        Box(2) = {2,0,0,1,1,1};
        BooleanDifference(newv) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
        Point(newp) = {5,5,5,1};
        """
    boolean_result=_execute_dynamic_tag_source(boolean_source)
    # Gmsh 4.15.2 parity: result volume 25; the disjoint tool leaves operand 1
    # unmodified so its corner Points 1-8 persist on the preserved result and
    # newp = 9.
    @test sort!(collect(keys(boolean_result.model.volumes)))==[25]
    @test sort!(collect(keys(boolean_result.model.points)))==[1:9;]

    boolean_keep_source=raw"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,1,1,1};
        Box(2) = {2,0,0,1,1,1};
        BooleanDifference(newv) = { Volume{1}; }{ Volume{2}; Delete; };
        Point(newp) = {5,5,5,1};
        """
    boolean_keep=_execute_dynamic_tag_source(boolean_keep_source)
    # Gmsh 4.15.2: the result IsSame the kept operand, so the model keeps
    # volume 1 only (outDimTags still reports the requested tag).
    @test sort!(collect(keys(boolean_keep.model.volumes)))==[1]
    @test sort!(collect(keys(boolean_keep.model.points)))==[1:9;]

    setmax_volume_source=raw"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,1,1,1};
        Box(2) = {2,0,0,1,1,1};
        SetMaxTag Volume(40);
        BooleanDifference(newv) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
        Point(newp) = {5,5,5,1};
        """
    setmax_volume=_execute_dynamic_tag_source(setmax_volume_source)
    # SetMaxTag floor survives the Boolean operand deletes; newp is unaffected
    # by the volume floor and reflects the live point maximum + 1.
    @test sort!(collect(keys(setmax_volume.model.volumes)))==[41]
    @test sort!(collect(keys(setmax_volume.model.points)))==[1:9;]

    # `Physical X(n) op= {..}` compound assignment (`modifyPhysicalGroup`
    # ops 1/2/3). `+=` appends unconditionally (`List_Add` — no dedup);
    # `-=` removes members and deletes a group left empty; `-=` on a
    # missing group is a silent no-op; `+=`/`*=`/`/=` on a missing group
    # and `*=`/`/=` generally are recorded errors.
    compound=_execute_dynamic_tag_source(raw"""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Physical Point(6) = {1};
        Physical Point(6) += {2};
        Physical Point(6) -= {1};
        Physical Point(9) -= {1};
        Physical Point("named") = {1,2};
        Physical Point("named") += {2};
        """)
    # `named` auto-assigns tag 7 (the explicit 6 already raised the
    # counter); `+= {2}` resolves the name, appends unconditionally, and —
    # Gmsh parity, verified against `newreg` — still bumps the physical
    # counter to 8 because `setMaxPhysicalTag` runs before the
    # `setPhysicalName` name resolution. The raw member list keeps the
    # duplicate `[1,2,2]`, but the observable view dedups entity memberships
    # (`std::find` in `GEO_Internals::synchronize`).
    @test compound.model.physical==Dict((0,6)=>[2],(0,7)=>[1,2])
    @test compound.model.physical_names==Dict((0,7)=>"named")
    @test compound.model.physical_tag_max==8

    # `-= {}` (and removing the last member) deletes the group.
    emptied=_execute_dynamic_tag_source(raw"""
        Point(1) = {0,0,0,1};
        Physical Point(4) = {1};
        Physical Point(4) -= {1};
        Physical Point(5) = {1};
        Physical Point(5) -= {};
        """)
    @test isempty(emptied.model.physical)
    @test isempty(emptied.model.physical_names)

    # Gmsh `.geo` definition-site tag semantics: `tag < 0` auto-assigns the
    # active kernel's per-dimension counter, `tag == 0` is a literal tag,
    # `tag > 0` is explicit. Curve/Surface Loop counters are independent
    # namespaces; point references look up `abs(tag)` (`CompareVertex`),
    # physical member lists look up `abs` too, and explicit physical tags
    # — including negatives — are literal.
    signed_tags=_execute_dynamic_tag_source(raw"""
        Point(0) = {0,0,0,1};
        Point(-3) = {1,0,0,1};
        Point(2) = {0,1,0,1};
        Line(-2) = {0,2};
        Line(3) = {2,1};
        Line(4) = {-1,0};
        Curve Loop(0) = {1,3,4};
        Curve Loop(-5) = {1,3,4};
        Surface Loop(-2) = {7};
        Physical Point(-4) = {-0,1};
        Physical Curve(7) = {-1,3};
        """)
    @test sort!(collect(keys(signed_tags.model.points)))==[0,1,2]
    @test sort!(collect(keys(signed_tags.model.curves)))==[1,3,4]
    @test sort!(collect(keys(signed_tags.model.loops)))==[0,1]
    @test sort!(collect(keys(signed_tags.model.surface_loops)))==[1]
    # `Physical Point(-4) = {-0,1}` — `orientedPhysicals` distributes the raw
    # group tag per member: member `-0` resolves entity `abs(0)=0` (which
    # exists here) and receives `gmsh_sign(0)*(-4)=0` → view group `0`;
    # member `1` receives `gmsh_sign(1)*(-4)=-4` → view group `abs(-4)=4`.
    # `Physical Curve(7) = {-1,3}` — member `-1` receives `-7`, member `3`
    # receives `7` — both land in view group 7.
    @test signed_tags.model.physical==Dict((0,0)=>[0],(0,4)=>[1],(1,7)=>[1,3])

    # Unknown members are skipped with a warning (Gmsh reports them when
    # the group is synchronized). The raw group is still created, but the
    # observable view derives from entity memberships — a group with no
    # resolvable member leaves no view entry at all.
    unknown_member=_execute_dynamic_tag_source(raw"""
        Point(1) = {0,0,0,1};
        Physical Point(5) = {99};
        """)
    @test isempty(unknown_member.model.physical)
    @test any(occursin("Skipping unknown point 99",w)
              for w in unknown_member.warnings)

    # A duplicate `=` declaration is a recoverable "already exists" error
    # in Gmsh — execution continues; `Msg::Error` detail lands on stderr
    # and the `yymsg` caller line becomes the thrown diagnostic.
    (dup,dup_stderr)=_dynamic_tag_stderr(raw"""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Physical Point(6) = {1};
        Physical Point(6) = {2};
        """)
    @test dup isa ArgumentError
    @test occursin("Could not modify physical point",sprint(showerror,dup))
    @test occursin("Physical point 6 already exists",dup_stderr)

    (compound_missing,missing_stderr)=_dynamic_tag_stderr(raw"""
        Point(1) = {0,0,0,1};
        Physical Point(6) += {1};
        """)
    @test compound_missing isa ArgumentError
    @test occursin("Could not modify physical point",
                   sprint(showerror,compound_missing))
    @test occursin("Physical point 6 does not exist",missing_stderr)

    for op in ("*=","/=")
        (unsupported,unsupported_stderr)=_dynamic_tag_stderr(
            "Point(1)={0,0,0,1}; Physical Point(6)={1}; " *
            "Physical Point(6) $op {1};")
        @test unsupported isa ArgumentError
        @test occursin("Could not modify physical point",
                       sprint(showerror,unsupported))
        @test occursin("Unsupported operation on physical point 6",
                       unsupported_stderr)
    end

    @testset ".geo Field tags are signed literals with a live newf counter" begin
        # `Field[i] = Kind` takes a raw `(int)` — zero and negatives are
        # literal ids — and `newf` is `FieldManager::maxId() + 1`, a live
        # maximum that shrinks when `Delete Field` removes the maximum and can
        # be nonpositive when every id is.
        fields_only=_execute_dynamic_tag_source("""
            Field[0] = Box;
            Field[-2] = MathEval;
            Field[newf] = Min;
            Field[5] = Attractor;
            Delete Field[5];
            Field[newf] = Threshold;
            """)
        @test sort!(collect(keys(fields_only.params.fields)))==[-2,0,1,2]

        negative=_execute_dynamic_tag_source("""
            Field[-2] = MathEval;
            Field[newf] = Box;
            """)
        # With every id negative `maxId()` is the maximum key — `newf` is -1.
        @test sort!(collect(keys(negative.params.fields)))==[-2,-1]

        refilled=_execute_dynamic_tag_source("""
            Field[3] = Box;
            Delete Field[3];
            Field[newf] = Min;
            """)
        @test sort!(collect(keys(refilled.params.fields)))==[1]
        @test refilled.params.fields[1].kind=="Min"

        # Recoverable diagnostics match `FieldManager`/`Gmsh.y` — the run ends
        # in an error only through the accumulated `yymsg(0)` parse errors;
        # `Msg::Error`-only statements just mark the run.
        duplicate=_dynamic_tag_stderr("Field[1]=Box; Field[1]=Min;")
        @test duplicate[1] isa ArgumentError
        @test occursin("Cannot create field 1 of type 'Min'",
                       sprint(showerror,duplicate[1]))
        @test occursin("Field id 1 is already defined",duplicate[2])

        function _params(source)
            path,io=mktemp()
            write(io,source);close(io)
            return Tessella.IO.read_geo_params(path)
        end
        @test _params("Field[1]=Box; Field[1]=Min;").fields[1].kind=="Box"

        unknown_kind=_dynamic_tag_stderr("Field[3]=Bogus;")
        @test unknown_kind[1] isa ArgumentError
        @test occursin("Cannot create field 3 of type 'Bogus'",
                       sprint(showerror,unknown_kind[1]))
        @test occursin("Unknown field type \"Bogus\"",unknown_kind[2])

        missing_write=_dynamic_tag_stderr("Field[9].VIn=0.1;")
        @test missing_write[1] isa ArgumentError
        @test occursin("No field with id 9",sprint(showerror,missing_write[1]))

        # An option write ahead of the declaration is dropped — the later
        # `Field[9]=Box` starts with no options, like Gmsh.
        ordering=_params("Field[9].VIn=0.1; Field[9]=Box;")
        @test ordering.scan_errors==["No field with id 9"]
        @test isempty(ordering.fields[9].options)

        missing_delete=_dynamic_tag_stderr("Field[1]=Box; Delete Field[9];")
        @test missing_delete[1]===nothing
        @test occursin("Cannot delete field id 9, it does not exist",
                       missing_delete[2])

        background_multi=_dynamic_tag_stderr(
            "Field[1]=Box; Background Field = {1,1};")
        @test background_multi[1] isa ArgumentError
        @test occursin("Only 1 field can be set as a background field.",
                       sprint(showerror,background_multi[1]))

        unknown_command=_dynamic_tag_stderr("Foo Field = {1};")
        @test unknown_command[1] isa ArgumentError
        @test occursin("Unknown command 'Foo Field'",
                       sprint(showerror,unknown_command[1]))

        # `BoundaryLayer Field` stores raw ids — an undeclared id is not a
        # parse error (it fails later at field build, like `get(id)`).
        boundary=_execute_dynamic_tag_source("BoundaryLayer Field = {7};")
        @test boundary.params.boundary_layer_fields==[7]
    end

    invalid_sources=(
        "newp = 2;"=>"read-only",
        "newreg[] = {2};"=>"read-only",
        "Point(newp[0]) = {0,0,0,1};"=>"scalar and cannot use []",
        "Point(1)={0,0,0,1}; Point(2)={1,0,0,1}; " *
        "Physical Point(\"same\")={1}; Physical Point(\"same\")={2};"=>
            "Could not modify physical point",
        "Point(1)={0,0,0,1}; Physical Point(\"\")={1};"=>
            "automatic Physical Point group requires a nonempty name",
        "Point(1)={0,0,0,1}; Physical Point={1};"=>
            "malformed Physical declaration",
        "Point(1)={0,0,0,1}; Physical Point()={1};"=>
            "malformed Physical declaration",
        "Point(1)={0,0,0,1}; Physical Point(\"bad\",)={1};"=>
            "malformed Physical declaration",
        "BooleanFragments{Volume{1};}{Volume{1};}; Mesh.MeshSizeMax = newv;"=>
            "topology-changing statement",
    )
    for (source,message) in invalid_sources
        err=_dynamic_tag_error(source)
        @test err isa ArgumentError
        @test occursin(message,sprint(showerror,err))
    end

    # `newp`/`newf` are C++ `int` increments — past INT32_MAX they wrap to
    # INT32_MIN — while `newreg` (`NEWREG()`) is the max over the *non-point*
    # dimensions plus physical, so a wrapped dimension loses to the other
    # counters and yields 1 (verified against Gmsh 4.15.2).
    wrapped_alloc=_execute_dynamic_tag_source(raw"""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(2147483647)={2,0,0,1};
        Point(newp)={3,0,0,1};
        Line(2147483647)={1,2};
        Line(newreg)={1,2};
        Field[2147483647]=Box; Field[newf]=Box;
        """)
    @test sort!(collect(keys(wrapped_alloc.model.points)))==
        [-2147483648,1,2,2147483647]
    @test sort!(collect(keys(wrapped_alloc.model.curves)))==[1,2147483647]
    @test haskey(wrapped_alloc.params.fields,2147483647)
    @test haskey(wrapped_alloc.params.fields,-2147483648)

    # `setMaxPhysicalTag(t + 1)` runs on a raw `int` — at `typemax(Int32)` the
    # name-only counter bump wraps to `typemin(Int32)` rather than erroring,
    # and the name binds to the wrapped tag (verified against Gmsh 4.15.2).
    wrapped=_execute_dynamic_tag_source(raw"""
        Point(1) = {0,0,0,1};
        Physical Point("last",2147483647) = {1};
        Physical Point("overflow") = {1};
        """)
    @test wrapped.model.physical_names==Dict(
        (0,2147483647)=>"last",(0,-2147483648)=>"overflow")
    @test wrapped.model.entity_physicals==Dict(
        (0,1)=>[2147483647,-2147483648])

    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end
