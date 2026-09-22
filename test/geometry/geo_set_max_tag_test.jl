using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, ntets
using Tessella.Elements: mixed_crc
using Tessella.IO: read_geo_params

const _GEO_SET_MAX_TAG_FIXTURE=normpath(joinpath(
    @__DIR__,"..","fixtures","geo_set_max_tags.geo"))

function _execute_set_max_tag_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _set_max_tag_error(source::AbstractString)
    try
        _execute_set_max_tag_source(source)
        return nothing
    catch err
        return err
    end
end

@testset "bounded .geo SetMaxTag counters" begin
    execution=execute_geo(_GEO_SET_MAX_TAG_FIXTURE)
    model=execution.model
    @test sort!(collect(keys(model.points)))==collect(101:104)
    @test sort!(collect(keys(model.curves)))==collect(201:206)
    @test sort!(collect(keys(model.loops)))==collect(207:210)
    @test sort!(collect(keys(model.surfaces)))==collect(401:404)
    @test sort!(collect(keys(model.surface_loops)))==[405]
    @test sort!(collect(keys(model.volumes)))==[601]
    @test model.physical_names==Dict(
        (0,603)=>"corners",(1,604)=>"edges",
        (2,605)=>"boundary",(3,606)=>"domain")
    @test execution.params.fields[1].options["PointsList"]==
          "{101, 102, 103, 104}"

    meshed=execute_geo(_GEO_SET_MAX_TAG_FIXTURE;mesh_dim=3)
    @test validate(meshed.mesh).ok
    @test nnodes(meshed.mesh)==30
    @test ntets(meshed.mesh)==60
    @test mesh_crc(meshed.mesh).sha==
          "ebb68b90154d9926a82e90d714a7eb4aecc3e70b5c51cb5cb88711a5ffce2851"
    projected=model_to_mixed(meshed.model,meshed.mesh,3,601)
    @test validate(projected).ok
    @test mixed_crc(projected).sha==
          "99234299a71ee9bd717638b770377e900973a6f0527ad554ecf8896178e04be6"

    lowered=_execute_set_max_tag_source(raw"""
        Point(10) = {0,0,0,1};
        SetMaxTag Point(20);
        SetMaxTag Point(5);
        Point(newp) = {1,0,0,1};
        Mesh.MeshSizeMin = newp / 100;
        """)
    @test sort!(collect(keys(lowered.model.points)))==[6,10]
    @test lowered.params.mesh_size_min==0.07

    fractional=_execute_set_max_tag_source(raw"""
        SetMaxTag Point(10.9);
        Point(newp) = {0,0,0,1};
        """)
    @test sort!(collect(keys(fractional.model.points)))==[11]

    allocator_expression=_execute_set_max_tag_source(raw"""
        SetMaxTag Point(10);
        SetMaxTag Point(newp + 4);
        Point(newp) = {0,0,0,1};
        """)
    @test sort!(collect(keys(allocator_expression.model.points)))==[16]

    namespaces=_execute_set_max_tag_source(raw"""
        SetMaxTag Point(10);
        SetMaxTag Curve(20);
        SetMaxTag Surface(30);
        SetMaxTag Volume(40);
        Mesh.MeshSizeMin = newp / 100;
        Mesh.MeshSizeMax = newreg / 100;
        """)
    @test namespaces.params.mesh_size_min==0.11
    @test namespaces.params.mesh_size_max==0.41

    auxiliary=_execute_set_max_tag_source(raw"""
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {0,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Line(3) = {3,1};
        Curve Loop(100) = {1,2,3};
        SetMaxTag Curve(5);
        SetMaxTag Surface(5);
        SetMaxTag Volume(5);
        Mesh.MeshSizeMin = newreg / 100;
        """)
    @test auxiliary.params.mesh_size_min==1.01

    primitive=_execute_set_max_tag_source(raw"""
        SetFactory("OpenCASCADE");
        Point(7) = {9,9,9,1};
        SetMaxTag Point(5);
        Box(1) = {0,0,0,1,1,1};
        Mesh.MeshSizeMax = newp / 10;
        """)
    # The Box materializes its eight corner Points on the lowest free tags
    # (past explicit Point 7), exactly like Gmsh's OCC `addBox`.
    @test sort!(collect(keys(primitive.model.points)))==[7:15;]
    @test sort!(collect(keys(primitive.model.volumes)))==[1]
    @test primitive.params.mesh_size_max==1.6

    # Gmsh 4.15.2 errors "OpenCASCADE point with tag 1 already exists"
    # (Msg::Error) + "Could not add point" (yymsg) when an explicit Point
    # collides with a materialized Box corner; the accumulated `Could not add`
    # diagnostic is what reaches the `execute_geo` exception.
    overlap_error=_set_max_tag_error(raw"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,1,1,1};
        Point(1) = {0.5,0.5,0.5,1};
        Mesh.MeshSizeMin = newp;
        """)
    @test overlap_error isa ArgumentError
    @test occursin("Could not add point",sprint(showerror,overlap_error))

    for kind in ("Point","Curve","Surface","Volume")
        allocator=kind=="Point" ? "newp" : "newreg"
        raised_only=_execute_set_max_tag_source(
            "SetFactory(\"OpenCASCADE\"); " *
            "SetMaxTag $kind(20); SetMaxTag $kind(5); " *
            "Mesh.MeshSizeMin = $allocator;")
        @test raised_only.params.mesh_size_min==21.0
    end
    switched=_execute_set_max_tag_source(raw"""
        SetFactory("OpenCASCADE");
        SetMaxTag Point(20);
        SetFactory("Built-in");
        SetMaxTag Point(5);
        Mesh.MeshSizeMin = newp;
        """)
    @test switched.params.mesh_size_min==21.0

    surface_primitive=_execute_set_max_tag_source(raw"""
        SetFactory("OpenCASCADE");
        Point(1) = {0,0,0,1};
        Point(2) = {1,0,0,1};
        Point(3) = {0,1,0,1};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Line(3) = {3,1};
        Curve Loop(4) = {1,2,3};
        Plane Surface(100) = {4};
        SetMaxTag Surface(5);
        Cone(newv) = {2,0,0,0,0,1,0.5,0};
        Mesh.MeshSizeMin = newreg / 100;
        """)
    @test sort!(collect(keys(surface_primitive.model.volumes)))==[101]
    @test surface_primitive.params.mesh_size_min==1.03

    negative=_execute_set_max_tag_source(
        "SetMaxTag Point(-1); Mesh.MeshSizeMin = newp;")
    @test negative.params.mesh_size_min==0.0
    negative_switched=_execute_set_max_tag_source(raw"""
        SetFactory("Built-in");
        SetMaxTag Point(-1);
        Mesh.MeshSizeMin = newp;
        SetFactory("OpenCASCADE");
        Mesh.MeshSizeMax = newp;
        """)
    @test negative_switched.params.mesh_size_min==0.0
    @test negative_switched.params.mesh_size_max==1.0

    recovered=_execute_set_max_tag_source(raw"""
        SetMaxTag Point(2147483647);
        SetMaxTag Point(0);
        Point(newp) = {0,0,0,1};
        """)
    @test sort!(collect(keys(recovered.model.points)))==[1]

    # `newp` is a C++ `int` increment — `SetMaxTag Point(2147483647)` wraps it
    # to INT32_MIN, and `Point(newp)` then auto-assigns the same wrapped tag
    # (verified against Gmsh 4.15.2: the point lands at tag -2147483648).
    wrapped=_execute_set_max_tag_source(raw"""
        SetMaxTag Point(2147483647);
        Point(newp) = {0,0,0,1};
        """)
    @test sort!(collect(keys(wrapped.model.points)))==[-2147483648]

    invalid_sources=(
        "SetMaxTag Point(2147483648);"=>"signed 32-bit integer range",
        "SetMaxTag Field(10);"=>"syntax error (Field)",
        "SetFactory(\"Unknown\");"=>
            "accepts only \"Built-in\" or \"OpenCASCADE\"",
        "SetFactory(factoryName);"=>"requires a literal",
    )
    for (source,message) in invalid_sources
        err=_set_max_tag_error(source)
        @test err isa ArgumentError
        @test occursin(message,sprint(showerror,err))
    end

    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end

function _read_geo_params_source(source::AbstractString)
    return mktemp() do path,io
        write(io,source)
        close(io)
        read_geo_params(path)
    end
end

@testset "A2 option writes, OldNewReg, lifecycle, factory gating" begin
    # NumberOption writes apply Gmsh's storage transforms: `Mesh.RandomSeed`
    # is an unsigned-int option — the float RHS truncates on write.
    seed=_read_geo_params_source(
        "Mesh.RandomSeed = 42.9; Mesh.MeshSizeMin = Mesh.RandomSeed;")
    @test seed.mesh_size_min==42.0
    @test seed.random_seed==42

    compound=_read_geo_params_source(
        "Mesh.RandomSeed = 10; Mesh.RandomSeed += 5; " *
        "Mesh.RandomSeed++; Mesh.MeshSizeMin = Mesh.RandomSeed;")
    @test compound.mesh_size_min==16.0

    # `Geometry.OldNewReg` (default 1): `news`/`newl`/`newv` read the shared
    # NEWREG allocator. SetMaxTag on the curve-loop dim raises NEWREG.
    old_new=_read_geo_params_source(raw"""
        SetMaxTag Surface(20);
        SetMaxTag GeoEntity{-1}(40);
        Mesh.MeshSizeMin = news;
        Mesh.MeshSizeMax = newll;
        """)
    @test old_new.mesh_size_min==41.0
    @test old_new.mesh_size_max==41.0
    # `GeoEntity{-1}` produces the "dim out of range" diagnostic but still
    # applies the loop counter — matching upstream's non-aborting yymsg.
    @test any(e->occursin("out of range",e),old_new.scan_errors)

    # Under `OldNewReg = 0` each symbol reads its own dimension counter.
    per_dim=_read_geo_params_source(raw"""
        Geometry.OldNewReg = 0;
        SetMaxTag Surface(20);
        SetMaxTag GeoEntity{-1}(40);
        Mesh.MeshSizeMin = news;
        Mesh.MeshSizeMax = newll;
        """)
    @test per_dim.mesh_size_min==21.0
    @test per_dim.mesh_size_max==41.0

    # Int options cast `(int)val` on write: 0.7 truncates to 0 → per-dim mode.
    truncated=_read_geo_params_source(raw"""
        Geometry.OldNewReg = 0.7;
        SetMaxTag Surface(20);
        SetMaxTag GeoEntity{-1}(40);
        Mesh.MeshSizeMin = news;
        """)
    @test truncated.mesh_size_min==21.0

    # `newreg` stays the shared allocator under `OldNewReg = 0`.
    reg=_read_geo_params_source(raw"""
        Geometry.OldNewReg = 0;
        SetMaxTag GeoEntity{-1}(40);
        Mesh.MeshSizeMin = newreg;
        """)
    @test reg.mesh_size_min==41.0

    # `GeoEntity{-2}` is the surface-loop counter; out-of-range dims error
    # without applying.
    loop2=_read_geo_params_source(
        "SetMaxTag GeoEntity{-2}(30); Mesh.MeshSizeMin = newsl;")
    @test loop2.mesh_size_min==31.0
    @test any(e->occursin("out of range",e),loop2.scan_errors)
    bad_dim=_read_geo_params_source(
        "SetMaxTag GeoEntity{9}(40); Mesh.MeshSizeMin = newll;")
    @test bad_dim.mesh_size_min==1.0
    @test any(e->occursin("out of range",e),bad_dim.scan_errors)

    # `SetFactory` synchronizes counters across dims -2..3: a built-in
    # curve-loop max carries into OCC reads.
    synced=_read_geo_params_source(raw"""
        SetMaxTag GeoEntity{-1}(40);
        SetFactory("OpenCASCADE");
        Mesh.MeshSizeMin = newll;
        """)
    @test synced.mesh_size_min==41.0

    # `Geometry.FirstEntityTag`/`FirstPhysicalTag` initialize fresh-model
    # counters (`NewModel`/`Delete All`/`Delete Model` all re-read them).
    fresh=_read_geo_params_source(
        "Geometry.FirstEntityTag = 5; NewModel; Mesh.MeshSizeMin = newp;")
    @test fresh.mesh_size_min==5.0
    clamped=_read_geo_params_source(
        "Geometry.FirstEntityTag = 0; NewModel; Mesh.MeshSizeMin = newp;")
    @test clamped.mesh_size_min==1.0
    fresh_phys=_execute_set_max_tag_source("""
        Geometry.FirstPhysicalTag = 7;
        NewModel;
        Point(1) = {0,0,0};
        Point(2) = {1,0,0};
        Physical Point("a") = {1};
        Physical Point("b") = {2};
        """)
    @test sort!(collect(keys(fresh_phys.model.physical)))==[(0,7),(0,8)]

    # `NewModel` resets every counter; `Delete Physicals` keeps the physical
    # counter AND the names (upstream `resetPhysicalGroups` semantics).
    reset=_read_geo_params_source(
        "Point(50) = {0,0,0,1}; NewModel; Mesh.MeshSizeMin = newp;")
    @test reset.mesh_size_min==1.0
    phys_kept=_execute_set_max_tag_source("""
        Point(1) = {0,0,0};
        Physical Point("a") = {1};
        Physical Point("b") = {1};
        Delete Physicals;
        Physical Point("c") = {1};
        """)
    @test haskey(phys_kept.model.physical,(0,3))

    # OCC-internals lifetime: `Delete Model` resets geometry but the kernel
    # stays alive (Box still works); `Delete All`/`NewModel` create a fresh
    # GModel — the OCC-gated primitives revert to the built-in diagnostic.
    occ_alive=_execute_set_max_tag_source("""
        SetFactory("OpenCASCADE");
        Delete Model;
        Box(1) = {0,0,0,1,1,1};
        """)
    @test sort!(collect(keys(occ_alive.model.volumes)))==[1]
    for lifecycle in ("Delete All;","NewModel;")
        err=_set_max_tag_error(
            "SetFactory(\"OpenCASCADE\"); $lifecycle " *
            "Box(1) = {0,0,0,1,1,1};")
        @test err isa ArgumentError
        @test occursin("Box only available with OpenCASCADE",sprint(showerror,err))
    end

    # Primitive gating under the built-in factory — recoverable diagnostic,
    # nothing created, counters untouched (so `newv` reads 1 afterwards).
    for (stmt,name) in (
            ("Box(1) = {0,0,0,1,1,1};","Box"),
            ("Cylinder(1) = {0,0,0,0,0,1,0.5};","Cylinder"),
            ("Sphere(1) = {0,0,0,0.5};","Sphere"),
            ("Cone(1) = {0,0,0,0,0,1,0.5,0.2};","Cone"),
            ("Torus(1) = {0,0,0,3,1};","Torus"))
        err=_set_max_tag_error(stmt)
        @test err isa ArgumentError
        @test occursin("$name only available with OpenCASCADE geometry kernel",
                       sprint(showerror,err))
        gated=_read_geo_params_source("$stmt Mesh.MeshSizeMin = newv;")
        @test gated.mesh_size_min==1.0
    end

    # The two-point `Sphere`/`PolarSphere` forms stay built-in under either
    # factory (upstream `newGeometrySphere`/`newGeometryPolarSphere`) — they
    # register parametric surfaces and create no volume entities.
    builtin_forms=_execute_set_max_tag_source("""
        Point(1) = {0,0,0};
        Point(2) = {1,0,0};
        Sphere(9) = {1,2};
        PolarSphere(10) = {1,2};
        """)
    @test isempty(builtin_forms.model.volumes)
    @test isempty(builtin_forms.model.surfaces)

    # Boolean gating asymmetry (Gmsh.y): the tagged `(t) =` form silently
    # no-ops under built-in; the standalone/list term reports the diagnostic.
    tagged_noop=_execute_set_max_tag_source("""
        Point(1) = {0,0,0};
        BooleanUnion(3) = {Volume{1}; Delete;}{Volume{2}; Delete;};
        """)
    @test isempty(tagged_noop.model.volumes)
    untagged=_set_max_tag_error(
        "BooleanUnion{Volume{1}; Delete;}{Volume{2}; Delete;};")
    @test untagged isa ArgumentError
    @test occursin("Boolean operators only available",sprint(showerror,untagged))
    term=_set_max_tag_error(
        "v() = BooleanUnion{Volume{1}; Delete;}{Volume{2}; Delete;};")
    @test term isa ArgumentError
    @test occursin("Boolean operators only available",sprint(showerror,term))
end
