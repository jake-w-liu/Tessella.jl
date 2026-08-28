using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, ntets
using Tessella.Elements: mixed_crc

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
    @test nnodes(meshed.mesh)==5
    @test ntets(meshed.mesh)==4
    @test mesh_crc(meshed.mesh).sha==
          "71ab10cf31fa64d469e1bc3985bd8c50bb240d1cdefaebbc17101bce22e7008b"
    projected=model_to_mixed(meshed.model,meshed.mesh,3,601)
    @test validate(projected).ok
    @test mixed_crc(projected).sha==
          "aa3127e5c1a302ebdb98890a4d7a62d5cb385e3cf73aa024372f809a7542a45d"

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
    @test sort!(collect(keys(primitive.model.points)))==[7]
    @test sort!(collect(keys(primitive.model.volumes)))==[1]
    @test primitive.params.mesh_size_max==1.6

    implicit_overlap=_execute_set_max_tag_source(raw"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,1,1,1};
        Point(1) = {0.5,0.5,0.5,1};
        Mesh.MeshSizeMin = newp;
        """)
    @test sort!(collect(keys(implicit_overlap.model.points)))==[1]
    @test implicit_overlap.params.mesh_size_min==9.0

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

    invalid_sources=(
        "SetMaxTag Point(2147483647); Point(newp)={0,0,0,1};"=>
            "no Point tags remain",
        "SetMaxTag Point(2147483648);"=>"signed 32-bit integer range",
        "SetMaxTag Field(10);"=>"unrecognized statement",
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
