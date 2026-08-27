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

@testset "bounded .geo dynamic tag allocators" begin
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
    @test nnodes(meshed.mesh)==11
    @test ntets(meshed.mesh)==16
    @test mesh_crc(meshed.mesh).sha==
          "2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808"
    @test length(model_periodic_nodes(meshed.model,meshed.mesh,2,22).slave_nodes)==5
    @test length(model_periodic_nodes(meshed.model,meshed.mesh,2,23).slave_nodes)==4
    projected=model_to_mixed(meshed.model,meshed.mesh,3,26)
    @test validate(projected).ok
    @test projected.physical_names==model.physical_names
    @test mixed_crc(projected).sha==
          "99aeefc2e269090b518f5896f2388bd419c8fb44d06e4743579d50d61aaedf81"

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
    @test sort!(collect(keys(primitives.model.points)))==[15,16]
    @test sort!(collect(keys(primitives.model.curves)))==[22]
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

    invalid_sources=(
        "newp = 2;"=>"read-only",
        "newreg[] = {2};"=>"read-only",
        "Point(newp[0]) = {0,0,0,1};"=>"scalar and cannot use []",
        "Point(2147483647)={0,0,0,1}; Point(newp)={1,0,0,1};"=>
            "no Point tags remain",
        "Point(1)={0,0,0,1}; Point(2)={1,0,0,1}; " *
        "Line(2147483647)={1,2}; Line(newreg)={1,2};"=>
            "no geometric region tags remain",
        "Field[2147483647]=Box; Field[newf]=Box;"=>
            "no Field tags remain",
        "Box(1)={0,0,0,1,1,1}; Box(13)={2,0,0,1,1,1}; " *
        "BooleanUnion(25)={Volume{1};Delete;}{Volume{13};Delete;}; " *
        "next = newv;"=>"topology-changing statement",
    )
    for (source,message) in invalid_sources
        err=_dynamic_tag_error(source)
        @test err isa ArgumentError
        @test occursin(message,sprint(showerror,err))
    end

    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end
