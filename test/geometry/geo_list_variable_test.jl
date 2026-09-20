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
          "7290d425e4b3e881889b8b3bb6661a077b390cce3f1870d26487c5c6fcca55c0"
    @test length(model_periodic_nodes(meshed.model,meshed.mesh,2,4).slave_nodes)==25
    @test length(model_periodic_nodes(meshed.model,meshed.mesh,2,5).slave_nodes)==25
    projected=model_to_mixed(meshed.model,meshed.mesh,3,1)
    @test validate(projected).ok
    @test projected.physical_names==model.physical_names
    @test mixed_crc(projected).sha==
          "1b3447d16ca7eea18b859b0ba8a0637a47bb859ce4ce327e113b9a1155dff7f5"

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
        "a[] = {}; Point(1)={0,0,0,1}; Line(1)=a[];"=>"entity list is empty",
        "Point(1)={0,0,0,1}; Point(2)={1,0,0,1}; Line(1)=1,2;"=>
            "unexpected token",
        "coords[] = {0,0,0}; Point(1) = {coords[], 1};"=>
            "requires exactly one scalar index",
    )
    for (source,message) in invalid_sources
        err=_list_variable_error(source)
        @test err isa ArgumentError
        @test occursin(message,sprint(showerror,err))
    end

    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end
