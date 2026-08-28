using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntris, triangle_area
using Tessella.Elements: mixed_crc

const _GEO_POINT_MESH_SIZE_FIXTURE=normpath(joinpath(
    @__DIR__,"..","fixtures","geo_point_mesh_sizes.geo"))

function _execute_point_mesh_size_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _point_mesh_size_error(source::AbstractString)
    try
        _execute_point_mesh_size_source(source)
        return nothing
    catch err
        return err
    end
end

@testset "atomic Point mesh-size constraints" begin
    model=GeoModel()
    @test add_point!(model,0,0,0;tag=1,mesh_size=0.8)==1
    @test add_point!(model,1,0,0;tag=2,mesh_size=0.6)==2
    @test set_point_mesh_size!(model,(1,2,1),0.4)===nothing
    @test model.point_size==Dict(1=>0.4,2=>0.4)

    @test set_point_mesh_size!(model,[1],1//2)===nothing
    @test model.point_size==Dict(1=>0.5,2=>0.4)
    stable=copy(model.point_size)
    for (points,size) in (
            ([1,99],0.25),([0],0.25),([-1],0.25),([true],0.25),
            ([1.0],0.25),([big(typemax(Int32))+1],0.25),
            (Int[],0.25),(Set([1]),0.25),([1],true),([1],0.0),
            ([1],-0.25),([1],Inf),([1],NaN),([1],"0.25"),
            ([1],big(10)^1000))
        @test_throws ArgumentError set_point_mesh_size!(model,points,size)
        @test model.point_size==stable
    end

    parsed=execute_geo(_GEO_POINT_MESH_SIZE_FIXTURE)
    @test parsed.model.point_size==Dict(
        1=>0.8,2=>0.4,3=>3*0.8/4,4=>0.4)
    @test parsed.model.physical==Dict((0,10)=>collect(1:4),(2,20)=>[1])
    @test parsed.model.physical_names==Dict(
        (0,10)=>"corners",(2,20)=>"domain")

    meshed=execute_geo(_GEO_POINT_MESH_SIZE_FIXTURE;mesh_dim=2)
    @test validate(meshed.mesh).ok
    @test nnodes(meshed.mesh)==23
    @test ntris(meshed.mesh)==28
    @test mesh_crc(meshed.mesh).sha==
          "f1bfa8a1cc61158cc6293540ad6ce6d7ce48054a616a3df57d93597857f1b089"
    area=sum(triangle_area(
        node(meshed.mesh,meshed.mesh.tris[1,triangle]),
        node(meshed.mesh,meshed.mesh.tris[2,triangle]),
        node(meshed.mesh,meshed.mesh.tris[3,triangle]))
        for triangle in 1:ntris(meshed.mesh))
    @test area≈1.0 atol=32eps(Float64)
    projected=model_to_mixed(meshed.model,meshed.mesh,2,1)
    @test validate(projected).ok
    @test projected.physical_names==parsed.model.physical_names
    @test mixed_crc(projected).sha==
          "1489fd244841d079350f33439a97b6f33205bcff32e4b99ac672e79a4387eaca"

    expressions=_execute_point_mesh_size_source(raw"""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        MeshSize {1:4}=0.9;
        selected[]={2,4,3}; MeshSize {selected[{0,1}]}=0.5;
        base=1; MeshSize {Max(1,base)+2}=0.6;
        Characteristic Length {base}=0.7;
        Point(newp)={2,2,0,1}; MeshSize {newp-1}=0.8;
        """)
    @test sort!(collect(keys(expressions.model.points)))==collect(1:5)
    @test expressions.model.point_size==Dict(
        1=>0.7,2=>0.5,3=>0.6,4=>0.5,5=>0.8)

    invalid_sources=(
        "Point(1)={0,0,0,1}; MeshSize {1,99}=0.5;"=>"unknown Point[99]",
        "MeshSize {1}=0.5; Point(1)={0,0,0,1};"=>"unknown Point[1]",
        "Point(1)={0,0,0,1}; MeshSize {1}=0;"=>"finite and positive",
        "Point(1)={0,0,0,1}; MeshSize {1}=-0.1;"=>"finite and positive",
        "Point(1)={0,0,0,1}; MeshSize {}=0.5;"=>"entity list is empty",
        "MeshSize {:}=0.5;"=>"matched no explicit modeled points",
        "Point(1)={0,0,0,1}; MeshSize {PointsOf{Surface{1};}}=0.5;"=>
            "topology queries are not supported",
    )
    for (source,message) in invalid_sources
        err=_point_mesh_size_error(source)
        @test err isa ArgumentError
        @test occursin(message,sprint(showerror,err))
    end

    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end
