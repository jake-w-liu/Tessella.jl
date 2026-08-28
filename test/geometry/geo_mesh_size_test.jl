using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntris, ntets, triangle_area
using Tessella.Elements: mixed_crc

const _GEO_POINT_MESH_SIZE_FIXTURE=normpath(joinpath(
    @__DIR__,"..","fixtures","geo_point_mesh_sizes.geo"))
const _GEO_POINT_MESH_SIZE_POINTS_OF_FIXTURE=normpath(joinpath(
    @__DIR__,"..","fixtures","geo_point_mesh_size_points_of.geo"))

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

function _point_mesh_size_square(sizes)
    model=GeoModel()
    for (tag,(x,y)) in enumerate(((0.0,0.0),(2.0,0.0),
                                  (2.0,2.0),(0.0,2.0)))
        add_point!(model,x,y,0;tag=tag,mesh_size=sizes[tag])
    end
    for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        add_line!(model,first,last;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_plane_surface!(model,[1];tag=1)
    return model
end

function _point_mesh_size_quadrants(mesh)
    counts=zeros(Int,4)
    for triangle in 1:ntris(mesh)
        nodes=mesh.tris[:,triangle]
        x=sum(mesh.coords[1,nodes])/3
        y=sum(mesh.coords[2,nodes])/3
        counts[(x<1 ? 0 : 1)+(y<1 ? 0 : 2)+1]+=1
    end
    return counts
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

    topology=execute_geo(_GEO_POINT_MESH_SIZE_POINTS_OF_FIXTURE)
    @test [topology.model.point_size[tag] for tag in 1:5]==
          [0.4,0.5,0.2,0.4,1.0]
    @test topology.model.physical==
          Dict((0,11)=>collect(1:4),(3,12)=>[1])
    @test topology.model.physical_names==
          Dict((0,11)=>"vertices",(3,12)=>"domain")
    for (entities,expected) in (
            ([(0,3)],[3]), ([(1,2)],[2,3]),
            ([(2,1)],[1,2,3]), ([(3,1)],collect(1:4)),
            ([(0,3),(1,4),(2,1),(3,1)],collect(1:4)))
        @test Tessella.Model._model_points_of(
            topology.model,entities,"test PointsOf")==expected
    end
    topology_stable=copy(topology.model.point_size)
    @test_throws ArgumentError Tessella.Model._model_points_of(
        topology.model,[(2,99)],"test PointsOf")
    @test topology.model.point_size==topology_stable

    topology_meshed=execute_geo(
        _GEO_POINT_MESH_SIZE_POINTS_OF_FIXTURE;mesh_dim=3)
    @test validate(topology_meshed.mesh).ok
    @test nnodes(topology_meshed.mesh)==6
    @test ntets(topology_meshed.mesh)==6
    @test mesh_crc(topology_meshed.mesh).sha==
          "cf091ac13ba325f5f68650b192598a5159679b40e19dbea744d66c58962357b5"

    holed_topology=_execute_point_mesh_size_source(raw"""
        Point(1)={0,0,0,1}; Point(2)={2,0,0,1};
        Point(3)={2,2,0,1}; Point(4)={0,2,0,1};
        Point(5)={0.5,0.5,0,1}; Point(6)={1.5,0.5,0,1};
        Point(7)={1.5,1.5,0,1}; Point(8)={0.5,1.5,0,1};
        Point(9)={1,1,0,1};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Line(5)={5,6}; Line(6)={6,7}; Line(7)={7,8}; Line(8)={8,5};
        Curve Loop(1)={1:4}; Curve Loop(2)={5:8};
        Plane Surface(1)={1,2}; Point{9} In Surface{1};
        MeshSize{PointsOf{Surface{:};}}=0.25;
        """)
    @test [holed_topology.model.point_size[tag] for tag in 1:9]==
          vcat(fill(0.25,8),1.0)

    meshed=execute_geo(_GEO_POINT_MESH_SIZE_FIXTURE;mesh_dim=2)
    @test validate(meshed.mesh).ok
    @test nnodes(meshed.mesh)==19
    @test ntris(meshed.mesh)==24
    @test mesh_crc(meshed.mesh).sha==
          "b3f1bf410e917d050eacceab998b0fdf7b4cd61d1d9f263805b5120c06f1f4df"
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
          "b7202dfa1cfb7469e7541c34e2b1bfae404c66f2462abc1953fa0b9374e5a010"

    asymmetric=_point_mesh_size_square([0.1,0.8,0.8,0.8])
    asymmetric_mesh=mesh_model_surface(asymmetric,1;min_angle_deg=20)
    @test validate(asymmetric_mesh).ok
    @test nnodes(asymmetric_mesh)==49
    @test ntris(asymmetric_mesh)==72
    @test _point_mesh_size_quadrants(asymmetric_mesh)==[45,8,11,8]
    @test mesh_crc(asymmetric_mesh).sha==
          "b36070c0c394727d037d35cd0d1944e4807208cc265df4c34ad96c6da82ea1e2"
    for triangle in 1:ntris(asymmetric_mesh)
        nodes=asymmetric_mesh.tris[:,triangle]
        points=ntuple(slot->node(asymmetric_mesh,nodes[slot]),3)
        centroid=(sum(point[1] for point in points)/3,
                  sum(point[2] for point in points)/3)
        target=min(0.8,muladd(0.35,centroid[1]+centroid[2],0.1))
        longest=maximum(hypot(
            points[mod1(slot+1,3)][1]-points[slot][1],
            points[mod1(slot+1,3)][2]-points[slot][2]) for slot in 1:3)
        @test longest<=target
    end

    xs=[0.0,2.0,2.0,0.0];ys=[0.0,0.0,2.0,2.0]
    initial=Tessella.Mesh2D.constrained_delaunay(
        xs,ys,[(1,2),(2,3),(3,4),(4,1)])
    interpolated=Tessella.Model._surface_point_size_field(
        initial,xs,ys,[0.1,0.8,0.8,0.8],1,"mesh_model_surface")
    @test interpolated(0.0,0.0)≈0.1
    @test interpolated(1.0,0.0)≈0.45
    @test interpolated(0.5,0.5)≈0.45
    @test interpolated(1.0,1.0)≈0.8
    extreme=Tessella.Model._surface_point_size_field(
        initial,xs,ys,[nextfloat(0.0),floatmax(Float64),
                       floatmax(Float64),floatmax(Float64)],1,
        "mesh_model_surface")
    @test extreme(0.0,0.0)==nextfloat(0.0)
    @test extreme(2.0,2.0)==floatmax(Float64)
    @test isfinite(extreme(1.0,0.0))
    interpolation_error=try
        interpolated(3.0,3.0)
        nothing
    catch err
        err
    end
    @test interpolation_error isa ErrorException
    @test occursin(
        "Point-size query (3.0,3.0) is outside the initial triangulation of Surface[1]",
        sprint(showerror,interpolation_error))

    generated=_point_mesh_size_square([0.1,0.8,0.8,0.8])
    generated_data=Tessella.Model._surface_pslg(
        generated,1,Dict(1=>[0.0,0.5,1.0]),"mesh_model_surface")
    generated_x,generated_y,generated_sizes=generated_data[1:3]
    midpoint=only(findall(eachindex(generated_x)) do index
        generated_x[index]==1.0 && generated_y[index]==0.0
    end)
    @test generated_sizes[midpoint]≈0.45

    add_point!(generated,1,1,0;tag=5,mesh_size=0.2)
    add_line!(generated,1,3;tag=5)
    embed!(generated,0,[5],2,1)
    embed!(generated,1,[5],2,1)
    coincident_data=Tessella.Model._surface_pslg(
        generated,1,Dict(5=>[0.0,0.5,1.0]),"mesh_model_surface")
    coincident_x,coincident_y,coincident_sizes=coincident_data[1:3]
    center=only(findall(eachindex(coincident_x)) do index
        coincident_x[index]==1.0 && coincident_y[index]==1.0
    end)
    @test coincident_sizes[center]==0.2

    merged=_point_mesh_size_square(fill(0.8,4))
    add_point!(merged,0,0,0;tag=5,mesh_size=0.2)
    add_point!(merged,1,1,0;tag=6,mesh_size=0.8)
    add_line!(merged,5,6;tag=5)
    embed!(merged,1,[5],2,1)
    merged_mesh=mesh_model_surface(merged,1;min_angle_deg=20)
    @test validate(merged_mesh).ok
    @test nnodes(merged_mesh)==39
    @test ntris(merged_mesh)==56
    @test mesh_crc(merged_mesh).sha==
          "0b48fd5481d111e3a2f4eb413b327dd3e31c83d6cb9fa29b7c6d37cd1bb102d5"

    uniform_mesh=mesh_model_surface(
        _point_mesh_size_square(fill(0.8,4)),1;min_angle_deg=20)
    @test validate(uniform_mesh).ok
    @test nnodes(uniform_mesh)==23
    @test ntris(uniform_mesh)==28
    @test mesh_crc(uniform_mesh).sha==
          "f1bfa8a1cc61158cc6293540ad6ce6d7ce48054a616a3df57d93597857f1b089"

    periodic=_point_mesh_size_square([0.4,0.8,0.8,0.4])
    set_periodic!(periodic,1,[2],[4],(
        1.0,0.0,0.0,2.0,
        0.0,1.0,0.0,0.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1.0))
    periodic_mesh=mesh_model_surface(periodic,1)
    @test validate(periodic_mesh).ok
    periodic_nodes=model_periodic_nodes(periodic,periodic_mesh,1,2)
    @test length(periodic_nodes.slave_nodes)==
          length(periodic_nodes.master_nodes)==9
    @test mesh_crc(periodic_mesh).sha==
          "43aa68464111f6c4c6e47faeed3ff94c503597b6bb5420d2916e3c2660215501"

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
        "Point(1)={0,0,0,1}; MeshSize{PointsOf{Surface{99};}}=0.5;"=>
            "unknown Surface[99]",
        "Box(1)={0,0,0,1,1,1}; MeshSize{PointsOf{Volume{1};}}=0.5;"=>
            "Volume[1] has no explicit surface-loop topology",
        "Point(1)={0,0,0,1}; MeshSize{PointsOf{Surface{1}}}=0.5;"=>
            "every PointsOf entity block must end with a semicolon",
        "Point(1)={0,0,0,1}; MeshSize{PointsOf{Surface{};}}=0.5;"=>
            "entity list is empty",
        "Point(1)={0,0,0,1}; MeshSize{PointsOf{}}=0.5;"=>
            "PointsOf must contain at least one entity block",
        "Point(1)={0,0,0,1}; MeshSize{PointsOf{Surface{-2147483648};}}=0.5;"=>
            "magnitude exceeds Int32",
        "Point(1)={0,0,0,1}; MeshSize{PointsOf{Point{1:40000}; Point{1:40000};}}=0.5;"=>
            "PointsOf expands beyond 65536 entities",
        "Point(1)={0,0,0,1}; MeshSize{PointsOf{Curve Loop{1};}}=0.5;"=>
            "supports Point, Curve/Line, Surface, and Volume blocks",
        "Point(1)={0,0,0,1}; MeshSize{Boundary{Point{1};}}=0.5;"=>
            "unsupported topology query",
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
