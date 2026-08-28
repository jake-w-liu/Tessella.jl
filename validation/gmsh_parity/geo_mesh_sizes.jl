#!/usr/bin/env julia
# P6: Point sizing and topology-derived Physical groups vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntris, ntets, triangle_area
using Tessella.Elements: mixed_crc

const GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures","geo_point_mesh_sizes.geo"))
const POINTS_OF_GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures",
    "geo_point_mesh_size_points_of.geo"))
const EXPECTED_SIZES=[0.8,0.4,3*0.8/4,0.4]
const EXPECTED_POINTS_OF_SIZES=[0.4,0.5,0.2,0.4,1.0]
const EXPECTED_PHYSICAL_NAMES=Dict(
    (0,10)=>"corners",(2,20)=>"domain")
const EXPECTED_TOPOLOGY_PHYSICAL=Dict(
    (0,11)=>collect(1:4),(3,12)=>[1],(0,21)=>[1,3],
    (1,22)=>[1,2,3],(2,23)=>collect(1:4),(1,25)=>[2,3,4,5],
    (1,26)=>collect(1:5),(0,27)=>collect(1:4))
const EXPECTED_TOPOLOGY_NAMES=Dict(
    (0,11)=>"vertices",(3,12)=>"domain",(0,21)=>"endpoints",
    (1,22)=>"face boundary",(2,23)=>"skin",(1,25)=>"two face rim",
    (1,26)=>"two face boundary",(0,27)=>"vertices query")

execution=execute_geo(GEO;mesh_dim=2)
mesh=execution.mesh
mesh===nothing && error("Tessella point-size fixture produced no mesh")
validate(mesh).ok || error("Tessella point-size mesh is invalid")
native_sizes=[execution.model.point_size[tag] for tag in 1:4]
native_sizes==EXPECTED_SIZES || error(
    "Tessella point-size constraints changed: $native_sizes")
execution.model.physical_names==EXPECTED_PHYSICAL_NAMES || error(
    "Tessella point-size physical names changed")
nnodes(mesh)==19 && ntris(mesh)==24 || error(
    "Tessella point-size mesh size changed")
mesh_crc(mesh).sha==
    "b3f1bf410e917d050eacceab998b0fdf7b4cd61d1d9f263805b5120c06f1f4df" ||
    error("Tessella point-size mesh CRC changed")
native_area=sum(triangle_area(
    node(mesh,mesh.tris[1,triangle]),node(mesh,mesh.tris[2,triangle]),
    node(mesh,mesh.tris[3,triangle])) for triangle in 1:ntris(mesh))
abs(native_area-1)<=32eps(Float64) || error(
    "Tessella point-size area is $native_area, expected 1")
projected=model_to_mixed(execution.model,mesh,2,1)
validate(projected).ok || error("Tessella point-size projection is invalid")
mixed_crc(projected).sha==
    "b7202dfa1cfb7469e7541c34e2b1bfae404c66f2462abc1953fa0b9374e5a010" ||
    error("Tessella point-size projection CRC changed")

points_of_execution=execute_geo(POINTS_OF_GEO)
points_of_sizes=[points_of_execution.model.point_size[tag] for tag in 1:5]
points_of_sizes==EXPECTED_POINTS_OF_SIZES || error(
    "Tessella PointsOf constraints changed: $points_of_sizes")
points_of_execution.model.physical==EXPECTED_TOPOLOGY_PHYSICAL || error(
    "Tessella topology-derived Physical memberships changed")
points_of_execution.model.physical_names==EXPECTED_TOPOLOGY_NAMES || error(
    "Tessella topology-derived Physical names changed")
Tessella.Model._model_points_of(
    points_of_execution.model,[(3,1)],"PointsOf differential")==collect(1:4) ||
    error("Tessella PointsOf recursive Volume boundary changed")
Tessella.Model._model_boundary(
    points_of_execution.model,[(2,1),(2,2)],"Boundary differential")==
        [1,2,3,1,5,4] || error("Tessella Boundary concatenation changed")
Tessella.Model._model_boundary(
    points_of_execution.model,[(2,1),(2,2)],"CombinedBoundary differential";
    combined=true)==[2,3,4,5] || error(
        "Tessella CombinedBoundary cancellation changed")
points_of_meshed=execute_geo(POINTS_OF_GEO;mesh_dim=3)
nnodes(points_of_meshed.mesh)==6 && ntets(points_of_meshed.mesh)==6 || error(
    "Tessella topology-derived Physical mesh size changed")
points_of_projected=model_to_mixed(
    points_of_meshed.model,points_of_meshed.mesh,3,1)
validate(points_of_projected).ok || error(
    "Tessella topology-derived Physical projection is invalid")
mixed_crc(points_of_projected).sha==
    "608dcd81b4ecab3138fd8da610ec109e1972c799ccb2c9d755917c9901250905" ||
    error("Tessella topology-derived Physical projection CRC changed")

function native_spatial_mesh()
    model=GeoModel()
    sizes=(0.1,0.8,0.8,0.8)
    for (tag,(x,y)) in enumerate(((0.0,0.0),(2.0,0.0),
                                  (2.0,2.0),(0.0,2.0)))
        add_point!(model,x,y,0;tag=tag,mesh_size=sizes[tag])
    end
    for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        add_line!(model,first,last;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_plane_surface!(model,[1];tag=1)
    return mesh_model_surface(model,1;min_angle_deg=20)
end

function native_quadrant_counts(mesh)
    counts=zeros(Int,4)
    for triangle in 1:ntris(mesh)
        nodes=mesh.tris[:,triangle]
        x=sum(mesh.coords[1,nodes])/3
        y=sum(mesh.coords[2,nodes])/3
        counts[(x<1 ? 0 : 1)+(y<1 ? 0 : 2)+1]+=1
    end
    return counts
end

spatial_mesh=native_spatial_mesh()
validate(spatial_mesh).ok || error("Tessella spatial Point-size mesh is invalid")
spatial_native_counts=native_quadrant_counts(spatial_mesh)
nnodes(spatial_mesh)==49 && ntris(spatial_mesh)==72 || error(
    "Tessella spatial Point-size mesh size changed")
spatial_native_counts==[45,8,11,8] || error(
    "Tessella spatial Point-size quadrant counts changed: $spatial_native_counts")
mesh_crc(spatial_mesh).sha==
    "b36070c0c394727d037d35cd0d1944e4807208cc265df4c34ad96c6da82ea1e2" ||
    error("Tessella spatial Point-size mesh CRC changed")
spatial_native_area=sum(triangle_area(
    node(spatial_mesh,spatial_mesh.tris[1,triangle]),
    node(spatial_mesh,spatial_mesh.tris[2,triangle]),
    node(spatial_mesh,spatial_mesh.tris[3,triangle]))
    for triangle in 1:ntris(spatial_mesh))
abs(spatial_native_area-4)<=128eps(Float64) || error(
    "Tessella spatial Point-size area is $spatial_native_area, expected 4")

Tessella.API.finalize()
try
    Tessella.API.initialize()
    Tessella.API.model.add_point(0,0,0;tag=1)
    Tessella.API.model.add_point(1,0,0;tag=2)
    Tessella.API.mesh.set_size([(0,1),(0,2)],0.35)
    [Tessella.API.CURRENT[].point_size[tag] for tag in 1:2]==[0.35,0.35] ||
        error("Tessella API point-size update changed")
finally
    Tessella.API.finalize()
end

function find_gmsh_api()
    explicit=get(ENV,"GMSH_JULIA_API","")
    !isempty(explicit) && isfile(explicit) && return explicit
    executable=Sys.which("gmsh")
    executable===nothing && error("gmsh is not on PATH")
    prefix=dirname(dirname(realpath(executable)))
    for path in (joinpath(prefix,"lib","gmsh.jl"),
                 "/opt/homebrew/opt/gmsh/lib/gmsh.jl")
        isfile(path) && return path
    end
    error("could not locate gmsh.jl")
end

include(find_gmsh_api())

function gmsh_coordinates()
    tags,values,_=gmsh.model.mesh.getNodes()
    return Dict(Int(tag)=>(values[3index-2],values[3index-1],values[3index])
                for (index,tag) in pairs(tags))
end

function gmsh_surface_area(coordinates,connectivity)
    area=0.0
    for offset in 1:3:length(connectivity)
        a=coordinates[Int(connectivity[offset])]
        b=coordinates[Int(connectivity[offset+1])]
        c=coordinates[Int(connectivity[offset+2])]
        ux=b[1]-a[1];uy=b[2]-a[2];uz=b[3]-a[3]
        vx=c[1]-a[1];vy=c[2]-a[2];vz=c[3]-a[3]
        area+=hypot(
            uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx)/2
    end
    return area
end

function gmsh_quadrant_counts(coordinates,connectivity)
    counts=zeros(Int,4)
    for offset in 1:3:length(connectivity)
        points=ntuple(
            slot->coordinates[Int(connectivity[offset+slot-1])],3)
        x=sum(point[1] for point in points)/3
        y=sum(point[2] for point in points)/3
        counts[(x<1 ? 0 : 1)+(y<1 ? 0 : 2)+1]+=1
    end
    return counts
end

gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "point-size differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    gmsh_sizes=gmsh.model.mesh.getSizes([(0,tag) for tag in 1:4])
    gmsh_sizes==EXPECTED_SIZES || error(
        "Gmsh point-size constraints changed: $gmsh_sizes")
    gmsh.parser.getNumber("finePoints")==[2.0,4.0] || error(
        "Gmsh point-size selector list changed")
    for (dim,expected) in (0=>Set(1:4),1=>Set(1:4),2=>Set([1]))
        actual=Set(Int(tag) for (_,tag) in gmsh.model.getEntities(dim))
        actual==expected || error(
            "Gmsh point-size dimension-$dim entities changed")
    end
    names=Dict(
        (Int(dim),Int(tag))=>gmsh.model.getPhysicalName(dim,tag)
        for (dim,tag) in gmsh.model.getPhysicalGroups())
    names==EXPECTED_PHYSICAL_NAMES || error(
        "Gmsh point-size physical names changed")

    gmsh.model.mesh.generate(2)
    coordinates=gmsh_coordinates()
    types,_,node_tags=gmsh.model.mesh.getElements(2,1)
    types==Int32[2] || error(
        "Gmsh point-size surface changed element type")
    connectivity=only(node_tags)
    (!isempty(connectivity) && length(connectivity)%3==0) || error(
        "Gmsh point-size triangle connectivity is malformed")
    gmsh_area=gmsh_surface_area(coordinates,connectivity)
    length(coordinates)==17 && length(connectivity)÷3==22 || error(
        "Gmsh point-size mesh size changed")
    abs(gmsh_area-1)<=32eps(Float64) || error(
        "Gmsh point-size area is $gmsh_area, expected 1")

    gmsh.clear()
    gmsh.open(POINTS_OF_GEO)
    points_of_gmsh_sizes=gmsh.model.mesh.getSizes([(0,tag) for tag in 1:5])
    points_of_gmsh_sizes==EXPECTED_POINTS_OF_SIZES || error(
        "Gmsh PointsOf constraints changed: $points_of_gmsh_sizes")
    gmsh.parser.getNumber("selectedSurfaces")==[-1.0,4.0] || error(
        "Gmsh PointsOf selector list changed")
    points_of_boundary=sort!([Int(tag) for (dim,tag) in
        gmsh.model.getBoundary([(3,1)],false,false,true) if dim==0])
    points_of_boundary==collect(1:4) || error(
        "Gmsh PointsOf recursive Volume boundary changed: $points_of_boundary")
    for ((dimension,tag),expected) in EXPECTED_TOPOLOGY_PHYSICAL
        membership=sort!(Int.(
            gmsh.model.getEntitiesForPhysicalGroup(dimension,tag)))
        membership==expected || error(
            "Gmsh Physical($dimension,$tag) membership changed: $membership")
        gmsh.model.getPhysicalName(dimension,tag)==
            EXPECTED_TOPOLOGY_NAMES[(dimension,tag)] || error(
                "Gmsh Physical($dimension,$tag) name changed")
    end

    primitive_points_of=mktempdir() do directory
        path=joinpath(directory,"primitive-points-of.geo")
        write(path,"""
            SetFactory("OpenCASCADE");
            Box(1) = {0,0,0,1,1,1};
            MeshSize {PointsOf{Volume{1};}} = 0.15;
            """)
        tessella_error=try
            execute_geo(path)
            nothing
        catch err
            err
        end
        tessella_error isa ArgumentError || error(
            "Tessella accepted PointsOf on implicit primitive topology")
        occursin("Volume[1] has no explicit surface-loop topology",
                 sprint(showerror,tessella_error)) || error(
            "Tessella primitive PointsOf blocker changed")
        gmsh.clear()
        gmsh.open(path)
        point_entities=gmsh.model.getEntities(0)
        length(point_entities)==8 || error(
            "Gmsh OCC Box PointsOf point count changed")
        sizes=gmsh.model.mesh.getSizes(point_entities)
        all(==(0.15),sizes) || error(
            "Gmsh OCC Box PointsOf sizes changed: $sizes")
        length(point_entities)
    end

    primitive_boundary=mktempdir() do directory
        path=joinpath(directory,"primitive-boundary.geo")
        write(path,"""
            SetFactory("OpenCASCADE");
            Box(1) = {0,0,0,1,1,1};
            Physical Surface("skin", 31) = CombinedBoundary{Volume{1};};
            """)
        tessella_error=try
            execute_geo(path)
            nothing
        catch err
            err
        end
        tessella_error isa ArgumentError || error(
            "Tessella accepted Boundary on implicit primitive topology")
        occursin("Volume[1] has no explicit surface-loop topology",
                 sprint(showerror,tessella_error)) || error(
            "Tessella primitive Boundary blocker changed")
        gmsh.clear()
        gmsh.open(path)
        surfaces=sort!(Int.(gmsh.model.getEntitiesForPhysicalGroup(2,31)))
        surfaces==collect(1:6) || error(
            "Gmsh OCC Box CombinedBoundary changed: $surfaces")
        length(surfaces)
    end

    gmsh.clear()
    gmsh.model.add("point_size_spatial")
    spatial_sizes=(0.1,0.8,0.8,0.8)
    for (tag,(x,y)) in enumerate(((0.0,0.0),(2.0,0.0),
                                  (2.0,2.0),(0.0,2.0)))
        gmsh.model.geo.addPoint(x,y,0,spatial_sizes[tag],tag)==tag ||
            error("Gmsh spatial Point tag changed")
    end
    for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        gmsh.model.geo.addLine(first,last,tag)==tag ||
            error("Gmsh spatial Curve tag changed")
    end
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)==1 ||
        error("Gmsh spatial Curve Loop tag changed")
    gmsh.model.geo.addPlaneSurface([1],1)==1 ||
        error("Gmsh spatial Surface tag changed")
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.RandomSeed",1)
    gmsh.model.mesh.generate(2)
    spatial_coordinates=gmsh_coordinates()
    spatial_types,_,spatial_node_tags=gmsh.model.mesh.getElements(2,1)
    spatial_types==Int32[2] || error(
        "Gmsh spatial Point-size surface changed element type")
    spatial_connectivity=only(spatial_node_tags)
    spatial_gmsh_counts=gmsh_quadrant_counts(
        spatial_coordinates,spatial_connectivity)
    length(spatial_coordinates)==45 && length(spatial_connectivity)÷3==68 ||
        error("Gmsh spatial Point-size mesh size changed")
    spatial_gmsh_counts==[41,11,8,8] || error(
        "Gmsh spatial Point-size quadrant counts changed: $spatial_gmsh_counts")
    spatial_gmsh_area=gmsh_surface_area(
        spatial_coordinates,spatial_connectivity)
    abs(spatial_gmsh_area-4)<=128eps(Float64) || error(
        "Gmsh spatial Point-size area is $spatial_gmsh_area, expected 4")
    spatial_native_counts[1]>3*maximum(spatial_native_counts[2:4]) || error(
        "Tessella did not localize the fine Point constraint")
    spatial_gmsh_counts[1]>3*maximum(spatial_gmsh_counts[2:4]) || error(
        "Gmsh did not localize the fine Point constraint")
    abs(ntris(spatial_mesh)-length(spatial_connectivity)÷3)<=4 || error(
        "Tessella and Gmsh spatial Point-size element counts diverged")

    gmsh.clear()
    gmsh.model.add("point_size_api")
    gmsh.model.geo.addPoint(0,0,0,1.0,1)==1 || error(
        "Gmsh API Point tag changed")
    gmsh.model.geo.addPoint(1,0,0,1.0,2)==2 || error(
        "Gmsh API Point tag changed")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setSize([(0,1),(0,2)],0.35)
    gmsh.model.mesh.getSizes([(0,1),(0,2)])==[0.35,0.35] || error(
        "Gmsh API point-size update changed")

    permissive_sizes=mktempdir() do directory
        path=joinpath(directory,"permissive.geo")
        write(path,"""
            MeshSize {2} = 0.25;
            Point(1) = {0,0,0,1};
            Point(2) = {1,0,0,1};
            Point(3) = {2,0,0,1};
            MeshSize {1} = -0.2;
            MeshSize {3} = 0;
            MeshSize {99} = 0.1;
            """)
        gmsh.clear()
        gmsh.open(path)
        gmsh.model.mesh.getSizes([(0,1),(0,2),(0,3)])
    end
    permissive_sizes==[-0.2,1.0,0.0] || error(
        "Gmsh permissive point-size behavior changed: $permissive_sizes")

    println("GMSH_PARITY_GEO_MESH_SIZES_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) tessella_nodes=$(nnodes(mesh)) " *
            "tessella_tris=$(ntris(mesh)) gmsh_nodes=$(length(coordinates)) " *
            "gmsh_tris=$(length(connectivity)÷3) area_error=" *
            "$(max(abs(native_area-1),abs(gmsh_area-1))) " *
            "spatial_tessella_tris=$(ntris(spatial_mesh)) " *
            "spatial_gmsh_tris=$(length(spatial_connectivity)÷3) " *
            "spatial_tessella_quadrants=$(join(spatial_native_counts,",")) " *
            "spatial_gmsh_quadrants=$(join(spatial_gmsh_counts,",")) " *
            "points_of_sizes=$(join(points_of_sizes,",")) " *
            "points_of_boundary=$(join(points_of_boundary,",")) " *
            "primitive_hidden_points=$primitive_points_of " *
            "physical_boundary=$(join(EXPECTED_TOPOLOGY_PHYSICAL[(2,23)],",")) " *
            "primitive_hidden_surfaces=$primitive_boundary " *
            "topology_projected_crc=$(mixed_crc(points_of_projected).sha) " *
            "mesh_crc=$(mesh_crc(mesh).sha) projected_crc=$(mixed_crc(projected).sha) " *
            "bounded_contract=positive_existing_points_piecewise_linear_surface_explicit_topology_queries")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
