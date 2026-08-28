#!/usr/bin/env julia
# P6: point-local .geo/API mesh-size constraints vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntris, triangle_area
using Tessella.Elements: mixed_crc

const GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures","geo_point_mesh_sizes.geo"))
const EXPECTED_SIZES=[0.8,0.4,3*0.8/4,0.4]
const EXPECTED_PHYSICAL_NAMES=Dict(
    (0,10)=>"corners",(2,20)=>"domain")

execution=execute_geo(GEO;mesh_dim=2)
mesh=execution.mesh
mesh===nothing && error("Tessella point-size fixture produced no mesh")
validate(mesh).ok || error("Tessella point-size mesh is invalid")
native_sizes=[execution.model.point_size[tag] for tag in 1:4]
native_sizes==EXPECTED_SIZES || error(
    "Tessella point-size constraints changed: $native_sizes")
execution.model.physical_names==EXPECTED_PHYSICAL_NAMES || error(
    "Tessella point-size physical names changed")
nnodes(mesh)==23 && ntris(mesh)==28 || error(
    "Tessella point-size mesh size changed")
mesh_crc(mesh).sha==
    "f1bfa8a1cc61158cc6293540ad6ce6d7ce48054a616a3df57d93597857f1b089" ||
    error("Tessella point-size mesh CRC changed")
native_area=sum(triangle_area(
    node(mesh,mesh.tris[1,triangle]),node(mesh,mesh.tris[2,triangle]),
    node(mesh,mesh.tris[3,triangle])) for triangle in 1:ntris(mesh))
abs(native_area-1)<=32eps(Float64) || error(
    "Tessella point-size area is $native_area, expected 1")
projected=model_to_mixed(execution.model,mesh,2,1)
validate(projected).ok || error("Tessella point-size projection is invalid")
mixed_crc(projected).sha==
    "1489fd244841d079350f33439a97b6f33205bcff32e4b99ac672e79a4387eaca" ||
    error("Tessella point-size projection CRC changed")

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
            "mesh_crc=$(mesh_crc(mesh).sha) projected_crc=$(mixed_crc(projected).sha) " *
            "bounded_contract=positive_existing_points")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
