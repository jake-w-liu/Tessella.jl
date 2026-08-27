#!/usr/bin/env julia
# P6: bounded geometry expressions and entity ranges vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntets, tet_volume
using Tessella.Elements: mixed_crc

const GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures","geo_geometry_expressions.geo"))
const EXPECTED_ENTITIES=Dict(
    0=>Set(vcat(collect(1:8),100)),
    1=>Set(10:21),
    2=>Set(40:45),
    3=>Set([60]))
const EXPECTED_PHYSICAL_NAMES=Dict(
    (0,70)=>"corners",(0,71)=>"probe",(1,72)=>"edges",
    (2,73)=>"boundary",(3,74)=>"domain")

execution=execute_geo(GEO;mesh_dim=3)
mesh=execution.mesh
mesh===nothing && error("Tessella geometry-expression fixture produced no mesh")
validate(mesh).ok || error("Tessella geometry-expression mesh is invalid")
nnodes(mesh)==9 && ntets(mesh)==12 || error(
    "Tessella geometry-expression mesh size changed")
mesh_crc(mesh).sha==
    "db4a080cdd8b4cdbd080d3ba42b798475d50a4590e67962c32edb8ac69205f24" ||
    error("Tessella geometry-expression mesh CRC changed")
tessella_volume=sum(tet_volume(
    node(mesh,mesh.tets[1,cell]),node(mesh,mesh.tets[2,cell]),
    node(mesh,mesh.tets[3,cell]),node(mesh,mesh.tets[4,cell]))
    for cell in 1:ntets(mesh))
abs(tessella_volume-1)<=1e-12 || error(
    "Tessella geometry-expression volume is $tessella_volume, expected 1")
for (dim,tags) in EXPECTED_ENTITIES
    entities=dim==0 ? keys(execution.model.points) :
             dim==1 ? keys(execution.model.curves) :
             dim==2 ? keys(execution.model.surfaces) :
                      keys(execution.model.volumes)
    Set(entities)==tags || error(
        "Tessella geometry-expression dimension-$dim tags changed")
end
execution.model.physical_names==EXPECTED_PHYSICAL_NAMES || error(
    "Tessella geometry-expression physical names changed")
projected=model_to_mixed(execution.model,mesh,3,60)
validate(projected).ok || error(
    "Tessella geometry-expression projection is invalid")
mixed_crc(projected).sha==
    "89ee7d39873b202e264917e98fce2756b038d3e6f2f06d1bdbce9c46f7e628cd" ||
    error("Tessella geometry-expression projection CRC changed")

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

function signed_tet_volume(a,b,c,d)
    bax=b[1]-a[1];bay=b[2]-a[2];baz=b[3]-a[3]
    cax=c[1]-a[1];cay=c[2]-a[2];caz=c[3]-a[3]
    dax=d[1]-a[1];day=d[2]-a[2];daz=d[3]-a[3]
    return (bax*(cay*daz-caz*day)-
            bay*(cax*daz-caz*dax)+
            baz*(cax*day-cay*dax))/6
end

gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "geometry-expression differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    for (name,value) in (
            "lc"=>0.45,"p"=>1.9,"c"=>10.9,"loop"=>30.9,
            "surface"=>40.9,"shell"=>50.9,"volume"=>60.9)
        gmsh.parser.getNumber(name)==[value] || error(
            "Gmsh parser variable $name changed")
    end
    for dim in 0:3
        tags=Set(Int(tag) for (_,tag) in gmsh.model.getEntities(dim))
        tags==EXPECTED_ENTITIES[dim] || error(
            "Gmsh geometry-expression dimension-$dim tags changed")
    end
    names=Dict(
        (Int(dim),Int(tag))=>gmsh.model.getPhysicalName(dim,tag)
        for (dim,tag) in gmsh.model.getPhysicalGroups())
    names==EXPECTED_PHYSICAL_NAMES || error(
        "Gmsh geometry-expression physical groups changed")
    boundary=Set((Int(dim),Int(tag)) for (dim,tag) in
        gmsh.model.getBoundary([(3,60)],false,true,false))
    boundary==Set((2,tag) for tag in 40:45) || error(
        "Gmsh geometry-expression volume boundary changed")

    gmsh.model.mesh.generate(3)
    coordinates=gmsh_coordinates()
    length(coordinates)>=9 || error(
        "Gmsh geometry-expression mesh omitted model points")
    types,_,node_tags=gmsh.model.mesh.getElements(3,60)
    types==Int32[4] || error(
        "Gmsh geometry-expression volume changed element type")
    connectivity=only(node_tags)
    (!isempty(connectivity) && length(connectivity)%4==0) || error(
        "Gmsh geometry-expression tetrahedron connectivity is malformed")
    gmsh_volume=0.0
    for offset in 1:4:length(connectivity)
        a=coordinates[Int(connectivity[offset])]
        b=coordinates[Int(connectivity[offset+1])]
        c=coordinates[Int(connectivity[offset+2])]
        d=coordinates[Int(connectivity[offset+3])]
        gmsh_volume+=abs(signed_tet_volume(a,b,c,d))
    end
    abs(gmsh_volume-1)<=2e-14 || error(
        "Gmsh geometry-expression volume is $gmsh_volume, expected 1")
    probe_tags,probe_values,_=gmsh.model.mesh.getNodes(0,100)
    length(probe_tags)==1 || error(
        "Gmsh geometry-expression embedded probe is absent")
    probe_error=hypot(
        probe_values[1]-0.5,probe_values[2]-0.5,probe_values[3]-0.5)
    probe_error<=1e-15 || error(
        "Gmsh geometry-expression embedded probe moved by $probe_error")

    println("GMSH_PARITY_GEO_GEOMETRY_EXPRESSIONS_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) tessella_nodes=$(nnodes(mesh)) " *
            "tessella_tets=$(ntets(mesh)) gmsh_nodes=$(length(coordinates)) " *
            "gmsh_tets=$(length(connectivity)÷4) volume_error=" *
            "$(max(abs(tessella_volume-1),abs(gmsh_volume-1))) " *
            "mesh_crc=$(mesh_crc(mesh).sha) projected_crc=$(mixed_crc(projected).sha)")
finally
    gmsh.finalize()
end
