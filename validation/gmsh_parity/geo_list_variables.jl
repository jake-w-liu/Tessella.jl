#!/usr/bin/env julia
# P6: bounded numeric list variables vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntets, tet_volume
using Tessella.Elements: mixed_crc

const GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures","geo_list_variables.geo"))
const EXPECTED_ENTITIES=Dict(
    0=>Set(vcat(collect(1:8),[101,102])),
    1=>Set(1:12),2=>Set(1:6),3=>Set([1]))
const EXPECTED_PHYSICAL_NAMES=Dict(
    (0,61)=>"corners",(0,65)=>"face probes",(1,62)=>"edges",
    (2,63)=>"boundary",(3,64)=>"domain")

execution=execute_geo(GEO;mesh_dim=3)
mesh=execution.mesh
mesh===nothing && error("Tessella list-variable fixture produced no mesh")
validate(mesh).ok || error("Tessella list-variable mesh is invalid")
nnodes(mesh)==11 && ntets(mesh)==16 || error(
    "Tessella list-variable mesh size changed")
mesh_crc(mesh).sha==
    "2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808" ||
    error("Tessella list-variable mesh CRC changed")
tessella_volume=sum(tet_volume(
    node(mesh,mesh.tets[1,cell]),node(mesh,mesh.tets[2,cell]),
    node(mesh,mesh.tets[3,cell]),node(mesh,mesh.tets[4,cell]))
    for cell in 1:ntets(mesh))
abs(tessella_volume-1)<=1e-12 || error(
    "Tessella list-variable volume is $tessella_volume, expected 1")
for (dim,tags) in EXPECTED_ENTITIES
    entities=dim==0 ? keys(execution.model.points) :
             dim==1 ? keys(execution.model.curves) :
             dim==2 ? keys(execution.model.surfaces) :
                      keys(execution.model.volumes)
    Set(entities)==tags || error(
        "Tessella list-variable dimension-$dim tags changed")
end
execution.model.physical_names==EXPECTED_PHYSICAL_NAMES || error(
    "Tessella list-variable physical names changed")
execution.params.fields[201].options["PointsList"]=="{101, 102}" || error(
    "Tessella list-variable field option changed")
[(constraint.dim,Int(constraint.slave_entity),Int(constraint.master_entity))
 for constraint in model_periodic_constraints(execution.model)]==
    [(2,4,6),(2,5,3)] || error(
        "Tessella list-variable periodic relations changed")
projected=model_to_mixed(execution.model,mesh,3,1)
validate(projected).ok || error(
    "Tessella list-variable projection is invalid")
mixed_crc(projected).sha==
    "27417f652cf93e0d6aad41c2f1b6c65af3751dfb3cb3166432d2e798f25a6493" ||
    error("Tessella list-variable projection CRC changed")

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
        "list-variable differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    expected_variables=Dict(
        "lc"=>[0.45],"pointTags"=>Float64.(vcat(collect(1:8),[101,102])),
        "edgeTags"=>Float64.(1:12),"loopTags"=>Float64.(1:6),
        "surfaceTags"=>Float64.(1:6),"probeTags"=>[101.0,102.0],
        "xSlave"=>[4.0],"xMaster"=>[6.0],"ySlave"=>[5.0],
        "yMaster"=>[3.0],"mutationOracle"=>[6.0,3.0,7.0],
        "mutationCopy"=>[6.0,3.0,7.0,-7.0,-3.0],
        "shellTags"=>[1.0],"volumeTags"=>[1.0],
        "physicalTags"=>Float64.(61:65),"cornerTags"=>Float64.(1:8),
        "fieldTags"=>[201.0])
    for (name,expected) in expected_variables
        gmsh.parser.getNumber(name)==expected || error(
            "Gmsh parser list $name changed")
    end
    for dim in 0:3
        tags=Set(Int(tag) for (_,tag) in gmsh.model.getEntities(dim))
        tags==EXPECTED_ENTITIES[dim] || error(
            "Gmsh list-variable dimension-$dim tags changed")
    end
    names=Dict(
        (Int(dim),Int(tag))=>gmsh.model.getPhysicalName(dim,tag)
        for (dim,tag) in gmsh.model.getPhysicalGroups())
    names==EXPECTED_PHYSICAL_NAMES || error(
        "Gmsh list-variable physical groups changed")
    gmsh.model.mesh.field.getNumbers(201,"PointsList")==[101.0,102.0] || error(
        "Gmsh list-variable field option changed")
    boundary=Set((Int(dim),Int(tag)) for (dim,tag) in
        gmsh.model.getBoundary([(3,1)],false,true,false))
    boundary==Set((2,tag) for tag in 1:6) || error(
        "Gmsh list-variable volume boundary changed")

    gmsh.model.mesh.generate(3)
    coordinates=gmsh_coordinates()
    types,_,node_tags=gmsh.model.mesh.getElements(3,1)
    types==Int32[4] || error(
        "Gmsh list-variable volume changed element type")
    connectivity=only(node_tags)
    (!isempty(connectivity) && length(connectivity)%4==0) || error(
        "Gmsh list-variable tetrahedron connectivity is malformed")
    gmsh_volume=0.0
    for offset in 1:4:length(connectivity)
        gmsh_volume+=abs(signed_tet_volume(
            coordinates[Int(connectivity[offset])],
            coordinates[Int(connectivity[offset+1])],
            coordinates[Int(connectivity[offset+2])],
            coordinates[Int(connectivity[offset+3])]))
    end
    abs(gmsh_volume-1)<=2e-14 || error(
        "Gmsh list-variable volume is $gmsh_volume, expected 1")
    for (slave,master) in ((4,6),(5,3))
        actual_master,slave_nodes,master_nodes,_=
            gmsh.model.mesh.getPeriodicNodes(2,slave)
        actual_master==master || error(
            "Gmsh list-variable Surface[$slave] master changed")
        length(slave_nodes)==length(master_nodes)>0 || error(
            "Gmsh list-variable Surface[$slave] has no node pairs")
    end

    println("GMSH_PARITY_GEO_LIST_VARIABLES_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) tessella_nodes=$(nnodes(mesh)) " *
            "tessella_tets=$(ntets(mesh)) gmsh_nodes=$(length(coordinates)) " *
            "gmsh_tets=$(length(connectivity)÷4) lists=$(length(expected_variables)) " *
            "volume_error=$(max(abs(tessella_volume-1),abs(gmsh_volume-1))) " *
            "mesh_crc=$(mesh_crc(mesh).sha) projected_crc=$(mixed_crc(projected).sha)")
finally
    gmsh.finalize()
end
