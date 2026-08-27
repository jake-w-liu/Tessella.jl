#!/usr/bin/env julia
# P6: bounded dynamic `.geo` tag allocators vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntets, tet_volume
using Tessella.Elements: mixed_crc

const GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures","geo_dynamic_tags.geo"))
const EXPECTED_ENTITIES=Dict(
    0=>Set(1:10),1=>Set(1:12),2=>Set(19:24),3=>Set([26]))
const EXPECTED_PHYSICAL_NAMES=Dict(
    (0,61)=>"corners",(0,65)=>"face probes",(1,62)=>"edges",
    (2,63)=>"boundary",(3,64)=>"domain")

execution=execute_geo(GEO;mesh_dim=3)
mesh=execution.mesh
mesh===nothing && error("Tessella dynamic-tag fixture produced no mesh")
validate(mesh).ok || error("Tessella dynamic-tag mesh is invalid")
nnodes(mesh)==11 && ntets(mesh)==16 || error(
    "Tessella dynamic-tag mesh size changed")
mesh_crc(mesh).sha==
    "2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808" ||
    error("Tessella dynamic-tag mesh CRC changed")
tessella_volume=sum(tet_volume(
    node(mesh,mesh.tets[1,cell]),node(mesh,mesh.tets[2,cell]),
    node(mesh,mesh.tets[3,cell]),node(mesh,mesh.tets[4,cell]))
    for cell in 1:ntets(mesh))
abs(tessella_volume-1)<=1e-12 || error(
    "Tessella dynamic-tag volume is $tessella_volume, expected 1")
for (dim,tags) in EXPECTED_ENTITIES
    entities=dim==0 ? keys(execution.model.points) :
             dim==1 ? keys(execution.model.curves) :
             dim==2 ? keys(execution.model.surfaces) :
                      keys(execution.model.volumes)
    Set(entities)==tags || error(
        "Tessella dynamic-tag dimension-$dim tags changed")
end
Set(keys(execution.model.loops))==Set(13:18) || error(
    "Tessella dynamic-tag Curve Loop tags changed")
Set(keys(execution.model.surface_loops))==Set([25]) || error(
    "Tessella dynamic-tag Surface Loop tag changed")
execution.model.physical_names==EXPECTED_PHYSICAL_NAMES || error(
    "Tessella dynamic-tag physical names changed")
execution.params.fields[1].options["PointsList"]=="{9, 10}" || error(
    "Tessella dynamic-tag field option changed")
[(constraint.dim,Int(constraint.slave_entity),Int(constraint.master_entity))
 for constraint in model_periodic_constraints(execution.model)]==
    [(2,22,24),(2,23,21)] || error(
        "Tessella dynamic-tag periodic relations changed")
projected=model_to_mixed(execution.model,mesh,3,26)
validate(projected).ok || error("Tessella dynamic-tag projection is invalid")
mixed_crc(projected).sha==
    "99aeefc2e269090b518f5896f2388bd419c8fb44d06e4743579d50d61aaedf81" ||
    error("Tessella dynamic-tag projection CRC changed")

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

function check_primitive_allocators()
    source=raw"""
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
        afterRegion = newreg;
        primitiveOracle[] = {
          first, firstAgain, second, third, fourth,
          point1, point2, curve, afterRegion, newp, newreg
        };
        """
    mktempdir() do directory
        path=joinpath(directory,"primitive_allocators.geo")
        write(path,source)
        gmsh.clear()
        gmsh.open(path)
        gmsh.parser.getNumber("primitiveOracle")==
            [1.0,1.0,13.0,16.0,19.0,15.0,16.0,22.0,23.0,17.0,23.0] ||
            error("Gmsh OCC primitive allocator sequence changed")
        Set(Int(tag) for (_,tag) in gmsh.model.getEntities(3))==
            Set([1,13,16,19]) || error("Gmsh OCC primitive volume tags changed")
        point_tags=Set(Int(tag) for (_,tag) in gmsh.model.getEntities(0))
        curve_tags=Set(Int(tag) for (_,tag) in gmsh.model.getEntities(1))
        all(tag in point_tags for tag in (15,16)) && 22 in curve_tags || error(
            "Gmsh OCC primitive explicit entity tags changed")

        tip_source=raw"""
            SetFactory("OpenCASCADE");
            Point(1) = {0,0,0,1};
            Point(2) = {1,0,0,1};
            Point(3) = {0,1,0,1};
            Line(1) = {1,2};
            Line(2) = {2,3};
            Line(3) = {3,1};
            Curve Loop(4) = {1,2,3};
            Plane Surface(100) = {4};
            tip = newv;
            Cone(tip) = {2,0,0,0,0,1,0.5,0};
            coneTipOracle[] = {tip, newreg, newp};
            """
        tip_path=joinpath(directory,"cone_tip_allocators.geo")
        write(tip_path,tip_source)
        gmsh.clear()
        gmsh.open(tip_path)
        gmsh.parser.getNumber("coneTipOracle")==[101.0,103.0,6.0] ||
            error("Gmsh cone-tip allocator sequence changed")
    end
    return nothing
end

gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "dynamic-tag differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    expected_variables=Dict(
        "point1"=>[1.0],"point1Again"=>[1.0],"point8"=>[8.0],
        "edge1"=>[1.0],"edge1Again"=>[1.0],"edge12"=>[12.0],
        "loop1"=>[13.0],"loop1Again"=>[13.0],"loop6"=>[18.0],
        "surface1"=>[19.0],"surface6"=>[24.0],
        "probe1"=>[9.0],"probe2"=>[10.0],
        "shell"=>[25.0],"volume"=>[26.0],"nextRegion"=>[27.0],
        "physicalBase"=>[60.0],"field"=>[1.0],"fieldAgain"=>[1.0],
        "nextField"=>[2.0],
        "allocatorSnapshot"=>
            [11.0,66.0,66.0,66.0,66.0,66.0,66.0,66.0,66.0,2.0])
    for (name,expected) in expected_variables
        gmsh.parser.getNumber(name)==expected || error(
            "Gmsh dynamic-tag parser variable $name changed")
    end
    for dim in 0:3
        tags=Set(Int(tag) for (_,tag) in gmsh.model.getEntities(dim))
        tags==EXPECTED_ENTITIES[dim] || error(
            "Gmsh dynamic-tag dimension-$dim tags changed")
    end
    names=Dict(
        (Int(dim),Int(tag))=>gmsh.model.getPhysicalName(dim,tag)
        for (dim,tag) in gmsh.model.getPhysicalGroups())
    names==EXPECTED_PHYSICAL_NAMES || error(
        "Gmsh dynamic-tag physical groups changed")
    gmsh.model.mesh.field.getNumbers(1,"PointsList")==[9.0,10.0] || error(
        "Gmsh dynamic-tag field option changed")
    boundary=Set((Int(dim),abs(Int(tag))) for (dim,tag) in
        gmsh.model.getBoundary([(3,26)],false,true,false))
    boundary==Set((2,tag) for tag in 19:24) || error(
        "Gmsh dynamic-tag volume boundary changed")

    gmsh.model.mesh.generate(3)
    coordinates=gmsh_coordinates()
    types,_,node_tags=gmsh.model.mesh.getElements(3,26)
    types==Int32[4] || error(
        "Gmsh dynamic-tag volume changed element type")
    connectivity=only(node_tags)
    (!isempty(connectivity) && length(connectivity)%4==0) || error(
        "Gmsh dynamic-tag tetrahedron connectivity is malformed")
    gmsh_volume=0.0
    for offset in 1:4:length(connectivity)
        gmsh_volume+=abs(signed_tet_volume(
            coordinates[Int(connectivity[offset])],
            coordinates[Int(connectivity[offset+1])],
            coordinates[Int(connectivity[offset+2])],
            coordinates[Int(connectivity[offset+3])]))
    end
    abs(gmsh_volume-1)<=2e-14 || error(
        "Gmsh dynamic-tag volume is $gmsh_volume, expected 1")
    for (slave,master) in ((22,24),(23,21))
        actual_master,slave_nodes,master_nodes,_=
            gmsh.model.mesh.getPeriodicNodes(2,slave)
        actual_master==master || error(
            "Gmsh dynamic-tag Surface[$slave] master changed")
        length(slave_nodes)==length(master_nodes)>0 || error(
            "Gmsh dynamic-tag Surface[$slave] has no node pairs")
    end
    gmsh_nodes=length(coordinates)
    gmsh_tets=length(connectivity)÷4
    volume_error=max(abs(tessella_volume-1),abs(gmsh_volume-1))

    check_primitive_allocators()
    println("GMSH_PARITY_GEO_DYNAMIC_TAGS_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) tessella_nodes=$(nnodes(mesh)) " *
            "tessella_tets=$(ntets(mesh)) gmsh_nodes=$gmsh_nodes " *
            "gmsh_tets=$gmsh_tets variables=$(length(expected_variables)) " *
            "volume_error=$volume_error mesh_crc=$(mesh_crc(mesh).sha) " *
            "projected_crc=$(mixed_crc(projected).sha)")
finally
    gmsh.finalize()
end
