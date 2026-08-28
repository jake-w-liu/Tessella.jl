#!/usr/bin/env julia
# P6: bounded dynamic `.geo` tag allocators vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntets, tet_volume
using Tessella.Elements: mixed_crc

const GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures","geo_dynamic_tags.geo"))
const SET_MAX_GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures","geo_set_max_tags.geo"))
const EXPECTED_ENTITIES=Dict(
    0=>Set(1:10),1=>Set(1:12),2=>Set(19:24),3=>Set([26]))
const EXPECTED_PHYSICAL_NAMES=Dict(
    (0,61)=>"corners",(0,65)=>"face probes",(1,62)=>"edges",
    (2,63)=>"boundary",(3,64)=>"domain")
const SET_MAX_ENTITIES=Dict(
    0=>Set(101:104),1=>Set(201:206),2=>Set(401:404),3=>Set([601]))
const SET_MAX_PHYSICAL_NAMES=Dict(
    (0,603)=>"corners",(1,604)=>"edges",
    (2,605)=>"boundary",(3,606)=>"domain")

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

set_max_execution=execute_geo(SET_MAX_GEO;mesh_dim=3)
set_max_mesh=set_max_execution.mesh
set_max_mesh===nothing && error("Tessella SetMaxTag fixture produced no mesh")
validate(set_max_mesh).ok || error("Tessella SetMaxTag mesh is invalid")
nnodes(set_max_mesh)==5 && ntets(set_max_mesh)==4 || error(
    "Tessella SetMaxTag mesh size changed")
mesh_crc(set_max_mesh).sha==
    "71ab10cf31fa64d469e1bc3985bd8c50bb240d1cdefaebbc17101bce22e7008b" ||
    error("Tessella SetMaxTag mesh CRC changed")
for (dim,tags) in SET_MAX_ENTITIES
    entities=dim==0 ? keys(set_max_execution.model.points) :
             dim==1 ? keys(set_max_execution.model.curves) :
             dim==2 ? keys(set_max_execution.model.surfaces) :
                      keys(set_max_execution.model.volumes)
    Set(entities)==tags || error("Tessella SetMaxTag dimension-$dim tags changed")
end
set_max_execution.model.physical_names==SET_MAX_PHYSICAL_NAMES || error(
    "Tessella SetMaxTag physical names changed")
set_max_execution.params.fields[1].options["PointsList"]==
    "{101, 102, 103, 104}" || error("Tessella SetMaxTag field option changed")
set_max_volume=sum(tet_volume(
    node(set_max_mesh,set_max_mesh.tets[1,cell]),
    node(set_max_mesh,set_max_mesh.tets[2,cell]),
    node(set_max_mesh,set_max_mesh.tets[3,cell]),
    node(set_max_mesh,set_max_mesh.tets[4,cell]))
    for cell in 1:ntets(set_max_mesh))
abs(set_max_volume-1/6)<=1e-12 || error(
    "Tessella SetMaxTag volume is $set_max_volume, expected 1/6")
set_max_projected=model_to_mixed(
    set_max_execution.model,set_max_mesh,3,601)
validate(set_max_projected).ok || error("Tessella SetMaxTag projection is invalid")
mixed_crc(set_max_projected).sha==
    "aa3127e5c1a302ebdb98890a4d7a62d5cb385e3cf73aa024372f809a7542a45d" ||
    error("Tessella SetMaxTag projection CRC changed")

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
            SetMaxTag Surface(5);
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

function check_set_max_tags()
    gmsh.clear()
    gmsh.open(SET_MAX_GEO)
    gmsh.parser.getNumber("allocatorSnapshot")==
        [105.0,607.0,607.0,607.0,607.0,607.0,607.0,607.0,607.0,2.0] ||
        error("Gmsh SetMaxTag allocator snapshot changed")
    for dim in 0:3
        tags=Set(Int(tag) for (_,tag) in gmsh.model.getEntities(dim))
        tags==SET_MAX_ENTITIES[dim] || error(
            "Gmsh SetMaxTag dimension-$dim tags changed")
    end
    names=Dict(
        (Int(dim),Int(tag))=>gmsh.model.getPhysicalName(dim,tag)
        for (dim,tag) in gmsh.model.getPhysicalGroups())
    names==SET_MAX_PHYSICAL_NAMES || error(
        "Gmsh SetMaxTag physical groups changed")
    gmsh.model.mesh.field.getNumbers(1,"PointsList")==
        [101.0,102.0,103.0,104.0] || error(
            "Gmsh SetMaxTag field option changed")
    boundary=Set((Int(dim),abs(Int(tag))) for (dim,tag) in
        gmsh.model.getBoundary([(3,601)],false,true,false))
    boundary==Set((2,tag) for tag in 401:404) || error(
        "Gmsh SetMaxTag volume boundary changed")

    gmsh.model.mesh.generate(3)
    coordinates=gmsh_coordinates()
    types,_,node_tags=gmsh.model.mesh.getElements(3,601)
    types==Int32[4] || error("Gmsh SetMaxTag volume changed element type")
    connectivity=only(node_tags)
    (!isempty(connectivity) && length(connectivity)%4==0) || error(
        "Gmsh SetMaxTag tetrahedron connectivity is malformed")
    volume=0.0
    for offset in 1:4:length(connectivity)
        volume+=abs(signed_tet_volume(
            coordinates[Int(connectivity[offset])],
            coordinates[Int(connectivity[offset+1])],
            coordinates[Int(connectivity[offset+2])],
            coordinates[Int(connectivity[offset+3])]))
    end
    abs(volume-1/6)<=2e-14 || error(
        "Gmsh SetMaxTag volume is $volume, expected 1/6")

    mktempdir() do directory
        lowering=raw"""
            SetFactory("Built-in");
            Point(10) = {0,0,0,1};
            SetMaxTag Point(20);
            SetMaxTag Point(5);
            lower = newp;
            Point(lower) = {1,0,0,1};
            lowerNext = newp;
            SetFactory("OpenCASCADE");
            occValue = newp;
            SetMaxTag Point(20);
            occRaised = newp;
            SetMaxTag Point(5);
            occPreserved = newp;
            SetFactory("Built-in");
            builtinRaised = newp;
            SetMaxTag Point(4);
            builtinLoweredAgain = newp;
            loweringOracle[] = {
              lower, lowerNext, occValue, occRaised, occPreserved,
              builtinRaised, builtinLoweredAgain
            };
            """
        lowering_path=joinpath(directory,"set_max_tag_lowering.geo")
        write(lowering_path,lowering)
        gmsh.clear()
        gmsh.open(lowering_path)
        lowering_actual=gmsh.parser.getNumber("loweringOracle")
        lowering_expected=[6.0,7.0,7.0,21.0,21.0,21.0,21.0]
        lowering_actual==lowering_expected || error(
            "Gmsh SetMaxTag lowering or factory sharing changed: expected " *
            "$lowering_expected, got $lowering_actual")

        negative=raw"""
            SetFactory("Built-in");
            SetMaxTag Point(-1);
            beforeFactoryActivation = newp;
            SetFactory("OpenCASCADE");
            afterFactoryActivation = newp;
            negativeOracle[] = {
              beforeFactoryActivation, afterFactoryActivation
            };
            """
        negative_path=joinpath(directory,"set_max_tag_negative.geo")
        write(negative_path,negative)
        gmsh.clear()
        gmsh.open(negative_path)
        negative_actual=gmsh.parser.getNumber("negativeOracle")
        negative_actual==[0.0,1.0] || error(
            "Gmsh SetMaxTag factory activation changed: expected [0.0, 1.0], " *
            "got $negative_actual")

        allocator_expression=raw"""
            SetFactory("Built-in");
            SetMaxTag Point(10);
            SetMaxTag Point(newp + 4);
            allocatorExpressionOracle[] = {newp};
            """
        allocator_expression_path=joinpath(
            directory,"set_max_tag_allocator_expression.geo")
        write(allocator_expression_path,allocator_expression)
        gmsh.clear()
        gmsh.open(allocator_expression_path)
        allocator_expression_actual=
            gmsh.parser.getNumber("allocatorExpressionOracle")
        allocator_expression_actual==[16.0] || error(
            "Gmsh SetMaxTag allocator expression changed: expected [16.0], " *
            "got $allocator_expression_actual")

        for (kind,allocator) in (("Curve","newreg"),("Surface","newreg"),
                                 ("Volume","newreg"))
            category_source="""
                SetFactory("OpenCASCADE");
                SetMaxTag $kind(20);
                SetMaxTag $kind(5);
                categoryOracle[] = {$allocator};
                """
            category_path=joinpath(
                directory,"set_max_tag_" * lowercase(kind) * ".geo")
            write(category_path,category_source)
            gmsh.clear()
            gmsh.open(category_path)
            category_actual=gmsh.parser.getNumber("categoryOracle")
            category_actual==[21.0] || error(
                "Gmsh OpenCASCADE SetMaxTag $kind lowering behavior changed: " *
                "expected [21.0], got $category_actual")
        end

        primitive_gap=raw"""
            SetFactory("OpenCASCADE");
            Point(7) = {9,9,9,1};
            SetMaxTag Point(5);
            Box(1) = {0,0,0,1,1,1};
            primitiveGapOracle[] = {newp, newreg};
            """
        primitive_gap_path=joinpath(directory,"set_max_tag_primitive_gap.geo")
        write(primitive_gap_path,primitive_gap)
        gmsh.clear()
        gmsh.open(primitive_gap_path)
        gmsh.parser.getNumber("primitiveGapOracle")==[16.0,13.0] || error(
            "Gmsh SetMaxTag primitive occupied-tag handling changed")
    end
    return (nodes=length(coordinates),tets=length(connectivity)÷4,
            volume_error=abs(volume-1/6))
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
    set_max_gmsh=check_set_max_tags()
    volume_error=max(volume_error,abs(set_max_volume-1/6),
                     set_max_gmsh.volume_error)
    println("GMSH_PARITY_GEO_DYNAMIC_TAGS_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) tessella_nodes=$(nnodes(mesh)) " *
            "tessella_tets=$(ntets(mesh)) gmsh_nodes=$gmsh_nodes " *
            "gmsh_tets=$gmsh_tets variables=$(length(expected_variables)) " *
            "setmax_nodes=$(set_max_gmsh.nodes) " *
            "setmax_tets=$(set_max_gmsh.tets) " *
            "volume_error=$volume_error mesh_crc=$(mesh_crc(mesh).sha) " *
            "projected_crc=$(mixed_crc(projected).sha) " *
            "setmax_crc=$(mixed_crc(set_max_projected).sha)")
finally
    gmsh.finalize()
end
