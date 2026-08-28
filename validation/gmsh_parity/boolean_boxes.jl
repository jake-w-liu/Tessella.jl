#!/usr/bin/env julia
# P6: BooleanDifference geometry, snapshot ownership, and Delete lifecycle vs Gmsh.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: bounding_box, ntets, tet_volume, node
using Tessella.Model: model_entities, model_entities_for_physical_group,
                      model_physical_groups, model_physical_name

const GEO=joinpath(@__DIR__,"boolean_boxes.geo")
result=execute_geo(GEO; mesh_dim=3)
mesh=result.mesh
mesh===nothing && error("Tessella boolean boxes produced no mesh")
ntets(mesh)>0 || error("Tessella boolean boxes have no tets")
length(result.model.volumes)==1 || error("Tessella boolean operands were not deleted")
function tessella_volume(mesh)
    return sum(tet_volume(node(mesh,mesh.tets[1,t]),node(mesh,mesh.tets[2,t]),
                          node(mesh,mesh.tets[3,t]),node(mesh,mesh.tets[4,t]))
               for t in 1:ntets(mesh))
end

V=tessella_volume(mesh)
abs(V-1)<=1e-12 || error("Tessella boolean-difference volume $V != 1")

retained=GeoModel()
add_box!(retained,0,0,0,2,1,1;tag=1)
add_box!(retained,0,0,0,1,1,1;tag=2)
boolean_volumes!(retained,:difference,1,2;tag=3)
retained_before=mesh_model_volume(retained,3)
translate_volume!(retained,1,(10,0,0))
retained_after=mesh_model_volume(retained,3)
tessella_retained=(before=tessella_volume(retained_before),
                   after=tessella_volume(retained_after),
                   before_bbox=bounding_box(retained_before),
                   after_bbox=bounding_box(retained_after))
tessella_retained.before==tessella_retained.after==1.0 || error(
    "Tessella retained-operand Boolean volume changed: $tessella_retained")
tessella_retained.before_bbox==tessella_retained.after_bbox==
    ((1.0,0.0,0.0),(2.0,1.0,1.0)) || error(
        "Tessella retained-operand Boolean bounds changed: $tessella_retained")

deleted=GeoModel()
add_box!(deleted,0,0,0,2,1,1;tag=1)
add_box!(deleted,0,0,0,1,1,1;tag=2)
add_box!(deleted,20,0,0,1,1,1;tag=4)
add_physical_group!(deleted,3,[1,4];tag=10,name="mixed")
add_physical_group!(deleted,3,[2];tag=11,name="tool")
boolean_volumes!(deleted,:difference,1,2;tag=3)
Tessella.Model._remove_volume_entity!(deleted,1)
Tessella.Model._remove_volume_entity!(deleted,2)
add_sphere!(deleted,30,0,0,1;tag=1)
add_box!(deleted,40,0,0,4,1,1;tag=2)
deleted_result=mesh_model_volume(deleted,3)
tessella_deleted=(volume=tessella_volume(deleted_result),
                  bbox=bounding_box(deleted_result),
                  groups=model_physical_groups(deleted),
                  members=model_entities_for_physical_group(deleted,3,10),
                  name=model_physical_name(deleted,3,10),
                  entities=model_entities(deleted,3))
tessella_deleted.volume==1.0 || error(
    "Tessella deleted-tag reuse changed Boolean volume: $tessella_deleted")
tessella_deleted.bbox==((1.0,0.0,0.0),(2.0,1.0,1.0)) || error(
    "Tessella deleted-tag reuse changed Boolean bounds: $tessella_deleted")
tessella_deleted.groups==[(3,10)] || error(
    "Tessella Boolean Delete Physical groups mismatch: $tessella_deleted")
tessella_deleted.members==[4] && tessella_deleted.name=="mixed" || error(
    "Tessella Boolean Delete Physical membership mismatch: $tessella_deleted")
tessella_deleted.entities==[(3,1),(3,2),(3,3),(3,4)] || error(
    "Tessella deleted volume tags were not reusable: $tessella_deleted")

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
gmsh.initialize(["gmsh","-v","0"])
try
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "Boolean parity requires Gmsh API 4.15.2, got $(gmsh.GMSH_API_VERSION)")
    gmsh.open(GEO)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    gmsh.model.mesh.generate(3)
    types, tags, _ = gmsh.model.mesh.getElements(3)
    isempty(types) && error("Gmsh boolean boxes have no volume elements")
    gmsh_tets=sum(length, tags; init=0)
    gmsh_tets>0 || error("Gmsh boolean boxes have no volume elements")

    gmsh.clear()
    gmsh.model.add("retained_operand_snapshot")
    gmsh.model.occ.addBox(0,0,0,2,1,1,1)
    gmsh.model.occ.addBox(0,0,0,1,1,1,2)
    out,_=gmsh.model.occ.cut([(3,1)],[(3,2)],3,false,false)
    out==[(3,3)] || error("Gmsh retained-operand Boolean returned $out")
    gmsh.model.occ.synchronize()
    gmsh_retained_before=gmsh.model.occ.getMass(3,3)
    gmsh_retained_bbox_before=gmsh.model.getBoundingBox(3,3)
    gmsh.model.occ.translate([(3,1)],10,0,0)
    gmsh.model.occ.synchronize()
    gmsh_retained_after=gmsh.model.occ.getMass(3,3)
    gmsh_retained_bbox_after=gmsh.model.getBoundingBox(3,3)
    isapprox(gmsh_retained_before,1.0;atol=1e-12,rtol=0) || error(
        "Gmsh retained-operand Boolean volume is $gmsh_retained_before")
    isapprox(gmsh_retained_after,gmsh_retained_before;atol=1e-12,rtol=0) || error(
        "Gmsh retained-operand Boolean changed to $gmsh_retained_after")
    gmsh_retained_bbox_after==gmsh_retained_bbox_before || error(
        "Gmsh retained-operand Boolean bounds changed")
    isapprox(tessella_retained.after,gmsh_retained_after;atol=1e-12,rtol=0) || error(
        "retained-operand volume mismatch: Tessella=$(tessella_retained.after), Gmsh=$gmsh_retained_after")

    gmsh.clear()
    gmsh.model.add("deleted_operand_snapshot")
    gmsh.model.occ.addBox(0,0,0,2,1,1,1)
    gmsh.model.occ.addBox(0,0,0,1,1,1,2)
    gmsh.model.occ.addBox(20,0,0,1,1,1,4)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(3,[1,4],10,"mixed")
    gmsh.model.addPhysicalGroup(3,[2],11,"tool")
    out,_=gmsh.model.occ.cut([(3,1)],[(3,2)],3,true,true)
    out==[(3,3)] || error("Gmsh deleting Boolean returned $out")
    gmsh.model.occ.synchronize()
    gmsh.model.occ.addSphere(30,0,0,1,1)
    gmsh.model.occ.addBox(40,0,0,4,1,1,2)
    gmsh.model.occ.synchronize()
    gmsh_deleted_volume=gmsh.model.occ.getMass(3,3)
    gmsh_deleted_bbox=gmsh.model.getBoundingBox(3,3)
    gmsh_groups=sort!(Tuple{Int,Int}.((Int(d),Int(t))
        for (d,t) in gmsh.model.getPhysicalGroups()))
    gmsh_members=sort!(Int.(gmsh.model.getEntitiesForPhysicalGroup(3,10)))
    gmsh_name=gmsh.model.getPhysicalName(3,10)
    gmsh_entities=sort!(Tuple{Int,Int}.((Int(d),Int(t))
        for (d,t) in gmsh.model.getEntities(3)))
    isapprox(gmsh_deleted_volume,1.0;atol=1e-12,rtol=0) || error(
        "Gmsh deleted-tag reuse changed Boolean volume to $gmsh_deleted_volume")
    gmsh_deleted_bbox==gmsh_retained_bbox_before || error(
        "Gmsh deleted-tag reuse changed Boolean bounds")
    gmsh_groups==tessella_deleted.groups || error(
        "Boolean Delete Physical-group mismatch: Tessella=$(tessella_deleted.groups), Gmsh=$gmsh_groups")
    gmsh_members==tessella_deleted.members && gmsh_name==tessella_deleted.name || error(
        "Boolean Delete Physical membership mismatch")
    gmsh_entities==tessella_deleted.entities || error(
        "deleted-tag reuse entity mismatch: Tessella=$(tessella_deleted.entities), Gmsh=$gmsh_entities")
    isapprox(tessella_deleted.volume,gmsh_deleted_volume;atol=1e-12,rtol=0) || error(
        "deleted-tag reuse volume mismatch")

    println("GMSH_PARITY_BOOLEAN_OK gmsh=$(gmsh.GMSH_API_VERSION) " *
            "tessella_volume=1 retained_snapshot=1 deleted_reuse=1 " *
            "physical_members=4 gmsh_tets=$gmsh_tets tessella_tets=$(ntets(mesh))")
finally
    gmsh.finalize()
end
