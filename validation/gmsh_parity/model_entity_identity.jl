#!/usr/bin/env julia
# P6: model entity names and live-reference retagging vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: bounding_box, mesh_crc, ntets, validate
using Tessella.Model: model_periodic_constraints

const API=Tessella.API

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

dim_tags(values)=sort!(Tuple{Int,Int}[
    (Int(dimension),Int(tag)) for (dimension,tag) in values])
tags(values)=sort!(Int.(values))

function tessella_surface_identity()
    API.finalize()
    return try
        API.initialize()
        for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                          (3,1.0,1.0),(4,0.0,1.0),
                          (5,0.25,0.25),(6,0.75,0.25))
            API.model.add_point(x,y,0.0;tag=tag)
        end
        for (tag,first_point,last_point) in
                ((1,1,2),(2,2,3),(3,3,4),(4,4,1),(5,5,6))
            API.model.add_line(first_point,last_point;tag=tag)
        end
        API.model.add_curve_loop([1,2,3,4];tag=1)
        API.model.add_plane_surface([1];tag=1)
        API.model.embed(0,[5],2,1)
        API.model.embed(1,[5],2,1)
        API.model.add_physical_group(0,[1,5];tag=11,name="vertices")
        API.model.add_physical_group(1,[1,5];tag=12,name="curves")
        API.model.add_physical_group(2,[1];tag=13,name="sheet")

        API.model.set_entity_name(0,1,"origin")
        API.model.set_entity_name(1,1,"boundary")
        API.model.set_entity_name(2,1,"surface")
        API.model.set_entity_name(0,2,"shared")
        API.model.set_entity_name(1,2,"shared")
        API.model.remove_entity_name("shared")
        API.model.set_entity_name(1,99,"preloaded")

        API.model.set_tag(0,1,101)
        API.model.set_tag(0,5,105)
        API.model.set_tag(1,1,201)
        API.model.set_tag(1,5,205)
        API.model.set_tag(2,1,301)
        model=API.CURRENT[]
        model===nothing && error("Tessella API lost its initialized model")

        positive_rejected=try
            API.model.set_tag(0,101,0)
            false
        catch err
            err isa ArgumentError || rethrow()
            true
        end
        dimension_rejected=try
            API.model.set_tag(4,101,102)
            false
        catch err
            err isa ArgumentError || rethrow()
            true
        end

        return (
            entities=API.model.get_entities(),
            line_boundary=dim_tags(API.model.get_boundary(
                [(1,201)],false,true,false)),
            surface_boundary=dim_tags(API.model.get_boundary(
                [(2,301)],false,true,false)),
            physical=(
                API.model.get_entities_for_physical_group(0,11),
                API.model.get_entities_for_physical_group(1,12),
                API.model.get_entities_for_physical_group(2,13)),
            embedded=sort!(copy(model.embeds[(2,301)])),
            removed_names=(API.model.get_entity_name(0,2),
                           API.model.get_entity_name(1,2)),
            names=(old=(API.model.get_entity_name(0,1),
                        API.model.get_entity_name(1,1),
                        API.model.get_entity_name(2,1)),
                   new=(API.model.get_entity_name(0,101),
                        API.model.get_entity_name(1,201),
                        API.model.get_entity_name(2,301)),
                   missing=API.model.get_entity_name(1,99)),
            positive_rejected=positive_rejected,
            dimension_rejected=dimension_rejected,
        )
    finally
        API.finalize()
    end
end

function gmsh_surface_identity()
    gmsh.clear()
    gmsh.model.add("model_entity_identity_surface")
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,1.0,1.0),(4,0.0,1.0),
                      (5,0.25,0.25),(6,0.75,0.25))
        gmsh.model.geo.addPoint(x,y,0.0,1.0,tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1),(5,5,6))
        gmsh.model.geo.addLine(first_point,last_point,tag)
    end
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0,[5],2,1)
    gmsh.model.mesh.embed(1,[5],2,1)
    gmsh.model.addPhysicalGroup(0,[1,5],11,"vertices")
    gmsh.model.addPhysicalGroup(1,[1,5],12,"curves")
    gmsh.model.addPhysicalGroup(2,[1],13,"sheet")

    gmsh.model.setEntityName(0,1,"origin")
    gmsh.model.setEntityName(1,1,"boundary")
    gmsh.model.setEntityName(2,1,"surface")
    gmsh.model.setEntityName(0,2,"shared")
    gmsh.model.setEntityName(1,2,"shared")
    gmsh.model.removeEntityName("shared")
    gmsh.model.setEntityName(1,99,"preloaded")

    gmsh.model.setTag(0,1,101)
    gmsh.model.setTag(0,5,105)
    gmsh.model.setTag(1,1,201)
    gmsh.model.setTag(1,5,205)
    gmsh.model.setTag(2,1,301)
    return (
        entities=dim_tags(gmsh.model.getEntities()),
        line_boundary=dim_tags(gmsh.model.getBoundary(
            [(1,201)],false,true,false)),
        surface_boundary=dim_tags(gmsh.model.getBoundary(
            [(2,301)],false,true,false)),
        physical=(
            tags(gmsh.model.getEntitiesForPhysicalGroup(0,11)),
            tags(gmsh.model.getEntitiesForPhysicalGroup(1,12)),
            tags(gmsh.model.getEntitiesForPhysicalGroup(2,13))),
        embedded=dim_tags(gmsh.model.mesh.getEmbedded(2,301)),
        removed_names=(gmsh.model.getEntityName(0,2),
                       gmsh.model.getEntityName(1,2)),
        names=(old=(gmsh.model.getEntityName(0,1),
                    gmsh.model.getEntityName(1,1),
                    gmsh.model.getEntityName(2,1)),
               new=(gmsh.model.getEntityName(0,101),
                    gmsh.model.getEntityName(1,201),
                    gmsh.model.getEntityName(2,301)),
               missing=gmsh.model.getEntityName(1,99)),
    )
end

const TRANSLATE_DOWN=(1.0,0.0,0.0,0.0,
                      0.0,1.0,0.0,-1.0,
                      0.0,0.0,1.0,0.0,
                      0.0,0.0,0.0,1.0)

function tessella_periodic_identity()
    API.finalize()
    return try
        API.initialize()
        for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                          (3,0.0,1.0),(4,1.0,1.0),
                          (5,0.0,2.0),(6,1.0,2.0))
            API.model.add_point(x,y,0.0;tag=tag)
        end
        for tag in 1:3
            API.model.add_line(2tag-1,2tag;tag=tag)
        end
        API.mesh.set_periodic(1,[1,2],[2,3],TRANSLATE_DOWN)
        API.model.set_tag(1,2,102)
        API.model.set_tag(1,1,101)
        API.model.set_tag(1,3,103)
        model=API.CURRENT[]
        model===nothing && error("Tessella API lost its periodic model")
        [(Int(constraint.slave_entity),Int(constraint.master_entity))
         for constraint in model_periodic_constraints(model)]
    finally
        API.finalize()
    end
end

function gmsh_periodic_identity()
    gmsh.clear()
    gmsh.model.add("model_entity_identity_periodic")
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,0.0,1.0),(4,1.0,1.0),
                      (5,0.0,2.0),(6,1.0,2.0))
        gmsh.model.geo.addPoint(x,y,0.0,1.0,tag)
    end
    for tag in 1:3
        gmsh.model.geo.addLine(2tag-1,2tag,tag)
    end
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setPeriodic(1,[1,2],[2,3],collect(TRANSLATE_DOWN))
    gmsh.model.setTag(1,2,102)
    gmsh.model.setTag(1,1,101)
    gmsh.model.setTag(1,3,103)
    masters=Int.(gmsh.model.mesh.getPeriodic(1,[101,102]))
    return collect(zip([101,102],masters))
end

function tessella_volume_identity()
    API.finalize()
    return try
        API.initialize()
        API.model.add_box(0,0,0,1,1,1;tag=1)
        API.model.add_physical_group(3,[1];tag=21,name="domain")
        API.model.set_entity_name(3,1,"box")
        before=API.mesh.generate(3)
        validate(before).ok || error("Tessella pre-retag box mesh is invalid")
        API.model.set_tag(3,1,401)
        after=API.mesh.generate(3)
        validate(after).ok || error("Tessella retagged box mesh is invalid")
        mesh_crc(before)==mesh_crc(after) || error(
            "Tessella box mesh changed after retagging")
        return (
            entities=API.model.get_entities(3),
            physical=API.model.get_entities_for_physical_group(3,21),
            bbox=bounding_box(after),
            name=(old=API.model.get_entity_name(3,1),
                  new=API.model.get_entity_name(3,401)),
            tets=ntets(after),
            crc=mesh_crc(after).sha,
        )
    finally
        API.finalize()
    end
end

function gmsh_volume_identity()
    gmsh.clear()
    gmsh.model.add("model_entity_identity_volume")
    gmsh.model.occ.addBox(0,0,0,1,1,1,1)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(3,[1],21,"domain")
    gmsh.model.setEntityName(3,1,"box")
    gmsh.model.setTag(3,1,401)
    gmsh.option.setNumber("Mesh.MeshSizeMax",0.4)
    gmsh.model.mesh.generate(3)
    types,element_tags,_=gmsh.model.mesh.getElements(3,401)
    isempty(types) && error("Gmsh retagged box has no volume elements")
    return (
        entities=dim_tags(gmsh.model.getEntities(3)),
        physical=tags(gmsh.model.getEntitiesForPhysicalGroup(3,21)),
        bbox=Tuple(Float64.(gmsh.model.getBoundingBox(3,401))),
        name=(old=gmsh.model.getEntityName(3,1),
              new=gmsh.model.getEntityName(3,401)),
        elements=sum(length,element_tags;init=0),
    )
end

function gmsh_tag_divergences()
    function point_model()
        gmsh.clear()
        gmsh.model.add("model_entity_identity_tag_domain")
        gmsh.model.geo.addPoint(0,0,0,1.0,1)
        gmsh.model.geo.synchronize()
    end
    point_model()
    gmsh.model.setTag(0,1,0)
    zero=dim_tags(gmsh.model.getEntities(0))
    point_model()
    gmsh.model.setTag(0,1,-1)
    negative=dim_tags(gmsh.model.getEntities(0))
    point_model()
    gmsh.model.setTag(4,1,2)
    invalid_dimension=dim_tags(gmsh.model.getEntities(0))

    gmsh.clear()
    gmsh.model.add("model_entity_identity_errors")
    gmsh.model.geo.addPoint(0,0,0,1.0,1)
    gmsh.model.geo.addPoint(1,0,0,1.0,2)
    gmsh.model.geo.synchronize()
    rejects=map(((1,1),(99,3),(1,2))) do (source,target)
        try
            gmsh.model.setTag(0,source,target)
            false
        catch
            true
        end
    end
    return (;zero,negative,invalid_dimension,rejects=Tuple(rejects))
end

tessella_surface=tessella_surface_identity()
tessella_periodic=tessella_periodic_identity()
tessella_volume=tessella_volume_identity()

gmsh.initialize(["gmsh","-v","0"])
try
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "model-identity differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh_surface=gmsh_surface_identity()
    gmsh_periodic=gmsh_periodic_identity()
    gmsh_volume=gmsh_volume_identity()
    gmsh_divergences=gmsh_tag_divergences()

    for field in (:entities,:line_boundary,:surface_boundary,:physical,
                  :embedded,:removed_names)
        getproperty(tessella_surface,field)==getproperty(gmsh_surface,field) ||
            error("surface identity $field differs: " *
                  "Tessella=$(getproperty(tessella_surface,field)) " *
                  "Gmsh=$(getproperty(gmsh_surface,field))")
    end
    tessella_surface.names.old==("","","") || error(
        "Tessella kept old name slots after retagging")
    tessella_surface.names.new==("origin","boundary","surface") || error(
        "Tessella did not migrate entity names")
    tessella_surface.names.missing=="" || error(
        "Tessella stored a name for a missing entity")
    gmsh_surface.names.old==("origin","boundary","surface") || error(
        "Gmsh 4.15.2 no longer keeps entity names on old tags")
    gmsh_surface.names.new==("","","") || error(
        "Gmsh 4.15.2 unexpectedly migrated entity names")
    gmsh_surface.names.missing=="preloaded" || error(
        "Gmsh 4.15.2 no longer stores missing-entity names")
    tessella_surface.positive_rejected || error(
        "Tessella accepted a nonpositive model tag")
    tessella_surface.dimension_rejected || error(
        "Tessella accepted an invalid model dimension")

    tessella_periodic==gmsh_periodic==[(101,102),(102,103)] || error(
        "periodic identity differs: Tessella=$tessella_periodic " *
        "Gmsh=$gmsh_periodic")
    tessella_volume.entities==gmsh_volume.entities==[(3,401)] || error(
        "volume entity retagging differs")
    tessella_volume.physical==gmsh_volume.physical==[401] || error(
        "volume Physical membership did not follow its entity")
    tessella_volume.name==(old="",new="box") || error(
        "Tessella volume name did not migrate")
    gmsh_volume.name==(old="box",new="") || error(
        "Gmsh 4.15.2 volume name ownership changed")
    tessella_bbox=(tessella_volume.bbox[1]...,tessella_volume.bbox[2]...)
    all(isapprox(tessella_value,gmsh_value;atol=2e-7,rtol=0)
        for (tessella_value,gmsh_value) in zip(tessella_bbox,gmsh_volume.bbox)) ||
        error("retagged box bounds differ: Tessella=$tessella_bbox " *
              "Gmsh=$(gmsh_volume.bbox)")
    tessella_volume.tets>0 && gmsh_volume.elements>0 || error(
        "retagged volume meshing produced no elements")

    gmsh_divergences.zero==[(0,0)] || error(
        "Gmsh 4.15.2 no longer accepts entity tag 0")
    gmsh_divergences.negative==[(0,-1)] || error(
        "Gmsh 4.15.2 no longer accepts negative entity tags")
    gmsh_divergences.invalid_dimension==[(0,1)] || error(
        "Gmsh 4.15.2 invalid-dimension setTag behavior changed")
    gmsh_divergences.rejects==(true,true,true) || error(
        "Gmsh source/collision setTag errors changed: $(gmsh_divergences.rejects)")

    println("GMSH_PARITY_MODEL_IDENTITY_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) topology_refs=18 periodic_pairs=2 " *
            "volume_crc=$(tessella_volume.crc) " *
            "bounded_divergences=positive_existing_names_migrate")
finally
    gmsh.finalize()
end
