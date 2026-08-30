#!/usr/bin/env julia
# P6: native entity type, plane-property, and nonpartition metadata vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.Model: model_entities, model_entity_type, model_entity_properties,
                      model_parent, model_number_of_partitions, model_partitions,
                      model_set_tag!

function find_gmsh_api()
    explicit=get(ENV,"GMSH_JULIA_API","")
    !isempty(explicit) && isfile(explicit) && return explicit
    executable=Sys.which("gmsh")
    executable===nothing && error("gmsh is not on PATH")
    prefix=dirname(dirname(realpath(executable)))
    for candidate in (joinpath(prefix,"lib","gmsh.jl"),
                      "/opt/homebrew/opt/gmsh/lib/gmsh.jl")
        isfile(candidate) && return candidate
    end
    error("could not locate gmsh.jl")
end

include(find_gmsh_api())

function tessella_metadata_tetrahedron()
    model=GeoModel()
    for (tag,x,y,z) in ((1,0.0,0.0,0.0),(2,2.0,0.0,0.0),
                        (3,0.0,3.0,0.0),(4,0.0,0.0,4.0))
        add_point!(model,x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,1),(4,1,4),(5,2,4),(6,3,4))
        add_line!(model,first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((1,[1,2,3]),(2,[1,5,-4]),
                         (3,[2,6,-5]),(4,[3,4,-6]))
        add_curve_loop!(model,curves;tag=tag)
        add_plane_surface!(model,[tag];tag=tag)
    end
    add_surface_loop!(model,[1,-2,3,-4];tag=1)
    add_volume!(model,[1];tag=1)
    return model
end

function gmsh_metadata_tetrahedron!()
    gmsh.clear()
    gmsh.model.add("model_entity_metadata_tetrahedron")
    for (tag,x,y,z) in ((1,0.0,0.0,0.0),(2,2.0,0.0,0.0),
                        (3,0.0,3.0,0.0),(4,0.0,0.0,4.0))
        gmsh.model.geo.addPoint(x,y,z,1.0,tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,1),(4,1,4),(5,2,4),(6,3,4))
        gmsh.model.geo.addLine(first_point,last_point,tag)
    end
    for (tag,curves) in ((1,[1,2,3]),(2,[1,5,-4]),
                         (3,[2,6,-5]),(4,[3,4,-6]))
        gmsh.model.geo.addCurveLoop(curves,tag)
        gmsh.model.geo.addPlaneSurface([tag],tag)
    end
    gmsh.model.geo.addSurfaceLoop([1,-2,3,-4],1)
    gmsh.model.geo.addVolume([1],1)
    gmsh.model.geo.synchronize()
    return nothing
end

function tessella_metadata_primitives()
    model=GeoModel()
    add_box!(model,-2,1,3,4,5,6;tag=1)
    add_cylinder!(model,10,20,30,2,3,6,4;tag=2)
    add_sphere!(model,-10,-20,-30,5;tag=3)
    add_cone!(model,1,2,3,-2,4,5,6,2;tag=4)
    return model
end

function gmsh_metadata_primitives!()
    gmsh.clear()
    gmsh.model.add("model_entity_metadata_primitives")
    gmsh.model.occ.addBox(-2,1,3,4,5,6,1)
    gmsh.model.occ.addCylinder(10,20,30,2,3,6,4,2)
    gmsh.model.occ.addSphere(-10,-20,-30,5,3)
    gmsh.model.occ.addCone(1,2,3,-2,4,5,6,2,4)
    gmsh.model.occ.synchronize()
    return nothing
end

function compare_metadata(model,entity,label)
    entity_type=model_entity_type(model,entity...)
    reference_type=gmsh.model.getEntityType(entity...)
    entity_type==reference_type || error(
        "$label type differs for $entity: Tessella=$entity_type Gmsh=$reference_type")
    entity_type==gmsh.model.getType(entity...) || error(
        "$label deprecated type alias differs for $entity")
    integers,reals=model_entity_properties(model,entity...)
    reference_integers,reference_reals=gmsh.model.getEntityProperties(entity...)
    integers==Int.(reference_integers) || error(
        "$label integer properties differ for $entity")
    if entity[1]==2
        length(reals)==4 && length(reference_reals)==4 || error(
            "$label Plane properties have the wrong length for $entity")
        reference_scale=hypot(reference_reals[1:3]...)
        reference_scale>0 || error(
            "$label Gmsh Plane normal is degenerate for $entity")
        reference=Float64[
            reference_reals[index]/reference_scale for index in 1:4]
        if sum(reals[index]*reference[index] for index in 1:3)<0
            reference .*= -1
        end
        all(isapprox.(reals,reference;rtol=64eps(Float64),atol=64eps(Float64))) ||
            error("$label Plane properties differ for $entity: " *
                  "Tessella=$reals normalized_Gmsh=$reference")
    else
        isempty(reals) && isempty(reference_reals) || error(
            "$label real properties differ for $entity")
    end
    parent=model_parent(model,entity...)
    reference_parent=Tuple(Int.(gmsh.model.getParent(entity...)))
    parent==reference_parent || error(
        "$label parent differs for $entity: Tessella=$parent Gmsh=$reference_parent")
    partitions=model_partitions(model,entity...)
    reference_partitions=Int.(gmsh.model.getPartitions(entity...))
    partitions==reference_partitions || error(
        "$label partitions differ for $entity: " *
        "Tessella=$partitions Gmsh=$reference_partitions")
    return nothing
end

gmsh.initialize(["gmsh","-v","0"])
try
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "model-metadata differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)

    tessella=tessella_metadata_tetrahedron()
    gmsh_metadata_tetrahedron!()
    entities=model_entities(tessella)
    entities==sort!(Tuple{Int,Int}[(Int(dim),Int(tag)) for (dim,tag) in
                                   gmsh.model.getEntities()]) ||
        error("explicit metadata entity sets differ")
    model_number_of_partitions(tessella)==gmsh.model.getNumberOfPartitions() ||
        error("explicit model partition count differs")
    for entity in entities
        compare_metadata(tessella,entity,"explicit model")
    end
    model_set_tag!(tessella,3,1,10)
    gmsh.model.setTag(3,1,10)
    compare_metadata(tessella,(3,10),"retagged model")

    primitives=tessella_metadata_primitives()
    gmsh_metadata_primitives!()
    model_number_of_partitions(primitives)==gmsh.model.getNumberOfPartitions() ||
        error("primitive model partition count differs")
    for tag in 1:4
        compare_metadata(primitives,(3,tag),"primitive model")
    end

    cases=length(entities)+1+4
    println("GMSH_PARITY_MODEL_METADATA_OK gmsh=",gmsh.GMSH_API_VERSION,
            " cases=",cases," type_queries=",2*cases,
            " property_queries=",cases,
            " parent_queries=",cases," partition_queries=",cases,
            " partition_counts=2",
            " bounded_divergences=implicit_primitive_subentities_no_partitions")
finally
    gmsh.finalize()
end
