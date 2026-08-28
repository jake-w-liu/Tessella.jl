#!/usr/bin/env julia
# P6: ordered dependency-safe entity removal vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Random
using Tessella
using Tessella.Model: model_entities, model_entities_for_physical_group,
                      model_entity_name, model_periodic_constraints,
                      model_physical_groups, model_physical_name,
                      remove_entities!, set_entity_name!

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

gmsh_entities()=sort!(Tuple{Int,Int}[
    (Int(dimension),Int(tag)) for (dimension,tag) in gmsh.model.getEntities()])
gmsh_groups()=sort!(Tuple{Int,Int}[
    (Int(dimension),Int(tag))
    for (dimension,tag) in gmsh.model.getPhysicalGroups()])

function tessella_two_triangles()
    model=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,0.0,1.0),(4,1.0,1.0))
        add_point!(model,x,y,0.0;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,1),(4,2,4),(5,4,3))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3];tag=1)
    add_plane_surface!(model,[1];tag=1)
    add_curve_loop!(model,[4,5,-2];tag=2)
    add_plane_surface!(model,[2];tag=2)
    return model
end

function gmsh_two_triangles!()
    gmsh.clear()
    gmsh.model.add("model_entity_removal_triangles")
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,0.0,1.0),(4,1.0,1.0))
        gmsh.model.geo.addPoint(x,y,0.0,1.0,tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,1),(4,2,4),(5,4,3))
        gmsh.model.geo.addLine(first_point,last_point,tag)
    end
    gmsh.model.geo.addCurveLoop([1,2,3],1)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.addCurveLoop([4,5,-2],2)
    gmsh.model.geo.addPlaneSurface([2],2)
    gmsh.model.geo.synchronize()
    return nothing
end

function tessella_topology_results()
    blocked=tessella_two_triangles()
    remove_entities!(blocked,[(0,1),(1,1)])

    nonrecursive=tessella_two_triangles()
    remove_entities!(nonrecursive,[(2,1)])
    remove_entities!(nonrecursive,[(1,1),(1,2)])

    recursive=tessella_two_triangles()
    remove_entities!(recursive,[(2,1)],true)

    low_first=tessella_two_triangles()
    remove_entities!(low_first,[(1,1),(2,1)])
    high_first=tessella_two_triangles()
    remove_entities!(high_first,[(2,1),(1,1)])
    return (
        blocked=model_entities(blocked),
        nonrecursive=model_entities(nonrecursive),
        recursive=model_entities(recursive),
        low_first=model_entities(low_first),
        high_first=model_entities(high_first),
    )
end

function gmsh_topology_results()
    gmsh_two_triangles!()
    gmsh.model.removeEntities([(0,1),(1,1)])
    blocked=gmsh_entities()

    gmsh_two_triangles!()
    gmsh.model.removeEntities([(2,1)])
    gmsh.model.removeEntities([(1,1),(1,2)])
    nonrecursive=gmsh_entities()

    gmsh_two_triangles!()
    gmsh.model.removeEntities([(2,1)],true)
    recursive=gmsh_entities()

    gmsh_two_triangles!()
    gmsh.model.removeEntities([(1,1),(2,1)])
    low_first=gmsh_entities()
    gmsh_two_triangles!()
    gmsh.model.removeEntities([(2,1),(1,1)])
    high_first=gmsh_entities()
    return (;blocked,nonrecursive,recursive,low_first,high_first)
end

function tessella_embedded_surface()
    model=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,1.0,1.0),(4,0.0,1.0),
                      (5,0.25,0.25),(6,0.75,0.25))
        add_point!(model,x,y,0.0;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1),(5,5,6))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_plane_surface!(model,[1];tag=1)
    embed!(model,0,[5],2,1)
    embed!(model,1,[5],2,1)
    return model
end

function gmsh_embedded_surface!()
    gmsh.clear()
    gmsh.model.add("model_entity_removal_embedded")
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
    return nothing
end

function tessella_embedding_results()
    blocked=tessella_embedded_surface()
    remove_entities!(blocked,[(0,5),(1,5)])
    recursive=tessella_embedded_surface()
    remove_entities!(recursive,[(2,1)],true)
    return (blocked=model_entities(blocked),
            recursive=model_entities(recursive))
end

function gmsh_embedding_results()
    gmsh_embedded_surface!()
    gmsh.model.removeEntities([(0,5),(1,5)])
    blocked=gmsh_entities()
    gmsh_embedded_surface!()
    gmsh.model.removeEntities([(2,1)],true)
    recursive=gmsh_entities()
    return (;blocked,recursive)
end

function tessella_tetrahedron()
    model=GeoModel()
    for (tag,x,y,z) in ((1,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (3,0.0,1.0,0.0),(4,0.0,0.0,1.0))
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

function gmsh_tetrahedron!()
    gmsh.clear()
    gmsh.model.add("model_entity_removal_volume")
    for (tag,x,y,z) in ((1,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (3,0.0,1.0,0.0),(4,0.0,0.0,1.0))
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

function tessella_metadata_results()
    model=GeoModel()
    add_point!(model,0,0,0;tag=1)
    add_point!(model,1,0,0;tag=2)
    set_entity_name!(model,0,1,"first")
    set_entity_name!(model,0,2,"second")
    add_physical_group!(model,0,[1,2];tag=10,name="nodes")
    add_physical_group!(model,0,[1];tag=11,name="first only")
    remove_entities!(model,[(0,1)])
    after_one=(
        entities=model_entities(model),
        groups=model_physical_groups(model),
        members=model_entities_for_physical_group(model,0,10),
        entity_name=model_entity_name(model,0,1),
        empty_group_name=model_physical_name(model,0,11),
    )
    remove_entities!(model,[(0,2)])
    after_all=(
        entities=model_entities(model),
        groups=model_physical_groups(model),
        entity_name=model_entity_name(model,0,2),
        group_names=(model_physical_name(model,0,10),
                     model_physical_name(model,0,11)),
    )
    return (;after_one,after_all)
end

function gmsh_metadata_results()
    gmsh.clear()
    gmsh.model.add("model_entity_removal_metadata")
    gmsh.model.geo.addPoint(0,0,0,1.0,1)
    gmsh.model.geo.addPoint(1,0,0,1.0,2)
    gmsh.model.geo.synchronize()
    gmsh.model.setEntityName(0,1,"first")
    gmsh.model.setEntityName(0,2,"second")
    gmsh.model.addPhysicalGroup(0,[1,2],10,"nodes")
    gmsh.model.addPhysicalGroup(0,[1],11,"first only")
    gmsh.model.removeEntities([(0,1)])
    after_one=(
        entities=gmsh_entities(),
        groups=gmsh_groups(),
        members=sort!(Int.(gmsh.model.getEntitiesForPhysicalGroup(0,10))),
        entity_name=gmsh.model.getEntityName(0,1),
        empty_group_name=gmsh.model.getPhysicalName(0,11),
    )
    gmsh.model.removeEntities([(0,2)])
    after_all=(
        entities=gmsh_entities(),
        groups=gmsh_groups(),
        entity_name=gmsh.model.getEntityName(0,2),
        group_names=(gmsh.model.getPhysicalName(0,10),
                     gmsh.model.getPhysicalName(0,11)),
    )
    return (;after_one,after_all)
end

function tessella_periodic_cleanup()
    model=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,0.0,1.0),(4,1.0,1.0))
        add_point!(model,x,y,0.0;tag=tag)
    end
    add_line!(model,1,2;tag=1)
    add_line!(model,3,4;tag=2)
    translate_y=(1.0,0.0,0.0,0.0,
                 0.0,1.0,0.0,1.0,
                 0.0,0.0,1.0,0.0,
                 0.0,0.0,0.0,1.0)
    set_periodic!(model,1,[2],[1],translate_y)
    remove_entities!(model,[(1,1)])
    return (entities=model_entities(model),
            constraints=length(model_periodic_constraints(model)))
end

tessella_topology=tessella_topology_results()
tessella_embedding=tessella_embedding_results()
tessella_metadata=tessella_metadata_results()
tessella_volume=tessella_tetrahedron()
remove_entities!(tessella_volume,[(3,1)],true)
tessella_volume_entities=model_entities(tessella_volume)
tessella_periodic=tessella_periodic_cleanup()

gmsh.initialize(["gmsh","-v","0"])
try
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "model-removal differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh_topology=gmsh_topology_results()
    gmsh_embedding=gmsh_embedding_results()
    gmsh_metadata=gmsh_metadata_results()
    gmsh_tetrahedron!()
    gmsh.model.removeEntities([(3,1)],true)
    gmsh_volume_entities=gmsh_entities()

    tessella_topology==gmsh_topology || error(
        "ordered topology removal differs: Tessella=$tessella_topology " *
        "Gmsh=$gmsh_topology")
    tessella_embedding==gmsh_embedding || error(
        "embedding-aware removal differs: Tessella=$tessella_embedding " *
        "Gmsh=$gmsh_embedding")
    tessella_volume_entities==gmsh_volume_entities==Tuple{Int,Int}[] || error(
        "recursive explicit-volume removal did not clear the model")

    for field in (:entities,:groups,:members)
        getproperty(tessella_metadata.after_one,field)==
            getproperty(gmsh_metadata.after_one,field) || error(
                "metadata removal $field differs")
    end
    tessella_metadata.after_all.entities==gmsh_metadata.after_all.entities ||
        error("final metadata entity set differs")
    tessella_metadata.after_all.groups==gmsh_metadata.after_all.groups ||
        error("final metadata Physical groups differ")
    tessella_metadata.after_one.entity_name=="" || error(
        "Tessella retained a removed entity name")
    tessella_metadata.after_one.empty_group_name=="" || error(
        "Tessella retained an empty Physical-group name")
    tessella_metadata.after_all.entity_name=="" || error(
        "Tessella retained the final removed entity name")
    tessella_metadata.after_all.group_names==("","") || error(
        "Tessella retained final empty Physical-group names")
    gmsh_metadata.after_one.entity_name=="first" || error(
        "Gmsh 4.15.2 entity-name removal behavior changed")
    gmsh_metadata.after_one.empty_group_name=="first only" || error(
        "Gmsh 4.15.2 empty Physical-name behavior changed")
    gmsh_metadata.after_all.entity_name=="second" || error(
        "Gmsh 4.15.2 final entity-name behavior changed")
    gmsh_metadata.after_all.group_names==("nodes","first only") || error(
        "Gmsh 4.15.2 final Physical-name behavior changed")

    tessella_periodic.constraints==0 || error(
        "Tessella retained a periodic relation to a removed master")
    gmsh_two_triangles!()
    gmsh.model.removeEntities([(2,99),(4,1),(-1,1)],true)
    gmsh_entities()==tessella_topology.blocked || error(
        "Gmsh missing/invalid removal no-op changed")

    candidates=model_entities(tessella_two_triangles())
    sequences=Vector{Vector{Tuple{Int,Int}}}()
    append!(sequences,([entity] for entity in candidates))
    append!(sequences,([first_entity,second_entity]
                       for first_entity in candidates
                       for second_entity in candidates))
    rng=Xoshiro(0x52454d4f56414c)
    for _ in 1:300
        push!(sequences,[rand(rng,candidates) for _ in 1:3])
    end
    for recursive in (false,true), sequence in sequences
        tessella=tessella_two_triangles()
        remove_entities!(tessella,sequence,recursive)
        gmsh_two_triangles!()
        gmsh.model.removeEntities(sequence,recursive)
        model_entities(tessella)==gmsh_entities() || error(
            "seeded removal differs for recursive=$recursive " *
            "sequence=$sequence")
    end
    seeded_cases=2length(sequences)

    println("GMSH_PARITY_MODEL_REMOVAL_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) topology_cases=5 embedding_cases=2 " *
            "seeded_cases=$seeded_cases recursive_volume=15 " *
            "periodic_cleanup=$(tessella_periodic.constraints) " *
            "bounded_divergences=positive_native_cleanup_names_groups_periodic")
finally
    gmsh.finalize()
end
