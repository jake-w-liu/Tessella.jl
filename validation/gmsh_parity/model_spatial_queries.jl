#!/usr/bin/env julia
# P6: analytical model bounding boxes and containment selection vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Random
using Tessella
using Tessella.Model: model_bounding_box, model_entities,
                      model_entities_in_bounding_box

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
gmsh_box(dimension,tag)=NTuple{6,Float64}(
    gmsh.model.getBoundingBox(dimension,tag))
gmsh_entities_in_box(box,dimension=-1)=sort!(Tuple{Int,Int}[
    (Int(dim),Int(tag)) for (dim,tag) in
    gmsh.model.getEntitiesInBoundingBox(box...,dimension)])

function tessella_tetrahedron()
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

function gmsh_tetrahedron!()
    gmsh.clear()
    gmsh.model.add("model_spatial_tetrahedron")
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

function containment_queries()
    result=NTuple{7,Float64}[
        (0.0,0.0,0.0,2.0,3.0,4.0,-1.0),
        (0.0,0.0,0.0,2.0,0.0,0.0,-1.0),
        (0.0,0.0,0.0,2.0,3.0,0.0,2.0),
        (1.0,1.0,1.0,0.0,0.0,0.0,-1.0),
        (-eps(),-eps(),-eps(),2.0,3.0,4.0,-1.0),
        (0.0,0.0,0.0,prevfloat(2.0),3.0,4.0,-1.0),
    ]
    rng=MersenneTwister(0x5350415449414c)
    for _ in 1:500
        values=ntuple(_->rand(rng)*7-2,6)
        dimension=Float64(rand(rng,-1:3))
        push!(result,(values...,dimension))
    end
    return result
end

function tessella_primitives()
    model=GeoModel()
    add_box!(model,-2,1,3,4,5,6;tag=1)
    add_cylinder!(model,10,20,30,2,3,6,4;tag=2)
    add_sphere!(model,-10,-20,-30,5;tag=3)
    add_cone!(model,1,2,3,-2,4,5,6,2;tag=4)
    return model
end

function gmsh_primitives!()
    gmsh.clear()
    gmsh.model.add("model_spatial_primitives")
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl",0)
    gmsh.model.occ.addBox(-2,1,3,4,5,6,1)
    gmsh.model.occ.addCylinder(10,20,30,2,3,6,4,2)
    gmsh.model.occ.addSphere(-10,-20,-30,5,3)
    gmsh.model.occ.addCone(1,2,3,-2,4,5,6,2,4)
    gmsh.model.occ.synchronize()
    return nothing
end

tessella=tessella_tetrahedron()
queries=containment_queries()

gmsh.initialize(["gmsh","-v","0"])
try
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "model-spatial differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh_tetrahedron!()
    entities=model_entities(tessella)
    entities==gmsh_entities() || error("explicit entity sets differ")
    for entity in entities
        tessella_box=model_bounding_box(tessella,entity...)
        reference_box=gmsh_box(entity...)
        tessella_box==reference_box || error(
            "explicit bounding box differs for $entity: " *
            "Tessella=$tessella_box Gmsh=$reference_box")
    end
    model_bounding_box(tessella,-1,-1)==gmsh_box(-1,-1) || error(
        "whole explicit-model bounding box differs")

    query_signatures=Set{Any}()
    nonempty_queries=0
    for query in queries
        box=NTuple{6,Float64}(query[1:6])
        dimension=Int(query[7])
        tessella_result=model_entities_in_bounding_box(tessella,box...,dimension)
        reference_result=gmsh_entities_in_box(box,dimension)
        tessella_result==reference_result || error(
            "containment query differs for box=$box dim=$dimension: " *
            "Tessella=$tessella_result Gmsh=$reference_result")
        push!(query_signatures,Tuple(tessella_result))
        isempty(tessella_result) || (nonempty_queries+=1)
    end

    primitives=tessella_primitives()
    gmsh_primitives!()
    maximum_padding=0.0
    for tag in 1:4
        tessella_box=model_bounding_box(primitives,3,tag)
        reference_box=gmsh_box(3,tag)
        for axis in 1:3
            reference_box[axis]<=tessella_box[axis]<=reference_box[axis+3] ||
                error("Gmsh OCC bounds do not enclose Tessella Volume[$tag]")
            reference_box[axis]<=tessella_box[axis+3]<=reference_box[axis+3] ||
                error("Gmsh OCC bounds do not enclose Tessella Volume[$tag]")
        end
        padding=maximum(abs(reference_box[index]-tessella_box[index])
                        for index in 1:6)
        padding<=2e-7 || error(
            "Volume[$tag] OCC padding $padding exceeds the bounded differential")
        maximum_padding=max(maximum_padding,padding)

        margin=2e-7
        expanded=(tessella_box[1]-margin,tessella_box[2]-margin,
                  tessella_box[3]-margin,tessella_box[4]+margin,
                  tessella_box[5]+margin,tessella_box[6]+margin)
        tessella_result=model_entities_in_bounding_box(
            primitives,expanded...,3)
        reference_result=gmsh_entities_in_box(expanded,3)
        tessella_result==reference_result || error(
            "primitive containment differs for expanded Volume[$tag] box")
    end

    println("GMSH_PARITY_MODEL_SPATIAL_OK gmsh=",gmsh.GMSH_API_VERSION,
            " exact_boxes=",length(entities)+1,
            " seeded_queries=",length(queries),
            " nonempty_queries=",nonempty_queries,
            " query_signatures=",length(query_signatures),
            " primitive_boxes=4 occ_padding_max=",maximum_padding,
            " bounded_divergences=finite_inputs_strict_dims_implicit_subentities")
finally
    gmsh.finalize()
end
