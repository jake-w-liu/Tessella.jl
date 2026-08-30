#!/usr/bin/env julia
# P6: native entity presentation, Point-coordinate, and model-attribute state.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.Model: set_entity_visibility!, model_entity_visibility,
                      set_entity_color!, model_entity_color,
                      set_point_coordinates!, set_model_attribute!,
                      model_attribute, model_attribute_names,
                      remove_model_attribute!, model_set_tag!, model_value,
                      model_parametrization_bounds

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

function assert_close(label,actual,reference;atol=1e-12)
    length(actual)==length(reference) || error(
        "$label length differs: Tessella=$(length(actual)) Gmsh=$(length(reference))")
    all(isapprox.(actual,reference;rtol=128eps(Float64),atol=atol)) || error(
        "$label differs: Tessella=$actual Gmsh=$reference")
    return nothing
end

function tessella_state_model()
    model=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,1.0,1.0),(4,0.0,1.0))
        add_point!(model,x,y,0;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_plane_surface!(model,[1];tag=1)
    return model
end

function gmsh_state_model!()
    gmsh.clear()
    gmsh.model.add("model_entity_state")
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,1.0,1.0),(4,0.0,1.0))
        gmsh.model.geo.addPoint(x,y,0,1.0,tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1))
        gmsh.model.geo.addLine(first_point,last_point,tag)
    end
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()
    return nothing
end

gmsh.initialize(["gmsh","-v","0"])
try
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "model-state differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    tessella=tessella_state_model()
    gmsh_state_model!()
    entities=Tuple{Int,Int}[
        (dimension,tag) for dimension in 0:2
        for tag in (dimension==2 ? (1,) : (1,2,3,4))]

    for (dimension,tag) in entities
        model_entity_visibility(tessella,dimension,tag)==
            gmsh.model.getVisibility(dimension,tag) || error(
                "default visibility differs for ($dimension,$tag)")
        model_entity_color(tessella,dimension,tag)==
            gmsh.model.getColor(dimension,tag) || error(
                "default color differs for ($dimension,$tag)")
    end

    set_entity_visibility!(tessella,[(2,1)],0)
    gmsh.model.setVisibility([(2,1)],0,false)
    for (dimension,tag) in entities
        model_entity_visibility(tessella,dimension,tag)==
            gmsh.model.getVisibility(dimension,tag) || error(
                "nonrecursive visibility differs for ($dimension,$tag)")
    end
    set_entity_visibility!(tessella,[(2,1)],-2,true)
    gmsh.model.setVisibility([(2,1)],-2,true)
    for (dimension,tag) in entities
        model_entity_visibility(tessella,dimension,tag)==
            gmsh.model.getVisibility(dimension,tag) || error(
                "recursive visibility differs for ($dimension,$tag)")
    end

    set_entity_color!(tessella,[(2,1)],10,20,30,40,true)
    gmsh.model.setColor([(2,1)],10,20,30,40,true)
    for (dimension,tag) in entities
        model_entity_color(tessella,dimension,tag)==
            gmsh.model.getColor(dimension,tag) || error(
                "recursive color differs for ($dimension,$tag)")
    end

    set_point_coordinates!(tessella,1,-0.25,0,0)
    gmsh.model.setCoordinates(1,-0.25,0,0)
    assert_close("updated Point value",model_value(tessella,0,1,[]),
                 gmsh.model.getValue(0,1,Float64[]))
    line_parameters=[0.0,0.5,1.0]
    assert_close("updated Line value",model_value(tessella,1,1,line_parameters),
                 gmsh.model.getValue(1,1,line_parameters))
    native_bounds=model_parametrization_bounds(tessella,2,1)
    reference_bounds=gmsh.model.getParametrizationBounds(2,1)
    native_bounds==([0.0,-0.25],[1.0,1.0]) || error(
        "native dependent Plane bounds did not follow the Point update")
    reference_bounds==([0.0,0.0],[1.0,1.0]) || error(
        "Gmsh's measured stale dependent Plane bounds changed")

    set_model_attribute!(tessella,"b",["2","x"])
    set_model_attribute!(tessella,"a",["1"])
    set_model_attribute!(tessella,"",String[])
    gmsh.model.setAttribute("b",["2","x"])
    gmsh.model.setAttribute("a",["1"])
    gmsh.model.setAttribute("",String[])
    model_attribute_names(tessella)==gmsh.model.getAttributeNames() ||
        error("attribute names differ")
    for name in model_attribute_names(tessella)
        model_attribute(tessella,name)==gmsh.model.getAttribute(name) ||
            error("attribute $name differs")
    end
    remove_model_attribute!(tessella,"a")
    gmsh.model.removeAttribute("a")
    model_attribute_names(tessella)==gmsh.model.getAttributeNames() ||
        error("attribute removal differs")

    model_set_tag!(tessella,2,1,10)
    gmsh.model.setTag(2,1,10)
    model_entity_visibility(tessella,2,10)==gmsh.model.getVisibility(2,10) ||
        error("retagged visibility differs")
    model_entity_color(tessella,2,10)==gmsh.model.getColor(2,10) ||
        error("retagged color differs")

    println("GMSH_PARITY_MODEL_STATE_OK gmsh=",gmsh.GMSH_API_VERSION,
            " entities=9 visibility_queries=27 color_queries=18 ",
            "coordinate_query_families=3 attributes=3 retagged_entities=1 ",
            "bounded_divergences=strict_dimensions_finite_coordinates_rgba_range_" *
            "nul_free_strings_implicit_boundaries_dynamic_dependent_plane_bounds")
finally
    gmsh.finalize()
end
