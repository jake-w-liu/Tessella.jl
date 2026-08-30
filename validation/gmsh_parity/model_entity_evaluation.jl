#!/usr/bin/env julia
# P6: explicit Point, straight-Line, and Plane evaluation vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.Model: model_value, model_derivative, model_second_derivative,
                      model_curvature, model_principal_curvatures, model_normal,
                      model_parametrization, model_parametrization_bounds,
                      model_is_inside, model_closest_point,
                      model_reparametrize_on_surface, model_set_tag!

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

function assert_close(label,actual,reference;rtol=128eps(Float64),atol=1e-12)
    length(actual)==length(reference) || error(
        "$label length differs: Tessella=$(length(actual)) Gmsh=$(length(reference))")
    all(isapprox.(actual,reference;rtol=rtol,atol=atol)) || error(
        "$label differs: Tessella=$actual Gmsh=$reference")
    return nothing
end

function tessella_plane_with_hole()
    model=GeoModel()
    for (tag,x,y,z) in ((1,1.0,2.0,3.0),(2,5.0,2.0,3.0),
                        (3,5.0,6.0,3.0),(4,1.0,6.0,3.0),
                        (5,2.0,3.0,3.0),(6,3.0,3.0,3.0),
                        (7,2.5,4.0,3.0))
        add_point!(model,x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1),
             (5,5,6),(6,6,7),(7,7,5))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_curve_loop!(model,[5,6,7];tag=2)
    add_plane_surface!(model,[1,2];tag=1)
    add_point!(model,0,0,4;tag=8)
    add_point!(model,1,0,4;tag=9)
    add_line!(model,8,9;tag=8)
    return model
end

function gmsh_plane_with_hole!()
    gmsh.clear()
    gmsh.model.add("model_entity_evaluation_plane")
    for (tag,x,y,z) in ((1,1.0,2.0,3.0),(2,5.0,2.0,3.0),
                        (3,5.0,6.0,3.0),(4,1.0,6.0,3.0),
                        (5,2.0,3.0,3.0),(6,3.0,3.0,3.0),
                        (7,2.5,4.0,3.0))
        gmsh.model.geo.addPoint(x,y,z,1.0,tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1),
             (5,5,6),(6,6,7),(7,7,5))
        gmsh.model.geo.addLine(first_point,last_point,tag)
    end
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    gmsh.model.geo.addCurveLoop([5,6,7],2)
    gmsh.model.geo.addPlaneSurface([1,2],1)
    gmsh.model.geo.addPoint(0,0,4,1.0,8)
    gmsh.model.geo.addPoint(1,0,4,1.0,9)
    gmsh.model.geo.addLine(8,9,8)
    gmsh.model.geo.synchronize()
    return nothing
end

function tessella_tilted_plane()
    model=GeoModel()
    for (tag,coordinate) in
            pairs(((1.0,2.0,3.0),(2.0,1.0,3.0),(2.0,2.0,2.0)))
        add_point!(model,coordinate...;tag=tag)
    end
    for (tag,first_point,last_point) in ((1,1,2),(2,2,3),(3,3,1))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3];tag=1)
    add_plane_surface!(model,[1];tag=1)
    return model
end

function gmsh_tilted_plane!()
    gmsh.clear()
    gmsh.model.add("model_entity_evaluation_tilted")
    for (tag,coordinate) in
            pairs(((1.0,2.0,3.0),(2.0,1.0,3.0),(2.0,2.0,2.0)))
        gmsh.model.geo.addPoint(coordinate...,1.0,tag)
    end
    for (tag,first_point,last_point) in ((1,1,2),(2,2,3),(3,3,1))
        gmsh.model.geo.addLine(first_point,last_point,tag)
    end
    gmsh.model.geo.addCurveLoop([1,2,3],1)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()
    return nothing
end

gmsh.initialize(["gmsh","-v","0"])
try
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "model-evaluation differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    tessella=tessella_plane_with_hole()
    gmsh_plane_with_hole!()

    assert_close("Point value",model_value(tessella,0,1,[]),
                 gmsh.model.getValue(0,1,Float64[]))
    model_parametrization_bounds(tessella,0,1)==
        gmsh.model.getParametrizationBounds(0,1) || error("Point bounds differ")
    model_parametrization(tessella,0,1,[1.0,2.0,3.0])==
        gmsh.model.getParametrization(0,1,[1.0,2.0,3.0]) ||
        error("Point parametrization differs")
    model_is_inside(tessella,0,1,[1.0,2.0,3.0])==
        gmsh.model.isInside(0,1,[1.0,2.0,3.0],false) ||
        error("Point containment differs")
    assert_close("Point surface reparametrization",
                 model_reparametrize_on_surface(tessella,0,1,[],1),
                 gmsh.model.reparametrizeOnSurface(0,1,Float64[],1))
    assert_close("off-Plane Point surface reparametrization",
                 model_reparametrize_on_surface(tessella,0,8,[],1),
                 gmsh.model.reparametrizeOnSurface(0,8,Float64[],1))

    line_parameters=[-1.0,0.0,0.5,1.0,2.0]
    assert_close("Line value",model_value(tessella,1,1,line_parameters),
                 gmsh.model.getValue(1,1,line_parameters);atol=1e-10)
    assert_close("Line derivative",
                 model_derivative(tessella,1,1,line_parameters),
                 gmsh.model.getDerivative(1,1,line_parameters);atol=1e-9)
    native_second=model_second_derivative(tessella,1,1,line_parameters)
    reference_second=gmsh.model.getSecondDerivative(1,1,line_parameters)
    all(iszero,native_second) || error("native straight-Line second derivative changed")
    maximum(abs,reference_second)<=1e-4 || error(
        "Gmsh straight-Line numerical second derivative exceeded the bound")
    assert_close("Line curvature",model_curvature(tessella,1,1,line_parameters),
                 gmsh.model.getCurvature(1,1,line_parameters);atol=1e-12)
    native_bounds=model_parametrization_bounds(tessella,1,1)
    reference_bounds=gmsh.model.getParametrizationBounds(1,1)
    native_bounds==reference_bounds || error("Line bounds differ")
    line_coordinates=[1.0,2.0,3.0,3.0,2.0,3.0,5.0,2.0,3.0]
    assert_close("Line parametrization",
                 model_parametrization(tessella,1,1,line_coordinates),
                 gmsh.model.getParametrization(1,1,line_coordinates);atol=1e-10)
    line_inside=[1.0,2.0,3.0,3.0,2.0,3.0,5.0,2.0,3.0,
                 7.0,2.0,3.0,3.0,2.1,3.0]
    model_is_inside(tessella,1,1,line_inside)==
        gmsh.model.isInside(1,1,line_inside,false) ||
        error("Line physical containment differs")
    model_is_inside(tessella,1,1,line_parameters,true)==
        gmsh.model.isInside(1,1,line_parameters,true) ||
        error("Line parametric containment differs")
    line_queries=[3.0,9.0,8.0,7.0,2.0,3.0]
    native_closest,native_closest_parameters=
        model_closest_point(tessella,1,1,line_queries)
    reference_closest,reference_closest_parameters=
        gmsh.model.getClosestPoint(1,1,line_queries)
    assert_close("Line closest points",native_closest,reference_closest;
                 rtol=1e-6,atol=1e-6)
    assert_close("Line closest parameters",native_closest_parameters,
                 reference_closest_parameters;rtol=1e-6,atol=1e-6)
    assert_close("Line surface reparametrization",
                 model_reparametrize_on_surface(
                     tessella,1,1,line_parameters,1,-2),
                 gmsh.model.reparametrizeOnSurface(
                     1,1,line_parameters,1,-2);atol=1e-12)
    off_plane_parameters=[0.0,0.5,1.0]
    assert_close("off-Plane Line surface reparametrization",
                 model_reparametrize_on_surface(
                     tessella,1,8,off_plane_parameters,1),
                 gmsh.model.reparametrizeOnSurface(
                     1,8,off_plane_parameters,1);atol=1e-12)

    plane_parameters=[0.0,0.0,1.0,2.0,-1.0,4.0,7.0,8.0]
    assert_close("Plane value",model_value(tessella,2,1,plane_parameters),
                 gmsh.model.getValue(2,1,plane_parameters);atol=1e-12)
    derivative_parameters=[0.0,0.0,1.0,2.0]
    assert_close("Plane derivative",
                 model_derivative(tessella,2,1,derivative_parameters),
                 gmsh.model.getDerivative(2,1,derivative_parameters);atol=1e-12)
    assert_close("Plane second derivative",
                 model_second_derivative(tessella,2,1,derivative_parameters),
                 gmsh.model.getSecondDerivative(2,1,derivative_parameters);
                 atol=1e-12)
    assert_close("Plane curvature",
                 model_curvature(tessella,2,1,derivative_parameters),
                 gmsh.model.getCurvature(2,1,derivative_parameters);atol=1e-12)
    native_principal=model_principal_curvatures(
        tessella,1,derivative_parameters)
    reference_principal=gmsh.model.getPrincipalCurvatures(
        1,derivative_parameters)
    for index in eachindex(native_principal)
        assert_close("Plane principal result $index",
                     native_principal[index],reference_principal[index];atol=1e-12)
    end
    assert_close("Plane normal",model_normal(tessella,1,derivative_parameters),
                 gmsh.model.getNormal(1,derivative_parameters);atol=1e-12)
    native_bounds=model_parametrization_bounds(tessella,2,1)
    reference_bounds=gmsh.model.getParametrizationBounds(2,1)
    for index in 1:2
        assert_close("Plane bound $index",native_bounds[index],
                     reference_bounds[index];atol=1e-12)
    end
    plane_coordinates=[1.0,2.0,3.0,3.0,4.0,3.0,9.0,9.0,8.0]
    assert_close("Plane parametrization",
                 model_parametrization(tessella,2,1,plane_coordinates),
                 gmsh.model.getParametrization(2,1,plane_coordinates);atol=1e-12)
    plane_inside=[1.0,2.0,3.0,1.1,2.1,3.0,3.0,4.0,3.0,
                  2.5,3.5,3.0,9.0,9.0,3.0,3.0,4.0,4.0]
    model_is_inside(tessella,2,1,plane_inside)==
        gmsh.model.isInside(2,1,plane_inside,false) ||
        error("Plane physical containment differs")
    plane_inside_parameters=[2.0,1.0,4.0,3.0,3.5,2.5,9.0,9.0]
    model_is_inside(tessella,2,1,plane_inside_parameters,true)==
        gmsh.model.isInside(2,1,plane_inside_parameters,true) ||
        error("Plane parametric containment differs")
    plane_queries=[3.0,4.0,8.0,9.0,9.0,8.0]
    native_closest,native_closest_parameters=
        model_closest_point(tessella,2,1,plane_queries)
    reference_closest,reference_closest_parameters=
        gmsh.model.getClosestPoint(2,1,plane_queries)
    assert_close("Plane closest points",native_closest,reference_closest;
                 atol=1e-12)
    assert_close("Plane closest parameters",native_closest_parameters,
                 reference_closest_parameters;atol=1e-12)

    tilted=tessella_tilted_plane()
    gmsh_tilted_plane!()
    tilted_bounds=model_parametrization_bounds(tilted,2,1)
    reference_tilted_bounds=gmsh.model.getParametrizationBounds(2,1)
    for index in 1:2
        assert_close("tilted Plane bound $index",tilted_bounds[index],
                     reference_tilted_bounds[index];atol=1e-12)
    end
    tilted_parameters=reduce(vcat,(
        tilted_bounds[1],tilted_bounds[2],[0.25,-0.5]))
    assert_close("tilted Plane value",
                 model_value(tilted,2,1,tilted_parameters),
                 gmsh.model.getValue(2,1,tilted_parameters);atol=2e-12)
    assert_close("tilted Plane derivative",
                 model_derivative(tilted,2,1,tilted_parameters),
                 gmsh.model.getDerivative(2,1,tilted_parameters);atol=2e-12)
    assert_close("tilted Plane normal",model_normal(tilted,1,tilted_parameters),
                 gmsh.model.getNormal(1,tilted_parameters);atol=2e-12)
    tilted_coordinates=[1.0,2.0,3.0,2.0,1.0,3.0,2.0,2.0,2.0]
    assert_close("tilted Plane parametrization",
                 model_parametrization(tilted,2,1,tilted_coordinates),
                 gmsh.model.getParametrization(2,1,tilted_coordinates);atol=2e-12)
    model_set_tag!(tilted,2,1,10)
    gmsh.model.setTag(2,1,10)
    assert_close("retagged Plane normal",model_normal(tilted,10,[0.0,0.0]),
                 gmsh.model.getNormal(10,[0.0,0.0]);atol=2e-12)

    println("GMSH_PARITY_MODEL_EVALUATION_OK gmsh=",gmsh.GMSH_API_VERSION,
            " point_queries=6 line_query_families=12 ",
            "plane_query_families=12 tilted_query_families=5 ",
            "bounded_divergences=strict_finite_shapes_exact_membership_" *
            "line_second_derivative_noise")
finally
    gmsh.finalize()
end
