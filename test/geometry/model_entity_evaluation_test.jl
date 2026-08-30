using Test
using Tessella
using Tessella.Model: model_value, model_derivative, model_second_derivative,
                      model_curvature, model_principal_curvatures, model_normal,
                      model_parametrization, model_parametrization_bounds,
                      model_is_inside, model_closest_point,
                      model_reparametrize_on_surface, model_set_tag!

function _evaluation_plane_with_hole()
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
    return model
end

@testset "native Point and Line evaluation" begin
    model=GeoModel()
    add_point!(model,1,2,3;tag=1)
    add_point!(model,5,2,3;tag=2)
    add_line!(model,1,2;tag=1)

    @test model_value(model,0,1,[])==[1.0,2.0,3.0]
    @test model_parametrization(model,0,1,[1,2,3])==Float64[]
    @test model_parametrization_bounds(model,0,1)==(Float64[],Float64[])
    @test model_is_inside(model,0,1,[1,2,3])==0

    parameters=[-1.0,0.0,0.5,1.0,2.0]
    @test model_value(model,1,1,parameters)==
          [-3.0,2.0,3.0,1.0,2.0,3.0,3.0,2.0,3.0,
           5.0,2.0,3.0,9.0,2.0,3.0]
    @test model_derivative(model,1,1,parameters)==
          repeat([4.0,0.0,0.0],length(parameters))
    @test model_second_derivative(model,1,1,parameters)==zeros(15)
    @test model_curvature(model,1,1,parameters)==zeros(5)
    lower,upper=model_parametrization_bounds(model,1,1)
    @test lower==[0.0]
    @test upper==[1.0]
    lower[1]=-9;upper[1]=9
    @test model_parametrization_bounds(model,1,1)==([0.0],[1.0])
    @test model_parametrization(
        model,1,1,[1,2,3,3,2,3,7,2,3,3,9,8])==[0.0,0.5,1.5,0.5]
    @test model_is_inside(
        model,1,1,[1,2,3,3,2,3,5,2,3,7,2,3,3,2.1,3])==3
    @test model_is_inside(model,1,1,parameters,true)==3
    closest,closest_parameters=model_closest_point(
        model,1,1,[3,9,8,7,2,3])
    @test closest==[3.0,2.0,3.0,5.0,2.0,3.0]
    @test closest_parameters==[0.5,1.0]
    fill!(closest,0);fill!(closest_parameters,0)
    @test model_closest_point(model,1,1,[3,9,8])==
          ([3.0,2.0,3.0],[0.5])
end

@testset "native Plane evaluation" begin
    model=_evaluation_plane_with_hole()
    @test model_parametrization_bounds(model,2,1)==([2.0,1.0],[6.0,5.0])
    parameters=[0.0,0.0,1.0,2.0,-1.0,4.0,7.0,8.0]
    @test model_value(model,2,1,parameters)==
          [0.0,0.0,3.0,2.0,1.0,3.0,4.0,-1.0,3.0,8.0,7.0,3.0]
    @test model_derivative(model,2,1,[0.0,0.0,1.0,2.0])==
          repeat([0.0,1.0,0.0,1.0,0.0,0.0],2)
    @test model_second_derivative(model,2,1,[0.0,0.0,1.0,2.0])==zeros(18)
    @test model_curvature(model,2,1,[0.0,0.0,1.0,2.0])==zeros(2)
    maximum_curvature,minimum_curvature,maximum_direction,minimum_direction=
        model_principal_curvatures(model,1,[0.0,0.0,1.0,2.0])
    @test maximum_curvature==minimum_curvature==zeros(2)
    @test maximum_direction==repeat([0.0,1.0,0.0],2)
    @test minimum_direction==repeat([1.0,0.0,0.0],2)
    @test model_normal(model,1,[0.0,0.0,1.0,2.0])==
          repeat([0.0,0.0,1.0],2)
    @test model_parametrization(
        model,2,1,[1,2,3,3,4,3,9,9,8])==[2.0,1.0,4.0,3.0,9.0,9.0]
    @test model_is_inside(
        model,2,1,[1,2,3,1.1,2.1,3,3,4,3,
                   2.5,3.5,3,9,9,3,3,4,4])==2
    @test model_is_inside(
        model,2,1,[2,1,4,3,3.5,2.5,9,9],true)==3
    @test model_closest_point(model,2,1,[3,4,8,9,9,8])==
          ([3.0,4.0,3.0,9.0,9.0,3.0],[4.0,3.0,9.0,9.0])
    @test model_reparametrize_on_surface(model,0,1,[],1)==[2.0,1.0]
    @test model_reparametrize_on_surface(
        model,1,1,[-1.0,0.0,0.5,1.0,2.0],1,-2)==
        [2.0,-3.0,2.0,1.0,2.0,3.0,2.0,5.0,2.0,9.0]

    add_point!(model,0,0,4;tag=8)
    add_point!(model,1,0,4;tag=9)
    add_line!(model,8,9;tag=8)
    @test model_reparametrize_on_surface(model,1,8,[0.0,0.5,1.0],1)==
          [0.0,0.0,0.0,0.5,0.0,1.0]

    model_set_tag!(model,2,1,10)
    @test model_normal(model,10,[4.0,3.0])==[0.0,0.0,1.0]
    @test model_reparametrize_on_surface(model,1,1,[0.5],10)==[2.0,3.0]
    @test_throws ArgumentError model_value(model,2,1,[4.0,3.0])

    tilted=GeoModel()
    for (tag,point) in pairs(((1.0,2.0,3.0),(2.0,1.0,3.0),(2.0,2.0,2.0)))
        add_point!(tilted,point...;tag=tag)
    end
    for (tag,first_point,last_point) in ((1,1,2),(2,2,3),(3,3,1))
        add_line!(tilted,first_point,last_point;tag=tag)
    end
    add_curve_loop!(tilted,[1,2,3];tag=1)
    add_plane_surface!(tilted,[1];tag=1)
    coordinates=[1.0,2.0,3.0,2.0,1.0,3.0,2.0,2.0,2.0]
    tilted_parameters=model_parametrization(tilted,2,1,coordinates)
    @test all(isapprox.(model_value(tilted,2,1,tilted_parameters),coordinates;
                       rtol=8eps(Float64),atol=8eps(Float64)))
    @test all(isapprox.(model_normal(tilted,1,[0.0,0.0]),
                       fill(inv(sqrt(3.0)),3);rtol=4eps(Float64)))
end

@testset "native evaluation validation" begin
    model=_evaluation_plane_with_hole()
    invalid_calls=(
        ()->model_value(model,3,1,[]),
        ()->model_value(model,0,1,[0.0]),
        ()->model_value(model,2,1,[0.0]),
        ()->model_derivative(model,0,1,[]),
        ()->model_derivative(model,1,99,[0.0]),
        ()->model_second_derivative(model,2,1,[0.0]),
        ()->model_curvature(model,1,1,[NaN]),
        ()->model_principal_curvatures(model,99,[0.0,0.0]),
        ()->model_normal(model,1,[0.0]),
        ()->model_parametrization(model,1,1,[0.0,0.0]),
        ()->model_parametrization_bounds(model,3,1),
        ()->model_is_inside(model,2,1,[0.0,0.0],false),
        ()->model_is_inside(model,2,1,[],1),
        ()->model_closest_point(model,0,1,[0.0,0.0,0.0]),
        ()->model_reparametrize_on_surface(model,2,1,[],1),
        ()->model_reparametrize_on_surface(model,0,1,[0.0],1),
        ()->model_reparametrize_on_surface(model,1,1,[0.5],99),
        ()->model_reparametrize_on_surface(model,1,1,[0.5],1,true),
        ()->model_reparametrize_on_surface(
            model,1,1,[0.5],1,big(typemax(Int32))+1),
        ()->model_reparametrize_on_surface(model,1,1,[NaN],1),
    )
    for call in invalid_calls
        @test_throws ArgumentError call()
    end

    degenerate=GeoModel()
    add_point!(degenerate,0,0,0;tag=1)
    add_point!(degenerate,0,0,0;tag=2)
    add_line!(degenerate,1,2;tag=1)
    @test_throws ArgumentError model_value(degenerate,1,1,[0.5])

    wide=GeoModel()
    add_point!(wide,-floatmax(Float64),0,0;tag=1)
    add_point!(wide,floatmax(Float64),0,0;tag=2)
    add_line!(wide,1,2;tag=1)
    @test model_value(wide,1,1,[0.5])==[0.0,0.0,0.0]
    @test_throws ArgumentError model_derivative(wide,1,1,[0.5])

    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end
