using Test
using Tessella

const _EVALUATION_API=Tessella.API

@testset "synchronized API entity evaluation" begin
    _EVALUATION_API.finalize()
    try
        _EVALUATION_API.initialize()
        for (tag,x,y,z) in ((1,1.0,2.0,0.0),(2,5.0,2.0,0.0),
                            (3,5.0,6.0,0.0),(4,1.0,6.0,0.0))
            _EVALUATION_API.model.add_point(x,y,z;tag=tag)
        end
        for (tag,first_point,last_point) in
                ((1,1,2),(2,2,3),(3,3,4),(4,4,1))
            _EVALUATION_API.model.add_line(first_point,last_point;tag=tag)
        end
        _EVALUATION_API.model.add_curve_loop([1,2,3,4];tag=1)
        _EVALUATION_API.model.add_plane_surface([1];tag=1)
        _EVALUATION_API.mesh.generate(2)
        cached=_EVALUATION_API.LAST_MESH[]

        @test _EVALUATION_API.model.get_value(0,1,[])==[1.0,2.0,0.0]
        @test _EVALUATION_API.model.get_value(1,1,[0.0,0.5,1.0])==
              [1.0,2.0,0.0,3.0,2.0,0.0,5.0,2.0,0.0]
        @test _EVALUATION_API.model.get_derivative(1,1,[0.5])==
              [4.0,0.0,0.0]
        @test _EVALUATION_API.model.get_second_derivative(2,1,[4.0,3.0])==
              zeros(9)
        @test _EVALUATION_API.model.get_curvature(2,1,[4.0,3.0])==[0.0]
        maximum,minimum,maximum_direction,minimum_direction=
            _EVALUATION_API.model.get_principal_curvatures(1,[4.0,3.0])
        @test maximum==minimum==[0.0]
        @test maximum_direction==[0.0,1.0,0.0]
        @test minimum_direction==[1.0,0.0,0.0]
        @test _EVALUATION_API.model.get_normal(1,[4.0,3.0])==
              [0.0,0.0,1.0]
        @test _EVALUATION_API.model.get_parametrization(
            2,1,[3.0,4.0,8.0])==[4.0,3.0]
        lower,upper=_EVALUATION_API.model.get_parametrization_bounds(2,1)
        @test lower==[2.0,1.0]
        @test upper==[6.0,5.0]
        fill!(lower,0);fill!(upper,0)
        @test _EVALUATION_API.model.get_parametrization_bounds(2,1)==
              ([2.0,1.0],[6.0,5.0])
        @test _EVALUATION_API.model.is_inside(
            2,1,[3.0,4.0,0.0,9.0,9.0,0.0])==1
        @test _EVALUATION_API.model.is_inside(
            2,1,[4.0,3.0,9.0,9.0],true)==1
        @test _EVALUATION_API.model.get_closest_point(
            2,1,[3.0,4.0,8.0])==([3.0,4.0,0.0],[4.0,3.0])
        @test _EVALUATION_API.LAST_MESH[]===cached

        for call in (
            ()->_EVALUATION_API.model.get_value(2,99,[0.0,0.0]),
            ()->_EVALUATION_API.model.get_derivative(0,1,[]),
            ()->_EVALUATION_API.model.get_normal(99,[0.0,0.0]),
            ()->_EVALUATION_API.model.is_inside(2,1,[NaN,0.0,0.0]),
            ()->_EVALUATION_API.model.get_closest_point(0,1,[1.0,2.0,3.0]),
        )
            @test_throws ArgumentError call()
            @test _EVALUATION_API.LAST_MESH[]===cached
        end

        @test (@doc _EVALUATION_API.model.get_value)!==nothing
        @test (@doc _EVALUATION_API.model.get_derivative)!==nothing
        @test (@doc _EVALUATION_API.model.get_second_derivative)!==nothing
        @test (@doc _EVALUATION_API.model.get_curvature)!==nothing
        @test (@doc _EVALUATION_API.model.get_principal_curvatures)!==nothing
        @test (@doc _EVALUATION_API.model.get_normal)!==nothing
        @test (@doc _EVALUATION_API.model.get_parametrization)!==nothing
        @test (@doc _EVALUATION_API.model.get_parametrization_bounds)!==nothing
        @test (@doc _EVALUATION_API.model.is_inside)!==nothing
        @test (@doc _EVALUATION_API.model.get_closest_point)!==nothing
        @test isempty(Docs.undocumented_names(Tessella.API;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.API;recursive=true))
    finally
        _EVALUATION_API.finalize()
    end
end
