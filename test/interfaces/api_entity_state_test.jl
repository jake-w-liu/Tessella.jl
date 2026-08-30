using Test
using Tessella

const _ENTITY_STATE_API=Tessella.API

@testset "synchronized API entity state" begin
    _ENTITY_STATE_API.finalize()
    try
        _ENTITY_STATE_API.initialize()
        for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                          (3,1.0,1.0),(4,0.0,1.0))
            _ENTITY_STATE_API.model.add_point(x,y,0;tag=tag)
        end
        for (tag,first_point,last_point) in
                ((1,1,2),(2,2,3),(3,3,4),(4,4,1))
            _ENTITY_STATE_API.model.add_line(first_point,last_point;tag=tag)
        end
        _ENTITY_STATE_API.model.add_curve_loop([1,2,3,4];tag=1)
        _ENTITY_STATE_API.model.add_plane_surface([1];tag=1)
        _ENTITY_STATE_API.mesh.generate(2)
        cached=_ENTITY_STATE_API.LAST_MESH[]

        @test _ENTITY_STATE_API.model.get_visibility(2,1)==1
        @test _ENTITY_STATE_API.model.get_color(2,1)==(0,0,255,0)
        @test _ENTITY_STATE_API.model.set_visibility([(2,1)],0,true)===nothing
        @test _ENTITY_STATE_API.model.get_visibility(0,1)==0
        @test _ENTITY_STATE_API.model.get_visibility(1,1)==0
        @test _ENTITY_STATE_API.model.get_visibility(2,1)==0
        @test _ENTITY_STATE_API.model.set_color(
            [(2,1)],10,20,30,40,true)===nothing
        @test _ENTITY_STATE_API.model.get_color(0,1)==(10,20,30,40)
        @test _ENTITY_STATE_API.model.get_color(1,1)==(10,20,30,40)
        @test _ENTITY_STATE_API.model.get_color(2,1)==(10,20,30,40)
        @test _ENTITY_STATE_API.LAST_MESH[]===cached

        @test _ENTITY_STATE_API.model.set_attribute("b",["2","x"])===nothing
        @test _ENTITY_STATE_API.model.set_attribute("a",["1"])===nothing
        @test _ENTITY_STATE_API.model.get_attribute_names()==["a","b"]
        values=_ENTITY_STATE_API.model.get_attribute("b")
        values[1]="changed"
        @test _ENTITY_STATE_API.model.get_attribute("b")==["2","x"]
        @test _ENTITY_STATE_API.model.get_attribute("missing")==String[]
        @test _ENTITY_STATE_API.model.remove_attribute("a")===nothing
        @test _ENTITY_STATE_API.model.get_attribute_names()==["b"]
        @test _ENTITY_STATE_API.LAST_MESH[]===cached

        @test_throws ArgumentError _ENTITY_STATE_API.model.set_coordinates(
            99,0,0,0)
        @test _ENTITY_STATE_API.LAST_MESH[]===cached
        @test _ENTITY_STATE_API.model.set_coordinates(1,-0.25,0,0)===nothing
        @test _ENTITY_STATE_API.LAST_MESH[]===nothing
        @test _ENTITY_STATE_API.model.get_value(0,1,[])==[-0.25,0.0,0.0]
        _ENTITY_STATE_API.mesh.generate(2)
        cached=_ENTITY_STATE_API.LAST_MESH[]

        _ENTITY_STATE_API.model.set_tag(2,1,10)
        @test _ENTITY_STATE_API.model.get_visibility(2,10)==0
        @test _ENTITY_STATE_API.model.get_color(2,10)==(10,20,30,40)
        @test_throws ArgumentError _ENTITY_STATE_API.model.get_visibility(2,1)
        @test _ENTITY_STATE_API.LAST_MESH[]===nothing

        for function_name in (
            :set_visibility,:get_visibility,:set_color,:get_color,
            :set_coordinates,:set_attribute,:get_attribute,
            :get_attribute_names,:remove_attribute)
            @test (@doc getproperty(_ENTITY_STATE_API.model,function_name))!==nothing
        end
        @test isempty(Docs.undocumented_names(Tessella.API;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.API;recursive=true))
    finally
        _ENTITY_STATE_API.finalize()
    end

    _ENTITY_STATE_API.initialize()
    try
        @test _ENTITY_STATE_API.model.get_attribute_names()==String[]
    finally
        _ENTITY_STATE_API.finalize()
    end
end
