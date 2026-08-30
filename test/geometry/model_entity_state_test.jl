using Test
using Tessella
using Tessella.Model: set_entity_visibility!, model_entity_visibility,
                      set_entity_color!, model_entity_color,
                      set_point_coordinates!, set_model_attribute!,
                      model_attribute, model_attribute_names,
                      remove_model_attribute!, model_set_tag!, remove_entities!,
                      model_value, model_bounding_box

function _entity_state_square()
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

@testset "native entity presentation state" begin
    model=_entity_state_square()
    @test [model_entity_visibility(model,dimension,1)
           for dimension in 0:2]==[1,1,1]
    @test [model_entity_color(model,dimension,1)
           for dimension in 0:2]==fill((0,0,255,0),3)

    @test set_entity_visibility!(model,[(2,1)],0)===nothing
    @test model_entity_visibility(model,2,1)==0
    @test model_entity_visibility(model,1,1)==1
    @test set_entity_visibility!(model,[(2,1)],-2,true)===nothing
    @test all(model_entity_visibility(model,dimension,tag)==-2
              for dimension in 0:2
              for tag in keys(Tessella.Model._model_entity_dictionary(
                  model,dimension)))
    @test set_entity_visibility!(model,[(2,99),(2,99)],7,true)===nothing

    @test set_entity_color!(model,[(2,1)],1,2,3,4)===nothing
    @test model_entity_color(model,2,1)==(1,2,3,4)
    @test model_entity_color(model,1,1)==(0,0,255,0)
    @test set_entity_color!(model,[(2,1)],10,20,30,40,true)===nothing
    @test all(model_entity_color(model,dimension,tag)==(10,20,30,40)
              for dimension in 0:2
              for tag in keys(Tessella.Model._model_entity_dictionary(
                  model,dimension)))
    @test_throws ArgumentError set_entity_color!(
        model,[(2,1),(1,1)],10,20,300,40,true)
    @test model_entity_color(model,2,1)==(10,20,30,40)

    model_set_tag!(model,2,1,10)
    @test model_entity_visibility(model,2,10)==-2
    @test model_entity_color(model,2,10)==(10,20,30,40)
    @test_throws ArgumentError model_entity_visibility(model,2,1)
    model.entity_visibility[(2,99)]=0
    @test_throws ArgumentError model_set_tag!(model,2,10,99)
    @test haskey(model.surfaces,10)
    @test model_entity_visibility(model,2,10)==-2
    delete!(model.entity_visibility,(2,99))

    @test remove_entities!(model,[(2,10)])==1
    add_plane_surface!(model,[1];tag=10)
    @test model_entity_visibility(model,2,10)==1
    @test model_entity_color(model,2,10)==(0,0,255,0)
end

@testset "native Point coordinates and model attributes" begin
    model=_entity_state_square()
    @test set_point_coordinates!(model,1,-1,0,0)===nothing
    @test model.points[1]==(-1.0,0.0,0.0)
    @test model.point_size[1]==1.0
    @test model_value(model,1,1,[0.0,0.5,1.0])==
          [-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0]
    @test model_bounding_box(model,2,1)==(-1.0,0.0,0.0,1.0,1.0,0.0)
    @test_throws ArgumentError set_point_coordinates!(model,1,NaN,0,0)
    @test_throws ArgumentError set_point_coordinates!(model,99,0,0,0)
    @test model.points[1]==(-1.0,0.0,0.0)

    values=["2","x"]
    @test set_model_attribute!(model,"b",values)===nothing
    values[1]="changed"
    set_model_attribute!(model,"a",["1"])
    set_model_attribute!(model,"",String[])
    @test model_attribute_names(model)==["","a","b"]
    @test model_attribute(model,"b")==["2","x"]
    detached=model_attribute(model,"b")
    detached[1]="changed"
    @test model_attribute(model,"b")==["2","x"]
    @test model_attribute(model,"missing")==String[]
    @test_throws ArgumentError set_model_attribute!(model,"bad\0name",["x"])
    @test_throws ArgumentError set_model_attribute!(model,"bad",["x\0y"])
    @test_throws ArgumentError set_model_attribute!(model,"bad",[1])
    @test model_attribute_names(model)==["","a","b"]
    @test remove_model_attribute!(model,"a")
    @test !remove_model_attribute!(model,"a")
    @test model_attribute_names(model)==["","b"]

    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end
