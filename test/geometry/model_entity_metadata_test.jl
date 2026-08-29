using Test
using Tessella
using Tessella.Model: model_entities, model_entity_type, model_parent,
                      model_number_of_partitions, model_partitions,
                      model_set_tag!, remove_entities!

function _metadata_tetrahedron()
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

@testset "native entity type metadata" begin
    model=_metadata_tetrahedron()
    expected=("Point","Line","Plane","Volume")
    @test model_number_of_partitions(model)==0
    for entity in model_entities(model)
        @test model_entity_type(model,entity...)==expected[entity[1]+1]
        @test model_parent(model,entity...)==(-1,-1)
        partitions=model_partitions(model,entity...)
        @test isempty(partitions)
        push!(partitions,99)
        @test isempty(model_partitions(model,entity...))
    end

    primitives=GeoModel()
    add_box!(primitives,0,0,0,1,1,1;tag=1)
    add_cylinder!(primitives,0,0,0,1,2,3,0.5;tag=2)
    add_sphere!(primitives,0,0,0,1;tag=3)
    add_cone!(primitives,0,0,0,1,2,3,1,0.25;tag=4)
    for tag in 1:4
        @test model_entity_type(primitives,3,tag)=="Volume"
    end

    model_set_tag!(model,3,1,10)
    @test model_entity_type(model,3,10)=="Volume"
    @test_throws ArgumentError model_entity_type(model,3,1)
    @test remove_entities!(model,[(3,10)])==1
    @test_throws ArgumentError model_entity_type(model,3,10)
end

@testset "entity metadata validation" begin
    model=_metadata_tetrahedron()
    for call in (
        ()->model_entity_type(model,-1,1),
        ()->model_entity_type(model,4,1),
        ()->model_entity_type(model,true,1),
        ()->model_entity_type(model,0,0),
        ()->model_entity_type(model,0,-1),
        ()->model_entity_type(model,0,true),
        ()->model_entity_type(model,0,99),
        ()->model_parent(model,2,99),
        ()->model_partitions(model,3,99),
    )
        @test_throws ArgumentError call()
    end
    @test model_number_of_partitions(GeoModel())==0
    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end
