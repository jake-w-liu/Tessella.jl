using Test
using Tessella
using Tessella.Model: model_bounding_box, model_entities,
                      model_entities_in_bounding_box, model_set_tag!,
                      remove_entities!, translate_volume!

_spatial_bounds_approx(first,second;rtol=4eps(Float64))=
    all(isapprox(first[index],second[index];rtol=rtol,atol=0.0)
        for index in 1:6)

function _spatial_tetrahedron()
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

@testset "explicit model bounding boxes" begin
    model=_spatial_tetrahedron()
    expected=Dict(
        (0,1)=>(0.0,0.0,0.0,0.0,0.0,0.0),
        (0,2)=>(2.0,0.0,0.0,2.0,0.0,0.0),
        (0,3)=>(0.0,3.0,0.0,0.0,3.0,0.0),
        (0,4)=>(0.0,0.0,4.0,0.0,0.0,4.0),
        (1,1)=>(0.0,0.0,0.0,2.0,0.0,0.0),
        (1,2)=>(0.0,0.0,0.0,2.0,3.0,0.0),
        (1,3)=>(0.0,0.0,0.0,0.0,3.0,0.0),
        (1,4)=>(0.0,0.0,0.0,0.0,0.0,4.0),
        (1,5)=>(0.0,0.0,0.0,2.0,0.0,4.0),
        (1,6)=>(0.0,0.0,0.0,0.0,3.0,4.0),
        (2,1)=>(0.0,0.0,0.0,2.0,3.0,0.0),
        (2,2)=>(0.0,0.0,0.0,2.0,0.0,4.0),
        (2,3)=>(0.0,0.0,0.0,2.0,3.0,4.0),
        (2,4)=>(0.0,0.0,0.0,0.0,3.0,4.0),
        (3,1)=>(0.0,0.0,0.0,2.0,3.0,4.0),
    )
    @test Set(keys(expected))==Set(model_entities(model))
    for entity in model_entities(model)
        @test model_bounding_box(model,entity...)==expected[entity]
    end
    @test model_bounding_box(model,-1,-1)==(0.0,0.0,0.0,2.0,3.0,4.0)

    add_point!(model,10,-2,7;tag=10)
    @test model_bounding_box(model,-1,-1)==(0.0,-2.0,0.0,10.0,3.0,7.0)
    model_set_tag!(model,0,10,20)
    @test model_bounding_box(model,0,20)==(10.0,-2.0,7.0,10.0,-2.0,7.0)
    @test_throws ArgumentError model_bounding_box(model,0,10)
end

@testset "bounding-box containment selection" begin
    model=_spatial_tetrahedron()
    @test model_entities_in_bounding_box(model,0,0,0,2,3,4)==
          model_entities(model)
    @test model_entities_in_bounding_box(model,0,0,0,2,0,0)==
          [(0,1),(0,2),(1,1)]
    @test model_entities_in_bounding_box(model,0,0,0,2,3,0,2)==
          [(2,1)]
    @test model_entities_in_bounding_box(model,0,0,0,2,3,4,2)==
          [(2,1),(2,2),(2,3),(2,4)]
    @test isempty(model_entities_in_bounding_box(model,1,1,1,0,0,0))
    detached=model_entities_in_bounding_box(model,0,0,0,2,3,4)
    empty!(detached)
    @test length(model_entities_in_bounding_box(model,0,0,0,2,3,4))==15

    embedded=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,1.0,1.0),(4,0.0,1.0),(9,3.0,3.0))
        add_point!(embedded,x,y,0;tag=tag)
    end
    for (tag,a,b) in ((1,1,2),(2,2,3),(3,3,4),(4,4,1))
        add_line!(embedded,a,b;tag=tag)
    end
    add_curve_loop!(embedded,[1,2,3,4];tag=1)
    add_plane_surface!(embedded,[1];tag=1)
    embed!(embedded,0,[9],2,1)
    @test model_bounding_box(embedded,2,1)==(0.0,0.0,0.0,1.0,1.0,0.0)
    @test model_bounding_box(embedded,-1,-1)==(0.0,0.0,0.0,3.0,3.0,0.0)
end

@testset "analytical primitive bounding boxes" begin
    model=GeoModel()
    add_box!(model,-2,1,3,4,5,6;tag=1)
    add_cylinder!(model,10,20,30,2,3,6,4;tag=2)
    add_sphere!(model,-10,-20,-30,5;tag=3)
    add_cone!(model,1,2,3,-2,4,5,6,2;tag=4)
    @test model_bounding_box(model,3,1)==(-2.0,1.0,3.0,2.0,6.0,9.0)
    @test model_bounding_box(model,3,3)==
          (-15.0,-25.0,-35.0,-5.0,-15.0,-25.0)

    setprecision(BigFloat,256) do
        cylinder_expected=Float64.((
            big"10"-big"4"*sqrt(big"45")/7,
            big"20"-big"4"*sqrt(big"40")/7,
            big"30"-big"4"*sqrt(big"13")/7,
            big"12"+big"4"*sqrt(big"45")/7,
            big"23"+big"4"*sqrt(big"40")/7,
            big"36"+big"4"*sqrt(big"13")/7,
        ))
        @test _spatial_bounds_approx(
            model_bounding_box(model,3,2),cylinder_expected)
        cone_expected=Float64.((
            min(big"1"-big"6"*sqrt(big"41")/sqrt(big"45"),
                big"-1"-big"2"*sqrt(big"41")/sqrt(big"45")),
            min(big"2"-big"6"*sqrt(big"29")/sqrt(big"45"),
                big"6"-big"2"*sqrt(big"29")/sqrt(big"45")),
            min(big"3"-big"6"*sqrt(big"20")/sqrt(big"45"),
                big"8"-big"2"*sqrt(big"20")/sqrt(big"45")),
            max(big"1"+big"6"*sqrt(big"41")/sqrt(big"45"),
                big"-1"+big"2"*sqrt(big"41")/sqrt(big"45")),
            max(big"2"+big"6"*sqrt(big"29")/sqrt(big"45"),
                big"6"+big"2"*sqrt(big"29")/sqrt(big"45")),
            max(big"3"+big"6"*sqrt(big"20")/sqrt(big"45"),
                big"8"+big"2"*sqrt(big"20")/sqrt(big"45")),
        ))
        @test _spatial_bounds_approx(
            model_bounding_box(model,3,4),cone_expected)
    end

    @test model_entities_in_bounding_box(
        model,-2,1,3,2,6,9,3)==[(3,1)]
    @test _spatial_bounds_approx(
        model_bounding_box(model,-1,-1),
        (-15.0,-25.0,-35.0,model_bounding_box(model,3,2)[4:6]...))
    translate_volume!(model,3,(100,200,300))
    @test model_bounding_box(model,3,3)==
          (85.0,175.0,265.0,95.0,185.0,275.0)
end

@testset "Boolean snapshot bounding boxes" begin
    model=GeoModel()
    add_box!(model,0,0,0,2,1,1;tag=1)
    add_box!(model,0,0,0,1,1,1;tag=2)
    boolean_volumes!(model,:difference,1,2;tag=3)
    expected=(1.0,0.0,0.0,2.0,1.0,1.0)
    @test model_bounding_box(model,3,3)==expected
    translate_volume!(model,1,(100,0,0))
    @test model_bounding_box(model,3,3)==expected
    @test remove_entities!(model,[(3,1),(3,2)])==2
    @test model_bounding_box(model,3,3)==expected
end

@testset "spatial query validation" begin
    empty_model=GeoModel()
    @test_throws ArgumentError model_bounding_box(empty_model,-1,-1)
    model=_spatial_tetrahedron()
    for call in (
        ()->model_bounding_box(model,-1,1),
        ()->model_bounding_box(model,0,-1),
        ()->model_bounding_box(model,4,1),
        ()->model_bounding_box(model,true,1),
        ()->model_bounding_box(model,0,true),
        ()->model_bounding_box(model,0,99),
        ()->model_entities_in_bounding_box(model,NaN,0,0,1,1,1),
        ()->model_entities_in_bounding_box(model,0,0,0,Inf,1,1),
        ()->model_entities_in_bounding_box(model,false,0,0,1,1,1),
        ()->model_entities_in_bounding_box(model,"0",0,0,1,1,1),
        ()->model_entities_in_bounding_box(model,0,0,0,1,1,1,4),
        ()->model_entities_in_bounding_box(model,0,0,0,1,1,1,true),
    )
        @test_throws ArgumentError call()
    end

    corrupt=GeoModel()
    add_point!(corrupt,0,0,0;tag=1)
    corrupt.points[1]=(NaN,0.0,0.0)
    @test_throws ArgumentError model_bounding_box(corrupt,0,1)
    overflow=GeoModel()
    add_box!(overflow,floatmax(Float64),0,0,floatmax(Float64),1,1;tag=1)
    @test_throws ArgumentError model_bounding_box(overflow,3,1)
    multiply_encoded=GeoModel()
    add_box!(multiply_encoded,0,0,0,1,1,1;tag=1)
    multiply_encoded.spheres[1]=(center=(0.0,0.0,0.0),radius=1.0)
    before=deepcopy(multiply_encoded)
    @test_throws ErrorException model_bounding_box(multiply_encoded,3,1)
    @test multiply_encoded.box_extents==before.box_extents
    @test multiply_encoded.spheres==before.spheres
    explicit_and_primitive=_spatial_tetrahedron()
    explicit_and_primitive.box_extents[1]=(0.0,0.0,0.0,1.0,1.0,1.0)
    @test_throws ErrorException model_bounding_box(explicit_and_primitive,3,1)
    explicit_and_primitive.points[1]=(NaN,0.0,0.0)
    @test isempty(model_entities_in_bounding_box(
        explicit_and_primitive,1,1,1,0,0,0))

    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end
