using Test
using Tessella
using Tessella.Model: model_bounding_box, model_entities,
                      model_entities_in_bounding_box, model_set_tag!,
                      remove_entities!, translate_volume!, model_distance,
                      add_discrete_entity!, add_discrete_nodes!,
                      add_discrete_elements!

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

@testset "near-axis cylinder and cone bounds retain radial extent" begin
    for dominant in 1:3, direction in (-1.0,1.0), scale in (1.0,1.0e5)
        axis=ntuple(i->scale*(i==dominant ? 3direction :
                              i==mod1(dominant+1,3) ? 3.0e-9 : 0.0),3)
        for (kind,r1,r2) in ((:cylinder,1.0,1.0),(:cone,1.0,0.5),
                              (:cone,0.5,1.0))
            model=GeoModel()
            tag=kind===:cylinder ? add_cylinder!(model,0,0,0,axis...,r1) :
                                  add_cone!(model,0,0,0,axis...,r1,r2)
            expected=setprecision(BigFloat,256) do
                a=BigFloat.(axis)
                magnitude=sqrt(sum(x->x^2,a))
                radial=ntuple(i->sqrt(sum(a[j]^2 for j in 1:3 if j!=i))/
                                  magnitude,3)
                lower=ntuple(i->min(-r1*radial[i],a[i]-r2*radial[i]),3)
                upper=ntuple(i->max(r1*radial[i],a[i]+r2*radial[i]),3)
                Float64.((lower...,upper...))
            end
            actual=model_bounding_box(model,3,tag)
            @test _spatial_bounds_approx(actual,expected)
            # The base circle extends across the zero coordinate on the
            # nearly aligned axis; a rounded-to-zero radial bound is wrong.
            @test direction>0 ? actual[dominant]<0 : actual[dominant+3]>0
        end
    end
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
    # Corners at floatmax+floatmax overflow to Inf — a materialized B-rep cannot
    # carry non-finite points, so construction itself fails (Gmsh's addBox
    # rejects degenerate boxes the same way).
    @test_throws ArgumentError add_box!(
        overflow,floatmax(Float64),0,0,floatmax(Float64),1,1;tag=1)
    @test isempty(overflow.volumes)
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

function _spatial_unit_square!(model,z,tbase)
    for (tag,x,y) in ((tbase+1,0.0,0.0),(tbase+2,1.0,0.0),
                      (tbase+3,1.0,1.0),(tbase+4,0.0,1.0))
        add_point!(model,x,y,z;tag=tag)
    end
    for (tag,a,b) in ((tbase+1,tbase+1,tbase+2),(tbase+2,tbase+2,tbase+3),
                      (tbase+3,tbase+3,tbase+4),(tbase+4,tbase+4,tbase+1))
        add_line!(model,a,b;tag=tag)
    end
    add_curve_loop!(model,[tbase+1,tbase+2,tbase+3,tbase+4];tag=tbase+1)
    add_plane_surface!(model,[tbase+1];tag=tbase+1)
    return tbase+1
end

# Assert the getDistance contract: the returned pair attains `d` and is
# symmetric under argument swap.
function _spatial_check_distance(model,d1,t1,d2,t2,expected;atol=1e-9)
    d,pa,pb=model_distance(model,d1,t1,d2,t2)
    @test d≈expected atol=atol*max(1.0,abs(expected))
    @test hypot(pa[1]-pb[1],pa[2]-pb[2],pa[3]-pb[3])≈d atol=1e-7*max(1.0,d)
    ds,pa_s,pb_s=model_distance(model,d2,t2,d1,t1)
    @test ds≈d atol=0.0 rtol=0.0
    @test pa_s==pb && pb_s==pa
    return d,pa,pb
end

@testset "entity distance — elementary pairs" begin
    # Oracles: gmsh 4.15.2 occ.getDistance on identical OCC geometry.
    model=GeoModel()
    add_point!(model,0,0,0;tag=1); add_point!(model,3,4,0;tag=2)
    d,pa,pb=_spatial_check_distance(model,0,1,0,2,5.0)
    @test pa==(0.0,0.0,0.0) && pb==(3.0,4.0,0.0)

    model=GeoModel()
    add_point!(model,0,0,0;tag=1); add_point!(model,2,0,0;tag=2)
    add_line!(model,1,2;tag=1); add_point!(model,1,4,0;tag=3)
    _,pa,pb=_spatial_check_distance(model,0,3,1,1,4.0)
    @test pb==(1.0,0.0,0.0)

    # Point on the query curve: distance 0 with the shared point.
    _,pa,pb=_spatial_check_distance(model,0,1,1,1,0.0)
    @test pa==pb==(0.0,0.0,0.0)

    # Quarter circle about the origin vs a point at (0,0,3): sqrt(1+9).
    model=GeoModel()
    add_point!(model,1,0,0;tag=1); add_point!(model,0,0,0;tag=2)
    add_point!(model,0,1,0;tag=3); add_point!(model,0,0,3;tag=4)
    add_circle_arc!(model,1,2,3;tag=1)
    _spatial_check_distance(model,0,4,1,1,sqrt(10.0))

    model=GeoModel()
    square=_spatial_unit_square!(model,0.0,0)
    add_point!(model,0.5,0.5,2;tag=9)
    _,pa,pb=_spatial_check_distance(model,0,9,2,square,2.0)
    @test all(pb.≈(0.5,0.5,0.0))

    # Projection outside the trim: answer is the boundary-curve distance.
    add_point!(model,3.0,0.5,0.0;tag=10)
    _spatial_check_distance(model,0,10,2,square,2.0)

    model=GeoModel()
    add_point!(model,-1,0,0;tag=1); add_point!(model,1,0,0;tag=2)
    add_point!(model,0,-1,0;tag=3); add_point!(model,0,1,0;tag=4)
    add_line!(model,1,2;tag=1); add_line!(model,3,4;tag=2)
    _spatial_check_distance(model,1,1,1,2,0.0)

    model=GeoModel()
    add_point!(model,0,0,0;tag=1); add_point!(model,1,0,0;tag=2)
    add_point!(model,0,1,3;tag=3); add_point!(model,0,1,4;tag=4)
    add_line!(model,1,2;tag=1); add_line!(model,3,4;tag=2)
    _spatial_check_distance(model,1,1,1,2,sqrt(10.0))

    model=GeoModel()
    square=_spatial_unit_square!(model,0.0,0)
    add_point!(model,0.5,0.5,2;tag=9); add_point!(model,0.5,0.6,2;tag=10)
    add_line!(model,9,10;tag=9)
    _spatial_check_distance(model,1,9,2,square,2.0)

    model=GeoModel()
    s1=_spatial_unit_square!(model,0.0,0)
    s2=_spatial_unit_square!(model,2.0,10)
    _spatial_check_distance(model,2,s1,2,s2,2.0)
    # Shared boundary edge: curve on s1 to s1 itself is 0.
    _spatial_check_distance(model,1,1,2,s1,0.0)
end

@testset "entity distance — solids and curved primitives" begin
    model=GeoModel()
    add_box!(model,0,0,0,1,1,1;tag=1)
    add_point!(model,3,0,0;tag=50)
    _,pa,pb=_spatial_check_distance(model,0,50,3,1,2.0)
    @test pa==(3.0,0.0,0.0) && pb==(1.0,0.0,0.0)

    model=GeoModel()
    add_box!(model,0,0,0,1,1,1;tag=1)
    add_point!(model,0.2,0.5,0.5;tag=50)
    _,pa,pb=_spatial_check_distance(model,0,50,3,1,0.0)
    # Containment returns the interior point on both entities (OCC convention).
    @test pa==pb==(0.2,0.5,0.5)

    model=GeoModel()
    add_box!(model,0,0,0,1,1,1;tag=1)
    add_box!(model,4,0,0,1,1,1;tag=2)
    _spatial_check_distance(model,3,1,3,2,3.0)

    model=GeoModel()
    add_box!(model,0,0,0,1,1,1;tag=1)
    add_sphere!(model,4,0.5,0.5,1.0;tag=2)
    _spatial_check_distance(model,3,1,3,2,2.0;atol=1e-8)

    model=GeoModel()
    add_sphere!(model,0,0,0,1.0;tag=1)
    add_sphere!(model,3,0,0,1.0;tag=2)
    _spatial_check_distance(model,3,1,3,2,1.0;atol=1e-8)

    model=GeoModel()
    add_box!(model,0,0,0,2,2,2;tag=1)
    add_sphere!(model,1,1,3,1.0;tag=2)
    _spatial_check_distance(model,3,1,3,2,0.0;atol=1e-8)

    model=GeoModel()
    add_cylinder!(model,0,0,0,0,0,3,1.0;tag=1)
    add_box!(model,4,0,0,1,1,1;tag=2)
    _spatial_check_distance(model,3,1,3,2,3.0;atol=1e-8)

    # Full torus (z-axis, R=3, r=1): the center point is 2 from the rim.
    model=GeoModel()
    add_torus!(model,0,0,0,3.0,1.0;tag=1)
    add_point!(model,0,0,0;tag=50)
    _spatial_check_distance(model,0,50,3,1,2.0;atol=1e-7)

    # A curve of the materialized tetrahedron to the volume is shared → 0.
    model=_spatial_tetrahedron()
    _spatial_check_distance(model,1,1,3,1,0.0)
    _spatial_check_distance(model,0,1,3,1,0.0)
    # Point (0.5,0.5,0.5) inside the tetrahedron (x/2+y/3+z/4<1) → 0.
    add_point!(model,0.5,0.5,0.5;tag=10)
    _spatial_check_distance(model,0,10,3,1,0.0)
    # Point (1,1,1) is just outside the slant face (13/12>1): distance
    # 1/√61 to the plane x/2+y/3+z/4=1 with the foot inside the triangle.
    add_point!(model,1,1,1;tag=12)
    _spatial_check_distance(model,0,12,3,1,1/sqrt(61.0);atol=1e-8)
    # Point (6,6,6) outside: distance to the tet boundary, matching
    # point↔surface minima.
    add_point!(model,6,6,6;tag=11)
    d,_,_=model_distance(model,0,11,3,1)
    dface=min(model_distance(model,0,11,2,1)[1],
              model_distance(model,0,11,2,2)[1],
              model_distance(model,0,11,2,3)[1],
              model_distance(model,0,11,2,4)[1])
    @test d≈dface atol=1e-9 rtol=0.0
end

@testset "entity distance — discrete entities" begin
    model=GeoModel()
    add_discrete_entity!(model,1,1)
    add_discrete_nodes!(model,1,1,[10,20],[0.,0,0,1,0,0])
    add_discrete_elements!(model,1,1,[1],[[50]],[[10,20]])
    add_point!(model,1,4,0;tag=1)
    d,pa,pb=model_distance(model,0,1,1,1)
    @test d≈4.0
    @test pb==(1.0,0.0,0.0)

    add_discrete_entity!(model,2,2)
    add_discrete_nodes!(model,2,2,[30,31,32],
                      [0.,0,3, 1,0,3, 0,1,3])
    add_discrete_elements!(model,2,2,[2],[[60]],[[30,31,32]])
    add_point!(model,0.25,0.25,0;tag=2)
    d,pa,pb=model_distance(model,0,2,2,2)
    @test d≈3.0
    @test all(pb.≈(0.25,0.25,3.0))
end

@testset "entity distance — validation and explicit blockers" begin
    model=_spatial_tetrahedron()
    for call in (
        ()->model_distance(model,-1,1,0,1),
        ()->model_distance(model,0,-1,0,1),
        ()->model_distance(model,4,1,0,1),
        ()->model_distance(model,0,1,4,1),
        ()->model_distance(model,true,1,0,1),
        ()->model_distance(model,0,1,0,true),
        ()->model_distance(model,0,99,0,1),
        ()->model_distance(model,0,1,0,99),
    )
        @test_throws ArgumentError call()
    end
    @test_throws ArgumentError model_distance(GeoModel(),0,1,0,2)

    # Unsupported surface kind: explicit ArgumentError, never a silent
    # approximate fallback.
    corrupt=_spatial_tetrahedron()
    corrupt.surface_types[1]=:nurbs
    err=try
        model_distance(corrupt,0,1,2,1); nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("unsupported kind",sprint(showerror,err))

    # A volume with no boundary and no primitive encoding is an explicit
    # blocker rather than a silent zero.
    hollow=GeoModel()
    hollow.volumes[1]=Int[]
    err=try
        model_distance(hollow,0,1,3,1); nothing
    catch e
        e
    end
    @test err isa ArgumentError

    # Distance between a surface with no loops and a point must raise an
    # explicit blocker (ErrorException for corrupt state, ArgumentError for
    # unsupported input) rather than a silent approximate answer.
    broken=_spatial_tetrahedron()
    empty!(broken.surfaces[1])
    err=try
        model_distance(broken,0,1,2,1); nothing
    catch e
        e
    end
    @test err isa Union{ArgumentError,ErrorException}
end
