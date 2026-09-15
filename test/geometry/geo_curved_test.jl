using Test
using Tessella
using Tessella.Model: model_value, model_derivative, model_second_derivative,
                      model_curvature, model_parametrization_bounds,
                      model_entity_type, model_bounding_box, model_set_tag!,
                      mesh_model_surface, transform_entities!, remove_entities!,
                      duplicate_entities!, coherence!, _affine_translation,
                      _affine_rotation

function _execute_curved_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _curved_error(source::AbstractString)
    try
        _execute_curved_source(source)
        return nothing
    catch err
        return err
    end
end

function _quarter_model()
    m=GeoModel()
    add_point!(m,1.0,0.0,0.0;tag=1)
    add_point!(m,0.0,0.0,0.0;tag=2)
    add_point!(m,0.0,1.0,0.0;tag=3)
    return m
end

@testset "circle arc construction and evaluation" begin
    m=_quarter_model()
    c=add_circle_arc!(m,1,2,3;tag=1)
    @test c==1
    @test m.curves[1]==(1,3)
    @test m.curve_control_points[1]==[1,2,3]
    @test m.curve_types[1]==:circle
    @test model_entity_type(m,1,c)=="Circle"
    @test model_parametrization_bounds(m,1,c)==([0.0],[1.0])
    @test model_value(m,1,c,[0.0])≈[1.0,0.0,0.0]
    @test model_value(m,1,c,[1.0])≈[0.0,1.0,0.0] atol=1e-15
    @test model_value(m,1,c,[0.5])≈[sqrt(0.5),sqrt(0.5),0.0]
    # Gmsh's `InterpolateCurve` derivatives are 1e-8-step finite differences,
    # not the analytic arc derivatives — these are the 4.15.2 binary's values.
    @test model_derivative(m,1,c,[0.5])==[-1.1107207320559809,1.110720737607096,0.0]
    @test model_curvature(m,1,c,[0.3])==[0.9805697128486208]
    @test model_curvature(m,1,c,[0.9])==[1.137437577809078]
    d2=model_second_derivative(m,1,c,[0.5])
    @test d2==[-2.220446049250313,-1.1102230246251565,0.0]
    @test collect(model_bounding_box(m,1,c))≈[0.0,0.0,0.0,1.0,1.0,0.0]
end

@testset "circle arcs in arbitrary planes" begin
    m=GeoModel()
    add_point!(m,0.0,0.0,1.0;tag=1)
    add_point!(m,0.0,0.0,0.0;tag=2)
    add_point!(m,0.0,1.0,0.0;tag=3)
    c=add_circle_arc!(m,1,2,3;tag=1)
    v=model_value(m,1,c,[0.5])
    @test v≈[0.0,sqrt(0.5),sqrt(0.5)]
    @test collect(model_bounding_box(m,1,c))≈[0.0,0.0,0.0,0.0,1.0,1.0]
    # A tilted arc whose plane normal comes only from the three points.
    m2=GeoModel()
    add_point!(m2,1.0,1.0,0.0;tag=1)
    add_point!(m2,0.0,0.0,0.0;tag=2)
    add_point!(m2,0.0,1.0,1.0;tag=3)
    c2=add_circle_arc!(m2,1,2,3;tag=1)
    @test model_value(m2,1,c2,[0.0])≈[1.0,1.0,0.0]
    @test model_value(m2,1,c2,[1.0])≈[0.0,1.0,1.0] atol=1e-15
    mid=model_value(m2,1,c2,[0.5])
    @test sqrt(sum(x->x^2,collect(mid)))≈sqrt(2.0)
end

@testset "semicircle uses the stored Plane normal" begin
    m=GeoModel()
    add_point!(m,1.0,0.0,0.0;tag=1)
    add_point!(m,0.0,0.0,0.0;tag=2)
    add_point!(m,-1.0,0.0,0.0;tag=3)
    # Collinear endpoints: the arc plane is defined by `plane_normal`, and a
    # (0,0,-1) normal sends the half arc through -y.
    c=add_circle_arc!(m,1,2,3;tag=1,plane_normal=(0.0,0.0,-1.0))
    @test model_value(m,1,c,[0.5])≈[0.0,-1.0,0.0] atol=1e-15
    @test model_value(m,1,c,[0.25])≈[sqrt(0.5),-sqrt(0.5),0.0]
    # The default stored normal is (0,0,1) like Gmsh.
    c2=add_circle_arc!(m,3,2,1;tag=2)
    @test model_value(m,1,c2,[0.5])≈[0.0,-1.0,0.0] atol=1e-15
end

@testset "ellipse arcs" begin
    m=GeoModel()
    add_point!(m,2.0,0.0,0.0;tag=1)
    add_point!(m,0.0,0.0,0.0;tag=2)
    add_point!(m,2.0,0.0,0.0;tag=5)
    add_point!(m,0.0,1.0,0.0;tag=3)
    # start=(2,0), center=(0,0), major direction +x, end=(0,1):
    # solving both endpoints on the axis-aligned ellipse gives a=2, b=1.
    e=add_ellipse_arc!(m,1,2,5,3;tag=1)
    @test m.curves[1]==(1,3)
    @test m.curve_control_points[1]==[1,2,5,3]
    @test m.curve_types[1]==:ellipse
    @test model_entity_type(m,1,e)=="Ellipse"
    @test model_value(m,1,e,[0.0])≈[2.0,0.0,0.0]
    @test model_value(m,1,e,[1.0])≈[0.0,1.0,0.0] atol=1e-15
    # Parametric angle 45°: (2 cos45, 1 sin45).
    @test model_value(m,1,e,[0.5])≈[2*sqrt(0.5),sqrt(0.5),0.0]
    @test collect(model_bounding_box(m,1,e))≈[0.0,0.0,0.0,2.0,1.0,0.0]
    # Gmsh's built-in 3-tag form maps major=start.
    m3=GeoModel()
    add_point!(m3,1.0,0.0,0.0;tag=1)
    add_point!(m3,0.0,0.0,0.0;tag=2)
    add_point!(m3,0.0,1.0,0.0;tag=3)
    e3=add_ellipse_arc!(m3,1,2,1,3;tag=1)
    @test m3.curve_control_points[1]==[1,2,1,3]
    @test model_value(m3,1,e3,[0.5])≈[sqrt(0.5),sqrt(0.5),0.0]
end

@testset "arc construction failures" begin
    m=_quarter_model()
    @test_throws ArgumentError add_circle_arc!(m,1,9,3)
    @test_throws ArgumentError add_ellipse_arc!(m,1,9,1,3)
    # Non-cocircular circle points fail like Gmsh's EndCurve check.
    add_point!(m,0.0,0.5,0.0;tag=4)
    @test_throws ArgumentError add_circle_arc!(m,1,2,4)
    # Zero radius: start == center.
    add_point!(m,1.0,0.0,0.0;tag=5)
    @test_throws ArgumentError add_circle_arc!(m,5,1,3)
    # Gmsh's "greater than Pi" gate is retained in `_arc_geometry` for
    # parity, but the start-normalized parametrization keeps every
    # solvable arc's span at or below Pi, so it is not reachable through
    # the public constructors.
    # A "wrong" ellipse: with the major point on the diagonal, the
    # antiparallel end point makes the sys2x2 axis solve singular, like
    # Gmsh's "Ellipse is wrong".
    m2=GeoModel()
    add_point!(m2,1.0,0.0,0.0;tag=1)
    add_point!(m2,0.0,0.0,0.0;tag=2)
    add_point!(m2,1.0,1.0,0.0;tag=5)
    add_point!(m2,-1.0,0.0,0.0;tag=3)
    @test_throws ArgumentError add_ellipse_arc!(m2,1,2,5,3)
end

@testset "curve loop sorting and closure" begin
    m=GeoModel()
    add_point!(m,0.0,0.0,0.0;tag=1)
    add_point!(m,1.0,0.0,0.0;tag=2)
    add_point!(m,1.0,1.0,0.0;tag=3)
    add_point!(m,0.0,1.0,0.0;tag=4)
    add_line!(m,1,2;tag=1); add_line!(m,2,3;tag=2)
    add_line!(m,3,4;tag=3); add_line!(m,4,1;tag=4)
    # Input order is irrelevant: SortEdgesInLoop chains forward from the
    # first tag by oriented endpoint connectivity (it does not rotate the
    # chain to the smallest tag).
    l=add_curve_loop!(m,[3,1,4,2];tag=1)
    @test m.loops[l]==[3,4,1,2]
    l2=add_curve_loop!(m,[-2,-4,-1,-3];tag=2)
    @test m.loops[l2]==[-2,-1,-4,-3]
    # A chain that dead-ends is rejected even though every tag is known:
    # curve 5 runs 1->3, so after 1->2->3 no member starts at 3's end.
    add_line!(m,1,3;tag=5)
    @test_throws ArgumentError add_curve_loop!(m,[1,2,5])
    # Loop tags live in their own namespace like Gmsh's _maxLineLoopNum:
    # allocating a loop must not advance the curve counter.
    m2=GeoModel()
    add_point!(m2,0.0,0.0,0.0;tag=1); add_point!(m2,1.0,0.0,0.0;tag=2)
    add_line!(m2,1,2;tag=1)
    add_curve_loop!(m2,[1])
    @test add_line!(m2,2,1)==2
end

@testset "ruled and triangular surfaces" begin
    m=GeoModel()
    add_point!(m,1.0,0.0,0.0;tag=1)
    add_point!(m,0.0,0.0,0.0;tag=2)
    add_point!(m,0.0,1.0,0.0;tag=3)
    add_point!(m,-1.0,0.0,0.0;tag=4)
    add_point!(m,0.0,-1.0,0.0;tag=5)
    add_circle_arc!(m,1,2,3;tag=1)
    add_circle_arc!(m,3,2,4;tag=2)
    add_circle_arc!(m,4,2,5;tag=3)
    add_circle_arc!(m,5,2,1;tag=4)
    l4=add_curve_loop!(m,[1,2,3,4])
    s=add_ruled_surface!(m,[l4];tag=1)
    @test m.surface_types[s]==:ruled
    @test model_entity_type(m,2,s)=="Surface"
    m3=GeoModel()
    for (tag,p) in ((1,(1.0,0.0,0.0)),(2,(0.0,0.0,0.0)),(3,(0.0,1.0,0.0)),
                    (4,(-1.0,0.0,0.0)))
        add_point!(m3,p...;tag=tag)
    end
    add_circle_arc!(m3,1,2,3;tag=1)
    add_circle_arc!(m3,3,2,4;tag=2)
    add_line!(m3,4,1;tag=3)
    l3=add_curve_loop!(m3,[1,2,3])
    t=add_ruled_surface!(m3,[l3];tag=1,sphere_center=2)
    @test m3.surface_types[t]==:tric
    @test model_entity_type(m3,2,t)=="Surface"
    @test m3.surface_geometry[t].sphere_center==2
    # 2-border and 5-border loops are rejected like Gmsh's addSurfaceFilling.
    m4=GeoModel()
    add_point!(m4,0.0,0.0,0.0;tag=1); add_point!(m4,1.0,0.0,0.0;tag=2)
    add_line!(m4,1,2;tag=1); add_line!(m4,2,1;tag=2)
    l2=add_curve_loop!(m4,[1,2])
    @test_throws ArgumentError add_ruled_surface!(m4,[l2])
    @test_throws ArgumentError add_ruled_surface!(m4,Int[];tag=9)
    @test_throws ArgumentError add_ruled_surface!(m4,[l3])
end

@testset "curved paths fail explicitly instead of degrading to chords" begin
    m=GeoModel()
    add_point!(m,1.0,0.0,0.0;tag=1)
    add_point!(m,0.0,0.0,0.0;tag=2)
    add_point!(m,0.0,1.0,0.0;tag=3)
    add_point!(m,-1.0,0.0,0.0;tag=4)
    add_circle_arc!(m,1,2,3;tag=1)
    add_circle_arc!(m,3,2,4;tag=2)
    add_line!(m,4,1;tag=3)
    l=add_curve_loop!(m,[1,2,3])
    s=add_ruled_surface!(m,[l];tag=1)
    err=try mesh_model_surface(m,s); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("Circle", sprint(showerror,err))
end

@testset "arc lifecycle" begin
    m=_quarter_model()
    add_point!(m,0.0,-1.0,0.0;tag=4)
    c=add_circle_arc!(m,1,2,3;tag=1)
    # Transforms move the center control point; geometry is re-derived.
    transform_entities!(m,_affine_translation((0.0,0.0,5.0),"test"),[(0,2)])
    mid=model_value(m,1,c,[0.5])
    @test mid[3] != 0.0
    # Retagging preserves type, control points, and geometry records.
    model_set_tag!(m,1,c,10)
    @test model_entity_type(m,1,10)=="Circle"
    @test m.curve_control_points[10]==[1,2,3]
    @test !haskey(m.curve_control_points,1)
    # An arc's center point is curve-owned: removing it is skipped like
    # Gmsh's boundary-owner protection.
    remove_entities!(m,[(0,2)])
    @test haskey(m.points,2) && haskey(m.curves,10)
    # Duplicata copies the arc with fresh control points.
    dup=duplicate_entities!(m,[(1,10)];caller="test")
    @test length(dup)==1
    @test model_entity_type(m,1,dup[1][2])=="Circle"
    # Coherence merges an identical reversed arc onto the kept record.
    m2=GeoModel()
    add_point!(m2,1.0,1.0,0.0;tag=5)
    add_point!(m2,0.0,1.0,0.0;tag=6)
    add_point!(m2,0.0,2.0,0.0;tag=7)
    add_circle_arc!(m2,5,6,7;tag=1)
    add_circle_arc!(m2,7,6,5;tag=2)
    coherence!(m2)
    @test sort!(collect(keys(m2.curves)))==[1]
    @test m2.curve_types==Dict(1=>:circle)
end

@testset ".geo curved entity execution" begin
    parsed=_execute_curved_source("""
        Point(1) = {1,0,0};
        Point(2) = {0,0,0};
        Point(3) = {0,1,0};
        Point(4) = {-1,0,0};
        Point(5) = {0,-1,0};
        Point(9) = {0,0,0.5};
        Circle(1) = {1,2,3};
        Circle(2) = {3,2,4} Plane {0,0,-1};
        Ellipse(3) = {1,2,5,3};
        Ellipse(4) = {1,2,3};
        Line(5) = {4,1};
        Line(6) = {4,3};
        Curve Loop(1) = {1,2,6,-4};
        Curve Loop(2) = {1,2,5};
        Surface(1) = {1} In Sphere {2};
        Ruled Surface(2) = {2} Using Point {9};
        Plane Surface(3) = {1};
        """)
    m=parsed.model
    @test model_entity_type(m,1,1)=="Circle"
    @test model_entity_type(m,1,2)=="Circle"
    @test model_entity_type(m,1,3)=="Ellipse"
    @test model_entity_type(m,1,4)=="Ellipse"
    @test model_entity_type(m,1,5)=="Line"
    @test m.loops[1]==[1,2,6,-4]
    @test m.loops[2]==[1,2,5]
    @test model_entity_type(m,2,1)=="Surface"
    @test m.surface_types[1]==:ruled
    @test m.surface_geometry[1].sphere_center==2
    @test model_entity_type(m,2,2)=="Surface"
    @test m.surface_types[2]==:tric
    @test m.surface_geometry[2].sphere_center==9
    @test model_entity_type(m,2,3)=="Plane"
    @test Tessella.Model._surface_type(m,3)==:plane
    # Ellipse(4)={1,2,3}: 3-tag form, major=start -> unit circle arc.
    @test m.curve_control_points[4]==[1,2,1,3]
end

@testset ".geo curved statement failures" begin
    @test_throws ArgumentError _execute_curved_source("""
        Point(1) = {1,0,0};
        Circle(1) = {1,1,1};
        """)
    @test_throws ArgumentError _execute_curved_source("""
        Point(1) = {1,0,0};
        Point(2) = {0,0,0};
        Circle(1) = {1,2};
        """)
    @test_throws ArgumentError _execute_curved_source("""
        Point(1) = {1,0,0};
        Point(2) = {0,0,0};
        Point(3) = {0,1,0};
        Circle(1) = {1,2,3} Sphere {0,0,1};
        """)
    @test_throws ArgumentError _execute_curved_source("""
        Point(1) = {0,0,0};
        Point(2) = {1,0,0};
        Line(1) = {1,2};
        Line(2) = {2,1};
        Curve Loop(1) = {1,2};
        Surface(1) = {1};
        """)
end

@testset "arc extrusion produces Gmsh surface kinds" begin
    parsed=_execute_curved_source("""
        Point(1) = {1,0,0};
        Point(2) = {0,0,0};
        Point(3) = {0,1,0};
        Circle(1) = {1,2,3};
        out[] = Extrude {0,0,1} { Curve{1}; };
        """)
    m=parsed.model
    # Gmsh: out = [top curve, lateral surface, generatrices...].
    @test parsed.lists["out"]==[2.0,5.0,4.0,-3.0]
    @test model_entity_type(m,1,2)=="Circle"
    @test model_entity_type(m,2,5)=="Surface"
    @test m.surface_types[5]==:ruled
end
