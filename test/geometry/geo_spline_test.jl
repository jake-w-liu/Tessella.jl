using Test
using Tessella
using Tessella.Model: model_value, model_derivative, model_second_derivative,
                      model_curvature, model_parametrization_bounds,
                      model_entity_type, model_bounding_box,
                      model_parametrization, model_closest_point,
                      model_is_inside, model_set_tag!, mesh_model_surface,
                      transform_entities!, remove_entities!, duplicate_entities!,
                      coherence!, _affine_translation

function _execute_spline_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _spline_error(source::AbstractString)
    try
        _execute_spline_source(source)
        return nothing
    catch err
        return err
    end
end

# The canonical five-point control polygon shared by most cases: a symmetric
# zig-zag on which Catmull-Rom, the UBS evaluator, De Casteljau, and a clamped
# cubic knot vector all produce distinct curves.
function _zigzag_model()
    m=GeoModel()
    for (tag,(x,y)) in enumerate(((0.0,0.0),(1.0,1.0),(2.0,0.0),
                                (3.0,1.0),(4.0,0.0)))
        add_point!(m,x,y,0.0;tag=tag)
    end
    return m
end

@testset "spline-family construction and records" begin
    m=_zigzag_model()
    sp=add_spline!(m,[1,2,3,4,5];tag=10)
    @test sp==10
    @test m.curves[10]==(1,5)
    @test m.curve_control_points[10]==[1,2,3,4,5]
    @test m.curve_types[10]==:spline
    @test model_entity_type(m,1,10)=="Nurb"
    bs=add_bspline!(m,[1,2,3,4,5])
    @test m.curve_types[bs]==:bspline
    bz=add_bezier!(m,[1,2,3,4])
    @test m.curve_types[bz]==:bezier
    @test m.curve_control_points[bz]==[1,2,3,4]
    nu=add_nurbs!(m,[1,2,3,4,5],
                  [0.0,0.0,0.0,0.0,1.0,2.0,2.0,2.0,2.0])
    @test m.curve_types[nu]==:nurbs
    @test m.curve_geometry[nu].deg==3
    @test (m.curve_geometry[nu].ubeg,m.curve_geometry[nu].uend)==(0.0,2.0)
    # knots are stored at float32 precision (Curve::k is `float*`)
    nu2=add_nurbs!(m,[1,2],[0.0,1.0/3.0,1.0])
    @test m.curve_geometry[nu2].knots==Float64.(Float32.([0.0,1.0/3.0,1.0]))
    @test m.curve_geometry[nu2].ubeg==0.0 && m.curve_geometry[nu2].uend==1.0
    @test m.curve_geometry[nu2].deg==0
end

@testset "spline-family construction errors" begin
    m=_zigzag_model()
    @test_throws ArgumentError add_spline!(m,[1])
    @test_throws ArgumentError add_bspline!(m,Int[])
    @test_throws ArgumentError add_bezier!(m,[1,99])
    @test_throws ArgumentError add_nurbs!(m,[1,2],[0.0,NaN])
    @test_throws ArgumentError add_nurbs!(m,[1,2],[0.0,-Inf])
    @test_throws ArgumentError add_nurbs!(m,[1,2],[1.0,0.0])
    @test_throws ArgumentError add_nurbs!(m,[1,2],[0.0,1.0])  # degree -1
    @test_throws ArgumentError add_nurbs!(m,[1,2],Bool[false,true])
    @test_throws ArgumentError add_spline!(m,[1,2];tag=-3)
    e=_spline_error("Point(1)={0,0,0};Point(2)={1,0,0};" *
                    "Spline(1)={1,2};Spline(1)={1,2};")
    @test e isa ArgumentError
end

@testset "Catmull-Rom spline evaluation" begin
    m=_zigzag_model()
    sp=add_spline!(m,[1,2,3,4,5])
    # uniform knots: the curve interpolates cp i+1 at u=i/(N-1)
    @test model_value(m,1,sp,[0.0])==[0.0,0.0,0.0]
    @test model_value(m,1,sp,[0.25])==[1.0,1.0,0.0]
    @test model_value(m,1,sp,[0.5])==[2.0,0.0,0.0]
    @test model_value(m,1,sp,[0.75])==[3.0,1.0,0.0]
    @test model_value(m,1,sp,[1.0])==[4.0,0.0,0.0]
    # mid-segment value — verified bit-for-bit against the pinned Gmsh 4.15.2
    @test model_value(m,1,sp,[0.3])==[1.1999999999999997,0.896,0.0]
    @test model_parametrization_bounds(m,1,sp)==([0.0],[1.0])
    @test model_entity_type(m,1,sp)=="Nurb"
    # `GEdge::bounds` is a 10-sample scan, not exact extrema
    @test collect(model_bounding_box(m,1,sp))==
        [0.0,0.0,0.0,4.0,0.9766803840877916,0.0]
    # two-point spline is a straight segment (ghost ends collapse)
    m2=GeoModel()
    add_point!(m2,0.0,0.0,0.0;tag=1); add_point!(m2,1.0,1.0,1.0;tag=2)
    sp2=add_spline!(m2,[1,2])
    @test model_value(m2,1,sp2,[0.5])==[0.5,0.5,0.5]
end

@testset "uniform BSpline evaluation" begin
    m=_zigzag_model()
    bs=add_bspline!(m,[1,2,3,4,5])
    @test model_value(m,1,bs,[0.0])==[0.0,0.0,0.0]
    @test model_value(m,1,bs,[1.0])==[4.0,0.0,0.0]
    # NbCurves=2 spans; the left-extremity matrix at t=0.5 lands mid-way
    @test model_value(m,1,bs,[0.25])==[1.1875,0.625,0.0]
    @test model_value(m,1,bs,[0.5])==[2.0,0.5,0.0]
    @test collect(model_bounding_box(m,1,bs))==
        [0.0,0.0,0.0,4.0,0.6200274348422496,0.0]
    # degenerate counts: 2 cps -> line, 3 -> quadratic Bezier, 4 -> cubic Bezier
    m2=GeoModel()
    add_point!(m2,0.0,0.0,0.0;tag=1); add_point!(m2,2.0,0.0,0.0;tag=2)
    @test model_value(m2,1,add_bspline!(m2,[1,2]),[0.25])==[0.5,0.0,0.0]
    m3=GeoModel()
    add_point!(m3,0.0,0.0,0.0;tag=1); add_point!(m3,1.0,2.0,0.0;tag=2)
    add_point!(m3,2.0,0.0,0.0;tag=3)
    @test model_value(m3,1,add_bspline!(m3,[1,2,3]),[0.5])==[1.0,1.0,0.0]
end

@testset "Bezier evaluation" begin
    m=GeoModel()
    for (tag,(x,y)) in enumerate(((0.0,0.0),(1.0,2.0),(3.0,2.0),(4.0,0.0)))
        add_point!(m,x,y,0.0;tag=tag)
    end
    bz=add_bezier!(m,[1,2,3,4])
    # exact cubic Bernstein values
    @test model_value(m,1,bz,[0.0])==[0.0,0.0,0.0]
    @test model_value(m,1,bz,[0.5])==[2.0,1.5,0.0]
    @test model_value(m,1,bz,[1.0])==[4.0,0.0,0.0]
    @test model_parametrization_bounds(m,1,bz)==([0.0],[1.0])
    @test model_entity_type(m,1,bz)=="Nurb"
end

@testset "Nurbs evaluation on its raw knot interval" begin
    m=_zigzag_model()
    nu=add_nurbs!(m,[1,2,3,4,5],
                  [0.0,0.0,0.0,0.0,1.0,2.0,2.0,2.0,2.0])
    # clamped cubic: interpolates the ends and follows the UBS segments
    @test model_value(m,1,nu,[0.0])==[0.0,0.0,0.0]
    @test model_value(m,1,nu,[0.5])==[1.1875,0.625,0.0]
    @test model_value(m,1,nu,[1.0])==[2.0,0.5,0.0]
    @test model_value(m,1,nu,[2.0])==[4.0,0.0,0.0]
    @test model_parametrization_bounds(m,1,nu)==([0.0],[2.0])
    @test collect(model_bounding_box(m,1,nu))==
        [0.0,0.0,0.0,4.0,0.6200274348422496,0.0]
    # degree-0 Nurbs is the piecewise-constant curve of `basisFuns`
    m2=GeoModel()
    add_point!(m2,0.0,0.0,0.0;tag=1); add_point!(m2,2.0,0.0,0.0;tag=2)
    n0=add_nurbs!(m2,[1,2],[0.0,0.5,1.0])
    @test model_value(m2,1,n0,[0.7])==[2.0,0.0,0.0]
    @test model_value(m2,1,n0,[0.2])==[0.0,0.0,0.0]
end

@testset "spline-family derivatives and curvature" begin
    m=_zigzag_model()
    sp=add_spline!(m,[1,2,3,4,5])
    # `InterpolateCurve` 1e-8 finite differences — pinned Gmsh 4.15.2 outputs
    @test model_derivative(m,1,sp,[0.3])==
        [3.9999999978945766,-3.8399999957583475,0.0]
    @test model_second_derivative(m,1,sp,[0.3])==
        [1.6653345369377348,-57.45404152435185,0.0]
    @test model_curvature(m,1,sp,[0.3])==[1.3105393943473993]
    # one-sided differences at the parameter bounds
    @test model_derivative(m,1,sp,[0.0])==[4.0,4.000000159999994,0.0]
    @test model_derivative(m,1,sp,[1.0])==[4.000000064507958,-4.000000169979145,0.0]
    nu=add_nurbs!(m,[1,2,3,4,5],
                  [0.0,0.0,0.0,0.0,1.0,2.0,2.0,2.0,2.0])
    @test model_derivative(m,1,nu,[0.5])==
        [1.8750000219114327,5.551115123125783e-9,0.0]
    @test model_curvature(m,1,nu,[0.5])==[0.8684411002522459]
end

@testset "spline projection, closest point, containment" begin
    m=_zigzag_model()
    sp=add_spline!(m,[1,2,3,4,5])
    # `parFromPoint` on the exact mid-segment evaluation point inverts to u
    @test model_parametrization(m,1,sp,[1.2,0.896,0.0])==
        [0.30000000000000004]
    closest,par=model_closest_point(m,1,sp,[2.0,0.5,0.0])
    @test closest==[1.6721139739491526,0.25202617931684557,0.0]
    @test par==[0.41802849348728816]
    @test model_is_inside(m,1,sp,[1.2,0.896,0.0])==1
    @test model_is_inside(m,1,sp,[0.0,5.0,0.0])==0
    @test model_is_inside(m,1,sp,[0.5],true)==1       # parametric: [0,1]
    @test model_is_inside(m,1,sp,[1.5],true)==0
    nu=add_nurbs!(m,[1,2,3,4,5],
                  [0.0,0.0,0.0,0.0,1.0,2.0,2.0,2.0,2.0])
    # parametric containment uses the raw knot interval
    @test model_is_inside(m,1,nu,[1.5],true)==1
    @test model_is_inside(m,1,nu,[2.5],true)==0
    @test model_is_inside(m,1,nu,[2.0,0.5,0.0])==1
    @test model_is_inside(m,1,nu,[9.0,9.0,0.0])==0
end

@testset "spline-family lifecycle" begin
    m=_zigzag_model()
    sp=add_spline!(m,[1,2,3,4,5])
    nu=add_nurbs!(m,[1,2,3,4,5],
                  [0.0,0.0,0.0,0.0,1.0,2.0,2.0,2.0,2.0])
    # transforms move the shared control-point vertices
    transform_entities!(m,_affine_translation((10.0,0.0,0.0),"test"),[(1,sp)])
    @test model_value(m,1,sp,[0.0])==[10.0,0.0,0.0]
    # the spline's +10x sweep already moved the shared control points
    transform_entities!(m,_affine_translation((0.0,5.0,0.0),"test"),[(1,nu)])
    @test model_value(m,1,nu,[1.0])==[12.0,5.5,0.0]
    # retagging preserves type, control points, and geometry records
    model_set_tag!(m,1,nu,77)
    @test m.curve_types[77]==:nurbs
    @test m.curve_geometry[77].deg==3
    @test model_value(m,1,77,[1.0])==[12.0,5.5,0.0]
    # removal drops the record and its control-point list
    remove_entities!(m,[(1,sp),(1,77)])
    @test !haskey(m.curves,sp) && !haskey(m.curve_control_points,sp)
    @test !haskey(m.curves,77) && !haskey(m.curve_geometry,77)
    coherence!(m)
    @test isempty(m.curves)
end

@testset "reversed Nurbs copies mirror knots and bounds" begin
    m=GeoModel()
    for (x,y) in ((0,0),(1,1),(2,0),(2,-1))
        add_point!(m,Float64(x),Float64(y),0.0)
    end
    n=add_nurbs!(m,[1,2,3],[0.0,0.0,1.0,2.0,2.0])
    l1=add_line!(m,3,4); l2=add_line!(m,4,1)
    # `Curve Loop = {-b,-a,-n}` traverses the Nurbs backwards; duplicating the
    # surface resolves the reversed record to a positive-tag copy, exactly
    # like Gmsh's `DuplicateSurface` on the `-n` record.
    lp=add_curve_loop!(m,[-l2,-l1,-n])
    sf=add_plane_surface!(m,[lp])
    dup=duplicate_entities!(m,[(2,sf)])
    newloop=m.loops[m.surfaces[dup[1][2]][1]]
    newcurve=last(newloop)
    @test newcurve>0
    @test m.curve_types[newcurve]==:nurbs
    g=m.curve_geometry[newcurve]
    @test g.knots==reverse(m.curve_geometry[n].knots)
    @test (g.ubeg,g.uend)==(1.0-2.0,1.0-0.0)
    # the copy's control-point positions run in reversed order
    @test [m.points[p] for p in m.curve_control_points[newcurve]]==
        reverse([m.points[p] for p in m.curve_control_points[n]])
end

@testset ".geo spline-family statements" begin
    parsed=_execute_spline_source("""
        Point(1) = {0,0,0};
        Point(2) = {1,1,0};
        Point(3) = {2,0,0};
        Point(4) = {3,1,0};
        Point(5) = {4,0,0};
        pts[] = {1,2,3,4,5};
        kk[] = {0,0,0,0,1,1,1,1};
        Spline(10) = pts[];
        BSpline(11) = {1,2,3,4,5};
        Bezier(12) = {1,2,3,4};
        Nurbs(13) = {1,2,3,4,5} Knots {0,0,0,0,1,2,2,2,2} Order 3;
        Nurbs(14) = {1,2,3} Knots {} Order 1;
        Nurbs(15) = pts[] Knots kk[] Order 99;
        """)
    m=parsed.model
    @test model_entity_type(m,1,10)=="Nurb"
    @test model_entity_type(m,1,11)=="Nurb"
    @test model_entity_type(m,1,12)=="Nurb"
    @test model_entity_type(m,1,13)=="Nurb"
    @test m.curve_types[10]==:spline
    @test m.curve_types[11]==:bspline
    @test m.curve_types[12]==:bezier
    @test m.curve_types[13]==:nurbs
    @test model_parametrization_bounds(m,1,13)==([0.0],[2.0])
    # `Nurbs ... Knots {}` builds the plain BSpline record
    @test m.curve_types[14]==:bspline
    # `ListOfDouble` positions accept list variables, and `Order` is parsed
    # but ignored by the built-in kernel (degree comes from the knot count)
    @test m.curve_types[15]==:nurbs
    @test m.curve_geometry[15].deg==2
    @test (m.curve_geometry[15].ubeg,m.curve_geometry[15].uend)==(0.0,1.0)
    @test model_value(m,1,10,[0.3])==[1.1999999999999997,0.896,0.0]
    @test model_value(m,1,11,[0.25])==[1.1875,0.625,0.0]
    @test model_value(m,1,12,[0.5])==[1.5,0.5,0.0]
    @test model_value(m,1,13,[1.0])==[2.0,0.5,0.0]
    @test model_value(m,1,15,[0.5])==[2.0,0.5,0.0]
end

@testset ".geo spline statement failures" begin
    @test _spline_error("""
        Point(1) = {0,0,0};
        Point(2) = {1,0,0};
        Nurbs(1) = {1,2} Order 1;
        """) isa ArgumentError
    @test _spline_error("""
        Point(1) = {0,0,0};
        Point(2) = {1,0,0};
        Nurbs(1) = {1,2} Knots {0,1};
        """) isa ArgumentError
    @test _spline_error("""
        Point(1) = {0,0,0};
        Point(2) = {1,0,0};
        Nurbs(1) = {1,2} Knots {0,1} Order;
        """) isa ArgumentError
    # Upstream stores a malformed Nurbs record without validation — an
    # over-sized degree/knot set parses silently (Gmsh 4.15.2 unrolls
    # `Nurbs(1) = {1,2,3} Knots {...} Order 4` unchanged).
    @test _execute_spline_source("""
        Point(1) = {0,0,0};
        Point(2) = {1,0,0};
        Point(3) = {2,0,0};
        Nurbs(1) = {1,2,3} Knots {0,0,0,0,0.5,1,1,1} Order 4;
        """) !== nothing
    @test _spline_error("""
        Point(1) = {0,0,0};
        Point(2) = {1,0,0};
        Spline(1) = {1,9};
        """) isa ArgumentError
    @test _spline_error("""
        Point(1) = {0,0,0};
        Spline(1) = {1};
        """) isa ArgumentError
end

@testset "spline-family curves reject straight-line meshing" begin
    m=_zigzag_model()
    sp=add_spline!(m,[1,2,3,4,5])
    close=add_line!(m,5,1)
    lp=add_curve_loop!(m,[sp,close])
    s=add_plane_surface!(m,[lp])
    # `_model_require_line_curve` — curved boundaries must not silently
    # degrade to chords in the straight-line surface kernel.
    @test_throws ArgumentError mesh_model_surface(m,s)
end
