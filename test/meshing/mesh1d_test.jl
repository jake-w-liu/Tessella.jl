# ── Stage-2 CRC suite: size fields + 1-D graded edge meshing ────────────────────
#
# Correctness  : analytic arc lengths; Gmsh 4.15.2 adaptive integration,
#                smoothing, count and placement oracles; graded spacing.
# Robustness   : resource caps, invalid controls, tiny/large sizes, closed curves.
# Completeness : scalar and BAMG policies, endpoint contexts, fixed-grid mode.

using Test
import Tessella
using Tessella.SizeField
using Tessella.Mesh1D
using Tessella.IO: GeoParams, GeoFieldSpec

dist(a,b) = sqrt(sum((a[i]-b[i])^2 for i in 1:3))

struct _EndpointMetricField <: AbstractAnisoField end

Tessella.SizeField.metric_at(::_EndpointMetricField, x, y, z) =
    isotropic_metric(0.25)
function Tessella.SizeField.metric_at(::_EndpointMetricField, x, y, z,
                                     entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}
    context = (Int(entity[1]), Int(entity[2]))
    h = context == (0, 11) ? 0.5 : context == (0, 12) ? 1.0 :
        context == (1, 7) ? 0.25 : 2.0
    return isotropic_metric(h)
end

@testset "SizeField + Mesh1D (Stage 2)" begin

    @testset "size fields evaluate correctly" begin
        @test size_at(ConstantSize(0.5), 1.0, 2.0, 3.0) == 0.5
        @test size_at(FunctionSize((x,y,z)->x+1.0), 4.0, 0.0, 0.0) == 5.0
        mn = MinSize([ConstantSize(0.3), FunctionSize((x,y,z)->0.1+x)])
        @test size_at(mn, 0.05, 0.0, 0.0) ≈ 0.15    # 0.1+0.05 < 0.3
        @test size_at(mn, 1.0, 0.0, 0.0) ≈ 0.3      # 0.3 < 1.1
        @test_throws ArgumentError ConstantSize(-1.0)
        @test_throws ArgumentError ConstantSize(Inf)
        @test_throws ArgumentError ConstantSize(NaN)
        @test_throws ArgumentError size_at(FunctionSize((x,y,z)->-1.0), 0.0,0.0,0.0)
        @test_throws ArgumentError size_at(FunctionSize((x,y,z)->"bad"), 0.0,0.0,0.0)
        @test_throws ArgumentError size_at(FunctionSize((x,y,z)->big(10)^1000), 0.0,0.0,0.0)
    end

    @testset "arc length vs analytic" begin
        seg(t) = (t*3.0, t*4.0, 0.0)                 # length 5
        @test curve_length(seg) ≈ 5.0 rtol=1e-9
        circ(t) = (2.0*cospi(2t), 2.0*sinpi(2t), 0.0)  # circumference 4π
        @test curve_length(circ; nsample=20000) ≈ 4π rtol=1e-6
    end

    @testset "curve contract and resource bounds" begin
        seg(t) = (t, 0.0, 0.0)
        @test_throws ArgumentError curve_length(seg; nsample=0)
        @test_throws ArgumentError curve_length(seg; nsample=big(typemax(Int))+1)
        @test curve_length(seg; nsample=big(2)) == 1.0
        @test_throws ArgumentError curve_length(seg; t0=1, t1=0)
        @test_throws ArgumentError curve_length(t->(NaN,0.0,0.0); nsample=2)
        @test_throws ArgumentError metric_length(seg, ConstantSize(1.0); nsample=-1)
        @test_throws ArgumentError mesh_curve(t->(0.0,0.0,0.0), ConstantSize(1.0))
        @test_throws ArgumentError mesh_curve(seg, ConstantSize(1.0); closed=true)
        @test_throws ArgumentError mesh_segment((0.0,0.0,0.0), (1.0,0.0,0.0),
                                                ConstantSize(1e-300); nsample=1)
        tiny(t)=(1e-320*t,0.0,0.0)
        @test metric_length(tiny,ConstantSize(1e-320);nsample=1) ≈ 1.0
        @test_throws ArgumentError metric_length(seg,ConstantSize(1.0);
            max_integration_points=128)
        @test_throws ArgumentError metric_length(seg,ConstantSize(1.0);
            min_integration_depth=8,max_integration_depth=7)
        @test_throws ArgumentError metric_length(seg,ConstantSize(1.0);
            smooth_ratio=0.9)
        @test_throws ArgumentError mesh_segment((0,0,0),(1,0,0),ConstantSize(.01);
            max_edges=10)
        @test_throws ArgumentError curve_length(seg;nsample=100,
            max_integration_points=100)
        @test_throws ArgumentError metric_length(seg,ConstantSize(1.0);
            derivative=t->(Inf,0,0))
    end

    @testset "constant size → uniform spacing, correct count" begin
        A=(0.0,0.0,0.0); B=(10.0,0.0,0.0)
        pts, par = mesh_segment(A, B, ConstantSize(1.0))
        @test length(pts) == 11                       # 10 edges
        @test pts[1] == A && pts[end] == B            # endpoints included
        gaps = [dist(pts[i],pts[i+1]) for i in 1:length(pts)-1]
        @test all(g -> isapprox(g, 1.0; atol=1e-6), gaps)   # uniform ≈ h
    end

    @testset "metric length = L/h" begin
        seg(t) = (t*12.0, 0.0, 0.0)
        @test metric_length(seg, ConstantSize(2.0)) ≈ 6.0 rtol=1e-6
    end

    @testset "curve and endpoint entity contexts" begin
        seg(t) = (t, 0.0, 0.0)
        point_size = ConstantField(; vin=0.2, vout=2.0, points=[11],
                                   include_boundary=false)
        curve_size = ConstantField(; vin=0.5, vout=2.0, curves=[7],
                                   include_boundary=false)
        scalar = MinSize((point_size, curve_size))

        # The begin vertex is finer than the curve; the unlisted end vertex is
        # coarser and cannot override the scalar curve constraint.
        @test metric_length(seg, scalar; nsample=2, entity=(1,7)) ≈ 2.0
        @test metric_length(seg, scalar; nsample=2, entity=(1,7),
            endpoint_entities=((0,11),(0,12)),smooth_ratio=nothing) ≈ 2.75
        @test metric_length(seg, scalar; nsample=2, entity=(1,7),
            endpoint_entities=((0,11),(0,12))) > 2.75
        @test metric_length(seg, scalar; nsample=2, entity=(1,7),
            endpoint_entities=(nothing,nothing)) ==
            metric_length(seg, scalar; nsample=2, entity=(1,7))

        segment_points, segment_params = mesh_segment(
            (0.0,0.0,0.0), (1.0,0.0,0.0), scalar; nsample=2, entity=(1,7),
            endpoint_entities=((0,11),(0,12)))
        curve_points, curve_params = mesh_curve(seg, scalar; nsample=2,
            derivative=t->(1.0,0.0,0.0), entity=(1,7),
            endpoint_entities=((0,11),(0,12)))
        @test length(segment_points) == length(segment_params) == 5
        @test segment_points == curve_points
        @test segment_params == curve_params

        # Gmsh's anisotropic edge integrand selects the vertex metric at each
        # endpoint instead of intersecting it with the curve metric.
        aniso = _EndpointMetricField()
        @test metric_length(seg, aniso; nsample=2, entity=(1,7)) ≈ 4.0
        @test metric_length(seg, aniso; nsample=2, entity=(1,7),
            endpoint_entities=((0,11),(0,12))) ≈ 4.0
        @test metric_length(seg, aniso; nsample=2, entity=(1,7),
            endpoint_entities=((0,11),(0,12)),anisotropic_metric=true) ≈ 2.75

        gmsh_policy=MathEvalAnisoField(m11="25",m22="25",m33="25")
        default_points,_=mesh_segment((0,0,0),(10,0,0),gmsh_policy;nsample=20)
        metric_points,_=mesh_segment((0,0,0),(10,0,0),gmsh_policy;
                                     nsample=20,anisotropic_metric=true)
        @test length(default_points)==2
        @test length(metric_points)==51

        @test_throws ArgumentError metric_length(seg, scalar;
            endpoint_entities=(0,11))
        @test_throws ArgumentError mesh_curve(seg, scalar;
            endpoint_entities=((1,11),nothing))
        @test_throws ArgumentError mesh_segment((0,0,0), (1,0,0), scalar;
            endpoint_entities=((0,0),nothing))
        @test_throws ArgumentError metric_length(seg, scalar;
            endpoint_entities=((0,big(typemax(Int))+1),nothing))
    end

    @testset "Gmsh 4.15.2 edge integration and placement parity" begin
        # Oracle coordinates below were generated with local Gmsh 4.15.2-git,
        # Mesh.{MeshSizeFromPoints,MeshSizeFromCurvature,
        # MeshSizeExtendFromBoundary}=0. They also trace meshGEdge.cpp's
        # Integration/smoothPrimitive/int(a+1.99) code paths.
        steep = FunctionSize((x,y,z)->0.05 + 0.8x)
        algorithm6_points, algorithm6_params = mesh_segment(
            (0,0,0),(1,0,0),steep;anisotropic_metric=false)
        algorithm7_points, algorithm7_params = mesh_segment(
            (0,0,0),(1,0,0),steep;anisotropic_metric=true)
        gmsh6 = [0.0, 0.05002998384178375, 0.1171624986079215,
                 0.2072439727060804, 0.3281193875826542,
                 0.4903155731930917, 0.7079577637354404, 1.0]
        gmsh7 = [0.0, 0.06440895201633101, 0.1951941025622427,
                 0.4607590133545323, 1.0]
        @test algorithm6_params ≈ gmsh6 atol=2e-11 rtol=0
        @test algorithm7_params ≈ gmsh7 atol=2e-11 rtol=0
        @test first.(algorithm6_points) ≈ gmsh6 atol=2e-11 rtol=0
        @test first.(algorithm7_points) ≈ gmsh7 atol=2e-11 rtol=0

        # The Gmsh count is ceil-like with a 0.01 tolerance, not round(a).
        threshold_points,_ = mesh_segment(
            (0,0,0),(10,0,0),ConstantSize(10 / 2.49))
        @test length(threshold_points) == 4

        point_size = ConstantField(;vin=0.2,vout=2.0,points=[1],
                                   include_boundary=false)
        curve_size = ConstantField(;vin=0.5,vout=2.0,curves=[1],
                                   include_boundary=false)
        endpoint_field = MinSize((point_size,curve_size))
        endpoint6, endpoint6_params = mesh_segment(
            (0,0,0),(1,0,0),endpoint_field;entity=(1,1),
            endpoint_entities=((0,1),(0,2)),anisotropic_metric=false)
        endpoint7, endpoint7_params = mesh_segment(
            (0,0,0),(1,0,0),endpoint_field;entity=(1,1),
            endpoint_entities=((0,1),(0,2)),anisotropic_metric=true)
        @test length(endpoint6) == 4
        @test endpoint6_params ≈
            [0.0,0.2314865831435398,0.5544886346204482,1.0] atol=2e-11 rtol=0
        @test length(endpoint7) == 3
        @test endpoint7_params ≈
            [0.0,0.499999991616791,1.0] atol=2e-11 rtol=0

        # End-to-end .geo background fields exercise Gmsh's model-size wrapper
        # and post-placement filterPoints pass. FinishUpBoundingBox pads a unit
        # line to model size √5; a raw ConstantField(VOut=10) is intentionally
        # not equivalent because it has no model bounding box.
        for (vin,gmsh_params) in (
                (0.02,[0.0,0.7087713251839173,1.0]),
                (0.2,[0.0,1.0]))
            params = GeoParams(NaN,NaN,1.0,0,
                Dict{Tuple{Int,Int},String}(),
                Dict(1=>GeoFieldSpec(1,"Constant",Dict(
                    "PointsList"=>"{1}","VIn"=>string(vin),"VOut"=>"10",
                    "IncludeBoundary"=>"0"))),1)
            point_only = build_geo_size_field(params,Dict();
                model_bbox=(0.0,1.0,0.0,0.0,0.0,0.0))
            @test size_at(point_only,0.5,0.0,0.0,(1,1)) ≈ sqrt(5.0)
            filtered, filtered_params = mesh_segment(
                (0,0,0),(1,0,0),point_only;entity=(1,1),
                endpoint_entities=((0,1),(0,2)),anisotropic_metric=false)
            @test length(filtered) == length(gmsh_params)
            @test filtered_params ≈ gmsh_params atol=2e-11 rtol=0
        end

        anisotropic = MathEvalAnisoField(m11="25",m22="4",m33="4")
        scalar_nodes,_ = mesh_segment((0,0,0),(10,0,0),anisotropic;
                                      anisotropic_metric=false)
        bamg_nodes,_ = mesh_segment((0,0,0),(10,0,0),anisotropic;
                                    anisotropic_metric=true)
        @test length(scalar_nodes) == 2
        @test length(bamg_nodes) == 51
    end

    @testset "graded size: local edge length tracks h" begin
        # size grows linearly along the segment; edge length should grow too, with
        # |edge| / h_mid ≈ constant (equal metric increments by construction).
        A=(0.0,0.0,0.0); B=(10.0,0.0,0.0)
        sf = FunctionSize((x,y,z) -> 0.2 + 0.3x)      # h from 0.2 to 3.2
        pts, par = mesh_segment(A, B, sf)
        gaps = [dist(pts[i],pts[i+1]) for i in 1:length(pts)-1]
        @test issorted(gaps; lt=(a,b)->a < b - 1e-9) || all(diff(gaps) .>= -1e-6)  # monotone increasing
        # ratio |edge|/h_mid nearly constant across the mesh (grading correctness)
        ratios = Float64[]
        for i in 1:length(pts)-1
            xm = 0.5*(pts[i][1]+pts[i+1][1]); h = 0.2+0.3xm
            push!(ratios, dist(pts[i],pts[i+1])/h)
        end
        @test maximum(ratios)/minimum(ratios) < 1.05  # within 5% (discretization)
        @test pts[1] == A && pts[end][1] ≈ 10.0
    end

    @testset "closed curve: nodes = edges, no duplicate" begin
        circ(t) = (5.0*cospi(2t), 5.0*sinpi(2t), 0.0)
        pts, par = mesh_curve(circ, ConstantSize(1.0); closed=true, nsample=8000)
        n = length(pts)
        @test n >= 3
        # perimeter recovered by summing edges incl. wrap-around
        per = sum(dist(pts[i], pts[mod1(i+1,n)]) for i in 1:n)
        @test per ≈ 2π*5.0 rtol=2e-3
        # first node not duplicated at the end
        @test dist(pts[1], pts[end]) > 1e-6
        # near-uniform spacing ≈ 1
        gaps = [dist(pts[i], pts[mod1(i+1,n)]) for i in 1:n]
        @test isapprox(sum(gaps)/n, 1.0; atol=0.05)
    end

    @testset "curve in 3-D (helix) meshes under size" begin
        helix(t) = (cos(2π*t), sin(2π*t), t)          # one turn, t∈[0,1]
        L = curve_length(helix; nsample=20000)
        @test L ≈ sqrt((2π)^2 + 1.0) rtol=1e-5        # analytic helix length
        pts, _ = mesh_curve(helix, ConstantSize(0.2))
        gaps = [dist(pts[i],pts[i+1]) for i in 1:length(pts)-1]
        @test isapprox(sum(gaps)/length(gaps), 0.2; atol=0.02)
    end
end
