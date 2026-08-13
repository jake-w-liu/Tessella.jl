# ── Stage-2 CRC suite: size fields + 1-D graded edge meshing ────────────────────
#
# Correctness  : arc length vs analytic (segment, circle); edge count = round(L/h);
#                uniform spacing under constant size; graded spacing tracks h.
# Robustness   : varying size fields, closed curves, tiny/large sizes.
# Completeness : MinSize combination, endpoint inclusion, closed-loop non-duplication.

using Test
using Tessella.SizeField
using Tessella.Mesh1D

dist(a,b) = sqrt(sum((a[i]-b[i])^2 for i in 1:3))

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
        @test_throws ArgumentError curve_length(seg; t0=1, t1=0)
        @test_throws ArgumentError curve_length(t->(NaN,0.0,0.0); nsample=2)
        @test_throws ArgumentError metric_length(seg, ConstantSize(1.0); nsample=-1)
        @test_throws ArgumentError mesh_curve(t->(0.0,0.0,0.0), ConstantSize(1.0))
        @test_throws ArgumentError mesh_curve(seg, ConstantSize(1.0); closed=true)
        @test_throws ArgumentError mesh_segment((0.0,0.0,0.0), (1.0,0.0,0.0),
                                                ConstantSize(1e-300); nsample=1)
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
