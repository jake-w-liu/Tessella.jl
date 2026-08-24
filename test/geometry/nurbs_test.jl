using Test
using Tessella

@testset "NURBS De Boor vs Bernstein and circle" begin
    # Quadratic Bézier = NURBS with knots [0,0,0,1,1,1].
    P=[(0.0,0.0,0.0),(1.0,1.0,0.0),(2.0,0.0,0.0)]
    c=NURBSCurve(2,[0,0,0,1,1,1],P)
    bernstein(t)=((1-t)^2 .* P[1] .+ 2(1-t)*t .* P[2] .+ t^2 .* P[3])
    for t in (0.0,0.25,0.5,0.75,1.0)
        q=nurbs_eval(c,t); b=bernstein(t)
        @test hypot(q[1]-b[1],q[2]-b[2],q[3]-b[3])<=1e-14
    end
    # Independent Cox–de Boor partition of unity on the same knot vector.
    U=Float64[0,0,0,1,1,1]
    for t in (0.1,0.5,0.9)
        s=sum(bspline_basis(i,2,U,t) for i in 1:3)
        @test s≈1.0 atol=1e-14
    end
    for t in (0.0,0.1,0.5,0.9,1.0)
        @test bspline_basis(1,2,U,t)≈(1-t)^2 atol=1e-14
        @test bspline_basis(2,2,U,t)≈2t*(1-t) atol=1e-14
        @test bspline_basis(3,2,U,t)≈t^2 atol=1e-14
    end
    @test bspline_basis(1,2,U,-1.0)==0.0
    @test bspline_basis(3,2,U,2.0)==0.0
    nonclamped=collect(0.0:6.0)
    for t in (2.0,2.5,3.5,4.0)
        @test sum(bspline_basis(i,2,nonclamped,t) for i in 1:4)≈1.0 atol=1e-14
    end
    @test bspline_basis(1,0,[0.0,1.0,2.0],0.5)==1.0
    @test bspline_basis(2,0,[0.0,1.0,2.0],1.5)==1.0
    @test bspline_basis(2,0,[0.0,1.0,2.0],2.0)==1.0
    @test nurbs_eval(c,-1.0)==P[1]
    @test nurbs_eval(c,2.0)==P[3]
    @test_throws ArgumentError NURBSCurve(2,[0,1],P)
    @test_throws ArgumentError NURBSCurve(true,[0,0,1,1],P[1:2])
    @test_throws ArgumentError NURBSCurve(1.0,[0,0,1,1],P[1:2])
    @test_throws ArgumentError NURBSCurve(big(2)^100,[0,0,1,1],P[1:2])
    @test_throws ArgumentError NURBSCurve(1,[0,0,0,0],P[1:2])
    @test_throws ArgumentError bspline_basis(1,-1,U,0.5)
    @test_throws ArgumentError bspline_basis(1,true,U,0.5)
    @test_throws ArgumentError bspline_basis(true,2,U,0.5)
    @test_throws ArgumentError bspline_basis(1,2,[0,1,0,1,1,1],0.5)
    @test_throws ArgumentError bspline_basis(1,2,[0,0,0,Inf,1,1],0.5)
    @test_throws ArgumentError bspline_basis(1,2,U,NaN)

    # Homogeneous weights are globally scale invariant. Normalizing them before
    # multiplication keeps a mathematically finite result finite at Float64 limits.
    xmax=floatmax(Float64)
    extreme=NURBSCurve(1,[0,0,1,1],[(xmax,0.0,0.0),(xmax,0.0,0.0)],[2.0,2.0])
    @test nurbs_eval(extreme,0.5)==(xmax,0.0,0.0)
    unrepresentable=NURBSCurve(1,[0,0,1,1],P[1:2],[nextfloat(0.0),xmax])
    @test_throws ArgumentError nurbs_eval(unrepresentable,0.5)

    # Quadratic NURBS circle of radius 1 in the xy-plane (9 controls, standard
    # square-to-circle weights). Evaluate at the four axis points.
    w=1/sqrt(2)
    ctrls=[(1.0,0.0,0.0),(1.0,1.0,0.0),(0.0,1.0,0.0),(-1.0,1.0,0.0),
           (-1.0,0.0,0.0),(-1.0,-1.0,0.0),(0.0,-1.0,0.0),(1.0,-1.0,0.0),(1.0,0.0,0.0)]
    knots=[0,0,0,1,1,2,2,3,3,4,4,4]./4
    weights=[1,w,1,w,1,w,1,w,1]
    circle=NURBSCurve(2,knots,ctrls,weights)
    for (u,target) in ((0.0,(1.0,0.0,0.0)),(0.25,(0.0,1.0,0.0)),
                       (0.5,(-1.0,0.0,0.0)),(0.75,(0.0,-1.0,0.0)),(1.0,(1.0,0.0,0.0)))
        q=nurbs_eval(circle,u)
        @test hypot(q[1]-target[1],q[2]-target[2],q[3]-target[3])<=1e-12
        @test hypot(q[1],q[2],q[3])≈1.0 atol=1e-12
    end

    # Bilinear NURBS patch equals bilinear interpolation of four corners.
    S=NURBSSurface(1,1,[0,0,1,1],[0,0,1,1],
                   [(0.0,0.0,0.0) (0.0,1.0,0.0); (1.0,0.0,0.0) (1.0,1.0,2.0)])
    q=nurbs_eval(S,0.5,0.5)
    @test q[1]≈0.5 && q[2]≈0.5 && q[3]≈0.5
    @test_throws ArgumentError NURBSSurface(true,1,[0,0,1,1],[0,0,1,1],S.controls)
    @test_throws ArgumentError NURBSSurface(1,1,[0,0,0,0],[0,0,1,1],S.controls)
end
