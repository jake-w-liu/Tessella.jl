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
    @test_throws ArgumentError NURBSCurve(2,[0,1],P)

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
end
