module ScalarFiniteInterpolationChecks
using Test,Tessella
const F=Tessella.SizeField

@testset "Finite threshold intervals and scalar interpolation" begin
    input=F.MathEvalField("x")
    field=F.ThresholdField(input;dist_min=-1e308,dist_max=1e308,size_min=1.,size_max=3.)
    for (x,expected) in ((-1e308,1.),(0.,2.),(1e308,3.))
        @test F.field_value(field,x,0.,0.)≈expected
    end
    sigmoid=F.ThresholdField(input;dist_min=-1e308,dist_max=1e308,
        size_min=1.,size_max=3.,sigmoid=true)
    @test F.field_value(sigmoid,0.,0.,0.)≈2.
    stopped=F.ThresholdField(input;dist_min=-1e308,dist_max=1e308,
        size_min=1.,size_max=3.,stop_at_dist_max=true)
    @test F.field_value(stopped,1e308,0.,0.)==F.GMSH_MAX_SIZE
    values=F.ThresholdField(input;dist_min=0.,dist_max=2.,size_min=-1e308,size_max=1e308)
    for (x,expected) in ((0.,-1e308),(1.,0.),(2.,1e308))
        @test F.field_value(values,x,0.,0.)==expected
    end
    @test_throws ArgumentError F.size_at(values,1.,0.,0.)
end

@testset "Finite primitive transition values" begin
    box=F.BoxField(0.,1.,0.,1.,0.,1.;vin=-1e308,vout=1e308,thickness=1.)
    @test F.field_value(box,1.5,.5,.5)==0.
    ball=F.BallField((0.,0.,0.),1.;vin=-1e308,vout=1e308,thickness=1.)
    @test F.field_value(ball,1.5,0.,0.)==0.
    frustum=F.FrustumField((0.,0.,0.),(0.,0.,1.);inner_v1=-1e308,outer_v1=-1e308,
        inner_v2=1e308,outer_v2=1e308)
    @test F.field_value(frustum,.5,0.,.5)==0.
    radii=F.FrustumField((0.,0.,0.),(0.,0.,1.);inner_r1=-1e308,outer_r1=1e308,
        inner_r2=-1e308,outer_r2=1e308,inner_v1=1.,inner_v2=1.,outer_v1=3.,outer_v2=3.)
    @test F.field_value(radii,0.,0.,.5)≈2.
end

@testset "Ball transitions retain physical displacement" begin
    large=F.BallField((-1e308,0.,0.),1e308;vin=1.,vout=4.,thickness=1.5e308)
    @test F.field_value(large,1e308,0.,0.)≈3.
    thin=F.BallField((0.,0.,0.),1.;vin=1.,vout=3.,thickness=1e-20)
    # sqrt(1+1e-20) rounds to1, but its shell coordinate is approximately1/2.
    @test F.field_value(thin,1.,1e-10,0.)≈2. rtol=1e-8
    @test F.field_value(thin,1.,0.,0.)==1.
    @test F.field_value(thin,1.,2e-10,0.)==3.
    sharp=F.BallField((0.,0.,0.),1.;vin=1.,vout=3.)
    @test F.field_value(sharp,1.,0.,0.)==3.
end

@testset "Subnormal corner distances retain a finite transition" begin
    for h in (nextfloat(0.),ldexp(1.,-1040),1e-200,1e-100,1.)
        box=F.BoxField(0.,0.,0.,0.,0.,0.;vin=1.,vout=3.,thickness=2h)
        ball=F.BallField((0.,0.,0.),0.;vin=1.,vout=3.,thickness=2h)
        @test F.field_value(box,h,h,0.)≈1+sqrt(2.) rtol=8eps(Float64)
        @test F.field_value(ball,h,h,0.)≈1+sqrt(2.) rtol=8eps(Float64)
    end
end

@testset "Finite interpolation retains small endpoint and midpoint values" begin
    R=Rational{BigInt}
    for (a,b,t) in ((-1e16,1e16+2,.5),(1e16,1.,prevfloat(1.)),
                    (1.,1e16,eps()/2),(-1e16,-1.,prevfloat(1.)))
        field=F.ThresholdField(F.MathEvalField("x");dist_min=0.,dist_max=1.,size_min=a,size_max=b)
        expected=Float64(R(a)+R(t)*(R(b)-R(a)))
        @test F.field_value(field,t,0.,0.)≈expected rtol=8eps(Float64) atol=0.
    end
    box=F.BoxField(0.,0.,0.,0.,0.,0.;vin=1e16,vout=1.,thickness=1.)
    t=prevfloat(1.);expected=Float64(R(1e16)+R(t)*(R(1.)-R(1e16)))
    @test F.field_value(box,t,0.,0.)≈expected rtol=8eps(Float64)
    ball=F.BallField((0.,0.,0.),1.;vin=-1e16,vout=1e16+2,thickness=1.)
    @test F.field_value(ball,1.5,0.,0.)==1.
    frustum=F.FrustumField((0.,0.,0.),(0.,0.,1.);inner_v1=-1e16,outer_v1=-1e16,
        inner_v2=1e16+2,outer_v2=1e16+2)
    @test F.field_value(frustum,.5,0.,.5)==1.
end
end
