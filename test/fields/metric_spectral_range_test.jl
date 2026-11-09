module MetricSpectralRangeChecks
using Test,Tessella
const F=Tessella.SizeField
const R=Rational{BigInt}

struct FixedMetric <: F.AbstractAnisoField
    metric::F.Metric3
end
F.metric_at(field::FixedMetric,x,y,z)=field.metric

physical_size(value)=setprecision(BigFloat,256) do
    Float64(inv(sqrt(BigFloat(value))))
end

@testset "Finite tensors retain sizes beyond the eigenvalue storage range" begin
    for a in (1e-200,1.,1e200,1e308)
        b=.9a;metric=F.Metric3(a,a,a,b,b,b)
        fine=physical_size(R(a)+2R(b));coarse=physical_size(R(a)-R(b))
        @test F.metric_size(metric)≈fine rtol=8eps(Float64)
        @test F.field_value(F.MaxSize((FixedMetric(metric),F.ConstantSize(1e-300))),0.,0.,0.)≈coarse rtol=8eps(Float64)
    end
    a=1e308;b=9e307;metric=F.Metric3(a,a,a,b,b,b);field=FixedMetric(metric)
    expected=physical_size(R(a)+2R(b))
    @test F.field_value(F.MinSize((field,F.ConstantSize(1.))),0.,0.,0.)≈expected rtol=8eps(Float64)
    @test F.field_value(F.MinAnisoField((field,)),0.,0.,0.)≈expected rtol=8eps(Float64)
    @test F.field_value(F.IntersectAnisoField((field,)),0.,0.,0.)≈expected rtol=8eps(Float64)
    @test F.metric_at(F.BoundedSize(field),0.,0.,0.)==metric
    @test F.metric_at(F.BoundedSize(field;factor=2.),0.,0.,0.).m11==a/4
    @test F.metric_at(F.BoundedSize(field;factor=2.),0.,0.,0.).m12==b/4
    bounded=F.BoundedSize(field;size_max=2e-154)
    @test F.directional_size(bounded,(0.,0.,0.),(1.,-1.,0.))≈2e-154 rtol=8eps(Float64)
    @test F.directional_size(bounded,(0.,0.,0.),(1.,1.,1.))≈expected rtol=8eps(Float64)
    coarse_bounded=F.BoundedSize(field;size_min=1e-154)
    @test F.directional_size(coarse_bounded,(0.,0.,0.),(1.,1.,1.))≈1e-154 rtol=8eps(Float64)
    smallest=R(a)-R(b);cap=smallest*R(4.5)^2
    capped=F.metric_aniso_cap(metric,4.5)
    @test capped.m11≈Float64((cap+2smallest)/3) rtol=8eps(Float64)
    @test F.directional_size(FixedMetric(capped),(0.,0.,0.),(1.,1.,1.))≈physical_size(cap) rtol=8eps(Float64)
    @test F.directional_size(FixedMetric(capped),(0.,0.,0.),(1.,-1.,0.))≈physical_size(smallest) rtol=8eps(Float64)
    @test_throws ArgumentError F.Metric3(a,a,a,1.1e308,1.1e308,1.1e308)
end

@testset "Physical sizes retain a positive eigenvalue below the Float64 grid" begin
    a=1e308;c=nextfloat(1e-308);metric=F.Metric3(a,c,1.,1.,0.,0.)
    determinant=R(a)*R(c)-1
    expected=setprecision(BigFloat,256) do
        trace=R(a)+R(c)
        largest=(BigFloat(trace)+sqrt(BigFloat(trace^2-4determinant)))/2
        Float64(sqrt(largest/BigFloat(determinant)))
    end
    field=FixedMetric(metric)
    @test F.field_value(F.MaxSize((field,F.ConstantSize(1.))),0.,0.,0.)≈expected rtol=8eps(Float64)
    @test F.metric_size(metric)≈1e-154 rtol=8eps(Float64)
    @test F.directional_size(F.BoundedSize(field),(0.,0.,0.),(-1e-308,1.,0.))≈F.GMSH_MAX_SIZE rtol=8eps(Float64)
    diagonal=FixedMetric(F.Metric3(nextfloat(0.),1.,1.,0.,0.,0.))
    @test F.field_value(F.MaxSize((diagonal,F.ConstantSize(1.))),0.,0.,0.)≈physical_size(R(nextfloat(0.))) rtol=8eps(Float64)
end
end
