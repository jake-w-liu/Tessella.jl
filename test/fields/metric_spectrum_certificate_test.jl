module MetricSpectrumCertificateChecks
using Test,Tessella
const F=Tessella.SizeField

struct FixedMetric <: F.AbstractAnisoField
    metric::F.Metric3
end
F.metric_at(field::FixedMetric,x,y,z)=field.metric

function analytic_spectrum(tiny,scale)
    setprecision(BigFloat,4096) do
        t=BigFloat(tiny);s=BigFloat(scale)
        discriminant=sqrt((3+t)^2-8t)
        (s*4t/(3+t+discriminant),s*(1-t),s*(3+t+discriminant)/2)
    end
end
eigen_query(metric)=F.metric_eigenvalues(metric)

@testset "Metric spectra certify tiny positive directions" begin
    for exponent in (100,500,1000), power in (-50,0,50), permutation in ((1,2,3),(3,1,2),(2,3,1))
        tiny=ldexp(1.,-exponent);scale=ldexp(1.,power)
        matrix=scale.*[1. tiny 1.;tiny 1. 1.;1. 1. 2.]
        matrix=matrix[collect(permutation),collect(permutation)]
        metric=F.Metric3(matrix[1,1],matrix[2,2],matrix[3,3],
            matrix[1,2],matrix[1,3],matrix[2,3])
        expected=analytic_spectrum(tiny,scale)
        @test eigen_query(metric)==map(Float64,expected)
        coarse=setprecision(BigFloat,4096) do
            Float64(inv(sqrt(expected[1])))
        end
        @test F.field_value(F.MaxSize((FixedMetric(metric),F.ConstantSize(1.))),0.,0.,0.)≈coarse rtol=2e-15
        fine=setprecision(BigFloat,4096) do
            Float64(inv(sqrt(expected[3])))
        end
        @test F.metric_size(metric)≈fine rtol=2e-15
    end
    metric=F.Metric3(1.,1.,2.,2.0^-1000,1.,1.)
    low_values,low_vectors=F._sym3_eigh_big(1.,1.,2.,2.0^-1000,1.,1.;bits=512)
    @test all(>(0),low_values)
    @test !F._metric_spectrum_certified(metric,low_values,low_vectors)
    capped=F.metric_aniso_cap(metric,10.)
    expected=analytic_spectrum(2.0^-1000,1.)
    @test collect(F.metric_eigenvalues(capped))≈collect(Float64.((expected[1],
        100expected[1],100expected[1]))) rtol=2e-13

    # The public Float64 spectrum rounds honestly, while physical size consumers
    # retain the positive direction below its smallest representable eigenvalue.
    underflow=F.Metric3(1.,nextfloat(0.),1.,2e-162,0.,0.)
    coarse=setprecision(BigFloat,4096) do
        R=Rational{BigInt}
        determinant=BigFloat(R(underflow.m11)*R(underflow.m22)-R(underflow.m12)^2)
        trace=BigFloat(R(underflow.m11)+R(underflow.m22))
        tiny=2determinant/(trace+sqrt(trace^2-4determinant))
        @test Float64(tiny)==0.
        Float64(inv(sqrt(tiny)))
    end
    @test eigen_query(underflow)==(0.,1.,1.)
    @test F.field_value(F.MaxSize((FixedMetric(underflow),F.ConstantSize(1.))),0.,0.,0.)≈coarse rtol=2e-15
    math=F.MathEvalAnisoField(m11="1",m22=string(nextfloat(0.)),m33="1",m12="2e-162")
    @test F.metric_at(math,0.,0.,0.)==underflow

    saved_precision=precision(BigFloat);saved_rounding=rounding(BigFloat)
    tasks=[Threads.@spawn begin
        bits=(128,192,256,320)[mod1(index,4)]
        mode=(RoundDown,RoundUp,RoundNearest,RoundToZero)[mod1(index,4)]
        setprecision(BigFloat,bits) do
            setrounding(BigFloat,mode) do
                eigen_query(metric)==map(Float64,expected) &&
                    precision(BigFloat)==bits && rounding(BigFloat)==mode
            end
        end
    end for index in 1:8]
    @test all(fetch,tasks)
    @test precision(BigFloat)==saved_precision
    @test rounding(BigFloat)==saved_rounding
    if Base.JLOptions().compile_enabled==1 && Base.JLOptions().opt_level>=1
        ordinary=F.Metric3(2.,3.,4.,.1,.2,.3)
        eigen_query(ordinary)
        @test @allocated(eigen_query(ordinary))==0
    end
end
end
