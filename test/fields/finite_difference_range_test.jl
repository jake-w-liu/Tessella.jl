module FiniteDifferenceRangeChecks
using Test,Tessella
const F=Tessella.SizeField

struct SampleField{T} <: F.AbstractField
    evaluate::T
end
F.field_value(field::SampleField,x,y,z)=field.evaluate(x,y,z)
query(field)=F.field_value(field,0.,0.,0.)
away_query(field)=F.field_value(field,.3,.4,.5)

@testset "Finite differences retain representable extreme-scale results" begin
    for value in (1e308,-1e308,1.,nextfloat(0.)), delta in (1e-200,1e-4,4.,1e200)
        field=SampleField((x,y,z)->value)
        @test query(F.MeanField(field;delta))==value
        @test query(F.GradientField(field;delta))==0.
        @test query(F.LaplacianField(field;delta))==0.
        @test query(F.MaxEigenHessianField(field;delta))==0.
    end
    # Intermediate differences and sums overflow, although the stencil divided
    # by its actual spacing is finite.
    linear=F.MathEvalField("5e307*x")
    quadratic=F.MathEvalField("1e307*(x^2+y^2+z^2)")
    @test query(F.GradientField(linear;kind=0,delta=4.))≈5e307
    @test query(F.GradientField(linear;kind=3,delta=4.))≈5e307
    @test query(F.LaplacianField(quadratic;delta=4.))≈6e307
    @test query(F.MaxEigenHessianField(quadratic;delta=4.))≈2e307
    for (delta,coefficient,expected) in ((1e-200,1e153,2e306),
                                       (1e-160,1e150,2e300),
                                       (1e200,1e-154,2e-308))
        scaled=SampleField((x,y,z)->(coefficient*x)^2)
        @test query(F.LaplacianField(scaled;delta))≈expected rtol=2e-15
        @test query(F.MaxEigenHessianField(scaled;delta))≈expected rtol=2e-15
    end

    calls=Ref(0)
    values=(1e308,1e308,-1e308,-1e308,1.,2.,4.)
    sampled=SampleField((x,y,z)->(calls[]+=1;values[calls[]]))
    @test query(F.MeanField(sampled))==1.
    @test calls[]==7

    # A monotone function of x+y+z has parallel planar level sets. Its unit
    # gradient remains representable when raw differences or their norm do not.
    planar=F.MathEvalField("1e308*tanh(x+y+z-2)")
    @test query(F.CurvatureField(planar;delta=4.))==0.

    mixed=F.MathEvalField("-1e308*x^2+1e308*x*y")
    expected=setprecision(BigFloat,256) do
        Float64((sqrt(BigFloat(2))-1)*BigFloat(1e308))
    end
    @test query(F.MaxEigenHessianField(mixed;delta=1.))≈expected rtol=2e-15
    concave=F.MathEvalField("-1e308*x^2-1e308*y^2+1e307*x*y")
    @test query(F.MaxEigenHessianField(concave;delta=1.))==0.

    # Real out-of-range derivatives still reject instead of returning a size.
    excessive=F.MathEvalField("1e308*(x^2+y^2+z^2)")
    @test_throws ArgumentError query(F.LaplacianField(excessive;delta=1.))
    @test_throws ArgumentError query(F.MaxEigenHessianField(excessive;delta=1.))

    if Base.JLOptions().compile_enabled==1 && Base.JLOptions().opt_level>=1
        for field in (F.MeanField(F.ConstantSize(1e308)),
                      F.LaplacianField(F.ConstantSize(1e308)),
                      F.GradientField(F.MathEvalField("x+y+z")),
                      F.MaxEigenHessianField(F.MathEvalField("x^2+y^2+z^2")))
            query(field)
            @test @allocated(query(field))==0
        end
        hessian=F.MaxEigenHessianField(F.MathEvalField("x^2+y^2+z^2"))
        away_query(hessian)
        @test @allocated(away_query(hessian))==0
    end
end
end
