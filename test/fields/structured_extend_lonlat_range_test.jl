module StructuredExtendLonLatRangeChecks
using Test,Tessella
const F=Tessella.SizeField
const R=Rational{BigInt}
query(field,point,entity=nothing)=F.field_value(field,point...,entity)

@testset "Structured interpolation retains original finite coordinates" begin
    data=[i+2j+3k for i in 0.:2.,j in 0.:1.,k in 0.:1.]
    field=F.StructuredField((0.,0.,0.),(1.,1.,1.),data)
    @test query(field,(.2,.3,.4))≈2.
    @test query(field,(1e308,-1e308,1e308))==5.
    @test query(field,(2.,1.,1.))==7.
    outside=F.StructuredField((0.,0.,0.),(1.,1.,1.),data;outside=17.)
    @test query(outside,(1e308,0.,0.))==17.
    # Preserve the established upper-endpoint OutsideValue convention.
    @test query(outside,(2.,1.,1.))==17.
    @test query(outside,(1.,.25,.25))≈2.25
    for (origin,delta,point) in ((-1e16,1e16,-.5),(1e16,-1e16,.5))
        near=F.StructuredField((origin,0.,0.),(delta,1.,1.),reshape([0.,1.],2,1,1);outside=17.)
        @test query(near,(point,0.,0.))==1.
        @test query(near,(0.,0.,0.))==17.
        @test query(near,(-point,0.,0.))==17.
    end
    lower=F.StructuredField((0.,0.,0.),(1e308,1.,1.),reshape([0.,1.],2,1,1);outside=17.)
    @test query(lower,(-nextfloat(0.),0.,0.))==17.
    singleton=F.StructuredField((0.,0.,0.),(0.,0.,0.),fill(.2,1,1,1);outside=17.)
    @test query(singleton,(1e308,-1e308,1e308))==.2
    reversed=F.StructuredField((1.,0.,0.),(-.5,1.,1.),reshape([1.,2.,3.],3,1,1))
    @test query(reversed,(.25,0.,0.))==2.5
    @test query(reversed,(1e308,0.,0.))==1.
    @test query(reversed,(-1e308,0.,0.))==3.

    shifted=F.StructuredField((1e16,0.,0.),(3.,1.,1.),reshape([0.,1.,2.],3,1,1))
    @test query(shifted,(1e16+4,0.,0.))==4/3
    wide=F.StructuredField((-1e308,0.,0.),(1e308,1.,1.),reshape([-1e308,0.,1e308],3,1,1))
    @test query(wide,(1e308,0.,0.))==1e308
    @test query(wide,(0.,0.,0.))==0.
    for (delta,point) in ((1e308,nextfloat(0.)),(1e150,1e-170),(1.,nextfloat(0.)))
        linear=F.StructuredField((0.,0.,0.),(delta,1.,1.),reshape([0.,delta],2,1,1))
        @test query(linear,(point,0.,0.))==point
    end
    subnormal=F.StructuredField((0.,0.,0.),(nextfloat(0.),1.,1.),reshape([1.,2.,3.],3,1,1))
    @test query(subnormal,(nextfloat(0.),0.,0.))==2.
    @test query(subnormal,(2nextfloat(0.),0.,0.))==3.
    for value in (floatmax(Float64),-floatmax(Float64))
        constant=F.StructuredField((0.,0.,0.),(1.,1.,1.),fill(value,2,2,2))
        @test query(constant,(.3,.7,.9))==value
    end
    cancelling=F.StructuredField((0.,0.,0.),(1.,1.,1.),
        reshape([floatmax(Float64),-floatmax(Float64),0.,0.,0.,0.,0.,0.],2,2,2))
    @test query(cancelling,(.5,.5,.5))==0.
    saved=copy(field.data)
    query(field,(.2,.3,.4))
    @test field.data==saved
    if Base.JLOptions().compile_enabled==1 && Base.JLOptions().opt_level>=1
        query(field,(.25,.4,.6))
        @test @allocated(query(field,(.25,.4,.6)))==0
    end
end

@testset "Extend scaling retains finite factored sizes" begin
    for (sample,entity_factor,global_factor) in ((1e308,1e308,1e308),
            (1e-308,1e-308,1e-308),(1e-20,1e-160,1e-160),
            (1e-200,1e-100,1e-100),(1.,1.,1.))
        factors=Dict((2,7)=>entity_factor)
        field=F.ExtendField([(0.,0.,0.)],[sample];size_max=.75,
            global_factor,entity_factors=factors)
        exact=R(sample)/(R(entity_factor)*R(global_factor))
        @test query(field,(0.,0.,0.),(2,7))≈Float64(exact) rtol=2e-15
        @test query(field,(.5,0.,0.),(2,7))≈Float64(exact/2+R(.75)/2) rtol=2e-15
        @test query(field,(2.,0.,0.),(2,7))==.75
        factors[(2,7)]=2.
        @test query(field,(0.,0.,0.),(2,7))≈Float64(exact) rtol=2e-15
        restricted=F.RestrictField(field;surfaces=[7])
        @test query(restricted,(0.,0.,0.),(2,7))≈Float64(exact) rtol=2e-15
    end
    excessive=F.ExtendField([(0.,0.,0.)],[1e308];global_factor=nextfloat(0.))
    @test_throws ArgumentError query(excessive,(0.,0.,0.),(2,7))
    if Base.JLOptions().compile_enabled==1 && Base.JLOptions().opt_level>=1
        ordinary=F.ExtendField([(0.,0.,0.)],[.2])
        query(ordinary,(.5,0.,0.),(2,7))
        @test @allocated(query(ordinary,(.5,0.,0.),(2,7)))==0
    end
end

@testset "Stereographic LonLat avoids physical squared-scale overflow" begin
    input=F.MathEvalField("1+x+2*y")
    for radius in (nextfloat(0.),1e-200,1.,1e200,1e308)
        constant=F.LonLatField(F.ConstantSize(.2);from_stereo=true,radius)
        @test query(constant,(0.,0.,0.))==.2
        field=F.LonLatField(input;from_stereo=true,radius)
        @test query(field,(0.,0.,0.))≈1+pi rtol=2e-15
        @test query(field,(radius,0.,0.))≈1+2asin(.6) rtol=2e-15
        @test query(field,(radius,-radius,0.))≈1-pi/4+2asin(1/3) rtol=2e-15
    end
    ordinary=F.LonLatField(input;radius=1.)
    @test query(ordinary,(.3,.4,.5))≈1+atan(.4,.3)+2asin(.5)
    @test query(ordinary,(0.,0.,2.))==F.GMSH_MAX_SIZE
    if Base.JLOptions().compile_enabled==1 && Base.JLOptions().opt_level>=1
        ordinary_stereo=F.LonLatField(input;from_stereo=true)
        query(ordinary_stereo,(1.,2.,0.))
        @test @allocated(query(ordinary_stereo,(1.,2.,0.)))==0
    end
end
end
