module ThresholdIntervalPrecisionTests
using Test,Tessella,Random
const F=Tessella.SizeField
const A=Tessella.API
const R=Rational{BigInt}

function expected_interval(d,l,u,a,b)
    ratio=clamp((R(d)-R(l))/(R(u)-R(l)),R(0),R(1))
    return Float64(R(a)+ratio*(R(b)-R(a)))
end

@testset "Threshold interpolation retains original finite interval operands" begin
    cases=((nextfloat(0.),0.,1e100,nextfloat(0.),1e308),
        (1e-300,0.,1e20,nextfloat(0.),1e308),
        (prevfloat(3.),-2.,3.,1e308,1.),(prevfloat(3.),0.,3.,1e308,1.),
        (1.,0.,3.,-1.,nextfloat(2.)))
    for (d,l,u,a,b) in cases
        field=F.ThresholdField(F.ConstantSize(d);dist_min=l,dist_max=u,size_min=a,size_max=b)
        expected=expected_interval(d,l,u,a,b)
        @test F.field_value(field,0.,0.,0.)≈expected rtol=8eps(Float64) atol=0.
        @test F.field_value(field,0.,0.,0.,(2,17))≈expected rtol=8eps(Float64) atol=0.
        @test F.size_at(field,0.,0.,0.)≈expected rtol=8eps(Float64) atol=0.
    end
    # The third/fourth controls retain a tiny interval complement after
    # normalizing near its upper endpoint; the fifth cancels to eps(2)/3.
    @test F.field_value(F.ThresholdField(F.ConstantSize(1.);dist_min=0.,dist_max=3.,
        size_min=-1.,size_max=nextfloat(2.)),0.,0.,0.)≈eps(2.)/3 rtol=8eps(Float64)
    for (d,l,u,a,b) in ((.5,0.,1.,.1,.5),(.9,.1,.9,1.,.1),(.1,.1,.9,1.,.1),
            (1.,1.,3.,.1,.5),(3.,1.,3.,.1,.5),(4.,1.,3.,.1,.5))
        field=F.ThresholdField(F.ConstantSize(d);dist_min=l,dist_max=u,size_min=a,size_max=b)
        @test F.field_value(field,0.,0.,0.)≈expected_interval(d,l,u,a,b) rtol=8eps(Float64)
    end
    stopped=F.ThresholdField(F.ConstantSize(3.);dist_min=1.,dist_max=3.,
        size_min=.1,size_max=.5,stop_at_dist_max=true)
    @test F.field_value(stopped,0.,0.,0.)==F.GMSH_MAX_SIZE
    for d in (.5,1.,2.,3.,4.)
        field=F.ThresholdField(F.ConstantSize(d);dist_min=1.,dist_max=3.,
            size_min=.1,size_max=.5,sigmoid=true)
        r=clamp((d-1)/2,0.,1.);t=exp(12r-6)/(1+exp(12r-6))
        @test F.field_value(field,0.,0.,0.)≈.1*(1-t)+.5t rtol=8eps(Float64)
    end
end

@testset "Threshold finite-scale controls use an exact independent reference" begin
    rng=MersenneTwister(492145)
    for trial in 1:256
        scale=ldexp(1.,rand(rng,-1000:1000));lower=scale;upper=3scale
        d=rand(rng,(prevfloat(lower),lower,nextfloat(lower),muladd(2scale,.1,lower),
            muladd(2scale,.5,lower),prevfloat(upper),upper,nextfloat(upper)))
        a=ldexp(.5+.5rand(rng),rand(rng,-1074:1023))*(rand(rng,Bool) ? 1. : -1.)
        b=ldexp(.5+.5rand(rng),rand(rng,-1074:1023))*(rand(rng,Bool) ? 1. : -1.)
        field=F.ThresholdField(F.ConstantSize(d);dist_min=lower,dist_max=upper,size_min=a,size_max=b)
        @test F.field_value(field,0.,0.,0.)≈expected_interval(d,lower,upper,a,b) rtol=8eps(Float64) atol=0.
    end
end

@testset "Distance view scaling retains small finite contributions" begin
    for far in (1e20,1e100)
        A.initialize()
        try
            A.option("General.Terminal",0)
            mktempdir() do dir
                path=joinpath(dir,"scaled.msh")
                open(path,"w") do io
                    println(io,raw"$MeshFormat");println(io,"2.2 0 8");println(io,raw"$EndMeshFormat")
                    println(io,raw"$Nodes");println(io,4)
                    println(io,"11 0 0 0\n21 1 0 0\n31 0.5 1e-300 0")
                    println(io,"41 0.5 ",-far," 0")
                    println(io,raw"$EndNodes");println(io,raw"$Elements");println(io,3)
                    println(io,"201 1 2 0 1 11 21\n301 2 2 0 7 11 21 31\n302 2 2 0 8 21 11 41")
                    println(io,raw"$EndElements")
                end
                A.open(path)
                before=(A.mesh.get_nodes(),A.mesh.get_elements())
                A.plugin.set_number("Distance","MinScale",1e-300)
                A.plugin.set_number("Distance","MaxScale",1e308)
                A.plugin.run("Distance")
                types,data=A.view.get_list_data(last(A.view.get_tags()))
                @test types==["ST"] && length(only(data))==24
                @test only(data)[12]≈expected_interval(1e-300,0.,far,1e-300,1e308) rtol=8eps(Float64)
                @test only(data)[[10,11,22,23]]==fill(1e-300,4)
                @test only(data)[24]==1e308
                @test (A.mesh.get_nodes(),A.mesh.get_elements())==before
            end
        finally
            A.finalize()
        end
    end
end

@noinline function sample_sum(field,n)
    result=0.
    for index in 1:n
        result+=F.field_value(field,Float64(index%37)/37,0.,0.)
    end
    return result
end

@testset "Conditioned interval ramps retain allocation-free ordinary queries" begin
    if Base.JLOptions().compile_enabled!=3
        fields=(F.ThresholdField(F.FunctionSize((x,y,z)->1+x/10);dist_min=1.,dist_max=3.,size_min=.1,size_max=.5),
            F.ThresholdField(F.FunctionSize((x,y,z)->1+2x);dist_min=1.,dist_max=3.,size_min=1.,size_max=.1),
            F.ThresholdField(F.FunctionSize((x,y,z)->.1+.8x);dist_min=.1,dist_max=.9,size_min=1.,size_max=.1),
            F.ThresholdField(F.ConstantSize(2.);dist_min=1.,dist_max=3.,size_min=nextfloat(0.),size_max=nextfloat(0.)))
        for field in fields
            sample_sum(field,1000)
            @test (@allocated sample_sum(field,1000))==0
        end
    end
end
end
