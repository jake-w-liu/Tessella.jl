module ModelVolumeBoundaryCacheChecks
using Test, Tessella
const M=Tessella.Model
const F=Tessella.SizeField

function cache_checksum(field,points,count)
    value=0.0
    for _ in 1:count,point in points
        value+=F.size_at(field,point...)
    end
    return value
end

@testset "boundary-size reuse preserves the owned view operator" begin
    view=PostViewField([0. 1. 0.;0. 0. 1.;0. 0. 0.],[.1,.2,.4];
        triangles=reshape(Int32[1,2,3],3,1),crop_negative=false)
    cache=M._VolumeBoundarySizeCache(
        view,Dict{NTuple{3,Float64},Float64}(),ReentrantLock())
    field=FunctionSize(cache)
    points=[(0.,0.,0.),(1.,0.,0.),(.25,.25,0.),(.25,.25,.001),
            (.25,.25,1.),(-2.,0.,0.),(2.,2.,2.)]
    expected=[F.field_value(view,p...) for p in points]
    @test [F.size_at(field,p...) for p in points]==expected
    @test [F.size_at(field,p...) for p in points]==expected
    @test F.size_at(field,0,0,0)==F.field_value(view,0,0,0)
    @test F.size_at(field,.25f0,.25f0,0f0)==F.field_value(view,.25f0,.25f0,0f0)
    @test_throws ArgumentError F.size_at(field,true,0.,0.)
    @test_throws ArgumentError F.size_at(field,NaN,0.,0.)

    # Exhaust more than one scratch generation, then revisit old queries.
    # Neither eviction nor interleaved readers may alter a field value.
    @test all(1:4200) do i
        p=(2.0+i/1000,2.,2.)
        F.size_at(field,p...)==F.field_value(view,p...)
    end
    @test length(cache.values)<=4096
    @test [F.size_at(field,p...) for p in points]==expected
    actual=Vector{Vector{Float64}}(undef,8)
    @sync for task in eachindex(actual)
        Threads.@spawn actual[task]=[F.size_at(field,p...) for _ in 1:8 for p in points]
    end
    @test all(values==repeat(expected,8) for values in actual)
    cache_checksum(field,points,2)
    if Base.JLOptions().compile_enabled!=3
        @test (@allocated cache_checksum(field,points,100))==0
    end
end

@testset "volume boundary extension installs its private reuse cache" begin
    model=GeoModel()
    M.add_box!(model,0.,0.,0.,1.,1.,1.;tag=1)
    surface=Tessella.Geometry.box_surface(0.,1.,0.,1.,0.,1.)
    field=M._volume_boundary_size_field(model,1,surface,"boundary-size reuse")
    @test field.f isa M._VolumeBoundarySizeCache
    @test isfinite(F.size_at(field,.5,.5,.5)) && F.size_at(field,.5,.5,.5)>0
    model.meshing.lc_extend_from_boundary=0
    @test M._volume_boundary_size_field(model,1,surface,"boundary-size reuse")===nothing
end
end
