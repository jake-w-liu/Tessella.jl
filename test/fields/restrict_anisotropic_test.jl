module RestrictAnisotropicChecks
using Test,Tessella
const F=Tessella.SizeField

@testset "Anisotropic restrictions preserve metrics and scalar operators" begin
    input=F.MathEvalAnisoField(;m11="100",m22="1",m33="1")
    field=F.RestrictField(input;surfaces=[7])
    @test field isa F.AbstractAnisoField
    # MathEvalAniso's scalar operator returns M11, not a mesh length.
    @test F.field_value(field,0.,0.,0.,(2,7))==100.
    @test F.metric_at(field,0.,0.,0.,(2,7)).m11==100.
    @test F.metric_at(field,0.,0.,0.,(2,7)).m22==1.
    @test F.directional_size(field,(0.,0.,0.),(1.,0.,0.);entity=(2,7))≈.1
    @test F.directional_size(field,(0.,0.,0.),(0.,1.,0.);entity=(2,7))≈1.
    @test F.field_value(field,0.,0.,0.,(2,8))==F.GMSH_MAX_SIZE
    @test F.directional_size(field,(0.,0.,0.),(0.,1.,0.);entity=(2,8))≈F.GMSH_MAX_SIZE
    @test F.directional_size(field,(0.,0.,0.),(0.,1.,0.))≈1.
    @test F.field_value(F.MinSize((field,F.ConstantSize(.5))),0.,0.,0.,(2,7))≈.1
    @test F.field_value(F.MaxSize((field,F.ConstantSize(.5))),0.,0.,0.,(2,7))≈1.
    bounded=F.BoundedSize(field;size_min=.2,size_max=.8)
    @test bounded isa F.AbstractAnisoField
    @test F.directional_size(bounded,(0.,0.,0.),(1.,0.,0.);entity=(2,7))≈.2
    @test F.directional_size(bounded,(0.,0.,0.),(0.,1.,0.);entity=(2,7))≈.8
    nested=F.RestrictField(field;surfaces=[7])
    @test F.metric_at(nested,0.,0.,0.,(2,7)).m22==1.
    scalar=F.RestrictField(F.ConstantSize(.25);surfaces=[7])
    @test !(scalar isa F.AbstractAnisoField)
    @test F.size_at(scalar,0.,0.,0.,(2,7))==.25
end

@testset "Anisotropic restriction expands boundary and embedded entities" begin
    input=F.MathEvalAnisoField(;m11="100",m22="1",m33="1")
    boundaries=Dict((3,9)=>[(2,7)],(2,7)=>[(1,3)],(1,3)=>[(0,2)])
    field=F.RestrictField(input;volumes=[9],entity_boundaries=boundaries)
    for entity in ((3,9),(2,7),(1,3),(0,2))
        @test F.metric_at(field,0.,0.,0.,entity).m22==1.
        @test F.field_value(field,0.,0.,0.,entity)==100.
    end
    embedded=F.RestrictField(input;volumes=[9],include_boundary=false,
        entity_boundaries=boundaries,entity_embedded=Dict((3,9)=>[(2,7)]))
    @test F.metric_at(embedded,0.,0.,0.,(0,2)).m22==1.
    excluded=F.RestrictField(input;volumes=[9],include_boundary=false,include_embedded=false,
        entity_boundaries=boundaries,entity_embedded=Dict((3,9)=>[(2,7)]))
    @test F.field_value(excluded,0.,0.,0.,(2,7))==F.GMSH_MAX_SIZE
end

struct EntityMetric <: F.AbstractAnisoField end
F.metric_at(::EntityMetric,x,y,z)=F.Metric3(4.,9.,16.,0.,0.,0.)
F.metric_at(::EntityMetric,x,y,z,entity::Tuple{T,U}) where {T<:Integer,U<:Integer}=
    F.Metric3(Float64(entity[2]),1.,1.,0.,0.,0.)

@testset "Restriction forwards metric entity context" begin
    field=F.RestrictField(EntityMetric();surfaces=[7])
    @test F.metric_at(field,0.,0.,0.).m11==4.
    @test F.metric_at(field,0.,0.,0.,(2,7)).m11==7.
end
end
