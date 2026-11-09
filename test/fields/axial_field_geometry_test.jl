module AxialFieldGeometryChecks
using Test,Tessella
const F=Tessella.SizeField

function thin_frustum(a,b,r)
    return F.FrustumField(a,b;inner_r1=0.,outer_r1=r,inner_r2=0.,outer_r2=r,
        inner_v1=1.,outer_v1=1.,inner_v2=3.,outer_v2=3.)
end

@testset "Cylinder predicates preserve the original finite axis" begin
    for axis in ((1.,1.,0.),(1.,2.,3.),(-1.,2.,-3.))
        field=F.CylinderField((0.,0.,0.),axis,1e-20;vin=2.,vout=9.)
        @test F.field_value(field,map(x->x/2,axis)...)==2.
        @test F.field_value(field,axis...)==9.
        @test F.field_value(field,map(-,axis)...)==9.
    end
    field=F.CylinderField((0.,0.,0.),(1.,1.,0.),1e-20;vin=2.,vout=9.)
    @test F.field_value(field,.5,.5,prevfloat(1e-20))==2.
    @test F.field_value(field,.5,.5,1e-20)==9.
    @test F.field_value(field,.5,.5,nextfloat(1e-20))==9.
    @test F.field_value(field,prevfloat(1.),prevfloat(1.),0.)==2.
    @test F.field_value(field,nextfloat(1.),nextfloat(1.),0.)==9.
    translated=F.CylinderField((1e16,-1e16,1e16),(4.,4.,0.),1e-20;vin=2.,vout=9.)
    @test F.field_value(translated,1e16+2,-1e16+2,1e16)==2.
    @test F.field_value(translated,1e16+4,-1e16+4,1e16)==9.
    legacy=F.CylinderField((0.,0.,0.),(0.,0.,1.),2.,1.,2.,9.)
    @test F.field_value(legacy,.5,0.,1.)==2.
    @test F.field_value(legacy,0.,0.,2.)==9.
end

@testset "Subnormal radial predicates retain their strict or closed boundary" begin
    for h in (nextfloat(0.),ldexp(1.,-1040),1e-200,1e-100)
        cylinder=F.CylinderField((0.,0.,0.),(0.,0.,h),3h;vin=2.,vout=9.)
        @test F.field_value(cylinder,2h,2h,0.)==2.
        @test F.field_value(cylinder,3h,0.,0.)==9.
        @test F.field_value(cylinder,0.,0.,h)==9.
        frustum=thin_frustum((0.,0.,0.),(0.,0.,1.),h)
        @test F.field_value(frustum,h,h,.5)==F.GMSH_MAX_SIZE
        @test F.field_value(frustum,h,0.,.5)==2.
        @test F.field_value(frustum,0.,0.,.5)==2.
    end
end

@testset "Frustum uses original endpoints and exact annular coordinates" begin
    field=thin_frustum((0.,0.,0.),(1.,1.,0.),1e-20)
    @test F.field_value(field,.5,.5,0.)==2.
    @test F.field_value(field,0.,0.,0.)==1.
    @test F.field_value(field,1.,1.,0.)==3.
    @test F.field_value(field,.5,.5,1e-20)==2.
    @test F.field_value(field,.5,.5,2e-20)==F.GMSH_MAX_SIZE
    translated=thin_frustum((1e16,-1e16,1e16),(1e16+4,-1e16+4,1e16),1e-20)
    @test F.field_value(translated,1e16+2,-1e16+2,1e16)==2.
    tapered=F.FrustumField((0.,0.,0.),(0.,0.,1.);inner_r1=1.,outer_r1=1.,
        inner_r2=1.,outer_r2=2.,inner_v1=2.,outer_v1=2.,inner_v2=2.,outer_v2=2.)
    @test F.field_value(tapered,1.,0.,1e-40)==2.
    @test F.field_value(tapered,1.,0.,0.)==F.GMSH_MAX_SIZE
    @test F.field_value(tapered,prevfloat(1.),0.,1e-40)==F.GMSH_MAX_SIZE
    legacy=F.FrustumField((1e16,0.,0.),(1.,0.,0.),1.,0.,1.,0.,1.,1.,1.,3.,3.)
    @test F.field_value(legacy,1e16,0.,0.)==1.
    @test F.field_value(legacy,nextfloat(1e16),0.,0.)==F.GMSH_MAX_SIZE
end

@testset "Finite endpoint geometry need not fit an intermediate length" begin
    axis=(1.2e308,1.2e308,1.2e308)
    cylinder=F.CylinderField((0.,0.,0.),axis,1.;vin=2.,vout=9.)
    @test all(isfinite,cylinder.unit_axis)
    @test F.field_value(cylinder,0.,0.,0.)==2.
    @test F.field_value(cylinder,axis...)==9.
    @test F.field_value(cylinder,1.,-1.,0.)==9.
    frustum=thin_frustum((-1e308,0.,0.),(1e308,0.,0.),1.)
    @test F.field_value(frustum,-1e308,0.,0.)==1.
    @test F.field_value(frustum,0.,0.,0.)==2.
    @test F.field_value(frustum,1e308,0.,0.)==3.
    @test F.field_value(frustum,0.,2.,0.)==F.GMSH_MAX_SIZE
end

@testset "Frustum retains a finite value below its axial parameter grid" begin
    field=F.FrustumField((0.,0.,0.),(1e308,0.,0.);inner_r1=-1.,inner_r2=-1.,
        outer_r1=1.,outer_r2=1.,inner_v1=0.,outer_v1=0.,inner_v2=1e308,outer_v2=1e308)
    for x in (1e-200,1e-100,1e-12)
        @test F.field_value(field,x,0.,0.)≈x rtol=8eps(Float64)
    end
    @test F.field_value(field,-1e-200,0.,0.)==F.GMSH_MAX_SIZE
    @test F.field_value(field,0.,0.,0.)==0.
    @test_throws ArgumentError F.size_at(field,0.,0.,0.)
end
end
