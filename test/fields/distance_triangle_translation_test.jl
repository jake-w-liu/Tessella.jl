using Test
using Tessella

@testset "Triangle distance in a local frame" begin
    # Plane x-y=offset: projection of (offset+2,1,1) has x=offset+1.5,
    # which is not representable near 1e16. Its distance is still sqrt(1/2).
    for offset in (0.0, 1e16, -1e16)
        triangle=((offset,0.0,0.0),(offset+4,4.0,0.0),(offset,0.0,4.0))
        field=DistanceField(triangles=[triangle])
        @test field_value(field,offset+2,1.0,1.0) ≈ sqrt(0.5) rtol=8eps(Float64)
        @test field_value(field,offset+2,2.0,1.0) == 0.0
        @test field_value(field,offset+2,3.0,0.5) ≈ sqrt(0.5) rtol=8eps(Float64)
        # Closest point on the edge/vertex still uses the bounded triangle.
        @test field_value(field,offset,0.0,-1.0) == 1.0
        @test field_value(field,offset,0.0,5.0) == 1.0
    end
    # Keep the triangle in a multi-leaf BVH: the public path must select its
    # distance rather than any of the unrelated point targets.
    field=DistanceField(points=[(1e16+100,Float64(i),0.0) for i in 1:12],
        triangles=[((1e16,0.0,0.0),(1e16+4,4.0,0.0),(1e16,0.0,4.0))])
    @test field_value(field,1e16+2,1.0,1.0) ≈ sqrt(0.5) rtol=8eps(Float64)
    # A near-plane offset can cancel in a Float64 dot product; use the exact
    # residual, rather than treating the rounded zero as an on-face query.
    field=DistanceField(triangles=[((0.,0.,0.),(1.,0.,1.),(0.,1.,1.))])
    @test field_value(field,0.25,nextfloat(0.25),0.5) ≈ eps(0.25)/sqrt(3.) rtol=8eps(Float64)
    # Huge, thin triangles can overflow the region determinants even when
    # their normal and the final bounded-triangle distance are representable.
    field=DistanceField(triangles=[((0.,0.,0.),(1e150,0.,0.),(1e150,1e-150,0.))])
    @test field_value(field,5e149,2e-150,0.) ≈ 1.5e-150 rtol=8eps(Float64)
end

module ClosestTriangleRegionChecks
using Test,Tessella
const F=Tessella.SizeField

@testset "Triangle region signs retain small interior projections" begin
    a=(-1e16,0.,0.);b=(0.,0.,0.);c=(0.,1.,0.);p=(-.25,.25,1.)
    # The projection is strictly interior: the tiny positive weight at a
    # is essential even though the region determinants nearly cancel.
    for triangle in ((a,b,c),(b,c,a),(c,a,b),(a,c,b),(c,b,a),(b,a,c))
        field=F.DistanceField(;triangles=[triangle])
        @test F.field_value(field,p...)≈1.
        distance,point=F._auto_closest_point_triangle(p,triangle...)
        @test distance≈1.
        @test all(isapprox.(point,(-.25,.25,0.);rtol=8eps(Float64),atol=0.))
    end
end

@testset "Automatic gap geometry retains a local approach vector" begin
    for offset in (0.,1e16,-1e16)
        triangle=((offset,0.,0.),(offset+4,4.,0.),(offset,0.,4.))
        p=(offset+2,1.,1.)
        distance,_,direction=F._auto_closest_triangle(p,triangle...)
        @test distance≈sqrt(.5) rtol=8eps(Float64)
        @test all(isapprox.(direction,(-.5,.5,0.);rtol=8eps(Float64),atol=0.))
        field=F.DistanceField(;triangles=[triangle])
        @test F._auto_lfs_search(field,1,p,Set{Int32}(),2.)≈sqrt(.5) rtol=8eps(Float64)
    end
    # A signed plane residual can cancel even without a large translation.
    triangle=((0.,0.,0.),(1.,0.,1.),(0.,1.,1.))
    p=(.25,nextfloat(.25),.5)
    distance,_,_=F._auto_closest_triangle(p,triangle...)
    @test distance≈eps(.25)/sqrt(3.) rtol=8eps(Float64)
    for scale in (1e-100,1e-170,1e100)
        triangle=((0.,0.,0.),(scale,0.,0.),(0.,scale,0.))
        distance,_,_=F._auto_closest_triangle((.25scale,.25scale,scale),triangle...)
        @test distance≈scale rtol=8eps(Float64)
        field=F.DistanceField(;triangles=[triangle])
        @test F._auto_lfs_search(field,1,(.25scale,.25scale,scale),Set{Int32}(),2scale)≈scale rtol=8eps(Float64)
        @test F._auto_lfs_search(field,1,(3scale,0.,scale),Set{Int32}(),4scale)==4scale
    end
    # The 60-degree boundary is inclusive, with exact zero separation
    # accepted only for a genuine triangle normal.
    triangle=((0.,0.,0.),(1.,-1.,0.),(0.,0.,1.))
    field=F.DistanceField(;triangles=[triangle])
    @test F._auto_lfs_search(field,1,(-1.,0.,-1.),Set{Int32}(),3.)≈sqrt(2.)
    @test F._auto_lfs_search(field,1,(-1.,0.,-2.),Set{Int32}(),3.)==3.
    @test F._auto_lfs_search(field,1,(0.,0.,0.),Set{Int32}(),3.)==0.
    degenerate=F.DistanceField(;triangles=[((0.,0.,0.),(1.,0.,0.),(2.,0.,0.))])
    @test F._auto_lfs_search(degenerate,1,(0.,0.,0.),Set{Int32}(),3.)==3.
end

@testset "Conditioned feature cones retain the original vertex geometry" begin
    triangle=((-1e16,-1e16,-1e16),(0.,0.,0.),(10001.,20000.,30000.))
    field=F.DistanceField(;triangles=[triangle])
    points=((3.0136355700966675,-1.0272697766363257,-.9863657934603419),
            (3.0136968094223087,-1.0273922491636756,-.9863045602586332),
            (3.0137036137918245,-1.02740585722227,-.9862977565695545),
            (3.0137104181613403,-1.0274194652808644,-.9862909528804757),
            (3.0137716574869815,-1.0275419378082145,-.986229719678767))
    R=Rational{BigInt};a,b,c=map(x->map(R,x),triangle)
    normal=(R(10000),R(-19999),R(9999))
    for point in points
        p=map(R,point)
        # Both edge derivatives at b are nonpositive, proving b is the
        # closest point. Its cone has this exact integer normal.
        @test sum(p.*a)<=0 && sum(p.*c)<=0
        accepted=4sum(p.*normal)^2>=sum(p.*p)*sum(normal.*normal)
        expected=accepted ? hypot(point...) : 10.
        @test F._auto_lfs_search(field,1,point,Set{Int32}(),10.)≈expected rtol=8eps(Float64)
    end
end

@testset "Automatic field construction consumes certified coarse gaps" begin
    quad(z)=Tessella.MeshTypes.Mesh([0. 1. 1. 0.;0. 0. 1. 1.;z z z z];
        tris=Int32[1 1;2 3;3 4])
    field=F.AutomaticMeshSizeField([quad(0.),quad(1.)];n_nodes_per_circle=0,
        n_points_per_gap=2,hmin=.25,hmax=1.,hbulk=1.,features=true,smoothing=false)
    @test all(≈(.5),field.h)
    @test F.field_value(field,0.,0.,0.)≈.5
end
end
