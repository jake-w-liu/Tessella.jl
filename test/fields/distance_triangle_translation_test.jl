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
