using Test
using Tessella
using LinearAlgebra: norm

@testset "mesh transformations preserve axis scale invariance" begin
    # An independent analytic oracle: a quarter-turn about (1,1,0) sends z
    # to (1,-1,0)/sqrt(2), and reflection normal to it exchanges -x and y.
    source=Mesh(Float64[0 1 0 0;0 0 1 0;0 0 0 1];
                tets=reshape(Int32[1,2,3,4],4,1),tet_tag=Int32[17])
    original=copy(source.coords)
    origin=(0.,0.,0.)
    reference_rotation=rotate_mesh(source,origin,(1.,1.,0.),pi/2)
    reference_reflection=mirror_mesh(source,origin,(1.,1.,0.))
    @test reference_rotation.coords[:,4]≈[inv(sqrt(2.)),-inv(sqrt(2.)),0.] atol=4eps()
    @test reference_reflection.coords[:,2]≈[0.,-1.,0.] atol=4eps()
    for scale in (nextfloat(0.),floatmin(Float64)/8,floatmax(Float64))
        axis=(scale,scale,0.)
        rotated=rotate_mesh(source,origin,axis,pi/2)
        reflected=mirror_mesh(source,origin,axis)
        @test rotated.coords≈reference_rotation.coords atol=4eps()
        @test reflected.coords≈reference_reflection.coords atol=4eps()
        @test norm(rotated.coords[:,4])≈1. atol=4eps()
        @test norm(reflected.coords[:,2])≈1. atol=4eps()
        @test validate(rotated).ok
        @test validate(reflected).ok
        @test rotated.tets==source.tets
        @test reflected.tets==reference_reflection.tets
        @test rotated.tet_tag==reflected.tet_tag==source.tet_tag
        @test mirror_mesh(reflected,origin,axis).coords≈source.coords atol=8eps()
    end
    @test source.coords==original
    @test_throws ArgumentError rotate_mesh(source,origin,(0.,0.,0.),pi/2)
    @test_throws ArgumentError mirror_mesh(source,origin,(Inf,1.,0.))
end
