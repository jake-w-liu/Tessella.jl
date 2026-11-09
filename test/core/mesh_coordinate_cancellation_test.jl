module MeshCoordinateCancellationChecks
using Test, Tessella, LinearAlgebra
const L=Tessella.MeshPointLocation
const R=Rational{BigInt}

function coordinate_checksum(mesh,point,count)
    value=0.0
    for _ in 1:count
        u,v,w=L.mesh_local_coordinates(mesh,1,point...)
        value+=u+v+w
    end
    return value
end

@testset "linear-simplex inversion certifies numerator cancellation" begin
    # These integer input coordinates are exactly orthogonal, even at the
    # largest translation along the segment's normal. A well-conditioned
    # denominator cannot certify a cancelled projection numerator.
    for edge in ((3.,5.,0.),(7.,11.,0.),(13.,19.,0.)), scale in (1e10,1e14,1e15)
        mesh=Mesh(hcat(zeros(3),collect(edge));segs=reshape(Int32[1,2],2,1))
        point=(edge[2]*scale,-edge[1]*scale,0.)
        @test sum(R(point[i])*R(edge[i]) for i in 1:3)==0
        @test L.mesh_local_coordinates(mesh,1,point...)==(-1.,0.,0.)
    end

    triangle=Mesh([0. 3. 7.;0. 5. 11.;0. 1. 2.];
        tris=reshape(Int32[1,2,3],3,1))
    for scale in (1e10,1e14,1e15)
        @test L.mesh_local_coordinates(triangle,1,-scale,scale,-2scale)==(0.,0.,0.)
    end

    coordinates=[0. 7. 3. 2.;0. 11. 5. 1.;0. 2. 1. 3.]
    tetrahedron=Mesh(coordinates;tets=reshape(Int32[1,2,3,4],4,1))
    @test validate(tetrahedron).ok
    for scale in (1e10,1e14), missing in 1:3
        reference=ntuple(i->i==missing ? 0. : scale,3)
        point=coordinates[:,2:4]*collect(reference)
        expected=Float64.(R.(coordinates[:,2:4])\R.(point))
        actual=L.mesh_local_coordinates(tetrahedron,1,point...)
        @test collect(actual)≈expected rtol=64eps(Float64) atol=0.
        @test actual[missing]==0.
    end
end

@testset "tiny solve numerators retain representable coordinates" begin
    triangle=Mesh([0. 1. 0.;0. 0. 1e-100;0. 0. 0.];
        tris=reshape(Int32[1,2,3],3,1))
    tetrahedron=Mesh([0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. 1e-200];
        tets=reshape(Int32[1,2,3,4],4,1))
    for (mesh,point) in ((triangle,(1e-200,2.5e-101,0.)),
                         (tetrahedron,(1e-200,.25,2.5e-201)))
        actual=L.mesh_local_coordinates(mesh,1,point...)
        @test actual[1]≈1e-200 rtol=64eps(Float64) atol=0.
        @test actual[2]≈.25 rtol=64eps(Float64)
        @test L.locate_elements(L.SimplexLocator(mesh),point...;strict=true)==UInt64[1]
    end
end

@testset "subnormal solve denominators use the exact ratio" begin
    triangle=Mesh([0. 1. 0.;0. 0. 1e-160;0. 0. 0.];
        tris=reshape(Int32[1,2,3],3,1))
    tetrahedron=Mesh([0. 1. 0. 0.;0. 0. 1e-160 0.;0. 0. 0. 1e-160];
        tets=reshape(Int32[1,2,3,4],4,1))
    expected=Float64(R(1e-60)/R(1e-160))
    @test L.mesh_local_coordinates(triangle,1,0.,1e-60,0.)[2]≈expected rtol=64eps(Float64)
    @test L.mesh_local_coordinates(tetrahedron,1,0.,0.,1e-60)[3]≈expected rtol=64eps(Float64)
end

@testset "ordinary affine coordinate queries retain their allocation bound" begin
    fixtures=(
        (Mesh([0. 3.;0. 5.;0. 0.];segs=reshape(Int32[1,2],2,1)),(.9,1.5,0.)),
        (Mesh([0. 2. .125;0. .25 1.;0. 0. 0.];
            tris=reshape(Int32[1,2,3],3,1)),(.4,.25,0.)),
        (Mesh([0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. 1.];
            tets=reshape(Int32[1,2,3,4],4,1)),(.2,.3,.1)))
    for (mesh,point) in fixtures
        value=coordinate_checksum(mesh,point,2)
        @test isfinite(value)
        if Base.JLOptions().compile_enabled!=3
            @test (@allocated coordinate_checksum(mesh,point,100))==0
        end
    end
end
end
