# P3 native mesh transformations. The expected coordinates below are analytic
# affine-map oracles; orientation and topology are checked independently.

using Test
using Tessella
using Tessella.MeshTypes: node, nnodes, nsegs, ntris, ntets, tet_signed_volume

function _transform_fixture()
    coordinates=Float64[0 1 0 0;
                        0 0 1 0;
                        0 0 0 1]
    segments=reshape(Int32[1,2],2,1)
    triangles=Int32[1 1 2 3;
                    3 2 3 1;
                    2 4 4 4]
    tetrahedra=reshape(Int32[1,2,3,4],4,1)
    return Mesh(coordinates;segs=segments,tris=triangles,tets=tetrahedra,
                seg_tag=Int32[7],tri_tag=Int32[11,12,13,14],tet_tag=Int32[19])
end

@inline function _transform_signed_tet(mesh)
    ids=mesh.tets[:,1]
    return tet_signed_volume(node(mesh,ids[1]),node(mesh,ids[2]),
                             node(mesh,ids[3]),node(mesh,ids[4]))
end

@noinline function _transform_allocated(n::Int)
    coordinates=Matrix{Float64}(undef,3,n)
    @inbounds for i in 1:n
        coordinates[1,i]=i/n
        coordinates[2,i]=(i%17)/17
        coordinates[3,i]=(i%29)/29
    end
    mesh=Mesh(coordinates)
    translate_mesh(mesh,(1.,2.,3.))
    GC.gc()
    return @allocated translate_mesh(mesh,(1.,2.,3.))
end

@testset "validated affine mesh transformations" begin
    source=_transform_fixture()
    @test validate(source).ok
    @test _transform_signed_tet(source)>0

    @testset "identity, tags, topology, and deterministic CRC" begin
        identity_map=affine_transform(source,[1 0 0;0 1 0;0 0 1])
        @test identity_map.coords==source.coords
        @test identity_map.segs==source.segs
        @test identity_map.tris==source.tris
        @test identity_map.tets==source.tets
        @test identity_map.seg_tag==source.seg_tag
        @test identity_map.tri_tag==source.tri_tag
        @test identity_map.tet_tag==source.tet_tag
        @test mesh_crc(identity_map)==mesh_crc(source)

        first=translate_mesh(source,(2.,-3.,4.))
        second=translate_mesh(source,(2.,-3.,4.))
        @test mesh_crc(first)==mesh_crc(second)
        @test first.coords==source.coords .+ [2.,-3.,4.]
        @test translate_mesh(first,(-2.,3.,-4.)).coords==source.coords
        @test first.segs==source.segs
        @test first.tris==source.tris
        @test first.tets==source.tets
        @test first.seg_tag==source.seg_tag
        @test first.tri_tag==source.tri_tag
        @test first.tet_tag==source.tet_tag
    end

    @testset "rotation and affine analytic oracles" begin
        rotated=rotate_mesh(source,(0.,0.,0.),(0.,0.,1.),pi/2)
        expected=Float64[0 0 -1 0;
                         0 1  0 0;
                         0 0  0 1]
        @test rotated.coords≈expected atol=4eps(Float64)
        @test _transform_signed_tet(rotated)>0

        matrix=[2.0 0.5 0.0;-0.25 1.5 0.0;0.0 0.0 3.0]
        origin=(1.,-2.,.5);shift=(-4.,1.,2.)
        mapped=affine_transform(source,matrix;origin=origin,translation=shift)
        oracle=Matrix{Float64}(undef,3,nnodes(source))
        @inbounds for i in axes(source.coords,2)
            p=source.coords[:,i]
            oracle[:,i]=collect(origin)+matrix*(p-collect(origin))+collect(shift)
        end
        @test mapped.coords≈oracle rtol=8eps(Float64) atol=8eps(Float64)
        @test validate(mapped).ok
        @test _transform_signed_tet(mapped)>0
    end

    @testset "dilation and reflection preserve orientation" begin
        expanded=dilate_mesh(source,(0.,0.,0.),(2.,3.,4.))
        @test expanded.coords==source.coords .* [2.,3.,4.]
        @test _transform_signed_tet(expanded)≈24*_transform_signed_tet(source)

        reflected=dilate_mesh(source,(0.,0.,0.),(-1.,1.,1.))
        mirrored=mirror_mesh(source,(0.,0.,0.),(1.,0.,0.))
        @test mirrored.coords≈reflected.coords atol=4eps(Float64)
        @test reflected.tris[1,:]==source.tris[1,:]
        @test reflected.tris[2,:]==source.tris[3,:]
        @test reflected.tris[3,:]==source.tris[2,:]
        @test reflected.tets[:,1]==source.tets[Int32[2,1,3,4],1]
        @test _transform_signed_tet(reflected)>0
        @test _transform_signed_tet(mirrored)>0
        @test validate(reflected).ok
        @test validate(mirrored).ok

        shifted_plane=mirror_mesh(source,(1.,0.,0.),(1.,0.,0.))
        @test shifted_plane.coords[1,:]≈2 .- source.coords[1,:]
        @test shifted_plane.coords[2:3,:]==source.coords[2:3,:]
    end

    @testset "exact certificates and extreme-coordinate fallback" begin
        # The 2×2 determinant is one Float64 ULP, so a tolerance-based singularity
        # test would reject a genuinely invertible map.
        near_singular=[1.0 1.0 0.0;1.0 nextfloat(1.0) 0.0;0.0 0.0 1.0]
        isolated=Mesh(reshape(Float64[1,2,3],3,1))
        @test validate(affine_transform(isolated,near_singular)).ok

        # `point-origin` overflows in Float64, while the exact affine result is zero.
        extreme=Mesh(reshape(Float64[-floatmax(Float64),0,0],3,1))
        cancelled=affine_transform(extreme,[.5 0 0;0 1 0;0 0 1];
                                   origin=(floatmax(Float64),0.,0.))
        @test node(cancelled,1)==(0.,0.,0.)

        @test_throws ArgumentError affine_transform(
            Mesh(reshape(Float64[floatmax(Float64),0,0],3,1)),
            [2.0 0 0;0 1 0;0 0 1])
        @test_throws ArgumentError affine_transform(isolated,[1 1 0;1 1 0;0 0 1])
    end

    @testset "input and error contracts" begin
        @test_throws ArgumentError affine_transform(source,[1 0;0 1])
        @test_throws ArgumentError affine_transform(source,[1 0 0;0 Inf 0;0 0 1])
        @test_throws ArgumentError affine_transform(source,Any[true 0 0;0 1 0;0 0 1])
        @test_throws ArgumentError affine_transform(source,ComplexF64[1 0 0;0 1 0;0 0 1])
        @test_throws ArgumentError affine_transform(source,[1 0 0;0 1 0;0 0 1];
                                                    origin=(0.,0.))
        @test_throws ArgumentError translate_mesh(source,(0.,NaN,0.))
        @test_throws ArgumentError rotate_mesh(source,(0.,0.,0.),(0.,0.,0.),1.)
        @test_throws ArgumentError rotate_mesh(source,(0.,0.,0.),(0.,0.,1.),true)
        @test_throws ArgumentError dilate_mesh(source,(0.,0.,0.),(1.,0.,1.))
        @test_throws ArgumentError mirror_mesh(source,(0.,0.,0.),(0.,0.,0.))

        inverted=Mesh(source.coords;tets=reshape(Int32[1,2,4,3],4,1))
        @test !validate(inverted).ok
        @test_throws ArgumentError translate_mesh(inverted,(1.,0.,0.))
    end

    @testset "resource growth" begin
        small=_transform_allocated(10_000)
        large=_transform_allocated(20_000)
        @test small>0
        @test large>small
        @test large<=2.25small+64_000
    end
end
