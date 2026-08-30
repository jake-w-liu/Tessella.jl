using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, node, tet_signed_volume, validate

const _MESH_TRANSFORM_API=Tessella.API

function _mesh_transform_positive_tets(mesh)
    return all(axes(mesh.tets,2)) do cell
        ids=mesh.tets[:,cell]
        tet_signed_volume(node(mesh,ids[1]),node(mesh,ids[2]),
                          node(mesh,ids[3]),node(mesh,ids[4]))>0
    end
end

function _mesh_transform_detached(first,second)
    return all((:coords,:segs,:tris,:tets,:seg_tag,:tri_tag,:tet_tag)) do field
        !Base.mightalias(getfield(first,field),getfield(second,field))
    end
end

@testset "atomic cached-mesh affine transformation through API" begin
    positive12=(2.0,0.0,0.0,1.0,
                0.0,3.0,0.0,-2.0,
                0.0,0.0,4.0,0.5)
    positive16=(positive12...,0.0,0.0,0.0,1.0)
    positive_matrix=[2.0 0.0 0.0 1.0;
                     0.0 3.0 0.0 -2.0;
                     0.0 0.0 4.0 0.5;
                     0.0 0.0 0.0 1.0]

    _MESH_TRANSFORM_API.finalize()
    @test_throws ArgumentError _MESH_TRANSFORM_API.mesh.affine_transform(positive12)

    try
        _MESH_TRANSFORM_API.initialize()
        @test_throws ArgumentError _MESH_TRANSFORM_API.mesh.affine_transform(positive12)
        @test_throws ArgumentError _MESH_TRANSFORM_API.mesh.affine_transform(
            positive12,[(3,1)])

        @test _MESH_TRANSFORM_API.model.add_box(
            0,0,0,1,1,1;tag=1)==1
        base=_MESH_TRANSFORM_API.mesh.generate(3)
        base_crc=mesh_crc(base)
        @test base_crc.sha==
              "e9f6cd048ad689d1566e9c6664824543863983b8df79d9c0fa50f1f35d31cf83"
        @test _MESH_TRANSFORM_API.model.get_bounding_box(-1,-1)==
              (0.0,0.0,0.0,1.0,1.0,1.0)

        transformed=_MESH_TRANSFORM_API.mesh.affine_transform(positive12)
        expected=similar(base.coords)
        expected[1,:].=2 .* base.coords[1,:] .+ 1
        expected[2,:].=3 .* base.coords[2,:] .- 2
        expected[3,:].=4 .* base.coords[3,:] .+ 0.5
        @test transformed.coords==expected
        @test transformed.segs==base.segs
        @test transformed.tris==base.tris
        @test transformed.tets==base.tets
        @test transformed.seg_tag==base.seg_tag
        @test transformed.tri_tag==base.tri_tag
        @test transformed.tet_tag==base.tet_tag
        @test validate(transformed).ok
        @test _mesh_transform_positive_tets(transformed)
        transformed_crc=mesh_crc(transformed)
        @test transformed_crc.bbox==((1.0,-2.0,0.5),(3.0,1.0,4.5))
        @test transformed_crc.sha==base_crc.sha

        cached=_MESH_TRANSFORM_API.mesh.get()
        @test cached.coords==expected
        @test _mesh_transform_detached(transformed,cached)
        transformed.coords[1,1]+=17
        @test _MESH_TRANSFORM_API.mesh.get().coords==expected

        failures=(
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(positive12[1:11]),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform((positive12...,0.0)),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(fill(1.0,3,3)),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(
                (0.0,0.0,0.0,0.0,
                 0.0,1.0,0.0,0.0,
                 0.0,0.0,1.0,0.0)),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(
                (Inf,positive12[2:end]...)),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(
                (true,positive12[2:end]...)),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(
                (floatmax(Float64),0.0,0.0,floatmax(Float64),
                 0.0,1.0,0.0,0.0,
                 0.0,0.0,1.0,0.0)),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(
                (positive12...,0.0,0.0,0.0,2.0)),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform("identity"),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(positive12,[(3,1)]),
            ()->_MESH_TRANSFORM_API.mesh.affine_transform(positive12,1),
        )
        for fail in failures
            @test_throws ArgumentError fail()
            after_failure=_MESH_TRANSFORM_API.mesh.get()
            @test after_failure.coords==expected
            @test mesh_crc(after_failure)==transformed_crc
        end

        @test _MESH_TRANSFORM_API.mesh.clear()===nothing
        regenerated=_MESH_TRANSFORM_API.mesh.generate(3)
        transformed16=_MESH_TRANSFORM_API.mesh.affine_transform(
            positive16,Int32[])
        @test transformed16.coords==expected
        @test mesh_crc(transformed16)==transformed_crc

        @test _MESH_TRANSFORM_API.mesh.clear()===nothing
        @test mesh_crc(_MESH_TRANSFORM_API.mesh.generate(3))==base_crc
        transformed_matrix=_MESH_TRANSFORM_API.mesh.affine_transform(positive_matrix)
        @test transformed_matrix.coords==expected
        @test mesh_crc(transformed_matrix)==transformed_crc

        @test _MESH_TRANSFORM_API.mesh.clear()===nothing
        original=_MESH_TRANSFORM_API.mesh.generate(3)
        reflection=(-1.0,0.0,0.0,2.0,
                    0.0,1.0,0.0,3.0,
                    0.0,0.0,1.0,4.0)
        reflected=_MESH_TRANSFORM_API.mesh.affine_transform(reflection)
        @test reflected.coords[1,:]==2 .- original.coords[1,:]
        @test reflected.coords[2,:]==original.coords[2,:] .+ 3
        @test reflected.coords[3,:]==original.coords[3,:] .+ 4
        @test reflected.tets==original.tets[[2,1,3,4],:]
        @test reflected.tris==original.tris[[1,3,2],:]
        @test validate(reflected).ok
        @test _mesh_transform_positive_tets(reflected)
        reflected_crc=mesh_crc(reflected)
        @test reflected_crc.bbox==((1.0,3.0,4.0),(2.0,4.0,5.0))
        @test reflected_crc.sha==base_crc.sha
        @test _MESH_TRANSFORM_API.model.get_bounding_box(-1,-1)==
              (0.0,0.0,0.0,1.0,1.0,1.0)

        @test _MESH_TRANSFORM_API.mesh.clear(())===nothing
        @test mesh_crc(_MESH_TRANSFORM_API.mesh.generate(3))==base_crc
    finally
        _MESH_TRANSFORM_API.finalize()
    end
end

@testset "mesh-only affine transforms retain periodic model relations" begin
    _MESH_TRANSFORM_API.finalize()
    try
        _MESH_TRANSFORM_API.initialize()
        for (tag,(x,y)) in enumerate(((0.0,0.0),(1.0,0.0),
                                      (1.0,1.0),(0.0,1.0)))
            @test _MESH_TRANSFORM_API.model.add_point(
                x,y,0;tag=tag,meshSize=0.5)==tag
        end
        for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            @test _MESH_TRANSFORM_API.model.add_line(
                first,last;tag=tag)==tag
        end
        @test _MESH_TRANSFORM_API.model.add_curve_loop(
            [1,2,3,4];tag=1)==1
        @test _MESH_TRANSFORM_API.model.add_plane_surface([1];tag=1)==1
        translate_x=(1.0,0.0,0.0,1.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0)
        @test _MESH_TRANSFORM_API.mesh.set_periodic(
            1,[2],[4],translate_x)===nothing
        generated=_MESH_TRANSFORM_API.mesh.generate(2)
        generated_crc=mesh_crc(generated)

        identity=(1.0,0.0,0.0,0.0,
                  0.0,1.0,0.0,0.0,
                  0.0,0.0,1.0,0.0)
        identical=_MESH_TRANSFORM_API.mesh.affine_transform(identity)
        @test mesh_crc(identical)==generated_crc
        @test length(_MESH_TRANSFORM_API.mesh.get_periodic_nodes(
            1,2).slave_nodes)==5

        rotate=(0.0,-1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,1.0,0.0)
        rotated=_MESH_TRANSFORM_API.mesh.affine_transform(rotate)
        @test validate(rotated).ok
        periodic_error=try
            _MESH_TRANSFORM_API.mesh.get_periodic_nodes(1,2)
            nothing
        catch err
            err
        end
        @test periodic_error isa ErrorException
        @test occursin("is not represented by a two-node mesh-edge chain",
                       sprint(showerror,periodic_error))

        @test _MESH_TRANSFORM_API.mesh.clear()===nothing
        regenerated=_MESH_TRANSFORM_API.mesh.generate(2)
        @test mesh_crc(regenerated)==generated_crc
        @test length(_MESH_TRANSFORM_API.mesh.get_periodic_nodes(
            1,2).slave_nodes)==5
    finally
        _MESH_TRANSFORM_API.finalize()
    end
end
