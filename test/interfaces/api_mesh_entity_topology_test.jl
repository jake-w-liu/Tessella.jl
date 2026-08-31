using Test
using SHA
using Tessella
using Tessella.MeshTypes: mesh_crc, unique_edges, unique_faces

const _MESH_ENTITY_API=Tessella.API

function _mesh_entity_api_fixture()
    return Mesh(
        Float64[0 1 0 0;0 0 1 0;0 0 0 1];
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1),
        seg_tag=Int32[11],tri_tag=Int32[22],tet_tag=Int32[33])
end

function _install_mesh_entity_fixture!(mesh)
    lock(_MESH_ENTITY_API.STATE_LOCK) do
        _MESH_ENTITY_API._replace_mesh_cache_locked!(
            _MESH_ENTITY_API._copy_mesh(mesh))
    end
    return nothing
end

function _mesh_entity_api_digest()
    arrays=(
        _MESH_ENTITY_API.mesh.get_all_edges()...,
        _MESH_ENTITY_API.mesh.get_edges(
            UInt64[1,2,2,1,1,4,4,1])...,
        _MESH_ENTITY_API.mesh.get_all_faces(3)...,
        _MESH_ENTITY_API.mesh.get_faces(
            3,UInt64[1,2,3,3,2,1,1,2,4])...)
    stream=IOBuffer()
    for values in arrays
        write(stream,htol(UInt64(length(values))))
        if eltype(values)==Int32
            foreach(value->write(
                stream,htol(reinterpret(UInt32,value))),values)
        else
            foreach(value->write(stream,htol(UInt64(value))),values)
        end
    end
    return bytes2hex(SHA.sha256(take!(stream)))
end

@testset "cached global edge and face topology through API" begin
    _MESH_ENTITY_API.finalize()
    for call in (
        ()->_MESH_ENTITY_API.mesh.create_edges(),
        ()->_MESH_ENTITY_API.mesh.create_faces(),
        ()->_MESH_ENTITY_API.mesh.get_edges(UInt64[]),
        ()->_MESH_ENTITY_API.mesh.get_faces(3,UInt64[]),
        ()->_MESH_ENTITY_API.mesh.get_all_edges(),
        ()->_MESH_ENTITY_API.mesh.get_all_faces(3),
    )
        @test_throws ArgumentError call()
    end

    try
        _MESH_ENTITY_API.initialize()
        for call in (
            ()->_MESH_ENTITY_API.mesh.create_edges(),
            ()->_MESH_ENTITY_API.mesh.create_faces(),
            ()->_MESH_ENTITY_API.mesh.get_edges(UInt64[]),
            ()->_MESH_ENTITY_API.mesh.get_faces(3,UInt64[]),
            ()->_MESH_ENTITY_API.mesh.get_all_edges(),
            ()->_MESH_ENTITY_API.mesh.get_all_faces(3),
        )
            @test_throws ArgumentError call()
        end

        fixture=_mesh_entity_api_fixture()
        fixture_crc=mesh_crc(fixture)
        _install_mesh_entity_fixture!(fixture)
        cached=_MESH_ENTITY_API.LAST_MESH[]

        @test _MESH_ENTITY_API.mesh.get_all_edges()==
              (UInt64[],UInt64[])
        @test _MESH_ENTITY_API.mesh.get_all_faces(3)==
              (UInt64[],UInt64[])
        @test _MESH_ENTITY_API.mesh.get_all_faces(4)==
              (UInt64[],UInt64[])
        @test _MESH_ENTITY_API.mesh.get_edges(UInt64[])==
              (UInt64[],Int32[])
        @test _MESH_ENTITY_API.mesh.get_faces(3,UInt64[])==
              (UInt64[],Int32[])
        @test_throws ArgumentError _MESH_ENTITY_API.mesh.get_edges(
            UInt64[1,2])
        @test_throws ArgumentError _MESH_ENTITY_API.mesh.get_faces(
            3,UInt64[1,2,3])

        for call in (
            ()->_MESH_ENTITY_API.mesh.create_edges(1),
            ()->_MESH_ENTITY_API.mesh.create_edges([(3,1)]),
            ()->_MESH_ENTITY_API.mesh.create_faces(missing),
            ()->_MESH_ENTITY_API.mesh.create_faces(((2,1),)),
        )
            @test_throws ArgumentError call()
        end
        @test _MESH_ENTITY_API.mesh.get_all_edges()==
              (UInt64[],UInt64[])
        @test _MESH_ENTITY_API.mesh.get_all_faces(3)==
              (UInt64[],UInt64[])

        @test _MESH_ENTITY_API.mesh.create_edges()===nothing
        @test _MESH_ENTITY_API.mesh.create_edges(())===nothing
        @test _MESH_ENTITY_API.LAST_MESH[]===cached
        @test _MESH_ENTITY_API.mesh.get_all_edges()==
              (UInt64.(1:6),UInt64[
                  1,2, 2,3, 3,1, 4,1, 4,3, 4,2])
        @test _MESH_ENTITY_API.mesh.get_edges(
            UInt64[1,2,2,1, 1,3,3,1, 1,4,4,1])==
            (UInt64[1,1,3,3,4,4],Int32[1,-1,1,-1,1,-1])

        @test _MESH_ENTITY_API.mesh.create_faces()===nothing
        @test _MESH_ENTITY_API.mesh.create_faces([])===nothing
        @test _MESH_ENTITY_API.LAST_MESH[]===cached
        @test _MESH_ENTITY_API.mesh.get_all_faces(3)==
              (UInt64.(1:4),UInt64[
                  1,2,3, 1,2,4, 1,4,3, 4,2,3])
        @test _MESH_ENTITY_API.mesh.get_all_faces(4)==
              (UInt64[],UInt64[])
        @test _MESH_ENTITY_API.mesh.get_faces(
            3,UInt64[1,2,3,3,2,1, 1,2,4])==
            (UInt64[1,1,2],Int32[0,0,0])
        @test mesh_crc(_MESH_ENTITY_API.mesh.get())==fixture_crc

        expected_edges=_MESH_ENTITY_API.mesh.get_all_edges()
        expected_faces=_MESH_ENTITY_API.mesh.get_all_faces(3)
        for call in (
            ()->_MESH_ENTITY_API.mesh.get_edges(1),
            ()->_MESH_ENTITY_API.mesh.get_edges(UInt64[1]),
            ()->_MESH_ENTITY_API.mesh.get_edges(Any[true,2]),
            ()->_MESH_ENTITY_API.mesh.get_edges(UInt64[0,1]),
            ()->_MESH_ENTITY_API.mesh.get_edges(UInt64[1,1]),
            ()->_MESH_ENTITY_API.mesh.get_faces(true,UInt64[]),
            ()->_MESH_ENTITY_API.mesh.get_faces(5,UInt64[]),
            ()->_MESH_ENTITY_API.mesh.get_faces(3,UInt64[1,2]),
            ()->_MESH_ENTITY_API.mesh.get_faces(3,UInt64[1,1,2]),
            ()->_MESH_ENTITY_API.mesh.get_faces(4,UInt64[1,2,3,4]),
            ()->_MESH_ENTITY_API.mesh.get_all_faces(false),
        )
            @test_throws ArgumentError call()
            @test mesh_crc(_MESH_ENTITY_API.mesh.get())==fixture_crc
            @test _MESH_ENTITY_API.mesh.get_all_edges()==expected_edges
            @test _MESH_ENTITY_API.mesh.get_all_faces(3)==expected_faces
        end

        detached_edges=_MESH_ENTITY_API.mesh.get_all_edges()
        detached_faces=_MESH_ENTITY_API.mesh.get_all_faces(3)
        detached_query=_MESH_ENTITY_API.mesh.get_edges(UInt64[1,2])
        detached_edges[1][1]=99;detached_edges[2][1]=99
        detached_faces[1][1]=99;detached_faces[2][1]=99
        detached_query[1][1]=99;detached_query[2][1]=99
        @test _MESH_ENTITY_API.mesh.get_all_edges()==expected_edges
        @test _MESH_ENTITY_API.mesh.get_all_faces(3)==expected_faces
        @test _MESH_ENTITY_API.mesh.get_edges(UInt64[1,2])==
              (UInt64[1],Int32[1])

        @test _mesh_entity_api_digest()==
              "25dc2bbdbc1100899eaf611a7b80d25cdd16ea19c0352c001dd29a3750affbda"

        translated=_MESH_ENTITY_API.mesh.affine_transform(
            (1.0,0.0,0.0,10.0,
             0.0,1.0,0.0,20.0,
             0.0,0.0,1.0,30.0))
        @test _MESH_ENTITY_API.mesh.get_all_edges()==
              (UInt64[],UInt64[])
        @test _MESH_ENTITY_API.mesh.get_all_faces(3)==
              (UInt64[],UInt64[])
        _MESH_ENTITY_API.mesh.create_edges()
        _MESH_ENTITY_API.mesh.create_faces()
        @test _MESH_ENTITY_API.mesh.get_all_edges()==expected_edges
        @test _MESH_ENTITY_API.mesh.get_all_faces(3)==expected_faces
        @test mesh_crc(_MESH_ENTITY_API.mesh.get())==mesh_crc(translated)

        _install_mesh_entity_fixture!(fixture)
        _MESH_ENTITY_API.mesh.create_edges()
        _MESH_ENTITY_API.mesh.create_faces()
        refined=_MESH_ENTITY_API.mesh.refine()
        @test _MESH_ENTITY_API.mesh.get_all_edges()==
              (UInt64[],UInt64[])
        @test _MESH_ENTITY_API.mesh.get_all_faces(3)==
              (UInt64[],UInt64[])
        _MESH_ENTITY_API.mesh.create_edges()
        _MESH_ENTITY_API.mesh.create_faces()
        refined_edges=_MESH_ENTITY_API.mesh.get_all_edges()
        refined_faces=_MESH_ENTITY_API.mesh.get_all_faces(3)
        @test length(refined_edges[1])==
              length(unique_edges(refined.tris,refined.tets))
        @test length(refined_faces[1])==
              length(unique_faces(refined.tris,refined.tets))

        @test _MESH_ENTITY_API.mesh.clear()===nothing
        for call in (
            ()->_MESH_ENTITY_API.mesh.get_all_edges(),
            ()->_MESH_ENTITY_API.mesh.get_all_faces(3),
            ()->_MESH_ENTITY_API.mesh.create_edges(),
            ()->_MESH_ENTITY_API.mesh.create_faces(),
        )
            @test_throws ArgumentError call()
        end

        @test isempty(Docs.undocumented_names(
            Tessella.API.mesh;private=false))
        @test isempty(Test.detect_ambiguities(
            Tessella.API.mesh;recursive=true))
    finally
        _MESH_ENTITY_API.finalize()
    end
end
