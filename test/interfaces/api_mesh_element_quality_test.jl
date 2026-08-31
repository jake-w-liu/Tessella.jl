using Test
using Tessella
using Tessella.MeshTypes: mesh_crc

const _MESH_QUALITY_API=Tessella.API

function _mesh_quality_api_fixture()
    return Mesh(
        Float64[0 2 0 0;
                0 0 3 0;
                0 0 0 4];
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1),
        seg_tag=Int32[11],tri_tag=Int32[22],tet_tag=Int32[33])
end

function _install_mesh_quality_fixture!(mesh)
    lock(_MESH_QUALITY_API.STATE_LOCK) do
        _MESH_QUALITY_API._replace_mesh_cache_locked!(
            _MESH_QUALITY_API._copy_mesh(mesh))
    end
    return nothing
end

@testset "cached linear-simplex element qualities through API" begin
    _MESH_QUALITY_API.finalize()
    @test_throws ArgumentError _MESH_QUALITY_API.mesh.get_element_qualities(
        UInt64[1])
    try
        _MESH_QUALITY_API.initialize()
        @test_throws ArgumentError _MESH_QUALITY_API.mesh.get_element_qualities(
            UInt64[1])

        fixture=_mesh_quality_api_fixture()
        fixture_crc=mesh_crc(fixture)
        _install_mesh_quality_fixture!(fixture)

        @test _MESH_QUALITY_API.mesh.get_element_qualities(
            UInt64[3,2,3],"volume")==[4.0,3.0,4.0]
        @test isapprox(
            _MESH_QUALITY_API.mesh.get_element_qualities(
                UInt64[2],"minSICN")[1],0.7994080650317894;
            atol=16eps(Float64),rtol=16eps(Float64))
        @test isapprox(
            _MESH_QUALITY_API.mesh.get_element_qualities(
                UInt64[3])[1],0.6988644589181429;
            atol=16eps(Float64),rtol=16eps(Float64))
        @test _MESH_QUALITY_API.mesh.get_element_qualities(
            UInt64[1],"innerRadius")==[1.0]
        @test isempty(_MESH_QUALITY_API.mesh.get_element_qualities(
            UInt64[],"minIsotropy"))

        detached=_MESH_QUALITY_API.mesh.get_element_qualities(
            UInt64[3,2],"volume")
        detached[1]=99.0
        @test _MESH_QUALITY_API.mesh.get_element_qualities(
            UInt64[3,2],"volume")==[4.0,3.0]

        for call in (
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(1),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(Bool[true]),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities([1.0]),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(UInt64[0]),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(UInt64[4]),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(
                UInt64[2,1],"minIsotropy"),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(
                UInt64[2],"unknown"),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(
                UInt64[2],"volume",true,1),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(
                UInt64[2],"volume",0,true),
            ()->_MESH_QUALITY_API.mesh.get_element_qualities(
                UInt64[2],"volume",1,2),
        )
            @test_throws ArgumentError call()
            @test mesh_crc(_MESH_QUALITY_API.mesh.get())==fixture_crc
        end

        translated=_MESH_QUALITY_API.mesh.affine_transform(
            (1.0,0.0,0.0,10.0,
             0.0,1.0,0.0,20.0,
             0.0,0.0,1.0,30.0))
        @test _MESH_QUALITY_API.mesh.get_element_qualities(
            UInt64[3,2,1],"volume")==[4.0,3.0,2.0]
        @test mesh_crc(_MESH_QUALITY_API.mesh.get())==mesh_crc(translated)

        @test _MESH_QUALITY_API.mesh.clear()===nothing
        @test_throws ArgumentError _MESH_QUALITY_API.mesh.get_element_qualities(
            UInt64[1])

        @test isempty(Docs.undocumented_names(
            Tessella.API.mesh;private=false))
        @test isempty(Test.detect_ambiguities(
            Tessella.API.mesh;recursive=true))
    finally
        _MESH_QUALITY_API.finalize()
    end
end
