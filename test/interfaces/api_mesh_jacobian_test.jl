using Test
using Tessella
using Tessella.MeshTypes: mesh_crc

const _MESH_JACOBIAN_API=Tessella.API

function _mesh_jacobian_api_fixture()
    return Mesh(
        Float64[0 2 0 0  1 0 0 0;
                0 0 3 0  0 0 1 0;
                0 0 0 4  0 1 0 4];
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=Int32[1 5;2 6;3 7;4 8],
        seg_tag=Int32[11],tri_tag=Int32[22],tet_tag=Int32[33,44])
end

function _install_mesh_jacobian_fixture!(mesh)
    lock(_MESH_JACOBIAN_API.STATE_LOCK) do
        _MESH_JACOBIAN_API._replace_mesh_cache_locked!(
            _MESH_JACOBIAN_API._copy_mesh(mesh))
    end
    return nothing
end

@testset "cached linear-simplex Jacobians through API" begin
    _MESH_JACOBIAN_API.finalize()
    @test_throws ArgumentError _MESH_JACOBIAN_API.mesh.get_jacobians(
        4,[0,0,0])
    @test_throws ArgumentError _MESH_JACOBIAN_API.mesh.get_jacobian(
        1,[0,0,0])
    try
        _MESH_JACOBIAN_API.initialize()
        @test_throws ArgumentError _MESH_JACOBIAN_API.mesh.get_jacobians(
            4,[0,0,0])
        @test_throws ArgumentError _MESH_JACOBIAN_API.mesh.get_jacobian(
            1,[0,0,0])

        fixture=_mesh_jacobian_api_fixture()
        fixture_crc=mesh_crc(fixture)
        _install_mesh_jacobian_fixture!(fixture)

        jacobians,determinants,coordinates=
            _MESH_JACOBIAN_API.mesh.get_jacobians(
                4,[0,0,0, 0.2,0.3,0.1])
        @test jacobians==[
            2.0,0,0,0,3.0,0,0,0,4.0,
            2.0,0,0,0,3.0,0,0,0,4.0,
            -1.0,0,1.0,-1.0,1.0,0,-1.0,0,4.0,
            -1.0,0,1.0,-1.0,1.0,0,-1.0,0,4.0]
        @test determinants==[24.0,24.0,-3.0,-3.0]
        @test isapprox(
            coordinates,
            [0.0,0.0,0.0, 0.4,0.9,0.4,
             1.0,0.0,0.0, 0.4,0.3,0.6];
            atol=4eps(Float64),rtol=4eps(Float64))
        positive=_MESH_JACOBIAN_API.mesh.get_jacobian(
            3,[0.2,0.3,0.1])
        @test positive[1]==[2.0,0,0,0,3.0,0,0,0,4.0]
        @test positive[2]==[24.0]
        @test isapprox(
            positive[3],[0.4,0.9,0.4];
            atol=4eps(Float64),rtol=4eps(Float64))
        inverted=_MESH_JACOBIAN_API.mesh.get_jacobian(
            4,[0.2,0.3,0.1])
        @test inverted[1]==[-1.0,0,1.0,-1.0,1.0,0,-1.0,0,4.0]
        @test inverted[2]==[-3.0]
        @test isapprox(
            inverted[3],[0.4,0.3,0.6];
            atol=4eps(Float64),rtol=4eps(Float64))
        @test _MESH_JACOBIAN_API.mesh.get_jacobians(
            4,[])==(Float64[],Float64[],Float64[])

        # Task partitioning selects contiguous Gmsh blocks: with two cached
        # tetrahedra, task 0 of 2 is the positive tet and task 1 of 2 the
        # inverted one; task>=num_tasks is the silently-empty Gmsh range.
        single_point=[0,0,0]
        @test _MESH_JACOBIAN_API.mesh.get_jacobians(
            4,single_point,-1,0,2)==(
            [2.0,0,0,0,3.0,0,0,0,4.0],[24.0],[0.0,0.0,0.0])
        @test _MESH_JACOBIAN_API.mesh.get_jacobians(
            4,single_point,-1,1,2)==(
            [-1.0,0,1.0,-1.0,1.0,0,-1.0,0,4.0],[-3.0],[1.0,0.0,0.0])
        @test _MESH_JACOBIAN_API.mesh.get_jacobians(
            4,single_point,-1,2,2)==(Float64[],Float64[],Float64[])
        @test _MESH_JACOBIAN_API.mesh.get_jacobians(
            3,single_point,-1,1,2)==(Float64[],Float64[],Float64[])
        two_point=(
            _MESH_JACOBIAN_API.mesh.get_jacobians(4,[0,0,0,0.2,0.3,0.1]))
        second_half=_MESH_JACOBIAN_API.mesh.get_jacobians(
            4,[0,0,0,0.2,0.3,0.1],-1,1,2)
        @test second_half[1]==two_point[1][19:36]
        @test second_half[2]==two_point[2][3:4]
        @test second_half[3]==two_point[3][7:12]
        first_half=_MESH_JACOBIAN_API.mesh.get_jacobians(
            4,[0,0,0,0.2,0.3,0.1],-1,0,2)
        @test vcat(first_half[1],second_half[1])==two_point[1]
        @test vcat(first_half[2],second_half[2])==two_point[2]
        @test vcat(first_half[3],second_half[3])==two_point[3]
        @test mesh_crc(_MESH_JACOBIAN_API.mesh.get())==fixture_crc

        jacobians[1]=99.0;determinants[1]=99.0;coordinates[1]=99.0
        @test _MESH_JACOBIAN_API.mesh.get_jacobian(
            3,[0,0,0])==
            ([2.0,0,0,0,3.0,0,0,0,4.0],[24.0],[0.0,0.0,0.0])
        @test mesh_crc(_MESH_JACOBIAN_API.mesh.get())==fixture_crc

        for call in (
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(true,[0,0,0]),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(999,[0,0,0]),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(4,[0,0]),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(4,[NaN,0,0]),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(
                4,Any[true,0,0]),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(4,[0,0,0],0),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(4,[0,0],-1,1,2),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(4,[0,0,0],-1,-1,1),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(4,[0,0,0],-1,0,0),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobians(4,[0,0,0],-1,true,1),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobian(true,[0,0,0]),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobian(0,[0,0,0]),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobian(5,[0,0,0]),
            ()->_MESH_JACOBIAN_API.mesh.get_jacobian(1,[Inf,0,0]),
        )
            @test_throws ArgumentError call()
            @test mesh_crc(_MESH_JACOBIAN_API.mesh.get())==fixture_crc
        end

        transform_fixture=Mesh(
            fixture.coords;
            segs=fixture.segs,tris=fixture.tris,
            tets=fixture.tets[:,1:1],seg_tag=fixture.seg_tag,
            tri_tag=fixture.tri_tag,tet_tag=fixture.tet_tag[1:1])
        _install_mesh_jacobian_fixture!(transform_fixture)
        translated=_MESH_JACOBIAN_API.mesh.affine_transform(
            (1.0,0.0,0.0,10.0,
             0.0,1.0,0.0,20.0,
             0.0,0.0,1.0,30.0))
        translated_result=_MESH_JACOBIAN_API.mesh.get_jacobian(
            3,[0.2,0.3,0.1])
        @test translated_result[1]==[2.0,0,0,0,3.0,0,0,0,4.0]
        @test translated_result[2]==[24.0]
        @test translated_result[3]==[10.4,20.9,30.4]
        @test mesh_crc(_MESH_JACOBIAN_API.mesh.get())==mesh_crc(translated)

        @test _MESH_JACOBIAN_API.mesh.clear()===nothing
        @test_throws ArgumentError _MESH_JACOBIAN_API.mesh.get_jacobians(
            4,[0,0,0])
        @test_throws ArgumentError _MESH_JACOBIAN_API.mesh.get_jacobian(
            1,[0,0,0])

        @test isempty(Docs.undocumented_names(
            Tessella.API.mesh;private=false))
        @test isempty(Test.detect_ambiguities(
            Tessella.API.mesh;recursive=true))
    finally
        _MESH_JACOBIAN_API.finalize()
    end
end
