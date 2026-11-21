using Test
using Tessella
using Tessella.MeshTypes: mesh_crc

const _MESH_LOCATION_API=Tessella.API

function _mesh_location_fixture()
    return Mesh(
        Float64[0 1 0 0;0 0 1 0;0 0 0 1];
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1),
        seg_tag=Int32[11],tri_tag=Int32[22],tet_tag=Int32[33])
end

function _mesh_location_install!(mesh;reset_locator=true)
    lock(_MESH_LOCATION_API.STATE_LOCK) do
        copied=_MESH_LOCATION_API._copy_mesh(mesh)
        if reset_locator
            _MESH_LOCATION_API._replace_mesh_cache_locked!(copied)
        else
            _MESH_LOCATION_API.LAST_MESH[]=copied
        end
        _MESH_LOCATION_API._native_lookup_invalidate_locked!()
    end
    return nothing
end

function _mesh_location_segment_fixture(count::Int)
    coordinates=zeros(3,count+1)
    segments=Matrix{Int32}(undef,2,count)
    for segment in 1:count
        coordinates[1,segment+1]=segment
        segments[:,segment].=Int32[segment,segment+1]
    end
    return Mesh(coordinates;segs=segments)
end

@noinline function _mesh_location_query_allocation(mesh)
    _mesh_location_install!(mesh)
    _MESH_LOCATION_API.mesh.get_element_by_coordinates(0.5,0.0,0.0,1,true)
    GC.gc()
    return @allocated _MESH_LOCATION_API.mesh.get_element_by_coordinates(
        0.5,0.0,0.0,1,true)
end

@testset "cached simplex point location through API" begin
    _MESH_LOCATION_API.finalize()
    @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_element_by_coordinates(
        0,0,0)
    @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
        0,0,0)
    @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
        1,0,0,0)

    try
        _MESH_LOCATION_API.initialize()
        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_element_by_coordinates(
            0,0,0)
        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0,0,0)
        # A tag that names no element is Gmsh's own diagnostic, raised under the
        # API default General.AbortOnError=2 (Gmsh 4.15.2: "Unknown element 1").
        @test_throws ErrorException("Unknown element 1") _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
            1,0,0,0)

        fixture=_mesh_location_fixture()
        fixture_crc=mesh_crc(fixture)
        _mesh_location_install!(fixture)
        # (This test used to follow the private `LAST_MESH_LOCATOR` cache through
        # these calls. Point queries have gone through the element octree since
        # c7cc213a and never fill that cache; the public results below are
        # what Gmsh 4.15.2 returns for the same mesh.)

        tags=_MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0.25,0.0,0.0,-1,true)
        @test tags==UInt64[3,2,1]
        first_result=_MESH_LOCATION_API.mesh.get_element_by_coordinates(
            0.25,0.0,0.0,-1,true)
        @test first_result==
              (UInt64(3),Int32(4),UInt64[1,2,3,4],0.25,0.0,0.0)
        @test _MESH_LOCATION_API.mesh.get_element_by_coordinates(
            0.25,0.0,0.0,2,true)==
              (UInt64(2),Int32(2),UInt64[1,2,3],0.25,0.0,0.0)
        @test _MESH_LOCATION_API.mesh.get_element_by_coordinates(
            0.25,0.0,0.0,1,true)==
              (UInt64(1),Int32(1),UInt64[1,2],-0.5,0.0,0.0)
        @test _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0.2,0.3,0.0,-1,true)==UInt64[3,2]
        @test _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0.1,0.2,0.3,-1,true)==UInt64[3]

        _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0.1,0.2,0.3,-1,true)

        tags[1]=99
        first_result[3][1]=99
        @test _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0.25,0.0,0.0,-1,true)==UInt64[3,2,1]
        @test _MESH_LOCATION_API.mesh.get_element_by_coordinates(
            0.25,0.0,0.0,-1,true)[3]==UInt64[1,2,3,4]

        # Off the line, the Newton inversion of a 1-D element extrapolates along
        # its completed frame: Gmsh 4.15.2 returns (-0.5, -38.0, -57.0) here.
        @test _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
            1,0.25,2.0,3.0)==(-0.5,-38.0,-57.0)
        @test _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
            1,0.25,0.0,0.0)==(-0.5,0.0,0.0)
        @test _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
            2,0.2,0.3,5.0)==(0.2,0.3,0.0)
        @test _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
            3,-0.25,0.5,1.25)==(-0.25,0.5,1.25)

        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            -5.0e-6,0.2,0.2,3,true)
        @test _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            -5.0e-6,0.2,0.2,3,false)==UInt64[3]
        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            -1.1,0.1,0.1,3,false)
        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0.25,0.0,0.0,0,true)

        for call in (
            ()->_MESH_LOCATION_API.mesh.get_element_by_coordinates(NaN,0,0),
            ()->_MESH_LOCATION_API.mesh.get_element_by_coordinates(0,Inf,0),
            ()->_MESH_LOCATION_API.mesh.get_element_by_coordinates(false,0,0),
            ()->_MESH_LOCATION_API.mesh.get_element_by_coordinates(0,0,0,false),
            ()->_MESH_LOCATION_API.mesh.get_element_by_coordinates(0,0,0,4),
            ()->_MESH_LOCATION_API.mesh.get_element_by_coordinates(
                0,0,0,-1,missing),
            ()->_MESH_LOCATION_API.mesh.get_elements_by_coordinates(
                0,0,0,big(typemax(Int))+1),
            ()->_MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
                true,0,0,0),
            ()->_MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
                1,missing,0,0),
        )
            @test_throws ArgumentError call()
            @test mesh_crc(_MESH_LOCATION_API.mesh.get())==fixture_crc
        end
        for unknown in (0,4)
            @test_throws ErrorException("Unknown element $unknown") _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
                unknown,0,0,0)
            @test mesh_crc(_MESH_LOCATION_API.mesh.get())==fixture_crc
        end

        # A direct internal cache replacement cannot make an already-built locator
        # stale: identity mismatch forces a rebuild on the next query.
        directly_moved=Mesh(
            fixture.coords .+ [10.0,20.0,30.0];
            segs=fixture.segs,tris=fixture.tris,tets=fixture.tets,
            seg_tag=fixture.seg_tag,tri_tag=fixture.tri_tag,
            tet_tag=fixture.tet_tag)
        _mesh_location_install!(directly_moved;reset_locator=false)
        @test _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            10.25,20.0,30.0,-1,true)==UInt64[3,2,1]

        @test _MESH_LOCATION_API.model.add_point(9,9,9;tag=1)==1
        @test _MESH_LOCATION_API.LAST_MESH[]===nothing

        _mesh_location_install!(fixture)
        _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0.25,0.0,0.0,-1,true)
        translation=(1.0,0.0,0.0,2.0,
                     0.0,1.0,0.0,3.0,
                     0.0,0.0,1.0,4.0)
        moved=_MESH_LOCATION_API.mesh.affine_transform(translation)
        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0.25,0.0,0.0,-1,true)
        @test _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            2.25,3.0,4.0,-1,true)==UInt64[3,2,1]
        @test moved.coords==fixture.coords .+ [2.0,3.0,4.0]

        refined=_MESH_LOCATION_API.mesh.refine()
        refined_tags=_MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            2.25,3.0,4.0,-1,true)
        @test !isempty(refined_tags)
        @test maximum(refined_tags)<=
              _MESH_LOCATION_API.mesh.get_max_element_tag()
        @test mesh_crc(_MESH_LOCATION_API.mesh.get())==mesh_crc(refined)

        @test _MESH_LOCATION_API.mesh.clear()===nothing
        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_element_by_coordinates(
            2.25,3.0,4.0)

        _mesh_location_install!(Mesh(zeros(3,0)))
        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_element_by_coordinates(
            0,0,0)
        @test_throws ArgumentError _MESH_LOCATION_API.mesh.get_elements_by_coordinates(
            0,0,0)
        # A tag that names no element is Gmsh's own diagnostic, raised under the
        # API default General.AbortOnError=2 (Gmsh 4.15.2: "Unknown element 1").
        @test_throws ErrorException("Unknown element 1") _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
            1,0,0,0)

        degenerate=Mesh(
            zeros(3,2);segs=reshape(Int32[1,2],2,1))
        _mesh_location_install!(degenerate)
        # A zero-length Line is located like any other element: Gmsh 4.15.2
        # returns element 1 with local coordinates (0,0,0) from both queries.
        @test _MESH_LOCATION_API.mesh.get_element_by_coordinates(0,0,0)==
              (UInt64(1),Int32(1),UInt64[1,2],0.0,0.0,0.0)
        @test _MESH_LOCATION_API.mesh.get_local_coordinates_in_element(
            1,0,0,0)==(0.0,0.0,0.0)

        allocation_small=_mesh_location_segment_fixture(5_000)
        allocation_large=_mesh_location_segment_fixture(10_000)
        small_bytes=_mesh_location_query_allocation(allocation_small)
        large_bytes=_mesh_location_query_allocation(allocation_large)
        @test small_bytes>0
        @test large_bytes<=small_bytes+4_096
    finally
        _MESH_LOCATION_API.finalize()
    end
end
