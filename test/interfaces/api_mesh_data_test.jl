using Test
using SHA
using Tessella
using Tessella.MeshTypes: mesh_crc

const _MESH_DATA_API=Tessella.API

function _mesh_data_fixture()
    coordinates=Float64[0 1 0 0;
                        0 0 1 0;
                        0 0 0 1]
    return Mesh(
        coordinates;
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1),
        seg_tag=Int32[11],tri_tag=Int32[22],tet_tag=Int32[33])
end

function _mesh_data_install!(mesh)
    lock(_MESH_DATA_API.STATE_LOCK) do
        _MESH_DATA_API.LAST_MESH[]=_MESH_DATA_API._copy_mesh(mesh)
    end
    return nothing
end

function _mesh_data_segment_fixture(count::Int)
    coordinates=zeros(3,count+1)
    segments=Matrix{Int32}(undef,2,count)
    for segment in 1:count
        coordinates[1,segment+1]=segment
        segments[:,segment].=Int32[segment,segment+1]
    end
    return Mesh(coordinates;segs=segments)
end

function _mesh_data_query_sha(groups)
    stream=IOBuffer()
    for values in groups
        marker=eltype(values)===UInt64 ? UInt8(1) : UInt8(2)
        write(stream,marker)
        write(stream,htol(UInt64(length(values))))
        for value in values
            bits=value isa Float64 ? reinterpret(UInt64,value) : UInt64(value)
            write(stream,htol(bits))
        end
    end
    return bytes2hex(SHA.sha256(take!(stream)))
end

@noinline function _mesh_data_type_query_allocation(mesh)
    _mesh_data_install!(mesh)
    _MESH_DATA_API.mesh.get_element_types()
    GC.gc()
    return @allocated _MESH_DATA_API.mesh.get_element_types()
end

@noinline function _mesh_data_element_query_allocation(mesh)
    _mesh_data_install!(mesh)
    _MESH_DATA_API.mesh.get_elements()
    GC.gc()
    return @allocated _MESH_DATA_API.mesh.get_elements()
end

@noinline function _mesh_data_nodes_by_type_allocation(mesh)
    _mesh_data_install!(mesh)
    _MESH_DATA_API.mesh.get_nodes_by_element_type(1)
    GC.gc()
    return @allocated _MESH_DATA_API.mesh.get_nodes_by_element_type(1)
end

@testset "detached bulk simplex data through API" begin
    _MESH_DATA_API.finalize()
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_nodes()
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_elements()
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_element_types()
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_elements_by_type(4)
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_nodes_by_element_type(4)
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_barycenters(4,-1,false,false)
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_element_edge_nodes(4)
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_element_face_nodes(4,3)
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_max_node_tag()
    @test_throws ArgumentError _MESH_DATA_API.mesh.get_max_element_tag()

    try
        _MESH_DATA_API.initialize()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_nodes()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_elements()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_element_types()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_elements_by_type(4)
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_nodes_by_element_type(4)
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_barycenters(4,-1,false,false)
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_element_edge_nodes(4)
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_element_face_nodes(4,3)
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_max_node_tag()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_max_element_tag()

        source=_mesh_data_fixture()
        source_crc=mesh_crc(source)
        @test source_crc.sha==
              "2c212be598e39946ad07859e38a85cec014d61839ba481009854bd5dd605c130"
        _mesh_data_install!(source)

        node_tags,coordinates,parameters=_MESH_DATA_API.mesh.get_nodes()
        @test node_tags==UInt64[1,2,3,4]
        @test coordinates==collect(vec(source.coords))
        @test parameters==Float64[]
        @test _MESH_DATA_API.mesh.get_nodes(-1,-2,true,false)==
              (UInt64[1,2,3,4],collect(vec(source.coords)),Float64[])
        node_tags[1]=99
        coordinates[1]=99
        push!(parameters,99)
        @test _MESH_DATA_API.mesh.get_nodes()==
              (UInt64[1,2,3,4],collect(vec(source.coords)),Float64[])

        for call in (
            ()->_MESH_DATA_API.mesh.get_nodes(false),
            ()->_MESH_DATA_API.mesh.get_nodes(-1,false),
            ()->_MESH_DATA_API.mesh.get_nodes(-2),
            ()->_MESH_DATA_API.mesh.get_nodes(4),
            ()->_MESH_DATA_API.mesh.get_nodes(3,-1),
            ()->_MESH_DATA_API.mesh.get_nodes(-1,1),
            ()->_MESH_DATA_API.mesh.get_nodes(-1,-1,0),
            ()->_MESH_DATA_API.mesh.get_nodes(-1,-1,false,missing),
            ()->_MESH_DATA_API.mesh.get_nodes(big(typemax(Int))+1),
        )
            @test_throws ArgumentError call()
            @test mesh_crc(_MESH_DATA_API.mesh.get())==source_crc
        end

        element_types,element_tags,element_nodes=
            _MESH_DATA_API.mesh.get_elements()
        @test element_types==Int32[1,2,4]
        @test element_tags==[UInt64[1],UInt64[2],UInt64[3]]
        @test element_nodes==[
            UInt64[1,2],UInt64[1,2,3],UInt64[1,2,3,4]]
        @test _MESH_DATA_API.mesh.get_elements(0,-1)==
              (Int32[],Vector{UInt64}[],Vector{UInt64}[])
        @test _MESH_DATA_API.mesh.get_elements(1,-2)==
              (Int32[1],[UInt64[1]],[UInt64[1,2]])
        @test _MESH_DATA_API.mesh.get_elements(2,-1)==
              (Int32[2],[UInt64[2]],[UInt64[1,2,3]])
        @test _MESH_DATA_API.mesh.get_elements(3,-1)==
              (Int32[4],[UInt64[3]],[UInt64[1,2,3,4]])

        element_types[1]=99
        element_tags[1][1]=99
        element_nodes[1][1]=99
        @test _MESH_DATA_API.mesh.get_elements()==
              (Int32[1,2,4],[UInt64[1],UInt64[2],UInt64[3]],
               [UInt64[1,2],UInt64[1,2,3],UInt64[1,2,3,4]])

        @test _MESH_DATA_API.mesh.get_element_types()==Int32[1,2,4]
        @test _MESH_DATA_API.mesh.get_element_types(0)==Int32[]
        @test _MESH_DATA_API.mesh.get_element_types(1)==Int32[1]
        @test _MESH_DATA_API.mesh.get_element_types(2)==Int32[2]
        @test _MESH_DATA_API.mesh.get_element_types(3)==Int32[4]
        @test _MESH_DATA_API.mesh.get_elements_by_type(1)==
              (UInt64[1],UInt64[1,2])
        @test _MESH_DATA_API.mesh.get_elements_by_type(2)==
              (UInt64[2],UInt64[1,2,3])
        @test _MESH_DATA_API.mesh.get_elements_by_type(4)==
              (UInt64[3],UInt64[1,2,3,4])
        @test _MESH_DATA_API.mesh.get_elements_by_type(3)==
              (UInt64[],UInt64[])
        @test _MESH_DATA_API.mesh.get_elements_by_type(15)==
              (UInt64[],UInt64[])

        segment_nodes=_MESH_DATA_API.mesh.get_nodes_by_element_type(1)
        triangle_nodes=_MESH_DATA_API.mesh.get_nodes_by_element_type(2)
        tetrahedron_nodes=_MESH_DATA_API.mesh.get_nodes_by_element_type(4)
        @test segment_nodes==
              (UInt64[1,2],Float64[0,0,0,1,0,0],Float64[])
        @test triangle_nodes==
              (UInt64[1,2,3],Float64[0,0,0,1,0,0,0,1,0],Float64[])
        @test tetrahedron_nodes==
              (UInt64[1,2,3,4],collect(vec(source.coords)),Float64[])
        @test _MESH_DATA_API.mesh.get_nodes_by_element_type(4,-2,false)==
              tetrahedron_nodes
        @test _MESH_DATA_API.mesh.get_nodes_by_element_type(3)==
              (UInt64[],Float64[],Float64[])

        segment_barycenter=_MESH_DATA_API.mesh.get_barycenters(1,-1,false,false)
        triangle_barycenter=_MESH_DATA_API.mesh.get_barycenters(2,-1,false,true)
        tetrahedron_barycenter=
            _MESH_DATA_API.mesh.get_barycenters(4,-1,false,false)
        @test segment_barycenter==Float64[0.5,0,0]
        @test triangle_barycenter==Float64[1/3,1/3,0]
        @test tetrahedron_barycenter==Float64[0.25,0.25,0.25]
        @test _MESH_DATA_API.mesh.get_barycenters(1,-1,true,true)==
              Float64[1,0,0]
        @test _MESH_DATA_API.mesh.get_barycenters(2,-1,true,false)==
              Float64[1,1,0]
        @test _MESH_DATA_API.mesh.get_barycenters(4,-1,true,true)==
              Float64[1,1,1]
        @test isempty(_MESH_DATA_API.mesh.get_barycenters(3,-1,false,false))

        segment_edges=_MESH_DATA_API.mesh.get_element_edge_nodes(1)
        triangle_edges=_MESH_DATA_API.mesh.get_element_edge_nodes(2,-2,true)
        tetrahedron_edges=_MESH_DATA_API.mesh.get_element_edge_nodes(4)
        @test segment_edges==UInt64[1,2]
        @test triangle_edges==UInt64[1,2,2,3,3,1]
        @test tetrahedron_edges==UInt64[1,2,2,3,3,1,4,1,4,3,4,2]
        @test isempty(_MESH_DATA_API.mesh.get_element_edge_nodes(3))

        triangle_faces=_MESH_DATA_API.mesh.get_element_face_nodes(2,3)
        tetrahedron_faces=_MESH_DATA_API.mesh.get_element_face_nodes(4,3,-2,true)
        @test triangle_faces==UInt64[1,2,3]
        @test tetrahedron_faces==UInt64[1,3,2,1,2,4,1,4,3,4,2,3]
        @test isempty(_MESH_DATA_API.mesh.get_element_face_nodes(1,3))
        @test isempty(_MESH_DATA_API.mesh.get_element_face_nodes(4,4))

        query_sha=_mesh_data_query_sha((
            segment_nodes[1],segment_nodes[2],triangle_nodes[1],triangle_nodes[2],
            tetrahedron_nodes[1],tetrahedron_nodes[2],segment_barycenter,
            triangle_barycenter,tetrahedron_barycenter,segment_edges,
            triangle_edges,tetrahedron_edges,triangle_faces,tetrahedron_faces))
        @test query_sha==
              "04e09b72ebf17bdc7ab2f9f96da2927c5a6e892e5c2d98c9ddbb8313dc4cab13"

        segment_nodes[1][1]=99
        segment_nodes[2][1]=99
        segment_barycenter[1]=99
        segment_edges[1]=99
        triangle_faces[1]=99
        @test _MESH_DATA_API.mesh.get_nodes_by_element_type(1)==
              (UInt64[1,2],Float64[0,0,0,1,0,0],Float64[])
        @test _MESH_DATA_API.mesh.get_barycenters(1,-1,false,false)==
              Float64[0.5,0,0]
        @test _MESH_DATA_API.mesh.get_element_edge_nodes(1)==UInt64[1,2]
        @test _MESH_DATA_API.mesh.get_element_face_nodes(2,3)==UInt64[1,2,3]
        @test _MESH_DATA_API.mesh.get_max_node_tag()==UInt64(4)
        @test _MESH_DATA_API.mesh.get_max_element_tag()==UInt64(3)

        for call in (
            ()->_MESH_DATA_API.mesh.get_elements(false),
            ()->_MESH_DATA_API.mesh.get_elements(-2),
            ()->_MESH_DATA_API.mesh.get_elements(4),
            ()->_MESH_DATA_API.mesh.get_elements(3,1),
            ()->_MESH_DATA_API.mesh.get_element_types(3,1),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(true),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(34),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(999),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,1),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,1,2),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,0,2),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,-1,1),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,0,0),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,true,1),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,0,false),
            ()->_MESH_DATA_API.mesh.get_nodes_by_element_type(true),
            ()->_MESH_DATA_API.mesh.get_nodes_by_element_type(34),
            ()->_MESH_DATA_API.mesh.get_nodes_by_element_type(999),
            ()->_MESH_DATA_API.mesh.get_nodes_by_element_type(4,1),
            ()->_MESH_DATA_API.mesh.get_nodes_by_element_type(4,-1,0),
            ()->_MESH_DATA_API.mesh.get_barycenters(4,1,false,false),
            ()->_MESH_DATA_API.mesh.get_barycenters(4,-1,0,false),
            ()->_MESH_DATA_API.mesh.get_barycenters(4,-1,false,0),
            ()->_MESH_DATA_API.mesh.get_barycenters(4,-1,false,false,1,2),
            ()->_MESH_DATA_API.mesh.get_element_edge_nodes(4,1),
            ()->_MESH_DATA_API.mesh.get_element_edge_nodes(4,-1,0),
            ()->_MESH_DATA_API.mesh.get_element_edge_nodes(4,-1,false,0,2),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,true),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,2),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,5),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,3,1),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,3,-1,0),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,3,-1,false,1,2),
        )
            @test_throws ArgumentError call()
            @test mesh_crc(_MESH_DATA_API.mesh.get())==source_crc
        end

        maximum=floatmax(Float64)
        overflow_fixture=Mesh(
            Float64[maximum maximum;0 0;0 0];
            segs=reshape(Int32[1,2],2,1))
        _mesh_data_install!(overflow_fixture)
        @test _MESH_DATA_API.mesh.get_barycenters(1,-1,false,false)==
              Float64[maximum,0,0]
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_barycenters(
            1,-1,true,false)
        @test _MESH_DATA_API.mesh.get_nodes_by_element_type(1)[2]==
              collect(vec(overflow_fixture.coords))
        _mesh_data_install!(source)

        translation=(1.0,0.0,0.0,2.0,
                     0.0,1.0,0.0,3.0,
                     0.0,0.0,1.0,4.0)
        moved=_MESH_DATA_API.mesh.affine_transform(translation)
        moved_nodes,moved_coordinates,_=_MESH_DATA_API.mesh.get_nodes()
        @test moved_nodes==UInt64[1,2,3,4]
        @test reshape(moved_coordinates,3,:)==moved.coords
        @test moved.coords==source.coords .+ [2.0,3.0,4.0]
        @test _MESH_DATA_API.mesh.get_elements()==
              (Int32[1,2,4],[UInt64[1],UInt64[2],UInt64[3]],
               [UInt64[1,2],UInt64[1,2,3],UInt64[1,2,3,4]])

        refined=_MESH_DATA_API.mesh.refine()
        refined_crc=mesh_crc(refined)
        @test refined_crc.sha==
              "db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031"
        @test refined_crc.bbox==((2.0,3.0,4.0),(3.0,4.0,5.0))
        @test _MESH_DATA_API.mesh.get_max_node_tag()==UInt64(10)
        @test _MESH_DATA_API.mesh.get_max_element_tag()==UInt64(14)
        refined_types,refined_tags,refined_nodes=
            _MESH_DATA_API.mesh.get_elements()
        @test refined_types==Int32[1,2,4]
        @test refined_tags==[
            UInt64[1,2],UInt64[3,4,5,6],UInt64[7,8,9,10,11,12,13,14]]
        @test length.(refined_nodes)==[4,12,32]
        for (element_type,cells) in
            ((1,refined.segs),(2,refined.tris),(4,refined.tets))
            type_tags,type_coordinates,_=
                _MESH_DATA_API.mesh.get_nodes_by_element_type(element_type)
            @test type_tags==UInt64.(vec(cells))
            @test reshape(type_coordinates,3,:)==
                  refined.coords[:,Int.(type_tags)]
        end
        @test length(_MESH_DATA_API.mesh.get_barycenters(4,-1,false,false))==24
        @test length(_MESH_DATA_API.mesh.get_element_edge_nodes(4))==96
        @test length(_MESH_DATA_API.mesh.get_element_face_nodes(4,3))==96

        @test _MESH_DATA_API.mesh.clear()===nothing
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_nodes()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_elements()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_max_node_tag()

        allocation_small=_mesh_data_segment_fixture(5_000)
        allocation_large=_mesh_data_segment_fixture(10_000)
        type_small=_mesh_data_type_query_allocation(allocation_small)
        type_large=_mesh_data_type_query_allocation(allocation_large)
        @test type_small>0
        @test type_large<=type_small+1_024
        elements_small=_mesh_data_element_query_allocation(allocation_small)
        elements_large=_mesh_data_element_query_allocation(allocation_large)
        @test elements_small>0
        @test elements_large>elements_small
        @test elements_large<=2.2elements_small+65_536
        nodes_by_type_small=
            _mesh_data_nodes_by_type_allocation(allocation_small)
        nodes_by_type_large=
            _mesh_data_nodes_by_type_allocation(allocation_large)
        @test nodes_by_type_small>0
        @test nodes_by_type_large>nodes_by_type_small
        @test nodes_by_type_large<=2.2nodes_by_type_small+65_536

        _mesh_data_install!(Mesh(zeros(3,0)))
        @test _MESH_DATA_API.mesh.get_nodes()==
              (UInt64[],Float64[],Float64[])
        @test _MESH_DATA_API.mesh.get_elements()==
              (Int32[],Vector{UInt64}[],Vector{UInt64}[])
        @test _MESH_DATA_API.mesh.get_nodes_by_element_type(4)==
              (UInt64[],Float64[],Float64[])
        @test isempty(_MESH_DATA_API.mesh.get_barycenters(4,-1,false,false))
        @test isempty(_MESH_DATA_API.mesh.get_element_edge_nodes(4))
        @test isempty(_MESH_DATA_API.mesh.get_element_face_nodes(4,3))
        @test _MESH_DATA_API.mesh.get_max_node_tag()==UInt64(0)
        @test _MESH_DATA_API.mesh.get_max_element_tag()==UInt64(0)
    finally
        _MESH_DATA_API.finalize()
    end
end
