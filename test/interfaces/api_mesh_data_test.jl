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

function _mesh_data_triangle_strip()
    coordinates=Float64[0 1 2 0 1 2;
                        0 0 0 1 1 1;
                        0 0 0 0 0 0]
    triangles=Int32[1 2 2 3;
                    2 5 3 6;
                    4 4 5 5]
    return Mesh(coordinates;tris=triangles)
end

@noinline function _mesh_data_partitioned_elements(mesh,task,num_tasks)
    _mesh_data_install!(mesh)
    GC.gc()
    return @allocated _MESH_DATA_API.mesh.get_elements_by_type(
        1,-1,task,num_tasks)
end

@noinline function _mesh_data_partitioned_barycenters(mesh,task,num_tasks)
    _mesh_data_install!(mesh)
    GC.gc()
    return @allocated _MESH_DATA_API.mesh.get_barycenters(
        1,-1,false,false,task,num_tasks)
end

@noinline function _mesh_data_partitioned_edge_nodes(mesh,task,num_tasks)
    _mesh_data_install!(mesh)
    GC.gc()
    return @allocated _MESH_DATA_API.mesh.get_element_edge_nodes(
        1,-1,false,task,num_tasks)
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

        # dim=-1 ignores `tag` entirely (verified against Gmsh 4.15.2).
        @test _MESH_DATA_API.mesh.get_nodes(-1,1)==
              _MESH_DATA_API.mesh.get_nodes(-1,-1)
        for call in (
            ()->_MESH_DATA_API.mesh.get_nodes(false),
            ()->_MESH_DATA_API.mesh.get_nodes(-1,false),
            ()->_MESH_DATA_API.mesh.get_nodes(-2),
            ()->_MESH_DATA_API.mesh.get_nodes(4),
            ()->_MESH_DATA_API.mesh.get_nodes(3,-1),
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
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,-1,1),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,0,0),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,0,-1),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(4,-1,1,0),
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
            ()->_MESH_DATA_API.mesh.get_barycenters(4,-1,false,false,-1,1),
            ()->_MESH_DATA_API.mesh.get_barycenters(4,-1,false,false,0,0),
            ()->_MESH_DATA_API.mesh.get_element_edge_nodes(4,1),
            ()->_MESH_DATA_API.mesh.get_element_edge_nodes(4,-1,0),
            ()->_MESH_DATA_API.mesh.get_element_edge_nodes(4,-1,false,-1,1),
            ()->_MESH_DATA_API.mesh.get_element_edge_nodes(4,-1,false,0,0),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,true),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,2),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,5),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,3,1),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,3,-1,0),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,3,-1,false,-1,1),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(4,3,-1,false,0,0),
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

@testset "deterministic task partitioning of detached queries" begin
    _MESH_DATA_API.finalize()
    try
        _MESH_DATA_API.initialize()
        chain=_mesh_data_segment_fixture(5)
        chain_crc=mesh_crc(chain)
        _mesh_data_install!(chain)

        # Gmsh contiguous blocks for five segments: (0,2)->1:2, (1,2)->3:5,
        # (0,3)->1:1, (1,3)->2:3, (2,3)->4:5. A strided i%num_tasks==task
        # scheme would return positions {2,5} for task 1 of 3 instead of
        # {2,3}, so these exact slices pin the documented formula.
        @test _MESH_DATA_API.mesh.get_elements_by_type(1,-1,0,2)==
              (UInt64[1,2],UInt64[1,2,2,3])
        @test _MESH_DATA_API.mesh.get_elements_by_type(1,-1,1,2)==
              (UInt64[3,4,5],UInt64[3,4,4,5,5,6])
        @test _MESH_DATA_API.mesh.get_elements_by_type(1,-1,0,3)==
              (UInt64[1],UInt64[1,2])
        @test _MESH_DATA_API.mesh.get_elements_by_type(1,-1,1,3)==
              (UInt64[2,3],UInt64[2,3,3,4])
        @test _MESH_DATA_API.mesh.get_elements_by_type(1,-1,2,3)==
              (UInt64[4,5],UInt64[4,5,5,6])
        # Adversarial task counts must not wrap before the truncating
        # division: (typemax(Int)-1)*5 overflows Int64, yet the block is 5:5.
        @test _MESH_DATA_API.mesh.get_elements_by_type(
            1,-1,typemax(Int)-1,typemax(Int))==
              (UInt64[5],UInt64[5,6])
        # task>=num_tasks is the silently-empty Gmsh range, never an error.
        for (task,num_tasks) in ((2,2),(3,3),(5,5),(0,7),(7,3))
            @test _MESH_DATA_API.mesh.get_elements_by_type(
                1,-1,task,num_tasks)==(UInt64[],UInt64[])
            @test isempty(_MESH_DATA_API.mesh.get_barycenters(
                1,-1,false,false,task,num_tasks))
            @test isempty(_MESH_DATA_API.mesh.get_element_edge_nodes(
                1,-1,false,task,num_tasks))
            @test isempty(_MESH_DATA_API.mesh.get_element_face_nodes(
                2,3,-1,false,task,num_tasks))
        end
        @test _MESH_DATA_API.mesh.get_elements_by_type(3,-1,1,2)==
              (UInt64[],UInt64[])

        @test _MESH_DATA_API.mesh.get_barycenters(1,-1,false,false,1,2)==
              Float64[2.5,0,0,3.5,0,0,4.5,0,0]
        @test _MESH_DATA_API.mesh.get_barycenters(1,-1,false,false,1,3)==
              Float64[1.5,0,0,2.5,0,0]
        @test _MESH_DATA_API.mesh.get_barycenters(1,-1,true,false,2,3)==
              Float64[7,0,0,9,0,0]
        @test _MESH_DATA_API.mesh.get_element_edge_nodes(1,-1,false,1,3)==
              UInt64[2,3,3,4]
        # Segments own no triangular faces: the slice stays empty, not an
        # error, exactly like the complete query.
        @test isempty(
            _MESH_DATA_API.mesh.get_element_face_nodes(1,3,-1,false,1,2))

        # Union over tasks reproduces every complete result bit-for-bit.
        full_tags,full_nodes=_MESH_DATA_API.mesh.get_elements_by_type(1)
        full_barycenters=_MESH_DATA_API.mesh.get_barycenters(
            1,-1,false,false)
        full_edges=_MESH_DATA_API.mesh.get_element_edge_nodes(1)
        for num_tasks in (1,2,3,5,6)
            union_tags=UInt64[]
            union_nodes=UInt64[]
            union_barycenters=Float64[]
            union_edges=UInt64[]
            for task in 0:num_tasks-1
                slice_tags,slice_nodes=
                    _MESH_DATA_API.mesh.get_elements_by_type(
                        1,-1,task,num_tasks)
                append!(union_tags,slice_tags)
                append!(union_nodes,slice_nodes)
                append!(union_barycenters,
                        _MESH_DATA_API.mesh.get_barycenters(
                            1,-1,false,false,task,num_tasks))
                append!(union_edges,
                        _MESH_DATA_API.mesh.get_element_edge_nodes(
                            1,-1,false,task,num_tasks))
            end
            @test union_tags==full_tags
            @test union_nodes==full_nodes
            @test union_barycenters==full_barycenters
            @test union_edges==full_edges
        end

        # Slices are detached: mutating one never leaks into the cache.
        slice_tags,slice_nodes=_MESH_DATA_API.mesh.get_elements_by_type(
            1,-1,1,2)
        slice_tags[1]=99
        slice_nodes[1]=99
        @test _MESH_DATA_API.mesh.get_elements_by_type(1,-1,1,2)==
              (UInt64[3,4,5],UInt64[3,4,4,5,5,6])
        @test mesh_crc(_MESH_DATA_API.mesh.get())==chain_crc

        strip=_mesh_data_triangle_strip()
        strip_crc=mesh_crc(strip)
        _mesh_data_install!(strip)
        # Four triangles: (0,3)->1:1, (1,3)->2:2, (2,3)->3:4.
        @test _MESH_DATA_API.mesh.get_elements_by_type(2,-1,2,3)==
              (UInt64[3,4],UInt64[2,3,5,3,6,5])
        @test _MESH_DATA_API.mesh.get_barycenters(2,-1,false,false,1,3)==
              Float64[2/3,2/3,0]
        @test _MESH_DATA_API.mesh.get_element_edge_nodes(2,-1,false,2,3)==
              UInt64[2,3,3,5,5,2,3,6,6,5,5,3]
        @test _MESH_DATA_API.mesh.get_element_face_nodes(2,3,-1,false,0,3)==
              UInt64[1,2,4]
        @test _MESH_DATA_API.mesh.get_element_face_nodes(2,3,-1,false,2,3)==
              UInt64[2,3,5,3,6,5]
        strip_tags,strip_nodes=_MESH_DATA_API.mesh.get_elements_by_type(2)
        strip_faces=_MESH_DATA_API.mesh.get_element_face_nodes(2,3)
        union_tags=UInt64[]
        union_nodes=UInt64[]
        union_faces=UInt64[]
        for task in 0:2
            slice_tags,slice_nodes=_MESH_DATA_API.mesh.get_elements_by_type(
                2,-1,task,3)
            append!(union_tags,slice_tags)
            append!(union_nodes,slice_nodes)
            append!(union_faces,_MESH_DATA_API.mesh.get_element_face_nodes(
                2,3,-1,false,task,3))
        end
        @test union_tags==strip_tags
        @test union_nodes==strip_nodes
        @test union_faces==strip_faces
        @test mesh_crc(_MESH_DATA_API.mesh.get())==strip_crc

        # Partitions of an empty cache stay empty and never allocate the
        # complete result.
        _mesh_data_install!(Mesh(zeros(3,0)))
        @test _MESH_DATA_API.mesh.get_elements_by_type(4,-1,1,2)==
              (UInt64[],UInt64[])
        @test isempty(
            _MESH_DATA_API.mesh.get_barycenters(4,-1,false,false,1,2))
        @test isempty(
            _MESH_DATA_API.mesh.get_element_edge_nodes(4,-1,false,1,2))
        @test isempty(_MESH_DATA_API.mesh.get_element_face_nodes(
            4,3,-1,false,1,2))

        # Allocation scales with the slice: a half partition must stay well
        # below the complete query. A full-compute-then-slice patch would
        # allocate at least the complete result and fail this ratchet.
        wide=_mesh_data_segment_fixture(10_000)
        full_elements=_mesh_data_partitioned_elements(wide,0,1)
        part_elements=_mesh_data_partitioned_elements(wide,1,2)
        @test full_elements>0
        @test 0<part_elements<full_elements
        @test 4part_elements<=3full_elements
        full_bary=_mesh_data_partitioned_barycenters(wide,0,1)
        part_bary=_mesh_data_partitioned_barycenters(wide,1,2)
        @test full_bary>0
        @test 0<part_bary<full_bary
        @test 4part_bary<=3full_bary
        full_edge=_mesh_data_partitioned_edge_nodes(wide,0,1)
        part_edge=_mesh_data_partitioned_edge_nodes(wide,1,2)
        @test full_edge>0
        @test 0<part_edge<full_edge
        @test 4part_edge<=3full_edge
    finally
        _MESH_DATA_API.finalize()
    end
end

@testset "entity-filtered mesh queries on generated caches" begin
    _MESH_DATA_API.initialize()
    try
        for (point,(x,y)) in enumerate(
                ((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)))
            _MESH_DATA_API.model.add_point(x,y,0.0;tag=point)
        end
        for (curve,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            _MESH_DATA_API.model.add_line(a,b;tag=curve)
        end
        _MESH_DATA_API.model.add_curve_loop([1,2,3,4];tag=1)
        _MESH_DATA_API.model.add_plane_surface([1];tag=1)
        _MESH_DATA_API.option("Mesh.MeshSizeMax",0.35)
        _MESH_DATA_API.mesh.generate(2)

        all_nodes=_MESH_DATA_API.mesh.get_nodes(-1,-1)[1]
        surface_nodes=_MESH_DATA_API.mesh.get_nodes(2,1)[1]
        surface_with_boundary=_MESH_DATA_API.mesh.get_nodes(2,1,true)[1]
        @test issorted(surface_nodes)
        @test Set(surface_nodes)⊆Set(all_nodes)
        @test surface_with_boundary[1:length(surface_nodes)]==surface_nodes
        # The include_boundary closure is the transitive boundary node union.
        expected=Set(surface_nodes)
        boundary_entities=_MESH_DATA_API.model.get_boundary(
            [(2,1)],false,false,false)
        queue=collect(boundary_entities)
        seen=Set{Tuple{Int,Int}}()
        while !isempty(queue)
            key=popfirst!(queue)
            key in seen && continue
            push!(seen,key)
            union!(expected,
                   _MESH_DATA_API.mesh.get_nodes(key[1],key[2])[1])
            append!(queue,_MESH_DATA_API.model.get_boundary(
                [key],false,false,false))
        end
        @test Set(surface_with_boundary)==expected
        # A curve query covers its own nodes; boundary adds its two points.
        curve_nodes=_MESH_DATA_API.mesh.get_nodes(1,1)[1]
        curve_with_boundary=_MESH_DATA_API.mesh.get_nodes(1,1,true)[1]
        @test curve_with_boundary[1:length(curve_nodes)]==curve_nodes
        @test Set(curve_with_boundary[length(curve_nodes)+1:end])==
              Set(_MESH_DATA_API.mesh.get_nodes(0,1)[1])∪
              Set(_MESH_DATA_API.mesh.get_nodes(0,2)[1])
        @test length(_MESH_DATA_API.mesh.get_nodes(0,1)[1])==1
        # Dimension-wide queries concatenate per-entity emissions.
        dim2=_MESH_DATA_API.mesh.get_nodes(2,-1)[1]
        @test Set(dim2)==Set(surface_nodes)
        dim2_with_boundary=_MESH_DATA_API.mesh.get_nodes(2,-1,true)[1]
        @test Set(dim2_with_boundary)==expected
        # Parametric coordinates ride the queried entity's parametrization.
        surf_tags,surf_coords,surf_par=
            _MESH_DATA_API.mesh.get_nodes(2,1,true,true)
        @test length(surf_par)==2length(surf_tags)
        for index in eachindex(surf_tags)
            evaluated=_MESH_DATA_API.model.get_value(
                2,1,surf_par[2index-1:2index])
            @test evaluated≈surf_coords[3index-2:3index]
        end
        curve_tags,curve_coords,curve_par=
            _MESH_DATA_API.mesh.get_nodes(1,1,true,true)
        @test length(curve_par)==length(curve_tags)
        for index in eachindex(curve_tags)
            evaluated=_MESH_DATA_API.model.get_value(
                1,1,[curve_par[index]])
            @test evaluated≈curve_coords[3index-2:3index]
        end
        # Points, all-dimension queries, and opt-outs emit no parameters.
        @test isempty(_MESH_DATA_API.mesh.get_nodes(0,1,true,true)[3])
        @test isempty(_MESH_DATA_API.mesh.get_nodes(-1,-1,false,true)[3])
        @test isempty(_MESH_DATA_API.mesh.get_nodes(2,1,true,false)[3])
        # Per-element queries pack each node's parameters on its owning entity.
        cell_tags,cell_coords,cell_par=
            _MESH_DATA_API.mesh.get_nodes_by_element_type(2,-1,true)
        owners=Dict{UInt64,Tuple{Int,Int}}()
        for (dim,entity_tags) in
                (0=>(1,2,3,4),1=>(1,2,3,4),2=>(1,))
            for entity in entity_tags,
                node in _MESH_DATA_API.mesh.get_nodes(dim,entity)[1]
                owners[node]=(dim,entity)
            end
        end
        expected_width=sum(
            node->owners[node][1] in (1,2) ? owners[node][1] : 0,cell_tags)
        @test length(cell_par)==expected_width
        position=1
        for (k,node) in enumerate(cell_tags)
            (dim,entity)=owners[node]
            width=dim in (1,2) ? dim : 0
            if width>0
                evaluated=_MESH_DATA_API.model.get_value(
                    dim,entity,cell_par[position:position+width-1])
                @test evaluated≈cell_coords[3k-2:3k]
            end
            position+=width
        end
        # Element queries filter cells onto the queried entity.
        types,tags,nodes=_MESH_DATA_API.mesh.get_elements(2,1)
        @test types==Int32[2]
        all_types,all_tags,all_nodes=_MESH_DATA_API.mesh.get_elements(2,-1)
        @test all_types==types && all_tags==tags && all_nodes==nodes
        @test _MESH_DATA_API.mesh.get_element_types(2,1)==Int32[2]
        @test _MESH_DATA_API.mesh.get_element_types(1,1)==Int32[]
        triangle_count=length(tags[1])
        # Element-by-tag classification resolves block position and owner.
        for element_tag in tags[1]
            element_type,element_nodes,entity_dim,entity_tag=
                _MESH_DATA_API.mesh.get_element(element_tag)
            @test element_type==Int32(2)
            @test (entity_dim,entity_tag)==(2,1)
            @test length(element_nodes)==3
            @test Set(element_nodes)⊆Set(surface_with_boundary)
        end
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_element(0)
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_element(
            length(tags[1])+1)
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_element(1.5)
        # Type-funnel queries resolve the tag in the type's own dimension.
        filtered_tags,filtered_nodes=
            _MESH_DATA_API.mesh.get_elements_by_type(2,1)
        @test sort!(filtered_tags)==sort!(tags[1])
        @test length(filtered_nodes)==3triangle_count
        @test length(
            _MESH_DATA_API.mesh.get_nodes_by_element_type(2,1)[1])==
            3triangle_count
        @test length(_MESH_DATA_API.mesh.get_barycenters(
            2,1,false,false))==3triangle_count
        @test length(_MESH_DATA_API.mesh.get_element_edge_nodes(2,1))==
            6triangle_count
        @test length(_MESH_DATA_API.mesh.get_element_face_nodes(2,3,1))==
            3triangle_count
        jacobians,determinants,_=_MESH_DATA_API.mesh.get_jacobians(
            2,[0.25,0.25,0.0],1)
        @test length(determinants)==triangle_count
        @test all(>(0),determinants)
        @test length(jacobians)==9triangle_count
        @test _MESH_DATA_API.mesh.get_basis_functions_orientation(
            2,"Lagrange",1)==zeros(Int32,triangle_count)
        _,key_entities,_=_MESH_DATA_API.mesh.get_keys(2,"Lagrange",1)
        @test sort!(unique!(key_entities))==
              sort!(unique!(surface_with_boundary))
        # Task partitioning subdivides the entity-filtered subset.
        union_tags=UInt64[]
        for task in 0:2
            slice,_=_MESH_DATA_API.mesh.get_elements_by_type(2,1,task,3)
            append!(union_tags,slice)
        end
        @test sort!(union_tags)==sort!(tags[1])
        @test _MESH_DATA_API.mesh.get_elements_by_type(
            2,1,triangle_count,triangle_count)==(UInt64[],UInt64[])
        # Unknown entities and dimension mismatches fail explicitly.
        for call in (
            ()->_MESH_DATA_API.mesh.get_nodes(2,77),
            ()->_MESH_DATA_API.mesh.get_nodes(1,77),
            ()->_MESH_DATA_API.mesh.get_elements(2,77),
            ()->_MESH_DATA_API.mesh.get_element_types(2,77),
            ()->_MESH_DATA_API.mesh.get_elements_by_type(2,77),
            ()->_MESH_DATA_API.mesh.get_nodes_by_element_type(2,77),
            ()->_MESH_DATA_API.mesh.get_barycenters(2,77,false,false),
            ()->_MESH_DATA_API.mesh.get_element_edge_nodes(2,77),
            ()->_MESH_DATA_API.mesh.get_element_face_nodes(2,3,77),
            ()->_MESH_DATA_API.mesh.get_jacobians(2,[0.25,0.25,0.0],77),
            ()->_MESH_DATA_API.mesh.get_basis_functions_orientation(
                2,"Lagrange",77),
            ()->_MESH_DATA_API.mesh.get_keys(2,"Lagrange",77),
            # The triangle type resolves its tag in dimension 2, so a curve
            # tag is rejected rather than silently returning curve cells.
            ()->_MESH_DATA_API.mesh.get_elements_by_type(2,1+77),
            ()->_MESH_DATA_API.mesh.get_nodes(3,1),
        )
            @test_throws ArgumentError call()
        end
        # A curve tag on a dim-2 query fails even though the entity exists.
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_nodes(2,3)
        # Refinement and affine transforms keep the classification live.
        _MESH_DATA_API.mesh.refine()
        refined_surface=_MESH_DATA_API.mesh.get_nodes(2,1)[1]
        @test length(refined_surface)>length(surface_nodes)
        refined_types,refined_tags,_=_MESH_DATA_API.mesh.get_elements(2,1)
        @test refined_types==Int32[2]
        @test length(refined_tags[1])==4triangle_count
        _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 1.5; 0 1 0 0; 0 0 1 0; 0 0 0 1])
        @test _MESH_DATA_API.mesh.get_nodes(2,1)[1]==refined_surface
        _,shifted,_=_MESH_DATA_API.mesh.get_nodes(2,1)
        @test all(>=(1.5),shifted[1:3:end])
        # Clearing drops the cache and its classification together.
        _MESH_DATA_API.mesh.clear()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_nodes(2,1)
    finally
        _MESH_DATA_API.finalize()
    end
end

@testset "entity-selective mesh mutations on generated caches" begin
    _MESH_DATA_API.initialize()
    try
        for (point,(x,y)) in enumerate(
                ((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)))
            _MESH_DATA_API.model.add_point(x,y,0.0;tag=point)
        end
        for (curve,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            _MESH_DATA_API.model.add_line(a,b;tag=curve)
        end
        _MESH_DATA_API.model.add_curve_loop([1,2,3,4];tag=1)
        _MESH_DATA_API.model.add_plane_surface([1];tag=1)
        _MESH_DATA_API.option("Mesh.MeshSizeMax",0.35)
        _MESH_DATA_API.mesh.generate(2)
        _MESH_DATA_API.mesh.refine()

        # --- selective affine_transform moves only entity-owned nodes ---
        curve_nodes_before=_MESH_DATA_API.mesh.get_nodes(1,1)
        surface_nodes_before=_MESH_DATA_API.mesh.get_nodes(2,1)
        point_nodes_before=_MESH_DATA_API.mesh.get_nodes(0,1)
        _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 0; 0 1 0 0; 0 0 1 3.0; 0 0 0 1],[(1,1)])
        curve_nodes_after=_MESH_DATA_API.mesh.get_nodes(1,1)
        @test curve_nodes_after[1]==curve_nodes_before[1]
        @test all(==(3.0),curve_nodes_after[2][3:3:end])
        @test all(==(0.0),curve_nodes_before[2][3:3:end])
        # Unlisted entities keep their coordinates exactly.
        @test _MESH_DATA_API.mesh.get_nodes(2,1)==surface_nodes_before
        @test _MESH_DATA_API.mesh.get_nodes(0,1)==point_nodes_before
        # Multiple entities select their combined owned nodes.
        _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 0; 0 1 0 0; 0 0 1 -3.0; 0 0 0 1],[(1,1),(1,2)])
        @test all(==(0.0),_MESH_DATA_API.mesh.get_nodes(1,1)[2][3:3:end])
        @test all(==(-3.0),_MESH_DATA_API.mesh.get_nodes(1,2)[2][3:3:end])
        # Selection errors and malformed entries fail explicitly.
        @test_throws ArgumentError _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],[(1,99)])
        @test_throws ArgumentError _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],[(4,1)])
        @test_throws ArgumentError _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],[(1,true)])
        @test_throws ArgumentError _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],[1.5])

        # --- selective create_edges / create_faces ---
        _MESH_DATA_API.mesh.generate(2)
        _MESH_DATA_API.mesh.refine()
        _MESH_DATA_API.mesh.create_edges([(2,1)])
        selective_edges,_=_MESH_DATA_API.mesh.get_all_edges()
        @test !isempty(selective_edges)
        # Curves own no cells in the surface cache, so they add nothing.
        _MESH_DATA_API.mesh.create_edges([(1,1),(1,2),(1,3),(1,4)])
        @test _MESH_DATA_API.mesh.get_all_edges()[1]==selective_edges
        _MESH_DATA_API.mesh.create_edges()
        all_edges,_=_MESH_DATA_API.mesh.get_all_edges()
        @test length(all_edges)>=length(selective_edges)
        @test selective_edges⊆all_edges
        @test_throws ArgumentError _MESH_DATA_API.mesh.create_edges([(1,99)])
        @test_throws ArgumentError _MESH_DATA_API.mesh.create_edges([(2,99)])
        _MESH_DATA_API.mesh.create_faces([(2,1)])
        @test !isempty(_MESH_DATA_API.mesh.get_all_faces(3)[1])
        # Selecting entities whose cells lack triangular faces is a no-op.
        _MESH_DATA_API.mesh.create_faces([(1,1),(0,1)])
        @test_throws ArgumentError _MESH_DATA_API.mesh.create_faces([(2,99)])

        # --- selective clear ---
        _MESH_DATA_API.mesh.generate(2)
        _MESH_DATA_API.mesh.refine()
        total=_MESH_DATA_API.mesh.get_nodes()[1]
        # Clearing an entity that owns no cells is a no-op, as in Gmsh.
        _MESH_DATA_API.mesh.clear([(1,1),(0,1)])
        @test _MESH_DATA_API.mesh.get_nodes()[1]==total
        @test _MESH_DATA_API.mesh.get_nodes(1,1)[1] != Int[]
        # Clearing the generating entity drops its cells and owned nodes while
        # boundary-owned nodes survive under their own classification.
        surface_node_count=
            length(_MESH_DATA_API.mesh.get_nodes(2,1)[1])
        boundary_count=length(total)-surface_node_count
        _MESH_DATA_API.mesh.clear([(2,1)])
        remaining=_MESH_DATA_API.mesh.get_nodes()[1]
        @test length(remaining)==boundary_count
        @test _MESH_DATA_API.mesh.get_nodes(2,1)[1]==Int[]
        @test _MESH_DATA_API.mesh.get_nodes(1,1)[1] != Int[]
        types,_,_=_MESH_DATA_API.mesh.get_elements(2,1)
        @test types==Int32[]
        @test_throws ArgumentError _MESH_DATA_API.mesh.clear([(2,99)])
        # Empty selection still clears the complete cache.
        _MESH_DATA_API.mesh.clear()
        @test_throws ArgumentError _MESH_DATA_API.mesh.get_nodes(2,1)
    finally
        _MESH_DATA_API.finalize()
    end
end

@testset "entity-selective mutations on a classified volume cache" begin
    _MESH_DATA_API.initialize()
    try
        source="""
Point(1)={0,0,0,0.25};Point(2)={1,0,0,0.25};Point(3)={1,1,0,0.25};Point(4)={0,1,0,0.25};
Point(5)={0,0,1,0.25};Point(6)={1,0,1,0.25};Point(7)={1,1,1,0.25};Point(8)={0,1,1,0.25};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};
Curve Loop(3)={1,10,-5,-9};Curve Loop(4)={3,12,-7,-11};
Curve Loop(5)={2,11,-6,-10};Curve Loop(6)={4,9,-8,-12};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
Plane Surface(4)={4};Plane Surface(5)={5};Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};
"""
        mktemp() do path,io
            write(io,source)
            close(io)
            _MESH_DATA_API.open_geo!(path)
        end
        _MESH_DATA_API.mesh.generate(3)
        _MESH_DATA_API.mesh.refine()
        total=_MESH_DATA_API.mesh.get_nodes()[1]
        surface_nodes=_MESH_DATA_API.mesh.get_nodes(2,1)
        curve_nodes=_MESH_DATA_API.mesh.get_nodes(1,1)
        @test !isempty(surface_nodes[1]) && !isempty(curve_nodes[1])
        # A boundary entity owns no cells in the tet-only cache, so clearing it
        # is a no-op — matching Gmsh 4.15.2's retention of the boundary mesh.
        _MESH_DATA_API.mesh.clear([(2,1),(1,1),(0,1)])
        @test _MESH_DATA_API.mesh.get_nodes()[1]==total
        @test _MESH_DATA_API.mesh.get_nodes(2,1)[1]==surface_nodes[1]
        # Clearing the volume drops every tet and volume-owned node; boundary
        # nodes survive under their own entity classification — with fresh
        # dense tags after compaction, so compare coordinates.
        _MESH_DATA_API.mesh.clear([(3,1)])
        remaining=_MESH_DATA_API.mesh.get_nodes()[1]
        volume_nodes=length(total)-length(remaining)
        @test volume_nodes>0
        @test _MESH_DATA_API.mesh.get_elements(3,1)[1]==Int32[]
        @test sort(_MESH_DATA_API.mesh.get_nodes(2,1)[2][1:3:end])==
              sort(surface_nodes[2][1:3:end])
        @test sort(_MESH_DATA_API.mesh.get_nodes(1,1)[2][1:3:end])==
              sort(curve_nodes[2][1:3:end])
        @test _MESH_DATA_API.mesh.get_nodes(3,1)[1]==Int[]
        # A selective transform that inverts tets is rejected atomically.
        _MESH_DATA_API.mesh.generate(3)
        _MESH_DATA_API.mesh.refine()
        before=_MESH_DATA_API.mesh.get_nodes()
        @test_throws ArgumentError _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 0; 0 1 0 0; 0 0 1 50.0; 0 0 0 1],[(2,1)])
        @test _MESH_DATA_API.mesh.get_nodes()==before
        # A gentle selective transform moves only the surface-owned nodes.
        _,coords_before,_=_MESH_DATA_API.mesh.get_nodes(2,1)
        _MESH_DATA_API.mesh.affine_transform(
            [1.0 0 0 0; 0 1 0 0; 0 0 1 0.05; 0 0 0 1],[(2,1)])
        _,coords_after,_=_MESH_DATA_API.mesh.get_nodes(2,1)
        @test coords_after[3:3:end]==coords_before[3:3:end].+0.05
        _,other,_=_MESH_DATA_API.mesh.get_nodes(2,3)
        @test all(<=(1.0),other[3:3:end])
        @test _MESH_DATA_API.mesh.get_nodes()[1]==before[1]
    finally
        _MESH_DATA_API.finalize()
    end
end
