using Test
using SHA
using Random
using Tessella
using Tessella.MeshTypes: unique_edges, unique_faces

const _ENTITY_TOPOLOGY=Tessella.MeshEntityTopology

function _entity_topology_fixture()
    return Mesh(
        Float64[0 1 0 0 1;0 0 1 0 1;0 0 0 1 1];
        segs=Int32[1 3;2 1],
        tris=Int32[1 3;2 2;3 4],
        tets=Int32[1 1;2 3;3 2;4 5])
end

function _edge_set(nodes)
    result=Set{NTuple{2,Int32}}()
    for index in 1:2:length(nodes)
        a=Int32(nodes[index]);b=Int32(nodes[index+1])
        push!(result,a<b ? (a,b) : (b,a))
    end
    return result
end

function _face_set(nodes)
    result=Set{NTuple{3,Int32}}()
    for index in 1:3:length(nodes)
        values=sort!(Int32[nodes[index],nodes[index+1],nodes[index+2]])
        push!(result,(values[1],values[2],values[3]))
    end
    return result
end

function _topology_digest(mesh)
    edges=_ENTITY_TOPOLOGY._mesh_edge_topology(mesh)
    faces=_ENTITY_TOPOLOGY._mesh_face_topology(mesh)
    arrays=(
        _ENTITY_TOPOLOGY._mesh_all_edges(edges)...,
        _ENTITY_TOPOLOGY._mesh_edges(
            edges,mesh,UInt64[1,2,2,1,1,5,5,1])...,
        _ENTITY_TOPOLOGY._mesh_all_faces(faces,3)...,
        _ENTITY_TOPOLOGY._mesh_faces(
            faces,mesh,3,UInt64[1,2,3,3,2,1,1,2,4])...)
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

@noinline function _entity_topology_allocations(count::Int)
    mesh=Mesh(
        zeros(3,4count);
        tets=reshape(Int32.(1:4count),4,count))
    _ENTITY_TOPOLOGY._mesh_edge_topology(mesh)
    _ENTITY_TOPOLOGY._mesh_face_topology(mesh)
    GC.gc()
    edge_bytes=@allocated _ENTITY_TOPOLOGY._mesh_edge_topology(mesh)
    face_bytes=@allocated _ENTITY_TOPOLOGY._mesh_face_topology(mesh)
    return edge_bytes,face_bytes
end

@testset "deterministic global simplex entity topology" begin
    fixture=_entity_topology_fixture()
    edge_topology=_ENTITY_TOPOLOGY._mesh_edge_topology(fixture)
    face_topology=_ENTITY_TOPOLOGY._mesh_face_topology(fixture)

    edge_tags,edge_nodes=
        _ENTITY_TOPOLOGY._mesh_all_edges(edge_topology)
    @test edge_tags==UInt64.(1:9)
    @test edge_nodes==UInt64[
        1,2, 3,1, 2,3, 2,4, 4,3, 4,1, 5,1, 5,2, 5,3]
    @test _ENTITY_TOPOLOGY._mesh_edges(
        edge_topology,fixture,
        UInt64[1,2,2,1, 1,3,3,1, 2,4,4,2, 1,5,5,1])==
        (UInt64[1,1,2,2,4,4,7,7],Int32[1,-1,1,-1,1,-1,1,-1])

    face_tags,face_nodes=
        _ENTITY_TOPOLOGY._mesh_all_faces(face_topology,3)
    @test face_tags==UInt64.(1:7)
    @test face_nodes==UInt64[
        1,2,3, 3,2,4, 1,2,4, 1,4,3,
        1,3,5, 1,5,2, 5,3,2]
    permutations=UInt64[
        1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1]
    @test _ENTITY_TOPOLOGY._mesh_faces(
        face_topology,fixture,3,permutations)==
        (fill(UInt64(1),6),zeros(Int32,6))
    @test _ENTITY_TOPOLOGY._mesh_all_faces(face_topology,4)==
          (UInt64[],UInt64[])

    @test _edge_set(edge_nodes)==Set(unique_edges(fixture.tris,fixture.tets))
    @test _face_set(face_nodes)==Set(unique_faces(fixture.tris,fixture.tets))

    rng=Xoshiro(0x327bc6341dd4f9a2)
    for _ in 1:64
        triangle_nodes=randperm(rng,8)[1:3]
        tetrahedron_nodes=randperm(rng,8)[1:4]
        mesh=Mesh(
            zeros(3,8);
            tris=reshape(Int32.(triangle_nodes),3,1),
            tets=reshape(Int32.(tetrahedron_nodes),4,1))
        random_edges=_ENTITY_TOPOLOGY._mesh_edge_topology(mesh)
        random_faces=_ENTITY_TOPOLOGY._mesh_face_topology(mesh)
        @test _edge_set(
            _ENTITY_TOPOLOGY._mesh_all_edges(random_edges)[2])==
            Set(unique_edges(mesh.tris,mesh.tets))
        @test _face_set(
            _ENTITY_TOPOLOGY._mesh_all_faces(random_faces,3)[2])==
            Set(unique_faces(mesh.tris,mesh.tets))
    end

    detached_edges=_ENTITY_TOPOLOGY._mesh_all_edges(edge_topology)
    detached_faces=_ENTITY_TOPOLOGY._mesh_all_faces(face_topology,3)
    detached_edges[1][1]=99;detached_edges[2][1]=99
    detached_faces[1][1]=99;detached_faces[2][1]=99
    @test _ENTITY_TOPOLOGY._mesh_all_edges(edge_topology)[1]==UInt64.(1:9)
    @test _ENTITY_TOPOLOGY._mesh_all_faces(face_topology,3)[1]==UInt64.(1:7)

    @test _ENTITY_TOPOLOGY._mesh_edges(nothing,fixture,UInt64[])==
          (UInt64[],Int32[])
    @test _ENTITY_TOPOLOGY._mesh_faces(nothing,fixture,3,UInt64[])==
          (UInt64[],Int32[])
    @test _ENTITY_TOPOLOGY._mesh_all_edges(nothing)==
          (UInt64[],UInt64[])
    @test _ENTITY_TOPOLOGY._mesh_all_faces(nothing,3)==
          (UInt64[],UInt64[])

    for call in (
        ()->_ENTITY_TOPOLOGY._mesh_edges(nothing,fixture,UInt64[1,2]),
        ()->_ENTITY_TOPOLOGY._mesh_faces(nothing,fixture,3,UInt64[1,2,3]),
        ()->_ENTITY_TOPOLOGY._mesh_edges(edge_topology,fixture,1),
        ()->_ENTITY_TOPOLOGY._mesh_edges(edge_topology,fixture,UInt64[1]),
        ()->_ENTITY_TOPOLOGY._mesh_edges(edge_topology,fixture,Any[true,2]),
        ()->_ENTITY_TOPOLOGY._mesh_edges(edge_topology,fixture,Any[missing,2]),
        ()->_ENTITY_TOPOLOGY._mesh_edges(
            edge_topology,fixture,Any[big(typemax(Int))+1,2]),
        ()->_ENTITY_TOPOLOGY._mesh_edges(edge_topology,fixture,UInt64[0,1]),
        ()->_ENTITY_TOPOLOGY._mesh_edges(edge_topology,fixture,UInt64[1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_edges(edge_topology,fixture,UInt64[1,1]),
        ()->_ENTITY_TOPOLOGY._mesh_edges(edge_topology,fixture,UInt64[4,5]),
        ()->_ENTITY_TOPOLOGY._mesh_faces(face_topology,fixture,true,UInt64[]),
        ()->_ENTITY_TOPOLOGY._mesh_faces(
            face_topology,fixture,big(typemax(Int))+1,UInt64[]),
        ()->_ENTITY_TOPOLOGY._mesh_faces(face_topology,fixture,2,UInt64[]),
        ()->_ENTITY_TOPOLOGY._mesh_faces(face_topology,fixture,3,1),
        ()->_ENTITY_TOPOLOGY._mesh_faces(face_topology,fixture,3,UInt64[1,2]),
        ()->_ENTITY_TOPOLOGY._mesh_faces(
            face_topology,fixture,3,UInt64[1,1,2]),
        ()->_ENTITY_TOPOLOGY._mesh_faces(
            face_topology,fixture,3,UInt64[1,4,5]),
        ()->_ENTITY_TOPOLOGY._mesh_faces(
            face_topology,fixture,4,UInt64[1,2,3,4]),
        ()->_ENTITY_TOPOLOGY._mesh_all_faces(face_topology,5),
    )
        @test_throws ArgumentError call()
    end

    repeated_edge=Mesh(
        zeros(3,1);segs=reshape(Int32[1,1],2,1))
    repeated_face=Mesh(
        zeros(3,2);tris=reshape(Int32[1,2,1],3,1))
    @test_throws ArgumentError _ENTITY_TOPOLOGY._mesh_edge_topology(
        repeated_edge)
    @test_throws ArgumentError _ENTITY_TOPOLOGY._mesh_face_topology(
        repeated_face)

    small_edges,small_faces=_entity_topology_allocations(2_000)
    large_edges,large_faces=_entity_topology_allocations(4_000)
    @test small_edges>0
    @test small_faces>0
    @test large_edges<=2small_edges+262_144
    @test large_faces<=2small_faces+262_144

    @test _topology_digest(fixture)==
          "cf3ece3db81a0e6820c166ece9ab2a4fc08947dacd515039969a3b69b260a964"
    @test isempty(Docs.undocumented_names(
        Tessella.MeshEntityTopology;private=false))
    @test isempty(Test.detect_ambiguities(
        Tessella.MeshEntityTopology;recursive=true))
end
