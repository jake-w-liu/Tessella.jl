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

function _entity_manual_topology_fixture()
    fixture=_entity_topology_fixture()
    coordinates=hcat(fixture.coords,Float64[2,1,0])
    return Mesh(
        coordinates;segs=fixture.segs,tris=fixture.tris,tets=fixture.tets)
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

function _manual_topology_digest(mesh)
    edges=_ENTITY_TOPOLOGY._mesh_add_edges(
        nothing,mesh,UInt64[100,200],UInt64[4,5,2,5])
    edges=_ENTITY_TOPOLOGY._mesh_edge_topology(mesh,edges)
    edges=_ENTITY_TOPOLOGY._mesh_add_edges(
        edges,mesh,UInt64[300],UInt64[5,6])

    faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        nothing,mesh,3,UInt64[100,300],UInt64[1,4,5,1,2,5])
    faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        faces,mesh,4,UInt64[200],UInt64[1,2,4,5])
    faces=_ENTITY_TOPOLOGY._mesh_face_topology(mesh,faces)
    faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        faces,mesh,3,UInt64[400],UInt64[2,5,6])

    arrays=(
        _ENTITY_TOPOLOGY._mesh_all_edges(edges)...,
        _ENTITY_TOPOLOGY._mesh_edges(
            edges,mesh,UInt64[4,5,5,4,5,6,6,5])...,
        _ENTITY_TOPOLOGY._mesh_all_faces(faces,3)...,
        _ENTITY_TOPOLOGY._mesh_all_faces(faces,4)...,
        _ENTITY_TOPOLOGY._mesh_faces(
            faces,mesh,3,UInt64[1,4,5,5,4,1,2,5,6])...,
        _ENTITY_TOPOLOGY._mesh_faces(
            faces,mesh,4,UInt64[1,2,4,5,5,4,2,1])...)
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

@noinline function _manual_entity_topology_allocations(count::Int)
    mesh=Mesh(zeros(3,4count))
    edge_tags=UInt64.(1:count)
    edge_nodes=UInt64.(1:2count)
    face_tags=UInt64.(1:count)
    face_nodes=UInt64.(1:4count)
    _ENTITY_TOPOLOGY._mesh_add_edges(
        nothing,mesh,edge_tags,edge_nodes)
    _ENTITY_TOPOLOGY._mesh_add_faces(
        nothing,mesh,4,face_tags,face_nodes)
    GC.gc()
    edge_bytes=@allocated _ENTITY_TOPOLOGY._mesh_add_edges(
        nothing,mesh,edge_tags,edge_nodes)
    face_bytes=@allocated _ENTITY_TOPOLOGY._mesh_add_faces(
        nothing,mesh,4,face_tags,face_nodes)
    return edge_bytes,face_bytes
end

function _quadrangle_permutations(nodes::NTuple{4,UInt64})
    result=UInt64[]
    for a in nodes,b in nodes,c in nodes,d in nodes
        allunique((a,b,c,d)) || continue
        append!(result,(a,b,c,d))
    end
    return result
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

@testset "atomic manual edge and face topology" begin
    fixture=_entity_manual_topology_fixture()

    manual_edges=_ENTITY_TOPOLOGY._mesh_add_edges(
        nothing,fixture,UInt64[100,200],UInt64[4,5,2,5])
    @test _ENTITY_TOPOLOGY._mesh_all_edges(manual_edges)==
          (UInt64[100,200],UInt64[4,5,2,5])
    @test _ENTITY_TOPOLOGY._mesh_edges(
        manual_edges,fixture,UInt64[4,5,5,4,2,5,5,2])==
        (UInt64[100,100,200,200],Int32[1,-1,1,-1])
    duplicate_edges=_ENTITY_TOPOLOGY._mesh_add_edges(
        manual_edges,fixture,UInt64[100,200],UInt64[5,4,5,2])
    @test _ENTITY_TOPOLOGY._mesh_all_edges(duplicate_edges)==
          _ENTITY_TOPOLOGY._mesh_all_edges(manual_edges)

    completed_edges=_ENTITY_TOPOLOGY._mesh_edge_topology(
        fixture,manual_edges)
    @test _ENTITY_TOPOLOGY._mesh_all_edges(completed_edges)==
          (UInt64[3,4,5,6,7,8,9,10,100,200],UInt64[
              1,2, 3,1, 2,3, 2,4, 4,3, 4,1, 5,1, 5,3,
              4,5, 2,5])
    @test _ENTITY_TOPOLOGY._mesh_edge_topology(
        fixture,completed_edges) |> _ENTITY_TOPOLOGY._mesh_all_edges ==
        _ENTITY_TOPOLOGY._mesh_all_edges(completed_edges)
    extended_edges=_ENTITY_TOPOLOGY._mesh_add_edges(
        completed_edges,fixture,UInt64[300],UInt64[5,6])
    @test last(_ENTITY_TOPOLOGY._mesh_all_edges(extended_edges)[1])==300
    @test _ENTITY_TOPOLOGY._mesh_edges(
        extended_edges,fixture,UInt64[5,6,6,5])==
        (UInt64[300,300],Int32[1,-1])

    edge_baseline=_ENTITY_TOPOLOGY._mesh_all_edges(extended_edges)
    for call in (
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,1,UInt64[1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,UInt64[400],1),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,UInt64[400],UInt64[1]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,Any[true],UInt64[1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,Any[missing],UInt64[1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,Any[0],UInt64[1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,Any[-1],UInt64[1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,
            Any[big(typemax(UInt64))+1],UInt64[1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,UInt64[400],UInt64[1,7]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,UInt64[400],UInt64[6,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,UInt64[400],UInt64[4,5]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,UInt64[100],UInt64[1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_edges(
            extended_edges,fixture,UInt64[400,100],UInt64[1,6,3,6]),
    )
        @test_throws ArgumentError call()
        @test _ENTITY_TOPOLOGY._mesh_all_edges(extended_edges)==edge_baseline
    end
    maximum_edge=_ENTITY_TOPOLOGY._mesh_add_edges(
        nothing,fixture,UInt64[typemax(UInt64)],UInt64[1,6])
    @test _ENTITY_TOPOLOGY._mesh_edges(
        maximum_edge,fixture,UInt64[1,6])==
        (UInt64[typemax(UInt64)],Int32[1])
    @test _ENTITY_TOPOLOGY._mesh_all_edges(
        _ENTITY_TOPOLOGY._mesh_add_edges(
            nothing,fixture,UInt64[],UInt64[]))==(UInt64[],UInt64[])

    colliding_generated_edges=_ENTITY_TOPOLOGY._mesh_edge_topology(
        fixture,_ENTITY_TOPOLOGY._mesh_add_edges(
            nothing,fixture,UInt64[3],UInt64[5,6]))
    @test _ENTITY_TOPOLOGY._mesh_all_edges(
        colliding_generated_edges)[1]==UInt64.(2:11)
    @test _ENTITY_TOPOLOGY._mesh_edges(
        colliding_generated_edges,fixture,UInt64[5,6])==
        (UInt64[3],Int32[1])

    manual_faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        nothing,fixture,3,UInt64[100,300],
        UInt64[1,4,5,1,2,5])
    manual_faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        manual_faces,fixture,4,UInt64[200],UInt64[1,2,4,5])
    @test _ENTITY_TOPOLOGY._mesh_all_faces(manual_faces,3)==
          (UInt64[100,300],UInt64[1,4,5,1,2,5])
    @test _ENTITY_TOPOLOGY._mesh_all_faces(manual_faces,4)==
          (UInt64[200],UInt64[1,2,4,5])
    quadrangle_permutations=_quadrangle_permutations(
        (UInt64(1),UInt64(2),UInt64(4),UInt64(5)))
    @test length(quadrangle_permutations)==96
    @test _ENTITY_TOPOLOGY._mesh_faces(
        manual_faces,fixture,4,quadrangle_permutations)==
        (fill(UInt64(200),24),zeros(Int32,24))
    duplicate_faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        manual_faces,fixture,3,UInt64[100,300],
        UInt64[5,4,1,5,2,1])
    duplicate_faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        duplicate_faces,fixture,4,UInt64[200],UInt64[5,4,2,1])
    @test _ENTITY_TOPOLOGY._mesh_all_faces(duplicate_faces,3)==
          _ENTITY_TOPOLOGY._mesh_all_faces(manual_faces,3)
    @test _ENTITY_TOPOLOGY._mesh_all_faces(duplicate_faces,4)==
          _ENTITY_TOPOLOGY._mesh_all_faces(manual_faces,4)

    completed_faces=_ENTITY_TOPOLOGY._mesh_face_topology(
        fixture,manual_faces)
    @test _ENTITY_TOPOLOGY._mesh_all_faces(completed_faces,3)==
          (UInt64[4,5,6,7,8,9,100,300],UInt64[
              1,2,3, 3,2,4, 1,2,4, 1,4,3, 1,3,5, 5,3,2,
              1,4,5, 1,2,5])
    @test _ENTITY_TOPOLOGY._mesh_all_faces(completed_faces,4)==
          (UInt64[200],UInt64[1,2,4,5])
    recompleted_faces=_ENTITY_TOPOLOGY._mesh_face_topology(
        fixture,completed_faces)
    @test (
        _ENTITY_TOPOLOGY._mesh_all_faces(recompleted_faces,3),
        _ENTITY_TOPOLOGY._mesh_all_faces(recompleted_faces,4))==
        (_ENTITY_TOPOLOGY._mesh_all_faces(completed_faces,3),
         _ENTITY_TOPOLOGY._mesh_all_faces(completed_faces,4))
    extended_faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        completed_faces,fixture,3,UInt64[400],UInt64[2,5,6])
    @test _ENTITY_TOPOLOGY._mesh_faces(
        extended_faces,fixture,3,UInt64[2,5,6,6,5,2])==
        (UInt64[400,400],Int32[0,0])

    face_baseline=(
        _ENTITY_TOPOLOGY._mesh_all_faces(extended_faces,3),
        _ENTITY_TOPOLOGY._mesh_all_faces(extended_faces,4))
    for call in (
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,5,UInt64[],UInt64[]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,1,UInt64[1,2,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,UInt64[500],1),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,UInt64[500],UInt64[1,2]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,Any[true],UInt64[1,2,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,Any[missing],UInt64[1,2,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,Any[0],UInt64[1,2,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,Any[-1],UInt64[1,2,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,
            Any[big(typemax(UInt64))+1],UInt64[1,2,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,UInt64[500],UInt64[1,2,7]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,UInt64[500],UInt64[1,1,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,UInt64[500],UInt64[1,4,5]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,4,UInt64[100],UInt64[1,2,3,6]),
        ()->_ENTITY_TOPOLOGY._mesh_add_faces(
            extended_faces,fixture,3,UInt64[500,100],
            UInt64[1,2,6,2,4,6]),
    )
        @test_throws ArgumentError call()
        @test (
            _ENTITY_TOPOLOGY._mesh_all_faces(extended_faces,3),
            _ENTITY_TOPOLOGY._mesh_all_faces(extended_faces,4))==face_baseline
    end
    maximum_face=_ENTITY_TOPOLOGY._mesh_add_faces(
        nothing,fixture,4,UInt64[typemax(UInt64)],UInt64[1,2,3,6])
    @test _ENTITY_TOPOLOGY._mesh_faces(
        maximum_face,fixture,4,UInt64[6,3,2,1])==
        (UInt64[typemax(UInt64)],Int32[0])
    @test _ENTITY_TOPOLOGY._mesh_all_faces(
        _ENTITY_TOPOLOGY._mesh_add_faces(
            nothing,fixture,4,UInt64[],UInt64[]),4)==
        (UInt64[],UInt64[])

    colliding_generated_faces=_ENTITY_TOPOLOGY._mesh_face_topology(
        fixture,_ENTITY_TOPOLOGY._mesh_add_faces(
            nothing,fixture,4,UInt64[3],UInt64[1,2,4,6]))
    @test sort!(vcat(
        _ENTITY_TOPOLOGY._mesh_all_faces(
            colliding_generated_faces,3)[1],
        _ENTITY_TOPOLOGY._mesh_all_faces(
            colliding_generated_faces,4)[1]))==UInt64.(2:9)
    @test _ENTITY_TOPOLOGY._mesh_faces(
        colliding_generated_faces,fixture,4,UInt64[6,4,2,1])==
        (UInt64[3],Int32[0])

    rng=Xoshiro(0x6f53d1b446fb2ac8)
    random_mesh=Mesh(zeros(3,12))
    edge_keys=shuffle!(rng,[(a,b) for a in 1:12 for b in (a+1):12])[1:40]
    edge_identifiers=UInt64.(shuffle!(rng,collect(1001:1040)))
    edge_input=UInt64[]
    for (a,b) in edge_keys
        append!(edge_input,rand(rng,Bool) ? (a,b) : (b,a))
    end
    random_edges=_ENTITY_TOPOLOGY._mesh_add_edges(
        nothing,random_mesh,edge_identifiers,edge_input)
    ascending_edges=UInt64[]
    foreach(pair->append!(ascending_edges,pair),edge_keys)
    @test _ENTITY_TOPOLOGY._mesh_edges(
        random_edges,random_mesh,ascending_edges)==
        (edge_identifiers,ones(Int32,length(edge_identifiers)))
    @test _ENTITY_TOPOLOGY._mesh_all_edges(random_edges)[1]==
          sort(edge_identifiers)

    quadrangle_keys=shuffle!(rng,[
        (a,b,c,d) for a in 1:12 for b in (a+1):12
        for c in (b+1):12 for d in (c+1):12])[1:48]
    face_identifiers=UInt64.(shuffle!(rng,collect(2001:2048)))
    quadrangle_input=UInt64[]
    for key in quadrangle_keys
        order=randperm(rng,4)
        append!(quadrangle_input,ntuple(i->key[order[i]],4))
    end
    random_faces=_ENTITY_TOPOLOGY._mesh_add_faces(
        nothing,random_mesh,4,face_identifiers,quadrangle_input)
    ascending_faces=UInt64[]
    foreach(face->append!(ascending_faces,face),quadrangle_keys)
    @test _ENTITY_TOPOLOGY._mesh_faces(
        random_faces,random_mesh,4,ascending_faces)==
        (face_identifiers,zeros(Int32,length(face_identifiers)))
    @test _ENTITY_TOPOLOGY._mesh_all_faces(random_faces,4)[1]==
          sort(face_identifiers)

    small_manual_edges,small_manual_faces=
        _manual_entity_topology_allocations(2_000)
    large_manual_edges,large_manual_faces=
        _manual_entity_topology_allocations(4_000)
    @test small_manual_edges>0
    @test small_manual_faces>0
    @test large_manual_edges<=2small_manual_edges+262_144
    @test large_manual_faces<=2small_manual_faces+262_144

    @test _manual_topology_digest(fixture)==
          "00110d2815d99ced60d72a8344958d0f0797fc5b85c938142dc4061f0abc8b06"
end
