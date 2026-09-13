using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, ntets, validate

const _TOPOLOGY_API=Tessella.API

function _add_api_topology_tetrahedron()
    for (tag,x,y,z) in ((10,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (7,0.0,1.0,0.0),(5,0.0,0.0,1.0))
        _TOPOLOGY_API.model.add_point(x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5))
        _TOPOLOGY_API.model.add_line(first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),
                         (23,[3,4,-12]),(24,[11,6,-4]))
        _TOPOLOGY_API.model.add_curve_loop(curves;tag=tag)
        _TOPOLOGY_API.model.add_plane_surface([tag];tag=tag)
    end
    _TOPOLOGY_API.model.add_surface_loop([21,-22,23,-24];tag=30)
    _TOPOLOGY_API.model.add_volume([30];tag=40)
    return nothing
end

@testset "cache-preserving explicit topology queries through API" begin
    _TOPOLOGY_API.finalize()
    @test_throws ArgumentError _TOPOLOGY_API.model.get_entities()
    @test_throws ArgumentError _TOPOLOGY_API.model.get_dimension()
    @test_throws ArgumentError _TOPOLOGY_API.model.get_boundary([])
    @test_throws ArgumentError _TOPOLOGY_API.model.get_adjacencies(0,1)

    try
        _TOPOLOGY_API.initialize()
        @test _TOPOLOGY_API.model.get_entities()==Tuple{Int,Int}[]
        @test _TOPOLOGY_API.model.get_dimension()==-1
        @test _TOPOLOGY_API.model.get_boundary([])==Tuple{Int,Int}[]
        _add_api_topology_tetrahedron()

        expected_entities=[
            (0,2),(0,5),(0,7),(0,10),
            (1,3),(1,4),(1,6),(1,8),(1,11),(1,12),
            (2,21),(2,22),(2,23),(2,24),(3,40)]
        @test _TOPOLOGY_API.model.get_entities()==expected_entities
        @test _TOPOLOGY_API.model.get_entities(2)==
              [(2,21),(2,22),(2,23),(2,24)]
        @test _TOPOLOGY_API.model.get_dimension()==3
        @test _TOPOLOGY_API.model.get_boundary(
            [(2,21),(2,22)],false,true,false)==
            [(1,8),(1,3),(1,11),(1,8),(1,12),(1,-6)]
        @test _TOPOLOGY_API.model.get_boundary(
            [(2,21),(2,22)],true,true,false)==
            [(1,3),(1,-6),(1,11),(1,12)]
        @test _TOPOLOGY_API.model.get_boundary(
            [(2,21),(2,22)],true,false,true)==[(0,5),(0,7)]
        @test _TOPOLOGY_API.model.get_boundary(
            [(3,40)],false,true,false)==
            [(2,21),(2,-22),(2,23),(2,-24)]
        @test _TOPOLOGY_API.model.get_adjacencies(0,10)==([6,8,11],Int[])
        @test _TOPOLOGY_API.model.get_adjacencies(1,8)==([21,22],[10,2])
        @test _TOPOLOGY_API.model.get_adjacencies(2,21)==([40],[8,3,11])
        @test _TOPOLOGY_API.model.get_adjacencies(3,40)==
              (Int[],[21,22,23,24])

        generated=_TOPOLOGY_API.mesh.generate(3)
        @test validate(generated).ok
        @test nnodes(generated)==9
        @test ntets(generated)==12
        generated_crc=mesh_crc(generated)
        @test generated_crc.sha==
              "979b12cba32c4e7e8317040d31636eca298b177460adec514374372613ba2f23"
        cached=_TOPOLOGY_API.LAST_MESH[]

        detached_entities=_TOPOLOGY_API.model.get_entities()
        detached_boundary=_TOPOLOGY_API.model.get_boundary(
            [(2,21),(2,22)],true,true,false)
        detached_upward,detached_downward=
            _TOPOLOGY_API.model.get_adjacencies(1,8)
        push!(detached_entities,(3,99))
        push!(detached_boundary,(1,99))
        push!(detached_upward,99)
        push!(detached_downward,99)
        @test _TOPOLOGY_API.model.get_entities()==expected_entities
        @test _TOPOLOGY_API.model.get_boundary(
            [(2,21),(2,22)],true,true,false)==
            [(1,3),(1,-6),(1,11),(1,12)]
        @test _TOPOLOGY_API.model.get_adjacencies(1,8)==([21,22],[10,2])
        @test _TOPOLOGY_API.LAST_MESH[]===cached
        @test mesh_crc(_TOPOLOGY_API.mesh.get())==generated_crc

        @test_throws ArgumentError _TOPOLOGY_API.model.get_entities(4)
        @test_throws ArgumentError _TOPOLOGY_API.model.get_boundary([(2,99)])
        @test_throws ArgumentError _TOPOLOGY_API.model.get_boundary(
            [(2,21)],true,1,false)
        @test_throws ArgumentError _TOPOLOGY_API.model.get_adjacencies(1,99)
        @test _TOPOLOGY_API.LAST_MESH[]===cached
        @test mesh_crc(_TOPOLOGY_API.mesh.get())==generated_crc

        names=(:get_entities,:get_dimension,:get_boundary,:get_adjacencies)
        metadata=Docs.meta(Tessella.API.model)
        @test all(name->haskey(
            metadata,Docs.Binding(Tessella.API.model,name)),names)
        @test isempty(Test.detect_ambiguities(Tessella.API;recursive=true))
    finally
        _TOPOLOGY_API.finalize()
    end

    @test_throws ArgumentError _TOPOLOGY_API.model.get_entities()
    @test_throws ArgumentError _TOPOLOGY_API.model.get_dimension()
    @test_throws ArgumentError _TOPOLOGY_API.model.get_boundary([])
    @test_throws ArgumentError _TOPOLOGY_API.model.get_adjacencies(0,1)
end

@testset "entity orphan queries through API" begin
    _TOPOLOGY_API.finalize()
    @test_throws ArgumentError _TOPOLOGY_API.model.is_entity_orphan(0,1)
    try
        _TOPOLOGY_API.initialize()
        _add_api_topology_tetrahedron()
        # Volume 40's transitive downward boundary closure is connected.
        for pair in ((3,40),(2,21),(2,24),(1,8),(1,4),(0,10),(0,5))
            @test _TOPOLOGY_API.model.is_entity_orphan(pair...)==false
        end
        # A floating point, a vertex-sharing curve, and an embedded point all
        # stay orphan: connectivity is downward-only and embeddings do not
        # connect, matching Gmsh 4.15.2's `isEntityOrphan`.
        _TOPOLOGY_API.model.add_point(9.0,9.0,9.0;tag=50)
        _TOPOLOGY_API.model.add_point(0.25,0.25,0.0;tag=63)
        _TOPOLOGY_API.model.add_line(10,63;tag=64)
        _TOPOLOGY_API.model.embed(0,[63],2,21)
        for pair in ((0,50),(0,63),(1,64))
            @test _TOPOLOGY_API.model.is_entity_orphan(pair...)==true
        end
        for call in (()->_TOPOLOGY_API.model.is_entity_orphan(0,99),
                     ()->_TOPOLOGY_API.model.is_entity_orphan(2,7),
                     ()->_TOPOLOGY_API.model.is_entity_orphan(0,-1),
                     ()->_TOPOLOGY_API.model.is_entity_orphan(0,0),
                     ()->_TOPOLOGY_API.model.is_entity_orphan(4,1),
                     ()->_TOPOLOGY_API.model.is_entity_orphan(0,1.5),
                     ()->_TOPOLOGY_API.model.is_entity_orphan(0,true))
            @test_throws ArgumentError call()
        end
        # Results are stable across mesh generation and its classification.
        generated=_TOPOLOGY_API.mesh.generate(3)
        @test validate(generated).ok
        for pair in ((3,40),(2,21),(0,10))
            @test _TOPOLOGY_API.model.is_entity_orphan(pair...)==false
        end
        for pair in ((0,50),(0,63),(1,64))
            @test _TOPOLOGY_API.model.is_entity_orphan(pair...)==true
        end
        metadata=Docs.meta(Tessella.API.model)
        @test haskey(
            metadata,Docs.Binding(Tessella.API.model,:is_entity_orphan))
        @test isempty(Test.detect_ambiguities(Tessella.API;recursive=true))
    finally
        _TOPOLOGY_API.finalize()
    end
    @test_throws ArgumentError _TOPOLOGY_API.model.is_entity_orphan(0,1)
end
