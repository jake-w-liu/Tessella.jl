using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, validate

const _MESH_LIFECYCLE_API=Tessella.API

@testset "atomic cached-mesh refinement and clearing through API" begin
    _MESH_LIFECYCLE_API.finalize()
    @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine()
    @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.clear()

    try
        _MESH_LIFECYCLE_API.initialize()
        @test _MESH_LIFECYCLE_API.mesh.clear()===nothing
        @test _MESH_LIFECYCLE_API.mesh.clear([])===nothing
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine()

        @test _MESH_LIFECYCLE_API.model.add_box(
            0,0,0,1,1,1;tag=1)==1
        generated=_MESH_LIFECYCLE_API.mesh.generate(3)
        generated_crc=mesh_crc(generated)
        @test generated_crc.sha==
              "e9f6cd048ad689d1566e9c6664824543863983b8df79d9c0fa50f1f35d31cf83"

        refined=_MESH_LIFECYCLE_API.mesh.refine()
        refined_crc=mesh_crc(refined)
        @test validate(refined).ok
        @test refined_crc.n_nodes==35
        @test refined_crc.n_tets==96
        @test refined_crc.sha==
              "6fb8a362968e08263e38a6c59444f7b5503ffdb74b7db54cb90dfd59b683cc69"

        stored=_MESH_LIFECYCLE_API.LAST_MESH[]
        @test stored!==nothing && stored!==refined
        for field in (:coords,:segs,:tris,:tets,:seg_tag,:tri_tag,:tet_tag)
            @test !Base.mightalias(getfield(refined,field),getfield(stored,field))
        end
        refined.coords[1,1]+=17.0
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==refined_crc

        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine(
            max_nodes=188)
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==refined_crc
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine(
            max_cells=767)
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==refined_crc
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine(
            max_nodes=true)
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==refined_crc

        twice=_MESH_LIFECYCLE_API.mesh.refine(
            max_nodes=189,max_cells=768)
        twice_crc=mesh_crc(twice)
        @test validate(twice).ok
        @test twice_crc.n_nodes==189
        @test twice_crc.n_tets==768
        @test twice_crc.sha==
              "09fd5ced56aba7a5b1b0380f8f9189dc95d3676793430f61fd01269baca1445c"

        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.clear([(3,1)])
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==twice_crc
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.clear(1)
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==twice_crc

        @test _MESH_LIFECYCLE_API.mesh.clear(())===nothing
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.get()
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine()
        @test _MESH_LIFECYCLE_API.model.get_entities()==[(3,1)]

        regenerated=_MESH_LIFECYCLE_API.mesh.generate(3)
        @test mesh_crc(regenerated)==generated_crc
        @test _MESH_LIFECYCLE_API.mesh.clear(Int32[])===nothing
        @test _MESH_LIFECYCLE_API.model.get_entities()==[(3,1)]
    finally
        _MESH_LIFECYCLE_API.finalize()
    end

    @test isempty(Docs.undocumented_names(
        Tessella.API.mesh;private=false))
    @test isempty(Test.detect_ambiguities(
        Tessella.API.mesh;recursive=true))
end

@testset "periodic maps follow cached refinement" begin
    _MESH_LIFECYCLE_API.finalize()
    try
        _MESH_LIFECYCLE_API.initialize()
        for (tag,(x,y)) in enumerate(((0.0,0.0),(1.0,0.0),
                                      (1.0,1.0),(0.0,1.0)))
            @test _MESH_LIFECYCLE_API.model.add_point(
                x,y,0;tag=tag,meshSize=0.5)==tag
        end
        for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            @test _MESH_LIFECYCLE_API.model.add_line(
                first,last;tag=tag)==tag
        end
        @test _MESH_LIFECYCLE_API.model.add_curve_loop(
            [1,2,3,4];tag=1)==1
        @test _MESH_LIFECYCLE_API.model.add_plane_surface([1];tag=1)==1
        translate_x=(1.0,0.0,0.0,1.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0)
        @test _MESH_LIFECYCLE_API.mesh.set_periodic(
            1,[2],[4],translate_x)===nothing

        generated=_MESH_LIFECYCLE_API.mesh.generate(2)
        @test validate(generated).ok
        @test length(_MESH_LIFECYCLE_API.mesh.get_periodic_nodes(
            1,2).slave_nodes)==5

        refined=_MESH_LIFECYCLE_API.mesh.refine()
        @test validate(refined).ok
        mapping=_MESH_LIFECYCLE_API.mesh.get_periodic_nodes(1,2)
        @test length(mapping.slave_nodes)==length(mapping.master_nodes)==9
        cached=_MESH_LIFECYCLE_API.mesh.get()
        for (slave,master) in zip(mapping.slave_nodes,mapping.master_nodes)
            @test Tuple(cached.coords[:,slave])==
                  (cached.coords[1,master]+1,cached.coords[2,master],
                   cached.coords[3,master])
        end

        @test _MESH_LIFECYCLE_API.mesh.clear()===nothing
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.get_periodic_nodes(1,2)
    finally
        _MESH_LIFECYCLE_API.finalize()
    end
end
