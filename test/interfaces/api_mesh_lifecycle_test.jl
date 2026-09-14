using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, validate, nnodes, ntets

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
              "eae8751b0dad3b89f2d7a4416ea079a352a3bd6b8eff31ef7eb7d8bb78d8509a"

        refined=_MESH_LIFECYCLE_API.mesh.refine()
        refined_crc=mesh_crc(refined)
        @test validate(refined).ok
        @test refined_crc.n_nodes==35
        @test refined_crc.n_tets==96
        @test refined_crc.sha==
              "83415c157f9b4daf2124e34562383d36be6fde037af214b0c34d14156dc93b15"

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
              "148948f0e3430d600b6b8e46214b047bf2074fe581fdd61b59c7eb962c1997e7"

        # `mesh.clear` removes only elements classified on the Volume; the
        # materialized corner Points keep their nodes, matching Gmsh's
        # `removeEntities`-style per-entity clearing. Unknown entities fail.
        @test _MESH_LIFECYCLE_API.mesh.clear([(3,1)])===nothing
        leftover=_MESH_LIFECYCLE_API.mesh.get()
        @test nnodes(leftover)==98 && ntets(leftover)==0
        _MESH_LIFECYCLE_API.mesh.generate(3)
        _MESH_LIFECYCLE_API.mesh.refine()
        _MESH_LIFECYCLE_API.mesh.refine()
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==twice_crc
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.clear([(3,99)])
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==twice_crc
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.clear(1)
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==twice_crc

        @test _MESH_LIFECYCLE_API.mesh.clear(())===nothing
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.get()
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine()
        @test _MESH_LIFECYCLE_API.model.get_entities()==vcat(
            [(0,i) for i in 1:8],[(1,i) for i in 1:12],
            [(2,i) for i in 1:6],[(3,1)])

        regenerated=_MESH_LIFECYCLE_API.mesh.generate(3)
        @test mesh_crc(regenerated)==generated_crc
        @test _MESH_LIFECYCLE_API.mesh.clear(Int32[])===nothing
        @test _MESH_LIFECYCLE_API.model.get_entities()==vcat(
            [(0,i) for i in 1:8],[(1,i) for i in 1:12],
            [(2,i) for i in 1:6],[(3,1)])
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
