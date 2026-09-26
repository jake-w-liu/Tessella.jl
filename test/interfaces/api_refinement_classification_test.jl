using Test
using Tessella

@testset "refinement classification uses node identity" begin
    api=Tessella.API
    function classified(mesh,owners;segs=Int32[],tris=Int32[],tets=Int32[])
        entity=(3,Int32(7))
        api._MeshClassification(mesh,entity,[entity],owners,
            Dict{Tuple{Int,Int32},Vector{Int32}}(),segs,tris,tets)
    end
    for scale in (1e-14,1.,1e8)
        # The unused first node coincides with a referenced point but has a
        # different identity and owner; it must disappear during compaction.
        mesh=Mesh(scale.*[0. 0 1 0;0 0 0 1;0 0 0 0];
                  tris=reshape(Int32[2,3,4],3,1))
        owners=Tuple{Int,Int32}[(0,99),(0,1),(1,2),(2,3)]
        class=classified(mesh,owners;tris=Int32[3])
        refined=Tessella.Refine.refine_uniform(mesh)
        result=api._inherit_refined_classification(class,refined,refined)
        @test result.node_entities==Tuple{Int,Int32}[
            (0,1),(1,2),(2,3),(1,2),(2,3),(2,3)]
        @test result.tri_entities==fill(Int32(3),4)
        @test result.mesh===refined
        @test class.node_entities==owners
    end
    mesh=Mesh([0. 1 0 1;0 0 0 0;0 0 0 0];segs=Int32[1 3;2 4])
    class=classified(mesh,Tuple{Int,Int32}[(1,1),(1,1),(1,2),(1,2)];
                     segs=Int32[1,2])
    refined=Tessella.Refine.refine_uniform(mesh)
    result=api._inherit_refined_classification(class,refined,refined)
    @test result.node_entities==Tuple{Int,Int32}[(1,1),(1,1),(1,2),(1,2),(1,1),(1,2)]
    @test result.seg_entities==Int32[1,1,2,2]

    mesh=Mesh([0. 1 0 0;0 0 1 0;0 0 0 1];
              segs=reshape(Int32[1,2],2,1),
              tris=reshape(Int32[1,2,3],3,1),
              tets=reshape(Int32[1,2,3,4],4,1))
    owners=Tuple{Int,Int32}[(0,i) for i in 1:4]
    class=classified(mesh,owners;segs=Int32[9],tris=Int32[8],tets=Int32[7])
    refined=Tessella.Refine.refine_uniform(mesh)
    result=api._inherit_refined_classification(class,refined,refined)
    @test result.node_entities[1:4]==owners
    # Lexicographic edges: 12, 13, 14, 23, 24, 34. Lower-dimensional
    # cells own their shared edge midpoints even when a tet also uses them.
    @test result.node_entities[5:end]==Tuple{Int,Int32}[(1,9),(2,8),(3,7),(2,8),(3,7),(3,7)]
    @test result.tet_entities==fill(Int32(7),8)

    empty_mesh=Mesh(zeros(3,0))
    class=classified(empty_mesh,Tuple{Int,Int32}[])
    refined=Tessella.Refine.refine_uniform(empty_mesh)
    @test isempty(api._inherit_refined_classification(class,refined,refined).node_entities)
    for higher in ((1,Int32(2)),(2,Int32(3)),(3,Int32(4)))
        @test api._merged_node_owner((0,Int32(1)),higher)==(0,Int32(1))
        @test api._merged_node_owner(higher,(0,Int32(1)))==(0,Int32(1))
        @test api._merged_node_owner((0,Int32(0)),higher)==higher
        @test api._merged_node_owner(higher,(0,Int32(0)))==higher
    end
end

@testset "multiple surface refinement retains classification" begin
    api=Tessella.API
    try
        api.initialize()
        for surface in 1:2
            offset=3(surface-1)
            for (i,(x,y)) in enumerate(((0.,0.),(1.,0.),(0.,1.)))
                api.model.add_point(x+2(surface-1),y,0;tag=offset+i)
            end
            for (i,(a,b)) in enumerate(((1,2),(2,3),(3,1)))
                api.model.add_line(offset+a,offset+b;tag=offset+i)
                api.mesh.set_transfinite_curve(offset+i,3)
            end
            api.model.add_curve_loop(collect(offset+1:offset+3);tag=surface)
            api.model.add_plane_surface([surface];tag=surface)
            api.mesh.set_transfinite_surface(surface)
        end
        api.mesh.generate(2)
        counts=[length(api.mesh.get_elements_by_type(2,tag)[1]) for tag in 1:2]
        for level in 1:2
            @test validate(api.mesh.refine()).ok
            for tag in 1:2
                @test length(api.mesh.get_elements_by_type(2,tag)[1])==4^level*counts[tag]
                @test !isempty(api.mesh.get_nodes(2,tag)[1])
            end
        end
    finally
        api.finalize()
    end
end
