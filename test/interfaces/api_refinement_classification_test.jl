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
                  segs=reshape(Int32[2,3],2,1),
                  tris=reshape(Int32[2,3,4],3,1))
        owners=Tuple{Int,Int32}[(0,99),(0,1),(1,2),(2,3)]
        # Curve2 owns an actual line support; endpoint labels alone cannot
        # establish the carrier of a quadratic midpoint.
        class=classified(mesh,owners;segs=Int32[2],tris=Int32[3])
        refined=Tessella.Refine.refine_uniform(mesh)
        result=api._inherit_refined_classification(class,refined,refined)
        @test result.node_entities==Tuple{Int,Int32}[
            (0,1),(1,2),(2,3),(1,2),(2,3),(2,3)]
        @test result.tri_entities==fill(Int32(3),4)
        @test result.seg_entities==fill(Int32(2),2)
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

@testset "single curve refinement retains classification and parameters" begin
    api=Tessella.API
    try
        api.initialize()
        api.model.add_point(0,0,0;tag=1)
        api.model.add_point(2,0,0;tag=2)
        api.model.add_line(1,2;tag=1)
        api.mesh.set_transfinite_curve(1,3)
        api.mesh.generate(1)
        for level in 1:2
            api.mesh.refine()
            segments=2^(level+1)
            @test api.mesh.get_element_types(1,1)==Int32[1]
            @test length(api.mesh.get_elements_by_type(1,1)[1])==segments
            @test length(api.mesh.get_nodes(1,1,true)[1])==segments+1
            tags,_,parameters=api.mesh.get_nodes(1,1)
            @test length(tags)==segments-1
            @test sort(parameters)==[index/segments for index in 1:segments-1]
            @test length(api.mesh.get_nodes(0,1)[1])==1
            @test length(api.mesh.get_nodes(0,2)[1])==1
        end
    finally
        api.finalize()
    end
end

@testset "single and multiple surface refinement retains classification" begin
    api=Tessella.API
    for surface_count in (1,2)
    try
        api.initialize()
        for surface in 1:surface_count
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
        counts=[length(api.mesh.get_elements_by_type(2,tag)[1]) for tag in 1:surface_count]
        for level in 1:2
            @test validate(api.mesh.refine()).ok
            for tag in 1:surface_count
                @test length(api.mesh.get_elements_by_type(2,tag)[1])==4^level*counts[tag]
                @test !isempty(api.mesh.get_nodes(2,tag)[1])
            end
        end
    finally
        api.finalize()
    end
    end
end

@testset "reversed volume refinement keeps entity classification" begin
    api=Tessella.API
    orientation(mesh)=Int[Tessella.Predicates.orient3(
        (Tessella.MeshTypes.node(mesh,mesh.tets[j,i]) for j in 1:4)...)
        for i in axes(mesh.tets,2)]
    try
        for reverse_all in (true,false)
            api.initialize()
            api.model.add_box(0,0,0,1,1,1)
            api.mesh.generate(3)
            reverse_all ? api.mesh.reverse() : api.mesh.reverse_elements([1])
            original=api.mesh.get()
            signs=orientation(original)
            projected=api._classification_skeleton(original,true)
            @test validate(projected).ok
            @test orientation(original)==signs
            @test projected.coords==original.coords
            @test projected.tet_tag==original.tet_tag
            refined=api.mesh.refine()
            @test orientation(refined)==repeat(signs;inner=8)
            @test api.mesh.get_element_types(3,1)==Int32[4]
            @test length(api.mesh.get_elements_by_type(4,1)[1])==size(refined.tets,2)
            @test !isempty(api.mesh.get_nodes(0,1)[1])
            @test validate(refined;require_positive_tets=false).ok
        end
    finally
        api.finalize()
    end
end

@testset "boundary size extension option accepts integer flags" begin
    api=Tessella.API
    try
        api.initialize()
        for (input,expected) in ((0.,0.),(.75,0.),(-1.9,-1.),(2.9,2.),(1.,1.))
            @test api.option("Mesh.MeshSizeExtendFromBoundary",input)==expected
            @test api.option("Mesh.MeshSizeExtendFromBoundary")==expected
        end
        before=api.option("Mesh.MeshSizeExtendFromBoundary")
        for value in (true,NaN,Inf,Float64(typemax(Int32))+1)
            @test_throws ArgumentError api.option("Mesh.MeshSizeExtendFromBoundary",value)
            @test api.option("Mesh.MeshSizeExtendFromBoundary")==before
        end
        api.option("Mesh.MeshSizeExtendFromBoundary",0)
        api.model.add_box(0,0,0,1,1,1)
        @test validate(api.mesh.generate(3)).ok
        @test api.CURRENT[].meshing.lc_extend_from_boundary==0
    finally
        api.finalize()
    end
end

@testset "boundary extension resolution honors 0 and negative modes" begin
    api=Tessella.API
    resolve=Tessella.Model._resolved_extend_from_boundary
    function square_model()
        api.initialize()
        api.model.add_point(0,0,0)
        api.model.add_point(1,0,0)
        api.model.add_point(1,1,0)
        api.model.add_point(0,1,0)
        api.model.add_line(1,2)
        api.model.add_line(2,3)
        api.model.add_line(3,4)
        api.model.add_line(4,1)
        api.model.add_curve_loop(collect(1:4))
        api.model.add_plane_surface([1])
        return api.CURRENT[]
    end
    # A stored per-entity flag resolves verbatim: 0 disables even while the
    # global default enables; a negative record defers to the global like
    # upstream's unset -1.
    try
        m=square_model()
        api.mesh.set_size_from_boundary(2,1,0)
        @test m.meshing.size_from_boundary[(2,1)]==0
        @test resolve(m,2,1)==0
        api.mesh.set_size_from_boundary(2,1,2)
        @test resolve(m,2,1)==2
        api.mesh.set_size_from_boundary(2,1,-1)
        @test resolve(m,2,1)==1
        @test_throws ArgumentError api.mesh.set_size_from_boundary(2,1,0.5)
        @test_throws ArgumentError api.mesh.set_size_from_boundary(3,1,1)
    finally
        api.finalize()
    end
    # With no per-entity record the surface falls back to the global — an
    # explicit per-entity 0 still overrides an enabled global. OPTIONS only
    # mirrors into `m.meshing.lc_extend_from_boundary` at generate time, so
    # the model field is set directly here (the generate-level sync is
    # covered by the mesh check below).
    for mode in (0,-2,-3)
        try
            m=square_model()
            m.meshing.lc_extend_from_boundary=mode
            @test resolve(m,2,1)==mode
            api.mesh.set_size_from_boundary(2,1,0)
            @test resolve(m,2,1)==0
            delete!(m.meshing.size_from_boundary,(2,1))
            @test resolve(m,2,1)==mode
        finally
            api.finalize()
        end
    end
    # Disabled extension still meshes — the global 0 option mirrors into the
    # model at generate time and reaches the surface path.
    try
        m=square_model()
        api.option("Mesh.MeshSizeExtendFromBoundary",0)
        api.mesh.set_size_from_boundary(2,1,0)
        @test validate(api.mesh.generate(2)).ok
        @test api.CURRENT[].meshing.lc_extend_from_boundary==0
    finally
        api.finalize()
    end
end
