module PeriodicLifecycleFunctionalChecks
using Test,Tessella
using Tessella.MeshTypes: validate
@testset "periodic high-order and refinement lifecycle" begin
    api = Tessella.API
    api.initialize()
    try
        api.option("Mesh.Renumber", 0)
        for (face, first, offset) in ((1, 1, 0.0), (2, 11, 2.0))
            for (i, p) in enumerate(((offset,0.,0.), (offset+1.,0.,0.),
                                     (offset+1.,1.,0.), (offset,1.,0.)))
                api.model.add_point(p...; tag=first+i-1)
            end
            for i in 0:3
                api.model.add_line(first+i, first+mod(i+1,4); tag=first+i)
            end
            api.model.add_curve_loop(collect(first:first+3); tag=face)
            api.model.add_plane_surface([face]; tag=face)
            tags = UInt64.(face*100 .+ (1:3))
            api.mesh.add_nodes(2, face, tags,
                [offset,0.,0., offset+1.,0.,0., offset,1.,0.])
            api.mesh.add_elements_by_type(face, 2, [face*300+1], tags)
        end
        affine = [1.,0.,0.,2., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
        api.mesh.set_periodic(2, [2], [1], affine)
        function correspondence(higher)
            relation = api.mesh.get_periodic_nodes(2, 2, higher)
            @test relation.master_entity == 1
            @test length(relation.slave_nodes) == length(relation.master_nodes)
            @test allunique(relation.slave_nodes) && allunique(relation.master_nodes)
            for (slave, master) in zip(relation.slave_nodes, relation.master_nodes)
                s = first(api.mesh.get_node(slave)); m = first(api.mesh.get_node(master))
                @test s ≈ m .+ [2.,0.,0.]
            end
            relation
        end
        @test length(correspondence(false).slave_nodes) == 3
        api.mesh.set_order(2)
        @test length(correspondence(false).slave_nodes) == 3
        @test length(correspondence(true).slave_nodes) == 6
        api.mesh.renumber_nodes([101], [0])
        @test 0 in correspondence(true).master_nodes
        @test first(api.mesh.get_node(0)) == [0.,0.,0.]
        api.mesh.refine()
        @test length(first(api.mesh.get_elements_by_type(2))) == 8
        @test length(correspondence(false).slave_nodes) == 6
        @test length(correspondence(true).slave_nodes) == 6
        @test validate(api.mesh.get()).ok
    finally
        api.finalize()
    end
end

end
