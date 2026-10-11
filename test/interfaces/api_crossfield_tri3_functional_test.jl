module CrossfieldTri3FunctionalChecks
using Test, Tessella
const A=Tessella.API

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
cache_snapshot()=let mesh=A.mesh.get()
    (copy(mesh.coords),[(block.msh,copy(block.nodes),copy(block.tags)) for block in mesh.blocks],
     deepcopy(mesh.entity_data))
end

function triangle_cache()
    A.option("General.Terminal",0)
    A.option("Mesh.Renumber",0)
    A.model.add_discrete_entity(1,4)
    A.model.add_discrete_entity(2,7)
    A.model.add_discrete_entity(3,9)
    A.mesh.add_nodes(3,9,[101,102,103,104],
                     [0.,0.,0.,2.,0.,0.,1.,1.,0.,0.,0.,1.])
    A.mesh.add_elements_by_type(4,1,[301],[101,102])
    A.mesh.add_elements_by_type(7,2,[401],[101,102,103])
    A.mesh.add_elements_by_type(9,4,[501],[101,102,103,104])
    # Publish actual classified data through a real public operation.
    A.mesh.generate(0)
end

function field_data(views)
    return [only(A.view.get_list_data(view)[2]) for view in views]
end

@testset "Tri3 cross fields use actual cache and caller identities" begin
    A.initialize()
    try
        triangle_cache()
        @test A.mesh.get() isa Tessella.Elements.MixedMesh
        @test A.mesh.get_element_types(2)==Int32[2]
        before=snapshot()
        cached=cache_snapshot()
        views=A.mesh.compute_cross_field()
        @test length(views)==3
        @test first.(A.view.get_list_data.(views))==[["ST"],["ST"],["VT"]]
        data=field_data(views)
        @test length.(data)==[12,12,18]
        @test data[1][1:9]==[0.,2.,1.,0.,0.,1.,0.,0.,0.]
        @test data[1][10:12]≈fill((2+2sqrt(2))/3,3)
        @test data[3][10:18]≈repeat([1.,0.,0.],3)
        @test snapshot()==before
        after=cache_snapshot()
        @test after[1:2]==cached[1:2]
        @test after[3].external_node_tags==cached[3].external_node_tags &&
              after[3].external_element_tags==cached[3].external_element_tags

        # Add a caller cell on the same authority-covered Surface, with
        # coincident coordinates but distinct UIDs and a different first edge.
        # It must be included once and remain disconnected from the cached cell.
        A.mesh.add_nodes(2,7,[201,202,203],[0.,0.,0.,2.,0.,0.,1.,1.,0.])
        A.mesh.add_elements_by_type(7,2,[402],[203,201,202])
        @test A.mesh.get_elements_by_type(2)==
            (UInt64[401,402],UInt64[101,102,103,203,201,202])
        before=snapshot()
        cached=cache_snapshot()
        views=A.mesh.compute_cross_field()
        data=field_data(views)
        @test length.(data)==[24,24,36]
        @test data[3][10:18]≈repeat([1.,0.,0.],3)
        @test data[3][28:36]≈repeat([-inv(sqrt(2)),-inv(sqrt(2)),0.],3)
        @test maximum(abs,data[2][vcat(10:12,22:24)])<=3e-8
        @test snapshot()==before
        @test cache_snapshot()[1:2]==cached[1:2]
        saved=copy(data[3])
        data[3].=99
        @test only(A.view.get_list_data(views[3])[2])==saved

        # A different caller triangle really shares the cache's edge UIDs.
        # It contributes to that component's smoothing while the coincident
        # UID-disjoint triangle retains its independent first-edge cross.
        A.mesh.add_nodes(2,7,[205],[1.,-1.,0.])
        A.mesh.add_elements_by_type(7,2,[403],[205,101,102])
        before=snapshot()
        views=A.mesh.compute_cross_field()
        data=field_data(views)
        @test length.(data)==[36,36,54]
        @test sum(abs2,data[3][10:12]-[1.,0.,0.])>1e-6
        @test data[3][28:36]≈repeat([-inv(sqrt(2)),-inv(sqrt(2)),0.],3)
        @test all(offset->isapprox(sum(abs2,data[3][offset:offset+2]),1.0;atol=1e-12),
                  (10,28,46))
        @test snapshot()==before
    finally
        A.finalize()
    end
end

@testset "unsupported Surface families reject before publishing views" begin
    for kind in (3,9)
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.model.add_discrete_entity(2,7)
            nodes=kind==3 ? [101,102,103,104] : [101,102,103,104,105,106]
            xyz=kind==3 ? [0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.] :
                [0.,0.,0.,1.,0.,0.,0.,1.,0.,.5,0.,0.,.5,.5,0.,0.,.5,0.]
            A.mesh.add_nodes(2,7,nodes,xyz)
            A.mesh.add_elements_by_type(7,kind,[401],nodes)
            before=snapshot()
            views=A.view.get_tags()
            @test_throws r"linear Tri3 Surface" A.mesh.compute_cross_field()
            @test A.view.get_tags()==views
            @test snapshot()==before
        finally
            A.finalize()
        end
    end
end
end
