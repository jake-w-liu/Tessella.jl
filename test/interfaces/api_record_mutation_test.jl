using Test
using Tessella

@testset "atomic sparse mesh record mutations" begin
    api=Tessella.API
    function records(;topology=false)
        api.initialize()
        api.model.add_discrete_entity(1,1)
        api.mesh.add_nodes(1,1,[10,20,30],[0.,0,0,1,0,0,2,0,0],[0.,0.5,1.])
        api.mesh.add_elements_by_type(1,1,[40,50],[10,20,20,30])
        topology && api.mesh.create_topology()
    end
    snapshot()=(api.mesh.get_nodes(),api.mesh.get_elements())
    try
        for topology in (false,true)
            records(;topology)
            before=snapshot()
            api.mesh.renumber_nodes([10,20,30],[20,30,10])
            @test api.mesh.get_node(20)[1]==[0.,0,0]
            @test api.mesh.get_node(30)[1]==[1.,0,0]
            @test api.mesh.get_node(10)[1]==[2.,0,0]
            @test api.mesh.get_elements_by_type(1)[2]==UInt64[20,30,30,10]
            # The inverse must restore every entity, including endpoint records.
            api.mesh.renumber_nodes([20,30,10],[10,20,30])
            @test snapshot()==before
            api.mesh.renumber_elements([40,50],[50,40])
            @test api.mesh.get_element(50)[2]==UInt64[10,20]
            @test api.mesh.get_element(40)[2]==UInt64[20,30]
            api.mesh.renumber_elements([50,40],[40,50])
            @test snapshot()==before
        end
        for operation in (:renumber_nodes,:renumber_elements)
            tags=operation==:renumber_nodes ? [10,20] : [40,50]
            for (old,new) in ((tags,[70]),(tags,[70,70]),
                              ([tags[1],tags[1]],[70,80]),(Int[],[70]),
                              ([tags[1]],[0]),([tags[1]],[big(typemax(Int32))+1]),
                              ([tags[1]],[tags[2]]))
                records()
                before=snapshot()
                @test_throws ArgumentError getproperty(api.mesh,operation)(old,new)
                @test snapshot()==before
            end
        end
        records()
        before=api.mesh.get_node(10)
        for params in ([0.,1.],[NaN],[Inf])
            @test_throws ArgumentError api.mesh.set_node(10,[9.,9,9],params)
            @test api.mesh.get_node(10)==before
        end
        api.mesh.set_node(10,[0.25,0,0],[0.125])
        @test api.mesh.get_node(10)==([0.25,0,0],[0.125],1,1)
        api.mesh.set_node(10,[0.5,0,0])
        @test api.mesh.get_node(10)==([0.5,0,0],[0.125],1,1)

        # Dense-cache validation must finish before a sparse map is committed.
        records()
        mesh=Mesh([0. 1 0;0 0 1;0 0 0];
                  segs=reshape(Int32[1,2],2,1),
                  tris=reshape(Int32[1,2,3],3,1))
        lock(api.STATE_LOCK) do
            api._replace_mesh_cache_locked!(mesh)
        end
        before=snapshot()
        @test_throws ArgumentError api.mesh.renumber_nodes([10,1],[70,99])
        @test snapshot()==before
        @test_throws ArgumentError api.mesh.renumber_elements([40,1,2],[70,2,1])
        @test snapshot()==before
        api.mesh.renumber_nodes([10,1,2],[70,2,1])
        @test api.mesh.get_node(70)[1]==[0.,0,0]
        @test api.mesh.get_elements_by_type(1)[2]==UInt64[2,1,70,20,20,30]
        @test api.mesh.get().segs==reshape(Int32[2,1],2,1)
        api.mesh.renumber_elements([40,1],[70,1])
        @test api.mesh.get_element(70)[2]==UInt64[70,20]
    finally
        api.finalize()
    end
end

@testset "global TransfiniteTri option lifecycle" begin
    api=Tessella.API
    function triangle()
        for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(0.,1.)))
            api.model.add_point(x,y,0;tag)
        end
        for (tag,(a,b)) in enumerate(((1,2),(2,3),(3,1)))
            api.model.add_line(a,b;tag)
            api.mesh.set_transfinite_curve(tag,5)
        end
        api.model.add_curve_loop([1,2,3];tag=1)
        api.model.add_plane_surface([1];tag=1)
        api.mesh.set_transfinite_surface(1)
    end
    try
        api.initialize()
        api.option("Mesh.TransfiniteTri",1)
        api.initialize()
        @test api.option("Mesh.TransfiniteTri")==0
        triangle()
        @test size(api.mesh.generate(2).coords,2)==21
        api.model.add("other")
        api.option("Mesh.TransfiniteTri",1)
        api.model.set_current("")
        @test api.option("Mesh.TransfiniteTri")==1
        @test size(api.mesh.generate(2).coords,2)==15
    finally
        api.finalize()
    end
end
