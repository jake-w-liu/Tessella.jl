@testset "quarter-turn periodic Surface copy" begin
    A=Tessella.API;A.initialize()
    try
        A.option("General.Terminal",0)
        points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                (0,0,1),(1,0,1),(1,1,1),(0,1,1))
        for (tag,p) in enumerate(points)
            A.model.add_point(p...;tag,meshSize=1.)
        end
        edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
               (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
        for (tag,ends) in enumerate(edges)
            A.model.add_line(ends...;tag)
            A.mesh.set_transfinite_curve(tag,2)
        end
        loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
               (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
        for (tag,lines) in enumerate(loops)
            A.model.add_curve_loop(collect(lines);tag)
            A.model.add_plane_surface([tag];tag)
            tag!=2 && A.mesh.set_transfinite_surface(tag)
            A.mesh.set_recombine(2,tag)
        end
        A.mesh.set_periodic(2,[2],[1],
            [0.,-1.,0.,1.,1.,0.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.])
        A.mesh.generate(2)
        kinds,tags,nodes=A.mesh.get_elements(2,2)
        @test kinds==Int32[3] && length(only(tags))==1
        owners=sort!([Tuple(sort(only(A.mesh.get_elements(2,surface)[3]))[1:4])
                      for surface in 1:6])
        @test length(unique(owners))==6
        link=A.mesh.get_periodic_nodes(2,2)
        @test link.master_entity==1 && length(link.slave_nodes)==4
        @test all(zip(link.slave_nodes,link.master_nodes)) do (slave,master)
            sm,pm=first(A.mesh.get_node(slave)),first(A.mesh.get_node(master))
            sm==[-pm[2]+1,pm[1],pm[3]+1]
        end
    finally;A.finalize();end
end
