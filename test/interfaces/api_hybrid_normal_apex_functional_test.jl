module HybridNormalApexFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements

function cube(;inward=true,embedded=true)
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
             (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,point) in enumerate(corners)
        A.model.add_point(point...;tag,meshSize=1.0)
    end
    embedded && A.model.add_point(0.25,0.3,0.4;tag=9,meshSize=1.0)
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(edges)
        A.model.add_line(ends...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,loop) in enumerate(loops)
        lines=inward && tag>1 ? [-line for line in reverse(loop)] : collect(loop)
        A.model.add_curve_loop(lines;tag)
        A.model.add_plane_surface([tag];tag)
        A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    embedded && A.mesh.embed(0,[9],3,1)
end

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())

@testset "OptimizePyramids retains its finite double value" begin
    A.initialize()
    try
        @test A.option("Mesh.OptimizePyramids")==0.0
        for factor in (-0.1,-1.7,0.7)
            A.option("Mesh.OptimizePyramids",factor)
            @test A.option("Mesh.OptimizePyramids")==factor
        end
        @test_throws ArgumentError A.option("Mesh.OptimizePyramids",NaN)
        @test A.option("Mesh.OptimizePyramids")==0.7
    finally
        A.finalize()
    end
end

@testset "negative factor constructs fixed interior apices with an embedded Point" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube()
        A.mesh.generate(2)
        point_tags,point_xyz,_=A.mesh.get_nodes(0,9)
        point=only(point_tags)
        lower_nodes=A.mesh.get_nodes()[1:2]
        lower=Dict((dim,tag)=>A.mesh.get_elements(dim,tag)
                   for (dim,tag) in A.model.get_entities() if dim<3)
        node_clock=A.mesh.get_max_node_tag()
        A.option("Mesh.MeshOnlyEmpty",1)
        A.option("Mesh.OptimizePyramids",-0.1)
        A.mesh.generate(3)
        @test A.option("Mesh.OptimizePyramids")==-0.1
        @test all(A.mesh.get_elements(owner...)==cells for (owner,cells) in lower)
        @test A.mesh.get_nodes(0,9)[1:2]==(point_tags,point_xyz)
        families,labels,vertices=A.mesh.get_elements(3,1)
        @test Set(families)==Set(Int32[4,7])
        cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
                   for (kind,ids,nodes) in zip(families,labels,vertices))
        @test size(cells[7],2)==6
        @test point in cells[4] && !(point in cells[7])
        ids,xyz,_=A.mesh.get_nodes()
        positions=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        @test all(positions[uid]==lower_nodes[2][3row-2:3row]
                  for (row,uid) in enumerate(lower_nodes[1]))
        apices=vec(cells[7][5,:])
        @test Set(apices)==Set(node_clock+UInt64(index) for index in 1:6)
        @test all(uid in first(A.mesh.get_nodes(3,1)) for uid in apices)
        offset=0.1sqrt(2.0)
        expected=([0.5,0.5,offset],[0.5,0.5,1-offset],
                  [0.5,offset,0.5],[1-offset,0.5,0.5],
                  [0.5,1-offset,0.5],[offset,0.5,0.5])
        @test all(isapprox(positions[node_clock+UInt64(index)],expected[index];atol=8eps(),rtol=0)
                  for index in 1:6)
        function tet_volume(nodes)
            a,b,c,d=(positions[uid] for uid in nodes)
            det(hcat(b-a,c-a,d-a))/6
        end
        tet_volumes=[tet_volume(cell) for cell in eachcol(cells[4])]
        pyramid_volumes=[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]])
                         for cell in eachcol(cells[7])]
        @test all(>(0),tet_volumes) && all(>(0),pyramid_volumes)
        @test sum(tet_volumes)+sum(pyramid_volumes)≈1.0 atol=2e-10
        @test all(isapprox(volume,offset/3;atol=8eps(),rtol=0) for volume in pyramid_volumes)
        incidence=Dict{Tuple,Int}()
        for (kind,nodes) in cells
            faces=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
                ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
            for cell in eachcol(nodes),face in faces
                key=Tuple(sort(cell[collect(face)]))
                incidence[key]=get(incidence,key,0)+1
            end
        end
        quads=Set(Tuple(sort(only(last(A.mesh.get_elements(2,surface))))) for surface in 1:6)
        @test Set(face for (face,count) in incidence if count==1)==quads
        @test all(count->count in (1,2),values(incidence))
        @test all(!(point in face) for (face,count) in incidence if count==1)
        @test Tessella.validate(A.mesh.get()).ok
        mktempdir() do directory
            path=joinpath(directory,"normal-apex.msh")
            view=A.view.add("normal-apex Point")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,
                                  [[positions[uid][1]] for uid in ids])
            A.view.write(view,path)
            recovered=E.read_mixed_msh(path);data=recovered.entity_data
            @test Set(data.external_node_tags)==Set(ids)
            @test Set(uid for tags in data.external_element_tags for uid in tags)==
                  Set(uid for tags in A.mesh.get_elements()[2] for uid in tags)
            row=findfirst(==(point),data.external_node_tags)
            @test row!==nothing && data.node_entities[row]==(0,Int32(9)) &&
                  recovered.coords[:,row]==point_xyz
            @test all(data.node_entities[findfirst(==(uid),data.external_node_tags)]==(3,Int32(1))
                      for uid in apices)
            @test Tessella.validate(recovered).ok
        end
        println("HYBRID_NORMAL_APEX factor=-0.1 Tet4=",size(cells[4],2),
                " Pyramid5=6 PointUID=",point," volume=",sum(tet_volumes)+sum(pyramid_volumes))
        flush(stdout)
    finally
        A.finalize()
    end
end

@testset "exterior normal births reject atomically" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube(;embedded=false)
        A.mesh.generate(2)
        A.option("Mesh.MeshOnlyEmpty",1)
        before=snapshot()
        for factor in (-1.0,)
            A.option("Mesh.OptimizePyramids",factor)
            @test_throws ArgumentError A.mesh.generate(3)
            @test snapshot()==before
        end
    finally
        A.finalize()
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube(;inward=false,embedded=false)
        A.mesh.generate(2)
        A.option("Mesh.MeshOnlyEmpty",1)
        A.option("Mesh.OptimizePyramids",-0.1)
        before=snapshot()
        @test_throws ArgumentError A.mesh.generate(3)
        @test snapshot()==before
    finally
        A.finalize()
    end
end
end
