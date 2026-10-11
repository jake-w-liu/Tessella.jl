module HybridEmbeddedPointFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements

function cube(point)
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
             (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(corners)
        A.model.add_point(p...;tag,meshSize=1.)
    end
    A.model.add_point(point...;tag=9,meshSize=1.)
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
        A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    A.mesh.embed(0,[9],3,1)
end

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())

@testset "actual embedded Point is a fixed conforming hybrid volume vertex" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube((0.25,0.3,0.4))
        A.mesh.generate(2)
        point_tags,point_xyz,_=A.mesh.get_nodes(0,9)
        point=only(point_tags)
        @test point_xyz==[0.25,0.3,0.4]
        lower=Dict((dim,tag)=>A.mesh.get_elements(dim,tag) for (dim,tag) in A.model.get_entities()
                   if dim<3)
        A.option("Mesh.MeshOnlyEmpty",1)
        A.mesh.generate(3)
        @test A.mesh.get_nodes(0,9)[1:2]==(point_tags,point_xyz)
        @test all(A.mesh.get_elements(owner...)==cells for (owner,cells) in lower)
        families,labels,vertices=A.mesh.get_elements(3,1)
        @test Set(families)==Set(Int32[4,7])
        cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
                   for (kind,ids,nodes) in zip(families,labels,vertices))
        @test point in cells[4] && !(point in cells[7])
        @test size(cells[7],2)==6
        ids,xyz,_=A.mesh.get_nodes()
        positions=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        function tet_volume(nodes)
            a,b,c,d=(positions[uid] for uid in nodes)
            det(hcat(b-a,c-a,d-a))/6
        end
        volumes=[tet_volume(cell) for cell in eachcol(cells[4])]
        append!(volumes,[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]])
                         for cell in eachcol(cells[7])])
        @test all(>(0),volumes)
        @test sum(volumes)≈1.0 atol=2e-10
        println("HYBRID_EMBEDDED_POINT Tet4=",size(cells[4],2),
                " Pyramid5=",size(cells[7],2)," PointUID=",point," volume=",sum(volumes));flush(stdout)
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
            path=joinpath(directory,"embedded-point.msh")
            view=A.view.add("embedded Point")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,
                                  [[positions[uid][1]] for uid in ids])
            A.view.write(view,path)
            recovered=E.read_mixed_msh(path)
            data=recovered.entity_data;row=findfirst(==(point),data.external_node_tags)
            @test row!==nothing && data.node_entities[row]==(0,Int32(9))
            @test recovered.coords[:,row]==point_xyz
            @test any(point in data.external_node_tags[block.nodes] for block in recovered.blocks
                      if block.msh==4)
            @test any(block.msh==15 && any(data.block_entities[index][column]==9 &&
                      data.external_node_tags[block.nodes[1,column]]==point
                      for column in axes(block.nodes,2))
                      for (index,block) in enumerate(recovered.blocks)
                      if block.msh==15 && 9 in data.block_entities[index])
            @test Tessella.validate(recovered).ok
        end
    finally
        A.finalize()
    end
end

@testset "invalid embedded Point constraints reject without mesh publication" begin
    for (point,missing_cell) in (((2.,0.3,0.4),false),((0.5,0.5,0.),false),
                                 ((0.25,0.3,0.4),true))
        A.initialize()
        try
            A.option("General.Terminal",0)
            cube(point)
            A.mesh.generate(2)
            A.option("Mesh.MeshOnlyEmpty",1)
            missing_cell && A.mesh.remove_elements(0,9)
            before=snapshot()
            @test_throws ArgumentError A.mesh.generate(3)
            @test snapshot()==before
        finally
            A.finalize()
        end
    end
end
end
