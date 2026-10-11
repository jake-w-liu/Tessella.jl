module HybridEmbeddedCurveFunctionalChecks
using Test, LinearAlgebra, Tessella
const A = Tessella.API
const E = Tessella.Elements

function cube(; ends=((0.25,0.3,0.4),(0.75,0.6,0.55)))
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
             (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(corners)
        A.model.add_point(p...;tag,meshSize=1.)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,pair) in enumerate(edges)
        A.model.add_line(pair...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        inward=tag>1 ? [-line for line in reverse(lines)] : collect(lines)
        A.model.add_curve_loop(inward;tag)
        A.model.add_plane_surface([tag];tag)
        A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    for (tag,p) in zip((101,102),ends)
        A.model.add_point(p...;tag,meshSize=1.)
    end
    A.model.add_line(101,102;tag=101)
    A.mesh.set_transfinite_curve(101,3)
    A.mesh.embed(1,[101],3,1)
end

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
cell_records(dim,tag)=begin
    kinds,labels,nodes=A.mesh.get_elements(dim,tag)
    Dict(uid=>(Int(kind),Tuple(matrix[:,column]))
         for (kind,uids,vertices) in zip(kinds,labels,nodes)
         for matrix in (reshape(vertices,E.msh_num_nodes(kind),length(uids)),)
         for (column,uid) in enumerate(uids))
end

function edge_link(cells,edge)
    links=Set(Tuple(sort(Int[uid for uid in cell if !(uid in edge)]))
              for cell in eachcol(cells) if all(uid->uid in cell,edge))
    degrees=Dict{Int,Int}()
    for link in links,node in link
        degrees[node]=get(degrees,node,0)+1
    end
    reached=Set{Int}()
    isempty(links) || push!(reached,first(first(links)))
    changed=true
    while changed
        changed=false
        for link in links
            any(node->node in reached,link) || continue
            for node in link
                node in reached || (push!(reached,node);changed=true)
            end
        end
    end
    return length(links)>=3 && all(==(2),values(degrees)) && reached==Set(keys(degrees))
end

@testset "retained embedded Curve segments remain conforming hybrid Tet edges" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.OptimizePyramids",-0.1)
        cube()
        A.mesh.generate(2)
        lower=cell_records(1,101)
        @test length(lower)==2 && all(first(record)==1 for record in values(lower))
        curve_uids=Set(uid for (_,vertices) in values(lower) for uid in vertices)
        @test length(curve_uids)==3
        source=Dict(owner=>A.mesh.get_nodes(owner...)[1:2] for owner in ((0,101),(0,102),(1,101)))
        @test length(source[(0,101)][1])==length(source[(0,102)][1])==length(source[(1,101)][1])==1
        @test source[(1,101)][2]≈[0.5,0.45,0.475] atol=2e-12
        A.option("Mesh.MeshOnlyEmpty",1)
        A.mesh.generate(3)
        @test cell_records(1,101)==lower
        @test all(A.mesh.get_nodes(owner...)[1:2]==record for (owner,record) in source)
        kinds,labels,vertices=A.mesh.get_elements(3,1)
        @test Set(kinds)==Set(Int32[4,7])
        cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
                   for (kind,ids,nodes) in zip(kinds,labels,vertices))
        @test size(cells[7],2)==6
        @test all(uid->uid in cells[4] && !(uid in cells[7]),curve_uids)
        @test all(edge_link(cells[4],record[2]) for record in values(lower))
        ids,xyz,_=A.mesh.get_nodes()
        positions=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        tet_volume(nodes)=begin
            a,b,c,d=(positions[uid] for uid in nodes)
            det(hcat(b-a,c-a,d-a))/6
        end
        volumes=[tet_volume(cell) for cell in eachcol(cells[4])]
        append!(volumes,[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]])
                         for cell in eachcol(cells[7])])
        @test all(>(0),volumes)
        @test sum(volumes)≈1.0 atol=2e-10
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
        @test all(isdisjoint(curve_uids,face) for (face,count) in incidence if count==1)
        @test Tessella.validate(A.mesh.get()).ok
        println("HYBRID_EMBEDDED_CURVE Line2=2 Tet4=",size(cells[4],2),
                " Pyramid5=",size(cells[7],2)," volume=",sum(volumes));flush(stdout)
        mktempdir() do directory
            path=joinpath(directory,"embedded-curve.msh")
            view=A.view.add("retained embedded Curve")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,
                                  [[positions[uid][1]] for uid in ids])
            A.view.write(view,path)
            recovered=E.read_mixed_msh(path)
            data=recovered.entity_data
            @test Set(data.external_node_tags)==Set(ids)
            for (owner,(tags,coordinates)) in source
                row=findfirst(==(only(tags)),data.external_node_tags)
                @test row!==nothing && data.node_entities[row]==(owner[1],Int32(owner[2]))
                @test recovered.coords[:,row]≈coordinates atol=8eps(Float64) rtol=0
            end
            actual=Dict(uid=>(Int(block.msh),Tuple(data.external_node_tags[block.nodes[:,column]]))
                        for (index,block) in enumerate(recovered.blocks) if block.msh==1
                        for (column,uid) in enumerate(data.external_element_tags[index])
                        if data.block_entities[index][column]==101)
            @test actual==lower
            @test all(uid->any(uid in data.external_node_tags[block.nodes] for block in recovered.blocks
                              if block.msh==4),curve_uids)
            @test Tessella.validate(recovered).ok
        end
    finally
        A.finalize()
    end
end

@testset "invalid embedded Curve constraints reject atomically" begin
    for (ends,missing_cell) in ((((0.25,0.3,0.4),(1.25,0.6,0.55)),false),
                                (((0.25,0.3,0.),(0.75,0.6,0.)),false),
                                (((0.25,0.3,0.4),(0.75,0.6,0.55)),true))
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.option("Mesh.OptimizePyramids",-0.1)
            cube(;ends)
            A.mesh.generate(2)
            A.option("Mesh.MeshOnlyEmpty",1)
            if missing_cell
                removed=first(first(A.mesh.get_elements(1,101)[2]))
                A.mesh.remove_elements(1,101,[removed])
            end
            before=snapshot()
            @test_throws ArgumentError A.mesh.generate(3)
            @test snapshot()==before
        finally
            A.finalize()
        end
    end
end
end
