module HybridEmbeddedSurfaceFunctionalChecks
using Test, LinearAlgebra, Tessella
const A = Tessella.API
const E = Tessella.Elements

function cube(;sheet=((0.35,0.35,0.45),(0.65,0.35,0.45),(0.5,0.65,0.45)))
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
    for (tag,p) in zip(101:103,sheet)
        A.model.add_point(p...;tag,meshSize=1.)
    end
    for (tag,pair) in zip(101:103,((101,102),(102,103),(103,101)))
        A.model.add_line(pair...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    A.model.add_curve_loop([101,102,103];tag=101)
    A.model.add_plane_surface([101];tag=101)
    A.mesh.set_transfinite_surface(101)
    A.mesh.embed(2,[101],3,1)
end

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
function cell_records(dim,tag)
    kinds,labels,nodes=A.mesh.get_elements(dim,tag)
    return Dict(uid=>(Int(kind),Tuple(matrix[:,column]))
                for (kind,uids,vertices) in zip(kinds,labels,nodes)
                for matrix in (reshape(vertices,E.msh_num_nodes(kind),length(uids)),)
                for (column,uid) in enumerate(uids))
end

function closed_edge_link(cells,edge)
    links=Set(Tuple(sort(Int[uid for uid in cell if !(uid in edge)]))
              for cell in eachcol(cells) if all(uid->uid in cell,edge))
    degrees=Dict{Int,Int}()
    for link in links,node in link;degrees[node]=get(degrees,node,0)+1;end
    reached=Set{Int}();isempty(links) || push!(reached,first(first(links)))
    changed=true
    while changed
        changed=false
        for link in links
            any(node->node in reached,link) || continue
            for node in link;node in reached || (push!(reached,node);changed=true);end
        end
    end
    return length(links)>=3 && all(==(2),values(degrees)) && reached==Set(keys(degrees))
end

@testset "retained interior Tri3 sheet is a protected two-sided hybrid Tet facet" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.OptimizePyramids",-0.1)
        cube()
        A.mesh.generate(2)
        owners=vcat([(0,tag) for tag in 101:103],[(1,tag) for tag in 101:103],[(2,101)])
        lower=Dict(owner=>cell_records(owner...) for owner in owners)
        @test length(lower[(2,101)])==1 && first(only(values(lower[(2,101)])))==2
        @test all(length(lower[(1,tag)])==1 && first(only(values(lower[(1,tag)])))==1 for tag in 101:103)
        face=last(only(values(lower[(2,101)])))
        @test length(Set(face))==3
        source=Dict(owner=>A.mesh.get_nodes(owner...)[1:2] for owner in ((0,101),(0,102),(0,103)))
        @test all(length(first(value))==1 for value in values(source))
        A.option("Mesh.MeshOnlyEmpty",1)
        A.mesh.generate(3)
        @test all(cell_records(owner...)==records for (owner,records) in lower)
        @test all(A.mesh.get_nodes(owner...)[1:2]==record for (owner,record) in source)
        kinds,labels,vertices=A.mesh.get_elements(3,1)
        @test Set(kinds)==Set(Int32[4,7])
        cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
                   for (kind,ids,nodes) in zip(kinds,labels,vertices))
        @test size(cells[7],2)==6
        @test all(uid->uid in cells[4] && !(uid in cells[7]),face)
        @test all(closed_edge_link(cells[4],last(only(values(lower[(1,tag)])))) for tag in 101:103)
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
        opposite=UInt64[only(uid for uid in cell if !(uid in face))
                        for cell in eachcol(cells[4]) if all(uid->uid in cell,face)]
        @test length(opposite)==2
        a,b,c=(positions[uid] for uid in face)
        sides=[dot(cross(b-a,c-a),positions[uid]-a) for uid in opposite]
        @test first(sides)*last(sides)<0
        incidence=Dict{Tuple,Int}()
        for (kind,nodes) in cells
            faces=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
                ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
            for cell in eachcol(nodes),roles in faces
                key=Tuple(sort(cell[collect(roles)]))
                incidence[key]=get(incidence,key,0)+1
            end
        end
        @test incidence[Tuple(sort(collect(face)))]==2
        quads=Set(Tuple(sort(only(last(A.mesh.get_elements(2,surface))))) for surface in 1:6)
        @test Set(panel for (panel,count) in incidence if count==1)==quads
        @test all(count->count in (1,2),values(incidence))
        @test all(isdisjoint(Set(face),panel) for (panel,count) in incidence if count==1)
        @test Tessella.validate(A.mesh.get()).ok
        println("HYBRID_EMBEDDED_SURFACE Tri3=1 Tet4=",size(cells[4],2),
                " Pyramid5=",size(cells[7],2)," volume=",sum(volumes)," sheet_incidence=2");flush(stdout)
        mktempdir() do directory
            path=joinpath(directory,"embedded-surface.msh")
            view=A.view.add("retained embedded Surface")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,
                                  [[positions[uid][1]] for uid in ids])
            A.view.write(view,path)
            recovered=E.read_mixed_msh(path);data=recovered.entity_data
            @test Set(data.external_node_tags)==Set(ids)
            for (owner,(tags,coordinates)) in source
                row=findfirst(==(only(tags)),data.external_node_tags)
                @test row!==nothing && data.node_entities[row]==(owner[1],Int32(owner[2]))
                @test recovered.coords[:,row]≈coordinates atol=8eps(Float64) rtol=0
            end
            for (owner,records) in lower
                actual=Dict(uid=>(Int(block.msh),Tuple(data.external_node_tags[block.nodes[:,column]]))
                            for (index,block) in enumerate(recovered.blocks)
                            if E.msh_dimension(block.msh)==owner[1]
                            for (column,uid) in enumerate(data.external_element_tags[index])
                            if data.block_entities[index][column]==owner[2])
                @test actual==records
            end
            @test Tessella.validate(recovered).ok
        end
    finally
        A.finalize()
    end
end

@testset "invalid retained Surface constraints reject without publication" begin
    for (sheet,missing_cell) in ((((0.35,0.35,0.45),(1.25,0.35,0.45),(0.5,0.65,0.45)),false),
                                 (((0.35,0.35,0.),(0.65,0.35,0.),(0.5,0.65,0.)),false),
                                 (((0.35,0.35,0.45),(0.65,0.35,0.45),(0.5,0.65,0.45)),true))
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.option("Mesh.OptimizePyramids",-0.1)
            cube(;sheet)
            A.mesh.generate(2)
            A.option("Mesh.MeshOnlyEmpty",1)
            missing_cell && A.mesh.remove_elements(2,101)
            before=snapshot()
            @test_throws ArgumentError A.mesh.generate(3)
            @test snapshot()==before
        finally
            A.finalize()
        end
    end
end
end
