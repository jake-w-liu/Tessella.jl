module HybridMultiConstraintsFunctionalChecks
using Test, LinearAlgebra, Tessella
const A = Tessella.API
const E = Tessella.Elements

const INTERIOR=((0.24,0.3,0.4),(0.64,0.3,0.45),
                (0.3,0.68,0.4),(0.72,0.7,0.58),
                (0.7,0.45,0.65),(0.35,0.5,0.65))

function cube(;interior=INTERIOR)
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
    for (tag,p) in zip(101:106,interior)
        A.model.add_point(p...;tag,meshSize=1.)
    end
    for (tag,pair) in zip((101,201),((101,102),(103,104)))
        A.model.add_line(pair...;tag)
        A.mesh.set_transfinite_curve(tag,3)
    end
    A.mesh.embed(1,[101,201],3,1)
    A.mesh.embed(0,[105,106],3,1)
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

@testset "complete retained two-Point two-Curve inventory reaches one hybrid Region" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.OptimizePyramids",-0.1)
        cube()
        @test Set(A.mesh.get_embedded(3,1))==Set(((0,105),(0,106),(1,101),(1,201)))
        A.mesh.generate(2)
        owners=vcat([(0,tag) for tag in 101:106],[(1,101),(1,201)])
        lower=Dict(owner=>cell_records(owner...) for owner in owners)
        source=Dict(owner=>A.mesh.get_nodes(owner...)[1:2] for owner in owners)
        @test all(length(lower[(0,tag)])==1 && first(only(values(lower[(0,tag)])))==15 for tag in 101:106)
        @test all(length(lower[(1,tag)])==2 && all(first(record)==1 for record in values(lower[(1,tag)]))
                  for tag in (101,201))
        @test all(length(first(record))==1 for record in values(source))
        @test source[(1,101)][2]≈[0.44,0.3,0.425] atol=2e-12
        @test source[(1,201)][2]≈[0.51,0.69,0.49] atol=2e-12
        constraints=Set(vcat(first(A.mesh.get_nodes(0,105)),first(A.mesh.get_nodes(0,106)),
                            [uid for tag in (101,201) for (_,nodes) in values(lower[(1,tag)]) for uid in nodes]))
        @test length(constraints)==8
        A.option("Mesh.MeshOnlyEmpty",1)
        A.mesh.generate(3)
        @test all(cell_records(owner...)==records for (owner,records) in lower)
        @test all(A.mesh.get_nodes(owner...)[1:2]==records for (owner,records) in source)
        kinds,labels,vertices=A.mesh.get_elements(3,1)
        @test Set(kinds)==Set(Int32[4,7])
        cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
                   for (kind,ids,nodes) in zip(kinds,labels,vertices))
        @test size(cells[7],2)==6
        @test all(uid->uid in cells[4] && !(uid in cells[7]),constraints)
        @test all(closed_edge_link(cells[4],last(record))
                  for tag in (101,201) for record in values(lower[(1,tag)]))
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
            for cell in eachcol(nodes),roles in faces
                key=Tuple(sort(cell[collect(roles)]))
                incidence[key]=get(incidence,key,0)+1
            end
        end
        quads=Set(Tuple(sort(only(last(A.mesh.get_elements(2,surface))))) for surface in 1:6)
        @test Set(face for (face,count) in incidence if count==1)==quads
        @test all(count->count in (1,2),values(incidence))
        @test all(isdisjoint(constraints,face) for (face,count) in incidence if count==1)
        @test Tessella.validate(A.mesh.get()).ok
        println("HYBRID_MULTI_CONSTRAINTS Point=2 Curve=2 Line2=4 Tet4=",size(cells[4],2),
                " Pyramid5=",size(cells[7],2)," volume=",sum(volumes));flush(stdout)
        mktempdir() do directory
            path=joinpath(directory,"hybrid-constraint-collection.msh")
            view=A.view.add("complete embedded inventory")
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

@testset "incomplete or outside multi-entity constraints reject atomically" begin
    for failure in (:outside_point,:boundary_curve,:missing_point,:missing_line)
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.option("Mesh.OptimizePyramids",-0.1)
            interior=collect(INTERIOR)
            failure==:outside_point && (interior[5]=(1.2,0.45,0.65))
            if failure==:boundary_curve
                interior[3]=(0.3,0.68,0.);interior[4]=(0.72,0.7,0.)
            end
            cube(;interior)
            A.mesh.generate(2)
            A.option("Mesh.MeshOnlyEmpty",1)
            failure==:missing_point && A.mesh.remove_elements(0,105)
            if failure==:missing_line
                removed=first(first(A.mesh.get_elements(1,201)[2]))
                A.mesh.remove_elements(1,201,[removed])
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
