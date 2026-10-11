module HybridPolygonSheetFunctionalChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const E=Tessella.Elements
const POLYGONS=(
    ((.25,.25,.45),(.625,.25,.45),(.75,.5,.45),(.625,.75,.45),(.25,.75,.45)),
    ((.25,.375,.45),(.375,.25,.45),(.625,.25,.45),(.75,.375,.45),(.625,.75,.45),(.375,.75,.45)))

function cube(sheet)
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,point) in enumerate(corners);A.model.add_point(point...;tag,meshSize=1.);end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,pair) in enumerate(edges)
        A.model.add_line(pair...;tag);A.mesh.set_transfinite_curve(tag,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        inward=tag>1 ? [-line for line in reverse(lines)] : collect(lines)
        A.model.add_curve_loop(inward;tag);A.model.add_plane_surface([tag];tag)
        A.mesh.set_transfinite_surface(tag);A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1);A.model.add_volume([1];tag=1)
    count=length(sheet)
    for (index,point) in enumerate(sheet);A.model.add_point(point...;tag=100+index,meshSize=1.);end
    for index in 1:count
        tag=100+index
        A.model.add_line(tag,100+mod1(index+1,count);tag);A.mesh.set_transfinite_curve(tag,2)
    end
    A.model.add_curve_loop(collect(101:100+count);tag=101)
    A.model.add_plane_surface([101];tag=101)
    A.mesh.set_algorithm(2,101,3)
    A.mesh.embed(2,[101],3,1)
    A.model.add_physical_group(2,[101];tag=61,name="sheet")
    A.model.add_physical_group(3,[1];tag=62,name="material")
end
function records(dim,tag)
    families,labels,nodes=A.mesh.get_elements(dim,tag)
    Dict(uid=>(Int(family),Tuple(matrix[:,column])) for (family,uids,vertices) in zip(families,labels,nodes)
        for matrix in (reshape(vertices,E.msh_num_nodes(family),length(uids)),)
        for (column,uid) in enumerate(uids))
end
snapshot()=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
function edge_link(cells,edge)
    links=Set(Tuple(sort(UInt64[uid for uid in cell if !(uid in edge)]))
        for cell in eachcol(cells) if all(uid->uid in cell,edge))
    graph=Dict{UInt64,Vector{UInt64}}()
    for (a,b) in links;push!(get!(graph,a,UInt64[]),b);push!(get!(graph,b,UInt64[]),a);end
    isempty(graph) && return false
    pending=UInt64[first(keys(graph))];reached=Set{UInt64}()
    while !isempty(pending)
        uid=pop!(pending);uid in reached && continue
        push!(reached,uid);append!(pending,graph[uid])
    end
    length(links)>=3 && all(neighbors->length(neighbors)==2,values(graph)) && reached==Set(keys(graph))
end
function session(body)
    A.initialize()
    try
        A.option("General.Terminal",0);A.option("Mesh.OptimizePyramids",-.1)
        A.option("Mesh.Renumber",0)
        body()
    finally;A.finalize();end
end

@testset "actual pentagon and hexagon Tri3 disks remain two-sided hybrid constraints" begin
    for polygon in POLYGONS
        session() do
            cube(polygon);A.mesh.generate(2)
            count=length(polygon)
            owners=vcat([(0,tag) for tag in 101:100+count],[(1,tag) for tag in 101:100+count],[(2,101)])
            lower=Dict(owner=>records(owner...) for owner in owners)
            faces=Tuple[last(record) for record in values(lower[(2,101)])]
            @test length(faces)==count-2 && all(record->first(record)==2,values(lower[(2,101)]))
            if count==5
                @test any(face->Tuple(sort(collect(face)))!=face,faces)
            end
            @test all(length(lower[(1,tag)])==1 && first(only(values(lower[(1,tag)])))==1 for tag in 101:100+count)
            vertices=Set(uid for face in faces for uid in face)
            edges=Set(minmax(face[i],face[mod1(i+1,3)]) for face in faces for i in 1:3)
            @test length(vertices)==count && length(edges)==2count-3
            source=Dict(owner=>A.mesh.get_nodes(owner...)[1:2] for owner in owners if owner[1]==0)
            @test all(length(first(value))==1 for value in values(source))
            A.option("Mesh.MeshOnlyEmpty",1);A.mesh.generate(3)
            published_lower=Dict(owner=>records(owner...) for owner in owners)
            lower_uids=UInt64[uid for cells in values(published_lower) for uid in keys(cells)]
            # Ordinary generation renews linear cell factory identities even
            # with renumber disabled. Literal cyclic owner incidence survives.
            @test all(Set(values(published_lower[owner]))==Set(values(cells)) for (owner,cells) in lower) &&
                allunique(lower_uids) && all(>(0),lower_uids) && maximum(lower_uids)<=A.mesh.get_max_element_tag()
            faces=Tuple[last(record) for record in values(published_lower[(2,101)])]
            @test all(A.mesh.get_nodes(owner...)[1:2]==row for (owner,row) in source)
            @test A.model.get_physical_groups_for_entity(2,101)==[61] &&
                A.model.get_physical_groups_for_entity(3,1)==[62]
            families,labels,nodes=A.mesh.get_elements(3,1)
            @test Set(families)==Set(Int32[4,7])
            cells=Dict(Int(family)=>reshape(connected,E.msh_num_nodes(family),length(uids))
                for (family,uids,connected) in zip(families,labels,nodes))
            @test size(cells[7],2)==6 && all(uid->uid in cells[4] && !(uid in cells[7]),vertices)
            @test all(edge_link(cells[4],edge) for edge in edges)
            ids,coordinates,_=A.mesh.get_nodes()
            positions=Dict(uid=>coordinates[3row-2:3row] for (row,uid) in enumerate(ids))
            tet_volume(cell)=begin
                a,b,c,d=(positions[uid] for uid in cell)
                det(hcat(b-a,c-a,d-a))/6
            end
            volumes=[tet_volume(cell) for cell in eachcol(cells[4])]
            append!(volumes,[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]]) for cell in eachcol(cells[7])])
            @test all(>(0),volumes) && isapprox(sum(volumes),1.;atol=2e-10)
            for face in faces
                opposite=UInt64[only(uid for uid in cell if !(uid in face))
                    for cell in eachcol(cells[4]) if all(uid->uid in cell,face)]
                @test length(opposite)==2
                a,b,c=(positions[uid] for uid in face)
                sides=[dot(cross(b-a,c-a),positions[uid]-a) for uid in opposite]
                @test first(sides)*last(sides)<0
            end
            incidence=Dict{Tuple,Int}()
            for (family,block) in cells
                roles=family==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
                    ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
                for cell in eachcol(block),slots in roles
                    key=Tuple(sort(cell[collect(slots)]));incidence[key]=get(incidence,key,0)+1
                end
            end
            quads=Set(Tuple(sort(only(last(A.mesh.get_elements(2,surface))))) for surface in 1:6)
            @test Set(face for (face,hits) in incidence if hits==1)==quads
            @test all(hits->hits in (1,2),values(incidence)) && all(incidence[Tuple(sort(collect(face)))]==2 for face in faces)
            @test all(isdisjoint(vertices,face) for (face,hits) in incidence if hits==1)
            @test Tessella.validate(A.mesh.get()).ok
            println("HYBRID_POLYGON_SHEET vertices=",count," Tri3=",length(faces)," Tet4=",size(cells[4],2),
                " Pyramid5=",size(cells[7],2)," volume=",sum(volumes));flush(stdout)
        end
    end
end
@testset "one reversed actual polygon Tri rejects atomically" begin
    session() do
        cube(first(POLYGONS));A.mesh.generate(2)
        tags=A.mesh.get_elements(2,101)[2][1]
        @test length(tags)==3
        A.mesh.reverse_elements([first(tags)]);A.option("Mesh.MeshOnlyEmpty",1)
        before=snapshot()
        failure=try;A.mesh.generate(3);nothing;catch error;error isa InterruptException && rethrow();error;end
        @test failure isa ArgumentError && occursin("embedded Surface lacks its owned actual Surface constructor",sprint(showerror,failure))
        @test snapshot()==before
    end
end
end
