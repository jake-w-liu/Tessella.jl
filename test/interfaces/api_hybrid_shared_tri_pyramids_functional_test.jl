module HybridSharedTriPyramidFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements

function adjacent_cubes(;shared_quad=false)
    for (tag,p) in enumerate(((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                              (0,0,1),(1,0,1),(1,1,1),(0,1,1)))
        A.model.add_point(p...;tag,meshSize=1.)
    end
    for (tag,ends) in enumerate(((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
                                 (7,8),(8,5),(1,5),(2,6),(3,7),(4,8)))
        A.model.add_line(ends...;tag);A.mesh.set_transfinite_curve(tag,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,loop) in enumerate(loops)
        A.model.add_curve_loop(tag==1 ? collect(loop) : [-line for line in reverse(loop)];tag)
        A.model.add_plane_surface([tag];tag);A.mesh.set_transfinite_surface(tag)
        (tag!=4 || shared_quad) && A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1);A.model.add_volume([1];tag=1)
    for (tag,p) in zip(101:104,((2,0,0),(2,1,0),(2,0,1),(2,1,1)))
        A.model.add_point(p...;tag,meshSize=1.)
    end
    for (tag,ends) in zip(101:108,((2,101),(101,102),(102,3),(6,103),(103,104),(104,7),(101,103),(102,104)))
        A.model.add_line(ends...;tag);A.mesh.set_transfinite_curve(tag,2)
    end
    right=((101,102,103,-2),(104,105,106,-6),(101,107,-104,-10),
           (102,108,-105,-107),(103,11,-106,-108))
    for (index,loop) in enumerate(right)
        tag=100+index
        A.model.add_curve_loop(index==1 ? collect(loop) : [-line for line in reverse(loop)];tag)
        A.model.add_plane_surface([tag];tag);A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop([101,102,103,104,105,-4];tag=101)
    A.model.add_volume([101];tag=2)
end

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
function records(cells)
    kinds,tags,vertices=cells
    sort!([(Int(kind),Tuple(cell)) for (kind,ids,nodes) in zip(kinds,tags,vertices)
           for cell in eachcol(reshape(nodes,E.msh_num_nodes(kind),length(ids)))])
end

function check_region(region,positions,factor,lower)
    kinds,tags,vertices=A.mesh.get_elements(3,region)
    cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
                for (kind,ids,nodes) in zip(kinds,tags,vertices))
    @test Set(keys(cells))==Set([4,7]) && size(cells[7],2)==5
    x=Float64(region-1);apices=vec(cells[7][5,:])
    @test allunique(apices) && all(uid in first(A.mesh.get_nodes(3,region)) for uid in apices)
    @test all(x<positions[uid][1]<x+1 && all(0<value<1 for value in positions[uid][2:3]) for uid in apices)
    function tet_volume(nodes)
        a,b,c,d=(positions[uid] for uid in nodes)
        det(hcat(b-a,c-a,d-a))/6
    end
    tets=[tet_volume(cell) for cell in eachcol(cells[4])]
    pyramids=[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]]) for cell in eachcol(cells[7])]
    @test all(>(0),tets) && all(>(0),pyramids)
    @test sum(tets)+sum(pyramids)≈1.0 atol=2e-10
    if factor<0
        delta=abs(factor)*sqrt(2.)
        expected=([x+.5,.5,delta],[x+.5,.5,1-delta],[x+.5,delta,.5],
                  [x+.5,1-delta,.5],region==1 ? [delta,.5,.5] : [2-delta,.5,.5])
        @test all(count(uid->isapprox(positions[uid],point;atol=32eps(),rtol=0),apices)==1 for point in expected)
        @test all(volume->isapprox(volume,delta/3;atol=32eps(),rtol=0),pyramids)
    end
    incidence=Dict{Tuple,Int}()
    for (kind,nodes) in cells
        faces=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
            ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
        for cell in eachcol(nodes),face in faces
            key=Tuple(sort(cell[collect(face)]));incidence[key]=get(incidence,key,0)+1
        end
    end
    surfaces=region==1 ? collect(1:6) : [101,102,103,104,105,4]
    boundary=Set{Tuple}()
    for surface in surfaces
        for (_,corners) in records(lower[(2,surface)])
            push!(boundary,Tuple(sort(collect(corners))))
        end
    end
    @test Set(face for (face,count) in incidence if count==1)==boundary
    @test all(count->count in (1,2),values(incidence))
    println("SHARED_TRI_PYRAMIDS factor=",factor," Region=",region," Tet4=",length(tets),
            " Pyramid5=5 volume=",sum(tets)+sum(pyramids));flush(stdout)
    return apices,incidence,Set(vec(cells[4]))∪Set(vec(cells[7]))
end

function run_shared(factor)
    A.initialize()
    try
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0);A.option("Mesh.MeshOnlyEmpty",1)
        adjacent_cubes();A.mesh.generate(2)
        lower=Dict(owner=>A.mesh.get_elements(owner...) for owner in A.model.get_entities() if owner[1]<3)
        lower_ids,lower_xyz,_=A.mesh.get_nodes()
        A.option("Mesh.OptimizePyramids",factor);A.mesh.generate(3)
        @test A.option("Mesh.OptimizePyramids")==factor
        @test all(records(A.mesh.get_elements(owner...))==records(cells) for (owner,cells) in lower)
        ids,xyz,_=A.mesh.get_nodes()
        positions=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        @test all(positions[uid]==lower_xyz[3row-2:3row] for (row,uid) in enumerate(lower_ids))
        left,incidence1,nodes1=check_region(1,positions,factor,lower)
        right,incidence2,nodes2=check_region(2,positions,factor,lower)
        shared=Set(Tuple(sort(collect(corners))) for (_,corners) in records(lower[(2,4)]))
        @test length(shared)==2 && all(get(incidence1,face,0)==get(incidence2,face,0)==1 for face in shared)
        @test intersect(nodes1,nodes2)==Set(UInt64[2,3,6,7])
        @test isempty(intersect(Set(left),Set(right)))
        @test all(uid in union(nodes1,nodes2) for uid in setdiff(Set(ids),Set(lower_ids)))
        @test allunique(ids) && allunique(vcat(A.mesh.get_elements()[2]...))
        @test Tessella.validate(A.mesh.get()).ok
        mktempdir() do directory
            path=joinpath(directory,"shared.msh");view=A.view.add("shared x")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,[[positions[uid][1]] for uid in ids])
            A.view.write(view,path)
            recovered=E.read_mixed_msh(path);data=recovered.entity_data
            @test Set(data.external_node_tags)==Set(ids)
            @test all(data.node_entities[findfirst(==(uid),data.external_node_tags)]==(3,Int32(region))
                       for (region,apices) in ((1,left),(2,right)) for uid in apices)
            @test all(owners[column] in (1,2) for (block,owners) in zip(recovered.blocks,data.block_entities)
                       if E.msh_dimension(Int(block.msh))==3 for column in axes(block.nodes,2))
            @test Tessella.validate(recovered).ok
        end
    finally
        A.finalize()
    end
end

@testset "shared-Tri Regions retain their actual negative apex witnesses" begin
    run_shared(-0.1)
end
@testset "shared-Tri Regions run separate real Winslow apex stars" begin
    run_shared(1.)
end
@testset "shared Quad still rejects without a two-sided Pyramid constructor" begin
    A.initialize()
    try
        A.option("General.Terminal",0);A.option("Mesh.MeshOnlyEmpty",1)
        adjacent_cubes(;shared_quad=true);A.mesh.generate(2);before=snapshot()
        A.option("Mesh.OptimizePyramids",1.)
        exception=try
            A.mesh.generate(3);nothing
        catch caught
            caught
        end
        @test exception isa ArgumentError
        @test exception isa ArgumentError && occursin("shared Quad",sprint(showerror,exception))
        @test snapshot()==before
    finally
        A.finalize()
    end
end
end
