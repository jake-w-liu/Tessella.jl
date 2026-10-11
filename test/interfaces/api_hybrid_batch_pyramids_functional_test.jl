module HybridBatchPyramidFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())

function cube(;offset=0,region=1,x=0.0,quad_faces=Tuple(1:6))
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (index,p) in enumerate(points)
        A.model.add_point(x+p[1],p[2],p[3];tag=offset+index,meshSize=1.0)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (index,ends) in enumerate(edges)
        A.model.add_line(offset+ends[1],offset+ends[2];tag=offset+index)
        A.mesh.set_transfinite_curve(offset+index,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (index,loop) in enumerate(loops)
        lines=index==1 ? collect(loop) : [-line for line in reverse(loop)]
        tags=[sign(line)*(offset+abs(line)) for line in lines]
        A.model.add_curve_loop(tags;tag=offset+index)
        A.model.add_plane_surface([offset+index];tag=offset+index)
        A.mesh.set_transfinite_surface(offset+index)
        index in quad_faces && A.mesh.set_recombine(2,offset+index)
    end
    A.model.add_surface_loop(collect(offset+1:offset+6);tag=offset+1)
    A.model.add_volume([offset+1];tag=region)
    return collect(offset+1:offset+6)
end

function ordered_cells(cells)
    kinds,tags,vertices=cells
    sort!([(Int(kind),Tuple(cell)) for (kind,ids,nodes) in zip(kinds,tags,vertices)
            for cell in eachcol(reshape(nodes,E.msh_num_nodes(kind),length(ids)))])
end

function public_cells()
    result=Dict{UInt64,Tuple}()
    for (dim,owner) in A.model.get_entities()
        kinds,tags,vertices=A.mesh.get_elements(dim,owner)
        for (kind,ids,nodes) in zip(kinds,tags,vertices)
            for (uid,cell) in zip(ids,eachcol(reshape(nodes,E.msh_num_nodes(kind),length(ids))))
                result[uid]=(Int(kind),Int(owner),Tuple(cell))
            end
        end
    end
    return result
end

function check_region(region,surfaces,x,positions,lower,factor)
    kinds,tags,vertices=A.mesh.get_elements(3,region)
    cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
               for (kind,ids,nodes) in zip(kinds,tags,vertices))
    @test Set(keys(cells))==Set([4,7])
    expected_count=factor==1.0 ? 1 : 6
    @test size(cells[4],2)>0 && size(cells[7],2)==expected_count
    apices=vec(cells[7][5,:])
    @test allunique(apices) && Set(apices)==Set(first(A.mesh.get_nodes(3,region)))
    @test all(x<positions[uid][1]<x+1 && all(0<value<1 for value in positions[uid][2:3])
              for uid in apices)
    function tet_volume(nodes)
        a,b,c,d=(positions[uid] for uid in nodes)
        det(hcat(b-a,c-a,d-a))/6
    end
    tets=[tet_volume(cell) for cell in eachcol(cells[4])]
    pyramids=[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]])
              for cell in eachcol(cells[7])]
    @test all(>(0),tets) && all(>(0),pyramids)
    @test sum(tets)+sum(pyramids)≈1.0 atol=2e-10
    if factor<0
        delta=abs(factor)*sqrt(2.0)
        expected=([x+.5,.5,delta],[x+.5,.5,1-delta],
                  [x+.5,delta,.5],[x+1-delta,.5,.5],
                  [x+.5,1-delta,.5],[x+delta,.5,.5])
        @test all(count(uid->isapprox(positions[uid],point;atol=32eps(),rtol=0),apices)==1
                  for point in expected)
        @test all(volume->isapprox(volume,delta/3;atol=32eps(),rtol=0),pyramids)
    elseif factor==1.0
        # Independent Gmsh4.15.2 one-apex observation on this translated
        # unit-cube context, with iterative termination tolerance.
        @test positions[only(apices)]≈[x+.5,.5,.3808238751460152] atol=5e-5
    end
    incidence=Dict{Tuple,Int}()
    for (kind,nodes) in cells
        faces=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
            ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
        for cell in eachcol(nodes),face in faces
            key=Tuple(sort(cell[collect(face)]))
            incidence[key]=get(incidence,key,0)+1
        end
    end
    exterior=Set{Tuple}();quads=Set{Tuple}()
    for surface in surfaces
        kinds,tags,vertices=lower[(2,surface)]
        for (kind,ids,nodes) in zip(kinds,tags,vertices),cell in eachcol(reshape(nodes,E.msh_num_nodes(kind),length(ids)))
            face=Tuple(sort(cell));push!(exterior,face)
            kind==3 && push!(quads,face)
        end
    end
    @test Set(face for (face,count) in incidence if count==1)==exterior
    @test all(count->count in (1,2),values(incidence))
    @test Set(Tuple(sort(cell[1:4])) for cell in eachcol(cells[7]))==quads
    println("HYBRID_BATCH factor=",factor," Region=",region," Tet4=",length(tets),
            " Pyramid5=",length(pyramids)," volume=",sum(tets)+sum(pyramids));flush(stdout)
    return Dict(uid=>positions[uid] for uid in apices)
end

function run_batch(factor)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.Renumber",0)
        A.option("Mesh.MeshOnlyEmpty",1)
        quad_faces=factor==1.0 ? (1,) : Tuple(1:6)
        faces1=cube(;quad_faces)
        faces2=cube(;offset=100,region=2,x=3.,quad_faces)
        A.mesh.generate(2)
        lower=Dict(owner=>A.mesh.get_elements(owner...) for owner in A.model.get_entities() if owner[1]<3)
        lower_ids,lower_xyz,_=A.mesh.get_nodes()
        A.option("Mesh.OptimizePyramids",factor)
        A.mesh.generate(3)
        @test A.option("Mesh.OptimizePyramids")==factor
        @test A.model.get_entities(3)==[(3,1),(3,2)]
        @test all(ordered_cells(A.mesh.get_elements(owner...))==ordered_cells(cells)
                  for (owner,cells) in lower)
        ids,xyz,_=A.mesh.get_nodes()
        positions=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        @test all(positions[uid]==lower_xyz[3row-2:3row] for (row,uid) in enumerate(lower_ids))
        left=check_region(1,faces1,0.,positions,lower,factor)
        right=check_region(2,faces2,3.,positions,lower,factor)
        @test isempty(intersect(keys(left),keys(right)))
        @test isempty(intersect(Set(lower_ids),union(keys(left),keys(right))))
        @test allunique(ids) && A.mesh.get_max_node_tag()>=maximum(ids)
        @test allunique(vcat(A.mesh.get_elements()[2]...))
        @test Tessella.validate(A.mesh.get()).ok
        mktempdir() do directory
            path=joinpath(directory,"batch.msh")
            view=A.view.add("batch x")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,
                                  [[positions[uid][1]] for uid in ids])
            A.view.write(view,path)
            recovered=E.read_mixed_msh(path);data=recovered.entity_data
            @test Set(data.external_node_tags)==Set(ids)
            @test all(data.node_entities[findfirst(==(uid),data.external_node_tags)]==(3,Int32(region))
                      for (region,apices) in ((1,left),(2,right)) for uid in keys(apices))
            @test Dict(uid=>(Int(block.msh),Int(owners[column]),
                            Tuple(data.external_node_tags[block.nodes[:,column]]))
                        for (block,tags,owners) in zip(recovered.blocks,data.external_element_tags,data.block_entities)
                        for (column,uid) in enumerate(tags))==public_cells()
            @test Tessella.validate(recovered).ok
        end
        return merge(left,right)
    finally
        A.finalize()
    end
end

@testset "separated Region batch uses real negative normal births" begin
    run_batch(-0.1)
end
@testset "separated Region batch runs Winslow apex optimization" begin
    run_batch(1.0)
end
function adjacent_cubes(;shared_quad=false)
    cube(;quad_faces=shared_quad ? Tuple(1:6) : (1,2,3,5,6))
    for (tag,p) in zip(101:104,((2,0,0),(2,1,0),(2,0,1),(2,1,1)))
        A.model.add_point(p...;tag,meshSize=1.)
    end
    for (tag,ends) in zip(101:108,((2,101),(101,102),(102,3),(6,103),(103,104),(104,7),(101,103),(102,104)))
        A.model.add_line(ends...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    for (tag,curves) in zip(101:105,((101,102,103,-2),(104,105,106,-6),
                                  (101,107,-104,-10),(102,108,-105,-107),(103,11,-106,-108)))
        A.model.add_curve_loop(collect(curves);tag)
        A.model.add_plane_surface([tag];tag)
        A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop([101,102,103,104,105,-4];tag=101)
    A.model.add_volume([101];tag=2)
end

@testset "shared-Quad batch guard and outside-apex failure are atomic" begin
    for shared in (false,true)
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.option("Mesh.MeshOnlyEmpty",1)
            if shared
                adjacent_cubes(;shared_quad=true)
            else
                cube();cube(offset=100,region=2,x=3.)
            end
            A.mesh.generate(2)
            before=snapshot()
            A.option("Mesh.OptimizePyramids",shared ? 1.0 : -1.0)
            exception=try
                A.mesh.generate(3)
                nothing
            catch caught
                caught
            end
            @test exception isa ArgumentError
            @test exception isa ArgumentError && occursin(shared ? "shared Quad" : "outside",sprint(showerror,exception))
            @test snapshot()==before
        finally
            A.finalize()
        end
    end
end
end
