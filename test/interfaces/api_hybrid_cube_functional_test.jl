module HybridCubeFunctionalChecks
using Test, LinearAlgebra, Tessella
const A=Tessella.API
const E=Tessella.Elements

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())

# Build actual public geometry, with one Quad4 per Face and no transfinite
# Volume attribute. The Region therefore needs the ordinary hybrid producer.
function cube(;offset=0,region=1,x=0.0)
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
    for (index,lines) in enumerate(loops)
        tags=[sign(line)*(offset+abs(line)) for line in lines]
        A.model.add_curve_loop(tags;tag=offset+index)
        A.model.add_plane_surface([offset+index];tag=offset+index)
        A.mesh.set_transfinite_surface(offset+index)
        A.mesh.set_recombine(2,offset+index)
    end
    A.model.add_surface_loop(collect(offset+1:offset+6);tag=offset+1)
    A.model.add_volume([offset+1];tag=region)
    return collect(offset+1:offset+6)
end

function volume_cells(region)
    families,labels,nodes=A.mesh.get_elements(3,region)
    result=Dict{Int,Matrix{UInt64}}()
    for (kind,ids,vertices) in zip(families,labels,nodes)
        width=kind==4 ? 4 : kind==7 ? 5 : error("unexpected volume family $kind")
        result[Int(kind)]=reshape(vertices,width,length(ids))
    end
    return result
end

function independent_volumes(cells,points)
    tets=Float64[];pyramids=Float64[]
    function determinant(nodes)
        a,b,c,d=(points[uid] for uid in nodes)
        return det(hcat(b-a,c-a,d-a))/6
    end
    for cell in eachcol(cells[4])
        push!(tets,determinant(cell))
    end
    for cell in eachcol(cells[7])
        # The two signed tetrahedra partition a cyclic Pyramid5 base.
        push!(pyramids,determinant(cell[[1,2,3,5]])+
                       determinant(cell[[1,3,4,5]]))
    end
    return tets,pyramids
end

function exterior_faces(cells)
    incidence=Dict{Tuple,Int}()
    for (kind,rows) in cells
        faces=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
                       ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
        for cell in eachcol(rows),face in faces
            key=Tuple(sort(cell[collect(face)]))
            incidence[key]=get(incidence,key,0)+1
        end
    end
    return incidence
end

function serialized_cells(mesh)
    data=mesh.entity_data
    return Dict(uid=>(Int(block.msh),Int(owner),
                     Tuple(data.external_node_tags[block.nodes[:,column]]))
                for (block,ids,owners) in zip(mesh.blocks,data.external_element_tags,data.block_entities)
                for (column,(uid,owner)) in enumerate(zip(ids,owners)))
end

function public_cells()
    result=Dict{UInt64,Tuple}()
    for (dim,owner) in A.model.get_entities()
        families,labels,nodes=A.mesh.get_elements(dim,owner)
        for (kind,uids,vertices) in zip(families,labels,nodes)
            cells=reshape(vertices,E.msh_num_nodes(kind),length(uids))
            for (column,uid) in enumerate(uids)
                result[uid]=(Int(kind),Int(owner),Tuple(cells[:,column]))
            end
        end
    end
    return result
end

@testset "OptimizePyramids public integer option" begin
    A.initialize()
    try
        @test A.option("Mesh.OptimizePyramids")==0.0
        A.option("Mesh.OptimizePyramids",0)
        @test A.option("Mesh.OptimizePyramids")==0.0
        A.option("Mesh.OptimizePyramids",-1.7)
        @test A.option("Mesh.OptimizePyramids")==-1.0
        @test_throws ArgumentError A.option("Mesh.OptimizePyramids",NaN)
        @test A.option("Mesh.OptimizePyramids")==-1.0
        A.option("Mesh.OptimizePyramids",0)
    finally
        A.finalize()
    end
end

@testset "public coarse cube produces conforming positive Tet4 and Pyramid5" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        surfaces=cube()
        A.mesh.generate(3)
        @test A.model.get_entities(3)==[(3,1)]
        cells=volume_cells(1)
        @test Set(keys(cells))==Set([4,7])
        @test size(cells[4],2)>0 && size(cells[7],2)==6
        ids,xyz,_=A.mesh.get_nodes()
        points=Dict(uid=>xyz[3index-2:3index] for (index,uid) in enumerate(ids))
        tets,pyramids=independent_volumes(cells,points)
        @test all(>(0),tets) && all(>(0),pyramids)
        @test sum(tets)+sum(pyramids)≈1.0 atol=2e-10
        println("HYBRID_CUBE Tet4=",length(tets)," Pyramid5=",length(pyramids),
                " volume=",sum(tets)+sum(pyramids));flush(stdout)
        incidence=exterior_faces(cells)
        @test all(count->count in (1,2),values(incidence))
        quads=Tuple[]
        for surface in surfaces
            kinds,tags,nodes=A.mesh.get_elements(2,surface)
            @test kinds==Int32[3] && length(only(tags))==1
            push!(quads,Tuple(sort(only(nodes))))
        end
        @test Set(face for (face,count) in incidence if count==1)==Set(quads)
        @test Set(Tuple(sort(cell[1:4])) for cell in eachcol(cells[7]))==Set(quads)

        before=snapshot()
        detached=A.mesh.get()
        @test detached isa E.MixedMesh && Tessella.validate(detached).ok
        pristine=copy(detached.coords)
        detached.coords[1,1]+=10
        @test A.mesh.get().coords==pristine
        @test snapshot()==before
        # The detached geometry carrier intentionally omits entity metadata.
        # A real NodeData view writer exports the classified session mesh;
        # compare that output with explicit public UID/owner queries.
        mktempdir() do directory
            path=joinpath(directory,"hybrid-cube.msh")
            view=A.view.add("hybrid cube x")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,
                                  [[points[uid][1]] for uid in ids])
            A.view.write(view,path)
            recovered=E.read_mixed_msh(path)
            @test Tessella.validate(recovered).ok
            @test recovered.entity_data!==nothing
            @test serialized_cells(recovered)==public_cells()
            original=Dict(uid=>Tuple(points[uid]) for uid in ids)
            restored=Dict(uid=>Tuple(recovered.coords[:,row])
                          for (row,uid) in enumerate(recovered.entity_data.external_node_tags))
            @test Set(keys(restored))==Set(keys(original))
            # The native view ASCII writer uses %.16g; UID and ownership
            # round trips are exact, while decimal XYZ may lose a few ULPs.
            error=maximum(abs(restored[uid][axis]-original[uid][axis])
                          for uid in ids for axis in 1:3)
            @test error<=8eps(Float64)
            println("HYBRID_MSH_MAX_COORD_ERROR ",error);flush(stdout)
        end
        @test snapshot()==before
    finally
        A.finalize()
    end
end

@testset "public hybrid batch rejects unsupported pyramid optimization atomically" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube()
        cube(offset=100,region=2,x=3.0)
        A.mesh.generate(2)
        before=snapshot()
        A.option("Mesh.OptimizePyramids",1)
        exception=try
            A.mesh.generate(3)
            nothing
        catch caught
            caught
        end
        @test exception isa ArgumentError
        @test exception isa ArgumentError && occursin("Mesh.OptimizePyramids=0",sprint(showerror,exception))
        @test snapshot()==before
        @test A.option("Mesh.OptimizePyramids")==1.0
        A.option("Mesh.OptimizePyramids",0)
        A.mesh.generate(3)
        @test A.model.get_entities(3)==[(3,1),(3,2)]
        @test all(region->Set(keys(volume_cells(region)))==Set([4,7]),(1,2))
        @test all(region->size(volume_cells(region)[7],2)==6,(1,2))
        println("HYBRID_BATCH_COUNTS ",[(region,size(volume_cells(region)[4],2),
                 size(volume_cells(region)[7],2)) for region in (1,2)]);flush(stdout)
    finally
        A.finalize()
    end
end
end
