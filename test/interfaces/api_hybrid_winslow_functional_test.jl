module HybridWinslowFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements
const W=Tessella.Model.ModelNativeHybridWinslow

@testset "regularized Winslow analytic gradient and constrained solve" begin
    points=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.))
    normals=((-1.,-1.,-1.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.))
    patch=[W.PatchTet(points,normals,4)]
    point=(0.2,0.3,0.4);epsilon=1.e-6
    function independent_energy(p)
        determinant=p[3]
        chi=(determinant+sqrt(determinant^2+epsilon^2))/2
        (2+sum(abs2,p))/chi^(2/3)+(determinant^2+1)/chi
    end
    value,gradient,_=W.energy_gradient(point,patch,epsilon)
    @test value≈independent_energy(point) rtol=2e-14
    for axis in 1:3
        h=1.e-6
        plus=ntuple(i->point[i]+(i==axis ? h : 0),3)
        minus=ntuple(i->point[i]-(i==axis ? h : 0),3)
        @test gradient[axis]≈(independent_energy(plus)-independent_energy(minus))/(2h) rtol=2e-8
    end
    inverted=W.energy_gradient((0.2,0.3,-0.4),patch,epsilon)
    @test all(isfinite,(inverted[1],inverted[2]...)) && inverted[3]<0
    function negative_energy(p)
        inverse_chi=2*(sqrt(p[3]^2+epsilon^2)-p[3])/epsilon^2
        (2+sum(abs2,p))*inverse_chi^(2/3)+(p[3]^2+1)*inverse_chi
    end
    far=(0.2,0.3,-100.0);h=1.e-3
    far_gradient=W.energy_gradient(far,patch,epsilon)[2][3]
    difference=(negative_energy((far[1],far[2],far[3]+h))-
                negative_energy((far[1],far[2],far[3]-h)))/(2h)
    @test far_gradient≈difference rtol=2e-8
    optimal,result=W.minimize((0.2,0.3,0.1),patch,epsilon)
    @test collect(optimal)≈[0.,0.,1.] atol=1e-4
    @test result[1]≈5.0 atol=1e-7
    @test norm(collect(result[2]))<1e-4
end

function cube(;quad_faces=(1,))
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
             (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,point) in enumerate(corners)
        A.model.add_point(point...;tag,meshSize=1.0)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(edges)
        A.model.add_line(ends...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,loop) in enumerate(loops)
        lines=tag==1 ? collect(loop) : [-line for line in reverse(loop)]
        A.model.add_curve_loop(lines;tag)
        A.model.add_plane_surface([tag];tag)
        A.mesh.set_transfinite_surface(tag)
        tag in quad_faces && A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
end

function run_cube(quad_faces)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.MeshOnlyEmpty",1)
        A.option("Mesh.Renumber",0)
        cube(;quad_faces)
        A.mesh.generate(2)
        lower_nodes=A.mesh.get_nodes()[1:2]
        lower=Dict((dim,tag)=>A.mesh.get_elements(dim,tag)
                   for (dim,tag) in A.model.get_entities() if dim<3)
        exterior=Set{Tuple}()
        for surface in 1:6
            families,labels,vertices=lower[(2,surface)]
            for (kind,ids,nodes) in zip(families,labels,vertices)
                cells=reshape(nodes,E.msh_num_nodes(kind),length(ids))
                union!(exterior,(Tuple(sort(cell)) for cell in eachcol(cells)))
            end
        end
        node_clock=A.mesh.get_max_node_tag()
        A.option("Mesh.MeshOnlyEmpty",1)
        A.option("Mesh.OptimizePyramids",1.0)
        A.mesh.generate(3)
        @test A.option("Mesh.OptimizePyramids")==1.0
        function records(cells)
            kinds,tags,connectivity=cells
            result=Tuple[]
            for (kind,ids,nodes) in zip(kinds,tags,connectivity)
                for corners in eachcol(reshape(nodes,E.msh_num_nodes(kind),length(ids)))
                    push!(result,(kind,Tuple(corners)))
                end
            end
            return sort!(result)
        end
        # Native generate(3) can recreate lower element tags even when only
        # empty entities are requested. Their owners and ordered corner UIDs
        # remain fixed; presentation order does not select a different cell.
        @test all(records(A.mesh.get_elements(owner...))==records(cells) for (owner,cells) in lower)
        families,labels,vertices=A.mesh.get_elements(3,1)
        @test Set(families)==Set(Int32[4,7])
        cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
                   for (kind,ids,nodes) in zip(families,labels,vertices))
        @test size(cells[7],2)==length(quad_faces)
        ids,xyz,_=A.mesh.get_nodes()
        positions=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        @test all(positions[uid]==lower_nodes[2][3row-2:3row]
                  for (row,uid) in enumerate(lower_nodes[1]))
        apices=vec(cells[7][5,:])
        @test Set(apices)==Set(node_clock+UInt64(index) for index in eachindex(quad_faces))
        @test all(uid in first(A.mesh.get_nodes(3,1)) for uid in apices)
        @test all(all(0.0<coordinate<1.0 for coordinate in positions[uid]) for uid in apices)
        if length(quad_faces)==1
            # Independent native 4.15.2 observation on this exact one-apex
            # context; this tolerance reflects iterative minimizer termination.
            @test positions[only(apices)]≈[0.5,0.5,0.3808238751460152] atol=5e-5
        end
        function tet_volume(nodes)
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
        @test Set(face for (face,count) in incidence if count==1)==exterior
        @test all(count->count in (1,2),values(incidence))
        @test Tessella.validate(A.mesh.get()).ok
        println("HYBRID_WINSLOW apices=",length(apices)," Tet4=",size(cells[4],2),
                " Pyramid5=",size(cells[7],2)," volume=",sum(volumes),
                " apexXYZ=",[positions[uid] for uid in apices]);flush(stdout)
        mktempdir() do directory
            path=joinpath(directory,"winslow.msh")
            view=A.view.add("Winslow")
            A.view.add_model_data(view,0,A.model.get_current(),"NodeData",ids,
                                  [[positions[uid][1]] for uid in ids])
            A.view.write(view,path)
            recovered=E.read_mixed_msh(path);data=recovered.entity_data
            @test Set(data.external_node_tags)==Set(ids)
            @test all(data.node_entities[findfirst(==(uid),data.external_node_tags)]==(3,Int32(1))
                      for uid in apices)
            @test Set(uid for tags in data.external_element_tags for uid in tags)==
                  Set(uid for tags in A.mesh.get_elements()[2] for uid in tags)
            @test Tessella.validate(recovered).ok
        end
    finally
        A.finalize()
    end
end

@testset "public one-Pyramid Winslow relocation preserves the conforming cube" begin
    run_cube((1,))
end
@testset "public six-Pyramid Winslow relocation preserves the conforming cube" begin
    run_cube(Tuple(1:6))
end
end
