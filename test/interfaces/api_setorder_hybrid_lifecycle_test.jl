module SetOrderHybridLifecycleChecks
# Regression for the order-pass constructor receipt transport: an order-2
# elevation used to drop the retained first-order Surface constructor
# receipts because the rebound authority only accepted linear product blocks.
# A subsequent `generate(3)` then lost the hybrid boundary's constructor
# vector and threw `mesh_model_volume: hybrid boundary lacks its owned
# Surface constructor vector`.
using Test, LinearAlgebra, Tessella
const A=Tessella.API
const E=Tessella.Elements

# One Quad4 per Face — the Region must take the ordinary hybrid producer,
# which consumes the retained first-order Surface constructor receipts.
function recombined_cube()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (index,p) in enumerate(points)
        A.model.add_point(Float64.(p)...;tag=index,meshSize=1.0)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (index,ends) in enumerate(edges)
        A.model.add_line(ends[1],ends[2];tag=index)
        A.mesh.set_transfinite_curve(index,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (index,lines) in enumerate(loops)
        A.model.add_curve_loop(collect(lines);tag=index)
        A.model.add_plane_surface([index];tag=index)
        A.mesh.set_transfinite_surface(index)
        A.mesh.set_recombine(2,index)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    return nothing
end

function constructor_count()
    authority=A.CURRENT[].meshing.mesh_authority
    return authority===nothing ? 0 : length(authority.surface_constructor_inputs)
end

function volume_signature()
    families,labels,nodes=A.mesh.get_elements(3,1)
    counts=Dict{Int,Int}()
    for (kind,ids) in zip(families,labels)
        counts[Int(kind)]=length(ids)
    end
    ids,xyz,_=A.mesh.get_nodes()
    points=Dict(uid=>xyz[3index-2:3index] for (index,uid) in enumerate(ids))
    volume=0.0
    for (kind,uids,vertices) in zip(families,labels,nodes)
        width=E.msh_num_nodes(kind)
        cells=reshape(vertices,width,length(uids))
        determinant(nodes)=det(hcat(points[nodes[2]]-points[nodes[1]],
                                   points[nodes[3]]-points[nodes[1]],
                                   points[nodes[4]]-points[nodes[1]]))/6
        for cell in eachcol(cells)
            cell_volume=if kind==4
                determinant(cell)
            elseif kind==7
                # The two signed tetrahedra partition a cyclic Pyramid5 base.
                determinant(cell[[1,2,3,5]])+determinant(cell[[1,3,4,5]])
            else
                error("unexpected volume family $kind")
            end
            cell_volume>0 || error("non-positive $kind cell volume $cell_volume")
            volume+=cell_volume
        end
    end
    return counts,volume
end

@testset "order elevation keeps retained Surface constructor receipts" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        recombined_cube()
        A.mesh.generate(2)
        @test constructor_count()==6
        A.mesh.set_order(2)
        # The order pass keeps every receipt bound through the elevated
        # products' corner skeletons.
        @test constructor_count()==6
        families,_,_=A.mesh.get_elements(2,1)
        @test families==Int32[10]
        A.mesh.generate(3)
        counts,volume=volume_signature()
        @test counts[4]>0 && counts[7]==6
        @test volume≈1.0 atol=2e-10
        mesh=A.mesh.get()
        @test mesh isa E.MixedMesh && Tessella.validate(mesh).ok
    finally
        A.finalize()
    end
end

@testset "recombined Region still meshes directly, de-elevated or remeshed" begin
    for sequence in (:direct,:de_elevate,:surface_remesh)
        A.initialize()
        try
            A.option("General.Terminal",0)
            recombined_cube()
            if sequence===:direct
                A.mesh.generate(3)
            elseif sequence===:de_elevate
                A.mesh.generate(2);A.mesh.set_order(2);A.mesh.set_order(1)
                A.mesh.generate(3)
            else
                A.mesh.generate(2);A.mesh.set_order(2)
                A.mesh.generate(2);A.mesh.generate(3)
            end
            counts,volume=volume_signature()
            @test counts[4]>0 && counts[7]==6
            @test volume≈1.0 atol=2e-10
            @test Tessella.validate(A.mesh.get()).ok
        finally
            A.finalize()
        end
    end
end

@testset "identity retagging keeps retained Surface constructor receipts" begin
    # `model.set_tag` rebuilds the mesh authority; retagging used to rebuild
    # through the nine-argument constructor and silently drop every Surface
    # constructor receipt, so the next `generate(3)` lost the hybrid boundary.
    for retag in (:surface,:point)
        A.initialize()
        try
            A.option("General.Terminal",0)
            recombined_cube()
            A.mesh.generate(2)
            @test constructor_count()==6
            if retag===:surface
                A.model.set_tag(2,1,10)
            else
                A.model.set_tag(0,1,10)
            end
            @test constructor_count()==6
            A.mesh.generate(3)
            counts,volume=volume_signature()
            @test counts[4]>0 && counts[7]==6
            @test volume≈1.0 atol=2e-10
            @test Tessella.validate(A.mesh.get()).ok
        finally
            A.finalize()
        end
    end
end

end
