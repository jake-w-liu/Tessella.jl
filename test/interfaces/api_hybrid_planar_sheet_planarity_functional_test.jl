module HybridSheetPlanarityFunctionalChecks
using Test, Tessella, LinearAlgebra
const A=Tessella.API
const M=Tessella.Model
function cube(;sheet=((.3,.3,.45),(.7,.3,.45),(.7,.7,.45),(.3,.7,.45)))
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
    for (tag,p) in zip(101:104,sheet)
        A.model.add_point(p...;tag,meshSize=1.)
    end
    for (tag,pair) in zip(101:104,((101,102),(102,103),(103,104),(104,101)))
        A.model.add_line(pair...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    A.model.add_curve_loop([101,102,103,104];tag=101)
    A.model.add_plane_surface([101];tag=101)
    A.mesh.set_transfinite_surface(101)
    A.mesh.embed(2,[101],3,1)
end

function actual_surface_input(sheet)
    cube(;sheet)
    A.mesh.generate(2)
    model=A.CURRENT[]
    authority=model.meshing.mesh_authority
    @test authority!==nothing && authority.mesh!==nothing
    product=M._model_mesh_product(authority,2,101;caller="planar sheet regression")
    @test product!==nothing && only(product.blocks).msh==2 && size(only(product.blocks).nodes)==(3,2)
    return model,authority,product
end

@testset "actual rotated two-Tri Surface binds its rounded planar coordinates" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        origin=(.35,.35,.35);u=(.8,.6,0.);v=(-.36,.48,.8);side=.25
        sheet=(origin,ntuple(i->origin[i]+side*u[i],3),
               ntuple(i->origin[i]+side*(u[i]+v[i]),3),ntuple(i->origin[i]+side*v[i],3))
        # The literal Float64 vertices are not exactly dyadically coplanar.
        @test Tessella.Predicates.orient3(sheet...)!=0
        model,authority,product=actual_surface_input(sheet)
        constructor=M._model_native_hybrid_volume_input(model,1,authority,"planar sheet regression")
        @test constructor!==nothing && length(constructor.embedded_surface.triangles)==2
        witness=constructor.embedded_surface
        @test Set(record.coordinate for record in witness.records)==Set(sheet)
        @test Set(record.node_owner for record in witness.records)==Set((0,Int32(tag)) for tag in 101:104)
        @test length(witness.lower_cells)==8 && length(Set(triangle.uid for triangle in witness.triangles))==2
        @test M._model_native_hybrid_surface_source_certificate(witness,authority,
            constructor.node_uids,constructor.coordinates,"planar sheet regression")===nothing
        constraints=M._model_native_hybrid_surface_constraints(constructor,constructor.node_uids,"planar sheet regression")
        @test length(constraints.facets)==2 && length(constraints.edges)==5
        @test M._model_native_hybrid_surface_containment_certificate(constructor,"planar sheet regression")===nothing
    finally
        A.finalize()
    end
end

@testset "actual lower sheet rejects a changed nonplanar corner without publication" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        model,authority,product=actual_surface_input(((.3,.3,.45),(.7,.3,.45),(.7,.7,.45),(.3,.7,.45)))
        A.model.set_coordinates(104,.3,.7,.46)
        A.option("Mesh.MeshOnlyEmpty",1)
        before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
        failure=try
            A.mesh.generate(3)
            nothing
        catch error
            error
        end
        @test failure isa ArgumentError && occursin("not a certified planar sheet",sprint(showerror,failure))
        @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==before
    finally
        A.finalize()
    end
end
end
