using Test
using Tessella
using Tessella.MeshTypes: mesh_crc

# Unrun isolated proposal: retain the original complete snapshot test as a
# mandatory neighbor. This adds the previously unobserved analytic B-rep
# invariant immediately before the nested Boolean consumes its pcurves.
@testset "OCC volume transform retains unrelated Boolean pcurve ownership" begin
    M=Tessella.Model
    for repeat in 1:3
        model=GeoModel()
        add_box!(model,0,0,0,2,1,1;tag=1)
        add_box!(model,0,0,0,1,1,1;tag=2)
        @test boolean_volumes!(model,:difference,1,2;tag=3)==3
        original=mesh_crc(mesh_model_volume(model,3))
        owned=sort!(unique!([abs(surface) for shell in model.volumes[3]
            for surface in model.surface_loops[shell]]))
        surface_records=deepcopy(Dict(surface=>model.surface_geometry[surface]
            for surface in owned))
        loop_records=deepcopy(Dict(loop=>model.loops[loop] for surface in owned
            for loop in model.surfaces[surface]))
        @test translate_volume!(model,1,(10,0,0))==1
        @test Dict(surface=>model.surface_geometry[surface] for surface in owned)==surface_records
        @test Dict(loop=>model.loops[loop] for surface in owned
            for loop in model.surfaces[surface])==loop_records
        @test mesh_crc(mesh_model_volume(model,3))==original
        for surface in owned, loop in model.surfaces[surface], signed_curve in model.loops[loop]
            curve=abs(signed_curve)
            @test haskey(model.surface_geometry[surface].pcurves,curve)
            geometry=model.curve_geometry[curve]
            parameters=[geometry.t0,(geometry.t0+geometry.t1)/2,geometry.t1]
            uv=M.model_reparametrize_on_surface(model,1,curve,parameters,surface)
            @test M.model_value(model,2,surface,uv)≈M.model_value(model,1,curve,parameters) atol=1e-12 rtol=0
        end
        add_box!(model,1,0,0,.5,1,1;tag=4)
        @test boolean_volumes!(model,:difference,3,4;tag=5)==5
        @test mesh_crc(mesh_model_volume(model,5)).bbox==((1.5,0.,0.),(2.,1.,1.))
    end
end
