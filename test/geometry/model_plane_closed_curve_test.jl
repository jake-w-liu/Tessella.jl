module ModelPlaneClosedCurveTests
using Test,Tessella,LinearAlgebra
const M=Tessella.Model

function closed_patch(;scale=1.0)
    model=GeoModel()
    points=[add_point!(model,scale*x,scale*y,2.0;mesh_size=scale) for
            (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
    curve=add_spline!(model,[points;points[1]])
    loop=add_curve_loop!(model,[curve])
    surface=add_plane_surface!(model,[loop])
    M.set_transfinite_curve!(model,curve,9)
    model.meshing.algorithm_2d=3
    return model,surface,curve
end

@testset "Direct closed-Curve Plane uses its constructor singular-vector fit" begin
    for scale in (1.0,1e-5)
        model,surface,curve=closed_patch(;scale)
        # A single closed Curve contributes only two fast-path samples; its
        # declared endpoint completes the same fit native GEO synchronization
        # already uses. No synchronization receipt exists for a direct model.
        @test isempty(model.meshing.native_geo_kernel.plane_frames)
        options=M._model_default_mesh1d_options(model,"closed Plane regression")
        fit=M._gmsh_native_surface_frame(model,surface,options,"closed Plane regression")
        @test abs(fit.normal[3])≈1.0
        @test fit.origin[3]≈2.0
        uv=M.model_parametrization(model,2,surface,[0.5scale,0.5scale,2.0])
        @test M.model_value(model,2,surface,uv)≈[0.5scale,0.5scale,2.0]
        @test M.model_normal(model,surface,uv)≈[0.,0.,1.]

        mesh=mesh_model_surface(model,surface)
        @test validate(mesh).ok && size(mesh.tris,2)>0
        @test all(z->isapprox(z,2.0;atol=8eps(2.0),rtol=0),mesh.coords[3,:])
        incidence=Dict{Tuple{Int32,Int32},Int}()
        for cell in eachcol(mesh.tris),(a,b) in ((1,2),(2,3),(3,1))
            edge=minmax(cell[a],cell[b])
            incidence[edge]=get(incidence,edge,0)+1
        end
        @test all(in((1,2)),values(incidence))
        boundary=[edge for (edge,count) in incidence if count==1]
        @test length(boundary)==8
        # Compare triangle area to the independently measured retained wire
        # polygon, rather than any mesher-specific triangle count or ordering.
        expected=0.0
        for index in 0:7
            p=M._model_curve_point(model,curve,index/8,"closed Plane regression")
            q=M._model_curve_point(model,curve,(index+1)/8,"closed Plane regression")
            expected+=p[1]*q[2]-p[2]*q[1]
        end
        actual=sum(norm(cross(mesh.coords[:,c[2]]-mesh.coords[:,c[1]],
                              mesh.coords[:,c[3]]-mesh.coords[:,c[1]]))/2
                   for c in eachcol(mesh.tris))
        @test actual≈abs(expected)/2 rtol=1e-12
    end
end
end
