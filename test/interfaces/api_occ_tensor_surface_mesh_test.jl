module APIOCCTensorSurfaceMeshChecks
using Test, Tessella
const A=Tessella.API
const M=Tessella.Model

# `occ.add_bspline_surface` and `occ.add_bezier_surface` build an untrimmed
# patch whose four edges are iso-parameter lines. The face record carried no
# pcurves for them, so the mesher did not recognise a curved OCC face, sent it
# to the plane mesher and `generate(2)` threw "plane geometry requires a
# Plane surface; Surface[1] is a BSpline surface". With the pcurves on the
# record the patch is meshed in its own (u,v) square like every other curved
# face (B-spline at size 0.25: 177 triangles, Gmsh 4.15.2: 177).
#
# Meshing a NURBS patch costs over a minute under the suite's interpreted
# flags whatever the size, so this checks what the mesher needs from the
# record instead: the face is admitted as a curved OCC face, and every
# boundary Curve is the lift of its pcurve.
@testset "untrimmed B-spline and Bezier surfaces carry their pcurves" begin
    # A 4x4 control grid over [0,2]^2 with height `height(x,y)`.
    function patch(add,height)
        tags=Int[]
        for j in 0:3,i in 0:3
            x=2i/3;y=2j/3
            push!(tags,A.model.occ.add_point(x,y,height(x,y)))
        end
        return add(tags,4)
    end
    for (add,height) in (
            (A.model.occ.add_bspline_surface,(x,y)->0.4*sin(2x)*cos(1.5y)),
            (A.model.occ.add_bezier_surface,(x,y)->0.5*(x-1)^2-0.3*(y-1)^2))
        A.initialize()
        try
            surface=patch(add,height);A.model.occ.synchronize()
            m=A.CURRENT[]
            geometry=m.surface_geometry[surface]
            @test M._model_occ_curved_surface(m,surface)
            @test sort!(collect(keys(geometry.pcurves)))==sort(geometry.boundary_curves)
            for curve in geometry.boundary_curves
                origin,direction=geometry.pcurves[curve].fwd.lin2d
                edge=m.curve_geometry[curve]
                for t in (edge.t0,0.3*edge.t0+0.7*edge.t1,edge.t1)
                    uv=origin.+t.*direction
                    lifted=M._occ_surface_point(geometry,uv[1],uv[2])
                    @test maximum(abs,lifted.-M._model_curve_point(m,curve,
                        Float64(t),"tensor pcurve"))<=1e-13
                end
                # The pcurve ends on the parameters of the edge's own vertices.
                first_vertex,last_vertex=m.curves[curve]
                @test origin.+edge.t0.*direction==geometry.vertex_uv[first_vertex]
                @test origin.+edge.t1.*direction==geometry.vertex_uv[last_vertex]
            end
            # The mesher asks for the point, or the point and first
            # derivatives, far more often than for second derivatives. Those
            # entries read one and three of the six nets; they must return
            # the very numbers of the full evaluation. (Meshing a filling
            # surface took 7 to 10 times Gmsh's time while every point paid
            # for all six nets through boxed intermediates.)
            nurbs=geometry.patch
            for (u,v) in ((0.0,0.0),(0.3,0.7),(1.0,0.45),(0.999,1.0))
                full=M._model_nurbs_surface_d2(geometry,u,v)
                @test M._occ_surface_point(geometry,u,v)===full[1]
                @test M._occ_surface_d1(geometry,u,v)===(full[1],full[2],full[3])
                @test M._occ_surface_d2(geometry,u,v)===full
                @test M._model_nurbs_homogeneous(nurbs,1,u,v)[4]==1.0
            end
            # Compiled, an evaluation owns two small buffers (the suite runs
            # interpreted, where every call allocates).
            if Base.JLOptions().compile_enabled==1
                M._model_nurbs_patch_point(nurbs,0.3,0.7)
                M._model_nurbs_patch_d2(nurbs,0.3,0.7)
                @test (@allocated M._model_nurbs_patch_point(nurbs,0.3,0.7))<=512
                @test (@allocated M._model_nurbs_patch_d2(nurbs,0.3,0.7))<=512
            end
        finally
            A.finalize()
        end
    end
end
end
