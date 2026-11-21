module OCCNURBSPeriodicPatchChecks
using Test, Tessella
const A=Tessella.API
const M=Tessella.Model

# A torus converted to a NURBS patch is closed in both directions, and its
# face is meshed on a sheet anchored where the seam pcurves lie, which can be
# one period away from the knot domain (v in [2pi,4pi] here). The patch
# evaluators extended the last span's polynomial there: still a point of the
# torus, but another one, with other derivatives, so the surface mesh folded
# over itself (at Mesh.MeshSizeMax 0.4: 312 triangles, 105 of them facing
# inward, 9 edges through facets, where the analytic torus has 638 and none)
# and `generate(3)` refused it. A parameter outside the knot domain of a
# closed direction is now taken one period back.
@testset "a closed NURBS patch is periodic outside its knot domain" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        torus=occ.add_torus(0,0,0,2,0.5)
        occ.convert_to_nurbs([(3,torus)])
        cylinder=occ.add_cylinder(5,0,0,0,0,2,0.5)
        occ.convert_to_nurbs([(3,cylinder)])
        occ.synchronize()
        model=A.CURRENT[]
        g=model.surface_geometry[1]
        @test M._model_nurbs_surface(g)
        @test M._model_nurbs_closed(g.patch.surface,1)
        @test M._model_nurbs_closed(g.patch.surface,2)
        close(a,b)=all(isapprox(x,y;atol=1e-11) for (p,q) in zip(a,b) for (x,y) in zip(p,q))
        for (u,v) in ((0.3,1.0),(2.5,3.0),(6.2,0.2),(4.0,5.5))
            inside=M._occ_surface_d2(g,u,v)
            # On the torus, as before, and now the same point and derivatives.
            @test abs(hypot(hypot(inside[1][1],inside[1][2])-2,inside[1][3])-0.5)<1e-12
            for (du,dv) in ((0.0,2π),(2π,0.0),(0.0,-2π),(-2π,2π),(4π,0.0))
                @test close(M._occ_surface_d2(g,u+du,v+dv),inside)
                @test close(M._occ_surface_d1(g,u+du,v+dv),inside[1:3])
                @test close((M._occ_surface_point(g,u+du,v+dv),),inside[1:1])
            end
        end
        # Within rounding of the domain the one-sided evaluation is kept.
        @test M._model_nurbs_closed_parameters(g.patch,2π+1e-12,-1e-12)==(2π+1e-12,-1e-12)
        @test M._model_nurbs_closed_parameters(g.patch,0.5,2π+1.0)[2]≈1.0
        # A cylinder's lateral patch closes around its axis only.
        lateral=first(geometry for (tag,geometry) in model.surface_geometry if tag>1 &&
            M._model_nurbs_surface(geometry) &&
            (M._model_nurbs_closed(geometry.patch.surface,1) ||
             M._model_nurbs_closed(geometry.patch.surface,2)))
        closed=(M._model_nurbs_closed(lateral.patch.surface,1),
                M._model_nurbs_closed(lateral.patch.surface,2))
        @test count(closed)==1
        lo,hi=M._model_nurbs_surface_bounds(lateral)
        open_axis=closed[1] ? 2 : 1
        outside=ntuple(axis->axis==open_axis ? hi[axis]+0.25*(hi[axis]-lo[axis]) :
                             0.5*(lo[axis]+hi[axis]),2)
        @test M._model_nurbs_closed_parameters(lateral.patch,outside...)==outside
        if Base.JLOptions().compile_enabled==1
            # Taking the parameter back costs no heap object: the bytes are
            # those of the boxed result, as for a parameter inside the domain.
            bytes(patch,u,v)=(M._model_nurbs_patch_d2(patch,u,v);
                              @allocated M._model_nurbs_patch_d2(patch,u,v))
            @test bytes(g.patch,0.3,1.0+2π)==bytes(g.patch,0.3,1.0)
        end
    finally
        A.finalize()
    end
end

end
