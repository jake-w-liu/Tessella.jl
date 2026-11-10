module OccNurbsSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # occ.convertToNURBS — per-entity exact NURBS conversion keeping tags
        # and topology while exposing the rational conic parameter.
        a=occ.add_point(0.0,0.0,0.0); b=occ.add_point(1.0,0.0,0.0)
        e=occ.add_line(a,b)
        c=occ.add_circle(0.0,0.0,0.0,1.0)
        occ.convert_to_nurbs([(1,e),(1,c),(0,a)])
        occ.synchronize()
        m=api._with_model(identity)
        @test m.curve_geometry[e].occ===:copied_curve
        @test m.curve_geometry[c].curve isa Tessella.Model._ModelOCCNURBSCurve
        chart=Tessella.Model._model_occ_conic_chart(0.,2π)
        rational=Tessella.Model._occ_conic_reparam(chart,1.234)[1]
        lo,hi=api.model.get_parametrization_bounds(1,c)
        @test lo≈[0.] && hi≈[2π]
        @test collect(api.model.get_value(1,c,[rational]))≈
            [cos(1.234),sin(1.234),0.0] atol=1e-13
        @test collect(Tessella.Model._model_occ_curve_d2(
            m.curve_geometry[c].curve,rational)[1])≈
            [cos(1.234),sin(1.234),0.0] atol=1e-13
        @test occ.get_mass(1,c)≈2π atol=1e-10
        # Surface conversion exposes the matching rational UV chart.
        cyl=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,4.0,1.5)
        occ.convert_to_nurbs([(3,cyl)])
        occ.synchronize()
        m=api._with_model(identity)
        lat=only([t for t in keys(m.surface_geometry) if
                  m.surface_geometry[t].occ===:cylinder])
        g=m.surface_geometry[lat]
        @test hasproperty(g,:patch) && hasproperty(g,:analytic_chart)
        u=Tessella.Model._occ_conic_reparam(g.analytic_chart.u,0.7)[1]
        p=Tessella.Model._occ_surface_point(g,u,2.0)
        @test collect(p)≈[1.5*cos(0.7),1.5*sin(0.7),2.0] atol=1e-13
        ulo,uhi=api.model.get_parametrization_bounds(2,lat)
        @test ulo≈[0.,0.] && uhi≈[2π,4.]
        @test collect(api.model.get_value(2,lat,[u,2.0]))≈
            [1.5*cos(0.7),1.5*sin(0.7),2.0] atol=1e-13
        @test occ.get_mass(3,cyl)≈π*1.5^2*4 atol=1e-9
    finally
        api.finalize()
    end
end
end
