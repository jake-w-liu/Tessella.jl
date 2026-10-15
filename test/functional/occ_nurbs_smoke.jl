module OccNurbsSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # occ.convertToNURBS — per-entity exact NURBS conversion keeping tags,
        # topology and the public parametrization (reparametrized conics), so
        # face pcurves stay valid; surface/solid conversion lands separately.
        a=occ.add_point(0.0,0.0,0.0); b=occ.add_point(1.0,0.0,0.0)
        e=occ.add_line(a,b)
        c=occ.add_circle(0.0,0.0,0.0,1.0)
        occ.convert_to_nurbs([(1,e),(1,c),(0,a)])
        occ.synchronize()
        m=api._with_model(identity)
        @test m.curve_geometry[e].occ===:copied_curve
        @test m.curve_geometry[c].curve isa Tessella.Model._ModelOCCReparamCurve
        @test collect(Tessella.Model._model_occ_curve_d2(
            m.curve_geometry[c].curve,1.234)[1])≈[cos(1.234),sin(1.234),0.0] atol=1e-13
        @test occ.get_mass(1,c)≈2π atol=1e-10
        # surface half: analytic faces gain a NURBS patch whose external
        # (u,v) parametrization is preserved exactly
        cyl=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,4.0,1.5)
        occ.convert_to_nurbs([(3,cyl)])
        occ.synchronize()
        m=api._with_model(identity)
        lat=only([t for t in keys(m.surface_geometry) if
                  m.surface_geometry[t].occ===:cylinder])
        g=m.surface_geometry[lat]
        @test hasproperty(g,:patch) && hasproperty(g,:reparam)
        p=Tessella.Model._occ_surface_point(g,0.7,2.0)
        @test collect(p)≈[1.5*cos(0.7),1.5*sin(0.7),2.0] atol=1e-13
        @test occ.get_mass(3,cyl)≈π*1.5^2*4 atol=1e-9
    finally
        api.finalize()
    end
end
end
