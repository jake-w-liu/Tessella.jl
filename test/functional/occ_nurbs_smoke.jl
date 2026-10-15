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
    finally
        api.finalize()
    end
end
end
