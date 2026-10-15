module OccOffsetSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # occ.offsetCurve — planar wire offsets; the interface test pins the
        # analytic trims, arc-join tangency, reflex corners, reversed members
        # and failure atomicity.
        p1=occ.add_point(0.0,0.0,0.0); p2=occ.add_point(1.0,0.0,0.0)
        p3=occ.add_point(1.0,1.0,0.0); p4=occ.add_point(0.0,1.0,0.0)
        es=[occ.add_line(p1,p2),occ.add_line(p2,p3),
            occ.add_line(p3,p4),occ.add_line(p4,p1)]
        lp=occ.add_curve_loop(es)
        inn=occ.offset_curve(lp,0.1)
        @test length(inn)==4
        out=occ.offset_curve(lp,-0.1)
        @test length(out)==8
        occ.synchronize()
        m=api._with_model(identity)
        @test all(m.curve_geometry[t].occ===:line for (_,t) in inn)
        @test collect(m.curve_geometry[inn[1][2]].origin)≈[0.0,0.1,0.0]
        @test m.curve_geometry[out[2][2]].occ===:circle
        @test m.curve_geometry[out[2][2]].r≈0.1
    finally
        api.finalize()
    end
end
end
