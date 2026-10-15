module OccCornerSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # OCC 2-D corner ops — topology publication only; the interface test
        # exercises analytic centers, tangency continuity, pointTag, reverse,
        # requested tags and failure atomicity.
        p1=occ.add_point(0.0,0.0,0.0); p2=occ.add_point(1.0,0.0,0.0)
        p3=occ.add_point(1.0,1.0,0.0)
        e1=occ.add_line(p1,p2); e2=occ.add_line(p2,p3)
        f=occ.fillet_2d(e1,e2,0.2)
        @test f==3
        q1=occ.add_point(0.0,10.0,0.0); q2=occ.add_point(1.0,10.0,0.0)
        q3=occ.add_point(1.0,11.0,0.0)
        f1=occ.add_line(q1,q2); f2=occ.add_line(q2,q3)
        c=occ.chamfer_2d(f1,f2,0.3,0.25)
        @test c==6
        occ.synchronize()
        m=api._with_model(identity)
        @test m.curve_geometry[f].occ===:circle
        @test collect(m.curve_geometry[f].center)≈[0.8,0.2,0.0]
        @test m.curve_geometry[c].occ===:line
        @test occ.get_mass(1,c)≈hypot(0.3,0.25)
    finally
        api.finalize()
    end
end
end
