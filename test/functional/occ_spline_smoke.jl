module OccSplineSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # OCC spline — GeomAPI_Interpolate-mirroring chord-length collocation
        # solve; the interface test exercises degree fallbacks, the independent
        # tridiagonal cross-check, periodic closure, Lagrange tangents, C2
        # continuity, atomicity, queries, transforms and wire meshing.
        p1=occ.add_point(0.0,0.0,0.0)
        p2=occ.add_point(1.0,1.0,0.0)
        p3=occ.add_point(2.0,0.5,0.0)
        p4=occ.add_point(3.0,0.0,0.0)
        s=occ.add_spline([p1,p2,p3,p4])
        @test s isa Integer && s>0
        @test occ.get_mass(1,s)>3.0
        line=occ.add_spline([p1,p4])
        @test occ.get_mass(1,line)≈3.0
        closed=occ.add_spline([p1,p2,p4,p1])
        @test occ.get_mass(1,closed)>4.0
        @test_throws ArgumentError occ.add_spline([p1])
        @test_throws ArgumentError occ.add_spline([p1,99])
        occ.synchronize()
        @test api.model.get_type(1,s)=="BSpline"
    finally
        api.finalize()
    end
end
end
