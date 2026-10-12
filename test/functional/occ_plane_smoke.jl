module OccPlaneSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # Lines-only square wire -> curve loop -> plane face; the center-form
        # arc paths (which need full compilation for their OCCT-replicated
        # arithmetic) are exercised in test/interfaces/api_occ_plane_test.jl.
        p=[occ.add_point(x,y,0.0) for (x,y) in [(0,0),(2,0),(2,1),(0,1)]]
        e=[occ.add_line(p[a],p[b]) for (a,b) in [(1,2),(2,3),(3,4),(4,1)]]
        w=occ.add_wire(e)
        f=occ.add_plane_surface([w])
        r=occ.add_rectangle(0.0,0.0,3.0,1.0,1.0)
        occ.synchronize()
        @test (2,f) in api.model.get_entities(2)
        @test (2,r) in api.model.get_entities(2)
        lo,hi=api.model.get_parametrization_bounds(2,f)
        @test lo≈[-1.0,-0.5] atol=1e-12
        @test hi≈[1.0,0.5] atol=1e-12
        @test api.model.get_value(2,f,[0.0,0.0])≈[1.0,0.5,0.0] atol=1e-12
        @test api.model.get_value(2,r,[0.0,0.0])≈[0.5,0.5,3.0] atol=1e-12
        @test sort!(api.model.get_boundary([(2,f)]))==[(1,1),(1,2),(1,3),(1,4)]
    finally
        api.finalize()
    end
end
end
