module OccCircleSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        s=occ.add_point(1.0,0.0,0.0)
        m=occ.add_point(0.0,1.0,0.0)
        e=occ.add_point(-1.0,0.0,0.0)
        arc=occ.add_circle_arc(s,m,e;center=false)
        full=occ.add_circle(5.0,0.0,0.0,2.0)
        occ.synchronize()
        # The through-point record evaluates through the generic `_gm_sincos`
        # path (the center form's OCCT-replicated record is exercised in
        # test/interfaces/api_occ_circle_test.jl).
        @test api.model.get_parametrization_bounds(1,arc)==([0.0],[Float64(π)])
        @test api.model.get_value(1,arc,[π/2])≈[0.0,1.0,0.0] atol=1e-12
        @test api.model.get_value(1,full,[0.0])≈[7.0,0.0,0.0] atol=1e-12
        @test (1,arc) in api.model.get_entities(1) && (1,full) in api.model.get_entities(1)
    finally
        api.finalize()
    end
end
end
