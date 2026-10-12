module OccEllipseSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # Generic `_gm_sincos`/`_gm_atan2` evaluation paths only — no
        # native-arithmetic records and no meshing (the interface test
        # exercises frames, pcurves, classification, and generate(2)).
        s=occ.add_point(4.0,0.0,0.0)
        c=occ.add_point(0.0,0.0,0.0)
        mj=occ.add_point(4.0,0.0,0.0)
        e=occ.add_point(0.0,2.0,0.0)
        arc=occ.add_ellipse_arc(s,c,mj,e)
        full=occ.add_ellipse(5.0,0.0,0.0,4.0,2.0)
        disk=occ.add_disk(0.0,0.0,3.0,4.0,2.0)
        occ.synchronize()
        @test (1,arc) in api.model.get_entities(1)
        @test (1,full) in api.model.get_entities(1)
        @test (2,disk) in api.model.get_entities(2)
        @test api.model.get_parametrization_bounds(1,full)==([0.0],[Float64(2π)])
        @test api.model.get_value(1,full,[π/2])≈[5.0,2.0,0.0] atol=1e-12
        @test api.model.get_value(1,arc,[π/2])≈[0.0,2.0,0.0] atol=1e-12
        @test api.model.get_value(2,disk,[0.0,0.0])≈[0.0,0.0,3.0] atol=1e-12
        @test api.model.get_value(2,disk,[4.0,0.0])≈[4.0,0.0,3.0] atol=1e-12
    finally
        api.finalize()
    end
end
end
