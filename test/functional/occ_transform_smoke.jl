module OccTransformSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # OCC transform family — topology publication and tag churn only;
        # the interface test exercises records, atomicity, detach semantics,
        # copy explosion, and generate(3) on a translated cylinder.
        v=occ.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
        p1=occ.add_point(10.0,0.0,0.0); p2=occ.add_point(11.0,0.0,0.0)
        l=occ.add_line(p1,p2)
        occ.synchronize()
        occ.translate([(3,v)],5.0,0.0,0.0)
        occ.synchronize()
        m=api._with_model(identity)
        @test collect(api.model.get_bounding_box(3,v))≈[5.0,0.0,0.0,6.0,1.0,1.0]
        @test sort(collect(keys(m.curves)))==[13;collect(14:25)]
        out=occ.copy([(3,v)])
        @test out==[(3,2)]
        occ.synchronize()
        @test sort(collect(keys(m.volumes)))==[1,2]
        occ.dilate([(3,2)],0.0,0.0,0.0, 2.0,2.0,2.0)
        occ.synchronize()
        @test collect(api.model.get_bounding_box(3,2))≈[10.0,0.0,0.0,12.0,2.0,2.0]
        occ.remove([(3,2)],true)
        occ.synchronize()
        @test sort(collect(keys(m.volumes)))==[1]
        @test length(m.surface_loops)==1
        @test m.curves[l]==(9,10)
    finally
        api.finalize()
    end
end
end
