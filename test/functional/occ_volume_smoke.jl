module OccVolumeSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # OCC shells/solids — topology publication and query paths only; the
        # interface test exercises wiring, cavity order, and generate(3).
        v=occ.add_box(0.0,0.0,0.0,1.0,2.0,3.0)
        occ.synchronize()
        m=api._with_model(identity)
        @test m.volumes==Dict(1=>[1])
        @test length(m.surface_loops)==1
        @test length(m.surfaces)==6
        @test length(m.curves)==12
        @test length(m.points)==8
        @test (3,v) in api.model.get_entities(3)
        api.initialize()
        occ=api.model.occ
        f=occ.add_rectangle(0.0,0.0,0.0,1.0,1.0)
        l=occ.add_surface_loop([f])
        vol=occ.add_volume([l])
        occ.synchronize()
        m=api._with_model(identity)
        @test m.volumes[vol]==[l]
        @test m.surface_loops[l]==[f]
    finally
        api.finalize()
    end
end
end
