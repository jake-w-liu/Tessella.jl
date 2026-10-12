module OccQueriesSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # OCC kernel-side queries — registry reads only; the interface
        # test exercises per-dimension boxes, curved-record bounds,
        # transform visibility, containment filters, loop traversal
        # order and set_size transport.
        v=occ.add_box(0.0,0.0,0.0, 2.0,1.0,3.0)
        p=occ.add_point(10.0,0.0,0.0)
        occ.synchronize()
        @test collect(occ.get_bounding_box(3,v))≈[0,0,0,2,1,3]
        @test collect(occ.get_bounding_box(0,1))≈[0,0,0,0,0,0]
        inside=occ.get_entities_in_bounding_box(-1.,-1.,-1.,3.,2.,4.)
        @test length(inside)==27
        @test (0,p) ∉ inside
        loops,surfaces=occ.get_surface_loops(v)
        @test loops==[1] && length(surfaces[1])==6
        loops2,curves=occ.get_curve_loops(1)
        @test loops2==[1] && length(curves[1])==4
        occ.mesh.set_size([(0,p)],0.25)
        occ.synchronize()
        m=api._with_model(identity)
        @test m.point_size[p]≈0.25
    finally
        api.finalize()
    end
end
end
