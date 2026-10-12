module OccMeasuresSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # OCC kernel-side BRepGProp measures — registry reads only;
        # the interface test exercises per-dim curves, conic-holed
        # planes, quadric faces, primitive and divergence-path solids,
        # the partial torus, degenerate edges, empty entities and the
        # conic-faced get_mass/get_center_of_mass split.
        v=occ.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
        @test occ.get_mass(3,v)≈1.0
        @test collect(occ.get_center_of_mass(3,v))≈[0.5,0.5,0.5]
        c=occ.add_cylinder(0.0,0.0,0.0, 0.0,0.0,2.0, 1.0)
        @test occ.get_mass(3,c)≈2.0π atol=1e-12
        @test collect(occ.get_center_of_mass(3,c))≈[0.0,0.0,1.0] atol=1e-12
        p1=occ.add_point(0.0,0.0,0.0); p2=occ.add_point(3.0,4.0,0.0)
        l=occ.add_line(p1,p2)
        @test occ.get_mass(1,l)≈5.0
        @test collect(occ.get_center_of_mass(1,l))≈[1.5,2.0,0.0]
        occ.synchronize()
        @test occ.get_mass(3,v)≈1.0
    finally
        api.finalize()
    end
end
end
