module OccWedgeSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # BRepPrimAPI_MakeWedge: xz trapezoid (0,0)→(dx,0)→(ltx,dz)→(0,dz)
        # extruded dy; ltx==0 collapses to the triangular ridge profile.
        # The interface test exercises the graph layout, signed extents,
        # zAxis frames, validation, lifecycle, transforms and generate(3).
        v=occ.add_wedge(0.0,0.0,0.0,2.0,3.0,4.0)
        occ.synchronize()
        m=api._with_model(identity)
        @test (3,v) in api.model.get_entities(3)
        @test length(m.surface_loops)==1
        @test occ.get_mass(3,v)≈12.0
        @test collect(occ.get_center_of_mass(3,v))≈collect((2/3,1.5,4/3))
        api.initialize()
        occ=api.model.occ
        vt=occ.add_wedge(0.0,0.0,0.0,3.0,2.0,4.0;ltx=1.0)
        @test occ.get_mass(3,vt)≈16.0
    finally
        api.finalize()
    end
end
end
