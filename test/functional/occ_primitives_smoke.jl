module OccPrimitivesSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # OCC primitive solids — topology publication and query paths only;
        # the interface test exercises records, pcurves, encodings, and
        # generate(3).
        v=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,2.0,1.0)
        vs=occ.add_sphere(0.0,0.0,0.0,2.0)
        occ.synchronize()
        m=api._with_model(identity)
        @test m.volumes==Dict(1=>[1],2=>[2])
        @test length(m.surface_loops)==2
        @test length(m.surfaces)==4
        @test length(m.curves)==6
        @test length(m.points)==4
        @test m.surface_types[1]==:cylinder
        @test m.surface_types[4]==:sphere
        @test haskey(m.cylinders,v)
        @test haskey(m.spheres,vs)
        @test (3,v) in api.model.get_entities(3)
        api.initialize()
        occ=api.model.occ
        vt=occ.add_torus(0.0,0.0,0.0,5.0,1.0)
        occ.synchronize()
        m=api._with_model(identity)
        @test m.volumes[vt]==[1]
        @test m.surface_types[1]==:torus
        @test m.surface_loops[1]==[1]
    finally
        api.finalize()
    end
end
end
