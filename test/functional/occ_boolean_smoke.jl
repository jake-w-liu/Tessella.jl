module OccBooleanSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # `occ.fuse`/`cut`/`intersect`/`fragment` — volume-operand
        # Booleans through the native analytic engine: the kernel
        # publishes preserve-numbering tags, per-operand
        # `outDimTagsMap` images, and atomic rejection.
        occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
        occ.add_box(0.5,0.5,0.5,1.0,1.0,1.0)
        occ.synchronize()
        out,map=occ.fuse([(3,1)],[(3,2)])
        @test out==[(3,1)] && map==[[(3,1)],[(3,1)]]
        @test occ.get_entities(3)==[(3,1)]
        @test occ.get_mass(3,1)≈1.875
        @test collect(occ.get_center_of_mass(3,1))≈collect((0.75,0.75,0.75))
        occ.add_box(0.5,0.5,0.5,1.0,1.0,1.0)
        occ.synchronize()
        t2=(t for (_,t) in occ.get_entities(3) if t!=1)|>only
        out,map=occ.cut([(3,1)],[(3,t2)])
        @test length(out)==1 && occ.get_mass(3,out[1][2])≈0.875
        # fresh overlapping pair for the corner intersection
        occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
        occ.add_box(0.5,0.5,0.5,1.0,1.0,1.0)
        occ.synchronize()
        ts=sort!([t for (_,t) in occ.get_entities(3) if t!=out[1][2]])
        out,_=occ.intersect([(3,ts[1])],[(3,ts[2])])
        @test length(out)==1 && occ.get_mass(3,out[1][2])≈0.125
        @test_throws ArgumentError occ.fuse([(3,999)],[(3,out[1][2])])
    finally
        api.finalize()
    end
end
end
