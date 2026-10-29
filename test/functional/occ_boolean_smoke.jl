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
    api.initialize()
    try
        occ=api.model.occ
        # torus × box — the doubly-periodic torus chart needs both u and
        # v seam unwrapping at the operand PLC (a raw v wrap retraces the
        # wire into an insert_segment! overlap), and the solid-side ray
        # parity needs a line × torus solver or every containment test
        # reports "outside" and kept-region classification drops the
        # section arcs off the plane wires. The box's vertical edges dip
        # inside the tube near its inner wall, so intersect returns the
        # four contact slivers; the kept torus-face wires bound a
        # non-canonical sub-domain, so their volumes integrate through
        # the trimmed divergence path. intersect+cut anchor the mass
        # path by inclusion–exclusion against the analytic torus volume.
        function build()
            t=occ.add_torus(0.0,0.0,0.0,1.0,0.3)
            b=occ.add_box(-0.5,-0.5,-0.5,1.0,1.0,1.0)
            occ.synchronize()
            return t,b
        end
        torus_vol=2*pi^2*1.0*0.09
        masses=Dict{Symbol,Float64}()
        counts=Dict{Symbol,Int}()
        for op in (:cut,:intersect)
            t,b=build()
            out,_=getfield(occ,op)([(3,t)],[(3,b)];
                                  removeObject=false,removeTool=false)
            occ.synchronize()
            counts[op]=count(e->e[1]==3,out)
            masses[op]=sum(occ.get_mass(3,t2) for (d,t2) in out if d==3;
                           init=0.0)
            occ.remove([out...,(3,t),(3,b)],true)
            occ.synchronize()
        end
        @test counts==Dict(:cut=>1,:intersect=>4)
        @test all(>(0.0),values(masses))
        # inclusion–exclusion: the cut remainder plus the four shared
        # slivers must reassemble the torus volume.
        @test masses[:cut]+masses[:intersect]≈torus_vol atol=1e-3
    finally
        api.finalize()
    end
end
end
