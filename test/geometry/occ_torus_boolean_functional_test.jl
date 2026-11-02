module OccTorusBooleanFunctionalChecks
using Test, Tessella
let api=Tessella.API
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
