module OccSweepSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # `occ.extrude` — bounded `MakePrism` sweeps: point→edge,
        # line→ruled face, all-line plane face (with holes)→prism,
        # full-circle disk along its normal→compact cylinder.
        # Topology/measure oracles live in the interface test.
        p=occ.add_point(0.0,0.0,0.0)
        out=occ.extrude([(0,p)],1.0,2.0,3.0)
        @test length(out)==1 && out[1][1]==1
        q=occ.add_point(2.0,0.0,0.0)
        e=occ.add_line(p,q)
        out=occ.extrude([(1,e)],0.0,0.0,4.0)
        @test [d for (d,_) in out]==[1,2,1,1]
        @test occ.get_mass(2,out[2][2])≈8.0
        pts=[occ.add_point(x,y,0.0) for (x,y) in
             [(0,0),(3,0),(3,2),(0,2)]]
        es=[occ.add_line(pts[i],pts[mod1(i+1,4)]) for i in 1:4]
        f=occ.add_plane_surface([occ.add_wire(es)])
        out=occ.extrude([(2,f)],0.0,0.0,5.0)
        @test [d for (d,_) in out]==[2,3,2,2,2,2]
        @test occ.get_mass(3,out[2][2])≈30.0
        @test collect(occ.get_center_of_mass(3,out[2][2]))≈
            collect((1.5,1.0,2.5))
        d=occ.add_disk(10.0,10.0,0.0,2.0,2.0)
        out=occ.extrude([(2,d)],0.0,0.0,3.0)
        @test [x for (x,_) in out]==[2,3,2]
        @test occ.get_mass(3,out[2][2])≈pi*4.0*3.0 rtol=1e-10
        @test_throws ArgumentError occ.extrude([(2,f)],0.0,0.0,0.0)
        @test_throws ArgumentError occ.extrude([(2,999)],0.0,0.0,1.0)
    finally
        api.finalize()
    end
end
end
