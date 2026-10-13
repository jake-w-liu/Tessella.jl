module OccRevolveSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # `occ.revolve` — bounded `MakeRevol` sweeps at angle 2·Pi:
        # point→rim circle, axial line→quadric/plane face, axial-plane
        # line-bounded face→revolved solid (compact cylinder/cone or
        # Pappus `:revolved`), axial-plane disk→torus. Topology/measure
        # oracles live in the interface test.
        p=occ.add_point(3.0,0.0,1.0)
        out=occ.revolve([(0,p)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
        @test out==[(0,p),(1,out[2][2])]
        a=occ.add_point(2.0,0.0,0.0)
        b=occ.add_point(2.0,0.0,3.0)
        l=occ.add_line(a,b)
        out=occ.revolve([(1,l)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
        @test [d for (d,_) in out]==[1,2,1,1]
        @test occ.get_mass(2,out[2][2])≈2pi*2*3
        ps=[occ.add_point(r,0.0,h) for (r,h) in
            [(0,0),(1.5,0),(1.5,2),(0,2)]]
        es=[occ.add_line(ps[i],ps[mod1(i+1,4)]) for i in 1:4]
        f=occ.add_plane_surface([occ.add_wire(es)])
        out=occ.revolve([(2,f)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
        @test [d for (d,_) in out]==[2,3,2,2,2]
        @test occ.get_mass(3,out[2][2])≈pi*1.5^2*2
        @test collect(occ.get_center_of_mass(3,out[2][2]))≈
            collect((0.0,0.0,1.0))
        d=occ.add_disk(5.0,0.0,0.0,1.0,1.0;zAxis=[0.0,1.0,0.0])
        out=occ.revolve([(2,d)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
        @test [x for (x,_) in out]==[2,3,2]
        @test occ.get_mass(3,out[2][2])≈2*pi*pi*5 rtol=1e-10
        @test_throws ArgumentError occ.revolve(
            [(0,p)],0.0,0.0,0.0,0.0,0.0,1.0,pi)
        @test_throws ArgumentError occ.revolve(
            [(2,999)],0.0,0.0,0.0,0.0,0.0,1.0,2pi)
    finally
        api.finalize()
    end
end
end
