module OCCBooleanTorusSeamsChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# A torus chart is periodic in both directions and its face is cut along two
# seams. The Boolean engine cut and folded sections, arcs and anchors only at
# the u seam: a torus poking through a box face lost the half of every
# section beyond the v seam, `fuse` then silently returned the box (8.0 for
# 8.1974) and `fragment` threw "Curve[23] must occur once on each of two
# distinct surfaces". The box face x=1 passes through the torus's axis here,
# so the sections are the two meridian circles, exact; half the torus,
# π²·R·r², lies on either side. Entity counts are Gmsh 4.15.2's.
@testset "OCC Booleans across both seams of a torus" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        ring=occ.add_torus(1,0,0,0.5,0.2)
        out,_=occ.fragment([(3,box)],[(3,ring)])
        occ.synchronize()
        half=π^2*0.5*0.2^2
        @test sort([occ.get_mass(3,t) for (_,t) in out])≈[half,half,8-half] rtol=1e-12
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(11,18,11,3)
    finally
        A.finalize()
    end

    # a plane normal to the axis cuts the tube in two parallels, through the
    # axis in two meridians: circles in space, straight lines in the chart
    ring=(occ=:torus,center=(0.0,0.0,0.0),axis=(0.0,0.0,1.0),
          X=(1.0,0.0,0.0),Y=(0.0,1.0,0.0),r1=2.0,r2=0.5)
    flat=(occ=:plane,center=(0.0,0.0,0.3),axis=(0.0,0.0,1.0),
          X=(1.0,0.0,0.0),Y=(0.0,1.0,0.0))
    parallels=M._brep_sec_plane_torus(flat,ring,true,"test")
    @test sort([s.curve.r for s in parallels])≈[1.6,2.4] rtol=1e-14
    @test all(s->s.curve.kind===:circle && haskey(s.pcB,:lin2d),parallels)
    upright=(occ=:plane,center=(0.0,0.0,0.0),axis=(0.0,1.0,0.0),
             X=(1.0,0.0,0.0),Y=(0.0,0.0,1.0))
    meridians=M._brep_sec_plane_torus(upright,ring,true,"test")
    @test [s.curve.r for s in meridians]==[0.5,0.5]
    @test sort([s.curve.c[1] for s in meridians])≈[-2.0,2.0] atol=1e-14
end
end
