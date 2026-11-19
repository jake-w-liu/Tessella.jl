module OCCBooleanClippedCylinderChecks
using Test,Tessella
const A=Tessella.API

# A cylinder through a box, off the box axis: the kept wall is the carrier's
# seam rectangle clipped to z∈[-1,1], so it starts at v≠0.
#   * its area was read off the carrier's whole chart (unclipped height);
#   * the wire-orientation sign of its divergence integral came from an
#     un-unwrapped plane shoelace, negative for such a rectangle, so the wall
#     entered the solid's volume with the wrong sign (the plug reported
#     −0.126 instead of 0.377).
@testset "OCC solid bounded by a clipped cylinder wall" begin
    r=0.245
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        cylinder=occ.add_cylinder(0.03,0.02,-1.6,0,0,3.2,r)
        out,_=occ.intersect([(3,box)],[(3,cylinder)])
        @test length(out)==1
        tag=out[1][2]
        @test occ.get_mass(3,tag)≈π*r^2*2 rtol=1e-10
        @test collect(occ.get_center_of_mass(3,tag))≈[0.03,0.02,0.0] atol=1e-10
        kernel=A.CURRENT[].occ_kernel
        walls=[t for ((d,t),id) in kernel.bindings if d==2 &&
            get(kernel.shapes[id].geometry.geometry,:occ,nothing)===:cylinder]
        @test length(walls)==1
        @test occ.get_mass(2,only(walls))≈2π*r*2 rtol=1e-12
    finally
        A.finalize()
    end
end
end
