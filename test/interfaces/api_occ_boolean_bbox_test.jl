module OCCBooleanBoundingBoxChecks
using Test,Tessella
const A=Tessella.API

# A sphere centred on the box face x=1: the intersection is the hemisphere
# x≤1, bounded by one great-circle disk in that face. The synchronized model
# boxed a shelled volume by its boundary curves only, so the cap lost its
# bulge (x∈[1,1]); the OCC registry boxed the cap by the whole ball. (The
# trimmed sphere face's own box is checked with the sphere-pole Boolean.)
@testset "OCC bounding box of a Boolean result with a sphere cap" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        ball=occ.add_sphere(1,0.05,0.1,0.45)
        out,_=occ.intersect([(3,box)],[(3,ball)])
        occ.synchronize()
        tag=out[1][2]
        expected=[0.55,-0.4,-0.35,1.0,0.5,0.55]
        @test collect(A.model.get_bounding_box(3,tag))≈expected atol=1e-9
        @test collect(occ.get_bounding_box(3,tag))≈expected atol=1e-9
    finally
        A.finalize()
    end
end
end
