module OCCDistanceEdgeEndChecks
using Test,Tessella
const A=Tessella.API

# The nearest point of a box to a point beyond one of its corners is that
# corner. The projection onto each edge has no foot inside the edge, and
# the search that takes over stopped 8.4e-10 short of the edge's end: the
# distance from Box(0.5,0.5,0.5,1,1,1) to (3,3,3) came out 4.8e-10 above
# sqrt(6.75), which is what Gmsh 4.15.2 returns to the last digit.
@testset "OCC distance to a point beyond the end of an edge" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0.5,0.5,0.5,1,1,1)
        far=occ.add_point(3,3,3)
        occ.synchronize()
        distance=occ.get_distance(3,box,0,far)
        @test distance[1]==sqrt(6.75)
        @test distance[2:4]==(1.5,1.5,1.5)
        # Curve 6 runs along y on xmax, zmax and Surface 2 is the xmax face
        for (dim,tag) in ((1,6),(2,2))
            @test occ.get_distance(dim,tag,0,far)[1]==sqrt(6.75)
        end
        # a foot inside the edge is still the foot
        beside=occ.add_point(3,1.0,1.5)
        @test occ.get_distance(1,6,0,beside)[1:4]==(1.5,1.5,1.0,1.5)
    finally
        A.finalize()
    end
end
end
