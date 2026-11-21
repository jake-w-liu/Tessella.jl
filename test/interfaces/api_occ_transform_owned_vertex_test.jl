module OCCTransformOwnedVertexChecks
using Test,Tessella
const A=Tessella.API

# `occ.translate` of a point that is a corner of a box dragged the corner:
# Point 1 moved, eight points remained, and Curve 1 no longer fitted its end
# points ("geometry is inconsistent with its endpoints" at meshing). Gmsh
# 4.15.2 moves a copy: Point 1 becomes a free vertex at the new place, the
# box is unchanged and its corner is bound again as Point 9, where Curves
# 1, 2 and 10 now end.
@testset "OCC translate of a vertex owned by a solid" begin
    for synchronized_first in (false,true)
        A.initialize()
        try
            A.option("General.Terminal",0)
            occ=A.model.occ
            box=occ.add_box(0,0,0,1,1,1)
            synchronized_first && occ.synchronize()
            occ.translate([(0,1)],0,0,1)
            occ.synchronize()
            @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(9,12,6,1)
            @test collect(A.model.get_value(0,1,Float64[]))==[0.0,0.0,2.0]
            @test collect(A.model.get_value(0,9,Float64[]))==[0.0,0.0,1.0]
            @test isempty(A.model.get_adjacencies(0,1)[1])
            @test sort(A.model.get_adjacencies(0,9)[1])==[1,2,10]
            @test A.model.get_boundary([(1,1)],false,true)==[(0,2),(0,9)]
            @test A.model.get_boundary([(1,2)],false,true)==[(0,9),(0,3)]
            @test A.model.get_boundary([(1,10)],false,true)==[(0,9),(0,5)]
            @test occ.get_mass(3,box)≈1 rtol=1e-12
            @test collect(A.model.get_bounding_box(3,box))≈[0,0,0,1,1,1] atol=1e-6
        finally
            A.finalize()
        end
    end
end
end
