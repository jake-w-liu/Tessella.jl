module OCCRevolveMixedDimensionsChecks
using Test,Tessella
const A=Tessella.API

# OpenCASCADE sweeps each input dimension as its own compound: the arc a
# point sweeps and the lateral arc of a line through that point are two
# curves, each with its own end point. The revolve shared them (the prism
# already kept them apart). Lists and counts are Gmsh 4.15.2's.
@testset "OCC revolve of a point and the line through it" begin
    for (angle,expected,counts) in (
            (π/2,[(0,3),(1,2),(1,5),(2,1),(1,3),(1,4)],(5,5,1,0)),
            (2π,[(0,1),(1,2),(1,1),(2,1),(1,3),(1,4)],(2,4,1,0)))
        A.initialize()
        try
            A.option("General.Terminal",0)
            occ=A.model.occ
            inner=occ.add_point(1,0,0)
            outer=occ.add_point(2,0,0)
            line=occ.add_line(inner,outer)
            @test occ.revolve([(0,inner),(1,line)],0,0,0,0,0,1,angle)==expected
            occ.synchronize()
            @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==counts
            # the point's own arc and the line's lateral arc at that point
            @test occ.get_mass(1,2)≈angle rtol=1e-12
            @test occ.get_mass(1,3)≈angle rtol=1e-12
            @test occ.get_mass(1,4)≈2angle rtol=1e-12
            @test occ.get_mass(2,1)≈angle*1.5 rtol=1e-12
            if angle<2π
                @test A.model.get_boundary([(1,2)],false,true)==[(0,1),(0,3)]
                @test A.model.get_boundary([(1,3)],false,true)==[(0,1),(0,4)]
                @test A.model.get_boundary([(1,4)],false,true)==[(0,2),(0,5)]
                @test A.model.get_boundary([(1,5)],false,true)==[(0,4),(0,5)]
            end
        finally
            A.finalize()
        end
    end
end
end
