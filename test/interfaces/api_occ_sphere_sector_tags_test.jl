# Tags of a partial OCC sphere, against Gmsh 4.15.2. `BRepPrim_Sphere` lists
# the top bound first (Gmsh: point 1 the north pole or top rim, curve 1 the
# top bound, then the meridian at the sweep's end, the bottom bound, the
# meridian at its start, the cap radii top then bottom, the axis; faces
# sphere, top cap, bottom cap, start side, end side). Tessella numbered the
# bottom first (examples/boolean/primitives.geo's spheres 2-5).
module ApiOccSphereSectorTagsTest
using Test, Tessella

const A = Tessella.API

@testset "partial sphere tags (Gmsh 4.15.2)" begin
    m = A.model; occ = m.occ
    bnd(d, t) = [b for (_, b) in m.get_boundary([(d, t)], false, true)]
    A.initialize()
    try
        occ.add_sphere(0, 0, 0, 0.3; angle1=-pi/4, angle2=pi/4, angle3=pi/2)
        occ.synchronize()
        pts = [round.(m.get_value(0, t, Float64[]); digits=4) .+ 0.0 for t in 1:6]
        @test pts == [[0, 0.2121, 0.2121], [0.2121, 0, 0.2121],
                      [0, 0.2121, -0.2121], [0.2121, 0, -0.2121],
                      [0, 0, 0.2121], [0, 0, -0.2121]]
        @test [bnd(1, t) for t in 1:9] == [[2, 1], [3, 1], [4, 3], [4, 2],
            [5, 2], [5, 1], [6, 3], [6, 4], [6, 5]]
        @test [bnd(2, t) for t in 1:5] == [[-1, -4, 3, 2], [1, -6, 5],
            [8, 3, -7], [8, 4, -5, -9], [2, -6, -9, 7]]
    finally
        A.finalize()
    end
    A.initialize()
    try
        occ.add_sphere(0, 0, 0, 0.3; angle3=pi/4); occ.synchronize()
        @test [round.(m.get_value(0, t, Float64[]); digits=4) .+ 0.0 for t in 1:2] ==
            [[0, 0, 0.3], [0, 0, -0.3]]
        @test [m.get_type(1, t) for t in 1:5] ==
            ["Unknown", "Circle", "Unknown", "Circle", "Line"]
        @test [bnd(2, t) for t in 1:3] == [[-1, -4, 3, 2], [-5, 4], [2, -5]]
    finally
        A.finalize()
    end
end

end
