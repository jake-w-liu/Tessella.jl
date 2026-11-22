module OptimizeExactOrientationTest
# The volume optimizer orients its output cells by the exact sign of their
# volume. Gmsh's `examples/simple_geo/cone.geo` left, in Volume 2, a cell on
# four almost coplanar boundary vertices (below) whose Float64 determinant
# has the opposite sign of its exact one; orienting by the Float64 value
# inverted it and `generate(3)` stopped with "1 negatively-oriented tets".
using Test, Tessella
import Tessella.Mesh3D as M3

const POINTS = [(1.1480684567129222e-17, 0.18750000000072575, 2.4374999999978226),
                (-0.2500000000009678, 3.0615158845677945e-17, 2.2499999999970965),
                (-0.18750000000072575, 2.2961369134258444e-17, 2.4374999999978226),
                (1.5307579422838973e-17, 0.2500000000009678, 2.2499999999970965)]

@testset "optimizer output orientation is exact" begin
    coords = reduce(hcat, [collect(p) for p in POINTS])
    a, b, c, d = POINTS
    # The Float64 determinant and the exact sign disagree on this cell.
    @test Tessella.MeshTypes.tet_signed_volume(a, b, c, d) < 0
    u = b .- a; v = c .- a; w = d .- a
    @test u[1] * (v[2] * w[3] - w[2] * v[3]) - v[1] * (u[2] * w[3] - w[2] * u[3]) +
          w[1] * (u[2] * v[3] - v[2] * u[3]) > 0
    for order in ([1, 2, 3, 4], [2, 1, 3, 4])
        mesh = Tessella.Mesh(coords; tets=reshape(Int32.(order), 4, 1))
        result = M3._gmsh_optimize_tetrahedra(mesh; threshold=0.0)
        @test result.mesh.tets[:, 1] == Int32[2, 1, 3, 4]
        @test Tessella.validate(result.mesh).ok
    end
end
end
