# Regression: refine! looked for the first bad triangle by rescanning slots from
# the lowest slot touched since the last scan; recycled freelist slots are low,
# so the scan was effectively quadratic (mesh_planar of a unit square with a hole
# at max_area=1e-5: ~150k triangles in 50 s; 2e-6 hit maxsteps).  Pending slots
# now live in a min-heap.  The selected triangle must still be the first bad
# triangle in slot order.
using Test
using Tessella
const M2 = Tessella.Mesh2D

function _first_bad_full_scan(T, interior, B, area)
    for t in eachindex(T.alive)
        (T.alive[t] && t <= length(interior) && interior[t] && !M2._is_ghost_tri(T, t)) || continue
        M2._needs_refine(T, t, B, area, nothing, nothing, nothing) && return Int32(t)
    end
    return Int32(0)
end

@testset "refine! bad-triangle scan" begin
    xs = [0.0, 1, 1, 0, 0.4, 0.6, 0.6, 0.4]; ys = [0.0, 0, 1, 1, 0.4, 0.4, 0.6, 0.6]
    segs = [(1, 2), (2, 3), (3, 4), (4, 1), (5, 8), (8, 7), (7, 6), (6, 5)]
    area = 2e-4; B = 1 / (2 * sind(25.0))
    for steps in (37, 400, 1500)
        T = M2.constrained_delaunay(xs, ys, segs)
        @test_throws ErrorException M2.refine!(T; max_area=area, maxsteps=steps)
        interior = Vector{Bool}(M2.classify_interior(T))
        M2._reset_bad_pending!(T)
        @test M2._find_bad(T, interior, B, area, nothing, nothing, nothing) ==
              _first_bad_full_scan(T, interior, B, area)
    end
    T = M2.constrained_delaunay(xs, ys, segs)
    interior = M2.refine!(T; max_area=area)
    @test _first_bad_full_scan(T, interior, B, area) == 0
    m = mesh_planar(xs, ys, segs; max_area=4e-4)
    @test validate(m).ok
    @test size(m.tris, 2) > 2_500
end
