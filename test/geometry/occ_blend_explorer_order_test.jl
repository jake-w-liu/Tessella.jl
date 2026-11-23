# Tags of a filleted or chamfered box, against Gmsh 4.15.2.
#
# Gmsh's `_multiBind` numbers the result's faces, curves and points in the
# order OpenCASCADE's `TopExp_Explorer` meets them, which is the order the
# fillet builder walks the rebuilt shell (breadth first from the first cut
# face, each face's edges breadth first through their vertices). Tessella
# numbered them in its own construction order: examples/boolean/chamfer.geo
# came out with the right shapes under other tags, and fillet_chamfer.geo
# filleted another edge than Gmsh (its `Fillet{1}{1}` names curve 1 of the
# chamfered box). Points (by position), curves (by end points) and surfaces
# (by curves) are Gmsh's, tag for tag.
module OccBlendExplorerOrderTest
using Test, Tessella

const A = Tessella.API

function tags_of()
    m = A.model
    P = Dict(t => round.(m.get_value(0, t, Float64[]); digits=4) .+ 0.0
             for (_, t) in m.get_entities(0))
    C = Dict(t => Set(abs(b) for (_, b) in m.get_boundary([(1, t)], false, false))
             for (_, t) in m.get_entities(1))
    S = Dict(t => Set(abs(b) for (_, b) in m.get_boundary([(2, t)], false, false))
             for (_, t) in m.get_entities(2))
    return P, C, S
end

function check(P, C, S, pts, crv, srf)
    @test length(P) == length(pts)
    @test all(i -> isapprox(P[i], collect(pts[i]); atol=1e-9), eachindex(pts))
    @test C == Dict(i => Set(c) for (i, c) in enumerate(crv))
    @test S == Dict(i => Set(s) for (i, s) in enumerate(srf))
end

@testset "fillet/chamfer result tags (Gmsh 4.15.2)" begin
    occ = A.model.occ
    # examples/boolean/chamfer.geo: two edges meeting at a corner
    chamfer_pts = ((0, 1, 0.95), (0, 1, 0), (0, 0.2, 0.95), (0, 0.2, 0),
        (0.05, 1, 1), (1, 1, 1), (1, 1, 0), (0.05, 0.15, 1), (0.2, 0, 0),
        (1, 0, 0), (0.2, 0, 1), (1, 0, 1))
    chamfer_crv = ((2, 1), (3, 1), (4, 2), (4, 3), (5, 6), (1, 5), (7, 6),
        (2, 7), (3, 8), (8, 5), (4, 9), (9, 10), (10, 7), (8, 11), (9, 11),
        (12, 6), (11, 12), (10, 12))
    chamfer_srf = ((1, 2, 4, 3), (5, 7, 8, 1, 6), (9, 2, 6, 10),
        (3, 11, 12, 13, 8), (9, 14, 15, 11, 4), (5, 10, 14, 17, 16),
        (18, 13, 7, 16), (12, 15, 17, 18))
    A.initialize()
    try
        occ.add_box(0, 0, 0, 1, 1, 1)
        occ.chamfer([1], [1, 2], [1, 1], [0.2, 0.05]); occ.synchronize()
        check(tags_of()..., chamfer_pts, chamfer_crv, chamfer_srf)
        m = A.model
        # the chamfer faces' end curves are trimmed lines: at least one
        # inner node each (`MinCurveNodes`), as in Gmsh
        @test [t for (_, t) in m.get_entities(1)
               if m.get_type(1, t) == "TrimmedCurve"] == [6, 9, 11, 14]
        # edges run Gmsh's way, faces keep the box's frames and orientation
        @test [b for (_, b) in m.get_boundary([(1, 1)], false, true)] == [2, 1]
        @test [b for (_, b) in m.get_boundary([(3, 1)], false, true)] ==
            [-1, 2, -3, -4, -5, 6, 7, -8]
        @test [b for (_, b) in m.get_boundary([(2, 3)], false, true)] ==
            [-9, 2, 6, -10]
        @test m.get_value(2, 3, [0.0, 0.0]) ≈ [0.025, 0, 0.975] atol=1e-12
        @test m.get_value(2, 3, [1.0, 0.0]) ≈ [0.025, 1, 0.975] atol=1e-12
        A.option("Mesh.MeshSizeMin", 0.1); A.option("Mesh.MeshSizeMax", 0.1)
        A.mesh.generate(1)
        @test [length(A.mesh.get_nodes(1, t, false, false)[1])
               for t in (6, 9, 11, 14)] == [1, 1, 2, 2]
    finally
        A.finalize()
    end
    # the same through a .geo Box, bound into the kernel by the Chamfer
    path = tempname() * ".geo"
    write(path, """
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0,1,1,1};
        Chamfer{1}{1,2}{1,1}{0.2,0.05}
        """)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        check(tags_of()..., chamfer_pts, chamfer_crv, chamfer_srf)
    finally
        A.finalize(); rm(path; force=true)
    end
    # two parallel edges of one face, filleted
    A.initialize()
    try
        occ.add_box(0, 0, 0, 1, 1, 1)
        occ.fillet([1], [1, 3], [0.2]); occ.synchronize()
        check(tags_of()...,
            ((0, 0.2, 1), (0, 0.8, 1), (0, 0.2, 0), (0, 0.8, 0), (0.2, 0, 1),
             (0.2, 1, 1), (1, 0, 1), (1, 1, 1), (0.2, 0, 0), (0.2, 1, 0),
             (1, 0, 0), (1, 1, 0)),
            ((1, 2), (3, 1), (4, 2), (3, 4), (1, 5), (2, 6), (5, 7), (6, 8),
             (7, 8), (3, 9), (9, 5), (4, 10), (10, 6), (9, 11), (10, 12),
             (11, 12), (11, 7), (12, 8)),
            ((1, 2, 4, 3), (1, 5, 7, 9, 8, 6), (10, 11, 5, 2), (12, 13, 6, 3),
             (4, 10, 14, 16, 15, 12), (14, 11, 7, 17), (15, 13, 8, 18),
             (17, 16, 18, 9)))
        @test occ.get_mass(3, 1) ≈ 1 - 2 * (0.04 - pi * 0.01) rtol=1e-12
    finally
        A.finalize()
    end
end

end
