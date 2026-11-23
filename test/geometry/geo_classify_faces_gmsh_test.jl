module GeoClassifyFacesGmshTest

using Test, Tessella

const API = Tessella.API

# `classifySurfaces` follows Gmsh 4.15.2's classifyFaces step by step: the
# surfaces, curves and points get Gmsh's tags (curves above the temporary
# ones the classification builds and drops), the points and curves get
# MPoint/MLine elements numbered after every temporary line Gmsh builds, and
# each entity lists its nodes by number and its elements in Gmsh's order;
# a surface's curves form Gmsh's curve loop (GEdgeLoop), signs included.
# Every expected value below was read from the Gmsh 4.15.2 Python wheel.

function entities()
    result = Any[]
    for (d, t) in API.model.get_entities()
        _, et, en = API.mesh.get_elements(d, t)
        bnd = d > 0 ? [b[2] for b in API.model.get_boundary([(d, t)], false, true)] : Int[]
        push!(result, (d, t, Int.(API.mesh.get_nodes(d, t)[1]), Int[x for e in et for x in e],
                       Int[x for e in en for x in e], Int.(bnd)))
    end
    return result
end

@testset "a cube's faces, edges and corners are Gmsh's" begin
    P = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    tags = [3, 7, 11, 2, 20, 5, 13, 9]
    conn = Int[]
    for (a, b, c, d) in ((1, 4, 3, 2), (5, 6, 7, 8), (1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8))
        append!(conn, (tags[a], tags[b], tags[c], tags[a], tags[c], tags[d]))
    end
    API.initialize()
    try
        API.option("General.Terminal", 0)
        API.model.add_discrete_entity(2, 1)
        API.mesh.add_nodes(2, 1, tags, Float64[c for p in P for c in p])
        API.mesh.add_elements(2, 1, [2], [collect(10:21)], [conn])
        API.mesh.classify_surfaces(40 * pi / 180, true, false, pi)
        @test entities() == [
            (0, 1, [2], [88], [2], Int[]), (0, 2, [3], [89], [3], Int[]),
            (0, 3, [9], [91], [9], Int[]), (0, 4, [11], [93], [11], Int[]),
            (0, 5, [7], [95], [7], Int[]), (0, 6, [20], [97], [20], Int[]),
            (0, 7, [5], [99], [5], Int[]), (0, 8, [13], [101], [13], Int[]),
            (1, 14, Int[], [90], [2, 3], [1, 2]), (1, 15, Int[], [92], [2, 9], [1, 3]),
            (1, 16, Int[], [94], [2, 11], [1, 4]), (1, 17, Int[], [96], [3, 7], [2, 5]),
            (1, 18, Int[], [98], [3, 20], [2, 6]), (1, 19, Int[], [100], [5, 7], [7, 5]),
            (1, 20, Int[], [102], [5, 13], [7, 8]), (1, 21, Int[], [103], [5, 20], [7, 6]),
            (1, 22, Int[], [104], [7, 11], [5, 4]), (1, 23, Int[], [105], [9, 13], [3, 8]),
            (1, 24, Int[], [106], [9, 20], [3, 6]), (1, 25, Int[], [107], [11, 13], [4, 8]),
            (2, 2, Int[], [10, 11], [3, 2, 11, 3, 11, 7], [14, 17, 22, -16]),
            (2, 3, Int[], [12, 13], [20, 5, 13, 20, 13, 9], [20, -23, 24, -21]),
            (2, 4, Int[], [14, 15], [3, 7, 5, 3, 5, 20], [17, -19, 21, -18]),
            (2, 5, Int[], [16, 17], [7, 11, 13, 7, 13, 5], [19, 22, 25, -20]),
            (2, 6, Int[], [18, 19], [11, 2, 9, 11, 9, 13], [15, 23, -25, -16]),
            (2, 7, Int[], [20, 21], [2, 3, 20, 2, 20, 9], [14, 18, -24, -15])]
    finally
        API.finalize()
    end
end

@testset "a split loop and a bent strip classify as Gmsh's" begin
    API.initialize()
    try
        API.option("General.Terminal", 0)
        API.model.add_discrete_entity(2, 1)
        API.mesh.add_nodes(2, 1, [1, 2, 3, 4], [0.0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0])
        API.mesh.add_elements(2, 1, [2], [[1, 2]], [[1, 2, 3, 1, 3, 4]])
        API.mesh.classify_surfaces(pi / 4, true, false, pi / 4)
        # Gmsh closes the last piece of the split loop with its first node
        # again: curve 6 holds the lines 4-1 and 1-1.
        @test [(e[1], e[2], e[4], e[5]) for e in entities()] == [
            (0, 1, [22], [1]), (0, 2, [23], [2]), (0, 3, [25], [3]), (0, 4, [27], [4]),
            (1, 3, [24], [1, 2]), (1, 4, [26], [2, 3]), (1, 5, [28], [3, 4]),
            (1, 6, [29, 30], [4, 1, 1, 1]), (2, 2, [1, 2], [1, 2, 3, 1, 3, 4])]
    finally
        API.finalize()
    end
    API.initialize()
    try
        API.option("General.Terminal", 0)
        API.model.add_discrete_entity(2, 1)
        h = [0.0, 0.0, tand(10.0), tand(10.0) + tand(20.0)]
        API.mesh.add_nodes(2, 1, collect(1:8), Float64[c for x in 0:3 for y in (0.0, 1.0) for c in (x, y, h[x+1])])
        conn = Int[]
        for k in 1:3
            a = 2k - 1
            append!(conn, (a, a + 2, a + 3, a, a + 3, a + 1))
        end
        API.mesh.add_elements(2, 1, [2], [collect(1:6)], [conn])
        API.mesh.classify_surfaces(15pi / 180, false, false, pi)
        # Boundary edges still bound the patch (computeNonManifoldEdges adds
        # them whatever includeBoundary says); the triangles are collected
        # depth first from the lowest-numbered one.
        @test [(e[1], e[2], e[4]) for e in entities()] ==
              [(0, 1, [46]), (1, 3, collect(47:54)), (2, 2, [1, 2, 4, 3, 6, 5])]
    finally
        API.finalize()
    end
end

# After CreateGeometry a discrete curve is remeshed along its geometry with
# no old polyline vertex kept: a unit-square plate of 4x4 cells, classified
# with its corners split and meshed at MeshSizeMax 0.4, gives 10 lines per
# side in Gmsh 4.15.2 (keeping the old vertices gave 12).
@testset "created curves are remeshed as Gmsh's discreteEdge" begin
    API.initialize()
    try
        API.option("General.Terminal", 0)
        API.model.add_discrete_entity(2, 1)
        n = 5
        API.mesh.add_nodes(2, 1, [1 + i + n * j for j in 0:n-1 for i in 0:n-1],
            Float64[c for j in 0:n-1 for i in 0:n-1 for c in (i / (n - 1), j / (n - 1), 0.0)])
        conn = Int[]
        for j in 0:n-2, i in 0:n-2
            a = 1 + i + n * j
            append!(conn, (a, a + 1, a + n + 1, a, a + n + 1, a + n))
        end
        API.mesh.add_elements(2, 1, [2], [collect(1:length(conn) ÷ 3)], [conn])
        API.mesh.classify_surfaces(pi / 4, true, false, pi / 4)
        API.mesh.create_geometry()
        API.option("Mesh.MeshSizeMax", 0.4)
        API.mesh.generate(1)
        lines(c) = length(API.mesh.get_elements(1, c)[2][1])
        @test [lines(c) for c in 3:5] == [10, 10, 10]
        # Curve 6 carries Gmsh's degenerate 1-1 line: its zero-length
        # segment adds nothing to meshGEdge's integral, so 10 lines.
        @test lines(6) == 10
    finally
        API.finalize()
    end
end

end
