@testset "bounded planar, spatial, structured and high-order meshing" begin
    square = mesh_planar([0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0],
                         [(1, 2), (2, 3), (3, 4), (4, 1)])
    @test validate(square).ok && MT.ntris(square) == 2
    @test surface_area(square) ≈ 1.0
    edges, incidence = MT.boundary_edges(square.tris)
    @test length(edges) == 4 && incidence == 2
    @test_throws ArgumentError mesh_planar([0.0, NaN], [0.0, 1.0], [(1, 2)])
    @test_throws ArgumentError mesh_planar([0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
                                         [(1, 2), (2, 3), (3, 1)]; max_area=0)

    spatial = Tessella.Mesh3D.delaunay3d([0.0, 1.0, 0.0, 0.0],
                                         [0.0, 0.0, 1.0, 0.0],
                                         [0.0, 0.0, 0.0, 1.0])
    tetrahedron = Tessella.Mesh3D.to_mesh3(spatial)
    @test validate(tetrahedron).ok && MT.ntets(tetrahedron) == 1
    @test signed_volume(tetrahedron) ≈ 1 / 6
    @test Tessella.Mesh3D.check_consistency3(spatial)[1]

    points, parameters = Tessella.Mesh1D.mesh_segment(
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), ConstantSize(0.5))
    @test first(points) == (0.0, 0.0, 0.0) && last(points) == (1.0, 0.0, 0.0)
    @test issorted(parameters) && all(>(0), diff(parameters))
    @test length(points) == length(parameters) == 3

    sides = ([(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
             [(1.0, 0.0, 0.0), (1.0, 0.5, 0.0), (1.0, 1.0, 0.0)],
             [(1.0, 1.0, 0.0), (0.5, 1.0, 0.0), (0.0, 1.0, 0.0)],
             [(0.0, 1.0, 0.0), (0.0, 0.5, 0.0), (0.0, 0.0, 0.0)])
    patch = mesh_transfinite_patch(sides...; face_tag=7, side_tags=(1, 2, 3, 4))
    @test validate(patch).ok && MT.nnodes(patch) == 9 && MT.ntris(patch) == 8
    @test surface_area(patch) ≈ 1.0 && all(==(7), patch.tri_tag)
    @test MT.node(patch, 5) == (0.5, 0.5, 0.0)
    @test_throws ArgumentError mesh_transfinite_patch(sides...; max_nodes=8)

    source = unit_square()
    refined = refine_uniform(source)
    @test validate(refined).ok && MT.ntris(refined) == 8
    @test surface_area(refined) ≈ surface_area(source)
    @test all(==(7), refined.tri_tag) && MT.ntris(source) == 2
    recombined = recombine_triangles(source)
    @test validate(recombined).ok
    quadrangles = filter(block -> block.msh == 3, recombined.blocks)
    @test length(quadrangles) == 1 && size(only(quadrangles).nodes) == (4, 1)
    @test Set(only(quadrangles).nodes) == Set(Int32[1, 2, 3, 4])

    quadratic = p2_tetmesh(unit_tetrahedron())
    @test size(quadratic.coords, 2) == 10 && size(quadratic.tet10) == (10, 1)
    @test validate(quadratic).ok && p2_min_jacobian(quadratic) > 0
    @test p2_volume(quadratic) ≈ 1 / 6
    @test quadratic.tet_tag == Int32[7]

    # A genuinely non-coplanar four-wire Surface meshes in its ruled
    # parameter domain (no planar projection exists) and every node lifts
    # onto the bilinear patch z = 0.5·(1−x)·y.
    warped = mktemp() do path, io
        write(io, """
            Point(1) = {0, 0, 0, 0.5};
            Point(2) = {1, 0, 0, 0.5};
            Point(3) = {1, 1, 0, 0.5};
            Point(4) = {0, 1, 0.5, 0.5};
            Line(1) = {1, 2};
            Line(2) = {2, 3};
            Line(3) = {3, 4};
            Line(4) = {4, 1};
            Curve Loop(1) = {1, 2, 3, 4};
            Surface(1) = {1};
            """)
        close(io)
        execute_geo(path; mesh_dim=2).mesh
    end
    @test validate(warped).ok && MT.ntris(warped) > 0
    for column in axes(warped.coords, 2)
        x, y, z = warped.coords[1, column], warped.coords[2, column],
                  warped.coords[3, column]
        u = clamp(x, 0.0, 1.0); v = clamp(y, 0.0, 1.0)
        @test abs(z - 0.5 * (1 - u) * v) < 1e-12
    end
end
