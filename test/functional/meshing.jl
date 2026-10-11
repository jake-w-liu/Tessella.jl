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
end
