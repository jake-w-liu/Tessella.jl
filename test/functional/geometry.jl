@testset "geometry entities, fields and geo parameters" begin
    model = GeoModel()
    @test add_point!(model, 0, 0, 0; tag=1) == 1
    @test add_point!(model, 1, 0, 0; tag=2) == 2
    @test add_point!(model, 0, 1, 0; tag=3) == 3
    @test add_line!(model, 1, 2; tag=1) == 1
    add_line!(model, 2, 3; tag=2)
    add_line!(model, 3, 1; tag=3)
    add_curve_loop!(model, [1, 2, 3]; tag=1)
    add_plane_surface!(model, [1]; tag=1)
    @test Tessella.Model.model_entities(model, 2) == [(2, 1)]
    @test Set(Tessella.Model.model_boundary(model, [(2, 1)], false, false, false)) ==
          Set([(1, 1), (1, 2), (1, 3)])
    @test_throws ArgumentError add_point!(model, 5, 5, 5; tag=1)
    @test_throws ArgumentError add_line!(model, 1, 99; tag=4)
    @test Tessella.Model.model_entities(model, 1) == [(1, 1), (1, 2), (1, 3)]

    curve = NURBSCurve(2, [0, 0, 0, 1, 1, 1],
                       [(0.0, 0.0, 0.0), (1.0, 1.0, 0.0), (2.0, 0.0, 0.0)])
    @test nurbs_eval(curve, 0.5) == (1.0, 0.5, 0.0)
    @test_throws ArgumentError NURBSCurve(2, [0, 1], [(0.0, 0.0, 0.0)])
    distance = DistanceField(points=[(0.0, 0.0, 0.0)])
    @test field_value(distance, 3.0, 4.0, 0.0) == 5.0
    field = MinSize([ConstantSize(0.5), FunctionSize((x, y, z) -> 0.1 + x)])
    @test size_at(field, 0.1, 0.0, 0.0) ≈ 0.2
    @test size_at(field, 1.0, 0.0, 0.0) == 0.5
    @test_throws ArgumentError ConstantSize(-1.0)
    @test_throws ArgumentError size_at(FunctionSize((x, y, z) -> NaN), 0.0, 0.0, 0.0)

    mktempdir() do directory
        path = joinpath(directory, "parameters.geo")
        write(path, """
            Mesh.MeshSizeMin = 0.125;
            Mesh.MeshSizeMax = 1.0;
            Mesh.RandomSeed = 3;
            """)
        parameters = Tessella.IO.read_geo_params(path)
        @test parameters.mesh_size_min == 0.125 && parameters.mesh_size_max == 1.0
        @test parameters.random_seed == 3
        @test_throws ArgumentError Tessella.IO.read_geo_params(path; max_file_bytes=5)
    end
end
