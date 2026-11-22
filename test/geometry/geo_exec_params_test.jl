module GeoExecParamsTest

using Test, Tessella

const GX = Tessella.GeoExec

# `execute_geo` reads the values `GeoParams` reports off its own run (one
# pass), so option values computed from the geometry are the executed ones.
@testset "execute_geo reports the run's options, names and fields" begin
    mktempdir() do dir
        path = joinpath(dir, "p.geo")
        write(path, """
            Point(1) = {0, 0, 0};
            Point(2) = {2, 0, 0};
            Line(1) = {1, 2};
            Physical Curve("edge") = {};
            Physical Curve("edge") += {1};
            Mesh.MeshSizeMax = Sqrt(4) / 4;
            Mesh.MeshSizeMin = Mesh.MeshSizeMax / 5;
            Mesh.MeshSizeFactor = 3;
            Mesh.RandomSeed = 7;
            Geometry.Tolerance = 1e-6;
            Field[1] = MathEval;
            Field[1].F = "0.1";
            Background Field = 1;
            """)
        p = cd(() -> GX.execute_geo(path; mesh_dim=0).params, dir)
        @test p.mesh_size_max == 0.5
        @test p.mesh_size_min == 0.1
        @test p.mesh_size_factor == 3.0
        @test p.random_seed == Tessella.IO._gmsh_random_seed(7.0)
        @test p.geometry_tolerance == 1e-6
        @test p.physical_groups == Dict((1, 1) => "edge")
        @test sort!(collect(keys(p.fields))) == [1]
        @test p.background_field == 1
        # `Delete Options` (ReInitOptions) forgets the option writes before it.
        write(path, "Mesh.MeshSizeMax = 4;\nDelete Options;\nMesh.MeshSizeMin = 1;\n")
        q = cd(() -> GX.execute_geo(path; mesh_dim=0).params, dir)
        @test isnan(q.mesh_size_max)
        @test q.mesh_size_min == 1.0
    end
end

end
