module IoGmshFormatParityTest

using Test, Tessella

# Each expectation below was checked against Gmsh 4.15.2 reading or writing the
# same bytes (see the commit message for the corpus).
const IO = Tessella.IO

function blocks_of(f)
    return [(b.msh, vec(b.nodes)) for b in f.mesh.blocks]
end

@testset "Gmsh-compatible foreign formats" begin
    mktempdir() do dir
        @testset "BDF: CQUAD9 and CPYRAM cards" begin
            grids = join(["GRID    $(i)       0       $(i % 3).0     $(i % 5).0     $(i % 7).0\n"
                          for i in 1:14])
            path = joinpath(dir, "g.bdf")
            # CQUAD9: six grids on the first line, three on the continuation.
            write(path, grids *
                "CQUAD9  1       4       1       2       3       4       5       6       +\n" *
                "+       7       8       9\n" *
                "CPYRAM  2       4       1       2       3       4       5\n" *
                "ENDDATA\n")
            f = IO.read_bdf(path)
            @test blocks_of(f) == [(10, Int32[1:9;]), (7, Int32[1:5;])]
            # The 14-node second-order pyramid Gmsh writes under CPYRAM.
            write(path, grids *
                "CPYRAM  2       4       1       2       3       4       5       6       +\n" *
                "+       7       8       9       10      11      12      13      14\n" *
                "ENDDATA\n")
            @test blocks_of(IO.read_bdf(path)) == [(14, Int32[1:14;])]
        end

        @testset "BDF writer layouts Gmsh can read" begin
            coords = [Float64(i % 3) + 0.1 * i for _ in 1:3, i in 1:6]
            mesh = Tessella.MixedMesh(coords,
                [Tessella.ElementBlock(9, reshape(Int32[1:6;], 6, 1), Int32[0])])
            path = joinpath(dir, "w.bdf")
            IO.write_bdf(path, mesh)
            line = only(filter(l -> startswith(l, "CTRIA6"), readlines(path)))
            # All six grids on the first line and an explicit property id: Gmsh
            # skips empty fields, so a blank PID or the MSC THETA/ZOFFS blanks
            # would shift the node list.
            fields = [strip(line[8k+1:min(end, 8k+8)]) for k in 1:8]
            @test fields == ["1", "0", "1", "2", "3", "4", "5", "6"]
            @test blocks_of(IO.read_bdf(path)) == [(9, Int32[1:6;])]
        end

        @testset "STL soup with a non-manifold edge" begin
            path = joinpath(dir, "fan.stl")
            facet(a, b, c) = "facet normal 0 0 1\nouter loop\n" *
                join(["vertex $(join(v, ' '))\n" for v in (a, b, c)]) * "endloop\nendfacet\n"
            write(path, "solid s\n" *
                facet((0, 0, 0), (1, 0, 0), (0, 1, 0)) *
                facet((0, 0, 0), (1, 0, 0), (0, -1, 0)) *
                facet((0, 0, 0), (1, 0, 0), (0, 0, 1)) * "endsolid s\n")
            mesh = IO.read_stl(path)
            @test size(mesh.coords, 2) == 5 && size(mesh.tris, 2) == 3
        end

        @testset "VTK: Gmsh's 14-node pyramid" begin
            path = joinpath(dir, "p.vtk")
            points = join(["$(i % 3) $(i % 5) $(i % 7)\n" for i in 1:14])
            write(path, "# vtk DataFile Version 2.0\nt\nASCII\nDATASET UNSTRUCTURED_GRID\n" *
                "POINTS 14 double\n" * points *
                "CELLS 1 15\n14 $(join(0:13, ' '))\nCELL_TYPES 1\n14\n")
            f = IO.read_vtk(path)
            @test blocks_of(f) == [(14, Int32[1:14;])]
        end
    end
end

end # module
