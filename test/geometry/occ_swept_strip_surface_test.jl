module OCCSweptStripSurfaceTests
# A rectangle and a disk swept along a helical spline, as in Gmsh's
# examples/boolean/pipe.geo (coarser sizes here). The rectangle's lateral
# faces are long BSpline strips; their parameter-plane triangulation fanned
# across the strip and the lifted cells, about one strip width long, crossed
# the swept surface, so the volume could not be recovered. Gmsh makes the
# initial mesh Delaunay in the first fundamental form first (delaunayizeBDS).
using Test, Tessella
const A = Tessella.API

const _PIPE = """
SetFactory("OpenCASCADE");
Mesh.MeshSizeMin = 0.2;
Mesh.MeshSizeMax = 0.2;
npts = 20; r = 1; rd = 0.2; h = 1;
For i In {0:npts-1}
  theta = i * 2*Pi/npts;
  Point(i + 1) = {r * Cos(theta), r * Sin(theta), i * h/npts};
EndFor
Spline(1) = {1:npts};
Wire(1) = {1};
Disk(1) = {1,0,0, rd};
Rectangle(2) = {1+2*rd,-rd,0, 2*rd,2*rd,rd/5};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi/2} { Surface{1,2}; }
Extrude { Surface{1,2}; } Using Wire {1}
Delete{ Surface{1,2}; }
"""

_count(dim, tag) = (blocks = A.mesh.get_elements(dim, tag)[2]; isempty(blocks) ? 0 : sum(length, blocks))
_count(dim) = (blocks = A.mesh.get_elements(dim)[2]; isempty(blocks) ? 0 : sum(length, blocks))

@testset "The lateral strips of a swept rectangle stay on their surface" begin
    A.initialize()
    try
        A.option("General.Terminal", 0)
        mktempdir() do directory
            path = joinpath(directory, "pipe.geo"); write(path, _PIPE)
            A.open(path)
            A.mesh.generate(2)
            # Gmsh 4.15.2: 193, 202, 170 and 162 triangles on the four strips of
            # the swept rectangle (its surfaces 7, 9, 11, 13).
            strips = [_count(2, tag) for tag in (6, 8, 10, 12)]
            @test all(abs.(strips .- [193, 202, 170, 162]) .<= 6)
            # Before, the strips' cells crossed and the volume stopped here.
            A.mesh.generate(3)
            @test 2000 <= _count(3) <= 2300        # Gmsh: 2120 tetrahedra
            @test validate(A.mesh.get()).ok
        end
    finally
        A.finalize()
    end
end
end
