module GmshMeshAdaptPlaneTests
# Mesh.Algorithm 1 (MeshAdapt) is Gmsh's refineMeshBDS: split, swap, collapse
# and smooth passes on the BDS mesh until every edge is between 0.7 and 1.4
# sizes long. Gmsh 4.15.2 on this unit square: 165 nodes, 288 triangles.
using Test, Tessella
const A = Tessella.API

const _SQUARE = """
Mesh.Algorithm = 1;
Point(1)={0,0,0,0.1};Point(2)={1,0,0,0.1};Point(3)={1,1,0,0.1};Point(4)={0,1,0,0.1};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
"""

# The first five face nodes Gmsh creates (node tag => x, y).
const _GMSH_NODES = Dict(41 => (0.2893905361706618, 0.840444501530364),
                         42 => (0.20463024340034047, 0.8426998206168165),
                         43 => (0.22757289436208716, 0.771199335204294),
                         44 => (0.14986999440936855, 0.8021521622030434),
                         45 => (0.8379665044099315, 0.7080796157376115))

@testset "MeshAdapt on a plane meshes as Gmsh does" begin
    A.initialize()
    try
        A.option("General.Terminal", 0)
        mktempdir() do directory
            path = joinpath(directory, "square.geo"); write(path, _SQUARE)
            A.open(path)
            A.mesh.generate(2)
            nodes = length(A.mesh.get_nodes()[1])
            triangles = sum(length, A.mesh.get_elements(2)[2])
            if Sys.ARCH === :x86_64
                @test (nodes, triangles) == (165, 288)
            else
                @test abs(nodes - 165) <= 5 && abs(triangles - 288) <= 10
            end
            tags, xyz, _ = A.mesh.get_nodes(2, 1)
            # Gmsh projects a smoothed node with an iterative closest-point
            # search; positions agree to about 1e-12, not to the last bit.
            for (k, tag) in enumerate(tags)
                haskey(_GMSH_NODES, Int(tag)) || continue
                x, y = _GMSH_NODES[Int(tag)]
                @test abs(xyz[3k-2] - x) < 1e-9 && abs(xyz[3k-1] - y) < 1e-9
            end
            @test validate(A.mesh.get()).ok
        end
    finally
        A.finalize()
    end
end
end
