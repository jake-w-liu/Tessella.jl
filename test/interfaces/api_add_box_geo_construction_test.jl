module ApiAddBoxGeoConstructionTest
# `model.add_box` builds a box from built-in kernel entities. Its mesh is the
# mesh of the equivalent `.geo` construction below (same Points, Line
# directions, Curve Loop signs and signed Surface Loop), node for node and cell
# for cell, and that construction's counts are Gmsh 4.15.2's:
#   default sizes: 343 nodes, 540 triangles, 1,138 tetrahedra;
#   size 0.25:     149 nodes, 264 triangles,   418 tetrahedra.
using Test, Tessella
const A = Tessella.API

const GEO = """
Point(1) = {0, 0, 1}; Point(2) = {0, 0, 0}; Point(3) = {0, 1, 1}; Point(4) = {0, 1, 0};
Point(5) = {1, 0, 1}; Point(6) = {1, 0, 0}; Point(7) = {1, 1, 1}; Point(8) = {1, 1, 0};
Line(1) = {2, 1}; Line(2) = {1, 3}; Line(3) = {4, 3}; Line(4) = {2, 4};
Line(5) = {6, 5}; Line(6) = {5, 7}; Line(7) = {8, 7}; Line(8) = {6, 8};
Line(9) = {2, 6}; Line(10) = {1, 5}; Line(11) = {4, 8}; Line(12) = {3, 7};
Curve Loop(1) = {4, 3, -2, -1}; Plane Surface(1) = {1};
Curve Loop(2) = {8, 7, -6, -5}; Plane Surface(2) = {2};
Curve Loop(3) = {1, 10, -5, -9}; Plane Surface(3) = {3};
Curve Loop(4) = {3, 12, -7, -11}; Plane Surface(4) = {4};
Curve Loop(5) = {9, 8, -11, -4}; Plane Surface(5) = {5};
Curve Loop(6) = {10, 6, -12, -2}; Plane Surface(6) = {6};
Surface Loop(1) = {-1, 2, -3, 4, -5, 6};
Volume(1) = {1};
"""

count_cells(dim) = sum(length, A.mesh.get_elements(dim)[2]; init=0)

function meshed(build, h; repeats=1)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        build()
        h === nothing || (A.option("Mesh.MeshSizeMax", h); A.option("Mesh.MeshSizeMin", h))
        snapshots = Any[]
        for _ in 1:repeats
            A.mesh.generate(3)
            tags, xyz, _ = A.mesh.get_nodes()
            types, elements, nodes = A.mesh.get_elements()
            push!(snapshots, (copy(tags), copy(xyz), copy(types), deepcopy(elements), deepcopy(nodes),
                              count_cells(2), count_cells(3)))
        end
        return snapshots
    finally
        A.finalize()
    end
end

@testset "model.add_box is the built-in .geo box Gmsh meshes" begin
    mktempdir() do directory
        path = joinpath(directory, "box.geo")
        write(path, GEO)
        for (h, nodes, triangles, tets) in ((nothing, 343, 540, 1138), (0.25, 149, 264, 418))
            geo = only(meshed(() -> A.open(path), h))
            box = meshed(() -> A.model.add_box(0, 0, 0, 1, 1, 1; tag=1), h; repeats=2)
            @test length(geo[1]) == nodes && geo[6] == triangles && geo[7] == tets
            @test box[1] == geo
            # Regenerating the unchanged model gives the same public mesh.
            @test box[2] == box[1]
        end
    end
end
end
