module GmshInsertionSheetSizesTest
# Gmsh's insertVerticesInRegion runs setLcs over the triangles of every GFace,
# embedded Surfaces included, before the tetrahedra: a vertex of an embedded
# sheet is a boundary vertex sized by its sheet triangles. Tessella sized it
# by its tetrahedron edges instead, so Gmsh's tutorial t15 (an embedded
# Surface in Volume 1) gave 431 tetrahedra where Gmsh meshes the same boundary
# with 460 (Tessella's 2-D mesh fed to Gmsh 4.15.2 as discrete surfaces); with
# the sheet sizes Tessella gives 461.
using Test, Tessella
import Tessella.Mesh3D as M3

@testset "embedded sheet vertices take their triangle sizes" begin
    corners = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
               (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
    faces = [(1, 3, 2), (1, 4, 3), (5, 6, 7), (5, 7, 8), (1, 2, 6), (1, 6, 5),
             (2, 3, 7), (2, 7, 6), (3, 4, 8), (3, 8, 7), (4, 1, 5), (4, 5, 8)]
    coords = hcat(reduce(hcat, [collect(p) for p in corners]), [0.3, 0.3, 0.3])
    center = Int32(9)
    tets = reduce(hcat, [Int32[f[1], f[2], f[3], center] for f in faces])
    mesh = Tessella.Mesh(coords; tets=tets)
    for column in axes(mesh.tets, 2)    # positive orientation
        a, b, c, d = (Tuple(mesh.coords[:, mesh.tets[k, column]]) for k in 1:4)
        if Tessella.MeshTypes.tet_signed_volume(a, b, c, d) < 0
            mesh.tets[1, column], mesh.tets[2, column] = mesh.tets[2, column], mesh.tets[1, column]
        end
    end
    @test Tessella.validate(mesh).ok
    T, _, exterior = M3._gmsh3_triangulation(mesh)
    plain, _ = M3._gmsh3_initial_sizes(T, mesh, exterior, false, Dict{Int32,Float64}())
    # Interior vertex: its longest tetrahedron edge, to the corner (1,1,1).
    @test plain[center] ≈ sqrt(3 * 0.7^2)
    sheet = Set([(Int32(1), Int32(2), center)])
    sized, _ = M3._gmsh3_initial_sizes(T, mesh, exterior, false, Dict{Int32,Float64}();
                                       sheet_faces=sheet)
    # A sheet vertex: its longest sheet-triangle edge, to the corner (1,0,0).
    @test sized[center] ≈ sqrt(0.7^2 + 2 * 0.3^2)
    @test sized[1:8] == plain[1:8]
end
end
