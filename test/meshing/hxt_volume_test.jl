module HxtVolumeTest
using Test, Tessella
const A = Tessella.API

# Mesh.Algorithm3D = 10 (HXT): counts are Gmsh 4.15.2's own (Python gmsh
# wheel, one thread), reproduced cell for cell by the port of HXT.
function hxt_counts(build, size)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        build(); A.model.occ.synchronize()
        A.option("Mesh.MeshSizeMin", size); A.option("Mesh.MeshSizeMax", size)
        A.option("Mesh.Algorithm3D", 10)
        A.mesh.generate(3)
        tets = sum(length, A.mesh.get_elements(3)[2]; init=0)
        return length(A.mesh.get_nodes()[1]), tets
    finally
        A.finalize()
    end
end

@testset "HXT volume mesher" begin
    # A box's surface mesh is Gmsh's, so the whole volume is: nodes, tags
    # and connectivity were checked identical (SHA) against the gmsh wheel.
    coarse = hxt_counts(() -> A.model.occ.add_box(0, 0, 0, 1, 1, 1), 0.25)
    fine = hxt_counts(() -> A.model.occ.add_box(0, 0, 0, 1, 1, 1), 0.2)
    if Sys.ARCH === :x86_64
        @test coarse == (145, 399)
        @test fine == (231, 699)
    else
        @test abs(coarse[2] - 399) <= 8
        @test abs(fine[2] - 699) <= 12
    end
    # A thin plate's boundary is not in the Delaunay tetrahedralization of its
    # vertices: Gmsh recovers 14 facets with TetGen/BR (hxt_boundary_recovery),
    # then refines and optimizes.
    plate = hxt_counts(() -> A.model.occ.add_box(0, 0, 0, 1, 1, 0.05), 0.3)
    if Sys.ARCH === :x86_64
        @test plate == (62, 132)
    else
        @test plate[1] == 62 && abs(plate[2] - 132) <= 6
    end
    # A Line embedded in the Volume: HXT constrains its edges (they are in no
    # triangle) and refines around them.
    embedded = hxt_counts(0.25) do
        A.model.occ.add_box(0, 0, 0, 1, 1, 1)
        a = A.model.occ.add_point(0.2, 0.5, 0.5); b = A.model.occ.add_point(0.8, 0.5, 0.5)
        l = A.model.occ.add_line(a, b); A.model.occ.synchronize()
        A.mesh.embed(1, [l], 3, 1)
    end
    if Sys.ARCH === :x86_64
        @test embedded == (146, 414)
    else
        @test embedded[1] >= 140 && abs(embedded[2] - 414) <= 10
    end
end
end
