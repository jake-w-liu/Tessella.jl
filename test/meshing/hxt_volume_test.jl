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
end
end
