module ApiRefineClearRegenerateTest

using Test, Tessella

# generate(3); refine(); clear([(3,1)]); generate(3) on an OCC box. The
# refined Surfaces keep no producer row order, so the second volume pass takes
# the surface mesh it finds, as for a loaded mesh. Gmsh 4.15.2 on the same
# sequence: 339 nodes / 1125 tetrahedra, then 2072 / 9000 after refine,
# 1082 nodes and 2160 triangles after clearing the volume, and 1800 nodes /
# 7573 tetrahedra after the second generate(3).
@testset "refine, clear the volume, generate again" begin
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.model.occ.add_box(0, 0, 0, 1, 1, 1; tag=1)
        A.model.occ.synchronize()
        counts() = (length(A.mesh.get_nodes()[1]),
                    Dict(Int(t) => length(e) for (t, e) in
                         zip(A.mesh.get_elements()[1], A.mesh.get_elements()[2])))
        A.mesh.generate(3)
        nodes, cells = counts()
        @test nodes == 339                  # Gmsh 4.15.2: 339
        @test cells[2] == 540               # Gmsh 4.15.2: 540
        A.mesh.refine()
        nodes, cells = counts()
        @test cells[2] == 2160              # Gmsh 4.15.2: 2160
        A.mesh.clear([(3, 1)])
        nodes, cells = counts()
        @test nodes == 1082                 # Gmsh 4.15.2: 1082
        @test cells[2] == 2160              # Gmsh 4.15.2: 2160
        @test !haskey(cells, 4)
        A.mesh.generate(3)
        nodes, cells = counts()
        @test cells[2] == 2160              # the refined surface is the input
        @test 1760 <= nodes <= 1840         # Gmsh 4.15.2: 1800 (Tessella 1795)
        @test 7350 <= cells[4] <= 7800      # Gmsh 4.15.2: 7573 (Tessella 7540)
    finally
        A.finalize()
    end
end

end # module
