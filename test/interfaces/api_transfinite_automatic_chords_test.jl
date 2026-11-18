using Test
using Tessella

# `setTransfiniteAutomatic` follows upstream `MeshSetTransfiniteFacesAutomatic`:
# chord node counts come from the current 1-D size settings (here
# Mesh.MeshSizeMax) and are rounded up to an even count; quad faces become
# `Right` transfinite surfaces, recombined by default, and the six-faced
# volume becomes transfinite. Gmsh 4.15.2 on this box: 320 nodes, 189
# hexahedra, chords of 4/8/10 nodes.
@testset "set_transfinite_automatic matches Gmsh chord counts" begin
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.model.occ.add_box(0, 0, 0, 1, 2, 3)
        A.model.occ.synchronize()
        A.option("Mesh.MeshSizeMax", 0.4)
        A.mesh.set_transfinite_automatic()
        A.mesh.generate(3)
        tags, xyz, _ = A.mesh.get_nodes()
        @test length(tags) == 320
        types, elements, _ = A.mesh.get_elements(3)
        @test types == [5]
        @test length(elements[1]) == 189
        counts = Set{Int}()
        for (_, curve) in A.model.get_entities(1)
            ctags, _, _ = A.mesh.get_nodes(1, curve, true)
            push!(counts, length(ctags))
        end
        @test counts == Set([4, 8, 10])
    finally
        A.finalize()
    end
end
