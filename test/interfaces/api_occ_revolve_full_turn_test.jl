module ApiOccRevolveFullTurnTest

using Test, Tessella

# `occ.revolve(..., 2pi; numElements, recombine)`: OCC_Internals drops the
# Layers attributes of a full turn ("Extruded meshes by revolution only for
# angle < 2*Pi"), so the solid meshes as an ordinary one. Tessella used to
# sweep it anyway and left 12 unwelded seam nodes. Gmsh 4.15.2: 169 nodes,
# 4 points, 84 lines, 314 triangles, 12 quadrangles (the transfinite
# source), 410 tetrahedra, no coincident nodes.
@testset "full-turn OCC revolve meshes without its Layers" begin
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.model.occ.add_rectangle(1, 0, 0, 1, 0.5)
        A.model.occ.rotate([(2, 1)], 0, 0, 0, 1, 0, 0, pi / 2)
        A.model.occ.revolve([(2, 1)], 0, 0, 0, 0, 0, 1, 2pi;
                            numElements=[12], recombine=true)
        A.model.occ.synchronize()
        for (curve, count) in ((1, 5), (2, 4), (3, 5), (4, 4))
            A.mesh.set_transfinite_curve(curve, count)
        end
        A.mesh.set_transfinite_surface(1)
        A.mesh.set_recombine(2, 1)
        A.mesh.generate(3)
        tags, xyz, _ = A.mesh.get_nodes()
        @test length(tags) == 169
        counts = Dict(Int(t) => length(e) for d in 0:3
                      for (t, e) in zip(A.mesh.get_elements(d)[1:2]...))
        @test counts == Dict(15 => 4, 1 => 84, 2 => 314, 3 => 12, 4 => 410)
        @test allunique(Tuple(xyz[3i-2:3i]) for i in eachindex(tags))
    finally
        A.finalize()
    end
end

end # module
