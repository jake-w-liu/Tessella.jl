module APIPartitionRegenerateTest

using Test, Tessella

const API = Tessella.API

# Generating again after PartitionMesh (tutorials/t21.geo through open +
# generate): Gmsh 4.15.2 meshes the original entities again, deletes the
# partition surfaces' cells in generate(1) and never remeshes them, and keeps
# the partition curves and points. The counts are Gmsh's for this model, with
# Surface 1 and its own Curves and Points in partition 1, the rest in 2.
function counts()
    surfaces = [(t, length(API.mesh.get_nodes(2, t)[1]),
                 sum(length, API.mesh.get_elements(2, t)[2]; init=0))
                for (_, t) in API.model.get_entities(2)]
    cells = [sum(length, API.mesh.get_elements(k)[2]; init=0) for k in 0:2]
    return length(API.mesh.get_nodes()[1]), cells, surfaces
end

@testset "generate after partitioning keeps Gmsh's partition semantics" begin
    API.initialize()
    try
        API.option("General.Terminal", 0)
        g = API.model.geo
        for (i, (x, y)) in enumerate(((0, 0), (1, 0), (1, 1), (0, 1), (2, 0), (2, 1)))
            g.add_point(x, y, 0, 0.25, i)
        end
        for (i, (a, b)) in enumerate(((1, 2), (2, 3), (3, 4), (4, 1), (2, 5), (5, 6), (6, 3)))
            g.add_line(a, b, i)
        end
        g.add_curve_loop([1, 2, 3, 4], 1); g.add_plane_surface([1], 1)
        g.add_curve_loop([5, 6, 7, -2], 2); g.add_plane_surface([2], 2)
        g.synchronize()
        API.mesh.generate(2)
        API.option("Mesh.PartitionCreateTopology", 1)
        tags = Int[]; parts = Int[]
        for (d, t) in API.model.get_entities()
            first = (d == 2 && t == 1) || (d == 1 && t in (1, 3, 4)) || (d == 0 && t in (1, 4))
            for block in API.mesh.get_elements(d, t)[2], x in block
                push!(tags, Int(x)); push!(parts, first ? 1 : 2)
            end
        end
        API.mesh.partition(2, tags, parts)
        @test counts()[3] == [(1, 0, 0), (2, 0, 0), (3, 14, 42), (4, 15, 44)]
        API.mesh.generate(1)
        @test counts()[3] == [(1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 0, 0)]
        API.mesh.generate(2)
        @test counts()[3] == [(1, 14, 42), (2, 15, 44), (3, 0, 0), (4, 0, 0)]
        @test counts()[2][3] == 86
        # Gmsh's partition topology adds 4 interface lines over Curve 2 and
        # the partition points; they survive the new generates.
        @test counts()[1:2] == (85, [17, 60, 86])
    finally
        API.finalize()
    end
end

end
