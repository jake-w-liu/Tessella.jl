module PeriodicHighOrderNodesTests
using Test, Tessella
const A = Tessella.API

pairs_of(link) = Dict(Int(s) => Int(m) for (s, m) in zip(link.slave_nodes, link.master_nodes))

# Two three-node transfinite Lines, the second periodic to the first, meshed at
# order 2. Gmsh 4.15.2 `getPeriodicNodes(1, 2)` returns the three primary pairs
# {3=>1, 4=>2, 8=>5}; with `includeHighOrderNodes` it adds the two mid-edge
# pairs {9=>6, 10=>7}. After the mesh is written as MSH 4.1 and read back,
# Gmsh returns all five pairs with or without the flag. The branch that answers
# for a loaded file used to answer for the generated mesh too, so the primary
# query returned five pairs.
@testset "getPeriodicNodes separates primary and high-order pairs of a generated mesh" begin
    A.finalize(); A.initialize()
    try
        A.option("General.Terminal", 0)
        A.option("Mesh.ElementOrder", 2)
        for (tag, x, y) in ((1, 0.0, 0.0), (2, 1.0, 0.0), (3, 0.0, 2.0), (4, 1.0, 2.0))
            A.model.add_point(x, y, 0; tag)
        end
        A.model.add_line(1, 2; tag=1); A.model.add_line(3, 4; tag=2)
        for curve in 1:2
            A.mesh.set_transfinite_curve(curve, 3)
        end
        A.mesh.set_periodic(1, [2], [1], [1.0, 0, 0, 0, 0, 1.0, 0, 2.0, 0, 0, 1.0, 0, 0, 0, 0, 1.0])
        A.mesh.generate(1)
        primary = A.mesh.get_periodic_nodes(1, 2)
        full = A.mesh.get_periodic_nodes(1, 2, true)
        @test primary.master_entity == 1 && full.master_entity == 1
        @test pairs_of(primary) == Dict(3 => 1, 4 => 2, 8 => 5)
        @test pairs_of(full) == Dict(3 => 1, 4 => 2, 8 => 5, 9 => 6, 10 => 7)
        @test pairs_of(A.mesh.get_periodic_nodes(0, 3)) == Dict(3 => 1)
        @test pairs_of(A.mesh.get_periodic_nodes(0, 4)) == Dict(4 => 2)
        mktempdir() do directory
            path = joinpath(directory, "periodic_p2.msh")
            A.write(path)
            A.finalize(); A.initialize()
            A.option("General.Terminal", 0)
            A.open(path)
            for high in (false, true)
                link = A.mesh.get_periodic_nodes(1, 2, high)
                @test link.master_entity == 1
                @test pairs_of(link) == Dict(3 => 1, 4 => 2, 8 => 5, 9 => 6, 10 => 7)
            end
        end
    finally
        A.finalize()
    end
end
end
