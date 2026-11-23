module APIRenumberNodesPhysicalChecks

using Test, Tessella

const A = Tessella.API

# Two squares, a physical group on the second only. Gmsh's generate and its
# renumberNodes (GModel::renumberMeshVertices for a potentially partial save)
# number the nodes of the physical entities' elements first, each entity in
# order, then the others: renumbering right after generate changes nothing.
# Tags per entity are Gmsh 4.15.2's.
const SOURCE = """
Point(1)={0,0,0,0.5}; Point(2)={1,0,0,0.5}; Point(3)={1,1,0,0.5}; Point(4)={0,1,0,0.5};
Point(5)={2,0,0,0.5}; Point(6)={2,1,0,0.5};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Line(5)={2,5}; Line(6)={5,6}; Line(7)={6,3};
Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
Curve Loop(2)={5,6,7,-2}; Plane Surface(2)={2};
Physical Surface(10)={2};
"""

@testset "renumber_nodes after generate keeps Gmsh's physical-first numbering" begin
    mktempdir() do directory
        path = joinpath(directory, "two_squares.geo")
        write(path, SOURCE)
        A.initialize()
        try
            A.option("General.Terminal", 0)
            A.open(path)
            A.mesh.generate(2)
            tags, coords, _ = A.mesh.get_nodes()
            @test length(tags) == 21
            expected = Dict((0, 1) => [13], (0, 2) => [1], (0, 4) => [14], (0, 5) => [3],
                            (1, 1) => [15], (1, 4) => [17], (1, 5) => [6],
                            (2, 1) => [18, 19, 20, 21], (2, 2) => [9, 10, 11, 12])
            for ((dim, tag), want) in expected
                @test Int.(A.mesh.get_nodes(dim, tag, false, false)[1]) == want
            end
            A.mesh.renumber_nodes()
            after, moved, _ = A.mesh.get_nodes()
            @test after == tags
            @test moved == coords
        finally
            A.finalize()
        end
    end
end

end
