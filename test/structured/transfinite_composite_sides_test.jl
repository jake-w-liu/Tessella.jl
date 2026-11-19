module TransfiniteCompositeSidesTest

using Test, Tessella

# Explicit transfinite corners may cut a boundary of more than four curves into
# sides that chain several curves (upstream `MeshTransfiniteSurface`; Gmsh
# tutorial t6). Here the bottom side chains Lines 1 and 2. Gmsh 4.15.2's
# triangles on the lattice of node coordinates (3x, 4y), indices i + 7j.
const _GMSH_TRIANGLES = Set([
    (0,1,7),(1,2,8),(1,7,8),(2,3,9),(2,8,9),(3,4,10),(3,9,10),(4,5,11),(4,10,11),
    (5,6,12),(5,11,12),(6,12,13),(7,8,14),(8,9,15),(8,14,15),(9,10,16),(9,15,16),
    (10,11,17),(10,16,17),(11,12,18),(11,17,18),(12,13,19),(12,18,19),(13,19,20),
    (14,15,21),(15,16,22),(15,21,22),(16,17,23),(16,22,23),(17,18,24),(17,23,24),
    (18,19,25),(18,24,25),(19,20,26),(19,25,26),(20,26,27),(21,22,28),(22,23,29),
    (22,28,29),(23,24,30),(23,29,30),(24,25,31),(24,30,31),(25,26,32),(25,31,32),
    (26,27,33),(26,32,33),(27,33,34)])

@testset "transfinite surface sides chaining several curves" begin
    path = joinpath(mktempdir(), "corners.geo")
    write(path, """
        Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={2,0,0}; Point(4)={2,1,0}; Point(5)={0,1,0};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,5}; Line(5)={5,1};
        Curve Loop(1)={1,2,3,4,5}; Plane Surface(1)={1};
        Transfinite Curve{1,2} = 4; Transfinite Curve{4} = 7; Transfinite Curve{3,5} = 5;
        Transfinite Surface{1} = {1,3,4,5};
        """)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(2)
        tags, xyz, _ = A.mesh.get_nodes()
        @test length(tags) == 35
        lattice = Dict(t => round(Int, 3xyz[3i-2]) + 7round(Int, 4xyz[3i-1])
                       for (i, t) in enumerate(tags))
        @test all(i -> abs(3xyz[3i-2] - round(3xyz[3i-2])) < 1e-9 &&
                       abs(4xyz[3i-1] - round(4xyz[3i-1])) < 1e-9, eachindex(tags))
        types, _, nodes = A.mesh.get_elements(2)
        @test types == [2]
        triangles = Set(Tuple(sort([lattice[n] for n in nodes[1][q:q+2]]))
                        for q in 1:3:length(nodes[1]))
        @test triangles == _GMSH_TRIANGLES
        # The boundary lines keep their own curves.
        @test length(A.mesh.get_elements(1, 1)[2][1]) == 3
        @test length(A.mesh.get_elements(1, 2)[2][1]) == 3
    finally
        A.finalize()
    end
end

end # module
