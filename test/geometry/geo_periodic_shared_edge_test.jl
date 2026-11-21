module GeoPeriodicSharedEdgeTest

using Test, Tessella

# Two side faces of an extruded square sharing the edge on the z axis, made
# periodic by the rotation about that edge: the shared points and the shared
# curve are their own images. Gmsh 4.15.2 sets only the three curve masters
# 1-4, 12-20 and 6-9, and meshes 132 nodes / 260 triangles in 2-D with 31
# periodic node pairs on Surface 13 (5 of them the shared edge's own nodes),
# and 141 nodes / 373 tetrahedra in 3-D.
const GEO = """
    Point(1)={0,0,0,0.25}; Point(2)={1,0,0,0.25}; Point(3)={1,1,0,0.25}; Point(4)={0,1,0,0.25};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    Extrude {0,0,1} { Surface{1}; }
    Periodic Surface{13} = {25} Rotate{{0,0,1},{0,0,0},-Pi/2};
    """

function run(dim)
    path = joinpath(mktempdir(), "shared.geo")
    write(path, GEO)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(dim)
        nodes = length(A.mesh.get_nodes()[1])
        types, tags, _ = A.mesh.get_elements(dim)
        cells = sum(length, tags)
        master, slaves, masters, _ = A.mesh.get_periodic_nodes(2, 13)
        return nodes, cells, master, length(slaves), count(slaves .== masters)
    finally
        A.finalize()
    end
end

@testset "periodic faces sharing an edge" begin
    nodes, triangles, master, pairs, fixed = run(2)
    @test (nodes, triangles) == (132, 260)   # Gmsh 4.15.2
    @test master == 25
    @test pairs == 31                        # Gmsh 4.15.2: 31
    @test fixed == 5                         # Gmsh 4.15.2: 5 self pairs
    nodes, tets, _, _, _ = run(3)
    @test nodes == 141                       # Gmsh 4.15.2: 141
    @test 360 <= tets <= 390                 # Gmsh 4.15.2: 373 (Tessella 375)
end

end # module
