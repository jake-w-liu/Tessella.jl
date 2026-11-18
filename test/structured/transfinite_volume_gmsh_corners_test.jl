using Test
using Tessella

# Unpinned transfinite Volume corners follow upstream
# `findTransfiniteCorners(GRegion*)`: the first bounding face's GEdgeLoop
# order (first curve with sign +1) fixes the (i,j,k) tab and hence the
# six-tetrahedron split. The source loop {3,4,1,2} starts at Point 3, not at
# the smallest Point tag. Expected tetrahedra are Gmsh 4.15.2's, as indices
# 9z+3y+x on the 3x3x3 lattice of node coordinates {0,0.5,1}.
const _GMSH_TETS = Set([
    (0,1,3,12),(0,1,10,12),(0,9,10,12),(1,2,4,13),(1,2,11,13),(1,3,4,13),
    (1,3,12,13),(1,10,11,13),(1,10,12,13),(2,4,5,14),(2,4,13,14),(2,11,13,14),
    (3,4,6,15),(3,4,13,15),(3,12,13,15),(4,5,7,16),(4,5,14,16),(4,6,7,16),
    (4,6,15,16),(4,13,14,16),(4,13,15,16),(5,7,8,17),(5,7,16,17),(5,14,16,17),
    (9,10,12,21),(9,10,19,21),(9,18,19,21),(10,11,13,22),(10,11,20,22),
    (10,12,13,22),(10,12,21,22),(10,19,20,22),(10,19,21,22),(11,13,14,23),
    (11,13,22,23),(11,20,22,23),(12,13,15,24),(12,13,22,24),(12,21,22,24),
    (13,14,16,25),(13,14,23,25),(13,15,16,25),(13,15,24,25),(13,22,23,25),
    (13,22,24,25),(14,16,17,26),(14,16,25,26),(14,23,25,26)])

@testset "transfinite Volume default corners follow Gmsh" begin
    A = Tessella.API
    path = joinpath(mktempdir(), "tv.geo")
    write(path, """
        Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Curve Loop(1)={3,4,1,2}; Plane Surface(1)={1};
        Extrude {0,0,1} { Surface{1}; }
        Transfinite Curve{:} = 3;
        Transfinite Surface{:};
        Transfinite Volume{1};
        """)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(3)
        tags, xyz, _ = A.mesh.get_nodes()
        lattice = Dict(t => round(Int, 2xyz[3i-2]) + 3round(Int, 2xyz[3i-1]) +
                            9round(Int, 2xyz[3i]) for (i, t) in enumerate(tags))
        types, _, nodes = A.mesh.get_elements(3)
        k = findfirst(==(4), types)
        @test k !== nothing
        tets = Set(Tuple(sort([lattice[n] for n in nodes[k][q:q+3]]))
                   for q in 1:4:length(nodes[k]))
        @test length(tets) == 48
        @test tets == _GMSH_TETS
    finally
        A.finalize()
    end
end
