module TransfiniteVolumeGmshCornersTest

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

# A first face whose GEdgeLoop runs clockwise seen from inside gives Gmsh a
# left-handed (i,j,k) tab; the kernel runs on the right-handed relabelling
# with the `CREATE_SIM` templates mirrored. Gmsh 4.15.2's tetrahedra.
const _GMSH_LEFT_HANDED_TETS = Set([
    (0,1,3,9),(1,2,4,10),(1,3,4,9),(1,4,9,10),(2,4,5,10),(2,5,10,11),
    (3,4,6,12),(3,4,9,12),(4,5,7,13),(4,5,10,13),(4,6,7,12),(4,7,12,13),
    (4,9,10,12),(4,10,12,13),(5,7,8,13),(5,8,13,14),(5,10,11,13),(5,11,13,14),
    (6,7,12,15),(7,8,13,16),(7,12,13,15),(7,13,15,16),(8,13,14,16),(8,14,16,17),
    (9,10,12,18),(10,11,13,19),(10,12,13,18),(10,13,18,19),(11,13,14,19),
    (11,14,19,20),(12,13,15,21),(12,13,18,21),(13,14,16,22),(13,14,19,22),
    (13,15,16,21),(13,16,21,22),(13,18,19,21),(13,19,21,22),(14,16,17,22),
    (14,17,22,23),(14,19,20,22),(14,20,22,23),(15,16,21,24),(16,17,22,25),
    (16,21,22,24),(16,22,24,25),(17,22,23,25),(17,23,25,26)])

const _RIGHT_HANDED_GEO = """
    Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={3,4,1,2}; Plane Surface(1)={1};
    Extrude {0,0,1} { Surface{1}; }
    Transfinite Curve{:} = 3;
    Transfinite Surface{:};
    Transfinite Volume{1};
    """

const _LEFT_HANDED_GEO = """
    Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
    Point(5)={0,0,1}; Point(6)={1,0,1}; Point(7)={1,1,1}; Point(8)={0,1,1};
    Line(1)={1,2}; Line(2)={3,2}; Line(3)={3,4}; Line(4)={4,1};
    Line(5)={5,6}; Line(6)={6,7}; Line(7)={8,7}; Line(8)={8,5};
    Line(9)={1,5}; Line(10)={2,6}; Line(11)={3,7}; Line(12)={4,8};
    Curve Loop(1)={4,1,-2,3}; Plane Surface(1)={1};
    Curve Loop(2)={5,6,-7,8}; Plane Surface(2)={2};
    Curve Loop(3)={1,10,-5,-9}; Plane Surface(3)={3};
    Curve Loop(4)={-2,11,-6,-10}; Plane Surface(4)={4};
    Curve Loop(5)={3,12,7,-11}; Plane Surface(5)={5};
    Curve Loop(6)={4,9,-8,-12}; Plane Surface(6)={6};
    Surface Loop(1)={3,1,2,4,5,6}; Volume(1)={1};
    Transfinite Curve{:} = 3;
    Transfinite Surface{:};
    Transfinite Volume{1};
    """

function lattice_tets(source)
    A = Tessella.API
    path = joinpath(mktempdir(), "tv.geo")
    write(path, source)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(3)
        tags, xyz, _ = A.mesh.get_nodes()
        lattice = Dict(t => round(Int, 2xyz[3i-2]) + 3round(Int, 2xyz[3i-1]) +
                            9round(Int, 2xyz[3i]) for (i, t) in enumerate(tags))
        types, _, nodes = A.mesh.get_elements(3)
        types == [4] || return nothing
        return Set(Tuple(sort([lattice[n] for n in nodes[1][q:q+3]]))
                   for q in 1:4:length(nodes[1]))
    finally
        A.finalize()
    end
end

@testset "transfinite Volume default corners follow Gmsh" begin
    @test lattice_tets(_RIGHT_HANDED_GEO) == _GMSH_TETS
    @test lattice_tets(_LEFT_HANDED_GEO) == _GMSH_LEFT_HANDED_TETS
end

end # module
