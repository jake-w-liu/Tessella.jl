module OCCBooleanPieceOrderChecks
using Test,Tessella
const A=Tessella.API

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        f(A.model.occ)
    finally
        A.finalize()
    end
end

# OpenCASCADE builds an operand's split solids from that operand's faces in
# order, the piece at the first vertex of a split face first; the result
# list follows the lowest operand's pieces, and each operand's image lists
# its own pieces in the same way. Values are Gmsh 4.15.2's. (A box with two
# spheres over its corners, out 2..6 and images [[2,3,4],[5,2],[4,6]] as in
# Gmsh, is checked outside the suite for time; geo_occ_boolean_piece_bounds
# covers the box∩sphere-first order.)
@testset "OCC fragment piece and image order" begin
    session() do occ
        occ.add_box(0.5,0.5,0.5,0.5,0.5,0.5)
        occ.add_box(0,0,0,1,1,1)
        occ.add_box(0.5,0.5,0.5,1,1,1)
        out,images=occ.fragment([(3,1),(3,2)],[(3,3)])
        @test out==[(3,1),(3,2),(3,3)]
        @test images==[[(3,1)],[(3,2),(3,1)],[(3,1),(3,3)]]
    end
end

# A sheet that cuts a box whole stays the same face, with its own edges and
# vertices (Gmsh 4.15.2: Surface 7 keeps Curves 13-16 and Points 9-12), and
# the fresh entities continue after them.
@testset "OCC fragment keeps whole walls" begin
    session() do occ
        occ.add_box(0,0,0,1,1,1)
        occ.add_rectangle(0,0,0.3,1,1)
        occ.add_rectangle(0,0,0.7,1,1)
        occ.synchronize()
        out,images=occ.fragment([(3,1)],[(2,7),(2,8)])
        occ.synchronize()
        @test images==[[(3,1),(3,2),(3,3)],[(2,7)],[(2,8)]]
        @test [t for (_,t) in A.model.get_entities(0)]==collect(9:24)
        @test [t for (_,t) in A.model.get_entities(1)]==collect(13:40)
        @test [t for (_,t) in A.model.get_entities(2)]==collect(7:22)
        @test [b[2] for b in A.model.get_boundary([(2,7)],false,true)]==[13,14,15,16]
        @test [b[2] for b in A.model.get_boundary([(2,8)],false,true)]==[17,18,19,20]
        # Gmsh's result list visits the walls 7,8,7,8,7,8: each slab's shell
        # grows breadth-first from its first face, a face's neighbours in
        # its stored edge order
        @test out==[(3,1),(3,2),(3,3),(2,7),(2,8),(2,7),(2,8),(2,7),(2,8)]
        B(d,t)=[b[2] for b in A.model.get_boundary([(d,t)],false,true)]
        @test B(3,1)==[-9,-10,7,11,-12,13]
        @test B(3,2)==[-14,-15,8,16,-7,17]
        @test B(3,3)==[-18,-19,20,21,-8,22]
        @test B(2,9)==[-21,23,22,16]
        @test B(2,14)==[-29,-16,30,20]
        @test B(2,18)==[-33,-20,35,-34]
        # a split face starts at the first piece of its first split source
        # edge in stored order: the y=0 face's x=1 edge, not its bottom
        @test B(2,10)==[-24,-25,21,13]
    end
end

# Two overlapping boxes. A rebuilt face keeps its source face's orientation
# (an OpenCASCADE box's planes all face +x, +y, +z) and starts at the first
# piece of its first split source edge; fresh edges are numbered in the
# stored order of each wire (its first edge, then the rest backwards). All
# boundaries are Gmsh 4.15.2's (`getBoundary(combined=false)`).
@testset "OCC cut and fuse keep OpenCASCADE's faces" begin
    B(d,t)=[b[2] for b in A.model.get_boundary([(d,t)],false,true)]
    session() do occ
        occ.add_box(0,0,0,1,1,1)
        occ.add_box(0.5,0.5,0.5,1,1,1)
        occ.cut([(3,1)],[(3,2)])
        occ.synchronize()
        @test B(3,1)==[-1,-2,3,4,-5,6,7,8,9]
        @test [B(2,s) for s in 1:9]==[[-1,4,3,-2],[-5,1,7,-6],
            [-8,-2,7,11,10,9],[-12,-15,3,8,14,13],[-4,5,16,-15],
            [-11,-6,16,12,18,17],[-19,20,-14,-9],[-21,19,-10,-17],
            [-20,21,-18,-13]]
    end
    session() do occ
        occ.add_box(0,0,0,1,1,1)
        occ.add_box(0.5,0.5,0.5,1,1,1)
        occ.fuse([(3,1)],[(3,2)])
        occ.synchronize()
        @test B(3,1)==[-1,-2,3,4,-5,6,-7,-8,-9,10,11,12]
        @test [B(2,s) for s in 7:12]==[[-19,9,14,22,21,-20],
            [-23,17,10,19,25,-24],[-22,13,18,23,27,-26],[-20,25,29,-28],
            [-26,21,28,-30],[-24,27,30,-29]]
        @test occ.get_mass(3,1)≈1.875
    end
end

# A sphere over a box corner. The piece of the sphere outside the box was
# measured 0.589 (9/8 of the octant too much): the containment rays wrapped
# hits on its 7/8-sphere face a period outside the face. Gmsh 4.15.2:
# 0.934550153, 0.065449847, 0.458148929.
@testset "OCC fragment of a sphere over a box corner" begin
    session() do occ
        occ.add_box(0,0,0,1,1,1)
        occ.add_sphere(1,1,1,0.5)
        out,_=occ.fragment([(3,1)],[(3,2)])
        @test out==[(3,1),(3,2),(3,3)]
        @test occ.get_mass(3,1)≈1-π/48 rtol=1e-9
        @test occ.get_mass(3,2)≈π/48 rtol=1e-9
        @test occ.get_mass(3,3)≈7π/48 rtol=1e-9
    end
end
end
