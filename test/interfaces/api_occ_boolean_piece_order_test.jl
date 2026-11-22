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
        # Gmsh's result list visits the walls 7,8,7,8,7,8; the order of a
        # rebuilt solid's faces is not reproduced yet
        @test_broken out==[(3,1),(3,2),(3,3),(2,7),(2,8),(2,7),(2,8),(2,7),(2,8)]
    end
end
end
