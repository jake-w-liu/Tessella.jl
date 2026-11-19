module OCCBooleanTagBaseChecks
using Test,Tessella
const A=Tessella.API

# OpenCASCADE unbinds the removed operands and recomputes every dimension's
# highest tag before it binds a Boolean's result, so the result's faces,
# curves and points continue from what is still bound. Tessella kept counting
# from the operands' tags: two boxes cut from one another left Surfaces 13–21
# where Gmsh 4.15.2 has 1–9 (and Curves 1–21, Points 1–14), so every tag a
# script used after a Boolean was off.
@testset "OCC Boolean results are numbered from the surviving tags" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_box(0,0,0,1,1,1)
        occ.add_box(0.5,0.5,0.5,1,1,1)
        out,_=occ.cut([(3,1)],[(3,2)])
        occ.synchronize()
        @test out==[(3,1)]
        @test [[t for (_,t) in A.model.get_entities(d)] for d in 0:2]==
            [collect(1:14),collect(1:21),collect(1:9)]
        @test occ.get_mass(3,1)≈0.875 rtol=1e-13
    finally
        A.finalize()
    end
end
end
