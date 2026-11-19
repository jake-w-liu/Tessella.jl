module OCCBooleanContainmentChecks
using Test,Tessella
const A=Tessella.API

# A tool lying strictly inside the object cuts none of the object's faces,
# yet the result is a different solid: the tool's boundary becomes an inner
# (cavity) shell. `occ.cut` used to return the object unchanged (6 faces,
# volume 8). Counts are Gmsh 4.15.2's; volume and area are closed forms.
# One planar Boolean is the cheapest input that shows it.
@testset "OCC cut with the tool fully inside the object" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        tool=occ.add_box(0.6,0.6,0.6,0.8,0.8,0.8)
        out,images=occ.cut([(3,box)],[(3,tool)])
        @test out==[(3,1)]
        @test images==[[(3,1)],Tuple{Int,Int}[]]
        @test Tuple(length(occ.get_entities(d)) for d in 0:3)==(16,24,12,1)
        @test occ.get_mass(3,1)≈8-0.512 rtol=1e-12
        @test sum(occ.get_mass(2,t) for (_,t) in occ.get_entities(2))≈
            24+6*0.64 rtol=1e-12
    finally
        A.finalize()
    end
end
end
