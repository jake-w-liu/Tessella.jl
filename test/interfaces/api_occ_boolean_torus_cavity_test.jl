module OCCBooleanTorusCavityChecks
using Test,Tessella
const A=Tessella.API

# A solid bounded by one torus face and plane faces was always measured as a
# torus sector (Pappus), so a box with a toroidal cavity reported the
# cavity's volume instead of its own.
@testset "OCC volume of a box with a toroidal cavity" begin
    torus_volume=2π^2*0.26*0.14^2
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        torus=occ.add_torus(0.26,0,0,0.26,0.14)
        out,_=occ.cut([(3,box)],[(3,torus)])
        @test out==[(3,1)]
        @test occ.get_mass(3,1)≈8-torus_volume rtol=1e-6
        @test collect(occ.get_center_of_mass(3,1))≈
            -torus_volume.*[0.26,0.0,0.0]./(8-torus_volume) atol=1e-6
        inertia=occ.get_matrix_of_inertia(3,1)
        @test all(>(4.0),(inertia[1],inertia[5],inertia[9]))
    finally
        A.finalize()
    end
    # sector solids keep the closed form
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_torus(0,0,0,1,0.3)
        occ.add_torus(4,0,0,1,0.3;angle=π/2)
        @test occ.get_mass(3,1)≈2π^2*0.09 rtol=1e-12
        @test occ.get_mass(3,2)≈π^2*0.09/2 rtol=1e-12
    finally
        A.finalize()
    end
end
end
