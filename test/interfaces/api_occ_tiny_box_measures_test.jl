module OCCTinyBoxMeasuresChecks
using Test,Tessella
const A=Tessella.API

@testset "OCC measures of a box with 1e-6 edges" begin
    # Faces of 1e-12 and a volume of 1e-18 were "zero": the measures used
    # `Precision::Confusion` (1e-7, a length) as their area and volume
    # threshold where OpenCASCADE's is 1e-20. Gmsh 4.15.2 reports this box's
    # mass 1.0000000002e-18 and face area 1.000000000279556e-12, refuses
    # edges of 1e-7, and puts the centre of a volume below 1e-20 at the
    # origin.
    A.initialize()
    try
        occ=A.model.occ
        box=occ.add_box(1,2,3,1e-6,1e-6,1e-6)
        @test occ.get_mass(3,box)≈1e-18 rtol=1e-8
        @test occ.get_mass(2,1)≈1.000000000279556e-12 rtol=1e-12
        @test collect(occ.get_center_of_mass(3,box))≈
            [1.0000005,2.0000005,3.0000005] atol=1e-9
        @test collect(occ.get_center_of_mass(2,1))≈
            [1.0,2.0000005,3.0000005] atol=1e-12
        @test_throws ArgumentError occ.add_box(0,0,0,1e-7,1,1)
        small=occ.add_box(0,0,0,2e-7,2e-7,2e-7)
        @test occ.get_mass(3,small)≈8e-21 rtol=1e-8
        @test occ.get_center_of_mass(3,small)==(0.0,0.0,0.0)
    finally
        A.finalize()
    end
end
end
