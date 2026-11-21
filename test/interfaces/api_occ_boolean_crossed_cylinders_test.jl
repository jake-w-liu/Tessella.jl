module OCCBooleanCrossedCylindersChecks
using Test,Tessella
const A=Tessella.API

# Two cylinders crossing at right angles, radii 0.5 and 0.3: the fuse once
# measured 0.63. The common part of two perpendicular cylinders through one
# axis point is 8 ∫₀^r √(r²−x²) √(R²−x²) dx = 0.2693703275 for these radii,
# so fuse = π/2 + 0.18π − 0.2693703275 and cut = π/2 − 0.2693703275. Gmsh
# 4.15.2 gives 1.8669200576 and 1.3014335187, off in the sixth digit. The
# section is a traced curve, minutes when interpreted: compiled runs only.
@testset "OCC Booleans of two crossed cylinders" begin
    if Base.JLOptions().compile_enabled==1
        common=0.2693703275
        for (operation,exact) in ((:fuse,π/2+0.18π-common),
                (:cut,π/2-common),(:intersect,common))
            A.initialize()
            try
                A.option("General.Terminal",0)
                occ=A.model.occ
                upright=occ.add_cylinder(0,0,0,0,0,2,0.5)
                across=occ.add_cylinder(-1,0,1,2,0,0,0.3)
                out,_=getproperty(occ,operation)([(3,upright)],[(3,across)])
                @test length(out)==1
                @test occ.get_mass(3,out[1][2])≈exact atol=2e-6
            finally
                A.finalize()
            end
        end
    end
end
end
