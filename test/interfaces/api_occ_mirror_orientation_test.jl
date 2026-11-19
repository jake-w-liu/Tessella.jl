module OCCMirrorOrientationChecks
using Test,Tessella
const A=Tessella.API

# A reflection maps a cylinder's (sphere's, cone's) frame vector by vector
# and leaves it left-handed, so the face's natural normal turns against the
# image of the old one. The shell kept using the face the same way round:
# mirroring a box with a cylindrical hole changed its volume from 1.7173 to
# 2.0942, and a mirrored cylinder cut by a box measured −0.2618. Both values
# are closed forms.
@testset "OCC mirror keeps curved faces' orientation in their shells" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,1,1)
        drill=occ.add_cylinder(0.5,0.5,-0.5,0,0,2,0.3)
        out,_=occ.cut([(3,box)],[(3,drill)])
        occ.mirror(out,1,0,0,0.5)
        @test occ.get_mass(3,out[1][2])≈2-π*0.09 rtol=1e-12
        @test collect(occ.get_center_of_mass(3,out[1][2]))≈
            [-1-(2-0.5*π*0.09)/(2-π*0.09),0.5,0.5] rtol=1e-12
    finally
        A.finalize()
    end

    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        rod=occ.add_cylinder(0,0,0,0,0,2,0.5)
        occ.mirror([(3,rod)],1,0,0,0.5)
        top=occ.add_box(-1.5,-1,1,2,2,2)
        out,_=occ.cut([(3,rod)],[(3,top)])
        @test occ.get_mass(3,out[1][2])≈π*0.25 rtol=1e-12
    finally
        A.finalize()
    end
end
end
