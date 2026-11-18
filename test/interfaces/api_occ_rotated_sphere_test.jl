module OCCRotatedSphereChecks
using Test,Tessella
const A=Tessella.API

# The sphere encoding stores a center and a radius, no polar axis. Its
# consistency check nevertheless required the pole vertices on the world z
# diameter, so any sphere rotated about another axis was reported as
# "inconsistent primitive encoding and shell topology" by
# `model.get_bounding_box` (and every query built on it).
@testset "OCC sphere stays consistent after a generic rotation" begin
    for transform in (
            occ->occ.rotate([(3,1)],0.1,0.2,0.3,1,1,0.5,0.7),
            occ->occ.rotate([(3,1)],0,0,0,1,0,0,π/2),
            occ->occ.mirror([(3,1)],1,1,0.5,0.25))
        A.initialize()
        try
            A.option("General.Terminal",0)
            occ=A.model.occ
            occ.add_sphere(0.3,-0.2,0.1,0.6)
            transform(occ)
            occ.synchronize()
            center=collect(occ.get_center_of_mass(3,1))
            box=collect(A.model.get_bounding_box(3,1))
            @test box≈vcat(center.-0.6,center.+0.6) atol=1e-9
            @test occ.get_mass(3,1)≈4π/3*0.6^3 rtol=1e-12
        finally
            A.finalize()
        end
    end
end
end
