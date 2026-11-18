module OCCDilateMeasureChecks
using Test,Tessella
const A=Tessella.API

total_area(occ)=sum(occ.get_mass(2,t) for (_,t) in A.model.get_entities(2))

# `occ.dilate` stretches a plane record's frame vectors and leaves its pcurves
# in the old uv units, so a uv shoelace area must be scaled by |X×Y| and the
# stored `axis` is no longer a unit normal. Plane-face areas used to stay at
# their pre-dilation values and every solid with planar faces reported a
# wrong mass (a 2×1×1 box dilated by 2 claimed volume 8 instead of 16).
@testset "OCC measures after dilation" begin
    cases=(
        ("box ×2",occ->occ.add_box(0,0,0,2,1,1),(2.,2.,2.),16.0,40.0),
        ("box anisotropic",occ->occ.add_box(0,0,0,2,1,1),(2.,1.,.5),2.0,13.0),
        ("wedge ×2",occ->occ.add_wedge(0,0,0,1,1,1),(2.,2.,2.),4.0,
            4*(3+sqrt(2))),
        ("cylinder ×2",occ->occ.add_cylinder(.5,.5,-.5,0,0,2,.3),(2.,2.,2.),
            π*.6^2*4,2π*.6*4+2π*.6^2),
        ("cone ×2",occ->occ.add_cone(0,0,0,0,0,1,.5,.2),(2.,2.,2.),
            π*2*(1+.4+.16)/3,π*1.4*hypot(2,.6)+π*(1+.16)),
        ("box ×0.5",occ->occ.add_box(0,0,0,2,1,1),(.5,.5,.5),0.25,2.5),
    )
    for (name,make,factors,volume,area) in cases
        @testset "$name" begin
            A.initialize()
            try
                A.option("General.Terminal",0)
                occ=A.model.occ
                tag=make(occ)
                occ.dilate([(3,tag)],0.1,0.2,0.3,factors...)
                occ.synchronize()
                @test occ.get_mass(3,tag)≈volume rtol=1e-12
                @test total_area(occ)≈area rtol=1e-12
                lo=A.model.get_bounding_box(3,tag)
                center=collect(occ.get_center_of_mass(3,tag))
                @test all(lo[i]<center[i]<lo[i+3] for i in 1:3)
            finally
                A.finalize()
            end
        end
    end
    # the centroid follows the dilation exactly
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,1,1)
        occ.dilate([(3,box)],0.1,0.2,0.3,2,2,2)
        @test collect(occ.get_center_of_mass(3,box))≈[1.9,0.8,0.7] atol=1e-12
        inertia=occ.get_matrix_of_inertia(3,box)
        # solid box 4×2×2 about its centroid, unit density
        @test [inertia[1],inertia[5],inertia[9]]≈
            16/12 .*[2^2+2^2,4^2+2^2,4^2+2^2] rtol=1e-12
    finally
        A.finalize()
    end
end
end
