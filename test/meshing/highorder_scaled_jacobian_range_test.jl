module HighOrderScaledJacobianRangeTests
using Test,Tessella
const H=Tessella.HighOrder
@testset "Dimensionless high-order Jacobians retain finite coordinate scales" begin
    for msh in (9,10,11,12,13,14)
        reference=Tessella.Elements.lagrange_nodes(msh)
        cells=reshape(Int32.(axes(reference,2)),size(reference,2),1)
        fam=H._p2opt_family(Int32(msh))
        for shear in (0.,.25)
            shape=copy(reference)
            shape[1,:].+=shear.*shape[2,:]
            baseline=H._p2opt_scaled_jac(shape,fam,cells,1,.2,.2,.2)
            @test isfinite(baseline) && baseline>0
            for scale in (1e-300,1e-150,1.,1e150,1e200)
                coords=shape.*scale
                value=H._p2opt_scaled_jac(coords,fam,cells,1,.2,.2,.2)
                @test value≈baseline atol=2e-14
                if fam.volume
                    coords[1,:].*=-1
                    @test H._p2opt_scaled_jac(coords,fam,cells,1,.2,.2,.2)≈-baseline atol=2e-14
                end
            end
        end
        @test H._p2opt_scaled_jac(zero(reference),fam,cells,1,.2,.2,.2)==-Inf
    end
end
end
