# `Mesh.AnisoMax` — the BGM metric aspect-ratio bound (`metric_aniso_cap`)
# and its transport through the model meshing attributes.

using Test
using Tessella

const _cap=Tessella.SizeField.metric_aniso_cap
const _from_axes=Tessella.SizeField.metric_from_axes

@testset "metric anisotropy cap" begin
    x=(1.0,0.0,0.0); y=(0.0,1.0,0.0); z=(0.0,0.0,1.0)
    # h_x=1, h_y=h_z=0.001 → λ=(1,1e6,1e6), physical aspect ratio 1000.
    m=_from_axes(1.0,1e6,1e6,x,y,z)

    capped=_cap(m,100.0)
    λ=sort(collect(metric_eigenvalues(capped)))
    @test sqrt(λ[3]/λ[1])<=100.0*(1+8eps(Float64))
    @test λ[1]≈1.0            # coarsest requested size never refines
    @test λ[2]≈1e4 && λ[3]≈1e4 # λmin·aniso_max² cap on both fine directions
    @test capped.m11≈1.0 && capped.m22≈1e4 && capped.m33≈1e4

    # A single over-fine direction clamps independently.
    stretched=_from_axes(1.0,4.0,1e6,x,y,z)
    capped2=_cap(stretched,10.0)
    λ2=sort(collect(metric_eigenvalues(capped2)))
    @test λ2≈[1.0,4.0,100.0]

    # Rotated eigenframe: the cap preserves the coarse direction exactly.
    u1=(0.6,0.8,0.0); u2=(-0.8,0.6,0.0); u3=(0.0,0.0,1.0)
    tilted=_from_axes(1.0,1e6,1e6,u1,u2,u3)
    tcapped=_cap(tilted,50.0)
    tλ=sort(collect(metric_eigenvalues(tcapped)))
    @test tλ≈[1.0,2500.0,2500.0]
    # u1 stays the coarse eigenaxis: xᵀMx on u1 is still λmin.
    @test Tessella.SizeField._metric_quad(tcapped,u1...)≈1.0

    # In-bound and unbounded inputs return the metric unchanged.
    iso=isotropic_metric(1.0)
    @test _cap(iso,10.0)===iso
    @test _cap(m,1e33)===m            # upstream unbounded default
    @test _cap(m,Inf)===m
    @test _cap(m,NaN)===m
    @test _cap(m,0.0)===m
    @test _cap(m,-2.0)===m

    # The model meshing attribute carries Mesh.AnisoMax into the 1-D options.
    model=GeoModel()
    @test model.meshing.aniso_max==1e33
    model.meshing.aniso_max=7.0
    options=Tessella.Model._model_default_mesh1d_options(model,"test")
    @test options.aniso_max==7.0
end
