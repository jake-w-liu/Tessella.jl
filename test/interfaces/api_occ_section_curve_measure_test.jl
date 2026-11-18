module OCCSectionCurveMeasureChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# A box face cutting a sphere off-center leaves a section circle that is
# neither a parallel nor a meridian of the sphere, so it is stored as a traced
# curve S(pc(t)) — e.g. `occ.cut(box(-1,-1,-1,2,2,2), sphere(.8,.55,.35,.6))`.
# `occ.get_mass(1, …)` on such an edge threw "OCC curve: quadrature did not
# converge": the centroid component normal to the cutting plane is pure
# evaluation noise and was held to a tolerance below rounding. The same curve
# is built here directly (a dense Hermite trace of the circle x=1/2 on the
# unit sphere) so the check does not have to run the section tracer.
@testset "OCC length of traced Boolean section curves" begin
    A.initialize()
    face=try
        A.option("General.Terminal",0)
        A.model.occ.add_sphere(0,0,0,1)
        kernel=A.CURRENT[].occ_kernel
        kernel.shapes[kernel.bindings[(2,1)]].geometry.geometry
    finally
        A.finalize()
    end
    ρ=sqrt(0.75)
    θs=collect(range(0.3,1.2,length=41))
    uv=NTuple{2,Float64}[(atan(ρ*cos(θ),0.5),asin(ρ*sin(θ))) for θ in θs]
    duv=NTuple{2,Float64}[(-0.5ρ*sin(θ)/(0.25+ρ^2*cos(θ)^2),
                           ρ*cos(θ)/sqrt(1-ρ^2*sin(θ)^2)) for θ in θs]
    knots,degree,poles=M._hermite_bspline_2d(θs,uv,duv)
    seeds=NTuple{3,Float64}[M._occ_surface_point(face,q...) for q in uv]
    polyline=M.NURBS.NURBSCurve(1,vcat(θs[1],θs,θs[end]),seeds)
    curve=M._ModelOCCSectionCurve(:bspline,face,
                                  (spl2d=(knots,degree,poles),),polyline)
    record=(occ=:copied_curve,curve=curve,t0=θs[1],t1=θs[end])

    # panels cut at the trace's own knots
    cuts=M._model_occ_curve_panels(curve,θs[1],θs[end])
    @test cuts==θs[2:end-1]
    @test M._model_occ_curve_panels(curve,θs[1],θs[2])==Float64[]

    stats=M._model_occ_curve_stats(record)
    @test stats.area≈ρ*0.9 rtol=1e-8
    centroid=stats.origin.+stats.centroid
    chord=2sin(0.45)/0.9          # arc centroid radius factor
    @test collect(centroid)≈[0.5,ρ*chord*cos(0.75),ρ*chord*sin(0.75)] atol=1e-7

    # the box follows the edge's own range, not the whole curve's hull
    half=(occ=:copied_curve,curve=curve,t0=θs[1],t1=θs[21])
    box=M._model_occ_copied_curve_bbox(half)
    @test box[2]≈ρ*cos(θs[21]) atol=1e-8
    @test box[5]≈ρ*cos(θs[1]) atol=1e-8
    @test box[3]≈ρ*sin(θs[1]) atol=1e-8
    @test box[6]≈ρ*sin(θs[21]) atol=1e-8
end
end
