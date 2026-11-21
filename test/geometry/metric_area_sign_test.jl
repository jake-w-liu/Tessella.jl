module MetricAreaSignChecks
using Test, Tessella
const M=Tessella.Model

# The metric producer's admissibility search asks, for every candidate split
# of an edge, whether the two child triangles have a finite positive area. It
# asked `triangle_area`, whose floating-point certificate fails on razor
# cells and then evaluates the area in `Rational{BigInt}`; on t17 a third of
# the producer's time went there. `_metric_area_sign` answers the two
# questions only: an area whose cross product clears its rounding bound is
# positive, anything closer is decided by the exact `orient2` on the three
# coordinate planes (the area is zero exactly when the rows are collinear),
# and only an input outside Float64's safe range still takes the exact area.
@testset "area sign of a child triangle agrees with the exact area" begin
    sign_of(a,b,c)=(area=M.triangle_area(a,b,c);(isfinite(area),area>0))
    ordinary=((0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.5))
    @test M._metric_area_sign(ordinary...)==(true,true)==sign_of(ordinary...)
    # Three collinear rows whose float cross product is not zero.
    a=(0.1,0.2,0.3);d=(0.7,-0.3,0.9)
    collinear=(a,ntuple(i->a[i]+d[i],3),ntuple(i->a[i]+2d[i],3))
    exact_line=((0.0,0.0,0.0),(1.0,2.0,3.0),(2.0,4.0,6.0))
    @test M._metric_area_sign(exact_line...)==(true,false)==sign_of(exact_line...)
    @test M._metric_area_sign(collinear...)==sign_of(collinear...)
    # Razor cells: the apex is a few ulps off the base line.
    for k in 0:40
        offset=k*eps(1.0)
        razor=((0.0,0.0,0.0),(1.0,1.0,1.0),(0.5,0.5,0.5+offset))
        @test M._metric_area_sign(razor...)==sign_of(razor...)
        tilted=((0.3,0.1,0.7),(1.3,1.1,1.7),(0.8+offset,0.6,1.2-offset))
        @test M._metric_area_sign(tilted...)==sign_of(tilted...)
    end
    @test M._metric_area_sign((0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,NaN,0.0))==(false,false)
    @test M._metric_area_sign((0.0,0.0,0.0),(Inf,0.0,0.0),(0.0,1.0,0.0))==(false,false)
    huge=((0.0,0.0,0.0),(1e200,0.0,0.0),(0.0,1e200,0.0))
    @test M._metric_area_sign(huge...)==sign_of(huge...)
    if Base.JLOptions().compile_enabled==1
        M._metric_area_sign(ordinary...)
        @test (@allocated M._metric_area_sign(ordinary...))==0
        razor=((0.0,0.0,0.0),(1.0,1.0,1.0),(0.5,0.5,0.5+4eps(1.0)))
        M._metric_area_sign(razor...);M.triangle_area(razor...)
        # The exact predicates cost a few words; the exact area costs kilobytes.
        @test (@allocated M._metric_area_sign(razor...))<=512
        @test (@allocated M.triangle_area(razor...))>2048
    end
end

end
