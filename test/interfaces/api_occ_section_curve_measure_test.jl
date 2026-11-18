module OCCSectionCurveMeasureChecks
using Test,Tessella
const A=Tessella.API

# A box face cutting a sphere off-center leaves a section circle that is
# neither a parallel nor a meridian of the sphere, so it is stored as a traced
# curve S(pc(t)). Its length quadrature never converged: the panels ignored
# the trace's own knots, and the centroid component normal to the cutting
# plane — pure evaluation noise — was held to a tolerance below rounding.
# `occ.get_mass(1, …)` threw "quadrature did not converge".
@testset "OCC length of traced Boolean section curves" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        ball=occ.add_sphere(0.8,0.55,0.35,0.6)
        out,_=occ.intersect([(3,box)],[(3,ball)])
        @test length(out)==1
        occ.synchronize()
        lengths=Float64[occ.get_mass(1,t) for (_,t) in A.model.get_entities(1)]
        @test all(isfinite,lengths)
        # the cap circles on x=1 and y=1: radii from the plane offsets
        rx=sqrt(0.6^2-0.2^2);ry=sqrt(0.6^2-0.45^2)
        # they cross at the two points of the edge x=y=1 inside the sphere,
        # so each plane keeps the arc on the box side of the other plane
        arc(r,d)=2r*acos(clamp(d/r,-1,1))
        expected=(2π*rx-arc(rx,0.45))+(2π*ry-arc(ry,0.2))
        # traced pieces only: box edges are Lines, sphere seam pieces Circles
        curved=sum(l for (l,(_,t)) in zip(lengths,A.model.get_entities(1))
                   if !(A.model.get_type(1,t) in ("Line","Circle")))
        @test curved≈expected rtol=1e-6
        for (_,t) in A.model.get_entities(1)
            center=collect(occ.get_center_of_mass(1,t))
            @test all(isfinite,center)
        end
    finally
        A.finalize()
    end
end
end
