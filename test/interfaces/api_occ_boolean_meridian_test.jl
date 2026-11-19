module OCCBooleanMeridianChecks
using Test,Tessella
const A=Tessella.API

# Planes through a sphere's polar axis cut it along great circles made of two
# meridian halves that share both endpoints (the poles).
@testset "OCC Boolean sections along sphere meridians" begin
    # Hemisphere: a box face through both poles. The half whose meridian
    # trace is valid one turn up (θ∈[2π,3π]) was evaluated on the first turn,
    # outside the sphere's chart, and lost; the two halves were then taken
    # for one edge because they share their endpoints. `occ.intersect` threw
    # "duplicate reverse pcurve on Curve[19]". Gmsh: volume 2π/3.
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        ball=occ.add_sphere(0,0,0,1)
        box=occ.add_box(0,-2,-2,4,4,4)
        out,_=occ.intersect([(3,ball)],[(3,box)])
        @test length(out)==1
        @test occ.get_mass(3,out[1][2])≈2π/3 rtol=1e-8
        # the cap is split at the sphere's seam, as in Gmsh: (2,7,3,1)
        @test Tuple(length(occ.get_entities(d)) for d in 0:3)==(2,7,3,1)
        @test sum(occ.get_mass(2,t) for (_,t) in occ.get_entities(2))≈3π rtol=1e-12
    finally
        A.finalize()
    end

    # Sphere centred on a box corner, the box in the octant away from the
    # seam: the plane y=1 cuts along the seam's great circle, but only the
    # half on the seam's side is the seam edge. The other half was dropped as
    # "coincident with an operand edge" and the cut's shell did not close
    # ("Curve[23] must occur once on each of two distinct surfaces").
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,1,1)
        ball=occ.add_sphere(2,1,1,0.5)
        out,_=occ.cut([(3,box)],[(3,ball)])
        @test out==[(3,1)]
        @test occ.get_mass(3,1)≈2-π*0.5^3/6 rtol=1e-8
        @test length(occ.get_entities(2))==7
    finally
        A.finalize()
    end
end
end
