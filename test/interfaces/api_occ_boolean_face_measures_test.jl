module OCCBooleanFaceMeasureChecks
using Test,Tessella
const A=Tessella.API

face_areas(occ,kind)=sort!(Float64[occ.get_mass(2,t)
    for (_,t) in A.model.get_entities(2) if A.model.get_type(2,t)==kind])

# A Boolean kept piece shares its carrier's surface record over a sub-domain
# wire. `occ.get_mass(2, …)` integrated the carrier's whole chart for
# sphere/cylinder/cone pieces — a sphere octant reported the full sphere's
# area, a cylinder wall clipped between two planes its unclipped height.
@testset "OCC area of trimmed quadric Boolean faces" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        ball=occ.add_sphere(0,0,0,1.5)
        occ.intersect([(3,box)],[(3,ball)])
        occ.synchronize()
        @test face_areas(occ,"Sphere")≈[π*1.5^2/2] rtol=1e-6
        @test face_areas(occ,"Plane")≈fill(π*1.5^2/4,3) rtol=1e-12
        sphere=only(t for (_,t) in A.model.get_entities(2)
                    if A.model.get_type(2,t)=="Sphere")
        @test collect(occ.get_center_of_mass(2,sphere))≈fill(0.75,3) atol=1e-5
        # centroidal inertia matrix of the octant: trace = 2·A·(R² − 3·(R/2)²)
        inertia=occ.get_matrix_of_inertia(2,sphere)
        @test inertia[1]+inertia[5]+inertia[9]≈
            2*π*1.5^2/2*(1.5^2-3*0.75^2) rtol=1e-5
    finally
        A.finalize()
    end

    r=0.245
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        cylinder=occ.add_cylinder(0.03,0.02,-1.6,0,0,3.2,r)
        occ.intersect([(3,box)],[(3,cylinder)])
        occ.synchronize()
        @test face_areas(occ,"Cylinder")≈[2π*r*2] rtol=1e-12
    finally
        A.finalize()
    end

    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        cone=occ.add_cone(0.03,0.02,-1.6,0,0,3.2,0.35,0.14)
        occ.intersect([(3,box)],[(3,cone)])
        occ.synchronize()
        radius(z)=0.35+(z+1.6)*(0.14-0.35)/3.2
        r1,r2=radius(-1.0),radius(1.0)
        @test face_areas(occ,"Cone")≈[π*(r1+r2)*hypot(2,r1-r2)] rtol=1e-6
    finally
        A.finalize()
    end

    # whole primitives keep their closed forms
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_sphere(0,0,0,1.5)
        occ.add_cylinder(3,0,0,0,0,2,0.5)
        occ.add_cone(6,0,0,0,0,1,0.5,0.2)
        occ.synchronize()
        @test face_areas(occ,"Sphere")==[4π*1.5^2]
        @test face_areas(occ,"Cylinder")≈[2π*0.5*2] rtol=1e-14
        @test face_areas(occ,"Cone")≈[π*0.7*hypot(1,0.3)] rtol=1e-14
    finally
        A.finalize()
    end
end
end
