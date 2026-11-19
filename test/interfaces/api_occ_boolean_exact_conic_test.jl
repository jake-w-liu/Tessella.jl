module OCCBooleanExactConicChecks
using Test,Tessella
const A=Tessella.API

curve_kinds(occ)=sort!(unique!(String[A.model.get_type(1,t)
    for (_,t) in A.model.get_entities(1) if occ.get_mass(1,t)>0]))

# A plane cuts a sphere in a circle, a cylinder obliquely in an ellipse, and
# two spheres meet in a circle — exact conics even when their trace on the
# curved face's chart is neither a parallel nor a meridian. These sections
# used to be traced numerically (1e-7): the edges were anonymous splines, the
# plane faces' areas were off by up to 1e-3, volumes by 1e-5, and the plane
# faces could not be meshed ("Point[N] is not coplanar with Surface[M]").
# The 3-D conic is kept exact now; only its chart trace is fitted. Expected
# values are closed forms (Gmsh 4.15.2 agrees to its own 1e-8).
@testset "OCC Boolean sections that are exact conics" begin
    # spherical cap: the cutting plane x=0.4 is off the centre
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        ball=occ.add_sphere(0,0,0,1)
        box=occ.add_box(0.4,-2,-2,4,4,4)
        out,_=occ.intersect([(3,ball)],[(3,box)])
        occ.synchronize()
        @test occ.get_mass(3,out[1][2])≈π*0.6^2*(3-0.6)/3 rtol=1e-9
        @test curve_kinds(occ)==["Circle"]
        areas=Dict{String,Float64}()
        for (_,t) in A.model.get_entities(2)
            kind=A.model.get_type(2,t)
            areas[kind]=get(areas,kind,0.0)+occ.get_mass(2,t)
        end
        @test areas["Plane"]≈π*(1-0.4^2) rtol=1e-12
        @test areas["Sphere"]≈2π*0.6 rtol=1e-9
    finally
        A.finalize()
    end

    # cylinder cut by an oblique half-space: an ellipse through the axis
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        cylinder=occ.add_cylinder(0,0,-1,0,0,2,0.5)
        box=occ.add_box(-2,-2,0,4,4,2)
        occ.rotate([(3,box)],0,0,0,0,1,0,0.5)
        out,_=occ.cut([(3,cylinder)],[(3,box)])
        occ.synchronize()
        @test occ.get_mass(3,out[1][2])≈π*0.25 rtol=1e-10
        @test "Ellipse" in curve_kinds(occ)
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,3,3,1)
    finally
        A.finalize()
    end

    # lens of two unit spheres 1.2 apart: V = π(4+d)(2−d)²/12
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        left=occ.add_sphere(0,0,0,1)
        right=occ.add_sphere(1.2,0,0,1)
        out,_=occ.intersect([(3,left)],[(3,right)])
        occ.synchronize()
        @test occ.get_mass(3,out[1][2])≈π*5.2*0.8^2/12 rtol=1e-9
        @test curve_kinds(occ)==["Circle"]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,3,3,1)
    finally
        A.finalize()
    end
end
end
