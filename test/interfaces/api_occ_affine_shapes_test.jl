module OCCAffineShapesChecks
using Test,Tessella
const A=Tessella.API

@testset "OCC dilate with unequal factors" begin
    # A dilation that is no similarity was refused for anything round ("it
    # requires an isotropic linear part"). The shape is now converted to
    # rational B-splines first, as OpenCASCADE does, and those take any
    # affine map exactly through their poles. Counts are Gmsh 4.15.2's; the
    # volumes are closed forms (Gmsh's own are 4.1869583 for the ellipsoid
    # and 1.5843294 for the cylinder: its conversion is approximate).
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        ball=occ.add_sphere(0,0,0,1)
        occ.dilate([(3,ball)],0,0,0,2.0,1.0,0.5)
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,3,1,1)
        @test A.model.get_type(2,1)=="BSpline surface"
        if Base.JLOptions().compile_enabled==1
            @test occ.get_mass(3,ball)≈4π/3 rtol=1e-10
        end
    finally
        A.finalize()
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        rod=occ.add_cylinder(0,0,0,0,0,2,0.5)
        occ.dilate([(3,rod)],0,0,0,2.0,1.0,0.5)
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,3,3,1)
        # the two caps are ellipses with semi-axes 1 and 0.5
        caps=[occ.get_mass(2,t) for (_,t) in A.model.get_entities(2)
              if A.model.get_type(2,t)=="Plane"]
        @test caps≈[π/2,π/2] rtol=1e-10
        if Base.JLOptions().compile_enabled==1
            @test occ.get_mass(3,rod)≈π/2 rtol=1e-10
        end
    finally
        A.finalize()
    end
end
end
