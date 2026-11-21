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
        # The two caps are ellipses with semi-axes 1 and 0.5. They are
        # degree-1 patches like Gmsh's, not planes whose frame took the
        # stretch: such a frame is not orthonormal, and the cap of a
        # dilated cone was meshed with half of Gmsh's triangles.
        @test all(A.model.get_type(2,t)=="BSpline surface"
                  for (_,t) in A.model.get_entities(2))
        # Surfaces 2 and 3 are the caps (the wall, Surface 1, is a rational
        # patch whose quadrature takes 18 s when interpreted)
        @test [occ.get_mass(2,t) for t in (2,3)]≈[π/2,π/2] rtol=1e-10
        if Base.JLOptions().compile_enabled==1
            @test occ.get_mass(3,rod)≈π/2 rtol=1e-10
        end
    finally
        A.finalize()
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,1,2,3)
        occ.dilate([(3,box)],0,0,0,2.0,1.0,0.5)
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(8,12,6,1)
        @test all(A.model.get_type(2,t)=="BSpline surface" for t in 1:6)
        @test sort!([occ.get_mass(2,t) for t in 1:6])≈[3,3,3,3,4,4] rtol=1e-12
        @test occ.get_mass(3,box)≈6 rtol=1e-12
        if Base.JLOptions().compile_enabled==1
            # Gmsh 4.15.2 at 0.2: 969 nodes, 1266 triangles, 3856
            # tetrahedra; the stretched plane frames gave 1098 triangles
            A.option("Mesh.MeshSizeMax",0.2)
            A.mesh.generate(2)
            triangles=sum(length(A.mesh.get_elements(2,t)[2][1]) for t in 1:6)
            @test 1200<triangles<1340
        end
    finally
        A.finalize()
    end
end
end
