module OCCBooleanContainmentChecks
using Test,Tessella
const A=Tessella.API

# A tool lying strictly inside the object cuts none of the object's faces,
# yet the result is a different solid: the tool's boundary becomes an inner
# (cavity) shell. Expected counts are Gmsh 4.15.2's; volumes and areas are
# closed forms.
@testset "OCC Booleans with a tool fully inside the object" begin
    sphere_volume=4π/3*0.5^3
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        ball=occ.add_sphere(1,1,1,0.5)
        out,images=occ.cut([(3,box)],[(3,ball)])
        @test out==[(3,1)]
        @test images==[[(3,1)],Tuple{Int,Int}[]]
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(10,15,7,1)
        @test occ.get_mass(3,1)≈8-sphere_volume rtol=1e-8
        @test collect(occ.get_center_of_mass(3,1))≈[1.0,1.0,1.0] atol=1e-7
        faces=A.model.get_boundary([(3,1)],false,false,false)
        @test length(faces)==7
        @test sum(occ.get_mass(2,t) for (_,t) in faces)≈24+π rtol=1e-12
    finally
        A.finalize()
    end

    # every primitive pair: the cut removes exactly the tool's volume
    tools=(
        (occ->occ.add_box(-0.4,-0.4,-0.4,0.8,0.8,0.8),0.512,12),
        (occ->occ.add_cylinder(0,0,-0.28,0,0,0.56,0.196),π*0.196^2*0.56,9),
        (occ->occ.add_cone(0,0,-0.4,0,0,0.8,0.4,0.16),
            π*0.8*(0.4^2+0.4*0.16+0.16^2)/3,9),
        (occ->occ.add_wedge(-0.4,-0.4,-0.4,0.8,0.8,0.8),0.256,11),
    )
    for (make,volume,nfaces) in tools
        A.initialize()
        try
            A.option("General.Terminal",0)
            occ=A.model.occ
            box=occ.add_box(-1,-1,-1,2,2,2)
            tool=make(occ)
            out,_=occ.cut([(3,box)],[(3,tool)])
            @test out==[(3,1)]
            occ.synchronize()
            @test occ.get_mass(3,1)≈8-volume rtol=1e-8
            @test length(A.model.get_boundary([(3,1)],false,false,false))==
                nfaces
        finally
            A.finalize()
        end
    end

    # the other operations on the same placement
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        ball=occ.add_sphere(1,1,1,0.5)
        out,images=occ.intersect([(3,box)],[(3,ball)])
        @test out==[(3,2)]
        @test images==[[(3,2)],[(3,2)]]
        @test occ.get_mass(3,2)≈sphere_volume rtol=1e-12
    finally
        A.finalize()
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        ball=occ.add_sphere(1,1,1,0.5)
        out,images=occ.fragment([(3,box)],[(3,ball)])
        @test out==[(3,3),(3,2)]
        @test images==[[(3,3),(3,2)],[(3,2)]]
        @test occ.get_mass(3,2)≈sphere_volume rtol=1e-12
        @test occ.get_mass(3,3)≈8-sphere_volume rtol=1e-8
    finally
        A.finalize()
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        ball=occ.add_sphere(1,1,1,0.5)
        out,_=occ.fuse([(3,box)],[(3,ball)])
        @test out==[(3,1)]
        @test occ.get_mass(3,1)≈8 rtol=1e-12
        occ.synchronize()
        @test length(A.model.get_boundary([(3,1)],false,false,false))==6
    finally
        A.finalize()
    end
end
end
