module OCCBooleanBoundingBoxChecks
using Test,Tessella
const A=Tessella.API

model_box(dim,tag)=collect(A.model.get_bounding_box(dim,tag))
kernel_box(dim,tag)=collect(A.model.occ.get_bounding_box(dim,tag))

# Bounding boxes of Boolean results with spherical faces. The synchronized
# model boxed a shelled volume by its boundary curves only (a spherical cap
# lost its bulge; a re-materialized sphere collapsed to its seam meridian),
# while both paths boxed any sphere piece by the whole ball and a traced
# section edge by the hull of the whole section curve.
@testset "OCC bounding boxes of Boolean results with sphere faces" begin
    # sphere centred on the box face x=1: the face cuts it along a great
    # circle through both poles, leaving a hemispherical cap on each side
    for (name,expected) in (
            ("fuse",[-1.,-1.,-1.,1.45,1.,1.]),
            ("intersect",[0.55,-0.4,-0.35,1.,0.5,0.55]),
            ("cut",[-1.,-1.,-1.,1.,1.,1.]))
        @testset "sphere through a box face: $name" begin
            A.initialize()
            try
                A.option("General.Terminal",0)
                occ=A.model.occ
                box=occ.add_box(-1,-1,-1,2,2,2)
                ball=occ.add_sphere(1,0.05,0.1,0.45)
                out,_=getfield(occ,Symbol(name))([(3,box)],[(3,ball)])
                occ.synchronize()
                tag=out[1][2]
                @test model_box(3,tag)≈expected atol=1e-6
                @test kernel_box(3,tag)≈expected atol=1e-6
            finally
                A.finalize()
            end
        end
    end

    # sphere octant: every face box, both paths
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        ball=occ.add_sphere(0,0,0,1.5)
        occ.intersect([(3,box)],[(3,ball)])
        occ.synchronize()
        @test model_box(3,1)≈[0,0,0,1.5,1.5,1.5] atol=1e-9
        for (_,t) in A.model.get_entities(2)
            A.model.get_type(2,t)=="Sphere" || continue
            @test model_box(2,t)≈[0,0,0,1.5,1.5,1.5] atol=1e-9
            @test kernel_box(2,t)≈[0,0,0,1.5,1.5,1.5] atol=1e-9
        end
    finally
        A.finalize()
    end

    # whole primitives are unchanged
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_sphere(0.3,-0.2,0.1,0.6)
        occ.synchronize()
        expected=[-0.3,-0.8,-0.5,0.9,0.4,0.7]
        @test model_box(3,1)≈expected atol=1e-12
        @test model_box(2,1)≈expected atol=1e-12
        @test kernel_box(2,1)≈expected atol=1e-12
    finally
        A.finalize()
    end
end
end
