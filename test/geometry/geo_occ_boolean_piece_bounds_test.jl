module GeoOCCBooleanPieceBoundsChecks
using Test,Tessella
const A=Tessella.API

function opened(f,text)
    path=tempname()*".geo"
    write(path,text)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        f()
    finally
        A.finalize()
        rm(path;force=true)
    end
end

volumes()=[t for (_,t) in A.model.get_entities(3)]

# t18's selection: pieces of a box fragmented by a sphere over its edge are
# picked by `Volume In BoundingBox`, which needs each piece's own bounds
# (the sphere's carrier box reaches far outside the box), and the sphere
# pieces keep evaluable pole and seam edges. The sphere sits on the box's
# y=0 face and crosses its x=2 face: two cutting planes, as t18's corner
# case but cheaper (that one took 22-37 s interpreted). Gmsh 4.15.2 with
# `Geometry.OCCBoundsUseStl`: 12/20/10/3 entities, the box rest 1, the cap
# inside the box 2 with bounds (2,0,0.3197,2.05,0.1803,0.6803), volumes 1
# and 2 selected.
@testset "Bounds of fragment pieces of a box and a sphere" begin
    opened("""
        SetFactory("OpenCASCADE");
        Box(10) = {2,0,0, 1,1,1};
        Sphere(11) = {1.7,0,0.5, 0.35};
        v() = BooleanFragments { Volume{10}; Delete; }{ Volume{11}; Delete; };
        Geometry.OCCBoundsUseStl = 1;
        eps = 1e-3;
        vin() = Volume In BoundingBox {2-eps,-eps,-eps, 3+eps,1+eps,1+eps};
        Physical Volume(1) = {vin()};
        """) do
        @test volumes()==[1,2,3]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(12,20,10,3)
        @test sort(A.model.get_entities_for_physical_group(3,1))==[1,2]
        r=sqrt(0.35^2-0.3^2)
        @test all(isapprox.(A.model.get_bounding_box(3,1),
                            (2.0,0.0,0.0,3.0,1.0,1.0);atol=1e-12))
        @test all(isapprox.(A.model.get_bounding_box(3,2),
                            (2.0,0.0,0.5-r,2.05,r,0.5+r);atol=1e-12))
        @test all(isapprox.(A.model.get_bounding_box(3,3),
                            (1.35,-0.35,0.15,2.05,0.35,0.85);atol=1e-12))
        for (_,c) in A.model.get_entities(1)
            @test all(isfinite,A.model.get_value(1,c,[0.5]))
        end
    end
end
end
