module GeoOCCBooleanCoplanarOperandChecks
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

# A `.geo` Boolean whose operand is itself a Boolean result resting on a
# coplanar box (Gmsh's examples/boolean/extend_field.geo). The operand
# snapshot for the analytic mesh path cannot be formed there: the
# triangle-soup Boolean refuses coplanar overlap. That must not stop the
# statement, since the conforming mesh never uses it. Gmsh 4.15.2: volumes 3, 10, 11.
@testset "Boolean of a Boolean operand on a coplanar box" begin
    opened("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 1,1,0.5};
        Sphere(2) = {1,1,0.5,0.4};
        BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
        Box(10) = {0.3,0.3,0.5, 0.1,0.1,0.05};
        BooleanFragments{ Volume{10}; Delete; }{ Volume{3}; Delete; }
        """) do
        @test volumes()==[3,10]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(18,28,13,2)
    end
end

# t18's selection: pieces of a box fragmented by a sphere at its corner are
# picked by `Volume In BoundingBox`, which needs each piece's own bounds
# (the sphere's carrier box reaches far outside the box), and the sphere
# piece below the pole keeps an evaluable pole edge. Gmsh 4.15.2 with
# `Geometry.OCCBoundsUseStl` selects volumes 1 and 2.
@testset "Bounds of fragment pieces of a box and a sphere" begin
    opened("""
        SetFactory("OpenCASCADE");
        Box(10) = {2,0,0, 1,1,1};
        Sphere(11) = {1.7,0,0, 0.35};
        v() = BooleanFragments { Volume{10}; Delete; }{ Volume{11}; Delete; };
        Geometry.OCCBoundsUseStl = 1;
        eps = 1e-3;
        vin() = Volume In BoundingBox {2-eps,-eps,-eps, 3+eps,1+eps,1+eps};
        Physical Volume(1) = {vin()};
        """) do
        @test volumes()==[1,2,3]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(13,22,11,3)
        # Gmsh numbers the box's corner piece 1 and the rest of the box 2
        boxes=Dict(t=>A.model.get_bounding_box(3,t) for t in 1:2)
        @test_broken boxes[1][4]<boxes[2][4]
        @test sort(A.model.get_entities_for_physical_group(3,1))==[1,2]
        small=argmin(t->boxes[t][4],[1,2])
        bb=boxes[small]
        @test all(isapprox.(boxes[3-small],(2.0,0.0,0.0,3.0,1.0,1.0);atol=1e-12))
        @test all(isapprox.(bb,(2.0,0.0,0.0,2.05,sqrt(0.35^2-0.3^2),sqrt(0.35^2-0.3^2));atol=1e-12))
        for (_,c) in A.model.get_entities(1)
            @test all(isfinite,A.model.get_value(1,c,[0.5]))
        end
    end
end
end
