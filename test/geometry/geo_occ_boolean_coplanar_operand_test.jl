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
end
