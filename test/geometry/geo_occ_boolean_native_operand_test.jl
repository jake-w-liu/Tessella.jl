module GeoOCCBooleanNativeOperandChecks
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

# `getMass` of a `.geo` volume: Tessella binds `.geo` entities into the OCC
# kernel lazily, so bind it first.
function mass(t)
    Tessella.API._with_model(m->Tessella.GeoExec._geo_occ_bind_entities!(m,
        [(3,t)],Tessella.GeoExec._geo_occ_prepare_kernel!(m),"test"))
    return A.model.occ.get_mass(3,t)
end

# Gmsh's examples/boolean/transform.geo, cut down: the translated box passes
# through the OCC kernel, the fragment then rebuilds it natively, and
# `Duplicata` must copy the fragment piece, not the kernel's earlier box.
# Gmsh 4.15.2: 36/57/27/4 entities, masses 6.5, 1.5, 6.5, 6.5.
@testset "Duplicata of a fragment piece of a transformed box" begin
    opened("""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 2,2,2};
        Box(2) = {1,1,1, 2,2,2};
        Translate{0.5,0,0}{ Volume{1}; }
        f() = BooleanFragments { Volume{1,2}; Delete; }{};
        Translate{4,0,0}{ Duplicata{ Volume{1}; } }
        """) do
        @test volumes()==[1,2,3,4]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(36,57,27,4)
        @test [mass(t) for t in 1:4]≈[6.5,1.5,6.5,6.5] rtol=1e-12
    end
end

# Gmsh's examples/boolean/simple4.geo with a square rod for the cylinder: a
# pyramid of Plane Surfaces whose Curve Loops wind both ways, all listed
# positive in its Surface Loop, fragmented with a rod through its apex.
# OpenCASCADE sews the shell into a solid of volume 1/3. Gmsh 4.15.2:
# 21/40/24/4 entities, masses 0.29866667, 0.03466667, 0.02533333, 0.02.
@testset "Fragment of an unoriented Plane Surface shell" begin
    opened("""
        SetFactory("OpenCASCADE");
        Point(1) = {0,0,0}; Point(2) = {1,0,0}; Point(3) = {1,1,0};
        Point(4) = {0,1,0}; Point(5) = {0.5,0.5,1};
        Line(1) = {1,2}; Line(2) = {3,2}; Line(3) = {4,3}; Line(4) = {1,4};
        Line(5) = {5,1}; Line(6) = {5,2}; Line(7) = {5,3}; Line(8) = {5,4};
        Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {1,5,6};
        Curve Loop(3) = {2,6,7}; Curve Loop(4) = {3,7,8};
        Curve Loop(5) = {4,5,8};
        Plane Surface(1) = {1}; Plane Surface(2) = {2}; Plane Surface(3) = {3};
        Plane Surface(4) = {4}; Plane Surface(5) = {5};
        Surface Loop(1) = {1,2,3,4,5};
        Volume(1) = {1};
        Box(2) = {0.4,0.4,-0.5, 0.2,0.2,2};
        BooleanFragments{ Volume{1,2}; Delete; }{}
        """) do
        @test volumes()==[1,2,3,4]
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(21,40,24,4)
        vols=[mass(t) for t in 1:4]
        @test sort(vols)≈sort([1/3-0.104/3,0.104/3,0.076/3,0.02]) rtol=1e-12
        # Gmsh orders the two rod pieces outside the pyramid top, bottom;
        # Tessella's split-face order (`_brep_regions_from_first_vertex`)
        # puts the bottom first here
        @test_broken vols≈[1/3-0.104/3,0.104/3,0.076/3,0.02]
    end
end
end
