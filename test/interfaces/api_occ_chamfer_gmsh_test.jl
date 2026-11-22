# `model.occ.chamfer` / `fillet` on `occ.add_box`, against Gmsh 4.15.2.
#
# The box's three minimum faces carry both a reversed shell use and a
# `reversed` record flag (the OpenCASCADE layout); reading both flipped their
# outward normals back, so a box chamfer failed ("kept faces do not close a
# solid") and a one-edge fillet removed the wrong prism. The Gmsh form of
# `chamfer` takes one surface per curve (the first distance is measured on
# it) and one distance, one per curve or two per curve; the result takes
# the lowest free tags (`_fillet` unbinds the volume first). Masses and
# counts are Gmsh's.
module ApiOccChamferGmshTest
using Test, Tessella

@testset "occ.chamfer / fillet on a box (Gmsh 4.15.2)" begin
    A=Tessella.API; occ=A.model.occ
    A.initialize()
    try
        occ.add_box(0,0,0,1,1,1)
        occ.fillet([1],[1],[0.25]); occ.synchronize()
        @test occ.get_mass(3,1)≈0.98658738521234 rtol=1e-12
        @test length(A.model.get_entities(2))==7
        @test length(A.model.get_entities(1))==15

    finally
        A.finalize()
    end
    A.initialize()
    try
        # examples/boolean/chamfer.geo: two edges meeting at a corner, one
        # distance each, measured on Surface 1
        occ.add_box(0,0,0,1,1,1)
        occ.chamfer([1],[1,2],[1,1],[0.2,0.05]); occ.synchronize()
        @test occ.get_mass(3,1)≈0.978979166666667 rtol=1e-12
        @test [length(A.model.get_entities(d)) for d in 0:3]==[12,18,8,1]
        @test sort!([t for (_,t) in A.model.get_entities(2)])==collect(1:8)
        @test sort!([occ.get_mass(2,t) for t in 1:8])≈sort!([0.76,0.99875,
            0.0583363094478901,0.98,0.281074945521653,0.93875,1.0,0.8]) rtol=1e-12
    finally
        A.finalize()
    end
    A.initialize()
    try
        # two distances: 0.1 on Surface 3 (ymin), 0.3 on the other face
        occ.add_box(0,0,0,1,1,1)
        occ.chamfer([1],[1],[3],[0.1,0.3]); occ.synchronize()
        @test occ.get_mass(3,1)≈0.985 rtol=1e-12
        @test sort!([occ.get_mass(2,t) for (_,t) in occ.get_entities(2)])≈
            sort!([0.7,0.985,0.316227766016838,1.0,0.985,0.9,1.0]) rtol=1e-12
    finally
        A.finalize()
    end
end
end
