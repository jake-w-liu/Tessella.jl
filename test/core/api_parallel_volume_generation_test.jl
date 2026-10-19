# PAR01 slice — parallel mesh *generation* of independent Regions under
# `General.NumThreads`.  `_higher_volume_phase!` detaches runs of eligible
# singleton Regions onto per-unit transaction copies under `Threads.@threads`
# and merges the receipts serially in Region order; stale or unsupported slots
# re-run the identical live serial call, so the produced mesh — its CRC, its
# serialized MSH bytes and its validation — must be identical for any
# requested worker count.  Run standalone (e.g.
# `julia -t4 test/core/api_parallel_volume_generation_test.jl`); the default
# suite does not register this file.  Under `-t1` both arms take the serial
# path and the checks still hold.
#
# The eligible set is deliberately narrow: implicit/analytic factory solids
# and factory-encoded shells without embedded constraints (the paths that
# provably never consume the shared `Mesh.RandomSeed` ledger — `Box`,
# `Sphere`, `Cylinder`, `Cone`).  Explicit PLC volumes whose boundary Surfaces
# all resolve to actual authority products take the native-initial
# Delaunay fill, which *does* thread the shared RNG; they keep the serial arm,
# as do transfinite, swept, compound and periodic Regions.

using Test
using Tessella
using Tessella.Elements: mixed_crc, write_mixed_msh

A=Tessella.API

# Two disjoint OCC Boxes — both factory-encoded, no embeds: both Regions are
# eligible and the run really fills under Threads.@threads.
const TWOBOX="""
SetFactory("OpenCASCADE");
Box(1)={0,0,0, 1,1,1};
Box(2)={2.5,0,0, 1,1.5,1};
"""

# Three disjoint Boxes — a run of three eligible Regions.
const THREEBOX="""
SetFactory("OpenCASCADE");
Box(1)={0,0,0, 1,1,1};
Box(2)={2.5,0,0, 1,1.5,1};
Box(3)={0,2.5,0, 1,1,1.5};
"""

# Two disjoint explicit PLC volumes — every boundary Surface resolves to an
# actual authority product, so each takes the native-initial fill that
# consumes the shared RNG ledger.  Both keep the serial arm; the mesh must
# still be byte-identical under -t1/-t4.
const TWOVOL="""
Point(1)={0,0,0,0.35};Point(2)={1,0,0,0.35};Point(3)={1,1,0,0.35};Point(4)={0,1,0,0.35};
Point(5)={0,0,1,0.35};Point(6)={1,0,1,0.35};Point(7)={1,1,1,0.35};Point(8)={0,1,1,0.35};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};
Curve Loop(3)={1,10,-5,-9};Curve Loop(4)={2,11,-6,-10};
Curve Loop(5)={3,12,-7,-11};Curve Loop(6)={4,9,-8,-12};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
Plane Surface(4)={4};Plane Surface(5)={5};Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};
Point(21)={2.5,0,0,0.35};Point(22)={3.5,0,0,0.35};Point(23)={3.5,1,0,0.35};Point(24)={2.5,1,0,0.35};
Point(25)={2.5,0,1,0.35};Point(26)={3.5,0,1,0.35};Point(27)={3.5,1,1,0.35};Point(28)={2.5,1,1,0.35};
Line(21)={21,22};Line(22)={22,23};Line(23)={23,24};Line(24)={24,21};
Line(25)={25,26};Line(26)={26,27};Line(27)={27,28};Line(28)={28,25};
Line(29)={21,25};Line(30)={22,26};Line(31)={23,27};Line(32)={24,28};
Curve Loop(21)={21,22,23,24};Curve Loop(22)={25,26,27,28};
Curve Loop(23)={21,30,-25,-29};Curve Loop(24)={22,31,-26,-30};
Curve Loop(25)={23,32,-27,-31};Curve Loop(26)={24,29,-28,-32};
Plane Surface(21)={21};Plane Surface(22)={22};Plane Surface(23)={23};
Plane Surface(24)={24};Plane Surface(25)={25};Plane Surface(26)={26};
Surface Loop(2)={21,22,23,24,25,26};Volume(2)={2};
"""

# Two adjacent explicit volumes sharing one Plane Surface — the shared shell
# keeps both Regions on the serial arm.
const SHAREDVOL="""
Point(1)={0,0,0,0.35};Point(2)={1,0,0,0.35};Point(3)={1,1,0,0.35};Point(4)={0,1,0,0.35};
Point(5)={0,0,1,0.35};Point(6)={1,0,1,0.35};Point(7)={1,1,1,0.35};Point(8)={0,1,1,0.35};
Point(9)={2,0,0,0.35};Point(10)={2,1,0,0.35};Point(11)={2,1,1,0.35};Point(12)={2,0,1,0.35};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Line(17)={2,9};Line(18)={9,10};Line(19)={10,3};Line(20)={9,12};
Line(21)={10,11};Line(22)={7,11};Line(23)={6,12};Line(24)={12,11};
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};
Curve Loop(3)={1,10,-5,-9};Curve Loop(4)={2,11,-6,-10};
Curve Loop(5)={3,12,-7,-11};Curve Loop(6)={4,9,-8,-12};
Curve Loop(7)={18,21,-24,-20};Curve Loop(8)={17,20,-23,-10};
Curve Loop(9)={11,22,-21,19};Curve Loop(10)={2,-19,-18,-17};
Curve Loop(11)={23,24,-22,-6};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
Plane Surface(4)={4};Plane Surface(5)={5};Plane Surface(6)={6};
Plane Surface(7)={7};Plane Surface(8)={8};Plane Surface(9)={9};
Plane Surface(10)={10};Plane Surface(11)={11};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};
Surface Loop(2)={7,8,9,10,11,-4};Volume(2)={2};
"""

# Mixed traversal: an explicit Region (serial arm) followed by two factory
# Boxes (a detached run) in one phase — run splitting, fill, replayed merge
# and serial re-run coexist.
const VOLBOX=first(split(TWOVOL,"Surface Loop(2)"))*
"""Surface Loop(2)={21,22,23,24,25,26};Volume(2)={2};
SetFactory("OpenCASCADE");
Box(3)={5,0,0, 1,1,1};
Box(4)={5,2.5,0, 1,1,1};
"""

const TFVOL="""
Point(1)={0,0,0,0.35};Point(2)={1,0,0,0.35};Point(3)={1,1,0,0.35};Point(4)={0,1,0,0.35};
Point(5)={0,0,1,0.35};Point(6)={1,0,1,0.35};Point(7)={1,1,1,0.35};Point(8)={0,1,1,0.35};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};
Curve Loop(3)={1,10,-5,-9};Curve Loop(4)={2,11,-6,-10};
Curve Loop(5)={3,12,-7,-11};Curve Loop(6)={4,9,-8,-12};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
Plane Surface(4)={4};Plane Surface(5)={5};Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};
Transfinite Curve{1,2,3,4,5,6,7,8,9,10,11,12}=3;
Transfinite Surface{1,2,3,4,5,6};
Transfinite Volume{1}={1,2,3,4,5,6,7,8};
SetFactory("OpenCASCADE");
Box(2)={2.5,0,0, 1,1,1};
"""

function _par_volume_generate(src,mesh_dim,workers;extra_geo="")
    path,tempio=mktemp()
    write(tempio,src*extra_geo*"\nMesh $mesh_dim;\n");close(tempio)
    A.finalize();A.initialize()
    A.option("General.Terminal",0)
    A.option("General.NumThreads",workers)
    ex=execute_geo(path;mesh_dim)
    rm(path)
    actual=ex.model.meshing.mesh_authority.mesh
    msh=joinpath(mktempdir(),"out.msh")
    write_mixed_msh(msh,actual;version=4.1)
    result=(mixed_crc(actual),read(msh),validate(actual).ok)
    A.finalize()
    return result
end

@testset "parallel detached Region generation stays deterministic" begin
    cases=(
        ("two factory boxes Mesh3",TWOBOX,3,""),
        ("three factory boxes Mesh3",THREEBOX,3,""),
        ("two explicit volumes Mesh3",TWOVOL,3,""),
        ("shared-boundary volumes Mesh3",SHAREDVOL,3,""),
        ("explicit+factory split run Mesh3",VOLBOX,3,""),
        ("transfinite+factory Mesh3",TFVOL,3,""),
        ("two boxes re-mesh authority Mesh3",TWOBOX,3,"Mesh 3;\n"),
        ("two boxes Mesh2 only",TWOBOX,2,""),
    )
    try
        for (name,src,dim,extra) in cases
            scrc,sbytes,sok=_par_volume_generate(src,dim,1;extra_geo=extra)
            tcrc,tbytes,tok=_par_volume_generate(src,dim,4;extra_geo=extra)
            @test sok
            @test tok
            @test scrc==tcrc
            @test sbytes==tbytes
            scrc==tcrc && sbytes==tbytes && sok && tok ||
                println("PARALLEL-DIFFERS: $name")
        end
    finally
        A.finalize()
    end
end
