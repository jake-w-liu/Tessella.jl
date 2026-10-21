# PAR01 slice — parallel mesh *generation* of independent Curves under
# `General.NumThreads`. `_dim01_curve_phase!` detaches runs of eligible
# singleton Curves onto per-unit transaction copies under `Threads.@threads`
# and merges the receipts serially in tag order; stale or unsupported slots
# re-run the identical live serial call and emission still mints every public
# UID on the live assembly, so the produced mesh — its CRC, its serialized
# MSH bytes and its validation — must be identical for any requested worker
# count. Run standalone (e.g.
# `julia -t4 test/core/api_parallel_curve_generation_test.jl`); the default
# suite does not register this file. Under `-t1` both arms take the serial
# path and the checks still hold.

using Test
using Tessella
using Tessella.Elements: mixed_crc, write_mixed_msh

A=Tessella.API

# Disjoint Lines with different per-vertex mesh sizes plus a Circle and a
# Spline: every Curve resolves to its own grading — a genuine detached run
# (and the native `_ModelNativeCurveGrade` receipt for the arc/spline arms).
const INDEPENDENT="""
Point(1)={0,0,0,0.10};Point(2)={1,0,0,0.30};
Point(3)={2,0,0,0.05};Point(4)={3,0,0,0.50};
Point(5)={4,0,0,0.20};Point(6)={5,0,0,0.20};
Point(7)={6,0,0,0.40};Point(8)={6,1,0,0.40};Point(9)={7,0.5,0,0.40};
Point(10)={8,0,0,0.15};Point(11)={8,1,0,0.15};Point(12)={9,1,0,0.15};Point(13)={9,0,0,0.15};
Line(1)={1,2};Line(2)={3,4};Line(3)={5,6};
Circle(4)={7,9,8};Spline(5)={10,11,12,13};
"""

# The wireframe of one box: twelve Lines sharing endpoints — shared boundary
# Points stay eligible (their state is immutable through the Curve phase) and
# `Mesh 1` emits a real detached run over connected geometry.
const BOXWIRE="""
Point(1)={0,0,0,0.35};Point(2)={1,0,0,0.35};Point(3)={1,1,0,0.35};Point(4)={0,1,0,0.35};
Point(5)={0,0,1,0.35};Point(6)={1,0,1,0.35};Point(7)={1,1,1,0.35};Point(8)={0,1,1,0.35};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
"""

# The same wireframe closed into a solid box: the Curve phase runs inside the
# full `Mesh 2`/`Mesh 3` pipeline before the Surface/Region arms consume its
# stored discretizations.
const BOX=BOXWIRE*"""
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};
Curve Loop(3)={1,10,-5,-9};Curve Loop(4)={2,11,-6,-10};
Curve Loop(5)={3,12,-7,-11};Curve Loop(6)={4,9,-8,-12};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
Plane Surface(4)={4};Plane Surface(5)={5};Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};
"""

# Serial-fallback cases mixed into eligible runs: Transfinite Curves, a
# Periodic slave/master pair, a `Degenerated` keep and a `Compound Curve`
# membership all keep the serial arm while the plain Lines around them still
# detach — run splitting, fill, replayed merge and the deferred periodic copy
# coexist in one pass.
const MIXED=BOXWIRE*"""
Point(21)={3,0,0,0.3};Point(22)={4,0,0,0.3};Point(23)={5,0,0,0.3};
Line(21)={21,22};Line(22)={22,23};
Transfinite Curve{1,3}=5;
Periodic Curve {5} = {1} Translate {0, 0, 1};
Degenerated Curve{7};
Compound Curve{21,22};
"""

# Curves first meshed under `Mesh 1`, then re-meshed by a second `Mesh 1`
# inside the same script — exercises the retained/keep arms next to detached
# production (via extra_geo below).
const REMESH=BOXWIRE

function _par_curve_generate(src,mesh_dim,workers;extra_geo="")
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

@testset "parallel detached Curve generation stays deterministic" begin
    cases=(
        ("independent graded curves Mesh1",INDEPENDENT,1,""),
        ("shared-endpoint box wire Mesh1",BOXWIRE,1,""),
        ("closed box Mesh2",BOX,2,""),
        ("closed box Mesh3",BOX,3,""),
        ("mixed TF+periodic+embedded Mesh1",MIXED,1,""),
        ("independent curves re-mesh Mesh1",REMESH,1,"Mesh 1;\n"),
    )
    try
        for (name,src,dim,extra) in cases
            scrc,sbytes,sok=_par_curve_generate(src,dim,1;extra_geo=extra)
            tcrc,tbytes,tok=_par_curve_generate(src,dim,4;extra_geo=extra)
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
