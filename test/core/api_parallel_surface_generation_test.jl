# PAR01 slice — parallel mesh *generation* of independent Surfaces under
# `General.NumThreads`. `_higher_surface_phase!` detaches runs of eligible
# singleton Surfaces onto per-unit transaction copies under `Threads.@threads`
# and merges the receipts serially in unit order; stale or unsupported slots
# re-run the identical live serial call, so the produced mesh — its CRC, its
# serialized MSH bytes and its validation — must be identical for any
# requested worker count. Run standalone (e.g.
# `julia -t4 test/core/api_parallel_surface_generation_test.jl`); the default
# suite does not register this file. Under `-t1` both arms take the serial
# path and the checks still hold.

using Test
using Tessella
using Tessella.Elements: mixed_crc, write_mixed_msh

A=Tessella.API

const BOX="""
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
"""

# Two disjoint surfaces sharing no Curve: the detached batch really runs in
# parallel instead of falling back on dirty footprints.
const INDEPENDENT="""
Point(1)={0,0,0,0.4};Point(2)={1,0,0,0.4};Point(3)={1,1,0,0.4};Point(4)={0,1,0,0.4};
Point(5)={2,0,0,0.25};Point(6)={3.5,0,0,0.25};Point(7)={3.5,1,0,0.25};Point(8)={2,1,0,0.25};
Point(9)={5,0,0,0.3};Point(10)={6,0,0,0.3};Point(11)={6,1,0,0.3};Point(12)={5,1,0,0.3};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={9,10};Line(10)={10,11};Line(11)={11,12};Line(12)={12,9};
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};Curve Loop(3)={9,10,11,12};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
"""

const TFBOX=BOX*"""
Transfinite Curve{1,2,3,4,5,6,7,8,9,10,11,12}=3;
Transfinite Surface{1,2,3,4,5,6};
Transfinite Volume{1}={1,2,3,4,5,6,7,8};
"""

const MIXED=BOX*"""
Point(21)={3,0,0,0.35};Point(22)={4,0,0,0.35};Point(23)={4,1,0,0.35};Point(24)={3,1,0,0.35};
Line(21)={21,22};Line(22)={22,23};Line(23)={23,24};Line(24)={24,21};
Curve Loop(21)={21,22,23,24};Plane Surface(21)={21};
Transfinite Curve{21,22,23,24}=4;Transfinite Surface{21}={21,22,23,24};
Recombine Surface{21};
"""

function _par_surface_generate(src,mesh_dim,workers;extra_geo="")
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

@testset "parallel detached Surface generation stays deterministic" begin
    cases=(
        ("independent surfaces Mesh2",INDEPENDENT,2,""),
        ("single box shared curves Mesh3",BOX,3,""),
        ("single box Mesh2",BOX,2,""),
        ("transfinite box Mesh3",TFBOX,3,""),
        ("mixed native+TF+recombine Mesh2",MIXED,2,""),
        ("re-mesh existing authority Mesh3",BOX,3,"Mesh 3;\n"),
        ("recombined box Mesh3",BOX,3,"Recombine Surface{1,2,3,4,5,6};\n"),
    )
    try
        for (name,src,dim,extra) in cases
            scrc,sbytes,sok=_par_surface_generate(src,dim,1;extra_geo=extra)
            tcrc,tbytes,tok=_par_surface_generate(src,dim,4;extra_geo=extra)
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
