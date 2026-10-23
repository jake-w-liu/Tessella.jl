# PAR01 tail — public-contract determinism of parallel mesh *generation*
# under `General.NumThreads`.  The Curve/Surface/Region phases detach runs of
# independent entities onto per-unit transaction copies (plus, since the
# eager scheduler, later eligible slots whose footprints are disjoint from
# every unit in between) and merge the receipts serially in entity order;
# stale or unsupported slots re-run the identical live serial call.  The
# produced mesh — per-entity node/element counts, the actual node and
# element UID assignment, its CRC, its serialized MSH4 bytes and its
# validation — must be identical for any requested worker count.  Run
# standalone (e.g. `julia -t4 test/core/api_parallel_determinism_test.jl`);
# the default suite does not register this file.  Under `-t1` both arms take
# the serial path and the checks still hold.
#
# Fixture notes: the `ESE_*` sources put a prepare-ineligible entity between
# two eligible, footprint-disjoint ones — the pattern the non-contiguous
# eager scheduler exists for.  `BOX` exercises the stale-slot serial fallback
# (all six faces share Curves, so only the first detached receipt can stay
# fresh).  `TWOBOX` exercises the Region arm's `Mesh.RandomSeed` ledger on
# genuinely detached OCC-factory fills.

using Test
using Tessella
using Tessella.Elements: mixed_crc, write_mixed_msh

A=Tessella.API

# Three disjoint Curves; the middle one is transfinite, so it keeps the
# serial arm while the outer two may be produced in one detached batch.
const ESE_CURVES="""
Point(1)={0,0,0,0.25};Point(2)={1,0,0,0.25};
Point(3)={2,0,0,0.15};Point(4)={3,0,0,0.15};
Point(5)={4,0,0,0.2};Point(6)={5,0,0,0.2};
Line(1)={1,2};Line(2)={3,4};Line(3)={5,6};
Transfinite Curve{2}=9;
"""

# A box volume with an embedded interior Surface (prepare-ineligible),
# followed by two standalone disjoint Surfaces — the embedded sheet separates
# the eligible runs.
const ESE_SURFACES="""
Point(1)={0,0,0,0.4};Point(2)={1,0,0,0.4};Point(3)={1,1,0,0.4};Point(4)={0,1,0,0.4};
Point(5)={0,0,1,0.4};Point(6)={1,0,1,0.4};Point(7)={1,1,1,0.4};Point(8)={0,1,1,0.4};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};
Curve Loop(3)={1,10,-5,-9};Curve Loop(4)={2,11,-6,-10};
Curve Loop(5)={3,12,-7,-11};Curve Loop(6)={4,9,-8,-12};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
Plane Surface(4)={4};Plane Surface(5)={5};Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};
Point(21)={0.25,0.25,0.5,0.15};Point(22)={0.75,0.25,0.5,0.15};
Point(23)={0.75,0.75,0.5,0.15};Point(24)={0.25,0.75,0.5,0.15};
Line(21)={21,22};Line(22)={22,23};Line(23)={23,24};Line(24)={24,21};
Curve Loop(21)={21,22,23,24};Plane Surface(7)={21};
Surface{7} In Volume{1};
Point(31)={3,0,0,0.3};Point(32)={4,0,0,0.3};Point(33)={4,1,0,0.3};Point(34)={3,1,0,0.3};
Line(31)={31,32};Line(32)={32,33};Line(33)={33,34};Line(34)={34,31};
Curve Loop(31)={31,32,33,34};Plane Surface(8)={31};
Point(41)={5,0,0,0.35};Point(42)={6,0,0,0.35};Point(43)={6,1,0,0.35};Point(44)={5,1,0,0.35};
Line(41)={41,42};Line(42)={42,43};Line(43)={43,44};Line(44)={44,41};
Curve Loop(41)={41,42,43,44};Plane Surface(9)={41};
"""

# An OCC factory Box, an explicit transfinite box (serial arm), then another
# OCC Box — two eligible Regions separated by an ineligible one.  The
# built-in tags start at 101 so they cannot collide with the OCC kernels'
# auto-assigned subentity tags.
const ESE_VOLUMES="""
SetFactory("OpenCASCADE");
Box(1)={0,0,0, 1,1,1};
SetFactory("Built-in");
Point(101)={3,0,0,0.4};Point(102)={4,0,0,0.4};Point(103)={4,1,0,0.4};Point(104)={3,1,0,0.4};
Point(105)={3,0,1,0.4};Point(106)={4,0,1,0.4};Point(107)={4,1,1,0.4};Point(108)={3,1,1,0.4};
Line(101)={101,102};Line(102)={102,103};Line(103)={103,104};Line(104)={104,101};
Line(105)={105,106};Line(106)={106,107};Line(107)={107,108};Line(108)={108,105};
Line(109)={101,105};Line(110)={102,106};Line(111)={103,107};Line(112)={104,108};
Curve Loop(101)={101,102,103,104};Curve Loop(102)={105,106,107,108};
Curve Loop(103)={101,110,-105,-109};Curve Loop(104)={102,111,-106,-110};
Curve Loop(105)={103,112,-107,-111};Curve Loop(106)={104,109,-108,-112};
Plane Surface(101)={101};Plane Surface(102)={102};Plane Surface(103)={103};
Plane Surface(104)={104};Plane Surface(105)={105};Plane Surface(106)={106};
Surface Loop(101)={101,102,103,104,105,106};Volume(102)={101};
Transfinite Curve{101,102,103,104,105,106,107,108,109,110,111,112}=3;
Transfinite Surface{101,102,103,104,105,106};
Transfinite Volume{102}={101,102,103,104,105,106,107,108};
SetFactory("OpenCASCADE");
Box(3)={6,0,0, 1,1,1};
"""

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

const TWOBOX="""
SetFactory("OpenCASCADE");
Box(1)={0,0,0, 1,1,1};
Box(2)={2.5,0,0, 1,1.5,1};
"""

# Per-entity node/element counts and actual UID assignment from the actual
# mesh's entity data — the public contract the parallel arm must preserve.
function _par_uid_table(mesh)
    data=mesh.entity_data
    nodes=Dict{Tuple{Int,Int32},Vector{UInt64}}()
    elements=Dict{Tuple{Int,Int32},Vector{UInt64}}()
    data===nothing && return (nodes,elements)
    for (row,owner) in enumerate(data.node_entities)
        push!(get!(()->UInt64[],nodes,owner),data.external_node_tags[row])
    end
    for (block,owners) in enumerate(data.block_entities),
        (slot,owner) in enumerate(owners)
        push!(get!(()->UInt64[],elements,
            (Tessella.Elements._block_dim(mesh.blocks[block]),Int32(owner))),
            data.external_element_tags[block][slot])
    end
    return (nodes,elements)
end

function _par_generate(src,mesh_dim,workers;extra_geo="")
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
    node_uids,element_uids=_par_uid_table(actual)
    result=(mixed_crc(actual),read(msh),validate(actual).ok,
            node_uids,element_uids)
    A.finalize()
    return result
end

@testset "parallel generation keeps counts and UID assignment deterministic" begin
    cases=(
        ("disjoint curves around transfinite Mesh1",ESE_CURVES,1,""),
        ("embedded sheet between surfaces Mesh2",ESE_SURFACES,2,""),
        ("embedded sheet between surfaces Mesh3",ESE_SURFACES,3,""),
        ("OCC box + transfinite + OCC box Mesh3",ESE_VOLUMES,3,""),
        ("shared-curve box Mesh3",BOX,3,""),
        ("shared-curve box Mesh2",BOX,2,""),
        ("two OCC boxes Mesh3",TWOBOX,3,""),
        ("two OCC boxes re-mesh Mesh3",TWOBOX,3,"Mesh 3;\n"),
    )
    try
        for (name,src,dim,extra) in cases
            scrc,sbytes,sok,snodes,selements=_par_generate(src,dim,1;extra_geo=extra)
            tcrc,tbytes,tnodes_ok,tnodes,telements=_par_generate(src,dim,4;extra_geo=extra)
            @test sok
            @test tnodes_ok
            # Node/element counts per owning entity.
            @test Dict(k=>length(v) for (k,v) in snodes)==
                Dict(k=>length(v) for (k,v) in tnodes)
            @test Dict(k=>length(v) for (k,v) in selements)==
                Dict(k=>length(v) for (k,v) in telements)
            # Actual UID assignment per owning entity.
            @test snodes==tnodes
            @test selements==telements
            @test scrc==tcrc
            @test sbytes==tbytes
            (scrc==tcrc && sbytes==tbytes && snodes==tnodes &&
                selements==telements && sok && tnodes_ok) ||
                println("PARALLEL-DIFFERS: $name")
        end
    finally
        A.finalize()
    end
end
