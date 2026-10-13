# `model.occ.fuse` / `cut` / `intersect` / `fragment` —
# `OCC_Internals::booleanOperator` on the pending-shape kernel. Bounded
# scope: volume (dim-3) operands materialize into a scratch GeoModel
# through the closest-entities materializer (kernel tags carry verbatim,
# so shared sub-shapes dedup through `seen`), the native analytic BREP
# engine runs there, and the result volumes publish back to the registry.
# Disposition mirrors `occBooleanPreserveNumbering`: preserved operands
# keep their binding, 1:1 images rebind the operand tag, deleted operands
# unbind recursively, and fresh pieces allocate tags after the surviving
# maximum in `out` order. `outDimTagsMap` reports each operand's image
# (itself when unmodified, the pieces it modified into otherwise, empty
# when deleted). Engine limitations (non-representable section pcurves,
# the nested containment + overlap fragment family) reject loudly before
# any kernel mutation — staging guarantees atomicity.

using Test
using Tessella
const A=Tessella.API
const M=Tessella.Model

function _occ_setup()
    A.initialize()
    return A.model.occ
end

_current_model()=A._with_model(identity)
_kernel()=A._with_model(m->m.occ_kernel)
_errtext(f)=(try;f();catch err;return sprint(showerror,err);end;"")

function _occ_counts()
    d=Dict{Int,Int}()
    for (dim,tag) in A.model.occ.get_entities()
        d[dim]=get(d,dim,0)+1
    end
    return d
end

@testset "occ.fuse — overlap union" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    out,map=occ.fuse([(3,1)],[(3,2)])
    @test out==[(3,1)]
    @test map==[[(3,1)],[(3,1)]]
    @test occ.get_entities(3)==[(3,1)]
    @test occ.get_mass(3,1)≈1.875
    @test collect(occ.get_center_of_mass(3,1))≈collect((0.75,0.75,0.75))
    @test collect(occ.get_bounding_box(3,1))≈
        collect((0.0,0.0,0.0,1.5,1.5,1.5))
    # a result entity behaves like any kernel shape
    cp=occ.copy([(3,1)])
    @test occ.get_mass(3,cp[1][2])≈1.875
    occ.translate(cp,10,0,0)
    bb=occ.get_bounding_box(3,cp[1][2])
    @test bb[1]≈10.0 && bb[4]≈11.5
    occ.remove(cp,true)
    @test occ.get_entities(3)==[(3,1)]
end

@testset "occ.fuse — disjoint operands preserve numbering" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(5,0,0,1,1,1)
    occ.synchronize()
    out,map=occ.fuse([(3,1)],[(3,2)])
    @test out==[(3,1),(3,2)]
    @test map==[[(3,1)],[(3,2)]]
    @test sort(occ.get_entities(3))==[(3,1),(3,2)]
    @test occ.get_mass(3,1)≈1.0 && occ.get_mass(3,2)≈1.0
end

@testset "occ.cut — single piece, tool removed" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    out,map=occ.cut([(3,1)],[(3,2)])
    @test out==[(3,1)]
    @test map==[[(3,1)],Tuple{Int,Int}[]]
    @test occ.get_entities(3)==[(3,1)]
    @test occ.get_mass(3,1)≈0.875
    # the result's full boundary publishes as bound sub-entities: three
    # unmodified box faces, three notched L-faces, three new notch faces
    @test length(occ.get_entities(2))==9
    @test length(occ.get_entities(1))==21
    @test length(occ.get_entities(0))==14
end

@testset "occ.cut — split produces two pieces" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.4,-0.5,-0.5,0.2,2,2)
    occ.synchronize()
    out,map=occ.cut([(3,1)],[(3,2)])
    @test out==[(3,1),(3,2)]
    @test map==[[(3,1),(3,2)],Tuple{Int,Int}[]]
    @test sort!(collect(occ.get_mass(3,t) for (_,t) in out))≈[0.4,0.4]
end

@testset "occ.cut — fresh tags follow the surviving maximum" begin
    occ=_occ_setup()
    # upstream allocates result tags as `maxTag+1` *after* operand removal:
    # deleting the high-tagged tool frees its claim, so the split pieces
    # renumber from the surviving maximum (empty here) — not from the
    # pre-operation maximum (upstream parity regression)
    occ.add_box(0,0,0,1,1,1;tag=5)
    occ.add_box(0.4,-0.5,-0.5,0.2,2,2;tag=9)
    occ.synchronize()
    out,map=occ.cut([(3,5)],[(3,9)])
    @test out==[(3,1),(3,2)]
    @test map==[[(3,1),(3,2)],Tuple{Int,Int}[]]
    @test sort(occ.get_entities(3))==[(3,1),(3,2)]
end

@testset "occ.intersect — corner overlap" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    out,map=occ.intersect([(3,1)],[(3,2)])
    @test out==[(3,1)]
    @test map==[[(3,1)],[(3,1)]]
    @test occ.get_entities(3)==[(3,1)]
    @test occ.get_mass(3,1)≈0.125
    @test collect(occ.get_center_of_mass(3,1))≈collect((0.75,0.75,0.75))
end

@testset "occ.fragment — two-box cells and images" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    out,map=occ.fragment([(3,1)],[(3,2)])
    @test out==[(3,1),(3,2),(3,3)]
    @test map==[[(3,1),(3,2)],[(3,2),(3,3)]]
    @test sort(occ.get_entities(3))==[(3,1),(3,2),(3,3)]
    masses=sort!(collect(occ.get_mass(3,t) for (_,t) in out))
    @test masses≈[0.125,0.875,0.875]
    @test sum(masses)≈1.875
end

@testset "occ.fragment — disjoint operand keeps its tag" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0,0,1,1,1)
    occ.add_box(5,0,0,1,1,1)
    occ.synchronize()
    out,map=occ.fragment([(3,1),(3,2)],[(3,3)])
    # the untouched tool rebinds its own tag and reports last
    @test out[end]==(3,3)
    @test map[3]==[(3,3)]
    @test length(out)==4 && length(out[1:end-1])==3
    @test sort(occ.get_entities(3))==[(3,3),(3,4),(3,5),(3,6)]
    @test occ.get_mass(3,3)≈1.0
    @test sum(occ.get_mass(3,t) for (_,t) in out)≈2.5
end

@testset "occ.fuse — objects-only operand list" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0,0,1,1,1)
    occ.synchronize()
    out,map=occ.fuse([(3,1),(3,2)],[])
    @test out==[(3,1)]
    @test map==[[(3,1)],[(3,1)]]
    @test occ.get_mass(3,1)≈1.5
end

@testset "occ booleans — remove flags" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    out,map=occ.cut([(3,1)],[(3,2)];removeObject=false,removeTool=false)
    # upstream keeps the originals bound; the difference piece binds fresh
    # and reports as the object's image
    @test out==[(3,3)]
    @test map==[[(3,3)],Tuple{Int,Int}[]]
    @test sort(occ.get_entities(3))==[(3,1),(3,2),(3,3)]
    @test occ.get_mass(3,1)≈1.0
    @test occ.get_mass(3,2)≈1.0
    @test occ.get_mass(3,3)≈0.875
end

@testset "occ booleans — explicit tag" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    out,map=occ.cut([(3,1)],[(3,2)];tag=42)
    @test out==[(3,42)]
    @test map==[[(3,42)],Tuple{Int,Int}[]]
    @test occ.get_entities(3)==[(3,42)]
    @test occ.get_mass(3,42)≈0.875
    # an explicit tag on a bound entity rejects
    occ.add_box(0,0,0,1,1,1)
    occ.synchronize()
    err=_errtext(()->occ.fuse([(3,42)],[(3,43)];tag=42))
    @test occursin("already exists",err)
end

@testset "occ booleans — explicit tag multi-piece rejection" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.4,-0.5,-0.5,0.2,2,2)
    occ.synchronize()
    before=_occ_counts()
    err=_errtext(()->occ.cut([(3,1)],[(3,2)];tag=77))
    @test occursin("cannot bind 2 result volumes",err)
    @test _occ_counts()==before
    @test sort(occ.get_entities(3))==[(3,1),(3,2)]
end

@testset "occ.cut — empty result" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(5,0,0,1,1,1)
    occ.synchronize()
    out,map=occ.intersect([(3,1)],[(3,2)])
    @test out==Tuple{Int,Int}[]
    @test map==[Tuple{Int,Int}[],Tuple{Int,Int}[]]
    @test occ.get_entities(3)==Tuple{Int,Int}[]
end

@testset "occ booleans — validation and atomicity" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    err=_errtext(()->occ.fuse([(3,99)],[(3,1)]))
    @test occursin("Unknown OpenCASCADE entity of dimension 3 with tag 99",
                   err)
    p=occ.add_point(0,0,0)
    err=_errtext(()->occ.fuse([(0,p)],[(3,1)]))
    @test occursin("volumes only",err)
    err=_errtext(()->occ.fuse([(3,1)],[(3,2)];tag=1.5))
    @test occursin("tag must be an integer",err)
    err=_errtext(()->occ.fuse([(3,1)],[(3,2)];removeObject=1))
    @test occursin("removeObject must be Bool",err)
    # empty object list is an upstream no-op
    @test occ.fuse([],[(3,1)])==(Tuple{Int,Int}[],Vector{Tuple{Int,Int}}[])
end

@testset "occ booleans — engine rejection stays atomic" begin
    occ=_occ_setup()
    # nested containment inside an overlapping third operand — the engine
    # does not yet materialize this fragment family; the rejection must
    # leave the registry byte-identical
    occ.add_box(0.5,0.5,0.5,0.5,0.5,0.5)
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    snap=occ.get_entities()
    counts=_occ_counts()
    err=_errtext(()->occ.fragment([(3,1),(3,2)],[(3,3)]))
    @test occursin("Curve",err) || occursin("shell",err)
    @test occ.get_entities()==snap
    @test _occ_counts()==counts
    # the model still answers afterwards
    out,_=occ.cut([(3,2)],[(3,1)])
    @test occ.get_mass(3,out[1][2])≈0.875
end

@testset "occ.fuse — post-sync queries and mesh" begin
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.5,0.5,0.5,1,1,1)
    occ.synchronize()
    out,_=occ.fuse([(3,1)],[(3,2)])
    occ.synchronize()
    @test occ.get_entities(3)==[(3,1)]
    @test occ.get_mass(3,1)≈1.875
    ents,dist,pts=occ.get_closest_entities(2.0,0.5,0.5,[(3,1)],1)
    @test ents==[(3,1)] && dist[1]≈0.5
    d,_=occ.get_distance(3,1,0,occ.add_point(3,3,3))
    @test d≈sqrt(6.75) rtol=1e-10
    A.mesh.generate(3)
    nt,_=A.mesh.get_nodes(3,1,true)
    @test length(nt)>0
    # every node lies inside the fused bounding box
    _,coords=A.mesh.get_nodes(3,1,true)
    @test all(-1e-9<=coords[3(i-1)+1]<=1.5+1e-9 for i in 1:length(nt))
end
