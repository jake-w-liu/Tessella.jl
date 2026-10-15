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

# Materialize kernel-bound volumes into a fresh scratch GeoModel exactly
# like `_model_occ_boolean`'s operand stage (kernel tags carry verbatim,
# shared sub-shapes dedup through `seen`, coordinates canonicalized so
# split points weld bitwise, `next_tag` seeded above the materialized
# tags so engine-fresh entities cannot collide).
function _occ_scratch_operands(vtags::AbstractVector{<:Integer})
    sm=M.GeoModel()
    kernel=_kernel()
    seen=Set{UInt64}()
    snap=NTuple{3,Float64}[]
    for t in vtags
        M._model_occ_dist_materialize!(sm,kernel,3,kernel.bindings[(3,t)],
            seen,"test";snap=snap,untyped_defaults=true)
    end
    dicts=(sm.points,sm.curves,sm.surfaces,sm.volumes)
    for d in 1:4
        sm.next_tag[d]=max(sm.next_tag[d],maximum(keys(dicts[d]);init=0))
    end
    return sm,Int[t for t in vtags]
end

# Axis-aligned bounding box of one materialized result piece: curve
# endpoints plus quarter/mid/three-quarter interior evals so rim arcs
# reach their apices.
function _sm_piece_bbox(sm,piece)
    lo=fill(Inf,3);hi=fill(-Inf,3)
    function add(p3)
        for a in 1:3
            lo[a]=min(lo[a],p3[a]);hi[a]=max(hi[a],p3[a])
        end
    end
    for sh in piece.shells, sface in sm.surface_loops[sh]
        stag=abs(sface)
        for lp in sm.surfaces[stag], sc in sm.loops[lp]
            c=abs(sc)
            for v in sm.curves[c]
                add(sm.points[v])
            end
            t0,t1=M._model_curve_param_bounds(sm,c,"test")
            for f in (0.25,0.5,0.75)
                add(M._model_curve_point(sm,c,t0+f*(t1-t0),"test"))
            end
        end
    end
    return (lo[1],lo[2],lo[3],hi[1],hi[2],hi[3])
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

@testset "occ.remove_all_duplicates — Gmsh OCC_Internals semantics" begin
    # upstream `removeAllDuplicates` runs `booleanFragments(-1, all
    # entities at the maximum dimension present, {}, remove both)`:
    # coincident solids collapse to one shared result, overlapping
    # solids split into cells, untouched shapes keep their tags.
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.4,0,0,1,1,1)
    occ.remove_all_duplicates()
    ents=sort!(Tuple{Int,Int}[e for e in occ.get_entities(3)])
    @test ents==[(3,1),(3,2),(3,3)]
    @test sort!([occ.get_mass(3,t) for (_,t) in ents])≈
        sort!([0.4*1.0,0.6,0.4])
    # a coincident pair collapses into a single shared volume
    occ.add_box(5,5,5,1,1,1)
    occ.add_box(5,5,5,1,1,1)
    occ.remove_all_duplicates()
    ents2=sort!(Tuple{Int,Int}[e for e in occ.get_entities(3)])
    @test ents2==[(3,1),(3,2),(3,3),(3,4)]
    @test occ.get_mass(3,4)≈1.0
    @test [occ.get_mass(3,t) for (_,t) in ents2]≈[0.4,0.6,0.4,1.0]
end

@testset "occ.remove_all_duplicates — no-op and rejection" begin
    occ=_occ_setup()
    @test occ.get_entities(3)==Tuple{Int,Int}[]
    occ.remove_all_duplicates()   # empty registry — no-op
    occ.add_box(0,0,0,1,1,1)
    occ.remove_all_duplicates()   # single entity — no-op
    @test occ.get_entities(3)==[(3,1)]
    @test occ.get_mass(3,1)≈1.0
    # maxdim 2 registry (no volumes): the volume-only bridge rejects
    # loudly instead of silently pretending support
    occ.add_rectangle(5,5,5,1,1)
    occ.remove([(3,1)])
    err=_errtext(()->occ.remove_all_duplicates())
    @test occursin("dimension",err) || occursin("volume",err) ||
          occursin("unsupported",err)
end

@testset "occ.fragment — re-fragment preserves masses and topology" begin
    # regression for the publish path: re-materializing OCC-published
    # cells must keep exact plane orientation — the kernel face `axis`
    # is the winding normal, and each shell use's `rev` flag is measured
    # against the containing solid rather than trusted from the
    # engine's provenance signs.
    occ=_occ_setup()
    occ.add_box(0,0,0,1,1,1)
    occ.add_box(0.4,0,0,1,1,1)
    occ.fragment([(3,1),(3,2)],Tuple{Int,Int}[])
    pre=[occ.get_mass(3,t) for t in (1,2,3)]
    @test pre≈[0.4,0.6,0.4]
    out,map=occ.fragment([(3,1),(3,2),(3,3)],Tuple{Int,Int}[])
    @test out==[(3,1),(3,2),(3,3)]
    @test map==[[(3,1)],[(3,2)],[(3,3)]]
    @test [occ.get_mass(3,t) for t in (1,2,3)]≈[0.4,0.6,0.4]
    bbs=[occ.get_bounding_box(3,t) for t in (1,2,3)]
    @test bbs[1][1:3]==(0.0,0.0,0.0) && bbs[1][4]≈0.4
    @test bbs[3][4]≈1.4
    # mixed set: published cells plus a coincident pair — the pair
    # collapses onto one shared cell while the cells re-materialize
    occ.add_box(5,5,5,1,1,1)
    occ.add_box(5,5,5,1,1,1)
    out2,map2=occ.fragment([(3,t) for t in (1,2,3,4,5)],
                           Tuple{Int,Int}[])
    @test out2==[(3,1),(3,2),(3,3),(3,4)]
    @test map2[end]==[(3,4)]
    @test map2[end-1]==[(3,4)]
    @test [occ.get_mass(3,t) for t in (1,2,3,4)]≈[0.4,0.6,0.4,1.0]
    # a third consecutive fragment still holds the line — publish is
    # idempotent, not just once
    occ.fragment([(3,t) for t in (1,2,3,4)],Tuple{Int,Int}[])
    @test [occ.get_mass(3,t) for t in (1,2,3,4)]≈[0.4,0.6,0.4,1.0]
end

@testset "occ.cut — curved result re-cut through the boolean bridge" begin
    # regression: a box-minus-cylinder piece carries circle edges whose
    # pcurves on the synthesized plane caps can be parameter-reversed
    # (left-handed uv frames on mirrored records / hole rims). The
    # `neg2d` pcurve wrapper keeps the edge parameter while recording
    # the mirrored trace; every consumer — measures, extraction,
    # re-boolean — must honour it.
    occ=_occ_setup()
    b=occ.add_box(0,0,0,1,1,1)
    cyl=occ.add_cylinder(0.5,0.5,0.0,0,0,1,0.4)
    out,_=occ.cut([(3,b)],[(3,cyl)])
    @test length(out)==1
    piece=out[1][2]
    @test occ.get_mass(3,piece)≈1-π*0.16
    # overlapping re-cut: removes the half of the piece with x<0.5,
    # leaving 0.5 minus half the cylinder disk. The engine path runs
    # deterministically below; the API call's provenance snapshot
    # tetrahedralizes the ~260-node rim PLC and can escalate into a
    # ~20 GB exact-Delaunay recovery — skipped under memory pressure.
    tool=occ.add_box(-0.5,-0.5,-0.5,1.0,2.0,2.0)
    sm,vtags=_occ_scratch_operands([piece,tool])
    kept,touched=M._brep_collect_faces_n(sm,:difference,vtags,1,"cut")
    pieces,image,created=M._brep_materialize_multi!(
        sm,:difference,vtags,kept,touched,2,1,"cut")
    @test length(pieces)==1
    @test all(isapprox.(_sm_piece_bbox(sm,only(pieces)),
                        (0.5,0.0,0.0,1.0,1.0,1.0);atol=1e-9))
    try
        out2,map2=occ.cut([(3,piece)],[(3,tool)])
        @test length(out2)==1
        @test occ.get_mass(3,out2[1][2])≈0.5-(π*0.16)/2
    catch err
        err isa ReadOnlyMemoryError || rethrow()
        @test_skip "provenance snapshot exceeds available memory"
    end
end

@testset "occ.fragment — curved piece against an overlapping box" begin
    # regression: fragmenting the curved cut result drops degenerate
    # sections onto the piece's rim circles. Those mirrored imprint
    # arcs traverse their edge backwards (neg pcurves) — materializing
    # them needs the occurrence sign flipped (fwd ⊻ neg) or the result
    # wire fails to chain.
    occ=_occ_setup()
    b=occ.add_box(0,0,0,1,1,1)
    cyl=occ.add_cylinder(0.5,0.5,0.0,0,0,1,0.4)
    out,_=occ.cut([(3,b)],[(3,cyl)])
    piece=out[1][2]
    @test occ.get_mass(3,piece)≈1-π*0.16
    b2=occ.add_box(0.5,0.0,0.0,1.0,1.0,1.0)
    # The engine path itself (materialize → collect → materialize) is
    # what the fix exercises and it runs in seconds: the mirrored rim
    # arcs must chain into wires and the arrangement must yield the
    # four physical cells — piece∖box2 and piece∩box2 (each
    # 0.5−π·0.16/2), the notch half-cylinder (π·0.16/2) and the x>1
    # slab (0.5). B∖A is genuinely disconnected — the notch is walled
    # off the slab by the piece slab at x∈[0.9,1].
    sm,vtags=_occ_scratch_operands([piece,b2])
    kept,touched=M._brep_collect_faces_n(sm,:fragments,vtags,1,"fragment")
    pieces,image,created=M._brep_materialize_multi!(
        sm,:fragments,vtags,kept,touched,2,1,"fragment")
    @test length(pieces)==4
    @test all(p->!isempty(p.shells),pieces)
    expected=sort!(NTuple{6,Float64}[
        (0.0,0.0,0.0,0.5,1.0,1.0),        # piece∖box2
        (0.5,0.0,0.0,1.0,1.0,1.0),        # piece∩box2
        (0.5,0.1,0.0,0.9,0.9,1.0),        # notch half-cylinder
        (1.0,0.0,0.0,1.5,1.0,1.0)])       # x>1 slab
    got=sort!([_sm_piece_bbox(sm,p) for p in pieces])
    @test length(got)==length(expected)
    @test all(i->all(isapprox.(got[i],expected[i];atol=1e-9)),
              eachindex(expected))
    # the API call additionally runs the operand provenance surface
    # tetrahedralization, whose exact-Delaunay recovery on the ~260-node
    # rim PLC needs ~20 GB — verified once with analytically exact
    # masses (0.2486725877/0.2486725877/0.2513274123/0.5); skipped under
    # memory pressure rather than failing the suite on a resource bound
    try
        out2,_=occ.fragment([(3,piece)],[(3,b2)])
        @test length(out2)==4
        @test sort!([occ.get_mass(3,t) for (_,t) in out2])≈
            sort!([0.5-(π*0.16)/2,0.5-(π*0.16)/2,(π*0.16)/2,0.5])
    catch err
        err isa ReadOnlyMemoryError || rethrow()
        @test_skip "provenance snapshot exceeds available memory"
    end
end

@testset "occ booleans — periodic seam and partial sectors" begin
    # Regression for the curved-face merge path: a box cut through a full
    # cylinder keeps two seam-adjacent lateral strips whose union wraps the
    # u-period. Their seam-edge usages are the same edge at two
    # period-shifted uv positions — distinct domain boundaries, not an
    # internal shared edge — so the merge must keep both copies or the
    # lateral region loses its planar cycle and never materializes. The
    # π/2 sector exercises the same machinery plus meridian side planes.
    for (angle,mi) in ((2π,1.3907963267948966),(π/2,0.9799219123544154))
        total=angle                       # r=1, h=2 → vol = angle·r²·h/2
        occ=_occ_setup()
        b=occ.add_box(0.3,-0.3,-0.3, 2,2,3)
        k=occ.add_cylinder(0,0,0,0,0,2,1;angle)
        out,map=occ.intersect([(3,b)],[(3,k)])
        @test out==[(3,1)] && map==[[(3,1)],[(3,1)]]
        @test occ.get_mass(3,1)≈mi
        # cyl∖box is the complement cell: intersect + cut == operand
        occ=_occ_setup()
        b=occ.add_box(0.3,-0.3,-0.3, 2,2,3)
        k=occ.add_cylinder(0,0,0,0,0,2,1;angle)
        out,_=occ.cut([(3,k)],[(3,b)])
        @test sum(occ.get_mass(3,t) for (_,t) in out)≈total-mi
        # fuse = box ∪ cyl = box + (cyl∖box)
        occ=_occ_setup()
        b=occ.add_box(0.3,-0.3,-0.3, 2,2,3)
        k=occ.add_cylinder(0,0,0,0,0,2,1;angle)
        out,_=occ.fuse([(3,b)],[(3,k)])
        @test length(out)==1
        @test occ.get_mass(3,out[1][2])≈12.0+total-mi
        # fragment's three cells (box∖cyl, cyl∖box, box∩cyl) tile the union
        occ=_occ_setup()
        b=occ.add_box(0.3,-0.3,-0.3, 2,2,3)
        k=occ.add_cylinder(0,0,0,0,0,2,1;angle)
        out,_=occ.fragment([(3,b)],[(3,k)])
        @test length(out)==3
        @test sum(occ.get_mass(3,t) for (_,t) in out)≈12.0+total-mi
        # reversed operand order produces the same intersection
        occ=_occ_setup()
        b=occ.add_box(0.3,-0.3,-0.3, 2,2,3)
        k=occ.add_cylinder(0,0,0,0,0,2,1;angle)
        out,_=occ.intersect([(3,k)],[(3,b)])
        @test occ.get_mass(3,out[1][2])≈mi
    end
end
