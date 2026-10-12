# `model.occ` shell/solid endpoints: `add_surface_loop` (`BRepPrim_Builder::
# MakeShell` + FORWARD `AddShellFace`, `sewing` gated), `add_volume`
# (`BRepBuilderAPI_MakeSolid`: first shell exterior, rest cavities), and
# `add_box` (`BRepPrimAPI_MakeBox` between opposite corners — signed extents
# allowed, zero extent rejects "Degenerate box"). The box publishes its whole
# sub-entity graph atomically — 8 vertices, 12 edges, 6 wires, 6 plane faces
# sharing the shell's edges, 1 shell, 1 solid — and `synchronize` imports
# shells as `m.surface_loops` (stable across repeated syncs) and solids as
# `m.volumes[tag] = [exterior loop, cavity loops...]`, matching
# `OCC_Internals::addSurfaceLoop`/`addVolume`/`addBox` semantics observed on
# Gmsh 4.15.2.

using Test
using Tessella
const A=Tessella.API
const M=Tessella.Model

function _occ_setup()
    A.initialize()
    return A.model.occ
end

_current_model()=A._with_model(identity)

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "occ.add_box" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0,1.0,2.0,3.0)
    @test v==1
    occ.synchronize()
    m=_current_model()
    # Complete sub-entity graph, all synchronized.
    @test length(m.points)==8
    @test length(m.curves)==12
    @test length(m.surfaces)==6
    @test length(m.surface_loops)==1
    @test m.volumes==Dict(1=>[1])
    # The shell carries the six faces FORWARD.
    @test sort(m.surface_loops[1])==collect(1:6)
    # Every surface is a single outer wire of four edges.
    for tag in keys(m.surfaces)
        loops=m.surfaces[tag]
        @test length(loops)==1
        @test length(m.loops[only(loops)])==4
    end
    # Bounding box is the analytic corner pair.
    xmin,ymin,zmin,xmax,ymax,zmax=A.model.get_bounding_box(3,v)
    @test collect((xmin,ymin,zmin))≈collect((0.0,0.0,0.0))
    @test collect((xmax,ymax,zmax))≈collect((1.0,2.0,3.0))
    # Face normals point outward on the six axis directions.
    normals=NTuple{3,Float64}[]
    for tag in 1:6
        g=m.surface_geometry[tag]
        push!(normals,g.axis)
    end
    @test sort!(normals)==sort!(collect((
        (0.0,0.0,-1.0),(0.0,0.0,1.0),(0.0,-1.0,0.0),
        (0.0,1.0,0.0),(-1.0,0.0,0.0),(1.0,0.0,0.0))))
    # Edges are shared between adjacent faces: each of the 12 edges appears
    # in exactly two distinct face wires.
    edge_faces=Dict{Int,Vector{Int}}()
    for (ftag,loops) in m.surfaces
        for ctag in m.loops[only(loops)]
            e=abs(ctag)
            push!(get!(edge_faces,e,Int[]),ftag)
        end
    end
    @test sort!(collect(keys(edge_faces)))==collect(1:12)
    @test all(length(unique(fs))==2 for fs in values(edge_faces))
    # get_entities reports the bound graph (wires/shells stay internal).
    ents=occ.get_entities()
    @test count(e->first(e)==0,ents)==8
    @test count(e->first(e)==1,ents)==12
    @test count(e->first(e)==2,ents)==6
    @test count(e->first(e)==3,ents)==1
end

@testset "occ.add_box variants" begin
    occ=_occ_setup()
    # Negative extents are opposite corners, matching BRepPrimAPI_MakeBox(P1,P2).
    v=occ.add_box(3.0,3.0,3.0,-1.0,-2.0,-3.0)
    occ.synchronize()
    xmin,ymin,zmin,xmax,ymax,zmax=A.model.get_bounding_box(3,v)
    @test collect((xmin,ymin,zmin))≈collect((2.0,1.0,0.0))
    @test collect((xmax,ymax,zmax))≈collect((3.0,3.0,3.0))
    # Explicit tag honored and collision rejected without mutation.
    tagged=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0;tag=9)
    @test tagged==9
    before=length(occ.get_entities())
    @test occursin("OCC Volume[9] already exists",
        _errtext(()->occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0;tag=9)))
    @test length(occ.get_entities())==before
    # Degenerate extents reject before any sub-entity allocation.
    for ext in ((0.0,1.0,1.0),(1.0,0.0,1.0),(1.0,1.0,0.0))
        @test occursin("Degenerate box",
            _errtext(()->occ.add_box(0.0,0.0,0.0,ext...)))
    end
    @test length(occ.get_entities())==before
    # Non-finite inputs reject on each corner/extent slot.
    @test occursin("must be finite",
        _errtext(()->occ.add_box(NaN,0.0,0.0,1.0,1.0,1.0)))
    @test occursin("must be finite",
        _errtext(()->occ.add_box(0.0,0.0,0.0,1.0,Inf,1.0)))
    @test length(occ.get_entities())==before
end

@testset "occ.add_surface_loop" begin
    occ=_occ_setup()
    # Two disjoint plane faces compose an open shell — the endpoint defers
    # manifold validity like upstream's FORWARD AddShellFace path.
    f1=occ.add_rectangle(0.0,0.0,0.0,1.0,1.0)
    f2=occ.add_rectangle(0.0,0.0,1.0,1.0,1.0)
    l=occ.add_surface_loop([f1,f2])
    @test l==1
    occ.synchronize()
    m=_current_model()
    @test m.surface_loops[l]==[f1,f2]
    # Re-synchronize is idempotent: no duplicate loops are materialized.
    occ.synchronize()
    @test length(m.surface_loops)==1
    # Face tags are used verbatim (no abs): negatives do not bind.
    A.initialize(); occ=_occ_setup()
    fa=occ.add_rectangle(0.0,0.0,0.0,1.0,1.0)
    @test occursin("unknown OCC Surface[-1]",
        _errtext(()->occ.add_surface_loop([-fa])))
    # Empty shells are accepted like upstream's deferred shell builder.
    le=occ.add_surface_loop(Int[])
    occ.synchronize()
    m=_current_model()
    @test m.surface_loops[le]==Int[]
    # sewing=true rejects loudly rather than approximating.
    @test occursin("sewing is not supported",
        _errtext(()->occ.add_surface_loop([fa];sewing=true)))
    @test occursin("sewing must be Bool",
        _errtext(()->occ.add_surface_loop([fa];sewing=1)))
    # Unknown surfaces reject without consuming a loop tag.
    @test occursin("unknown OCC Surface[77]",
        _errtext(()->occ.add_surface_loop([fa,77])))
    @test occ.add_surface_loop([fa])==2
    # Explicit tag honored; collision rejected.
    @test occ.add_surface_loop([fa];tag=5)==5
    @test occursin("OCC Surface Loop[5] already exists",
        _errtext(()->occ.add_surface_loop([fa];tag=5)))
end

@testset "occ.add_volume" begin
    occ=_occ_setup()
    # First shell is the exterior; subsequent shells are cavities.
    outer=occ.add_box(0.0,0.0,0.0,4.0,4.0,4.0)
    inner=occ.add_box(1.0,1.0,1.0,2.0,2.0,2.0)
    occ.synchronize()
    m=_current_model()
    loops=sort!(collect(keys(m.surface_loops)))
    @test length(loops)==2
    v=occ.add_volume(loops)
    occ.synchronize()
    m=_current_model()
    @test m.volumes[v]==loops
    # Swapped order yields a different boundary list (order is preserved).
    v2=occ.add_volume(reverse(loops))
    occ.synchronize()
    @test m.volumes[v2]==reverse(loops)
    # Empty volumes are accepted like upstream's deferred solid builder.
    ve=occ.add_volume(Int[])
    occ.synchronize()
    @test m.volumes[ve]==Int[]
    # Unknown shells reject.
    @test occursin("unknown OCC Surface Loop[42]",
        _errtext(()->occ.add_volume([loops[1],42])))
    # Explicit tag honored; collision rejected.
    @test occ.add_volume([loops[1]];tag=8)==8
    @test occursin("OCC Volume[8] already exists",
        _errtext(()->occ.add_volume([loops[1]];tag=8)))
end

@testset "occ.add_box removal" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    occ.synchronize()
    m=_current_model()
    @test haskey(m.volumes,v)
    # A face inside the bound shell is owned: direct removal is a no-op and
    # the face stays bound (matching the existing owned-shape convention).
    occ.remove([(2,1)])
    @test (2,1) in occ.get_entities()
    occ.remove([(3,v)],true)
    occ.synchronize()
    m=_current_model()
    @test !haskey(m.volumes,v)
    # The recursive removal carries the shell's GModel loop away too.
    @test isempty(m.surface_loops)
end

@testset "occ.add_box volume meshing" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    occ.synchronize()
    A.mesh.generate(3)
    node_tags,coords=A.mesh.get_nodes(3,v,true)
    @test length(node_tags)>0
    for i in 1:3:3length(node_tags)
        x,y,z=coords[i],coords[i+1],coords[i+2]
        @test -1e-9<=x<=1.0+1e-9
        @test -1e-9<=y<=1.0+1e-9
        @test -1e-9<=z<=1.0+1e-9
    end
    element_types,element_tags,_=A.mesh.get_elements(3,v)
    @test element_types==Int32[4]
    @test length(element_tags[1])>0
end
