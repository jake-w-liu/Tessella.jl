# `model.occ` wire and planar-face endpoints: `add_wire` (signed head-to-tail
# chaining, first-tag reversal, `checkClosed`), `add_curve_loop`, and the
# `BRepBuilderAPI_MakeFace` pair `add_plane_surface` / `add_rectangle`.
# Wires bind as dimension -1 shape records; faces carry a fitted `:plane`
# record with per-edge `lin2d`/`circ2d` pcurves that `synchronize` imports
# as `m.surfaces[tag] = [outer loop, hole loops...]` — matching
# `OCC_Internals::addWire`/`addCurveLoop`/`addPlaneSurface`/`addRectangle`
# semantics observed on Gmsh 4.15.2.

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

# A unit-square point/edge set on z=0.
function _occ_square(occ,origin=(0.0,0.0))
    ox,oy=origin
    p=[occ.add_point(ox+x,oy+y,0.0) for (x,y) in [(0,0),(1,0),(1,1),(0,1)]]
    e=[occ.add_line(p[a],p[b]) for (a,b) in [(1,2),(2,3),(3,4),(4,1)]]
    return p,e
end

# A square of side `side` with lower-left corner at (ox,oy).
function _occ_box(occ,ox,oy,side)
    p=[occ.add_point(ox+x,oy+y,0.0) for (x,y) in
       [(0,0),(side,0),(side,side),(0,side)]]
    e=[occ.add_line(p[a],p[b]) for (a,b) in [(1,2),(2,3),(3,4),(4,1)]]
    return p,e
end

@testset "OCC add_wire chaining and tags" begin
    occ=_occ_setup()
    p,e=_occ_square(occ)
    w=occ.add_wire(e)
    @test w==1
    @test occ.get_max_tag(-1)==1
    # Wires are kernel shapes, not public dim entities.
    @test sort!(collect(occ.get_entities(1)))==[(1,1),(1,2),(1,3),(1,4)]
    @test sort!(collect(occ.get_entities()))==[(0,1),(0,2),(0,3),(0,4),
                                              (1,1),(1,2),(1,3),(1,4)]
    # Explicit tag honored; a second wire gets the next -1 tag.
    @test occ.add_wire(e;tag=9)==9
    @test occ.add_wire(e)==10
    @test occursin("already exists",_errtext(()->occ.add_wire(e;tag=9)))
    # Negative first tag: the synced wire shape must carry the reversed
    # use list (order and orientation flipped together).
    p2,e2=_occ_square(occ,(2.0,0.0))
    wr=occ.add_wire([-e2[1],e2[2],e2[3],e2[4]])
    fr=occ.add_plane_surface([wr])
    occ.synchronize()
    m=_current_model()
    @test m.loops[only(m.surfaces[fr])]==[-8,-7,-6,-5]
end

@testset "OCC add_wire mixed sign chaining" begin
    occ=_occ_setup()
    p1=occ.add_point(0.0,0.0,0.0)
    p2=occ.add_point(1.0,0.0,0.0)
    p3=occ.add_point(1.0,1.0,0.0)
    e1=occ.add_line(p1,p2)
    e2r=occ.add_line(p3,p2)   # defined v3->v2: must chain REVERSED after e1
    e3r=occ.add_line(p3,p1)   # defined v3->v1: must chain FORWARD after e2r
    w=occ.add_wire([e1,e2r,e3r])
    f=occ.add_plane_surface([w])
    occ.synchronize()
    m=_current_model()
    @test m.loops[only(m.surfaces[f])]==[1,-2,3]
end

@testset "OCC add_wire guards" begin
    occ=_occ_setup()
    p,e=_occ_square(occ)
    @test occursin("Could not create wire",_errtext(()->occ.add_wire(Int[])))
    @test occursin("unknown OCC Curve[99]",_errtext(()->occ.add_wire([99])))
    @test occursin("checkClosed must be a boolean",
                   _errtext(()->occ.add_wire(e;checkClosed=1)))
    # Disconnected chain (e1 then e4 both start at a junction but e3 cannot
    # continue from the tail): MakeWire semantics reject it.
    @test occursin("Could not create wire",_errtext(()->occ.add_wire([e[1],e[3]])))
    @test occursin("Curve loop is not closed",
                   _errtext(()->occ.add_wire([e[1],e[2]];checkClosed=true)))
    @test occursin("Curve loop is not closed",
                   _errtext(()->occ.add_curve_loop([e[1],e[2]])))
    @test occursin("tag must be an integer",_errtext(()->occ.add_wire(e;tag="w")))
    @test occursin("tag must be an integer",_errtext(()->occ.add_wire([e[1],1.5])))
    # Failures bind nothing.
    @test occ.get_max_tag(-1)==0
end

@testset "OCC add_curve_loop" begin
    occ=_occ_setup()
    p,e=_occ_square(occ)
    @test occ.add_curve_loop(e)==1
    @test occ.get_max_tag(-1)==1
    # Closed check passes on a proper ring and the wire binds.
    f=occ.add_plane_surface([1])
    @test f==1
end

@testset "OCC add_plane_surface square" begin
    occ=_occ_setup()
    p,e=_occ_square(occ)
    w=occ.add_wire(e)
    f=occ.add_plane_surface([w])
    @test f==1
    @test sort!(collect(occ.get_entities(2)))==[(2,1)]
    occ.synchronize()
    m=_current_model()
    @test m.loops[only(m.surfaces[f])]==e
    @test m.surface_types[f]==:plane
    g=m.surface_geometry[f]
    @test g.axis==(0.0,0.0,1.0)
    @test collect(g.X)≈[1.0,0.0,0.0] atol=1e-12
    @test collect(g.Y)≈[0.0,1.0,0.0] atol=1e-12
    # UV bounds cover the unit square about its center.
    lo,hi=A.model.get_parametrization_bounds(2,f)
    @test lo≈[-0.5,-0.5] atol=1e-12
    @test hi≈[0.5,0.5] atol=1e-12
    @test A.model.get_value(2,f,[0.0,0.0])≈[0.5,0.5,0.0] atol=1e-12
    @test A.model.get_value(2,f,[0.5,0.5])≈[1.0,1.0,0.0] atol=1e-12
    @test sort!(A.model.get_boundary([(2,f)]))==[(1,1),(1,2),(1,3),(1,4)]
end

@testset "OCC add_plane_surface hole" begin
    occ=_occ_setup()
    po,eo=_occ_box(occ,0.0,0.0,2.0)
    pi,ei=_occ_box(occ,0.5,0.5,1.0)
    w_o=occ.add_wire(eo)
    w_i=occ.add_wire(ei)
    f=occ.add_plane_surface([w_o,w_i])
    occ.synchronize()
    m=_current_model()
    @test length(m.surfaces[f])==2
    # Hole wire stored REVERSED: reversed order with flipped signs.
    @test m.loops[m.surfaces[f][1]]==eo
    @test m.loops[m.surfaces[f][2]]==[-ei[4],-ei[3],-ei[2],-ei[1]]
    g=m.surface_geometry[f]
    # Every boundary edge (outer and hole) carries a pcurve.
    for tag in (eo...,ei...)
        @test haskey(g.pcurves,tag)
        pc=g.pcurves[tag].fwd
        @test hasproperty(pc,:lin2d)
    end
end

@testset "OCC add_plane_surface guards" begin
    occ=_occ_setup()
    p,e=_occ_square(occ)
    w=occ.add_wire(e)
    @test occursin("Could not create surface",_errtext(()->occ.add_plane_surface(Int[])))
    @test occursin("unknown OCC Wire[99]",_errtext(()->occ.add_plane_surface([99])))
    f=occ.add_plane_surface([w])
    @test occursin("already exists",_errtext(()->occ.add_plane_surface([w];tag=f)))
    # Non-planar wire: last edge leaves the z=0 plane.
    p5=occ.add_point(0.0,0.0,1.0)
    e5=occ.add_line(p[4],p5)
    wnp=occ.add_wire([e[1],e[2],e[3],e5])
    @test occursin("Could not create surface",_errtext(()->occ.add_plane_surface([wnp])))
    # A Bezier edge cannot carry an analytic plane pcurve (closed wire
    # needed so the check reaches the kind test).
    pb=[occ.add_point(x,5.0,0.0) for x in (0.0,0.5,1.0)]
    bz=occ.add_bezier(pb)
    bl=occ.add_line(pb[3],pb[1])
    wb=occ.add_wire([bz,bl])
    @test occursin("only line and circle",_errtext(()->occ.add_plane_surface([wb])))
    # An open wire is not a face boundary upstream either.
    wopen=occ.add_wire([e[1],e[2]])
    @test occursin("Could not create surface",_errtext(()->occ.add_plane_surface([wopen])))
    # Wire tag allocation is not consumed by the failures above.
    @test occ.get_max_tag(-1)==4
    @test occ.get_max_tag(2)==1
end

@testset "OCC add_plane_surface annulus" begin
    occ=_occ_setup()
    c1=occ.add_circle(0.0,0.0,0.0,2.0)
    c2=occ.add_circle(0.0,0.0,0.0,1.0)
    wd=occ.add_wire([c1]); wh=occ.add_wire([c2])
    fd=occ.add_plane_surface([wd])
    fa=occ.add_plane_surface([wd,wh])
    occ.synchronize()
    m=_current_model()
    g=m.surface_geometry[fa]
    @test g.axis==(0.0,0.0,1.0)
    @test collect(g.center)≈[0.0,0.0,0.0] atol=1e-12
    @test m.loops[m.surfaces[fa][1]]==[c1]
    @test m.loops[m.surfaces[fa][2]]==[-c2]
    @test g.pcurves[c1].fwd.circ2d[3]≈2.0
    @test g.pcurves[c2].fwd.circ2d[3]≈1.0
    lo,hi=A.model.get_parametrization_bounds(2,fa)
    @test lo≈[-2.0,-2.0] atol=1e-12
    @test hi≈[2.0,2.0] atol=1e-12
    # Disk on the single circle wire.
    gd=m.surface_geometry[fd]
    @test gd.axis==(0.0,0.0,1.0)
    # Pcurve round-trip: edge point -> UV -> surface point recovers the edge.
    for tag in (c1,c2)
        rec=m.curve_geometry[tag]; pc=g.pcurves[tag].fwd
        for t in range(rec.t0,rec.t1,length=7)
            p3=A.model.get_value(1,tag,[t])
            uv=M._occ_pcurve_d0(pc,t)
            p3b=A.model.get_value(2,fa,[uv[1],uv[2]])
            @test p3b≈p3 atol=1e-10
        end
    end
end

@testset "OCC add_rectangle" begin
    occ=_occ_setup()
    r=occ.add_rectangle(0.0,0.0,0.0,2.0,1.0)
    @test r==1
    occ.synchronize()
    m=_current_model()
    # Recursive upstream binding exposes all sub-entities.
    @test sort!(collect(occ.get_entities(0)))==[(0,t) for t in 1:4]
    @test sort!(collect(occ.get_entities(1)))==[(1,t) for t in 1:4]
    @test m.loops[only(m.surfaces[r])]==[1,2,3,4]
    @test m.surface_types[r]==:plane
    g=m.surface_geometry[r]
    @test g.axis==(0.0,0.0,1.0)
    lo,hi=A.model.get_parametrization_bounds(2,r)
    @test lo≈[-1.0,-0.5] atol=1e-12
    @test hi≈[1.0,0.5] atol=1e-12
    @test A.model.get_value(2,r,[0.0,0.0])≈[1.0,0.5,0.0] atol=1e-12
    @test A.model.get_normal(r,[0.0,0.0])≈[0.0,0.0,1.0] atol=1e-12
    @test sort!(A.model.get_boundary([(2,r)]))==[(1,1),(1,2),(1,3),(1,4)]
    # Negative dx/dy normalize to the same CCW cell (x∈[2,4], y∈[2,3]).
    rn=occ.add_rectangle(4.0,3.0,0.0,-2.0,-1.0)
    occ.synchronize()
    @test A.model.get_value(2,rn,[0.0,0.0])≈[3.0,2.5,0.0] atol=1e-12
end

@testset "OCC add_rectangle guards" begin
    occ=_occ_setup()
    @test occursin("zero width or height",_errtext(()->occ.add_rectangle(0,0,0,0,1)))
    @test occursin("zero width or height",_errtext(()->occ.add_rectangle(0,0,0,1,0)))
    @test occursin("must be finite",_errtext(()->occ.add_rectangle(0,0,0,NaN,1)))
    @test occursin("must be finite",_errtext(()->occ.add_rectangle(0,0,0,1,1;roundedRadius=NaN)))
    @test occ.get_entities(2)==[]
    @test occ.add_rectangle(0,0,0,1,1;tag=7)==7
    @test occursin("already exists",_errtext(()->occ.add_rectangle(0,0,0,1,1;tag=7)))
end

@testset "OCC add_rectangle rounded corners" begin
    occ=_occ_setup()
    r=occ.add_rectangle(0.0,0.0,0.0,4.0,2.0;roundedRadius=0.5)
    occ.synchronize()
    m=_current_model()
    @test length(occ.get_entities(0))==8
    @test length(occ.get_entities(1))==8
    loop=m.loops[only(m.surfaces[r])]
    @test loop==[1,2,3,4,5,6,7,8]
    @test [m.curve_types[t] for t in loop]==
        [:line,:circle,:line,:circle,:line,:circle,:line,:circle]
    g=m.surface_geometry[r]
    @test g.axis==(0.0,0.0,1.0)
    lo,hi=A.model.get_parametrization_bounds(2,r)
    @test lo≈[-2.0,-1.0] atol=1e-9
    @test hi≈[2.0,1.0] atol=1e-9
    # Corner arcs land their endpoints on the corner centers' rims and
    # their pcurves round-trip to the same 3-D points.
    for tag in loop
        rec=m.curve_geometry[tag]; pc=g.pcurves[tag].fwd
        for t in range(rec.t0,rec.t1,length=4)
            p3=A.model.get_value(1,tag,[t])
            uv=M._occ_pcurve_d0(pc,t)
            p3b=A.model.get_value(2,r,[uv[1],uv[2]])
            @test p3b≈p3 atol=1e-9
        end
    end
end

@testset "OCC wire/surface failure atomicity" begin
    occ=_occ_setup()
    p,e=_occ_square(occ)
    w=occ.add_wire(e)
    f=occ.add_plane_surface([w])
    before=((occ.get_max_tag(-1),occ.get_max_tag(0),
             occ.get_max_tag(1),occ.get_max_tag(2)),
            sort!(collect(occ.get_entities())))
    # Unknown edges, bad chains, open loops, non-planar wires, duplicate
    # tags, bad rectangles — none may bind or bump a maximum.
    for call in (()->occ.add_wire([99]),
                 ()->occ.add_wire([e[1],e[3]]),
                 ()->occ.add_wire(e;tag=w),
                 ()->occ.add_curve_loop([e[1],e[2]]),
                 ()->occ.add_plane_surface([99]),
                 ()->occ.add_plane_surface(Int[]),
                 ()->occ.add_plane_surface([w];tag=f),
                 ()->occ.add_rectangle(0,0,0,0,1),
                 ()->occ.add_rectangle(0,0,0,1,1;tag=f))
        @test !isempty(_errtext(call))
    end
    after=((occ.get_max_tag(-1),occ.get_max_tag(0),
            occ.get_max_tag(1),occ.get_max_tag(2)),
           sort!(collect(occ.get_entities())))
    @test before==after
end

@testset "OCC plane surface end-to-end mesh" begin
    occ=_occ_setup()
    # Square-with-square-hole face, lines only.
    p=[occ.add_point(x,y,0.0) for (x,y) in
       [(0,0),(4,0),(4,4),(0,4),(1,1),(3,1),(3,3),(1,3)]]
    e=[occ.add_line(p[a],p[b]) for (a,b) in
       [(1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5)]]
    wo=occ.add_wire(e[1:4]); wh=occ.add_wire(e[5:8])
    f=occ.add_plane_surface([wo,wh])
    occ.synchronize()
    A.mesh.generate(2)
    _,_,ntags=A.mesh.get_elements(2,f)
    face_nodes=unique(vcat(ntags...))
    all_tags,coords,_=A.mesh.get_nodes()
    pos=Dict(all_tags[i]=>(coords[3i-2],coords[3i-1]) for i in eachindex(all_tags))
    inside=count(n->begin q=pos[n]; 1.0<q[1]<3.0 && 1.0<q[2]<3.0 end,face_nodes)
    @test length(face_nodes)>10
    @test inside==0
end
