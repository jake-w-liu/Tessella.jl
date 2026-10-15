# `model.occ.offsetCurve` — `BRepOffsetAPI_MakeOffset` with `GeomAbs_Arc`
# joins. Every member of a planar wire is offset `offset` to the left of its
# travel direction (the wire normal follows the right-hand winding); adjacent
# offsets trim to a mutual intersection or bridge through a |offset|-radius
# arc centered on the shared vertex. Upstream `_multiBind` returns only the
# offset edges (dim-1 tags); the wire itself stays unbound.

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

function _unit_square(occ)
    p1=occ.add_point(0.0,0.0,0.0); p2=occ.add_point(1.0,0.0,0.0)
    p3=occ.add_point(1.0,1.0,0.0); p4=occ.add_point(0.0,1.0,0.0)
    es=[occ.add_line(p1,p2),occ.add_line(p2,p3),
        occ.add_line(p3,p4),occ.add_line(p4,p1)]
    return es,occ.add_curve_loop(es)
end

@testset "occ.offsetCurve closed polygon (trim joins)" begin
    occ=_occ_setup()
    es,lp=_unit_square(occ)
    out=occ.offset_curve(lp,0.1)
    # CCW wire: positive offset goes inward — inner square [0.1,0.9]^2.
    @test length(out)==4
    @test all(d==1 for (d,_) in out)
    occ.synchronize(); m=_current_model()
    spans=[(t0=m.curve_geometry[t].t0,t1=m.curve_geometry[t].t1)
        for (_,t) in out]
    @test all(s.t1-s.t0≈0.8 for s in spans)
    # Inward corners: (0.1,0.1) (0.9,0.1) (0.9,0.9) (0.1,0.9)
    g=m.curve_geometry[out[1][2]]
    @test g.occ===:line && collect(g.dir)≈[1.0,0.0,0.0]
    @test collect(g.origin)≈[0.0,0.1,0.0]
    # Edges chain head-to-tail through shared fresh vertices.
    for k in 1:4
        t=out[k][2]; tn=out[mod1(k+1,4)][2]
        @test m.curves[t][2]==m.curves[tn][1]
    end
    # All four new vertices are the inner corners.
    pts=[m.points[m.curves[t][1]] for (_,t) in out]
    for c in ((0.1,0.1,0.0),(0.9,0.1,0.0),(0.9,0.9,0.0),(0.1,0.9,0.0))
        @test any(p->isapprox(collect(p),collect(c);atol=1e-12),pts)
    end
end

@testset "occ.offsetCurve outward polygon (arc joins)" begin
    occ=_occ_setup()
    es,lp=_unit_square(occ)
    out=occ.offset_curve(lp,-0.1)
    @test length(out)==8                     # 4 lines + 4 corner arcs
    occ.synchronize(); m=_current_model()
    kinds=[m.curve_geometry[t].occ for (_,t) in out]
    @test kinds==[:line,:circle,:line,:circle,:line,:circle,:line,:circle]
    # Arc joins are centered on the original corners, radius |o|, quarter.
    for (_,t) in out[2:2:8]
        g=m.curve_geometry[t]
        @test g.r≈0.1 atol=1e-12
        @test g.t1-g.t0≈pi/2 atol=1e-12
    end
    @test collect(m.curve_geometry[out[2][2]].center)≈[1.0,0.0,0.0]
    @test collect(m.curve_geometry[out[8][2]].center)≈[0.0,0.0,0.0]
    # Tangent continuity: line end tangent == arc start tangent.
    for k in (2,4,6,8)
        ga=m.curve_geometry[out[k][2]]; gl=m.curve_geometry[out[k-1][2]]
        tl=M._occ_circle_derivative(ga,ga.t0)
        @test abs(tl[3])<1e-12 && abs(abs(tl[1])+abs(tl[2])-0.1)<1e-9
        @test gl.occ===:line
    end
end

@testset "occ.offsetCurve closed circle" begin
    occ=_occ_setup()
    c=occ.add_circle(0.0,0.0,0.0,2.0)
    lp=occ.add_curve_loop([c])
    out=occ.offset_curve(lp,0.5)
    @test length(out)==1
    occ.synchronize(); m=_current_model()
    g=m.curve_geometry[out[1][2]]
    @test g.occ===:circle
    @test g.r≈1.5
    @test g.t1-g.t0≈2pi
    # Vertices of a closed result coincide (single vertex both ends).
    @test m.curves[out[1][2]][1]==m.curves[out[1][2]][2]
end

@testset "occ.offsetCurve open wire and mixed members" begin
    occ=_occ_setup()
    # Open L: inward (+0.1) trims the corner — 2 lines, no arcs.
    p1=occ.add_point(0.0,0.0,0.0); p2=occ.add_point(1.0,0.0,0.0)
    p3=occ.add_point(1.0,1.0,0.0)
    e1=occ.add_line(p1,p2); e2=occ.add_line(p2,p3)
    w=occ.add_wire([e1,e2])
    out=occ.offset_curve(w,0.1)
    @test length(out)==2
    occ.synchronize(); m=_current_model()
    g1,g2=m.curve_geometry[out[1][2]],m.curve_geometry[out[2][2]]
    @test g1.occ===:line && collect(g1.origin)≈[0.0,0.1,0.0]
    @test g1.t1≈0.9
    @test g2.occ===:line && collect(g2.origin)≈[0.9,0.0,0.0]
    @test g2.t0≈0.1
    # D-shape: line + semicircle (center origin, bulge +x), inward offset
    # trims the arc, outward (-0.1) adds quarter-round caps at both ends.
    occ=_occ_setup()
    q1=occ.add_point(0.0,-1.0,0.0); q2=occ.add_point(0.0,1.0,0.0)
    qC=occ.add_point(0.0,0.0,0.0)
    ln=occ.add_line(q1,q2)
    arc=occ.add_circle_arc(q2,qC,q1)
    lp=occ.add_curve_loop([ln,arc])
    out=occ.offset_curve(lp,-0.1)
    @test length(out)==4
    occ.synchronize(); m=_current_model()
    kinds=[m.curve_geometry[t].occ for (_,t) in out]
    @test kinds==[:line,:circle,:circle,:circle]
    gm=m.curve_geometry[out[3][2]]
    @test gm.r≈1.1
    @test collect(m.curve_geometry[out[2][2]].center)≈[0.0,1.0,0.0]
    @test collect(m.curve_geometry[out[4][2]].center)≈[0.0,-1.0,0.0]
end

@testset "occ.offsetCurve reflex corner" begin
    occ=_occ_setup()
    # CCW L-shape with a reflex vertex at (0.5,0.5): inward offset rounds
    # the notch with a vertex-centered arc instead of a trim.
    ps=[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,0.0),
        (0.5,1.0,0.0),(0.5,0.5,0.0),(0.0,0.5,0.0)]
    pt=[occ.add_point(p...) for p in ps]
    es=[occ.add_line(pt[i],pt[mod1(i+1,6)]) for i in 1:6]
    lp=occ.add_curve_loop(es)
    out=occ.offset_curve(lp,0.1)
    @test length(out)==7
    occ.synchronize(); m=_current_model()
    kinds=[m.curve_geometry[t].occ for (_,t) in out]
    @test count(==(:circle),kinds)==1
    gi=findfirst(==(:circle),kinds)
    g=m.curve_geometry[out[gi][2]]
    @test collect(g.center)≈[0.5,0.5,0.0]
    @test g.r≈0.1
    @test g.t1-g.t0≈pi/2 atol=1e-12
end

@testset "occ.offsetCurve guards and atomicity" begin
    occ=_occ_setup()
    @test occursin("Unknown OpenCASCADE curve loop",
        _errtext(()->occ.offset_curve(999,0.1)))
    # Unsupported member kind.
    p1=occ.add_point(-2.0,0.0,0.0); p2=occ.add_point(0.0,0.0,0.0)
    pm=occ.add_point(2.0,0.0,0.0); p3=occ.add_point(0.0,1.0,0.0)
    el=occ.add_ellipse_arc(p1,p2,pm,p3)
    w=occ.add_wire([el])
    @test occursin("only line and circle",_errtext(()->occ.offset_curve(w,0.1)))
    # Non-planar wire.
    a=occ.add_point(0.0,0.0,0.0); b=occ.add_point(1.0,0.0,0.0)
    c=occ.add_point(1.0,0.0,1.0); d=occ.add_point(1.0,1.0,1.0)
    n1=occ.add_line(a,b); n2=occ.add_line(b,c); n3=occ.add_line(c,d)
    nw=occ.add_wire([n1,n2,n3])
    @test occursin("not planar",_errtext(()->occ.offset_curve(nw,0.1)))
    # Disconnected members (fillet separates their vertices afterward).
    x=occ.add_point(0.0,5.0,0.0); y=occ.add_point(1.0,5.0,0.0)
    z=occ.add_point(1.0,6.0,0.0)
    fa=occ.add_line(x,y); fb=occ.add_line(y,z)
    dw=occ.add_wire([fa,fb])
    occ.fillet_2d(fa,fb,0.2)
    @test occursin("not connected",_errtext(()->occ.offset_curve(dw,0.1)))
    # Offset swallowing a circle.
    sc=occ.add_circle(0.0,0.0,0.0,0.05)
    sl=occ.add_curve_loop([sc])
    @test occursin("offset swallows",_errtext(()->occ.offset_curve(sl,0.1)))
    # Failed calls bind nothing.
    occ.synchronize(); m=_current_model()
    kernel=m.occ_kernel
    @test !isempty(m.curves)
    nb=length(kernel.bindings)
    _errtext(()->occ.offset_curve(999,0.1))
    _errtext(()->occ.offset_curve(w,0.1))
    @test length(kernel.bindings)==nb
end
