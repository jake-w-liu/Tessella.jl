# `model.occ.fillet2D`/`chamfer2D` — `ChFi2d_FilletAPI`/`ChFi2d_ChamferAPI`
# corner trims. Two edges sharing a coincident endpoint are cut back to the
# fillet tangency / chamfer distance points and a fresh joining edge binds
# under `tag` (auto otherwise); the trimmed edges keep their tags and the
# orphaned corner vertex stays bound, matching `OCC_Internals::fillet2D`/
# `chamfer2D` on Gmsh 4.15.2 (temporary wire/face gates, `pointTag` result
# selection, `reverse` side flip, recursive rebind conventions).

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

@testset "occ.fillet2D line-line" begin
    occ=_occ_setup()
    p1=occ.add_point(0,0,0); p2=occ.add_point(1,0,0); p3=occ.add_point(1,1,0)
    e1=occ.add_line(p1,p2); e2=occ.add_line(p2,p3)
    ftag=occ.fillet_2d(e1,e2,0.2)
    @test ftag==3
    occ.synchronize()
    m=_current_model()
    g1,g2,gf=m.curve_geometry[e1],m.curve_geometry[e2],m.curve_geometry[ftag]
    @test g1.t1≈0.8; @test g2.t0≈0.2
    @test gf.occ===:circle
    @test collect(gf.center)≈[0.8,0.2,0.0]
    @test gf.r≈0.2
    @test gf.t1-gf.t0≈pi/2
    @test collect(M._occ_circle_point(gf,gf.t0))≈[0.8,0,0]
    @test collect(M._occ_circle_point(gf,gf.t1))≈[1.0,0.2,0]
    @test occ.get_mass(1,ftag)≈0.2*pi/2
    # tangent continuity with e1's travel direction
    dT=M._occ_circle_derivative(gf,gf.t0)
    @test dT[1]>0 && abs(dT[2])<1e-12 && abs(dT[3])<1e-12
    # shared vertices: trimmed e1 ends where the fillet starts; fillet ends
    # where trimmed e2 starts; the orphaned corner vertex stays bound.
    @test m.curves[e1][2]==m.curves[ftag][1]
    @test m.curves[ftag][2]==m.curves[e2][1]
    @test haskey(m.points,p2)
    @test length(m.points)==5
end

@testset "occ.fillet2D concave and end-end corners" begin
    occ=_occ_setup()
    # Right turn: wedge is the lower-left sector.
    a=occ.add_point(0,0,0); v=occ.add_point(1,0,0); b=occ.add_point(1,-1,0)
    e1=occ.add_line(a,v); e2=occ.add_line(v,b)
    ftag=occ.fillet_2d(e1,e2,0.2)
    occ.synchronize(); m=_current_model()
    @test collect(m.curve_geometry[ftag].center)≈[0.8,-0.2,0.0]
    @test m.curve_geometry[ftag].t1-m.curve_geometry[ftag].t0≈pi/2
    # Both edges ENDING at the corner: e4 travels -y into V.
    v2=occ.add_point(1,0,0); d=occ.add_point(1,1,0); a2=occ.add_point(0,0,0)
    e3=occ.add_line(a2,v2); e4=occ.add_line(d,v2)
    ft2=occ.fillet_2d(e3,e4,0.2)
    occ.synchronize(); m=_current_model()
    @test collect(m.curve_geometry[ft2].center)≈[0.8,0.2,0.0]
    @test m.curve_geometry[e4].t1≈0.8   # e4 trimmed at its END
end

@testset "occ.fillet2D line-circle (native record)" begin
    occ=_occ_setup()
    pa=occ.add_point(0,0,0); pv=occ.add_point(1,0,0)
    pc=occ.add_point(0,0,0); pend=occ.add_point(0,1,0)
    e1=occ.add_line(pa,pv)
    e2=occ.add_circle_arc(pv,pc,pend)   # arc (1,0)->(0,1) about origin, r=1
    r=0.2
    ftag=occ.fillet_2d(e1,e2,r)
    occ.synchronize(); m=_current_model()
    gf=m.curve_geometry[ftag]
    # Analytic solution: C on y=r offset of the line, |C|=R-r internally.
    @test collect(gf.center)≈[sqrt(0.6),0.2,0.0]
    @test collect(M._occ_circle_point(gf,gf.t0))≈[sqrt(0.6),0.0,0.0]
    T2=M._occ_circle_point(gf,gf.t1)
    @test abs(M._occ_modulus(M._arc_sub(T2,gf.center))-r)<1e-12
    @test collect(T2)≈[sqrt(0.6)/0.8,0.25,0.0]
    # The trimmed circle arc still evaluates through its native frame.
    g2=m.curve_geometry[e2]
    @test collect(M._occ_circle_point(g2,g2.t0))≈collect(T2)
    @test collect(M._occ_circle_point(g2,g2.t1))≈[0,1,0]
end

@testset "occ.fillet2D circle-circle with pointTag" begin
    occ=_occ_setup()
    # circle1 c=(0,0) r=1 and circle2 c=(0.8,0) r=0.5 meet transversally at
    # V±=(0.86875,±0.49525); arcs run V+ -> V- through the right of each.
    x1=0.86875; y1=sqrt(1-x1^2)
    v1=occ.add_point(x1,y1,0); v2=occ.add_point(x1,-y1,0)
    mid1=occ.add_point(1,0,0); mid2=occ.add_point(1.3,0,0)
    arc1=occ.add_circle_arc(v1,mid1,v2;center=false)
    arc2=occ.add_circle_arc(v1,mid2,v2;center=false)
    r=0.1
    ftag=occ.fillet_2d(arc1,arc2,r;pointTag=v1)
    occ.synchronize(); m=_current_model()
    gf=m.curve_geometry[ftag]
    T1=M._occ_circle_point(gf,gf.t0); T2=M._occ_circle_point(gf,gf.t1)
    @test abs(M._occ_modulus(M._arc_sub(T1,gf.center))-r)<1e-12
    @test abs(M._occ_modulus(M._arc_sub(T2,gf.center))-r)<1e-12
    # The v1 corner was cut (v1 side trimmed, v2 side kept).
    g1=m.curve_geometry[arc1]
    dT1v1=M._occ_modulus(M._arc_sub(T1,(x1,y1,0.0)))
    dT1v2=M._occ_modulus(M._arc_sub(T1,(x1,-y1,0.0)))
    @test dT1v1<dT1v2
    @test collect(M._occ_circle_point(g1,g1.t1))≈[x1,-y1,0]
    # Tangent continuity on both sides.
    tT1=M._occ_2df_param(m.curve_geometry[arc1],T1)
    d1=M._occ_circle_derivative(m.curve_geometry[arc1],tT1)
    dF=M._occ_circle_derivative(gf,gf.t0)
    @test M._occ_modulus(M._occ_cross(d1,dF))<1e-10*max(M._occ_modulus(d1),1.0)
end

@testset "occ.fillet2D closed-circle seam corner" begin
    occ=_occ_setup()
    circ=occ.add_circle(0.5,0,0,0.5)
    occ.synchronize(); m=_current_model()
    gc=m.curve_geometry[circ]
    P0=M._occ_circle_point(gc,gc.t0)
    svtag=only(filter(t->M._points_close(m.points[t],P0,1e-7),
                      [t for (d,t) in occ.get_entities(0)]))
    ln=occ.add_line(svtag,occ.add_point(P0[1]-1,P0[2],P0[3]))
    r=0.1
    ftag=occ.fillet_2d(circ,ln,r)
    occ.synchronize(); m=_current_model()
    gf=m.curve_geometry[ftag]
    T1=M._occ_circle_point(gf,gf.t0)
    @test abs(M._occ_modulus(M._arc_sub(T1,gf.center))-r)<1e-12
    # Tangency point is on the circle and near the seam corner.
    gcl=m.curve_geometry[circ]
    tT=M._occ_2df_param(gcl,T1)
    @test M._points_close(M._occ_circle_point(gcl,tT),T1,1e-9)
    @test M._occ_modulus(M._arc_sub(T1,P0))<0.5
    # The circle kept its tag with a narrowed span.
    @test m.curve_geometry[circ].t1-m.curve_geometry[circ].t0 < 2pi-0.05
end

@testset "occ.chamfer2D" begin
    occ=_occ_setup()
    q1=occ.add_point(0,10,0); q2=occ.add_point(1,10,0); q3=occ.add_point(1,11,0)
    f1=occ.add_line(q1,q2); f2=occ.add_line(q2,q3)
    ctag=occ.chamfer_2d(f1,f2,0.3,0.25)
    @test ctag==3
    occ.synchronize(); m=_current_model()
    gc=m.curve_geometry[ctag]
    @test gc.occ===:line
    @test collect(M._occ_line_point(gc,gc.t0))≈[0.7,10,0]
    @test collect(M._occ_line_point(gc,gc.t1))≈[1.0,10.25,0]
    @test occ.get_mass(1,ctag)≈hypot(0.3,0.25)
    @test m.curves[f1][2]==m.curves[ctag][1]
    @test m.curves[ctag][2]==m.curves[f2][1]
    # Chamfer on a circle edge: distance converts to angle/R.
    p1=occ.add_point(0,20,0); pv=occ.add_point(1,20,0)
    pc=occ.add_point(0,20,0); pe=occ.add_point(0,21,0)
    h1=occ.add_line(p1,pv)
    h2=occ.add_circle_arc(pv,pc,pe)   # arc (1,20)->(0,21) about (0,20), r=1
    ct2=occ.chamfer_2d(h1,h2,0.4,0.5)
    occ.synchronize(); m=_current_model()
    @test collect(M._occ_line_point(m.curve_geometry[h1],
        m.curve_geometry[h1].t1))≈[0.6,20,0]
    @test collect(M._occ_circle_point(m.curve_geometry[h2],
        m.curve_geometry[h2].t0))≈[cos(0.5),20+sin(0.5),0]
end

@testset "occ.fillet2D/chamfer2D guards and atomicity" begin
    occ=_occ_setup()
    a=occ.add_point(0,0,0); v=occ.add_point(1,0,0); b=occ.add_point(1,1,0)
    e1=occ.add_line(a,v); e2=occ.add_line(v,b)
    # A disconnected edge and a bspline (ChFi2d scope: lines/circles).
    ua=occ.add_point(0,5,0); ub=occ.add_point(1,5,0)
    e3=occ.add_line(ua,ub)
    bs=occ.add_bspline([v,occ.add_point(1.5,0.5,0),
                        occ.add_point(2,1,0)])
    nb_curve=length(occ.get_entities(1)); nb_vert=length(occ.get_entities(0))
    for (name,f) in (
        ("radius zero",()->occ.fillet_2d(e1,e2,0.0)),
        ("radius negative",()->occ.fillet_2d(e1,e2,-1.0)),
        ("radius NaN",()->occ.fillet_2d(e1,e2,NaN)),
        ("radius exceeds",()->occ.fillet_2d(e1,e2,5.0)),
        ("unconnected",()->occ.fillet_2d(e1,e3,0.2)),
        ("unknown edge",()->occ.fillet_2d(e1,999,0.2)),
        ("tag clash",()->occ.fillet_2d(e1,e2,0.2;tag=e1)),
        ("bspline edge",()->occ.fillet_2d(bs,e2,0.2)),
        ("unknown pointTag",()->occ.fillet_2d(e1,e2,0.2;pointTag=999)),
        ("reverse reflex",()->occ.fillet_2d(e1,e2,0.2;reverse=true)),
        ("chamfer exceeds",()->occ.chamfer_2d(e1,e2,2.0,0.5)),
        ("chamfer zero",()->occ.chamfer_2d(e1,e2,0.0,0.5)),
        ("chamfer NaN",()->occ.chamfer_2d(e1,e2,0.2,NaN)))
        msg=_errtext(f)
        @test msg!=""
    end
    @test length(occ.get_entities(1))==nb_curve
    @test length(occ.get_entities(0))==nb_vert
    # Requested tags and return values.
    @test occ.fillet_2d(e1,e2,0.2;tag=42)==42
    @test (1,42) in occ.get_entities(1)
    c=occ.add_point(3,0,0); d2=occ.add_point(3,1,0)
    e5=occ.add_line(c,occ.add_point(4,0,0)); e6=occ.add_line(occ.add_point(4,0,0),d2)
    # e5 took the next auto tag (43) after the explicit 42 bumped the maximum.
    @test e5==43
    @test occ.chamfer_2d(e5,e6,0.1,0.1;tag=77)==77
    @test (1,77) in occ.get_entities(1)
    occ.synchronize()
    m=_current_model()
    @test haskey(m.curve_geometry,42)
    @test haskey(m.curve_geometry,77)
end
