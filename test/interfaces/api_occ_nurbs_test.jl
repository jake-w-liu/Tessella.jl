# `model.occ.convertToNURBS` — `BRepBuilderAPI_NurbsConvert` per bound
# entity. Edge geometry is replaced by an exact NURBS form under the same
# tag and topology; lines become degree-1 B-splines on the identical
# arc-length parametrization, and conics become quadratic-rational
# B-splines whose `_ModelOCCConicReparam` half-angle map preserves the
# public angle parametrization exactly — every stored pcurve stays valid
# bitwise. Vertices, wires, degenerate and `surface_iso` records pass
# through; surface/solid conversion is still staged separately.

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

_eval(g,t)=first(M._model_occ_curve_d2(g.curve,t))

@testset "occ.convertToNURBS line (exact arc-length)" begin
    occ=_occ_setup()
    a=occ.add_point(0.0,0.0,0.0); b=occ.add_point(3.0,4.0,0.0)
    e=occ.add_line(a,b)
    occ.convert_to_nurbs([(1,e)])
    occ.synchronize(); m=_current_model()
    g=m.curve_geometry[e]
    @test g.occ===:copied_curve
    @test m.curve_types[e]===:tensor_bspline
    @test g.t0==0.0 && g.t1≈5.0
    for t in (0.0,1.7,3.3,5.0)
        @test collect(_eval(g,t))≈[0.6t,0.8t,0.0] atol=1e-14
        _,d1,_=M._model_occ_curve_d2(g.curve,t)
        @test collect(d1)≈[0.6,0.8,0.0] atol=1e-14
    end
end

@testset "occ.convertToNURBS circle (reparam-exact)" begin
    occ=_occ_setup()
    c=occ.add_circle(1.0,2.0,0.0,2.0)
    occ.convert_to_nurbs([(1,c)])
    occ.synchronize(); m=_current_model()
    g=m.curve_geometry[c]
    @test g.occ===:copied_curve
    @test g.curve isa M._ModelOCCReparamCurve
    @test g.curve.curve.degree==2
    @test length(g.curve.curve.controls)==9   # 4 quarter spans
    for t in range(0.02,6.28,25)
        p,d1,d2=M._model_occ_curve_d2(g.curve,t)
        @test collect(p)≈[1+2cos(t),2+2sin(t),0.0] atol=1e-13
        @test collect(d1)≈[-2sin(t),2cos(t),0.0] atol=1e-12
        @test collect(d2)≈[-2cos(t),-2sin(t),0.0] atol=1e-11
    end
end

@testset "occ.convertToNURBS trimmed arc + ellipse" begin
    occ=_occ_setup()
    e=occ.add_ellipse(0.0,0.0,0.0,3.0,1.0)
    occ.convert_to_nurbs([(1,e)])
    occ.synchronize(); m=_current_model()
    g=m.curve_geometry[e]
    for t in range(0.05,6.2,40)
        p=_eval(g,t)
        @test (p[1]/3)^2+p[2]^2≈1.0 atol=1e-13
    end
    # trimmed arc — partial span, native-backed record
    occ2=A.model.occ
    p1=occ2.add_point(1.0,0.0,0.0); pc=occ2.add_point(0.0,0.0,0.0)
    p2=occ2.add_point(0.0,1.0,0.0)
    arc=occ2.add_circle_arc(p1,pc,p2)
    ga=let k=A._with_model(identity); k.occ_kernel.shapes[k.occ_kernel.bindings[(1,arc)]].geometry.geometry end
    occ2.convert_to_nurbs([(1,arc)])
    occ2.synchronize(); m2=_current_model()
    g2=m2.curve_geometry[arc]
    @test g2.occ===:copied_curve
    @test g2.t0≈ga.t0 && g2.t1≈ga.t1   # public parametrization preserved
    @test collect(_eval(g2,g2.t0))≈collect(M._occ_circle_point(ga,ga.t0)) atol=1e-13
    @test collect(_eval(g2,g2.t1))≈collect(M._occ_circle_point(ga,ga.t1)) atol=1e-13
    mid=(g2.t0+g2.t1)/2
    @test collect(_eval(g2,mid))≈collect(M._occ_circle_point(ga,mid)) atol=1e-13
end

@testset "occ.convertToNURBS bspline passthrough + bezier" begin
    occ=_occ_setup()
    ps=[occ.add_point(0.,0.,0.),occ.add_point(1.,0.,0.),
        occ.add_point(1.,1.,0.),occ.add_point(0.,1.,0.)]
    b=occ.add_bspline(ps)
    z=occ.add_bezier(ps)
    occ.convert_to_nurbs([(1,b),(1,z)])
    occ.synchronize(); m=_current_model()
    @test m.curve_geometry[b].occ===:copied_curve
    @test m.curve_geometry[b].curve.family===:bspline
    @test m.curve_types[z]===:tensor_bspline
    # bspline family re-normalizes :bezier → :bspline
    @test m.curve_geometry[z].curve.family===:bspline
end

@testset "occ.convertToNURBS keeps pcurves valid on faces" begin
    occ=_occ_setup()
    pts=[occ.add_point(0.,0.,0.),occ.add_point(1.,0.,0.),
         occ.add_point(1.,1.,0.),occ.add_point(0.,1.,0.)]
    es=[occ.add_line(pts[1],pts[2]),occ.add_line(pts[2],pts[3]),
        occ.add_line(pts[3],pts[4]),occ.add_line(pts[4],pts[1])]
    lp=occ.add_curve_loop(es)
    f=occ.add_plane_surface(lp)
    occ.convert_to_nurbs([(1,t) for t in es])
    occ.convert_to_nurbs([(-1,lp)])      # wire: member edges already converted
    occ.synchronize(); m=_current_model()
    g=m.surface_geometry[f]
    for (etag,ent) in g.pcurves
        pc=ent.fwd; ge=m.curve_geometry[etag]
        for t in range(ge.t0+1e-9,ge.t1-1e-9,5)
            uv=M._occ_pcurve_d0(pc,t)
            ps=M._occ_surface_point(g,uv...)
            pe=_eval(ge,t)
            @test collect(ps)≈collect(pe) atol=1e-12
        end
    end
end

@testset "occ.convertToNURBS downstream queries" begin
    occ=_occ_setup()
    c=occ.add_circle(0.0,0.0,0.0,1.0)
    occ.convert_to_nurbs([(1,c)])
    @test occ.get_mass(1,c)≈2π atol=1e-10
    occ.translate([(1,c)],1.0,2.0,3.0)
    occ.synchronize(); m=_current_model()
    g=m.curve_geometry[c]
    @test g.curve isa M._ModelOCCReparamCurve
    @test collect(_eval(g,0.0))≈[2.0,2.0,3.0] atol=1e-13
    # get_closest_entities materializes through the converted evaluator
    pts,_=occ.get_closest_entities(1.0,2.0,3.0,[(1,c)],1)
    @test pts[1]==(1,c)
end

@testset "occ.convertToNURBS atomicity + rejections" begin
    occ=_occ_setup()
    a=occ.add_point(0.0,0.0,0.0); b=occ.add_point(1.0,0.0,0.0)
    e=occ.add_line(a,b)
    kernel=_current_model().occ_kernel
    nshapes=length(kernel.shapes)
    @test occursin(
        "OCC.convertToNURBS: Unknown OpenCASCADE entity of dimension 1 with tag 999",
        _errtext(()->occ.convert_to_nurbs([(1,999)])))
    @test length(kernel.shapes)==nshapes
    # dim 0 vertex: accepted no-op (upstream converts trivially)
    occ.convert_to_nurbs([(0,a)])
    # unimplemented dims reject loudly before any mutation
    c=occ.add_point(0.0,1.0,0.0); d=occ.add_point(1.0,1.0,0.0)
    w=occ.add_curve_loop([occ.add_line(a,b),occ.add_line(b,d),
                          occ.add_line(d,c),occ.add_line(c,a)])
    f=occ.add_plane_surface(w)
    @test occursin("not implemented",
        _errtext(()->occ.convert_to_nurbs([(2,f)])))
    @test length(kernel.shapes)>nshapes   # face publication succeeded
    n2=length(kernel.shapes)
    _errtext(()->occ.convert_to_nurbs([(2,f)]))
    @test length(kernel.shapes)==n2     # rejected conversion adds nothing
    # reconversion is idempotent
    occ.convert_to_nurbs([(1,e)])
    occ.convert_to_nurbs([(1,e)])
    occ.synchronize(); m=_current_model()
    @test m.curve_geometry[e].occ===:copied_curve
    @test collect(_eval(m.curve_geometry[e],0.5))≈[0.5,0.0,0.0] atol=1e-14
end
