# `model.occ.convertToNURBS` — `BRepBuilderAPI_NurbsConvert` per bound
# entity. Edge geometry is replaced by an exact NURBS form under the same
# tag and topology; lines become degree-1 B-splines on the identical
# arc-length parametrization, and conics become quadratic-rational
# B-splines whose `_ModelOCCConicReparam` half-angle map preserves the
# public angle parametrization exactly — every stored pcurve stays valid
# bitwise. Faces keep their analytic record as metadata and gain a NURBS
# `patch` plus per-axis `reparam`; evaluation funnels through
# `_model_nurbs_surface_d2` while external (u,v) stay identical. Vertices,
# wires, degenerate and `surface_iso` records pass through; solids convert
# recursively preserving their primitive encodings.

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
    # reconversion is idempotent
    occ.convert_to_nurbs([(1,e)])
    occ.convert_to_nurbs([(1,e)])
    occ.synchronize(); m=_current_model()
    @test m.curve_geometry[e].occ===:copied_curve
    @test collect(_eval(m.curve_geometry[e],0.5))≈[0.5,0.0,0.0] atol=1e-14
end

_surf_eval(g,u,v)=M._occ_surface_point(g,u,v)
_surf_analytic(g)=Base.structdiff(g,(patch=nothing,reparam=nothing))

@testset "occ.convertToNURBS plane + quadrics (exact reparam)" begin
    occ=_occ_setup()
    # plane face → bilinear degree-1 patch
    a=occ.add_point(0.0,0.0,0.0); b=occ.add_point(2.0,0.0,0.0)
    c=occ.add_point(2.0,3.0,0.0); d=occ.add_point(0.0,3.0,0.0)
    w=occ.add_curve_loop([occ.add_line(a,b),occ.add_line(b,c),
                          occ.add_line(c,d),occ.add_line(d,a)])
    f=occ.add_plane_surface(w)
    occ.convert_to_nurbs([(2,f)])
    occ.synchronize(); m=_current_model()
    g=m.surface_geometry[f]
    @test hasproperty(g,:patch)
    @test g.patch.surface.degree_u==1 && g.patch.surface.degree_v==1
    ga=_surf_analytic(g)
    lo,hi=M._occ_surface_bounds(g)
    for i in 1:4, j in 1:4
        u=Float64(lo[1]+(hi[1]-lo[1])*(i-0.5)/4)
        v=Float64(lo[2]+(hi[2]-lo[2])*(j-0.5)/4)
        @test collect(_surf_eval(g,u,v))≈
              collect(M._occ_surface_point(ga,u,v)) atol=1e-13
    end
    # cylinder lateral: quadratic-rational u-axis + linear v
    cyl=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,4.0,1.5)
    occ.convert_to_nurbs([(3,cyl)])
    occ.synchronize(); m=_current_model()
    lat=only([t for t in keys(m.surface_geometry) if
              m.surface_geometry[t].occ===:cylinder])
    g=m.surface_geometry[lat]
    @test hasproperty(g,:patch)
    @test g.patch.surface.degree_u==2
    ga=_surf_analytic(g)
    for u in (0.37,1.11,2.83,4.4,5.9), v in (0.4,1.7,3.2)
        pu=Float64(u);pv=Float64(v)
        p=_surf_eval(g,pu,pv)
        pa=M._occ_surface_point(ga,pu,pv)
        @test collect(p)≈collect(pa) atol=1e-12
        _,du,dv=M._occ_surface_d1(g,pu,pv)
        _,dua,dva=M._occ_surface_d1(ga,pu,pv)
        @test collect(du)≈collect(dua) atol=1e-10
        @test collect(dv)≈collect(dva) atol=1e-10
    end
    # u knot locations must also be exact (left-limit at shared knots)
    for u in (π/2,π,3π/2), v in (1.0,2.5)
        pu=Float64(u);pv=Float64(v)
        _,_,_,duu,dvv,duv=M._occ_surface_d2(g,pu,pv)
        _,_,_,duua,dvva,duva=M._occ_surface_d2(ga,pu,pv)
        @test collect(duu)≈collect(duua) atol=1e-10
        @test collect(dvv)≈collect(dvva) atol=1e-10
        @test collect(duv)≈collect(duva) atol=1e-10
    end
    # volume measure through the materialized encoding path
    @test occ.get_mass(3,cyl)≈π*1.5^2*4 atol=1e-9
end

@testset "occ.convertToNURBS sphere/torus/cone surfaces" begin
    occ=_occ_setup()
    sph=occ.add_sphere(0.0,0.0,0.0,2.0)
    occ.synchronize(); m=_current_model()
    sf=only([t for t in keys(m.surface_geometry) if
             m.surface_geometry[t].occ===:sphere])
    occ.convert_to_nurbs([(2,sf)])
    occ.synchronize(); m=_current_model()
    g=m.surface_geometry[sf]
    @test hasproperty(g,:patch) && g.reparam.u!==nothing &&
          g.reparam.v!==nothing
    ga=_surf_analytic(g)
    for u in (0.4,1.9,3.3,5.1), v in (-1.2,-0.4,0.3,1.25)
        pu=Float64(u);pv=Float64(v)
        @test collect(_surf_eval(g,pu,pv))≈
              collect(M._occ_surface_point(ga,pu,pv)) atol=1e-12
        # pole-adjacent v stays on the sphere
        @test abs(M._occ_modulus(_surf_eval(g,pu,pv))-2.0)<1e-12
    end
    # pcurve lift: every boundary edge's uv-trace still lands on the seam
    for loop in m.surfaces[sf], etag in m.loops[loop]
        ent=M._occ_face_pcurve(m,g,etag)
        ent===nothing && continue
        pc,t0,t1=ent
        for k in 1:5
            t=t0+(t1-t0)*(k-0.5)/5
            e3=M._periodic_curve_point(m,abs(etag),t,"test")
            uv=M._occ_pcurve_d0(pc,t)
            @test collect(M._occ_surface_point(g,uv[1],uv[2]))≈
                  collect(e3) atol=1e-10
        end
    end
    # torus — product of two circle nets
    occ=_occ_setup()
    tor=occ.add_torus(0.0,0.0,0.0,3.0,1.0)
    occ.synchronize(); m=_current_model()
    tf=only([t for t in keys(m.surface_geometry) if
             m.surface_geometry[t].occ===:torus])
    occ.convert_to_nurbs([(2,tf)])
    occ.synchronize(); m=_current_model()
    g=m.surface_geometry[tf];ga=_surf_analytic(g)
    for u in (0.7,2.2,4.0,5.5), v in (0.3,1.8,3.9,5.7)
        pu=Float64(u);pv=Float64(v)
        @test collect(_surf_eval(g,pu,pv))≈
              collect(M._occ_surface_point(ga,pu,pv)) atol=1e-12
    end
    # cone — slant-length v parametrization preserved
    occ=_occ_setup()
    cn=occ.add_cone(0.0,0.0,0.0,0.0,0.0,3.0,2.0,1.0)
    occ.synchronize(); m=_current_model()
    cf=only([t for t in keys(m.surface_geometry) if
             m.surface_geometry[t].occ===:cone])
    occ.convert_to_nurbs([(2,cf)])
    occ.synchronize(); m=_current_model()
    g=m.surface_geometry[cf];ga=_surf_analytic(g)
    lo,hi=M._occ_surface_bounds(g)
    for u in (0.5,2.6,4.8), v in (0.4,1.5,2.8)
        pu=Float64(u);pv=Float64(v)
        @test collect(_surf_eval(g,pu,pv))≈
              collect(M._occ_surface_point(ga,pu,pv)) atol=1e-12
    end
    @test hi[2]≈sqrt(9.0+1.0) atol=1e-12
end

@testset "occ.convertToNURBS transformed + projected converted face" begin
    occ=_occ_setup()
    cyl=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,4.0,1.5)
    occ.synchronize(); m=_current_model()
    lat=only([t for t in keys(m.surface_geometry) if
              m.surface_geometry[t].occ===:cylinder])
    occ.convert_to_nurbs([(3,cyl)])
    occ.rotate([(3,cyl)],0.0,0.0,0.0,0.0,0.0,1.0,π/6)
    occ.translate([(3,cyl)],1.0,2.0,3.0)
    occ.synchronize(); m=_current_model()
    g=m.surface_geometry[lat]
    lo,hi=M._occ_surface_bounds(g)
    for i in 1:5, j in 1:5
        u=Float64(lo[1]+(hi[1]-lo[1])*(i-0.5)/5)
        v=Float64(lo[2]+(hi[2]-lo[2])*(j-0.5)/5)
        p=M._occ_surface_point(g,u,v)
        @test hypot(p[1]-1.0,p[2]-2.0)≈1.5 atol=1e-11
    end
    @test occ.get_mass(3,cyl)≈π*1.5^2*4 atol=1e-9
    # projection recovers the external (angle, height) parametrization
    tgt=collect(M._occ_surface_point(g,1.1,2.5))
    u,v,xyz=M._model_nurbs_surface_project(g,tuple(tgt...))
    @test u≈1.1 atol=1e-9
    @test v≈2.5 atol=1e-9
    # reflection rejects loudly
    @test occursin("not implemented",
        _errtext(()->occ.mirror([(3,cyl)],1.0,0.0,0.0,0.0)))
end

@testset "occ.convertToNURBS converted face still meshes" begin
    occ=_occ_setup()
    cyl=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,4.0,1.5)
    occ.convert_to_nurbs([(3,cyl)])
    occ.synchronize()
    A.mesh.generate(3)
    ets,etags,_=A.mesh.get_elements(3,cyl)
    @test ets==Int32[4]
    @test length(etags[1])>0
    m=_current_model()
    lat=only([t for t in keys(m.surface_geometry) if
              m.surface_geometry[t].occ===:cylinder])
    @test hasproperty(m.surface_geometry[lat],:patch)
    # the primitive encoding survives — members still satisfy it exactly
    @test haskey(m.cylinders,cyl)
end
