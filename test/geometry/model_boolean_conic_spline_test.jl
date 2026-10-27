# Boolean engine support for ellipse and B-spline edge records —
# `_brep_edge_curve` accepts the `occ=:ellipse` family and imported
# `occ=:copied_curve` NURBS evaluators, pcurve machinery covers
# `ellips2d`/`spl2d`, and the mixed-dimensional engine materializes
# ellipse-arc result edges. Genuinely unrepresentable records and
# unsupported surface pairs still reject atomically by name.

using Test
using Tessella
const A=Tessella.API
const M=Tessella.Model

const _FIX=joinpath(@__DIR__,"..","fixtures")

_errtext(f)=(try;f();catch err;return sprint(showerror,err);end;"")
_ueq(a,b,rt=1e-12)=all(abs.(collect(a).-collect(b)).<=rt*max(1.0,maximum(abs,collect(b))))

@testset "Boolean edge extraction — ellipse and copied curves" begin
    A.initialize()
    try
        occ=A.model.occ
        e=occ.add_ellipse(0.0,0.0,0.0,1.4,0.8)
        pa=[occ.add_point(x,y,0.0) for (x,y) in
            ((0.0,0.0),(1.0,0.3),(2.0,-0.2),(3.0,0.5))]
        bs=occ.add_bspline(pa)
        occ.synchronize()
        m=A._with_model(identity)

        ce,t0,t1=M._brep_edge_curve(m,e)
        @test ce.kind===:ellipse
        @test ce.rx≈1.4 && ce.ry≈0.8
        @test _ueq(ce.c,(0.0,0.0,0.0))
        @test t0≈0.0 && t1≈2π
        @test _ueq(M._brep_curve_eval(ce,0.0),(1.4,0.0,0.0))
        @test _ueq(M._brep_curve_eval(ce,π/2),(0.0,0.8,0.0))

        cb,bt0,bt1=M._brep_edge_curve(m,bs)
        @test cb.kind===:bspline
        @test M._brep_curve_eval(cb,0.5*(bt0+bt1)) ==
            M._model_occ_curve_d2(cb.curve,0.5*(bt0+bt1))[1]

        # a genuinely unrepresentable record keeps the named rejection
        ctag=maximum(keys(m.curve_geometry))+1
        m.curves[ctag]=m.curves[e]
        m.curve_geometry[ctag]=(occ=:bezier,)
        err=_errtext(()->M._brep_edge_curve(m,ctag))
        @test occursin("_brep_edge_curve: unsupported curve record",err)
    finally
        A.finalize()
    end
end

@testset "Boolean pcurve machinery — ellips2d and spl2d" begin
    pc=(ellips2d=((0.0,0.0),(1.0,0.0),1.4,0.8),)
    tol=1e-9
    ev(t)=M._brep_pc_eval(pc,t)
    @test _ueq(ev(0.0),(1.4,0.0))
    @test ev(π/2)[2]≈0.8

    # ray parity: inside odd, outside even (an outside ray hits a closed
    # conic twice)
    for (pu,pv,inside) in ((0.0,0.0,true),(2.0,0.0,false),(0.0,0.5,true),
                           (0.0,0.9,false),(-2.0,0.3,false))
        n=M._brep_pc_ray_crossings(pc,0.0,2π,pu,pv,tol)
        @test isodd(n)==inside
    end

    # section-line × ellipse intersections invert onto the line's parameter
    sec=(lin2d=((1.0,0.0),(0.0,1.0)),)
    its=M._brep_pc_intersect(sec,pc,0.0,2π,tol)
    ey=0.8*sqrt(1-1/1.96)
    @test length(its)==2 && all(abs.(sort(abs.(its)).-ey).<1e-6)

    # Green's-theorem arc area: the full period is the ellipse area
    @test M._brep_pc_arc_area(pc,0.0,2π)≈π*1.4*0.8
    @test _ueq(M._brep_pc_urange(pc,0.0,2π),(-1.4,1.4))

    # inversion, flip, affine reparam
    s=M._brep_pc_param(pc,(1.0,ey),0.0,2π,tol)
    @test s!==nothing && abs(ev(s)[1]-1.0)<1e-7 && abs(ev(s)[2]-ey)<1e-7
    # near/param must use the ABSOLUTE conic level — a point deep inside
    # the disk is neither near the locus nor invertible onto it (the
    # signed level is negative there and once passed the band check)
    circ=(circ2d=((0.0,0.0),(1.0,0.0),1.0),)
    @test !M._brep_pc_near(circ,0.0,2π,(0.0,0.0),1e-7)
    @test M._brep_pc_param(circ,(0.0,0.0),0.0,2π,tol)===nothing
    @test !M._brep_pc_near(pc,0.0,2π,(0.0,0.0),1e-7)
    @test M._brep_pc_param(pc,(0.0,0.0),0.0,2π,tol)===nothing
    pf=M._brep_flip_pc(pc,"test")
    @test _ueq(M._brep_pc_eval(pf,-0.7),ev(0.7))
    pq,neg=M._brep_pc_reparam(pc,1.0,π/2,tol,"test")
    @test !neg && _ueq(M._brep_pc_eval(pq,0.4),ev(0.4+π/2))
    # an axis-mismatched phase map falls back to a spl2d carrier
    pb,negb=M._brep_pc_reparam(pc,1.0,0.3,tol,"test")
    @test !negb && hasproperty(pb,:spl2d) &&
          _ueq(M._brep_pc_eval(pb,0.4),ev(0.7),3e-4)

    # spl2d pcurve through the same machinery (PCHIP fit of the ellipse)
    ts=[2π*i/16 for i in 0:16]
    sp=(spl2d=M._pchip_bspline_2d(ts,[ev(t) for t in ts]),)
    @test isodd(M._brep_pc_ray_crossings(sp,0.0,2π,0.0,0.0,tol))
    @test M._brep_pc_ray_crossings(sp,0.0,2π,2.0,0.0,tol)==0
    @test abs(M._brep_pc_arc_area(sp,0.0,2π)-π*1.4*0.8)<1e-2
    @test M._brep_pc_urange(sp,0.0,2π)[2]≈1.4 atol=1e-6
    spf=M._brep_flip_pc(sp,"test")
    @test _ueq(M._brep_pc_eval(spf,-1.0),M._brep_pc_eval(sp,1.0))
    spp=M._brep_pc_param(sp,M._brep_pc_eval(sp,1.0),0.0,2π,tol)
    @test spp≈1.0 atol=1e-7
end

@testset "occ.intersect — elliptical disk clipped by a box slab" begin
    A.initialize()
    try
        occ=A.model.occ
        # 1.4×0.8 disk at z=0 ∩ slab |x|≤1: the result boundary carries
        # ellipse arcs; the enclosed area is analytic.
        d=occ.add_disk(0.0,0.0,0.0,1.4,0.8)
        occ.add_box(-1.0,-2.0,-1.0,2.0,4.0,2.0)
        occ.synchronize()
        out,_=occ.intersect([(2,d)],[(3,1)])
        occ.synchronize()
        @test out==[(2,d)]
        rx,ry=1.4,0.8
        expected=2*ry*(sqrt(1-1/rx^2)+rx*asin(1/rx))
        @test occ.get_mass(2,d)≈expected rtol=1e-9
        # the result face wires still reference ellipse edge pieces
        m=A._with_model(identity)
        stags=A.model.get_entities(2)
        @test length(stags)==1
        sg=m.surface_geometry[stags[1][2]]
        kinds=Set{Symbol}()
        for pcs in values(sg.pcurves), ent in (pcs.fwd,pcs.rev)
            ent===nothing && continue
            inner=hasproperty(ent,:neg2d) ? ent.neg2d[1] : ent
            push!(kinds,first(propertynames(inner)))
        end
        @test :ellips2d in kinds && :lin2d in kinds
    finally
        A.finalize()
    end
end

@testset "B-spline edge pcurve synthesis on a plane record" begin
    A.initialize()
    try
        occ=A.model.occ
        # `_brep_pc_on_surface`'s `_brep_synth_pcurve` branch: a planar
        # non-rational B-spline projects its poles into the plane frame,
        # giving an exact spl2d pcurve evaluated against the 3-D trace
        pa=[occ.add_point(x,y,0.0) for (x,y) in
            ((0.2,0.1),(0.6,0.8),(1.1,0.2),(1.4,0.6))]
        bs=occ.add_bspline(pa)
        occ.synchronize()
        m=A._with_model(identity)
        crv,t0,t1=M._brep_edge_curve(m,bs)
        sg=(occ=:plane,center=(0.0,0.0,0.0),axis=(0.0,0.0,1.0),
            X=(1.0,0.0,0.0),Y=(0.0,1.0,0.0),reversed=false,
            pcurves=Dict{Int,NamedTuple}())
        pc,neg=M._brep_pc_on_surface(sg,crv,t0,t1,1e-9,"test")
        @test !neg && hasproperty(pc,:spl2d)
        for s in range(t0,t1,length=9)
            p3=M._brep_curve_eval(crv,s)
            uv=M._brep_pc_eval(pc,s)
            @test uv[1]≈p3[1] atol=1e-12
            @test uv[2]≈p3[2] atol=1e-12
        end
    finally
        A.finalize()
    end
end

@testset "occ.intersect — coincident sheets: ellipse disk vs rectangle" begin
    A.initialize()
    try
        occ=A.model.occ
        # coincident z=0 sheets: the 2×2 rectangle's corners sit outside
        # the 1.4×0.8 ellipse, so the overlap is the slab-clipped lens —
        # 2 ellipse arcs joined by 2 chords, analytic area
        d=occ.add_disk(0.0,0.0,0.0,1.4,0.8)
        r=occ.add_rectangle(-1.0,-1.0,0.0,2.0,2.0)
        occ.synchronize()
        out,_=occ.intersect([(2,d)],[(2,r)])
        occ.synchronize()
        @test length(out)==1
        @test occ.get_mass(2,out[1][2])≈
            2*0.8*(sqrt(1-1/1.96)+1.4*asin(1/1.4)) rtol=1e-9
    finally
        A.finalize()
    end
    A.initialize()
    try
        occ=A.model.occ
        # and the reverse: a wide rectangle swallows the ellipse whole —
        # the intersection is the full disk with its ellipse boundary
        d=occ.add_disk(0.0,0.0,0.0,1.4,0.8)
        r=occ.add_rectangle(-3.0,-3.0,0.0,6.0,6.0)
        occ.synchronize()
        out,_=occ.intersect([(2,d)],[(2,r)])
        occ.synchronize()
        @test length(out)==1
        @test occ.get_mass(2,out[1][2])≈π*1.4*0.8 rtol=1e-9
    finally
        A.finalize()
    end
end

@testset "occ.fragment — tensor B-spline face fragmented by box" begin
    A.initialize()
    try
        occ=A.model.occ
        # imported trimmed B-spline surface × a box's six plane faces —
        # the tensor_bspline×plane pair sections the face on each plane's
        # zero level set and emits same-parameter section curves.
        v=occ.import_shapes(joinpath(_FIX,"iges_occ_trimmed_bspline.igs"))
        occ.synchronize()
        (x0,y0,z0,x1,y1,z1)=occ.get_bounding_box(v[1][1],v[1][2])
        pad=max(x1-x0,y1-y0)/4
        occ.add_box(x0+0.1*(x1-x0),y0-pad,z0-pad,
                    0.4*(x1-x0),y1-y0+2pad,z1-z0+2pad)
        occ.synchronize()
        out,_=occ.fragment([(v[1][1],v[1][2])],[(3,1)])
        occ.synchronize()
        # the pair no longer rejects — both operands survive as fragments
        @test (3,1) in out
        ents=occ.get_entities()
        @test count(e->e[1]==2,ents)>=7
    finally
        A.finalize()
    end
end

@testset "Boolean face sections — tensor patch × analytic partners" begin
    A.initialize()
    try
        occ=A.model.occ
        # mildly wavy 4×4 B-spline sheet over [0,2.1]² × z≤0.15, pierced
        # by each analytic primitive — the marching-squares contour and
        # the Hermite-fitted same-parameter pcurves must ride both
        # surfaces within the fit tolerance.
        tags=Int[]
        for j in 0:3, i in 0:3
            x,y=i*0.7,j*0.7
            push!(tags,occ.add_point(x,y,0.15*sin(π*x)*sin(π*y)))
        end
        occ.add_bspline_surface(tags,4)
        occ.add_sphere(1.05,1.05,0.0,0.9)
        occ.add_cylinder(0.4,0.4,-1.0,0.0,0.0,3.0,0.5)
        occ.add_cone(1.5,0.5,-0.5,0.0,0.0,2.0,0.7,0.2)
        occ.add_torus(1.05,1.05,0.0,0.8,0.25)
        occ.synchronize()
        m=A._with_model(identity)
        tsg=nothing
        sgs=Dict{Symbol,NamedTuple}()
        for sg in values(m.surface_geometry)
            k=get(sg,:occ,nothing)
            k in (:tensor_bspline,:tensor_bezier) && (tsg=sg)
            k in (:sphere,:cylinder,:cone,:torus) &&
                !haskey(sgs,k) && (sgs[k]=sg)
        end
        @test tsg!==nothing
        @test length(sgs)==4
        for kind in (:sphere,:cylinder,:cone,:torus)
            sgB=sgs[kind]
            fld,_=M._brep_partner_field(sgB)
            secs=M._brep_face_sections((sg=tsg,),(sg=sgB,),"test")
            @test !isempty(secs)
            for s in secs
                kn=s.curve.g.curve.curve.knots
                t0,t1=kn[1],kn[end]
                md=0.0
                for i in 0:20
                    tt=fma(t1-t0,i/20,t0)
                    p=M._brep_curve_eval(s.curve,tt)
                    md=max(md,abs(fld(p)))
                end
                # the section curve stays on the partner surface —
                # the fitted pcurves share the section parameter, so the
                # mid-span lifts name the same contour point
                @test md<5e-6
                tm=(t0+t1)/2
                PA=M._brep_surface_eval(tsg,M._brep_pc_eval(s.pcA,tm)...)
                PB=M._brep_surface_eval(sgB,M._brep_pc_eval(s.pcB,tm)...)
                @test M._b3dist(PA,PB)<5e-6
            end
        end
    finally
        A.finalize()
    end
end

@testset "Boolean face sections — analytic × analytic pairs" begin
    A.initialize()
    try
        occ=A.model.occ
        # mutually piercing analytic primitives covering the formerly
        # unsupported pairs — the generic level-set section traces the
        # contour on the better-conditioned chart, unwraps periodic
        # parameters, and emits shared-parameter fitted pcurves.
        occ.add_sphere(0.0,0.0,0.0,1.2)
        occ.add_cone(0.3,0.0,-0.7,0.0,0.0,1.6,0.5,0.95)
        occ.add_cylinder(0.0,0.0,-1.0,0.0,0.0,2.4,0.55)
        occ.add_torus(0.55,0.0,0.0,0.9,0.3)
        occ.synchronize()
        m=A._with_model(identity)
        sgs=NamedTuple[]
        for sg in values(m.surface_geometry)
            get(sg,:occ,nothing) in (:sphere,:cone,:cylinder,:torus) &&
                push!(sgs,sg)
        end
        @test length(sgs)>=4
        by=Dict(s.occ=>s for s in sgs)
        for (ka,kb) in ((:sphere,:cone),(:cone,:cylinder),
                        (:torus,:cylinder),(:sphere,:torus))
            sa,sb=by[ka],by[kb]
            fld,_=M._brep_partner_field(sb)
            secs=M._brep_face_sections((sg=sa,),(sg=sb,),"test")
            @test !isempty(secs)
            for s in secs
                kn=s.curve.g.curve.curve.knots
                t0,t1=kn[1],kn[end]
                md=0.0
                for i in 0:20
                    tt=fma(t1-t0,i/20,t0)
                    p=M._brep_curve_eval(s.curve,tt)
                    md=max(md,abs(fld(p)))
                end
                @test md<5e-6
                tm=(t0+t1)/2
                PA=M._brep_surface_eval(sa,M._brep_pc_eval(s.pcA,tm)...)
                PB=M._brep_surface_eval(sb,M._brep_pc_eval(s.pcB,tm)...)
                @test M._b3dist(PA,PB)<5e-6
            end
        end
    finally
        A.finalize()
    end
end

@testset "occ.boolean — sphere × cone all four operations materialize" begin
    # Upstream Gmsh 4.15.2 on this fixture: cut/fuse/intersect each emit
    # one volume, fragment emits three. The generic level-set section
    # crosses the sphere's periodic-u seam; the fitted pcurve endpoints
    # must lift into each piece's u band (not the canonical band) and the
    # closed-carrier merge must not join across the evaluator's domain
    # wrap — both regressions produced dangling/one-edge wires here.
    A.initialize()
    try
        occ=A.model.occ
        occ.add_sphere(0.0,0.0,0.0,1.2)
        occ.add_cone(0.3,0.0,-0.7,0.0,0.0,1.6,0.5,0.95)
        occ.synchronize()
        for (op,want) in ((occ.cut,1),(occ.fuse,1),
                          (occ.intersect,1),(occ.fragment,3))
            out,_=op([(3,1)],[(3,2)];removeObject=false,removeTool=false)
            @test length(out)==want
            for (_,t) in out
                # materialized volume evaluates a finite bounding box —
                # a dropped/merged section would leave a degenerate loop
                bb=collect(occ.get_bounding_box(3,t))
                @test all(isfinite,bb)
                @test bb[4]-bb[1]>0.1 && bb[5]-bb[2]>0.1 &&
                      bb[6]-bb[3]>0.1
            end
        end
    finally
        A.finalize()
    end
end
