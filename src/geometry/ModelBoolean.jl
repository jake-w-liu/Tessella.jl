# ── Analytic B-rep Booleans over materialized OCC-style solids ──────────────
#
# `boolean_volumes!` materializes the result boundary like Gmsh's OCC kernel
# (`BRepAlgoAPI_*` over `BOPAlgo`): operand faces are intersected pairwise to
# section curves, operand edges split where they pierce the other solid, each
# face's parameter domain is re-wired by directed half-edge traversal, kept
# regions are classified against the other solid, and the surviving pieces
# are materialized as ordinary points/curves/surfaces/loops/volume entities
# with OCC geometry records. The operand triangle snapshots remain for the
# volume meshing path (`mesh_boolean`), which is unchanged.
#
# Supported sections are restricted to what the entity model can represent —
# line and circle curves on plane/cylinder/sphere/cone faces. Anything else
# (ellipses, general surface-surface curves, torus sections, coincident
# same-domain surfaces) fails explicitly rather than approximating.
#
# Internal representation:
#   solid = (pts, ptags, edges, ctags, faces)
#     pts/ptags — vertex coordinates and source point tags (0 when new)
#     edges/ctags — (curve, t0, t1, v1, v2) and source curve tags (0 when new);
#         curve = (kind=:line,o,d) | (kind=:circle,c,n,x,y,r) |
#                 (kind=:degenerate,xyz)
#     faces — (sg, stag, orev, wires)
#         sg    — OCC surface record (:plane/:cylinder/:sphere/:cone/:torus)
#         orev  — ±1 outward orientation flag (the surface's sign in the shell)
#         wires — Vector of wires; wires[1] is the outer wire. A wire is a
#                 Vector of wedges (eidx, sign, pc); pc is the pcurve record
#                 for that direction ((lin2d=(o,d))/(circ2d=(c,x,r)) evaluated
#                 at the edge parameter, matching `m.surface_geometry` pcurves)

const _BREP_TOL = 1e-9

# ── small 3-D helpers ────────────────────────────────────────────────────────

@inline _b3(a,b) = (a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _b3add(a,b) = (a[1]+b[1],a[2]+b[2],a[3]+b[3])
@inline _b3mul(a,s) = (a[1]*s,a[2]*s,a[3]*s)
@inline _b3dot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _b3cross(a,b) = (a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
@inline _b3norm(a) = sqrt(_b3dot(a,a))
@inline _b3unit(a) = (l=_b3norm(a); (a[1]/l,a[2]/l,a[3]/l))
@inline _b3dist(a,b) = _b3norm(_b3(a,b))

function _brep_curve_eval(curve::NamedTuple,t::Float64)
    if curve.kind===:line
        return _b3add(curve.o,_b3mul(curve.d,t))
    elseif curve.kind===:circle
        s,c=_gm_sincos(t)
        return _b3add(curve.c,_b3add(_b3mul(curve.x,c*curve.r),
                                   _b3mul(curve.y,s*curve.r)))
    elseif curve.kind===:ellipse
        return _occ_ellipse_point(curve.g,t)
    elseif curve.kind===:bspline
        # exact evaluator of the stored/copied OCC curve — homogeneous
        # de Boor through `_model_occ_curve_d2`
        return _model_occ_curve_d2(curve.curve,t)[1]
    end
    throw(ArgumentError(
        "_brep_curve_eval: $(curve.kind) edge has no 3-D curve evaluator"))
end

# ── surface evaluation / inversion ──────────────────────────────────────────

_brep_surface_eval(sg::NamedTuple,u::Float64,v::Float64) =
    _occ_surface_point(sg,u,v)

# Geometric unit normal at (u,v), oriented by the record's `reversed` flag
# like `_occ_surface_normal` — mesh winding aligns to the face normal.
function _brep_surface_normal(sg::NamedTuple,u::Float64,v::Float64)
    n=if sg.occ===:plane
        sg.axis
    elseif sg.occ===:cylinder
        s,c=_gm_sincos(u)
        _b3unit((fma(c,sg.X[1],s*sg.Y[1]),fma(c,sg.X[2],s*sg.Y[2]),
                 fma(c,sg.X[3],s*sg.Y[3])))
    elseif sg.occ===:sphere
        su,cu=_gm_sincos(u); sv,cv=_gm_sincos(v)
        _b3unit((fma(cv*cu,sg.X[1],fma(cv*su,sg.Y[1],sv*sg.axis[1])),
                 fma(cv*cu,sg.X[2],fma(cv*su,sg.Y[2],sv*sg.axis[2])),
                 fma(cv*cu,sg.X[3],fma(cv*su,sg.Y[3],sv*sg.axis[3]))))
    elseif sg.occ===:cone
        s,c=_gm_sincos(u)
        slope=(sg.r2-sg.r1)/sg.height
        l=sqrt(1.0+slope*slope)
        radial=(fma(c,sg.X[1],s*sg.Y[1]),fma(c,sg.X[2],s*sg.Y[2]),
                fma(c,sg.X[3],s*sg.Y[3]))
        ((radial[1]-slope*sg.axis[1])/l,
         (radial[2]-slope*sg.axis[2])/l,
         (radial[3]-slope*sg.axis[3])/l)
    elseif sg.occ===:torus
        # outward tube normal: cos v·(radial dir) + sin v·axis
        s,c=_gm_sincos(u); sv,cv=_gm_sincos(v)
        _b3unit((fma(cv*c,sg.X[1],fma(cv*s,sg.Y[1],sv*sg.axis[1])),
                 fma(cv*c,sg.X[2],fma(cv*s,sg.Y[2],sv*sg.axis[2])),
                 fma(cv*c,sg.X[3],fma(cv*s,sg.Y[3],sv*sg.axis[3]))))
    elseif sg.occ in (:tensor_bspline,:tensor_bezier)
        _,du,dv=_occ_surface_d1(sg,u,v)
        _b3unit(_b3cross(du,dv))
    else
        throw(ArgumentError("_brep_surface_normal: unsupported surface $(sg.occ)"))
    end
    return get(sg,:reversed,false) ? (-n[1],-n[2],-n[3]) : n
end

# Analytic inversion — (u,v) of a point known to lie on the surface.
function _brep_surface_uv(sg::NamedTuple,p::NTuple{3,Float64})
    if sg.occ===:plane
        d=_b3(p,sg.center)
        return (_b3dot(d,sg.X),_b3dot(d,sg.Y))
    elseif sg.occ===:cylinder || sg.occ===:cone
        d=_b3(p,sg.center)
        u=atan(_b3dot(d,sg.Y),_b3dot(d,sg.X))
        u<0 && (u+=_OCC_TWO_PI)
        v=_b3dot(d,sg.axis)
        if sg.occ===:cone
            # the record's v is slant length along the generatrix
            # (wire pcurves put the top rim at v=slant); axial z is
            # v·cos(half-angle)
            _,csa=_gm_sincos(_gm_atan((sg.r2-sg.r1)/sg.height))
            v/=csa
        end
        return (u,v)
    elseif sg.occ===:sphere
        d=_b3(p,sg.center)
        z=_b3dot(d,sg.axis); rr=_b3dot(d,sg.X); ri=_b3dot(d,sg.Y)
        v=atan(z,sqrt(rr*rr+ri*ri))
        u=atan(ri,rr); u<0 && (u+=_OCC_TWO_PI)
        return (u,v)
    elseif sg.occ===:torus
        # P = C + (r1+r2·cos v)·(cos u·X + sin u·Y) + r2·sin v·axis —
        # u is the spine angle, v the tube angle measured from the
        # outward radial toward +axis.
        d=_b3(p,sg.center)
        rr=_b3dot(d,sg.X); ri=_b3dot(d,sg.Y)
        u=atan(ri,rr); u<0 && (u+=_OCC_TWO_PI)
        v=atan(_b3dot(d,sg.axis),sqrt(rr*rr+ri*ri)-sg.r1)
        v<0 && (v+=_OCC_TWO_PI)
        return (u,v)
    elseif sg.occ in (:tensor_bspline,:tensor_bezier)
        u,v,_=_model_nurbs_surface_project(sg,p)
        return (u,v)
    end
    throw(ArgumentError("_brep_surface_uv: unsupported surface $(sg.occ)"))
end

# u-period of a parametrization: analytic revolutions wrap at 2π; a
# tensor patch's u-period is its knot-domain width when its seam is
# geometrically closed (STEP/IGES B-spline bounded faces on periodic
# supports) — projection pcurves unwrap across the seam with it.
_brep_uperiod(sg::NamedTuple) =
    if sg.occ in (:cylinder,:sphere,:cone,:torus)
        _OCC_TWO_PI
    elseif sg.occ in (:tensor_bspline,:tensor_bezier)
        _model_nurbs_patch_seam_period(sg,1)
    else
        0.0
    end
# v-chart period of a parametrization: a torus's tube angle and a
# tensor patch's v-seam period; other charts are non-periodic in v.
_brep_vperiod(sg::NamedTuple) =
    if sg.occ===:torus
        _OCC_TWO_PI
    elseif sg.occ in (:tensor_bspline,:tensor_bezier)
        _model_nurbs_patch_seam_period(sg,2)
    else
        0.0
    end

# ── pcurve evaluation ────────────────────────────────────────────────────────

# `neg2d` pcurve wrapper: evaluated at the edge parameter t as `inner(−t)` —
# the form a left-handed uv trace takes, e.g. a circle rim whose axis is
# anti-parallel to the face's synthesized plane frame (the trace is the same
# circle run backwards, which no right-handed `circ2d` can write down).
_brep_pc_unneg(pc) =
    hasproperty(pc,:neg2d) ? (pc.neg2d[1],true) : (pc,false)

# 2π-periodic pcurves — intersection and inversion params report mod 2π, so
# spanning and lift logic treats circles and ellipses identically (`neg2d`
# keeps the same periodic inner record).
_brep_pc_iscirc(pc) = hasproperty(pc,:circ2d) ||
    hasproperty(pc,:ellips2d) ||
    (hasproperty(pc,:neg2d) && _brep_pc_iscirc(pc.neg2d[1]))

# Every stored pcurve kind the OCC face model produces — `lin2d`, `circ2d`,
# `ellips2d`, `spl2d`, the loft `tan2d`/`tan2dv` conic-reparam traces and the
# `neg2d` mirror — evaluates through the shared `Geom2d` adaptors.
@inline _brep_pc_eval(pc,t) = _occ_pcurve_d0(pc,Float64(t))

# Tangent vector of a pcurve at t (not normalized).
@inline _brep_pc_tangent(pc,t) = _occ_pcurve_d1(pc,Float64(t))[2]

# ── generic pcurve machinery ────────────────────────────────────────────────
#
# The analytic lin2d/circ2d fast paths below each carry a generic sibling
# covering `ellips2d`, `spl2d`, `tan2d`/`tan2dv` and any future stored record:
# structure-aware sample grids + bisected sign-changes for level roots, and a
# grid-seeded Newton projection for parameter inversion. These are the same
# primitives the meshing evaluators use — the records already carry exact
# `_occ_pcurve_d0`/`d1`/`d2` semantics.

# parameter-space sample grid over [lo,hi]: B-splines densify per knot span,
# smooth records take a uniform grid
function _brep_pc_grid(pc,lo,hi;subs::Int=8)
    local bk
    if hasproperty(pc,:spl2d)
        bk=unique(pc.spl2d[1])
    elseif hasproperty(pc,:tan2d)
        bk=pc.tan2d[1].knots
    elseif hasproperty(pc,:tan2dv)
        bk=pc.tan2dv[1].knots
    else
        bk=nothing
    end
    ts=Float64[]
    if bk!==nothing
        bounds=vcat(Float64[lo],
                    Float64[k for k in bk if lo-1e-12<k<hi+1e-12],
                    Float64[hi])
        sort!(bounds)
        for i in 1:length(bounds)-1
            a,b=bounds[i],bounds[i+1]
            for j in 0:subs
                push!(ts,fma((b-a),j/subs,a))
            end
        end
    else
        n=64
        for i in 0:n
            push!(ts,fma(hi-lo,i/n,lo))
        end
    end
    return ts
end

# all t in the sorted sample grid `ts` where the smooth scalar field `f`
# crosses zero — sign changes bisect to parameter accuracy; interior grid
# points whose |f| sits below `tol`·scale report as tangent contacts (the
# analytic pair solvers' disc≈0 convention)
function _brep_scalar_roots(f,ts::AbstractVector{<:Real},tol;
                            tangents::Bool=true)
    n=length(ts)
    fs=[f(Float64(t)) for t in ts]
    fscl=max(1.0,maximum(abs,fs))
    ztol=tol*fscl
    out=Float64[]
    for i in 1:n-1
        a,b=ts[i],ts[i+1];fa,fb=fs[i],fs[i+1]
        if fa==0.0
            (isempty(out) || a-out[end]>1e-9*max(1.0,b-a)) &&
                push!(out,Float64(a))
        elseif signbit(fa)!=signbit(fb)
            flo=fa;lo_,hi_=a,b
            for _ in 1:80
                mid=(lo_+hi_)/2;fm=f(mid)
                fm==0.0 && (lo_=hi_=mid;break)
                signbit(fm)==signbit(flo) ? (lo_=mid;flo=fm) : (hi_=mid)
            end
            r=(lo_+hi_)/2
            (isempty(out) || r-out[end]>1e-9*max(1.0,b-a)) &&
                push!(out,r)
        elseif tangents && i>1 && i<n-1 && abs(fa)<=ztol &&
                abs(fa)<=abs(fs[i-1]) && abs(fa)<=abs(fb)
            # tangent touch at a grid node — mirror the analytic
            # disc==0 → two-coincident-roots convention as one entry
            (isempty(out) || ts[i]-out[end]>1e-9*max(1.0,b-a)) &&
                push!(out,Float64(ts[i]))
        end
    end
    fs[n]==0.0 && (isempty(out) || ts[n]-out[end]>1e-9*max(1.0,ts[n]-ts[n-1])) &&
        push!(out,Float64(ts[n]))
    return out
end

# implicit level of an elementary pcurve — ~signed uv distance (negative
# inside the bounded side of conics). `nothing` for spline/reparam traces.
# Generated spl2d section traces carry `lvl2d` — the level set the trace was
# extracted from — which is stronger than distance-field minima where the
# trace merely brushes a wire without crossing it.
function _brep_pc_level(pc)
    hasproperty(pc,:neg2d) &&
        return _brep_pc_level(pc.neg2d[1])
    if hasproperty(pc,:lvl2d)
        return pc.lvl2d
    elseif hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        nn=hypot(d[1],d[2])
        nn==0.0 && return nothing
        return uv->((uv[1]-o[1])*d[2]-(uv[2]-o[2])*d[1])/nn
    elseif hasproperty(pc,:circ2d)
        c,x,r=pc.circ2d
        return uv->(hypot(uv[1]-c[1],uv[2]-c[2])-r)
    elseif hasproperty(pc,:ellips2d)
        o,x,a,b=pc.ellips2d
        return uv->(hypot(((uv[1]-o[1])*x[1]+(uv[2]-o[2])*x[2])/a,
                          ((uv[1]-o[1])*x[2]-(uv[2]-o[2])*x[1])/b)-1.0)*
                    min(a,b)
    end
    return nothing
end

# natural parameter domain of a pcurve record — splines and conic-reparam
# traces are bounded by their knot tables; elementary curves are unbounded
function _brep_pc_domain(pc)
    if hasproperty(pc,:spl2d)
        k=pc.spl2d[1]
        return (k[1],k[end])
    elseif hasproperty(pc,:tan2d)
        k=pc.tan2d[1].knots
        return (k[1],k[end])
    elseif hasproperty(pc,:tan2dv)
        k=pc.tan2dv[1].knots
        return (k[1],k[end])
    end
    return (-Inf,Inf)
end

# Newton refinement of a point-on-pcurve parameter over [lo,hi] — the
# g(t)=(p(t)−uv)·p′(t) normal-equation root, seeded on the curve's grid.
# `nothing` when no in-range root matches within `tol`.
function _brep_pc_project(pc,uv,lo,hi,tol)
    r=_brep_pc_project_res(pc,uv,lo,hi)
    r===nothing && return nothing
    t,res=r
    scl=max(1.0,abs(uv[1]),abs(uv[2]))
    res<10tol*scl || return nothing
    return t
end

# project uv onto pc over [lo,hi] — returns (t, uv residual) at the
# converged nearest point regardless of tolerance
function _brep_pc_project_res(pc,uv,lo,hi)
    ts=_brep_pc_grid(pc,lo,hi;subs=4)
    best=ts[1];bd=Inf
    for t in ts
        p=_brep_pc_eval(pc,t)
        d2=(p[1]-uv[1])*(p[1]-uv[1])+(p[2]-uv[2])*(p[2]-uv[2])
        d2<bd && (bd=d2;best=t)
    end
    t=best
    for _ in 1:40
        p,d1,d2=_occ_pcurve_d2(pc,t)
        r=(p[1]-uv[1],p[2]-uv[2])
        g=r[1]*d1[1]+r[2]*d1[2]
        gp=d1[1]*d1[1]+d1[2]*d1[2]+r[1]*d2[1]+r[2]*d2[2]
        gp==0.0 && break
        step=g/gp
        t=clamp(t-step,lo,hi)
        abs(step)<=2eps(Float64)*max(1.0,abs(t)) && break
    end
    p=_brep_pc_eval(pc,t)
    return t,hypot(p[1]-uv[1],p[2]-uv[2])
end

# ── operand extraction ──────────────────────────────────────────────────────

@inline _brep_plane_uv(sg,p) =
    (d=_b3(p,sg.center); (_b3dot(d,sg.X),_b3dot(d,sg.Y)))

function _brep_edge_curve(m::GeoModel,ctag::Int)
    rec=get(m.curve_geometry,ctag,nothing)
    a,b=m.curves[ctag]
    if rec===nothing || !hasproperty(rec,:occ)
        pa=m.points[a];pb=m.points[b]
        d=_b3(pb,pa);len=_b3norm(d)
        return (kind=:line,o=pa,d=(d[1]/len,d[2]/len,d[3]/len)),0.0,len
    elseif rec.occ===:line ||
            (rec.occ===:copied_curve && get(rec,:converted_from,nothing)===:line)
        # NURBS-converted records keep the analytic parametrization as inert
        # metadata — the materializer consumes the preserved shadow, whose
        # locus is bitwise identical.
        return (kind=:line,o=rec.origin,d=rec.dir),rec.t0,rec.t1
    elseif rec.occ===:circle ||
            (rec.occ===:copied_curve && get(rec,:converted_from,nothing)===:circle)
        return (kind=:circle,c=rec.center,n=rec.n,x=rec.X,y=rec.Y,r=rec.r),
               rec.t0,rec.t1
    elseif rec.occ===:ellipse ||
            (rec.occ===:copied_curve && get(rec,:converted_from,nothing)===:ellipse)
        # like the converted line/circle shadows above: a NURBS-converted
        # ellipse keeps its analytic record as inert metadata — the same
        # evaluator `_occ_ellipse_point` serves both forms bitwise.
        return (kind=:ellipse,c=rec.center,n=rec.n,x=rec.X,y=rec.Y,
                rx=rec.rx,ry=rec.ry,g=rec),rec.t0,rec.t1
    elseif rec.occ===:copied_curve && hasproperty(rec,:curve)
        # imported/swept B-spline or NURBS edge — the record owns the exact
        # `_ModelOCCNURBSCurve`/`_ModelOCCReparamCurve` evaluator; the engine
        # consumes it through `_brep_curve_eval`/`_brep_curve_param`.
        return (kind=:bspline,curve=rec.curve,g=rec,
                ctype=get(m.curve_types,ctag,:copied_curve)),rec.t0,rec.t1
    elseif rec.occ===:surface_iso && hasproperty(rec,:face) &&
            _model_nurbs_surface(rec.face)
        # isoparametric edge of a tensor patch — the composed
        # surface∘lin2d evaluator rides the patch exactly; the seed
        # polyline supplies knots/controls for grids and projections
        lo,hi=minmax(rec.t0,rec.t1)
        n=max(9,2*length(unique(rec.axis==1 ?
            rec.face.patch.surface.knots_u :
            rec.face.patch.surface.knots_v))+1)
        pts=[_model_surface_iso_point(rec,fma(hi-lo,i/(n-1),lo))
             for i in 0:n-1]
        pknots=Float64[lo,lo]
        append!(pknots,[fma(hi-lo,i/(n-1),lo) for i in 1:n-2])
        push!(pknots,hi,hi)
        poly=NURBS.NURBSCurve(1,pknots,pts)
        pc=rec.axis==1 ? (lin2d=((0.0,rec.fixed),(1.0,0.0)),) :
                         (lin2d=((rec.fixed,0.0),(0.0,1.0)),)
        ev=_ModelOCCSectionCurve(:bspline,rec.face,pc,poly)
        return (kind=:bspline,curve=ev,g=rec,
                ctype=get(m.curve_types,ctag,:surface_iso)),rec.t0,rec.t1
    elseif rec.occ===:degenerate
        return (kind=:degenerate,xyz=m.points[a]),rec.t0,rec.t1
    end
    throw(ArgumentError("_brep_edge_curve: unsupported curve record"))
end

# 3-D sample points of a stored curve loop's traversal: edge start points
# plus arc midpoints (an arc's winding contribution needs a sample inside
# the span).
function _brep_wire_pts(m::GeoModel,loop_tag::Int)
    pts=NTuple{3,Float64}[]
    for ce in m.loops[loop_tag]
        a,b=m.curves[abs(ce)]
        push!(pts,m.points[ce>0 ? a : b])
        curve,t0,t1=_brep_edge_curve(m,abs(ce))
        # an arc's winding contribution needs a sample inside the span —
        # every non-straight evaluable curve contributes its midpoint
        curve.kind in (:line,:degenerate) && continue
        if a==b
            # a closed edge's two endpoint samples coincide — the winding
            # comes from the span samples alone, so take several
            for i in 1:3
                push!(pts,_brep_curve_eval(curve,t0+(t1-t0)*i/4))
            end
        else
            push!(pts,_brep_curve_eval(curve,(t0+t1)/2))
        end
    end
    return pts
end

# Newell normal of a stored curve loop's traversal. The surface-loop sign
# convention ties a face's reversal flag to its outer wire's winding
# normal, so synthesized face records take their axis from the winding —
# `_model_plane_frame`'s certified normal is a canonical fit whose sign is
# unrelated to traversal.
function _brep_wire_newell(m::GeoModel,loop_tag::Int)
    pts=_brep_wire_pts(m,loop_tag)
    n=(0.0,0.0,0.0)
    for i in eachindex(pts)
        p=pts[i];q=pts[mod1(i+1,length(pts))]
        n=(n[1]+(p[2]-q[2])*(p[3]+q[3]),
           n[2]+(p[3]-q[3])*(p[1]+q[1]),
           n[3]+(p[1]-q[1])*(p[2]+q[2]))
    end
    L=sqrt(_b3dot(n,n))
    return L==0.0 ? (0.0,0.0,1.0) : (n[1]/L,n[2]/L,n[3]/L)
end

# Synthesize a :plane record for a native plane surface so every operand face
# carries an occ-style record. The record's `axis` must match the outer
# wire's winding normal — the surface-loop sign then reads the use reversal
# uniformly across operands (`orev·axis` = outward); the frame's exact
# fitted normal supplies the magnitude, the winding supplies the sign.
function _brep_plane_record(m::GeoModel,stag::Int,caller::AbstractString)
    frame=_model_plane_frame(m,stag,caller)
    n=_brep_wire_newell(m,m.surfaces[stag][1])
    if _b3dot(n,frame.normal)<0
        return (occ=:plane,center=frame.origin,axis=.-frame.normal,
                X=frame.first_direction,Y=.-frame.second_direction,
                reversed=false,pcurves=Dict{Int,NamedTuple}())
    end
    return (occ=:plane,center=frame.origin,axis=frame.normal,
            X=frame.first_direction,Y=frame.second_direction,
            reversed=false,pcurves=Dict{Int,NamedTuple}())
end

# The pcurve of an edge on a synthesized plane record. A line projects its
# endpoints over its own parameter range; a circle projects its center and
# X direction into the face frame. An anti-aligned circle (e.g. a hole rim
# wound against the face normal) traces a left-handed uv circle, which a
# right-handed `circ2d` cannot represent — the `neg2d` wrapper evaluates
# the same record at −t, exactly the mirrored trace.
function _brep_synth_pcurve(sg,curve::NamedTuple,t0,t1,caller::AbstractString)
    if curve.kind===:circle
        xu=(_b3dot(curve.x,sg.X),_b3dot(curve.x,sg.Y))
        pc=(circ2d=(_brep_plane_uv(sg,curve.c),xu,curve.r),)
        return _b3dot(curve.n,_b3cross(sg.X,sg.Y))>=0 ? pc : (neg2d=(pc,),)
    elseif curve.kind===:ellipse
        # exact frame projection — same handedness rule as the circle
        xu=(_b3dot(curve.x,sg.X),_b3dot(curve.x,sg.Y))
        pc=(ellips2d=(_brep_plane_uv(sg,curve.c),xu,curve.rx,curve.ry),)
        return _b3dot(curve.n,_b3cross(sg.X,sg.Y))>=0 ? pc : (neg2d=(pc,),)
    elseif curve.kind===:bspline
        cv=curve.curve
        if cv isa _ModelOCCNURBSCurve &&
                all(w->w==1.0,cv.curve.weights)
            # exact: the plane's orthonormal projection acts on the poles
            # affinely — a non-rational B-spline's uv trace is the B-spline
            # of its projected control points over the same knots
            cps=NTuple{2,Float64}[
                _brep_plane_uv(sg,c) for c in cv.curve.controls]
            return (spl2d=(cv.curve.knots,cv.curve.degree,cps),)
        end
        # rational / reparametrized trace — the importer's own convention:
        # a monotone C1 cubic through edge-parameter samples
        ts=[t0+(t1-t0)*i/16 for i in 0:16]
        uvs=[_brep_plane_uv(sg,_brep_curve_eval(curve,t)) for t in ts]
        return (spl2d=_pchip_bspline_2d(ts,uvs),)
    end
    curve.kind in (:line,:degenerate) || throw(ArgumentError(
        "$caller: unsupported $(curve.kind) edge on an untyped plane"))
    o=_brep_plane_uv(sg,_brep_curve_eval(curve,t0))
    e=_brep_plane_uv(sg,_brep_curve_eval(curve,t1))
    dt=t1-t0
    return (lin2d=(o,((e[1]-o[1])/dt,(e[2]-o[2])/dt)),)
end

function _brep_extract(m::GeoModel,tag::Int,caller::AbstractString)
    pts=NTuple{3,Float64}[];ptags=Int[]
    edges=NamedTuple[];ctags=Int[]
    faces=NamedTuple[]
    pmap=Dict{Int,Int}();emap=Dict{Int,Int}()
    function vertex(p)
        get!(pmap,p) do
            push!(pts,m.points[p]);push!(ptags,p);length(pts)
        end
    end
    for shell in m.volumes[tag]
        for ssigned in m.surface_loops[shell]
            stag=abs(ssigned);orev=ssigned>0 ? 1 : -1
            sg=get(m.surface_geometry,stag,nothing)
            occ_face=sg!==nothing && hasproperty(sg,:occ)
            occ_face || (sg=_brep_plane_record(m,stag,caller))
            wires=Vector{Tuple{Int,Int,NamedTuple}}[]
            for loop_tag in m.surfaces[stag]
                wire=Tuple{Int,Int,NamedTuple}[]
                for sc in m.loops[loop_tag]
                    ctag=abs(sc);sgn=sc>0 ? 1 : -1
                    eidx=get!(emap,ctag) do
                        curve,t0,t1=_brep_edge_curve(m,ctag)
                        a,b=m.curves[ctag]
                        push!(edges,(curve=curve,t0=t0,t1=t1,
                                     v1=vertex(a),v2=vertex(b)))
                        push!(ctags,ctag);length(edges)
                    end
                    local pc
                    if occ_face
                        ent=_brep_occ_face_pc(m,sg,ctag,stag,caller)
                        pc=(sgn<0 && ent.rev!==nothing) ? ent.rev : ent.fwd
                        pc===nothing && throw(ArgumentError(
                            "$caller: Surface[$stag] lacks a pcurve for " *
                            "Curve[$ctag]"))
                        erec=get(m.curve_geometry,ctag,nothing)
                        if erec===nothing || !hasproperty(erec,:occ)
                            # A record-free edge (an untyped materialized
                            # line) stores its pcurve on the `[0,1]`
                            # segment parameter; the engine runs it on the
                            # arc-length `[t0,t1]` convention — restretch
                            # the trace to this edge's range.
                            e=edges[eidx]
                            Δ=e.t1-e.t0
                            pc=_occ_pcurve_reparam(pc,-e.t0/Δ,1.0/Δ)
                        end
                    else
                        e=edges[eidx]
                        pc=_brep_synth_pcurve(sg,e.curve,e.t0,e.t1,caller)
                    end
                    push!(wire,(eidx,sgn,pc))
                end
                push!(wires,wire)
            end
            isempty(wires) && throw(ArgumentError(
                "$caller: Surface[$stag] has no boundary wires"))
            push!(faces,(sg=sg,stag=stag,orev=orev,wires=wires))
        end
    end
    isempty(faces) && throw(ArgumentError(
        "$caller: Volume[$tag] has no materialized boundary"))
    return (pts=pts,ptags=ptags,edges=edges,ctags=ctags,faces=faces)
end

# pcurve entry `(fwd,rev)` for an edge on an OCC-geometry face — the stored
# pcurve table, or the domain-axis lin2d for a `surface_iso` boundary edge
# of an untrimmed tensor patch (which carries no pcurve table).
function _brep_occ_face_pc(m,sg,ctag,stag,caller)
    pcdict=hasproperty(sg,:pcurves) ? sg.pcurves : nothing
    ent=pcdict===nothing ? nothing : get(pcdict,ctag,nothing)
    ent===nothing || return ent
    erec0=get(m.curve_geometry,ctag,nothing)
    (erec0!==nothing && hasproperty(erec0,:occ) &&
     erec0.occ===:surface_iso && hasproperty(erec0,:face) &&
     hasproperty(erec0.face,:patch) && hasproperty(sg,:patch) &&
     erec0.face.patch===sg.patch) || throw(ArgumentError(
        "$caller: Surface[$stag] lacks a pcurve for Curve[$ctag]"))
    return (fwd=erec0.axis==1 ?
            (lin2d=((0.0,erec0.fixed),(1.0,0.0)),) :
            (lin2d=((erec0.fixed,0.0),(0.0,1.0)),),
            rev=nothing)
end

# ── 2-D wire containment (even-odd parity) ──────────────────────────────────

# Rightward crossings of a horizontal ray from (pu,pv) against one pcurve over
# the edge's parameter range. Half-open convention: a vertex-level crossing
# counts only when the edge enters the half-plane from below.
function _brep_pc_ray_crossings(pc,t0,t1,pu,pv,tol)
    if hasproperty(pc,:neg2d)
        return _brep_pc_ray_crossings(pc.neg2d[1],-t1,-t0,pu,pv,tol)
    end
    lo,hi=minmax(t0,t1)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        abs(d[2])<tol && return 0
        t=(pv-o[2])/d[2]
        (t<lo-tol || t>hi+tol) && return 0
        x=o[1]+t*d[1]
        x<=pu+tol && return 0
        (t<=lo+tol && d[2]<0) && return 0
        (t>=hi-tol && d[2]>0) && return 0
        return 1
    elseif hasproperty(pc,:circ2d)
        c,x,r=pc.circ2d
        dy=pv-c[2]
        dy*dy>r*r && return 0
        # u(t)=c1+r·cos(t+α), v(t)=c2+r·sin(t+α), α=atan2(x2,x1)
        α=atan(x[2],x[1])
        β=asin(clamp(dy/r,-1.0,1.0))
        hits=0
        for γ in (β,π-β)
            t=mod(γ-α,_OCC_TWO_PI)
            px=c[1]+r*cos(γ)
            px<=pu+tol && continue
            # v'(t)=r·cos(t+α)=px−c1 — the endpoint rule mirrors lin2d: a
            # crossing at an arc's v-lower endpoint counts, at its v-upper
            # endpoint does not (vertex crossings count exactly once).
            dv=r*cos(γ)
            ok=false
            for k in -2:2
                tk=t+k*_OCC_TWO_PI
                (tk<lo-tol || tk>hi+tol) && continue
                if (tk>lo+tol && tk<hi-tol) ||
                   (tk<=lo+tol && dv>0) || (tk>=hi-tol && dv<0)
                    ok=true;break
                end
            end
            ok && (hits+=1)
        end
        return hits
    elseif hasproperty(pc,:ellips2d)
        o,x,a,b=pc.ellips2d
        # v(t)=o2+a·x2·cos+b·x1·sin — A·cos+B·sin=pv−o2
        A=a*x[2];B=b*x[1];C=pv-o[2]
        amp=hypot(A,B)
        amp<tol && return 0
        abs(C)>amp && return 0
        φ=atan(B,A);δ=acos(clamp(C/amp,-1.0,1.0))
        hits=0
        for γ in (φ-δ,φ+δ)
            px=o[1]+a*x[1]*cos(γ)-b*x[2]*sin(γ)
            px<=pu+tol && continue
            dv=-a*x[2]*sin(γ)+b*x[1]*cos(γ)
            ok=false
            for k in -2:2
                tk=γ+k*_OCC_TWO_PI
                (tk<lo-tol || tk>hi+tol) && continue
                if (tk>lo+tol && tk<hi-tol) ||
                   (tk<=lo+tol && dv>0) || (tk>=hi-tol && dv<0)
                    ok=true;break
                end
            end
            ok && (hits+=1)
        end
        return hits
    end
    # generic records (spl2d, tan2d/tan2dv): sign changes of v(t)−pv over
    # the curve's own sample grid; tangency on the ray adds an even count,
    # so parity follows the interior sign-change roots alone
    ts=_brep_pc_grid(pc,lo,hi)
    hits=0
    for t in _brep_scalar_roots(t->_brep_pc_eval(pc,t)[2]-pv,ts,tol;
                                tangents=false)
        px=_brep_pc_eval(pc,t)[1]
        px<=pu+tol && continue
        dv=_brep_pc_tangent(pc,t)[2]
        if (t>lo+tol && t<hi-tol) ||
           (t<=lo+tol && dv>0) || (t>=hi-tol && dv<0)
            hits+=1
        end
    end
    return hits
end

function _brep_uv_in_wire(edges,wire,pu,pv,tol)
    crossings=0
    for (eidx,sign,pc) in wire
        e=edges[eidx]
        crossings+=_brep_pc_ray_crossings(pc,e.t0,e.t1,pu,pv,tol)
    end
    return isodd(crossings)
end

# uv containment against a face (inside outer wire, outside hole wires).
function _brep_uv_in_face(edges,face,uv,tol)
    _brep_uv_in_wire(edges,face.wires[1],uv[1],uv[2],tol) || return false
    for k in 2:length(face.wires)
        _brep_uv_in_wire(edges,face.wires[k],uv[1],uv[2],tol) && return false
    end
    return true
end

# ── surface × surface sections ──────────────────────────────────────────────

_brep_pc_line_on_plane(sg,o,d) =
    (lin2d=(_brep_plane_uv(sg,o),(_b3dot(d,sg.X),_b3dot(d,sg.Y))),)

# A section circle's uv trace is a right-handed circ2d only when the
# circle's normal agrees with the uv frame's handedness; an anti-aligned
# trace mirrors the parameter — the neg2d wrapper evals the same record at
# −t so pcurve params keep tracking the edge parameter. Synthesized plane
# records may carry a left-handed frame (axis mirrored out of X×Y), so the
# test runs against X×Y rather than axis.
function _brep_pc_circ_on_plane(sg,c,n,x,r)
    pc=(circ2d=(_brep_plane_uv(sg,c),(_b3dot(x,sg.X),_b3dot(x,sg.Y)),r),)
    return _b3dot(n,_b3cross(sg.X,sg.Y))>=0 ? pc : (neg2d=(pc,),)
end

# circle at constant height on a cylinder/cone-like uv domain (u=angle,v=z).
function _brep_pc_circ_on_cylinder(sg,center,n,x,r)
    sgn=_b3dot(n,sg.axis)>=0 ? 1.0 : -1.0
    ca=_b3dot(x,sg.X);cb=_b3dot(x,sg.Y)
    u0=atan(cb,ca);u0<0 && (u0+=_OCC_TWO_PI)
    sgn<0 && (u0=_OCC_TWO_PI-u0)
    v=_b3dot(_b3(center,sg.center),sg.axis)
    if sg.occ===:cone
        # cone v is slant length along the generatrix, not axial height
        _,csa=_gm_sincos(_gm_atan((sg.r2-sg.r1)/sg.height))
        v/=csa
    end
    return (lin2d=((u0,v),(sgn,0.0)),)
end

# Pcurve(s) of a 3-D section circle on a sphere uv domain (u=lon,v=lat).
# A sphere section trace is representable as lin2d only in two cases:
#   latitude circle — the circle's plane ⊥ the polar axis (n∥axis): the
#     trace is a horizontal line u=u0+t·sgn at constant v;
#   meridian great circle — the circle's plane contains the axis and the
#     sphere center: the trace is two vertical half-lines at antipodal
#     longitudes, split at the poles. On each pole-to-pole half the latitude
#     is piecewise-linear in the circle parameter: v=π/2+δ−θ over
#     θ∈[δ,δ+π] and v=θ−δ−3π/2 over θ∈[δ+π,δ+2π], with
#     δ=atan2(y·axis,x·axis) the north-pole parameter.
# Every other section has a sinusoidal uv trace the model cannot represent —
# it must fail rather than approximate.
function _brep_sphere_section_pcs(sg,c,n,x,y,r,caller)
    axn=_b3dot(n,sg.axis)
    dc=_b3(c,sg.center)
    z=_b3dot(dc,sg.axis)
    rho=_b3norm(_b3(dc,_b3mul(sg.axis,z)))
    if abs(abs(axn)-1.0)<=_BREP_TOL && rho<=_BREP_TOL
        sgn=axn>=0 ? 1.0 : -1.0
        u0=_brep_surface_uv(sg,_b3add(c,_b3mul(x,r)))[1]
        sgn<0 && (u0=_OCC_TWO_PI-u0)
        v=atan(z,r)
        return NamedTuple[(lin2d=((u0,v),(sgn,0.0)),)]
    end
    if abs(axn)<=_BREP_TOL && rho<=_BREP_TOL &&
            abs(r-sg.radius)<=_BREP_TOL*max(1.0,sg.radius)
        xa=_b3dot(x,sg.axis);ya=_b3dot(y,sg.axis)
        δ=atan(ya,xa)
        out=NamedTuple[]
        for (tm,vo,vd) in ((δ+π/2,π/2+δ,-1.0),(δ+3π/2,-δ-3π/2,1.0))
            p=_b3add(c,_b3add(_b3mul(x,r*cos(tm)),_b3mul(y,r*sin(tm))))
            u=_brep_surface_uv(sg,p)[1]
            push!(out,(lin2d=((u,vo),(0.0,vd)),))
        end
        return out
    end
    throw(ArgumentError(
        "$caller: sphere section pcurve is not representable " *
        "(non-latitude, non-meridian circle)"))
end

# Is a section circle's sphere-domain pcurve lin2d-representable? — the
# same latitude/meridian conditions `_brep_sphere_section_pcs` enforces.
function _brep_sphere_pc_ok(sg,c,n,r)
    axn=_b3dot(n,sg.axis)
    dc=_b3(c,sg.center)
    z=_b3dot(dc,sg.axis)
    rho=_b3norm(_b3(dc,_b3mul(sg.axis,z)))
    return (abs(abs(axn)-1.0)<=_BREP_TOL && rho<=_BREP_TOL) ||
           (abs(axn)<=_BREP_TOL && rho<=_BREP_TOL &&
            abs(r-sg.radius)<=_BREP_TOL*max(1.0,sg.radius))
end

# line (parallel to axis) at angle u0 on a cylinder/cone domain.
function _brep_pc_line_on_cylinder(sg,o,d)
    p0=_brep_surface_uv(sg,o)
    sgn=_b3dot(d,sg.axis)>=0 ? 1.0 : -1.0
    return (lin2d=(p0,(0.0,sgn)),)
end

# ── tensor-patch × analytic-surface sections ────────────────────────────────
#
# A tensor B-spline/Bezier face meets an analytic face (plane, sphere,
# cylinder, cone, torus) in the zero contour of the partner's signed-
# distance field over the patch's (u,v) domain. The contour is traced by
# marching squares on a knot-aware grid (each crossing Newton-refined onto
# the level set), chained into components and fit by a cubic `spl2d`
# interpolant whose domain doubles as the section curve's parameter. The
# 3-D section curve evaluates S∘pc — it rides the patch exactly, which is
# what welds its endpoints onto the faces' wire vertices.

# Signed-distance field of an analytic face record: returns
# (fld(P)→distance, grd(P)→unit ∇fld). ∇ of a signed distance is the
# partner's geometric normal at the closest point — the direction the
# projected-gradient Newton correction needs. The zero set is the
# (extended) analytic surface; the fields stay continuous wherever a
# traced contour can reach, with an arbitrary unit fallback at the
# measure-zero gradient singularities (axis/centre/apex) where |fld|
# is bounded away from zero for non-degenerate records.
function _brep_partner_field(sg::NamedTuple)
    if sg.occ===:plane
        c=sg.center;n=sg.axis
        return (P->_b3dot(_b3(P,c),n), P->n)
    elseif sg.occ===:sphere
        c=sg.center;r=sg.radius
        return (P->(w=_b3(P,c);_b3norm(w)-r),
                P->(w=_b3(P,c);l=_b3norm(w);
                    l>eps() ? (w[1]/l,w[2]/l,w[3]/l) : sg.X))
    elseif sg.occ===:cylinder
        c=sg.center;a=sg.axis;r=sg.radius
        rad=P->(w=_b3(P,c);_b3(w,_b3mul(a,_b3dot(w,a))))
        return (P->(ρ=rad(P);_b3norm(ρ)-r),
                P->(ρ=rad(P);l=_b3norm(ρ);
                    l>eps() ? (ρ[1]/l,ρ[2]/l,ρ[3]/l) : sg.X))
    elseif sg.occ===:cone
        # In the (axial z, radial ρ) half-plane the cone's lateral
        # surface is the line through (0,r1)–(h,r2); the signed
        # distance to that line lifts to the signed distance to the
        # single nappe (ρ≥0 keeps the second nappe off the zero set).
        c=sg.center;a=sg.axis
        h=sg.height;dr=sg.r2-sg.r1
        L=sqrt(h*h+dr*dr)
        cz=h/L;cr=dr/L
        return (P->(w=_b3(P,c);z=_b3dot(w,a);
                    ρ=_b3norm(_b3(w,_b3mul(a,z)));
                    (ρ-sg.r1)*cz-z*cr),
                P->(w=_b3(P,c);z=_b3dot(w,a);
                    ρv=_b3(w,_b3mul(a,z));l=_b3norm(ρv);
                    ρ̂=l>eps() ? (ρv[1]/l,ρv[2]/l,ρv[3]/l) : sg.X;
                    (fma(cz,ρ̂[1],-cr*a[1]),
                     fma(cz,ρ̂[2],-cr*a[2]),
                     fma(cz,ρ̂[3],-cr*a[3]))))
    elseif sg.occ===:torus
        # distance to the spine circle minus the tube radius — a true
        # signed distance to the torus
        c=sg.center;a=sg.axis;r1=sg.r1;r2=sg.r2
        return (P->(w=_b3(P,c);z=_b3dot(w,a);
                    ρ=_b3norm(_b3(w,_b3mul(a,z)));
                    sqrt(fma(ρ-r1,ρ-r1,z*z))-r2),
                P->(w=_b3(P,c);z=_b3dot(w,a);
                    ρv=_b3(w,_b3mul(a,z));l=_b3norm(ρv);
                    if l<=eps()
                        return a
                    end;
                    s=(l-r1)/l;
                    _b3unit((fma(s,ρv[1],z*a[1]),
                             fma(s,ρv[2],z*a[2]),
                             fma(s,ρv[3],z*a[3])))))
    elseif _model_nurbs_surface(sg)
        # generic surface partner: project P onto the patch and take the
        # signed distance along the foot normal. The zero set is exactly
        # the patch (P on S ⇒ foot≈P ⇒ f≈0); off-surface the field is a
        # faithful signed distance, which is all the level-set tracer
        # needs.
        fld=P->begin
            u,v,_=_model_nurbs_surface_project(sg,P)
            F=_brep_surface_eval(sg,u,v)
            _b3dot(_b3(P,F),_brep_surface_normal(sg,u,v))
        end
        grd=P->begin
            u,v,_=_model_nurbs_surface_project(sg,P)
            _brep_surface_normal(sg,u,v)
        end
        return (fld,grd)
    end
    throw(ArgumentError(
        "_brep_partner_field: unsupported surface $(sg.occ)"))
end

# Move (u,v) onto the partner level set by projected-gradient Newton —
# the contour's correction direction is ∇f with ∇f=⟨(S_u,S_v),∇dist⟩.
# Returns (u,v,on_level).
function _brep_level_refine(sg,fld,grd,u,v,ftol,ulo,uhi,vlo,vhi)
    ok=false
    for _ in 1:50
        p,du,dv=_occ_surface_d1(sg,u,v)
        f=fld(p)
        if abs(f)<=ftol
            ok=true;break
        end
        n=grd(p)
        gu=_b3dot(du,n);gv=_b3dot(dv,n)
        g2=gu*gu+gv*gv
        g2<=0.0 && break
        λ=f/g2
        du_step=λ*gu;dv_step=λ*gv
        u=clamp(u-du_step,ulo,uhi);v=clamp(v-dv_step,vlo,vhi)
        hypot(du_step,dv_step)<=1e-14*max(1.0,hypot(u,v)) && break
    end
    p=_brep_surface_eval(sg,u,v)
    return (u,v,abs(fld(p))<=ftol)
end

# Marching-squares trace of f=0 on a surface domain. `us`/`vs` are the
# grid node coordinates (knot-densified); F[i,j]=f(us[i],vs[j]). `fld2`
# is the partner level field on the domain and `refine(u,v)` the Newton
# corrector. Nodes within `ztol` of zero snap to node-keyed crossings so
# adjacent cells weld through the shared vertex. `uper`/`vper` mark a
# periodic domain direction: the last grid line then coincides with the
# first and the chainer welds the seam instead of ending the contour
# there. Returns a list of components (verts, closed).
function _brep_msq_trace(sg,fld2,refine,us,vs,F,ztol,ftol,caller;
                       uper=0.0,vper=0.0)
    nu,nv=length(us),length(vs)
    ci(i)=uper>0.0 && i==nu ? 1 : i      # canonical node column
    cj(j)=vper>0.0 && j==nv ? 1 : j      # canonical node row
    # canonicalize a refined point sitting on a periodic seam line so
    # both adjacent cells record bitwise-identical crossing vertices
    function seam_uv(u,v)
        uper>0.0 && abs(u-us[nu])<=ztol && (u=us[1])
        vper>0.0 && abs(v-vs[nv])<=ztol && (v=vs[1])
        return (u,v)
    end
    # crossing record per grid edge: h-edges (i,j)–(i+1,j), v-edges (i,j)–(i,j+1)
    hcross=Matrix{Any}(nothing,nu-1,nv)
    vcross=Matrix{Any}(nothing,nu,nv-1)
    function edge_refine(ua,va,fa,ub,vb,fb)
        # bisect on the edge param, then Newton onto the level
        a=0.0;b=1.0;flo=fa
        for _ in 1:60
            m=(a+b)/2
            fm=fld2(fma(ua,1-m,ub*m),fma(va,1-m,vb*m))
            fm==0.0 && (a=b=m;break)
            signbit(fm)==signbit(flo) ? (a=m;flo=fm) : (b=m)
        end
        m=(a+b)/2
        u2,v2,_=refine(fma(ua,1-m,ub*m),fma(va,1-m,vb*m))
        return seam_uv(u2,v2)
    end
    for j in 1:(vper>0.0 ? nv-1 : nv), i in 1:nu-1
        fa=F[i,j];fb=F[i+1,j]
        if abs(fa)<=ztol
            hcross[i,j]=((:n,ci(i),cj(j)),(us[ci(i)],vs[cj(j)]))
        elseif abs(fb)<=ztol
            hcross[i,j]=((:n,ci(i+1),cj(j)),(us[ci(i+1)],vs[cj(j)]))
        elseif signbit(fa)!=signbit(fb)
            hcross[i,j]=((:h,i,cj(j)),
                edge_refine(us[i],vs[j],fa,us[i+1],vs[j],fb))
        end
    end
    for i in 1:(uper>0.0 ? nu-1 : nu), j in 1:nv-1
        fa=F[i,j];fb=F[i,j+1]
        if abs(fa)<=ztol
            vcross[i,j]=((:n,ci(i),cj(j)),(us[ci(i)],vs[cj(j)]))
        elseif abs(fb)<=ztol
            vcross[i,j]=((:n,ci(i),cj(j+1)),(us[ci(i)],vs[cj(j+1)]))
        elseif signbit(fa)!=signbit(fb)
            vcross[i,j]=((:v,ci(i),j),
                edge_refine(us[i],vs[j],fa,us[i],vs[j+1],fb))
        end
    end
    # per-cell boundary crossings, cyclic order B,R,T,L; pair into segments
    segs=Tuple{Any,NTuple{2,Float64},Any,NTuple{2,Float64}}[]
    for j in 1:nv-1, i in 1:nu-1
        cr=NTuple{3,Any}[]
        for e in (hcross[i,j],vcross[ci(i+1),j],hcross[i,cj(j+1)],vcross[i,j])
            e===nothing || push!(cr,(e[1],e[2],length(cr)))
        end
        length(cr)==0 && continue
        if length(cr)==3
            # a zero-node contact contributes two incidences — merge the
            # duplicate node entry so the count is even again
            di=findfirst(x->cr[x][1] isa Tuple && cr[x][1][1]===:n &&
                           any(y->y!=x && cr[y][1]==cr[x][1],
                               eachindex(cr)),eachindex(cr))
            di===nothing && throw(ErrorException(
                "$caller: tensor×analytic trace hit an odd crossing count"))
            cr=NTuple{3,Any}[cr[k] for k in eachindex(cr) if k!=di]
        end
        if length(cr)==2
            push!(segs,(cr[1][1],cr[1][2],cr[2][1],cr[2][2]))
        elseif length(cr)==4
            uc=(us[i]+us[i+1])/2;vc=(vs[j]+vs[j+1])/2
            fc=fld2(uc,vc)
            # cyclic entries (B,R,T,L): same sign as the BL corner pairs
            # (B,R)+(T,L), opposite pairs (B,L)+(R,T)
            same=signbit(fc)==signbit(F[i,j])
            if same
                push!(segs,(cr[1][1],cr[1][2],cr[2][1],cr[2][2]),
                           (cr[3][1],cr[3][2],cr[4][1],cr[4][2]))
            else
                push!(segs,(cr[1][1],cr[1][2],cr[4][1],cr[4][2]),
                           (cr[2][1],cr[2][2],cr[3][1],cr[3][2]))
            end
        else
            throw(ErrorException(
                "$caller: tensor×analytic trace hit an odd crossing count"))
        end
    end
    # A level set running exactly along a shared grid edge emits the same
    # segment (same node/edge key pair, same endpoints) from both adjacent
    # cells — dedup, else the chainer welds the coincident pair into an
    # out-and-back hairpin that later trips the zero-area wire guard.
    if !isempty(segs)
        seen=Dict{Tuple{Any,Any},Int}()
        keep=trues(length(segs))
        for (i,s) in enumerate(segs)
            kk=isless(s[1],s[3]) ? (s[1],s[3]) : (s[3],s[1])
            j=get(seen,kk,0)
            if j!=0
                o=segs[j]
                same_uv=(hypot(s[2][1]-o[2][1],s[2][2]-o[2][2])<=ztol &&
                         hypot(s[4][1]-o[4][1],s[4][2]-o[4][2])<=ztol) ||
                        (hypot(s[2][1]-o[4][1],s[2][2]-o[4][2])<=ztol &&
                         hypot(s[4][1]-o[2][1],s[4][2]-o[2][2])<=ztol)
                same_uv && (keep[i]=false;continue)
            end
            seen[kk]=i
        end
        segs=segs[keep]
    end
    return _brep_msq_chains(segs,caller)
end

# Chain cell segments through shared end keys (edge or node) into ordered
# components. Edge keys have degree ≤2; a node key can reach degree 4 when
# the level set runs through a shared grid vertex (saddle contact) — the
# four incidences then pair by opposite departure direction so both smooth
# branches continue through the vertex rather than turning.
function _brep_msq_chains(segs,caller)
    ends=Dict{Any,Vector{Tuple{Int,Int}}}()
    for (i,s) in enumerate(segs)
        push!(get!(Vector{Tuple{Int,Int}},ends,s[1]),(i,1))
        push!(get!(Vector{Tuple{Int,Int}},ends,s[3]),(i,2))
    end
    partner=Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for v in values(ends)
        n=length(v)
        n<=2 && continue
        n>4 && throw(ErrorException(
            "$caller: tensor×analytic section self-intersects at a grid vertex"))
        vuv=nothing;angs=Float64[]
        for (si,sd) in v
            s=segs[si]
            if sd==1
                vuv===nothing && (vuv=s[2])
                push!(angs,atan(s[4][2]-vuv[2],s[4][1]-vuv[1]))
            else
                vuv===nothing && (vuv=s[4])
                push!(angs,atan(s[2][2]-vuv[2],s[2][1]-vuv[1]))
            end
        end
        ord=sortperm(angs)
        if n==4
            for k in 1:2
                a=v[ord[k]];b=v[ord[k+2]]
                partner[a]=b;partner[b]=a
            end
        else
            # degree 3: a through-branch plus a tangent contact — pair the
            # most-opposite incidences; the odd one stays an open endpoint
            bi,bj,bd=1,2,0.0
            for a in 1:n, b in a+1:n
                d=abs(angs[a]-angs[b]);d=min(d,2π-d)
                d>bd && (bd=d;bi=a;bj=b)
            end
            partner[v[bi]]=v[bj];partner[v[bj]]=v[bi]
        end
    end
    used=falses(length(segs))
    comps=Tuple{Vector{NTuple{2,Float64}},Bool}[]
    for s0 in eachindex(segs)
        used[s0] && continue
        s=segs[s0]
        used[s0]=true
        verts=NTuple{2,Float64}[s[2],s[4]]
        closed=false
        # forward: exit through end-2's key; the last vertex equals the first
        # on a closed loop since shared keys carry identical uv
        key=s[3];cur=s0;curside=2
        while true
            nxt=get(partner,(cur,curside),nothing)
            if nxt===nothing
                for (si,sd) in get(ends,key,Tuple{Int,Int}[])
                    # a used segment belongs to an already-chained arc
                    # (except s0 — returning to the start seg IS the wrap
                    # that closes a ring); an incidence holding a planned
                    # pairing belongs to this node's through-branch — only
                    # an unpaired incidence (degree-2 node, or a degree-3
                    # spur) is a free continuation, and the spur correctly
                    # ends the arc
                    (si==cur || (used[si] && si!=s0) ||
                     haskey(partner,(si,sd))) && continue
                    nxt=(si,sd);break
                end
            end
            nxt===nothing && break
            nn,nd=nxt
            if nn==s0
                closed=true;break
            end
            used[nn] && throw(ErrorException(
                "$caller: tensor×analytic section chaining is inconsistent"))
            used[nn]=true
            ns=segs[nn]
            okey,ouv=nd==1 ? (ns[3],ns[4]) : (ns[1],ns[2])
            push!(verts,ouv);key=okey;cur=nn;curside=nd==1 ? 2 : 1
        end
        if !closed
            key=s[1];cur=s0;curside=1
            while true
                nxt=get(partner,(cur,curside),nothing)
                if nxt===nothing
                    for (si,sd) in get(ends,key,Tuple{Int,Int}[])
                        (si==cur || (used[si] && si!=s0) ||
                         haskey(partner,(si,sd))) && continue
                        nxt=(si,sd);break
                    end
                end
                nxt===nothing && break
                nn,nd=nxt
                if nn==s0
                    closed=true;break
                end
                used[nn] && throw(ErrorException(
                    "$caller: tensor×analytic section chaining is inconsistent"))
                used[nn]=true
                ns=segs[nn]
                okey,ouv=nd==1 ? (ns[3],ns[4]) : (ns[1],ns[2])
                pushfirst!(verts,ouv);key=okey;cur=nn;curside=nd==1 ? 2 : 1
            end
        end
        push!(comps,(verts,closed))
    end
    return comps
end

# Gauss–Newton projection of a 3-D point onto a tensor patch, seeded by the
# nearest traced sample — used only while the section is built (the stored
# lvl2d closure for the plane-side pcurve consults it per evaluation).
function _brep_nurbs_project_seed(sg,P,seeds_uv,seeds_pt,ulo,uhi,vlo,vhi)
    bi=1;bd=Inf
    for i in eachindex(seeds_pt)
        d=_b3dist(seeds_pt[i],P)
        d<bd && (bd=d;bi=i)
    end
    u,v=seeds_uv[bi]
    for _ in 1:60
        p,du,dv=_occ_surface_d1(sg,u,v)
        r=_b3(p,P)
        g1=_b3dot(r,du);g2=_b3dot(r,dv)
        a=_b3dot(du,du);b=_b3dot(du,dv);c=_b3dot(dv,dv)
        det=a*c-b*b
        det<=0.0 && break
        Δu=(g1*c-g2*b)/det;Δv=(g2*a-g1*b)/det
        u=clamp(u-Δu,ulo,uhi);v=clamp(v-Δv,vlo,vhi)
        hypot(Δu,Δv)<=1e-13*max(1.0,hypot(u,v)) && break
    end
    return u,v
end

# implicit level on the PARTNER's uv domain for a tensor×analytic section:
# signed distance from the lifted partner-surface point to the tensor
# surface. Returns (value, partner-uv gradient) — the gradient of a
# signed distance is the tensor normal at the foot, contracted with the
# partner's parametric partials.
function _brep_partner_lvl(sg,sgB,seeds_uv,seeds_pt,ulo,uhi,vlo,vhi,x,y)
    P,bu,bv=_occ_surface_d1(sgB,x,y)
    u,v=_brep_nurbs_project_seed(sg,P,seeds_uv,seeds_pt,ulo,uhi,vlo,vhi)
    p,du,dv=_occ_surface_d1(sg,u,v)
    n=_b3cross(du,dv)
    nn=_b3norm(n)
    nn==0.0 && return (Inf,(1.0,0.0))
    n=_b3mul(n,1/nn)
    r=_b3(P,p)
    g=_b3dot(r,n)
    return g,(_b3dot(bu,n),_b3dot(bv,n))
end

# Newton-refine a partner-uv point onto the zero level of
# `_brep_partner_lvl` — the partner-domain analog of `_brep_level_refine`.
function _brep_partner_lvl_refine(sg,sgB,seeds_uv,seeds_pt,ulo,uhi,vlo,vhi,
                                  x,y,ftol)
    for _ in 1:50
        g,(gx,gy)=_brep_partner_lvl(sg,sgB,seeds_uv,seeds_pt,
                                    ulo,uhi,vlo,vhi,x,y)
        abs(g)<=ftol && return (x,y,true)
        g2=gx*gx+gy*gy
        g2<=0.0 && break
        λ=g/g2
        x-=λ*gx;y-=λ*gy
        abs(λ)*hypot(gx,gy)<=1e-14*max(1.0,hypot(x,y)) && break
    end
    g,_=_brep_partner_lvl(sg,sgB,seeds_uv,seeds_pt,ulo,uhi,vlo,vhi,x,y)
    return (x,y,abs(g)<=ftol)
end

# Shift each uv's first coordinate by whole `uper` periods so consecutive
# samples stay within half a period — the partner-side pcurve of a
# contour crossing the u-seam is then a continuous spline in the
# unwrapped chart (periodic evaluators identify u mod `uper`, and the
# seam-aware machinery already unwraps projection pcurves the same way).
function _brep_unwrap_axis!(pts,per,ax)
    per<=0.0 && return pts
    for i in 2:length(pts)
        k=round((pts[i][ax]-pts[i-1][ax])/per)
        k!=0 && (pts[i]=ax==1 ?
            (pts[i][1]-k*per,pts[i][2]) :
            (pts[i][1],pts[i][2]-k*per))
    end
    return pts
end
_brep_unwrap_u!(pts,uper)=_brep_unwrap_axis!(pts,uper,1)
_brep_unwrap_v!(pts,vper)=_brep_unwrap_axis!(pts,vper,2)

# Contour tangent direction on the tensor domain at (u,v): the level
# field's parameter-space gradient is ⟨S_u,n⟩·∂u + ⟨S_v,n⟩·∂v with n the
# partner field gradient at S(u,v) — the contour runs perpendicular.
# Returns a NaN pair where the gradient vanishes (tangent intersection)
# so the caller's Hermite fit falls back to a secant slope there.
function _brep_contour_tangent(sg,grd,u,v)
    p,du,dv=_occ_surface_d1(sg,u,v)
    n=grd(p)
    fu=_b3dot(du,n);fv=_b3dot(dv,n)
    l=hypot(fu,fv)
    l<=eps() && return (NaN,NaN)
    return (-fv/l,fu/l)
end

# Joint cubic `spl2d` fit of one traced component in BOTH operand domains.
# The tensor-side trace is chord-timed in (u,v); each refined tensor-side
# sample lifts to the partner-domain uv exactly — the two interpolants
# share the section parameter, so pcA(t) and pcB(t) always name the same
# contour point (required for the engine's shared-parameter section
# convention). Both sides fit Hermite cubics on EXACT contour tangents —
# the estimated pchip slopes only converge ~O(h²) and cannot reach
# `tolfit` inside the point cap on a strongly curved trace. The tensor
# side's tangent comes from `_brep_contour_tangent`; the partner side's
# is the same 3-D contour velocity pushed through the partner
# parametrization (J_B⁺·J_A·dA), which keeps the shared-parameter
# contract exact at every knot. Intervals are split where the
# interpolants' midpoints fail to refine back onto the level set within
# `tolfit`. On u-periodic partners the pcB samples are unwrapped by
# `_brep_unwrap_u!` so the fit never interpolates across the seam.
# Returns (ts, pta, ptb, fitA, fitB) with fitA/fitB the (knots,deg,cps)
# triples. `tuper`/`tvper` are the TRACE domain's chart periods — the
# pca samples are unwrapped by them so a contour crossing the trace
# domain's own seam interpolates in the unwrapped chart.
function _brep_fit_section_pcs(sg,sgB,fld2,refine,grd,verts,closed,
                              ftol,tolfit,ulo,uhi,vlo,vhi,caller;
                              tuper=0.0,tvper=0.0)
    # drop consecutive near-duplicate vertices (below parameter resolution)
    pta=NTuple{2,Float64}[verts[1]]
    for q in verts
        hypot(q[1]-pta[end][1],q[2]-pta[end][2])>1e-14 &&
            push!(pta,q)
    end
    closed || length(pta)>=2 || return nothing
    closed && length(pta)<4 && return nothing
    _brep_unwrap_u!(pta,tuper)
    _brep_unwrap_v!(pta,tvper)
    function chords(ws)
        out=Float64[0.0]
        for i in 2:length(ws)
            push!(out,out[end]+hypot(ws[i][1]-ws[i-1][1],
                                     ws[i][2]-ws[i-1][2]))
        end
        return out
    end
    ts=chords(pta)
    ts[end]<=0.0 && return nothing
    uper=_brep_uperiod(sgB)
    vper=_brep_vperiod(sgB)
    ptb=_brep_unwrap_v!(_brep_unwrap_u!(NTuple{2,Float64}[
        _brep_surface_uv(sgB,_brep_surface_eval(sg,q...)) for q in pta],
        uper),vper)
    for it in 1:30
        n=length(pta)
        n>=2 || return nothing
        # exact contour tangents on the tensor side, oriented along the
        # chain and scaled by the local chord speed so they read as
        # d(uv)/dt under the chord parameter
        duv=Vector{NTuple{2,Float64}}(undef,n)
        prev=(NaN,NaN)
        for i in 1:n
            dir=_brep_contour_tangent(sg,grd,pta[i][1],pta[i][2])
            if isnan(dir[1])
                duv[i]=dir;continue
            end
            if i>1 && !isnan(prev[1]) &&
               dir[1]*prev[1]+dir[2]*prev[2]<0.0
                dir=(-dir[1],-dir[2])
            elseif i==1 && n>1
                sx,sy=pta[2][1]-pta[1][1],pta[2][2]-pta[1][2]
                dir[1]*sx+dir[2]*sy<0.0 && (dir=(-dir[1],-dir[2]))
            end
            prev=dir
            lo=i==1 ? 1 : i-1; hi=i==n ? n : i+1
            h=ts[hi]-ts[lo]
            mag=h>0.0 ?
                hypot(pta[hi][1]-pta[lo][1],pta[hi][2]-pta[lo][2])/h : 1.0
            duv[i]=(dir[1]*mag,dir[2]*mag)
        end
        if closed && n>=2
            # the ring's two ends are the same contour point — force the
            # shared tangent so the clamped spline closes C1 instead of
            # bulging off the level set at the last interval
            duv[n]=duv[1]
        end
        # partner-side tangents: push each tensor-side velocity through
        # the partner parametrization — duvB solves J_B·duvB = J_A·duvA
        # (the same 3-D contour velocity), so pcB'(t) is exact wherever
        # the partner chart is regular
        duvB=Vector{NTuple{2,Float64}}(undef,n)
        for i in 1:n
            dA=duv[i]
            if isnan(dA[1])
                duvB[i]=dA;continue
            end
            _,au,av=_occ_surface_d1(sg,pta[i][1],pta[i][2])
            V=(fma(dA[1],au[1],dA[2]*av[1]),
               fma(dA[1],au[2],dA[2]*av[2]),
               fma(dA[1],au[3],dA[2]*av[3]))
            _,bu,bv=_occ_surface_d1(sgB,ptb[i][1],ptb[i][2])
            a=_b3dot(bu,bu);b=_b3dot(bu,bv);c=_b3dot(bv,bv)
            det=a*c-b*b
            if !(det>eps()*a*c)  # singular partner chart (pole/apex)
                duvB[i]=(NaN,NaN);continue
            end
            r1=_b3dot(V,bu);r2=_b3dot(V,bv)
            duvB[i]=((c*r1-b*r2)/det,(a*r2-b*r1)/det)
        end
        closed && n>=2 && (duvB[n]=duvB[1])
        pca=_hermite_bspline_2d(ts,pta,duv)
        pcb=_hermite_bspline_2d(ts,ptb,duvB)
        bad=Int[]
        badt=Float64[]
        wf=wg=0.0;wfi=wgi=0
        for i in 1:n-1
            for frac in (0.3,0.5,0.7)
                tm=fma(ts[i+1]-ts[i],frac,ts[i])
                # tensor side: the interpolant's midpoint must land on the
                # partner within `tolfit` — |f| IS the 3-D distance to it
                uvm=_occ_spl2d_eval(pca[1],pca[2],pca[3],tm)
                rA=abs(fld2(uvm[1],uvm[2]))
                rA>wf && (wf=rA;wfi=i)
                rA>tolfit && (push!(bad,i);push!(badt,tm);break)
                # partner side: pcB(t) must name the SAME contour point as
                # pcA(t) (the shared-parameter contract) — compare the
                # lifted partner point to the tensor-side eval directly
                bvm=_occ_spl2d_eval(pcb[1],pcb[2],pcb[3],tm)
                Q=_brep_surface_eval(sgB,bvm[1],bvm[2])
                rB=_b3dist(Q,_brep_surface_eval(sg,uvm[1],uvm[2]))
                rB>wg && (wg=rB;wgi=i)
                rB>tolfit && (push!(bad,i);push!(badt,tm);break)
            end
        end
        isempty(bad) && return ts,pta,ptb,pca,pcb
        length(pta)+length(bad)>6000 && throw(ArgumentError(
            "$caller: tensor×analytic section trace did not converge to " *
            "the pcurve fit tolerance"))
        for k in length(bad):-1:1
            i=bad[k];tm=badt[k]
            uvm=_occ_spl2d_eval(pca[1],pca[2],pca[3],tm)
            u2,v2,ok=refine(uvm[1],uvm[2])
            ok || throw(ArgumentError(
                "$caller: tensor×analytic section trace is singular " *
                "(no local level-set point)"))
            splice!(pta,i+1:i,[(u2,v2)])
            splice!(ptb,i+1:i,
                [_brep_surface_uv(sgB,_brep_surface_eval(sg,u2,v2))])
            splice!(ts,i+1:i,[tm])
        end
        _brep_unwrap_v!(_brep_unwrap_u!(ptb,uper),vper)
        _brep_unwrap_v!(_brep_unwrap_u!(pta,tuper),tvper)
        ts=chords(pta)
    end
    throw(ArgumentError(
        "$caller: tensor×analytic section trace did not converge to " *
        "the pcurve fit tolerance"))
end

# Generic section-builder for a level-set trace: the zero contour of the
# partner's signed-distance field is marched over the trace face's uv
# domain, Newton-refined, and fitted into shared-parameter spl2d pcurves
# on both domains. `ts` is the trace-domain face record — a tensor
# patch or an analytic surface; `sgB` is the partner face record
# (plane/sphere/cylinder/cone/torus or another tensor patch). Emits
# (curve,closed,pcA,pcB) records like the analytic solvers;
# `sgB_is_b` selects which operand side the partner belongs to.
function _brep_sec_msq(ts,sgB,sgB_is_b,caller)
    isT=_model_nurbs_surface(ts)
    isT && hasproperty(ts,:reparam) && ts.reparam!==nothing &&
        throw(ArgumentError(
            "$caller: level-set section on a reparametrized patch is " *
            "unsupported"))
    fld,grd=_brep_partner_field(sgB)
    # trace-domain window + node grid: knot-densified for tensor patches,
    # a uniform subdivision for analytic charts
    if isT
        lo,hi=_model_nurbs_surface_bounds(ts)
        ulo,uhi,vlo,vhi=lo[1],hi[1],lo[2],hi[2]
        s=ts.patch.surface
        ku=unique(s.knots_u[s.degree_u+1:end-s.degree_u])
        kv=unique(s.knots_v[s.degree_v+1:end-s.degree_v])
        subu=clamp(div(160,length(ku)-1),6,24)
        subv=clamp(div(160,length(kv)-1),6,24)
        us=Float64[]
        for i in 1:length(ku)-1, j in 0:subu-1
            push!(us,fma(ku[i+1]-ku[i],j/subu,ku[i]))
        end
        push!(us,ku[end])
        vs=Float64[]
        for i in 1:length(kv)-1, j in 0:subv-1
            push!(vs,fma(kv[i+1]-kv[i],j/subv,kv[i]))
        end
        push!(vs,kv[end])
    else
        lo,hi=_occ_surface_bounds(ts)
        ulo,uhi,vlo,vhi=lo[1],hi[1],lo[2],hi[2]
        us=collect(range(ulo,uhi,97))
        vs=collect(range(vlo,vhi,97))
    end
    # domain wrap: a parametrization period only welds the domain's grid
    # lines when the face window actually spans the whole period — a
    # partial-revolution face's chart seam is a real boundary
    uper=_brep_uperiod(ts)
    uper = uper>0.0 && uhi-ulo>=uper-1e-9*max(1.0,uper) ? uper : 0.0
    vper=_brep_vperiod(ts)
    vper = vper>0.0 && vhi-vlo>=vper-1e-9*max(1.0,vper) ? vper : 0.0
    fld2=(u,v)->fld(_brep_surface_eval(ts,u,v))
    F=[fld2(u,v) for u in us, v in vs]
    fscl=max(1.0,maximum(abs,F))
    ztol=1e-9*fscl
    # a field below the zero tolerance on EVERY grid node is a coincident
    # surface pair (same-domain contact) — every cell would flood with
    # node crossings; reject like the other same-domain carriers
    maximum(abs,F)<=ztol && throw(ArgumentError(
        "$caller: coincident surfaces (same-domain) are not yet supported"))
    # model scale for the 3-D residual (Newton target): control-pole
    # bound on tensor patches, sampled-lift bound on analytic charts
    mscale=1.0
    if isT
        for c in ts.patch.surface.controls
            mscale=max(mscale,abs(c[1]),abs(c[2]),abs(c[3]))
        end
    else
        su=max(1,div(length(us),8));sv=max(1,div(length(vs),8))
        for i in 1:su:length(us), j in 1:sv:length(vs)
            p=_brep_surface_eval(ts,us[i],vs[j])
            mscale=max(mscale,abs(p[1]),abs(p[2]),abs(p[3]))
        end
    end
    ftol=1e-9*mscale
    # canonicalize unwrapped-chart seeds into the window before Newton —
    # a seed past the seam clamps to the domain boundary and loses the
    # level entirely (the field is periodic, the box clamp is not)
    function refine(u,v)
        uper>0.0 && (u=ulo+mod(u-ulo,uper))
        vper>0.0 && (v=vlo+mod(v-vlo,vper))
        return _brep_level_refine(ts,fld,grd,u,v,ftol,
                                  ulo,uhi,vlo,vhi)
    end
    comps=_brep_msq_trace(ts,fld2,refine,us,vs,F,ztol,ftol,caller;
                          uper=uper,vper=vper)
    out=NamedTuple[]
    tolfit=1e-7*mscale
    for (verts,closed) in comps
        length(verts)<2 && continue
        res=_brep_fit_section_pcs(ts,sgB,fld2,refine,grd,verts,closed,
                                  ftol,tolfit,ulo,uhi,vlo,vhi,caller;
                                  tuper=uper,tvper=vper)
        res===nothing && continue
        ts_,pta,ptb,(knotsA,pdegA,cpsA),(knotsB,pdegB,cpsB)=res
        seeds_uv=copy(pta)
        seeds_pt=NTuple{3,Float64}[_brep_surface_eval(ts,q...) for q in pta]
        lvlA=uv->fld2(uv[1],uv[2])
        lvlB=uv->_brep_partner_lvl(ts,sgB,seeds_uv,seeds_pt,
                                   ulo,uhi,vlo,vhi,uv[1],uv[2])[1]
        pca=(spl2d=(knotsA,pdegA,cpsA),lvl2d=lvlA)
        pcb=(spl2d=(knotsB,pdegB,cpsB),lvl2d=lvlB)
        # degree-1 seed polyline through the refined contour points: supplies
        # the NURBSCurve hull (knots/controls for panels, bounds, bbox) while
        # evaluation itself composes surface∘pcurve exactly
        pknots=Float64[ts_[1],ts_[1]]
        append!(pknots,ts_[2:end-1])
        push!(pknots,ts_[end],ts_[end])
        poly=NURBS.NURBSCurve(1,pknots,seeds_pt)
        ev=_ModelOCCSectionCurve(:bspline,ts,
                                 (spl2d=(knotsA,pdegA,cpsA),),poly)
        curve=(kind=:bspline,curve=ev,
               g=(occ=:copied_curve,curve=ev),ctype=:copied_curve)
        push!(out,(curve=curve,closed=closed,
                   pcA=sgB_is_b ? pca : pcb,
                   pcB=sgB_is_b ? pcb : pca))
    end
    return out
end
const _brep_sec_tensor_analytic=_brep_sec_msq

# Vertex anchoring for numerically-fitted section endpoints. A traced
# contour's fitted pcurve terminates a fit-residual off the wire edge it
# meets; projected per face that endpoint disagrees with the true pierce
# vertex by more than the uv node tolerance, so section arcs dangle and
# faces never split. Resolve each span end once in 3-D and consistently:
# prefer the already-computed edge pierce params (one authoritative vertex
# per edge crossing — every section landing there then shares bit-identical
# uv on the partner face), falling back to a projected wire point when the
# endpoint sits where no pierce was registered.
function _brep_sec_vertex(solids,oi,oj,esA,esB,sec,s,tol,mscale)
    p3=_brep_curve_eval(sec.curve,s)
    snap3=1e-4*max(1.0,mscale)
    recs=Any[nothing,nothing]
    for (side,(S,fidx,spc,es)) in enumerate((
            (solids[oi],sec.fa,sec.pcA,esA),
            (solids[oj],sec.fb,sec.pcB,esB)))
        f=S.faces[fidx]
        uv=_brep_pc_eval(spc,s)
        # canonicalize on this face's chart — a seam-crossing trace's
        # endpoint may sit a whole period off the wire pcurves' coordinates
        uper=_brep_uperiod(f.sg);vper=_brep_vperiod(f.sg)
        if uper>0 || vper>0
            uw1=uper>0 ?
                _brep_face_ubase(S.edges,f)+
                    mod(uv[1]-_brep_face_ubase(S.edges,f),uper) : uv[1]
            uv=(uw1,vper>0 ? mod(uv[2],vper) : uv[2])
        end
        uvtol=1e-4*max(1.0,hypot(uv[1],uv[2]))
        # per side: the authoritative vertex is this face's own wire-edge
        # pierce — its pcurve eval at the pierce param is bitwise the edge
        # split's node uv, where a projected candidate can sit the STEP
        # pcurve/edge inconsistency off and the section arc dangles
        bestd=Inf;res=nothing;proj=nothing;projd=Inf
        for wire in f.wires,(eidx,_,epc) in wire
            e=S.edges[eidx]
            e.curve.kind===:degenerate && continue
            for tp in get(es,eidx,Float64[])
                cand=_brep_curve_eval(e.curve,tp)
                d=_b3dist(p3,cand)
                d<bestd && (bestd=d;res=(cand,eidx,tp))
            end
            elo,ehi=minmax(e.t0,e.t1)
            te=_brep_pc_param(epc,uv,elo,ehi,1e-3)
            te===nothing && continue
            # uv-domain residual of the endpoint onto this edge — fitted
            # section endpoints are level-exact on the boundary in uv
            # while their 3-D image carries the trace's fit slack, so the
            # uv distance (not the 3-D one) decides the anchor
            qe=_brep_pc_eval(epc,te)
            duv=hypot(qe[1]-uv[1],qe[2]-uv[2])
            duv<projd || continue
            projd=duv
            proj=(_brep_curve_eval(e.curve,te),eidx,te,duv)
        end
        if get(ENV,"BREP_DEBUG_ARCS","")!=""
            println("DBG secvert s=$(round(s,digits=3)) side=$side " *
                    "uv=$(round.(uv,digits=3)) bestd=$(round(bestd,sigdigits=3)) " *
                    "resedge=$(res===nothing ? -1 : res[2]) " *
                    "projd=$(round(projd,sigdigits=3)) " *
                    "projedge=$(proj===nothing ? -1 : proj[2]) uvtol=$uvtol")
        end
        if res!==nothing && bestd<snap3
            recs[side]=res
        elseif proj!==nothing && projd<=uvtol
            # the endpoint sits on proj's edge in uv; if it also lies
            # within the section's 3-D slack of a registered pierce on
            # the same edge, the canonical pierce param keeps every
            # referencing face bitwise aligned
            if res!==nothing && res[2]==proj[2] &&
                    bestd<5e-3*mscale
                recs[side]=res
            else
                recs[side]=proj[1:3]
            end
        elseif proj!==nothing && _b3dist(p3,proj[1])<snap3
            recs[side]=proj[1:3]
        end
    end
    (recs[1]===nothing && recs[2]===nothing) && return nothing
    return (A=recs[1],B=recs[2])
end

# An anchored section span end's authoritative 3-D vertex — the operand-
# edge pierce evaluation the anchored side's wire pieces materialize
# bitwise. Unanchored endpoints keep the section curve's own evaluation.
_brep_sec_vfix(a)=a===nothing ? nothing :
    (a.A!==nothing ? a.A[1] : a.B[1])

# An anchored section endpoint's node uv on face `f` (which belongs to
# operand side `myside`): when the anchor is this operand's own wire-edge
# pierce, the edge pcurve eval at the pierce param reproduces the edge
# split's node uv bitwise — a surface projection of the 3-D vertex lands
# the pcurve/edge inconsistency off and the arc dangles. When only the
# partner side anchored (its edge pierces this face's interior), fall back
# to projecting its 3-D vertex.
function _brep_anchor_uv(f,myside,a)
    a===nothing && return nothing
    r=myside==1 ? a.A : a.B
    if r!==nothing
        for wire in f.wires,(eidx,_,epc) in wire
            eidx==r[2] || continue
            return _brep_pc_eval(epc,r[3])
        end
    end
    ro=myside==1 ? a.B : a.A
    ro===nothing && return nothing
    return _brep_surface_uv(f.sg,ro[1])
end

# did this endpoint anchor on `myside`'s own wire edge? — then the pierce
# param is already in that side's split table and no projected duplicate
# may be added
_brep_anchor_own(myside,a) = a!==nothing &&
    (myside==1 ? a.A : a.B)!==nothing

function _brep_face_sections(fa,fb,caller)
    ka,kb=fa.sg.occ,fb.sg.occ
    # identical curved carriers are a same-domain contact — their wires
    # coincide on one physical surface, so no section curve exists.
    # Planes go through the coincident-imprint path upstream; cylinders
    # and spheres have their own coaxial/concentric branches below.
    if ka===kb && ka in (:cone,:torus,:tensor_bspline,:tensor_bezier) &&
            _brep_sg_identical(fa.sg,fb.sg,_BREP_TOL)
        throw(ArgumentError(
            "$caller: coincident $ka faces (same-domain) are not yet supported"))
    end
    if ka===:plane && kb===:plane
        return _brep_sec_plane_plane(fa.sg,fb.sg,caller)
    elseif ka===:plane && kb===:cylinder
        return _brep_sec_plane_cylinder(fa.sg,fb.sg,true,caller)
    elseif ka===:cylinder && kb===:plane
        return _brep_sec_plane_cylinder(fb.sg,fa.sg,false,caller)
    elseif ka===:plane && kb===:sphere
        return _brep_sec_plane_sphere(fa.sg,fb.sg,true,caller)
    elseif ka===:sphere && kb===:plane
        return _brep_sec_plane_sphere(fb.sg,fa.sg,false,caller)
    elseif ka===:plane && kb===:cone
        return _brep_sec_plane_cone(fa.sg,fb.sg,true,caller)
    elseif ka===:cone && kb===:plane
        return _brep_sec_plane_cone(fb.sg,fa.sg,false,caller)
    elseif ka===:sphere && kb===:sphere
        return _brep_sec_sphere_sphere(fa.sg,fb.sg,caller)
    elseif ka===:cylinder && kb===:sphere
        return _brep_sec_cylinder_sphere(fa.sg,fb.sg,true,caller)
    elseif ka===:sphere && kb===:cylinder
        return _brep_sec_cylinder_sphere(fb.sg,fa.sg,false,caller)
    elseif ka===:cylinder && kb===:cylinder
        return _brep_sec_cylinder_cylinder(fa.sg,fb.sg,caller)
    end
    # generic level-set section: march the zero contour of one face's
    # signed-distance field over the other face's uv domain. The trace
    # domain is picked by chart quality — a tensor patch's knot-adapted
    # grid first, then the least-singular analytic chart; a plane's
    # unbounded uvb window never hosts a trace.
    _msq_cap(x)=x===:plane || x===:sphere || x===:cylinder ||
                x===:cone || x===:torus ||
                x===:tensor_bspline || x===:tensor_bezier
    if _msq_cap(ka) && _msq_cap(kb)
        _msq_rank(x)=x===:tensor_bspline||x===:tensor_bezier ? 0 :
                     x===:sphere ? 1 : x===:cylinder ? 2 :
                     x===:cone ? 3 : x===:torus ? 4 : 5
        if _msq_rank(ka)<=_msq_rank(kb)
            return _brep_sec_msq(fa.sg,fb.sg,true,caller)
        else
            return _brep_sec_msq(fb.sg,fa.sg,false,caller)
        end
    end
    throw(ArgumentError(
        "$caller: unsupported Boolean face pair $(ka)×$(kb)"))
end

function _brep_sec_plane_plane(pa,pb,caller)
    n=_b3cross(pa.axis,pb.axis);l=_b3norm(n)
    if l<_BREP_TOL
        d=_b3dot(_b3(pb.center,pa.center),pa.axis)
        abs(d)<_BREP_TOL && throw(ArgumentError(
            "$caller: coincident faces (same-domain) are not yet supported"))
        return NamedTuple[]
    end
    d=_b3unit(n)
    ca,cb=pa.center,pb.center
    na,nb=pa.axis,pb.axis
    # point on both planes: o = ca + λ1·na + λ2·nb with
    # na·(o-ca)=0 and nb·(o-cb)=0 → λ1 + a12·λ2 = 0, a12·λ1 + λ2 = nb·(cb-ca)
    a12=_b3dot(na,nb)
    det=1.0-a12*a12
    b2=_b3dot(_b3(cb,ca),nb)
    λ2=b2/det;λ1=-a12*λ2
    o=_b3add(ca,_b3add(_b3mul(na,λ1),_b3mul(nb,λ2)))
    curve=(kind=:line,o=o,d=d)
    return [(curve=curve,closed=false,
             pcA=_brep_pc_line_on_plane(pa,o,d),
             pcB=_brep_pc_line_on_plane(pb,o,d))]
end

function _brep_sec_plane_cylinder(pl,cy,plane_is_a,caller)
    dot=_b3dot(pl.axis,cy.axis)
    if abs(abs(dot)-1.0)<_BREP_TOL
        # plane ⊥ axis → circle
        axpt=_b3add(cy.center,
            _b3mul(cy.axis,_b3dot(_b3(pl.center,cy.center),cy.axis)))
        c=_b3(axpt,_b3mul(pl.axis,_b3dot(_b3(axpt,pl.center),pl.axis)))
        n=dot>0 ? pl.axis : (.-pl.axis)
        x,y=_occ_ax2_setx(n,cy.X)
        curve=(kind=:circle,c=c,n=n,x=x,y=y,r=cy.radius)
        pc_pl=_brep_pc_circ_on_plane(pl,c,n,x,cy.radius)
        pc_cy=_brep_pc_circ_on_cylinder(cy,c,n,x,cy.radius)
        return [(curve=curve,closed=true,
                 pcA=plane_is_a ? pc_pl : pc_cy,
                 pcB=plane_is_a ? pc_cy : pc_pl)]
    elseif abs(dot)<_BREP_TOL
        # plane ∥ axis → 0/1/2 axial lines
        dist=_b3dot(_b3(pl.center,cy.center),pl.axis)
        ad=abs(dist)
        ad>cy.radius+_BREP_TOL && return NamedTuple[]
        axpt=_b3add(cy.center,
            _b3mul(cy.axis,_b3dot(_b3(pl.center,cy.center),cy.axis)))
        foot=_b3add(axpt,_b3mul(pl.axis,dist))
        out=NamedTuple[]
        off=_b3unit(_b3cross(pl.axis,cy.axis))
        w=sqrt(max(0.0,cy.radius*cy.radius-dist*dist))
        for sgn in (w<=_BREP_TOL ? (1.0,) : (1.0,-1.0))
            o=_b3add(foot,_b3mul(off,sgn*w))
            pcA=plane_is_a ? _brep_pc_line_on_plane(pl,o,cy.axis) :
                             _brep_pc_line_on_cylinder(cy,o,cy.axis)
            pcB=plane_is_a ? _brep_pc_line_on_cylinder(cy,o,cy.axis) :
                             _brep_pc_line_on_plane(pl,o,cy.axis)
            push!(out,(curve=(kind=:line,o=o,d=cy.axis),closed=false,
                       pcA=pcA,pcB=pcB))
        end
        return out
    end
    # oblique section — an ellipse the analytic path can't emit; the
    # generic level-set section traces it on the cylinder domain
    return _brep_sec_msq(cy,pl,!plane_is_a,caller)
end

# Cartesian product of per-face pcurve lists into section entries.
function _brep_sec_entries(curve,closed,pcs_a,pcs_b,plane_is_a)
    out=NamedTuple[]
    for pa in pcs_a, pb in pcs_b
        push!(out,(curve=curve,closed=closed,
                   pcA=plane_is_a ? pa : pb,
                   pcB=plane_is_a ? pb : pa))
    end
    return out
end

function _brep_sec_plane_sphere(pl,sp,plane_is_a,caller)
    dist=_b3dot(_b3(pl.center,sp.center),pl.axis)
    ad=abs(dist)
    ad>sp.radius+_BREP_TOL && return NamedTuple[]
    ad>=sp.radius-_BREP_TOL && throw(ArgumentError(
        "$caller: tangent plane×sphere contact is unsupported"))
    r=sqrt(sp.radius*sp.radius-dist*dist)
    c=_b3add(sp.center,_b3mul(pl.axis,dist))
    n=pl.axis
    x,y=_occ_ax2_frame(n)
    # a tilted-axis sphere (imported/rotated record) sees a non-latitude,
    # non-meridian circle — its uv trace isn't lin2d, so trace on the
    # sphere chart instead of rejecting
    _brep_sphere_pc_ok(sp,c,n,r) ||
        return _brep_sec_msq(sp,pl,!plane_is_a,caller)
    curve=(kind=:circle,c=c,n=n,x=x,y=y,r=r)
    pc_pl=NamedTuple[_brep_pc_circ_on_plane(pl,c,n,x,r)]
    pc_sp=_brep_sphere_section_pcs(sp,c,n,x,y,r,caller)
    return _brep_sec_entries(curve,true,pc_pl,pc_sp,plane_is_a)
end

function _brep_sec_plane_cone(pl,cn,plane_is_a,caller)
    dot=_b3dot(pl.axis,cn.axis)
    if abs(abs(dot)-1.0)<_BREP_TOL
        v=_b3dot(_b3(pl.center,cn.center),cn.axis)
        rv=cn.r1+v*(cn.r2-cn.r1)/cn.height
        rv<-_BREP_TOL && return NamedTuple[]
        rv<=_BREP_TOL && throw(ArgumentError(
            "$caller: plane×cone apex contact is unsupported"))
        axpt=_b3add(cn.center,_b3mul(cn.axis,v))
        c=_b3(axpt,_b3mul(pl.axis,_b3dot(_b3(axpt,pl.center),pl.axis)))
        n=dot>0 ? pl.axis : (.-pl.axis)
        x,y=_occ_ax2_setx(n,cn.X)
        curve=(kind=:circle,c=c,n=n,x=x,y=y,r=rv)
        pc_pl=_brep_pc_circ_on_plane(pl,c,n,x,rv)
        pc_cn=_brep_pc_circ_on_cylinder(cn,c,n,x,rv)
        return [(curve=curve,closed=true,
                 pcA=plane_is_a ? pc_pl : pc_cn,
                 pcB=plane_is_a ? pc_cn : pc_pl)]
    end
    # non-axial section — a conic the analytic path can't emit; the
    # generic level-set section traces it on the cone domain
    return _brep_sec_msq(cn,pl,!plane_is_a,caller)
end

function _brep_sec_sphere_sphere(sa,sb,caller)
    d=_b3(sb.center,sa.center);dist=_b3norm(d)
    if dist<_BREP_TOL
        # concentric: equal radii → coincident surfaces (same-domain
        # contact, no section); different radii never touch
        abs(sa.radius-sb.radius)<=_BREP_TOL && throw(ArgumentError(
            "$caller: coincident spheres (same-domain) are not yet supported"))
        return NamedTuple[]
    end
    (dist>sa.radius+sb.radius-_BREP_TOL ||
     dist<abs(sa.radius-sb.radius)+_BREP_TOL) && return NamedTuple[]
    a=(sa.radius^2-sb.radius^2+dist^2)/(2dist)
    h2=sa.radius^2-a^2
    h2<=_BREP_TOL && throw(ArgumentError(
        "$caller: tangent sphere×sphere contact is unsupported"))
    n=_b3mul(d,1/dist)
    c=_b3add(sa.center,_b3mul(n,a))
    r=sqrt(h2)
    x,y=_occ_ax2_frame(n)
    if !(_brep_sphere_pc_ok(sa,c,n,r) && _brep_sphere_pc_ok(sb,c,n,r))
        # the section circle's pcurve isn't lin2d-representable on at
        # least one sphere — trace it on the first sphere's domain
        return _brep_sec_msq(sa,sb,true,caller)
    end
    curve=(kind=:circle,c=c,n=n,x=x,y=y,r=r)
    pcA=_brep_sphere_section_pcs(sa,c,n,x,y,r,caller)
    pcB=_brep_sphere_section_pcs(sb,c,n,x,y,r,caller)
    out=NamedTuple[]
    for pa in pcA,pb in pcB
        push!(out,(curve=curve,closed=true,pcA=pa,pcB=pb))
    end
    return out
end

function _brep_sec_cylinder_sphere(cy,sp,cy_is_a,caller)
    # the analytic path covers the coaxial case → circles; anything else
    # traces on the cylinder domain
    d=_b3(sp.center,cy.center)
    axial=_b3dot(d,cy.axis)
    radial=_b3(d,_b3mul(cy.axis,axial))
    _b3norm(radial)>_BREP_TOL &&
        return _brep_sec_msq(cy,sp,cy_is_a,caller)
    dz2=sp.radius^2-cy.radius^2
    # the sphere reaches the lateral wall only when r > R — a smaller
    # on-axis sphere never meets it; r ≈ R is tangent contact
    dz2<-_BREP_TOL && return NamedTuple[]
    dz2<=_BREP_TOL && throw(ArgumentError(
        "$caller: tangent cylinder×sphere contact is unsupported"))
    dz=sqrt(dz2)
    # the sphere-side pcurve is lin2d only when the sphere's own axis
    # lines up with the circle plane — otherwise trace on the cylinder
    for sgn in (1.0,-1.0)
        c=_b3add(cy.center,_b3mul(cy.axis,axial+sgn*dz))
        _brep_sphere_pc_ok(sp,c,cy.axis,cy.radius) ||
            return _brep_sec_msq(cy,sp,cy_is_a,caller)
    end
    out=NamedTuple[]
    for sgn in (1.0,-1.0)
        z0=axial+sgn*dz
        c=_b3add(cy.center,_b3mul(cy.axis,z0))
        n=cy.axis
        x,y=cy.X,cy.Y
        curve=(kind=:circle,c=c,n=n,x=x,y=y,r=cy.radius)
        pc_cy=_brep_pc_circ_on_cylinder(cy,c,n,x,cy.radius)
        for pc_sp in _brep_sphere_section_pcs(sp,c,n,x,y,cy.radius,caller)
            push!(out,(curve=curve,closed=true,
                       pcA=cy_is_a ? pc_cy : pc_sp,
                       pcB=cy_is_a ? pc_sp : pc_cy))
        end
    end
    return out
end

function _brep_sec_cylinder_cylinder(ca,cb,caller)
    if _b3norm(_b3cross(ca.axis,cb.axis))<_BREP_TOL
        d=_b3(cb.center,ca.center)
        proj=_b3(d,_b3mul(ca.axis,_b3dot(d,ca.axis)))
        rd=_b3norm(proj)
        if rd<_BREP_TOL
            # coaxial: equal radii → coincident surfaces (same-domain
            # contact, no section); different radii never touch
            abs(ca.radius-cb.radius)<=_BREP_TOL && throw(ArgumentError(
                "$caller: coaxial cylinders (same-domain) are not yet " *
                "supported"))
            return NamedTuple[]
        end
        # parallel off-axis: the section is axial lines where the
        # lateral fields meet — the generic trace finds them; a
        # single-line tangent contact is a degenerate root the trace
        # cannot certify
        (abs(rd-(ca.radius+cb.radius))<=_BREP_TOL ||
         abs(rd-abs(ca.radius-cb.radius))<=_BREP_TOL) &&
            throw(ArgumentError(
                "$caller: tangent cylinder×cylinder contact is unsupported"))
        return _brep_sec_msq(ca,cb,true,caller)
    end
    # non-parallel axes — the generic level-set section traces the
    # curve on the first cylinder's domain
    return _brep_sec_msq(ca,cb,true,caller)
end

# ── 2-D pcurve × pcurve intersections ───────────────────────────────────────
#
# Intersect a section pcurve (infinite line or full circle) with a boundary
# edge's pcurve (over the edge parameter range [et0,et1]); return
# section-parameter values.

function _brep_pc_intersect(sec_pc,edge_pc,et0,et1,tol;
                            sec_crv=nothing,edge_crv=nothing,
                            pu=0.0,pv=0.0)
    # Unwrap neg2d pcurves: the edge parameter mirrors the inner parameter,
    # so the reported span reflects to (−et1,−et0) before dispatch. A negated
    # section pcurve reports mirrored section parameters back to the caller.
    sec_pc,nsec=_brep_pc_unneg(sec_pc)
    edge_pc,nedge=_brep_pc_unneg(edge_pc)
    nedge && ((et0,et1)=(-et1,-et0))
    out=if hasproperty(sec_pc,:lin2d) && hasproperty(edge_pc,:lin2d)
        _brep_ix_lin_lin(sec_pc.lin2d,edge_pc.lin2d,et0,et1,tol;pu=pu,pv=pv)
    elseif hasproperty(sec_pc,:lin2d) && hasproperty(edge_pc,:circ2d)
        _brep_ix_lin_circ(sec_pc.lin2d,edge_pc.circ2d,et0,et1,tol)
    elseif hasproperty(sec_pc,:circ2d) && hasproperty(edge_pc,:lin2d)
        _brep_ix_circ_lin(sec_pc.circ2d,edge_pc.lin2d,et0,et1,tol)
    elseif hasproperty(sec_pc,:circ2d) && hasproperty(edge_pc,:circ2d)
        _brep_ix_circ_circ(sec_pc.circ2d,edge_pc.circ2d,et0,et1,tol)
    elseif hasproperty(edge_pc,:lin2d)
        # generic trace × line-pcurve edge: solve on the section's own
        # chart — more robust than the level-field + inversion route for
        # fitted traces whose partner-projection degrades near clamps.
        # Elementary carriers (ellips2d) report an unbounded domain —
        # their exact implicit level keeps them on the generic solver.
        dlo,dhi=_brep_pc_domain(sec_pc)
        if isfinite(dlo)&&isfinite(dhi)
            _brep_ix_sec_lin(sec_pc,edge_pc.lin2d,et0,et1,tol;pu=pu,pv=pv)
        else
            _brep_ix_generic(sec_pc,edge_pc,et0,et1,tol;
                             sec_crv=sec_crv,edge_crv=edge_crv,pu=pu,pv=pv)
        end
    else
        _brep_ix_generic(sec_pc,edge_pc,et0,et1,tol;
                         sec_crv=sec_crv,edge_crv=edge_crv,pu=pu,pv=pv)
    end
    return nsec ? .-out : out
end

# snap a circle parameter into [lo,hi] modulo 2π; nothing when out of range.
function _brep_snap_theta(t,lo,hi,tol)
    base=mod(t,_OCC_TWO_PI)
    for k in -2:2
        tk=base+k*_OCC_TWO_PI
        (tk>=lo-tol && tk<=hi+tol) && return tk
    end
    return nothing
end

function _brep_ix_lin_lin(sl,el,et0,et1,tol;pu=0.0,pv=0.0)
    so,sd=sl;eo,ed=el
    det=sd[1]*ed[2]-sd[2]*ed[1]
    abs(det)<tol && return Float64[]
    lo,hi=minmax(et0,et1)
    out=Float64[]
    # the section line may run unwrapped on a periodic domain — try whole-
    # period translates of its origin so a canonical edge span still matches
    for ku in (pu>0 ? (-1:1) : (0:0)), kv in (pv>0 ? (-1:1) : (0:0))
        o2=(so[1]+ku*pu,so[2]+kv*pv)
        w=(eo[1]-o2[1],eo[2]-o2[2])
        t=(w[1]*ed[2]-w[2]*ed[1])/det
        s=(w[1]*sd[2]-w[2]*sd[1])/det
        (s<lo-tol || s>hi+tol) && continue
        isempty(out) || abs(t-out[end])<=tol || push!(out,t)
        isempty(out) && push!(out,t)
    end
    return out
end

# params on a generic section trace where it crosses a boundary edge whose
# pcurve is a straight uv line — solved on the SECTION's own chart as roots
# of the line's implicit equation, so no uv inversion or level field is
# needed. Whole-period translates cover a seam-crossing trace's unwrapped
# coordinates: copies parallel to the edge leave the line invariant (they
# only matter to the edge-span check), copies normal to it shift the roots.
function _brep_ix_sec_lin(sec_pc,el,et0,et1,tol;pu=0.0,pv=0.0)
    eo,ed=el
    dd=ed[1]*ed[1]+ed[2]*ed[2]
    dd<=0.0 && return Float64[]
    lo,hi=minmax(et0,et1)
    dlo,dhi=_brep_pc_domain(sec_pc)
    (isfinite(dlo)&&isfinite(dhi)) || return Float64[]
    ts=_brep_pc_grid(sec_pc,dlo,dhi)
    umin,umax,vmin,vmax=Inf,-Inf,Inf,-Inf
    for tt in ts
        q=_brep_pc_eval(sec_pc,tt)
        umin=min(umin,q[1]);umax=max(umax,q[1])
        vmin=min(vmin,q[2]);vmax=max(vmax,q[2])
    end
    # edge uv extent over its span (linear → ends suffice)
    e1=(eo[1]+lo*ed[1],eo[2]+lo*ed[2])
    e2=(eo[1]+hi*ed[1],eo[2]+hi*ed[2])
    eu_lo,eu_hi=minmax(e1[1],e2[1]);ev_lo,ev_hi=minmax(e1[2],e2[2])
    # copies of the line in u: distinct only when the edge is not u-parallel
    un=abs(ed[2])>1e-12*sqrt(dd)          # u-shift moves the line
    vn=abs(ed[1])>1e-12*sqrt(dd)          # v-shift moves the line
    kus=(pu>0 && un) ?
        (floor(Int,(umin-eu_hi)/pu)-1:ceil(Int,(umax-eu_lo)/pu)+1) : (0:0)
    kvs=(pv>0 && vn) ?
        (floor(Int,(vmin-ev_hi)/pv)-1:ceil(Int,(vmax-ev_lo)/pv)+1) : (0:0)
    # validation copies cover both axes fully (a u-parallel edge's copies
    # coincide as lines but its span still lives on the canonical band)
    kusV=pu>0 ?
        (floor(Int,(umin-eu_hi)/pu)-1:ceil(Int,(umax-eu_lo)/pu)+1) : (0:0)
    kvsV=pv>0 ?
        (floor(Int,(vmin-ev_hi)/pv)-1:ceil(Int,(vmax-ev_lo)/pv)+1) : (0:0)
    ns=sqrt(dd)
    out=Float64[]
    for ku in kus, kv in kvs
        δu,δv=ku*pu,kv*pv
        hits=_brep_scalar_roots(
            t->(q=_brep_pc_eval(sec_pc,t);
                ((q[1]-δu-eo[1])*ed[2]-(q[2]-δv-eo[2])*ed[1])/ns),
            ts,tol)
        for t in hits
            for ku2 in kusV, kv2 in kvsV
                uw=(let q=_brep_pc_eval(sec_pc,t);
                    (q[1]-ku2*pu,q[2]-kv2*pv);end)
                s=((uw[1]-eo[1])*ed[1]+(uw[2]-eo[2])*ed[2])/dd
                (s<lo-tol || s>hi+tol) && continue
                r=(uw[1]-fma(s,ed[1],eo[1]),uw[2]-fma(s,ed[2],eo[2]))
                hypot(r[1],r[2])<=
                    1000tol*max(1.0,abs(uw[1]),abs(uw[2])) || continue
                push!(out,t);break
            end
        end
    end
    sort!(out)
    ded=Float64[]
    for t in out
        (isempty(ded) || t-ded[end]>tol) && push!(ded,t)
    end
    return ded
end

# params on the LINE where it meets the arc.
function _brep_ix_lin_circ(sl,ec,et0,et1,tol)
    so,sd=sl;c,x,r=ec
    rel=(so[1]-c[1],so[2]-c[2])
    a=sd[1]*x[1]+sd[2]*x[2]
    b=sd[2]*x[1]-sd[1]*x[2]
    e=rel[1]*x[1]+rel[2]*x[2]
    f=rel[2]*x[1]-rel[1]*x[2]
    A=a*a+b*b;B=2(e*a+f*b);C=e*e+f*f-r*r
    A<tol && return Float64[]
    disc=B*B-4A*C
    disc<0 && return Float64[]
    out=Float64[];sq=sqrt(max(0.0,disc))
    lo,hi=minmax(et0,et1)
    for t in ((-B+sq)/(2A),(-B-sq)/(2A))
        uv=(so[1]+t*sd[1],so[2]+t*sd[2])
        cx=(uv[1]-c[1])*x[1]+(uv[2]-c[2])*x[2]
        cy=(uv[2]-c[2])*x[1]-(uv[1]-c[1])*x[2]
        θ=_brep_snap_theta(atan(cy,cx),lo,hi,tol)
        θ===nothing || push!(out,t)
    end
    return out
end

# params on the CIRCLE (section) where it meets the line segment.
function _brep_ix_circ_lin(sc,el,et0,et1,tol)
    c,x,r=sc;eo,ed=el
    rel=(eo[1]-c[1],eo[2]-c[2])
    a=ed[1]*x[1]+ed[2]*x[2]
    b=ed[2]*x[1]-ed[1]*x[2]
    e=rel[1]*x[1]+rel[2]*x[2]
    f=rel[2]*x[1]-rel[1]*x[2]
    A=a*a+b*b;B=2(e*a+f*b);C=e*e+f*f-r*r
    A<tol && return Float64[]
    disc=B*B-4A*C
    disc<0 && return Float64[]
    out=Float64[];sq=sqrt(max(0.0,disc))
    lo,hi=minmax(et0,et1)
    for s in ((-B+sq)/(2A),(-B-sq)/(2A))
        (s<lo-tol || s>hi+tol) && continue
        uv=(eo[1]+s*ed[1],eo[2]+s*ed[2])
        cx=(uv[1]-c[1])*x[1]+(uv[2]-c[2])*x[2]
        cy=(uv[2]-c[2])*x[1]-(uv[1]-c[1])*x[2]
        push!(out,mod(atan(cy,cx),_OCC_TWO_PI))
    end
    return out
end

# params on the FIRST (section) circle where it meets the arc.
function _brep_ix_circ_circ(sc,ec,et0,et1,tol)
    c1,x1,r1=sc;c2,x2,r2=ec
    d=(c2[1]-c1[1],c2[2]-c1[2]);dist=sqrt(d[1]^2+d[2]^2)
    (dist>r1+r2-tol || dist<abs(r1-r2)+tol || dist<tol) && return Float64[]
    a=(r1*r1-r2*r2+dist*dist)/(2dist)
    h2=r1*r1-a*a
    h2<0 && return Float64[]
    h=sqrt(max(0.0,h2))
    ux=(d[1]*x1[1]+d[2]*x1[2])/dist
    uy=(d[2]*x1[1]-d[1]*x1[2])/dist
    out=Float64[]
    lo,hi=minmax(et0,et1)
    for sgn in (1.0,-1.0)
        bx=a*ux-sgn*h*uy; by=a*uy+sgn*h*ux
        t=atan(by,bx)
        uv=_occ_circ2d_eval(c1,x1,r1,t)
        cx=(uv[1]-c2[1])*x2[1]+(uv[2]-c2[2])*x2[2]
        cy=(uv[2]-c2[2])*x2[1]-(uv[1]-c2[1])*x2[2]
        θ=_brep_snap_theta(atan(cy,cx),lo,hi,tol)
        θ===nothing || push!(out,mod(t,_OCC_TWO_PI))
    end
    return out
end

# parameter on `pc` of the point nearest `uv`, over the natural domain or
# [lo,hi] when supplied — analytic for elementary records, the seeded Newton
# projection for splines and conic-reparam traces.
function _brep_pc_uv_to_t(pc,uv,tol,lo=-Inf,hi=Inf)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        dd=d[1]*d[1]+d[2]*d[2]
        dd==0.0 && return nothing
        t=((uv[1]-o[1])*d[1]+(uv[2]-o[2])*d[2])/dd
        (t<lo-tol || t>hi+tol) && return nothing
        return t
    elseif hasproperty(pc,:circ2d)
        c,x,r=pc.circ2d
        cx=(uv[1]-c[1])*x[1]+(uv[2]-c[2])*x[2]
        cy=(uv[2]-c[2])*x[1]-(uv[1]-c[1])*x[2]
        return _brep_snap_theta(atan(cy,cx),lo,hi,tol)
    elseif hasproperty(pc,:ellips2d)
        o,x,a,b=pc.ellips2d
        cx=((uv[1]-o[1])*x[1]+(uv[2]-o[2])*x[2])/a
        cy=((uv[2]-o[2])*x[1]-(uv[1]-o[1])*x[2])/b
        return _brep_snap_theta(atan(cy,cx),lo,hi,tol)
    end
    dlo,dhi=_brep_pc_domain(pc)
    _brep_pc_project(pc,uv,max(lo,dlo),min(hi,dhi),tol)
end

# hits of an edge pcurve against a generic section pcurve — the distance
# field edge(t)→dist(t,sec) has no sign change, so local minima below
# tolerance on the sample grid get a golden-section refinement instead.
function _brep_pc_curve_hits(sec_pc,edge_pc,ts,tol)
    slo,shi=_brep_pc_domain(sec_pc)
    slo=isinf(slo) ? -4*_OCC_TWO_PI : slo
    shi=isinf(shi) ?  4*_OCC_TWO_PI : shi
    function dist(t)
        uv=_brep_pc_eval(edge_pc,t)
        s=_brep_pc_project(sec_pc,uv,slo,shi,tol)
        s===nothing && return Inf
        p=_brep_pc_eval(sec_pc,s)
        return hypot(p[1]-uv[1],p[2]-uv[2])
    end
    fs=[dist(Float64(t)) for t in ts]
    scl=max(1.0,maximum(x->isfinite(x) ? x : 0.0,fs))
    hits=Float64[];φ=(sqrt(5)-1)/2
    for i in 2:length(ts)-1
        (isfinite(fs[i]) && fs[i]<=fs[i-1] && fs[i]<=fs[i+1] &&
         fs[i]<10tol*scl) || continue
        a,b=ts[i-1],ts[i+1]
        c=b-φ*(b-a);d=a+φ*(b-a);fc=dist(c);fd=dist(d)
        for _ in 1:60
            if fc<fd
                b=d;d=c;fd=fc;c=b-φ*(b-a);fc=dist(c)
            else
                a=c;c=d;fc=fd;d=a+φ*(b-a);fd=dist(d)
            end
        end
        tm=(a+b)/2
        if dist(tm)<10tol*scl &&
                (isempty(hits) || tm-hits[end]>1e-9*max(1.0,b-a))
            push!(hits,tm)
        end
    end
    return hits
end

# generic pair solver — the edge pcurve stays the bounded operand: march
# its own grid for crossings of the section's implicit level (or its raw
# distance field when the section has none), then invert each crossing uv
# back onto the section. `pu`/`pv` are the face domain's chart periods: a
# seam-crossing trace keeps unwrapped coordinates while every wire pcurve
# stays canonical, so an inversion that fails by a whole period is retried
# on the matching uv translate.
function _brep_ix_generic(sec_pc,edge_pc,et0,et1,tol;
                          sec_crv=nothing,edge_crv=nothing,
                          pu=0.0,pv=0.0)
    lo,hi=minmax(et0,et1)
    ts=_brep_pc_grid(edge_pc,lo,hi)
    lvl=_brep_pc_level(sec_pc)
    if lvl!==nothing
        hits=_brep_scalar_roots(t->lvl(_brep_pc_eval(edge_pc,t)),ts,tol)
    else
        hits=_brep_pc_curve_hits(sec_pc,edge_pc,ts,tol)
    end
    out=Float64[]
    # fitted traces (spl2d carrying lvl2d) deviate from the true contour by
    # their fit residual — in a distorted uv metric the gap to a wire
    # crossing defeats the nominal inversion tolerance, so such sections
    # get a second, slack-bounded inversion pass
    fitted=lvl!==nothing && hasproperty(sec_pc,:spl2d)
    ext=NaN
    ubox=nothing
    for te in hits
        uv=_brep_pc_eval(edge_pc,te)
        st=_brep_pc_uv_to_t(sec_pc,uv,tol)
        if st===nothing && (fitted || pu>0 || pv>0)
            dlo,dhi=_brep_pc_domain(sec_pc)
            if isfinite(dlo)&&isfinite(dhi)
                # candidate uv translates: the canonical hit plus whole-
                # period shifts that reach the trace's unwrapped band
                if ubox===nothing
                    umin,umax,vmin,vmax=Inf,-Inf,Inf,-Inf
                    for tt in _brep_pc_grid(sec_pc,dlo,dhi;subs=2)
                        q=_brep_pc_eval(sec_pc,tt)
                        umin=min(umin,q[1]);umax=max(umax,q[1])
                        vmin=min(vmin,q[2]);vmax=max(vmax,q[2])
                    end
                    ubox=(umin,umax,vmin,vmax)
                end
                umin,umax,vmin,vmax=ubox
                kulo=pu>0 ? floor(Int,(uv[1]-umax)/pu)-1 : 0
                kuhi=pu>0 ? ceil(Int,(uv[1]-umin)/pu)+1 : 0
                kvlo=pv>0 ? floor(Int,(uv[2]-vmax)/pv)-1 : 0
                kvhi=pv>0 ? ceil(Int,(uv[2]-vmin)/pv)+1 : 0
                for ku in kulo:kuhi, kv in kvlo:kvhi
                    (ku==0 && kv==0) && continue
                    st=_brep_pc_uv_to_t(
                        sec_pc,(uv[1]+ku*pu,uv[2]+kv*pv),tol)
                    st===nothing || break
                end
                if st===nothing && fitted
                    isnan(ext) &&
                        (ext=_brep_pc_uv_extent(sec_pc,dlo,dhi))
                    best=nothing;bd=Inf
                    for ku in kulo:kuhi, kv in kvlo:kvhi
                        pr=_brep_pc_project_res(
                            sec_pc,(uv[1]+ku*pu,uv[2]+kv*pv),dlo,dhi)
                        pr===nothing && continue
                        pr[2]<bd && (bd=pr[2];best=pr[1])
                    end
                    if best!==nothing && bd<=max(0.05*ext,1000tol)
                        st=best
                    elseif sec_crv!==nothing && edge_crv!==nothing
                        # near-degenerate uv metric or a residual just past
                        # the slack bound: verify in 3-D — the wire pierce
                        # lies on the traced curve within its fit residual,
                        # a different component stays macroscopically far
                        V=_brep_curve_eval(edge_crv,te)
                        vtol=2e-3*max(1.0,_b3norm(V))
                        if best!==nothing &&
                                _b3dist(_brep_curve_eval(sec_crv,best),V)<=vtol
                            st=best
                        else
                            t2=_brep_curve_param(sec_crv,V,0.0)
                            if t2!==nothing &&
                                    _b3dist(_brep_curve_eval(sec_crv,t2),V)<=
                                    vtol
                                st=t2
                            end
                        end
                    end
                elseif st===nothing &&
                        sec_crv!==nothing && edge_crv!==nothing
                    # unfitted trace on a periodic domain: no level slack to
                    # spend, but a shifted projection may still be nearest —
                    # admit it only on a strict 3-D check
                    V=_brep_curve_eval(edge_crv,te)
                    vtol=2e-3*max(1.0,_b3norm(V))
                    best=nothing;bd=Inf
                    for ku in kulo:kuhi, kv in kvlo:kvhi
                        pr=_brep_pc_project_res(
                            sec_pc,(uv[1]+ku*pu,uv[2]+kv*pv),dlo,dhi)
                        pr===nothing && continue
                        pr[2]<bd && (bd=pr[2];best=pr[1])
                    end
                    if best!==nothing &&
                            _b3dist(_brep_curve_eval(sec_crv,best),V)<=vtol
                        st=best
                    end
                end
            end
        end
        st===nothing || push!(out,st)
    end
    return out
end

# uv bounding-box diagonal of a pcurve over [lo,hi] — the scale a fitted
# trace's fit residual is measured against
function _brep_pc_uv_extent(pc,lo,hi)
    umin=umax=vmin=vmax=nothing
    for t in _brep_pc_grid(pc,lo,hi;subs=2)
        p=_brep_pc_eval(pc,t)
        if umin===nothing
            umin=umax=p[1];vmin=vmax=p[2]
        else
            umin=min(umin,p[1]);umax=max(umax,p[1])
            vmin=min(vmin,p[2]);vmax=max(vmax,p[2])
        end
    end
    umin===nothing && return 0.0
    return hypot(umax-umin,vmax-vmin)
end

# ── section trimming ────────────────────────────────────────────────────────
#
# Restrict a section to the parameter intervals lying inside one face's
# domain: gather crossings of the section pcurve with the face's wires, then
# keep intervals whose midpoints classify inside. Returns (a,b) parameter
# spans in the section curve's natural parameter (circle: angle in [0,2π)
# wraps; line: unbounded → spans end at wire crossings only, so an inside
# line section always terminates at crossings).

# leftmost u of a face's own wire arcs — origin of the canonical cut domain
# for U-periodic surfaces (0 for non-periodic).
function _brep_face_ubase(edges,face)
    _brep_uperiod(face.sg)>0 || return 0.0
    ulo=Inf
    for wire in face.wires,(eidx,sign,pc) in wire
        e=edges[eidx]
        lo,_=_brep_pc_urange(pc,e.t0,e.t1)
        ulo=min(ulo,lo)
    end
    return isfinite(ulo) ? ulo : 0.0
end

# wrap a uv point's u into [ulo, ulo+period) for periodic faces
@inline function _brep_wrap_u(uv,ulo,period)
    period>0 || return uv
    u=ulo+mod(uv[1]-ulo,period)
    return (u,uv[2])
end

# Does the section pcurve run along a domain boundary at `uv`? A span whose
# trace coincides with a wire edge is in-domain even though the parity probe
# is ambiguous on the boundary — the section's 3-D curve then coincides with
# that operand edge and the edge-section conversion claims it downstream.
function _brep_uv_on_wire(edges,face,spc,uv,tol)
    for wire in face.wires,(eidx,_,epc) in wire
        e=edges[eidx]
        e.curve.kind===:degenerate && continue
        _brep_pc_same_carrier(spc,epc,tol) || continue
        elo,ehi=minmax(e.t0,e.t1)
        _brep_pc_param(epc,uv,elo,ehi,tol)===nothing || return true
    end
    return false
end

function _brep_trim_section(solid,face,pc,closed,tol;scrv=nothing)
    uper=_brep_uperiod(face.sg)
    vper=_brep_vperiod(face.sg)
    ulo=_brep_face_ubase(solid.edges,face)
    if uper>0
        # fold the section pcurve into the face's canonical u band so its
        # crossings with the wires compare in the same coordinates — a trace
        # that still lands a whole period off a wire edge is reconciled by
        # the periodic-translate retries inside `_brep_pc_intersect`
        uc=_brep_pc_eval(pc,0.0)[1]
        kshift=uc-(ulo+mod(uc-ulo,uper))
        kshift!=0 && (pc=_brep_shift_pc(pc,kshift))
    end
    params=Float64[]
    for wire in face.wires
        for (eidx,sign,epc) in wire
            e=solid.edges[eidx]
            # degenerate pole/apex edges carry real domain-boundary pcurves
            # (e.g. v=±π/2 on a sphere) — their crossings delimit sections
            append!(params,_brep_pc_intersect(pc,epc,e.t0,e.t1,tol;
                                             sec_crv=scrv,
                                             edge_crv=e.curve,
                                             pu=uper,pv=vper))
        end
    end
    sort!(params)
    uniq=Float64[]
    for p in params
        (isempty(uniq) || p-uniq[end]>tol) && push!(uniq,p)
    end
    out=Tuple{Float64,Float64}[]
    # a section whose uv trace closes on itself (circ2d) wraps its parameter;
    # a lin2d trace (even a latitude circle, closed only by periodicity) is
    # handled as an open curve — the seam cut normalizes the wrap
    if _brep_pc_iscirc(pc)
        isempty(uniq) && begin
            uv0=_brep_wrap_u(_brep_pc_eval(pc,0.0),ulo,uper)
            (_brep_uv_in_face(solid.edges,face,uv0,tol) ||
             _brep_uv_on_wire(solid.edges,face,pc,uv0,tol)) &&
                push!(out,(0.0,_OCC_TWO_PI))
            return out
        end
        u0=mod(uniq[1],_OCC_TWO_PI)
        # crossings arrive in any lift — negated pcurves report mirrored
        # (often negative) params — rebase every crossing into [u0,u0+2π)
        # before adjoining the wrap endpoint; an interval that wrapped the
        # whole circle otherwise gets a midpoint on the wrong arc
        uniq=sort!([u0+mod(p-u0,_OCC_TWO_PI) for p in uniq])
        push!(uniq,u0+_OCC_TWO_PI)
        for i in 1:length(uniq)-1
            a,b=uniq[i],uniq[i+1]
            b-a<tol && continue
            uvm=_brep_wrap_u(_brep_pc_eval(pc,(a+b)/2),ulo,uper)
            (_brep_uv_in_face(solid.edges,face,uvm,tol) ||
             _brep_uv_on_wire(solid.edges,face,pc,uvm,tol)) || continue
            push!(out,(a,b))
        end
        return out
    end
    dlo,dhi=_brep_pc_domain(pc)
    if isfinite(dlo) && isfinite(dhi)
        # bounded traces (spl2d/tan2d): the section exists only on its own
        # knot domain — wire crossings delimit sub-intervals of it, and a
        # crossing-free contour is kept or dropped by a single probe
        cuts=Float64[dlo]
        for p in uniq
            dlo+tol<p<dhi-tol && push!(cuts,p)
        end
        push!(cuts,dhi)
        for i in 1:length(cuts)-1
            a,b=cuts[i],cuts[i+1]
            b-a<tol && continue
            uvm=_brep_wrap_u(_brep_pc_eval(pc,(a+b)/2),ulo,uper)
            (_brep_uv_in_face(solid.edges,face,uvm,tol) ||
             _brep_uv_on_wire(solid.edges,face,pc,uvm,tol)) || continue
            push!(out,(a,b))
        end
        return out
    end
    isempty(uniq) && return out
    n=length(uniq)
    for i in 0:n
        if i==0
            a=-Inf;b=uniq[1];probe=b-max(1.0,uniq[min(n,2)]-uniq[1])
        elseif i==n
            a=uniq[n];b=Inf;probe=a+max(1.0,uniq[n]-uniq[max(1,n-1)])
        else
            a,b=uniq[i],uniq[i+1];probe=(a+b)/2
        end
        uvm=_brep_wrap_u(_brep_pc_eval(pc,probe),ulo,uper)
        (_brep_uv_in_face(solid.edges,face,uvm,tol) ||
         _brep_uv_on_wire(solid.edges,face,pc,uvm,tol)) || continue
        push!(out,(a,b))
    end
    return out
end

# ── edge × face piercing (operand edge splitting) ───────────────────────────

function _brep_edge_face_params(solid,e,face,tol)
    curve=e.curve
    curve.kind===:degenerate && return Float64[]
    sg=face.sg
    params=Float64[]
    if sg.occ in (:tensor_bspline,:tensor_bezier)
        # signed closest-foot distance roots — see `_brep_edge_pierce_tensor`
        append!(params,_brep_edge_pierce_tensor(sg,curve,
                                                minmax(e.t0,e.t1)...,tol))
    elseif curve.kind===:line
        append!(params,_brep_line_surface(sg,curve))
    elseif curve.kind===:circle
        append!(params,_brep_circle_surface(sg,curve,callerless=true))
    elseif curve.kind===:ellipse && sg.occ===:plane
        # same A·cos+B·sin=C normal-distance structure as the circle,
        # axis-weighted by the two radii
        A=curve.rx*_b3dot(curve.x,sg.axis)
        B=curve.ry*_b3dot(curve.y,sg.axis)
        C=_b3dot(_b3(sg.center,curve.c),sg.axis)
        amp=hypot(A,B)
        if amp>=_BREP_TOL && abs(C)<=amp+_BREP_TOL
            φ=atan(B,A);δ=acos(clamp(C/amp,-1.0,1.0))
            for g in (φ-δ,φ+δ), k in -1:1
                push!(params,g+k*_OCC_TWO_PI)
            end
        end
    else
        # generic records — level roots on the curve's own knot grid
        lvl=_brep_surface_level(sg)
        lvl===nothing ||
            append!(params,_brep_scalar_roots(
                t->lvl(_brep_curve_eval(curve,t)),
                _brep_curve_grid(curve,minmax(e.t0,e.t1)...),
                tol;tangents=false))
    end
    out=Float64[]
    lo,hi=minmax(e.t0,e.t1)
    for t in params
        (t<lo-tol || t>hi+tol) && continue
        p=_brep_curve_eval(curve,t)
        uv=_brep_surface_uv(sg,p)
        _brep_uv_in_face(solid.edges,face,uv,tol) || continue
        push!(out,t)
    end
    # a grazing contact reports two roots straddling the touch point; their
    # images coincide in 3-D — collapse them to one cut so the wire doesn't
    # grow a sliver arc around a single vertex
    length(out)<2 && return out
    return _brep_edge_splits(curve,out,e.t0,e.t1,tol)
end

# sort + dedup an edge's split parameters — crossings whose 3-D images
# coincide (grazing contacts, two faces meeting at one point) or that sit
# within the parameter range's own resolution collapse to a single cut; a
# sliver arc around one vertex is not a region
function _brep_edge_splits(curve,params,t0,t1,tol)
    ps=sort!(collect(params))
    ptol=tol*max(1.0,abs(t0),abs(t1))
    out=Float64[]
    for t in ps
        if !isempty(out)
            if t-out[end]<ptol
                continue
            end
            pa=_brep_curve_eval(curve,out[end])
            pb=_brep_curve_eval(curve,t)
            _b3dist(pa,pb)<100tol && continue
        end
        push!(out,t)
    end
    return out
end

# Edge × tensor-patch piercing: the unsigned distance to a free-form surface
# has no scalar level to root, so track the closest-foot point instead —
# g(t)=⟨C(t)−S(u(t),v(t)),n(u(t),v(t))⟩ with (u(t),v(t)) the projection foot.
# A transversal piercing changes g's sign; bisect each sign-change bracket.
# A discontinuous foot jump can mimic a sign change at a large residual, so
# a candidate is kept only when g re-evaluates ~0 at the root. Tangential
# grazes (g touches 0 without crossing) are not piercings.
function _brep_edge_pierce_tensor(sg,crv,lo,hi,tol)
    mscale=1.0
    for c in sg.patch.surface.controls
        mscale=max(mscale,abs(c[1]),abs(c[2]),abs(c[3]))
    end
    ftol=50tol*mscale
    function footlvl(t)
        P=_brep_curve_eval(crv,t)
        u,v,xyz=_model_nurbs_surface_project(sg,P)
        _,du,dv=_occ_surface_d1(sg,u,v)
        n=_b3cross(du,dv)
        nn=_b3norm(n)
        nn<=eps(Float64) && return Inf
        return _b3dot(_b3(P,xyz),n)/nn
    end
    ts=_brep_curve_grid(crv,lo,hi)
    out=Float64[]
    gb=footlvl(ts[1])
    for i in 2:length(ts)
        a,b=ts[i-1],ts[i]
        ga=gb;gb=footlvl(b)
        (isfinite(ga) && isfinite(gb)) || continue
        signbit(ga)==signbit(gb) && continue
        fa=ga
        for _ in 1:60
            m=(a+b)/2
            fm=footlvl(m)
            isfinite(fm) || break
            fm==0.0 && (a=b=m;break)
            signbit(fm)==signbit(fa) ? (a=m;fa=fm) : (b=m)
        end
        r=(a+b)/2
        fr=footlvl(r)
        isfinite(fr) && abs(fr)<=ftol && push!(out,r)
    end
    return out
end

# implicit level of an elementary surface record (~signed distance to the
# underlying surface) — nothing for NURBS/revolution records, whose generic
# edge crossings are not enumerated
function _brep_surface_level(sg)
    if sg.occ===:plane
        return p->_b3dot(_b3(p,sg.center),sg.axis)
    elseif sg.occ===:sphere
        return p->(_b3norm(_b3(p,sg.center))-sg.radius)
    elseif sg.occ===:cylinder || sg.occ===:cone
        return p->begin
            rel=_b3(p,sg.center);oa=_b3dot(rel,sg.axis)
            rp=_b3(rel,_b3mul(sg.axis,oa))
            if sg.occ===:cylinder
                sqrt(_b3dot(rp,rp))-sg.radius
            else
                k=(sg.r2-sg.r1)/sg.height
                sqrt(_b3dot(rp,rp))-(sg.r1+k*oa)
            end
        end
    end
    return nothing
end

# parameter-space sample grid for a 3-D curve record over [lo,hi] —
# B-splines densify at their knot-span breakpoints like the pcurve grid
function _brep_curve_grid(crv,lo,hi;subs::Int=8)
    local bk
    if crv.kind===:bspline && crv.curve isa _ModelOCCNURBSCurve
        bk=unique(crv.curve.curve.knots)
    elseif crv.kind===:bspline && crv.curve isa _ModelOCCReparamCurve
        bk=crv.curve.reparam.knots
    elseif crv.kind===:bspline && crv.curve isa _ModelOCCSectionCurve
        bk=unique(crv.curve.curve.knots)
    else
        bk=nothing
    end
    ts=Float64[]
    if bk!==nothing
        bounds=vcat(Float64[lo],
                    Float64[k for k in bk if lo-1e-12<k<hi+1e-12],
                    Float64[hi])
        sort!(bounds)
        for i in 1:length(bounds)-1
            a,b=bounds[i],bounds[i+1]
            for j in 0:subs
                push!(ts,fma(b-a,j/subs,a))
            end
        end
    else
        for i in 0:64
            push!(ts,fma(hi-lo,i/64,lo))
        end
    end
    return ts
end

# line × surface → curve params (unclamped).
function _brep_line_surface(sg,curve)
    o,d=curve.o,curve.d
    if sg.occ===:plane
        den=_b3dot(d,sg.axis)
        abs(den)<_BREP_TOL && return Float64[]
        return [_b3dot(_b3(sg.center,o),sg.axis)/den]
    elseif sg.occ===:cylinder || sg.occ===:cone
        rel=_b3(o,sg.center)
        da=_b3dot(d,sg.axis);oa=_b3dot(rel,sg.axis)
        dp=_b3(d,_b3mul(sg.axis,da));rp=_b3(rel,_b3mul(sg.axis,oa))
        if sg.occ===:cylinder
            A=_b3dot(dp,dp);B=2*_b3dot(dp,rp);C=_b3dot(rp,rp)-sg.radius^2
        else
            k=(sg.r2-sg.r1)/sg.height
            A=_b3dot(dp,dp)-(k*da)^2
            B=2*(_b3dot(dp,rp)-k*k*oa*da-k*sg.r1*da)
            C=_b3dot(rp,rp)-(sg.r1+k*oa)^2
        end
        abs(A)<_BREP_TOL && return Float64[]
        disc=B*B-4A*C
        disc<0 && return Float64[]
        sq=sqrt(max(0.0,disc))
        return [(-B+sq)/(2A),(-B-sq)/(2A)]
    elseif sg.occ===:sphere
        rel=_b3(o,sg.center)
        A=_b3dot(d,d);B=2*_b3dot(d,rel);C=_b3dot(rel,rel)-sg.radius^2
        disc=B*B-4A*C
        disc<0 && return Float64[]
        sq=sqrt(max(0.0,disc))
        return [(-B+sq)/(2A),(-B-sq)/(2A)]
    end
    return Float64[]
end

# circle × surface → angle params.
function _brep_circle_surface(sg,curve;callerless::Bool=false)
    if sg.occ===:plane
        A=curve.r*_b3dot(curve.x,sg.axis)
        B=curve.r*_b3dot(curve.y,sg.axis)
        C=_b3dot(_b3(sg.center,curve.c),sg.axis)
        amp=sqrt(A*A+B*B)
        amp<_BREP_TOL && return Float64[]
        abs(C)>amp+_BREP_TOL && return Float64[]
        φ=atan(B,A)
        δ=acos(clamp(C/amp,-1.0,1.0))
        return [mod(φ-δ,_OCC_TWO_PI),mod(φ+δ,_OCC_TWO_PI)]
    end
    callerless && return Float64[]
    throw(ArgumentError("_brep_circle_surface: circle×$(sg.occ) unsupported"))
end

# ── uv arcs and the face arrangement ────────────────────────────────────────
#
# An arc is a directed-capable uv segment: (a, b, pc, t0, t1, key, eidx).
# `a`/`b` are node indices with `a` at parameter t0 and `b` at t1. `key`
# identifies the source for 3-D edge deduplication:
#   (:edge, operand, edge_idx, sub) — piece of an operand boundary edge
#   (:sec, operand, sec_idx, sub)   — piece of a section curve
#   (:seam, operand, u_side, sub)   — periodic-domain cut arc
# `eidx` links back to the operand edge record for :edge arcs, or to a
# pending-edge slot for :sec/:seam arcs.

# Node table with tolerance merging.
function _brep_uv_nodes(points,tol)
    nodes=Tuple{Float64,Float64}[]
    idx=Int[]
    for p in points
        found=0
        for (i,q) in enumerate(nodes)
            if abs(p[1]-q[1])<tol && abs(p[2]-q[2])<tol
                found=i;break
            end
        end
        if found==0
            push!(nodes,p);found=length(nodes)
        end
        push!(idx,found)
    end
    return nodes,idx
end

# Departure tangent of arc at a node it leaves from (a-leaving or b-leaving).
function _brep_arc_dep_tangent(ar,from_a::Bool)
    d=_brep_pc_tangent(ar.pc,from_a ? ar.t0 : ar.t1)
    return from_a ? d : (.-d)
end

# Arrival tangent of arc at a node it enters (into a means b→a traversal).
function _brep_arc_arr_tangent(ar,into_b::Bool)
    d=_brep_pc_tangent(ar.pc,into_b ? ar.t1 : ar.t0)
    return into_b ? d : (.-d)
end

# u-extent of a pcurve arc over [t0,t1].
function _brep_pc_urange(pc,t0,t1)
    hasproperty(pc,:neg2d) &&
        return _brep_pc_urange(pc.neg2d[1],-t1,-t0)
    lo,hi=minmax(t0,t1)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return minmax(o[1]+d[1]*t0,o[1]+d[1]*t1)
    elseif hasproperty(pc,:circ2d)
        c,x,r=pc.circ2d
        umin,umax=Inf,-Inf
        for t in (t0,t1)
            u=_brep_pc_eval(pc,t)[1]
            umin=min(umin,u);umax=max(umax,u)
        end
        α=atan(x[2],x[1])
        for (g,s) in ((0.0,r),(π,-r))   # u extremum where t+α≡g
            t=mod(g-α,_OCC_TWO_PI)
            for k in -3:3
                tk=t+k*_OCC_TWO_PI
                lo<tk<hi && (umin=min(umin,c[1]+s);umax=max(umax,c[1]+s))
            end
        end
        return umin,umax
    elseif hasproperty(pc,:ellips2d)
        # u(t)=o1+a·x1·cos−b·x2·sin — extrema where du/dt=0
        o,x,a,b=pc.ellips2d
        umin,umax=Inf,-Inf
        for t in (t0,t1)
            u=_brep_pc_eval(pc,t)[1]
            umin=min(umin,u);umax=max(umax,u)
        end
        φ=atan(-b*x[2],a*x[1])          # u=A·cos+B·sin, max at t=atan(B,A)
        for g in (φ,φ+π)
            for k in -3:3
                tk=g+k*_OCC_TWO_PI
                if lo<tk<hi
                    u=_brep_pc_eval(pc,tk)[1]
                    umin=min(umin,u);umax=max(umax,u)
                end
            end
        end
        return umin,umax
    end
    # generic — du/dt=0 roots on the curve's own knot grid are the exact
    # extrema candidates; endpoints close the bracket
    umin,umax=Inf,-Inf
    for t in (t0,t1)
        u=_brep_pc_eval(pc,t)[1]
        umin=min(umin,u);umax=max(umax,u)
    end
    for t in _brep_scalar_roots(t->_brep_pc_tangent(pc,t)[1],
                                _brep_pc_grid(pc,lo,hi),1e-9;tangents=false)
        u=_brep_pc_eval(pc,t)[1]
        umin=min(umin,u);umax=max(umax,u)
    end
    return umin,umax
end

# v-extent of a pcurve arc over [t0,t1] — the generic mirror of
# `_brep_pc_urange` (endpoints + dv/dt=0 roots on the knot grid).
function _brep_pc_vrange(pc,t0,t1)
    hasproperty(pc,:neg2d) &&
        return _brep_pc_vrange(pc.neg2d[1],-t1,-t0)
    lo,hi=minmax(t0,t1)
    vmin,vmax=Inf,-Inf
    for t in (t0,t1)
        v=_brep_pc_eval(pc,t)[2]
        vmin=min(vmin,v);vmax=max(vmax,v)
    end
    for t in _brep_scalar_roots(t->_brep_pc_tangent(pc,t)[2],
                                _brep_pc_grid(pc,lo,hi),1e-9;tangents=false)
        v=_brep_pc_eval(pc,t)[2]
        vmin=min(vmin,v);vmax=max(vmax,v)
    end
    return vmin,vmax
end

# Split a pcurve segment at u=u_seam0+k·period crossings (periodic cut-domain
# seam). Returns (t_ranges) covering [t0,t1] with no interior seam crossing.
function _brep_pc_split_at_useam(pc,t0,t1,u_seam0,period,tol)
    if hasproperty(pc,:neg2d)
        inner=_brep_pc_split_at_useam(
            pc.neg2d[1],-t1,-t0,u_seam0,period,tol)
        return Tuple{Float64,Float64}[(-b,-a) for (a,b) in reverse(inner)]
    end
    cuts=Float64[]
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        if abs(d[1])>tol
            ua=o[1]+d[1]*t0;ub=o[1]+d[1]*t1
            kmin=ceil(Int,(min(ua,ub)-u_seam0)/period)-1
            kmax=floor(Int,(max(ua,ub)-u_seam0)/period)+1
            for k in kmin:kmax
                t=(u_seam0+k*period-o[1])/d[1]
                (t>t0+tol && t<t1-tol) && push!(cuts,t)
            end
        end
    elseif hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d)
        # u(t)=A·cos+B·sin+C=u_seam0+k·period for both conic kinds
        if hasproperty(pc,:circ2d)
            c,x,r=pc.circ2d
            A=r*x[1];B=-r*x[2];C=c[1]
        else
            o,x,a,b=pc.ellips2d
            A=a*x[1];B=-b*x[2];C=o[1]
        end
        amp=hypot(A,B)
        amp<tol && return [(t0,t1)]
        φ=atan(B,A)
        lo,hi=minmax(t0,t1)
        kmin=floor(Int,(C-amp-u_seam0)/period)-1
        kmax=ceil(Int,(C+amp-u_seam0)/period)+1
        for k in kmin:kmax
            du=(u_seam0+k*period-C)/amp
            abs(du)<=1.0 || continue
            δ=acos(clamp(du,-1.0,1.0))
            for γ in (φ-δ,φ+δ)
                for j in -2:2
                    tk=γ+j*_OCC_TWO_PI
                    (tk>t0+tol && tk<t1-tol) && push!(cuts,tk)
                end
            end
        end
    else
        # generic records — u(t)=u_seam0+k·period roots on the knot grid
        umin,umax=_brep_pc_urange(pc,t0,t1)
        kmin=ceil(Int,(umin-u_seam0)/period)-1
        kmax=floor(Int,(umax-u_seam0)/period)+1
        lo,hi=minmax(t0,t1)
        ts=_brep_pc_grid(pc,lo,hi)
        for k in kmin:kmax
            append!(cuts,_brep_scalar_roots(
                t->_brep_pc_eval(pc,t)[1]-u_seam0-k*period,ts,tol;
                tangents=false))
        end
        filter!(t->t>t0+tol && t<t1-tol,cuts)
    end
    sort!(cuts)
    unique!(cuts)
    out=Tuple{Float64,Float64}[]
    a=t0
    for t in cuts
        push!(out,(a,t));a=t
    end
    push!(out,(a,t1))
    return out
end

# Directed half-edge traversal. `arcs` is a Vector of NamedTuples
# (a,b,pc,t0,t1,key,eidx). Returns cycles: each a Vector of (arc_idx,fwd).
function _brep_uv_arrange(arcs,nnodes::Int,tol,caller)
    n=length(arcs)
    used=falses(n,2)
    adj=[Tuple{Int,Bool}[] for _ in 1:nnodes]
    for (i,ar) in enumerate(arcs)
        push!(adj[ar.a],(i,true))
        push!(adj[ar.b],(i,false))
    end
    cycles=Vector{Tuple{Int,Bool}}[]
    for i0 in 1:n, fwd0 in (true,false)
        used[i0,fwd0 ? 1 : 2] && continue
        cycle=Tuple{Int,Bool}[]
        cur=i0;fwd=fwd0
        closed=false
        steps=0
        while true
            steps+=1
            steps>4n+8 && throw(ErrorException(
                "$caller: face wiring traversal did not terminate"))
            used[cur,fwd ? 1 : 2]=true
            push!(cycle,(cur,fwd))
            ar=arcs[cur]
            node=fwd ? ar.b : ar.a
            tin=_brep_arc_arr_tangent(ar,fwd)   # travel direction INTO node
            rev=(-tin[1],-tin[2])
            angrev=atan(rev[2],rev[1])
            best=0;bestfwd=true;bestang=Inf
            for (j,departs_a) in adj[node]
                jfwd=departs_a
                if j==i0 && jfwd==fwd0
                    # closing candidate allowed
                else
                    used[j,jfwd ? 1 : 2] && continue
                end
                tdep=_brep_arc_dep_tangent(arcs[j],jfwd)
                angdep=atan(tdep[2],tdep[1])
                # smallest clockwise turn from reversed-arrival
                da=mod(angrev-angdep,_OCC_TWO_PI)
                da<tol && (da=_OCC_TWO_PI)  # straight back — last resort
                if da<bestang
                    bestang=da;best=j;bestfwd=jfwd
                end
            end
            if best==0
                break
            elseif best==i0 && bestfwd==fwd0
                closed=true
                break
            end
            cur=best;fwd=bestfwd
        end
        closed || throw(ErrorException(
            "$caller: face wiring traversal hit a dead end"))
        push!(cycles,cycle)
    end
    return cycles
end

# Signed uv area of a cycle (∮u dv − v du)/2 via Green's theorem on arcs.
function _brep_cycle_area(arcs,cycle)
    total=0.0
    for (i,fwd) in cycle
        ar=arcs[i]
        ta,tb=fwd ? (ar.t0,ar.t1) : (ar.t1,ar.t0)
        total+=_brep_pc_arc_area(ar.pc,ta,tb)
    end
    return total
end

# Green's theorem contribution of one pcurve arc over [ta,tb] (signed by
# direction of evaluation): ∫(u·v' − v·u')/2 dt.
function _brep_pc_arc_area(pc,ta,tb)
    # ∫(u v'−v u')dt with u(t)=inner(−t): substituting s=−t yields
    # −∫[−tb,−ta] G(s)ds — the mirrored traversal flips the sign.
    hasproperty(pc,:neg2d) &&
        return -_brep_pc_arc_area(pc.neg2d[1],-tb,-ta)
    if hasproperty(pc,:lin2d)
        u0,v0=_brep_pc_eval(pc,ta)
        u1,v1=_brep_pc_eval(pc,tb)
        return (u0*v1-u1*v0)/2
    elseif hasproperty(pc,:circ2d)
        c,x,r=pc.circ2d
        # u = c1 + r(cos t·x1 − sin t·x2), v = c2 + r(cos t·x2 + sin t·x1)
        # (Y = rot90(X)). u v'−v u' = r² + r·[cos t·(c1x1+c2x2) + sin t·(c2x1−c1x2)]
        # — the concentric term is constant; the off-centre term integrates the
        # frame components against ds = sin tb−sin ta, dc = cos tb−cos ta.
        s0,co0=_gm_sincos(ta);s1,co1=_gm_sincos(tb)
        dc=co1-co0;ds=s1-s0
        return (r*r*(tb-ta) +
                r*(c[1]*(x[1]*ds + x[2]*dc) + c[2]*(x[2]*ds - x[1]*dc)))/2
    elseif hasproperty(pc,:ellips2d)
        # u=o1+A cos+B sin, v=o2+C cos+D sin with A=a·x1,B=−b·x2,C=a·x2,
        # D=b·x1 — the concentric term collapses to a·b(x1²+x2²)=a·b.
        o,x,a,b=pc.ellips2d
        A=a*x[1];B=-b*x[2];C=a*x[2];D=b*x[1]
        s0,co0=_gm_sincos(ta);s1,co1=_gm_sincos(tb)
        dc=co1-co0;ds=s1-s0
        return (a*b*(tb-ta) + o[1]*(C*dc+D*ds) - o[2]*(A*dc+B*ds))/2
    end
    # generic records — 5-point Gauss–Legendre over each knot span
    gn=(-0.9061798459386640,-0.5384693101056831,0.0,
         0.5384693101056831,0.9061798459386640)
    gw=(0.2369268850561891,0.4786286704993665,0.5688888888888889,
        0.4786286704993665,0.2369268850561891)
    sgn=ta<=tb ? 1.0 : -1.0
    breaks=_brep_pc_grid(pc,min(ta,tb),max(ta,tb);subs=1)
    total=0.0
    for i in 1:length(breaks)-1
        x0,x1=breaks[i],breaks[i+1]
        h=(x1-x0)/2
        h<=0 && continue
        m=(x0+x1)/2;s=0.0
        for j in 1:5
            t=m+h*gn[j]
            uv=_brep_pc_eval(pc,t);d=_brep_pc_tangent(pc,t)
            s+=gw[j]*(uv[1]*d[2]-uv[2]*d[1])
        end
        total+=h*s
    end
    return sgn*total/2
end

# Point-in-cycle test (even-odd) for arc lists.
function _brep_cycle_inside(arcs,cycle,uv,tol)
    crossings=0
    for (i,fwd) in cycle
        ar=arcs[i]
        t0,t1=fwd ? (ar.t0,ar.t1) : (ar.t1,ar.t0)
        crossings+=_brep_pc_ray_crossings(ar.pc,t0,t1,uv[1],uv[2],tol)
    end
    return isodd(crossings)
end

# Interior (bounded-side) sample point of a single traversal cycle: midpoint
# of its longest arc nudged left of travel for positive-area (CCW) cycles,
# right for negative ones.
function _brep_cycle_sample(arcs,cycle,area,tol,caller)
    order=sortperm([abs(arcs[i].t1-arcs[i].t0) for (i,f) in cycle],rev=true)
    for oi in order
        i,fwd=cycle[oi]
        ar=arcs[i]
        tm=(ar.t0+ar.t1)/2
        uv=_brep_pc_eval(ar.pc,tm)
        d=_brep_pc_tangent(ar.pc,tm)
        fwd || (d=(.-d))
        l=sqrt(d[1]*d[1]+d[2]*d[2])
        l<eps() && continue
        s=area>0 ? 1.0 : -1.0
        left=(s*-d[2]/l,s*d[1]/l)
        for scale in (1e-3,1e-4,1e-5,1e-6,1e-7)
            cand=(uv[1]+left[1]*scale,uv[2]+left[2]*scale)
            ins=_brep_cycle_inside(arcs,cycle,cand,tol)
            ins && return cand
        end
    end
    throw(ErrorException(
        "$caller: could not sample a Boolean wiring region"))
end

# Group traversal cycles into face regions. Each bounded partition region is
# a positive-area (CCW) cycle; the complement boundary is negative. A negative
# cycle h is a hole of the region that contains the region h encloses: the
# smallest positive cycle d containing h's interior sample, then the smallest
# positive cycle strictly containing d's sample is h's region. Negative cycles
# with no such parent (complement boundaries) are discarded.
function _brep_group_cycles(arcs,cycles,edges,f,tol,caller)
    areas=[_brep_cycle_area(arcs,cy) for cy in cycles]
    samples=Tuple{Float64,Float64}[]
    in_domain=Bool[]
    for (i,a) in enumerate(areas)
        if abs(a)<tol
            throw(ErrorException(
                "$caller: degenerate zero-area wire in Boolean face split"))
        end
        s=_brep_cycle_sample(arcs,cycles[i],a,tol,caller)
        push!(samples,s)
        push!(in_domain,_brep_uv_in_face(edges,f,s,tol))
    end
    pos=[i for i in eachindex(cycles) if areas[i]>0]
    neg=[i for i in eachindex(cycles) if areas[i]<0]
    # smallest positive cycle (by |area|) containing uv, excluding `excl`
    smallest_containing(uv,excl)=begin
        best=0;ba=Inf
        for p in pos
            p==excl && continue
            abs(areas[p])<ba || continue
            _brep_cycle_inside(arcs,cycles[p],uv,tol) || continue
            best=p;ba=abs(areas[p])
        end
        best
    end
    regions=Vector{Int}[]
    for c in pos
        in_domain[c] || continue
        wires=[c]
        for h in neg
            d=smallest_containing(samples[h],0)
            d==0 && continue              # h bounds no region → complement
            parent=smallest_containing(samples[d],d)
            parent==c && push!(wires,h)
        end
        push!(regions,wires)
    end
    return regions,areas
end

# ── solid classification ────────────────────────────────────────────────────
#
# Point-in-solid by ray parity: intersect a ray from P with every face
# surface, keep hits whose uv lands inside the trimmed face, count. Ambiguous
# (grazing/boundary) directions fall through to the next candidate.

# u-coordinate of a face wire's period image center — the mean of its
# arc midpoints' u. Union-merged faces on periodic carriers live on an
# unwrapped chart (u off the canonical [0,2π) window), so containment
# tests must first re-center the probe uv onto the wire's period image.
function _brep_face_wire_ucenter(solid,face)
    s=0.0;n=0
    for wire in face.wires,(eidx,_,pc) in wire
        e=solid.edges[eidx]
        e.curve.kind===:degenerate && continue
        s+=_brep_pc_eval(pc,(e.t0+e.t1)/2)[1];n+=1
    end
    n==0 && return 0.0
    return s/n
end

function _brep_ray_face_hits(solid,face,P,D,tol)
    sg=face.sg
    params=_brep_line_surface(sg,(o=P,d=D))
    count=0;ambiguous=false
    for t in params
        t<=tol && continue
        p=_b3add(P,_b3mul(D,t))
        uv=_brep_surface_uv(sg,p)
        uper=_brep_uperiod(sg)
        if uper>0
            cu=_brep_face_wire_ucenter(solid,face)
            uv=(uv[1]+uper*round((cu-uv[1])/uper),uv[2])
        end
        # inside the trimmed face?
        inside=_brep_uv_in_face(solid.edges,face,uv,tol)
        # near a wire boundary → ambiguous direction
        near=false
        for wire in face.wires,(eidx,sign,pc) in wire
            e=solid.edges[eidx]
            e.curve.kind===:degenerate && continue
            # coarse proximity: evaluate edge midpoint? use cheap endpoint
            # distance + interval projection for lin2d only
            if _brep_pc_near(pc,e.t0,e.t1,uv,20tol)
                near=true;break
            end
        end
        near && (ambiguous=true)
        inside && (count+=1)
    end
    return count,ambiguous
end

# Cheap uv proximity of a point to a pcurve arc (endpoint + line distance).
function _brep_pc_near(pc,t0,t1,uv,tol)
    hasproperty(pc,:neg2d) &&
        return _brep_pc_near(pc.neg2d[1],-t1,-t0,uv,tol)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        n2=d[1]*d[1]+d[2]*d[2]
        n2<tol*tol && begin
            return abs(uv[1]-o[1])<tol && abs(uv[2]-o[2])<tol
        end
        s=((uv[1]-o[1])*d[1]+(uv[2]-o[2])*d[2])/n2
        s=clamp(s,min(t0,t1),max(t0,t1))
        dx=uv[1]-(o[1]+s*d[1]);dy=uv[2]-(o[2]+s*d[2])
        return dx*dx+dy*dy<tol*tol
    elseif hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d)
        # coarse radial band (absolute — the level is signed inside the
        # conic's disk), then the exact angular range check
        lvl=_brep_pc_level(pc)
        abs(lvl(uv))<tol || return false
        return _brep_pc_uv_to_t(pc,uv,tol,min(t0,t1),max(t0,t1))!==nothing
    end
    # generic records — seeded Newton projection gives the closest point;
    # the caller's `tol` is the near band itself (typically ~20·tol)
    s=_brep_pc_project(pc,uv,min(t0,t1),max(t0,t1),tol)
    s===nothing && return false
    p=_brep_pc_eval(pc,s)
    return hypot(p[1]-uv[1],p[2]-uv[2])<tol
end

const _BREP_RAY_DIRS = (
    (0.890120899537191, 0.342847862075555, 0.304948994628101),
    (0.334354466399804, 0.887197286487756, 0.318721580386145),
    (0.312217584669563, 0.301256145227518, 0.900553430549812),
    (0.577350269189626, 0.577350269189626, 0.577350269189626))

function _brep_point_in_solid(solid,P,tol,caller)
    for D in _BREP_RAY_DIRS
        count=0;ambiguous=false
        for face in solid.faces
            c,amb=_brep_ray_face_hits(solid,face,P,D,tol)
            count+=c;ambiguous|=amb
        end
        ambiguous && continue
        return isodd(count)
    end
    throw(ErrorException(
        "$caller: solid containment query is degenerate on all ray directions"))
end
# ── per-face arrangement ────────────────────────────────────────────────────

# shift a pcurve record's u-coordinate by -shift
function _brep_shift_pc(pc,shift)
    hasproperty(pc,:neg2d) &&
        return (neg2d=(_brep_shift_pc(pc.neg2d[1],shift),),)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=((o[1]-shift,o[2]),d),)
    elseif hasproperty(pc,:circ2d)
        c,x,r=pc.circ2d
        return (circ2d=((c[1]-shift,c[2]),x,r),)
    elseif hasproperty(pc,:ellips2d)
        o,x,a,b=pc.ellips2d
        return (ellips2d=((o[1]-shift,o[2]),x,a,b),)
    elseif hasproperty(pc,:spl2d)
        knots,p,cps=pc.spl2d
        out=(spl2d=(knots,p,
            NTuple{2,Float64}[(c[1]-shift,c[2]) for c in cps]),)
        if hasproperty(pc,:lvl2d)
            l=pc.lvl2d
            out=merge(out,(lvl2d=(uv->l((uv[1]+shift,uv[2]))),))
        end
        return out
    elseif hasproperty(pc,:tan2d)
        inner,u0,us,vv=pc.tan2d
        return (tan2d=(inner,u0-shift,us,vv),)
    elseif hasproperty(pc,:tan2dv)
        inner,uu,v0,vs=pc.tan2dv
        return (tan2dv=(inner,uu-shift,v0,vs),)
    end
    throw(ArgumentError("_brep_shift_pc: unsupported pcurve record"))
end

# every parameter in (t0,t1) where the pcurve's u crosses a cut-domain seam
# (u = ulo + k·period) for U-periodic surfaces
function _brep_pc_split_all_useams(pc,t0,t1,period,ulo,tol)
    if hasproperty(pc,:neg2d)
        inner=_brep_pc_split_all_useams(
            pc.neg2d[1],-t1,-t0,period,ulo,tol)
        return Tuple{Float64,Float64}[(-b,-a) for (a,b) in reverse(inner)]
    end
    cuts=Float64[]
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        if abs(d[1])>tol
            ua=o[1]+d[1]*t0; ub=o[1]+d[1]*t1
            kmin=ceil(Int,(min(ua,ub)-ulo)/period)-1
            kmax=floor(Int,(max(ua,ub)-ulo)/period)+1
            for k in kmin:kmax
                t=(ulo+k*period-o[1])/d[1]
                (t>t0+tol && t<t1-tol) && push!(cuts,t)
            end
        end
    elseif hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d)
        if hasproperty(pc,:circ2d)
            c,x,r=pc.circ2d
            A=r*x[1];B=-r*x[2];C=c[1]
        else
            o,x,a,b=pc.ellips2d
            A=a*x[1];B=-b*x[2];C=o[1]
        end
        amp=hypot(A,B)
        amp<tol && return [(t0,t1)]
        φ=atan(B,A)
        kmin=floor(Int,(C-amp-ulo)/period)-1
        kmax=ceil(Int,(C+amp-ulo)/period)+1
        for k in kmin:kmax
            v=(ulo+k*period-C)/amp
            abs(v)<=1.0 || continue
            δ=acos(clamp(v,-1.0,1.0))
            for γ in (φ-δ,φ+δ)
                for j in -2:2
                    tk=γ+j*_OCC_TWO_PI
                    (tk>t0+tol && tk<t1-tol) && push!(cuts,tk)
                end
            end
        end
    else
        # generic records — seam crossings found on the curve's knot grid
        umin,umax=_brep_pc_urange(pc,t0,t1)
        kmin=ceil(Int,(umin-ulo)/period)-1
        kmax=floor(Int,(umax-ulo)/period)+1
        lo,hi=minmax(t0,t1)
        ts=_brep_pc_grid(pc,lo,hi)
        for k in kmin:kmax
            append!(cuts,_brep_scalar_roots(
                t->_brep_pc_eval(pc,t)[1]-ulo-k*period,ts,tol;
                tangents=false))
        end
        filter!(t->t>t0+tol && t<t1-tol,cuts)
    end
    sort!(cuts)
    unique!(cuts)
    out=Tuple{Float64,Float64}[]
    a=t0
    for t in cuts
        push!(out,(a,t));a=t
    end
    push!(out,(a,t1))
    return out
end

# split [t0,t1] at strictly-interior params (sorted, deduplicated)
function _brep_split_ranges(t0,t1,params,tol)
    uniq=Float64[]
    for p in sort!(collect(params))
        (p>t0+tol && p<t1-tol) || continue
        (isempty(uniq) || p-uniq[end]>tol) && push!(uniq,p)
    end
    out=Tuple{Float64,Float64}[]
    a=t0
    for p in uniq
        push!(out,(a,p));a=p
    end
    push!(out,(a,t1))
    return out
end

# invert a pcurve: parameter t∈[t0,t1] with eval(pc,t)≈uv, or nothing
function _brep_pc_param(pc,uv,t0,t1,tol)
    if hasproperty(pc,:neg2d)
        s=_brep_pc_param(pc.neg2d[1],uv,-t1,-t0,tol)
        return s===nothing ? nothing : -s
    end
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        n2=d[1]*d[1]+d[2]*d[2]
        n2<tol*tol && return nothing
        s=((uv[1]-o[1])*d[1]+(uv[2]-o[2])*d[2])/n2
        e=_brep_pc_eval(pc,s)
        (abs(e[1]-uv[1])<10tol && abs(e[2]-uv[2])<10tol) || return nothing
        (s>=t0-tol && s<=t1+tol) || return nothing
        return s
    elseif hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d)
        lvl=_brep_pc_level(pc)
        rs=hasproperty(pc,:circ2d) ? pc.circ2d[3] :
           min(pc.ellips2d[3],pc.ellips2d[4])
        # absolute band — the conic level is signed (negative inside)
        abs(lvl(uv))>10tol*max(1.0,rs) && return nothing
        return _brep_pc_uv_to_t(pc,uv,tol,t0,t1)
    end
    # generic records — grid-seeded projection restricted to [t0,t1]
    return _brep_pc_project(pc,uv,min(t0,t1),max(t0,t1),tol)
end

# Build the face's uv arrangement (boundary sub-edges + section pieces, cut at
# the periodic seam) and extract kept-region wires by leftmost-turn traversal.
# Returns (regions, arcs, cycles) — each region is a Vector of cycles (outer
# first); a cycle is a Vector of (arc_index, forward).
function _brep_face_regions(solid,f,side,edge_splits,secs,extras,tol,caller)
    uper=_brep_uperiod(f.sg)
    if get(ENV,"BREP_DEBUG_ARCS","")!=""
        println("DBG face sg=$(f.sg.occ) stag=$(get(f,:stag,0)) secs=",
                [(se=sf.se,span=sf.span) for sf in secs],
                " esplits=",Dict(k=>sort(round.(v,digits=3)) for (k,v) in edge_splits))
    end
    raw=NamedTuple[]
    # a section arc coincident with an earlier section arc is the same
    # geometric edge (e.g. an operand edge lying in this face's plane gives
    # the section of both adjacent faces) — keep the first; a wire piece
    # coincident with a section drops in favour of the section's shared
    # (:s,se) key so both faces reference one result edge
    kept_secs=NamedTuple[]
    for sf in secs
        uvm=_brep_pc_eval(sf.pc,(sf.span[1]+sf.span[2])/2)
        dup=any(ks->_brep_pc_same_carrier(sf.pc,ks.pc,tol) &&
                _brep_pc_param(ks.pc,uvm,ks.span[1],ks.span[2],tol)!==nothing,
                kept_secs)
        dup || push!(kept_secs,sf)
    end
    for wire in f.wires,(eidx,sign,pc) in wire
        e=solid.edges[eidx]
        # every writer upstream (pierce roots, section endpoints, imprint
        # contacts) feeds one vertex set — collapse params whose 3-D images
        # coincide so a grazing or doubly-reported contact can't grow a
        # sliver arc around a single vertex
        extra=_brep_edge_splits(e.curve,get(edge_splits,eidx,Float64[]),
                                e.t0,e.t1,tol)
        for (k,(ta,tb)) in enumerate(_brep_split_ranges(e.t0,e.t1,extra,tol))
            uvm=_brep_pc_eval(pc,(ta+tb)/2)
            dup=any(ks->_brep_pc_same_carrier(pc,ks.pc,tol) &&
                    _brep_pc_param(ks.pc,uvm,ks.span[1],ks.span[2],tol)!==
                    nothing, kept_secs)
            # a degenerate-section imprint carries the shared edge's own key
            # — the coincident wire piece yields to it
            if !dup
                dup=any(ex->get(ex,:esec,false) &&
                        _brep_pc_same_carrier(pc,ex.pc,tol) &&
                        _brep_pc_param(ex.pc,uvm,ex.t0,ex.t1,tol)!==
                        nothing, extras)
            end
            dup && continue
            # the canonical edge identity is the model curve tag plus the
            # edge-parameter span — two operand solids sharing the same
            # Curve (e.g. republished cells on a common wall) resolve to
            # one result edge; pieces split differently per operand keep
            # distinct keys rather than union into a broken loop
            push!(raw,(pc=pc,t0=ta,t1=tb,
                       key=(:e,solid.ctags[eidx],ta,tb),crv=e.curve))
        end
    end
    for sf in kept_secs
        push!(raw,(pc=sf.pc,t0=sf.span[1],t1=sf.span[2],
                   key=(:s,sf.se),crv=sf.crv,
                   fix0=get(sf,:fix0,nothing),fix1=get(sf,:fix1,nothing)))
    end
    append!(raw,extras)
    # canonical cut domain: [ulo, ulo+uper], ulo = leftmost u of the face's
    # own wires. Stray section pieces fold into it; arcs ON the right seam
    # (u≈ulo+uper) stay — the two seam wedges of a periodic edge are distinct
    # arcs and must not collapse onto each other.
    ulo=0.0
    if uper>0
        ulo=Inf
        for ar in raw
            ar.key[1]===:e || continue
            lo,_=_brep_pc_urange(ar.pc,ar.t0,ar.t1)
            ulo=min(ulo,lo)
        end
        isfinite(ulo) || (ulo=0.0)
    end
    arcs=NamedTuple[]
    for ar in raw
        # `neg` arcs are imprinted pcurves evaluated at s = −t_edge — the
        # result edge's own parameter range is recovered as (−t1,−t0)
        neg=get(ar,:neg,false)
        if uper>0
            for (k,(ta,tb)) in enumerate(
                    _brep_pc_split_all_useams(ar.pc,ar.t0,ar.t1,uper,ulo,tol))
                umid=_brep_pc_eval(ar.pc,(ta+tb)/2)[1]
                band=floor(Int,(umid-ulo)/uper)
                band>0 && umid-ulo-band*uper<tol && (band-=1)
                pc2=band==0 ? ar.pc : _brep_shift_pc(ar.pc,band*uper)
                # section arcs keep (:s,se) — spans were already split at
                # both faces' seams, so each se maps to one result edge on
                # either side; wire arcs carry their own piece key
                va,vb=neg ? (-tb,-ta) : (ta,tb)
                fx0=ta==ar.t0 ? get(ar,:fix0,nothing) : nothing
                fx1=tb==ar.t1 ? get(ar,:fix1,nothing) : nothing
                # anchor fixes are expressed in the face's canonical u
                # chart while a fitted trace's pcurve can sit whole
                # periods off (unwrapped contour, right-seam wedge):
                # lift each fix onto the band its own piece's endpoint
                # eval occupies so the shared vertex stays bitwise
                fx0!==nothing &&
                    (fx0=(fx0[1]+round(Int,
                          (_brep_pc_eval(pc2,ta)[1]-fx0[1])/uper)*uper,
                          fx0[2]))
                fx1!==nothing &&
                    (fx1=(fx1[1]+round(Int,
                          (_brep_pc_eval(pc2,tb)[1]-fx1[1])/uper)*uper,
                          fx1[2]))
                push!(arcs,(pc=pc2,t0=ta,t1=tb,key=ar.key,crv=ar.crv,
                            vt0=va,vt1=vb,neg=neg,a=0,b=0,
                            fix0=fx0,fix1=fx1))
            end
        else
            va,vb=neg ? (-ar.t1,-ar.t0) : (ar.t0,ar.t1)
            push!(arcs,(pc=ar.pc,t0=ar.t0,t1=ar.t1,key=ar.key,crv=ar.crv,
                        vt0=va,vt1=vb,neg=neg,a=0,b=0,
                        fix0=get(ar,:fix0,nothing),
                        fix1=get(ar,:fix1,nothing)))
        end
    end
    isempty(arcs) && return NamedTuple[],arcs,Vector{Tuple{Int,Bool}}[]
    pts=Tuple{Float64,Float64}[]
    for ar in arcs
        # anchored section ends place the shared vertex exactly — the
        # fitted pcurve's own endpoint can sit a residual off it
        push!(pts,ar.fix0===nothing ? _brep_pc_eval(ar.pc,ar.t0) : ar.fix0)
        push!(pts,ar.fix1===nothing ? _brep_pc_eval(ar.pc,ar.t1) : ar.fix1)
    end
    nodes,idx=_brep_uv_nodes(pts,tol)
    arcs=[merge(ar,(a=idx[2i-1],b=idx[2i])) for (i,ar) in enumerate(arcs)]
    # coincident duplicate arcs — e.g. two operands' edges imprinted onto the
    # same span of this face — only pair into zero-area lens cycles; keep the
    # first (face-wire arcs precede section/imprint extras)
    keep_arc=trues(length(arcs))
    for i in eachindex(arcs), j in i+1:length(arcs)
        keep_arc[i] && keep_arc[j] || continue
        ai,aj=arcs[i],arcs[j]
        (ai.a==aj.a && ai.b==aj.b || ai.a==aj.b && ai.b==aj.a) || continue
        mi=_brep_pc_eval(ai.pc,(ai.t0+ai.t1)/2)
        mj=_brep_pc_eval(aj.pc,(aj.t0+aj.t1)/2)
        hypot(mi[1]-mj[1],mi[2]-mj[2])<tol && (keep_arc[j]=false)
    end
    arcs=arcs[keep_arc]
    # collapsed self-loops: an arc whose endpoints merged to one node only
    # bounds a region when the loop has real extent (a closed conic edge).
    # A sliver piece whose whole loop fits inside tol — e.g. a grazing
    # contact split the wire twice around one vertex — carries no region;
    # drop it before traversal pairs it into a zero-area lens cycle
    keep_arc=trues(length(arcs))
    for (i,ar) in enumerate(arcs)
        ar.a==ar.b || continue
        p0=_brep_pc_eval(ar.pc,ar.t0)
        dev=0.0
        for t in _brep_pc_grid(ar.pc,ar.t0,ar.t1;subs=2)
            q=_brep_pc_eval(ar.pc,t)
            dev=max(dev,hypot(q[1]-p0[1],q[2]-p0[2]))
            dev>=tol && break
        end
        dev<tol && (keep_arc[i]=false)
    end
    arcs=arcs[keep_arc]
    # dangling chords: a section piece may end strictly inside this face
    # where the OTHER operand's boundary pierces it (e.g. a bounded sheet
    # stabbed through a solid's wall). Such pieces split no region — prune
    # them iteratively like graph leaves so wiring traversal never pairs a
    # chord with itself into a zero-area lens.
    while true
        deg=zeros(Int,length(nodes))
        for ar in arcs
            deg[ar.a]+=1;deg[ar.b]+=1
        end
        drop=findall(i->deg[arcs[i].a]<2 || deg[arcs[i].b]<2,
                     eachindex(arcs))
        isempty(drop) && break
        if get(ENV,"BREP_DEBUG_ARCS","")!=""
            for i in drop
                ar=arcs[i]
                println("DBG   DROP key=$(ar.key) t=[$(ar.t0),$(ar.t1)] deg" *
                        " $(deg[ar.a]),$(deg[ar.b])")
            end
        end
        deleteat!(arcs,drop)
    end
    isempty(arcs) && return NamedTuple[],arcs,Vector{Tuple{Int,Bool}}[]
    if get(ENV,"BREP_DEBUG_ARCS","")!=""
        println("DBG regions face sg=$(f.sg.occ) stag=$(get(f,:stag,0)) arcs:")
        for (i,ar) in enumerate(arcs)
            u0=_brep_pc_eval(ar.pc,ar.t0);u1=_brep_pc_eval(ar.pc,ar.t1)
            println("   arc$i key=$(ar.key) t=[$(ar.t0),$(ar.t1)] " *
                    "uv $(round.(u0,digits=3))->$(round.(u1,digits=3)) " *
                    "nodes $(ar.a)->$(ar.b)")
        end
    end
    cycles=_brep_uv_arrange(arcs,length(nodes),tol,caller)
    regions,_=_brep_group_cycles(arcs,cycles,solid.edges,f,tol,caller)
    if get(ENV,"BREP_DEBUG_ARCS","")!=""
        println("DBG   cycles=$(cycles) regions=$regions")
        for (ci,cy) in enumerate(cycles)
            println("DBG   cycle$ci area=",
                    round(_brep_cycle_area(arcs,cy),sigdigits=6))
        end
    end
    return regions,arcs,cycles
end

# A uv point strictly inside the region's outer cycle and outside its holes:
# midpoint of the longest outer arc, nudged toward the interior (left of
# travel for a CCW outer cycle).
function _brep_region_sample(arcs,cycles,region,tol,caller)
    outer=cycles[region[1]]
    # domain scale from arc endpoints
    umin=vmin=Inf;umax=vmax=-Inf
    for cy in region,(i,fwd) in cycles[cy]
        ar=arcs[i]
        for t in (ar.t0,ar.t1)
            uv=_brep_pc_eval(ar.pc,t)
            umin=min(umin,uv[1]);vmin=min(vmin,uv[2])
            umax=max(umax,uv[1]);vmax=max(vmax,uv[2])
        end
    end
    ds=max(1e-3,umax-umin,vmax-vmin)
    order=sortperm([abs(arcs[i].t1-arcs[i].t0) for (i,fwd) in outer],rev=true)
    for oi in order
        i,fwd=outer[oi]
        ar=arcs[i]
        tm=(ar.t0+ar.t1)/2
        uv=_brep_pc_eval(ar.pc,tm)
        d=_brep_pc_tangent(ar.pc,tm)
        fwd || (d=(.-d))
        l=sqrt(d[1]*d[1]+d[2]*d[2])
        l<eps() && continue
        left=(-d[2]/l,d[1]/l)
        for scale in (1e-3,1e-4,1e-5,1e-6,1e-7)
            cand=(uv[1]+left[1]*scale*ds,uv[2]+left[2]*scale*ds)
            _brep_cycle_inside(arcs,outer,cand,tol) || continue
            ok=true
            for c in region[2:end]
                _brep_cycle_inside(arcs,cycles[c],cand,tol) && (ok=false;break)
            end
            ok && return cand
        end
    end
    throw(ErrorException(
        "$caller: could not place a sample point inside a Boolean face region"))
end

# Intersect two kept-span sets on a section curve. Open curves use plain
# interval intersection; closed (circ2d-trace) curves decompose the period
# [0,2π) into cells at all span endpoints so differently-lifted spans
# intersect correctly. Returns merged (lo,hi) ranges.
# natural parameter period of a closed section pcurve: conic records are
# natively 2π-periodic; a bounded fitted trace closes over its own knot
# domain. Returns 0 for unbounded/non-periodic parametrizations.
function _brep_pc_param_period(pc)
    hasproperty(pc,:neg2d) && return _brep_pc_param_period(pc.neg2d[1])
    (hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d)) &&
        return _OCC_TWO_PI
    dlo,dhi=_brep_pc_domain(pc)
    (isfinite(dlo)&&isfinite(dhi)) || return 0.0
    return dhi-dlo
end

# whether the pcurve evaluator itself wraps parameters — conic records
# evaluate sin/cos of t so an out-of-domain t stays a valid trace point;
# bounded spline/reparam traces do not (their wrap point must become an
# ordinary span boundary instead).
_brep_pc_native_periodic(pc) =
    hasproperty(pc,:neg2d) ? _brep_pc_native_periodic(pc.neg2d[1]) :
    (hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d))

function _brep_intersect_spans(spansA,spansB,closed,tol,caller;pc=nothing)
    if !closed
        ovl=Tuple{Float64,Float64}[]
        for (a0,a1) in spansA,(b0,b1) in spansB
            lo=max(a0,b0);hi=min(a1,b1)
            lo<hi-tol || continue
            push!(ovl,(lo,hi))
        end
    else
        # the closed parameter circle's period: conic traces wrap mod 2π,
        # bounded fitted traces close over their own domain — and since
        # their evaluator does not wrap, spans must stay inside it
        per=_OCC_TWO_PI;base=0.0;native=true
        if pc!==nothing
            per=_brep_pc_param_period(pc)
            dlo,_=_brep_pc_domain(pc)
            isfinite(dlo) && (base=dlo)
            native=_brep_pc_native_periodic(pc)
        end
        per>0 || (per=_OCC_TWO_PI)
        cuts=Float64[base]
        for s in Iterators.flatten((spansA,spansB)), t in s
            isfinite(t) || continue
            c=base+mod(t-base,per)
            any(q->abs(q-c)<tol,cuts) || push!(cuts,c)
        end
        sort!(cuts)
        incell(spans,m)=any(
            ((a,b),)->begin
                for k in -1:1
                    mk=m+k*per
                    (mk>a+tol && mk<b-tol) && return true
                end
                false
            end, spans)
        ovl=Tuple{Float64,Float64}[]
        n=length(cuts)
        for i in 1:n
            a=cuts[i];b=i<n ? cuts[i+1] : cuts[1]+per
            m=(a+b)/2
            (incell(spansA,m) && incell(spansB,m)) || continue
            push!(ovl,(a,b))
        end
        # merge across the wrap point when both end cells are kept — only
        # natively periodic traces may carry the wrapped (negative-t) span;
        # a bounded trace keeps two in-domain pieces welded at the wrap node
        while native && length(ovl)>1 &&
                ovl[end][2]>=ovl[1][1]+per-tol
            ovl[1]=(ovl[end][1]-per,ovl[1][2])
            pop!(ovl)
        end
    end
    merged=Tuple{Float64,Float64}[]
    for (lo,hi) in sort!(ovl)
        (isfinite(lo) && isfinite(hi)) || throw(ErrorException(
            "$caller: an unbounded section piece survived trimming"))
        if !isempty(merged) && lo<=merged[end][2]+tol
            merged[end]=(merged[end][1],max(merged[end][2],hi))
        else
            push!(merged,(lo,hi))
        end
    end
    return merged
end

# ── coincident (same-domain) plane faces ────────────────────────────────────
#
# Two faces sharing a carrier plane produce no section curve — their overlap
# is resolved by imprinting: each face's boundary wires are mapped into the
# other's uv domain and join the arrangement as ordinary arcs, so kept
# regions split along the coincident partner's outline. Every crossing and
# collinear-overlap endpoint becomes a global operand-edge split or a section
# cut so the implied result vertices exist on every face that needs them.

_brep_planes_coincident(pa,pb) =
    _b3norm(_b3cross(pa.axis,pb.axis))<_BREP_TOL &&
    abs(_b3dot(_b3(pb.center,pa.center),pa.axis))<_BREP_TOL

# affine map src-uv → dst-uv between coincident planes (exact: the columns
# are the images of src.X/src.Y in dst's frame)
function _brep_plane_map(dst,src)
    o=_brep_plane_uv(dst,_brep_surface_eval(src,0.0,0.0))
    u=_brep_plane_uv(dst,_brep_surface_eval(src,1.0,0.0))
    v=_brep_plane_uv(dst,_brep_surface_eval(src,0.0,1.0))
    return (c1=(u[1]-o[1],u[2]-o[2]),c2=(v[1]-o[1],v[2]-o[2]),b=o)
end

@inline _brep_map_uv(mp,uv) =
    (fma(uv[1],mp.c1[1],fma(uv[2],mp.c2[1],mp.b[1])),
     fma(uv[1],mp.c1[2],fma(uv[2],mp.c2[2],mp.b[2])))

@inline _brep_map_dir(mp,d) =
    (mp.c1[1]*d[1]+mp.c2[1]*d[2],mp.c1[2]*d[1]+mp.c2[2]*d[2])

# pcurve under an orthogonal affine map. A reflection (det M<0) cannot be a
# `circ2d` evaluated at the same parameter — the returned pcurve evaluated at
# -t gives the mapped point and `neg` reports the reversal (lines are
# unaffected by orientation).
function _brep_map_pc(pc,mp)
    inner,wrapped=_brep_pc_unneg(pc)
    if hasproperty(inner,:lin2d)
        o,d=inner.lin2d
        return (lin2d=(_brep_map_uv(mp,o),_brep_map_dir(mp,d)),),false
    elseif hasproperty(inner,:circ2d) || hasproperty(inner,:ellips2d)
        neg=(mp.c1[1]*mp.c2[2]-mp.c1[2]*mp.c2[1])<0
        mapped=if hasproperty(inner,:circ2d)
            c,x,r=inner.circ2d
            (circ2d=(_brep_map_uv(mp,c),_brep_map_dir(mp,x),r),)
        else
            o,x,a,b=inner.ellips2d
            (ellips2d=(_brep_map_uv(mp,o),_brep_map_dir(mp,x),a,b),)
        end
        # a reflection inside the map composes with an existing neg2d
        # wrapper: mapped(inner(−t)) needs −t evaluated on the mapped conic.
        return (wrapped ? (neg2d=(mapped,),) : mapped),neg
    elseif hasproperty(inner,:spl2d)
        # affine pole map is exact for B-splines and carries no handedness
        knots,p,cps=inner.spl2d
        mapped=(spl2d=(knots,p,
            NTuple{2,Float64}[_brep_map_uv(mp,c) for c in cps]),)
        return (wrapped ? (neg2d=(mapped,),) : mapped),false
    end
    throw(ArgumentError(
        "_brep_map_pc: $(first(propertynames(inner))) pcurve has no " *
        "plane-map representable image"))
end

# pcurve with its parameter reversed (t → −t): used to store an imprinted
# arc's pcurve under the operand edge's own parametrization. A reversed
# line negates its direction; a reversed circle's left-handed trace carries
# the `neg2d` wrapper.
function _brep_flip_pc(pc,caller)
    hasproperty(pc,:neg2d) && return pc.neg2d[1]
    hasproperty(pc,:lin2d) &&
        return (lin2d=(pc.lin2d[1],(-pc.lin2d[2][1],-pc.lin2d[2][2])),)
    (hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d)) &&
        return (neg2d=(pc,),)
    if hasproperty(pc,:spl2d)
        # exact parameter reversal: p(−t) is the same spline over mirrored
        # knots with mirrored poles
        knots,p,cps=pc.spl2d
        out=(spl2d=(Float64[-k for k in Iterators.reverse(knots)],p,
                    NTuple{2,Float64}[cp for cp in
                                      Iterators.reverse(cps)]),)
        hasproperty(pc,:lvl2d) && (out=merge(out,(lvl2d=pc.lvl2d,)))
        return out
    end
    throw(ArgumentError(
        "$caller: unrepresentable parameter-reversed pcurve"))
end

# same carrier test: collinear lines or coincident conics — generic records
# decide by sampling one operand's trace and projecting onto the other
function _brep_pc_same_carrier(a,b,tol)
    a,_=_brep_pc_unneg(a);b,_=_brep_pc_unneg(b)
    if hasproperty(a,:lin2d) && hasproperty(b,:lin2d)
        (o1,d1),(o2,d2)=a.lin2d,b.lin2d
        n1=d1[1]*d1[1]+d1[2]*d1[2];n2=d2[1]*d2[1]+d2[2]*d2[2]
        (n1<tol*tol || n2<tol*tol) && return false
        abs(d1[1]*d2[2]-d1[2]*d2[1])>100tol*sqrt(n1*n2) && return false
        w1=o2[1]-o1[1];w2=o2[2]-o1[2]
        return abs(w1*d1[2]-w2*d1[1])<100tol*sqrt(n1)
    elseif hasproperty(a,:circ2d) && hasproperty(b,:circ2d)
        c1,x1,r1=a.circ2d;c2,x2,r2=b.circ2d
        return abs(c1[1]-c2[1])<100tol && abs(c1[2]-c2[2])<100tol &&
               abs(r1-r2)<100tol
    elseif hasproperty(a,:ellips2d) && hasproperty(b,:ellips2d)
        o1,x1,a1,b1=a.ellips2d;o2,x2,a2,b2=b.ellips2d
        (abs(o1[1]-o2[1])<100tol && abs(o1[2]-o2[2])<100tol &&
         abs(a1-a2)<100tol && abs(b1-b2)<100tol) || return false
        # a circle's x direction is arbitrary; a true ellipse needs the
        # first semi-axis directions ±parallel
        abs(a1-b1)<100tol &&
            return true
        return abs(x1[1]*x2[1]+x1[2]*x2[2])>1-100tol
    end
    # mismatched kinds or splines: seven samples of a's locus must each lie
    # on b — periodic records sample one period, splines their knot domain
    dlo,dhi=_brep_pc_domain(a)
    if isinf(dlo)
        dlo,dhi=0.0,_OCC_TWO_PI
    end
    blo,bhi=_brep_pc_domain(b)
    if isinf(blo)
        blo,bhi=-4*_OCC_TWO_PI,4*_OCC_TWO_PI
    end
    for i in 0:6
        t=fma(dhi-dlo,i/6,dlo)
        uv=_brep_pc_eval(a,t)
        _brep_pc_project(b,uv,blo,bhi,tol)===nothing && return false
    end
    return true
end

# same 3-D carrier test for result curves (collinear lines / same circle)
function _brep_same_curve3d(c1,c2,tol)
    c1.kind!==c2.kind && return false
    if c1.kind===:line
        _b3norm(_b3cross(c1.d,c2.d))>100tol && return false
        return _b3norm(_b3cross(_b3(c2.o,c1.o),c1.d))<100tol
    elseif c1.kind===:circle
        return _b3dist(c1.c,c2.c)<100tol && abs(c1.r-c2.r)<100tol &&
               _b3norm(_b3cross(c1.n,c2.n))<100tol
    elseif c1.kind===:ellipse
        _b3dist(c1.c,c2.c)<100tol || return false
        _b3norm(_b3cross(c1.n,c2.n))<100tol || return false
        # the conic metric xxᵀ/rx²+yyᵀ/ry² is frame-rotation invariant
        for i in 1:3, j in i:3
            m1=(c1.x[i]*c1.x[j])/c1.rx^2+(c1.y[i]*c1.y[j])/c1.ry^2
            m2=(c2.x[i]*c2.x[j])/c2.rx^2+(c2.y[i]*c2.y[j])/c2.ry^2
            abs(m1-m2)>100tol*max(1.0,abs(m1),abs(m2)) && return false
        end
        return true
    elseif c1.kind===:bspline
        c1.curve===c2.curve && return true
        (c1.curve isa _ModelOCCNURBSCurve &&
         c2.curve isa _ModelOCCNURBSCurve) || return false
        k1=c1.curve.curve;k2=c2.curve.curve
        (k1.degree==k2.degree && k1.knots==k2.knots &&
         k1.weights==k2.weights &&
         length(k1.controls)==length(k2.controls)) || return false
        return all(i->_b3dist(k1.controls[i],k2.controls[i])<100tol,
                   eachindex(k1.controls))
    end
    return false
end

# parameter of a 3-D point on a result curve, lifted near tref (conics)
function _brep_curve_param(crv,p,tref)
    if crv.kind===:line
        return _b3dot(_b3(p,crv.o),crv.d)
    elseif crv.kind===:circle
        rel=_b3(p,crv.c)
        θ=atan(_b3dot(rel,crv.y),_b3dot(rel,crv.x))
        return tref+mod(θ-tref+π,_OCC_TWO_PI)-π
    elseif crv.kind===:ellipse
        rel=_b3(p,crv.c)
        θ=atan(_b3dot(rel,crv.y)/crv.ry,_b3dot(rel,crv.x)/crv.rx)
        return tref+mod(θ-tref+π,_OCC_TWO_PI)-π
    elseif crv.kind===:bspline
        lo,hi=_model_occ_curve_bounds(crv.curve)
        res=_model_occ_curve_project(crv.curve,p,lo,hi)
        res===nothing && return nothing
        t,pp=res
        _b3dist(pp,p)<1e-6*max(1.0,_b3norm(p)) || return nothing
        return t
    end
    return nothing
end

# ── sections degenerating onto operand edges ────────────────────────────────
#
# When a section curve coincides with an operand edge the result edge IS that
# operand edge — the owning faces keep their wire piece keys, and every face
# of the other operand carrying the line receives an imprint arc keyed by the
# operand edge's piece so all references resolve to one result edge.

# pcurve under the affine parameter change t_old = α·t_new + β (α = ±1 for
# circles, arbitrary nonzero for lines). Returns (pcurve, neg): when a circle
# reparametrization reverses orientation the returned pcurve is evaluated at
# s = −t_new (neg=true) since a left-handed 2-D frame is not a `circ2d`.
function _brep_pc_reparam(pc,α,β,tol,caller)
    # a negated pcurve mirrors the affine map: t_old = −α·t_new − β — the
    # recursive neg flag composes into the returned (pcurve,neg) contract.
    hasproperty(pc,:neg2d) &&
        return _brep_pc_reparam(pc.neg2d[1],-α,-β,tol,caller)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=((fma(β,d[1],o[1]),fma(β,d[2],o[2])),
                       (α*d[1],α*d[2])),),false
    elseif hasproperty(pc,:circ2d)
        abs(abs(α)-1)>100tol && throw(ArgumentError(
            "$caller: coincident conic section with a non-angular " *
            "parameter map has no representable pcurve"))
        c,x=pc.circ2d[1],pc.circ2d[2]
        sb,cb=_gm_sincos(β)
        # frame rotated by β: eval(t) = pc(t+β); for α<0 evaluated at
        # s=−t it yields pc(−t+β) — exactly the reversed-parameter map
        xr=(fma(cb,x[1],-sb*x[2]),fma(cb,x[2],sb*x[1]))
        return (circ2d=(c,xr,pc.circ2d[3]),),α<0
    elseif hasproperty(pc,:ellips2d)
        abs(abs(α)-1)>100tol && throw(ArgumentError(
            "$caller: coincident ellipse section with a non-angular " *
            "parameter map has no representable pcurve"))
        c,x,rx,ry=pc.ellips2d
        y=(-x[2],x[1])
        sb,cb=_gm_sincos(β)
        # t_old = α·s + β expanded into conjugate diameters of the shifted
        # trace:  X' = rx·cosβ·x + α·ry·sinβ·y,
        #         Y' = −rx·sinβ·x + α·ry·cosβ·y.
        # ellips2d stores an axis-aligned frame, so this is representable
        # only when X'⊥Y' — i.e. rx≈ry or β≡0 (mod π/2). Otherwise the
        # mapped trace is sampled into spl2d below.
        Xp=(fma(rx*cb,x[1],α*ry*sb*y[1]),fma(rx*cb,x[2],α*ry*sb*y[2]))
        Yp=(fma(-rx*sb,x[1],α*ry*cb*y[1]),fma(-rx*sb,x[2],α*ry*cb*y[2]))
        dot=Xp[1]*Yp[1]+Xp[2]*Yp[2]
        nXp=hypot(Xp[1],Xp[2]); nYp=hypot(Yp[1],Yp[2])
        if abs(dot)<=1e-9*nXp*nYp && nXp>tol && nYp>tol
            xr=(Xp[1]/nXp,Xp[2]/nXp)
            # perp(xr)=(−xr[2],xr[1]); Yp = ry'·perp(xr) fixes the signed
            # minor radius (negative bakes the α<0 reflection into ry')
            ryp=Yp[1]*(-xr[2])+Yp[2]*xr[1]
            abs(ryp)>tol ||
                throw(ArgumentError(
                    "$caller: coincident ellipse section's parameter " *
                    "map collapses a semi-axis"))
            return (ellips2d=(c,xr,nXp,ryp),),false
        end
        # genuine axis-mismatch: sample the mapped conic into a spl2d
        # carrier so coincident-arc merging still succeeds. s∈[0,2π] maps
        # to t_old over the full period regardless of the sign of α.
        nS=48
        ss=collect(range(0.0,2π,nS+1))
        uvs=[let t=fma(α,s,β)
                 (c[1]+rx*cos(t)*x[1]+ry*sin(t)*y[1],
                  c[2]+rx*cos(t)*x[2]+ry*sin(t)*y[2])
             end for s in ss]
        return (spl2d=_pchip_bspline_2d(ss,uvs),),false
    elseif hasproperty(pc,:spl2d)
        # exact affine reparametrization of a B-spline: knots map through
        # (k−β)/α; α<0 mirrors the knot vector and pole order
        knots,p,cps=pc.spl2d
        out=if α>0
            (spl2d=(Float64[(k-β)/α for k in knots],p,cps),)
        else
            (spl2d=(Float64[(k-β)/α for k in Iterators.reverse(knots)],
                    p,NTuple{2,Float64}[cp for cp in
                                        Iterators.reverse(cps)]),)
        end
        hasproperty(pc,:lvl2d) && (out=merge(out,(lvl2d=pc.lvl2d,)))
        return out,false
    end
    throw(ArgumentError(
        "$caller: coincident section's pcurve has no representable " *
        "affine reparametrization"))
end

# pcurve of a 3-D curve known to lie in the target surface, parametrized by
# the curve's own parameter. The uv trace is classified as a line or a
# (possibly reversed) circle; any other trace has no representation and
# fails explicitly.
function _brep_pc_on_surface(sg,crv,t0,t1,tol,caller)
    ws=NTuple{2,Float64}[]
    us=NTuple{2,Float64}[]
    ts=(t0,t0+0.25(t1-t0),0.5(t0+t1),t0+0.75(t1-t0),t1)
    for t in ts
        push!(us,_brep_surface_uv(sg,_brep_curve_eval(crv,t)))
    end
    # an ellipse or B-spline trace on an untyped plane is exactly its frame
    # projection — use the synthesized record when it reproduces the samples
    if sg.occ===:plane && crv.kind in (:ellipse,:bspline)
        pc=_brep_synth_pcurve(sg,crv,t0,t1,caller)
        if all(i->(pe=_brep_pc_eval(pc,ts[i]);
                   abs(pe[1]-us[i][1])<100tol &&
                   abs(pe[2]-us[i][2])<100tol),eachindex(us))
            return pc,false
        end
    end
    uper=_brep_uperiod(sg)
    u0=us[1][1]
    push!(ws,us[1])
    for i in 2:length(us)
        u=uper>0 ? us[i][1]+uper*round((u0-us[i][1])/uper) : us[i][1]
        push!(ws,(u,us[i][2]))
    end
    # straight uv trace → lin2d (parametrization must stay affine in t)
    d=((ws[5][1]-ws[1][1])/(t1-t0),(ws[5][2]-ws[1][2])/(t1-t0))
    o=(fma(-t0,d[1],ws[1][1]),fma(-t0,d[2],ws[1][2]))
    if hypot(d[1],d[2])>tol && all(2:4) do i
            pe=_brep_pc_eval((lin2d=(o,d),),ts[i])
            abs(pe[1]-ws[i][1])<100tol && abs(pe[2]-ws[i][2])<100tol
        end
        return (lin2d=(o,d),),false
    end
    # circular trace: circumcenter of three samples, then frame and
    # handedness from the first — a reversed trace needs the negated
    # parameter since `circ2d` frames are right-handed
    (a,b,c)=ws[1],ws[3],ws[5]
    d1=(b[1]-a[1],b[2]-a[2]);d2=(c[1]-b[1],c[2]-b[2])
    m1=((a[1]+b[1])/2,(a[2]+b[2])/2);m2=((b[1]+c[1])/2,(b[2]+c[2])/2)
    det=d1[1]*d2[2]-d1[2]*d2[1]
    abs(det)<tol*hypot(d1[1],d1[2])*hypot(d2[1],d2[2]) && throw(ArgumentError(
        "$caller: coincident edge's uv trace is neither line nor circle " *
        "— unrepresentable pcurve"))
    r1=m1[1]*d1[1]+m1[2]*d1[2];r2=m2[1]*d2[1]+m2[2]*d2[2]
    ctr=((r1*d2[2]-r2*d1[2])/det,(d1[1]*r2-d2[1]*r1)/det)
    rr=hypot(ws[1][1]-ctr[1],ws[1][2]-ctr[2])
    rr<tol && throw(ArgumentError(
        "$caller: coincident edge degenerates to a point in the face " *
        "domain — unrepresentable pcurve"))
    for i in 2:length(ws)
        abs(hypot(ws[i][1]-ctr[1],ws[i][2]-ctr[2])-rr)<100tol ||
            throw(ArgumentError(
                "$caller: coincident edge's uv trace is neither line nor " *
                "circle — unrepresentable pcurve"))
    end
    w0=((ws[1][1]-ctr[1])/rr,(ws[1][2]-ctr[2])/rr)
    t0s=ts[1]
    for neg in (false,true)
        sgn=neg ? -1.0 : 1.0
        # frame mapping param s=sgn·t onto the trace: x = cos(s·t0)·w0 −
        # sin(s·t0)·rot90(w0) with rot90(u,v)=(−v,u)
        s0,c0=_gm_sincos(sgn*t0s)
        x=(fma(c0,w0[1],s0*w0[2]),fma(c0,w0[2],-s0*w0[1]))
        pc=(circ2d=(ctr,x,rr),)
        all(2:length(ws)) do i
            pe=_brep_pc_eval(pc,sgn*ts[i])
            abs(pe[1]-ws[i][1])<100tol && abs(pe[2]-ws[i][2])<100tol
        end && return pc,neg
    end
    throw(ArgumentError(
        "$caller: coincident edge's circular uv trace has a " *
        "non-angular parametrization — unrepresentable pcurve"))
end

# Contact pass for one edge trace `mpc` imprinted onto S's face `fidx`
# (O's edge `eOidx`, parametrized by s with edge param t = neg ? −s : s over
# [slo,shi]): crossings and collinear contacts with the face's wires and
# section arcs push splits onto both operands' edges (global — every face
# owning the edge sees them) and cuts onto the shared section spans.
# `skip_sec` excludes the degenerate section the imprint itself replaces.
function _brep_imprint_contacts!(S,O,fidx,eOidx,mpc,neg,slo,shi,es_S,es_O,
                                 sec_cuts,sections,S_is_A,tol,
                                 skip_sec::Int=0)
    f=S.faces[fidx]
    s2t(s)=neg ? -s : s
    # snap a first-arg pcurve param into [slo,shi] (circ2d intersections
    # report mod-2π; lin2d report unbounded)
    function inrange(s)
        if _brep_pc_iscirc(mpc)
            return _brep_snap_theta(s,slo,shi,tol)
        end
        return (s>=slo-tol && s<=shi+tol) ? s : nothing
    end
    for wire2 in f.wires,(eidx2,_,epc2) in wire2
        e2=S.edges[eidx2]
        e2.curve.kind===:degenerate && continue
        lo2,hi2=minmax(e2.t0,e2.t1)
        for t in _brep_pc_intersect(epc2,mpc,slo,shi,tol)
            push!(get!(Vector{Float64},es_S,eidx2),t)
        end
        for s in _brep_pc_intersect(mpc,epc2,lo2,hi2,tol)
            ss=inrange(s)
            ss===nothing ||
                push!(get!(Vector{Float64},es_O,eOidx),s2t(ss))
        end
        _brep_pc_same_carrier(epc2,mpc,tol) || continue
        # collinear contact: only endpoint projections delimit the
        # overlap — the transverse solve misses them entirely
        for te in (e2.t0,e2.t1)
            uv=_brep_pc_eval(epc2,te)
            sp=_brep_pc_param(mpc,uv,slo,shi,tol)
            sp===nothing ||
                push!(get!(Vector{Float64},es_O,eOidx),s2t(sp))
        end
        for se in (slo,shi)
            uv=_brep_pc_eval(mpc,se)
            tp=_brep_pc_param(epc2,uv,lo2,hi2,tol)
            tp===nothing ||
                push!(get!(Vector{Float64},es_S,eidx2),tp)
        end
    end
    for (si,sec) in enumerate(sections)
        si==skip_sec && continue
        (S_is_A ? sec.fa : sec.fb)==fidx || continue
        spc=S_is_A ? sec.pcA : sec.pcB
        # circ2d intersections report mod-2π — push every lift so spans at
        # other lifts are still cut
        nlifts=_brep_pc_iscirc(spc) ? (-1,0,1) : (0,)
        for t in _brep_pc_intersect(spc,mpc,slo,shi,tol), k in nlifts
            push!(get!(Vector{Float64},sec_cuts,si),t+k*_OCC_TWO_PI)
        end
        for span in sec.spans
            for s in _brep_pc_intersect(mpc,spc,span[1],span[2],tol)
                ss=inrange(s)
                ss===nothing ||
                    push!(get!(Vector{Float64},es_O,eOidx),s2t(ss))
            end
            _brep_pc_same_carrier(spc,mpc,tol) || continue
            for te in span
                uv=_brep_pc_eval(spc,te)
                sp=_brep_pc_param(mpc,uv,slo,shi,tol)
                sp===nothing ||
                    push!(get!(Vector{Float64},es_O,eOidx),s2t(sp))
            end
            for se in (slo,shi)
                uv=_brep_pc_eval(mpc,se)
                tp=_brep_pc_param(spc,uv,span[1],span[2],tol)
                tp===nothing ||
                    push!(get!(Vector{Float64},sec_cuts,si),tp)
            end
        end
    end
    return nothing
end

# Register one direction of a coincident pair: O's face `oidx` is imprinted
# onto S's face `fidx`. The mapped pcurve per edge is recorded for arc
# construction once all splits settle.
function _brep_coinc_imprints!(S,O,fidx,oidx,es_S,es_O,sec_cuts,imps,
                               sections,S_is_A,tol)
    f=S.faces[fidx];of=O.faces[oidx]
    mp=_brep_plane_map(f.sg,of.sg)
    seen=Set{Int}()
    for wire in of.wires,(eidx,sign,epc) in wire
        eidx in seen && continue
        push!(seen,eidx)
        e=O.edges[eidx]
        e.curve.kind===:degenerate && continue
        mpc,neg=_brep_map_pc(epc,mp)
        elo,ehi=minmax(e.t0,e.t1)
        slo,shi=neg ? (-ehi,-elo) : (elo,ehi)
        _brep_imprint_contacts!(S,O,fidx,eidx,mpc,neg,slo,shi,
                                es_S,es_O,sec_cuts,sections,S_is_A,tol)
        push!(imps,(eidx=eidx,pc=mpc,neg=neg,crv=e.curve))
    end
    return nothing
end

# Winner-edge parameter sub-spans for a degenerate-section span (lo,hi),
# returned as (σa,σb,ta,tb) tuples pairing each directed section-parameter
# bound with the matching winner-edge parameter in canonical edge
# coordinates. On a closed (full-period) edge a span can cross the
# parametric seam: endpoint projections alone are then ambiguous (antipodal
# endpoints bound two complementary arcs), so the section midpoint
# disambiguates and the wrapped span splits into canonical pieces. nothing
# when an endpoint does not project onto the edge.
function _brep_esec_tspans(we,scrv,lo,hi,tol)
    p1=_brep_curve_eval(scrv,lo)
    p2=_brep_curve_eval(scrv,hi)
    wmid=(we.t0+we.t1)/2
    t1=_brep_curve_param(we.curve,p1,wmid)
    t2=_brep_curve_param(we.curve,p2,wmid)
    (t1===nothing || t2===nothing) && return nothing
    per=(we.curve.kind in (:circle,:ellipse) &&
         abs(we.t1-we.t0-_OCC_TWO_PI)<1e-6) ? _OCC_TWO_PI : 0.0
    a,b=minmax(t1,t2)
    if per>0
        pm=_brep_curve_eval(scrv,0.5(lo+hi))
        tm=_brep_curve_param(we.curve,pm,wmid)
        if tm!==nothing && !any(k->a-tol<tm+k*per<b+tol,-1:1)
            # the midpoint projects onto the complementary arc: the true
            # directed span t1→t2 crosses the seam. Unwrap t2 onto the
            # branch containing the midpoint, then cut the span at every
            # domain boundary it crosses into canonical pieces.
            T2=t1>t2 ? t2+per : t2-per
            Tlo,Thi=minmax(t1,T2)
            if !any(k->Tlo-tol<tm+k*per<Thi+tol,-1:1)
                T2=t1>t2 ? t2-per : t2+per
                abs(T2-t1)<per-tol || return [(lo,hi,t1,t2)]
                Tlo,Thi=minmax(t1,T2)
                any(k->Tlo-tol<tm+k*per<Thi+tol,-1:1) ||
                    return [(lo,hi,t1,t2)]
            end
            σof=T->lo+(hi-lo)*(T-t1)/(T2-t1)
            subs=NTuple{4,Float64}[]
            klo=floor(Int,(Tlo-we.t0)/per)-1
            khi=ceil(Int,(Thi-we.t1)/per)+1
            for k in klo:khi
                plo=max(Tlo,we.t0+k*per)
                phi=min(Thi,we.t1+k*per)
                phi-plo>tol || continue
                push!(subs,(σof(plo),σof(phi),plo-k*per,phi-k*per))
            end
            return subs
        end
    end
    return [(lo,hi,t1,t2)]
end

# Imprint arcs for sections degenerating onto operand edges: the winner
# edge's pieces covering each recorded span, keyed like the winner's own
# wire pieces so every face's references resolve to one result edge. Each
# spec carries the target-face pcurve (parametrized by the winner edge's
# parameter, or its negation for `neg`).
function _brep_esec_arcs(W,L,f,esecs,es_W,wside,tol)
    arcs=NamedTuple[]
    emitted=Set{Tuple{Int,Int}}()
    for es in esecs
        e=W.edges[es.weidx]
        pieces=_brep_split_ranges(e.t0,e.t1,get(es_W,es.weidx,Float64[]),tol)
        for (k,(ta,tb)) in enumerate(pieces)
            any(sp->ta>=sp[1]-tol && tb<=sp[2]+tol,es.spans) || continue
            sa,sb=es.neg ? (-tb,-ta) : (ta,tb)
            _brep_uv_in_face(L.edges,f,_brep_pc_eval(es.pc,(sa+sb)/2),tol) ||
                continue
            # two degenerate sections can project the same winner-edge piece
            # onto this face — emit it once
            (es.weidx,k) in emitted && continue
            push!(emitted,(es.weidx,k))
            push!(arcs,(pc=es.pc,t0=sa,t1=sb,
                        key=(:e,W.ctags[es.weidx],ta,tb),crv=es.crv,
                        neg=es.neg,esec=true))
        end
    end
    return arcs
end

# Build the imprint arcs on face f: O's boundary-edge pieces (split by the
# finalized global edge splits) mapped into f's domain, kept where they lie
# strictly inside it and are not already carried by f's own wires or by a
# section arc. Arc keys carry O's canonical curve tag plus the edge span,
# so they match O's own wire pieces for the same edge piece.
function _brep_imprint_arcs(S,O,f,imps,es_O,secs_f,oside,also,tol)
    arcs=NamedTuple[]
    for im in imps
        e=O.edges[im.eidx]
        pieces=_brep_split_ranges(e.t0,e.t1,get(es_O,im.eidx,Float64[]),tol)
        for (k,(ta,tb)) in enumerate(pieces)
            sa,sb=im.neg ? (-tb,-ta) : (ta,tb)
            uvm=_brep_pc_eval(im.pc,(sa+sb)/2)
            covered=false
            for wire in f.wires,(eidx2,_,epc2) in wire
                e2=S.edges[eidx2]
                _brep_pc_same_carrier(im.pc,epc2,tol) || continue
                _brep_pc_param(epc2,uvm,e2.t0,e2.t1,tol)===nothing &&
                    continue
                covered=true;break
            end
            if !covered
                for sf in secs_f
                    _brep_pc_same_carrier(im.pc,sf.pc,tol) || continue
                    _brep_pc_param(sf.pc,uvm,sf.span[1],sf.span[2],tol)===
                        nothing && continue
                    covered=true;break
                end
            end
            # degenerate-section imprints and earlier imprint arcs already
            # carry the coincident trace — a duplicate would fold a region
            # to zero area
            if !covered
                covered=any(ex->_brep_pc_same_carrier(im.pc,ex.pc,tol) &&
                            _brep_pc_param(ex.pc,uvm,ex.t0,ex.t1,tol)!==
                            nothing, also)
            end
            if !covered
                covered=any(ex->_brep_pc_same_carrier(im.pc,ex.pc,tol) &&
                            _brep_pc_param(ex.pc,uvm,ex.t0,ex.t1,tol)!==
                            nothing, arcs)
            end
            covered && continue
            _brep_uv_in_face(S.edges,f,uvm,tol) || continue
            push!(arcs,(pc=im.pc,t0=sa,t1=sb,
                        key=(:e,O.ctags[im.eidx],ta,tb),crv=im.crv,neg=im.neg))
        end
    end
    return arcs
end

# ── union same-domain unification ───────────────────────────────────────────
#
# A fuse glues coplanar faces across operands: kept regions that share a
# result edge on one geometric surface merge into a single face spanning the
# combined domain (two side-by-side boxes yield one spanning face per plane,
# not split pieces). Applied to :union and :intersection — a difference keeps
# its imprint splits (OCC cut does not unify same-domain pieces).

# Can two kept specs merge into one face? Same geometric carrier, same
# outward normal direction, and a representable frame map. Plane pairs map
# affinely; other carriers merge only under identical uv frames (conservative
# — an unmapped merge is left as separate valid faces rather than
# approximated).
function _brep_unifiable(a,b,tol)
    a.sg.occ===b.sg.occ || return false
    if a.sg.occ===:plane
        na=_b3mul(a.sg.axis,Float64(a.sgn))
        nb=_b3mul(b.sg.axis,Float64(b.sgn))
        _b3dist(na,nb)>10tol && return false
        return _brep_planes_coincident(a.sg,b.sg)
    end
    a.sgn==b.sgn || return false
    for k in (:axis,:center,:X,:Y)
        hasproperty(a.sg,k) || continue
        _b3dist(getproperty(a.sg,k),getproperty(b.sg,k))>100tol &&
            return false
    end
    for k in (:radius,:radius1,:radius2,:half_angle,:r1,:r2)
        hasproperty(a.sg,k) || continue
        abs(getproperty(a.sg,k)-getproperty(b.sg,k))>100tol && return false
    end
    if a.sg.occ===:cone
        abs(a.sg.height-b.sg.height)>100tol && return false
    end
    # tensor carriers merge only across one shared patch record — two
    # distinct patches' uv domains are unrelated even when their
    # surfaces coincide
    if hasproperty(a.sg,:patch)
        (hasproperty(b.sg,:patch) && a.sg.patch===b.sg.patch) ||
            return false
    elseif hasproperty(b.sg,:patch)
        return false
    end
    return true
end

# Same-physical-surface certificate for a cross-operand face pair: same
# carrier kind with every geometric parameter equal — the uv frame fields
# (:X,:Y) are chart choices, not surface identity, so rotated frames of one
# cylinder still count; tensor carriers must share one patch record. A
# coincident pair is a same-domain contact, not a Boolean section.
function _brep_sg_identical(a,b,tol)
    a.occ===b.occ || return false
    for k in (:axis,:center)
        hasproperty(a,k) || continue
        _b3dist(getproperty(a,k),getproperty(b,k))>100tol &&
            return false
    end
    # carrier-identity scalar fields: :radius for sphere/cylinder,
    # :r1/:r2 for cone (frustum radii — with the axial :height they fix
    # apex and slope) and torus (spine and tube radii). A torus's :angle
    # and a cylinder's :height are chart extents, not carrier identity.
    for k in (:radius,:radius1,:radius2,:half_angle,:r1,:r2)
        hasproperty(a,k) || continue
        abs(getproperty(a,k)-getproperty(b,k))>100tol && return false
    end
    if a.occ===:cone
        # same apex+frame but a different slant extent is a different cone
        abs(a.height-b.height)>100tol && return false
    end
    if hasproperty(a,:patch)
        (hasproperty(b,:patch) && a.patch===b.patch) || return false
    elseif hasproperty(b,:patch)
        return false
    end
    return true
end

# Merge the kept specs of one edge-connected component into unified faces.
# Every member wire's arcs map into the reference member's uv frame; arcs
# whose edge key is used by more than one member are internal to the merged
# domain and drop out. The remainder is re-traversed and each resulting
# bounded face that overlaps a member region becomes a merged face spec.
function _brep_merge_union_faces(kept,comp,tol,caller)
    ref=kept[comp[1]]
    maps=Dict{Int,Any}(mi => (kept[mi].sg.occ===:plane ?
        _brep_plane_map(ref.sg,kept[mi].sg) : nothing) for mi in comp)
    # Periodic-carrier unwrap: two member pieces of one periodic patch are
    # split by the canonical chart cut — their shared result edge's two
    # pcurve copies sit whole periods apart (the face's own seam edge).
    # OCC unifies them into one face, so place every member's chart at the
    # whole-period offset that lands its shared arcs on an already-placed
    # member's position; the coincident copies then drop as internal below
    # and the merged domain reads contiguously across the seam. Members
    # sharing no key with a placed member keep the canonical chart — any
    # non-coincident shared arcs simply stay boundary (previous behavior).
    uper=_brep_uperiod(ref.sg)
    # interior side of each arc at its stored direction: an outer (CCW)
    # wire puts the region left of a forward traversal, a hole (CW) wire
    # puts it right — used by the sub-span certificate to require both
    # regions on the same side of the shared edge, and by the coincident
    # drop below to certify a wall carried by distinct edges is internal
    # to the merged domain (opposite sides)
    member_side=Dict{Int,Vector{Int}}(
        mi=>Int[(wi==1)==ar.fwd ? 1 : -1
                for (wi,w) in enumerate(kept[mi].wires) for ar in w]
        for mi in comp)
    ushift=Dict{Int,Float64}()
    if uper>0
        ushift[comp[1]]=0.0
        member_arcs=Dict{Int,Vector{NamedTuple}}(
            mi=>NamedTuple[ar for w in kept[mi].wires for ar in w]
            for mi in comp)
        # constant-u ("slit") arcs per member — `(marc base index is the
        # flattened position below)`, u position and v extent in the
        # member's own chart. Two such arcs a whole u-period apart trace
        # the same physical slit line — the members bound the chart cut
        # through distinct edges (a periodic face whose window ends are
        # carried by separate edges, e.g. an imported trimmed periodic
        # face or an imprint on the cut) and share no result-edge key.
        member_cu=Dict{Int,
            Vector{Tuple{Int,NamedTuple,Float64,Float64,Float64,Int}}}(
            mi=>Tuple{Int,NamedTuple,Float64,Float64,Float64,Int}[]
            for mi in comp)
        let flat=0
            for mi in comp
                ai=0
                for ar in member_arcs[mi]
                    flat+=1;ai+=1
                    # lin2d-only: the residual v→t cut below is exact
                    # for a straight uv trace — a curved slit trace
                    # stays conservative (no partner, previous behavior)
                    hasproperty(ar.pc,:lin2d) || continue
                    u0,u1=_brep_pc_urange(ar.pc,ar.t0,ar.t1)
                    abs(u1-u0)<tol || continue
                    v0,v1=_brep_pc_vrange(ar.pc,ar.t0,ar.t1)
                    push!(member_cu[mi],
                          (flat,ar,u0,min(v0,v1),max(v0,v1),
                           member_side[mi][ai]))
                end
            end
        end
        grown=true
        while grown
            grown=false
            for mi in comp
                haskey(ushift,mi) || continue
                si=ushift[mi]
                for (ai,ar) in enumerate(member_arcs[mi])
                    mu=_brep_pc_eval(ar.pc,(ar.t0+ar.t1)/2)
                    for mj in comp
                        mj==mi && continue
                        haskey(ushift,mj) && continue
                        for (bi,br) in enumerate(member_arcs[mj])
                            br.key==ar.key || continue
                            if abs((br.t1-br.t0)-(ar.t1-ar.t0))<tol
                                bu=_brep_pc_eval(br.pc,(br.t0+br.t1)/2)
                                # the two copies of one periodic-edge span
                                # carry identical v; a whole-period u
                                # offset aligns them. When the copies are
                                # not coincident the pair may still be
                                # adjacent pieces of one seam-split span,
                                # so fall through to the junction
                                # certificate below instead of skipping.
                                if abs(bu[2]-mu[2])<=10tol
                                    k=round(Int,(mu[1]+si-bu[1])/uper)
                                    sj=k*uper
                                    if abs(mu[1]+si-(bu[1]+sj))<=10tol
                                        ushift[mj]=sj;grown=true;break
                                    end
                                end
                            end
                            # Sub-span shared key: the canonical u-seam
                            # split cuts one edge piece into adjacent
                            # spans that keep the parent's key and land
                            # in different members — a periodic carrier
                            # whose face wires hold no seam edge (tensor
                            # patch, imported trimmed periodic face)
                            # gives no coincident-span copy to place the
                            # far member through. Two adjacent same-key
                            # sub-spans still certify the period shift:
                            # the only producer of such pairs is the
                            # seam split, so their shared endpoint sits
                            # exactly on a u-seam line — evaluate both
                            # pcurves at the junction and require the
                            # same v at a u a whole number of periods
                            # apart (non-adjacent or overlapping
                            # sub-spans prove nothing and stay
                            # conservative). The members must also bound
                            # the shared edge from the same side — the
                            # two pieces continue through the junction
                            # along one curve, so equal interior sides
                            # put both regions on one side of it and the
                            # union tiles the seam wedge; opposite
                            # senses only touch at the point and a shift
                            # would weld a pinched domain.
                            tj=abs(ar.t1-br.t0)<=tol ? (ar.t1,br.t0) :
                               (abs(br.t1-ar.t0)<=tol ?
                                (ar.t0,br.t1) : nothing)
                            tj===nothing && continue
                            member_side[mi][ai]==
                                member_side[mj][bi] || continue
                            muj=_brep_pc_eval(ar.pc,tj[1])
                            buj=_brep_pc_eval(br.pc,tj[2])
                            abs(buj[2]-muj[2])>10tol && continue
                            k=round(Int,(muj[1]+si-buj[1])/uper)
                            sj=k*uper
                            abs(muj[1]+si-(buj[1]+sj))>10tol &&
                                continue
                            ushift[mj]=sj;grown=true;break
                        end
                    end
                end
                # Slit partners: two members may bound the chart cut
                # through distinct constant-u edges a whole period
                # apart whose v extents overlap — the same physical
                # slit line carried by different edges (an imported
                # trimmed periodic face stores the two window uses as
                # separate edges, or an imprint lies on the cut) — so
                # their domains are contiguous across it. Certify at
                # the overlap midpoint: the two traces must land on one
                # uv point under the candidate whole-period shift. A
                # same-key slit edge is already matched by the loops
                # above.
                for mj in comp
                    mj==mi && continue
                    haskey(ushift,mj) && continue
                    placed=false
                    for (_,sa,ua,va0,va1,sidea) in member_cu[mi],
                        (_,sb,ub,vb0,vb1,sideb) in member_cu[mj]
                        abs(abs(ua-ub)-uper)>10tol && continue
                        # a real slit wall is internal to the union —
                        # the members must bound it from opposite sides
                        sidea==sideb && continue
                        vl=max(va0,vb0);vh=min(va1,vb1)
                        vh-vl>tol || continue
                        vm=(vl+vh)/2
                        ta=_brep_pc_param(sa.pc,(ua,vm),sa.t0,sa.t1,tol)
                        tb=_brep_pc_param(sb.pc,(ub,vm),sb.t0,sb.t1,tol)
                        (ta===nothing || tb===nothing) && continue
                        muj=_brep_pc_eval(sa.pc,ta)
                        buj=_brep_pc_eval(sb.pc,tb)
                        abs(muj[2]-buj[2])>10tol && continue
                        k=round(Int,(muj[1]+si-buj[1])/uper)
                        sj=k*uper
                        abs(muj[1]+si-(buj[1]+sj))>10tol && continue
                        ushift[mj]=sj;grown=true;placed=true;break
                    end
                    placed && break
                end
            end
        end
    end
    if get(ENV,"BREP_DEBUG_UNWRAP","")!="" && uper>0
        unplaced=[mi for mi in comp if !haskey(ushift,mi)]
        println("DBGUNWRAP comp=$comp uper=$uper ushift=$ushift unplaced=$unplaced")
        for mi in comp, (wi,w) in enumerate(kept[mi].wires), ar in w
            mu=_brep_pc_eval(ar.pc,(ar.t0+ar.t1)/2)
            println("DBGUNWRAP   m$mi w$wi key=$(ar.key) t=($(ar.t0),$(ar.t1)) " *
                    "vt=($(ar.vt0),$(ar.vt1)) uvmid=$(round.(mu,digits=4)) neg=$(ar.neg)")
        end
    end
    marcs=NamedTuple[]
    test_arcs=NamedTuple[]
    member_cycles=Vector{Vector{Vector{Tuple{Int,Bool}}}}()
    keylocs=Dict{Any,Vector{Tuple{Int,Int}}}()   # key -> [(member, marcs idx)]
    marc_member=Int[]                            # marcs idx -> member
    for mi in comp
        mp=maps[mi]
        sj=get(ushift,mi,0.0)
        cycs=Vector{Tuple{Int,Bool}}[]
        for w in kept[mi].wires
            cy=Tuple{Int,Bool}[]
            for ar in w
                rec=ar
                if mp!==nothing
                    mpc,flip=_brep_map_pc(ar.pc,mp)
                    a0=ar.t0;a1=ar.t1
                    if flip
                        a0=-ar.t1;a1=-ar.t0
                    end
                    n2=xor(ar.neg,flip)
                    rec=(key=ar.key,fwd=ar.fwd,pc=mpc,t0=a0,t1=a1,
                         crv=ar.crv,vt0=(n2 ? -a1 : a0),
                         vt1=(n2 ? -a0 : a1),neg=n2)
                elseif sj!=0.0
                    # `_brep_shift_pc` subtracts its argument from the
                    # evaluated u — negate to land the arc at +sj
                    rec=(key=ar.key,fwd=ar.fwd,
                         pc=_brep_shift_pc(ar.pc,-sj),
                         t0=ar.t0,t1=ar.t1,crv=ar.crv,
                         vt0=ar.vt0,vt1=ar.vt1,neg=ar.neg)
                end
                push!(test_arcs,rec)
                push!(cy,(length(test_arcs),ar.fwd))
                push!(marcs,merge(rec,(a=0,b=0)))
                push!(marc_member,mi)
                push!(get!(Vector{Tuple{Int,Int}},keylocs,ar.key),
                      (mi,length(marcs)))
            end
            push!(cycs,cy)
        end
        push!(member_cycles,cycs)
    end
    isempty(marcs) && return [kept[mi] for mi in comp]
    # An arc is internal to the merged domain only when another member uses
    # the same edge over the same uv segment. Period-shifted copies of a
    # periodic seam edge share the key at different uv positions; they are
    # distinct domain boundaries and must be kept, or the merged region would
    # wrap the seam and lose its planar cycle representation — unless the
    # members were unwrapped onto one contiguous domain above, which lands
    # the copies at the same uv so they cancel like any shared edge.
    keeparc=trues(length(marcs))
    for lst in values(keylocs)
        length(lst)>1 || continue
        for x in 1:length(lst), y in x+1:length(lst)
            mi,ix=lst[x];mj,iy=lst[y]
            mi==mj && continue
            ax=marcs[ix];ay=marcs[iy]
            abs((ax.t1-ax.t0)-(ay.t1-ay.t0))<tol || continue
            mx=_brep_pc_eval(ax.pc,(ax.t0+ax.t1)/2)
            my=_brep_pc_eval(ay.pc,(ay.t0+ay.t1)/2)
            hypot(mx[1]-my[1],mx[2]-my[2])<tol || continue
            keeparc[ix]=keeparc[iy]=false
        end
    end
    # The same physical wall may reach the merge through distinct edges —
    # a periodic face whose window ends are stored as separate edges, an
    # imprint lying on a boundary, or a seam slit carried by non-lin2d
    # pcurves the overlap pass below cannot cut. Two such arcs bound one
    # uv segment as an internal wall only when the members' regions lie
    # on opposite sides of it — sample each member's domain on both
    # sides of the coincident midpoint (the member_cycles polygon is
    # already in this reference chart, so no map bookkeeping is needed);
    # equal-side or uncovered coincidences stay boundary.
    cpos=Dict{Int,Int}(mi=>k for (k,mi) in enumerate(comp))
    in_member=function(mpos,uv)
        cycs=member_cycles[mpos]
        _brep_cycle_inside(test_arcs,cycs[1],uv,tol) || return false
        any(h->_brep_cycle_inside(test_arcs,h,uv,tol),cycs[2:end]) &&
            return false
        return true
    end
    side_of=function(mpos,uv,n)
        off=100tol
        pin=(uv[1]+n[1]*off,uv[2]+n[2]*off)
        pout=(uv[1]-n[1]*off,uv[2]-n[2]*off)
        in_member(mpos,pin) && return 1
        in_member(mpos,pout) && return -1
        return 0
    end
    for i in 1:length(marcs)-1, j in i+1:length(marcs)
        mi=marc_member[i];mj=marc_member[j]
        mi==mj && continue
        ax=marcs[i];ay=marcs[j]
        ax.key==ay.key && continue
        abs((ax.t1-ax.t0)-(ay.t1-ay.t0))<tol || continue
        mx=_brep_pc_eval(ax.pc,(ax.t0+ax.t1)/2)
        my=_brep_pc_eval(ay.pc,(ay.t0+ay.t1)/2)
        hypot(mx[1]-my[1],mx[2]-my[2])<tol || continue
        d=_brep_pc_tangent(ax.pc,(ax.t0+ax.t1)/2)
        L=hypot(d[1],d[2]);L<tol && continue
        n=(-d[2]/L,d[1]/L)
        sa=side_of(cpos[mi],mx,n)
        sb=side_of(cpos[mj],mx,n)
        (sa==0 || sb==0 || sa==sb) && continue
        keeparc[i]=keeparc[j]=false
    end
    # Different-edge slit partners: two uses of one physical slit line
    # carried by distinct edges land coincident once the members are
    # unwrapped — the shared part is internal to the unified domain.
    # The overlap drops; each use's residual ends stay boundary (the
    # slit arcs are lin2d by the partner certificate, so the residual
    # v→parameter cut is exact).
    extra_marcs=NamedTuple[]
    if uper>0 && !isempty(member_cu)
        drop_iv=Dict{Int,Vector{Tuple{Float64,Float64}}}()
        for (ix,mi) in enumerate(comp),jx in ix+1:length(comp)
            mj=comp[jx]
            si=get(ushift,mi,0.0);sjb=get(ushift,mj,0.0)
            for (fi,sa,ua,va0,va1,sidea) in member_cu[mi],
                (fj,sb,ub,vb0,vb1,sideb) in member_cu[mj]
                sa.key==sb.key && continue
                # only a wall bounded by the two members from opposite
                # sides is internal to the merged domain
                sidea==sideb && continue
                vl=max(va0,vb0);vh=min(va1,vb1)
                vh-vl>tol || continue
                # the partner copies land coincident only when both
                # members were placed consistently — a member that
                # failed to place leaves its arc's u a whole period
                # off, where the pair can stay boundary as before
                abs(ua+si-(ub+sjb))>10tol && continue
                push!(get!(Vector{Tuple{Float64,Float64}},drop_iv,fi),
                      (vl,vh))
                push!(get!(Vector{Tuple{Float64,Float64}},drop_iv,fj),
                      (vl,vh))
            end
        end
        for (ix,ivs) in drop_iv
            keeparc[ix] || continue
            ar=marcs[ix]
            o,d=ar.pc.lin2d
            v0,v1=_brep_pc_vrange(ar.pc,ar.t0,ar.t1)
            sort!(ivs)
            res=Tuple{Float64,Float64}[]
            cur=v0
            for (dl,dh) in ivs
                dl-cur>tol && push!(res,(cur,dl))
                cur=max(cur,dh)
            end
            v1-cur>tol && push!(res,(cur,v1))
            keeparc[ix]=false
            for (rl,rh) in res
                rh-rl>tol || continue
                ta=(rl-o[2])/d[2];tb=(rh-o[2])/d[2]
                q0,q1=minmax(ta,tb)
                (q0>=ar.t0-tol && q1<=ar.t1+tol) || continue
                f0=(q0-ar.t0)/(ar.t1-ar.t0)
                f1=(q1-ar.t0)/(ar.t1-ar.t0)
                push!(extra_marcs,
                      (key=ar.key,fwd=ar.fwd,pc=ar.pc,t0=q0,t1=q1,
                       crv=ar.crv,
                       vt0=ar.vt0+f0*(ar.vt1-ar.vt0),
                       vt1=ar.vt0+f1*(ar.vt1-ar.vt0),
                       neg=ar.neg,a=0,b=0))
            end
        end
    end
    marcs=marcs[keeparc]
    append!(marcs,extra_marcs)
    if isempty(marcs)
        # every member arc was claimed by another member at the same uv —
        # the members cover one identical domain (a duplicate split), so
        # the unified result is a single face carrying the union of the
        # members' cell memberships, like OCC's same-domain merge
        spec=(sg=ref.sg,sgn=ref.sgn,wires=ref.wires,stag=ref.stag)
        hasproperty(ref,:sev) && (spec=merge(spec,(sev=ref.sev,)))
        if hasproperty(ref,:uses)
            uses=Tuple{UInt64,Int}[]
            for mi in comp, u in kept[mi].uses
                u in uses || push!(uses,u)
            end
            spec=merge(spec,(uses=uses,))
        end
        return [spec]
    end
    pts=Tuple{Float64,Float64}[]
    for ar in marcs
        push!(pts,_brep_pc_eval(ar.pc,ar.t0))
        push!(pts,_brep_pc_eval(ar.pc,ar.t1))
    end
    nodes,idx=_brep_uv_nodes(pts,tol)
    narcs=[merge(ar,(a=idx[2i-1],b=idx[2i])) for (i,ar) in enumerate(marcs)]
    cycles=_brep_uv_arrange(narcs,length(nodes),tol,caller)
    areas=[_brep_cycle_area(narcs,cy) for cy in cycles]
    # merged-domain containment: inside a member's outer wire (cycs[1]) and
    # outside its hole wires
    in_merged=function(uv)
        for cycs in member_cycles
            isempty(cycs) && continue
            _brep_cycle_inside(test_arcs,cycs[1],uv,tol) || continue
            any(h->_brep_cycle_inside(test_arcs,h,uv,tol),cycs[2:end]) ||
                return true
        end
        return false
    end
    pos=[i for i in eachindex(cycles) if areas[i]>tol]
    neg=[i for i in eachindex(cycles) if areas[i]<-tol]
    for i in eachindex(cycles)
        abs(areas[i])<=tol && throw(ErrorException(
            "$caller: degenerate zero-area wire in Boolean face merge"))
    end
    samples=[_brep_cycle_sample(narcs,cycles[i],areas[i],tol,caller)
             for i in eachindex(cycles)]
    smallest_containing=function(uv,excl)
        best=0;ba=Inf
        for p in pos
            p==excl && continue
            abs(areas[p])<ba || continue
            _brep_cycle_inside(narcs,cycles[p],uv,tol) || continue
            best=p;ba=abs(areas[p])
        end
        best
    end
    out=NamedTuple[]
    for c in pos
        in_merged(samples[c]) || continue
        wires=[ [(key=narcs[i].key,fwd=fwd,pc=narcs[i].pc,
                  t0=narcs[i].t0,t1=narcs[i].t1,crv=narcs[i].crv,
                  vt0=narcs[i].vt0,vt1=narcs[i].vt1,neg=narcs[i].neg)
                for (i,fwd) in cycles[c]] ]
        for h in neg
            d=smallest_containing(samples[h],0)
            d==0 && continue
            parent=smallest_containing(samples[d],d)
            parent==c || continue
            push!(wires,[(key=narcs[i].key,fwd=fwd,pc=narcs[i].pc,
                          t0=narcs[i].t0,t1=narcs[i].t1,crv=narcs[i].crv,
                          vt0=narcs[i].vt0,vt1=narcs[i].vt1,neg=narcs[i].neg)
                         for (i,fwd) in cycles[h]])
        end
        spec=(sg=ref.sg,sgn=ref.sgn,wires=wires,stag=ref.stag)
        hasproperty(ref,:sev) && (spec=merge(spec,(sev=ref.sev,)))
        if hasproperty(ref,:uses)
            # N-way specs keep their cell memberships — the merged face
            # carries the union of its members' uses
            uses=Tuple{UInt64,Int}[]
            for mi in comp, u in kept[mi].uses
                u in uses || push!(uses,u)
            end
            spec=merge(spec,(uses=uses,))
        end
        push!(out,spec)
    end
    return out
end

function _brep_unify_union(kept,tol,caller;same_uses::Bool=false)
    reps=Int[];gid=zeros(Int,length(kept))
    for (i,s) in enumerate(kept)
        j=findfirst(reps) do r
            _brep_unifiable(kept[r],s,tol) &&
                (!same_uses || Set(kept[r].uses)==Set(s.uses))
        end
        j===nothing && (push!(reps,i);j=length(reps))
        gid[i]=j
    end
    out=NamedTuple[]
    for g in eachindex(reps)
        mem=findall(==(g),gid)
        if length(mem)==1
            push!(out,kept[mem[1]]);continue
        end
        # connected components via shared result-edge keys
        uf=Dict{Int,Int}(mi=>mi for mi in mem)
        findr=function(x)
            while uf[x]!=x; uf[x]=uf[uf[x]];x=uf[x] end
            x
        end
        keyowner=Dict{Any,Int}()
        for mi in mem,w in kept[mi].wires,ar in w
            o=get(keyowner,ar.key,0)
            if o==0
                keyowner[ar.key]=mi
            else
                ra,rb=findr(mi),findr(o)
                ra!=rb && (uf[ra]=rb)
            end
        end
        # slit partners carry no shared key: two members may bound the
        # periodic chart cut through distinct constant-u edges a whole
        # u-period apart (a periodic face whose window ends are stored
        # as separate edges — e.g. an imported trimmed periodic face —
        # or an imprint lying on the cut). The traces are the same
        # physical slit line, so whenever their v extents overlap AND
        # the members bound them from opposite sides — a real slit wall
        # is internal — the members' domains are contiguous across the
        # slit: group them. Same-side traces are duplicate-domain
        # artifacts the merge certificates reject anyway.
        uperg=_brep_uperiod(kept[reps[g]].sg)
        if uperg>0
            for (ix,mi) in enumerate(mem),jx in ix+1:length(mem)
                mj=mem[jx]
                findr(mi)==findr(mj) && continue
                hit=false
                for (wi,wa) in enumerate(kept[mi].wires),aa in wa
                    hasproperty(aa.pc,:lin2d) || continue
                    ua0,ua1=_brep_pc_urange(aa.pc,aa.t0,aa.t1)
                    abs(ua1-ua0)<tol || continue
                    va0,va1=_brep_pc_vrange(aa.pc,aa.t0,aa.t1)
                    sa=(wi==1)==aa.fwd ? 1 : -1
                    for (wj,wb) in enumerate(kept[mj].wires),bb in wb
                        hasproperty(bb.pc,:lin2d) || continue
                        ub0,ub1=_brep_pc_urange(bb.pc,bb.t0,bb.t1)
                        abs(ub1-ub0)<tol || continue
                        abs(abs(ua0-ub0)-uperg)>10tol && continue
                        vb0,vb1=_brep_pc_vrange(bb.pc,bb.t0,bb.t1)
                        min(va1,vb1)-max(va0,vb0)>tol || continue
                        sb=(wj==1)==bb.fwd ? 1 : -1
                        sa==sb && continue
                        hit=true;break
                    end
                    hit && break
                end
                hit && (uf[findr(mi)]=findr(mj))
            end
        end
        comps=Dict{Int,Vector{Int}}()
        for mi in mem
            push!(get!(Vector{Int},comps,findr(mi)),mi)
        end
        # components emit in first-face order — `Dict` iteration order
        # would make result order depend on the runtime's hashing
        for comp in sort!(collect(values(comps));by=minimum)
            if length(comp)==1
                push!(out,kept[comp[1]])
            else
                append!(out,_brep_merge_union_faces(kept,comp,tol,caller))
            end
        end
    end
    out=_brep_merge_edge_pieces(out,tol,caller)
    return out
end

# affine map v_new → v_old between two result curves on one carrier, or
# nothing when the carrier kinds do not afford one (degenerate/bspline)
function _brep_canonical_map(cm,cc,vm0,tol)
    if cc.kind===:line
        cm.kind===:line || return nothing
        dc2=_b3dot(cc.d,cc.d)
        α=_b3dot(cm.d,cc.d)/dc2
        β=_b3dot(_b3(cm.o,cc.o),cc.d)/dc2
        return α,β
    elseif cc.kind===:circle && cm.kind===:circle
        _b3dist(cm.c,cc.c)<100tol || return nothing
        abs(cm.r-cc.r)<100tol || return nothing
        nn=_b3dot(cm.n,cc.n)
        (abs(abs(nn)-1)>100tol) && return nothing
        α=nn>0 ? 1.0 : -1.0
        pA=_brep_curve_eval(cm,vm0)
        θA=_brep_curve_param(cc,pA,α*vm0)
        θA===nothing && return nothing
        return α,θA-α*vm0
    elseif cc.kind===:ellipse && cm.kind===:ellipse
        # same angular convention as circles — the metric check gates the
        # shared carrier; matched radii keep the map affine
        _brep_same_curve3d(cm,cc,tol) || return nothing
        (abs(cm.rx-cc.rx)<100tol && abs(cm.ry-cc.ry)<100tol) ||
            return nothing
        nn=_b3dot(cm.n,cc.n)
        (abs(abs(nn)-1)>100tol) && return nothing
        α=nn>0 ? 1.0 : -1.0
        pA=_brep_curve_eval(cm,vm0)
        θA=_brep_curve_param(cc,pA,α*vm0)
        θA===nothing && return nothing
        return α,θA-α*vm0
    elseif cc.kind===:bspline && cm.kind===:bspline
        # one shared evaluator → the same parameter is the canonical map;
        # distinct evaluators stay unmerged (conservative)
        cm.curve===cc.curve || return nothing
        return 1.0,0.0
    end
    return nothing
end

# Fuse piece-edges: result-edge keys on one carrier sharing a 3-D endpoint
# and bounded by the same pair of faces merge into a single edge (OCC's
# same-domain edge gluing — the vertex between two collinear pieces of one
# boundary is not a result vertex). Every arc keyed to a member is relabeled
# into the canonical edge's parametrization, then consecutive arcs of one
# key inside a wire join.
function _brep_merge_edge_pieces(kept,tol,caller)
    scale=1.0
    for s in kept,w in s.wires,ar in w
        for t in (ar.vt0,ar.vt1)
            p=_brep_vertex_eval(ar.crv,t)
            scale=max(scale,abs(p[1]),abs(p[2]),abs(p[3]))
        end
    end
    tol3=1e-9*scale
    keyfaces=Dict{Any,Set{Int}}()
    keyrep=Dict{Any,NamedTuple}()
    # every key whose span ends at a junction point — a junction used by a
    # third edge key is a real degree-3+ vertex, not a removable piece split
    junction=Dict{NTuple{3,Float64},Set{Any}}()
    jpts=NTuple{3,Float64}[]
    for (fi,s) in enumerate(kept),w in s.wires,ar in w
        push!(get!(Set{Int},keyfaces,ar.key),fi)
        haskey(keyrep,ar.key) ||
            (keyrep[ar.key]=(crv=ar.crv,vt0=ar.vt0,vt1=ar.vt1))
        for tt in (ar.vt0,ar.vt1)
            p=_brep_vertex_eval(ar.crv,tt)
            i=findfirst(q->_b3dist(q,p)<tol3,jpts)
            i===nothing && (push!(jpts,p);i=length(jpts))
            push!(get!(Set{Any},junction,jpts[i]),ar.key)
        end
    end
    ks=collect(keys(keyrep))
    uf=Dict{Any,Any}(k=>k for k in ks)
    function findk(k)
        r=k
        while uf[r]!=r; r=uf[r] end
        uf[k]=r
        return r
    end
    for i in eachindex(ks), j in i+1:length(ks)
        k1,k2=ks[i],ks[j]
        keyfaces[k1]==keyfaces[k2] || continue
        r1,r2=keyrep[k1],keyrep[k2]
        _brep_same_curve3d(r1.crv,r2.crv,tol) || continue
        # a junction pairing one piece's domain-start to the other's
        # domain-end on the same bounded evaluator is the carrier's own
        # periodic wrap, not a removable vertex: the merged edge would
        # span the whole closed domain and appear twice in one wire
        wrap=false
        if r1.crv.kind===:bspline && r1.crv.curve===r2.crv.curve
            blo,bhi=_model_occ_curve_bounds(r1.crv.curve)
            if bhi>blo
                for s1 in (r1.vt0,r1.vt1), s2 in (r2.vt0,r2.vt1)
                    _b3dist(_brep_vertex_eval(r1.crv,s1),
                            _brep_vertex_eval(r2.crv,s2))<tol3 || continue
                    ((abs(s1-blo)<tol3 && abs(s2-bhi)<tol3) ||
                     (abs(s1-bhi)<tol3 && abs(s2-blo)<tol3)) &&
                        (wrap=true;break)
                end
            end
        end
        wrap && continue
        adj=false
        for s1 in (r1.vt0,r1.vt1), s2 in (r2.vt0,r2.vt1)
            _b3dist(_brep_vertex_eval(r1.crv,s1),
                    _brep_vertex_eval(r2.crv,s2))<tol3 || continue
            p=_brep_vertex_eval(r1.crv,s1)
            i=findfirst(q->_b3dist(q,p)<tol3,jpts)
            # merge through this junction only when no third key ends there
            if i!==nothing && issubset(junction[jpts[i]],Set{Any}((k1,k2)))
                adj=true;break
            end
        end
        adj || continue
        # every member must afford an affine map onto the root's carrier —
        # otherwise the pieces stay distinct rather than misparametrized
        findk(k1)==findk(k2) && continue
        uf[findk(k1)]=findk(k2)
    end
    # components → canonical (root) key; verify affine maps for all members
    comps=Dict{Any,Vector{Any}}()
    for k in ks
        push!(get!(Vector{Any},comps,findk(k)),k)
    end
    canon=Dict{Any,Tuple{Any,Float64,Float64}}()   # key → (root,α,β)
    for (root,mem) in comps
        length(mem)>1 || continue
        cc=keyrep[root].crv
        maps=NamedTuple[]
        ok=true
        for k in mem
            r=keyrep[k]
            ab=_brep_canonical_map(r.crv,cc,r.vt0,tol)
            ab===nothing && (ok=false;break)
            push!(maps,(k=k,α=ab[1],β=ab[2]))
        end
        ok || continue
        # a closed conic's canonical map is unique only mod 2π — members
        # meeting the root through the parametrization seam must unwrap onto
        # the covered span or their canonical intervals stay disjoint
        per=cc.kind in (:circle,:ellipse) ? _OCC_TWO_PI : 0.0
        covered=Tuple{Float64,Float64}[]
        iroot=findfirst(mp->mp.k==root,maps)
        if iroot===nothing
            mp=maps[1];r=keyrep[mp.k]
            push!(covered,minmax(mp.α*r.vt0+mp.β,mp.α*r.vt1+mp.β))
        else
            mp=maps[iroot];r=keyrep[root]
            push!(covered,minmax(mp.α*r.vt0+mp.β,mp.α*r.vt1+mp.β))
        end
        for mp in maps
            α,β=mp.α,mp.β
            if per>0 && mp.k!=root
                r=keyrep[mp.k]
                lo,hi=minmax(α*r.vt0+β,α*r.vt1+β)
                bg=Inf;bk=0
                for ks in (-2,-1,0,1,2)
                    l2=lo+ks*per;h2=hi+ks*per
                    g=Inf
                    for (cl,ch) in covered
                        g=min(g,l2>ch ? l2-ch : (h2<cl ? cl-h2 : 0.0))
                    end
                    g<bg-tol && (bg=g;bk=ks)
                end
                β+=bk*per
                push!(covered,(lo+bk*per,hi+bk*per))
            end
            canon[mp.k]=(root,α,β)
        end
    end
    isempty(canon) && return kept
    out=NamedTuple[]
    for s in kept
        # join consecutive arc-uses of one merged edge inside a wire: the
        # chain links exit parameter to entry parameter (each is t0 or t1
        # depending on traversal direction)
        exitp(ar)=ar.fwd ? ar.t1 : ar.t0
        entryp(ar)=ar.fwd ? ar.t0 : ar.t1
        joinable(last,rec)=last.key==rec.key && last.neg==rec.neg &&
            last.fwd==rec.fwd && _brep_pc_same_carrier(last.pc,rec.pc,tol) &&
            (abs(exitp(last)-entryp(rec))<tol ||
             # seam-unwrapped canonical pieces meet end-to-start in stored
             # parameters independent of the traversal flag convention
             abs(last.t1-rec.t0)<tol || abs(rec.t1-last.t0)<tol)
        join(last,rec)=(key=last.key,fwd=last.fwd,pc=last.pc,
                        t0=min(last.t0,rec.t0),t1=max(last.t1,rec.t1),
                        crv=last.crv,vt0=min(last.vt0,rec.vt0),
                        vt1=max(last.vt1,rec.vt1),neg=last.neg)
        wires=Vector{NamedTuple}[]
        for w in s.wires
            nw=NamedTuple[]
            for ar in w
                c=get(canon,ar.key,nothing)
                rec=ar
                if c!==nothing
                    root,α,β=c
                    # canonical param u = α·t_member + β, so the merged
                    # pcurve in u composes the member pcurve with the
                    # inverse map t = (u−β)/α — negated for a neg arc,
                    # whose stored pcurve takes s = −t_member
                    γ,δ=ar.neg ? (-1/α,β/α) : (1/α,-β/α)
                    vlo,vhi=minmax(α*ar.vt0+β,α*ar.vt1+β)
                    pc2,f=_brep_pc_reparam(ar.pc,γ,δ,tol,caller)
                    # `fwd` is a traversal flag in the arc's stored
                    # parameters — the affine map from member to canonical
                    # stored params changes sign whenever the neg flags or
                    # α differ, so the traversal flag composes with it
                    fwd2=xor(ar.fwd,α<0,ar.neg,f)
                    rec= f ? (key=root,fwd=fwd2,pc=pc2,t0=-vhi,t1=-vlo,
                              crv=keyrep[root].crv,vt0=vlo,vt1=vhi,neg=true) :
                             (key=root,fwd=fwd2,pc=pc2,t0=vlo,t1=vhi,
                              crv=keyrep[root].crv,vt0=vlo,vt1=vhi,neg=false)
                end
                if !isempty(nw) && joinable(nw[end],rec)
                    nw[end]=join(nw[end],rec)
                else
                    push!(nw,rec)
                end
            end
            # the cycle may wrap: the last piece joins the first
            if length(nw)>1 && joinable(nw[end],nw[1])
                nw[1]=join(nw[end],nw[1])
                pop!(nw)
            end
            push!(wires,nw)
        end
        spec=(sg=s.sg,sgn=s.sgn,wires=wires,stag=s.stag)
        hasproperty(s,:uses) && (spec=merge(spec,(uses=s.uses,)))
        hasproperty(s,:sev) && (spec=merge(spec,(sev=s.sev,)))
        push!(out,spec)
    end
    return out
end

# ── the driver ──────────────────────────────────────────────────────────────

# Sections + trims + edge splits + imprints + per-face kept regions, for both
# operands. Returns the kept face specs: (sg, sgn, wires) — wires are
# traversal-ordered cycles of arc-uses (key,fwd,pc,t0,t1,crv,vt0,vt1,neg).
function _brep_collect_faces(m::GeoModel,op::Symbol,ta::Int,tb::Int,caller)
    tol=_BREP_TOL
    A=_brep_extract(m,ta,caller)
    B=_brep_extract(m,tb,caller)
    # coincident (same-domain) plane pairs produce no section — their wires
    # imprint into each other's domain instead
    coinc=Tuple{Int,Int}[]
    sections=NamedTuple[]
    for (ia,fa) in enumerate(A.faces),(ib,fb) in enumerate(B.faces)
        if fa.sg.occ===:plane && fb.sg.occ===:plane &&
                _brep_planes_coincident(fa.sg,fb.sg)
            push!(coinc,(ia,ib));continue
        end
        for s in _brep_face_sections(fa,fb,caller)
            # trim each section on both faces; intersect the kept span sets
            spansA=_brep_trim_section(A,fa,s.pcA,s.closed,tol;
                                      scrv=s.curve)
            spansB=_brep_trim_section(B,fb,s.pcB,s.closed,tol;
                                      scrv=s.curve)
            merged=_brep_intersect_spans(spansA,spansB,s.closed,tol,
                                         caller;pc=s.pcA)
            isempty(merged) && continue
            push!(sections,(curve=s.curve,closed=s.closed,pcA=s.pcA,
                            pcB=s.pcB,fa=ia,fb=ib,spans=merged))
        end
    end
    # operand edge splits — 3-D pierce params first; a split is a 3-D vertex
    # shared by every face owning the edge, so all are collected before any
    # arrangement is built
    esplits=(Dict{Int,Vector{Float64}}(),Dict{Int,Vector{Float64}}())
    for (sd,S,other) in ((1,A,B),(2,B,A))
        es=esplits[sd]
        for (eidx,e) in enumerate(S.edges)
            e.curve.kind===:degenerate && continue
            ps=Float64[]
            for of in other.faces
                append!(ps,_brep_edge_face_params(other,e,of,tol))
            end
            isempty(ps) ||
                (es[eidx]=_brep_edge_splits(e.curve,ps,e.t0,e.t1,tol))
        end
    end
    # coincident-face imprints: crossings and collinear contacts become edge
    # splits and section cuts
    sec_cuts=Dict{Int,Vector{Float64}}()
    imps=([NamedTuple[] for _ in A.faces],[NamedTuple[] for _ in B.faces])
    for (ia,ib) in coinc
        _brep_coinc_imprints!(A,B,ia,ib,esplits[1],esplits[2],sec_cuts,
                              imps[1][ia],sections,true,tol)
        _brep_coinc_imprints!(B,A,ib,ia,esplits[2],esplits[1],sec_cuts,
                              imps[2][ib],sections,false,tol)
    end
    # coincident section curves (an operand edge lying in a face's plane
    # gives the section of both adjacent faces) share their span endpoints as
    # cuts so the pieces align on every face carrying them
    for i in eachindex(sections), j in i+1:length(sections)
        si,sj=sections[i],sections[j]
        _brep_same_curve3d(si.curve,sj.curve,tol) || continue
        for span in si.spans, te in span, spanj in sj.spans
            p=_brep_curve_eval(si.curve,te)
            tj=_brep_curve_param(sj.curve,p,(spanj[1]+spanj[2])/2)
            (tj===nothing || tj<=spanj[1]+tol || tj>=spanj[2]-tol) &&
                continue
            push!(get!(Vector{Float64},sec_cuts,j),tj)
        end
        for span in sj.spans, te in span, spani in si.spans
            p=_brep_curve_eval(sj.curve,te)
            ti=_brep_curve_param(si.curve,p,(spani[1]+spani[2])/2)
            (ti===nothing || ti<=spani[1]+tol || ti>=spani[2]-tol) &&
                continue
            push!(get!(Vector{Float64},sec_cuts,i),ti)
        end
    end
    # sections whose curve coincides with an operand edge degenerate: the
    # result edge IS that edge — the owning faces keep their wire pieces,
    # the other operand's faces receive imprint arcs with the same key
    _brep_wire_uses_edge(wires,e)=any(w->any(x->x[1]==e,w),wires)
    esec_imps=([NamedTuple[] for _ in A.faces],[NamedTuple[] for _ in B.faces])
    esec_skip=falses(length(sections))
    for (si,sec) in enumerate(sections)
        fa=A.faces[sec.fa];fb=B.faces[sec.fb]
        eA=0;eB=0
        for wire in fa.wires,(eidx,_,_) in wire
            _brep_same_curve3d(sec.curve,A.edges[eidx].curve,tol) &&
                (eA=eidx;break)
        end
        for wire in fb.wires,(eidx,_,_) in wire
            _brep_same_curve3d(sec.curve,B.edges[eidx].curve,tol) &&
                (eB=eidx;break)
        end
        (eA==0 && eB==0) && continue
        esec_skip[si]=true
        wside,weidx = eA!=0 ? (1,eA) : (2,eB)
        lside=3-wside
        W=wside==1 ? A : B
        L=wside==1 ? B : A
        we=W.edges[weidx]
        leidx=wside==1 ? eB : eA
        # target faces of the loser operand: the section's other face, plus
        # (4-way coincidence) every loser face adjacent to the loser edge
        targets=Tuple{Int,Int,Any}[]
        if wside==1
            push!(targets,(2,sec.fb,sec.pcB))
            if eB!=0
                for (tf,tface) in enumerate(B.faces)
                    tf==sec.fb && continue
                    _brep_wire_uses_edge(tface.wires,eB) || continue
                    push!(targets,(2,tf,nothing))
                end
            end
        else
            push!(targets,(1,sec.fa,sec.pcA))
            if eA!=0
                for (tf,tface) in enumerate(A.faces)
                    tf==sec.fa && continue
                    _brep_wire_uses_edge(tface.wires,eA) || continue
                    push!(targets,(1,tf,nothing))
                end
            end
        end
        for (lo,hi) in sec.spans
            subs=_brep_esec_tspans(we,sec.curve,lo,hi,tol)
            subs===nothing && throw(ArgumentError(
                "$caller: a coincident section endpoint does not project " *
                "onto its operand edge"))
            for (sa,sb,ta,tb) in subs
                p1=_brep_curve_eval(sec.curve,sa)
                p2=_brep_curve_eval(sec.curve,sb)
                wlo,whi=minmax(ta,tb)
                push!(get!(Vector{Float64},esplits[wside],weidx),wlo,whi)
                if leidx!=0
                    le=L.edges[leidx]
                    lmid=(le.t0+le.t1)/2
                    l1=_brep_curve_param(le.curve,p1,lmid)
                    l2=_brep_curve_param(le.curve,p2,lmid)
                    l1===nothing || l2===nothing ||
                        push!(get!(Vector{Float64},esplits[lside],leidx),
                              min(l1,l2),max(l1,l2))
                end
                for (ts,tf,spc) in targets
                    tface=L.faces[tf]
                    pc,neg = if spc!==nothing
                        α=(sb-sa)/(tb-ta);β=sa-α*ta
                        _brep_pc_reparam(spc,α,β,tol,caller)
                    else
                        _brep_pc_on_surface(tface.sg,we.curve,wlo,whi,
                                            tol,caller)
                    end
                    sslo,sshi=neg ? (-whi,-wlo) : (wlo,whi)
                    _brep_imprint_contacts!(L,W,tf,weidx,pc,neg,sslo,sshi,
                                            esplits[ts],esplits[wside],
                                            sec_cuts,sections,ts==1,tol,si)
                    # span endpoints landing on the target face's wires
                    # split those edges too
                    for te in (sslo,sshi)
                        uv=_brep_pc_eval(pc,te)
                        for wire in tface.wires,(e2,_,epc2) in wire
                            e2r=L.edges[e2]
                            lo2,hi2=minmax(e2r.t0,e2r.t1)
                            tp=_brep_pc_param(epc2,uv,lo2,hi2,10tol)
                            (tp===nothing || tp<=lo2+tol || tp>=hi2-tol) &&
                                continue
                            push!(get!(Vector{Float64},esplits[ts],e2),tp)
                        end
                    end
                    # the same 3-D edge can be reached through several
                    # section records (each face pair of the 4-way
                    # coincidence emits one) — collect the spans under one
                    # spec rather than duplicating the carrier. A seam-
                    # wrapped span produces one spec per canonical piece:
                    # their pcurves differ by a period shift and must not
                    # share a record.
                    ex_i=findfirst(es->es.weidx==weidx && es.neg==neg &&
                                       es.pc==pc,esec_imps[ts][tf])
                    if ex_i===nothing
                        push!(esec_imps[ts][tf],
                              (pc=pc,spans=[(wlo,whi)],neg=neg,
                               wside=wside,weidx=weidx,crv=we.curve))
                    else
                        es=esec_imps[ts][tf][ex_i]
                        any(sp->abs(sp[1]-wlo)<tol && abs(sp[2]-whi)<tol,
                            es.spans) || push!(es.spans,(wlo,whi))
                    end
                end
            end
        end
    end
    # section sub-edges, cut at both faces' cut-domain seams (a seam crossing
    # is a real result vertex shared by both faces' arcs, so the same
    # sub-edge keys must appear on each side) and at coincident contacts
    sec_edges=NamedTuple[]
    secs_A=[NamedTuple[] for _ in A.faces]
    secs_B=[NamedTuple[] for _ in B.faces]
    # fitted-section endpoint anchors — see `_brep_sec_vertex`
    mscale=1.0
    for S in (A,B), p in S.pts
        mscale=max(mscale,abs(p[1]),abs(p[2]),abs(p[3]))
    end
    for (si,sec) in enumerate(sections)
        esec_skip[si] && continue
        for (lo,hi) in sec.spans
            va=_brep_sec_vertex((A,B),1,2,esplits[1],esplits[2],
                                sec,lo,tol,mscale)
            vb=_brep_sec_vertex((A,B),1,2,esplits[1],esplits[2],
                                sec,hi,tol,mscale)
            for a in (va,vb)
                a===nothing && continue
                a.A!==nothing && push!(
                    get!(Vector{Float64},esplits[1],a.A[2]),a.A[3])
                a.B!==nothing && push!(
                    get!(Vector{Float64},esplits[2],a.B[2]),a.B[3])
            end
            cuts=copy(get(sec_cuts,si,Float64[]))
            for (S,face,spc) in ((A,A.faces[sec.fa],sec.pcA),
                                 (B,B.faces[sec.fb],sec.pcB))
                uper=_brep_uperiod(face.sg)
                uper>0 || continue
                ulo=_brep_face_ubase(S.edges,face)
                pieces=_brep_pc_split_all_useams(spc,lo,hi,uper,ulo,tol)
                for p in pieces[2:end]
                    push!(cuts,p[1])
                end
            end
            for (s0,s1) in _brep_split_ranges(lo,hi,cuts,tol)
                a0=s0-lo<tol ? va : nothing
                a1=s1-hi>-tol ? vb : nothing
                push!(sec_edges,(curve=sec.curve,t0=s0,t1=s1,
                                 v0=_brep_sec_vfix(a0),
                                 v1=_brep_sec_vfix(a1)))
                se=length(sec_edges)
                fA=_brep_anchor_uv(A.faces[sec.fa],1,a0)
                fA1=_brep_anchor_uv(A.faces[sec.fa],1,a1)
                fB=_brep_anchor_uv(B.faces[sec.fb],2,a0)
                fB1=_brep_anchor_uv(B.faces[sec.fb],2,a1)
                if get(ENV,"BREP_DEBUG_ARCS","")!=""
                    println("DBG se$se span=($(round(s0,digits=3))," *
                            "$(round(s1,digits=3))) fa=$(sec.fa) fb=$(sec.fb)" *
                            " fixA0=$(fA===nothing ? nothing : round.(fA,digits=3))" *
                            " fixA1=$(fA1===nothing ? nothing : round.(fA1,digits=3))" *
                            " fixB0=$(fB===nothing ? nothing : round.(fB,digits=3))" *
                            " fixB1=$(fB1===nothing ? nothing : round.(fB1,digits=3))")
                end
                push!(secs_A[sec.fa],(se=se,span=(s0,s1),pc=sec.pcA,
                                      crv=sec.curve,anc0=a0,anc1=a1,
                                      fix0=fA,fix1=fA1))
                push!(secs_B[sec.fb],(se=se,span=(s0,s1),pc=sec.pcB,
                                      crv=sec.curve,anc0=a0,anc1=a1,
                                      fix0=fB,fix1=fB1))
            end
        end
    end
    # section endpoints landing on a wedge split that edge globally
    for (sd,S,secs) in ((1,A,secs_A),(2,B,secs_B))
        es=esplits[sd]
        for (fidx,f) in enumerate(S.faces), sf in secs[fidx]
            uper=_brep_uperiod(f.sg)
            ulo=_brep_face_ubase(S.edges,f)
            for (tend,anc) in ((sf.span[1],get(sf,:anc0,nothing)),
                               (sf.span[2],get(sf,:anc1,nothing)))
                # an endpoint anchored on this operand's own wire edge
                # already pushed its authoritative pierce param —
                # re-projecting the fitted uv would add a residual-off
                # duplicate split
                _brep_anchor_own(sd,anc) && continue
                uv=_brep_wrap_u(_brep_pc_eval(sf.pc,tend),ulo,uper)
                for wire in f.wires,(eidx,sign,epc) in wire
                    e=S.edges[eidx]
                    lo,hi=minmax(e.t0,e.t1)
                    t=_brep_pc_param(epc,uv,lo,hi,10tol)
                    t===nothing && continue
                    (t>lo+tol && t<hi-tol) || continue
                    push!(get!(Vector{Float64},es,eidx),t)
                end
            end
        end
    end
    # per-face coincident partner lists for the region keep-table
    coinc_by_face=(Dict{Int,Vector{Int}}(),Dict{Int,Vector{Int}}())
    for (ia,ib) in coinc
        push!(get!(Vector{Int},coinc_by_face[1],ia),ib)
        push!(get!(Vector{Int},coinc_by_face[2],ib),ia)
    end
    kept=NamedTuple[]
    for (side,S,other,secs) in ((1,A,B,secs_A),(2,B,A,secs_B))
        es=esplits[side];oes=esplits[3-side]
        # imprint arcs are built only after every split settles — the pieces
        # are indexed by the other operand's own split table; degenerate-
        # section arcs come first so coincident-face imprints covered by
        # them drop instead of duplicating a carrier
        extras=[begin
                    ea=_brep_esec_arcs(other,S,S.faces[fidx],
                                     esec_imps[side][fidx],oes,3-side,tol)
                    append!(_brep_imprint_arcs(S,other,S.faces[fidx],
                                               imps[side][fidx],oes,
                                               secs[fidx],3-side,ea,tol),ea)
                end for fidx in eachindex(S.faces)]
        for (fidx,f) in enumerate(S.faces)
            regions,arcs,cycles=_brep_face_regions(
                S,f,side,es,secs[fidx],extras[fidx],tol,caller)
            cofaces=get(coinc_by_face[side],fidx,Int[])
            for region in regions
                uv=_brep_region_sample(arcs,cycles,region,tol,caller)
                P=_brep_surface_eval(f.sg,uv[1],uv[2])
                # a region lying inside a coincident partner face sits on the
                # other solid's boundary — containment cannot classify it;
                # the op's same-domain rule resolves by outward-normal sign
                ci=0
                for oi in cofaces
                    of=other.faces[oi]
                    uv2=_brep_surface_uv(of.sg,P)
                    _brep_uv_in_face(other.edges,of,uv2,tol) || continue
                    ci=oi;break
                end
                keep,flip=if ci>0
                    of=other.faces[ci]
                    same=(f.orev*of.orev*_b3dot(f.sg.axis,of.sg.axis))>0
                    # like-oriented coincident faces share material on one
                    # side: union/intersection keep the A copy (the B copy is
                    # dropped, deduplicating), difference drops both;
                    # opposite-oriented faces sandwich material — difference
                    # keeps the A copy, union/intersection drop both
                    (side==1 ? (op===:difference ? !same : same) : false),
                    false
                else
                    inside=_brep_point_in_solid(other,P,tol,caller)
                    if side==1
                        (op===:intersection ? inside : !inside),false
                    else
                        if op===:difference
                            inside,true
                        else
                            (op===:intersection ? inside : !inside),false
                        end
                    end
                end
                if get(ENV,"BREP_DEBUG_ARCS","")!=""
                    println("DBG   side=$side fidx=$fidx region=$region " *
                            "uv=$(round.(uv,digits=3)) P=$(round.(P,digits=3)) " *
                            "inside=$(ci>0 ? "coinc" : _brep_point_in_solid(other,P,tol,caller)) " *
                            "keep=$keep flip=$flip")
                end
                keep || continue
                wires=[ [ (key=arcs[i].key,fwd=fwd,pc=arcs[i].pc,
                           t0=arcs[i].t0,t1=arcs[i].t1,crv=arcs[i].crv,
                           vt0=arcs[i].vt0,vt1=arcs[i].vt1,neg=arcs[i].neg)
                         for (i,fwd) in cycles[c] ] for c in region ]
                sgn=flip ? -f.orev : f.orev
                spec=(sg=f.sg,sgn=sgn,wires=wires,stag=f.stag)
                isempty(sec_edges) ||
                    (spec=merge(spec,(sev=sec_edges,)))
                push!(kept,spec)
            end
        end
    end
    # fuse/intersection glue same-domain pieces on one carrier into single
    # faces; a difference keeps its imprint splits
    (op===:union || op===:intersection) &&
        (kept=_brep_unify_union(kept,tol,caller))
    return kept
end

# ── multi-operand arrangement ───────────────────────────────────────────────
#
# N-way generalization of the pairwise pipeline for `.geo` Boolean lists —
# Gmsh's OCC `booleanOperator` over object/tool operand vectors. Every operand
# pair contributes sections and coincident-face imprints exactly like the
# binary pass; each resulting face region is then classified by the SET of
# operands containing its inside and outside sides (a membership bitmask over
# the operand order, objects in the low bits). The kept cells of that
# arrangement each materialize into one result volume, matching OCC's
# BOPAlgo/BuilderAlgo cell decomposition: union fuses all nonempty cells into
# connected solids, difference keeps object cells clear of every tool,
# intersection keeps cells inside at least one object and one tool, and
# fragments keeps every nonempty cell with the partition walls shared.

# Operand-membership masks are UInt64 bit sets over the operand order.
_brep_bit(i::Int) = UInt64(1) << (i-1)

# Is a cell with membership `label` part of the result? Objects occupy the
# low `nobj` bits.
function _brep_label_keep(op::Symbol,label::UInt64,nobj::Int)
    label==0 && return false
    (op===:union || op===:fragments) && return true
    objmask=(UInt64(1)<<nobj)-UInt64(1)
    if op===:difference
        return (label & objmask)!=0 && (label & ~objmask)==0
    end
    return (label & objmask)!=0 && (label & ~objmask)!=0
end

# `_brep_imprint_contacts!` for the N-way sections table: a section matches
# the imprinted face when either endpoint operand is the face's operand.
function _brep_imprint_contacts_n!(solids,si,fidx,oi,eOidx,mpc,neg,slo,shi,
                                   esplits,sec_cuts,sections,tol,
                                   skip_sec::Int=0)
    S=solids[si];O=solids[oi]
    es_S=esplits[si];es_O=esplits[oi]
    f=S.faces[fidx]
    s2t(s)=neg ? -s : s
    function inrange(s)
        if _brep_pc_iscirc(mpc)
            return _brep_snap_theta(s,slo,shi,tol)
        end
        return (s>=slo-tol && s<=shi+tol) ? s : nothing
    end
    for wire2 in f.wires,(eidx2,_,epc2) in wire2
        e2=S.edges[eidx2]
        e2.curve.kind===:degenerate && continue
        lo2,hi2=minmax(e2.t0,e2.t1)
        for t in _brep_pc_intersect(epc2,mpc,slo,shi,tol)
            push!(get!(Vector{Float64},es_S,eidx2),t)
        end
        for s in _brep_pc_intersect(mpc,epc2,lo2,hi2,tol)
            ss=inrange(s)
            ss===nothing ||
                push!(get!(Vector{Float64},es_O,eOidx),s2t(ss))
        end
        _brep_pc_same_carrier(epc2,mpc,tol) || continue
        for te in (e2.t0,e2.t1)
            uv=_brep_pc_eval(epc2,te)
            sp=_brep_pc_param(mpc,uv,slo,shi,tol)
            sp===nothing ||
                push!(get!(Vector{Float64},es_O,eOidx),s2t(sp))
        end
        for se in (slo,shi)
            uv=_brep_pc_eval(mpc,se)
            tp=_brep_pc_param(epc2,uv,lo2,hi2,tol)
            tp===nothing ||
                push!(get!(Vector{Float64},es_S,eidx2),tp)
        end
    end
    for (nsi,sec) in enumerate(sections)
        nsi==skip_sec && continue
        local spc
        if sec.oi==si && sec.fa==fidx
            spc=sec.pcA
        elseif sec.oj==si && sec.fb==fidx
            spc=sec.pcB
        else
            continue
        end
        nlifts=_brep_pc_iscirc(spc) ? (-1,0,1) : (0,)
        for t in _brep_pc_intersect(spc,mpc,slo,shi,tol), k in nlifts
            push!(get!(Vector{Float64},sec_cuts,nsi),t+k*_OCC_TWO_PI)
        end
        for span in sec.spans
            for s in _brep_pc_intersect(mpc,spc,span[1],span[2],tol)
                ss=inrange(s)
                ss===nothing ||
                    push!(get!(Vector{Float64},es_O,eOidx),s2t(ss))
            end
            _brep_pc_same_carrier(spc,mpc,tol) || continue
            for te in span
                uv=_brep_pc_eval(spc,te)
                sp=_brep_pc_param(mpc,uv,slo,shi,tol)
                sp===nothing ||
                    push!(get!(Vector{Float64},es_O,eOidx),s2t(sp))
            end
            for se in (slo,shi)
                uv=_brep_pc_eval(mpc,se)
                tp=_brep_pc_param(spc,uv,span[1],span[2],tol)
                tp===nothing ||
                    push!(get!(Vector{Float64},sec_cuts,nsi),tp)
            end
        end
    end
    return nothing
end

# One direction of a coincident pair: operand `oj`'s face `fj` is imprinted
# onto operand `si`'s face `fi`.
function _brep_coinc_imprints_n!(solids,si,fi,oj,fj,esplits,sec_cuts,imps_out,
                                 sections,tol)
    S=solids[si];O=solids[oj]
    f=S.faces[fi];of=O.faces[fj]
    mp=_brep_plane_map(f.sg,of.sg)
    seen=Set{Int}()
    for wire in of.wires,(eidx,sign,epc) in wire
        eidx in seen && continue
        push!(seen,eidx)
        e=O.edges[eidx]
        e.curve.kind===:degenerate && continue
        mpc,neg=_brep_map_pc(epc,mp)
        elo,ehi=minmax(e.t0,e.t1)
        slo,shi=neg ? (-ehi,-elo) : (elo,ehi)
        _brep_imprint_contacts_n!(solids,si,fi,oj,eidx,mpc,neg,slo,shi,
                                  esplits,sec_cuts,sections,tol)
        push!(imps_out,(eidx=eidx,pc=mpc,neg=neg,crv=e.curve))
    end
    return nothing
end

# Collect the kept face regions of an N-way Boolean. `tags` lists the operand
# volumes in order — objects first, then tools; `nobj` is the object count.
# Returns `(kept, touched)` where each kept spec carries
# `uses::Vector{Tuple{UInt64,Int}}` — the adjacent result cells (membership
# mask) plus the face sign inside that cell's shell — and `touched` marks
# operands participating in any geometric interaction (needed for OCC's
# modified-vs-unmodified tag semantics downstream).
function _brep_collect_faces_n(m::GeoModel,op::Symbol,tags::Vector{Int},
                               nobj::Int,caller)
    tol=_BREP_TOL
    N=length(tags)
    N<=62 || throw(ArgumentError(
        "$caller: at most 62 Boolean operands per operation"))
    solids=[_brep_extract(m,t,caller) for t in tags]
    coinc=NTuple{4,Int}[]
    sections=NamedTuple[]
    touched=falses(N)
    for i in 1:N, j in i+1:N
        Si=solids[i];Sj=solids[j]
        for (fi,fI) in enumerate(Si.faces), (fj,fJ) in enumerate(Sj.faces)
            if fI.sg.occ===:plane && fJ.sg.occ===:plane &&
                    _brep_planes_coincident(fI.sg,fJ.sg)
                push!(coinc,(i,fi,j,fj))
                continue
            end
            for s in _brep_face_sections(fI,fJ,caller)
                spansI=_brep_trim_section(Si,fI,s.pcA,s.closed,tol;
                                          scrv=s.curve)
                spansJ=_brep_trim_section(Sj,fJ,s.pcB,s.closed,tol;
                                          scrv=s.curve)
                merged=_brep_intersect_spans(spansI,spansJ,s.closed,tol,
                                             caller;pc=s.pcA)
                isempty(merged) && continue
                push!(sections,(curve=s.curve,closed=s.closed,pcA=s.pcA,
                                pcB=s.pcB,oi=i,fa=fi,oj=j,fb=fj,
                                spans=merged))
                touched[i]=touched[j]=true
            end
        end
    end
    # operand edges piercing other operands' faces split those edges
    esplits=[Dict{Int,Vector{Float64}}() for _ in 1:N]
    for i in 1:N
        Si=solids[i];es=esplits[i]
        for (eidx,e) in enumerate(Si.edges)
            e.curve.kind===:degenerate && continue
            ps=Float64[]
            for j in 1:N
                j==i && continue
                for of in solids[j].faces
                    append!(ps,_brep_edge_face_params(solids[j],e,of,tol))
                end
            end
            if !isempty(ps)
                es[eidx]=sort!(ps);touched[i]=true
            end
        end
    end
    # coincident face pairs imprint each partner's boundary onto the other
    sec_cuts=Dict{Int,Vector{Float64}}()
    imps=[[Dict{Int,Vector{NamedTuple}}() for _ in solids[i].faces]
          for i in 1:N]
    for (i,fi,j,fj) in coinc
        _brep_coinc_imprints_n!(solids,i,fi,j,fj,esplits,sec_cuts,
            get!(Vector{NamedTuple},imps[i][fi],j),sections,tol)
        _brep_coinc_imprints_n!(solids,j,fj,i,fi,esplits,sec_cuts,
            get!(Vector{NamedTuple},imps[j][fj],i),sections,tol)
    end
    # sections whose curve coincides with an operand edge degenerate: the
    # lowest-indexed owning operand's edge becomes the result edge and the
    # other operands' faces receive imprint arcs carrying its key
    esec_imps=[[NamedTuple[] for _ in solids[i].faces] for i in 1:N]
    esec_skip=falses(length(sections))
    for (si,sec) in enumerate(sections)
        owners=Int[];owned=Int[]
        for k in (sec.oi,sec.oj)
            S=solids[k]
            efound=0
            for wire in S.faces[k==sec.oi ? sec.fa : sec.fb].wires,
                (eidx,_,_) in wire
                if _brep_same_curve3d(sec.curve,S.edges[eidx].curve,tol)
                    efound=eidx;break
                end
            end
            efound==0 || (push!(owners,k);push!(owned,efound))
        end
        isempty(owners) && continue
        esec_skip[si]=true
        # a third operand's edge can coincide with the section line over the
        # same span — it is the same result edge, so it joins the owners and
        # the globally lowest-indexed owner wins; without this both it and
        # the section's own winner would imprint duplicate arcs
        for k in 1:N
            (k==sec.oi || k==sec.oj) && continue
            S=solids[k];efound=0
            for tface in S.faces, wire in tface.wires,(eidx,_,_) in wire
                e=S.edges[eidx]
                _brep_same_curve3d(sec.curve,e.curve,tol) || continue
                covers=any(sec.spans) do (lo,hi)
                    pm=_brep_curve_eval(sec.curve,(lo+hi)/2)
                    tp=_brep_curve_param(e.curve,pm,(e.t0+e.t1)/2)
                    tp!==nothing && tp>=min(e.t0,e.t1)-tol &&
                        tp<=max(e.t0,e.t1)+tol
                end
                covers && (efound=eidx;break)
            end
            efound==0 || (push!(owners,k);push!(owned,efound))
        end
        wpos=argmin(owners);w=owners[wpos];weidx=owned[wpos]
        we=solids[w].edges[weidx]
        targets=Tuple{Int,Int,Any}[]
        for k in (sec.oi,sec.oj)
            k==w && continue
            kf=k==sec.oi ? sec.fa : sec.fb
            kpc=k==sec.oi ? sec.pcA : sec.pcB
            push!(targets,(k,kf,kpc))
            kp=findfirst(==(k),owners)
            if kp!==nothing
                ke=owned[kp]
                for (tf,tface) in enumerate(solids[k].faces)
                    tf==kf && continue
                    any(wr->any(x->x[1]==ke,wr),tface.wires) || continue
                    push!(targets,(k,tf,nothing))
                end
            end
        end
        # non-endpoint owners are never the section face; every face of a
        # losing non-endpoint owner carrying the edge needs the imprint
        for (k2,oe) in zip(owners,owned)
            (k2==w || k2==sec.oi || k2==sec.oj) && continue
            for (tf,tface) in enumerate(solids[k2].faces)
                any(wr->any(x->x[1]==oe,wr),tface.wires) || continue
                push!(targets,(k2,tf,nothing))
            end
        end
        for (lo,hi) in sec.spans
            subs=_brep_esec_tspans(we,sec.curve,lo,hi,tol)
            subs===nothing && throw(ArgumentError(
                "$caller: a coincident section endpoint does not project " *
                "onto its operand edge"))
            for (sa,sb,ta,tb) in subs
                p1=_brep_curve_eval(sec.curve,sa)
                p2=_brep_curve_eval(sec.curve,sb)
                wlo,whi=minmax(ta,tb)
                push!(get!(Vector{Float64},esplits[w],weidx),wlo,whi)
                touched[w]=true
                for (k2,oe) in zip(owners,owned)
                    k2==w && continue
                    oer=solids[k2].edges[oe]
                    omid=(oer.t0+oer.t1)/2
                    l1=_brep_curve_param(oer.curve,p1,omid)
                    l2=_brep_curve_param(oer.curve,p2,omid)
                    l1===nothing || l2===nothing || (push!(get!(
                        Vector{Float64},esplits[k2],oe),min(l1,l2),
                        max(l1,l2));touched[k2]=true)
                end
                for (tk,tf,spc) in targets
                    tface=solids[tk].faces[tf]
                    pc,neg = if spc!==nothing
                        α=(sb-sa)/(tb-ta);β=sa-α*ta
                        _brep_pc_reparam(spc,α,β,tol,caller)
                    else
                        _brep_pc_on_surface(tface.sg,we.curve,wlo,whi,
                                            tol,caller)
                    end
                    sslo,sshi=neg ? (-whi,-wlo) : (wlo,whi)
                    _brep_imprint_contacts_n!(solids,tk,tf,w,weidx,pc,neg,
                                              sslo,sshi,esplits,sec_cuts,
                                              sections,tol,si)
                    for te in (sslo,sshi)
                        uv=_brep_pc_eval(pc,te)
                        for wire in tface.wires,(e2,_,epc2) in wire
                            e2r=solids[tk].edges[e2]
                            lo2,hi2=minmax(e2r.t0,e2r.t1)
                            tp=_brep_pc_param(epc2,uv,lo2,hi2,10tol)
                            (tp===nothing || tp<=lo2+tol || tp>=hi2-tol) &&
                                continue
                            push!(get!(Vector{Float64},esplits[tk],e2),tp)
                            touched[tk]=true
                        end
                    end
                    ex_i=findfirst(es->es.wside==w && es.weidx==weidx &&
                                       es.neg==neg && es.pc==pc,
                                   esec_imps[tk][tf])
                    if ex_i===nothing
                        push!(esec_imps[tk][tf],
                              (pc=pc,spans=[(wlo,whi)],neg=neg,
                               wside=w,weidx=weidx,crv=we.curve))
                    else
                        es=esec_imps[tk][tf][ex_i]
                        any(sp->abs(sp[1]-wlo)<tol && abs(sp[2]-whi)<tol,
                            es.spans) || push!(es.spans,(wlo,whi))
                    end
                end
            end
        end
    end
    # section sub-edges, cut at both faces' seam lifts and contact cuts
    sec_edges=NamedTuple[]
    secs_pf=[[NamedTuple[] for _ in solids[i].faces] for i in 1:N]
    # fitted-section endpoint anchors — see `_brep_sec_vertex`
    mscale=1.0
    for S in solids, p in S.pts
        mscale=max(mscale,abs(p[1]),abs(p[2]),abs(p[3]))
    end
    for (si,sec) in enumerate(sections)
        esec_skip[si] && continue
        for (lo,hi) in sec.spans
            va=_brep_sec_vertex(solids,sec.oi,sec.oj,
                                esplits[sec.oi],esplits[sec.oj],
                                sec,lo,tol,mscale)
            vb=_brep_sec_vertex(solids,sec.oi,sec.oj,
                                esplits[sec.oi],esplits[sec.oj],
                                sec,hi,tol,mscale)
            for a in (va,vb)
                a===nothing && continue
                a.A!==nothing && push!(get!(Vector{Float64},
                    esplits[sec.oi],a.A[2]),a.A[3])
                a.B!==nothing && push!(get!(Vector{Float64},
                    esplits[sec.oj],a.B[2]),a.B[3])
            end
            cuts=copy(get(sec_cuts,si,Float64[]))
            for k in (sec.oi,sec.oj)
                S=solids[k]
                face=k==sec.oi ? S.faces[sec.fa] : S.faces[sec.fb]
                spc=k==sec.oi ? sec.pcA : sec.pcB
                uper=_brep_uperiod(face.sg)
                uper>0 || continue
                ulo=_brep_face_ubase(S.edges,face)
                pieces=_brep_pc_split_all_useams(spc,lo,hi,uper,ulo,tol)
                for p in pieces[2:end]
                    push!(cuts,p[1])
                end
            end
            for (s0,s1) in _brep_split_ranges(lo,hi,cuts,tol)
                a0=s0-lo<tol ? va : nothing
                a1=s1-hi>-tol ? vb : nothing
                push!(sec_edges,(curve=sec.curve,t0=s0,t1=s1,
                                 v0=_brep_sec_vfix(a0),
                                 v1=_brep_sec_vfix(a1)))
                se=length(sec_edges)
                fa=solids[sec.oi].faces[sec.fa]
                fb=solids[sec.oj].faces[sec.fb]
                push!(secs_pf[sec.oi][sec.fa],
                      (se=se,span=(s0,s1),pc=sec.pcA,crv=sec.curve,
                       anc0=a0,anc1=a1,myside=1,
                       fix0=_brep_anchor_uv(fa,1,a0),
                       fix1=_brep_anchor_uv(fa,1,a1)))
                push!(secs_pf[sec.oj][sec.fb],
                      (se=se,span=(s0,s1),pc=sec.pcB,crv=sec.curve,
                       anc0=a0,anc1=a1,myside=2,
                       fix0=_brep_anchor_uv(fb,2,a0),
                       fix1=_brep_anchor_uv(fb,2,a1)))
            end
        end
    end
    # section endpoints landing on a wire edge split that edge globally
    for i in 1:N
        es=esplits[i];S=solids[i]
        for (fidx,f) in enumerate(S.faces), sf in secs_pf[i][fidx]
            uper=_brep_uperiod(f.sg)
            ulo=_brep_face_ubase(S.edges,f)
            for (tend,anc) in ((sf.span[1],get(sf,:anc0,nothing)),
                               (sf.span[2],get(sf,:anc1,nothing)))
                # an endpoint anchored on this operand's own wire edge
                # already pushed its pierce param — re-projecting the
                # fitted uv would add a residual-off duplicate split
                _brep_anchor_own(get(sf,:myside,0),anc) && continue
                uv=_brep_wrap_u(_brep_pc_eval(sf.pc,tend),ulo,uper)
                for wire in f.wires,(eidx,sign,epc) in wire
                    e=S.edges[eidx]
                    lo,hi=minmax(e.t0,e.t1)
                    t=_brep_pc_param(epc,uv,lo,hi,10tol)
                    t===nothing && continue
                    (t>lo+tol && t<hi-tol) || continue
                    push!(get!(Vector{Float64},es,eidx),t)
                    touched[i]=true
                end
            end
        end
    end
    # per-face coincident partner lists for the region keep-table
    coinc_by_face=[Dict{Int,Vector{Tuple{Int,Int}}}() for _ in 1:N]
    for (i,fi,j,fj) in coinc
        push!(get!(Vector{Tuple{Int,Int}},coinc_by_face[i],fi),(j,fj))
        push!(get!(Vector{Tuple{Int,Int}},coinc_by_face[j],fj),(i,fi))
    end
    kept=NamedTuple[]
    for i in 1:N
        S=solids[i];es=esplits[i]
        # imprint arcs settle only after every split does — the pieces are
        # indexed by the partner operands' own split tables; degenerate-
        # section arcs come first so covered imprints drop instead of
        # duplicating a carrier
        extras=[begin
                    f=S.faces[fidx]
                    ea=NamedTuple[]
                    by_winner=Dict{Int,Vector{NamedTuple}}()
                    for spec in esec_imps[i][fidx]
                        push!(get!(Vector{NamedTuple},by_winner,spec.wside),
                              spec)
                    end
                    for w in sort!(collect(keys(by_winner)))
                        append!(ea,_brep_esec_arcs(solids[w],S,f,
                                                 by_winner[w],esplits[w],
                                                 w,tol))
                    end
                    arcs=copy(ea)
                    for j in sort!(collect(keys(imps[i][fidx])))
                        append!(arcs,_brep_imprint_arcs(
                            S,solids[j],f,imps[i][fidx][j],esplits[j],
                            secs_pf[i][fidx],j,arcs,tol))
                    end
                    arcs
                end for fidx in eachindex(S.faces)]
        for (fidx,f) in enumerate(S.faces)
            # an operand is geometrically modified when a face carries
            # interior imprint/section arcs or loses a region — a coincident
            # plane without overlapping imprint does not modify it
            isempty(extras[fidx]) || (touched[i]=true)
            regions,arcs,cycles=_brep_face_regions(
                S,f,i,es,secs_pf[i][fidx],extras[fidx],tol,caller)
            cofaces=get(coinc_by_face[i],fidx,Tuple{Int,Int}[])
            for region in regions
                uv=_brep_region_sample(arcs,cycles,region,tol,caller)
                P=_brep_surface_eval(f.sg,uv[1],uv[2])
                if get(ENV,"BREP_DEBUG_ARCS","")!=""
                    println("DBG   region i=$i fidx=$fidx stag=" *
                            "$(get(f,:stag,0)) uv=$(round.(uv,digits=3)) " *
                            "P=$(round.(P,digits=3)) area=",
                            round(sum(abs,_brep_cycle_area(arcs,cycles[c])
                                     for c in region),sigdigits=6))
                end
                # a region sitting on coincident partner faces is on those
                # solids' boundary — point-in-solid cannot classify them; the
                # partners join the membership labels by normal direction
                same_p=Int[];opp_p=Int[]
                for (oj,ofidx) in cofaces
                    of=solids[oj].faces[ofidx]
                    uv2=_brep_surface_uv(of.sg,P)
                    _brep_uv_in_face(solids[oj].edges,of,uv2,tol) || continue
                    sameside=(f.orev*of.orev*_b3dot(f.sg.axis,of.sg.axis))>0
                    push!(sameside ? same_p : opp_p,oj)
                end
                coincident=!isempty(same_p)||!isempty(opp_p)
                lin=_brep_bit(i);lout=UInt64(0)
                for oj in same_p
                    lin|=_brep_bit(oj)
                end
                for oj in opp_p
                    lout|=_brep_bit(oj)
                end
                for k in 1:N
                    k==i && continue
                    (k in same_p || k in opp_p) && continue
                    _brep_point_in_solid(solids[k],P,tol,caller) || continue
                    b=_brep_bit(k);lin|=b;lout|=b
                end
                keep_in=_brep_label_keep(op,lin,nobj)
                keep_out=_brep_label_keep(op,lout,nobj)
                uses=Tuple{UInt64,Int}[]
                if coincident
                    # the lowest-indexed coincident copy carries the wall —
                    # the same computation on that operand's own face emits
                    # the identical uses, so non-owners emit nothing
                    if minimum(vcat(i,same_p,opp_p))!=i
                        touched[i]=true;continue
                    end
                    if !isempty(opp_p) &&
                            (op===:union || op===:intersection)
                        # opposite-normal copies sandwich distinct material —
                        # under fuse/common that seam is internal and drops
                        touched[i]=true;continue
                    end
                    keep_in && push!(uses,(lin,f.orev))
                    keep_out && push!(uses,(lout,-f.orev))
                else
                    keep_in && push!(uses,(lin,f.orev))
                    keep_out && push!(uses,(lout,-f.orev))
                end
                # a wall between two kept cells is internal to the fuse — it
                # carries no result boundary. The cell ops keep both uses:
                # the wall separates two different result pieces
                if op===:union && length(uses)==2
                    empty!(uses)
                end
                isempty(uses) && (touched[i]=true;continue)
                wires=[ [ (key=arcs[k].key,fwd=fwd,pc=arcs[k].pc,
                           t0=arcs[k].t0,t1=arcs[k].t1,crv=arcs[k].crv,
                           vt0=arcs[k].vt0,vt1=arcs[k].vt1,neg=arcs[k].neg)
                         for (k,fwd) in cycles[c] ] for c in region ]
                spec=(sg=f.sg,stag=f.stag,wires=wires,uses=uses,
                      sgn=uses[1][2])
                isempty(sec_edges) ||
                    (spec=merge(spec,(sev=sec_edges,)))
                push!(kept,spec)
            end
        end
    end
    # fuse glues same-domain pieces on one carrier into single faces;
    # fragments/difference keep their imprint splits. Common glues only
    # regions serving the same cell set — a cross-cell merge would erase the
    # shared wall edges that bound each piece's shell
    op===:union && (kept=_brep_unify_union(kept,tol,caller))
    op===:intersection &&
        (kept=_brep_unify_union(kept,tol,caller;same_uses=true))
    return kept,touched
end

# ── result materialization ──────────────────────────────────────────────────

function _brep_vertex_eval(crv,t)
    crv.kind===:degenerate && return crv.xyz
    return _brep_curve_eval(crv,t)
end

# Solid-shaped view of a kept-face subset for `_brep_point_in_solid` — the
# containment tests that group result shells into volumes.
function _brep_spec_solid_view(kept,fis)
    faces=NamedTuple[];vedges=NamedTuple[]
    for fi in fis
        spec=kept[fi]
        wires=Vector{Tuple{Int,Int,NamedTuple}}[]
        for wire in spec.wires
            w=Tuple{Int,Int,NamedTuple}[]
            for ar in wire
                # a per-arc edge view — the arc's t0/t1 are the pcurve's own
                # parameter range, which a mirrored imprint negates
                push!(vedges,(curve=ar.crv,t0=ar.t0,t1=ar.t1))
                push!(w,(length(vedges),ar.fwd ? 1 : -1,ar.pc))
            end
            push!(wires,w)
        end
        push!(faces,(sg=spec.sg,stag=spec.stag,orev=1,wires=wires))
    end
    return (edges=vedges,faces=faces)
end

# One boundary vertex of a kept-face subset — the representative point used
# for shell containment tests.
_brep_spec_rep(kept,fis) =
    _brep_vertex_eval(kept[fis[1]].wires[1][1].crv,
                      kept[fis[1]].wires[1][1].vt0)

# Containment forms a forest: a component's parent is the component containing
# it one nesting level up. Even-depth components are disjoint solids — each
# becomes its own result volume (Gmsh's OCC Boolean returns one volume per
# solid); odd-depth components are void shells attached to their parent
# volume. Returns the volume groups (each a list of component indices, first
# entry the outer shell).
function _brep_volume_groups(inside,inside_count,caller)
    ncomp=length(inside_count)
    parent=zeros(Int,ncomp)
    for i in 1:ncomp
        inside_count[i]==0 && continue
        for j in 1:ncomp
            if inside[j,i] && inside_count[j]==inside_count[i]-1
                parent[i]=j;break
            end
        end
        parent[i]==0 && throw(ErrorException(
            "$caller: nested Boolean result shell has no parent component"))
    end
    groups=Vector{Int}[]
    vol_of=zeros(Int,ncomp)
    for i in sortperm(inside_count)
        if iseven(inside_count[i])
            push!(groups,[i]);vol_of[i]=length(groups)
        else
            push!(groups[vol_of[parent[i]]],i)
            vol_of[i]=vol_of[parent[i]]
        end
    end
    return groups
end

function _brep_materialize_boolean!(m::GeoModel,t::Int,op::Symbol,
                                    ta::Int,tb::Int,caller)
    kept=_brep_collect_faces(m,op,ta,tb,caller)
    isempty(kept) && return nothing   # empty result — the caller drops the
                                      # preallocated records and binds nothing
    tol=_BREP_TOL
    # global result edge table + 3-D vertex dedupe, first-appearance order
    scale=1.0
    for p in values(m.points)
        scale=max(scale,abs(p[1]),abs(p[2]),abs(p[3]))
    end
    tol3=1e-9*scale
    res_edges=NamedTuple[]            # (curve,t0,t1) — same shape as operands
    res_verts=NTuple{3,Float64}[]
    emap=Dict{Any,Int}()
    function vertex_of(p)
        for (i,q) in enumerate(res_verts)
            _b3dist(p,q)<tol3 && return i
        end
        push!(res_verts,p);return length(res_verts)
    end
    # anchored section span ends materialize at the operand-edge pierce
    # eval — every fitted trace's endpoint sits a few ulps off it, which
    # would split one result vertex into two and dead-end the wire
    sev=nothing
    for spec in kept
        hasproperty(spec,:sev) && spec.sev!==nothing &&
            (sev=spec.sev;break)
    end
    vfix=Dict{Int,NamedTuple}()
    for spec in kept, wire in spec.wires, ar in wire
        ei=get!(emap,ar.key) do
            # vt0/vt1 are the 3-D curve's own parameters (imprint arcs carry
            # a reversed pcurve parametrization under mirrored frame maps)
            push!(res_edges,(curve=ar.crv,t0=ar.vt0,t1=ar.vt1))
            length(res_edges)
        end
        # a merged key spans several source pieces — widen to their union
        e=res_edges[ei]
        lo=min(e.t0,ar.vt0,ar.vt1);hi=max(e.t1,ar.vt0,ar.vt1)
        (lo==e.t0 && hi==e.t1) ||
            (res_edges[ei]=(curve=e.curve,t0=lo,t1=hi))
        if sev!==nothing && !haskey(vfix,ei) && ar.key isa Tuple &&
                length(ar.key)>=2 && ar.key[1]===:s &&
                1<=ar.key[2]<=length(sev)
            vfix[ei]=sev[ar.key[2]]
        end
    end
    # endpoints → vertices (deduped, first appearance in wire order); the
    # endpoints of a merged edge are its span extrema — interior piece
    # junctions are not result vertices
    edge_verts=Vector{Tuple{Int,Int}}(undef,length(res_edges))
    for (ei,e) in enumerate(res_edges)
        fx=get(vfix,ei,nothing)
        p0=fx!==nothing && e.t0==fx.t0 && fx.v0!==nothing ?
            fx.v0 : _brep_vertex_eval(e.curve,e.t0)
        p1=fx!==nothing && e.t1==fx.t1 && fx.v1!==nothing ?
            fx.v1 : _brep_vertex_eval(e.curve,e.t1)
        va=vertex_of(p0)
        vb=vertex_of(p1)
        edge_verts[ei]=(va,vb)
    end
    # shell components: faces adjacent through shared result edges
    parent=collect(1:length(kept))
    function findroot(uf,x)
        while uf[x]!=x
            uf[x]=uf[uf[x]];x=uf[x]
        end
        x
    end
    edge_faces=Dict{Any,Int}()
    for (fi,spec) in enumerate(kept), wire in spec.wires, ar in wire
        ei=emap[ar.key]
        j=get(edge_faces,ei,0)
        if j==0
            edge_faces[ei]=fi
        else
            ri,rj=findroot(parent,fi),findroot(parent,j)
            ri!=rj && (parent[ri]=rj)
        end
    end
    comps=Dict{Int,Vector{Int}}()
    for fi in 1:length(kept)
        push!(get!(Vector{Int},comps,findroot(parent,fi)),fi)
    end
    # operand order: A's kept faces precede B's, so sorting components by
    # first face index puts the A-derived component's solid first
    comp_faces=sort!(collect(values(comps));by=minimum)
    # exterior components first: a component inside another component's shell
    # is a cavity
    ncomp=length(comp_faces)
    inside=falses(ncomp,ncomp)          # inside[j,i] — comp i sits inside comp j
    inside_count=zeros(Int,ncomp)
    for i in 1:ncomp, j in 1:ncomp
        i==j && continue
        p=_brep_spec_rep(kept,comp_faces[i])
        if _brep_point_in_solid(
                _brep_spec_solid_view(kept,comp_faces[j]),p,tol,caller)
            inside[j,i]=true;inside_count[i]+=1
        end
    end
    groups=_brep_volume_groups(inside,inside_count,caller)
    # materialize — every new entity is tracked for rollback
    points=Int[];curves=Int[];loops=Int[];surfaces=Int[];shells=Int[]
    extra_volumes=Int[]
    try
        ptag=Vector{Int}(undef,length(res_verts))
        for (i,p) in enumerate(res_verts)
            ptag[i]=_add_occ_point!(m,p)
            push!(points,ptag[i])
        end
        # like the OCC primitive materializers, result vertices carry no
        # explicit mesh-size constraint — `_volume_boundary_size_field`
        # derives them from incident edges
        for p in points
            delete!(m.point_size,p)
            delete!(m.point_w,p)
        end
        ctag=Vector{Int}(undef,length(res_edges))
        for (i,e) in enumerate(res_edges)
            va,vb=edge_verts[i]
            crv=e.curve
            local c
            if crv.kind===:degenerate
                c=_add_occ_degenerate!(m,ptag[va])
            elseif crv.kind===:line
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_geometry[c]=(occ=:line,origin=crv.o,dir=crv.d,
                                     t0=e.t0,t1=e.t1)
            elseif crv.kind===:circle
                c=_add_occ_circle!(m,ptag[va],ptag[vb],crv.c,crv.n,
                                   crv.x,crv.y,crv.r,e.t0,e.t1)
            elseif crv.kind===:ellipse
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_types[c]=:ellipse
                m.curve_geometry[c]=(occ=:ellipse,center=crv.c,n=crv.n,
                                     X=crv.x,Y=crv.y,rx=crv.rx,ry=crv.ry,
                                     t0=e.t0,t1=e.t1)
            else
                # :bspline — re-emit the evaluator record restamped to the
                # result span; `crv.ctype` preserves the source curve kind
                # (:copied_curve/:tensor_bspline/:tensor_bezier)
                crv.kind===:bspline || throw(ArgumentError(
                    "$caller: unsupported result curve kind $(crv.kind)"))
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_types[c]=crv.ctype
                m.curve_geometry[c]=merge(crv.g,(t0=e.t0,t1=e.t1))
            end
            ctag[i]=c;push!(curves,c)
        end
        # a multi-component result extracts its own shells' components from the
        # shared snapshot mesh — the primary tag links to itself
        length(groups)>1 && (m.boolean_components[t]=t)
        for (gi,group) in enumerate(groups)
            volume_shells=Int[]
            for ci in group
                stags=Int[]
                for fi in comp_faces[ci]
                    spec=kept[fi]
                    pcdict=Dict{Int,NamedTuple}()
                    loop_tags=Int[]
                    for wire in spec.wires
                        ids=Int[]
                        for ar in wire
                            c=ctag[emap[ar.key]]
                            # a neg arc's pcurve params run opposite the edge
                            # parameter — uv-forward traversal is then
                            # edge-backward, so the OCC occurrence
                            # orientation (signed tag + pcurve slot) flips
                            efwd=ar.fwd ⊻ ar.neg
                            push!(ids,efwd ? c : -c)
                            # stored pcurves are edge-parametrized — mirrored
                            # imprints negate their pcurve's parameter
                            pcst=ar.neg ? _brep_flip_pc(ar.pc,caller) : ar.pc
                            ent=get(pcdict,c,nothing)
                            if efwd
                                if !(ent===nothing || ent.fwd===nothing)
                                    throw(ErrorException(
                                        "$caller: duplicate forward pcurve " *
                                        "on Curve[$c]"))
                                end
                                pcdict[c]=(fwd=pcst,
                                           rev=ent===nothing ? nothing :
                                                                  ent.rev)
                            else
                                ent===nothing || ent.rev===nothing ||
                                    throw(ErrorException(
                                        "$caller: duplicate reverse pcurve " *
                                        "on Curve[$c]"))
                                pcdict[c]=(fwd=ent===nothing ? nothing :
                                                                  ent.fwd,
                                           rev=pcst)
                            end
                        end
                        lt=add_curve_loop!(m,ids)
                        push!(loops,lt);push!(loop_tags,lt)
                    end
                    sg=merge(spec.sg,(pcurves=pcdict,))
                    st=_add_occ_surface!(m,spec.sg.occ,loop_tags,sg)
                    push!(surfaces,st);push!(stags,spec.sgn*st)
                    face=m.surface_geometry[st]
                    for wire in spec.wires, ar in wire
                        c=ctag[emap[ar.key]]
                        ar.crv.kind===:degenerate &&
                            _occ_bind_degenerate!(m,c,face)
                    end
                end
                sh=add_surface_loop!(m,stags)
                push!(shells,sh);push!(volume_shells,sh)
            end
            if gi==1
                m.volumes[t]=volume_shells
            else
                et=_alloc_tag!(m,3,0,caller)
                m.volumes[et]=volume_shells
                m.booleans[et]=(op=op,a=ta,b=tb)
                m.boolean_operands[et]=m.boolean_operands[t]
                m.boolean_components[et]=t
                push!(extra_volumes,et)
            end
        end
    catch
        for et in extra_volumes
            delete!(m.volumes,et);delete!(m.booleans,et)
            delete!(m.boolean_operands,et);delete!(m.boolean_components,et)
        end
        for sh in shells
            delete!(m.surface_loops,sh)
        end
        _occ_materialize_rollback!(m,points,curves,loops,surfaces,0)
        m.volumes[t]=Int[]
        delete!(m.boolean_components,t)
        rethrow()
    end
    return nothing
end

# ── multi-operand materialization ───────────────────────────────────────────
#
# Group an N-way Boolean's kept face regions into result solids. A fuse
# merges every kept cell through shared result edges; the other ops keep one
# piece per membership cell — a face between two kept cells enters both
# pieces' shells with opposite signs (OCC's shared partition wall). Each
# group's boundary components then decompose under the same containment
# forest as the binary path.
#
# Returns `(pieces, image, created)`:
#   pieces — sorted `(srcs, pseudo, shells)`; `srcs` is the union of the
#       piece's cell membership bits, `pseudo` is the operand index the piece
#       coincides with when it is exactly that operand's own unmodified solid
#       (0 otherwise), `shells` its materialized surface loops
#   image  — per operand, the indices of `pieces` forming its modified image
#       (union/difference credit objects only, matching OCC's GF maps)
#   created — every materialized entity tag, for caller-side rollback

function _brep_materialize_multi!(m::GeoModel,op::Symbol,tags,
                                  kept,touched,nops::Int,nobj::Int,caller)
    tol=_BREP_TOL
    scale=1.0
    for p in values(m.points)
        scale=max(scale,abs(p[1]),abs(p[2]),abs(p[3]))
    end
    tol3=1e-9*scale
    # piece grouping: fuse merges all cells; the cell ops group by mask
    if op===:union
        gspecs=[collect(eachindex(kept))]
        gmasks=[UInt64(0)]
    else
        gmasks=sort!(collect(Set{UInt64}(
            u[1] for spec in kept for u in spec.uses)))
        gspecs=[[fi for (fi,spec) in enumerate(kept)
                 if any(u->u[1]==gm,spec.uses)] for gm in gmasks]
    end
    pieces=NamedTuple[]     # (srcs,pseudo,comps,order)
    for (gi,fis) in enumerate(gspecs)
        n=length(fis)
        pos=Dict{Int,Int}(fi=>k for (k,fi) in enumerate(fis))
        parent=collect(1:n)
        function findroot(uf,x)
            while uf[x]!=x
                uf[x]=uf[uf[x]];x=uf[x]
            end
            x
        end
        edge_faces=Dict{Any,Int}()
        for fi in fis, wire in kept[fi].wires, ar in wire
            j=get(edge_faces,ar.key,0)
            if j==0
                edge_faces[ar.key]=fi
            else
                ri,rj=findroot(parent,pos[fi]),findroot(parent,pos[j])
                ri!=rj && (parent[ri]=rj)
            end
        end
        compsd=Dict{Int,Vector{Int}}()
        for fi in fis
            push!(get!(Vector{Int},compsd,findroot(parent,pos[fi])),fi)
        end
        comp_faces=sort!(collect(values(compsd));by=minimum)
        ncomp=length(comp_faces)
        inside=falses(ncomp,ncomp);inside_count=zeros(Int,ncomp)
        for i in 1:ncomp, j in 1:ncomp
            i==j && continue
            if _brep_point_in_solid(
                    _brep_spec_solid_view(kept,comp_faces[j]),
                    _brep_spec_rep(kept,comp_faces[i]),tol,caller)
                inside[j,i]=true;inside_count[i]+=1
            end
        end
        for vg in _brep_volume_groups(inside,inside_count,caller)
            comps=[comp_faces[ci] for ci in vg]
            srcs=UInt64(0)
            if op===:union
                for cf in comps, fi in cf, u in kept[fi].uses
                    srcs|=u[1]
                end
            else
                srcs=gmasks[gi]
            end
            push!(pieces,(srcs=srcs,pseudo=0,comps=comps,
                          order=minimum(minimum,comps)))
        end
    end
    # operand images — a piece credits every operand in its membership mask
    # (OCC's modified/generated maps list the piece under each member)
    image=[Int[] for _ in 1:nops]
    for (pi,p) in enumerate(pieces)
        s=p.srcs
        for i in 1:nops
            (s&_brep_bit(i))!=0 && push!(image[i],pi)
        end
    end
    # a piece is exactly operand i (OCC Extent==0 / IsSame) when i never
    # interacted and this piece is i's only image — all of i's faces then
    # bound it whole. Membership-mask cardinality is not the criterion: a
    # contained operand's face labels include the container's bit
    for (pi,p) in enumerate(pieces)
        for i in 1:nops
            (p.srcs&_brep_bit(i))!=0 || continue
            if !touched[i] && image[i]==[pi]
                pieces[pi]=merge(p,(pseudo=i,));break
            end
        end
    end
    # result order: fuse follows first-appearance (operand) order; the cell
    # ops order by membership — lowest member index, then mask value —
    # matching OCC's traversal of the GF decomposition
    if op===:union
        sort!(pieces,by=p->p.order)
    else
        sort!(pieces,by=p->(trailing_zeros(p.srcs),p.srcs,p.order))
    end
    # the image indices were computed pre-sort — rebuild them
    image=[Int[] for _ in 1:nops]
    for (pi,p) in enumerate(pieces)
        s=p.srcs
        for i in 1:nops
            (s&_brep_bit(i))!=0 && push!(image[i],pi)
        end
    end
    emit=falses(length(kept))
    for p in pieces
        p.pseudo==0 || continue
        for cf in p.comps, fi in cf
            emit[fi]=true
        end
    end
    # faces shared between a materialized piece and a preserved (pseudo)
    # operand reuse the operand's existing surface — OCC shares the face
    # TShape, so the neighbor piece's shell references the same Surface
    piece_of=Dict{UInt64,Int}(p.srcs=>pi for (pi,p) in enumerate(pieces))
    pseudo_surf=Dict{Int,Set{Int}}()
    pseudo_ctag=Dict{Int,Vector{Int}}()
    for p in pieces
        p.pseudo==0 && continue
        ss=Set{Int}()
        for sh in m.volumes[tags[p.pseudo]], sc in m.surface_loops[sh]
            push!(ss,abs(sc))
        end
        pseudo_surf[p.pseudo]=ss
        pseudo_ctag[p.pseudo]=_brep_extract(m,tags[p.pseudo],caller).ctags
    end
    reuse=falses(length(kept))
    reuse_ctag=Dict{Any,Int}()
    for (fi,spec) in enumerate(kept)
        emit[fi] || continue
        spec.stag>0 || continue
        i=0
        for u in spec.uses
            pi=get(piece_of,u[1],0)
            if pi!=0 && pieces[pi].pseudo!=0
                i=pieces[pi].pseudo;break
            end
        end
        i==0 && continue
        spec.stag in pseudo_surf[i] || continue
        ok=true
        pset=Set(pseudo_ctag[i])
        for w in spec.wires, ar in w
            (ar.key isa Tuple && length(ar.key)>=2 &&
             ar.key[1]===:e && ar.key[2] in pset) || (ok=false;break)
        end
        ok || continue
        reuse[fi]=true;emit[fi]=false
        for w in spec.wires, ar in w
            # canonical keys already carry the operand's own curve tag
            reuse_ctag[ar.key]=ar.key[2]
        end
    end
    # global result edge table + 3-D vertex dedupe, first-appearance order
    res_edges=NamedTuple[]
    res_verts=NTuple{3,Float64}[]
    emap=Dict{Any,Int}()
    function vertex_of(p)
        for (i,q) in enumerate(res_verts)
            _b3dist(p,q)<tol3 && return i
        end
        push!(res_verts,p);return length(res_verts)
    end
    # anchored section span ends materialize at the operand-edge pierce
    # eval — every fitted trace's endpoint sits a few ulps off it, which
    # would split one result vertex into two and dead-end the wire
    sev=nothing
    for spec in kept
        hasproperty(spec,:sev) && spec.sev!==nothing &&
            (sev=spec.sev;break)
    end
    vfix=Dict{Int,NamedTuple}()
    for (fi,spec) in enumerate(kept)
        emit[fi] || continue
        for wire in spec.wires, ar in wire
            ei=get!(emap,ar.key) do
                push!(res_edges,(curve=ar.crv,t0=ar.vt0,t1=ar.vt1))
                length(res_edges)
            end
            e=res_edges[ei]
            lo=min(e.t0,ar.vt0,ar.vt1);hi=max(e.t1,ar.vt0,ar.vt1)
            (lo==e.t0 && hi==e.t1) ||
                (res_edges[ei]=(curve=e.curve,t0=lo,t1=hi))
            if sev!==nothing && !haskey(vfix,ei) && ar.key isa Tuple &&
                    length(ar.key)>=2 && ar.key[1]===:s &&
                    1<=ar.key[2]<=length(sev)
                vfix[ei]=sev[ar.key[2]]
            end
        end
    end
    edge_verts=Vector{Tuple{Int,Int}}(undef,length(res_edges))
    for (ei,e) in enumerate(res_edges)
        fx=get(vfix,ei,nothing)
        p0=fx!==nothing && e.t0==fx.t0 && fx.v0!==nothing ?
            fx.v0 : _brep_vertex_eval(e.curve,e.t0)
        p1=fx!==nothing && e.t1==fx.t1 && fx.v1!==nothing ?
            fx.v1 : _brep_vertex_eval(e.curve,e.t1)
        va=vertex_of(p0)
        vb=vertex_of(p1)
        edge_verts[ei]=(va,vb)
    end
    points=Int[];curves=Int[];loops=Int[];surfaces=Int[];shells=Int[]
    pieces_out=Vector{NamedTuple}(undef,length(pieces))
    try
        ptag=Vector{Int}(undef,length(res_verts))
        for (i,p) in enumerate(res_verts)
            ptag[i]=_add_occ_point!(m,p)
            push!(points,ptag[i])
        end
        # same convention as the OCC primitive materializers: result
        # vertices carry no explicit mesh-size constraint
        for p in points
            delete!(m.point_size,p)
            delete!(m.point_w,p)
        end
        ctag=Vector{Int}(undef,length(res_edges))
        for (i,e) in enumerate(res_edges)
            va,vb=edge_verts[i]
            crv=e.curve
            local c
            if crv.kind===:degenerate
                c=_add_occ_degenerate!(m,ptag[va])
            elseif crv.kind===:line
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_geometry[c]=(occ=:line,origin=crv.o,dir=crv.d,
                                     t0=e.t0,t1=e.t1)
            elseif crv.kind===:circle
                c=_add_occ_circle!(m,ptag[va],ptag[vb],crv.c,crv.n,
                                   crv.x,crv.y,crv.r,e.t0,e.t1)
            elseif crv.kind===:ellipse
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_types[c]=:ellipse
                m.curve_geometry[c]=(occ=:ellipse,center=crv.c,n=crv.n,
                                     X=crv.x,Y=crv.y,rx=crv.rx,ry=crv.ry,
                                     t0=e.t0,t1=e.t1)
            else
                crv.kind===:bspline || throw(ArgumentError(
                    "$caller: unsupported result curve kind $(crv.kind)"))
                c=_alloc_tag!(m,1,0,caller)
                m.curves[c]=(ptag[va],ptag[vb])
                m.curve_types[c]=crv.ctype
                m.curve_geometry[c]=merge(crv.g,(t0=e.t0,t1=e.t1))
            end
            ctag[i]=c;push!(curves,c)
        end
        # one surface per emitted spec — a face shared by two cells is a
        # single surface referenced by both pieces' shells (opposite signs)
        stag=zeros(Int,length(kept))
        for (fi,spec) in enumerate(kept)
            if reuse[fi]
                stag[fi]=spec.stag
                continue
            end
            emit[fi] || continue
            pcdict=Dict{Int,NamedTuple}()
            loop_tags=Int[]
            for wire in spec.wires
                ids=Int[]
                for ar in wire
                    c=haskey(reuse_ctag,ar.key) ? reuse_ctag[ar.key] :
                        ctag[emap[ar.key]]
                    # a neg arc's pcurve params run opposite the edge
                    # parameter — uv-forward traversal is then edge-backward,
                    # so the OCC occurrence orientation flips
                    efwd=ar.fwd ⊻ ar.neg
                    push!(ids,efwd ? c : -c)
                    pcst=ar.neg ? _brep_flip_pc(ar.pc,caller) : ar.pc
                    ent=get(pcdict,c,nothing)
                    if efwd
                        if !(ent===nothing || ent.fwd===nothing)
                            throw(ErrorException(
                                "$caller: duplicate forward pcurve " *
                                "on Curve[$c]"))
                        end
                        pcdict[c]=(fwd=pcst,
                                   rev=ent===nothing ? nothing : ent.rev)
                    else
                        ent===nothing || ent.rev===nothing ||
                            throw(ErrorException(
                                "$caller: duplicate reverse pcurve " *
                                "on Curve[$c]"))
                        pcdict[c]=(fwd=ent===nothing ? nothing : ent.fwd,
                                   rev=pcst)
                    end
                end
                lt=try
                    add_curve_loop!(m,ids)
                catch err
                    if get(ENV,"BREP_DEBUG_LOOP","")!=""
                        println("DBGLOOP fi=$fi sg=$(spec.sg.occ) " *
                                "nwire=$(length(spec.wires)) ids=$ids")
                        for (k,id) in enumerate(ids)
                            a,b=m.curves[abs(id)]
                            pa=m.points[a];pb=m.points[b]
                            dprev=0.0
                            if k>1
                                a0,b0=m.curves[abs(ids[k-1])]
                                pend=m.points[ids[k-1]>0 ? b0 : a0]
                                dprev=_b3dist(pend,id>0 ? pa : pb)
                            end
                            println("  id=$id ends=($a,$b) " *
                                    "pa=$(round.(pa,digits=10)) " *
                                    "pb=$(round.(pb,digits=10)) " *
                                    "dgap=$(round(dprev,sigdigits=3))")
                        end
                        for (wi,ww) in enumerate(spec.wires), ar in ww
                            println("  w$wi key=$(ar.key) vt=($(ar.vt0)," *
                                    "$(ar.vt1)) fwd=$(ar.fwd) neg=$(ar.neg)")
                        end
                    end
                    rethrow()
                end
                push!(loops,lt);push!(loop_tags,lt)
            end
            sg=merge(spec.sg,(pcurves=pcdict,))
            st=_add_occ_surface!(m,spec.sg.occ,loop_tags,sg)
            push!(surfaces,st);stag[fi]=st
            face=m.surface_geometry[st]
            for wire in spec.wires, ar in wire
                c=haskey(reuse_ctag,ar.key) ? reuse_ctag[ar.key] :
                    ctag[emap[ar.key]]
                ar.crv.kind===:degenerate &&
                    _occ_bind_degenerate!(m,c,face)
            end
        end
        # shells per piece component — the shell entry sign is the spec's use
        # sign for this piece's cell
        for (pi,p) in enumerate(pieces)
            p.pseudo!=0 && continue
            pshells=Int[]
            for cf in p.comps
                stags=Int[]
                for fi in cf
                    spec=kept[fi]
                    sgn=if op===:union
                        spec.uses[1][2]
                    else
                        ui=findfirst(u->u[1]==p.srcs,spec.uses)
                        ui===nothing && throw(ErrorException(
                            "$caller: face piece lost its cell membership"))
                        spec.uses[ui][2]
                    end
                    push!(stags,sgn*stag[fi])
                end
                sh=add_surface_loop!(m,stags)
                push!(shells,sh);push!(pshells,sh)
            end
            pieces_out[pi]=(srcs=p.srcs,pseudo=p.pseudo,shells=pshells)
        end
        for (pi,p) in enumerate(pieces)
            p.pseudo==0 || (pieces_out[pi]=(srcs=p.srcs,pseudo=p.pseudo,
                                            shells=Int[]))
        end
    catch
        for sh in shells
            delete!(m.surface_loops,sh)
        end
        _occ_materialize_rollback!(m,points,curves,loops,surfaces,0)
        rethrow()
    end
    return pieces_out,image,(points=points,curves=curves,loops=loops,
                             surfaces=surfaces,shells=shells)
end

# ── multi-component result meshes ────────────────────────────────────────────
#
# A multi-component Boolean result owns one volume per disjoint solid. Extra
# volumes carry the same (op, operand-snapshot) records as the primary tag
# plus a `boolean_components` link; their boundary surface is the subset of
# the shared result mesh containing the volume's own materialized shells.

# Connected tri components of a surface mesh, unioned over shared nodes.
function _mesh_surface_components(mesh::Mesh)
    nt=ntris(mesh)
    parent=collect(1:nt)
    function root(i)
        while parent[i]!=i
            parent[i]=parent[parent[i]];i=parent[i]
        end
        i
    end
    first=Dict{Int32,Int}()
    for ti in 1:nt, k in 1:3
        v=mesh.tris[k,ti]
        j=get(first,v,0)
        if j==0
            first[v]=ti
        else
            ri,rj=root(ti),root(j)
            ri!=rj && (parent[ri]=rj)
        end
    end
    comps=Dict{Int,Vector{Int}}()
    for ti in 1:nt
        push!(get!(Vector{Int},comps,root(ti)),ti)
    end
    return collect(values(comps))
end

function _mesh_tris_subset(mesh::Mesh,trilist::Vector{Int})
    coords=NTuple{3,Float64}[]
    nmap=Dict{Int32,Int32}()
    tris=NTuple{3,Int32}[]
    for ti in trilist
        push!(tris,ntuple(3) do k
            v=mesh.tris[k,ti]
            get!(nmap,v) do
                push!(coords,(mesh.coords[1,v],mesh.coords[2,v],
                              mesh.coords[3,v]))
                Int32(length(coords))
            end
        end)
    end
    cm=Matrix{Float64}(undef,3,length(coords))
    for (i,c) in enumerate(coords)
        cm[:,i].=c
    end
    tm=Matrix{Int32}(undef,3,length(tris))
    for (i,tr) in enumerate(tris)
        tm[:,i].=tr
    end
    return Mesh(cm;tris=tm)
end

# A boundary vertex of a materialized shell — any point on its first face.
function _boolean_shell_rep(m::GeoModel,sh::Int,caller::AbstractString)
    for st in m.surface_loops[sh], lt in m.surfaces[abs(st)]
        ids=m.loops[lt]
        isempty(ids) && continue
        p=m.curves[abs(ids[1])][1]
        return m.points[p]
    end
    throw(ErrorException(
        "$caller: Boolean result shell $sh has no boundary vertices"))
end

# Boundary mesh of one result volume: the connected components of the full
# Boolean result mesh containing the volume's materialized shells.
function _boolean_component_surface(m::GeoModel,t::Int,full::Mesh,
                                    caller::AbstractString)
    reps=NTuple{3,Float64}[
        _boolean_shell_rep(m,sh,caller) for sh in m.volumes[t]]
    isempty(reps) && throw(ErrorException(
        "$caller: Boolean Volume[$t] has no materialized result shells"))
    comps=_mesh_surface_components(full)
    scale=1.0
    for i in axes(full.coords,2)
        scale=max(scale,abs(full.coords[1,i]),abs(full.coords[2,i]),
                  abs(full.coords[3,i]))
    end
    tol=1e-7*scale
    keep=falses(length(comps))
    for rep in reps
        best=0;bestd=Inf;secondd=Inf
        for (ci,tris) in enumerate(comps)
            d=Inf
            for ti in tris, k in 1:3
                v=full.tris[k,ti]
                dd=_b3dist(rep,(full.coords[1,v],full.coords[2,v],
                                full.coords[3,v]))
                dd<d && (d=dd)
            end
            if d<bestd
                secondd=bestd;bestd=d;best=ci
            elseif d<secondd
                secondd=d
            end
        end
        (best!=0 && bestd<=tol) || throw(ErrorException(
            "$caller: Boolean result component is missing from the result " *
            "surface mesh"))
        secondd>tol || throw(ErrorException(
            "$caller: Boolean result components are not separated enough to " *
            "assign mesh components unambiguously"))
        keep[best]=true
    end
    trilist=Int[ti for (ci,tris) in enumerate(comps) if keep[ci]
                for ti in tris]
    return _mesh_tris_subset(full,trilist)
end


# ── parametric surface tessellation ─────────────────────────────────────────
# `mesh_model_surface` meshes planar (and ruled/tric) faces; OCC curved faces
# mesh through their (u,v) domain instead. The boundary reuses the planar
# path's curve discretization — the same `_surface_curve_parameters` grading
# and canonical 3-D coordinates — so neighbouring faces weld shared-edge
# nodes bitwise. The wire's pcurves bound the face's uv domain; on periodic
# surfaces `u` unwraps along the traversal so the domain is an ordinary PSLG
# for `constrained_delaunay`. Interior vertices evaluate through
# `_occ_surface_point`; boundary vertices keep their canonical coordinates.
function _model_occ_curved_surface_mesh_tracked(m::GeoModel,t::Int,
                                        caller::AbstractString;
                                        initial_only::Bool=false,
                                        _mesh1d_options=nothing)
    sg=get(m.surface_geometry,t,nothing)
    occ=sg===nothing ? nothing : get(sg,:occ,nothing)
    occ in (:cylinder,:sphere,:cone,:torus,:tensor_bspline,:tensor_bezier) ||
        throw(ArgumentError(
            "$caller: Surface[$t] has no parametric OCC record"))
    pcurves=get(sg,:pcurves,nothing)
    pcurves===nothing && throw(ArgumentError(
        "$caller: Surface[$t] carries no pcurves"))
    period=_brep_uperiod(sg)
    forced=Dict{Int,Vector{Float64}}()
    _attribute_forced_parameters(m,t,forced,caller;
        _mesh1d_options=_mesh1d_options)
    xs=Float64[];ys=Float64[]
    lift=NTuple{3,Float64}[]          # canonical 3-D per PSLG vertex
    segs=Tuple{Int,Int}[]
    index=Dict{NTuple{2,Float64},Int}()
    # (owner, parametric) provenance per emitted uv vertex — a boundary
    # sample's Point/Curve identity and curve parameter, or the face itself
    # plus its (u,v) for interior vertices. Same `emit` dedup rule as `lift`:
    # the first writer of a uv position owns the metadata too.
    vertex_meta=Dict{NTuple{2,Float64},
        Tuple{Tuple{Int,Int32},Union{Nothing,Vector{Float64}}}}()
    function emit(uv,p3,owner,param)
        return get!(index,uv) do
            push!(xs,uv[1]);push!(ys,uv[2]);push!(lift,p3)
            vertex_meta[uv]=(owner,param)
            return length(xs)
        end
    end
    for loop_id in m.surfaces[t]
        _verify_loop_closed(m,loop_id,caller,"Surface[$t]")
        loop_idx=Int[]
        u_prev=nothing
        for signed in m.loops[loop_id]
            curve=abs(signed)
            a,b=m.curves[curve]
            t0,t1=_model_curve_param_bounds(m,curve,caller)
            ent=get(pcurves,curve,nothing)
            ent===nothing && throw(ArgumentError(
                "$caller: Surface[$t] has no pcurve for Curve[$curve]"))
            crec=get(m.curve_geometry,curve,nothing)
            pc=signed<0 && ent.rev!==nothing ? ent.rev : ent.fwd
            if (crec!==nothing && get(crec,:occ,nothing)===:degenerate) ||
                    get(m.curve_types,curve,nothing)===:degenerate
                # A degenerate edge is a single 3-D vertex, but its pcurve
                # still bounds the uv domain: a point-like trace simply
                # deduplicates into a neighbour vertex, while a real span
                # collapses an entire domain row onto the pole. Emit the
                # collapsed row as an ordinary boundary chain — the pole fan.
                # Every sample lifts to the canonical pole vertex, so the
                # bitwise weld below folds the row into the fan apex and each
                # CDT cell touching two row vertices collapses out of the
                # lifted mesh. The row is deliberately NOT u-periodic
                # unwrapped: its two endpoint representatives legitimately
                # differ by the domain's full period, and anchoring the row to
                # `u_prev` would fold the whole wire onto one seam column.
                u0,v0=_occ_pcurve_d0(pc,t0)
                u1,v1=_occ_pcurve_d0(pc,t1)
                tol=1e-9*max(1.0,abs(u0),abs(u1),abs(v0),abs(v1))
                if abs(u1-u0)>tol || abs(v1-v0)>tol
                    # A genuine degenerate edge collapses onto ONE vertex;
                    # distinct endpoints on a row-spanning pcurve is a
                    # corrupt record, not a fan boundary.
                    a==b || throw(ArgumentError(
                        "$caller: Surface[$t] degenerate Curve[$curve] " *
                        "has distinct endpoints — not a collapsed " *
                        "pole/apex"))
                    pole=m.points[a]
                    pscale=max(1.0,abs(pole[1]),abs(pole[2]),abs(pole[3]))
                    ptol=1e-8*pscale
                    uvtol=1e-9*max(1.0,hypot(u1-u0,v1-v0))
                    plist=get(forced,curve,nothing)
                    ordered=plist===nothing ? Float64[t0,t1] :
                                             collect(Float64,plist)
                    isempty(ordered) && (ordered=Float64[t0,t1])
                    # `_surface_curve_parameters` endpoint snapping, so a
                    # stored parameter an ulp off a bound still emits the
                    # exact corner representative.
                    snap=128eps(Float64)*max(1.0,abs(t1-t0))
                    for i in eachindex(ordered)
                        ordered[i]-t0<=snap && (ordered[i]=t0)
                        t1-ordered[i]<=snap && (ordered[i]=t1)
                    end
                    signed<0 && reverse!(ordered)
                    tfirst=signed>0 ? t0 : t1
                    tlast=signed>0 ? t1 : t0
                    first(ordered)==tfirst || pushfirst!(ordered,tfirst)
                    last(ordered)==tlast || push!(ordered,tlast)
                    # Flatten the collapsed trace: a straight lin2d row keeps
                    # just its endpoint representatives; a curved trace
                    # subdivides until the polyline follows the pcurve.
                    seq=Float64[]
                    function _pole_row_flat(pa,pb,depth)
                        pm=(pa+pb)/2
                        if depth<10
                            ua,va=_occ_pcurve_d0(pc,pa)
                            ub,vb=_occ_pcurve_d0(pc,pb)
                            um,vm=_occ_pcurve_d0(pc,pm)
                            if hypot(um-(ua+ub)/2,vm-(va+vb)/2)>uvtol
                                _pole_row_flat(pa,pm,depth+1)
                                _pole_row_flat(pm,pb,depth+1)
                                return nothing
                            end
                        end
                        push!(seq,pa)
                        return nothing
                    end
                    for k in 2:length(ordered)
                        _pole_row_flat(ordered[k-1],ordered[k],0)
                    end
                    push!(seq,last(ordered))
                    for parameter in seq
                        u,v=_occ_pcurve_d0(pc,parameter)
                        # The collapsed row must lift to the pole — a pcurve
                        # that does not is a corrupt degenerate record, not a
                        # meshable fan boundary.
                        lifted=_occ_surface_point(sg,u,v)
                        _b3dist(lifted,pole)<=ptol || throw(ArgumentError(
                            "$caller: Surface[$t] degenerate Curve[$curve] " *
                            "pcurve does not collapse onto its pole vertex"))
                        u_prev=u
                        # Every sample on the collapsed row is the same 3-D
                        # vertex — the degenerate edge's own Point entity.
                        push!(loop_idx,emit((u,v),pole,(0,Int32(a)),
                            Float64[]))
                    end
                    continue
                end
            end
            for parameter in _surface_curve_parameters(
                    m,forced,curve,signed,caller)
                if parameter==t0
                    p3=m.points[a];owner=(0,Int32(a));param=Float64[]
                elseif parameter==t1
                    p3=m.points[b];owner=(0,Int32(b));param=Float64[]
                else
                    p3=_periodic_curve_point(m,curve,parameter,caller)
                    owner=(1,Int32(curve));param=Float64[parameter]
                end
                u,v=_occ_pcurve_d0(pc,parameter)
                if period>0.0 && u_prev!==nothing
                    u+=period*round((u_prev-u)/period)
                end
                u_prev=u
                push!(loop_idx,emit((u,v),p3,owner,param))
            end
        end
        nloop=length(loop_idx)
        nloop>=3 || throw(ArgumentError(
            "$caller: Loop[$loop_id] needs at least three points"))
        for k in 1:nloop
            i,j=loop_idx[k],loop_idx[mod1(k+1,nloop)]
            i==j || (xs[i]==xs[j] && ys[i]==ys[j]) || push!(segs,(i,j))
        end
    end
    isempty(segs) && throw(ArgumentError(
        "$caller: Surface[$t] has an empty parametric boundary"))
    # Interior sampling: a boundary-only CDT of a sparse uv domain produces
    # chords between distant angles that dive deep into the solid (most
    # visible when the two rims of a lofted/cylindrical band carry different
    # node counts).  Seed a uv lattice at the boundary's 3-D target spacing
    # mapped through the local patch metric so lifted triangles stay local.
    umin=minimum(xs);umax=maximum(xs);vmin=minimum(ys);vmax=maximum(ys)
    uspan=umax-umin;vspan=vmax-vmin
    # Upstream `ALGO_2D_INITIAL_ONLY` stops at the recovered 1-D boundary
    # mesh, so initial-only suppresses every interior (u,v) lattice seed.
    if !initial_only && uspan>0.0 && vspan>0.0
        # Direction-aware targets: a swept band's boundary mixes long
        # v-direction seams with short u-direction rim chords, so a single
        # mean spacing starves the short direction of lattice columns and
        # the CDT bridges whole profile arcs — the lifted triangles then fold
        # through the sweep (a self-intersecting patch).  Classify each
        # boundary edge by its dominant uv direction and target each lattice
        # direction at its own rim spacing.
        h3u=0.0;nu_e=0;h3v=0.0;nv_e=0
        for (i,j) in segs
            a=lift[i];b=lift[j]
            l3=sqrt((b[1]-a[1])^2+(b[2]-a[2])^2+(b[3]-a[3])^2)
            if abs(xs[j]-xs[i])>=abs(ys[j]-ys[i])
                h3u+=l3;nu_e+=1
            else
                h3v+=l3;nv_e+=1
            end
        end
        h3=(h3u+h3v)/max(nu_e+nv_e,1)
        # A direction whose ONLY boundary edges collapse to zero 3-D length —
        # a sphere's pole rows lift to a single apex vertex — has no rim
        # spacing of its own; fall back to the global mean rather than
        # seeding at a zero (Inf uv) step.
        h3u=nu_e>0 && h3u>0.0 ? h3u/nu_e : h3
        h3v=nv_e>0 && h3v>0.0 ? h3v/nv_e : h3
        # The metric of a swept/lofted patch varies across the domain
        # (a pipe's inner bend, a cone's apex); take the densest of a few
        # stations so no region is seeded coarser than its own curvature.
        su=0.0;sv=0.0
        for fu in (0.25,0.5,0.75), fv in (0.25,0.5,0.75)
            _,du_c,dv_c=_occ_surface_d1(sg,umin+fu*uspan,vmin+fv*vspan)
            su=max(su,sqrt(_b3dot(du_c,du_c)))
            sv=max(sv,sqrt(_b3dot(dv_c,dv_c)))
        end
        du_t=h3u/max(su,1e-12);dv_t=h3v/max(sv,1e-12)
        # Skip seeding in a direction whose uv step is not a positive finite
        # length (a singular metric, or a fully collapsed boundary) — an
        # empty lattice degrades to boundary-only CDT, never a hard error.
        nu=du_t>0.0 && isfinite(uspan/du_t) ?
            clamp(floor(Int,uspan/du_t)-1,0,64) : 0
        nv=dv_t>0.0 && isfinite(vspan/dv_t) ?
            clamp(floor(Int,vspan/dv_t)-1,0,64) : 0
        # Degrade proportionally under the seed cap rather than dropping the
        # lattice wholesale — an empty interior is exactly what seeds prevent.
        while nu>1 && nv>1 && nu*nv>4096
            nu>=nv ? (nu=(nu+1)÷2) : (nv=(nv+1)÷2)
        end
        if nu>0 && nv>0 && nu*nv<=4096
            du=uspan/(nu+1);dv=vspan/(nv+1)
            tol2=(0.35*min(du,dv))^2
            for i in 1:nu,j in 1:nv
                u=umin+i*du;v=vmin+j*dv
                # skip lattice points that would create slivers against the
                # boundary discretization
                near=false
                for k in eachindex(xs)
                    (xs[k]-u)^2+(ys[k]-v)^2<tol2 && (near=true;break)
                end
                near && continue
                p3=_occ_surface_point(sg,u,v)
                emit((u,v),p3,(2,Int32(t)),[u,v])
            end
        end
    end
    T=constrained_delaunay(xs,ys,segs)
    interior=classify_interior(T)
    mesh=to_mesh(T;interior)
    # Lift: boundary vertices keep the canonical 3-D they share with the
    # neighbouring faces' curve nodes; interior vertices evaluate the
    # analytic patch.
    canonical=Dict{NTuple{2,Float64},NTuple{3,Float64}}()
    for k in eachindex(lift)
        canonical[(xs[k],ys[k])]=lift[k]
    end
    uvcoords=mesh.coords
    coords=Matrix{Float64}(undef,3,nnodes(mesh))
    for k in 1:nnodes(mesh)
        u=uvcoords[1,k];v=uvcoords[2,k]
        p3=get(canonical,(u,v),nothing)
        p3===nothing && (p3=_occ_surface_point(sg,u,v))
        coords[1,k]=p3[1];coords[2,k]=p3[2];coords[3,k]=p3[3]
    end
    # Degenerate-face collapse: the two uv boundary representatives of a pole
    # or apex (the degenerate edge's endpoints at u=0 and u=a) lift to the same
    # 3-D vertex, so a uv triangle touching both collapses to an edge. Weld
    # bitwise-identical lifted nodes and drop collapsed cells — each surviving
    # edge still pairs the same two neighbours, so incidence is preserved.
    welded=Dict{NTuple{3,Float64},Int32}()
    remap=Vector{Int32}(undef,nnodes(mesh))
    kept=NTuple{3,Float64}[]
    kept_uv=NTuple{2,Float64}[]
    for k in 1:nnodes(mesh)
        p3=(coords[1,k],coords[2,k],coords[3,k])
        remap[k]=get!(welded,p3) do
            push!(kept,p3)
            push!(kept_uv,(uvcoords[1,k],uvcoords[2,k]))
            return Int32(length(kept))
        end
    end
    # Periodic-seam repair: on a closed-section face the seam appears as two uv
    # boundary columns whose nodes weld pairwise in 3-D. An interior uv edge
    # (seamL_i, X) and its twin (seamR_i, X) then weld to the same 3-D edge and
    # would carry four incident triangles. Whenever a welded edge collects more
    # than two incident uv triangles, subdivide surplus interior preimage edges
    # at their uv midpoint; the split halves weld uniquely, so every repair
    # strictly reduces the over-use count and cannot create a new collision.
    uxs=[uvcoords[1,k] for k in 1:nnodes(mesh)]
    uys=[uvcoords[2,k] for k in 1:nnodes(mesh)]
    # Preimage provenance parallel to uxs/uys: emitted vertices carry their
    # recorded owner/parameter; CDT Steiner vertices are face-owned with
    # their uv as the parametrization.
    pre_owner=Vector{Tuple{Int,Int32}}(undef,length(uxs))
    pre_param=Vector{Union{Nothing,Vector{Float64}}}(undef,length(uxs))
    for k in eachindex(uxs)
        meta=get(vertex_meta,(uxs[k],uys[k]),nothing)
        if meta===nothing
            pre_owner[k]=(2,Int32(t));pre_param[k]=[uxs[k],uys[k]]
        else
            pre_owner[k]=meta[1];pre_param[k]=meta[2]
        end
    end
    tris_uv=NTuple{3,Int32}[(mesh.tris[1,c],mesh.tris[2,c],mesh.tris[3,c])
                            for c in axes(mesh.tris,2)]
    cells=NTuple{3,Int32}[]
    seen_cells=Set{NTuple{3,Int32}}()
    for _ in 1:(8*length(tris_uv)+8)
        # uv edge incidence (true triangulation degrees) and welded-edge
        # incidence over surviving, deduplicated cells.
        uvcnt=Dict{NTuple{2,Int32},Int}()
        for (a,b,c) in tris_uv, (x,y) in ((a,b),(b,c),(c,a))
            uvk=x<y ? (x,y) : (y,x)
            uvcnt[uvk]=get(uvcnt,uvk,0)+1
        end
        einc=Dict{NTuple{2,Int32},Int}()
        preim=Dict{NTuple{2,Int32},Vector{NTuple{2,Int32}}}()
        empty!(cells);empty!(seen_cells)
        for (a,b,c) in tris_uv
            wa,wb,wc=remap[a],remap[b],remap[c]
            (wa==wb || wb==wc || wc==wa) && continue
            tri=(wa,wb,wc)
            key=tuple(sort!(Int[tri[1],tri[2],tri[3]])...)
            key in seen_cells && continue
            push!(seen_cells,key);push!(cells,tri)
            for (x,y) in ((a,b),(b,c),(c,a))
                wx,wy=remap[x],remap[y]
                wx==wy && continue
                we=wx<wy ? (wx,wy) : (wy,wx)
                einc[we]=get(einc,we,0)+1
                uvk=x<y ? (x,y) : (y,x)
                lst=get!(preim,we,NTuple{2,Int32}[])
                uvk in lst || push!(lst,uvk)
            end
        end
        bad_e=nothing
        for (e,c) in einc
            c>2 && (bad_e===nothing || e<bad_e) && (bad_e=e)
        end
        bad_e===nothing && break
        cands=[e for e in preim[bad_e] if get(uvcnt,e,0)==2]
        isempty(cands) && throw(ArgumentError(
            "$caller: Surface[$t] weld repair found no interior edge to split"))
        (a,b)=minimum(cands)
        mu=(uxs[a]+uxs[b])/2; mv=(uys[a]+uys[b])/2
        p3=_occ_surface_point(sg,mu,mv)
        push!(uxs,mu);push!(uys,mv)
        push!(pre_owner,(2,Int32(t)));push!(pre_param,[mu,mv])
        push!(remap,get!(welded,p3) do
            push!(kept,p3);push!(kept_uv,(mu,mv))
            return Int32(length(kept))
        end)
        mn=Int32(length(uxs))
        split=NTuple{3,Int32}[]
        for (x,y,z) in tris_uv
            if x==a && y==b
                push!(split,(a,mn,z),(mn,b,z))
            elseif x==b && y==a
                push!(split,(b,mn,z),(mn,a,z))
            elseif y==a && z==b
                push!(split,(a,mn,x),(mn,b,x))
            elseif y==b && z==a
                push!(split,(b,mn,x),(mn,a,x))
            elseif z==a && x==b
                push!(split,(a,mn,y),(mn,b,y))
            elseif z==b && x==a
                push!(split,(b,mn,y),(mn,a,y))
            else
                push!(split,(x,y,z))
            end
        end
        tris_uv=split
    end
    # Vertex-link repair: on a closed-section face the two uv seam columns
    # weld onto the same 3-D nodes, but CDT chords that reach both columns
    # asymmetrically can leave a welded vertex with a disconnected triangle
    # link (a pinch — two fans meeting only at the vertex).  Rebuild the
    # vertex's whole star on one covering sheet — shifting each incident
    # preimage's neighbourhood back by the domain's u period — so the star
    # becomes a single polygon that a fresh CDT fans around the vertex.
    umin_r=minimum(uxs);umax_r=maximum(uxs);P_r=umax_r-umin_r
    if P_r>0.0
        for _ in 1:16
            empty!(cells);empty!(seen_cells)
            for (a,b,c) in tris_uv
                wa,wb,wc=remap[a],remap[b],remap[c]
                (wa==wb || wb==wc || wc==wa) && continue
                tri=(wa,wb,wc)
                key=tuple(sort!(Int[tri[1],tri[2],tri[3]])...)
                key in seen_cells && continue
                push!(seen_cells,key);push!(cells,tri)
            end
            pre_of=Dict{Int32,Vector{Int32}}()
            for (i,r) in enumerate(remap)
                push!(get!(pre_of,r,Int32[]),Int32(i))
            end
            badw=Int32(0)
            for (a,b,c) in cells, w in (a,b,c)
                pre=get(pre_of,w,nothing)
                (pre===nothing || length(pre)<2) && continue
                lk=Dict{Int32,Vector{Int32}}()
                for (x,y,z) in cells
                    p1,p2=x==w ? (y,z) : y==w ? (z,x) :
                          z==w ? (x,y) : continue
                    push!(get!(lk,p1,Int32[]),p2)
                    push!(get!(lk,p2,Int32[]),p1)
                end
                seen_v=Set{Int32}();ncomp=0
                for s in keys(lk)
                    s in seen_v && continue
                    ncomp+=1;stack=Int32[s]
                    while !isempty(stack)
                        z=pop!(stack);z in seen_v && continue
                        push!(seen_v,z)
                        for q2 in lk[z]
                            q2 in seen_v || push!(stack,q2)
                        end
                    end
                end
                ncomp>1 && (badw=w;break)
            end
            badw==0 && break
            pre=pre_of[badw]
            qmin=pre[argmin([uxs[q] for q in pre])]
            up=uxs[qmin];wv=uys[qmin]
            # Every incident cell contains the welded vertex, so its union on
            # the covering sheet is star-shaped about w: shifting each cell so
            # that its badw preimage lands at (up,wv) places every neighbour
            # at a unique position — including the two seam preimages of a
            # welded neighbour, which coincide there.  The star polygon is
            # then the neighbours in angular order around w.
            inc=Set{Int}()
            npmap=Dict{NTuple{2,Float64},Int32}()
            for (ci,(a,b,c)) in enumerate(tris_uv)
                np=(remap[a]==badw)+(remap[b]==badw)+(remap[c]==badw)
                np==1 || continue
                q=remap[a]==badw ? a : remap[b]==badw ? b : c
                d=up-uxs[q]
                push!(inc,ci)
                for v in (a,b,c)
                    v==q && continue
                    pu=uxs[v]+d
                    pos=(pu==0.0 ? 0.0 : pu, uys[v]==0.0 ? 0.0 : uys[v])
                    old=get(npmap,pos,Int32(0))
                    if old==0
                        npmap[pos]=v
                    elseif pu>up
                        # right of w: keep the umin-side preimage
                        uxs[v]<uxs[old] && (npmap[pos]=v)
                    elseif pu<up
                        # left of w: keep the umax-side preimage
                        uxs[v]>uxs[old] && (npmap[pos]=v)
                    end
                end
            end
            length(npmap)>=3 || break
            order=sort!(collect(keys(npmap));
                        by=p->atan(p[2]-wv,p[1]-up))
            # A welded neighbour may occupy several covering positions (one
            # per sheet): the cylinder link visits it once, so pick one
            # representative position per welded vertex — the nearest to w —
            # for the fan; the skipped copies form rim "ears" that are
            # retriangulated separately without w.
            chosen=Dict{Int32,NTuple{2,Float64}}()
            for p in order
                W=remap[npmap[p]]
                old=get(chosen,W,nothing)
                if old===nothing ||
                        (p[1]-up)^2+(p[2]-wv)^2 <
                        (old[1]-up)^2+(old[2]-wv)^2
                    chosen[W]=p
                end
            end
            repset=Set(values(chosen))
            fan=NTuple{2,Float64}[p for p in order if p in repset]
            length(fan)>=3 || break
            newt=NTuple{3,Int32}[]
            # fan triangles around w over the representative positions
            fi=indexin(fan,order)
            for k in eachindex(fan)
                a=npmap[fan[k]]
                b=npmap[fan[mod1(k+1,length(fan))]]
                push!(newt,(qmin,a,b))
            end
            # ears: runs of skipped positions between consecutive fan verts
            for (k,ia) in enumerate(fi)
                ib=fi[mod1(k+1,length(fi))]
                skips=Int[]
                j=ia
                while true
                    j=mod1(j+1,length(order))
                    j==ib && break
                    push!(skips,j)
                end
                isempty(skips) && continue
                epoly=NTuple{2,Float64}[fan[k]]
                append!(epoly,(order[j] for j in skips))
                push!(epoly,fan[mod1(k+1,length(fan))])
                exs=Float64[p[1] for p in epoly]
                eys=Float64[p[2] for p in epoly]
                en=length(epoly)
                esegs=Tuple{Int,Int}[(i,mod1(i+1,en)) for i in 1:en]
                T3=constrained_delaunay(exs,eys,esegs)
                interior3=classify_interior(T3)
                em=to_mesh(T3;interior=interior3)
                for c in axes(em.tris,2)
                    vs=Int32[]
                    ok=true
                    for r in (em.tris[1,c],em.tris[2,c],em.tris[3,c])
                        pos=(em.coords[1,r]==0.0 ? 0.0 : em.coords[1,r],
                             em.coords[2,r]==0.0 ? 0.0 : em.coords[2,r])
                        vi=get(npmap,pos,Int32(0))
                        if vi==0
                            # CDT Steiner point in the ear
                            ub=umin_r+mod(pos[1]-umin_r,P_r)
                            ub>=umax_r && (ub=umin_r)
                            p3=_occ_surface_point(sg,ub,pos[2])
                            push!(uxs,ub);push!(uys,pos[2])
                            push!(pre_owner,(2,Int32(t)))
                            push!(pre_param,[ub,pos[2]])
                            push!(remap,get!(welded,p3) do
                                push!(kept,p3);push!(kept_uv,(ub,pos[2]))
                                return Int32(length(kept))
                            end)
                            vi=Int32(length(uxs))
                            npmap[pos]=vi
                        end
                        push!(vs,vi)
                    end
                    ok && push!(newt,(vs[1],vs[2],vs[3]))
                end
            end
            isempty(newt) && break
            tris_uv=NTuple{3,Int32}[
                t for (i,t) in enumerate(tris_uv) if !(i in inc)]
            append!(tris_uv,newt)
        end
    end
    post=Dict{NTuple{2,Int32},Int}()
    for tri in cells, (x,y) in ((tri[1],tri[2]),(tri[2],tri[3]),(tri[3],tri[1]))
        we=x<y ? (x,y) : (y,x)
        post[we]=get(post,we,0)+1
    end
    for (we,c) in post
        c<=2 || throw(ArgumentError(
            "$caller: Surface[$t] seam weld did not converge to a manifold mesh"))
    end
    # Resolve each welded node's ownership from its preimages: the
    # lowest-dimensional owner wins (a Point beats a Curve beats the face),
    # matching the precedence `model_to_mixed` applies on planar faces, and
    # ties keep the first preimage — deterministic by emit order. The stored
    # parameter follows the winning preimage, so a curve-owned boundary
    # vertex keeps its curve parameter and a face-owned vertex keeps a uv.
    kept_owner=Vector{Tuple{Int,Int32}}(undef,length(kept))
    kept_param=Vector{Union{Nothing,Vector{Float64}}}(undef,length(kept))
    resolved=falses(length(kept))
    for (k,r) in enumerate(remap)
        owner=pre_owner[k]
        (!resolved[r] || owner[1]<kept_owner[r][1]) || continue
        resolved[r]=true;kept_owner[r]=owner;kept_param[r]=pre_param[k]
    end
    # Periodic faces unwrap u across the domain seam while traversing the
    # wire; fold stored surface parameters back into the analytic interval so
    # `model_parametrization`/`model_value` round-trip on the native bounds.
    if period>0.0
        ulo=_occ_surface_bounds(sg)[1][1]
        for r in eachindex(kept_param)
            kept_owner[r][1]==2 || continue
            p=kept_param[r]
            p===nothing && continue
            p[1]=ulo+mod(p[1]-ulo,period)
        end
    end
    coords=Matrix{Float64}(undef,3,length(kept))
    for k in eachindex(kept)
        coords[1,k]=kept[k][1];coords[2,k]=kept[k][2];coords[3,k]=kept[k][3]
    end
    tris=Matrix{Int32}(undef,3,length(cells))
    for k in eachindex(cells)
        tris[1,k]=cells[k][1];tris[2,k]=cells[k][2];tris[3,k]=cells[k][3]
    end
    mesh=Mesh(coords;tris=tris)
    isempty(cells) && throw(ArgumentError(
        "$caller: Surface[$t] parametric triangulation collapsed entirely"))
    # The unwrapped (u,v) triangulation behind the welded product: the
    # algorithm-7 metric producer refines this preimage representation —
    # seam columns and collapsed degenerate rows stay distinct parametric
    # vertices sharing world positions — so its certified local operations
    # run on an ordinary manifold disk before the same weld folds the
    # refined sheet back onto the closed surface. `world` carries the
    # canonical lifted position per row (`kept[remap]`, so boundary rows
    # keep the curve product's bitwise coordinates and every row of a weld
    # class shares one position); `owners`/`params` are the per-row
    # provenance `pre_owner`/`pre_param`; `cells` index the unwrapped rows.
    unwrapped=(u=uxs,v=uys,cells=tris_uv,owners=pre_owner,params=pre_param,
        world=NTuple{3,Float64}[kept[remap[k]] for k in eachindex(uxs)],
        period=period)
    # The uv-CCW winding lifts to either side of the face depending on the
    # parametrization's handedness; align the patch with the face's
    # geometric normal like a planar patch follows its `axis`.
    cell=mesh.tris[:,1]
    um=(kept_uv[cell[1]][1]+kept_uv[cell[2]][1]+kept_uv[cell[3]][1])/3
    vm=(kept_uv[cell[1]][2]+kept_uv[cell[2]][2]+kept_uv[cell[3]][2])/3
    nn=_brep_surface_normal(sg,um,vm)
    pa=coords[:,cell[1]];pb=coords[:,cell[2]];pc_=coords[:,cell[3]]
    n3=_b3cross(_b3(pb,pa),_b3(pc_,pa))
    if _b3dot(n3,nn)<0
        flipped=Matrix{Int32}(undef,3,size(mesh.tris,2))
        flipped[1,:].=mesh.tris[1,:]
        flipped[2,:].=mesh.tris[3,:]
        flipped[3,:].=mesh.tris[2,:]
        mesh=Mesh(coords;tris=flipped)
    end
    return mesh,kept_owner,kept_param,unwrapped
end

# The lifted tessellation alone — the volume-PLC consumer classifies boundary
# vertices geometrically and does not need per-node provenance.
function _model_occ_curved_surface_mesh(m::GeoModel,t::Int,
                                        caller::AbstractString;kwargs...)
    mesh,_,_,_=_model_occ_curved_surface_mesh_tracked(m,t,caller;kwargs...)
    return mesh
end

# Wrap the lifted (u,v) tessellation in the native producer's identity record,
# the same contract `_gmsh_native_surface_export` emits so `generate(2)`
# classifies the patch like any other surface product: boundary vertices that
# coincide bitwise with an actual Point/Curve product node keep that node's
# UID (retained), every other vertex allocates above the consumed high-water
# mark, and each node carries its owner plus parametric coordinates — the
# curve parameter on boundary edges and the face's own (u,v) inside.
function _occ_curved_surface_production(m::GeoModel,t::Int,mesh::Mesh,
                                        owners::Vector{Tuple{Int,Int32}},
                                        params,input,limit,caller::AbstractString)
    tags=Vector{UInt64}(undef,nnodes(mesh))
    retained=falses(nnodes(mesh))
    lookup=Dict{Tuple{Int,Int32,NTuple{3,Float64}},UInt64}()
    source=input===nothing ? nothing : input.mesh
    if source!==nothing && source.entity_data!==nothing
        data=source.entity_data
        for row in 1:nnodes(source)
            entity=data.node_entities[row]
            key=(entity[1],entity[2],
                 (source.coords[1,row],source.coords[2,row],
                  source.coords[3,row]))
            haskey(lookup,key) || (lookup[key]=data.external_node_tags[row])
        end
    end
    next=UInt64(limit)
    for k in 1:nnodes(mesh)
        owner=owners[k]
        point=(mesh.coords[1,k],mesh.coords[2,k],mesh.coords[3,k])
        tag=get(lookup,(owner[1],owner[2],point),UInt64(0))
        if tag!=0
            tags[k]=tag;retained[k]=true
        else
            next+=UInt64(1)
            iszero(next) && throw(ArgumentError(
                "$caller: Surface[$t] vertex births exhaust UInt64 labels"))
            tags[k]=next
        end
    end
    return _ModelSurfaceProduction(mesh,tags,owners,params,retained,next)
end

# ── result surfaces (any Boolean record) ────────────────────────────────────
#
# `m.boolean_operands` carries either a binary snapshot `(A,B)` or a
# multi-operand record `(meshes, cell)`. A multi piece's boundary mesh
# composes its arrangement cell — `∩ members ∖ ∪ non-members` — from binary
# `mesh_boolean` calls on the operation-time snapshots; a fuse piece takes the
# full fused mesh and lets the component extraction pick its own shells.

# Boundary mesh of one arrangement cell over operand snapshot `meshes`:
# intersect the members, subtract every other operand.
function _boolean_cell_mesh(meshes::Vector{Mesh},cell::UInt64,caller)
    n=length(meshes)
    members=[i for i in 1:n if (cell & _brep_bit(i))!=0]
    isempty(members) && throw(ErrorException(
        "$caller: Boolean result cell has no member operands"))
    acc=meshes[members[1]]
    for i in members[2:end]
        acc=mesh_boolean(acc,meshes[i],:intersection)
    end
    for i in 1:n
        (cell & _brep_bit(i))==0 || continue
        acc=mesh_boolean(acc,meshes[i],:difference)
    end
    return acc
end

# The result surface of any Boolean record: the binary snapshot path, or the
# per-piece composition for multi-operand results. Multi-component fuse
# pieces extract their own connected components via `boolean_components`.
function _boolean_result_surface(m::GeoModel,t::Int,
                                 caller::AbstractString)
    spec=m.booleans[t]
    rec=m.boolean_operands[t]
    surface=if rec isa Tuple{Mesh,Mesh}
        mesh_boolean(rec[1],rec[2],spec.op)
    else
        if rec.cell===nothing
            foldl((a,b)->mesh_boolean(a,b,:union),rec.meshes)
        else
            _boolean_cell_mesh(rec.meshes,rec.cell,caller)
        end
    end
    haskey(m.boolean_components,t) || return surface
    return _boolean_component_surface(m,t,surface,caller)
end
