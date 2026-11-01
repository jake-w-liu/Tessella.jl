# ── trimmed NURBS faces: uv-chart constrained triangulation + on-surface lift ─
#
# A closed shell carrying at least one ADVANCED_FACE on a
# B_SPLINE_SURFACE_WITH_KNOTS (direct or complex-rational form) with the
# remaining faces PLANEs covers a "NURBS-trimmed" solid — a free-form cap or
# wall bounded by sampled trims — when:
#   * every geometric edge is shared by exactly two face uses — LINE edges by
#     canonical endpoint pair, degree-1 B-splines by their control sub-chain,
#     higher-degree open B-splines by the whole-domain sampled chain of
#     `_step_spline_chain`, and closed B-spline rim loops by the cyclic
#     whole-domain chain;
#   * every PLANE face's sampled boundary points lie in its plane;
#   * each NURBS face's boundary chains project onto its patch — every chain
#     point lifts back onto the surface within tolerance — and its unwrapped
#     uv wires form simple polygons (the nonlinear uv chart).
#
# Each NURBS face is meshed on its own uv chart: boundary samples seed a
# constrained Delaunay domain, a metric-scaled uv lattice fills the interior,
# boundary vertices keep the canonical 3-D positions they share with the
# neighbouring faces' sampled chains and interior/Steiner vertices lift
# through nurbs_eval. A geometrically periodic parametric direction is
# unwrapped along the wire so a tube wall triangulates as one uv strip whose
# two seam columns weld positionally — the torus-unwrapping precedent from
# the OCC curved-surface producer. Positional dedup welds every shared chain
# bitwise, then the directed shared-edge audit, a positive signed volume and
# `validate` certify the closed shell.

const _BREP_NURBS_TRIM_SUBDIV=8        # sampled segments per nonzero knot span
const _BREP_NURBS_TRIM_MAXSEED=4096    # interior uv lattice cap (OCC parity)

# ── NURBS patch with exact first derivatives ─────────────────────────────────
# Homogeneous nets let one De Boor evaluator cover both position and the two
# partial-derivative control nets of a (possibly rational) tensor surface.
struct _BrepNurbsPatch
    s::NURBSSurface
    net::Matrix{NTuple{4,Float64}}      # w-scaled homogeneous controls n_u×n_v
    dunet::Matrix{NTuple{4,Float64}}    # u-derivative net (n_u-1)×n_v, degree pu-1
    dvnet::Matrix{NTuple{4,Float64}}    # v-derivative net n_u×(n_v-1), degree pv-1
    u_lo::Float64; u_hi::Float64
    v_lo::Float64; v_hi::Float64
end

function _brep_nurbs_patch(s::NURBSSurface)
    wscale=maximum(s.weights)
    nu,nv=size(s.controls); pu,pv=s.degree_u,s.degree_v
    net=Matrix{NTuple{4,Float64}}(undef,nu,nv)
    @inbounds for j in 1:nv,i in 1:nu
        w=s.weights[i,j]/wscale;c=s.controls[i,j]
        net[i,j]=(w*c[1],w*c[2],w*c[3],w)
    end
    U,V=s.knots_u,s.knots_v
    dunet=Matrix{NTuple{4,Float64}}(undef,nu-1,nv)
    @inbounds for j in 1:nv,i in 1:nu-1
        den=U[i+pu+1]-U[i+1];f=pu/den
        a=net[i,j];b=net[i+1,j]
        dunet[i,j]=(f*(b[1]-a[1]),f*(b[2]-a[2]),f*(b[3]-a[3]),f*(b[4]-a[4]))
    end
    dvnet=Matrix{NTuple{4,Float64}}(undef,nu,nv-1)
    @inbounds for j in 1:nv-1,i in 1:nu
        den=V[j+pv+1]-V[j+1];f=pv/den
        a=net[i,j];b=net[i,j+1]
        dvnet[i,j]=(f*(b[1]-a[1]),f*(b[2]-a[2]),f*(b[3]-a[3]),f*(b[4]-a[4]))
    end
    return _BrepNurbsPatch(s,net,dunet,dvnet,
                           U[pu+1],U[end-pu],V[pv+1],V[end-pv])
end

# 4D de Boor on a homogeneous control row — same recursion as NURBS._deboor_h
# with an inline span lookup so degree-0 rows (derivative nets) also work.
function _brep_deboor4(p::Int,U::Vector{Float64},
                       Pw::AbstractVector{NTuple{4,Float64}},u::Float64)
    n=length(U)-p-2
    i=u>=U[end-p] ? n : begin
        low=p;high=length(U)-p-1;mid=(low+high)÷2
        while u<U[mid+1] || u>=U[mid+2]
            u<U[mid+1] ? (high=mid) : (low=mid)
            mid=(low+high)÷2
        end
        mid
    end
    d=Vector{NTuple{4,Float64}}(undef,p+1)
    @inbounds for j in 0:p
        d[j+1]=Pw[i-p+j+1]
    end
    @inbounds for r in 1:p, j in p:-1:r
        ulo=U[i-p+j+1];uhi=U[i+j-r+2];den=uhi-ulo
        den>0 || return d[p+1]
        α=(u-ulo)/den;a=d[j];b=d[j+1]
        d[j+1]=((1-α)*a[1]+α*b[1],(1-α)*a[2]+α*b[2],
                (1-α)*a[3]+α*b[3],(1-α)*a[4]+α*b[4])
    end
    return d[p+1]
end

# Point + both partials of the tensor patch at (u,v); parameters clamp into the
# active knot interval. Rational patches differentiate by the quotient rule.
function _brep_nurbs_d1(patch::_BrepNurbsPatch,u::Float64,v::Float64)
    s=patch.s;pu,pv=s.degree_u,s.degree_v
    U,V=s.knots_u,s.knots_v
    nu,nv=size(s.controls)
    uu=clamp(u,patch.u_lo,patch.u_hi);vv=clamp(v,patch.v_lo,patch.v_hi)
    tmp=Vector{NTuple{4,Float64}}(undef,nu)
    @inbounds for i in 1:nu
        tmp[i]=_brep_deboor4(pv,V,[patch.net[i,j] for j in 1:nv],vv)
    end
    H=_brep_deboor4(pu,U,tmp,uu)
    Ud=U[2:end-1];Vd=V[2:end-1]
    tmpu=Vector{NTuple{4,Float64}}(undef,nu-1)
    @inbounds for i in 1:nu-1
        tmpu[i]=_brep_deboor4(pv,V,[patch.dunet[i,j] for j in 1:nv],vv)
    end
    Hu=_brep_deboor4(pu-1,Ud,tmpu,uu)
    tmpv=Vector{NTuple{4,Float64}}(undef,nu)
    @inbounds for i in 1:nu
        tmpv[i]=_brep_deboor4(pv-1,Vd,[patch.dvnet[i,j] for j in 1:nv-1],vv)
    end
    Hv=_brep_deboor4(pu,U,tmpv,uu)
    w=H[4]
    (isfinite(w) && w>0) || throw(ArgumentError(
        "trimmed NURBS face: non-positive homogeneous weight at parameter"))
    P=(H[1]/w,H[2]/w,H[3]/w)
    du=((Hu[1]-P[1]*Hu[4])/w,(Hu[2]-P[2]*Hu[4])/w,(Hu[3]-P[3]*Hu[4])/w)
    dv=((Hv[1]-P[1]*Hv[4])/w,(Hv[2]-P[2]*Hv[4])/w,(Hv[3]-P[3]*Hv[4])/w)
    return P,du,dv
end

_brep_lift(patch,u,v)=nurbs_eval(patch.s,u,v)

# Gauss–Newton on the normal equations of the closest-point problem. The
# boundary chain points are authored on the surface, so the residual is tiny
# and convergence is effectively quadratic. Returns (u,v,converged).
function _brep_nurbs_newton(patch::_BrepNurbsPatch,
                            P::NTuple{3,Float64},u::Float64,v::Float64)
    uspan=patch.u_hi-patch.u_lo;vspan=patch.v_hi-patch.v_lo
    for _ in 1:30
        S,du,dv=_brep_nurbs_d1(patch,u,v)
        r=_brep_vsub(S,P)
        f1=_brep_vdot(r,du);f2=_brep_vdot(r,dv)
        g11=_brep_vdot(du,du);g12=_brep_vdot(du,dv);g22=_brep_vdot(dv,dv)
        det=g11*g22-g12*g12
        (isfinite(det) && det>0) || return (u,v,false)
        d1=(f2*g12-f1*g22)/det;d2=(f1*g12-f2*g11)/det
        (isfinite(d1) && isfinite(d2)) || return (u,v,false)
        u2=clamp(u+d1,patch.u_lo,patch.u_hi)
        v2=clamp(v+d2,patch.v_lo,patch.v_hi)
        step=hypot(u2-u,v2-v)
        u,v=u2,v2
        step<=1e-13*(uspan+vspan) && return (u,v,true)
    end
    return (u,v,false)
end

# Project one 3-D point onto the patch: hint first (wire-continuation), then a
# knot-span midpoint grid; keep the converged root nearest the target. The
# caller still enforces the on-surface residual — nothing returned here is
# trusted without it.
function _brep_nurbs_project(patch::_BrepNurbsPatch,P::NTuple{3,Float64},
                             hint)::Union{NTuple{2,Float64},Nothing}
    seeds=NTuple{2,Float64}[]
    hint!==nothing && push!(seeds,hint)
    pu,pv=patch.s.degree_u,patch.s.degree_v
    U,V=patch.s.knots_u,patch.s.knots_v
    function mids(knots,p)
        out=Float64[]
        for i in p+1:length(knots)-p-1
            knots[i+1]>knots[i] && push!(out,(knots[i]+knots[i+1])/2)
        end
        return out
    end
    for a in mids(U,pu), b in mids(V,pv)
        push!(seeds,(a,b))
    end
    for a in (patch.u_lo,patch.u_hi), b in (patch.v_lo,patch.v_hi)
        push!(seeds,(a,b))
    end
    best=nothing;bestd=Inf
    for (a,b) in seeds
        u,v,ok=_brep_nurbs_newton(patch,P,a,b)
        ok || continue
        q=_brep_lift(patch,u,v)
        d=_brep_vdot(_brep_vsub(q,P),_brep_vsub(q,P))
        d<bestd && (best=(u,v);bestd=d)
    end
    return best
end

# Geometric seam: opposite parametric rows coincide pointwise → period, else 0.
function _brep_nurbs_seam(patch::_BrepNurbsPatch,tol::Float64)
    uper=vper=0.0
    nv=24
    seam=true
    for k in 0:nv
        v=patch.v_lo+(patch.v_hi-patch.v_lo)*k/nv
        _nurbs_pt_eq(_brep_lift(patch,patch.u_lo,v),
                     _brep_lift(patch,patch.u_hi,v),tol) || (seam=false;break)
    end
    seam && (uper=patch.u_hi-patch.u_lo)
    nu=24
    seam=true
    for k in 0:nu
        u=patch.u_lo+(patch.u_hi-patch.u_lo)*k/nu
        _nurbs_pt_eq(_brep_lift(patch,u,patch.v_lo),
                     _brep_lift(patch,u,patch.v_hi),tol) || (seam=false;break)
    end
    seam && (vper=patch.v_hi-patch.v_lo)
    return uper,vper
end

# ── parametric face charts: NURBS patch, cylinder, sphere ────────────────────
# A "chart" is one face's parametric domain: the trimmed-NURBS certificate
# triangulates every curved face in its own (u,v) coordinates with the same
# boundary-project → period-unwrap → CDT → on-surface-lift pipeline. The
# NURBS chart wraps the Gauss–Newton machinery above; the analytic charts
# project in closed form — a cylinder maps u to the axial angle and v to the
# axial coordinate, a sphere maps (u,v) to (longitude, latitude) in the
# placement's own (x,y,z) frame. Both are u-periodic with period 2π; the
# sphere additionally carries collapsed poles at v=±π/2.
struct _BrepNurbsChart
    pk::_BrepNurbsPatch
end
struct _BrepCylChart
    c0::NTuple{3,Float64}; z::NTuple{3,Float64}
    x::NTuple{3,Float64};  y::NTuple{3,Float64}; R::Float64
end
struct _BrepSphChart
    c::NTuple{3,Float64};  z::NTuple{3,Float64}
    x::NTuple{3,Float64};  y::NTuple{3,Float64}; R::Float64
end

_chart_lift(c::_BrepNurbsChart,u,v)=_brep_lift(c.pk,u,v)
function _chart_lift(c::_BrepCylChart,u,v)
    cu,su=cos(u),sin(u)
    return ntuple(ax->c.c0[ax]+v*c.z[ax]+c.R*(cu*c.x[ax]+su*c.y[ax]),3)
end
function _chart_lift(c::_BrepSphChart,u,v)
    cu,su=cos(u),sin(u);cv,sv=cos(v),sin(v)
    return ntuple(ax->c.c[ax]+c.R*(cv*(cu*c.x[ax]+su*c.y[ax])+sv*c.z[ax]),3)
end

# (point, ∂S/∂u, ∂S/∂v) — the interior-lattice metric samples.
_chart_d1(c::_BrepNurbsChart,u,v)=_brep_nurbs_d1(c.pk,u,v)
function _chart_d1(c::_BrepCylChart,u,v)
    cu,su=cos(u),sin(u)
    du=ntuple(ax->c.R*(-su*c.x[ax]+cu*c.y[ax]),3)
    return _chart_lift(c,u,v),du,c.z
end
function _chart_d1(c::_BrepSphChart,u,v)
    cu,su=cos(u),sin(u);cv,sv=cos(v),sin(v)
    du=ntuple(ax->c.R*cv*(-su*c.x[ax]+cu*c.y[ax]),3)
    dv=ntuple(ax->c.R*(-sv*(cu*c.x[ax]+su*c.y[ax])+cv*c.z[ax]),3)
    return _chart_lift(c,u,v),du,dv
end

# Closed-form projections for the analytic charts return exact (u,v); the
# caller's lift check supplies the on-surface certificate.
_chart_project(c::_BrepNurbsChart,P,hint)=_brep_nurbs_project(c.pk,P,hint)
function _chart_project(c::_BrepCylChart,P,hint)
    d=_brep_vsub(P,c.c0)
    return (atan(_brep_vdot(d,c.y),_brep_vdot(d,c.x)),_brep_vdot(d,c.z))
end
function _chart_project(c::_BrepSphChart,P,hint)
    d=_brep_vsub(P,c.c)
    return (atan(_brep_vdot(d,c.y),_brep_vdot(d,c.x)),
            atan(_brep_vdot(d,c.z),
                 hypot(_brep_vdot(d,c.x),_brep_vdot(d,c.y))))
end

_chart_periods(c::_BrepNurbsChart,stol)=_brep_nurbs_seam(c.pk,stol)
_chart_periods(c::_BrepCylChart,stol)=(2π,0.0)
_chart_periods(c::_BrepSphChart,stol)=(2π,0.0)

# uv length scale for coincidence tolerances.
_chart_uscale(c::_BrepNurbsChart)=
    max(c.pk.u_hi-c.pk.u_lo,c.pk.v_hi-c.pk.v_lo)
_chart_uscale(c::_BrepCylChart)=2π
_chart_uscale(c::_BrepSphChart)=2π

# A pole vertex maps to a chart singularity where the u coordinate is free:
# snap v to the canonical pole and continue the wire's u coordinate so the
# degenerate edge collapses to a zero-length uv segment.
function _chart_pole_snap(c::_BrepSphChart,p,u,v,prev,stol)
    d=_brep_vsub(p,c.c)
    ρ=hypot(_brep_vdot(d,c.x),_brep_vdot(d,c.y))
    ρ<=stol || return (u,v)
    vv=_brep_vdot(d,c.z)>=0.0 ? π/2 : -π/2
    uu=prev===nothing ? u : prev[1]
    return (uu,vv)
end
_chart_pole_snap(c,p,u,v,prev,stol)=(u,v)

# Collapsed-boundary v parameters for a bare closed wire wrapping the u
# period: (v_lo_pole, v_hi_pole) — nothing where the chart has no pole on
# that side. A wire spanning the full period without returning to its start
# is closed through the pole on its left side (the STEP outer-bound
# convention: the face domain lies left of the traversed wire). For the
# NURBS chart a pole side is a boundary control row collapsed to a single
# physical position (clamped knots evaluate it to that point).
_chart_pole_ends(c::_BrepSphChart)=(-π/2,π/2)
function _chart_pole_ends(c::_BrepNurbsChart)
    ctrls=c.pk.s.controls
    tol=1e-9*max(_poly_scale([p for p in ctrls]),1.0)
    lo=hi=nothing
    all(i->_nurbs_pt_eq(ctrls[i,1],ctrls[1,1],tol),
        axes(ctrls,1)) && (lo=c.pk.v_lo)
    all(i->_nurbs_pt_eq(ctrls[i,end],ctrls[1,end],tol),
        axes(ctrls,1)) && (hi=c.pk.v_hi)
    return (lo,hi)
end
_chart_pole_ends(c)=(nothing,nothing)

# Base-domain coordinate of an unwrapped sheet coordinate: a periodic
# direction is reduced modulo its period before lifting — the NURBS de Boor
# clamps out-of-domain parameters, which would fold every unwrapped sheet
# copy onto the base boundary instead of wrapping onto the seam's twin. The
# analytic charts lift periodically anyway, so their wrap is a no-op
# positionally but keeps the same call shape.
_chart_base_u(c::_BrepNurbsChart,u,per)=
    per>0.0 ? c.pk.u_lo+mod(u-c.pk.u_lo,per) : u
_chart_base_v(c::_BrepNurbsChart,v,per)=
    per>0.0 ? c.pk.v_lo+mod(v-c.pk.v_lo,per) : v
_chart_base_u(c,u,per)=per>0.0 ? mod(u,per) : u
_chart_base_v(c,v,per)=per>0.0 ? mod(v,per) : v
_chart_lift_unwrapped(c,u,v,uper,vper)=_chart_lift(c,
    _chart_base_u(c,u,uper),_chart_base_v(c,v,vper))
_chart_d1_unwrapped(c,u,v,uper,vper)=_chart_d1(c,
    _chart_base_u(c,u,uper),_chart_base_v(c,v,vper))

# Canonical whole-domain cyclic chain of a closed B-spline edge: the curve's
# parameter range, `_BREP_NURBS_TRIM_SUBDIV` segments per nonzero span, with
# both evaluated endpoints pinned on the shared vertex. Returns the sampled
# positions in parameter order (no closing duplicate) or nothing for open
# edges, non-splines or partial-domain uses.
function _step_spline_chain_closed(entities,use)
    c=_step_nurbs_curve(use.base,entities)
    c===nothing && return nothing
    ec=use.ec
    length(ec.args)>=3 || throw(ArgumentError(
        "trimmed NURBS shell: EDGE_CURVE #$(ec.id) is truncated"))
    vs=_step_vertex_point(entities,ec.args[2])
    ve=_step_vertex_point(entities,ec.args[3])
    (vs===nothing || ve===nothing) && throw(ArgumentError(
        "trimmed NURBS shell: EDGE_CURVE #$(ec.id) has unresolvable "*
        "vertex geometry"))
    vs==ve || return nothing
    U=c.knots
    u0=U[c.degree+1];u1=U[end-c.degree]
    u1-u0>0 || return nothing
    p0=nurbs_eval(c,u0);p1=nurbs_eval(c,u1)
    tol=1e-9*max(_poly_scale(c.controls),1.0)
    (_posnear(vs,p0,tol) && _posnear(vs,p1,tol)) || return nothing
    spans=count(k->U[k+1]>U[k],c.degree+1:length(U)-c.degree-1)
    m=max(1,_BREP_NURBS_TRIM_SUBDIV*spans)
    cyc=NTuple{3,Float64}[nurbs_eval(c,u0+(u1-u0)*k/m) for k in 0:m-1]
    cyc[1]=vs
    return cyc
end

# Canonical sampled chain of a circular arc edge: the use's own [start,stop]
# window on `_BREP_PRISM_ARC_SEGMENTS` segments per full turn — the same
# sampling rule the prismatic certificate applies, so two uses of the same
# edge entity emit bitwise-identical interior nodes and weld positionally.
# Parameter order runs start→stop along the circle frame; traversal reversal
# happens in `_nurbs_trim_chain`.
function _step_trim_arc_chain(a)
    lo=atan(_brep_vdot(_brep_vsub(a.start,a.center),a.y),
            _brep_vdot(_brep_vsub(a.start,a.center),a.x))
    m=max(1,ceil(Int,_BREP_PRISM_ARC_SEGMENTS*a.span/(2π)))
    pts=NTuple{3,Float64}[ntuple(ax->a.center[ax]+
        a.r*(cos(lo+a.span*k/m)*a.x[ax]+sin(lo+a.span*k/m)*a.y[ax]),3)
        for k in 0:m]
    # Pin the junction ends on the shared vertex positions: evaluation can
    # sit ~1e-16 off the authored vertex, and a later weld tolerance would
    # still collapse the gap — pinning makes the shared endpoints bitwise.
    pts[1]=a.start;pts[end]=a.stop
    return pts
end

# Canonical cyclic chain of a whole-circle edge, anchored at its shared
# vertex so both uses of the same EDGE_CURVE produce identical samples.
function _step_trim_circle_cyc(cf,v1)
    θ1=atan(_brep_vdot(_brep_vsub(v1,cf.center),cf.y),
            _brep_vdot(_brep_vsub(v1,cf.center),cf.x))
    m=_BREP_PRISM_ARC_SEGMENTS
    cyc=NTuple{3,Float64}[ntuple(ax->cf.center[ax]+
        cf.r*(cos(θ1+2π*k/m)*cf.x[ax]+sin(θ1+2π*k/m)*cf.y[ax]),3)
        for k in 0:m-1]
    cyc[1]=v1
    return cyc
end

# Traversal-ordered canonical 3-D chain of one classified use (for wire rings,
# uv projection and segment incidence). Open chains run junction→junction; a
# closed-spline or whole-circle chain closes back onto its vertex —
# direction follows the ORIENTED_EDGE flag folded with the EDGE_CURVE
# same_sense flag, captured at classification into aux.pts. A degenerate
# (pole) edge contributes its single vertex and no segments.
function _nurbs_trim_chain(ur,ptol)
    if ur.kind==:line
        return NTuple{3,Float64}[ur.u.v1,ur.u.v2]
    elseif ur.kind==:closedspline || ur.kind==:circle
        return ur.aux.pts
    elseif ur.kind==:degenerate
        return NTuple{3,Float64}[ur.u.v1]
    end
    ch=ur.aux.pts
    _posnear(ur.u.v1,first(ch),ptol) && return collect(ch)
    return reverse(collect(ch))
end

# Certify and emit a closed shell that mixes trimmed B-spline faces with
# planar faces and analytic cylindrical/spherical faces. Every curved face —
# NURBS or analytic — is triangulated in its own parametric chart through the
# shared chart abstraction; a shared edge's canonical 3-D chain is emitted
# identically by both adjacent faces so segments weld positionally. Returns
# a Mesh or nothing when the inventory is not this certificate's family;
# throws on certified mismatches.
function _step_nurbs_trimmed_shell(entities,faces,curvekinds,caller)
    length(faces)>=2 || return nothing
    all(f->f.role==1,faces) || return nothing
    all(f->f.surf!==nothing,faces) || return nothing
    nurbs_i=Int[];chart_i=Int[];chart_kind=Dict{Int,Symbol}();plane_i=Int[]
    for (i,f) in enumerate(faces)
        if _step_is_nurbs_surface(f.surf)
            push!(nurbs_i,i);push!(chart_i,i);chart_kind[i]=:nurbs
        elseif f.surf.kind=="PLANE"
            push!(plane_i,i)
        elseif f.surf.kind=="CYLINDRICAL_SURFACE"
            push!(chart_i,i);chart_kind[i]=:cyl
        elseif f.surf.kind=="SPHERICAL_SURFACE"
            push!(chart_i,i);chart_kind[i]=:sph
        else
            return nothing
        end
    end
    isempty(nurbs_i) && return nothing

    # ── per-face classified uses ──────────────────────────────────────────
    fuse=Vector{Vector{Vector{NamedTuple}}}(undef,length(faces))
    for (i,f) in enumerate(faces)
        fe=_deref_entity(entities,(:ref,f.id))
        fe===nothing && throw(ArgumentError(
            "$caller: face #$(f.id) record is missing"))
        bounds=_step_face_bound_uses(entities,fe,curvekinds,caller)
        fused=Vector{NamedTuple}[]
        for uses in bounds
            out=NamedTuple[]
            for u in uses
                b=u.base.kind
                if b=="LINE"
                    if u.v1==u.v2
                        # A degenerate (pole) edge is admissible only on a
                        # chart whose domain actually collapses there —
                        # sphere/cylinder faces certify it through the pole
                        # snap below.
                        get(chart_kind,i,:none) in (:cyl,:sph) ||
                            (push!(curvekinds,
                                "degenerate bound edge — pole NURBS face");
                             return nothing)
                        push!(out,(kind=:degenerate,u=u,aux=nothing))
                    else
                        push!(out,(kind=:line,u=u,aux=nothing))
                    end
                elseif b=="CIRCLE"
                    c,nr,xr,yr,cr=_step_circle_frame_full(entities,u.base,
                                                          caller)
                    cf=(center=c,nr=nr,x=xr,y=yr,r=cr)
                    if u.v1==u.v2
                        cyc=_step_trim_circle_cyc(cf,u.v1)
                        sense=length(u.ec.args)>=5 ?
                            _brep_enum(u.ec.args[5]) : nothing
                        sense===nothing && (sense=true)
                        ch=u.fwd==sense ?
                            [cyc;cyc[1:1]] :
                            [cyc[1];reverse(cyc[2:end]);cyc[1]]
                        push!(out,(kind=:circle,u=u,
                                   aux=(cf=cf,cyc=cyc,pts=ch)))
                    else
                        a=_step_arc_span(entities,u,caller)
                        push!(out,(kind=:arc,u=u,
                                   aux=(span=a,cf=(center=a.center,nr=a.nr,
                                        x=a.x,y=a.y,r=a.r),
                                        pts=_step_trim_arc_chain(a))))
                    end
                elseif b=="B_SPLINE_CURVE_WITH_KNOTS" || b=="COMPLEX"
                    ch=_step_polyline_subchain(entities,u)
                    if ch===nothing
                        if u.v1==u.v2
                            cyc=_step_spline_chain_closed(entities,u)
                            cyc===nothing &&
                                (push!(curvekinds,b);return nothing)
                            # Same_sense .F. reverses the parameter direction
                            # of the edge relative to its curve; the wire
                            # traversal is fwd==sense.
                            sense=length(u.ec.args)>=5 ?
                                _brep_enum(u.ec.args[5]) : nothing
                            sense===nothing && (sense=true)
                            ch=u.fwd==sense ?
                                [cyc;cyc[1:1]] :
                                [cyc[1];reverse(cyc[2:end]);cyc[1]]
                            push!(out,(kind=:closedspline,u=u,
                                       aux=(cyc=cyc,pts=ch)))
                        else
                            ch=_step_spline_chain(entities,u)
                            ch===nothing &&
                                (push!(curvekinds,b);return nothing)
                            push!(out,(kind=:spline,u=u,aux=(pts=ch,)))
                        end
                    else
                        push!(out,(kind=:polyline,u=u,aux=(pts=ch,)))
                    end
                else
                    push!(curvekinds,b);return nothing
                end
            end
            push!(fused,out)
        end
        fuse[i]=fused
    end

    allv=NTuple{3,Float64}[]
    for bound in Iterators.flatten(fuse),ur in bound
        push!(allv,ur.u.v1);push!(allv,ur.u.v2)
    end
    pscale=max(_poly_scale(allv),1.0);ptol=1e-9*pscale

    # ── straight-segment incidence: every sampled edge exactly twice ──────
    vpos=NTuple{3,Float64}[];vranges=Dict{Tuple{Int,Int,Int},UnitRange{Int}}()
    for (i,bound) in enumerate(fuse),(bi,b) in enumerate(bound),
            (ui,ur) in enumerate(b)
        ch=_nurbs_trim_chain(ur,ptol)
        r0=length(vpos)+1;append!(vpos,ch)
        vranges[(i,bi,ui)]=r0:length(vpos)
    end
    vnodes,vcanon=_dedup_positions(vpos,caller)
    segcount=Dict{Tuple{Int32,Int32},Int}()
    for (key,rng) in vranges
        ids=[vcanon[j] for j in rng]
        for k in 1:length(ids)-1
            a=ids[k];b=ids[k+1]
            a==b && throw(ArgumentError(
                "$caller: face #$(faces[key[1]].id) carries a degenerate "*
                "edge segment"))
            e=a<b ? (a,b) : (b,a)
            segcount[e]=get(segcount,e,0)+1
        end
    end
    for (e,c) in segcount
        c==2 || return nothing
    end

    # ── plane faces: every sampled boundary point in-plane ────────────────
    for i in plane_i
        p0,pn=_step_plane_frame(entities,faces[i].surf,caller)
        atol=1e-9
        for bound in fuse[i],ur in bound
            ch=_nurbs_trim_chain(ur,ptol)
            for v in ch
                abs(_brep_vdot(_brep_vsub(v,p0),pn))<=ptol || throw(
                    ArgumentError("$caller: sampled edge #$(ur.u.ec.id) "*
                        "leaves PLANE #$(faces[i].surf.id)"))
            end
            if ur.kind==:arc || ur.kind==:circle
                cf=ur.aux.cf
                abs(_brep_vdot(cf.nr,pn))>=1.0-atol || throw(ArgumentError(
                    "$caller: CIRCLE #$(ur.u.base.id) is not contained in "*
                    "PLANE #$(faces[i].surf.id)"))
                abs(_brep_vdot(_brep_vsub(cf.center,p0),pn))<=ptol ||
                    throw(ArgumentError(
                        "$caller: CIRCLE #$(ur.u.base.id) center is off "*
                        "PLANE #$(faces[i].surf.id)"))
            end
        end
    end

    # ── chart faces: uv domain polygons + on-surface lift certificate ─────
    # nuv[i] = per-bound unwrapped uv vertex lists; ncan[i] = parallel
    # canonical 3-D positions. Shared-edge lifts agree bitwise: the same
    # canonical chain projects on the NURBS patch's uv chart and on the
    # analytic neighbour's parametric domain — both lifts must hit the same
    # 3-D positions within the surface tolerance.
    nuv=Dict{Int,Vector{Vector{NTuple{2,Float64}}}}()
    ncan=Dict{Int,Vector{Vector{NTuple{3,Float64}}}}()
    charts=Dict{Int,Any}()
    periods=Dict{Int,Tuple{Float64,Float64}}()
    for i in chart_i
        f=faces[i]
        if chart_kind[i]==:nurbs
            nurbs=try
                _step_nurbs_surface(f.surf,entities,caller)
            catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError(
                    "$caller: face #$(f.id) B-spline surface is malformed: "*
                    sprint(showerror,err)))
            end
            nurbs isa NURBSSurface || (push!(curvekinds,
                "face surface $(_step_surface_name(f.surf))");return nothing)
            pk=_brep_nurbs_patch(nurbs)
            charts[i]=_BrepNurbsChart(pk)
            stol=1e-9*max(_poly_scale([nurbs.controls...]),1.0)
        elseif chart_kind[i]==:cyl
            c0,z,x,y,R=_step_cylinder_frame(entities,f.surf,caller)
            charts[i]=_BrepCylChart(c0,z,x,y,R)
            stol=1e-9*max(R,1.0)
        else
            c,z,R=_step_sphere_params(entities,f.surf,caller)
            x=_brep_default_x(z,caller);y=_brep_vcross(z,x)
            charts[i]=_BrepSphChart(c,z,x,y,R)
            stol=1e-9*max(R,1.0)
        end
        chart=charts[i]
        uper,vper=_chart_periods(chart,stol)
        periods[i]=(uper,vper)
        bounds_uv=Vector{NTuple{2,Float64}}[]
        bounds_p3=Vector{NTuple{3,Float64}}[]
        uvtol=1e-7*_chart_uscale(chart)
        for bound in fuse[i]
            uvs=NTuple{2,Float64}[]
            p3s=NTuple{3,Float64}[]
            prev=nothing
            prev3=nothing
            for (ui,ur) in enumerate(bound)
                ch=_nurbs_trim_chain(ur,ptol)
                # Wire chaining: the use's first chain point is the previous
                # use's junction vertex and is only re-emitted as its own
                # last point; the first use contributes its start vertex.
                prev3===nothing ||
                    _posnear(prev3,first(ch),ptol) || throw(ArgumentError(
                        "$caller: face #$(f.id) boundary edges do not "*
                        "form an ordered closed chain"))
                # An authored PCURVE on this face's surface supplies the uv
                # wire directly — the 3-D chain is only inverted onto the
                # authored parametric curve, never re-projected onto the
                # surface (see `_step_use_pcurve`).
                pc=_step_use_pcurve(entities,ur.u,f.surf.id,caller)
                pc_t=nothing
                for k in (isempty(uvs) ? 1 : 2):length(ch)
                    p=ch[k]
                    uv=_chart_project(chart,p,prev===nothing ?
                        nothing : (prev[1],prev[2]))
                    uv===nothing && throw(ArgumentError(
                        "$caller: edge #$(ur.u.ec.id) of face #$(f.id) "*
                        "does not project onto its "*
                        "$(chart_kind[i]==:nurbs ? "B-spline" :
                            f.surf.kind=="CYLINDRICAL_SURFACE" ?
                            "cylindrical" : "spherical") surface"))
                    u0,v0=uv
                    if pc!==nothing
                        pc_t=_pc2_param(pc,(u0,v0),pc_t,caller,
                                        ur.u.ec.id,f.id)
                        pc_t===nothing && throw(ArgumentError(
                            "$caller: edge #$(ur.u.ec.id) does not follow "*
                            "its PCURVE on face #$(f.id)"))
                        u0,v0=_pc2_eval(pc,pc_t)
                    end
                    # On-surface certificate on the in-domain projection
                    # first — the unwrapped sheet coordinate below is only
                    # for the uv polygon (clamped lift would alias copies).
                    q=_chart_lift(chart,u0,v0)
                    _nurbs_pt_eq(q,p,stol) || throw(ArgumentError(
                        "$caller: edge #$(ur.u.ec.id) leaves the "*
                        "$(chart_kind[i]==:nurbs ? "B-spline" :
                            f.surf.kind=="CYLINDRICAL_SURFACE" ?
                            "cylindrical" : "spherical") surface of "*
                        "face #$(f.id)"))
                    u0,v0=_chart_pole_snap(chart,p,u0,v0,prev,stol)
                    uper>0.0 && prev!==nothing &&
                        (u0+=uper*round((prev[1]-u0)/uper))
                    vper>0.0 && prev!==nothing &&
                        (v0+=vper*round((prev[2]-v0)/vper))
                    uvc=(u0,v0)
                    prev=uvc;prev3=last(ch)
                    push!(uvs,uvc);push!(p3s,p)
                end
            end
            # The wire closes on its first vertex — possibly one period away
            # in the unwrapped chart; the modulo-segment below closes it.
            n=length(uvs)
            n>=3 || throw(ArgumentError(
                "$caller: face #$(f.id) wire has fewer than three "*
                "projected vertices"))
            du0=uvs[n][1]-uvs[1][1];dv0=abs(uvs[n][2]-uvs[1][2])
            lo_p,hi_p=_chart_pole_ends(chart)
            # A bare single-use closed wire ending a full u period from its
            # start (a latitude rim on a sphere face, or the equator rim of a
            # NURBS dome whose v boundary row collapses to a pole) is closed
            # through the collapsed pole on its left side: the kept period
            # endpoint makes the polygon span the whole period, and the two
            # appended uv vertices carry the pole's 3-D position so the
            # boundary degenerates exactly onto the physical pole.
            period_dup=uper>0.0 && dv0<=uvtol &&
                abs(abs(du0)-round(abs(du0)/uper)*uper)<=uvtol
            bare=period_dup && length(bound)==1 &&
                bound[1].kind in (:circle,:closedspline)
            vp=bare ? (du0>0.0 ? hi_p : lo_p) : nothing
            if vp!==nothing
                p_pole=_chart_lift(chart,0.0,vp)
                push!(uvs,(uvs[n][1],vp));push!(p3s,p_pole)
                push!(uvs,(uvs[1][1],vp));push!(p3s,p_pole)
                n+=2
            elseif abs(du0)<=uvtol && dv0<=uvtol
                pop!(uvs);pop!(p3s);n-=1
            elseif period_dup &&
                all(ur->ur.kind!=:degenerate,bound)
                # The last vertex may sit a full period from the first —
                # identical uv position, distinct sheet. With an authored
                # degenerate pole edge in the bound, though, the closing
                # segment is that edge's uv image and the last point must be
                # kept as the real closure.
                pop!(uvs);pop!(p3s);n-=1
            end
            n>=3 || throw(ArgumentError(
                "$caller: face #$(f.id) wire has fewer than three "*
                "projected vertices"))
            push!(bounds_uv,uvs);push!(bounds_p3,p3s)
        end
        nuv[i]=bounds_uv;ncan[i]=bounds_p3
    end

    # ── uv polygon audit: simple closed wires, one unique outer ──────────
    for i in chart_i
        wires=nuv[i]
        length(wires)>=1 || return nothing
        if length(wires)>1
            # A unique containment certificate picks the outer wire, as in
            # the planar-hole path; holes must oppose its winding.
            candidates=[k for k in eachindex(wires) if
                all(l->l==k || _brep_point_in_wire(wires[l][1],wires[k]),
                    eachindex(wires))]
            length(candidates)==1 || return nothing
            order=vcat(candidates[1],
                       setdiff(eachindex(wires),candidates))
            wires=wires[order]
            ncan[i]=ncan[i][order];nuv[i]=wires
            a1=_signed_area2(first.(wires[1]),last.(wires[1]))
            for k in 2:length(wires)
                ak=_signed_area2(first.(wires[k]),last.(wires[k]))
                signbit(ak)==signbit(a1) && (reverse!(wires[k]);
                                            reverse!(ncan[i][k]))
            end
        end
        for (k,wire) in enumerate(wires)
            nw=length(wire)
            for a in 1:nw, b in 1:nw
                a==b && continue
                # non-consecutive segment pairs must not cross or touch
                (a+1==b || b+1==a || (a==1 && b==nw) || (b==1 && a==nw)) &&
                    continue
                _brep_wires_intersect(wire[a],wire[mod1(a+1,nw)],
                    wire[b],wire[mod1(b+1,nw)]) && throw(ArgumentError(
                    "$caller: NURBS face #$(faces[i].id) uv boundary "*
                    "self-intersects"))
            end
        end
    end

    # ── emission ──────────────────────────────────────────────────────────
    pool=NTuple{3,Float64}[]
    items=Tuple{Symbol,UnitRange{Int},Int}[]
    emit!(kind,pts,fi)=(push!(items,(kind,
        length(pool)+1:length(pool)+length(pts),fi));append!(pool,pts);nothing)

    # Traversal-ordered ring for one wire of a plane face.
    function trim_ring(bound)
        ring=NTuple{3,Float64}[]
        if length(bound)==1 && bound[1].kind in (:closedspline,:circle)
            append!(ring,bound[1].aux.cyc)
            return ring
        end
        for ur in bound
            ur.kind in (:closedspline,:circle) && throw(ArgumentError(
                "$caller: closed $(ur.u.base.kind) edge #$(ur.u.ec.id) "*
                "shares a wire with other edges"))
            ch=_nurbs_trim_chain(ur,ptol)
            append!(ring,@view ch[2:end])
        end
        return ring
    end

    for i in plane_i
        rings=[trim_ring(bound) for bound in fuse[i]]
        any(r->length(r)<3,rings) && throw(ArgumentError(
            "$caller: plane face #$(faces[i].id) has a degenerate wire"))
        if length(rings)>1
            for tri in _brep_planar_face_holes(rings,caller)
                emit!(:tri,tri,i)
            end
        else
            emit!(:ring,only(rings),i)
        end
    end

    for i in chart_i
        chart=charts[i]
        uper_e,vper_e=periods[i]
        xs=Float64[];ys=Float64[]
        lift=NTuple{3,Float64}[]
        segs=Tuple{Int,Int}[]
        index=Dict{NTuple{2,Float64},Int}()
        function emituv(uv,p3)
            get!(index,uv) do
                push!(xs,uv[1]);push!(ys,uv[2]);push!(lift,p3)
                return length(xs)
            end
        end
        for (uvs,p3s) in zip(nuv[i],ncan[i])
            loop_idx=[emituv(uvs[k],p3s[k]) for k in eachindex(uvs)]
            nloop=length(loop_idx)
            for k in 1:nloop
                a,b=loop_idx[k],loop_idx[mod1(k+1,nloop)]
                a==b || (xs[a]==xs[b] && ys[a]==ys[b]) ||
                    push!(segs,(a,b))
            end
        end
        isempty(segs) && throw(ArgumentError(
            "$caller: face #$(faces[i].id) has an empty uv boundary"))
        # Metric-scaled interior lattice at the boundary's 3-D spacing —
        # the OCC curved-surface seeding rule, on the same footprint caps.
        umin=minimum(xs);umax=maximum(xs);vmin=minimum(ys);vmax=maximum(ys)
        uspan=umax-umin;vspan=vmax-vmin
        if uspan>0.0 && vspan>0.0
            h3u=0.0;nu_e=0;h3v=0.0;nv_e=0
            for (a,b) in segs
                pa=lift[a];pb=lift[b]
                l3=sqrt((pb[1]-pa[1])^2+(pb[2]-pa[2])^2+(pb[3]-pa[3])^2)
                if abs(xs[b]-xs[a])>=abs(ys[b]-ys[a])
                    h3u+=l3;nu_e+=1
                else
                    h3v+=l3;nv_e+=1
                end
            end
            h3=(h3u+h3v)/max(nu_e+nv_e,1)
            h3u=nu_e>0 && h3u>0.0 ? h3u/nu_e : h3
            h3v=nv_e>0 && h3v>0.0 ? h3v/nv_e : h3
            su=0.0;sv=0.0
            for fu in (0.25,0.5,0.75),fv in (0.25,0.5,0.75)
                _,du_c,dv_c=_chart_d1_unwrapped(chart,
                    umin+fu*uspan,vmin+fv*vspan,uper_e,vper_e)
                su=max(su,sqrt(_brep_vdot(du_c,du_c)))
                sv=max(sv,sqrt(_brep_vdot(dv_c,dv_c)))
            end
            du_t=h3u/max(su,1e-12);dv_t=h3v/max(sv,1e-12)
            nu=du_t>0.0 && isfinite(uspan/du_t) ?
                clamp(floor(Int,uspan/du_t)-1,0,64) : 0
            nv=dv_t>0.0 && isfinite(vspan/dv_t) ?
                clamp(floor(Int,vspan/dv_t)-1,0,64) : 0
            while nu>1 && nv>1 && nu*nv>_BREP_NURBS_TRIM_MAXSEED
                nu>=nv ? (nu=(nu+1)÷2) : (nv=(nv+1)÷2)
            end
            if nu>0 && nv>0 && nu*nv<=_BREP_NURBS_TRIM_MAXSEED
                du=uspan/(nu+1);dv=vspan/(nv+1)
                tol2=(0.35*min(du,dv))^2
                outer=nuv[i][1]
                holes=nuv[i][2:end]
                for iu in 1:nu,jv in 1:nv
                    u=umin+iu*du;v=vmin+jv*dv
                    _brep_point_in_wire((u,v),outer) || continue
                    any(w->_brep_point_in_wire((u,v),w),holes) && continue
                    near=false
                    for k in eachindex(xs)
                        (xs[k]-u)^2+(ys[k]-v)^2<tol2 && (near=true;break)
                    end
                    near && continue
                    emituv((u,v),
                        _chart_lift_unwrapped(chart,u,v,uper_e,vper_e))
                end
            end
        end
        T=Mesh2D.constrained_delaunay(xs,ys,segs)
        interior=Mesh2D.classify_interior(T)
        mesh=Mesh2D.to_mesh(T;interior)
        canonical=Dict{NTuple{2,Float64},NTuple{3,Float64}}()
        for k in eachindex(lift)
            canonical[(xs[k],ys[k])]=lift[k]
        end
        uvcoords=mesh.coords
        cells=NTuple{3,NTuple{3,Float64}}[]
        for c in axes(mesh.tris,2)
            tri=NTuple{3,Float64}[]
            for row in 1:3
                u=uvcoords[1,mesh.tris[row,c]];v=uvcoords[2,mesh.tris[row,c]]
                p3=get(canonical,(u,v),nothing)
                p3===nothing && (p3=_chart_lift_unwrapped(chart,u,v,
                    uper_e,vper_e))
                push!(tri,p3)
            end
            # Lifted weld-collapse: a cell whose two vertices lift within
            # the positional weld tolerance of one 3-D node (pole fan,
            # degenerate-edge image or a folded domain corner) drops out —
            # the pool dedup below would weld them into a self-edge anyway.
            _posnear(tri[1],tri[2],ptol) || _posnear(tri[2],tri[3],ptol) ||
                _posnear(tri[3],tri[1],ptol) ||
                push!(cells,(tri[1],tri[2],tri[3]))
        end
        isempty(cells) && throw(ArgumentError(
            "$caller: face #$(faces[i].id) triangulates to nothing"))
        for tri in cells
            emit!(:tri,tri,i)
        end
    end

    nodes,canon=_dedup_positions(pool,caller)
    itris=NTuple{3,Int32}[];tface=Int[]
    for (kind,rng,fi) in items
        if kind==:ring
            ir=Int32[];prev=Int32(0)
            for j in rng
                v=canon[j]
                v!=prev && push!(ir,v);prev=v
            end
            length(ir)>1 && first(ir)==last(ir) && pop!(ir)
            n0=length(itris)
            _ear_clip!(itris,ir,nodes,caller)
            append!(tface,fill(fi,length(itris)-n0))
        else
            a,b,c=(Int32(canon[j]) for j in rng)
            push!(itris,(a,b,c));push!(tface,fi)
        end
    end
    _checked_brep_mul(3,length(itris),caller,"NURBS-trimmed triangle count")
    length(itris)<=_BREP_MAX_SURFACE_TRIS || throw(ArgumentError(
        "$caller: NURBS-trimmed shell exceeds the surface triangle limit"))
    tris,edge_use=_step_shell_flip_parity(itris,tface,length(faces),
        "NURBS-trimmed shell",caller)
    return _step_shell_audit(tris,nodes,tface,length(faces),edge_use,
        "NURBS-trimmed shell",caller)
end

# ── STEP PCURVE trims ────────────────────────────────────────────────────────
# An EDGE_CURVE's geometry may be authored as a SURFACE_CURVE carrying one
# PCURVE per adjacent face (the curve-on-surface representation). The
# PCURVE's DEFINITIONAL_REPRESENTATION holds a 2-D LINE or B-spline curve in
# the face surface's own uv coordinates. When a pcurve is authored for this
# face's surface, its evaluated uv positions supply the trim wire directly —
# the canonical 3-D chain is inverted onto the authored parametric curve
# (seeded by the face chart's projection) instead of being re-projected onto
# the surface; the shared on-surface lift certificate still verifies the
# authored curve against the canonical 3-D positions, so a pcurve that
# drifts off the surface fails loudly rather than silently falling back.
struct _Pcurve2d
    kind::Symbol            # :line | :bspline
    p0::NTuple{2,Float64}   # :line origin
    d::NTuple{2,Float64}    # :line direction
    degree::Int             # :bspline degree (0 for :line)
    knots::Vector{Float64}  # :bspline expanded knot vector
    controls::Vector{NTuple{2,Float64}}
    weights::Vector{Float64}
end

# A 2-D Cartesian point (exactly two coordinates).
function _step_point2d(entities,arg)
    ent=_deref_entity(entities,arg)
    ent===nothing && return nothing
    ent.kind=="CARTESIAN_POINT" || return nothing
    c=ent.args[end]
    c isa Vector && length(c)==2 || return nothing
    all(v->v isa Real && isfinite(Float64(v)),c) || return nothing
    return (Float64(c[1]),Float64(c[2]))
end

# A 2-D direction (exactly two components).
function _step_dir2d(entities,arg)
    ent=_deref_entity(entities,arg)
    ent===nothing && return nothing
    ent.kind=="DIRECTION" || return nothing
    c=ent.args[end]
    c isa Vector && length(c)==2 || return nothing
    all(v->v isa Real && isfinite(Float64(v)),c) || return nothing
    return (Float64(c[1]),Float64(c[2]))
end

# 2-D control-point list for a parametric-space B-spline.
function _ref_points2d(entities,arg)
    arg isa Vector || return nothing
    pts=NTuple{2,Float64}[]
    for item in arg
        p=_step_point2d(entities,item)
        p===nothing && return nothing
        push!(pts,p)
    end
    return pts
end

# Parse the DEFINITIONAL_REPRESENTATION's 2-D geometry of a PCURVE entity:
# a LINE on 2-D point+vector or a (possibly rational) 2-D B-spline. Returns
# nothing when the representation is missing or carries another curve kind.
function _pcurve2d_curve(entities,pc,caller)
    i=_skip_name(pc.args)
    length(pc.args)>=i+1 || return nothing
    dr=_deref_entity(entities,pc.args[i+1])
    dr===nothing && return nothing
    dr.kind=="DEFINITIONAL_REPRESENTATION" || return nothing
    j=_skip_name(dr.args)
    items=dr.args[j]
    items isa Vector && !isempty(items) || return nothing
    cv=_deref_entity(entities,first(items))
    cv===nothing && return nothing
    if cv.kind=="LINE"
        p=_step_point2d(entities,_first_ref(cv.args))
        refs=[a for a in cv.args if a isa Tuple && a[1]===:ref]
        vec=length(refs)>=2 ? _deref_entity(entities,refs[2]) : nothing
        (p===nothing || vec===nothing || vec.kind!="VECTOR") && return nothing
        i2=_skip_name(vec.args)
        length(vec.args)>=i2+1 || return nothing
        d=_step_dir2d(entities,vec.args[i2])
        d===nothing && return nothing
        # The VECTOR magnitude is the line's parametric speed.
        mag=try Float64(vec.args[end]) catch; return nothing end
        (isfinite(mag) && mag!=0.0) || return nothing
        return _Pcurve2d(:line,p,(d[1]*mag,d[2]*mag),0,Float64[],
                         NTuple{2,Float64}[],Float64[])
    elseif cv.kind=="B_SPLINE_CURVE_WITH_KNOTS"
        k=_skip_name(cv.args)
        length(cv.args)>=k+6 || return nothing
        deg=_as_int(cv.args[k],caller,"pcurve degree")
        pts=_ref_points2d(entities,cv.args[k+1])
        pts===nothing && return nothing
        knots=_expand_knots(cv.args[k+5],cv.args[k+6],caller)
        return _Pcurve2d(:bspline,(0.0,0.0),(0.0,0.0),deg,knots,pts,
                         ones(length(pts)))
    elseif cv.kind=="COMPLEX"
        parts=_complex_parts(cv,caller)
        (haskey(parts,"B_SPLINE_CURVE") &&
         haskey(parts,"B_SPLINE_CURVE_WITH_KNOTS")) || return nothing
        bc=parts["B_SPLINE_CURVE"];bk=parts["B_SPLINE_CURVE_WITH_KNOTS"]
        length(bc.args)>=2 && length(bk.args)>=2 || return nothing
        deg=_as_int(bc.args[1],caller,"pcurve degree")
        pts=_ref_points2d(entities,bc.args[2])
        pts===nothing && return nothing
        knots=_expand_knots(bk.args[1],bk.args[2],caller)
        weights=haskey(parts,"RATIONAL_B_SPLINE_CURVE") &&
            !isempty(parts["RATIONAL_B_SPLINE_CURVE"].args) ?
            _float_vec(parts["RATIONAL_B_SPLINE_CURVE"].args[1],caller) :
            ones(length(pts))
        length(weights)==length(pts) || return nothing
        return _Pcurve2d(:bspline,(0.0,0.0),(0.0,0.0),deg,knots,pts,weights)
    end
    return nothing
end

# The PCURVE authored for `use` on the surface `surf_id`: walks the
# EDGE_CURVE's geometry — a SURFACE_CURVE wrapper carries a list of
# curve-on-surface representations; a PCURVE selects its basis surface. A
# pcurve whose basis surface is this face's returns its parsed 2-D curve;
# anything else (plain 3-D geometry, pcurve on another surface, an
# unparseable representation) yields nothing and the caller falls back to
# chart projection.
function _step_use_pcurve(entities,use,surf_id,caller)
    geom=_deref_entity(entities,use.ec.args[4])
    geom===nothing && return nothing
    geom.kind=="SURFACE_CURVE" || return nothing
    i=_skip_name(geom.args)
    length(geom.args)>=i+1 || return nothing
    reps=geom.args[i+1]
    reps isa Vector || (reps=[reps])
    for r in reps
        pc=_deref_entity(entities,r)
        (pc===nothing || pc.kind!="PCURVE") && continue
        j=_skip_name(pc.args)
        length(pc.args)>=j || continue
        basis=_deref_entity(entities,pc.args[j])
        (basis===nothing || basis.id!=surf_id) && continue
        return _pcurve2d_curve(entities,pc,caller)
    end
    return nothing
end

# Evaluate a 2-D pcurve: a LINE extends origin+t·direction; a B-spline runs
# the same clamped De Boor as `_brep_lift` in homogeneous (x·w, y·w, w).
function _pc2_eval(c::_Pcurve2d,t)
    c.kind==:line && return (c.p0[1]+t*c.d[1],c.p0[2]+t*c.d[2])
    U=c.knots;p=c.degree
    uu=clamp(t,U[p+1],U[end-p])
    wmax=maximum(c.weights)
    Pw=NTuple{4,Float64}[(c.controls[k][1]*c.weights[k]/wmax,
                          c.controls[k][2]*c.weights[k]/wmax,
                          0.0,c.weights[k]/wmax)
                         for k in eachindex(c.controls)]
    h=_brep_deboor4(p,U,Pw,uu)
    h[4]>0.0 || return (h[1],h[2])
    return (h[1]/h[4],h[2]/h[4])
end

# Invert a 2-D point onto a pcurve: a LINE projects in closed form; a
# B-spline runs damped Gauss–Newton on |c(t)-q|² inside the active domain,
# seeded by the previous chain parameter (or, for the first point, by a
# coarse knot-span scan). Returns the converged parameter or nothing when
# no seed converges — the caller turns that into the PCURVE mismatch
# certificate.
function _pc2_param(c::_Pcurve2d,q,seed,caller,ecid,fid)
    if c.kind==:line
        dn=c.d[1]*c.d[1]+c.d[2]*c.d[2]
        dn>0.0 || return nothing
        return ((q[1]-c.p0[1])*c.d[1]+(q[2]-c.p0[2])*c.d[2])/dn
    end
    U=c.knots;p=c.degree
    u0=U[p+1];u1=U[end-p]
    u1-u0>0.0 || return nothing
    seeds=Float64[]
    seed!==nothing && push!(seeds,seed)
    spans=[k for k in p+1:length(U)-p-1 if U[k+1]>U[k]]
    for k in spans
        push!(seeds,(U[k]+U[k+1])/2)
    end
    best=Inf;bt=NaN
    for s in seeds
        t=s
        ok=false
        for _ in 1:60
            hfd=max(1e-7*(u1-u0),1e-12)
            ta=clamp(t-hfd,u0,u1);tb=clamp(t+hfd,u0,u1)
            tb-ta>0 || break
            ca=_pc2_eval(c,ta);cb=_pc2_eval(c,tb)
            ds=(cb[1]-ca[1],cb[2]-ca[2])
            dn=ds[1]*ds[1]+ds[2]*ds[2]
            dn>0.0 || break
            ct=_pc2_eval(c,t)
            r=(ct[1]-q[1],ct[2]-q[2])
            step=(r[1]*ds[1]+r[2]*ds[2])/dn*(tb-ta)
            t2=clamp(t-step,u0,u1)
            if abs(t2-t)<=1e-12*max(u1-u0,1.0)
                c2=_pc2_eval(c,t2)
                res=hypot(c2[1]-q[1],c2[2]-q[2])
                res<best && (best=res;bt=t2;ok=true)
                break
            end
            t=t2
        end
        if !ok
            ct=_pc2_eval(c,t)
            res=hypot(ct[1]-q[1],ct[2]-q[2])
            res<best && (best=res;bt=t)
        end
    end
    isnan(bt) && return nothing
    return bt
end
