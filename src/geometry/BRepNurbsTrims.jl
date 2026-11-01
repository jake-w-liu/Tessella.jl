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
    return NTuple{3,Float64}[nurbs_eval(c,u0+(u1-u0)*k/m) for k in 0:m-1]
end

# Traversal-ordered canonical 3-D chain of one classified use (for wire rings,
# uv projection and segment incidence). Open chains run junction→junction; a
# closed-spline chain closes back onto its vertex — direction follows the
# ORIENTED_EDGE flag folded with the EDGE_CURVE same_sense flag, captured at
# classification into aux.pts.
function _nurbs_trim_chain(ur,ptol)
    if ur.kind==:line
        return NTuple{3,Float64}[ur.u.v1,ur.u.v2]
    elseif ur.kind==:closedspline
        return ur.aux.pts
    end
    ch=ur.aux.pts
    _posnear(ur.u.v1,first(ch),ptol) && return collect(ch)
    return reverse(collect(ch))
end

# Certify and emit a closed shell that mixes trimmed B-spline faces with
# planar faces. Returns a Mesh or nothing when the inventory is not this
# certificate's family; throws on certified mismatches.
function _step_nurbs_trimmed_shell(entities,faces,curvekinds,caller)
    length(faces)>=2 || return nothing
    all(f->f.role==1,faces) || return nothing
    all(f->f.surf!==nothing,faces) || return nothing
    nurbs_i=Int[]
    for (i,f) in enumerate(faces)
        if _step_is_nurbs_surface(f.surf)
            push!(nurbs_i,i)
        elseif f.surf.kind!="PLANE"
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
                    u.v1!=u.v2 || (push!(curvekinds,
                        "degenerate bound edge — pole NURBS face");
                        return nothing)
                    push!(out,(kind=:line,u=u,aux=nothing))
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
    plane_i=setdiff(eachindex(faces),nurbs_i)
    for i in plane_i
        p0,pn=_step_plane_frame(entities,faces[i].surf,caller)
        for bound in fuse[i],ur in bound
            ch=_nurbs_trim_chain(ur,ptol)
            for v in ch
                abs(_brep_vdot(_brep_vsub(v,p0),pn))<=ptol || throw(
                    ArgumentError("$caller: sampled edge #$(ur.u.ec.id) "*
                        "leaves PLANE #$(faces[i].surf.id)"))
            end
        end
    end

    # ── NURBS faces: uv domain polygons + on-surface lift certificate ─────
    # nuv[i] = per-bound unwrapped uv vertex lists; ncan[i] = parallel
    # canonical 3-D positions.
    nuv=Dict{Int,Vector{Vector{NTuple{2,Float64}}}}()
    ncan=Dict{Int,Vector{Vector{NTuple{3,Float64}}}}()
    patch=Dict{Int,_BrepNurbsPatch}()
    periods=Dict{Int,Tuple{Float64,Float64}}()
    for i in nurbs_i
        f=faces[i]
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
        patch[i]=pk
        stol=1e-9*max(_poly_scale([nurbs.controls...]),1.0)
        uper,vper=_brep_nurbs_seam(pk,stol)
        periods[i]=(uper,vper)
        bounds_uv=Vector{NTuple{2,Float64}}[]
        bounds_p3=Vector{NTuple{3,Float64}}[]
        uvtol=1e-7*max(pk.u_hi-pk.u_lo,pk.v_hi-pk.v_lo)
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
                for k in (isempty(uvs) ? 1 : 2):length(ch)
                    p=ch[k]
                    uv=_brep_nurbs_project(pk,p,prev===nothing ?
                        nothing : (prev[1],prev[2]))
                    uv===nothing && throw(ArgumentError(
                        "$caller: edge #$(ur.u.ec.id) of face #$(f.id) "*
                        "does not project onto its B-spline surface"))
                    u0,v0=uv
                    # On-surface certificate on the in-domain projection
                    # first — the unwrapped sheet coordinate below is only
                    # for the uv polygon (clamped lift would alias copies).
                    q=_brep_lift(pk,u0,v0)
                    _nurbs_pt_eq(q,p,stol) || throw(ArgumentError(
                        "$caller: edge #$(ur.u.ec.id) leaves the B-spline "*
                        "surface of face #$(f.id)"))
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
                "$caller: NURBS face #$(f.id) wire has fewer than three "*
                "projected vertices"))
            if abs(uvs[n][1]-uvs[1][1])<=uvtol && abs(uvs[n][2]-uvs[1][2])<=uvtol
                pop!(uvs);pop!(p3s);n-=1
            elseif uper>0.0
                # The last vertex may sit a full period from the first —
                # identical uv position, distinct sheet.
                du=abs(uvs[n][1]-uvs[1][1]);dv=abs(uvs[n][2]-uvs[1][2])
                (abs(du-round(du/uper)*uper)<=uvtol && dv<=uvtol) &&
                    (pop!(uvs);pop!(p3s);n-=1)
            end
            n>=3 || throw(ArgumentError(
                "$caller: NURBS face #$(f.id) wire has fewer than three "*
                "projected vertices"))
            push!(bounds_uv,uvs);push!(bounds_p3,p3s)
        end
        nuv[i]=bounds_uv;ncan[i]=bounds_p3
    end

    # ── uv polygon audit: simple closed wires, one unique outer ──────────
    for i in nurbs_i
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
        if length(bound)==1 && bound[1].kind==:closedspline
            append!(ring,bound[1].aux.cyc)
            return ring
        end
        for ur in bound
            ur.kind==:closedspline && throw(ArgumentError(
                "$caller: closed B-spline edge #$(ur.u.ec.id) shares a "*
                "wire with other edges"))
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

    for i in nurbs_i
        pk=patch[i]
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
            "$caller: NURBS face #$(faces[i].id) has an empty uv boundary"))
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
                _,du_c,dv_c=_brep_nurbs_d1(pk,umin+fu*uspan,vmin+fv*vspan)
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
                    emituv((u,v),_brep_lift(pk,u,v))
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
                p3===nothing && (p3=_brep_lift(pk,u,v))
                push!(tri,p3)
            end
            # Lifted weld-collapse: a cell whose two vertices lift to the
            # same 3-D node (pole fan or a folded domain corner) drops out.
            tri[1]==tri[2] || tri[2]==tri[3] || tri[3]==tri[1] ||
                push!(cells,(tri[1],tri[2],tri[3]))
        end
        isempty(cells) && throw(ArgumentError(
            "$caller: NURBS face #$(faces[i].id) triangulates to nothing"))
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
