# ── single-face NURBS closed shells (B_SPLINE_SURFACE_WITH_KNOTS) ────────────
#
# OCC exports a closed free-form solid as one ADVANCED_FACE over a
# B_SPLINE_SURFACE_WITH_KNOTS whose boundary loop is only seam edges — the
# same layout the torus importer handles: every EDGE_CURVE self-loops on a
# vertex that sits on the surface's parametric boundary (param corners in
# practice), and each seam curve is traversed once in each direction. Pole-
# collapsed sides (sphere-style) appear as degenerate self-loop edges whose
# vertex coincides with the collapsed parametric row.
#
# The certificate does not trust the declared u_closed/v_closed flags: it
# measures the surface itself. Each of the four parametric boundary rows is
# classified by evaluation — :seam when the two opposite rows coincide within
# tolerance at every sampled partner coordinate, :poles when both opposite
# rows independently collapse to a constant point, :open otherwise (genuine
# trim — reject). The declared bound edges are then cross-checked: every
# bound vertex must land on an emitted parametric-boundary vertex, every
# edge used twice or more is a seam traversal pair, and every edge used once
# must self-loop on a pole vertex — so an interior trim (hole circle, slit)
# cannot masquerade as a seam.
#
# Emission welds seams and poles by construction — the u_hi row is never
# minted as separate vertices; cells wrap modulo the grid — so the triangle
# complex is watertight exactly and the directed shared-edge audit, a
# positive signed volume and `validate` certify the result. The emitted grid
# is a straight-line approximation of the curved surface (48x24 cells, same
# density as the torus producer); topology is certified, geometric accuracy
# is a sampling bound, not exact.

const _BREP_NURBS_NU=48
const _BREP_NURBS_NV=24

# Surface entity is a supported B-spline record: either the direct
# B_SPLINE_SURFACE_WITH_KNOTS form or the complex rational form.
function _step_is_nurbs_surface(ent)
    ent===nothing && return false
    ent.kind=="B_SPLINE_SURFACE_WITH_KNOTS" && return true
    if ent.kind=="COMPLEX"
        for part in ent.args
            part isa StepEntity && part.kind=="B_SPLINE_SURFACE" && return true
        end
    end
    return false
end

# Short readable name for a surface entity, expanding complex parts for the
# unsupported-topology diagnostic.
function _step_surface_name(ent)
    ent===nothing && return "(missing surface)"
    ent.kind=="COMPLEX" || return ent.kind
    kinds=String[]
    for part in ent.args
        part isa StepEntity && push!(kinds,part.kind)
    end
    return isempty(kinds) ? "COMPLEX" : "COMPLEX($(join(kinds,"+")))"
end

_nurbs_pt_eq(a,b,tol)=abs(a[1]-b[1])<=tol && abs(a[2]-b[2])<=tol && abs(a[3]-b[3])<=tol

# All-equal within tol test for a sampled row.
function _nurbs_row_degenerate(row,tol)
    a=first(row)
    for p in row
        _nurbs_pt_eq(p,a,tol) || return false
    end
    return true
end

# Classify one pair of opposite parametric rows:
#   :seam   rows coincide pointwise within tol (closed by welding)
#   :poles  each row independently collapses to a constant point
#   :open   genuine boundary — not a closed surface direction
function _nurbs_edge_mode(row_lo,row_hi,tol)
    all(j->_nurbs_pt_eq(row_lo[j],row_hi[j],tol),eachindex(row_lo)) &&
        return :seam
    (_nurbs_row_degenerate(row_lo,tol) && _nurbs_row_degenerate(row_hi,tol)) &&
        return :poles
    return :open
end

# Certify and emit a closed shell over exactly one B-spline surface face.
# Returns a Mesh or nothing when the face is not a certifiably closed
# parametric domain (genuine trims keep the classified reject).
function _step_nurbs_face_shell(entities,face,curvekinds,caller)
    surf_ent=face.surf
    nurbs=try
        _step_nurbs_surface(surf_ent,entities,caller)
    catch err
        err isa InterruptException && rethrow()
        push!(curvekinds,"malformed B-spline surface: $(sprint(showerror,err))")
        return nothing
    end
    nurbs isa NURBSSurface ||
        (push!(curvekinds,"face surface $(_step_surface_name(surf_ent))");
         return nothing)

    nu,nv=_BREP_NURBS_NU,_BREP_NURBS_NV
    pu,pv=nurbs.degree_u,nurbs.degree_v
    U,V=nurbs.knots_u,nurbs.knots_v
    u_lo,u_hi=U[pu+1],U[end-pu]
    v_lo,v_hi=V[pv+1],V[end-pv]

    # Scale-aware tolerance from the control-polygon extent.
    cx0=cx1=cy0=cy1=cz0=cz1=0.0
    for row in eachrow(nurbs.controls), p in row
        cx0=min(cx0,p[1]);cx1=max(cx1,p[1])
        cy0=min(cy0,p[2]);cy1=max(cy1,p[2])
        cz0=min(cz0,p[3]);cz1=max(cz1,p[3])
    end
    scale=max(cx1-cx0,cy1-cy0,cz1-cz0,1.0)
    tol=1e-9*scale

    row_u_lo=[nurbs_eval(nurbs,u_lo,v_lo+(v_hi-v_lo)*k/nv) for k in 0:nv]
    row_u_hi=[nurbs_eval(nurbs,u_hi,v_lo+(v_hi-v_lo)*k/nv) for k in 0:nv]
    row_v_lo=[nurbs_eval(nurbs,u_lo+(u_hi-u_lo)*k/nu,v_lo) for k in 0:nu]
    row_v_hi=[nurbs_eval(nurbs,u_lo+(u_hi-u_lo)*k/nu,v_hi) for k in 0:nu]

    u_mode=_nurbs_edge_mode(row_u_lo,row_u_hi,tol)
    u_mode==:open &&
        (push!(curvekinds,"B-spline surface open in u — trimmed NURBS face");
         return nothing)
    v_mode=_nurbs_edge_mode(row_v_lo,row_v_hi,tol)
    v_mode==:open &&
        (push!(curvekinds,"B-spline surface open in v — trimmed NURBS face");
         return nothing)

    pole_lo_u=u_mode==:poles ? row_u_lo[1] : nothing
    pole_hi_u=u_mode==:poles ? row_u_hi[1] : nothing
    pole_lo_v=v_mode==:poles ? row_v_lo[1] : nothing
    pole_hi_v=v_mode==:poles ? row_v_hi[1] : nothing

    # ── emit grid vertices with weld keys ────────────────────────────────
    # Logical node (i,j) maps onto an emitted vertex: seams wrap modulo the
    # cell count, poles collapse the whole row to a single vertex. u-poles
    # take precedence at corners — a doubly-poled corner is the same point.
    verts=NTuple{3,Float64}[]
    keys=Dict{Any,Int32}()
    function vid(i,j)
        ii=u_mode==:seam ? mod(i,nu) : i
        jj=v_mode==:seam ? mod(j,nv) : j
        key=if u_mode==:poles && ii==0
            (:poleu,0)
        elseif u_mode==:poles && ii==nu
            (:poleu,1)
        elseif v_mode==:poles && jj==0
            (:polev,0)
        elseif v_mode==:poles && jj==nv
            (:polev,1)
        else
            (:g,ii,jj)
        end
        idx=get(keys,key,0)
        idx!=0 && return Int32(idx)
        p=key[1]==:poleu ? (key[2]==0 ? pole_lo_u : pole_hi_u) :
          key[1]==:polev ? (key[2]==0 ? pole_lo_v : pole_hi_v) :
          nurbs_eval(nurbs,u_lo+(u_hi-u_lo)*ii/nu,v_lo+(v_hi-v_lo)*jj/nv)
        push!(verts,p)
        keys[key]=length(verts)
        return Int32(length(verts))
    end

    itris=NTuple{3,Int32}[]
    for i in 0:nu-1, j in 0:nv-1
        quad=(vid(i,j),vid(i+1,j),vid(i+1,j+1),vid(i,j+1))
        out=Int32[];prev=Int32(0)
        for q in quad
            q!=prev && push!(out,q);prev=q
        end
        length(out)>1 && first(out)==last(out) && pop!(out)
        length(out)>=3 || continue
        for k in 2:length(out)-1
            push!(itris,(out[1],out[k],out[k+1]))
        end
    end
    _checked_brep_mul(3,length(itris),caller,"NURBS shell triangle count")
    length(verts)<=_BREP_MAX_SURFACE_NODES || throw(ArgumentError(
        "$caller: NURBS shell exceeds the surface node limit"))
    length(itris)<=_BREP_MAX_SURFACE_TRIS || throw(ArgumentError(
        "$caller: NURBS shell exceeds the surface triangle limit"))

    # Parametric-boundary vertex set (emitted positions): poles plus the
    # perimeter rows — seam sides contribute their lo row only.
    boundary=NTuple{3,Float64}[]
    if u_mode==:poles
        push!(boundary,pole_lo_u,pole_hi_u)
    else
        append!(boundary,@view row_u_lo[1:end-1])
    end
    if v_mode==:poles
        push!(boundary,pole_lo_v,pole_hi_v)
    else
        append!(boundary,@view row_v_lo[1:end-1])
    end
    on_boundary(v)=any(p->_nurbs_pt_eq(p,v,tol),boundary)
    poles=NTuple{3,Float64}[]
    u_mode==:poles && append!(poles,(pole_lo_u,pole_hi_u))
    v_mode==:poles && append!(poles,(pole_lo_v,pole_hi_v))
    on_pole(v)=any(p->_nurbs_pt_eq(p,v,tol),poles)

    # ── declared-bound cross-check ───────────────────────────────────────
    fe=_deref_entity(entities,(:ref,face.id))
    fe===nothing && throw(ArgumentError(
        "$caller: face #$(face.id) record is missing"))
    edge_use_count=Dict{Int,Int}()
    seen_uses=Vector{_BREP_EDGE_USE}[]
    for uses in _step_face_bound_uses(entities,fe,curvekinds,caller)
        push!(seen_uses,uses)
        for u in uses
            edge_use_count[u.ec.id]=get(edge_use_count,u.ec.id,0)+1
            for v in (u.v1,u.v2)
                on_boundary(v) || (push!(curvekinds,
                    "bound vertex is interior — trimmed NURBS face");
                    return nothing)
            end
        end
    end
    for uses in seen_uses, u in uses
        if edge_use_count[u.ec.id]==1
            # A single-use edge must be a pole-degenerate self-loop (sphere
            # layout). Any other once-traversed edge — an interior hole
            # circle or a boundary-reaching slit — is a genuine trim.
            (u.v1==u.v2 && on_pole(u.v1)) ||
                (push!(curvekinds,"single-use non-pole bound edge — trimmed NURBS face");
                 return nothing)
        end
    end

    # ── closed-complex audit + orientation + volume ──────────────────────
    edge_use=Dict{Tuple{Int32,Int32},Int}()
    edge_dir=Dict{Tuple{Int32,Int32},Int}()
    for t in itris,(a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
        key=a<b ? (a,b) : (b,a)
        edge_use[key]=get(edge_use,key,0)+1
        edge_dir[key]=get(edge_dir,key,0)+(a==key[1] ? 1 : -1)
    end
    for (key,n) in edge_use
        n==2 || throw(ArgumentError(
            "$caller: NURBS shell is not closed — edge $key has $n incident triangles"))
        edge_dir[key]==0 || throw(ArgumentError(
            "$caller: NURBS shell is inconsistently oriented across edge $key"))
    end
    vscale=max(_poly_scale(verts),1.0)
    V=sum(_brep_vdot(verts[t[1]],
        _brep_vcross(verts[t[2]],verts[t[3]]))/6 for t in itris)
    abs(V)>1e-12*vscale*vscale*vscale || throw(ArgumentError(
        "$caller: NURBS shell bounds zero volume"))
    if V<0
        itris=map(t->(t[1],t[3],t[2]),itris)
    end
    C=Matrix{Float64}(undef,3,length(verts))
    for (i,p) in enumerate(verts); C[:,i]=[p...]; end
    T=Matrix{Int32}(undef,3,length(itris))
    for (i,t) in enumerate(itris); T[:,i]=[t...]; end
    surface=Mesh(C;tris=T)
    diag=validate(surface)
    diag.ok || throw(ErrorException(
        "$caller: NURBS shell boundary mesh is invalid — "*
        join(diag.messages,"; ")))
    return surface
end
