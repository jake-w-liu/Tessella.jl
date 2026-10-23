# `occ.importShapes` / `OCC_Internals::importShapes` — ISO-10303-21 STEP
# solids imported as native OCC records (points/lines/wires/plane
# faces/shells/volumes) rather than the triangulated `import_step`
# meshes. The BRep module supplies the parser and the face/bound/edge
# graph walk; this file converts each closed shell's rings into shared
# vertex/edge entities — identical vertices and edges publish once, so
# the shell validator's two-surface-per-curve rule holds by
# construction.
#
# Two face paths share one per-solid dedupe state: fully straight-edged
# faces on `PLANE` or untyped `FACE` surfaces keep the vertex-ring
# assembly (identical output to earlier revisions); every other face
# takes the structural walk below — ORIENTED_EDGE uses on shared
# EDGE_CURVE entities materialize as native OCC records (line, circle,
# ellipse, degenerate, NURBS `copied_curve`), STEP surfaces become
# plane/cylinder/cone/sphere/torus/tensor_bspline records, and each
# bound wires through `_occ_wire_from_uses!` so the loop-validity
# certificate's <3-edge rejection still applies to degenerate straight
# loops while a legal two-edge curved loop imports intact. Pcurves are
# projected onto the stored face record and shifted into the wire's
# unwrapped uv window; STEP `SEAM_CURVE`/`SURFACE_CURVE` wrappers resolve
# to their innermost geometry.

# Solid-level dedupe state: vertex positions → point tags, unordered
# vertex-tag pairs → (line id, stored stop tag).
struct _OCCImportDedupe
    ppos::Vector{NTuple{3,Float64}}
    ptag::Vector{Int}
    edge::Dict{NTuple{2,Int},Tuple{UInt64,Int}}
end
_OCCImportDedupe()=_OCCImportDedupe(NTuple{3,Float64}[],Int[],
                                   Dict{NTuple{2,Int},Tuple{UInt64,Int}}())

function _occ_import_point!(m,d:: _OCCImportDedupe,p)
    for i in eachindex(d.ppos)
        _model_occ_vertices_coincide(d.ppos[i],p) && return d.ptag[i]
    end
    t=_model_occ_add_point!(m,p[1],p[2],p[3])
    push!(d.ppos,p); push!(d.ptag,t)
    return t
end

# One ring (ordered positions, no duplicated closing point) → wire tag.
# `flip` winds the loop backwards (hole wires are consumed CW relative
# to the face normal by `add_plane_surface` convention checks — STEP
# bounds already carry the correct orientation so `flip` stays false).
function _occ_import_wire!(m,d,kernel,ring,flip,caller)
    verts=[_occ_import_point!(m,d,p) for p in ring]
    uses=Tuple{UInt64,Bool}[]
    K=length(verts)
    for j in 1:K
        a,b=verts[j],verts[mod1(j+1,K)]
        key=a<b ? (a,b) : (b,a)
        eu=get(d.edge,key,nothing)
        if eu===nothing
            et=_model_occ_add_line!(m,a,b)
            eu=(kernel.bindings[(1,et)],b)
            d.edge[key]=eu
        end
        rev=eu[2]!=b
        push!(uses,flip ? (eu[1],!rev) : (eu[1],rev))
    end
    flip && (uses=reverse(uses))
    return _occ_wire_from_uses!(m,uses,caller)
end

# ── structural STEP→OCC conversion ──────────────────────────────────────────
#
# Per-solid state for the structural path. STEP topology identities map
# straight onto OCC entities — `vert_map`/`edge_map` hold the shared
# shape ids the shell validator's two-use rule needs; `d.edge` is the
# ring path's line-edge dedupe keyed by vertex tags so a straight edge
# shared between a ring-routed face and a structural face still binds
# once.
struct _OCCImportStructural
    vert_map::Dict{Int,UInt64}
    edge_map::Dict{Int,Tuple{UInt64,Bool}}
end
_OCCImportStructural()=_OCCImportStructural(Dict{Int,UInt64}(),
                                            Dict{Int,Tuple{UInt64,Bool}}())

# STEP VERTEX_POINT → OCC vertex shape id. Dedupe by entity id first,
# then by position against the ring path's table so unwelded-but-
# coincident files still share topology.
function _occ_import_vertex_id!(m,d,st,entities,ref,caller)
    ent=BRep._deref_entity(entities,ref)
    (ent!==nothing &&
     (ent.kind=="VERTEX_POINT" || ent.kind=="CARTESIAN_POINT")) ||
        throw(ArgumentError(
            "$caller: $(ref) does not resolve to a STEP vertex"))
    hit=get(st.vert_map,ent.id,nothing)
    hit!==nothing && return hit
    p=BRep._step_vertex_point(entities,ref)
    p===nothing && throw(ArgumentError(
        "$caller: VERTEX_POINT #$(ent.id) has unresolvable geometry"))
    tag=_occ_import_point!(m,d,p)
    id=_occ_import_bound_vertex!(m,tag,caller)
    st.vert_map[ent.id]=id
    return id
end

# Tag → bound OCC vertex id (points may already exist in the kernel).
function _occ_import_bound_vertex!(m,tag,caller)
    id=get(m.occ_kernel.bindings,(0,tag),nothing)
    id===nothing && throw(ArgumentError(
        "$caller: internal — OCC Point[$tag] missing after add_point"))
    return id
end

# AXIS2_PLACEMENT_3D → (location, axis, X, Y) — `$` axis/ref fall back
# to the frame defaults the STEP spec gives an incomplete placement.
function _occ_import_placement(entities,ent,caller)
    refs=[a for a in ent.args if a isa Tuple && a[1]===:ref]
    loc=length(refs)>=1 ? BRep._deref_point(entities,refs[1]) : nothing
    loc===nothing && (loc=(0.0,0.0,0.0))
    z=length(refs)>=2 ? BRep._deref_direction(entities,refs[2]) : nothing
    z===nothing && (z=(0.0,0.0,1.0))
    n=BRep._brep_vunit(z,caller)
    x=length(refs)>=3 ? BRep._deref_direction(entities,refs[3]) : nothing
    X=x===nothing ? BRep._brep_default_x(n,caller) :
        BRep._brep_ortho(x,n,caller)
    Y=BRep._brep_vcross(n,X)
    return loc,n,X,Y
end

# Reverse a NURBSCurve so its forward param runs against the stored
# edge direction (STEP `same_sense` `.F.`): knot mirror + pole/weight
# reversal, same transform `_occ_heal_pcurve_reverse` uses for spl2d.
function _occ_import_reverse_nurbs(curve)
    T=curve.knots[1]+curve.knots[end]
    return NURBS.NURBSCurve(curve.degree,
        T .- reverse(curve.knots),reverse(curve.controls),
        reverse(curve.weights))
end

# Azimuth of `p` about the conic frame — used to seed circle/ellipse
# edge ranges from their endpoint vertices. `rx`/`ry` scale the frame
# for the ellipse's eccentric angle (circles pass 1.0).
function _occ_import_theta(g,p;rx=1.0,ry=1.0)
    d=BRep._brep_vsub(p,g.center)
    return _gm_atan2(BRep._brep_vdot(d,g.Y)/ry,
                     BRep._brep_vdot(d,g.X)/rx)
end

# EDGE_CURVE → (edge shape id, align). `st.edge_map` dedupes shared
# STEP entities; `align` reports whether the bound edge's stored
# direction agrees with STEP's `vertex_start → vertex_end` — only
# edges materialized fresh here carry STEP-ordered children, so a line
# edge reused from the ring path may need its use flipped. `same_sense`
# `.F.` reverses the parametrization so `eval(t0)` lands on v1.
function _occ_import_step_edge!(m,d,st,entities,eref,curvekinds,caller)
    kernel=m.occ_kernel
    ec=BRep._deref_entity(entities,eref)
    (ec!==nothing && ec.kind=="EDGE_CURVE") || throw(ArgumentError(
        "$caller: $(eref) does not resolve to an EDGE_CURVE"))
    hit=get(st.edge_map,ec.id,nothing)
    hit!==nothing && return hit
    length(ec.args)>=4 || throw(ArgumentError(
        "$caller: EDGE_CURVE #$(ec.id) is truncated"))
    pa=BRep._step_vertex_point(entities,ec.args[2])
    pb=BRep._step_vertex_point(entities,ec.args[3])
    (pa===nothing || pb===nothing) && throw(ArgumentError(
        "$caller: EDGE_CURVE #$(ec.id) has unresolvable vertex geometry"))
    sense=length(ec.args)>=5 ? BRep._brep_enum(ec.args[5]) : nothing
    sense===nothing && (sense=true)
    va=_occ_import_vertex_id!(m,d,st,entities,ec.args[2],caller)
    vb=_occ_import_vertex_id!(m,d,st,entities,ec.args[3],caller)
    closed=_model_occ_vertices_coincide(pa,pb)
    curve=BRep._step_curve_base(entities,
        BRep._deref_entity(entities,ec.args[4]),0)
    curve===nothing && throw(ArgumentError(
        "$caller: EDGE_CURVE #$(ec.id) has unresolvable geometry"))
    geometry=kind=nothing
    if closed && curve.kind=="LINE"
        # A LINE cannot be a closed edge — this is a degenerate (pole)
        # edge; the uv connector lands at face-bind time.
        geometry=(occ=:degenerate,t0=0.0,t1=_OCC_TWO_PI)
        kind=:degenerate
    elseif curve.kind=="LINE"
        atag=_occ_import_point!(m,d,pa); btag=_occ_import_point!(m,d,pb)
        key=atag<btag ? (atag,btag) : (btag,atag)
        prev=get(d.edge,key,nothing)
        if prev!==nothing
            # A straight edge shared with a ring-routed face — reuse it
            # so the shell validator sees one edge, two faces. The ring
            # path's stored stop vertex decides `align`.
            entry=(prev[1],prev[2]==btag)
            st.edge_map[ec.id]=entry
            return entry
        end
        delta=BRep._brep_vsub(pb,pa); len=BRep._brep_vdot(delta,delta)
        len>0 || throw(ArgumentError(
            "$caller: EDGE_CURVE #$(ec.id) has coincident endpoints"))
        len=sqrt(len)
        geometry=(occ=:line,origin=pa,
                  dir=BRep._brep_vunit(delta,caller),
                  t0=0.0,t1=len)
        kind=:line
    elseif curve.kind=="CIRCLE"
        C,n,X,Y,r=_occ_import_conic(entities,curve,caller)
        if !sense
            n=(-n[1],-n[2],-n[3]); Y=(-Y[1],-Y[2],-Y[3])
        end
        g=(center=C,X=X,Y=Y)
        ta=_occ_import_theta(g,pa); tb=_occ_import_theta(g,pb)
        if closed
            tb=ta+_OCC_TWO_PI
        else
            dth=mod(tb-ta,_OCC_TWO_PI)
            dth<1e-9 && (dth=_OCC_TWO_PI)
            tb=ta+dth
        end
        geometry=(occ=:circle,center=C,X=X,Y=Y,n=n,r=r,t0=ta,t1=tb)
        kind=:circle
    elseif curve.kind=="ELLIPSE"
        C,n,X,Y,r1r2=_occ_import_conic(entities,curve,caller;r2=true)
        rx,ry=r1r2
        if !sense
            n=(-n[1],-n[2],-n[3]); Y=(-Y[1],-Y[2],-Y[3])
        end
        g=(center=C,X=X,Y=Y)
        ta=_occ_import_theta(g,pa;rx=rx,ry=ry)
        tb=_occ_import_theta(g,pb;rx=rx,ry=ry)
        if closed
            tb=ta+_OCC_TWO_PI
        else
            dth=mod(tb-ta,_OCC_TWO_PI)
            dth<1e-9 && (dth=_OCC_TWO_PI)
            tb=ta+dth
        end
        geometry=(occ=:ellipse,center=C,X=X,Y=Y,n=n,rx=rx,ry=ry,
                  t0=ta,t1=tb)
        kind=:ellipse
    else
        nc=BRep._step_nurbs_curve(curve,entities)
        nc===nothing && throw(ArgumentError(
            "$caller: EDGE_CURVE #$(ec.id) uses $(curve.kind) — only " *
            "LINE/CIRCLE/ELLIPSE/B_SPLINE curves import currently"))
        sense || (nc=_occ_import_reverse_nurbs(nc))
        cv=_ModelOCCNURBSCurve(:bspline,nc)
        lo,hi=nc.knots[1],nc.knots[end]
        if closed
            t0,t1=lo,hi
        else
            t0=_model_occ_curve_project(cv,pa,lo,hi)[1]
            t1=_model_occ_curve_project(cv,pb,lo,hi)[1]
            t1>t0 || throw(ArgumentError(
                "$caller: EDGE_CURVE #$(ec.id) vertex parameters do not " *
                "order along its B-spline"))
        end
        geometry=(occ=:copied_curve,curve=cv,t0=t0,t1=t1)
        kind=:copied_curve
    end
    allocated=_model_occ_next_tag(kernel,1,-1,caller)
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,1,allocated,id,
        _ModelOCCShape(1,(kind=kind,geometry=geometry),UInt64[va,vb]))
    kernel.next_shape=id+UInt64(1)
    st.edge_map[ec.id]=(id,true)
    if kind===:line
        atag=_model_occ_shape_tag(kernel,va,caller)
        btag=_model_occ_shape_tag(kernel,vb,caller)
        key=atag<btag ? (atag,btag) : (btag,atag)
        d.edge[key]=(id,btag)
    end
    return (id,true)
end

# CIRCLE/ELLIPSE record → (center, axis, X, Y, radius(es)). `r2=true`
# reads the ellipse's two radii.
function _occ_import_conic(entities,curve,caller;r2=false)
    refs=[a for a in curve.args if a isa Tuple && a[1]===:ref]
    plc=isempty(refs) ? nothing : BRep._deref_entity(entities,refs[1])
    (plc!==nothing && plc.kind=="AXIS2_PLACEMENT_3D") || throw(
        ArgumentError(
            "$caller: $(curve.kind) #$(curve.id) lacks an " *
            "AXIS2_PLACEMENT_3D"))
    C,n,X,Y=_occ_import_placement(entities,plc,caller)
    nums=[Float64(a) for a in curve.args if a isa Real]
    r2 ? length(nums)>=2 : length(nums)>=1 || throw(ArgumentError(
        "$caller: $(curve.kind) #$(curve.id) is missing its radii"))
    r2 && return (C,n,X,Y,(nums[1],nums[2]))
    return C,n,X,Y,nums[1]
end

# STEP EDGE_LOOP → `uses` in traversal order. `flip` inverts the whole
# bound (`.F.` bound flag — traversed backwards with each use flipped),
# matching `_step_loop_ring`'s `reverse!(ring)` on the vertex path.
function _occ_import_loop_uses(m,d,st,entities,loop,sense,curvekinds,
                               caller)
    kernel=m.occ_kernel
    loop.kind=="EDGE_LOOP" || throw(ArgumentError(
        "$caller: $(loop.kind) #$(loop.id) — only EDGE_LOOP bounds " *
        "import on curved faces"))
    i=BRep._skip_name(loop.args)
    length(loop.args)>=i || throw(ArgumentError(
        "$caller: EDGE_LOOP #$(loop.id) is truncated"))
    refs=loop.args[i]
    refs isa Vector || throw(ArgumentError(
        "$caller: EDGE_LOOP #$(loop.id) member list is malformed"))
    uses=Tuple{UInt64,Bool}[]
    allstraight=true
    for r in refs
        oe=BRep._deref_entity(entities,r)
        (oe!==nothing && oe.kind=="ORIENTED_EDGE") || throw(ArgumentError(
            "$caller: EDGE_LOOP #$(loop.id) member is not an ORIENTED_EDGE"))
        ecref=nothing; eflag=true
        for a in oe.args
            a isa Tuple && a[1]===:ref && (ecref=a)
            flag=BRep._brep_enum(a); flag===nothing || (eflag=flag)
        end
        ecref===nothing && throw(ArgumentError(
            "$caller: ORIENTED_EDGE #$(oe.id) has no edge reference"))
        eid,align=_occ_import_step_edge!(m,d,st,entities,ecref,
                                         curvekinds,caller)
        ec=BRep._deref_entity(entities,ecref)
        cv=BRep._deref_entity(entities,ec.args[4])
        allstraight &= BRep._step_straight_edge(entities,cv,0,curvekinds)
        push!(uses,(eid,eflag!=align))
    end
    sense || (uses=reverse([(eid,!rev) for (eid,rev) in uses]))
    # The vertex-ring certificate, kept verbatim for straight loops: a
    # boundary of fewer than three straight, open edges is invalid. A
    # curved loop of one or two edges is legal STEP topology — the
    # loop-closure check below still applies.
    if length(uses)<3
        allself=all(uses) do (eid,_)
            a,b=kernel.shapes[eid].children; a==b
        end
        (allstraight && !allself) && throw(ArgumentError(
            "$caller: face boundary loop has fewer than three edges"))
    end
    return uses
end

# Endpoint-continuity check on one wired `uses` list: each traversal
# must end where the next begins, within the confusion tolerance.
function _occ_import_check_chain(kernel,uses,caller)
    pos(id)=kernel.shapes[id].geometry.point
    n=length(uses)
    n>0 || return
    for i in 1:n
        eid,rev=uses[i]; a,b=kernel.shapes[eid].children
        stop=rev ? a : b
        eid2,rev2=uses[mod1(i+1,n)]; a2,b2=kernel.shapes[eid2].children
        start=rev2 ? b2 : a2
        _occ_modulus(_occ3_sub(pos(stop),pos(start)))<=_OCC_CONFUSION*10 ||
            throw(ArgumentError(
                "$caller: face boundary edges do not form an ordered " *
                "closed chain"))
    end
    return nothing
end

# Shift a pcurve record by (du,dv) — the seam-side copies of an edge
# used twice on a periodic face are the same trace translated by one
# period.
function _occ_import_pc_shift(pc,du::Float64,dv::Float64)
    (du==0.0 && dv==0.0) && return pc
    if hasproperty(pc,:neg2d)
        return (neg2d=(_occ_import_pc_shift(pc.neg2d[1],du,dv),),)
    elseif hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=((o[1]+du,o[2]+dv),d),)
    elseif hasproperty(pc,:circ2d)
        o,x,r=pc.circ2d
        return (circ2d=((o[1]+du,o[2]+dv),x,r),)
    elseif hasproperty(pc,:ellips2d)
        o,x,a,b=pc.ellips2d
        return (ellips2d=((o[1]+du,o[2]+dv),x,a,b),)
    elseif hasproperty(pc,:spl2d)
        knots,p,cps=pc.spl2d
        return (spl2d=(knots,p,[(c[1]+du,c[2]+dv) for c in cps]),)
    end
    throw(ArgumentError(
        "_occ_import_pc_shift: unsupported pcurve kind"))
end

# A pcurve reversed in the edge parameter: `pc'(t)=pc(t0+t1−t)`. A closed
# edge on a periodic face wires down the same uv strip the wire climbed —
# its copy must run backwards to close the polygon (the rim circles of a
# slit lateral face travel opposite senses around the window). Lines and
# splines reverse natively; conics wrap in `neg2d` — a right-handed 2-D
# circle cannot write the reversed trace, so the record carries the
# period-rolled forward copy and `_occ_pcurve_d0` evaluates `inner(−t)`.
function _occ_import_pc_reverse(pc,t0::Float64,t1::Float64,caller)
    T=t0+t1
    if hasproperty(pc,:neg2d)
        # `pc(t)=inner(−t)` reverses to `pc'(t)=inner(t−T)` — a pure
        # parameter shift of the inner record, no wrapper needed.
        inner=pc.neg2d[1]
        return _occ_import_pc_paramshift(inner,-T,caller)
    end
    if hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d)
        o,x,rest=if hasproperty(pc,:circ2d)
            pc.circ2d[1],pc.circ2d[2],pc.circ2d[3:end]
        else
            pc.ellips2d[1],pc.ellips2d[2],pc.ellips2d[3:end]
        end
        y=(-x[2],x[1])
        co,si=cos(T),sin(T)
        x2=(co*x[1]+si*y[1],co*x[2]+si*y[2])
        inner=hasproperty(pc,:circ2d) ?
            (circ2d=(o,x2,rest[1])) : (ellips2d=(o,x2,rest[1],rest[2]))
        return (neg2d=(inner,),)
    end
    return _occ_heal_pcurve_reverse(pc,T,nothing,caller)
end

# Advance a pcurve's parameter: `pc'(t)=pc(t+s)`. Circles/ellipses roll
# their phase frame; splines shift their knots back so `pc'` evaluates
# over the same `[t0,t1]` domain.
function _occ_import_pc_paramshift(pc,s::Float64,caller)
    if hasproperty(pc,:neg2d)
        inner=_occ_import_pc_paramshift(pc.neg2d[1],-s,caller)
        return (neg2d=(inner,),)
    elseif hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=((o[1]+s*d[1],o[2]+s*d[2]),d),)
    elseif hasproperty(pc,:circ2d) || hasproperty(pc,:ellips2d)
        o,x,rest=if hasproperty(pc,:circ2d)
            pc.circ2d[1],pc.circ2d[2],pc.circ2d[3:end]
        else
            pc.ellips2d[1],pc.ellips2d[2],pc.ellips2d[3:end]
        end
        y=(-x[2],x[1])
        co,si=cos(s),sin(s)
        x2=(co*x[1]+si*y[1],co*x[2]+si*y[2])
        return hasproperty(pc,:circ2d) ?
            (circ2d=(o,x2,rest[1])) : (ellips2d=(o,x2,rest[1],rest[2]))
    elseif hasproperty(pc,:spl2d)
        knots,p,cps=pc.spl2d
        return (spl2d=(knots.-s,p,cps),)
    end
    throw(ArgumentError(
        "_occ_import_pc_paramshift: unsupported pcurve kind"))
end

# ---- hint-chained NURBS projection --------------------------------------
#
# `_model_nurbs_projection_cached` returns every converged normal-equation
# root. On a periodic or symmetric tensor patch several roots are
# geometrically equivalent (same 3D point at uvs a symmetry step apart)
# and the global pick lands in an arbitrary copy — exactly what broke
# wire continuity on t20's closed B-spline faces. The importer instead
# picks, among the near-best-3D candidates, the root whose uv sits
# nearest a continuity hint (the previous wire sample's uv).
function _occ_import_project_hinted(face_g,P,hint)
    result=_model_nurbs_projection_cached(face_g.patch,P)
    states=result.states
    isempty(states) && return nothing
    dmin=minimum(_model_nurbs_root_distance(s[2],P) for s in states)
    tol=dmin+max(1e-12,dmin*1e-6)
    best=nothing; bd=Inf
    for s in states
        _model_nurbs_root_distance(s[2],P)<=tol || continue
        d=hint===nothing ? 0.0 : _model_nurbs_root_distance(s[1],hint)
        (best===nothing || d<bd) && (best=s; bd=d)
    end
    best===nothing && return nothing
    u,v=best[1][1],best[1][2]
    if hasproperty(face_g,:reparam) && face_g.reparam!==nothing
        ru,rv=face_g.reparam.u,face_g.reparam.v
        u=ru===nothing ? u : _occ_conic_reparam_inv(ru,u)
        v=rv===nothing ? v : _occ_conic_reparam_inv(rv,v)
    end
    return u,v
end

# `_trim_project_pcurve`'s NURBS twin: the projection order follows the
# wire traversal (param t0→t1 forward, t1→t0 reversed) and each sample
# hints the next, so the spl2d stays on one uv copy even when the
# surface's geometric period makes copies indistinguishable.
function _occ_import_project_pcurve(m,kernel,face_g,eid,hint,rev,caller)
    g=kernel.shapes[eid].geometry.geometry
    etag=_model_occ_shape_tag(kernel,eid,caller)
    t0,t1=g.t0,g.t1
    ns=17
    uvs=Vector{NTuple{2,Float64}}(undef,ns)
    prev=hint
    for k in 0:ns-1
        i=rev ? ns-k : k+1
        t=rev ? t1-(t1-t0)*k/(ns-1) : t0+(t1-t0)*k/(ns-1)
        P=_pipe_edge_point(m,g,etag,t)
        uv=_occ_import_project_hinted(face_g,P,prev)
        uv===nothing && throw(ArgumentError(
            "$caller: Curve[$etag] does not project onto the surface"))
        prev=uv; uvs[i]=uv
    end
    a0,aN=uvs[1],uvs[end]
    du=(aN[1]-a0[1])/(t1-t0),(aN[2]-a0[2])/(t1-t0)
    scale=max(1.0,maximum(abs(u[1]-a0[1]) for u in uvs),
              maximum(abs(u[2]-a0[2]) for u in uvs))
    lin=true
    for i in 2:ns-1
        t=t0+(t1-t0)*(i-1)/(ns-1)
        pu=a0[1]+du[1]*(t-t0); pv=a0[2]+du[2]*(t-t0)
        (abs(uvs[i][1]-pu)>1e-9*scale || abs(uvs[i][2]-pv)>1e-9*scale) &&
            (lin=false; break)
    end
    # `lin2d` evaluates `o+d·t` at the ABSOLUTE edge parameter — store
    # the uv at t=0, not at t0 (edges like #282 start at t0<0).
    lin && return (lin2d=((a0[1]-du[1]*t0,a0[2]-du[2]*t0),du),)
    return (spl2d=_pchip_bspline_2d(
        Float64[t0+(t1-t0)*i/(ns-1) for i in 0:ns-1],uvs),)
end

# The uv trace of a degenerate (pole/slit) edge in edge-param order:
# the straight connector between its predecessor's end and successor's
# start uv positions (spanning one period across a periodic cut). `rev`
# flips the traversal direction.
function _occ_import_degenerate_pcurve(start_uv,end_uv,rev)
    delta=(end_uv[1]-start_uv[1],end_uv[2]-start_uv[2])
    L=sqrt(delta[1]*delta[1]+delta[2]*delta[2])
    if L<1e-12
        L=1e-12; dir=(1.0,0.0)
    else
        dir=(delta[1]/L,delta[2]/L)
    end
    o=rev ? end_uv : start_uv
    d=rev ? (-dir[1],-dir[2]) : dir
    return (lin2d=(o,d),),L
end

# The provisional face geometry for one STEP surface entity, evaluated
# per face — `height`/`r1`-style extent fields are fitted to this face's
# own boundary vertices. `ring` is the outer wire's ordered traversal
# positions (untyped `FACE` records fit their frame to it). Returns the
# record plus the surface kind tag.
function _occ_import_surface_base(entities,surf,sense,verts,ring,caller)
    if surf===nothing
        # Untyped `FACE` — the traversal's own Newell normal is the
        # stored frame (the structural path only sees untyped faces
        # inside curved solids; a straight-edged one routes to rings).
        n=BRep._brep_vunit(BRep._newell_normal(ring),caller)
        X=BRep._brep_default_x(n,caller); Y=BRep._brep_vcross(n,X)
        return ((occ=:plane,center=ring[1],axis=n,X=X,Y=Y,
                 reversed=false),:plane)
    end
    if surf.kind=="PLANE"
        refs=[a for a in surf.args if a isa Tuple && a[1]===:ref]
        plc=isempty(refs) ? nothing : BRep._deref_entity(entities,refs[1])
        (plc!==nothing && plc.kind=="AXIS2_PLACEMENT_3D") || throw(
            ArgumentError(
                "$caller: PLANE #$(surf.id) lacks an AXIS2_PLACEMENT_3D"))
        loc,n,X,Y=_occ_import_placement(entities,plc,caller)
        return ((occ=:plane,center=loc,axis=n,X=X,Y=Y,
                 reversed=!sense),:plane)
    elseif surf.kind=="CYLINDRICAL_SURFACE"
        loc,n,X,Y,R=_occ_import_conic(entities,surf,caller)
        hs=[BRep._brep_vdot(BRep._brep_vsub(v,loc),n) for v in verts]
        v0,v1=extrema(hs)
        return ((occ=:cylinder,center=_occ3_add(loc,_occ3_scale(n,v0)),
                 axis=n,X=X,Y=Y,radius=R,height=v1-v0,
                 reversed=!sense),:cylinder)
    elseif surf.kind=="CONICAL_SURFACE"
        loc,n,X,Y,(R,alpha)=_occ_import_conic(entities,surf,caller;r2=true)
        tana=tan(alpha)
        abs(tana)<1e-12 && throw(ArgumentError(
            "$caller: CONICAL_SURFACE #$(surf.id) has ~zero semi-angle"))
        hs=[BRep._brep_vdot(BRep._brep_vsub(v,loc),n) for v in verts]
        h0,h1=extrema(hs)
        r1=R+h0*tana; r2=R+h1*tana
        (r1>=-_OCC_CONFUSION && r2>=-_OCC_CONFUSION &&
         h1>h0+_OCC_CONFUSION) || throw(ArgumentError(
            "$caller: CONICAL_SURFACE #$(surf.id) does not bound a " *
            "finite slant segment"))
        return ((occ=:cone,
                 center=_occ3_add(loc,_occ3_scale(n,h0)),axis=n,X=X,Y=Y,
                 r1=max(r1,0.0),r2=r2,height=h1-h0,
                 reversed=!sense),:cone)
    elseif surf.kind=="SPHERICAL_SURFACE"
        loc,n,X,Y,R=_occ_import_conic(entities,surf,caller)
        return ((occ=:sphere,center=loc,radius=R,axis=n,X=X,Y=Y,
                 reversed=!sense),:sphere)
    elseif surf.kind=="TOROIDAL_SURFACE"
        loc,n,X,Y,(r1,r2)=_occ_import_conic(entities,surf,caller;r2=true)
        return ((occ=:torus,center=loc,axis=n,X=X,Y=Y,r1=r1,r2=r2,
                 angle=_OCC_TWO_PI,reversed=!sense),:torus)
    else
        ns=BRep._step_nurbs_surface(surf,entities)
        ns===nothing && throw(ArgumentError(
            "$caller: $(surf.kind) #$(surf.id) — only PLANE/" *
            "CYLINDRICAL_SURFACE/CONICAL_SURFACE/SPHERICAL_SURFACE/" *
            "TOROIDAL_SURFACE/B-spline surfaces import currently"))
        return ((occ=:tensor_bspline,
                 patch=_ModelNURBSSurface(:bspline,ns),
                 reversed=!sense),:tensor_bspline)
    end
end

# Position of a bound vertex shape id.
_occ_import_vertex_pos(kernel,vid)=kernel.shapes[vid].geometry.point

# The structural face: STEP bounds → wires → provisional surface →
# pcurves/vertex_uv → bound face shape. Returns (ftag, shell_rev) where
# `shell_rev` is the use sign the solid's shells must carry (`true` when
# the STEP `same_sense` flips the face normal off the surface's).
function _occ_import_step_face_structural(m,d,st,entities,fe,bounds,surf,
                                          sense,curvekinds,caller)
    kernel=m.occ_kernel
    wire_uses=Vector{Tuple{UInt64,Bool}}[]
    outers=Int[]
    for bref in bounds
        be=BRep._deref_entity(entities,bref)
        (be!==nothing && be.kind in BRep._STEP_BREP_BOUND_KINDS) || throw(
            ArgumentError("$caller: $(fe.kind) #$(fe.id) boundary is " *
                          "not a FACE_BOUND/FACE_OUTER_BOUND"))
        length(be.args)>=3 || throw(ArgumentError(
            "$caller: $(be.kind) #$(be.id) is truncated"))
        loop=BRep._deref_entity(entities,be.args[2])
        loop===nothing && throw(ArgumentError(
            "$caller: $(be.kind) #$(be.id) references a missing loop"))
        flag=BRep._brep_enum(be.args[3]); flag===nothing && (flag=true)
        uses=_occ_import_loop_uses(m,d,st,entities,loop,flag,curvekinds,
                                   caller)
        _occ_import_check_chain(kernel,uses,caller)
        push!(wire_uses,uses)
        be.kind=="FACE_OUTER_BOUND" && push!(outers,length(wire_uses))
    end
    if length(wire_uses)>1
        length(outers)<=1 || throw(ArgumentError(
            "$caller: holed face #$(fe.id) declares multiple " *
            "FACE_OUTER_BOUND wires"))
        isempty(outers) || (wire_uses=vcat(wire_uses[outers[1]:outers[1]],
            wire_uses[1:outers[1]-1],wire_uses[outers[1]+1:end]))
    end
    # Boundary vertices feed the analytic surface extent fit and the
    # vertex_uv record; the outer wire's ordered traversal ring frames
    # untyped `FACE` records.
    vpos=NTuple{3,Float64}[]
    vids=UInt64[]
    ring=NTuple{3,Float64}[]
    for uses in wire_uses,(eid,rev) in uses
        a,b=kernel.shapes[eid].children
        first_vid,last_vid=rev ? (b,a) : (a,b)
        for vid in (first_vid,last_vid)
            vid in vids && continue
            push!(vids,vid)
            push!(vpos,_occ_import_vertex_pos(kernel,vid))
        end
    end
    for (eid,rev) in wire_uses[1]
        a,b=kernel.shapes[eid].children
        push!(ring,_occ_import_vertex_pos(kernel,rev ? b : a))
    end
    base_g,fkind=_occ_import_surface_base(entities,surf,sense,vpos,ring,
                                          caller)
    # Plane records carry `uvb` — `_occ_surface_closest` reads it during
    # pcurve projection, so seed it from the boundary vertices now; the
    # pcurve-derived bounds overwrite it below.
    if fkind===:plane
        umin=vmin=Inf; umax=vmax=-Inf
        for p in vpos
            u=BRep._brep_vdot(BRep._brep_vsub(p,base_g.center),base_g.X)
            v=BRep._brep_vdot(BRep._brep_vsub(p,base_g.center),base_g.Y)
            umin=min(umin,u); umax=max(umax,u)
            vmin=min(vmin,v); vmax=max(vmax,v)
        end
        pad=1e-6*max(1.0,umax-umin,vmax-vmin)
        base_g=merge(base_g,(uvb=([umin-pad,vmin-pad],
                                  [umax+pad,vmax+pad]),))
    end
    # pcurves — each non-degenerate edge projects once (canonical uv);
    # the chain pass then shifts every use's copy into the wire's
    # unwrapped window, and degenerate edges become straight connectors.
    period_u=_brep_uperiod(base_g)
    period_v=_occ_surface_isvperiodic(base_g) ? _OCC_TWO_PI : 0.0
    nurbs=_model_nurbs_surface(base_g)
    base_pc=Dict{UInt64,NamedTuple}()
    nurbs || for uses in wire_uses,(eid,_) in uses
        haskey(base_pc,eid) && continue
        eg=kernel.shapes[eid].geometry.geometry
        eg.occ===:degenerate && continue
        base_pc[eid]=if fkind===:plane &&
            (eg.occ===:line ||
             (eg.occ in (:circle,:ellipse) &&
              BRep._brep_vdot(eg.n,base_g.axis)>0.0))
            _model_occ_plane_pcurve(eg,base_g.center,base_g.X,base_g.Y)
        else
            _trim_project_pcurve(m,kernel,base_g,eid,caller)
        end
    end
    pcurves=Dict{Int,NamedTuple}()
    trange=Dict{Int,Tuple{Float64,Float64}}()
    deg_spans=Dict{UInt64,Float64}()
    use_pc=Dict{Tuple{UInt64,Bool},NamedTuple}()
    for uses in wire_uses
        slot=Dict{Tuple{UInt64,Bool},NamedTuple}()
        uprev=nothing
        pending=Int[]          # degenerate-use indices awaiting a successor
        pending_start=nothing  # uv where the pending run began
        function flush_pending!(estart)
            prev=pending_start===nothing ? estart : pending_start
            for di in pending
                did,drev=uses[di]
                dpc,span=_occ_import_degenerate_pcurve(prev,estart,drev)
                slot[(did,drev)]=dpc
                use_pc[(did,drev)]=dpc
                deg_spans[did]=span
                trange[_model_occ_shape_tag(kernel,did,caller)]=(0.0,span)
                prev=estart
            end
            empty!(pending)
            return nothing
        end
        # One wire pass. For NURBS faces each use's pcurve is projected
        # hint-chained in traversal order; for analytic bases the
        # canonical projection is period-shifted into the running window.
        # `shift_extra` adds a fixed uv offset to one use's copy — the
        # seam-pair repair below uses it to move a retraversed seam copy
        # a full period. `rev_extra` swaps a use's copy for its
        # t-reversed twin — a closed edge wired through the second rim
        # of a slit face traverses its window backwards.
        # Returns (wire-start uv, wire-end uv).
        function chain_pass!(anchor;shift_extra=Dict{Int,Float64}(),
                             rev_extra=Set{Int}())
            empty!(slot); empty!(pending); pending_start=nothing
            uprev=anchor; first_start=nothing
            for (ui,(eid,rev)) in enumerate(uses)
                eg=kernel.shapes[eid].geometry.geometry
                if eg.occ===:degenerate
                    isempty(pending) && (pending_start=uprev)
                    push!(pending,ui)
                    continue
                end
                pc=if nurbs
                    pcp=_occ_import_project_pcurve(m,kernel,base_g,eid,
                                                 uprev,rev,caller)
                    ui in rev_extra ?
                        _occ_import_pc_reverse(pcp,eg.t0,eg.t1,caller) :
                        pcp
                else
                    pc0=base_pc[eid]
                    ui in rev_extra && (pc0=_occ_import_pc_reverse(
                        pc0,eg.t0,eg.t1,caller))
                    s0=rev ? _occ_pcurve_d0(pc0,eg.t1) :
                        _occ_pcurve_d0(pc0,eg.t0)
                    du=period_u>0.0 && uprev!==nothing ?
                        period_u*round((uprev[1]-s0[1])/period_u) : 0.0
                    dv=period_v>0.0 && uprev!==nothing ?
                        period_v*round((uprev[2]-s0[2])/period_v) : 0.0
                    _occ_import_pc_shift(pc0,du,dv)
                end
                extra=get(shift_extra,ui,0.0)
                extra!=0.0 && (pc=_occ_import_pc_shift(pc,extra,0.0))
                st0=rev ? _occ_pcurve_d0(pc,eg.t1) :
                    _occ_pcurve_d0(pc,eg.t0)
                slot[(eid,rev)]=pc
                use_pc[(eid,rev)]=pc
                isempty(pending) ||
                    (flush_pending!(st0); pending_start=nothing)
                first_start===nothing && (first_start=st0)
                eend=rev ? _occ_pcurve_d0(pc,eg.t0) :
                    _occ_pcurve_d0(pc,eg.t1)
                uprev=eend
                trange[_model_occ_shape_tag(kernel,eid,caller)]=
                    (eg.t0,eg.t1)
            end
            # A degenerate run reaching the wire end closes onto the
            # first use's trace start.
            if !isempty(pending)
                eid,rev=uses[1]
                eg=kernel.shapes[eid].geometry.geometry
                estart=if haskey(slot,(eid,rev))
                    pc=slot[(eid,rev)]
                    rev ? _occ_pcurve_d0(pc,eg.t1) :
                        _occ_pcurve_d0(pc,eg.t0)
                else
                    uprev
                end
                estart===nothing && throw(ArgumentError(
                    "$caller: cannot anchor a degenerate edge's uv " *
                    "trace"))
                flush_pending!(estart)
            end
            return first_start,uprev
        end
        chain_anchor=nothing
        if nurbs
            first_start,uprev=chain_pass!(nothing)
            # Wire-closure re-anchor: when the first use's canonical
            # projection picked a different symmetry copy than the
            # chain-end, re-project it anchored at the wire's end uv —
            # the closing vertex's copy — and re-chain the wire. (A
            # degenerate last use legitimately leaves `uprev` short of
            # the start — its connector bridges the pole gap.)
            last_deg=kernel.shapes[uses[end][1]].
                geometry.geometry.occ===:degenerate
            if !last_deg && first_start!==nothing && uprev!==nothing
                sc=max(1.0,abs(first_start[1]),abs(first_start[2]),
                       abs(uprev[1]),abs(uprev[2]))
                gap=max(abs(uprev[1]-first_start[1]),
                        abs(uprev[2]-first_start[2]))
                if gap>1e-5*sc
                    first_start,uprev=chain_pass!(uprev)
                    chain_anchor=uprev
                    gap=max(abs(uprev[1]-first_start[1]),
                            abs(uprev[2]-first_start[2]))
                    gap>1e-4*sc && throw(ArgumentError(
                        "$caller: face #$(fe.id) wire pcurves do not " *
                        "close on the NURBS surface (start=" *
                        "$(first_start) end=$(uprev) gap=$(gap))"))
                end
            end
        else
            chain_pass!(nothing)
        end
        # Periodic-seam repair — an edge used in both orientations in one
        # wire is the face's seam: its two use copies bound the face's uv
        # window at opposite period ends. The greedy nearest-copy anchor
        # above can collapse both uses onto the same window (the copies
        # then retrace one another and the sampled wire self-overlaps —
        # a slit wire like a full-cone lateral reports duplicate CDT
        # constraints, a sphere's seam-meridian pair crosses them).
        # Re-chain with the second use shifted one u-period, keeping the
        # sign that leaves the smallest closure gap.
        if period_u>0.0
            first_use=Dict{UInt64,Int}()
            shift_extra=Dict{Int,Float64}()
            for (ui,(eid,rev)) in enumerate(uses)
                eg=kernel.shapes[eid].geometry.geometry
                eg.occ===:degenerate && continue
                prev=get(first_use,eid,nothing)
                prev===nothing && (first_use[eid]=ui; continue)
                uses[prev][2]==rev && continue
                pc_a=get(slot,(eid,!rev),nothing)
                pc_b=get(slot,(eid,rev),nothing)
                (pc_a===nothing || pc_b===nothing) && continue
                qa=_occ_pcurve_d0(pc_a,eg.t0); qb=_occ_pcurve_d0(pc_b,eg.t0)
                ra=_occ_pcurve_d0(pc_a,eg.t1); rb=_occ_pcurve_d0(pc_b,eg.t1)
                vs=max(1.0,abs(qa[2]),abs(ra[2]),abs(qb[2]),abs(rb[2]))
                (abs(qb[1]-qa[1])<0.5*period_u &&
                 abs(rb[1]-ra[1])<0.5*period_u &&
                 abs(qb[2]-qa[2])<=1e-6*vs &&
                 abs(rb[2]-ra[2])<=1e-6*vs) || continue
                best=0.0;bestgap=Inf
                for sgn in (1.0,-1.0)
                    shift_extra[ui]=sgn*period_u
                    wstart,wend=chain_pass!(chain_anchor;shift_extra)
                    if wstart!==nothing && wend!==nothing
                        sc=max(1.0,abs(wstart[1]),abs(wstart[2]),
                               abs(wend[1]),abs(wend[2]))
                        gap=max(abs(wend[1]-wstart[1]),
                                abs(wend[2]-wstart[2]))/sc
                        gap<bestgap && (best=sgn*period_u;bestgap=gap)
                    end
                end
                best==0.0 ? delete!(shift_extra,ui) : (shift_extra[ui]=best)
            end
            isempty(shift_extra) ||
                chain_pass!(chain_anchor;shift_extra)
            # Wrapped-window repair — a closed edge (both endpoints the
            # same vertex, like a rim circle) on a periodic face can run
            # its uv trace either direction around the window; the
            # greedy anchor keeps the forward sense, so a slit wire like
            # a frustum's lateral face wraps the window twice — the
            # second rim circle retraverses the first's strip and the
            # wire's closing edge overlaps the sampled trace (duplicate
            # CDT constraints, doubled face mass). When the wire's
            # u-extent exceeds one period, try each closed use's
            # t-reversed copy and keep the reversal that best shrinks
            # the span; repeat while a single reversal still helps.
            function wire_uspan()
                lo=Inf; hi=-Inf
                for (eid2,rev2) in uses
                    pc=get(slot,(eid2,rev2),nothing)
                    pc===nothing && continue
                    eg2=kernel.shapes[eid2].geometry.geometry
                    for i in 0:8
                        u=_occ_pcurve_d0(pc,
                            eg2.t0+(eg2.t1-eg2.t0)*i/8)[1]
                        u<lo && (lo=u); u>hi && (hi=u)
                    end
                end
                return hi-lo
            end
            uspan=wire_uspan()
            if uspan>1.5*period_u
                closed=[ui for (ui,(eid,rev)) in enumerate(uses)
                        if kernel.shapes[eid].geometry.geometry.occ!==
                           :degenerate &&
                           length(kernel.shapes[eid].children)>=2 &&
                           kernel.shapes[eid].children[1]==
                           kernel.shapes[eid].children[2]]
                rev_extra=Set{Int}()
                while uspan>1.5*period_u && !isempty(closed)
                    best_ui=0; best=uspan
                    for ui in closed
                        chain_pass!(chain_anchor;shift_extra,
                                    rev_extra=union(rev_extra,(ui,)))
                        s2=wire_uspan()
                        s2<best && (best=s2;best_ui=ui)
                    end
                    best_ui==0 && break
                    push!(rev_extra,best_ui); uspan=best
                    filter!(!=(best_ui),closed)
                end
                chain_pass!(chain_anchor;shift_extra,rev_extra)
            end
        end
        # Slot bookkeeping: both-direction uses split into fwd/rev;
        # single-direction uses always publish `fwd` so the shared
        # wire-polygon passes see the trace.
        dirs=Dict{UInt64,Vector{Bool}}()
        for (eid,rev) in uses
            push!(get!(dirs,eid,Bool[]),rev)
        end
        for (eid,_) in dirs
            etag=_model_occ_shape_tag(kernel,eid,caller)
            fwd_use=get(slot,(eid,false),nothing)
            rev_use=get(slot,(eid,true),nothing)
            pcurves[etag]=if fwd_use!==nothing && rev_use!==nothing
                (fwd=fwd_use,rev=rev_use)
            elseif fwd_use!==nothing
                (fwd=fwd_use,rev=nothing)
            else
                (fwd=rev_use,rev=nothing)
            end
        end
    end
    # vertex_uv for the patch path (`_occ_vertex_reparam_on_face` reads
    # it only on `_model_nurbs_surface` records).
    extra=NamedTuple()
    if _model_nurbs_surface(base_g)
        # Vertex uvs come from the use copies so they sit in the wire's
        # own uv window — a fresh closest-point solve could return an
        # equivalent point in another symmetry copy.
        vuv=Dict{Int,NTuple{2,Float64}}()
        for uses in wire_uses,(eid,rev) in uses
            pc=get(use_pc,(eid,rev),nothing)
            pc===nothing && continue
            eg=kernel.shapes[eid].geometry.geometry
            a,b=kernel.shapes[eid].children
            for (vid,tp) in ((a,eg.t0),(b,eg.t1))
                vtag=_model_occ_shape_tag(kernel,vid,caller)
                haskey(vuv,vtag) && continue
                vuv[vtag]=_occ_pcurve_d0(pc,tp)
            end
        end
        for vid in vids
            vtag=_model_occ_shape_tag(kernel,vid,caller)
            haskey(vuv,vtag) && continue
            pt=kernel.shapes[vid].geometry.point
            pr=_occ_surface_closest(base_g,pt,1.0)
            pr===nothing && throw(ArgumentError(
                "$caller: wire vertex does not project onto the surface"))
            vuv[vtag]=(pr[1],pr[2])
        end
        extra=(vertex_uv=vuv,)
    end
    # Plane records carry `uvb` (their canonical sampling domain); other
    # kinds pull bounds from `_occ_surface_bounds`.
    if fkind===:plane
        extra=merge(extra,(uvb=_model_occ_prim_uvb(pcurves,trange),))
    elseif fkind!==:plane
        umin=vmin=Inf; umax=vmax=-Inf
        for (etag,e) in pcurves,(pc) in (e.fwd,e.rev)
            pc===nothing && continue
            t0,t1=trange[abs(etag)]
            lo,hi=_occ_pcurve_fillbndbox(pc,t0,t1)
            umin=min(umin,lo[1]); vmin=min(vmin,lo[2])
            umax=max(umax,hi[1]); vmax=max(vmax,hi[2])
        end
        lo_b,hi_b=_occ_surface_bounds(base_g)
        covers=umin<=lo_b[1]+1e-9 && umax>=hi_b[1]-1e-9 &&
               vmin<=lo_b[2]+1e-9 && vmax>=hi_b[2]-1e-9
        # Periodic faces whose wire only touches a period end still
        # cover the domain — compare modulo the period.
        if !covers && period_u>0.0 && umax-umin>=period_u-1e-9
            covers=vmin<=lo_b[2]+1e-9 && vmax>=hi_b[2]-1e-9
        end
        trimmed=length(wire_uses)>1 || !covers
        trimmed && (extra=merge(extra,(trimmed=true,)))
    end
    geometry=merge(base_g,(pcurves=pcurves,),extra)
    # `_occ_degenerate_point` reads the edge's `face` record for its
    # connector pcurve — point it at this face's final geometry.
    for (did,span) in deg_spans
        _occ_import_patch_degenerate!(kernel,did,span,geometry,caller)
    end
    wire_tags=[_occ_wire_from_uses!(m,uses,caller) for uses in wire_uses]
    wire_ids=map(wire_tags) do t
        id=get(kernel.bindings,(-1,t),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: internal — OCC Wire[$t] unbound"))
        id
    end
    ftag=_model_occ_next_tag(kernel,2,-1,caller)
    fid=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,2,ftag,fid,
        _ModelOCCShape(2,(kind=fkind,geometry=geometry),wire_ids))
    kernel.next_shape=fid+UInt64(1)
    return ftag,!sense
end

# Give a degenerate edge its uv span and the face record it lives on —
# `_occ_degenerate_point` reads `g.face.pcurves[tag].fwd` and evaluates
# the connector at the edge parameter.
function _occ_import_patch_degenerate!(kernel,eid,span,face_g,caller)
    sh=kernel.shapes[eid]
    eg=sh.geometry.geometry
    kernel.shapes[eid]=_ModelOCCShape(1,
        merge(sh.geometry,(geometry=merge(eg,(t1=span,face=face_g)),)),
        sh.children)
    return nothing
end

# Face routing preamble — the ring path's `BRep._step_face_parts`
# cannot be probed first: its `_chain_ring` certificate throws on a
# legal two-edge curved loop. Resolve (bounds, surf, sense) and whether
# every bound is poly/straight here, without walking the chain.
function _occ_import_face_bounds(entities,fe,caller)
    i=BRep._skip_name(fe.args)
    length(fe.args)>=i || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) has no boundary list"))
    bounds=fe.args[i]
    bounds isa Vector || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) boundary list is malformed"))
    isempty(bounds) && throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) has no face boundaries"))
    surf=nothing; sense=true
    if fe.kind!="FACE"
        length(fe.args)>=i+2 || throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) is missing its surface or " *
            "sense flag"))
        surf=BRep._deref_entity(entities,fe.args[i+1])
        surf===nothing && throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) references a missing surface"))
        flag=BRep._brep_enum(fe.args[i+2])
        flag===nothing || (sense=flag)
    end
    return bounds,surf,sense
end

# True when every bound is a POLY_LOOP or an EDGE_LOOP over straight
# edges — the condition under which the vertex-ring path applies.
function _occ_import_face_all_straight(entities,bounds,curvekinds,caller)
    for bref in bounds
        be=BRep._deref_entity(entities,bref)
        (be!==nothing && be.kind in BRep._STEP_BREP_BOUND_KINDS) || throw(
            ArgumentError("$caller: face boundary is not a " *
                          "FACE_BOUND/FACE_OUTER_BOUND"))
        length(be.args)>=3 || throw(ArgumentError(
            "$caller: $(be.kind) #$(be.id) is truncated"))
        loop=BRep._deref_entity(entities,be.args[2])
        loop===nothing && throw(ArgumentError(
            "$caller: $(be.kind) #$(be.id) references a missing loop"))
        loop.kind=="POLY_LOOP" && continue
        loop.kind=="EDGE_LOOP" || return false
        i=BRep._skip_name(loop.args)
        i>length(loop.args) && return false
        refs=loop.args[i]
        refs isa Vector || return false
        for r in refs
            oe=BRep._deref_entity(entities,r)
            (oe!==nothing && oe.kind=="ORIENTED_EDGE") || return false
            ecref=nothing
            for a in oe.args
                a isa Tuple && a[1]===:ref && (ecref=a)
            end
            ec=BRep._deref_entity(entities,ecref)
            (ec!==nothing && ec.kind=="EDGE_CURVE" &&
             length(ec.args)>=4) || return false
            cv=BRep._deref_entity(entities,ec.args[4])
            BRep._step_straight_edge(entities,cv,0,curvekinds) ||
                return false
        end
    end
    return true
end

# One STEP face → (ftag, shell_rev). Straight-edged planar faces keep
# the ring path (identical to earlier revisions); everything else walks
# the structural conversion.
function _occ_import_step_face(m,d,st,kernel,entities,fe,curvekinds,caller)
    bounds,surf,sense=_occ_import_face_bounds(entities,fe,caller)
    if (surf===nothing || surf.kind=="PLANE") &&
        _occ_import_face_all_straight(entities,bounds,curvekinds,caller)
        rings,_,_,_=BRep._step_face_parts(entities,fe,curvekinds,caller)
        rings===nothing && throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) uses an unsupported loop kind"))
        wires=[_occ_import_wire!(m,d,kernel,r,false,caller) for r in rings]
        return _model_occ_add_plane_surface!(m,wires),false
    end
    return _occ_import_step_face_structural(m,d,st,entities,fe,bounds,
                                            surf,sense,curvekinds,caller)
end

# Bind one shell with explicit per-face use signs — `rev` carries the
# STEP `same_sense` flip the divergence/Green paths fold into the face
# normal.
function _occ_import_bind_shell!(m,uses,caller)
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,-2,-1,caller)
    face_ids=UInt64[fid for (fid,_) in uses]
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,-2,allocated,id,
        _ModelOCCShape(-2,(uses=uses,),face_ids))
    kernel.next_shape=id+UInt64(1)
    return allocated
end

# One solid root → [(3, tag)]. Each CLOSED_SHELL becomes a surface
# loop; `BREP_WITH_VOIDS`-style later shells carry role −1 (cavities).
function _occ_import_step_solid(m,brep,entities,caller)
    kernel=m.occ_kernel
    length(brep.args)>=2 || throw(ArgumentError(
        "$caller: $(brep.kind) #$(brep.id) has no shell"))
    d=_OCCImportDedupe(); st=_OCCImportStructural()
    curvekinds=Set{String}()
    shell_tags=Int[]
    for (si,sarg) in enumerate(brep.args[2:end])
        items=sarg isa Vector ? sarg : Any[sarg]
        for item in items
            se=BRep._deref_entity(entities,item)
            se===nothing && throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references a missing shell"))
            se.kind=="CLOSED_SHELL" || throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references $(se.kind) " *
                "#$(se.id) where a CLOSED_SHELL was expected"))
            i=BRep._skip_name(se.args)
            flist=se.args[i]
            flist isa Vector || throw(ArgumentError(
                "$caller: CLOSED_SHELL #$(se.id) face list is malformed"))
            uses=Tuple{UInt64,Bool}[]
            for fref in flist
                fe=BRep._deref_entity(entities,fref)
                (fe!==nothing &&
                 fe.kind in BRep._STEP_BREP_FACE_KINDS) || throw(
                    ArgumentError(
                        "$caller: CLOSED_SHELL #$(se.id) member is not a face"))
                ftag,rev=_occ_import_step_face(m,d,st,kernel,entities,fe,
                                               curvekinds,caller)
                fid=get(kernel.bindings,(2,ftag),nothing)
                fid===nothing && throw(ArgumentError(
                    "$caller: internal — face tag $ftag unbound"))
                push!(uses,(fid,rev))
            end
            isempty(uses) && throw(ArgumentError(
                "$caller: CLOSED_SHELL #$(se.id) is empty"))
            push!(shell_tags,_occ_import_bind_shell!(m,uses,caller))
        end
    end
    isempty(shell_tags) && throw(ArgumentError(
        "$caller: $(brep.kind) #$(brep.id) has no shells"))
    return [(3,_model_occ_add_volume!(m,shell_tags))]
end

function _model_occ_import_shapes(m::GeoModel,path;highest_dim_only=true,
                                  caller="OCC.importShapes")
    isfile(path) || throw(ArgumentError("$caller: missing file $path"))
    ext=lowercase(splitext(path)[2])
    if ext in (".igs",".iges")
        return _model_occ_import_iges(m,path;
            highest_dim_only=highest_dim_only,caller=caller)
    end
    source=read(path,String)
    occursin("ISO-10303-21",source) || throw(ArgumentError(
        "$caller: $path is not an ISO-10303-21 STEP file"))
    kernel=m.occ_kernel
    # Native entities already present — the GEO OCC route runs this
    # mid-script — push the kernel's tag maxima past the model's so
    # imported tags cannot collide.
    for d in 0:3
        kernel.maxima[d+3]=max(kernel.maxima[d+3],m.next_tag[d+1])
    end
    isempty(m.loops) ||
        (kernel.maxima[2]=max(kernel.maxima[2],maximum(keys(m.loops))))
    isempty(m.surface_loops) ||
        (kernel.maxima[1]=max(kernel.maxima[1],
                              maximum(keys(m.surface_loops))))
    entities=BRep.parse_step_entities(source)
    breps=sort!([e for e in values(entities)
                 if e.kind in BRep._STEP_BREP_SOLID_KINDS]; by=e->e.id)
    isempty(breps) && throw(ArgumentError(
        "$caller: no MANIFOLD_SOLID_BREP/FACETED_BREP/BREP_WITH_VOIDS " *
        "found in $path"))
    out=Tuple{Int,Int}[]
    for brep in breps
        append!(out,_occ_import_step_solid(m,brep,entities,caller))
    end
    return out
end

# ── `healShapes` ──────────────────────────────────────────────────────────────
#
# `OCC_Internals::healShapes`: weld coincident-but-distinct vertices and
# edges in the target entities' topology within `tolerance`
# (`Precision::Confusion` default), then — when `sewFaces` — regroup
# healed faces into shells and bind solids. Healing works on the stored
# graph: vertex entities cluster by position, edge children rewire to
# the cluster's canonical vertex, and line edges sharing both endpoints
# collapse to one curve entity. Faces/curves are not re-fit — tolerance
# healing moves topology, not geometry.
# Reverse a pcurve onto an oppositely-directed replacement edge:
# `pc'(s)=pc(T−s)` with `T=t0+t1` over the shared `[0,L]` line range —
# lin2d via affine reparam, spl2d by knot-mirror + pole reversal.
function _occ_heal_pcurve_reverse(pc,T::Float64,kernel,caller)
    if hasproperty(pc,:neg2d)
        inner=pc.neg2d[1]
        hasproperty(inner,:lin2d) || throw(ArgumentError(caller*
            ": cannot weld edge — unsupported negated pcurve kind"))
        o,d=inner.lin2d
        return (neg2d=((lin2d=((o[1]-T*d[1],o[2]-T*d[2]),d),),))
    end
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=((o[1]+T*d[1],o[2]+T*d[2]),(-d[1],-d[2])),)
    elseif hasproperty(pc,:spl2d)
        knots,p,cps=pc.spl2d
        return (spl2d=(T .- reverse(knots),p,reverse(cps)),)
    end
    throw(ArgumentError(caller*
        ": cannot weld edge — unsupported pcurve kind for reversal"))
end

function _model_occ_heal_shapes(m::GeoModel,dim_tags;tolerance=-1.0,
                                fixDegenerated=true,fixSmallEdges=true,
                                fixSmallFaces=true,sewFaces=true,
                                makeSolids=true,
                                caller="OCC.healShapes")
    kernel=m.occ_kernel
    tol=tolerance>0.0 ? tolerance : _OCC_CONFUSION
    roots=Tuple{Int,Int}[]
    if isempty(dim_tags)
        for (d,t) in keys(kernel.bindings)
            d>=0 && push!(roots,(d,t))
        end
    else
        for (d,t) in dim_tags
            d=Int(d); t=abs(Int(t))
            haskey(kernel.bindings,(d,t)) || throw(ArgumentError(
                "$caller: unknown OCC entity ($d,$t)"))
            push!(roots,(d,t))
        end
    end
    # Reachable-vertex weld: gather every vertex id under the roots
    # (shapes' children), cluster positions within `tol`, rewire edges.
    vids=UInt64[]
    stack=UInt64[kernel.bindings[r] for r in roots]
    seen=Set{UInt64}()
    while !isempty(stack)
        id=pop!(stack)
        id in seen && continue; push!(seen,id)
        sh=kernel.shapes[id]
        if sh.dimension==0
            push!(vids,id); continue
        end
        append!(stack,sh.children)
        hasproperty(sh.geometry,:uses) || continue
        for (uid,_) in sh.geometry.uses
            uid in seen || push!(stack,uid)
        end
    end
    # Cluster by position — canonical = first member.
    canon=Dict{UInt64,UInt64}()
    reps=UInt64[]
    for id in vids
        p=kernel.shapes[id].geometry.point
        hit=0
        for r in reps
            _occ_modulus(_occ3_sub(kernel.shapes[r].geometry.point,p))<=tol &&
                (hit=r; break)
        end
        if hit==0
            push!(reps,id); canon[id]=id
        else
            canon[id]=hit
        end
    end
    # Rewire edge children to canonical vertices.
    for id in seen
        sh=kernel.shapes[id]
        sh.dimension==1 || continue
        ch=sh.children
        nc=UInt64[get(canon,c,c) for c in ch]
        nc==ch && continue
        kernel.shapes[id]=_ModelOCCShape(1,sh.geometry,nc)
    end
    # Edge dedupe: two lines sharing both endpoints collapse — faces
    # then share the curve entity (needed for shell validation).
    ecanon=Dict{UInt64,UInt64}()
    if fixSmallEdges || sewFaces
        sig=Dict{NTuple{2,UInt64},UInt64}()
        for id in seen
            sh=kernel.shapes[id]
            sh.dimension==1 || continue
            get(sh.geometry.geometry,:occ,nothing)==:line || continue
            a,b=sh.children; key=a<b ? (a,b) : (b,a)
            prev=get(sig,key,0)
            if prev==0
                sig[key]=id
            else
                ecanon[id]=prev
            end
        end
        if !isempty(ecanon)
            # Re-key face pcurves while wires still reference the
            # face's own edges: each pcurve's key moves from the old
            # edge tag to the canonical edge's tag. When the canonical
            # edge runs opposite the dropped one, `fwd` — parametrized
            # along the edge's own direction — must be reversed.
            for id in seen
                sh=kernel.shapes[id]
                sh.dimension==2 || continue
                fr=sh.geometry.geometry
                hasproperty(fr,:pcurves) || continue
                pc=Dict{Int,NamedTuple}()
                for cid in sh.children
                    wsh=kernel.shapes[cid]
                    hasproperty(wsh.geometry,:uses) || continue
                    for (u,_) in wsh.geometry.uses
                        ot=kernel.shape_tags[u][2]
                        nu=get(ecanon,u,u)
                        nt=kernel.shape_tags[nu][2]
                        ent=fr.pcurves[ot]
                        if nu!=u && kernel.shapes[u].children[1]!=
                                kernel.shapes[nu].children[1]
                            eg=kernel.shapes[u].geometry.geometry
                            ent=(fwd=_occ_heal_pcurve_reverse(ent.fwd,
                                     eg.t0+eg.t1,kernel,caller),
                                 rev=nothing)
                        end
                        pc[nt]=ent
                    end
                end
                kernel.shapes[id]=_ModelOCCShape(2,
                    merge(sh.geometry,(geometry=merge(fr,(pcurves=pc,)),)),
                    sh.children)
            end
            # Rewire wire/shell `uses` to canonical edges — flipping
            # `rev` when the dropped edge's stored direction opposes
            # the canonical edge's.
            for id in seen
                sh=kernel.shapes[id]
                hasproperty(sh.geometry,:uses) || continue
                uses=Tuple{UInt64,Bool}[]
                for (u,r) in sh.geometry.uses
                    nu=get(ecanon,u,u)
                    if nu==u
                        push!(uses,(u,r))
                    else
                        och=kernel.shapes[u].children
                        nch=kernel.shapes[nu].children
                        push!(uses,(nu,och[1]==nch[1] ? r : !r))
                    end
                end
                uses==sh.geometry.uses && continue
                kernel.shapes[id]=_ModelOCCShape(sh.dimension,
                    merge(sh.geometry,(uses=uses,)),
                    UInt64[u for (u,_) in uses])
            end
        end
    end
    # Unbind welded-away entities: duplicates no longer exist as
    # distinct geometry once every reference points at a canonical.
    for (id,rep) in Iterators.flatten((pairs(canon),pairs(ecanon)))
        id==rep && continue
        entry=get(kernel.shape_tags,id,nothing)
        entry===nothing && continue
        _model_occ_entity_map_unbind!(
            kernel.entity_maps[first(entry)+3],last(entry))
        delete!(kernel.bindings,entry);delete!(kernel.shape_tags,id)
        push!(kernel.removed,entry);kernel.changed=true
    end
    return nothing
end
