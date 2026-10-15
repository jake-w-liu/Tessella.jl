# `gmsh.model.occ.fillet2D`/`chamfer2D` — `ChFi2d_FilletAPI`/`ChFi2d_ChamferAPI`:
# two edges connected at a coincident corner vertex in one plane are trimmed
# back to the fillet tangency/chamfer points and a fresh joining edge is bound.
# Upstream replaces the two edge objects in place (their tags are rebound) and
# binds the joining edge under `tag`; vertices orphaned by the trim stay bound.

# Unit tangent at an edge endpoint pointing INTO the edge interior (`at_end`
# picks the t1 side). ChFi2d's analytic scope is lines and circles.
function _occ_2df_inward(g,at_end::Bool,caller)
    if g.occ===:line
        d=_occ_dir(g.dir)
        return at_end ? (-d[1],-d[2],-d[3]) : d
    elseif g.occ===:circle
        T=_occ_circle_derivative(g,at_end ? g.t1 : g.t0)
        u=_occ_dir(T)
        return at_end ? (-u[1],-u[2],-u[3]) : u
    end
    throw(ArgumentError(
        "$caller: only line and circle edges are supported — got $(g.occ)"))
end

function _occ_2df_eval(g,t)
    g.occ===:line && return _occ_line_point(g,t)
    g.occ===:circle && return _occ_circle_point(g,t)
    throw(ArgumentError("_occ_2df_eval: unsupported curve kind $(g.occ)"))
end

# The record eval dispatches to: native-backed edges evaluate through their
# nested `native_geometry` frame, whose X/Y may differ from the outer fields —
# parameter work must use the same frame.
_occ_2df_frame(g)=hasproperty(g,:native_geometry) ? g.native_geometry : g

# Parameter of a point on the record — line projection, circle angle
# branch-adjusted into the stored span.
function _occ_2df_param(g,p)
    if g.occ===:line
        d=g.dir; rel=_arc_sub(p,g.origin)
        return (rel[1]*d[1]+rel[2]*d[2]+rel[3]*d[3])/(d[1]*d[1]+d[2]*d[2]+d[3]*d[3])
    end
    f=_occ_2df_frame(g)
    rel=_arc_sub(p,f.center); Y=_occ_circle_y(f)
    θ=atan(_dot3(rel,Y),_dot3(rel,f.X))
    t=θ+_OCC_TWO_PI*round((g.t0-θ)/_OCC_TWO_PI)
    t<g.t0-1e-9 && (t+=_OCC_TWO_PI)
    t>g.t1+1e-9 && (t-=_OCC_TWO_PI)
    return t
end

# The cut must leave a positive-length remainder on the far side of the edge —
# `step` converts the 3-D tolerance into the record's parameter units.
function _occ_2df_trim_ok(g,at_end,tcut,tol)
    g.occ===:line && (tol/=max(_occ_modulus(g.dir),eps(Float64)))
    g.occ===:circle && (tol/=max(g.r,eps(Float64)))
    return at_end ? (g.t0+tol<tcut<g.t1) : (g.t0<tcut<g.t1-tol)
end

# Every coincident endpoint pair between the two edges — the corner candidates.
function _occ_2df_corners(shapes,id1,id2,tol)
    out=NamedTuple[]
    a1,b1=shapes[id1].children; a2,b2=shapes[id2].children
    for (end1,v1) in ((false,a1),(true,b1)),
        (end2,v2) in ((false,a2),(true,b2))
        p1=shapes[v1].geometry.point; p2=shapes[v2].geometry.point
        _model_occ_vertices_coincide(p1,p2) &&
            push!(out,(point=p1,end1=end1,end2=end2))
    end
    return out
end

_occ_2df_norm(v)=(m=_occ_modulus(v); (v[1]/m,v[2]/m,v[3]/m))

# Offset-line ∩ circle centers: C = A + λ·d̂ with |C−c| = ρ.
function _occ_2df_line_circle(A,d,radius,c,ρ)
    out=NTuple{3,Float64}[]
    rel=_arc_sub(A,c)
    λc=_dot3(rel,d)
    q2=_dot3(rel,rel)-λc*λc
    q2>ρ*ρ+1e-12*max(1.0,ρ*ρ) && return out
    h=sqrt(max(0.0,ρ*ρ-q2))
    for sgn in (1.0,-1.0)
        push!(out,_occ3_add(A,_occ3_scale(d,sgn*h-λc)))
    end
    return out
end

# Circle ∩ circle inside the corner plane.
function _occ_2df_circle_circle(c1,ρ1,c2,ρ2,n)
    out=NTuple{3,Float64}[]
    dvec=_arc_sub(c2,c1); d=_occ_modulus(dvec)
    d<1e-14 && return out                    # concentric: no isolated roots
    e=_occ3_scale(dvec,1.0/d)
    a=(d*d+ρ1*ρ1-ρ2*ρ2)/(2d)
    h2=ρ1*ρ1-a*a
    h2<-1e-12*max(1.0,ρ1*ρ1) && return out
    h=sqrt(max(0.0,h2))
    base=_occ3_add(c1,_occ3_scale(e,a))
    m=_occ_cross(n,e)
    for sgn in (1.0,-1.0)
        push!(out,_occ3_add(base,_occ3_scale(m,sgn*h)))
    end
    return out
end

# Tangency point of the fillet circle (center C) on the edge: the line foot,
# or the radial projection onto a circle; also its trim parameter.
function _occ_2df_tangent(g,C)
    if g.occ===:line
        t=_occ_2df_param(g,C)
        return _occ_line_point(g,t),t
    end
    rel=_arc_sub(C,g.center)
    T=_occ3_add(g.center,_occ3_scale(rel,g.r/_occ_modulus(rel)))
    return T,_occ_2df_param(g,T)
end

# Fillet solve: intersect the two tangent-center loci — a line contributes its
# offsets `X·m̂ = V·m̂ ± r`, a circle the rings `|X−c| = R±r`. A candidate is
# the corner fillet when its tangencies land strictly inside both edges on the
# corner sides AND the center sits in the corner wedge — the same side of each
# edge's corner tangent as the other edge's inward direction (`mi = ûi×n̂`
# is the in-plane normal of edge i's tangent, `u_j·mi` marks the side holding
# the wedge; `side` flips to the reflex sector like upstream's `reverse`).
# Nearest to the corner point wins (upstream iSolution = -1).
function _occ_2df_solve_fillet(g1,g2,u1,u2,n,V,corner,r,side,tol,caller)
    (isfinite(r) && r>0) || throw(ArgumentError(
        "$caller: radius must be positive"))
    m1=_occ_2df_norm(_occ_cross(u1,n)); m2=_occ_2df_norm(_occ_cross(u2,n))
    s1=side*_dot3(u2,m1); s2=side*_dot3(u1,m2)
    candidates=NTuple{3,Float64}[]
    isline1=g1.occ===:line; isline2=g2.occ===:line
    if isline1 && isline2
        ln1=_occ_2df_norm(_occ_cross(n,_occ_dir(g1.dir)))
        ln2=_occ_2df_norm(_occ_cross(n,_occ_dir(g2.dir)))
        for s1 in (1.0,-1.0), s2 in (1.0,-1.0)
            b=(_dot3(V,ln1)+s1*r,_dot3(V,ln2)+s2*r,_dot3(V,n))
            det=_dot3(ln1,_occ_cross(ln2,n))
            abs(det)<1e-14 && continue
            push!(candidates,
                _occ3_add(_occ3_scale(_occ_cross(ln2,n),b[1]/det),
                    _occ3_add(_occ3_scale(_occ_cross(n,ln1),b[2]/det),
                              _occ3_scale(_occ_cross(ln1,ln2),b[3]/det))))
        end
    elseif isline1 && !isline2
        d=_occ_dir(g1.dir); m=_occ_2df_norm(_occ_cross(n,d))
        for s in (1.0,-1.0), ρ in (g2.r+r,abs(g2.r-r))
            ρ<tol && continue
            append!(candidates,_occ_2df_line_circle(
                _occ3_add(V,_occ3_scale(m,s*r)),d,r,g2.center,ρ))
        end
    elseif !isline1 && isline2
        d=_occ_dir(g2.dir); m=_occ_2df_norm(_occ_cross(n,d))
        for s in (1.0,-1.0), ρ in (g1.r+r,abs(g1.r-r))
            ρ<tol && continue
            append!(candidates,_occ_2df_line_circle(
                _occ3_add(V,_occ3_scale(m,s*r)),d,r,g1.center,ρ))
        end
    else
        for ρ1 in (g1.r+r,abs(g1.r-r)), ρ2 in (g2.r+r,abs(g2.r-r))
            (ρ1<tol || ρ2<tol) && continue
            append!(candidates,
                _occ_2df_circle_circle(g1.center,ρ1,g2.center,ρ2,n))
        end
    end
    best=nothing; bestd=Inf
    for C in candidates
        # interior wedge: same side of each corner tangent as the other ray
        _dot3(_arc_sub(C,V),m1)*s1<=tol && continue
        _dot3(_arc_sub(C,V),m2)*s2<=tol && continue
        T1,tT1=_occ_2df_tangent(g1,C); T2,tT2=_occ_2df_tangent(g2,C)
        _occ_2df_trim_ok(g1,corner.end1,tT1,tol) || continue
        _occ_2df_trim_ok(g2,corner.end2,tT2,tol) || continue
        _dot3(_arc_sub(T1,V),u1)<=tol && continue      # cut the corner side
        _dot3(_arc_sub(T2,V),u2)<=tol && continue
        d=_sqdist(C,V)
        d<bestd && (best=(C,T1,tT1,T2,tT2); bestd=d)
    end
    best===nothing && throw(ArgumentError("$caller: Could not compute fillet"))
    C,T1,tT1,T2,tT2=best
    # The arc runs T1→T2 in the path's traversal sense: forward tangent at T1
    # equals the incoming travel direction −u1; the exit tangent must match
    # +u2 (guards against a same-tangent candidate looping the wrong way).
    ρ1=_occ3_scale(_arc_sub(T1,C),1.0/r); ρ2=_occ3_scale(_arc_sub(T2,C),1.0/r)
    nc=_occ_2df_norm(_occ_cross(u1,ρ1))
    tend=_occ_cross(nc,ρ2)
    _dot3(tend,u2)>0 || throw(ArgumentError("$caller: Could not compute fillet"))
    span=atan(_dot3(_occ_cross(ρ1,ρ2),nc),_dot3(ρ1,ρ2))
    span<=0 && (span+=_OCC_TWO_PI)
    arc=(occ=:circle,center=C,n=nc,X=ρ1,Y=_occ_cross(nc,ρ1),r=r,
         t0=0.0,t1=span)
    return (T1=T1,T2=T2,tT1=tT1,tT2=tT2,edge=arc,kind=:circle)
end

# Chamfer solve: the cut point sits `dist` along the edge from the corner —
# arc length for lines, angle·R for circles.
function _occ_2df_solve_chamfer(g1,g2,V,corner,d1,d2,tol,caller)
    (d1 isa Real && isfinite(d1) && d1>0) || throw(ArgumentError(
        "$caller: chamfer distances must be positive"))
    (d2 isa Real && isfinite(d2) && d2>0) || throw(ArgumentError(
        "$caller: chamfer distances must be positive"))
    cut(g,at_end,d)=begin
        step=g.occ===:line ? d/_occ_modulus(g.dir) : d/g.r
        at_end ? g.t1-step : g.t0+step
    end
    tT1=cut(g1,corner.end1,Float64(d1)); tT2=cut(g2,corner.end2,Float64(d2))
    _occ_2df_trim_ok(g1,corner.end1,tT1,tol) || throw(ArgumentError(
        "$caller: chamfer distance exceeds the first edge"))
    _occ_2df_trim_ok(g2,corner.end2,tT2,tol) || throw(ArgumentError(
        "$caller: chamfer distance exceeds the second edge"))
    T1=_occ_2df_eval(g1,tT1); T2=_occ_2df_eval(g2,tT2)
    diff=_arc_sub(T2,T1); len=_occ_modulus(diff)
    len<=tol && throw(ArgumentError("$caller: degenerate chamfer edge"))
    line=(occ=:line,origin=T1,dir=_occ3_scale(diff,1.0/len),t0=0.0,t1=len)
    return (T1=T1,T2=T2,tT1=tT1,tT2=tT2,edge=line,kind=:line)
end

# Publish the plan: bind the two fresh tangency vertices, replace the two
# trimmed edge records in place under their tags (record fields otherwise
# preserved), then the joining edge under `allocated` sharing the tangency
# vertices — upstream's `unbind`+`bind(…, recursive)` collapse.
function _occ_2df_apply!(kernel,id1,id2,g1,g2,corner,allocated,sol,caller)
    vT1,vT2,eid=_model_occ_reserve_shapes(kernel,3,caller)
    kernel.maxima[3]<=typemax(Int32)-2 || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    # Capture both tags up front — each bind bumps the vertex maximum.
    vtag1,vtag2=kernel.maxima[3]+1,kernel.maxima[3]+2
    _model_occ_bind_shape!(kernel,0,vtag1,vT1,
        _ModelOCCShape(0,(point=sol.T1,mesh_size=0.0),UInt64[]))
    _model_occ_bind_shape!(kernel,0,vtag2,vT2,
        _ModelOCCShape(0,(point=sol.T2,mesh_size=0.0),UInt64[]))
    for (id,g,at_end,tT,vnew) in
        ((id1,g1,corner.end1,sol.tT1,vT1),
         (id2,g2,corner.end2,sol.tT2,vT2))
        a,b=kernel.shapes[id].children
        rec=at_end ? merge(g,(t1=tT,)) : merge(g,(t0=tT,))
        # Trim the nested native frame's bounds too so evaluators reading it
        # stay consistent with the outer record.
        hasproperty(rec,:native_geometry) && (
            rec=merge(rec,(native_geometry=at_end ?
                merge(rec.native_geometry,(t1=tT,)) :
                merge(rec.native_geometry,(t0=tT,)),)))
        kids=at_end ? UInt64[a,vnew] : UInt64[vnew,b]
        kernel.shapes[id]=_ModelOCCShape(1,
            (kind=kernel.shapes[id].geometry.kind,geometry=rec),kids)
    end
    _model_occ_bind_shape!(kernel,1,allocated,eid,
        _ModelOCCShape(1,(kind=sol.kind,geometry=sol.edge),UInt64[vT1,vT2]))
    kernel.next_shape=eid+UInt64(1)
    return nothing
end

# Shared entry: resolve edges, locate the corner (nearest `point` — upstream's
# FirstVertex(ed1) default), fit the corner plane, solve, publish. Every
# fallible step precedes the first mutation.
function _model_occ_corner_op!(m::GeoModel,e1t,e2t;tag,point_tag,reverse,
                               caller,radius,dist1,dist2)
    kernel=m.occ_kernel; shapes=kernel.shapes; tol=_OCC_CONFUSION
    t1=_model_occ_requested_tag(e1t,caller)
    t2=_model_occ_requested_tag(e2t,caller)
    allocated=_model_occ_next_tag(kernel,1,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(1,allocated)) && throw(ArgumentError(
        "$caller: OCC curve with tag $allocated already exists"))
    id1=get(kernel.bindings,(1,t1),nothing); id2=get(kernel.bindings,(1,t2),nothing)
    id1===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE curve with tag $t1"))
    id2===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE curve with tag $t2"))
    g1=shapes[id1].geometry.geometry; g2=shapes[id2].geometry.geometry
    g1.occ in (:line,:circle) || throw(ArgumentError(
        "$caller: only line and circle edges are supported — got $(g1.occ)"))
    g2.occ in (:line,:circle) || throw(ArgumentError(
        "$caller: only line and circle edges are supported — got $(g2.occ)"))
    # Native-backed records still carry center/frame/t-bounds on the outer
    # record — trimming narrows both levels in `_occ_2df_apply!`.
    corners=_occ_2df_corners(shapes,id1,id2,tol)
    isempty(corners) && throw(ArgumentError(
        "$caller: Could not create temporary wire from given edges. " *
        "Are edges connected?"))
    # Upstream resolves `pointTag` only after the temporary wire/face gates.
    pt=if point_tag<0
        shapes[shapes[id1].children[1]].geometry.point   # FirstVertex(ed1)
    else
        ptv=_model_occ_requested_tag(point_tag,caller)
        vid=get(kernel.bindings,(0,ptv),nothing)
        vid===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE point with tag $ptv"))
        shapes[vid].geometry.point
    end
    # Edges may share several vertices (closed/looped pairs) — upstream's
    # `Result(point,…,-1)` picks the nearest valid corner, so try each
    # coincident corner in distance order until the solve succeeds.
    sort!(corners,by=c->_sqdist(c.point,pt))
    sol=nothing; corner=nothing; lasterr=nothing
    for c in corners
        V=c.point
        u1=_occ_2df_inward(g1,c.end1,caller)
        u2=_occ_2df_inward(g2,c.end2,caller)
        cr=_occ_cross(u1,u2)
        _occ_modulus(cr)<=1e3*eps(Float64) && (lasterr=
            "$caller: edges do not form a corner at the shared vertex";
            continue)
        n=_occ_2df_norm(cr)
        # The corner plane must actually contain both curves — the temp-face
        # build upstream fails otherwise.
        (_occ_2df_planar(g1,V,n,tol) && _occ_2df_planar(g2,V,n,tol)) ||
            (lasterr="$caller: Could not create temporary face from given " *
                "edges. Are edges planar?"; continue)
        try
            sol=radius===nothing ?
                _occ_2df_solve_chamfer(g1,g2,V,c,dist1,dist2,tol,caller) :
                _occ_2df_solve_fillet(g1,g2,u1,u2,n,V,c,Float64(radius),
                    reverse ? -1.0 : 1.0,tol,caller)
            corner=c
            break
        catch err
            err isa ArgumentError || rethrow()
            lasterr=err.msg
        end
    end
    sol===nothing && throw(ArgumentError(lasterr===nothing ?
        "$caller: Could not compute the corner operation" : lasterr))
    _occ_2df_apply!(kernel,id1,id2,g1,g2,corner,allocated,sol,caller)
    # The trimmed edges keep their tags with new geometry and children —
    # upstream rebinds them, so drop any previously synchronized model
    # entries to force re-import (untouched keys are no-ops).
    for t in (t1,t2)
        delete!(m.curves,t); delete!(m.curve_types,t)
        delete!(m.curve_geometry,t)
    end
    return allocated
end

# Both endpoints and (for circles) the mid-span must sit in the corner plane.
function _occ_2df_planar(g,V,n,tol)
    ts=g.occ===:line ? (g.t0,g.t1) : (g.t0,(g.t0+g.t1)/2,g.t1)
    for t in ts
        abs(_dot3(_arc_sub(_occ_2df_eval(g,t),V),n))>10tol && return false
    end
    return true
end

function _model_occ_fillet_2d!(m::GeoModel,e1t,e2t,radius;
                              tag=-1,point_tag=-1,reverse=false)
    radius isa Real || throw(ArgumentError("OCC.fillet_2d: radius must be real"))
    _model_occ_corner_op!(m,e1t,e2t;tag,point_tag,reverse,
        radius=Float64(radius),dist1=0.0,dist2=0.0,caller="OCC.fillet_2d")
end

function _model_occ_chamfer_2d!(m::GeoModel,e1t,e2t,d1,d2;
                                tag=-1)
    _model_occ_corner_op!(m,e1t,e2t;tag,point_tag=-1,reverse=false,
        radius=nothing,dist1=d1,dist2=d2,caller="OCC.chamfer_2d")
end

# --- occ.offsetCurve (upstream OCC_Internals::offsetCurve) -------------------
# BRepOffsetAPI_MakeOffset(wire, GeomAbs_Arc): every member edge of a planar
# wire is offset by `offset` along its left-of-travel normal n̂×t̂ and
# consecutive offsets are re-joined — trimmed to a mutual intersection when
# the offset curves meet inside their spans, otherwise bridged by a
# vertex-centered arc of radius |offset| (the offset image of each original
# endpoint is at distance |o| from its vertex, so that arc is tangent to both
# neighbors). Upstream `_multiBind` publishes the result's edges and vertices
# recursively — the wire itself stays unbound — and returns the (1,tag) pairs.

# Map a member edge to its traversal-directed record — reversal folds the
# original end-to-start evaluation into a same-span record running forward.
# Reversed circles flip n and carry X̂=radial(t0+t1), Ŷ=X sin(t0+t1)−Y cos(t0+t1)
# so eval'(s)=eval(t0+t1−s); the offset then uses the directed frame directly.
function _occ_off_directed(g,rev)
    rev || return g
    g.occ===:line && return (occ=:line,origin=_occ_2df_eval(g,g.t1),
        dir=_occ3_scale(_occ_dir(g.dir),-1.0),t0=0.0,t1=g.t1-g.t0)
    # Native-backed circles evaluate through the nested frame — reverse that.
    f=_occ_2df_frame(g)
    a=g.t0+g.t1; ca=cos(a); sa=sin(a)
    Xp=_occ3_add(_occ3_scale(f.X,ca),_occ3_scale(f.Y,sa))
    Yp=_occ3_add(_occ3_scale(f.X,sa),_occ3_scale(f.Y,-ca))
    np=_occ3_scale(f.n,-1.0)
    ng=(occ=:circle,center=f.center,n=np,X=Xp,Y=Yp,r=f.r,t0=g.t0,t1=g.t1)
    hasproperty(g,:native_geometry) && (
        ng=merge(ng,(native_geometry=merge(g.native_geometry,
            (n=np,X=Xp,Y=Yp)),)))
    return ng
end

# Shift a directed record by o·(n̂×t̂): translate a line along its left normal;
# shrink a circle traveling CCW about n̂ (interior on the left), grow CW ones.
function _occ_off_shift(g,o,n,caller)
    if g.occ===:line
        d=_occ_dir(g.dir)
        L=_occ_2df_norm(_occ_cross(n,d))
        return (occ=:line,origin=_occ3_add(g.origin,_occ3_scale(L,o)),dir=d,
            t0=g.t0,t1=g.t1)
    end
    # Emit through the frame evaluators actually use — a native-backed record
    # can carry an outer frame that differs from its nested `native_geometry`.
    f=_occ_2df_frame(g)
    s=_dot3(f.n,n)>=0.0 ? 1.0 : -1.0
    r=f.r-o*s
    r>_OCC_CONFUSION || throw(ArgumentError(
        "$caller: Could not compute offset curve (offset swallows a circle)"))
    return (occ=:circle,center=f.center,n=f.n,X=f.X,Y=f.Y,r=r,
        t0=g.t0,t1=g.t1)
end

# All intersections of two directed offset records, as (point,paramA,paramB).
function _occ_off_intersections(gA,gB,n)
    out=Tuple{NTuple{3,Float64},Float64,Float64}[]
    if gA.occ===:line && gB.occ===:line
        dA=_occ_dir(gA.dir); dB=_occ_dir(gB.dir)
        cr=_occ_cross(dA,dB); m2=_dot3(cr,cr)
        m2<=1e-14 && return Tuple{NTuple{3,Float64},Float64,Float64}[]
        w=_arc_sub(gB.origin,gA.origin)
        tA=_dot3(_occ_cross(w,dB),cr)/m2
        p=_occ3_add(gA.origin,_occ3_scale(dA,tA))
        return [(p,_occ_2df_param(gA,p),_occ_2df_param(gB,p))]
    elseif gA.occ===:line && gB.occ===:circle
        for C in _occ_2df_line_circle(gA.origin,_occ_dir(gA.dir),1.0,gB.center,gB.r)
            push!(out,(C,_occ_2df_param(gA,C),_occ_2df_param(gB,C)))
        end
    elseif gA.occ===:circle && gB.occ===:line
        for C in _occ_2df_line_circle(gB.origin,_occ_dir(gB.dir),1.0,gA.center,gA.r)
            push!(out,(C,_occ_2df_param(gA,C),_occ_2df_param(gB,C)))
        end
    else
        for C in _occ_2df_circle_circle(gA.center,gA.r,gB.center,gB.r,n)
            push!(out,(C,_occ_2df_param(gA,C),_occ_2df_param(gB,C)))
        end
    end
    return out
end

# Pick the join between directed offset records gA (ending near V) and gB
# (starting near V): contiguous pass-through, interior-intersection trim, or
# a vertex-centered arc of radius |o| (each offset image of V is |o| away, so
# that arc is tangent to both neighbors — GeomAbs_Arc).
function _occ_off_join(gA,gB,V,n,o,tol,caller)
    Aprev=_occ_2df_eval(gA,gA.t1); Bnext=_occ_2df_eval(gB,gB.t0)
    cands=Tuple{NTuple{3,Float64},Float64,Float64}[]
    for (p,tA,tB) in _occ_off_intersections(gA,gB,n)
        # A valid cut trims each span on the corner side: inside A before its
        # end, inside B after its start. Crossings on the far side of both
        # spans (past A's start AND past B's end, or into either extension)
        # belong to the infinite curves, not the offset wire.
        tA>gA.t1-tol && continue          # past A's end — extension
        tB<gB.t0+tol && continue          # before B's start — extension
        if tA<=gA.t0+tol || tB>=gB.t1-tol
            # Past the far ends of both — an irrelevant crossing.
            tA<=gA.t0+tol && tB>=gB.t1-tol && continue
            # …but a crossing inside one span and past the other eats a
            # member: upstream's MakeOffset absorbs it; report loudly.
            throw(ArgumentError(
                "$caller: offset consumes an edge (self-intersecting " *
                "result not supported)"))
        end
        push!(cands,(p,tA,tB))
    end
    if !isempty(cands)
        # The cut nearest the vertex wins.
        i=argmin(_sqdist(c[1],V)+abs(gA.t1-c[2])+abs(c[3]-gB.t0)
            for c in cands)
        p,tA,tB=cands[i]
        return (tA=tA,tB=tB,arc=nothing,pa=Aprev,pb=Bnext)
    end
    # Tangent-continuous members: the offsets already meet — shared vertex.
    _sqdist(Aprev,Bnext)<=(10tol)^2 && return
        (tA=gA.t1,tB=gB.t0,arc=nothing,pa=Aprev,pb=Bnext)
    # Arc join centered at the corner vertex.
    r=abs(o)
    r<=tol && throw(ArgumentError("$caller: degenerate offset join"))
    ρ1=_occ3_scale(_arc_sub(Aprev,V),1.0/r); ρ2=_occ3_scale(_arc_sub(Bnext,V),1.0/r)
    uA=_occ_2df_inward(gA,true,caller); uB=_occ_2df_inward(gB,false,caller)
    nc=_occ_2df_norm(_occ_cross(uA,ρ1))
    tend=_occ_cross(nc,ρ2)
    _dot3(tend,uB)>0 || throw(ArgumentError(
        "$caller: Could not compute offset curve"))
    span=atan(_dot3(_occ_cross(ρ1,ρ2),nc),_dot3(ρ1,ρ2))
    span<=0 && (span+=_OCC_TWO_PI)
    return (tA=gA.t1,tB=gB.t0,
        arc=(occ=:circle,center=V,n=nc,X=ρ1,Y=_occ_cross(nc,ρ1),r=r,
            t0=0.0,t1=span),pa=Aprev,pb=Bnext)
end

# Offset a planar wire (upstream OCC_Internals::offsetCurve). Every fallible
# step precedes the first mutation; the result keeps the wire's traversal
# order and publishes each offset edge (+vertices) under fresh tags.
function _model_occ_offset_wire!(m::GeoModel,loop_tag,offset;caller)
    kernel=m.occ_kernel; shapes=kernel.shapes; tol=_OCC_CONFUSION
    offset isa Real && isfinite(offset) || throw(ArgumentError(
        "$caller: offset must be a finite number"))
    o=Float64(offset)
    lt=_model_occ_requested_tag(loop_tag,caller)
    wid=get(kernel.bindings,(-1,lt),nothing)
    wid===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE curve loop with tag $lt"))
    uses=shapes[wid].geometry.uses
    isempty(uses) && throw(ArgumentError("$caller: empty curve loop"))
    gs=[_occ_off_directed(shapes[id].geometry.geometry,rev) for (id,rev) in uses]
    for g in gs
        g.occ in (:line,:circle) || throw(ArgumentError(
            "$caller: only line and circle edges are supported — got $(g.occ)"))
    end
    N=length(gs)
    V0=shapes[shapes[uses[1][1]].children[uses[1][2] ? 2 : 1]].geometry.point
    # Wire plane normal via Newell accumulation over member sample points —
    # consistent for CCW order at convex AND reflex corners (a corner-tangent
    # cross flips sign at reflex corners).
    nacc=(0.0,0.0,0.0)
    for g in gs
        ts=g.occ===:line ? (g.t0,g.t1) : (g.t0,(g.t0+g.t1)/2,g.t1)
        P=[_occ_2df_eval(g,t) for t in ts]
        for k in 1:length(P)-1
            nacc=_occ3_add(nacc,_occ_cross(P[k],P[k+1]))
        end
    end
    if _occ_modulus(nacc)<=1e3*eps(Float64)
        # Degenerate accumulation — smooth wire (closed circle, tangent chain).
        i=findfirst(g->g.occ===:circle,gs)
        i===nothing && throw(ArgumentError(
            "$caller: Could not compute offset curve " *
            "(wire plane is degenerate)"))
        n=_occ_dir(gs[i].n)
    else
        n=_occ_2df_norm(nacc)
    end
    for g in gs
        _occ_2df_planar(g,V0,n,tol) || throw(ArgumentError(
            "$caller: Could not compute offset curve (wire is not planar)"))
    end
    h1=shapes[uses[1][1]].children[uses[1][2] ? 2 : 1]
    tN=shapes[uses[N][1]].children[uses[N][2] ? 1 : 2]
    closed=_model_occ_vertices_coincide(
        shapes[h1].geometry.point,shapes[tN].geometry.point)
    offs=[_occ_off_shift(g,o,n,caller) for g in gs]
    joins=Vector{Any}(undef,N)
    npairs=closed ? N : N-1
    for i in 1:npairs
        j=mod1(i+1,N)
        ai=uses[i][1]; aj=uses[j][1]
        Vi=shapes[ai].children[uses[i][2] ? 1 : 2]
        Vj=shapes[aj].children[uses[j][2] ? 2 : 1]
        Vi==Vj || _model_occ_vertices_coincide(
            shapes[Vi].geometry.point,shapes[Vj].geometry.point) ||
            throw(ArgumentError("$caller: curve loop edges are not connected"))
        joins[i]=_occ_off_join(offs[i],offs[j],
            shapes[Vi].geometry.point,n,o,tol,caller)
    end
    # Trimmed member spans — join i−1 trims the start, join i the end.
    t0s=Vector{Float64}(undef,N); t1s=Vector{Float64}(undef,N)
    emitted=falses(N)
    for i in 1:N
        jp=closed ? joins[mod1(i-1,N)] : (i>1 ? joins[i-1] : nothing)
        jn=i<=npairs ? joins[i] : nothing
        t0s[i]=jp===nothing ? offs[i].t0 : jp.tB
        t1s[i]=jn===nothing ? offs[i].t1 : jn.tA
        t1s[i]-t0s[i]<=tol && throw(ArgumentError(
            "$caller: offset consumes an edge (self-intersecting result " *
            "not supported)"))
        emitted[i]=true
    end
    # Items in traversal order: member i's trimmed edge, then join i's arc.
    items=Tuple{Any,Symbol,NTuple{3,Float64},NTuple{3,Float64}}[]
    for i in 1:N
        emitted[i] || continue
        rec=merge(offs[i],(t0=t0s[i],t1=t1s[i]))
        push!(items,(rec,offs[i].occ,_occ_2df_eval(offs[i],t0s[i]),
            _occ_2df_eval(offs[i],t1s[i])))
        jn=i<=npairs ? joins[i] : nothing
        if jn!==nothing && jn.arc!==nothing && emitted[mod1(i+1,N)]
            push!(items,(jn.arc,:circle,jn.pa,jn.pb))
        end
    end
    isempty(items) && throw(ArgumentError(
        "$caller: Could not compute offset curve"))
    # Non-adjacent result pieces must not cross — upstream resolves the
    # arising loops; report loudly.
    for a in 1:length(items), b in a+1:length(items)
        b==a+1 && continue
        closed && a==1 && b==length(items) && continue
        for (_,tA,tB) in _occ_off_intersections(items[a][1],items[b][1],n)
            gA=items[a][1]; gB=items[b][1]
            gA.t0+tol<tA<gA.t1-tol && gB.t0+tol<tB<gB.t1-tol &&
                throw(ArgumentError("$caller: offset result self-intersects"))
        end
    end
    # Merge coincident consecutive endpoints into shared vertices.
    verts=NTuple{3,Float64}[]; spans=Tuple{Int,Int}[]
    for (_,_,ps,pe) in items
        a=if isempty(verts) || _sqdist(verts[end],ps)>(10tol)^2
            push!(verts,ps); length(verts)
        else
            length(verts)
        end
        push!(verts,pe); push!(spans,(a,length(verts)))
    end
    if closed
        _sqdist(verts[end],verts[1])<=(10tol)^2 || throw(ArgumentError(
            "$caller: offset result is not closed"))
        pop!(verts); spans[end]=(spans[end][1],1)
    end
    nv=length(verts); ne=length(items)
    kernel.maxima[3]+nv<=typemax(Int32) || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]+ne<=typemax(Int32) || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    # One contiguous reservation — `reserve_shapes` does not advance
    # `next_shape`, so separate calls would hand back overlapping ids.
    ids=_model_occ_reserve_shapes(kernel,nv+ne,caller)
    vids=ids[1:nv]; eids=ids[nv+1:nv+ne]
    vtag0=kernel.maxima[3]
    for (i,p) in enumerate(verts)
        _model_occ_bind_shape!(kernel,0,vtag0+i,vids[i],
            _ModelOCCShape(0,(point=p,mesh_size=0.0),UInt64[]))
    end
    etag0=kernel.maxima[4]
    out=Tuple{Int,Int}[]
    for (i,(rec,kind,ps,pe)) in enumerate(items)
        a,b=spans[i]
        _model_occ_bind_shape!(kernel,1,etag0+i,eids[i],
            _ModelOCCShape(1,(kind=kind,geometry=rec),
                UInt64[vids[a],vids[b]]))
        push!(out,(1,etag0+i))
    end
    kernel.next_shape=eids[ne]+UInt64(1)
    return out
end
