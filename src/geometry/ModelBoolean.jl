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
    end
    throw(ArgumentError("_brep_curve_eval: degenerate edge has no 3-D curve"))
end

# ── surface evaluation / inversion ──────────────────────────────────────────

_brep_surface_eval(sg::NamedTuple,u::Float64,v::Float64) =
    _occ_surface_point(sg,u,v)

# Geometric (unoriented) unit normal at (u,v).
function _brep_surface_normal(sg::NamedTuple,u::Float64,v::Float64)
    if sg.occ===:plane
        return sg.axis
    elseif sg.occ===:cylinder
        s,c=_gm_sincos(u)
        return _b3unit((fma(c,sg.X[1],s*sg.Y[1]),fma(c,sg.X[2],s*sg.Y[2]),
                        fma(c,sg.X[3],s*sg.Y[3])))
    elseif sg.occ===:sphere
        su,cu=_gm_sincos(u); sv,cv=_gm_sincos(v)
        return _b3unit((fma(cv*cu,sg.X[1],fma(cv*su,sg.Y[1],sv*sg.axis[1])),
                        fma(cv*cu,sg.X[2],fma(cv*su,sg.Y[2],sv*sg.axis[2])),
                        fma(cv*cu,sg.X[3],fma(cv*su,sg.Y[3],sv*sg.axis[3]))))
    elseif sg.occ===:cone
        s,c=_gm_sincos(u)
        slope=(sg.r2-sg.r1)/sg.height
        l=sqrt(1.0+slope*slope)
        radial=(fma(c,sg.X[1],s*sg.Y[1]),fma(c,sg.X[2],s*sg.Y[2]),
                fma(c,sg.X[3],s*sg.Y[3]))
        return ((radial[1]-slope*sg.axis[1])/l,
                (radial[2]-slope*sg.axis[2])/l,
                (radial[3]-slope*sg.axis[3])/l)
    end
    throw(ArgumentError("_brep_surface_normal: unsupported surface $(sg.occ)"))
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
        return (u,_b3dot(d,sg.axis))
    elseif sg.occ===:sphere
        d=_b3(p,sg.center)
        z=_b3dot(d,sg.axis); rr=_b3dot(d,sg.X); ri=_b3dot(d,sg.Y)
        v=atan(z,sqrt(rr*rr+ri*ri))
        u=atan(ri,rr); u<0 && (u+=_OCC_TWO_PI)
        return (u,v)
    end
    throw(ArgumentError("_brep_surface_uv: unsupported surface $(sg.occ)"))
end

_brep_uperiod(sg::NamedTuple) =
    sg.occ in (:cylinder,:sphere,:cone,:torus) ? _OCC_TWO_PI : 0.0

# ── pcurve evaluation ────────────────────────────────────────────────────────

@inline _brep_pc_eval(pc,t) =
    hasproperty(pc,:lin2d) ?
        (fma(t,pc.lin2d[2][1],pc.lin2d[1][1]),
         fma(t,pc.lin2d[2][2],pc.lin2d[1][2])) :
        _occ_circ2d_eval(pc.circ2d[1],pc.circ2d[2],pc.circ2d[3],t)

# Tangent vector of a pcurve at t (not normalized).
@inline function _brep_pc_tangent(pc,t)
    if hasproperty(pc,:lin2d)
        return pc.lin2d[2]
    end
    c,x,r=pc.circ2d;s,co=_gm_sincos(t)
    return (r*(-s*x[1]-co*x[2]),r*(-s*x[2]+co*x[1]))
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
    elseif rec.occ===:line
        return (kind=:line,o=rec.origin,d=rec.dir),rec.t0,rec.t1
    elseif rec.occ===:circle
        return (kind=:circle,c=rec.center,n=rec.n,x=rec.X,y=rec.Y,r=rec.r),
               rec.t0,rec.t1
    elseif rec.occ===:degenerate
        return (kind=:degenerate,xyz=m.points[a]),rec.t0,rec.t1
    end
    throw(ArgumentError("_brep_edge_curve: unsupported curve record"))
end

# Synthesize a :plane record for a native plane surface so every operand face
# carries an occ-style record.
function _brep_plane_record(m::GeoModel,stag::Int,caller::AbstractString)
    frame=_model_plane_frame(m,stag,caller)
    return (occ=:plane,center=frame.origin,axis=frame.normal,
            X=frame.first_direction,Y=frame.second_direction,
            reversed=false,pcurves=Dict{Int,NamedTuple}())
end

# The pcurve of an edge on a synthesized plane record: project the curve's
# endpoints over its own parameter range.
function _brep_synth_pcurve(sg,curve::NamedTuple,t0,t1)
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
                        ent=get(sg.pcurves,ctag,nothing)
                        ent===nothing && throw(ArgumentError(
                            "$caller: Surface[$stag] lacks a pcurve for " *
                            "Curve[$ctag]"))
                        pc=(sgn<0 && ent.rev!==nothing) ? ent.rev : ent.fwd
                        pc===nothing && throw(ArgumentError(
                            "$caller: Surface[$stag] lacks a pcurve for " *
                            "Curve[$ctag]"))
                    else
                        e=edges[eidx]
                        pc=_brep_synth_pcurve(sg,e.curve,e.t0,e.t1)
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

# ── 2-D wire containment (even-odd parity) ──────────────────────────────────

# Rightward crossings of a horizontal ray from (pu,pv) against one pcurve over
# the edge's parameter range. Half-open convention: a vertex-level crossing
# counts only when the edge enters the half-plane from below.
function _brep_pc_ray_crossings(pc,t0,t1,pu,pv,tol)
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
    else
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
    end
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

_brep_pc_circ_on_plane(sg,c,x,r) =
    (circ2d=(_brep_plane_uv(sg,c),(_b3dot(x,sg.X),_b3dot(x,sg.Y)),r),)

# circle at constant height on a cylinder/cone-like uv domain (u=angle,v=z).
function _brep_pc_circ_on_cylinder(sg,center,n,x,r)
    sgn=_b3dot(n,sg.axis)>=0 ? 1.0 : -1.0
    ca=_b3dot(x,sg.X);cb=_b3dot(x,sg.Y)
    u0=atan(cb,ca);u0<0 && (u0+=_OCC_TWO_PI)
    sgn<0 && (u0=_OCC_TWO_PI-u0)
    v=_b3dot(_b3(center,sg.center),sg.axis)
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

# line (parallel to axis) at angle u0 on a cylinder/cone domain.
function _brep_pc_line_on_cylinder(sg,o,d)
    p0=_brep_surface_uv(sg,o)
    sgn=_b3dot(d,sg.axis)>=0 ? 1.0 : -1.0
    return (lin2d=(p0,(0.0,sgn)),)
end

function _brep_face_sections(fa,fb,caller)
    ka,kb=fa.sg.occ,fb.sg.occ
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
        pc_pl=_brep_pc_circ_on_plane(pl,c,x,cy.radius)
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
        foot=_b3(axpt,_b3mul(pl.axis,dist))
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
    throw(ArgumentError(
        "$caller: oblique plane×cylinder sections (ellipses) are unsupported"))
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
    curve=(kind=:circle,c=c,n=n,x=x,y=y,r=r)
    pc_pl=NamedTuple[_brep_pc_circ_on_plane(pl,c,x,r)]
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
        pc_pl=_brep_pc_circ_on_plane(pl,c,x,rv)
        pc_cn=_brep_pc_circ_on_cylinder(cn,c,n,x,rv)
        return [(curve=curve,closed=true,
                 pcA=plane_is_a ? pc_pl : pc_cn,
                 pcB=plane_is_a ? pc_cn : pc_pl)]
    end
    throw(ArgumentError(
        "$caller: non-axial plane×cone sections (conics) are unsupported"))
end

function _brep_sec_sphere_sphere(sa,sb,caller)
    d=_b3(sb.center,sa.center);dist=_b3norm(d)
    dist<_BREP_TOL && return NamedTuple[]
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
    # supported when the sphere center lies on the cylinder axis → circles
    d=_b3(sp.center,cy.center)
    axial=_b3dot(d,cy.axis)
    radial=_b3(d,_b3mul(cy.axis,axial))
    _b3norm(radial)>_BREP_TOL && throw(ArgumentError(
        "$caller: off-axis cylinder×sphere sections are unsupported"))
    dz2=sp.radius^2-cy.radius^2
    # the sphere reaches the lateral wall only when r > R — a smaller
    # on-axis sphere never meets it; r ≈ R is tangent contact
    dz2<-_BREP_TOL && return NamedTuple[]
    dz2<=_BREP_TOL && throw(ArgumentError(
        "$caller: tangent cylinder×sphere contact is unsupported"))
    dz=sqrt(dz2)
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
        _b3norm(proj)<_BREP_TOL && throw(ArgumentError(
            "$caller: coaxial cylinders (same-domain) are not yet supported"))
        return NamedTuple[]
    end
    throw(ArgumentError(
        "$caller: non-coaxial cylinder×cylinder sections are unsupported"))
end

# ── 2-D pcurve × pcurve intersections ───────────────────────────────────────
#
# Intersect a section pcurve (infinite line or full circle) with a boundary
# edge's pcurve (over the edge parameter range [et0,et1]); return
# section-parameter values.

function _brep_pc_intersect(sec_pc,edge_pc,et0,et1,tol)
    if hasproperty(sec_pc,:lin2d) && hasproperty(edge_pc,:lin2d)
        return _brep_ix_lin_lin(sec_pc.lin2d,edge_pc.lin2d,et0,et1,tol)
    elseif hasproperty(sec_pc,:lin2d) && hasproperty(edge_pc,:circ2d)
        return _brep_ix_lin_circ(sec_pc.lin2d,edge_pc.circ2d,et0,et1,tol)
    elseif hasproperty(sec_pc,:circ2d) && hasproperty(edge_pc,:lin2d)
        return _brep_ix_circ_lin(sec_pc.circ2d,edge_pc.lin2d,et0,et1,tol)
    else
        return _brep_ix_circ_circ(sec_pc.circ2d,edge_pc.circ2d,et0,et1,tol)
    end
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

function _brep_ix_lin_lin(sl,el,et0,et1,tol)
    so,sd=sl;eo,ed=el
    det=sd[1]*ed[2]-sd[2]*ed[1]
    abs(det)<tol && return Float64[]
    w=(eo[1]-so[1],eo[2]-so[2])
    t=(w[1]*ed[2]-w[2]*ed[1])/det
    s=(w[1]*sd[2]-w[2]*sd[1])/det
    lo,hi=minmax(et0,et1)
    (s<lo-tol || s>hi+tol) && return Float64[]
    return [t]
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

function _brep_trim_section(solid,face,pc,closed,tol)
    uper=_brep_uperiod(face.sg)
    ulo=_brep_face_ubase(solid.edges,face)
    if uper>0
        # fold the section pcurve into the face's canonical u band so its
        # crossings with the wires compare in the same coordinates
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
            append!(params,_brep_pc_intersect(pc,epc,e.t0,e.t1,tol))
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
    if hasproperty(pc,:circ2d)
        isempty(uniq) && begin
            uv0=_brep_wrap_u(_brep_pc_eval(pc,0.0),ulo,uper)
            (_brep_uv_in_face(solid.edges,face,uv0,tol) ||
             _brep_uv_on_wire(solid.edges,face,pc,uv0,tol)) &&
                push!(out,(0.0,_OCC_TWO_PI))
            return out
        end
        u0=mod(uniq[1],_OCC_TWO_PI)
        while !isempty(uniq) && uniq[end]>=u0+_OCC_TWO_PI-tol
            pop!(uniq)
        end
        pushfirst!(uniq,u0)
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
    if curve.kind===:line
        append!(params,_brep_line_surface(sg,curve))
    elseif curve.kind===:circle
        append!(params,_brep_circle_surface(sg,curve,callerless=true))
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
    return out
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
    lo,hi=minmax(t0,t1)
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return minmax(o[1]+d[1]*t0,o[1]+d[1]*t1)
    end
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
end

# Split a pcurve segment at u=u_seam0+k·period crossings (periodic cut-domain
# seam). Returns (t_ranges) covering [t0,t1] with no interior seam crossing.
function _brep_pc_split_at_useam(pc,t0,t1,u_seam0,period,tol)
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
    else
        c,x,r=pc.circ2d
        # u(t)=c1+r·cos(t+α)=u_seam0+k·period
        α=atan(x[2],x[1])
        kmin=floor(Int,(c[1]-r-u_seam0)/period)-1
        kmax=ceil(Int,(c[1]+r-u_seam0)/period)+1
        for k in kmin:kmax
            du=(u_seam0+k*period-c[1])/r
            abs(du)<=1.0 || continue
            β=acos(clamp(du,-1.0,1.0))
            for γ in (β,-β)
                tt=mod(γ-α,_OCC_TWO_PI)
                for j in -2:2
                    tk=tt+j*_OCC_TWO_PI
                    (tk>t0+tol && tk<t1-tol) && push!(cuts,tk)
                end
            end
        end
    end
    sort!(cuts)
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
    if hasproperty(pc,:lin2d)
        u0,v0=_brep_pc_eval(pc,ta)
        u1,v1=_brep_pc_eval(pc,tb)
        return (u0*v1-u1*v0)/2
    end
    c,x,r=pc.circ2d
    # u·v'−v·u' = r² + r·(c1·y'(t) − c2·x'(t))·(in x,y coords):
    # rel(t)=r(cos t·x + sin t·y); rel' = r(−sin·x + cos·y)
    # u v'−v u' = r² + c·rel' − rel·c' where c'=0 → r² + c1·rv' − c2·ru'
    # integrate: r²·Δt + c1·r·(sin x·(−cos t) + cos x·(sin t))... computed:
    s0,co0=_gm_sincos(ta);s1,co1=_gm_sincos(tb)
    # ∫(u v' − v u')dt = r²·(tb−ta) + r·[c1·(x2(cos tb−cos ta)+y2(sin tb−sin ta))
    #                      − c2·(x1(cos tb−cos ta)+y1(sin tb−sin ta))]
    dc=co1-co0;ds=s1-s0
    return (r*r*(tb-ta) +
            r*(c[1]*(x[2]*dc - x[1]*ds) - c[2]*(x[1]*dc + x[2]*ds)))/2
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
            _brep_cycle_inside(arcs,cycle,cand,tol) && return cand
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
        abs(a)<tol && throw(ErrorException(
            "$caller: degenerate zero-area wire in Boolean face split"))
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

function _brep_ray_face_hits(solid,face,P,D,tol)
    sg=face.sg
    params=_brep_line_surface(sg,(o=P,d=D))
    count=0;ambiguous=false
    for t in params
        t<=tol && continue
        p=_b3add(P,_b3mul(D,t))
        uv=_brep_surface_uv(sg,p)
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
    end
    c,x,r=pc.circ2d
    dx=uv[1]-c[1];dy=uv[2]-c[2]
    dist=abs(sqrt(dx*dx+dy*dy)-r)
    dist<tol || return false
    θ=atan((dy*x[1]-dx*x[2]),(dx*x[1]+dy*x[2]))
    return _brep_snap_theta(θ,min(t0,t1),max(t0,t1),tol)!==nothing
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
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=((o[1]-shift,o[2]),d),)
    end
    c,x,r=pc.circ2d
    return (circ2d=((c[1]-shift,c[2]),x,r),)
end

# every parameter in (t0,t1) where the pcurve's u crosses a cut-domain seam
# (u = ulo + k·period) for U-periodic surfaces
function _brep_pc_split_all_useams(pc,t0,t1,period,ulo,tol)
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
    else
        c,x,r=pc.circ2d
        α=atan(x[2],x[1])
        kmin=floor(Int,(c[1]-r-ulo)/period)-1
        kmax=ceil(Int,(c[1]+r-ulo)/period)+1
        for k in kmin:kmax
            v=(ulo+k*period-c[1])/r
            abs(v)<=1.0 || continue
            β=acos(clamp(v,-1.0,1.0))
            for γ in (β,-β)
                tt=mod(γ-α,_OCC_TWO_PI)
                for j in -2:2
                    tk=tt+j*_OCC_TWO_PI
                    (tk>t0+tol && tk<t1-tol) && push!(cuts,tk)
                end
            end
        end
    end
    sort!(cuts)
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
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        n2=d[1]*d[1]+d[2]*d[2]
        n2<tol*tol && return nothing
        s=((uv[1]-o[1])*d[1]+(uv[2]-o[2])*d[2])/n2
        e=_brep_pc_eval(pc,s)
        (abs(e[1]-uv[1])<10tol && abs(e[2]-uv[2])<10tol) || return nothing
        (s>=t0-tol && s<=t1+tol) || return nothing
        return s
    end
    c,x,r=pc.circ2d
    dx,dy=uv[1]-c[1],uv[2]-c[2]
    cx=dx*x[1]+dy*x[2];cy=dy*x[1]-dx*x[2]
    abs(sqrt(cx*cx+cy*cy)-r)>10tol*max(1.0,r) && return nothing
    return _brep_snap_theta(atan(cy,cx),t0,t1,tol)
end

# Build the face's uv arrangement (boundary sub-edges + section pieces, cut at
# the periodic seam) and extract kept-region wires by leftmost-turn traversal.
# Returns (regions, arcs, cycles) — each region is a Vector of cycles (outer
# first); a cycle is a Vector of (arc_index, forward).
function _brep_face_regions(solid,f,side,edge_splits,secs,extras,tol,caller)
    uper=_brep_uperiod(f.sg)
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
        extra=get(edge_splits,eidx,Float64[])
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
            push!(raw,(pc=pc,t0=ta,t1=tb,key=(:e,side,eidx,k),crv=e.curve))
        end
    end
    for sf in kept_secs
        push!(raw,(pc=sf.pc,t0=sf.span[1],t1=sf.span[2],
                   key=(:s,sf.se),crv=sf.crv))
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
                push!(arcs,(pc=pc2,t0=ta,t1=tb,key=ar.key,crv=ar.crv,
                            vt0=va,vt1=vb,neg=neg,a=0,b=0))
            end
        else
            va,vb=neg ? (-ar.t1,-ar.t0) : (ar.t0,ar.t1)
            push!(arcs,(pc=ar.pc,t0=ar.t0,t1=ar.t1,key=ar.key,crv=ar.crv,
                        vt0=va,vt1=vb,neg=neg,a=0,b=0))
        end
    end
    isempty(arcs) && return NamedTuple[],arcs,Vector{Tuple{Int,Bool}}[]
    pts=Tuple{Float64,Float64}[]
    for ar in arcs
        push!(pts,_brep_pc_eval(ar.pc,ar.t0))
        push!(pts,_brep_pc_eval(ar.pc,ar.t1))
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
    cycles=_brep_uv_arrange(arcs,length(nodes),tol,caller)
    regions,_=_brep_group_cycles(arcs,cycles,solid.edges,f,tol,caller)
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
function _brep_intersect_spans(spansA,spansB,closed,tol,caller)
    if !closed
        ovl=Tuple{Float64,Float64}[]
        for (a0,a1) in spansA,(b0,b1) in spansB
            lo=max(a0,b0);hi=min(a1,b1)
            lo<hi-tol || continue
            push!(ovl,(lo,hi))
        end
    else
        cuts=Float64[0.0]
        for s in Iterators.flatten((spansA,spansB)), t in s
            isfinite(t) || continue
            c=mod(t,_OCC_TWO_PI)
            any(q->abs(q-c)<tol,cuts) || push!(cuts,c)
        end
        sort!(cuts)
        incell(spans,m)=any(
            ((a,b),)->begin
                for k in -1:1
                    mk=m+k*_OCC_TWO_PI
                    (mk>a+tol && mk<b-tol) && return true
                end
                false
            end, spans)
        ovl=Tuple{Float64,Float64}[]
        n=length(cuts)
        for i in 1:n
            a=cuts[i];b=i<n ? cuts[i+1] : cuts[1]+_OCC_TWO_PI
            m=(a+b)/2
            (incell(spansA,m) && incell(spansB,m)) || continue
            push!(ovl,(a,b))
        end
        # merge across the wrap point when both end cells are kept
        while length(ovl)>1 && ovl[end][2]>=ovl[1][1]+_OCC_TWO_PI-tol
            ovl[1]=(ovl[end][1]-_OCC_TWO_PI,ovl[1][2])
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
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=(_brep_map_uv(mp,o),_brep_map_dir(mp,d)),),false
    end
    c,x,r=pc.circ2d
    neg=(mp.c1[1]*mp.c2[2]-mp.c1[2]*mp.c2[1])<0
    return (circ2d=(_brep_map_uv(mp,c),_brep_map_dir(mp,x),r),),neg
end

# pcurve with its parameter reversed (t → −t): used to store an imprinted
# arc's pcurve under the operand edge's own parametrization. Reversible for
# lines; a reversed circle is not a `circ2d` and cannot be represented.
function _brep_flip_pc(pc,caller)
    hasproperty(pc,:lin2d) &&
        return (lin2d=(pc.lin2d[1],(-pc.lin2d[2][1],-pc.lin2d[2][2])),)
    throw(ArgumentError(
        "$caller: coincident-face imprint of a circle under a mirrored " *
        "frame has no representable pcurve"))
end

# same carrier test: collinear lines or coincident circles
function _brep_pc_same_carrier(a,b,tol)
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
    end
    return false
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
    end
    return false
end

# parameter of a 3-D point on a result curve, lifted near tref (circles)
function _brep_curve_param(crv,p,tref)
    if crv.kind===:line
        return _b3dot(_b3(p,crv.o),crv.d)
    elseif crv.kind===:circle
        rel=_b3(p,crv.c)
        θ=atan(_b3dot(rel,crv.y),_b3dot(rel,crv.x))
        return tref+mod(θ-tref+π,_OCC_TWO_PI)-π
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
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=((fma(β,d[1],o[1]),fma(β,d[2],o[2])),
                       (α*d[1],α*d[2])),),false
    end
    abs(abs(α)-1)>100tol && throw(ArgumentError(
        "$caller: coincident circle section with a non-angular parameter " *
        "map has no representable pcurve"))
    c,x,r=pc.circ2d
    cb,sb=_gm_sincos(β)
    # frame rotated by β: eval(t) = pc(t+β); for α<0 evaluated at s=−t it
    # yields pc(−t+β) — exactly the reversed-parameter map
    xr=(fma(cb,x[1],-sb*x[2]),fma(cb,x[2],sb*x[1]))
    return (circ2d=(c,xr,r),),α<0
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
        c0,s0=_gm_sincos(sgn*t0s)
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
        if hasproperty(mpc,:circ2d)
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
        nlifts=hasproperty(spc,:circ2d) ? (-1,0,1) : (0,)
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
                        key=(:e,wside,es.weidx,k),crv=es.crv,
                        neg=es.neg,esec=true))
        end
    end
    return arcs
end

# Build the imprint arcs on face f: O's boundary-edge pieces (split by the
# finalized global edge splits) mapped into f's domain, kept where they lie
# strictly inside it and are not already carried by f's own wires or by a
# section arc. `oside` is O's operand index so the arc keys match O's own
# wire pieces for the same edge piece.
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
                        key=(:e,oside,im.eidx,k),crv=im.crv,neg=im.neg))
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
    for k in (:radius,:radius1,:radius2,:half_angle)
        hasproperty(a.sg,k) || continue
        abs(getproperty(a.sg,k)-getproperty(b.sg,k))>100tol && return false
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
    keymem=Dict{Any,Set{Int}}()
    for mi in comp,w in kept[mi].wires,ar in w
        push!(get!(Set{Int},keymem,ar.key),mi)
    end
    marcs=NamedTuple[]
    test_arcs=NamedTuple[]
    member_cycles=Vector{Vector{Vector{Tuple{Int,Bool}}}}()
    for mi in comp
        mp=maps[mi]
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
                end
                push!(test_arcs,rec)
                push!(cy,(length(test_arcs),ar.fwd))
                if length(keymem[ar.key])>1
                    continue   # internal to the merged domain
                end
                push!(marcs,merge(rec,(a=0,b=0)))
            end
            push!(cycs,cy)
        end
        push!(member_cycles,cycs)
    end
    isempty(marcs) && return [kept[mi] for mi in comp]
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
        tol<abs(areas[i]) || throw(ErrorException(
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
        comps=Dict{Int,Vector{Int}}()
        for mi in mem
            push!(get!(Vector{Int},comps,findr(mi)),mi)
        end
        for comp in values(comps)
            if length(comp)==1
                push!(out,kept[comp[1]])
            else
                append!(out,_brep_merge_union_faces(kept,comp,tol,caller))
            end
        end
    end
    return _brep_merge_edge_pieces(out,tol,caller)
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
        for mp in maps
            canon[mp.k]=(root,mp.α,mp.β)
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
            last.fwd==rec.fwd && abs(exitp(last)-entryp(rec))<tol &&
            _brep_pc_same_carrier(last.pc,rec.pc,tol)
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
                    rec= f ? (key=root,fwd=ar.fwd,pc=pc2,t0=-vhi,t1=-vlo,
                              crv=keyrep[root].crv,vt0=vlo,vt1=vhi,neg=true) :
                             (key=root,fwd=ar.fwd,pc=pc2,t0=vlo,t1=vhi,
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
            spansA=_brep_trim_section(A,fa,s.pcA,s.closed,tol)
            spansB=_brep_trim_section(B,fb,s.pcB,s.closed,tol)
            merged=_brep_intersect_spans(spansA,spansB,s.closed,tol,caller)
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
            isempty(ps) || (es[eidx]=sort!(ps))
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
            p1=_brep_curve_eval(sec.curve,lo)
            p2=_brep_curve_eval(sec.curve,hi)
            wmid=(we.t0+we.t1)/2
            t1=_brep_curve_param(we.curve,p1,wmid)
            t2=_brep_curve_param(we.curve,p2,wmid)
            (t1===nothing || t2===nothing) && throw(ArgumentError(
                "$caller: a coincident section endpoint does not project " *
                "onto its operand edge"))
            wlo,whi=minmax(t1,t2)
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
                    α=(hi-lo)/(t2-t1);β=lo-α*t1
                    _brep_pc_reparam(spc,α,β,tol,caller)
                else
                    _brep_pc_on_surface(tface.sg,we.curve,wlo,whi,tol,caller)
                end
                sslo,sshi=neg ? (-whi,-wlo) : (wlo,whi)
                _brep_imprint_contacts!(L,W,tf,weidx,pc,neg,sslo,sshi,
                                        esplits[ts],esplits[wside],sec_cuts,
                                        sections,ts==1,tol,si)
                # span endpoints landing on the target face's wires split
                # those edges too
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
                # the same 3-D edge can be reached through several section
                # records (each face pair of the 4-way coincidence emits
                # one) — collect the spans under one spec rather than
                # duplicating the carrier
                ex_i=findfirst(es->es.weidx==weidx,esec_imps[ts][tf])
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
    # section sub-edges, cut at both faces' cut-domain seams (a seam crossing
    # is a real result vertex shared by both faces' arcs, so the same
    # sub-edge keys must appear on each side) and at coincident contacts
    sec_edges=NamedTuple[]
    secs_A=[NamedTuple[] for _ in A.faces]
    secs_B=[NamedTuple[] for _ in B.faces]
    for (si,sec) in enumerate(sections)
        esec_skip[si] && continue
        for (lo,hi) in sec.spans
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
                push!(sec_edges,(curve=sec.curve,t0=s0,t1=s1))
                se=length(sec_edges)
                push!(secs_A[sec.fa],(se=se,span=(s0,s1),pc=sec.pcA,
                                      crv=sec.curve))
                push!(secs_B[sec.fb],(se=se,span=(s0,s1),pc=sec.pcB,
                                      crv=sec.curve))
            end
        end
    end
    # section endpoints landing on a wedge split that edge globally
    for (sd,S,secs) in ((1,A,secs_A),(2,B,secs_B))
        es=esplits[sd]
        for (fidx,f) in enumerate(S.faces), sf in secs[fidx]
            uper=_brep_uperiod(f.sg)
            ulo=_brep_face_ubase(S.edges,f)
            for tend in (sf.span[1],sf.span[2])
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
                keep || continue
                wires=[ [ (key=arcs[i].key,fwd=fwd,pc=arcs[i].pc,
                           t0=arcs[i].t0,t1=arcs[i].t1,crv=arcs[i].crv,
                           vt0=arcs[i].vt0,vt1=arcs[i].vt1,neg=arcs[i].neg)
                         for (i,fwd) in cycles[c] ] for c in region ]
                sgn=flip ? -f.orev : f.orev
                push!(kept,(sg=f.sg,sgn=sgn,wires=wires,stag=f.stag))
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
        if hasproperty(mpc,:circ2d)
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
        nlifts=hasproperty(spc,:circ2d) ? (-1,0,1) : (0,)
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
                spansI=_brep_trim_section(Si,fI,s.pcA,s.closed,tol)
                spansJ=_brep_trim_section(Sj,fJ,s.pcB,s.closed,tol)
                merged=_brep_intersect_spans(spansI,spansJ,s.closed,tol,caller)
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
            p1=_brep_curve_eval(sec.curve,lo)
            p2=_brep_curve_eval(sec.curve,hi)
            wmid=(we.t0+we.t1)/2
            t1=_brep_curve_param(we.curve,p1,wmid)
            t2=_brep_curve_param(we.curve,p2,wmid)
            (t1===nothing || t2===nothing) && throw(ArgumentError(
                "$caller: a coincident section endpoint does not project " *
                "onto its operand edge"))
            wlo,whi=minmax(t1,t2)
            push!(get!(Vector{Float64},esplits[w],weidx),wlo,whi)
            touched[w]=true
            for (k2,oe) in zip(owners,owned)
                k2==w && continue
                oer=solids[k2].edges[oe]
                omid=(oer.t0+oer.t1)/2
                l1=_brep_curve_param(oer.curve,p1,omid)
                l2=_brep_curve_param(oer.curve,p2,omid)
                l1===nothing || l2===nothing || (push!(get!(
                    Vector{Float64},esplits[k2],oe),min(l1,l2),max(l1,l2));
                    touched[k2]=true)
            end
            for (tk,tf,spc) in targets
                tface=solids[tk].faces[tf]
                pc,neg = if spc!==nothing
                    α=(hi-lo)/(t2-t1);β=lo-α*t1
                    _brep_pc_reparam(spc,α,β,tol,caller)
                else
                    _brep_pc_on_surface(tface.sg,we.curve,wlo,whi,tol,caller)
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
                ex_i=findfirst(es->es.wside==w && es.weidx==weidx,
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
    # section sub-edges, cut at both faces' seam lifts and contact cuts
    sec_edges=NamedTuple[]
    secs_pf=[[NamedTuple[] for _ in solids[i].faces] for i in 1:N]
    for (si,sec) in enumerate(sections)
        esec_skip[si] && continue
        for (lo,hi) in sec.spans
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
                push!(sec_edges,(curve=sec.curve,t0=s0,t1=s1))
                se=length(sec_edges)
                push!(secs_pf[sec.oi][sec.fa],
                      (se=se,span=(s0,s1),pc=sec.pcA,crv=sec.curve))
                push!(secs_pf[sec.oj][sec.fb],
                      (se=se,span=(s0,s1),pc=sec.pcB,crv=sec.curve))
            end
        end
    end
    # section endpoints landing on a wire edge split that edge globally
    for i in 1:N
        es=esplits[i];S=solids[i]
        for (fidx,f) in enumerate(S.faces), sf in secs_pf[i][fidx]
            uper=_brep_uperiod(f.sg)
            ulo=_brep_face_ubase(S.edges,f)
            for tend in (sf.span[1],sf.span[2])
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
                push!(kept,(sg=f.sg,stag=f.stag,wires=wires,uses=uses,
                            sgn=uses[1][2]))
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
    end
    # endpoints → vertices (deduped, first appearance in wire order); the
    # endpoints of a merged edge are its span extrema — interior piece
    # junctions are not result vertices
    edge_verts=Vector{Tuple{Int,Int}}(undef,length(res_edges))
    for (ei,e) in enumerate(res_edges)
        va=vertex_of(_brep_vertex_eval(e.curve,e.t0))
        vb=vertex_of(_brep_vertex_eval(e.curve,e.t1))
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
            else
                c=_add_occ_circle!(m,ptag[va],ptag[vb],crv.c,crv.n,
                                   crv.x,crv.y,crv.r,e.t0,e.t1)
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
                            push!(ids,ar.fwd ? c : -c)
                            # stored pcurves are edge-parametrized — mirrored
                            # imprints negate their pcurve's parameter
                            pcst=ar.neg ? _brep_flip_pc(ar.pc,caller) : ar.pc
                            ent=get(pcdict,c,nothing)
                            if ar.fwd
                                ent===nothing || ent.fwd===nothing ||
                                    throw(ErrorException(
                                        "$caller: duplicate forward pcurve " *
                                        "on Curve[$c]"))
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
        for w in spec.wires, ar in w
            (ar.key isa Tuple && length(ar.key)>=3 &&
             ar.key[1]===:e && ar.key[2]==i) || (ok=false;break)
        end
        ok || continue
        reuse[fi]=true;emit[fi]=false
        for w in spec.wires, ar in w
            reuse_ctag[ar.key]=pseudo_ctag[i][ar.key[3]]
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
        end
    end
    edge_verts=Vector{Tuple{Int,Int}}(undef,length(res_edges))
    for (ei,e) in enumerate(res_edges)
        va=vertex_of(_brep_vertex_eval(e.curve,e.t0))
        vb=vertex_of(_brep_vertex_eval(e.curve,e.t1))
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
            else
                c=_add_occ_circle!(m,ptag[va],ptag[vb],crv.c,crv.n,
                                   crv.x,crv.y,crv.r,e.t0,e.t1)
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
                    push!(ids,ar.fwd ? c : -c)
                    pcst=ar.neg ? _brep_flip_pc(ar.pc,caller) : ar.pc
                    ent=get(pcdict,c,nothing)
                    if ar.fwd
                        ent===nothing || ent.fwd===nothing ||
                            throw(ErrorException(
                                "$caller: duplicate forward pcurve " *
                                "on Curve[$c]"))
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
                lt=add_curve_loop!(m,ids)
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
