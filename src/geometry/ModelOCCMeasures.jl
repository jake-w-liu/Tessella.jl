# `gmsh.model.occ` BRepGProp measures — `getMass` and `getCenterOfMass`
# (LinearProperties / SurfaceProperties / VolumeProperties) on the bound
# kernel graph, matching upstream's unsynchronized registry queries.
# Curves measure per-record: line/circle in closed form; ellipse, NURBS
# (`copied_curve`) and `surface_iso` by Gauss-Legendre quadrature;
# `degenerate` collapses to (0, point). Faces measure per-record: the
# plane's pcurve shoelace with conic-segment corrections and the
# outer-minus-holes convention, closed-form `∫P·dA` over canonical
# cylinder/cone/sphere/torus domains, and full-knot-domain tensor
# quadrature for NURBS patches (OCC tensor faces always cover their
# whole patch). Solids measure through the compact `primitive` record
# where present, the torus-sector Pappus integrals, or the divergence
# theorem over all-plane shells — anything else rejects loudly rather
# than approximating, per the codebase's reject-before-wrong-answer
# convention. Entities outside dims 1:3 follow upstream: the switch
# never fires, so the empty GProp reports mass 0 and centre (0,0,0).

# 8-point Gauss-Legendre abscissae/weights on [-1,1] — exact for
# polynomials through degree 15, ample for the conic/spline measures.
const _OCC_GL_T=(0.9602898564975363,0.7966664774136267,
                 0.5255324099163290,0.1834346424956498)
const _OCC_GL_W=(0.1012285362903763,0.2223810344533745,
                 0.3137066458778873,0.3626837833783620)

@inline _occ2_cross(a,b)=a[1]*b[2]-a[2]*b[1]
@inline _occ2_sub(a,b)=(a[1]-b[1],a[2]-b[2])
@inline _occ3_add(a,b)=(a[1]+b[1],a[2]+b[2],a[3]+b[3])
@inline _occ3_scale(a,s)=(a[1]*s,a[2]*s,a[3]*s)

# GL8 on each panel of [t0,t1] — `f` maps a parameter to (value, weight)
# and the accumulators ride the weighted sum. Panel subdivision keeps
# oscillatory speed integrands (ellipse eccentricity, spline knot
# boundaries) at machine-precision rather than single-panel accuracy.
function _model_occ_quad1d(visit,t0,t1,panels)
    panels<1 && return 0.0,(0.0,0.0,0.0)
    L=0.0; M=(0.0,0.0,0.0)
    span=(t1-t0)/panels
    for k in 0:panels-1
        a=t0+k*span; h=span/2; mid=a+h
        for i in 1:4, sgn in (-1,1)
            t=mid+sgn*h*_OCC_GL_T[i]; w=_OCC_GL_W[i]*h
            p,speed=visit(t)
            dw=w*speed
            L+=dw; M=_occ3_add(M,_occ3_scale(p,dw))
        end
    end
    return L,M
end

# Distinct interior breakpoints of a NURBS curve's knot range — the
# natural smoothness panels for the measure quadrature.
function _model_occ_curve_panels(curve,t0,t1)
    cuts=Float64[]
    for k in curve.curve.knots
        t0<k<t1 && (isempty(cuts) || k>cuts[end]) && push!(cuts,k)
    end
    return cuts
end

# Nonzero knot spans inside [t0,t1] as (a,b) pairs — the surface analog
# of `_model_occ_curve_panels`.
function _model_occ_span_cuts(knots,t0,t1)
    cuts=Float64[t0]
    for k in knots
        t0<k<t1 && k>cuts[end] && push!(cuts,k)
    end
    push!(cuts,t1)
    return Tuple{Float64,Float64}[(cuts[i],cuts[i+1])
                                  for i in 1:length(cuts)-1]
end

# ---- curves ----------------------------------------------------------

function _model_occ_curve_measure(kernel,id)
    rec=kernel.shapes[id].geometry.geometry
    if rec.occ===:line
        L=abs(rec.t1-rec.t0)*_occ_modulus(rec.dir)
        return L,_occ_line_point(rec,(rec.t0+rec.t1)/2)
    elseif rec.occ===:circle
        dθ=rec.t1-rec.t0
        L=rec.r*abs(dθ)
        # Arc centroid rides the bisector at r·sinc(Δθ/2) from center.
        f=dθ==0.0 ? 0.0 : sin(dθ/2)/(dθ/2)
        s,c=_gm_sincos((rec.t0+rec.t1)/2)
        return L,(fma(rec.r*f*c,rec.X[1],fma(rec.r*f*s,rec.Y[1],rec.center[1])),
                  fma(rec.r*f*c,rec.X[2],fma(rec.r*f*s,rec.Y[2],rec.center[2])),
                  fma(rec.r*f*c,rec.X[3],fma(rec.r*f*s,rec.Y[3],rec.center[3])))
    elseif rec.occ===:degenerate
        tag=kernel.shape_tags[id][2]
        return 0.0,_occ_degenerate_point(rec,tag,rec.t0)
    end
    evaluator=rec.occ===:ellipse ?
        t->(_occ_ellipse_point(rec,t),_occ_ellipse_derivative(rec,t)) :
        rec.occ===:copied_curve ?
        t->begin values=_model_occ_curve_d2(rec.curve,t)
            (values[1],values[2]) end :
        rec.occ===:surface_iso ?
        t->(_model_surface_iso_point(rec,t),_model_surface_iso_d1(rec,t)) :
        throw(ErrorException(
            "OCC.get_mass: unsupported OCC curve kind $(rec.occ)"))
    # GL8 per knot span for splines; 8 uniform panels for the conic/iso
    # families — panel subdivision keeps oscillatory speed integrands at
    # machine precision.
    cuts=rec.occ===:copied_curve ?
        _model_occ_curve_panels(rec.curve,rec.t0,rec.t1) : Float64[]
    bounds=Float64[rec.t0]
    append!(bounds,cuts)
    push!(bounds,rec.t1)
    L=0.0; M=(0.0,0.0,0.0)
    for k in 1:length(bounds)-1
        l,m=_model_occ_quad1d(t->begin p,d=evaluator(t)
            (p,_occ_modulus(d)) end,bounds[k],bounds[k+1],
            isempty(cuts) ? 8 : 1)
        L+=l; M=_occ3_add(M,m)
    end
    c=L>_OCC_CONFUSION ? _occ3_scale(M,1.0/L) :
        evaluator((rec.t0+rec.t1)/2)[1]
    return L,c
end

# ---- plane-face shoelace ---------------------------------------------

# One wire's raw signed shoelace area and first moments in the plane's
# physical uv frame. Reversed uses traverse the `fwd` pcurve backward;
# conic edges add their signed sector-minus-chord segment.
function _model_occ_wire_shoelace(kernel,face_rec,uses)
    area=0.0; Mx=0.0; My=0.0
    for (eid,rev) in uses
        etag=kernel.shape_tags[eid][2]
        eg=kernel.shapes[eid].geometry.geometry
        pc=face_rec.pcurves[etag].fwd
        t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
        if haskey(pc,:lin2d)
            a=_occ_pcurve_d0(pc,t0); b=_occ_pcurve_d0(pc,t1)
            cr=_occ2_cross(a,b)
            area+=cr/2
            Mx+=(a[1]+b[1])*cr/6
            My+=(a[2]+b[2])*cr/6
            continue
        end
        # Conic edge (circ2d/ellips2d). Region between chord and arc =
        # conic sector − chord triangle through the conic center.
        o,x,ra,rb=haskey(pc,:circ2d) ?
            (pc.circ2d[1],pc.circ2d[2],pc.circ2d[3],pc.circ2d[3]) :
            (pc.ellips2d[1],pc.ellips2d[2],pc.ellips2d[3],pc.ellips2d[4])
        a=_occ_pcurve_d0(pc,t0); b=_occ_pcurve_d0(pc,t1)
        Δ=t1-t0
        cr=_occ2_cross(a,b)
        area+=cr/2
        Mx+=(a[1]+b[1])*cr/6
        My+=(a[2]+b[2])*cr/6
        sec_area=ra*rb*Δ/2
        tri_area=_occ2_cross(_occ2_sub(a,o),_occ2_sub(b,o))/2
        seg=sec_area-tri_area
        area+=seg
        if abs(seg)>0.0
            θm=(t0+t1)/2
            # Sector centroid distance 4·sin(Δθ/2)/(3·Δθ) along the
            # mid-parameter direction, affine-scaled per semi-axis; the
            # 2-D frame's derived YDirection is (−x.y, x.x).
            du=abs(Δ)>1e-7 ? 4.0*sin(Δ/2)/(3*Δ) : 2.0/3
            s,c=_gm_sincos(θm)
            sec_c=(o[1]+du*(ra*c*x[1]-rb*s*x[2]),
                   o[2]+du*(ra*c*x[2]+rb*s*x[1]))
            tri_c=((o[1]+a[1]+b[1])/3,(o[2]+a[2]+b[2])/3)
            Mx+=sec_area*sec_c[1]-tri_area*tri_c[1]
            My+=sec_area*sec_c[2]-tri_area*tri_c[2]
        end
    end
    return area,Mx,My
end

@inline _model_occ_plane_eval(rec,u,v)=
    (fma(u,rec.X[1],v*rec.Y[1])+rec.center[1],
     fma(u,rec.X[2],v*rec.Y[2])+rec.center[2],
     fma(u,rec.X[3],v*rec.Y[3])+rec.center[3])

# Plane-face area, 3-D centroid and second-moment tensor T_F = ∬_F r⊗r
# dA — the tensor via the signed fan triangulation from the face's own
# shoelace centroid (exact for any line-bounded polygon, holes included;
# the midpoint rule integrates quadratics exactly). Conic-edged solid
# faces fall back to (area, centroid, nothing).
function _model_occ_plane_stats(kernel,id;want_tensor=false)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    A=0.0
    raw_per_wire=NamedTuple[]
    for wire_id in shape.children
        uses=kernel.shapes[wire_id].geometry.uses
        push!(raw_per_wire,(uses=uses,
            raw=_model_occ_wire_shoelace(kernel,rec,uses)))
    end
    # Pass 1: face area + uv centroid under the |outer|−Σ|holes| rule.
    for (index,wire) in enumerate(raw_per_wire)
        area=abs(wire.raw[1])
        A+=index==1 ? area : -area
    end
    A>_OCC_CONFUSION || throw(ErrorException(
        "OCC.get_mass: plane face has zero area"))
    cu=cv=0.0
    for (index,wire) in enumerate(raw_per_wire)
        area=abs(wire.raw[1])
        area>_OCC_CONFUSION || continue
        w=(index==1 ? area : -area)/A
        cu+=w*wire.raw[2]/area; cv+=w*wire.raw[3]/area
    end
    c3=_model_occ_plane_eval(rec,cu,cv)
    want_tensor || return A,c3,nothing
    T=zeros(3,3)
    # Pass 2: per-wire signed fan from (cu,cv); each wire's contribution
    # is its raw (winding-signed) tensor scaled by sign(raw area).
    for (index,wire) in enumerate(raw_per_wire)
        sgn=wire.raw[1]>=0 ? 1.0 : -1.0
        w=(index==1 ? 1.0 : -1.0)*sgn
        c2=(cu,cv); a3=c3
        for (eid,rev) in wire.uses
            etag=kernel.shape_tags[eid][2]
            eg=kernel.shapes[eid].geometry.geometry
            pc=rec.pcurves[etag].fwd
            haskey(pc,:lin2d) || throw(ErrorException(
                "OCC.get_center_of_mass: solid faces with conic " *
                "boundaries are not supported"))
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            p1=_occ_pcurve_d0(pc,t0); p2=_occ_pcurve_d0(pc,t1)
            d1=_occ2_sub(p1,c2); d2=_occ2_sub(p2,c2)
            A2=_occ2_cross(d1,d2)/2
            b3=_model_occ_plane_eval(rec,p1[1],p1[2])
            e3=_model_occ_plane_eval(rec,p2[1],p2[2])
            mab=_occ3_scale(_occ3_add(a3,b3),0.5)
            mbc=_occ3_scale(_occ3_add(b3,e3),0.5)
            mca=_occ3_scale(_occ3_add(e3,a3),0.5)
            f=w*A2/3
            for m in (mab,mbc,mca)
                for i in 1:3, j in 1:3
                    T[i,j]+=f*m[i]*m[j]
                end
            end
        end
    end
    return A,c3,T
end

# ---- faces -----------------------------------------------------------

function _model_occ_face_measure(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    occ=get(rec,:occ,nothing)
    if occ===:plane
        A,c3,_=_model_occ_plane_stats(kernel,id)
        return A,c3
    elseif haskey(rec,:patch)
        # NURBS tensor face — OCC tensor faces cover the full knot
        # domain. GL8 per nonzero knot span in each direction keeps the
        # piecewise-smooth patch integrand at machine precision.
        surface=rec.patch.surface
        (u0,v0),(u1,v1)=_model_nurbs_surface_bounds(rec)
        uspans=_model_occ_span_cuts(surface.knots_u,u0,u1)
        vspans=_model_occ_span_cuts(surface.knots_v,v0,v1)
        A=0.0; M=(0.0,0.0,0.0)
        for (ua,ub) in uspans, (va,vb) in vspans
            hu=(ub-ua)/2; mu=(ua+ub)/2; hv=(vb-va)/2; mv=(va+vb)/2
            for i in 1:4, si in (-1,1), j in 1:4, sj in (-1,1)
                u=mu+si*hu*_OCC_GL_T[i]; v=mv+sj*hv*_OCC_GL_T[j]
                w=_OCC_GL_W[i]*_OCC_GL_W[j]*hu*hv
                p,du,dv,_,_,_=_model_nurbs_surface_d2(rec,u,v)
                dA=_occ_modulus(_occ_cross(du,dv))
                A+=w*dA; M=_occ3_add(M,_occ3_scale(p,w*dA))
            end
        end
        A>_OCC_CONFUSION || throw(ErrorException(
            "OCC.get_mass: NURBS face has zero area"))
        return A,_occ3_scale(M,1.0/A)
    end
    # Quadric patches over their canonical domains — the closed-form
    # ∫P·dA integrals (u-mean cos/sin cover the partial torus sweep).
    u0,u1,v0,v1=if occ===:cylinder
        (0.0,_OCC_TWO_PI,0.0,rec.height)
    elseif occ===:sphere
        (0.0,_OCC_TWO_PI,-π/2,π/2)
    elseif occ===:cone
        (0.0,_OCC_TWO_PI,0.0,
         sqrt(rec.height*rec.height+(rec.r1-rec.r2)*(rec.r1-rec.r2)))
    elseif occ===:torus
        (0.0,rec.angle,0.0,_OCC_TWO_PI)
    else
        throw(ErrorException(
            "OCC.get_mass: unsupported OCC surface kind $occ"))
    end
    du=u1-u0; dv=v1-v0
    cu=du==0.0 ? 0.0 : (sin(u1)-sin(u0))/du
    su=du==0.0 ? 0.0 : (cos(u0)-cos(u1))/du
    X,Y,Z,C=rec.X,_occ_frame_y(rec),rec.axis,rec.center
    radial=k->(X[1]*k*cu+Y[1]*k*su,X[2]*k*cu+Y[2]*k*su,X[3]*k*cu+Y[3]*k*su)
    Sv=sin(v1)-sin(v0)
    if occ===:cylinder
        A=rec.radius*du*dv
        ax=Z.*((v0+v1)/2); rad=radial(rec.radius)
        return A,_occ3_add(C,_occ3_add(ax,rad))
    elseif occ===:cone
        sa=_gm_atan((rec.r2-rec.r1)/rec.height)
        ssa,csa=_gm_sincos(sa)
        Ir=rec.r1*dv+ssa*(v1*v1-v0*v0)/2
        Ir2=rec.r1*rec.r1*dv+rec.r1*ssa*(v1*v1-v0*v0)+
            ssa*ssa*(v1*v1*v1-v0*v0*v0)/3
        Iv1=rec.r1*(v1*v1-v0*v0)/2+ssa*(v1*v1*v1-v0*v0*v0)/3
        A=du*Ir
        ax=Z.*(csa*Iv1/Ir); rad=radial(Ir2/Ir)
        return A,_occ3_add(C,_occ3_add(ax,rad))
    elseif occ===:sphere
        A=rec.radius*rec.radius*du*Sv
        Jv=(dv+sin(v1)*cos(v1)-sin(v0)*cos(v0))/2
        ax=Z.*(rec.radius*(sin(v1)+sin(v0))/2)
        rad=radial(rec.radius*Jv/Sv)
        return A,_occ3_add(C,_occ3_add(ax,rad))
    else # torus
        Jv=(dv+sin(v1)*cos(v1)-sin(v0)*cos(v0))/2
        Ir=rec.r1*dv+rec.r2*Sv
        Ir2=rec.r1*rec.r1*dv+2*rec.r1*rec.r2*Sv+rec.r2*rec.r2*Jv
        Jz=rec.r1*(cos(v0)-cos(v1))+
           rec.r2*(sin(v1)^2-sin(v0)^2)/2
        A=rec.r2*du*Ir
        ax=Z.*(rec.r2*Jz/Ir); rad=radial(Ir2/Ir)
        return A,_occ3_add(C,_occ3_add(ax,rad))
    end
end

# ---- solids ----------------------------------------------------------

# ∬_F r⊗r dA · n contraction for one plane face: T·n̂ accumulated during
# the tensor pass is what the divergence centroid needs.
_occ_tensor_apply(T,n)=ntuple(i->T[i,1]*n[1]+T[i,2]*n[2]+T[i,3]*n[3],3)

function _model_occ_solid_measure(kernel,id;want_centroid=true)
    shape=kernel.shapes[id]; g=shape.geometry
    if haskey(g,:primitive)
        r=g.primitive.record
        if g.primitive.occ===:cylinder
            n=_occ_dir(r.axis)
            return π*r.radius*r.radius*r.height,
                _occ3_add(r.center,_occ3_scale(n,r.height/2))
        elseif g.primitive.occ===:sphere
            return 4π/3*r.radius^3,r.center
        elseif g.primitive.occ===:cone
            q=r.r1*r.r1+r.r1*r.r2+r.r2*r.r2
            frac=(r.r1*r.r1+2*r.r1*r.r2+3*r.r2*r.r2)/(4*q)
            n=_occ_dir(r.axis)
            return π*r.height*q/3,
                _occ3_add(r.center,_occ3_scale(n,r.height*frac))
        end
        throw(ErrorException(
            "OCC.get_mass: unsupported primitive kind $(g.primitive.occ)"))
    end
    # Collect (face_id, use_reversed) across all shells. An empty solid
    # accumulates nothing — the empty-GProp (0, origin) upstream reports.
    uses=Tuple{UInt64,Bool}[]
    for shell_id in shape.children
        for (face_id,rev) in kernel.shapes[shell_id].geometry.uses
            push!(uses,(face_id,rev))
        end
    end
    isempty(uses) && return 0.0,(0.0,0.0,0.0)
    face_kind=u->get(kernel.shapes[u[1]].geometry.geometry,:occ,nothing)
    torus=count(u->face_kind(u)===:torus,uses)
    if torus>0
        # Torus sector — the swept profile's Pappus integrals. V =
        # φ·πr2²·r1; the centroid rides the sector bisector at
        # (r1+r2²/4r1)·sinc(φ/2) from center.
        (torus==1 && all(u->face_kind(u)===:torus ||
            face_kind(u)===:plane,uses)) ||
            throw(ErrorException(
                "OCC.get_mass: torus solid with non-canonical boundary"))
        rec=only(kernel.shapes[fid].geometry.geometry
                 for (fid,_) in uses
                 if get(kernel.shapes[fid].geometry.geometry,:occ,
                        nothing)===:torus)
        φ=rec.angle
        V=φ*π*rec.r2*rec.r2*rec.r1
        mean_r=rec.r1+rec.r2*rec.r2/(4*rec.r1)
        arm=φ==0.0 ? 0.0 : mean_r*sin(φ/2)/(φ/2)
        s,c=_gm_sincos(φ/2)
        return V,_occ3_add(rec.center,
            (arm*c*rec.X[1]+arm*s*rec.Y[1],
             arm*c*rec.X[2]+arm*s*rec.Y[2],
             arm*c*rec.X[3]+arm*s*rec.Y[3]))
    end
    all(u->face_kind(u)===:plane,uses) ||
        throw(ErrorException("OCC.get_mass: unsupported solid boundary"))
    # Divergence theorem over outward-oriented plane faces:
    #   V = (1/3)Σ A·(n̂·c_F),  ∭r dV = (1/4)Σ n̂·∬_F r⊗r dA.
    V=0.0; M=(0.0,0.0,0.0)
    for (fid,rev) in uses
        rec=kernel.shapes[fid].geometry.geometry
        σ=rev ? -1.0 : 1.0
        n=(σ*rec.axis[1],σ*rec.axis[2],σ*rec.axis[3])
        # The tensor pass (conic-free wires only) runs only for the
        # centroid — get_mass on a conic-faced solid still reports V.
        A,c3,T=_model_occ_plane_stats(kernel,fid;
                                      want_tensor=want_centroid)
        V+=A*_dot3(n,c3)/3
        want_centroid && (M=_occ3_add(M,_occ_tensor_apply(T,n)))
    end
    # Upstream reports the signed divergence sum (inward-wound shells go
    # negative); a degenerate one leaves the empty-GProp centroid.
    abs(V)>_OCC_CONFUSION || return V,(0.0,0.0,0.0)
    return V,_occ3_scale(M,0.25/V)
end

# ---- public endpoints -------------------------------------------------

function _model_occ_measure_entity(m::GeoModel,dim,tag,caller;
                                   want_centroid=true)
    d=_model_occ_dimension(dim,caller)
    t=_model_occ_requested_tag(tag,caller)
    kernel=m.occ_kernel
    id=get(kernel.bindings,(d,t),nothing)
    id===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE entity of dimension $d with tag $t"))
    if d==1
        return _model_occ_curve_measure(kernel,id)
    elseif d==2
        return _model_occ_face_measure(kernel,id)
    elseif d==3
        return _model_occ_solid_measure(kernel,id;
                                        want_centroid=want_centroid)
    end
    return 0.0,(0.0,0.0,0.0)
end

_model_occ_get_mass(m,dim,tag)=first(
    _model_occ_measure_entity(m,dim,tag,"OCC.get_mass";want_centroid=false))
_model_occ_get_center_of_mass(m,dim,tag)=last(
    _model_occ_measure_entity(m,dim,tag,"OCC.get_center_of_mass"))
