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
        if haskey(pc,:neg2d)
            pc=pc.neg2d[1];t0,t1=-t0,-t1
        end
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
function _model_occ_plane_stats(kernel,id;want_tensor=false,
                                caller="OCC.get_center_of_mass")
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
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            if haskey(pc,:neg2d)
                pc=pc.neg2d[1];t0,t1=-t0,-t1
            end
            haskey(pc,:lin2d) || throw(ErrorException(
                "$caller: faces with conic " *
                "boundaries are not supported"))
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

# ---- cylinder patches ------------------------------------------------
#
# A trimmed cylinder face p(u,v)=C+R·u_r+v·Z has surface element
# dA=R·du·dv and unit normal ±u_r (the frame handedness h=X·(Y×Z)
# picks the sign). The divergence needs ∬_D(p·u_r)dA with
# p·u_r=Cx·cos u+Cy·sin u+R, Cx=C·X, Cy=C·Y — every domain moment
# reduces by Green's theorem to a boundary line integral over the
# face's pcurves:
#   ∬_D dudv=∮u dv,  ∬_D cos u dudv=∮sin u dv,  ∬_D sin u dudv=∮−cos u dv.
# Each wire's signed integral is normalized by its raw shoelace sign so
# the outer wire counts positive and holes negative, matching the
# plane-face convention.
function _model_occ_cyl_face_div(kernel,fid)
    shape=kernel.shapes[fid]; rec=shape.geometry.geometry
    I1=Ic=Is=0.0
    for (wi,wire_id) in enumerate(shape.children)
        uses=kernel.shapes[wire_id].geometry.uses
        raw=_model_occ_wire_shoelace(kernel,rec,uses)[1]
        j1=jc=js=0.0
        u_prev=nothing
        for (eid,rev) in uses
            etag=kernel.shape_tags[eid][2]
            eg=kernel.shapes[eid].geometry.geometry
            pc=rec.pcurves[etag].fwd
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            # u is periodic: a stored pcurve picks an arbitrary
            # representative mod 2π — unwrap the wire so each edge
            # continues the previous one's u, or ∮u dv telescopes to
            # zero instead of the domain area
            us=_occ_pcurve_d0(pc,t0)[1]
            shift=u_prev===nothing ? 0.0 :
                _OCC_TWO_PI*round((u_prev-us)/_OCC_TWO_PI)
            panels=max(1,ceil(Int,abs(t1-t0)/(π/2)))
            span=(t1-t0)/panels
            for k in 0:panels-1
                a=t0+k*span; h=span/2; mid=a+h
                for i in 1:4, s in (-1,1)
                    t=mid+s*h*_OCC_GL_T[i]; w=_OCC_GL_W[i]*h
                    (u,v),(du,dv)=_occ_pcurve_d1(pc,t)
                    u+=shift
                    j1+=w*u*dv; jc+=w*sin(u)*dv; js+=w*(-cos(u)*dv)
                end
            end
            u_prev=_occ_pcurve_d0(pc,t1)[1]+shift
        end
        w=(wi==1 ? 1.0 : -1.0)*(raw>=0 ? 1.0 : -1.0)
        I1+=w*j1; Ic+=w*jc; Is+=w*js
    end
    return I1,Ic,Is
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
    # Partial sweeps carry `angle` (azimuth span) and `lat1`/`lat2`
    # (sphere latitude range) in the record.
    u0,u1,v0,v1=if occ===:cylinder
        (0.0,get(rec,:angle,_OCC_TWO_PI),0.0,rec.height)
    elseif occ===:sphere
        (0.0,get(rec,:angle,_OCC_TWO_PI),
         get(rec,:lat1,-π/2),get(rec,:lat2,π/2))
    elseif occ===:cone
        (0.0,get(rec,:angle,_OCC_TWO_PI),0.0,
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

# ---- sector solids ----------------------------------------------------
#
# Partial primitive sweeps carry a `sector` record instead of a compact
# `primitive` — the solid is the wedge
#   cylinder/cone: {c + ρ·E(u) + z·n : u∈[0,a], ρ∈[0,R(z)], z∈[0,h]}
#   sphere:        {c + ρ·d(u,v) : u∈[0,a], v∈[v1,v2], ρ∈[0,R]}
# with E(u)=cos u·X+sin u·Y and d=cos v·E+sin v·n. V, M=∭p dV and
# Q=∭p⊗p dV come from closed-form azimuth integrals; the cone's axial
# moments ∫z^m R(z)^k dz are polynomial in z (degree ≤ 4), so a 4-point
# Gauss rule is exact.
#
# Azimuth moments over [0,a]: U0=∫du, Uc=∫cos u, Us=∫sin u and the
# second-order Uc2/Us2/Ucs. `Uv` collects Uc·X+Us·Y as a 3-vector, `UE`
# the symmetric tensor ∫E⊗E du.
function _occ_sector_u(a,X,Y)
    s,c=_gm_sincos(a)
    U0=a; Uc=s; Us=1.0-c
    s2,c2=_gm_sincos(2a)
    Uc2=a/2+s2/4; Us2=a/2-s2/4; Ucs=s*s/2
    Uv=(Uc*X[1]+Us*Y[1],Uc*X[2]+Us*Y[2],Uc*X[3]+Us*Y[3])
    UE=zeros(3,3)
    for i in 1:3,j in 1:3
        UE[i,j]=Uc2*X[i]*X[j]+Ucs*(X[i]*Y[j]+Y[i]*X[j])+Us2*Y[i]*Y[j]
    end
    return U0,Uv,UE
end

# Axisymmetric wedge (cylinder/cone sector): R(z)=r1+dr·z/h. Returns
# V, first-moment vector M=∭p dV and second-moment tensor Q=∭p⊗p dV
# in the local frame anchored at the world center `c` (p measured from
# the origin; `c` folds in through the caller's shift).
function _occ_sector_stats_axis(rec,U0,Uv,UE)
    n=_occ_dir(rec.axis)
    r1,r2,h=rec.occ===:cylinder ?
        (rec.radius,rec.radius,rec.height) : (rec.r1,rec.r2,rec.height)
    dr=r2-r1
    V=0.0
    M=(0.0,0.0,0.0)
    Q=zeros(3,3)
    for i in 1:4,s in (-1,1)
        z=h*(0.5+0.5*s*_OCC_GL_T[i]); w=_OCC_GL_W[i]*h/2
        R=r1+dr*(z/h)
        R2=R*R; R3=R2*R; R4=R3*R
        V+=w*U0*R2/2
        mz=w*z*U0*R2/2
        mr=w*R3/3
        M=(M[1]+mr*Uv[1]+mz*n[1],M[2]+mr*Uv[2]+mz*n[2],
           M[3]+mr*Uv[3]+mz*n[3])
        qE=w*R4/4; qx=w*z*R3/3; qn=w*z*z*U0*R2/2
        for i2 in 1:3,j2 in 1:3
            Q[i2,j2]+=qE*UE[i2,j2]+qx*(Uv[i2]*n[j2]+n[i2]*Uv[j2])+
                qn*n[i2]*n[j2]
        end
    end
    return V,M,Q
end

# Spherical sector solid — a *cylindrical* wedge: the flat latitude-cap
# planes bound {z1≤z≤z2, ρ≤√(R²−z²), u∈[0,a]}, not a spherical-coordinate
# ρ-wedge. All needed integrals ∫z^m·ρ(z)^k dz are closed-form:
#   ρ²=R²−z² and ρ⁴ are polynomial; ρ³ uses the √(a²−x²)^{3/2}
#   antiderivative; z·ρ³ is the exact power form −(R²−z²)^{5/2}/5.
function _occ_sector_stats_sph(rec,U0,Uv,UE)
    n=_occ_dir(rec.axis)
    R=rec.radius
    z1=R*_gm_sincos(rec.lat1)[1]; z2=R*_gm_sincos(rec.lat2)[1]
    R2=R*R; R4=R2*R2
    rho(z)=sqrt(max(R2-z*z,0.0))
    I02=(R2*z2-z2*z2*z2/3)-(R2*z1-z1*z1*z1/3)
    I12=(R2*z2*z2/2-z2^4/4)-(R2*z1*z1/2-z1^4/4)
    I03=(z2/8*(5*R2-2*z2*z2)*rho(z2)+3*R4/8*asin(clamp(z2/R,-1,1)))-
        (z1/8*(5*R2-2*z1*z1)*rho(z1)+3*R4/8*asin(clamp(z1/R,-1,1)))
    I13=(rho(z1)^5-rho(z2)^5)/5
    I22=(R2*z2*z2*z2/3-z2^5/5)-(R2*z1*z1*z1/3-z1^5/5)
    I04=(R4*z2-2*R2*z2*z2*z2/3+z2^5/5)-
        (R4*z1-2*R2*z1*z1*z1/3+z1^5/5)
    V=U0*I02/2
    M=(Uv[1]*I03/3+U0*n[1]*I12/2,
       Uv[2]*I03/3+U0*n[2]*I12/2,
       Uv[3]*I03/3+U0*n[3]*I12/2)
    Q=zeros(3,3)
    for i in 1:3,j in 1:3
        Q[i,j]=UE[i,j]*I04/4+(Uv[i]*n[j]+n[i]*Uv[j])*I13/3+
               U0*n[i]*n[j]*I22/2
    end
    return V,M,Q
end

# Sector solid (V, centroid, Q) — centroid measured in world coords.
function _occ_sector_stats(rec)
    U0,Uv,UE=_occ_sector_u(rec.angle,rec.X,rec.Y)
    V,M,Q=rec.occ===:sphere ? _occ_sector_stats_sph(rec,U0,Uv,UE) :
        _occ_sector_stats_axis(rec,U0,Uv,UE)
    c=rec.center
    # Shift p→c+p: M→M+V·c, Q→Q+c⊗M+M⊗cᵀ+V·c⊗cᵀ.
    M=(M[1]+V*c[1],M[2]+V*c[2],M[3]+V*c[3])
    for i in 1:3,j in 1:3
        Q[i,j]+=c[i]*M[j]+M[i]*c[j]-V*c[i]*c[j]
    end
    abs(V)>_OCC_CONFUSION || return V,(0.0,0.0,0.0),Q
    return V,_occ3_scale(M,1.0/V),Q
end

function _model_occ_solid_measure(kernel,id;want_centroid=true)
    shape=kernel.shapes[id]; g=shape.geometry
    if haskey(g,:sector)
        V,c,_=_occ_sector_stats(g.sector)
        return V,want_centroid ? c : (0.0,0.0,0.0)
    end
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
        elseif g.primitive.occ===:revolved
            haskey(r,:angle) || return r.V,
                _occ3_add(r.point,_occ3_scale(r.axis,r.h_c))
            # Partial turn — the centroid rides the sector bisector at
            # `arm` off the axis in addition to h_c along it.
            return r.V,_occ3_add(r.point,_occ3_add(
                _occ3_scale(r.axis,r.h_c),
                _occ3_scale(r.bisector,r.arm)))
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
    all(u->face_kind(u)===:plane || face_kind(u)===:cylinder,uses) ||
        throw(ErrorException("OCC.get_mass: unsupported solid boundary"))
    # Divergence theorem over outward-oriented faces:
    #   V = (1/3)Σ ∬_F r·n̂ dA,  ∭r dV = (1/4)Σ ∬_F r(r·n̂) dA.
    # Plane faces integrate exactly through `_model_occ_plane_stats`;
    # cylinder patches reduce to wire-boundary line integrals.
    V=0.0; M=(0.0,0.0,0.0)
    for (fid,rev) in uses
        rec=kernel.shapes[fid].geometry.geometry
        σ=rev ? -1.0 : 1.0
        if get(rec,:occ,nothing)===:cylinder
            want_centroid && throw(ErrorException(
                "OCC.get_center_of_mass: cylinder-bounded solids are " *
                "not supported"))
            I1,Ic,Is=_model_occ_cyl_face_div(kernel,fid)
            h=_dot3(_occ_cross(rec.X,_occ_frame_y(rec)),rec.axis)>=0 ?
                1.0 : -1.0
            Cx=_dot3(rec.center,rec.X); Cy=_dot3(rec.center,rec.Y)
            V+=σ*h*rec.radius*(Cx*Ic+Cy*Is+rec.radius*I1)/3
            continue
        end
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

# ---- second moments: get_matrix_of_inertia -----------------------------
#
# Upstream's `GProps::MatrixOfInertia` is the inertia tensor about the
# origin, I = tr(Q)I₃ − Q with Q = ∫r⊗r dμ under the same per-dim
# Linear/Surface/VolumeProperties switch. Curves and faces integrate
# r⊗r per record — closed forms for line/circle and the plane tensor
# path, GL8 panels for the quadrature families (the canonical quadric
# domain, NURBS knot spans). Solids reuse the measure split: compact
# `primitive` records in closed form (sphere 2mr²/5·I, cylinder/cone
# axial forms through the ∭x²dV/∭z²dV layer integrals), and all-plane
# shells through the divergence theorem — ∭x_i x_j dV = n_k∬x_i x_j x_k
# dA/5 with a signed-fan cubic moment exact for line-bounded polygons
# (the degree-3 Stroud rule on the same fan). Conic wires, quadric-
# bounded shells and other unsupported boundaries reject like
# `get_center_of_mass`; non-1:3 dims and empty shells report the
# empty-GProp zero matrix.

function _model_occ_curve_moment2(kernel,id)
    rec=kernel.shapes[id].geometry.geometry
    m=zeros(3,3)
    if rec.occ===:line
        o,d=rec.origin,rec.dir
        md=_occ_modulus(d)
        lo,hi=minmax(rec.t0,rec.t1)
        a=md*(hi-lo); b=md*(hi*hi-lo*lo)/2; c=md*(hi*hi*hi-lo*lo*lo)/3
        for i in 1:3,j in 1:3
            m[i,j]=a*o[i]*o[j]+b*(o[i]*d[j]+d[i]*o[j])+c*d[i]*d[j]
        end
        return m
    elseif rec.occ===:circle
        lo,hi=minmax(rec.t0,rec.t1); Δ=hi-lo
        s0,c0=_gm_sincos(lo); s1,c1=_gm_sincos(hi)
        Su=s1-s0; Cv=c0-c1
        C2=(Δ+s1*c1-s0*c0)/2; S2=(Δ-s1*c1+s0*c0)/2
        CS=(s1*s1-s0*s0)/2
        C,X,Y,r=rec.center,rec.X,rec.Y,rec.r
        for i in 1:3,j in 1:3
            m[i,j]=r*(C[i]*C[j]*Δ +
                r*(C[i]*(X[j]*Su+Y[j]*Cv)+(X[i]*Su+Y[i]*Cv)*C[j]) +
                r*r*(C2*X[i]*X[j]+CS*(X[i]*Y[j]+Y[i]*X[j])+S2*Y[i]*Y[j]))
        end
        return m
    elseif rec.occ===:degenerate
        return m
    end
    evaluator=rec.occ===:ellipse ?
        t->(_occ_ellipse_point(rec,t),_occ_ellipse_derivative(rec,t)) :
        rec.occ===:copied_curve ?
        t->begin values=_model_occ_curve_d2(rec.curve,t)
            (values[1],values[2]) end :
        rec.occ===:surface_iso ?
        t->(_model_surface_iso_point(rec,t),_model_surface_iso_d1(rec,t)) :
        throw(ErrorException("OCC.get_matrix_of_inertia: unsupported " *
            "OCC curve kind $(rec.occ)"))
    cuts=rec.occ===:copied_curve ?
        _model_occ_curve_panels(rec.curve,rec.t0,rec.t1) : Float64[]
    bounds=Float64[rec.t0]
    append!(bounds,cuts)
    push!(bounds,rec.t1)
    for k in 1:length(bounds)-1
        a,bk=bounds[k],bounds[k+1]
        panels=isempty(cuts) ? 8 : 1
        span=(bk-a)/panels
        for p in 0:panels-1
            aa=a+p*span; h=span/2; mid=aa+h
            for i in 1:4, sgn in (-1,1)
                t=mid+sgn*h*_OCC_GL_T[i]; w=_OCC_GL_W[i]*h
                pt,d=evaluator(t)
                dw=w*_occ_modulus(d)
                for i2 in 1:3,j2 in 1:3
                    m[i2,j2]+=dw*pt[i2]*pt[j2]
                end
            end
        end
    end
    return m
end

# Cubic boundary moment U_F[i,j,k] = ∬_F x_i x_j x_k dA for a plane
# face — the divergence theorem needs ∭x_i x_j dV = n_k U_F[i,j,k]/5.
# Same signed fan as `_model_occ_plane_stats`, evaluated with the
# degree-3 Stroud rule (exact for cubics, hence also for the fan's
# per-triangle integrand).
function _model_occ_plane_moment3(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    raw_per_wire=NamedTuple[]
    for wire_id in shape.children
        uses=kernel.shapes[wire_id].geometry.uses
        push!(raw_per_wire,(uses=uses,
            raw=_model_occ_wire_shoelace(kernel,rec,uses)))
    end
    A=0.0
    for (index,wire) in enumerate(raw_per_wire)
        A+=index==1 ? abs(wire.raw[1]) : -abs(wire.raw[1])
    end
    A>_OCC_CONFUSION || throw(ErrorException(
        "OCC.get_matrix_of_inertia: plane face has zero area"))
    cu=cv=0.0
    for (index,wire) in enumerate(raw_per_wire)
        area=abs(wire.raw[1])
        area>_OCC_CONFUSION || continue
        w=(index==1 ? area : -area)/A
        cu+=w*wire.raw[2]/area; cv+=w*wire.raw[3]/area
    end
    c2=(cu,cv); a3=_model_occ_plane_eval(rec,cu,cv)
    U=zeros(3,3,3)
    for (index,wire) in enumerate(raw_per_wire)
        sgn=wire.raw[1]>=0 ? 1.0 : -1.0
        w=(index==1 ? 1.0 : -1.0)*sgn
        for (eid,rev) in wire.uses
            etag=kernel.shape_tags[eid][2]
            eg=kernel.shapes[eid].geometry.geometry
            pc=rec.pcurves[etag].fwd
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            if haskey(pc,:neg2d)
                pc=pc.neg2d[1];t0,t1=-t0,-t1
            end
            haskey(pc,:lin2d) || throw(ErrorException(
                "OCC.get_matrix_of_inertia: faces with conic " *
                "boundaries are not supported"))
            p1=_occ_pcurve_d0(pc,t0); p2=_occ_pcurve_d0(pc,t1)
            A2=_occ2_cross(_occ2_sub(p1,c2),_occ2_sub(p2,c2))/2
            b3=_model_occ_plane_eval(rec,p1[1],p1[2])
            e3=_model_occ_plane_eval(rec,p2[1],p2[2])
            ctr=_occ3_scale(_occ3_add(a3,_occ3_add(b3,e3)),1/3)
            f=w*A2
            for i in 1:3,j in 1:3,k in 1:3
                acc=9/20*ctr[i]*ctr[j]*ctr[k]+
                    1/20*(a3[i]*a3[j]*a3[k]+b3[i]*b3[j]*b3[k]+
                          e3[i]*e3[j]*e3[k])
                mab=(a3[i]+b3[i])*(a3[j]+b3[j])*(a3[k]+b3[k])
                mbc=(b3[i]+e3[i])*(b3[j]+e3[j])*(b3[k]+e3[k])
                mca=(e3[i]+a3[i])*(e3[j]+a3[j])*(e3[k]+a3[k])
                U[i,j,k]+=f*(acc+(mab+mbc+mca)/60)
            end
        end
    end
    return U
end

function _model_occ_face_moment2(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    occ=get(rec,:occ,nothing)
    m=zeros(3,3)
    if occ===:plane
        _,_,T=_model_occ_plane_stats(kernel,id;want_tensor=true,
            caller="OCC.get_matrix_of_inertia")
        return T
    elseif haskey(rec,:patch)
        surface=rec.patch.surface
        (u0,v0),(u1,v1)=_model_nurbs_surface_bounds(rec)
        uspans=_model_occ_span_cuts(surface.knots_u,u0,u1)
        vspans=_model_occ_span_cuts(surface.knots_v,v0,v1)
        for (ua,ub) in uspans, (va,vb) in vspans
            hu=(ub-ua)/2; mu=(ua+ub)/2; hv=(vb-va)/2; mv=(va+vb)/2
            for i in 1:4, si in (-1,1), j in 1:4, sj in (-1,1)
                u=mu+si*hu*_OCC_GL_T[i]; v=mv+sj*hv*_OCC_GL_T[j]
                w=_OCC_GL_W[i]*_OCC_GL_W[j]*hu*hv
                p,du,dv,_,_,_=_model_nurbs_surface_d2(rec,u,v)
                dA=_occ_modulus(_occ_cross(du,dv))
                for i2 in 1:3,j2 in 1:3
                    m[i2,j2]+=w*dA*p[i2]*p[j2]
                end
            end
        end
        return m
    end
    # Quadric patches over their canonical domains — 8×8 GL8 panels,
    # matching the curve-panel resolution for these smooth integrands.
    u0,u1,v0,v1=if occ===:cylinder
        (0.0,get(rec,:angle,_OCC_TWO_PI),0.0,rec.height)
    elseif occ===:sphere
        (0.0,get(rec,:angle,_OCC_TWO_PI),
         get(rec,:lat1,-π/2),get(rec,:lat2,π/2))
    elseif occ===:cone
        (0.0,get(rec,:angle,_OCC_TWO_PI),0.0,
         sqrt(rec.height*rec.height+(rec.r1-rec.r2)*(rec.r1-rec.r2)))
    elseif occ===:torus
        (0.0,rec.angle,0.0,_OCC_TWO_PI)
    else
        throw(ErrorException("OCC.get_matrix_of_inertia: unsupported " *
            "OCC surface kind $occ"))
    end
    for pu in 0:7, pv in 0:7
        ua=u0+pu*(u1-u0)/8; ub=ua+(u1-u0)/8
        va=v0+pv*(v1-v0)/8; vb=va+(v1-v0)/8
        hu=(ub-ua)/2; mu=(ua+ub)/2; hv=(vb-va)/2; mv=(va+vb)/2
        for i in 1:4, si in (-1,1), j in 1:4, sj in (-1,1)
            u=mu+si*hu*_OCC_GL_T[i]; v=mv+sj*hv*_OCC_GL_T[j]
            w=_OCC_GL_W[i]*_OCC_GL_W[j]*hu*hv
            p,du,dv=_occ_surface_d1(rec,u,v)
            dA=_occ_modulus(_occ_cross(du,dv))
            for i2 in 1:3,j2 in 1:3
                m[i2,j2]+=w*dA*p[i2]*p[j2]
            end
        end
    end
    return m
end

function _model_occ_solid_moment2(kernel,id)
    shape=kernel.shapes[id]; g=shape.geometry
    m=zeros(3,3)
    if haskey(g,:sector)
        _,_,Q=_occ_sector_stats(g.sector)
        return Q
    end
    if haskey(g,:primitive)
        r=g.primitive.record
        if g.primitive.occ===:cylinder
            n=_occ_dir(r.axis)
            V=π*r.radius*r.radius*r.height
            cb=_occ3_add(r.center,_occ3_scale(n,r.height/2))
            qp=V*r.radius*r.radius/4; qa=V*r.height*r.height/12
            for i in 1:3,j in 1:3
                m[i,j]=qp*(Float64(i==j)-n[i]*n[j])+qa*n[i]*n[j]+
                    V*cb[i]*cb[j]
            end
            return m
        elseif g.primitive.occ===:sphere
            V=4π/3*r.radius^3; c=r.center
            q=V*r.radius*r.radius/5
            for i in 1:3,j in 1:3
                m[i,j]=q*Float64(i==j)+V*c[i]*c[j]
            end
            return m
        elseif g.primitive.occ===:cone
            n=_occ_dir(r.axis)
            q2=r.r1*r.r1+r.r1*r.r2+r.r2*r.r2
            V=π*r.height*q2/3
            zb=r.height*(r.r1*r.r1+2*r.r1*r.r2+3*r.r2*r.r2)/(4q2)
            cb=_occ3_add(r.center,_occ3_scale(n,zb))
            S4=r.r1^4+r.r1^3*r.r2+r.r1^2*r.r2^2+r.r1*r.r2^3+r.r2^4
            qp=π*r.height*S4/20
            qa=π*r.height^3*(r.r1^2+3*r.r1*r.r2+6*r.r2^2)/30-V*zb*zb
            for i in 1:3,j in 1:3
                m[i,j]=qp*(Float64(i==j)-n[i]*n[j])+qa*n[i]*n[j]+
                    V*cb[i]*cb[j]
            end
            return m
        end
        throw(ErrorException("OCC.get_matrix_of_inertia: unsupported " *
            "primitive kind $(g.primitive.occ)"))
    end
    uses=Tuple{UInt64,Bool}[]
    for shell_id in shape.children
        for (face_id,rev) in kernel.shapes[shell_id].geometry.uses
            push!(uses,(face_id,rev))
        end
    end
    isempty(uses) && return m
    face_kind=u->get(kernel.shapes[u[1]].geometry.geometry,:occ,nothing)
    all(u->face_kind(u)===:plane,uses) || throw(ErrorException(
        "OCC.get_matrix_of_inertia: unsupported solid boundary"))
    for (fid,rev) in uses
        rec=kernel.shapes[fid].geometry.geometry
        σ=rev ? -1.0 : 1.0
        n=(σ*rec.axis[1],σ*rec.axis[2],σ*rec.axis[3])
        U=_model_occ_plane_moment3(kernel,fid)
        for i in 1:3,j in 1:3,k in 1:3
            m[i,j]+=n[k]*U[i,j,k]/5
        end
    end
    return m
end

function _model_occ_get_matrix_of_inertia(m,dim,tag)
    d=_model_occ_dimension(dim,"OCC.get_matrix_of_inertia")
    t=_model_occ_requested_tag(tag,"OCC.get_matrix_of_inertia")
    kernel=m.occ_kernel
    id=get(kernel.bindings,(d,t),nothing)
    id===nothing && throw(ArgumentError(
        "OCC.get_matrix_of_inertia: Unknown OpenCASCADE entity of " *
        "dimension $d with tag $t"))
    Q=if d==1
        _model_occ_curve_moment2(kernel,id)
    elseif d==2
        _model_occ_face_moment2(kernel,id)
    elseif d==3
        _model_occ_solid_moment2(kernel,id)
    else
        zeros(3,3)
    end
    tr=Q[1,1]+Q[2,2]+Q[3,3]
    return Float64[tr-Q[1,1],-Q[1,2],-Q[1,3],
                   -Q[2,1],tr-Q[2,2],-Q[2,3],
                   -Q[3,1],-Q[3,2],tr-Q[3,3]]
end
