# `gmsh.model.occ` BRepGProp measures — `getMass` and `getCenterOfMass`
# (LinearProperties / SurfaceProperties / VolumeProperties) on the bound
# kernel graph, matching upstream's unsynchronized registry queries.
# Curves measure per-record: line/circle in closed form; ellipse, NURBS
# (`copied_curve`) and `surface_iso` by Gauss-Legendre quadrature;
# `degenerate` collapses to (0, point). Faces measure per-record: the
# plane's pcurve shoelace with conic-segment corrections and the
# outer-minus-holes convention, closed-form `∫P·dA` over canonical
# cylinder/cone/sphere/torus domains, and tensor or trimmed-chart
# quadrature for NURBS patches. Solids use compact primitive records,
# torus-sector Pappus integrals, plane/cylinder boundary moments, or
# divergence quadrature on finite natural domains. Natural-solid moments
# retain a nearby reference frame through their centroidal integration.
# Entities outside dims 1:3 follow upstream: the switch
# never fires, so the empty GProp reports mass 0 and centre (0,0,0).

# `GProp`'s `EPS_DIM`: a length, area or volume below it is no measure, and
# its centre is reported at the origin. The measures used
# `Precision::Confusion` (1e-7, a length) for this, so a box with 1e-6
# edges — faces of 1e-12, volume 1e-18 — had "zero area". Gmsh 4.15.2
# reports that box's mass, centre and face areas, and only below 1e-20
# (edges of 2e-7) puts the centre at the origin.
const _OCC_EPS_DIM=1e-20

# 8-point Gauss-Legendre abscissae/weights on [-1,1] — exact for
# polynomials through degree 15. Natural NURBS faces select higher orders
# and bounded refinement below according to their degree and weights.
const _OCC_GL_T=(0.9602898564975363,0.7966664774136267,
                 0.5255324099163290,0.1834346424956498)
const _OCC_GL_W=(0.1012285362903763,0.2223810344533745,
                 0.3137066458778873,0.3626837833783620)

@inline _occ2_cross(a,b)=a[1]*b[2]-a[2]*b[1]
@inline _occ2_sub(a,b)=(a[1]-b[1],a[2]-b[2])
@inline _occ3_add(a,b)=(a[1]+b[1],a[2]+b[2],a[3]+b[3])
@inline _occ3_sub(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _occ3_scale(a,s)=(a[1]*s,a[2]*s,a[3]*s)

# Distinct interior breakpoints of a NURBS curve's knot range — the
# natural smoothness panels for the measure quadrature.
function _model_occ_curve_panels(curve,t0,t1)
    cuts=Float64[]
    for k in curve.curve.knots
        t0<k<t1 && (isempty(cuts) || k>cuts[end]) && push!(cuts,k)
    end
    return cuts
end

# A Boolean section curve is the composition S(pc(t)): its smoothness breaks
# sit at the trace pcurve's own knots, which the seed polyline's knot vector
# need not contain. Panels cut only at the polyline knots leave derivative
# kinks inside a panel and the adaptive rule never reaches its tolerance.
function _model_occ_curve_panels(curve::_ModelOCCSectionCurve,t0,t1)
    cuts=Float64[k for k in curve.curve.knots if t0<k<t1]
    pc=curve.pc
    sign=1.0
    while hasproperty(pc,:neg2d)
        pc=pc.neg2d[1];sign=-sign
    end
    hasproperty(pc,:spl2d) &&
        append!(cuts,Float64[sign*k for k in pc.spl2d[1] if t0<sign*k<t1])
    sort!(cuts)
    out=Float64[]
    for k in cuts
        (isempty(out) || k>out[end]) && push!(out,k)
    end
    return out
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

include("ModelOCCCurveMeasures.jl")

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
    stats=_model_occ_curve_stats(rec)
    return stats.area,_occ3_add(stats.origin,stats.centroid)
end

# ---- plane-face shoelace ---------------------------------------------

# One wire's raw signed shoelace area and first moments in the plane's
# physical uv frame. Reversed uses traverse the `fwd` pcurve backward;
# conic edges add their signed sector-minus-chord segment.
function _model_occ_wire_shoelace(kernel,face_rec,uses)
    area=0.0; Mx=0.0; My=0.0
    first_pt=nothing; last_pt=nothing
    starts=NTuple{2,Float64}[]; ends=NTuple{2,Float64}[]
    for (eid,rev) in uses
        etag=kernel.shape_tags[eid][2]
        eg=kernel.shapes[eid].geometry.geometry
        pc=face_rec.pcurves[etag].fwd
        t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
        if haskey(pc,:neg2d)
            pc=pc.neg2d[1];t0,t1=-t0,-t1
        end
        pt_a=_occ_pcurve_d0(pc,t0); pt_b=_occ_pcurve_d0(pc,t1)
        first_pt===nothing && (first_pt=pt_a)
        last_pt=pt_b
        push!(starts,pt_a); push!(ends,pt_b)
        if haskey(pc,:lin2d)
            a=_occ_pcurve_d0(pc,t0); b=_occ_pcurve_d0(pc,t1)
            cr=_occ2_cross(a,b)
            area+=cr/2
            Mx+=(a[1]+b[1])*cr/6
            My+=(a[2]+b[2])*cr/6
            continue
        end
        if !(haskey(pc,:circ2d) || haskey(pc,:ellips2d))
            # Non-elementary trace (spl2d): shoelace moments by Green's
            # theorem in the SAME symmetric convention as the elementary
            # branches — A=∮(u dv−v du)/2, ∬u=∮(u²dv−uv du)/3,
            # ∬v=∮(uv dv−v² du)/3 — GL8 per B-spline span. A mixed
            # convention (say ∮u dv here) leaves each such edge's boundary
            # term (u1v1−u0v0)/2 uncancelled inside a closed wire.
            # panels end at the knots in either direction of traversal: a
            # wire running the edge backwards (t0 > t1) got one panel
            # across every knot, which left its area good to five digits
            cuts=_model_occ_pcurve_measure_panels(pc,t0,t1)
            for k in 1:length(cuts)-1
                a0=cuts[k]; h=(cuts[k+1]-a0)/2; mid=a0+h
                for i in 1:4,s in (-1,1)
                    t=mid+s*h*_OCC_GL_T[i]; w=_OCC_GL_W[i]*h
                    (u,v),(du,dv)=_occ_pcurve_d1(pc,t)
                    area+=w*(u*dv-v*du)/2
                    Mx+=w*(u*u*dv-u*v*du)/3
                    My+=w*(u*v*dv-v*v*du)/3
                end
            end
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
    # Periodic surfaces: a wire wrapping the period closes only modulo
    # 2π, so the cut-uv boundary includes the chord last→first — add its
    # shoelace term or seam-wrapping domains lose the closure wedge.
    # A wire whose edges are stored out of traversal order (OpenCASCADE's
    # partial-cylinder side faces) still closes: every start point is some
    # edge's end point, and no chord is missing.
    closed=!isempty(starts) && all(p->any(q->abs(q[1]-p[1])<=1e-9 &&
        abs(q[2]-p[2])<=1e-9,ends),starts)
    if first_pt!==nothing && !closed &&
            (first_pt[1]!=last_pt[1] || first_pt[2]!=last_pt[2])
        cr=_occ2_cross(last_pt,first_pt)
        area+=cr/2
        Mx+=(last_pt[1]+first_pt[1])*cr/6
        My+=(last_pt[2]+first_pt[2])*cr/6
    end
    return area,Mx,My
end

# uv→physical area factor |X×Y| and unit normal of a plane record. A
# transform maps the stored frame vectors without renormalizing them — a
# dilation stretches `X`/`Y` (and `axis`) while the pcurves keep their uv
# coordinates — so a uv shoelace area is not the physical area and `axis`
# is not a unit normal. Orthonormal frames return `(1.0, axis)` verbatim.
function _occ_plane_metric(rec)
    X=(Float64(rec.X[1]),Float64(rec.X[2]),Float64(rec.X[3]))
    Y=(Float64(rec.Y[1]),Float64(rec.Y[2]),Float64(rec.Y[3]))
    c=_occ_cross(X,Y)
    J=_occ_modulus(c)
    (abs(J-1.0)<=1e-14 && abs(_dot3(rec.axis,rec.axis)-1.0)<=1e-14) &&
        return 1.0,rec.axis
    (isfinite(J) && J>0) || throw(ErrorException(
        "OCC.get_mass: plane face has a degenerate frame"))
    s=_dot3(c,rec.axis)>=0 ? 1.0/J : -1.0/J
    return J,(c[1]*s,c[2]*s,c[3]*s)
end

@inline _model_occ_plane_eval(rec,u,v)=
    (fma(u,rec.X[1],v*rec.Y[1])+rec.center[1],
     fma(u,rec.X[2],v*rec.Y[2])+rec.center[2],
     fma(u,rec.X[3],v*rec.Y[3])+rec.center[3])

# Green's-theorem planar moments through the actual pcurves (including conic
# arcs and spline traces). Raw per-wire integrals stay in traversal
# orientation; the aggregator normalizes outer/hole signs exactly like
# `_model_occ_wire_shoelace`.
const _OCC_PLANE_MOMENT_KEYS=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2),
    (3,0),(2,1),(1,2),(0,3))

function _model_occ_plane_wire_moments(kernel,rec,uses)
    m=Dict{Tuple{Int,Int},Float64}(key=>0.0 for key in _OCC_PLANE_MOMENT_KEYS)
    for (eid,rev) in uses
        etag=kernel.shape_tags[eid][2]
        eg=kernel.shapes[eid].geometry.geometry
        pc=rec.pcurves[etag].fwd
        t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
        if haskey(pc,:neg2d)
            pc=pc.neg2d[1];t0,t1=-t0,-t1
        end
        cuts=if haskey(pc,:lin2d)
            Float64[t0,t1]
        else
            base=_model_occ_pcurve_measure_panels(pc,t0,t1)
            refined=Float64[base[1]]
            for i in 1:length(base)-1
                a,b=base[i],base[i+1]
                for j in 1:8
                    push!(refined,a+(b-a)*j/8)
                end
            end
            refined
        end
        for k in 1:length(cuts)-1
            a=cuts[k];h=(cuts[k+1]-a)/2;mid=a+h
            for i in 1:4,s in (-1,1)
                t=mid+s*h*_OCC_GL_T[i];w=_OCC_GL_W[i]*h
                (u,v),(du,dv)=_occ_pcurve_d1(pc,t)
                m[(0,0)]+=w*u*dv
                m[(1,0)]+=w*(u*u/2)*dv
                m[(0,1)]+=w*(-v*v/2)*du
                m[(2,0)]+=w*(u^3/3)*dv
                m[(1,1)]+=w*(u*u*v/2)*dv
                m[(0,2)]+=w*(-v^3/3)*du
                m[(3,0)]+=w*(u^4/4)*dv
                m[(2,1)]+=w*(u^3*v/3)*dv
                m[(1,2)]+=w*(u*u*v*v/2)*dv
                m[(0,3)]+=w*(-v^4/4)*du
            end
        end
    end
    return m
end

function _model_occ_plane_boundary_moments(kernel,rec,wires)
    A=0.0
    total=Dict{Tuple{Int,Int},Float64}(key=>0.0 for key in _OCC_PLANE_MOMENT_KEYS)
    for (index,uses) in enumerate(wires)
        raw=_model_occ_plane_wire_moments(kernel,rec,uses)
        area=raw[(0,0)]
        abs(area)<=_OCC_EPS_DIM && continue
        sgn=area>=0 ? 1.0 : -1.0
        weight=(index==1 ? 1.0 : -1.0)*sgn
        A+=(index==1 ? 1.0 : -1.0)*abs(area)
        for key in _OCC_PLANE_MOMENT_KEYS
            total[key]+=weight*raw[key]
        end
    end
    return A,total
end

function _model_occ_plane_curved(kernel,shape,rec)
    for wire_id in shape.children,
            (eid,rev) in kernel.shapes[wire_id].geometry.uses
        pc=rec.pcurves[kernel.shape_tags[eid][2]].fwd
        haskey(pc,:neg2d) && (pc=pc.neg2d[1])
        haskey(pc,:lin2d) || return true
    end
    return false
end

# Expand ∬(base + du·X + dv·Y)_i(...)_j dA from the centered 2D moments.
function _model_occ_plane_tensor(base,X,Y,m)
    T=zeros(3,3)
    for i in 1:3,j in 1:3
        first=((base[i],0,0),(X[i],1,0),(Y[i],0,1))
        second=((base[j],0,0),(X[j],1,0),(Y[j],0,1))
        for (ci,pu,qu) in first,(cj,pv,qv) in second
            T[i,j]+=ci*cj*m[(pu+pv,qu+qv)]
        end
    end
    return T
end

function _model_occ_plane_cubic(base,X,Y,m)
    U=zeros(3,3,3)
    factors=((base[1],0,0),(X[1],1,0),(Y[1],0,1))
    for i in 1:3,j in 1:3,k in 1:3
        fi=((base[i],0,0),(X[i],1,0),(Y[i],0,1))
        fj=((base[j],0,0),(X[j],1,0),(Y[j],0,1))
        fk=((base[k],0,0),(X[k],1,0),(Y[k],0,1))
        for (ci,pu,qu) in fi,(cj,pv,qv) in fj,(ck,pw,qw) in fk
            U[i,j,k]+=ci*cj*ck*m[(pu+pv+pw,qu+qv+qw)]
        end
    end
    return U
end

function _model_occ_centered_plane_moments(raw,cu,cv)
    m=Dict{Tuple{Int,Int},Float64}()
    A=raw[(0,0)]
    m[(0,0)]=A;m[(1,0)]=0.0;m[(0,1)]=0.0
    m[(2,0)]=raw[(2,0)]-2cu*raw[(1,0)]+cu*cu*A
    m[(1,1)]=raw[(1,1)]-cu*raw[(0,1)]-cv*raw[(1,0)]+cu*cv*A
    m[(0,2)]=raw[(0,2)]-2cv*raw[(0,1)]+cv*cv*A
    m[(3,0)]=raw[(3,0)]-3cu*raw[(2,0)]+3cu*cu*raw[(1,0)]-cu^3*A
    m[(2,1)]=raw[(2,1)]-cv*raw[(2,0)]-2cu*raw[(1,1)]+
        2cu*cv*raw[(1,0)]+cu*cu*raw[(0,1)]-cu*cu*cv*A
    m[(1,2)]=raw[(1,2)]-cu*raw[(0,2)]-2cv*raw[(1,1)]+
        2cu*cv*raw[(0,1)]+cv*cv*raw[(1,0)]-cu*cv*cv*A
    m[(0,3)]=raw[(0,3)]-3cv*raw[(0,2)]+3cv*cv*raw[(0,1)]-cv^3*A
    return m
end

# Plane-face area, 3-D centroid and second-moment tensor T_F = ∬_F r⊗r
# dA — the tensor via the signed fan triangulation from the face's own
# shoelace centroid (exact for any line-bounded polygon, holes included;
# the midpoint rule integrates quadratics exactly). `centered=true`
# evaluates ∬_F(r-c)⊗(r-c)dA. Conic/spline wires use boundary 1-forms.
function _model_occ_plane_stats(kernel,id;want_tensor=false,
                                caller="OCC.get_center_of_mass",
                                face_rec=nothing,centered=false)
    shape=kernel.shapes[id]
    rec=face_rec===nothing ? shape.geometry.geometry : face_rec
    J,_=_occ_plane_metric(rec)
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
    A>_OCC_EPS_DIM || throw(ErrorException(
        "OCC.get_mass: plane face has zero area"))
    cu=cv=0.0
    for (index,wire) in enumerate(raw_per_wire)
        area=abs(wire.raw[1])
        area>_OCC_EPS_DIM || continue
        w=(index==1 ? area : -area)/A
        # the wire's own centroid is moment over area with their common
        # winding sign — dividing by |area| negated it for a wire wound
        # clockwise (a mirrored face's circular hole)
        cu+=w*wire.raw[2]/wire.raw[1]; cv+=w*wire.raw[3]/wire.raw[1]
    end
    c3=_model_occ_plane_eval(rec,cu,cv)
    want_tensor || return A*J,c3,nothing
    # Conic and spline traces cannot use the straight-chord fan; integrate
    # the actual boundary 1-forms instead (this also covers disk faces).
    curved=_model_occ_plane_curved(kernel,shape,rec)
    if curved
        wires=[wire.uses for wire in raw_per_wire]
        _,raw=_model_occ_plane_boundary_moments(kernel,rec,wires)
        raw[(0,0)]>_OCC_EPS_DIM || throw(ErrorException(
            "OCC.get_mass: plane face has zero area"))
        cu_b=raw[(1,0)]/raw[(0,0)];cv_b=raw[(0,1)]/raw[(0,0)]
        if centered
            mc=_model_occ_centered_plane_moments(raw,cu_b,cv_b)
            T=_model_occ_plane_tensor((0.0,0.0,0.0),rec.X,rec.Y,mc)
        else
            T=_model_occ_plane_tensor(_model_occ_plane_eval(rec,cu_b,cv_b),
                rec.X,rec.Y,raw)
        end
        return A*J,c3,T.*J
    end
    T=zeros(3,3)
    # Pass 2: per-wire signed fan from (cu,cv); each wire's contribution
    # is its raw (winding-signed) tensor scaled by sign(raw area). `centered`
    # evaluates ∬(r-c3)⊗(r-c3)dA directly, avoiding cancellation for faces far
    # from the origin.
    origin=centered ? c3 : (0.0,0.0,0.0)
    for (index,wire) in enumerate(raw_per_wire)
        sgn=wire.raw[1]>=0 ? 1.0 : -1.0
        w=(index==1 ? 1.0 : -1.0)*sgn
        c2=(cu,cv); a3=_occ3_sub(c3,origin)
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
            b3=_occ3_sub(_model_occ_plane_eval(rec,p1[1],p1[2]),origin)
            e3=_occ3_sub(_model_occ_plane_eval(rec,p2[1],p2[2]),origin)
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
    return A*J,c3,T.*J
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
function _model_occ_pcurve_measure_panels(pc,t0,t1)
    panels=max(1,ceil(Int,abs(t1-t0)/(π/2)))
    cuts=Float64[t0+(t1-t0)*i/panels for i in 0:panels]
    inner,neg=_brep_pc_unneg(pc)
    if hasproperty(inner,:spl2d)
        lo,hi=minmax(t0,t1)
        for knot in inner.spl2d[1]
            t=neg ? -knot : knot
            lo<t<hi && push!(cuts,t)
        end
    elseif hasproperty(inner,:brep_conversion2d)
        # the trace of an edge converted to a rational curve runs through
        # that curve's span-wise reparametrization, smooth only between its
        # breakpoints: a NURBS-converted disk measured one panel across
        # them was 2e-4 off
        input=inner.brep_conversion2d[2]
        if input isa _ModelOCCConicReparam
            lo,hi=minmax(t0,t1)
            for knot in input.knots
                t=neg ? -knot : knot
                lo<t<hi && push!(cuts,t)
            end
        end
    end
    sort!(unique!(cuts))
    t1<t0 && reverse!(cuts)
    return cuts
end

function _model_occ_cyl_face_div(kernel,fid;want_vmoment=false)
    shape=kernel.shapes[fid]; rec=shape.geometry.geometry
    I1=Ic=Is=Iv=0.0
    for (wi,wire_id) in enumerate(shape.children)
        uses=kernel.shapes[wire_id].geometry.uses
        j1=jc=js=jv=0.0
        u_prev=nothing; uv_first=nothing; uv_last=nothing
        for (eid,rev) in uses
            etag=kernel.shape_tags[eid][2]
            eg=kernel.shapes[eid].geometry.geometry
            entry=rec.pcurves[etag]
            pc=rev && entry.rev!==nothing ? entry.rev : entry.fwd
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            # u is periodic: a stored pcurve picks an arbitrary
            # representative mod 2π — unwrap the wire so each edge
            # continues the previous one's u, or ∮u dv telescopes to
            # zero instead of the domain area
            us,vs0=_occ_pcurve_d0(pc,t0)
            shift=u_prev===nothing ? 0.0 :
                _OCC_TWO_PI*round((u_prev-us)/_OCC_TWO_PI)
            uv_first===nothing && (uv_first=(us+shift,vs0))
            # A fitted spline is smooth within each knot span; a single
            # angular panel can cross many derivative discontinuities.
            cuts=_model_occ_pcurve_measure_panels(pc,t0,t1)
            for k in 1:length(cuts)-1
                a=cuts[k];h=(cuts[k+1]-a)/2;mid=a+h
                for i in 1:4, s in (-1,1)
                    t=mid+s*h*_OCC_GL_T[i]; w=_OCC_GL_W[i]*h
                    (u,v),(du,dv)=_occ_pcurve_d1(pc,t)
                    u+=shift
                    j1+=w*u*dv; jc+=w*sin(u)*dv; js+=w*(-cos(u)*dv)
                    want_vmoment && (jv+=w*(-(v*v/2))*du)
                end
            end
            ue,ve=_occ_pcurve_d0(pc,t1)
            uv_last=(ue+shift,ve)
            u_prev=ue+shift
        end
        # Seam-wrapping wire: the uv loop closes through the periodic
        # identification — append the last→first chord so Green's
        # theorem sees the complete contour (degenerates to zero for
        # wires closed exactly in uv).
        if uv_first!==nothing &&
                (uv_first[1]!=uv_last[1] || uv_first[2]!=uv_last[2])
            du0=uv_first[1]-uv_last[1]; dv0=uv_first[2]-uv_last[2]
            um=(uv_last[1]+uv_first[1])/2
            j1+=um*dv0
            mdu=du0==0.0 ? 0.0 : dv0/du0
            jc+=mdu*(cos(uv_last[1])-cos(uv_first[1]))
            js+=-mdu*(sin(uv_first[1])-sin(uv_last[1]))
            va,vb=uv_last[2],uv_first[2]
            jv+=-du0*(va*va+va*vb+vb*vb)/6
        end
        # Orientation from the wire's own period-unwrapped lift: ∮u dv over
        # the chord-closed contour is its signed uv area. The plane shoelace
        # (`_model_occ_wire_shoelace`) reads un-unwrapped `.fwd` pcurves, so
        # on a seam rectangle whose lower rim is not at v=0 — a wall clipped
        # between two planes — it returned the wrong sign and the patch
        # entered the divergence sum negated.
        w=(wi==1 ? 1.0 : -1.0)*(j1>=0 ? 1.0 : -1.0)
        I1+=w*j1; Ic+=w*jc; Is+=w*js; want_vmoment && (Iv+=w*jv)
    end
    return I1,Ic,Is,Iv
end

# Complete Green's-theorem domain integrals for a trimmed cylinder patch,
# used by the solid first/second boundary moments. Each integrand reduces to
# ∮Pdu+Qdv over the stored pcurves; the artificial seam-closing chord is
# integrated numerically with many GL4 panels (smooth polynomial-trig
# integrands) while each real pcurve edge keeps the same knot-panel resolve
# as `_model_occ_cyl_face_div`.
const _OCC_CYL_INTEGRALS=(:I1,:Ic,:Is,:Icc,:Ics,:Iss,:Iccc,:Iccs,:Icss,
    :Isss,:Iv,:Ivc,:Ivs,:Iv2,:Ivc2,:Ivcs,:Ivs2,:Ivvc,:Ivvs)

function _occ_cyl_chord_integrals!(acc::AbstractDict,a,b)
    du=b[1]-a[1];dv=b[2]-a[2]
    (du==0.0 && dv==0.0) && return nothing
    panels=256
    for panel in 0:panels-1
        mid=(panel+0.5)/panels;h=0.5/panels
        for i in 1:4,s in (-1,1)
            t=mid+s*h*_OCC_GL_T[i];w=_OCC_GL_W[i]*h
            u=a[1]+t*du;v=a[2]+t*dv
            su,cu=_gm_sincos(u)
            v2=v*v
            acc[:I1]+=w*u*dv
            acc[:Ic]+=w*su*dv
            acc[:Is]+=w*(-cu)*dv
            acc[:Icc]+=w*((u+su*cu)/2)*dv
            acc[:Ics]+=w*(-cu*cu/2)*dv
            acc[:Iss]+=w*((u-su*cu)/2)*dv
            acc[:Iccc]+=w*(su-su^3/3)*dv
            acc[:Iccs]+=w*(-cu^3/3)*dv
            acc[:Icss]+=w*(su^3/3)*dv
            acc[:Isss]+=w*(-cu+cu^3/3)*dv
            acc[:Iv]+=w*(-v2/2)*du
            acc[:Ivc]+=w*(-v2/2*cu)*du
            acc[:Ivs]+=w*(-v2/2*su)*du
            acc[:Iv2]+=w*(u*v2)*dv
            acc[:Ivc2]+=w*(-v2/2*cu*cu)*du
            acc[:Ivcs]+=w*(-v2/2*cu*su)*du
            acc[:Ivs2]+=w*(-v2/2*su*su)*du
            acc[:Ivvc]+=w*(v2*su)*dv
            acc[:Ivvs]+=w*(-v2*cu)*dv
        end
    end
    return nothing
end

function _model_occ_cyl_face_integrals(kernel,fid)
    shape=kernel.shapes[fid]; rec=_model_occ_analytic_surface_record(shape.geometry.geometry)
    totals=Dict(name=>0.0 for name in _OCC_CYL_INTEGRALS)
    for (wi,wire_id) in enumerate(shape.children)
        uses=kernel.shapes[wire_id].geometry.uses
        local_acc=Dict(name=>0.0 for name in _OCC_CYL_INTEGRALS)
        u_prev=nothing;uv_first=nothing;uv_last=nothing
        for (eid,rev) in uses
            etag=kernel.shape_tags[eid][2]
            eg=kernel.shapes[eid].geometry.geometry
            entry=rec.pcurves[etag]
            pc=rev && entry.rev!==nothing ? entry.rev : entry.fwd
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            us,vs0=_occ_pcurve_d0(pc,t0)
            shift=u_prev===nothing ? 0.0 :
                _OCC_TWO_PI*round((u_prev-us)/_OCC_TWO_PI)
            uv_first===nothing && (uv_first=(us+shift,vs0))
            cuts=_model_occ_pcurve_measure_panels(pc,t0,t1)
            for k in 1:length(cuts)-1
                a=cuts[k];h=(cuts[k+1]-a)/2;mid=a+h
                for i in 1:4, s in (-1,1)
                    t=mid+s*h*_OCC_GL_T[i];w=_OCC_GL_W[i]*h
                    (u,v),(du,dv)=_occ_pcurve_d1(pc,t)
                    u+=shift
                    su,cu=_gm_sincos(u)
                    v2=v*v
                    local_acc[:I1]+=w*u*dv
                    local_acc[:Ic]+=w*su*dv
                    local_acc[:Is]+=w*(-cu)*dv
                    local_acc[:Icc]+=w*((u+su*cu)/2)*dv
                    local_acc[:Ics]+=w*(-cu*cu/2)*dv
                    local_acc[:Iss]+=w*((u-su*cu)/2)*dv
                    local_acc[:Iccc]+=w*(su-su^3/3)*dv
                    local_acc[:Iccs]+=w*(-cu^3/3)*dv
                    local_acc[:Icss]+=w*(su^3/3)*dv
                    local_acc[:Isss]+=w*(-cu+cu^3/3)*dv
                    local_acc[:Iv]+=w*(-v2/2)*du
                    local_acc[:Ivc]+=w*(-v2/2*cu)*du
                    local_acc[:Ivs]+=w*(-v2/2*su)*du
                    local_acc[:Iv2]+=w*(u*v2)*dv
                    local_acc[:Ivc2]+=w*(-v2/2*cu*cu)*du
                    local_acc[:Ivcs]+=w*(-v2/2*cu*su)*du
                    local_acc[:Ivs2]+=w*(-v2/2*su*su)*du
                    local_acc[:Ivvc]+=w*(v2*su)*dv
                    local_acc[:Ivvs]+=w*(-v2*cu)*dv
                end
            end
            ue,ve=_occ_pcurve_d0(pc,t1)
            uv_last=(ue+shift,ve)
            u_prev=ue+shift
        end
        if uv_first!==nothing &&
                (uv_first[1]!=uv_last[1] || uv_first[2]!=uv_last[2])
            _occ_cyl_chord_integrals!(local_acc,uv_last,uv_first)
        end
        # signed uv area of the unwrapped contour — see `_model_occ_cyl_face_div`
        w=(wi==1 ? 1.0 : -1.0)*(local_acc[:I1]>=0 ? 1.0 : -1.0)
        for (index,name) in enumerate(_OCC_CYL_INTEGRALS)
            totals[name]+=w*local_acc[name]
        end
    end
    return NamedTuple{_OCC_CYL_INTEGRALS}(Tuple(totals[name] for name in _OCC_CYL_INTEGRALS))
end

# ∬ p (p·u_r+R) dudv, ∬ (p-offset)⊗(p-offset)((p-offset)·u_r+R) dudv and
# J0=∬(p·u_r+R)dudv for one cylinder patch. The caller multiplies by the
# frame handedness, shell orientation and radius.
function _model_occ_cyl_face_moment_terms(kernel,fid,offset)
    J=_model_occ_cyl_face_integrals(kernel,fid)
    rec=_model_occ_analytic_surface_record(kernel.shapes[fid].geometry.geometry)
    R=rec.radius;X=rec.X;Y=rec.Y;Z=rec.axis
    C=_occ3_sub(rec.center,offset)
    Cx=_dot3(C,X);Cy=_dot3(C,Y)
    J0=Cx*J.Ic+Cy*J.Is+R*J.I1
    Wx=Cx*J.Icc+Cy*J.Ics+R*J.Ic
    Wy=Cx*J.Ics+Cy*J.Iss+R*J.Is
    W=_occ3_add(_occ3_scale(X,Wx),_occ3_scale(Y,Wy))
    K=Cx*J.Ivc+Cy*J.Ivs+R*J.Iv
    M=_occ3_add(_occ3_scale(C,J0),
        _occ3_add(_occ3_scale(W,R),_occ3_scale(Z,K)))
    Axx=Cx*J.Iccc+Cy*J.Iccs+R*J.Icc
    Axy=Cx*J.Iccs+Cy*J.Icss+R*J.Ics
    Ayy=Cx*J.Icss+Cy*J.Isss+R*J.Iss
    Lx=Cx*J.Ivc2+Cy*J.Ivcs+R*J.Ivc
    Ly=Cx*J.Ivcs+Cy*J.Ivs2+R*J.Ivs
    L=_occ3_add(_occ3_scale(X,Lx),_occ3_scale(Y,Ly))
    K2=Cx*J.Ivvc+Cy*J.Ivvs+R*J.Iv2
    Q=zeros(3,3)
    for i in 1:3,j in 1:3
        Q[i,j]+=C[i]*C[j]*J0
        Q[i,j]+=R*(C[i]*W[j]+W[i]*C[j])
        Q[i,j]+=R*R*(X[i]*X[j]*Axx+
            (X[i]*Y[j]+Y[i]*X[j])*Axy+Y[i]*Y[j]*Ayy)
        Q[i,j]+=K*(Z[i]*C[j]+C[i]*Z[j])
        Q[i,j]+=R*(Z[i]*L[j]+L[i]*Z[j])
        Q[i,j]+=K2*Z[i]*Z[j]
    end
    return J0,M,Q
end

# ---- trimmed faces ----------------------------------------------------
#
# `addTrimmedSurface` faces carry `trimmed=true` and a boundary that is
# not the base's canonical domain, so measures integrate over the stored
# pcurves instead. Planes reuse the shoelace; cylinders extend the
# Green's-theorem divergence with the v-moment ∬v dudv = ∮(−v²/2)du —
# giving the exact patch centroid
#   ∬p dA = R·(C·I1 + R·(X·Ic+Y·Is) + Z·Iv),  A = R·I1.
# Other bases fall back to a CDT of the sampled uv boundary with a
# per-triangle Gauss rule on dA=|p_u×p_v| — the only approximation is the
# pcurve's own resolution, which the `spl2d` knot spacing keeps tight.
function _model_occ_trimmed_face_measure(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    occ=get(rec,:occ,nothing)
    if occ===:plane
        A,c3,_=_model_occ_plane_stats(kernel,id)
        return A,c3
    elseif occ===:cylinder
        I1,Ic,Is,Iv=_model_occ_cyl_face_div(kernel,id;want_vmoment=true)
        # Trim wires wind either way — normalize the moment tuple to the
        # counter-clockwise orientation before turning it into a measure.
        s=I1>=0.0 ? 1.0 : -1.0
        I1*=s; Ic*=s; Is*=s; Iv*=s
        A=rec.radius*I1
        A>_OCC_EPS_DIM || throw(ErrorException(
            "OCC.get_mass: trimmed cylinder face has zero area"))
        X,Y,Z,C=rec.X,rec.Y,rec.axis,rec.center
        M=_occ3_scale(
            _occ3_add(_occ3_scale(C,I1),
                      _occ3_add(_occ3_scale(X,rec.radius*Ic),
                                _occ3_add(_occ3_scale(Y,rec.radius*Is),
                                          _occ3_scale(Z,Iv)))),
            rec.radius)
        return A,_occ3_scale(M,1.0/A)
    end
    if occ===:sphere
        J=_model_occ_sphere_face_integrals(kernel,id)
        J===nothing || return _model_occ_sphere_face_measure(rec,J)
    end
    green=_model_occ_green_face_integrals(kernel,id)
    green===nothing || return green[1],_occ3_scale(green[2],1.0/green[1])
    # Generic parametric base: sample the pcurve boundary at each edge's
    # resolution, triangulate the uv domain, integrate |p_u×p_v| with a
    # three-point rule over 4×-split subtriangles.
    xs,ys,T,interior=_model_occ_trimmed_face_cdt(kernel,id,shape,rec)
    A=0.0; M=(0.0,0.0,0.0)
    # reference-triangle rule on (0,0),(1,0),(0,1): weights 1/3 at
    # (1/6,1/6),(2/3,1/6),(1/6,2/3), applied on 4 uniform subtriangles.
    nodes=((1/6,1/6),(2/3,1/6),(1/6,2/3))
    for t in eachindex(interior)
        interior[t] || continue
        p1=(xs[T.tv[3t-2]],ys[T.tv[3t-2]])
        p2=(xs[T.tv[3t-1]],ys[T.tv[3t-1]])
        p3=(xs[T.tv[3t]],ys[T.tv[3t]])
        m12=_occ2_mid(p1,p2); m23=_occ2_mid(p2,p3); m31=_occ2_mid(p3,p1)
        for (q1,q2,q3) in ((p1,m12,m31),(m12,p2,m23),(m31,m23,p3),
                          (m12,m23,m31))
            jac=abs((q2[1]-q1[1])*(q3[2]-q1[2])-
                    (q3[1]-q1[1])*(q2[2]-q1[2]))
            for (s,tt) in nodes
                u=q1[1]+s*(q2[1]-q1[1])+tt*(q3[1]-q1[1])
                v=q1[2]+s*(q2[2]-q1[2])+tt*(q3[2]-q1[2])
                p,du,dv=_occ_surface_d1(rec,u,v)
                dA=_occ_modulus(_occ_cross(du,dv))
                w=jac/2/3
                A+=w*dA; M=_occ3_add(M,_occ3_scale(p,w*dA))
            end
        end
    end
    A>_OCC_EPS_DIM || throw(ErrorException(
        "OCC.get_mass: trimmed face has zero area"))
    centroid=_occ3_scale(M,1.0/A)
    if occ===:cone
        # the area has an exact boundary form; the centroid keeps the rule
        exact=_model_occ_cone_face_integrals(kernel,id)
        exact===nothing || (A=exact[1])
    end
    return A,centroid
end

# ---- cone patches ----------------------------------------------------
#
# A trimmed cone face p(u,v)=C+z·a+ρ·e(u), z=v·cos γ, ρ=r1+k·z (v the slant
# length, k the slope, γ=atan k) has dA=ρ·du·dv exactly and
#   p·(p_u×p_v) = ρ·cos γ·[C·e(u) + r1 − k·(C·a)]
# for a right-handed frame, so its area and its divergence volume are the
# boundary integrals ∮−W(v)du and ∮−cos γ·W(v)·g(u)du with W=r1·v+k·cos γ·v²/2
# — exact to the boundary quadrature, where the triangulated-domain rule
# left a box cut by a cone 1e-5 off. Returns (area, ∬p·(p_u×p_v)/3) for the
# right-handed normal, or `nothing` for a wire that winds around the axis
# without a seam or lacks a pcurve.
function _model_occ_cone_face_integrals(kernel,fid)
    shape=kernel.shapes[fid]; rec=shape.geometry.geometry
    haskey(rec,:pcurves) || return nothing
    k=(rec.r2-rec.r1)/rec.height
    cγ=1/sqrt(1+k*k)
    Yc=_occ_frame_y(rec)
    Cx=_dot3(rec.center,rec.X);Cy=_dot3(rec.center,Yc)
    c0=rec.r1-k*_dot3(rec.center,rec.axis)
    area=0.0;vol=0.0
    for (wi,wire_id) in enumerate(shape.children)
        uses=kernel.shapes[wire_id].geometry.uses
        a=0.0;w3=0.0;turn=0.0
        for (eid,rev) in uses
            entry=get(rec.pcurves,kernel.shape_tags[eid][2],nothing)
            entry===nothing && return nothing
            eg=kernel.shapes[eid].geometry.geometry
            pc=rev && entry.rev!==nothing ? entry.rev : entry.fwd
            pc===nothing && return nothing
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            cuts=_model_occ_pcurve_measure_panels(pc,t0,t1)
            for j in 1:length(cuts)-1
                lo=cuts[j];h=(cuts[j+1]-lo)/2;mid=lo+h
                for i in 1:4, sgn in (-1,1)
                    t=mid+sgn*h*_OCC_GL_T[i]
                    (u,v),(du,_)=_occ_pcurve_d1(pc,t)
                    w=_OCC_GL_W[i]*h*du
                    su,cu=_gm_sincos(u)
                    W=v*(rec.r1+k*cγ*v/2)
                    turn+=w
                    a-=w*W
                    w3-=w*cγ*W*(Cx*cu+Cy*su+c0)
                end
            end
        end
        abs(turn)>1.0 && return nothing
        # outer wire positive, holes negative, whichever way each winds
        sgn=(wi==1 ? 1.0 : -1.0)*(a>=0 ? 1.0 : -1.0)
        area+=sgn*a;vol+=sgn*w3
    end
    area>0 || return nothing
    return area,vol/3
end

# ---- sphere patches --------------------------------------------------
#
# A trimmed sphere face p(u,v)=C+R·n̂, n̂=cos v(cos u·X+sin u·Y)+sin v·Z, has
# dA=R²cos v·du·dv. Every domain moment of n̂ that the area, the divergence
# volume and their first moments need reduces by Green's theorem to a
# boundary line integral ∬_D f du dv = ∮_{∂D} P du with ∂P/∂v = −f:
#   1·cos v          → −sin v
#   n̂ₓ,n̂ᵧ ·cos v     → −(cos u, sin u)·(v/2+sin 2v/4)
#   n̂_z·cos v        → −sin²v/2
#   n̂ᵢn̂ⱼ·cos v       → products of (cos u, sin u) with −(sin v−sin³v/3),
#                       cos³v/3, −sin³v/3
# so the measures are exact to the boundary quadrature — the pole's
# degenerate edge contributes through its own du — instead of a
# triangulated-domain rule. Returns `nothing` for a wire that winds around
# the axis without a seam (no planar contour to orient) or lacks a pcurve;
# the caller then keeps the triangulated rule.
const _OCC_SPH_INTEGRALS=(:I0,:Ix,:Iy,:Iz,:Txx,:Txy,:Tyy,:Txz,:Tyz,:Tzz)

_occ_sphere_conic_panels(t0,t1)=(n=max(1,ceil(Int,abs(t1-t0)/(π/16)));
    Float64[t0+(t1-t0)*i/n for i in 0:n])

# Chart data of the sphere `rec` at parameter `t` of a circle or ellipse
# edge record lying on it: (sin u, cos u, sin v, cos v, v, du/dt).
function _occ_sphere_conic_chart(rec,eg,t)
    s,c=_gm_sincos(t)
    rx,ry=eg.occ===:circle ? (eg.r,eg.r) : (eg.rx,eg.ry)
    Yc=_occ_circle_y(eg)
    p=ntuple(k->eg.center[k]+rx*c*eg.X[k]+ry*s*Yc[k]-rec.center[k],3)
    dp=ntuple(k->-rx*s*eg.X[k]+ry*c*Yc[k],3)
    Ys=_occ_frame_y(rec)
    x=_dot3(p,rec.X);y=_dot3(p,Ys);z=_dot3(p,rec.axis)
    dx=_dot3(dp,rec.X);dy=_dot3(dp,Ys)
    ρ2=x*x+y*y;ρ=sqrt(ρ2);R=sqrt(ρ2+z*z)
    return y/ρ,x/ρ,z/R,ρ/R,atan(z,ρ),(x*dy-y*dx)/ρ2
end

function _model_occ_sphere_face_integrals(kernel,fid)
    shape=kernel.shapes[fid]; rec=shape.geometry.geometry
    haskey(rec,:pcurves) || return nothing
    total=zeros(10)
    for (wi,wire_id) in enumerate(shape.children)
        uses=kernel.shapes[wire_id].geometry.uses
        acc=zeros(10); turn=0.0
        for (eid,rev) in uses
            entry=get(rec.pcurves,kernel.shape_tags[eid][2],nothing)
            entry===nothing && return nothing
            eg=kernel.shapes[eid].geometry.geometry
            pc=rev && entry.rev!==nothing ? entry.rev : entry.fwd
            pc===nothing && return nothing
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            # An exact conic edge whose chart trace is only a fit (a circle
            # that is no parallel or meridian) is integrated on the conic
            # itself: chart point and du come from the 3-D point and
            # velocity, so the measure does not inherit the fit tolerance.
            exact=hasproperty(_brep_pc_unneg(pc)[1],:spl2d) &&
                get(eg,:occ,nothing) in (:circle,:ellipse)
            cuts=exact ? _occ_sphere_conic_panels(t0,t1) :
                _model_occ_pcurve_measure_panels(pc,t0,t1)
            for k in 1:length(cuts)-1
                a=cuts[k];h=(cuts[k+1]-a)/2;mid=a+h
                for i in 1:4, sgn in (-1,1)
                    t=mid+sgn*h*_OCC_GL_T[i]
                    local w,su,cu,sv,cv,v
                    if exact
                        su,cu,sv,cv,v,du=_occ_sphere_conic_chart(rec,eg,t)
                        w=_OCC_GL_W[i]*h*du
                    else
                        (u,v),(du,_)=_occ_pcurve_d1(pc,t)
                        w=_OCC_GL_W[i]*h*du
                        su,cu=_gm_sincos(u); sv,cv=_gm_sincos(v)
                    end
                    q=v/2+sv*cv/2          # ∫cos²v dv
                    c3=sv-sv*sv*sv/3       # ∫cos³v dv
                    turn+=w
                    acc[1]+=w*(-sv)
                    acc[2]+=w*(-cu*q)
                    acc[3]+=w*(-su*q)
                    acc[4]+=w*(-sv*sv/2)
                    acc[5]+=w*(-cu*cu*c3)
                    acc[6]+=w*(-su*cu*c3)
                    acc[7]+=w*(-su*su*c3)
                    acc[8]+=w*(cu*cv*cv*cv/3)
                    acc[9]+=w*(su*cv*cv*cv/3)
                    acc[10]+=w*(-sv*sv*sv/3)
                end
            end
        end
        abs(turn)>1.0 && return nothing
        # outer wire positive, holes negative, whichever way each winds
        sgn=(wi==1 ? 1.0 : -1.0)*(acc[1]>=0 ? 1.0 : -1.0)
        total.+=sgn.*acc
    end
    total[1]>0 || return nothing
    return NamedTuple{_OCC_SPH_INTEGRALS}(Tuple(total))
end

# Area and centroid of a trimmed sphere face from its boundary integrals.
function _model_occ_sphere_face_measure(rec,J)
    R=rec.radius;X,Y,Z,C=rec.X,_occ_frame_y(rec),rec.axis,rec.center
    A=R*R*J.I0
    n1=ntuple(k->X[k]*J.Ix+Y[k]*J.Iy+Z[k]*J.Iz,3)
    return A,ntuple(k->C[k]+R*n1[k]/J.I0,3)
end

# Divergence contributions (V, M) = (∬p·n̂ dA/3, ∬p(p·n̂)dA) of a trimmed
# sphere face for the outward radial normal.
function _model_occ_sphere_face_div(rec,J)
    R=rec.radius;X,Y,Z,C=rec.X,_occ_frame_y(rec),rec.axis,rec.center
    cl=(_dot3(C,X),_dot3(C,Y),_dot3(C,Z))
    n1=(J.Ix,J.Iy,J.Iz)
    cn=cl[1]*n1[1]+cl[2]*n1[2]+cl[3]*n1[3]          # ∬C·n̂ cos v
    V=R*R*(cn+R*J.I0)/3
    tc=(J.Txx*cl[1]+J.Txy*cl[2]+J.Txz*cl[3],
        J.Txy*cl[1]+J.Tyy*cl[2]+J.Tyz*cl[3],
        J.Txz*cl[1]+J.Tyz*cl[2]+J.Tzz*cl[3])       # ∬n̂(C·n̂) cos v
    world(a)=ntuple(k->X[k]*a[1]+Y[k]*a[2]+Z[k]*a[3],3)
    wn=world(n1);wt=world(tc)
    R2=R*R
    M=ntuple(k->R2*(C[k]*cn+R*C[k]*J.I0+R*wt[k]+R2*wn[k]),3)
    return V,M
end

# CDT of the sampled pcurve boundary — shared scaffold for the trimmed-face
# measure and the divergence integral below. Each wire's uses walk `pcurves
# [etag].fwd` over their orientation's parameter range; consecutive samples
# unwrap by whole periods so periodic faces triangulate a contiguous uv
# polygon. Returns (xs, ys, T, interior) — CDT vertex coords, the
# constrained triangulation, and `classify_interior`'s inside mask.
function _model_occ_trimmed_face_cdt(kernel,id,shape,rec)
    period=_brep_uperiod(rec)
    vperiod=_brep_vperiod(rec)
    xs=Float64[]; ys=Float64[]; segs=Tuple{Int,Int}[]
    index=Dict{NTuple{2,Float64},Int}()
    emit(uv)=get!(index,uv) do
        push!(xs,uv[1]); push!(ys,uv[2]); return length(xs)
    end
    lo=[Inf,Inf]; hi=[-Inf,-Inf]
    for wire_id in shape.children
        uses=kernel.shapes[wire_id].geometry.uses
        loop_idx=Int[]; u_prev=nothing; v_prev=nothing
        for (eid,rev) in uses
            etag=kernel.shape_tags[eid][2]
            eg=kernel.shapes[eid].geometry.geometry
            ent=rec.pcurves[etag]
            pc=(rev && ent.rev!==nothing) ? ent.rev : ent.fwd
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            ruv=0.002*(isfinite(hi[1]-lo[1]) ? max(hi[1]-lo[1],1.0) : 1.0)
            step=_occ_pcurve_resolution(pc,ruv)
            ns=max(9,ceil(Int,abs(t1-t0)/max(step,eps())))
            ns=min(ns,2048)
            for i in 0:ns
                t=t0+(t1-t0)*i/ns
                u,v=_occ_pcurve_d0(pc,t)
                if period>0.0 && u_prev!==nothing
                    u+=period*round((u_prev-u)/period)
                end
                if vperiod>0.0 && v_prev!==nothing
                    v+=vperiod*round((v_prev-v)/vperiod)
                end
                u_prev=u;v_prev=v
                lo[1]=min(lo[1],u); hi[1]=max(hi[1],u)
                lo[2]=min(lo[2],v); hi[2]=max(hi[2],v)
                # Consecutive uses meet at the wire's shared vertices —
                # each edge's own pcurve evaluates that vertex ~1e-12
                # apart, which the CDT reads as a constraint crossing.
                # Weld the use's first sample onto the previous use's
                # terminal index instead (the chain certificate already
                # guarantees the two ends are one vertex).
                i==0 && !isempty(loop_idx) &&
                    (push!(loop_idx,loop_idx[end]);continue)
                push!(loop_idx,emit((u,v)))
            end
        end
        n=length(loop_idx)
        n>=3 || throw(ArgumentError(
            "OCC.get_mass: trimmed face wire needs at least three points"))
        # The wire's closing vertex welds the same way: the last use's
        # terminal sample and the first use's start are one vertex —
        # reuse the first index so the polygon's final real edge ends
        # there and the degenerate closing edge drops out below.
        scale=max(1.0,hi[1]-lo[1],hi[2]-lo[2])
        hypot(xs[loop_idx[end]]-xs[loop_idx[1]],
              ys[loop_idx[end]]-ys[loop_idx[1]])<=1e-9*scale &&
            (loop_idx[end]=loop_idx[1])
        for k in 1:n
            i,j=loop_idx[k],loop_idx[mod1(k+1,n)]
            i==j || (xs[i]==xs[j] && ys[i]==ys[j]) || push!(segs,(i,j))
        end
    end
    isempty(segs) && throw(ErrorException(
        "OCC.get_mass: trimmed face has an empty uv boundary"))
    # Interior lattice: the CDT's natural vertices are the boundary
    # samples, so interior triangles span nearly the whole domain and
    # the fixed-order per-triangle quadrature loses the surface's
    # curvature (the implied area vector stopped matching the wire at
    # ~1e-3 on curved lofts). A half-cell-offset lattice over the uv
    # bbox refines the interior cells; the parity classification is
    # unaffected by extra vertices and the points are deduped against
    # the boundary samples.
    Wu=hi[1]-lo[1]; Wv=hi[2]-lo[2]
    pitch=max(Wu,Wv)/32
    if isfinite(pitch) && pitch>0.0
        nx=clamp(ceil(Int,Wu/pitch),2,64)
        ny=clamp(ceil(Int,Wv/pitch),2,64)
        for i in 0:nx-1, j in 0:ny-1
            emit((lo[1]+Wu*(i+0.5)/nx,lo[2]+Wv*(j+0.5)/ny))
        end
    end
    T=try
        constrained_delaunay(xs,ys,segs)
    catch err
        if get(ENV,"TESS_CDT_DEBUG","")!=""
            io=IOBuffer()
            println(io,"trimmed face CDT failed: ",
                    sprint(showerror,err))
            println(io,"uperiod=$period vperiod=$vperiod")
            println(io,"bbox u=[$(lo[1]),$(hi[1])] v=[$(lo[2]),$(hi[2])] " *
                    "nv=$(length(xs)) nseg=$(length(segs))")
            for (i,(a,b)) in enumerate(segs)
                println(io,"seg$i ($a)->($b) " *
                        "($(xs[a]),$(ys[a])) -> ($(xs[b]),$(ys[b]))")
            end
            path=joinpath(pwd(),"_tcdt_dump.txt")
            open(io2->write(io2,take!(io)),path,"w")
            println(stderr,"trimmed face CDT debug dump -> $path")
        end
        rethrow()
    end
    interior=classify_interior(T)
    return xs,ys,T,interior
end

# Divergence integral ∬_F p·(p_u×p_v) over a trimmed parametric face —
# the CDT-restricted counterpart of `_model_occ_patch_face_div` for
# faces whose wire bounds a sub-region of the surface's parameter
# domain (STEP-imported patches, quadrics with pcurve wires). Returns
# (V, M) with V=∬rn/3, M=∬p·rn — the same signed contributions
# `get_mass`/`get_center_of_mass` accumulate per shell use.
function _model_occ_trimmed_face_div(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    green=_model_occ_green_face_integrals(kernel,id)
    green===nothing || return green[3],green[4]
    xs,ys,T,interior=_model_occ_trimmed_face_cdt(kernel,id,shape,rec)
    nodes=((1/6,1/6),(2/3,1/6),(1/6,2/3))
    V=0.0; M=(0.0,0.0,0.0)
    for t in eachindex(interior)
        interior[t] || continue
        p1=(xs[T.tv[3t-2]],ys[T.tv[3t-2]])
        p2=(xs[T.tv[3t-1]],ys[T.tv[3t-1]])
        p3=(xs[T.tv[3t]],ys[T.tv[3t]])
        m12=_occ2_mid(p1,p2); m23=_occ2_mid(p2,p3); m31=_occ2_mid(p3,p1)
        for (q1,q2,q3) in ((p1,m12,m31),(m12,p2,m23),(m31,m23,p3),
                          (m12,m23,m31))
            jac=abs((q2[1]-q1[1])*(q3[2]-q1[2])-
                    (q3[1]-q1[1])*(q2[2]-q1[2]))
            for (s,tt) in nodes
                u=q1[1]+s*(q2[1]-q1[1])+tt*(q3[1]-q1[1])
                v=q1[2]+s*(q2[2]-q1[2])+tt*(q3[2]-q1[2])
                p,du,dv=_occ_surface_d1(rec,u,v)
                w=jac/2/3
                rn=_dot3(p,_occ_cross(du,dv))
                V+=w*rn/3; M=_occ3_add(M,_occ3_scale(p,w*rn))
            end
        end
    end
    return V,M
end

@inline _occ2_mid(a,b)=((a[1]+b[1])/2,(a[2]+b[2])/2)

# Trimmed tensor faces by Green's theorem on the exact boundary: with
# F(u,v) = ∫ f(s,v) ds from the patch's lower u bound, ∬_D f = ∮_∂D F dv.
# F is a Gauss sum per u knot span, the boundary integral adaptive Gauss
# along each pcurve between its own knots. The triangulated domain of
# `_model_occ_trimmed_face_cdt` (a polygonal boundary) left the faces of an
# imported STEP part 1e-4 off. f collects |p_u×p_v|, p|p_u×p_v| and
# p·(p_u×p_v). Returns (area, ∬p dA, ∬p·n/3, ∬p (p·n)) for the counter-
# clockwise domain, or `nothing` where the method does not apply (a
# periodic patch, a missing pcurve).
const _OCC_GREEN_INNER=16
function _model_occ_green_face_integrals(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    (hasproperty(rec,:patch) && hasproperty(rec,:pcurves)) || return nothing
    (_brep_uperiod(rec)>0.0 || _brep_vperiod(rec)>0.0) && return nothing
    surface=rec.patch.surface
    ku=surface.knots_u; pu=surface.degree_u
    ulo=ku[pu+1]; uhi=ku[end-pu]
    ubreaks=unique!(Float64[k for k in ku if ulo<=k<=uhi])
    nodes,weights=_occ_natural_gauss_rule(_OCC_GREEN_INNER)
    function inner(u,v)
        acc=zeros(8)   # area, first(3), volume, divfirst(3)
        u==ulo && return acc
        a0=min(u,ulo); a1=max(u,ulo); sgn=u>=ulo ? 1.0 : -1.0
        cuts=Float64[a0]
        for k in ubreaks
            a0<k<a1 && push!(cuts,k)
        end
        push!(cuts,a1)
        for c in 1:length(cuts)-1
            h=(cuts[c+1]-cuts[c])/2; mid=(cuts[c+1]+cuts[c])/2
            for i in eachindex(nodes),si in (-1,1)
                s=mid+si*h*nodes[i]; w=sgn*weights[i]*h
                p,du,dv=_occ_surface_d1(rec,s,v)
                n=_occ_cross(du,dv); a=_occ_modulus(n); rn=_dot3(p,n)
                acc[1]+=w*a
                acc[2]+=w*a*p[1]; acc[3]+=w*a*p[2]; acc[4]+=w*a*p[3]
                acc[5]+=w*rn/3
                acc[6]+=w*rn*p[1]; acc[7]+=w*rn*p[2]; acc[8]+=w*rn*p[3]
            end
        end
        return acc
    end
    g16=_occ_natural_gauss_rule(16)
    # The boundary point is the edge's own 3-D point taken back onto the
    # patch (Gauss-Newton from the pcurve's estimate): a fitted pcurve is
    # only as good as its samples, the edge is exact. Its uv velocity is
    # the edge tangent in the patch's tangent frame.
    function boundary_uv(eg,pc,t)
        (u,v),(dut,dvt)=_occ_pcurve_d1(pc,t)
        kind=get(eg,:occ,nothing)
        if kind===:copied_curve && hasproperty(eg,:curve)
            C,Ct,_=_model_occ_curve_d2(eg.curve,t)
        elseif kind===:circle && !hasproperty(eg,:native_geometry) &&
                !hasproperty(eg,:native_arithmetic)
            sn,cs=_gm_sincos(t); Y=_occ_circle_y(eg)
            C=_occ_circle_point(eg,t)
            Ct=_occ3_add(_occ3_scale(eg.X,-eg.r*sn),_occ3_scale(Y,eg.r*cs))
        elseif kind===:ellipse
            sn,cs=_gm_sincos(t)
            C=_occ_ellipse_point(eg,t)
            Ct=_occ3_add(_occ3_scale(eg.X,-eg.rx*sn),_occ3_scale(eg.Y,eg.ry*cs))
        elseif kind===:line
            C=_occ3_add(eg.origin,_occ3_scale(eg.dir,t)); Ct=eg.dir
        else
            return u,v,dvt
        end
        for _ in 1:8
            p,Su,Sv=_occ_surface_d1(rec,u,v)
            r=_occ3_sub(C,p)
            a=_dot3(Su,Su); b=_dot3(Su,Sv); c=_dot3(Sv,Sv)
            det=a*c-b*b; det>0.0 || return u,v,dvt
            ru=_dot3(Su,r); rv=_dot3(Sv,r)
            du=(c*ru-b*rv)/det; dv=(a*rv-b*ru)/det
            u+=du; v+=dv
            abs(du)+abs(dv)<=1e-15*(1.0+abs(u)+abs(v)) && break
        end
        _,Su,Sv=_occ_surface_d1(rec,u,v)
        a=_dot3(Su,Su); b=_dot3(Su,Sv); c=_dot3(Sv,Sv)
        det=a*c-b*b
        det>0.0 || return u,v,dvt
        cu=_dot3(Su,Ct); cv=_dot3(Sv,Ct)
        return u,v,(a*cv-b*cu)/det
    end
    function panel(eg,pc,a,b)
        acc=zeros(8)
        h=(b-a)/2; mid=(a+b)/2
        for i in eachindex(g16[1]),si in (-1,1)
            t=mid+si*h*g16[1][i]
            u,v,dvt=boundary_uv(eg,pc,t)
            acc.+=(g16[2][i]*h*dvt).*inner(u,v)
        end
        return acc
    end
    function adapt(eg,pc,a,b,whole,depth)
        m=(a+b)/2
        left=panel(eg,pc,a,m); right=panel(eg,pc,m,b)
        both=left.+right
        scale=max(abs(both[1]),1e-300)
        (depth>=10 || abs(both[1]-whole[1])<=1e-14*max(scale,1.0)) &&
            return both
        return adapt(eg,pc,a,m,left,depth+1).+adapt(eg,pc,m,b,right,depth+1)
    end
    total=zeros(8)
    for wire_id in shape.children
        for (eid,rev) in kernel.shapes[wire_id].geometry.uses
            etag=get(kernel.shape_tags,eid,nothing)
            etag===nothing && return nothing
            ent=get(rec.pcurves,etag[2],nothing)
            ent===nothing && return nothing
            pc=(rev && ent.rev!==nothing) ? ent.rev : ent.fwd
            pc===nothing && return nothing
            eg=kernel.shapes[eid].geometry.geometry
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            cuts=_model_occ_pcurve_measure_panels(pc,t0,t1)
            for c in 1:length(cuts)-1
                a,b=cuts[c],cuts[c+1]
                total.+=adapt(eg,pc,a,b,panel(eg,pc,a,b),0)
            end
        end
    end
    s=total[1]>=0.0 ? 1.0 : -1.0
    total.*=s
    total[1]>_OCC_EPS_DIM || return nothing
    return (total[1],(total[2],total[3],total[4]),total[5],
            (total[6],total[7],total[8]))
end

# ---- faces -----------------------------------------------------------

# Does a quadric face's wire bound exactly the rectangular chart window
# [u0,u1]×[v0,v1]? Primitive faces do (four iso-line uses running corner to
# corner); a Boolean kept piece shares the carrier's record over a sub-domain
# wire — a sphere octant, a cylinder wall clipped between two planes — and
# must integrate over its pcurves instead of the carrier's whole chart.
# Records without a complete pcurve table keep the canonical reading.
function _occ_quadric_window_wire(kernel,shape,rec,u0,u1,v0,v1)
    haskey(rec,:pcurves) || return true
    for wire_id in shape.children, (eid,_) in kernel.shapes[wire_id].geometry.uses
        haskey(rec.pcurves,kernel.shape_tags[eid][2]) || return true
    end
    length(shape.children)==1 || return false
    uses=kernel.shapes[shape.children[1]].geometry.uses
    length(uses)==4 || return false
    tol=1e-9*max(1.0,abs(u1-u0),abs(v1-v0))
    corner(uv)=(abs(uv[1]-u0)<=tol || abs(uv[1]-u1)<=tol) &&
               (abs(uv[2]-v0)<=tol || abs(uv[2]-v1)<=tol)
    for (eid,rev) in uses
        ent=rec.pcurves[kernel.shape_tags[eid][2]]
        pc=(rev && ent.rev!==nothing) ? ent.rev : ent.fwd
        (pc!==nothing && hasproperty(pc,:lin2d)) || return false
        eg=kernel.shapes[eid].geometry.geometry
        (haskey(eg,:t0) && haskey(eg,:t1)) || return false
        a=_occ_pcurve_d0(pc,eg.t0); b=_occ_pcurve_d0(pc,eg.t1)
        (corner(a) && corner(b)) || return false
        (abs(a[1]-b[1])<=tol || abs(a[2]-b[2])<=tol) || return false
    end
    return true
end

# Canonical chart window of a quadric record — partial sweeps carry `angle`
# (azimuth span) and `lat1`/`lat2` (sphere latitude range).
function _occ_quadric_window(rec)
    occ=rec.occ
    occ===:cylinder && return (0.0,get(rec,:angle,_OCC_TWO_PI),0.0,rec.height)
    occ===:sphere && return (0.0,get(rec,:angle,_OCC_TWO_PI),
                             get(rec,:lat1,-π/2),get(rec,:lat2,π/2))
    occ===:cone && return (0.0,get(rec,:angle,_OCC_TWO_PI),0.0,
        sqrt(rec.height*rec.height+(rec.r1-rec.r2)*(rec.r1-rec.r2)))
    return nothing
end

function _model_occ_face_measure(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    if hasproperty(rec,:natural_bounds)
        integrals=_occ_natural_surface_integrals(rec;fields=(:area,:first))
        integrals.area>0 || throw(ArgumentError("OCC.get_mass: natural face has zero area"))
        return integrals.area,_occ3_scale(integrals.first,1/integrals.area)
    end
    occ=get(rec,:occ,nothing)
    # A tensor face reclassified from a genuine plane (bitwise-coplanarity
    # admission for the parametric mesher) carries its exact plane record —
    # integrate that rather than the trimmed-CDT approximation.
    occ!==:plane && haskey(rec,:plane_rec) &&
        return _model_occ_plane_stats(kernel,id;
            face_rec=rec.plane_rec)[1:2]
    if occ===:plane
        A,c3,_=_model_occ_plane_stats(kernel,id)
        return A,c3
    elseif haskey(rec,:patch)
        # Trimmed patch — the wire bounds a sub-region of the knot
        # domain (STEP-imported tensor faces); integrate the CDT of
        # the sampled uv boundary like the other trimmed bases.
        (get(rec,:trimmed,false) && haskey(rec,:pcurves)) &&
            return _model_occ_trimmed_face_measure(kernel,id)
        integrals=_model_occ_full_patch_integrals(rec;fields=(:area,:first))
        integrals.area>0 || throw(ArgumentError("OCC.get_mass: NURBS face has zero area"))
        return integrals.area,_occ3_scale(integrals.first,1/integrals.area)
    end
    if get(rec,:trimmed,false)
        return _model_occ_trimmed_face_measure(kernel,id)
    end
    # A torus face whose wire bounds a sub-domain of the tube chart
    # (boolean kept pieces, STEP trims) is no canonical sector — the
    # quadric closed form below would integrate the carrier's whole
    # domain. Route it through the trimmed-CDT measure like the other
    # pcurve-bounded faces.
    if occ===:torus && haskey(rec,:pcurves) &&
            _occ_torus_sector_window(kernel,id,rec)===nothing
        return _model_occ_trimmed_face_measure(kernel,id)
    end
    # The same holds for cylinder/sphere/cone pieces: only a wire tracing
    # the record's own chart window may use the closed forms.
    if occ===:cylinder || occ===:sphere || occ===:cone
        _occ_quadric_window_wire(kernel,shape,rec,
                                 _occ_quadric_window(rec)...) ||
            return _model_occ_trimmed_face_measure(kernel,id)
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
        # a Boolean can leave a genuine sector narrower than the carrier's
        # own sweep (a torus halved by a plane through its axis): the
        # closed form runs over the wire's window, not `rec.angle`
        win=haskey(rec,:pcurves) ?
            _occ_torus_sector_window(kernel,id,rec) : nothing
        win===nothing ? (0.0,rec.angle,0.0,_OCC_TWO_PI) :
            (win[1],win[2],0.0,_OCC_TWO_PI)
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

# NURBS patch face divergence integrals: `∬ r·(∂u×∂v) dudv` over the knot
# spans (GL8 per span like `_model_occ_face_measure`); returns the face's
# V-contribution `I/3` and the raw first moment `∬ p·(p·n)dudv` — the
# caller folds in the use sign and the global 1/4.
function _model_occ_patch_face_div(kernel,id;want_centroid=true)
    rec=kernel.shapes[id].geometry.geometry
    fields=want_centroid ? (:volume,:div_first) : (:volume,)
    integrals=_model_occ_full_patch_integrals(rec;fields)
    return integrals.volume,integrals.div_first
end

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

# Sector solid (V, centroid, Qc) — centroid measured in world coords and
# Qc the centroidal second-moment tensor ∭(p-c̄)⊗(p-c̄)dV.
function _occ_sector_stats(rec)
    U0,Uv,UE=_occ_sector_u(rec.angle,rec.X,rec.Y)
    V,M,Q=rec.occ===:sphere ? _occ_sector_stats_sph(rec,U0,Uv,UE) :
        _occ_sector_stats_axis(rec,U0,Uv,UE)
    c=rec.center
    abs(V)>_OCC_EPS_DIM || return V,(0.0,0.0,0.0),Q
    local_centroid=_occ3_scale(M,1.0/V)
    # Q is about the local frame origin at the geometry center c; shift the
    # centroidal tensor only by the (bounded) local centroid offset.
    for i in 1:3,j in 1:3
        Q[i,j]-=V*local_centroid[i]*local_centroid[j]
    end
    return V,_occ3_add(c,local_centroid),Q
end

# Divergence integral ∬_F p·(p_u×p_v) over a trimmed torus face — the
# boundary-integral counterpart of `_model_occ_trimmed_face_div`, the
# same Green's-theorem reduction `_model_occ_cyl_face_div` uses for
# cylinders. For p(u,v)=C+R(v)·(cos u·X+sin u·Y)+r2·sin v·Z with
# R(v)=r1+r2·cos v, the surface element is p_u×p_v =
# R·r2·(cos v·q̂ + sin v·Z) (q̂ = cos u·X+sin u·Y the radial unit
# direction — the parametrization normal is outward), so
#   p·(p_u×p_v) = r2·R·[R cos v + r2 sin²v
#                      + cos v·(Cx cos u + Cy sin u) + Cz sin v]
# with Ci = C·Xi. Every term is trig-polynomial in u, so the domain
# integral reduces to the boundary line integral ∬p·n̂ dA = ∮ I(u,v) dv,
#   I = r2·R·[ u·(R cos v + r2 sin²v + Cz sin v)
#             + cos v·(Cx sin u − Cy cos u) ].
# The stored pcurve walk unwraps u across the 2π azimuth seam (the bare
# u factor is rep-sensitive) and v across the tube seam so the closure
# chord measures the true geometric gap; a residual period-folded gap
# is a chart-identified seam, not a chord, so the closing leg only
# fires on a real sub-tolerance break. Wire orientation is normalized
# by shoelace sign exactly like the cylinder path — outer positive,
# holes negative — and the face use's σ applies the BRep sense as
# usual. Returns ∬rn/3 to match `_model_occ_trimmed_face_div`.
function _model_occ_torus_face_div(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    r1,r2=rec.r1,rec.r2
    X,Y,Z=rec.X,rec.Y,rec.axis
    C=rec.center
    Cx=_dot3(C,X); Cy=_dot3(C,Y); Cz=_dot3(C,Z)
    J=0.0
    @inline function Igv(u,v)
        su,cu=_gm_sincos(u); sv,cv=_gm_sincos(v)
        R=fma(r2,cv,r1)
        return r2*R*(u*(R*cv+r2*sv*sv+Cz*sv)+cv*(Cx*su-Cy*cu))
    end
    for (wi,wire_id) in enumerate(shape.children)
        uses=kernel.shapes[wire_id].geometry.uses
        # Orientation: classify the wire on its OWN period-unwrapped
        # lift (accumulate the same Green's shoelace ∮(u dv−v du)/2 the
        # line integral walks). `_model_occ_wire_shoelace` cannot be
        # reused — it reads un-unwrapped `.fwd` pcurves, so a seam-
        # rectangle outer wire (the doubly-periodic chart's boundary —
        # its two copies of each seam edge cancel exactly) collapses to
        # raw≈0 and leaves the sign to rounding noise.
        raw=0.0
        j=0.0; u_prev=nothing; v_prev=nothing
        uv_first=nothing; uv_last=nothing
        for (eid,rev) in uses
            etag=kernel.shape_tags[eid][2]
            eg=kernel.shapes[eid].geometry.geometry
            ent=rec.pcurves[etag]
            pc=(rev && ent.rev!==nothing) ? ent.rev : ent.fwd
            t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
            if haskey(pc,:neg2d)
                pc=pc.neg2d[1];t0,t1=-t0,-t1
            end
            # u and v are 2π-periodic: stored pcurves pick arbitrary
            # representatives mod the period — unwrap the wire so each
            # edge continues the previous one in both coordinates,
            # like the cylinder path's u-unwrap.
            us,vs0=_occ_pcurve_d0(pc,t0)
            ushift=u_prev===nothing ? 0.0 :
                _OCC_TWO_PI*round((u_prev-us)/_OCC_TWO_PI)
            vshift=v_prev===nothing ? 0.0 :
                _OCC_TWO_PI*round((v_prev-vs0)/_OCC_TWO_PI)
            uv_first===nothing &&
                (uv_first=(us+ushift,vs0+vshift))
            # panels end at the trace's knots: a fitted section trace is
            # only C1 there, and quarter-turn panels across hundreds of
            # knots left the volume of a torus piece good to five digits
            cuts=_model_occ_pcurve_measure_panels(pc,t0,t1)
            for k in 1:length(cuts)-1
                a=cuts[k]; h=(cuts[k+1]-a)/2; mid=a+h
                for i in 1:4, s in (-1,1)
                    t=mid+s*h*_OCC_GL_T[i]; w=_OCC_GL_W[i]*h
                    (u,v),(du,dv)=_occ_pcurve_d1(pc,t)
                    u+=ushift; v+=vshift
                    j+=w*Igv(u,v)*dv
                    raw+=w*(u*dv-v*du)/2
                end
            end
            ue,ve=_occ_pcurve_d0(pc,t1)
            uv_last=(ue+ushift,ve+vshift)
            u_prev=ue+ushift; v_prev=ve+vshift
        end
        # A wire closed only through the chart's periodic
        # identification unwraps to a period-sized gap — that is a
        # seam, not a chord to integrate. A residual sub-tolerance
        # break (weld mismatch between stored pcurves) still needs the
        # closing leg like the cylinder path's chord.
        if uv_first!==nothing
            du0=uv_first[1]-uv_last[1]
            dv0=uv_first[2]-uv_last[2]
            du0-=2π*round(du0/2π); dv0-=2π*round(dv0/2π)
            if du0!=0.0 || dv0!=0.0
                raw+=_occ2_cross(uv_last,uv_first)/2
                for k in 0:7
                    a=(k+0.5)/8
                    for i in 1:4, s in (-1,1)
                        t=a+s*(1/16)*_OCC_GL_T[i]
                        w=_OCC_GL_W[i]*(1/16)
                        u=fma(t,du0,uv_last[1]); v=fma(t,dv0,uv_last[2])
                        j+=w*Igv(u,v)*dv0
                    end
                end
            end
        end
        w=(wi==1 ? 1.0 : -1.0)*(raw>=0 ? 1.0 : -1.0)
        J+=w*j
    end
    return J/3
end

# The Pappus sector shortcut in `_model_occ_solid_measure` is exact only
# when the torus face is the carrier's canonical sector boundary: a
# single wire whose pcurves tile one u-window × the full tube-circle
# rectangle. Boolean kept pieces and STEP trims reuse the same carrier
# record over sub-domain wires — extra hole wires, curved section
# windows, or a v-window all fail the tiling test and integrate through
# the generic divergence path instead. Returns the sector's (umin,umax)
# window — its azimuth span is the Pappus angle and its midpoint the
# centroid's bisector — or `nothing` for a non-sector wire.
function _occ_torus_sector_window(kernel,fid,rec)
    shape=kernel.shapes[fid]
    length(shape.children)==1 || return nothing
    haskey(rec,:pcurves) || return nothing
    uses=kernel.shapes[first(shape.children)].geometry.uses
    lo_b,hi_b=_occ_surface_bounds(rec)
    vtol=1e-9*max(1.0,hi_b[1]-lo_b[1],hi_b[2]-lo_b[2])
    boxes=NTuple{2,Float64}[]
    umin=vmin=Inf; umax=vmax=-Inf
    for (eid,rev) in uses
        etag=kernel.shape_tags[eid][2]
        eg=kernel.shapes[eid].geometry.geometry
        ent=get(rec.pcurves,etag,nothing)
        ent===nothing && return nothing
        pc=(rev && ent.rev!==nothing) ? ent.rev : ent.fwd
        t0,t1=rev ? (eg.t1,eg.t0) : (eg.t0,eg.t1)
        if haskey(pc,:neg2d)
            pc=pc.neg2d[1];t0,t1=-t0,-t1
        end
        lo,hi=_occ_pcurve_fillbndbox(pc,min(t0,t1),max(t0,t1))
        push!(boxes,lo);push!(boxes,hi)
        umin=min(umin,lo[1]); vmin=min(vmin,lo[2])
        umax=max(umax,hi[1]); vmax=max(vmax,hi[2])
    end
    isfinite(umin) || return nothing
    # The sector sweeps the full tube circle: the wire's v extent must
    # cover one full v period. Pcurves can sit on any period-shifted
    # branch (seam unwrapping anchors the window wherever the profile
    # rim landed), so compare the SPAN, not the canonical anchor points.
    vspan=hi_b[2]-lo_b[2]
    vmax-vmin>=vspan-vtol || return nothing
    # ... and every use must lie on a rectangle side: a v-degenerate
    # rim spanning the u window (the shifted seam can sit at any v), or
    # a u-degenerate meridian spanning the wire's full v extent at a
    # window end. A curved section pcurve fails both degeneracy tests.
    for k in 1:2:length(boxes)
        lo,hi=boxes[k],boxes[k+1]
        if hi[2]-lo[2]<=vtol
            lo[1]<=umin+vtol && hi[1]>=umax-vtol || return nothing
        elseif hi[1]-lo[1]<=vtol
            (abs(lo[1]-umin)<=vtol || abs(lo[1]-umax)<=vtol) ||
                return nothing
            lo[2]<=vmin+vtol && hi[2]>=vmax-vtol || return nothing
        else
            return nothing
        end
    end
    return (umin,umax)
end

# The carrier record re-framed onto a sector window: X at azimuth `umin`,
# sweep `umax-umin` — the form the canonical-sector helpers expect.
function _occ_torus_window_record(rec,win)
    u0,u1=win
    (u0==0.0 && u1==rec.angle) && return rec
    s,c=_gm_sincos(u0)
    Y=_occ_frame_y(rec)
    return merge(rec,(X=ntuple(i->c*rec.X[i]+s*Y[i],3),
                      Y=ntuple(i->c*Y[i]-s*rec.X[i],3),angle=u1-u0))
end

# Divergence moments belong in a nearby coordinate frame. Integrating
# world-coordinate first moments and then dividing cancels terms of size
# |origin|²; a rounded erroneous centroid also corrupts centered inertia.
# Keep the local centroid separate until the public world-coordinate return.
function _model_occ_natural_solid_stats(kernel,uses;want_centroid=true,want_second=false)
    bounds=_MODEL_EMPTY_BOUNDS
    for (fid,_) in uses
        rec=kernel.shapes[fid].geometry.geometry
        bounds=_model_bounds_union(bounds,_occ_natural_surface_bbox(rec))
    end
    all(isfinite,bounds) || throw(ArgumentError(
        "OCC measures: natural solid has no finite bounds"))
    finite_bounds=bounds
    origin=ntuple(axis->_convex_coordinate(finite_bounds[axis],finite_bounds[axis+3],.5),3)
    identity=(1.,0.,0.,0.,1.,0.,0.,0.,1.)
    offset=.-origin
    transform=_AffineTransform(identity,offset,NTuple{16,Float64}[],
        (linear=identity,loc=offset))
    local_faces=Tuple{Any,Float64}[]
    volume=0.;first=(0.,0.,0.)
    for (fid,reversed) in uses
        rec=kernel.shapes[fid].geometry.geometry
        geometry=_transform_occ_surface_record(rec,transform,
            "OCC measures","natural face")
        sign=reversed ? -1. : 1.
        fields=(want_centroid || want_second) ? (:volume,:div_first) : (:volume,)
        integrals=_occ_natural_surface_integrals(geometry;fields)
        volume+=sign*integrals.volume
        (want_centroid || want_second) &&
            (first=_occ3_add(first,_occ3_scale(integrals.div_first,sign)))
        want_second && push!(local_faces,(geometry,sign))
    end
    # A length tolerance cannot classify a nonzero volume as empty. Small
    # representable solids retain their own centroid and nonzero inertia.
    total_first=first;total_volume=volume
    centroid=(want_centroid || want_second) && !iszero(total_volume) ?
        ntuple(axis->(total_first[axis]/total_volume)/4,3) : (0.,0.,0.)
    second=zeros(3,3)
    if want_second && !iszero(volume)
        for (geometry,sign) in local_faces
            integrals=_occ_natural_surface_integrals(geometry;
                centered=true,center=centroid,fields=(:centered_div_second,))
            for i in 1:3,j in 1:3
                second[i,j]+=sign*integrals.centered_div_second[i,j]
            end
        end
    end
    return (volume=total_volume,origin=origin,centroid=centroid,second=second)
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
    if all(use->hasproperty(kernel.shapes[use[1]].geometry.geometry,:natural_bounds),uses)
        stats=_model_occ_natural_solid_stats(kernel,uses;want_centroid)
        center=want_centroid ? _occ3_add(stats.origin,stats.centroid) : (0.,0.,0.)
        return stats.volume,center
    end
    face_kind=u->get(kernel.shapes[u[1]].geometry.geometry,:occ,nothing)
    torus=count(u->face_kind(u)===:torus,uses)
    # A sector solid is the torus face alone (full turn) or with its two
    # plane caps — a box with a toroidal cavity is torus+planes too, and
    # must not be taken for the cavity's own sector.
    if torus==1 && length(uses)<=3 &&
            all(u->face_kind(u)===:torus ||
                   face_kind(u)===:plane,uses)
        # A purely torus+plane boundary may be a torus sector — the
        # swept profile's Pappus integrals. V = φ·πr2²·r1; the centroid
        # rides the sector bisector at (r1+r2²/4r1)·sinc(φ/2) from
        # center. Verified against the actual wire below: trimmed
        # pieces share the carrier record over sub-domain wires and
        # fall through to the generic divergence path instead.
        tfid=only(fid for (fid,_) in uses
                  if get(kernel.shapes[fid].geometry.geometry,:occ,
                         nothing)===:torus)
        rec=kernel.shapes[tfid].geometry.geometry
        win=_occ_torus_sector_window(kernel,tfid,rec)
        if win!==nothing
            # Torus sector [umin,umax] — the swept profile's Pappus
            # integrals. V = φ·πr2²·r1; the centroid rides the sector
            # bisector at (r1+r2²/4r1)·sinc(φ/2) from center.
            φ=win[2]-win[1]
            V=φ*π*rec.r2*rec.r2*rec.r1
            mean_r=rec.r1+rec.r2*rec.r2/(4*rec.r1)
            arm=φ==0.0 ? 0.0 : mean_r*sin(φ/2)/(φ/2)
            bm=(win[1]+win[2])/2
            sb,cb=_gm_sincos(bm)
            return V,_occ3_add(rec.center,
                (arm*(cb*rec.X[1]+sb*rec.Y[1]),
                 arm*(cb*rec.X[2]+sb*rec.Y[2]),
                 arm*(cb*rec.X[3]+sb*rec.Y[3])))
        end
        # A trimmed torus piece shares the carrier's record over a
        # sub-domain wire (boolean kept faces, STEP trims, section
        # windows) — it is no sector, so it integrates through the
        # generic divergence path below via its pcurves.
    end
    all(u->face_kind(u)===:plane || face_kind(u)===:cylinder ||
           haskey(kernel.shapes[u[1]].geometry.geometry,:patch) ||
           haskey(kernel.shapes[u[1]].geometry.geometry,:pcurves),uses) ||
        throw(ErrorException("OCC.get_mass: unsupported solid boundary"))
    # Divergence theorem over outward-oriented faces:
    #   V = (1/3)Σ ∬_F r·n̂ dA,  ∭r dV = (1/4)Σ ∬_F r(r·n̂) dA.
    # Plane faces integrate exactly through `_model_occ_plane_stats`;
    # cylinder patches reduce to wire-boundary line integrals; NURBS
    # patches over the full knot domain integrate `r·(∂u×∂v)` per span.
    # Faces whose wire bounds only a sub-region — STEP-imported trimmed
    # patches and quadrics carrying pcurves — integrate the same
    # integrand over a CDT of the sampled uv boundary.
    V=0.0; M=(0.0,0.0,0.0)
    for (fid,rev) in uses
        rec=kernel.shapes[fid].geometry.geometry
        σ=rev ? -1.0 : 1.0
        if hasproperty(rec,:natural_bounds)
            integrals=_occ_natural_surface_integrals(rec)
            V+=σ*integrals.volume
            want_centroid && (M=_occ3_add(M,_occ3_scale(integrals.div_first,σ)))
            continue
        end
        if haskey(rec,:plane_rec)
            # Reclassified-from-plane tensor cap: the exact plane record
            # keeps shoelace precision on conic wires.
            prec=rec.plane_rec
            nh=_occ_plane_metric(prec)[2]
            n=(σ*nh[1],σ*nh[2],σ*nh[3])
            A,c3,_=_model_occ_plane_stats(kernel,fid;face_rec=prec)
            d=A*_dot3(n,c3)
            V+=d/3
            want_centroid && (M=_occ3_add(M,_occ3_scale(c3,d)))
            continue
        end
        if haskey(rec,:patch)
            if haskey(rec,:pcurves) && get(rec,:trimmed,false)
                Vp,Mp=_model_occ_trimmed_face_div(kernel,fid)
            else
                Vp,Mp=_model_occ_patch_face_div(kernel,fid;
                                              want_centroid=want_centroid)
            end
            V+=σ*Vp
            want_centroid && (M=_occ3_add(M,_occ3_scale(Mp,σ)))
            continue
        end
        occ=get(rec,:occ,nothing)
        if occ===:torus && haskey(rec,:pcurves)
            if want_centroid
                # The boundary reduction carries only the volume
                # integral — the centroid moment still walks the CDT
                # quadrature for trimmed torus pieces.
                Vp,Mp=_model_occ_trimmed_face_div(kernel,fid)
                V+=σ*Vp
                M=_occ3_add(M,_occ3_scale(Mp,σ))
            else
                V+=σ*_model_occ_torus_face_div(kernel,fid)
            end
            continue
        end
        if occ===:sphere && haskey(rec,:pcurves)
            J=_model_occ_sphere_face_integrals(kernel,fid)
            if J!==nothing
                # boundary integrals for the outward radial normal; a
                # left-handed frame's natural normal is the inward one
                h=_dot3(_occ_cross(rec.X,_occ_frame_y(rec)),rec.axis)>=0 ?
                    1.0 : -1.0
                Vp,Mp=_model_occ_sphere_face_div(rec,J)
                V+=σ*h*Vp
                want_centroid && (M=_occ3_add(M,_occ3_scale(Mp,σ*h)))
                continue
            end
        end
        if occ===:cone && haskey(rec,:pcurves) && !want_centroid
            exact=_model_occ_cone_face_integrals(kernel,fid)
            if exact!==nothing
                h=_dot3(_occ_cross(rec.X,_occ_frame_y(rec)),rec.axis)>=0 ?
                    1.0 : -1.0
                V+=σ*h*exact[2]
                continue
            end
        end
        if haskey(rec,:pcurves) &&
            occ!==:cylinder && occ!==:plane
            Vp,Mp=_model_occ_trimmed_face_div(kernel,fid)
            V+=σ*Vp
            want_centroid && (M=_occ3_add(M,_occ3_scale(Mp,σ)))
            continue
        end
        if get(rec,:occ,nothing)===:cylinder
            I1,Ic,Is,=_model_occ_cyl_face_div(kernel,fid)
            h=_dot3(_occ_cross(rec.X,_occ_frame_y(rec)),rec.axis)>=0 ?
                1.0 : -1.0
            Cx=_dot3(rec.center,rec.X); Cy=_dot3(rec.center,rec.Y)
            J0=Cx*Ic+Cy*Is+rec.radius*I1
            V+=σ*h*rec.radius*J0/3
            if want_centroid
                _,Mf,_=_model_occ_cyl_face_moment_terms(kernel,fid,
                    (0.0,0.0,0.0))
                M=_occ3_add(M,_occ3_scale(Mf,σ*h*rec.radius))
            end
            continue
        end
        nh=_occ_plane_metric(rec)[2]
        n=(σ*nh[1],σ*nh[2],σ*nh[3])
        # On a plane face r·n̂ is the constant face offset, so the centroid
        # moment ∬_F r(r·n̂)dA = (n̂·c̄)·A·c̄ — no per-edge tensor pass, and
        # conic-bounded plane faces stay supported.
        A,c3,_=_model_occ_plane_stats(kernel,fid)
        d=A*_dot3(n,c3)
        V+=d/3
        want_centroid && (M=_occ3_add(M,_occ3_scale(c3,d)))
    end
    # Upstream reports the signed divergence sum (inward-wound shells go
    # negative); a degenerate one leaves the empty-GProp centroid.
    abs(V)>_OCC_EPS_DIM || return V,(0.0,0.0,0.0)
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
# Gmsh's `GProps::MatrixOfInertia` is the centroidal inertia tensor,
# I = tr(Qc)I₃ − Qc with Qc = ∫(r-c)⊗(r-c)dμ (public 4.15.2 control: a unit
# box is diag(1/6), not the origin-moment diag(2/3) with −1/4 products).
# Curves and faces integrate the centered tensor per record — closed forms
# for line/circle, the plane fan or conic boundary 1-forms, and GL8 panels
# for the quadrature families (the canonical quadric domain, NURBS knot
# spans). Solids reuse the measure split: compact `primitive` records in
# closed form (sphere 2mr²/5·I, cylinder/cone axial forms through the
# ∭x²dV/∭z²dV layer integrals), and plane/cylinder/natural shells through
# the centered divergence theorem — ∭(r-c)_i(r-c)_j dV =
# n_k∬(r-c)_i(r-c)_j(r-c)_k dA/5. Unsupported boundaries reject like
# `get_center_of_mass`; non-1:3 dims and empty shells report the
# empty-GProp zero matrix.

function _model_occ_curve_moment2(kernel,id)
    rec=kernel.shapes[id].geometry.geometry
    rec.occ in (:line,:circle,:degenerate) ||
        return _model_occ_curve_stats(rec;want_second=true).second
    m=zeros(3,3)
    _,centroid=_model_occ_curve_measure(kernel,id)
    if rec.occ===:line
        o=_occ3_sub(rec.origin,centroid);d=rec.dir
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
        C=_occ3_sub(rec.center,centroid)
        X,Y,r=rec.X,rec.Y,rec.r
        for i in 1:3,j in 1:3
            m[i,j]=r*(C[i]*C[j]*Δ +
                r*(C[i]*(X[j]*Su+Y[j]*Cv)+(X[i]*Su+Y[i]*Cv)*C[j]) +
                r*r*(C2*X[i]*X[j]+CS*(X[i]*Y[j]+Y[i]*X[j])+S2*Y[i]*Y[j]))
        end
        return m
    elseif rec.occ===:degenerate
        return m
    end
    return m
end

# Cubic boundary moment U_F[i,j,k] = ∬_F (x-center)_i(x-center)_j
# (x-center)_k dA for a plane face — the centered divergence theorem needs
# ∭(r-c)_i(r-c)_j dV = n_k U_F[i,j,k]/5. Same signed fan as
# `_model_occ_plane_stats`, evaluated with the degree-3 Stroud rule (exact
# for cubics, hence also for the fan's per-triangle integrand); conic/spline
# wires use boundary 1-forms.
function _model_occ_plane_moment3(kernel,id;center=(0.0,0.0,0.0))
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
    A>_OCC_EPS_DIM || throw(ErrorException(
        "OCC.get_matrix_of_inertia: plane face has zero area"))
    cu=cv=0.0
    for (index,wire) in enumerate(raw_per_wire)
        area=abs(wire.raw[1])
        area>_OCC_EPS_DIM || continue
        w=(index==1 ? area : -area)/A
        # the wire's own centroid is moment over area with their common
        # winding sign — dividing by |area| negated it for a wire wound
        # clockwise (a mirrored face's circular hole)
        cu+=w*wire.raw[2]/wire.raw[1]; cv+=w*wire.raw[3]/wire.raw[1]
    end
    c2=(cu,cv); a3=_occ3_sub(_model_occ_plane_eval(rec,cu,cv),center)
    if _model_occ_plane_curved(kernel,shape,rec)
        wires=[kernel.shapes[w].geometry.uses for w in shape.children]
        _,raw=_model_occ_plane_boundary_moments(kernel,rec,wires)
        raw[(0,0)]>_OCC_EPS_DIM || throw(ErrorException(
            "OCC.get_matrix_of_inertia: plane face has zero area"))
        cu_b=raw[(1,0)]/raw[(0,0)];cv_b=raw[(0,1)]/raw[(0,0)]
        mc=_model_occ_centered_plane_moments(raw,cu_b,cv_b)
        base=_occ3_sub(_model_occ_plane_eval(rec,cu_b,cv_b),center)
        return _model_occ_plane_cubic(base,rec.X,rec.Y,mc).*
            _occ_plane_metric(rec)[1]
    end
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
            b3=_occ3_sub(_model_occ_plane_eval(rec,p1[1],p1[2]),center)
            e3=_occ3_sub(_model_occ_plane_eval(rec,p2[1],p2[2]),center)
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
    return U.*_occ_plane_metric(rec)[1]
end

# Centered face second-moment tensor T_F = ∬_F(r-c)⊗(r-c)dA over a trimmed
# parametric face — the CDT quadrature of `_model_occ_trimmed_face_measure`
# carrying the tensor integrand (the same 4×-split three-point rule). Used for
# any face whose wire bounds a sub-domain of its base — the canonical-domain
# quadric arms below would otherwise integrate the carrier's whole chart.
function _model_occ_trimmed_face_moment2(kernel,id;centroid=nothing)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    if centroid===nothing
        area,measured_centroid=_model_occ_trimmed_face_measure(kernel,id)
        centroid=measured_centroid
    end
    xs,ys,T,interior=_model_occ_trimmed_face_cdt(kernel,id,shape,rec)
    nodes=((1/6,1/6),(2/3,1/6),(1/6,2/3))
    m=zeros(3,3)
    for t in eachindex(interior)
        interior[t] || continue
        p1=(xs[T.tv[3t-2]],ys[T.tv[3t-2]])
        p2=(xs[T.tv[3t-1]],ys[T.tv[3t-1]])
        p3=(xs[T.tv[3t]],ys[T.tv[3t]])
        m12=_occ2_mid(p1,p2); m23=_occ2_mid(p2,p3); m31=_occ2_mid(p3,p1)
        for (q1,q2,q3) in ((p1,m12,m31),(m12,p2,m23),(m31,m23,p3),
                          (m12,m23,m31))
            jac=abs((q2[1]-q1[1])*(q3[2]-q1[2])-
                    (q3[1]-q1[1])*(q2[2]-q1[2]))
            for (s,tt) in nodes
                u=q1[1]+s*(q2[1]-q1[1])+tt*(q3[1]-q1[1])
                v=q1[2]+s*(q2[2]-q1[2])+tt*(q3[2]-q1[2])
                p,du,dv=_occ_surface_d1(rec,u,v)
                w=jac/2/3*_occ_modulus(_occ_cross(du,dv))
                d=_occ3_sub(p,centroid)
                for i in 1:3,j in 1:3
                    m[i,j]+=w*d[i]*d[j]
                end
            end
        end
    end
    return m
end

# Centered divergence face tensor ∬(p-c)⊗(p-c)((p-c)·n)dA over a trimmed
# parametric face. The same CDT and 4×-split rule as the mass/centroid
# divergence path; used for quadric and patch faces outside the exact
# plane/cylinder/torus-sector arms.
function _model_occ_trimmed_face_div_moment2(kernel,id,center)
    shape=kernel.shapes[id];rec=shape.geometry.geometry
    xs,ys,T,interior=_model_occ_trimmed_face_cdt(kernel,id,shape,rec)
    nodes=((1/6,1/6),(2/3,1/6),(1/6,2/3))
    Q=zeros(3,3)
    for t in eachindex(interior)
        interior[t] || continue
        p1=(xs[T.tv[3t-2]],ys[T.tv[3t-2]])
        p2=(xs[T.tv[3t-1]],ys[T.tv[3t-1]])
        p3=(xs[T.tv[3t]],ys[T.tv[3t]])
        m12=_occ2_mid(p1,p2);m23=_occ2_mid(p2,p3);m31=_occ2_mid(p3,p1)
        for (q1,q2,q3) in ((p1,m12,m31),(m12,p2,m23),(m31,m23,p3),
                          (m12,m23,m31))
            jac=abs((q2[1]-q1[1])*(q3[2]-q1[2])-
                    (q3[1]-q1[1])*(q2[2]-q1[2]))
            for (s,tt) in nodes
                u=q1[1]+s*(q2[1]-q1[1])+tt*(q3[1]-q1[1])
                v=q1[2]+s*(q2[2]-q1[2])+tt*(q3[2]-q1[2])
                p,du,dv=_occ_surface_d1(rec,u,v)
                w=jac/2/3
                normal=_occ_cross(du,dv)
                d=_occ3_sub(p,center)
                rn=_dot3(d,normal)
                for i in 1:3,j in 1:3
                    Q[i,j]+=w*rn*d[i]*d[j]/5
                end
            end
        end
    end
    return Q
end

function _model_occ_face_moment2(kernel,id)
    shape=kernel.shapes[id]; rec=shape.geometry.geometry
    hasproperty(rec,:natural_bounds) &&
        return _occ_natural_surface_integrals(rec;centered=true,
            fields=(:centered_second,)).centered_second
    occ=get(rec,:occ,nothing)
    if occ===:plane
        _,_,T=_model_occ_plane_stats(kernel,id;want_tensor=true,
            caller="OCC.get_matrix_of_inertia",centered=true)
        return T
    end
    # Public Gmsh reports the inertia about the entity's centroid; every
    # remaining arm integrates (r-c)⊗(r-c) directly.
    _,centroid=_model_occ_face_measure(kernel,id)
    m=zeros(3,3)
    # A wire that bounds only a sub-domain of the base chart — STEP trims
    # and Boolean kept pieces — cannot use the canonical-domain arms
    # below: integrate the CDT of the sampled uv boundary instead. The
    # sector-window check keeps genuine canonical torus sectors on the
    # exact arm.
    if (get(rec,:trimmed,false) && haskey(rec,:pcurves)) ||
         (occ===:torus && haskey(rec,:pcurves) &&
          _occ_torus_sector_window(kernel,id,rec)===nothing) ||
         ((occ===:cylinder || occ===:sphere || occ===:cone) &&
          !_occ_quadric_window_wire(kernel,shape,rec,
                                    _occ_quadric_window(rec)...))
        return _model_occ_trimmed_face_moment2(kernel,id;centroid)
    end
    if haskey(rec,:patch)
        return _model_occ_full_patch_integrals(rec;centered=true,
            center=centroid,fields=(:centered_second,)).centered_second
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
        # a Boolean can leave a genuine sector narrower than the carrier's
        # own sweep (a torus halved by a plane through its axis): the
        # closed form runs over the wire's window, not `rec.angle`
        win=haskey(rec,:pcurves) ?
            _occ_torus_sector_window(kernel,id,rec) : nothing
        win===nothing ? (0.0,rec.angle,0.0,_OCC_TWO_PI) :
            (win[1],win[2],0.0,_OCC_TWO_PI)
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
            d=_occ3_sub(p,centroid)
            for i2 in 1:3,j2 in 1:3
                m[i2,j2]+=w*dA*d[i2]*d[j2]
            end
        end
    end
    return m
end

# Torus sector (Pappus) centroidal second-moment tensor. The solid is the
# swept tube over u∈[u0,u1]; analytic radial/azimuth integrals avoid a
# curved-CDT approximation. Q is returned in the torus local frame and then
# mapped through the stored X/Y/axis basis.
function _occ_torus_sector_moment2(rec,win)
    u0,u1=win
    φ=u1-u0
    R=rec.r1;r=rec.r2
    s0,c0=_gm_sincos(u0);s1,c1=_gm_sincos(u1)
    s20,c20=_gm_sincos(2u0);s21,c21=_gm_sincos(2u1)
    Su=s1-s0; Cu=c0-c1
    U2=φ/2+(s21-s20)/4
    U2y=φ/2-(s21-s20)/4
    UV=(s1*s1-s0*s0)/2
    V=φ*π*R*r*r
    V>_OCC_EPS_DIM || return zeros(3,3)
    # ∫ρ(R+ρcosv)³ dρdv (the x²/y² weight) and the z² weight.
    qcoef=π*R*r*r*(R*R+3r*r/4)
    qzz=π*φ*R*r^4/4
    mcoef=π*r*r*(R*R+r*r/4)
    mX=mcoef*Su
    mY=mcoef*Cu
    cX=mX/V; cY=mY/V
    Q=zeros(3,3)
    Q[1,1]=qcoef*U2-V*cX*cX
    Q[2,2]=qcoef*U2y-V*cY*cY
    Q[1,2]=Q[2,1]=qcoef*UV-V*cX*cY
    Q[3,3]=qzz
    X=rec.X;Y=rec.Y;Z=rec.axis
    basis=(X,Y,Z)
    Qw=zeros(3,3)
    for a in 1:3,b in 1:3
        Q[a,b]==0.0 && continue
        for i in 1:3,j in 1:3
            Qw[i,j]+=Q[a,b]*basis[a][i]*basis[b][j]
        end
    end
    return Qw
end

# Centered divergence tensor over a full NURBS patch's knot domain.
function _model_occ_patch_face_div_moment2(kernel,id,center)
    rec=kernel.shapes[id].geometry.geometry
    surface=rec.patch.surface
    (u0,v0),(u1,v1)=_model_nurbs_surface_bounds(rec)
    uspans=_model_occ_span_cuts(surface.knots_u,u0,u1)
    vspans=_model_occ_span_cuts(surface.knots_v,v0,v1)
    Q=zeros(3,3)
    for (ua,ub) in uspans, (va,vb) in vspans
        hu=(ub-ua)/2;mu=(ua+ub)/2;hv=(vb-va)/2;mv=(va+vb)/2
        for i in 1:4,si in (-1,1),j in 1:4,sj in (-1,1)
            u=mu+si*hu*_OCC_GL_T[i];v=mv+sj*hv*_OCC_GL_T[j]
            w=_OCC_GL_W[i]*_OCC_GL_W[j]*hu*hv
            p,du,dv,_,_,_=_model_nurbs_surface_d2(rec,u,v)
            d=_occ3_sub(p,center)
            rn=_dot3(d,_occ_cross(du,dv))
            for i2 in 1:3,j2 in 1:3
                Q[i2,j2]+=w*rn*d[i2]*d[j2]/5
            end
        end
    end
    return Q
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
            qp=V*r.radius*r.radius/4; qa=V*r.height*r.height/12
            for i in 1:3,j in 1:3
                m[i,j]=qp*(Float64(i==j)-n[i]*n[j])+qa*n[i]*n[j]
            end
            return m
        elseif g.primitive.occ===:sphere
            V=4π/3*r.radius^3
            q=V*r.radius*r.radius/5
            for i in 1:3,j in 1:3
                m[i,j]=q*Float64(i==j)
            end
            return m
        elseif g.primitive.occ===:cone
            n=_occ_dir(r.axis)
            q2=r.r1*r.r1+r.r1*r.r2+r.r2*r.r2
            V=π*r.height*q2/3
            zb=r.height*(r.r1*r.r1+2*r.r1*r.r2+3*r.r2*r.r2)/(4q2)
            S4=r.r1^4+r.r1^3*r.r2+r.r1^2*r.r2^2+r.r1*r.r2^3+r.r2^4
            qp=π*r.height*S4/20
            qa=π*r.height^3*(r.r1^2+3*r.r1*r.r2+6*r.r2^2)/30-V*zb*zb
            for i in 1:3,j in 1:3
                m[i,j]=qp*(Float64(i==j)-n[i]*n[j])+qa*n[i]*n[j]
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
    if all(use->hasproperty(kernel.shapes[use[1]].geometry.geometry,:natural_bounds),uses)
        return _model_occ_natural_solid_stats(kernel,uses;want_second=true).second
    end
    face_kind=u->get(kernel.shapes[u[1]].geometry.geometry,:occ,nothing)
    torus_faces=[fid for (fid,_) in uses if face_kind((fid,false))===:torus]
    if length(torus_faces)==1 && length(uses)<=3 &&
            all(u->face_kind(u) in (:torus,:plane),uses)
        trec=kernel.shapes[only(torus_faces)].geometry.geometry
        win=_occ_torus_sector_window(kernel,only(torus_faces),trec)
        win===nothing || return _occ_torus_sector_moment2(trec,win)
    end
    face_ok=u->begin
        rec=kernel.shapes[u[1]].geometry.geometry
        face_kind(u) in (:plane,:cylinder) ||
            hasproperty(rec,:natural_bounds) ||
            haskey(rec,:pcurves) || haskey(rec,:patch)
    end
    all(face_ok,uses) || throw(ErrorException(
        "OCC.get_matrix_of_inertia: unsupported solid boundary"))
    V,centroid=_model_occ_solid_measure(kernel,id)
    abs(V)>_OCC_EPS_DIM || return m
    for (fid,rev) in uses
        rec=kernel.shapes[fid].geometry.geometry
        σ=rev ? -1.0 : 1.0
        if hasproperty(rec,:natural_bounds)
            integrals=_occ_natural_surface_integrals(rec;centered=true,
                center=centroid)
            for i in 1:3,j in 1:3
                m[i,j]+=σ*integrals.centered_div_second[i,j]
            end
        elseif rec.occ===:plane
            nh=_occ_plane_metric(rec)[2]
            n=(σ*nh[1],σ*nh[2],σ*nh[3])
            U=_model_occ_plane_moment3(kernel,fid;center=centroid)
            for i in 1:3,j in 1:3,k in 1:3
                m[i,j]+=n[k]*U[i,j,k]/5
            end
        elseif rec.occ===:cylinder
            h=_dot3(_occ_cross(rec.X,_occ_frame_y(rec)),rec.axis)>=0 ?
                1.0 : -1.0
            _,_,Q=_model_occ_cyl_face_moment_terms(kernel,fid,centroid)
            scale=σ*h*rec.radius/5
            for i in 1:3,j in 1:3
                m[i,j]+=scale*Q[i,j]
            end
        elseif haskey(rec,:pcurves)
            Q=_model_occ_trimmed_face_div_moment2(kernel,fid,centroid)
            for i in 1:3,j in 1:3
                m[i,j]+=σ*Q[i,j]
            end
        elseif haskey(rec,:patch)
            Q=_model_occ_patch_face_div_moment2(kernel,fid,centroid)
            for i in 1:3,j in 1:3
                m[i,j]+=σ*Q[i,j]
            end
        else
            throw(ErrorException(
                "OCC.get_matrix_of_inertia: unsupported solid boundary"))
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
