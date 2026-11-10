# Native edge/surface coordinate maps for implicit planar pcurves and
# BREP edges whose stored pcurve has a different parameterization.
function _occ_brep_edge_d2(g,t::Float64)
    g.occ===:line && return (_occ_line_point(g,t),g.dir,(0.0,0.0,0.0))
    g.occ===:circle && return (_occ_circle_point(g,t),
        _occ_circle_derivative(g,t),_occ_circle_second_derivative(g,t))
    g.occ===:ellipse && return (_occ_ellipse_point(g,t),
        _occ_ellipse_derivative(g,t),_occ_ellipse_second_derivative(g,t))
    g.occ===:copied_curve && return _model_occ_curve_d2(g.curve,t)
    throw(ArgumentError("BREP pcurve: edge has no regular native evaluator"))
end

function _occ_brep_plane_uv(g,point;vector=false)
    q=vector ? point : point.-g.center
    sx=maximum(abs,g.X); sy=maximum(abs,g.Y)
    sx>0 && sy>0 || throw(ArgumentError("BREP plane has a singular frame"))
    X=g.X./sx; Y=g.Y./sy
    xx=_dot3(X,X); xy=_dot3(X,Y); yy=_dot3(Y,Y)
    determinant=xx*yy-xy*xy
    determinant>0 || throw(ArgumentError("BREP plane has a singular frame"))
    x=_dot3(q,X); y=_dot3(q,Y)
    uv=((yy*x-xy*y)/determinant/sx,(xx*y-xy*x)/determinant/sy)
    all(isfinite,uv) || throw(ArgumentError("BREP plane parameters exceed Float64"))
    return uv
end

function _occ_brep_plane_pcurve(edge,plane)
    if edge.occ===:line
        return (lin2d=(_occ_brep_plane_uv(plane,edge.origin),
                       _occ_brep_plane_uv(plane,edge.dir;vector=true)),)
    elseif edge.occ in (:circle,:ellipse)
        center=_occ_brep_plane_uv(plane,edge.center)
        X=_occ_brep_plane_uv(plane,edge.X;vector=true)
        Y=_occ_brep_plane_uv(plane,edge.Y;vector=true)
        sx=sqrt(sum(abs2,X)); sy=sqrt(sum(abs2,Y))
        sx>0 && sy>0 || throw(ArgumentError("BREP conic does not span its plane"))
        x=X./sx
        a=(edge.occ===:circle ? edge.r : edge.rx)*sx
        b=(edge.occ===:circle ? edge.r : edge.ry)*sy
        inner=abs(a-b)<=32eps(Float64)*max(a,b) ?
            (circ2d=(center,x,a),) : (ellips2d=(center,x,a,b),)
        return X[1]*Y[2]-X[2]*Y[1]<0 ? (neg2d=(inner,)) : inner
    elseif edge.occ===:copied_curve && edge.curve isa _ModelOCCNURBSCurve
        cv=edge.curve
        poles=[_occ_brep_plane_uv(plane,p) for p in cv.curve.controls]
        if all(==(cv.curve.weights[1]),cv.curve.weights)
            return (spl2d=(cv.curve.knots,cv.curve.degree,poles),)
        end
        nc=NURBS.NURBSCurve(cv.curve.degree,cv.curve.knots,
            [(p[1],p[2],0.0) for p in poles],cv.curve.weights)
        return (brep_nurbs2d=_ModelOCCNURBSCurve(cv.family,nc),)
    end
    return (brep_plane2d=(edge,plane),)
end

function _occ_pcurve_d2(pc::NamedTuple{(:brep_plane2d,)},t::Float64)
    edge,plane=pc.brep_plane2d
    p,d,dd=_occ_brep_edge_d2(edge,t)
    return _occ_brep_plane_uv(plane,p),
        _occ_brep_plane_uv(plane,d;vector=true),
        _occ_brep_plane_uv(plane,dd;vector=true)
end
_occ_pcurve_d0(pc::NamedTuple{(:brep_plane2d,)},t::Float64)=first(_occ_pcurve_d2(pc,t))
function _occ_pcurve_d1(pc::NamedTuple{(:brep_plane2d,)},t::Float64)
    p,d,_=_occ_pcurve_d2(pc,t)
    return p,d
end

function _occ_brep_inverse_trace(data,t)
    edge,face,source,lo,hi,tolerance=data
    point=_occ_brep_edge_d2(edge,t)[1]
    fraction=(t-edge.t0)/(edge.t1-edge.t0)
    seed=clamp(lo+fraction*(hi-lo),lo,hi)
    sample=q->begin
        sample_uv=_occ_pcurve_d0(source,Float64(q))
        sample_point=_occ_surface_point(face,Float64(sample_uv[1]),Float64(sample_uv[2]))
        return sample_uv,sample_point,sum(abs2,sample_point.-point)
    end
    scale=max(1.,maximum(abs,point))
    epsilon=max(64eps(Float64)*scale,1e-13*scale)
    best=seed
    uv,mapped,distance=sample(best)
    candidates=Tuple{Float64,Float64}[(distance,seed)]
    # A closed edge has two coincident endpoint images. Its stored domain
    # and direction select which UV endpoint belongs to each native end.
    if (t==edge.t0 || t==edge.t1) && sqrt(distance)<=max(tolerance,epsilon)
        return uv
    end
    for index in 0:64
        q=lo+(hi-lo)*index/64
        candidate,position,squared=sample(q)
        push!(candidates,(squared,q))
        if squared<distance-epsilon^2 ||
                abs(squared-distance)<=epsilon^2 && abs(q-seed)<abs(best-seed)
            best=q; uv=candidate; mapped=position; distance=squared
        end
    end
    function refine(start)
        local q=start
        local candidate,position,squared=sample(q)
        for _ in 1:96
            squared<=epsilon^2 && break
            _,direction=_occ_pcurve_d1(source,Float64(q))
            _,du,dv,_,_,_=_occ_surface_d2(face,Float64(candidate[1]),Float64(candidate[2]))
            tangent=ntuple(i->du[i]*direction[1]+dv[i]*direction[2],3)
            denominator=_dot3(tangent,tangent)
            denominator>0 || break
            step=_dot3(position.-point,tangent)/denominator
            isfinite(step) || break
            accepted=false
            for power in 0:24
                trial=clamp(q-step*2.0^(-power),lo,hi)
                trial==q && continue
                next_uv,next_point,next_squared=sample(trial)
                if next_squared<squared
                    q=trial; candidate=next_uv; position=next_point; squared=next_squared
                    accepted=true; break
                end
            end
            accepted || break
        end
        return q,candidate,position,squared
    end
    sort!(candidates;by=entry->(entry[1],abs(entry[2]-seed)))
    starts=unique(vcat(best,lo,hi,last.(candidates[1:min(8,length(candidates))])))
    for start in starts
        q,candidate,position,squared=refine(start)
        if squared<distance-epsilon^2 ||
                abs(squared-distance)<=epsilon^2 && abs(q-seed)<abs(best-seed)
            best=q; uv=candidate; mapped=position; distance=squared
        end
    end
    # Invert along the actual source pcurve. A nearest-angle hint can jump
    # by one period under a nonlinear speed map even when XYZ is exact.
    residual=sqrt(distance)
    residual<=max(tolerance,64eps(Float64)*maximum(abs,point)) || throw(ArgumentError(
        "BREP pcurve: 3-D edge and surface disagree by $residual"))
    return uv
end

_occ_pcurve_d0(pc::NamedTuple{(:brep_inverse2d,)},t::Float64)=
    _occ_brep_inverse_trace(pc.brep_inverse2d,t)

function _occ_pcurve_d2(pc::NamedTuple{(:brep_inverse2d,)},t::Float64)
    data=pc.brep_inverse2d
    edge,face,source,lo,hi,_=data
    uv=_occ_brep_inverse_trace(data,t)
    _,tangent,acceleration=_occ_brep_edge_d2(edge,t)
    _,du,dv,duu,dvv,duv=_occ_surface_d2(face,uv[1],uv[2])
    a=_dot3(du,du); b=_dot3(du,dv); c=_dot3(dv,dv)
    determinant=a*c-b*b
    regular=determinant>64eps(Float64)*a*c
    solve(vector)=begin
        x=_dot3(du,vector); y=_dot3(dv,vector)
        regular ? ((c*x-b*y)/determinant,(a*y-b*x)/determinant) :
            a>=c && a>0 ? (x/a,0.0) : c>0 ? (0.0,y/c) : (0.0,0.0)
    end
    tangent_uv=solve(tangent)
    correction=ntuple(i->acceleration[i]-duu[i]*tangent_uv[1]^2-
        dvv[i]*tangent_uv[2]^2-2duv[i]*tangent_uv[1]*tangent_uv[2],3)
    first=tangent_uv
    second=solve(correction)
    if !regular
        # Resolve the invisible pole coordinate by the one-sided regular
        # limit of this exact same surface inverse, retaining its UV branch.
        span=edge.t1-edge.t0
        h=min(span/16,max(span*1e-6,64eps(Float64)*max(1.,abs(t))))
        direction=t+h<=edge.t1 ? 1.0 : -1.0
        t1=clamp(t+direction*h,edge.t0,edge.t1)
        t2=clamp(t+direction*2h,edge.t0,edge.t1)
        if t1!=t && t2!=t1
            p1=_occ_brep_inverse_trace(data,t1)
            p2=_occ_brep_inverse_trace(data,t2)
            first=ntuple(i->direction*(-3uv[i]+4p1[i]-p2[i])/(2h),2)
            second=ntuple(i->(uv[i]-2p1[i]+p2[i])/(h*h),2)
        end
    end
    return uv,first,second
end

function _occ_pcurve_d1(pc::NamedTuple{(:brep_inverse2d,)},t::Float64)
    p,d,_=_occ_pcurve_d2(pc,t)
    return p,d
end
