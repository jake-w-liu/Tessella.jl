# A BREP edge can be defined entirely by a curve on a surface. Keep that
# composition as its native geometry; a sampled 3-D fit cannot preserve its
# derivatives, parameters or subsequent affine images.
struct _ModelOCCSurfaceTraceCurve{F,P}
    family::Symbol
    face::F
    pc::P
    curve::NURBS.NURBSCurve # projection seeds and quadrature panels only
    transforms::Vector{@NamedTuple{linear::NTuple{9,Float64},loc::NTuple{3,Float64}}}
    bbox::NTuple{6,Float64}
end

function _model_occ_curve_d2(g::_ModelOCCSurfaceTraceCurve,parameter)
    uv,d,dd=_occ_pcurve_d2(g.pc,Float64(parameter))
    p,Su,Sv,Suu,Svv,Suv=_occ_surface_d2(g.face,uv...)
    tangent=ntuple(i->Su[i]*d[1]+Sv[i]*d[2],3)
    acceleration=ntuple(i->Suu[i]*d[1]^2+2Suv[i]*d[1]*d[2]+
        Svv[i]*d[2]^2+Su[i]*dd[1]+Sv[i]*dd[2],3)
    # Apply finite affine images in order. Multiplying their coefficients
    # can overflow even when every transformed point and derivative is finite.
    point=p
    for transform in g.transforms
        point=_occ_matvec(transform.linear,point).+transform.loc
        tangent=_occ_matvec(transform.linear,tangent)
        acceleration=_occ_matvec(transform.linear,acceleration)
    end
    result=(point,tangent,acceleration)
    all(vector->all(isfinite,vector),result) || throw(ArgumentError(
        "OCC Curve: surface trace derivatives exceed Float64"))
    return result
end

_model_occ_curve_bounds(g::_ModelOCCSurfaceTraceCurve)=
    (g.curve.knots[g.curve.degree+1],g.curve.knots[end-g.curve.degree])
_model_occ_curve_bbox(g::_ModelOCCSurfaceTraceCurve)=g.bbox

function _model_occ_curve_project(g::_ModelOCCSurfaceTraceCurve,target,lower,upper)
    lo,hi=minmax(lower,upper)
    cuts=sort!(unique!(Float64[lo;[k for k in g.curve.knots if lo<k<hi];hi]))
    states=Tuple{Float64,NTuple{3,Float64},Float64}[]
    tolerance=(hi-lo)*1e-9+eps(Float64)
    for index in 1:length(cuts)-1
        _model_occ_curve_interval_extrema!(states,g,target,
            cuts[index],cuts[index+1],4,tolerance)
    end
    for parameter in (lo,hi)
        point,_,_,_,distance=_model_occ_curve_extrema_state(g,target,parameter)
        push!(states,(parameter,point,distance))
    end
    best=argmin(state->state[3],states)
    return best[1],best[2]
end

function _model_occ_transform_nurbs_curve(g::_ModelOCCSurfaceTraceCurve,t,caller)
    poles=NTuple{3,Float64}[_finite_result(_occ_trsf_apply(t,p),caller)
        for p in g.curve.controls]
    seed=NURBS.NURBSCurve(g.curve.degree,g.curve.knots,poles,g.curve.weights)
    box=g.bbox
    bounds=_model_occ_box_points(_occ_trsf_apply(t,(x,y,z))
        for x in (box[1],box[4]),y in (box[2],box[5]),z in (box[3],box[6]))
    all(isfinite,bounds) || throw(ArgumentError(
        "$caller: surface trace transform exceeds Float64"))
    transforms=copy(g.transforms);push!(transforms,t.occ)
    return _ModelOCCSurfaceTraceCurve(g.family,g.face,g.pc,seed,transforms,bounds)
end

# Conservative source-coordinate bounds. Polynomial/rational control hulls
# and full conic bounds cover the complete trace, including inter-sample
# extrema. These bounds describe the support, never the seed polyline.
function _occ_brep_curve2_bounds(curve,lower,upper,caller)
    base=_occ_brep_untrim(curve)
    if base.kind==1
        points=(base.origin.+lower.*base.dir,base.origin.+upper.*base.dir)
        lo=ntuple(i->min(points[1][i],points[2][i]),2)
        hi=ntuple(i->max(points[1][i],points[2][i]),2)
    elseif base.kind in (2,3)
        ry=base.kind==2 ? base.a : base.b
        radii=ntuple(i->abs(base.a*base.X[i])+abs(ry*base.Y[i]),2)
        lo=base.center.-radii;hi=base.center.+radii
    elseif base.kind in (6,7)
        native=_occ_brep_nurbs_curve(base,caller)
        lo=ntuple(i->minimum(p[i] for p in native.controls),2)
        hi=ntuple(i->maximum(p[i] for p in native.controls),2)
    elseif base.kind in (4,5)
        bound=max(abs(lower),abs(upper))
        base.a>0 || throw(ArgumentError("$caller: invalid BREP conic radius"))
        x=base.kind==4 ? bound^2/(4base.a) : abs(base.a)*cosh(bound)
        y=base.kind==4 ? bound : abs(base.b)*sinh(bound)
        radii=ntuple(i->x*abs(base.X[i])+y*abs(base.Y[i]),2)
        lo=base.center.-radii;hi=base.center.+radii
    else
        # The evaluator reports the precise unsupported family before this
        # bound is used; a loose or sampled bound would misstate geometry.
        _occ_brep_pcurve(base,caller)
        throw(ArgumentError("$caller: BREP pcurve has no finite native bound"))
    end
    all(isfinite,(lo...,hi...)) || throw(ArgumentError(
        "$caller: BREP pcurve bounds exceed Float64"))
    return lo,hi
end

function _occ_brep_trace_curve(surface,curve,lower,upper,caller)
    lower<upper || throw(ArgumentError("$caller: BREP edge range must increase"))
    pc=_occ_brep_pcurve(curve,caller)
    lo,hi=_occ_brep_curve2_bounds(curve,lower,upper,caller)
    # A constant coordinate still needs a regular two-dimensional carrier
    # for D2. Extending an extrusion/conic's analytic support does not alter
    # any point of the edge, and no extra Face entity is published.
    supports=ntuple(2) do axis
        lo[axis]!=hi[axis] && return (lo[axis],hi[axis])
        lower_support=lo[axis]-.5;upper_support=hi[axis]+.5
        if lower_support==upper_support
            lower_support=isfinite(prevfloat(lo[axis])) ? prevfloat(lo[axis]) : lo[axis]
            upper_support=isfinite(nextfloat(hi[axis])) ? nextfloat(hi[axis]) : hi[axis]
        end
        lower_support<upper_support || throw(ArgumentError(
            "$caller: BREP trace has no representable surface carrier"))
        return lower_support,upper_support
    end
    support_lo=first.(supports);support_hi=last.(supports)
    face,scale,shift=_occ_brep_surface_geometry(surface,(support_lo,support_hi),caller)
    actual_pc=_occ_brep_pcurve_map(pc,1.,0.,scale,shift)
    mapped_lo=ntuple(i->scale[i]*lo[i]+shift[i],2)
    mapped_hi=ntuple(i->scale[i]*hi[i]+shift[i],2)
    box_domain=face.occ in (:cylinder,:cone) ?
        ((0.,mapped_lo[2]),(_OCC_TWO_PI,mapped_hi[2])) : (mapped_lo,mapped_hi)
    box=_occ_natural_surface_bbox(merge(face,(natural_bounds=box_domain,)))
    # Keep enough parameter panels to distinguish repeated periodic turns.
    periods=(_brep_uperiod(face),_brep_vperiod(face))
    count=32
    for axis in 1:2
        period=periods[axis]
        period>0 || continue
        turns=(hi[axis]-lo[axis])/period
        isfinite(turns) && turns<=8192 || throw(ArgumentError(
            "$caller: BREP surface trace exceeds periodic panel capacity"))
        count=max(count,ceil(Int,8turns))
    end
    count<=65536 || throw(ArgumentError("$caller: BREP surface trace exceeds panel capacity"))
    cuts=sort!(unique!(Float64[
        _convex_coordinate(lower,upper,index/count) for index in 0:count]))
    points=NTuple{3,Float64}[_occ_surface_point(face,_occ_pcurve_d0(actual_pc,t)...)
        for t in cuts]
    all(point->all(isfinite,point),points) || throw(ArgumentError(
        "$caller: BREP surface trace exceeds Float64"))
    seed=NURBS.NURBSCurve(1,[first(cuts);cuts;last(cuts)],points)
    native=_ModelOCCSurfaceTraceCurve(:parametric,face,actual_pc,seed,
        @NamedTuple{linear::NTuple{9,Float64},loc::NTuple{3,Float64}}[],box)
    return (occ=:copied_curve,curve=native,t0=lower,t1=upper)
end
