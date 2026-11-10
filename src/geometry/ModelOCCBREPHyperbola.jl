# A hyperbola's native parameter is hyperbolic, not a rational B-spline
# parameter. Its exact quadratic rational image supplies projection and
# swept-surface geometry; tanh transports the original parameter.
# https://dev.opencascade.org/doc/refman/html/class_geom___hyperbola.html
struct _ModelOCCHyperbolicReparam
    knots::Vector{Float64}
    tanhq::Vector{Float64}
end

function _occ_conic_reparam(r::_ModelOCCHyperbolicReparam,t)
    index=clamp(searchsortedlast(r.knots,t),1,length(r.tanhq))
    a,b=r.knots[index],r.knots[index+1]
    span=b-a;middle=a/2+b/2
    tangent=tanh((t-middle)/2)
    parameter=clamp(a+span/2*(1+tangent/r.tanhq[index]),a,b)
    derivative=span/(4r.tanhq[index])*(1-tangent^2)
    return parameter,derivative,-derivative*tangent
end

function _occ_conic_reparam_inv(r::_ModelOCCHyperbolicReparam,s)
    index=clamp(searchsortedlast(r.knots,s),1,length(r.tanhq))
    a,b=r.knots[index],r.knots[index+1]
    tangent=(2*(s-a)/(b-a)-1)*r.tanhq[index]
    abs(tangent)<1 || throw(ArgumentError("OCC hyperbola: parameter exceeds finite branch"))
    return a/2+b/2+2atanh(tangent)
end

function _occ_hyperbolic_inverse_d2(r::_ModelOCCHyperbolicReparam,s)
    index=clamp(searchsortedlast(r.knots,s),1,length(r.tanhq))
    a,b=r.knots[index],r.knots[index+1]
    q=r.tanhq[index];z=(2*(s-a)/(b-a)-1)*q
    abs(z)<1 || throw(ArgumentError("OCC hyperbola: rational parameter exceeds finite branch"))
    first=4q/((b-a)*(1-z*z))
    return a/2+b/2+2atanh(z),first,first*first*z
end

function _occ_hyperbolic_full_map(mapping::_ModelOCCHyperbolicReparam)
    lo,hi=first(mapping.knots),last(mapping.knots)
    return _ModelOCCHyperbolicReparam(Float64[lo,hi],Float64[tanh((hi-lo)/4)])
end

function _occ_hyperbolic_surface_maps(geometry)
    mappings=get(geometry,:reparam,(u=nothing,v=nothing))
    return (u=mappings.u isa _ModelOCCHyperbolicReparam ? _occ_hyperbolic_full_map(mappings.u) : nothing,
        v=mappings.v isa _ModelOCCHyperbolicReparam ? _occ_hyperbolic_full_map(mappings.v) : nothing)
end

function _occ_hyperbolic_surface_axis(surface,axis,mapping,caller)
    lo,hi=first(mapping.knots),last(mapping.knots)
    parameter=first(_occ_conic_reparam(mapping,lo/2+hi/2))
    weight=cosh((hi-lo)/2)
    count=size(surface.controls,3-axis)
    controls=Matrix{NTuple{3,Float64}}(undef,axis==1 ? 3 : count,axis==1 ? count : 3)
    weights=Matrix{Float64}(undef,size(controls))
    knots=axis==1 ? surface.knots_u : surface.knots_v
    degree=axis==1 ? surface.degree_u : surface.degree_v
    for index in 1:count
        poles=axis==1 ? collect(@view surface.controls[:,index]) : collect(@view surface.controls[index,:])
        old_weights=axis==1 ? collect(@view surface.weights[:,index]) : collect(@view surface.weights[index,:])
        first(old_weights)==last(old_weights) || throw(ArgumentError(
            "$caller: hyperbolic surface has inconsistent endpoint weights"))
        curve=_ModelOCCNURBSCurve(:bspline,NURBS.NURBSCurve(degree,knots,poles,old_weights))
        middle=first(_model_occ_curve_d2(curve,parameter))
        p0,p2=first(poles),last(poles)
        pole=ntuple(i->middle[i]+(middle[i]-(p0[i]/2+p2[i]/2))/weight,3)
        for local_index in 1:3
            row,column=axis==1 ? (local_index,index) : (index,local_index)
            controls[row,column]=local_index==1 ? p0 : local_index==3 ? p2 : pole
            weights[row,column]=first(old_weights)*(local_index==2 ? weight : 1.)
        end
    end
    full_knots=Float64[lo,lo,lo,hi,hi,hi]
    return axis==1 ? NURBS.NURBSSurface(2,surface.degree_v,full_knots,
        surface.knots_v,controls,weights) : NURBS.NURBSSurface(surface.degree_u,2,
        surface.knots_u,full_knots,controls,weights)
end

function _model_occ_hyperbolic_surface_nurbs(geometry,caller)
    mappings=geometry.reparam
    surface=geometry.patch.surface
    mappings.u isa _ModelOCCHyperbolicReparam &&
        (surface=_occ_hyperbolic_surface_axis(surface,1,mappings.u,caller))
    mappings.v isa _ModelOCCHyperbolicReparam &&
        (surface=_occ_hyperbolic_surface_axis(surface,2,mappings.v,caller))
    remaining=(u=mappings.u isa _ModelOCCHyperbolicReparam ? nothing : mappings.u,
        v=mappings.v isa _ModelOCCHyperbolicReparam ? nothing : mappings.v)
    out=merge(geometry,(patch=_ModelNURBSSurface(:bspline,surface),))
    return remaining==(u=nothing,v=nothing) ? Base.structdiff(out,(reparam=nothing,)) :
        merge(out,(reparam=remaining,))
end

function _occ_pcurve_d2(pc::NamedTuple{(:brep_conversion2d,)},t::Float64)
    base,input,output=pc.brep_conversion2d
    parameter,dt,ddt=input===nothing ? (t,1.,0.) : _occ_hyperbolic_inverse_d2(input,t)
    point,first,second=_occ_pcurve_d2(base,parameter)
    tangent=ntuple(i->dt*first[i],2)
    acceleration=ntuple(i->ddt*first[i]+dt^2*second[i],2)
    states=ntuple(2) do axis
        mapping=axis==1 ? output.u : output.v
        mapping===nothing && return (point[axis],tangent[axis],acceleration[axis])
        value,d,dd=_occ_conic_reparam(mapping,point[axis])
        return (value,d*tangent[axis],dd*tangent[axis]^2+d*acceleration[axis])
    end
    return ntuple(i->states[i][1],2),ntuple(i->states[i][2],2),ntuple(i->states[i][3],2)
end

_occ_pcurve_d0(pc::NamedTuple{(:brep_conversion2d,)},t::Float64)=first(_occ_pcurve_d2(pc,t))

struct _ModelOCCHyperbolaCurve{G}
    family::Symbol
    source::G
    curve::NURBS.NURBSCurve
    nets::NTuple{3,Union{Nothing,Vector{NTuple{4,Float64}}}}
    reparam::_ModelOCCHyperbolicReparam
    transforms::Vector{@NamedTuple{linear::NTuple{9,Float64},loc::NTuple{3,Float64}}}
end

_model_occ_native_curve_type(g::_ModelOCCHyperbolaCurve)="Hyperbola"

function _occ_brep_hyperbola_d2(curve,t)
    sh,ch=sinh(t),cosh(t)
    point=ntuple(i->curve.center[i]+curve.a*ch*curve.X[i]+curve.b*sh*curve.Y[i],3)
    first=ntuple(i->curve.a*sh*curve.X[i]+curve.b*ch*curve.Y[i],3)
    second=ntuple(i->curve.a*ch*curve.X[i]+curve.b*sh*curve.Y[i],3)
    all(v->all(isfinite,v),(point,first,second)) || throw(ArgumentError(
        "OCC hyperbola: derivatives exceed Float64"))
    return point,first,second
end

function _occ_brep_hyperbola(curve,lo,hi,caller;max_span=.5)
    curve.a>=0 && curve.b>=0 || throw(ArgumentError(
        "$caller: hyperbola radii must be nonnegative"))
    width=hi-lo
    isfinite(width) && width>0 || throw(ArgumentError(
        "$caller: hyperbola range must be finite and increasing"))
    # Short hyperbolic spans keep the rational tangent-intersection weights
    # near one and supply smooth local measure/projection intervals.
    isfinite(max_span) && max_span>0 || throw(ArgumentError(
        "$caller: hyperbola span width must be finite and positive"))
    count=ceil(Int,min(width/max_span,16385.))
    count<=16384 || throw(ArgumentError("$caller: hyperbola span capacity exceeded"))
    cuts=Float64[_convex_coordinate(lo,hi,i/count) for i in 0:count]
    all(i->cuts[i]<cuts[i+1],1:count) || throw(ArgumentError(
        "$caller: hyperbola spans are not representable"))
    poles=NTuple{3,Float64}[_occ_brep_hyperbola_d2(curve,lo)[1]]
    weights=Float64[1.];tanhq=Float64[]
    for index in 1:count
        a,b=cuts[index],cuts[index+1]
        middle=a/2+b/2;weight=cosh((b-a)/2)
        sh,ch=sinh(middle),cosh(middle)
        pole=ntuple(i->curve.center[i]+(curve.a*ch/weight)*curve.X[i]+
            (curve.b*sh/weight)*curve.Y[i],3)
        all(isfinite,pole) || throw(ArgumentError(
            "$caller: hyperbola rational controls exceed Float64"))
        push!(poles,pole,_occ_brep_hyperbola_d2(curve,b)[1])
        push!(weights,weight,1.);push!(tanhq,tanh((b-a)/4))
    end
    knots=Float64[lo,lo,lo]
    for knot in cuts[2:end-1];push!(knots,knot,knot);end
    append!(knots,(hi,hi,hi))
    native=NURBS.NURBSCurve(2,knots,poles,weights)
    inner=_ModelOCCNURBSCurve(:bspline,native)
    return _ModelOCCHyperbolaCurve(:hyperbola,curve,native,inner.nets,
        _ModelOCCHyperbolicReparam(cuts,tanhq),
        @NamedTuple{linear::NTuple{9,Float64},loc::NTuple{3,Float64}}[])
end

# Explicit conversion exposes the rational parameter. OCCT's hyperbola
# conversion uses one quadratic span over the trimmed edge range; retaining
# the source hyperbolic evaluator would leave its interior parameters wrong.
function _model_occ_hyperbola_nurbs(g::_ModelOCCHyperbolaCurve,caller)
    lo,hi=_model_occ_curve_bounds(g)
    source=_occ_brep_hyperbola(g.source,lo,hi,caller;max_span=hi-lo)
    poles=copy(source.curve.controls)
    for index in eachindex(poles)
        point=poles[index]
        for transform in g.transforms
            point=_finite_result(_occ_matvec(transform.linear,point).+transform.loc,caller)
        end
        poles[index]=point
    end
    native=NURBS.NURBSCurve(2,source.curve.knots,poles,source.curve.weights)
    return _ModelOCCNURBSCurve(:bspline,native)
end

function _model_occ_curve_d2(g::_ModelOCCHyperbolaCurve,t)
    point,first,second=_occ_brep_hyperbola_d2(g.source,Float64(t))
    for transform in g.transforms
        point=_occ_matvec(transform.linear,point).+transform.loc
        first=_occ_matvec(transform.linear,first)
        second=_occ_matvec(transform.linear,second)
    end
    all(v->all(isfinite,v),(point,first,second)) || throw(ArgumentError(
        "OCC hyperbola: transformed derivatives exceed Float64"))
    return point,first,second
end

_model_occ_curve_bounds(g::_ModelOCCHyperbolaCurve)=(first(g.reparam.knots),last(g.reparam.knots))

function _model_occ_curve_bbox(g::_ModelOCCHyperbolaCurve)
    axis_a=g.source.a.*g.source.X;axis_b=g.source.b.*g.source.Y
    for transform in g.transforms
        axis_a=_occ_matvec(transform.linear,axis_a)
        axis_b=_occ_matvec(transform.linear,axis_b)
    end
    all(isfinite,(axis_a...,axis_b...)) || return _model_occ_controls_box(g.curve.controls)
    lo,hi=_model_occ_curve_bounds(g)
    parameters=Float64[lo,hi]
    for axis in 1:3
        abs(axis_b[axis])<abs(axis_a[axis]) || continue
        parameter=atanh(-axis_b[axis]/axis_a[axis])
        lo<parameter<hi && push!(parameters,parameter)
    end
    bounds=_model_occ_box_points(first(_model_occ_curve_d2(g,t)) for t in parameters)
    return (bounds[1]-1e-7,bounds[2]-1e-7,bounds[3]-1e-7,
        bounds[4]+1e-7,bounds[5]+1e-7,bounds[6]+1e-7)
end

function _model_occ_curve_resolution(g::_ModelOCCHyperbolaCurve)
    inner=_ModelOCCNURBSCurve(:bspline,g.curve,g.nets)
    derivative_scale=maximum((g.reparam.knots[i+1]-g.reparam.knots[i])/
        (4g.reparam.tanhq[i]) for i in eachindex(g.reparam.tanhq))
    return _model_occ_curve_resolution(inner)/derivative_scale
end

function _model_occ_curve_project(g::_ModelOCCHyperbolaCurve,target,lower,upper)
    inner=_ModelOCCNURBSCurve(:bspline,g.curve,g.nets)
    result=_model_occ_curve_project(inner,target,
        first(_occ_conic_reparam(g.reparam,lower)),
        first(_occ_conic_reparam(g.reparam,upper)))
    result===nothing && return nothing
    return _occ_conic_reparam_inv(g.reparam,result[1]),result[2]
end

function _model_occ_transform_nurbs_curve(g::_ModelOCCHyperbolaCurve,t,caller)
    poles=NTuple{3,Float64}[_finite_result(_occ_trsf_apply(t,p),caller)
        for p in g.curve.controls]
    native=NURBS.NURBSCurve(2,g.curve.knots,poles,g.curve.weights)
    inner=_ModelOCCNURBSCurve(:bspline,native)
    transforms=copy(g.transforms);push!(transforms,t.occ)
    return _ModelOCCHyperbolaCurve(g.family,g.source,native,inner.nets,
        g.reparam,transforms)
end
