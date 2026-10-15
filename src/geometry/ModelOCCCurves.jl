# Standalone OCC curve geometry owns copied poles. The pending shape registry
# retains endpoint shape identities independently of this evaluator's controls.
struct _ModelOCCNURBSCurve
    family::Symbol
    curve::NURBS.NURBSCurve
    nets::NTuple{3,Union{Nothing,Vector{NTuple{4,Float64}}}}
end

function _model_occ_curve_derivative_net(net,degree,knots)
    degree==0 && return nothing
    result=Vector{NTuple{4,Float64}}(undef,length(net)-1)
    for index in eachindex(result)
        denominator=knots[index+degree+1]-knots[index+1]
        factor=denominator==0.0 ? 0.0 : degree/denominator
        first,last=net[index],net[index+1]
        result[index]=ntuple(axis->factor*(last[axis]-first[axis]),4)
    end
    return result
end

function _ModelOCCNURBSCurve(family,curve::NURBS.NURBSCurve)
    net=NTuple{4,Float64}[(weight*point[1],weight*point[2],weight*point[3],weight)
                        for (point,weight) in zip(curve.controls,curve.weights)]
    all(point->all(isfinite,point),net) || throw(ArgumentError(
        "OCC Curve: homogeneous controls must be finite"))
    first=_model_occ_curve_derivative_net(net,curve.degree,curve.knots)
    second=curve.degree<2 ? nothing :
        _model_occ_curve_derivative_net(first,curve.degree-1,curve.knots[2:end-1])
    return _ModelOCCNURBSCurve(family,curve,(net,first,second))
end

function _model_occ_curve_homogeneous(geometry,order,parameter)
    net=geometry.nets[order+1]
    net===nothing && return (0.0,0.0,0.0,0.0)
    knots=geometry.curve.knots[1+order:end-order]
    return _model_nurbs_deboor_h(net,geometry.curve.degree-order,knots,parameter)
end

function _model_occ_curve_d2(geometry::_ModelOCCNURBSCurve,parameter)
    value=_model_occ_curve_homogeneous(geometry,0,parameter)
    first=_model_occ_curve_homogeneous(geometry,1,parameter)
    second=_model_occ_curve_homogeneous(geometry,2,parameter)
    isfinite(value[4]) && value[4]!=0.0 || throw(ArgumentError(
        "OCC Curve: homogeneous weight must be finite and nonzero"))
    point=ntuple(axis->value[axis]/value[4],3)
    derivative=ntuple(axis->(first[axis]-first[4]*point[axis])/value[4],3)
    second_derivative=ntuple(axis->(second[axis]-second[4]*point[axis]-
        2.0*first[4]*derivative[axis])/value[4],3)
    result=(point,derivative,second_derivative)
    all(vector->all(isfinite,vector),result) || throw(ArgumentError(
        "OCC Curve: evaluated derivatives must be finite"))
    return result
end

function _model_occ_curve_bounds(geometry::_ModelOCCNURBSCurve)
    curve=geometry.curve
    return curve.knots[curve.degree+1],curve.knots[end-curve.degree]
end

function _model_occ_curve_geometry(family,controls;degree=3,weights=Float64[],
                                   knots=Float64[],multiplicities=Int[],
                                   caller="OCC Curve")
    family in (:bezier,:bspline) || throw(ArgumentError("$caller: unknown curve family"))
    count=length(controls)
    count>=2 || throw(ArgumentError("$caller: at least two controls are required"))
    degree isa Integer && !(degree isa Bool) || throw(ArgumentError(
        "$caller: degree must be an integer"))
    actual_degree=family===:bezier ? count-1 : min(degree<=0 ? 3 : degree,count-1)
    actual_degree<=25 || throw(ArgumentError("$caller: curve degree exceeds 25"))
    expanded=family===:bezier ? vcat(zeros(count),ones(count)) :
        _model_nurbs_knots(count,actual_degree,knots,multiplicities,caller,"Curve")
    coefficients=isempty(weights) ? ones(Float64,count) :
        Float64[_finite_scalar(weight,caller,"weight") for weight in weights]
    length(coefficients)==count || throw(ArgumentError("$caller: weight count mismatch"))
    curve=NURBS.NURBSCurve(actual_degree,expanded,controls,coefficients)
    return _ModelOCCNURBSCurve(family,curve)
end

# Geom_{Bezier,BSpline}Curve::Resolution / BSplCLib::Resolution.
function _model_occ_curve_resolution(geometry::_ModelOCCNURBSCurve)
    curve=geometry.curve
    controls,weights,knots,degree=curve.controls,curve.weights,curve.knots,curve.degree
    count=length(controls)
    rational=any(weight->weight!=weights[1],weights)
    maximum_derivative=0.0
    for i in 1:count-1
        inverse=1.0/(knots[i+degree+1]-knots[i+1])
        current,previous=controls[i+1],controls[i]
        if rational
            for j in max(i-degree-1,0):min(2degree+1+i,count)-1
                pole=controls[j+1]
                value=0.0
                for axis in 1:3
                    factor=(pole[axis]-current[axis])*weights[i+1]-
                        (pole[axis]-previous[axis])*weights[i]
                    value+=abs(factor)
                end
                value*=inverse
                maximum_derivative<value && (maximum_derivative=value)
            end
        else
            value=0.0
            for axis in 1:3
                value+=abs(current[axis]-previous[axis])
            end
            value*=inverse
            maximum_derivative<value && (maximum_derivative=value)
        end
    end
    rational && (maximum_derivative/=minimum(weights))
    maximum_derivative*=degree
    denominator=maximum_derivative>floatmin(Float64) ? maximum_derivative : floatmin(Float64)
    return 1e-7*(1.0/denominator)
end

function _model_occ_curve_extrema_state(geometry,target,parameter)
    point,derivative,acceleration=_model_occ_curve_d2(geometry,parameter)
    residual=ntuple(axis->point[axis]-target[axis],3)
    value=(residual[1]*derivative[1]+residual[2]*derivative[2])+residual[3]*derivative[3]
    slope=((derivative[1]*derivative[1]+residual[1]*acceleration[1])+
        (derivative[2]*derivative[2]+residual[2]*acceleration[2]))+
        (derivative[3]*derivative[3]+residual[3]*acceleration[3])
    distance=(residual[1]*residual[1]+residual[2]*residual[2])+residual[3]*residual[3]
    return point,residual,value,slope,distance
end

# The interval solver retains the existing Newton computation. The surrounding
# Extrema_GExtPC admission/order is literal; math_FunctionRoots arithmetic is
# still a separate native parity requirement.
function _model_occ_curve_interval_extrema!(states,geometry,target,first,last,samples,tolerance)
    for seed in 0:samples
        parameter=first+(last-first)*(seed/samples)
        converged=false
        for iteration in 1:100
            _,_,value,slope,_=_model_occ_curve_extrema_state(geometry,target,parameter)
            slope==0.0 && break
            proposal=clamp(parameter-value/slope,first,last)
            if abs(proposal-parameter)<=2eps(Float64)*max(1.0,abs(parameter))
                parameter=proposal
                converged=true
                break
            end
            parameter=proposal
        end
        converged || continue
        point,residual,normal,_,distance=_model_occ_curve_extrema_state(geometry,target,parameter)
        derivative=_model_occ_curve_d2(geometry,parameter)[2]
        scale=max(1.0,sqrt(sum(value*value for value in residual)),
            sqrt(sum(value*value for value in derivative)))
        abs(normal)<=1e-10*scale*scale || continue
        finished_parameter=parameter
        any(state->abs(state[1]-finished_parameter)<=tolerance,states) && continue
        push!(states,(parameter,point,distance))
    end
    return nothing
end

function _model_occ_curve_project(geometry::_ModelOCCNURBSCurve,target,lower,upper)
    curve=geometry.curve
    tolerance=_model_occ_curve_resolution(geometry)
    states=Tuple{Float64,NTuple{3,Float64},Float64}[]
    if geometry.family===:bezier
        _model_occ_curve_interval_extrema!(states,geometry,target,lower,upper,2length(curve.controls),tolerance)
    else
        actual_lower,actual_upper=_model_occ_curve_bounds(geometry)
        knots=unique(Float64[knot for knot in curve.knots if actual_lower<=knot<=actual_upper])
        coefficient=(upper-lower)*1e-9
        first_used=1
        for i in eachindex(knots)
            lower>=knots[i]-coefficient || break
            first_used=i
        end
        last_used=length(knots)
        for i in length(knots):-1:1
            upper<=knots[i]+coefficient || break
            last_used=i
        end
        if first_used==last_used
            first_used=1
            last_used=2
        end
        first_state=_model_occ_curve_extrema_state(geometry,target,lower)
        last_state=_model_occ_curve_extrema_state(geometry,target,upper)
        square_confusion=1e-7*1e-7
        samples=curve.degree+1
        if curve.degree==1
            distance_min=floatmax(Float64)
            minimum_state=(0.0,(0.0,0.0,0.0),distance_min)
            previous_distance=0.0
            for i in first_used:last_used-1
                first,last=knots[i],knots[i+1]
                if i==first_used
                    first=lower
                elseif i==last_used-1
                    last=upper
                end
                first_point,first_residual,_,_,first_distance=
                    _model_occ_curve_extrema_state(geometry,target,first)
                last_point,last_residual,_,_,last_distance=
                    _model_occ_curve_extrema_state(geometry,target,last)
                chord=ntuple(axis->first_point[axis]-last_point[axis],3)
                value1=(chord[1]*first_residual[1]+chord[2]*first_residual[2])+chord[3]*first_residual[3]
                value2=(chord[1]*last_residual[1]+chord[2]*last_residual[2])+chord[3]*last_residual[3]
                if i!=first_used && distance_min>previous_distance
                    distance_min=previous_distance
                    minimum_state=(first,first_point,distance_min)
                end
                previous_distance=last_distance
                minimum_distance=min(first_distance,last_distance)
                minimum_derivative=min(abs(value1),abs(value2))
                if abs(value1)<1e100 && abs(value2)<1e100 &&
                    (value1*value2<=0.0 || minimum_distance<100.0*square_confusion ||
                     2.0*minimum_derivative<1e-7)
                    _model_occ_curve_interval_extrema!(states,geometry,target,first,last,samples,tolerance)
                end
            end
            finished_distance=distance_min
            finished_minimum=minimum_state
            if finished_distance<1e100 && all(state->finished_distance<state[3] &&
                abs(state[1]-finished_minimum[1])>tolerance,states)
                push!(states,minimum_state)
            end
        else
            parameters=Float64[]
            distances=Float64[]
            for i in first_used:last_used-1
                first=i==first_used ? lower : knots[i]
                last=i==last_used-1 ? upper : knots[i+1]
                step=(last-first)/samples
                for sample in 0:samples-1
                    parameter=first+step*sample
                    push!(parameters,parameter)
                    push!(distances,_model_occ_curve_extrema_state(geometry,target,parameter)[5])
                end
            end
            push!(parameters,upper)
            push!(distances,last_state[5])
            for i in 2:length(parameters)-1
                if distances[i]<=square_confusion
                    point=_model_occ_curve_d2(geometry,parameters[i])[1]
                    push!(states,(parameters[i],point,distances[i]))
                end
                if (distances[i]>=distances[i+1] && distances[i]>=distances[i-1]) ||
                    (distances[i]<=distances[i+1] && distances[i]<=distances[i-1])
                    _model_occ_curve_interval_extrema!(states,geometry,target,
                        parameters[i-1],parameters[i+1],samples,tolerance)
                end
            end
            for endpoint in (1,length(parameters))
                endpoint_state=endpoint==1 ? first_state : last_state
                distance=endpoint_state[5]
                distance>square_confusion && distance<1e100 || continue
                first_index=endpoint==1 ? 1 : length(parameters)-1
                first,last=parameters[first_index],parameters[first_index+1]
                _,residual1,value1,_,_=_model_occ_curve_extrema_state(geometry,target,first)
                _,residual2,value2,_,_=_model_occ_curve_extrema_state(geometry,target,last)
                product=(residual1[1]*residual2[1]+residual1[2]*residual2[2])+residual1[3]*residual2[3]
                value=endpoint==1 ? value1 : value2
                if abs(value1)<1e100 && abs(value2)<1e100 &&
                    (value1*value2<=0.0 || product<=0.0 || 2.0*abs(value)<1e-7)
                    _model_occ_curve_interval_extrema!(states,geometry,target,first,last,samples,tolerance)
                end
            end
        end
        if first_state[5]<square_confusion && !any(state->abs(state[1]-lower)<tolerance,states)
            pushfirst!(states,(lower,first_state[1],first_state[5]))
        end
        last_added=any(state->abs(state[1]-lower)>=tolerance &&
            abs(upper-state[1])<tolerance,states)
        if last_state[5]<square_confusion && !last_added
            push!(states,(upper,last_state[1],last_state[5]))
        end
    end
    best=nothing
    for state in states
        (best===nothing || state[3]<best[3]) && (best=state)
    end
    return best===nothing ? nothing : (best[1],best[2])
end

# ── converted-curve reparametrization ────────────────────────────────────────
# `convertToNURBS` replaces a conic record by a quadratic-rational NURBS
# whose public parameter map is preserved exactly: internal knots sit at the
# angular breakpoints and the tangent half-angle map
#   s = a + Δ/2·(1 + tan((t−m)/2)/tan(Δ/4)),   m = (a+b)/2
# carries the external angle t onto the internal B-spline parameter s.
# Because the external parametrization is unchanged, every stored pcurve and
# upstream/downstream parameter reference stays bitwise valid — stricter than
# the upstream behaviour, which refits 2-D curves approximately.
struct _ModelOCCConicReparam
    knots::Vector{Float64}              # span breakpoints (knot values)
    tanq::Vector{Float64}               # tan(Δ_i/4) per span
end

# Same field layout as `_ModelOCCNURBSCurve` plus the reparametrization data,
# so generic field access (`curve.curve.knots`, `curve.nets`, `curve.family`)
# keeps working for every consumer that does not evaluate.
struct _ModelOCCReparamCurve
    family::Symbol
    curve::NURBS.NURBSCurve
    nets::NTuple{3,Union{Nothing,Vector{NTuple{4,Float64}}}}
    reparam::_ModelOCCConicReparam
end

function _occ_conic_reparam_lookup(r::_ModelOCCConicReparam,t)
    index=clamp(searchsortedlast(r.knots,t),1,length(r.tanq))
    a=r.knots[index]
    span=r.knots[index+1]-a
    return a,span,a+span/2,r.tanq[index]
end

# t → (s, ds/dt, d²s/dt²) on the interior of each span; at a shared knot
# either adjacent span yields the same s and the composed curve stays C¹.
function _occ_conic_reparam(r::_ModelOCCConicReparam,t)
    a,span,mid,factor=_occ_conic_reparam_lookup(r,Float64(t))
    sx,cx=_gm_sincos((t-mid)/2)
    tangent=sx/cx
    ds=span/(4*factor*cx*cx)
    # At a span endpoint the `1+tangent/factor` cancellation is not exact in
    # floating point; pinning s inside its own span keeps de Boor's
    # span selection (`searchsortedlast` → the span STARTING at the knot)
    # consistent with the ds/dds factors computed here.
    s=clamp(a+span/2*(1+tangent/factor),a,a+span)
    return s,ds,ds*tangent
end

# s → t: the angle the internal parameter reaches, θ = m + 2·atan(u·tan(Δ/4)).
function _occ_conic_reparam_inv(r::_ModelOCCConicReparam,s)
    a,span,mid,factor=_occ_conic_reparam_lookup(r,Float64(s))
    return mid+2.0*_gm_atan((2*(s-a)/span-1)*factor)
end

function _model_occ_curve_d2(geometry::_ModelOCCReparamCurve,parameter)
    s,ds,dds=_occ_conic_reparam(geometry.reparam,Float64(parameter))
    point,first,second=_model_occ_curve_d2(
        _ModelOCCNURBSCurve(geometry.family,geometry.curve,geometry.nets),s)
    return (point,
            ntuple(axis->first[axis]*ds,3),
            ntuple(axis->first[axis]*dds+second[axis]*ds*ds,3))
end

# The internal knot range equals the external parameter range numerically.
_model_occ_curve_bounds(geometry::_ModelOCCReparamCurve)=
    _model_occ_curve_bounds(
        _ModelOCCNURBSCurve(geometry.family,geometry.curve,geometry.nets))
_model_occ_curve_resolution(geometry::_ModelOCCReparamCurve)=
    _model_occ_curve_resolution(
        _ModelOCCNURBSCurve(geometry.family,geometry.curve,geometry.nets))

function _model_occ_curve_project(geometry::_ModelOCCReparamCurve,target,
                                  lower,upper)
    inner=_ModelOCCNURBSCurve(geometry.family,geometry.curve,geometry.nets)
    result=_model_occ_curve_project(inner,target,
        _occ_conic_reparam(geometry.reparam,lower)[1],
        _occ_conic_reparam(geometry.reparam,upper)[1])
    result===nothing && return nothing
    return _occ_conic_reparam_inv(geometry.reparam,result[1]),result[2]
end
