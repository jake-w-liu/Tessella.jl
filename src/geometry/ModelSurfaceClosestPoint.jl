# Two-variable L-BFGS with a safeguarded cubic line search, independently
# expressed from the published numerical methods (Liu/Nocedal; More/Thuente,
# DOI 10.1145/192115.192132). Native public controls bind the arithmetic,
# compatibility constants and stopping order. No restricted solver is linked.
@inline _model_closest_inner2(a,b)=(0.0+a[1]*b[1])+a[2]*b[2]
@inline _model_closest_length2(a)=sqrt(_model_closest_inner2(a,a))
@inline _model_closest_scale2(a,factor)=(a[1]*factor,a[2]*factor)

function _model_closest_interval2(best,other,trial,bracketed,minimum,maximum)
    bx,bf,bg=best;ox,of,og=other;tx,tf,tg=trial
    invalid=maximum<minimum || bg*(tx-bx)>=0.0
    if bracketed
        invalid |= tx<=min(bx,ox) || tx>=max(bx,ox)
    end
    invalid && return best,other,tx,bracketed,false
    opposed=tg*(bg/abs(bg))<0.0
    limit_to_interval=false
    if tf>bf
        theta=3.0*(bf-tf)/(tx-bx)+bg+tg
        magnitude=max(abs(theta),max(abs(bg),abs(tg)))
        root=magnitude*sqrt((theta/magnitude)*(theta/magnitude)-bg/magnitude*(tg/magnitude))
        tx<bx && (root=-root)
        fraction=(root-bg+theta)/(root-bg+root+tg)
        cubic=bx+fraction*(tx-bx)
        quadratic=bx+bg/((bf-tf)/(tx-bx)+bg)/2.0*(tx-bx)
        chosen=abs(cubic-bx)<abs(quadratic-bx) ? cubic : cubic+(quadratic-cubic)/2.0
        bracketed=true;limit_to_interval=true
    elseif opposed
        theta=3.0*(bf-tf)/(tx-bx)+bg+tg
        magnitude=max(abs(theta),max(abs(bg),abs(tg)))
        root=magnitude*sqrt((theta/magnitude)*(theta/magnitude)-bg/magnitude*(tg/magnitude))
        tx>bx && (root=-root)
        fraction=(root-tg+theta)/(root-tg+root+bg)
        cubic=tx+fraction*(bx-tx)
        secant=tx+tg/(tg-bg)*(bx-tx)
        chosen=abs(cubic-tx)>abs(secant-tx) ? cubic : secant
        bracketed=true
    elseif abs(tg)<abs(bg)
        theta=3.0*(bf-tf)/(tx-bx)+bg+tg
        magnitude=max(abs(theta),max(abs(bg),abs(tg)))
        discriminant=max(0.0,(theta/magnitude)*(theta/magnitude)-bg/magnitude*(tg/magnitude))
        root=magnitude*sqrt(discriminant)
        tx>bx && (root=-root)
        fraction=(root-tg+theta)/(root+(bg-tg)+root)
        cubic=fraction<0.0 && root!=0.0 ? tx+fraction*(bx-tx) : (tx>bx ? maximum : minimum)
        secant=tx+tg/(tg-bg)*(bx-tx)
        chosen=if bracketed
            abs(tx-cubic)<abs(tx-secant) ? cubic : secant
        else
            abs(tx-cubic)>abs(tx-secant) ? cubic : secant
        end
        limit_to_interval=true
    elseif bracketed
        theta=3.0*(tf-of)/(ox-tx)+og+tg
        magnitude=max(abs(theta),max(abs(og),abs(tg)))
        root=magnitude*sqrt((theta/magnitude)*(theta/magnitude)-og/magnitude*(tg/magnitude))
        tx>ox && (root=-root)
        fraction=(root-tg+theta)/(root-tg+root+og)
        chosen=tx+fraction*(ox-tx)
    else
        chosen=tx>bx ? maximum : minimum
    end
    if tf>bf
        other=trial
    else
        opposed && (other=best)
        best=trial
    end
    step=max(minimum,min(maximum,chosen))
    if bracketed && limit_to_interval
        guard=best[1]+0.66*(other[1]-best[1])
        step=other[1]>best[1] ? min(guard,step) : max(guard,step)
    end
    return best,other,step,bracketed,true
end

function _model_closest_search2(callback::F,position,value,gradient,direction,step,ceiling) where F
    lower_limit,upper_limit=1.0e-50,1.0e50
    step=max(lower_limit,min(upper_limit,step))
    initial_slope=_model_closest_inner2(gradient,direction)
    (initial_slope>=0.0 || step<=0.0) && return position,value,gradient,step,false,0
    best=(0.0,value,initial_slope);other=best
    bracketed=false;modified=true;valid_interpolation=true
    interval_width=upper_limit-lower_limit;prior_width=interval_width/0.5
    sufficient_slope=0.001*initial_slope
    for evaluation in 1:20
        lower,upper=bracketed ? (min(best[1],other[1]),max(best[1],other[1])) :
            (best[1],step+4.0*(step-best[1]))
        step=max(lower_limit,min(upper_limit,step))
        outside=bracketed && (step<=lower || step>=upper)
        too_narrow=bracketed && upper-lower<=5.0e-14*upper
        if outside || evaluation==20 || !valid_interpolation || too_narrow
            step=best[1]
        end
        candidate=(position[1]+step*direction[1],position[2]+step*direction[2])
        candidate_value,candidate_gradient=callback(candidate)
        if candidate_value>=ceiling
            candidate_value=ceiling;candidate_gradient=(0.0,0.0)
        end
        slope=_model_closest_inner2(candidate_gradient,direction)
        sufficient_value=value+step*sufficient_slope
        success=candidate_value<=sufficient_value && abs(slope)<=-0.4*initial_slope
        warning=(bracketed && (step<=lower || step>=upper)) || !valid_interpolation
        warning |= step==upper_limit && candidate_value<=sufficient_value && slope<=sufficient_slope
        warning |= step==lower_limit && (candidate_value>sufficient_value || slope>=sufficient_slope)
        warning |= evaluation==20 || too_narrow
        if success || warning
            return candidate,candidate_value,candidate_gradient,step,success,evaluation
        end
        if modified && candidate_value<=sufficient_value && slope>=0.001*initial_slope
            modified=false
        end
        trial=(step,candidate_value,slope)
        if modified && candidate_value<=best[2] && candidate_value>sufficient_value
            shifted_best=(best[1],best[2]-best[1]*sufficient_slope,best[3]-sufficient_slope)
            shifted_other=(other[1],other[2]-other[1]*sufficient_slope,other[3]-sufficient_slope)
            shifted_trial=(step,candidate_value-step*sufficient_slope,slope-sufficient_slope)
            shifted_best,shifted_other,step,bracketed,valid_interpolation=
                _model_closest_interval2(shifted_best,shifted_other,shifted_trial,bracketed,lower,upper)
            best=(shifted_best[1],shifted_best[2]+shifted_best[1]*sufficient_slope,shifted_best[3]+sufficient_slope)
            other=(shifted_other[1],shifted_other[2]+shifted_other[1]*sufficient_slope,shifted_other[3]+sufficient_slope)
        else
            best,other,step,bracketed,valid_interpolation=
                _model_closest_interval2(best,other,trial,bracketed,lower,upper)
        end
        if bracketed
            if abs(other[1]-best[1])>=0.66*prior_width
                step=best[1]+0.5*(other[1]-best[1])
            end
            prior_width,interval_width=interval_width,abs(other[1]-best[1])
        end
    end
    error("internal closest-point line search did not terminate")
end

function _model_closest_minimize2(callback::F,initial::NTuple{2,Float64}) where F
    position=initial;value,gradient=callback(position)
    ceiling=10.0*(abs(value)+1.0);evaluations=1
    _model_closest_length2(gradient)<=1.0e-12 && return position,4,0,evaluations
    direction=(-gradient[1],-gradient[2]);step=min(1.0/_model_closest_length2(gradient),1.0)
    newest=nothing;older=nothing
    for iteration in 1:500
        old_position,old_gradient,old_value=position,gradient,value
        newest===nothing || (step=1.0)
        largest=max(abs(direction[1]),abs(direction[2]))
        if largest!=0.0
            factor=1.0/largest;direction=_model_closest_scale2(direction,factor);step/=factor
            factor=1.0/_model_closest_length2(direction);direction=_model_closest_scale2(direction,factor);step/=factor
        end
        position,value,gradient,step,accepted,count=
            _model_closest_search2(callback,position,value,gradient,direction,step,ceiling)
        evaluations+=count
        displacement=(-old_position[1]+position[1],-old_position[2]+position[2])
        gradient_change=(-old_gradient[1]+gradient[1],-old_gradient[2]+gradient[2])
        iteration==500 && return position,5,iteration,evaluations
        _model_closest_length2(gradient)<=1.0e-12 && return position,4,iteration,evaluations
        old_value-value<=0.0*max(abs(old_value),max(abs(value),1.0)) && return position,1,iteration,evaluations
        _model_closest_length2(displacement)<=0.0 && return position,2,iteration,evaluations
        if !accepted
            direction=(-gradient[1],-gradient[2]);continue
        end
        curvature=_model_closest_inner2(gradient_change,displacement)
        gradient_squared=_model_closest_inner2(gradient_change,gradient_change)
        (curvature==0.0 || gradient_squared==0.0) && return position,-2,iteration,evaluations
        older=newest;newest=(displacement,gradient_change,1.0/curvature)
        work=gradient
        theta_new=_model_closest_inner2(newest[1],work)
        weight=theta_new*newest[3]
        work=(work[1]-weight*newest[2][1],work[2]-weight*newest[2][2])
        theta_old=0.0
        if older!==nothing
            theta_old=_model_closest_inner2(older[1],work);weight=theta_old*older[3]
            work=(work[1]-weight*older[2][1],work[2]-weight*older[2][2])
        end
        work=_model_closest_scale2(work,curvature/gradient_squared)
        if older!==nothing
            weight=older[3]*(-_model_closest_inner2(older[2],work)+theta_old)
            work=(work[1]+weight*older[1][1],work[2]+weight*older[1][2])
        end
        weight=newest[3]*(-_model_closest_inner2(newest[2],work)+theta_new)
        work=(work[1]+weight*newest[1][1],work[2]+weight*newest[1][2])
        direction=(-work[1],-work[2])
    end
    error("internal closest-point minimizer did not terminate")
end

@inline _model_closest_plane_value(plane,u,v)=ntuple(axis->
    plane.origin[axis]+plane.first_direction[axis]*u+plane.second_direction[axis]*v,3)

function _model_closest_plane_native(plane,coordinate,lower,upper)
    evaluate(u,v)=_model_closest_plane_value(plane,u,v)
    function distance(u,v)
        point=evaluate(u,v)
        dx=coordinate[1]-point[1];dy=coordinate[2]-point[2];dz=coordinate[3]-point[3]
        return sqrt(dx*dx+dy*dy+dz*dz)
    end
    guess=(0.0,0.0);best=distance(guess...)
    du=(upper[1]-lower[1])/10.0;dv=(upper[2]-lower[2])/10.0
    epsilon_u=1.0e-5*(upper[1]-lower[1]);epsilon_v=1.0e-5*(upper[2]-lower[2])
    # A nonzero finite interval is required for the native repeated-addition grid.
    all(isfinite,(du,dv,epsilon_u,epsilon_v)) && du>0.0 && dv>0.0 ||
        throw(ArgumentError("native Plane closest-point sampling bounds are degenerate"))
    u=lower[1]
    while u<=upper[1]+epsilon_u
        v=lower[2]
        while v<=upper[2]+epsilon_v
            candidate=distance(u,v)
            if candidate<best
                best=candidate;guess=(u,v)
            end
            v+=dv
        end
        u+=du
    end
    function callback(uv)
        point=evaluate(uv...)
        dx=coordinate[1]-point[1];dy=coordinate[2]-point[2];dz=coordinate[3]-point[3]
        # These public native objective and gradient contracts have different
        # coordinate weights. Replacing them with a consistent quadratic alters
        # the observed line-search decisions and final Float64 parameters.
        value=0.5*dx*dx+dy*dy+dz*dz
        first=plane.first_direction;second=plane.second_direction
        gradient=(-dx*first[1]-dy*first[2]-dz*first[3],
                  -dx*second[1]-dy*second[2]-dz*second[3])
        return value,gradient
    end
    uv,termination,iterations,evaluations=_model_closest_minimize2(callback,guess)
    return evaluate(uv...),uv,(termination=termination,iterations=iterations,evaluations=evaluations)
end

# The native built-in Surface provider evaluates its signed generatrices.
# Admission depends on the live CAD topology and existing curve evaluators.
function _model_high_order_ruled_surface(m,tag,caller)
    haskey(m.surfaces,tag) && !haskey(m.discrete,(2,tag)) || return false
    _surface_type(m,tag) in (:ruled,:tric) || return false
    loops=m.surfaces[tag]
    length(loops)==1 && haskey(m.loops,only(loops)) || return false
    generatrices=m.loops[only(loops)]
    length(generatrices) in (3,4) || return false
    return all(curve->_model_high_order_curve_geometry(m,abs(curve),caller)!==nothing,
               generatrices)
end

function _model_closest_ruled_native(m,tag,coordinate,caller;initial=nothing)
    _model_high_order_ruled_surface(m,tag,caller) || throw(ArgumentError(
        "$caller: Surface[$tag] lacks a supported built-in geometry provider"))
    evaluate(u,v)=_ruled_surface_point(m,tag,u,v,caller)
    function distance(u,v)
        p=evaluate(u,v)
        dx=coordinate[1]-p[1];dy=coordinate[2]-p[2];dz=coordinate[3]-p[3]
        return sqrt(dx*dx+dy*dy+dz*dz)
    end
    if initial===nothing
        guess=(0.0,0.0)
    else
        length(initial)==2 && all(isfinite,initial) || throw(ArgumentError(
            "$caller: initial Surface parameters must contain two finite values"))
        guess=(Float64(initial[1]),Float64(initial[2]))
    end
    best=distance(guess...)
    u=0.0
    while u<=1.0+1.0e-5
        v=0.0
        while v<=1.0+1.0e-5
            candidate=distance(u,v)
            if candidate<best
                best=candidate;guess=(u,v)
            end
            v+=0.1
        end
        u+=0.1
    end
    function callback(uv)
        p=evaluate(uv...)
        dx=coordinate[1]-p[1];dy=coordinate[2]-p[2];dz=coordinate[3]-p[3]
        first=_ruled_surface_d1(m,tag,uv[1],uv[2],1,caller)
        second=_ruled_surface_d1(m,tag,uv[1],uv[2],2,caller)
        return 0.5*dx*dx+dy*dy+dz*dz,
            (-dx*first[1]-dy*first[2]-dz*first[3],
             -dx*second[1]-dy*second[2]-dz*second[3])
    end
    uv,termination,iterations,evaluations=_model_closest_minimize2(callback,guess)
    xyz=evaluate(uv...)
    all(isfinite,xyz) && all(isfinite,uv) || throw(ArgumentError(
        "$caller: built-in Surface projection is non-finite"))
    return xyz,uv,(termination=termination,iterations=iterations,evaluations=evaluations)
end
