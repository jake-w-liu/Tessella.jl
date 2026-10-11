# Normal-equation roots for an operation-owned copied tensor patch, using the
# pinned OCCT traversal and Julia-native arithmetic.
mutable struct _ModelNURBSProjectionFunction
    cache::_ModelNURBSProjectionCache
    target::NTuple{3,Float64}
    last_uv::Vector{Float64}
    last_xyz::NTuple{3,Float64}
end
function _ModelNURBSProjectionFunction(cache,target)
    return _ModelNURBSProjectionFunction(cache,Tuple(target),Float64[],(NaN,NaN,NaN))
end
function _model_nurbs_root_dot(a,b)
    value=0.0
    for i in eachindex(a)
        value+=a[i]*b[i]
    end
    return value
end
_model_nurbs_root_distance(a,b)=_model_nurbs_root_dot(a.-b,a.-b)

function _model_nurbs_root_value!(function_value,uv)
    point,du,dv=_model_nurbs_cache_evaluate(function_value.cache,uv[1],uv[2],1)
    residual=point.-function_value.target
    function_value.last_uv=copy(uv);function_value.last_xyz=point
    return [_model_nurbs_root_dot(residual,du),_model_nurbs_root_dot(residual,dv)]
end
function _model_nurbs_root_values!(function_value,uv)
    point,du,dv,duu,dvv,duv=_model_nurbs_cache_evaluate(function_value.cache,uv[1],uv[2],2)
    residual=point.-function_value.target
    function_value.last_uv=copy(uv);function_value.last_xyz=point
    ff=[_model_nurbs_root_dot(residual,du),_model_nurbs_root_dot(residual,dv)]
    any(value->value<=-1e100 || value>=1e100,ff) && return nothing
    a=_model_nurbs_root_dot(du,du)+_model_nurbs_root_dot(residual,duu)
    b=_model_nurbs_root_dot(dv,du)+_model_nurbs_root_dot(residual,duv)
    c=_model_nurbs_root_dot(dv,dv)+_model_nurbs_root_dot(residual,dvv)
    jacobian=[a b;b c]
    gh=zeros(2)
    for i in 1:2,j in 1:2
        gh[i]+=jacobian[j,i]*ff[j]
    end
    any(value->value<=-2e100 || value>=2e100,gh) && return nothing
    return ff,jacobian,gh,0.5*_model_nurbs_root_dot(ff,ff),_model_nurbs_root_dot(gh,gh)
end

struct _ModelNURBSProjectionLine
    function_value::_ModelNURBSProjectionFunction
    origin::Vector{Float64}
    direction::Vector{Float64}
end
function (line::_ModelNURBSProjectionLine)(parameter)
    uv=[line.direction[i]*parameter+line.origin[i] for i in 1:2]
    ff=_model_nurbs_root_value!(line.function_value,uv)
    any(value->value<=-1e100 || value>=1e100,ff) && return nothing
    return 0.5*_model_nurbs_root_dot(ff,ff)
end

function _model_nurbs_root_brent(function_value,ax,bx,cx,tolerance)
    a=min(ax,cx);b=max(ax,cx);x=bx;w=bx;v=bx
    initial=function_value(x);initial===nothing && return nothing
    fx=initial;fw=initial;fv=initial;e=0.0;d=floatmax(Float64)
    for _ in 1:100
        middle=0.5*(a+b);tol1=tolerance*abs(x)+tolerance;tol2=2.0*tol1
        x<=tol2+a && x>=b-tol2 && return (x,fx)
        if abs(e)>tol1
            r=(x-w)*(fx-fv);q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r;q=2.0*(q-r)
            q>0.0 && (p=-p)
            q=abs(q);previous_e=e;e=d
            if abs(p)>=abs(0.5*q*previous_e) || p<=q*(a-x) || p>=q*(b-x)
                e=x>=middle ? a-x : b-x;d=0.3819660*e
            else
                d=p/q;u=x+d
                if u-a<tol2 || b-u<tol2
                    d=_model_nurbs_root_sign(tol1,middle-x)
                end
            end
        else
            e=x>=middle ? a-x : b-x;d=0.3819660*e
        end
        u=abs(d)>=tol1 ? x+d : x+_model_nurbs_root_sign(tol1,d)
        evaluated=function_value(u);evaluated===nothing && return nothing
        fu=evaluated
        if fu<=fx
            u>=x ? (a=x) : (b=x)
            v,w,x=w,x,u;fv,fw,fx=fw,fx,fu
        else
            u<x ? (a=u) : (b=u)
            if fu<=fw || w==x
                v,w=w,u;fv,fw=fw,fu
            elseif fu<=fv || v==x || v==w
                v=u;fv=fu
            end
        end
    end
    return nothing
end

function _model_nurbs_root_minimize2(function_value,point,direction,
                                     point_value,direction_value,gradient,dgradient,tolerance)
    (point_value<=-2e100 || point_value>=2e100 ||
     direction_value<=-2e100 || direction_value>=2e100) && return false,direction
    eps=1e-20;tol1d=1.1
    for d in direction
        abs(d)>eps && (tol1d=min(tol1d,tolerance/abs(d)))
    end
    tol1d>0.9 && return false,direction
    df1=_model_nurbs_root_dot(gradient,direction)
    df2=_model_nurbs_root_dot(dgradient,direction)
    if df1 < -eps && df2>eps
        solution=-df1/(df2-df1)
    else
        c=point_value;b=df1;a=direction_value-(b+c)
        if abs(a)<=eps
            solution=abs(b)>=eps ? -c/b : 0.0
        else
            delta=b*b-4.0*a*c
            if delta>1e-9
                delta=sqrt(delta);solution=-(b+delta);other=delta-b
                abs(other)<abs(solution) && (solution=other)
                solution/=2.0*a
            else
                solution=-(0.5*b)/a
            end
        end
    end
    abs(solution)>=1.0 && return false,direction
    line=_ModelNURBSProjectionLine(function_value,copy(point),copy(direction))
    evaluated=line(solution);evaluated===nothing && return false,direction
    value=evaluated;good=value<point_value;result=good ? value : point_value
    if value>0.2*point_value && tol1d<0.5
        ax,bx,cx=solution<0.0 ? (solution,0.0,1.0) : (0.0,solution,1.0)
        answer=_model_nurbs_root_brent(line,ax,bx,cx,tol1d)
        if answer!==nothing && answer[2]<=result
            solution,result=answer;good=true
            if _model_nurbs_root_dot(gradient,gradient)>1e14 && ax<solution<cx
                answer=_model_nurbs_root_brent(line,ax,0.5*(ax+solution),solution,tol1d)
                answer!==nothing && answer[2]<=result && ((solution,result)=answer)
                answer=_model_nurbs_root_brent(line,solution,0.5*(cx+solution),cx,tol1d)
                answer!==nothing && answer[2]<=result && ((solution,result)=answer)
            end
        end
    end
    return good,good ? direction.*solution : direction
end

function _model_nurbs_root_minimize3(function_value,point0,point1,point2,
                                     value1,direction,tolerance)
    tol1d=2.1
    for d in direction
        abs(d)>1e-16 && (tol1d=min(tol1d,tolerance/abs(d)))
    end
    tol1d>1.9 && return false,direction
    tol1d/=3.0;new_direction=point1.-point0
    length=sqrt(_model_nurbs_root_dot(new_direction,new_direction))
    length<=1e-16 && return false,new_direction
    inverse=1.0/length
    difference=point2.-point1
    cx=sqrt(_model_nurbs_root_dot(difference,difference))*inverse
    cx<1e-2 && return false,new_direction
    line=_ModelNURBSProjectionLine(function_value,copy(point1),new_direction)
    answer=_model_nurbs_root_brent(line,-1.0,0.0,cx,tol1d)
    return answer!==nothing && answer[2]<value1 ?
        (true,new_direction.*answer[1]) : (false,new_direction)
end

function _model_nurbs_root_direction(jacobian,gh,ff,gradient,inverse,previous,constraints=nothing)
    free=Int[i for i in 1:2 if constraints===nothing || iszero(constraints[i])]
    direction=zeros(2)
    isempty(free) && return direction,0.0
    if !gradient
        if length(free)==1
            i=free[1]
            denominator=0.0;rhs=0.0
            for j in 1:2
                denominator+=jacobian[j,i]*jacobian[j,i]
                rhs+=jacobian[j,i]*(-ff[j])
            end
            if abs(denominator)>1e-20
                direction[i]=rhs/denominator
            else
                gradient=true
                direction[i]=previous[i]
            end
        else
            solution=_model_nurbs_root_gauss(jacobian,-ff,1e-9)
            solution===nothing && (solution=_model_nurbs_root_svd(jacobian,-ff))
            if solution===nothing
                gradient=true;direction=copy(previous)
            else
                direction=solution
            end
        end
    else
        direction=copy(previous)
    end
    ratio=abs(direction[free[1]]*inverse[free[1]])
    for i in free[2:end]
        ratio=max(ratio,abs(direction[i]*inverse[i]))
    end
    ratio>1.0 && (direction./=ratio)
    dy=0.0
    for i in free
        dy+=direction[i]*gh[i]
    end
    dy>=-1e-32 && (gradient=true)
    if gradient
        dy=0.0
        for i in free
            direction[i]=-gh[i];dy+=gh[i]*gh[i]
        end
        dy=-dy
    end
    if constraints!==nothing
        for i in 1:2
            !iszero(constraints[i]) && (direction[i]=0.0)
        end
    end
    return direction,dy
end

function _model_nurbs_root_bounds(solution,saved,lower,upper,tolerance)
    constraints=zeros(Int,2);delta=solution.-saved;ratio=1.0;out=false
    for i in 1:2
        if lower[i]==upper[i]
            constraints[i]=1;out=true
        elseif solution[i]<lower[i]
            constraints[i]=1;out=true
            -delta[i]>tolerance && (ratio=min(ratio,(lower[i]-saved[i])/delta[i]))
        elseif solution[i]>upper[i]
            constraints[i]=1;out=true
            delta[i]>tolerance && (ratio=min(ratio,(upper[i]-saved[i])/delta[i]))
        end
    end
    new=true
    if out
        if ratio==0.0
            new=false;solution=copy(saved);delta=zeros(2)
        else
            delta=delta.*ratio;solution=saved.+delta
            for i in 1:2
                if solution[i]<lower[i]
                    solution[i]=lower[i];delta[i]=solution[i]-saved[i]
                elseif solution[i]>upper[i]
                    solution[i]=upper[i];delta[i]=solution[i]-saved[i]
                end
            end
        end
    end
    return solution,delta,constraints,out,new
end

function _model_nurbs_projection_root!(function_value,start,lower,upper;tolerance=1e-9)
    inverse=[1.0/max((min(upper[i],2e100)-max(lower[i],-2e100))/4.0,1e-9) for i in 1:2]
    solution=[max(lower[i],min(upper[i],start[i])) for i in 1:2]
    initial=_model_nurbs_root_values!(function_value,solution)
    # Extrema records the last evaluated state even when the initial normal
    # equations exceed the native guard; that state still reaches its collector.
    initial===nothing && return function_value
    ff,jacobian,gh,f2,gnr=initial
    lambda2=gnr;save=[max(f2,1e-16)];tol_function=nextfloat(f2)-f2;gradient=false
    (f2<=1e-32 || gnr<=1e-64) && return function_value
    for count in 1:100
        minimum=f2;old_gradient=gnr;previous=copy(solution);saved=copy(solution)
        direction,dy=_model_nurbs_root_direction(jacobian,gh,ff,gradient,inverse,zeros(2))
        if abs(dy)<=1e-32
            _model_nurbs_root_value!(function_value,solution);return function_value
        end
        if gradient
            lambda=min(lambda2/sqrt(abs(dy)),1.0)
        else
            lambda=1.0;lambda2=0.5/sqrt(_model_nurbs_root_dot(direction,direction))
        end
        solution=solution.+direction.*lambda
        solution,delta,constraints,out,new=_model_nurbs_root_bounds(solution,saved,lower,upper,tolerance)
        gradient_save=copy(gh)
        if new
            values=_model_nurbs_root_values!(function_value,solution)
            values===nothing && return function_value
            ff,jacobian,gh,f2,gnr=values
        end
        if f2<=1e-32 || gnr<=1e-64
            _model_nurbs_root_value!(function_value,solution);return function_value
        end
        if out || f2/minimum>0.005
            dy=_model_nurbs_root_dot(gh,direction);old_f=minimum;stop=false;descent=0
            if !out
                while f2/minimum>0.005 && !stop
                    if f2<old_f && dy<0.0
                        descent+=1;saved=copy(solution);old_f=f2
                        solution=solution.+direction.*lambda
                        solution,delta,constraints,stop,new=_model_nurbs_root_bounds(solution,saved,lower,upper,tolerance)
                        lambda*=1.7
                    else
                        if f2>=old_f || f2>=minimum
                            good=false
                            if descent==0
                                descent+=1
                                good,delta=_model_nurbs_root_minimize2(function_value,saved,delta,old_f,f2,gradient_save,gh,tolerance)
                            elseif gradient || descent>1 || old_f>minimum
                                descent+=1
                                good,delta=_model_nurbs_root_minimize3(function_value,previous,saved,solution,old_f,delta,tolerance)
                            end
                            if !good
                                solution=copy(saved);f2=old_f
                            else
                                solution=saved.+delta
                                solution,delta,constraints,out,new=_model_nurbs_root_bounds(solution,saved,lower,upper,tolerance)
                            end
                            out=false
                        end
                        stop=true
                    end
                    gradient_save=copy(gh)
                    if new
                        values=_model_nurbs_root_values!(function_value,solution)
                        values===nothing && return function_value
                        ff,jacobian,gh,f2,gnr=values
                    end
                    dy=_model_nurbs_root_dot(gh,direction)
                    if abs(dy)<=1e-32
                        f2>old_f && (solution=copy(saved))
                        _model_nurbs_root_value!(function_value,solution);return function_value
                    end
                    descent>=100 && (stop=true)
                end
            end
            if out
                stop=f2>1.001*old_f;out_bis=out;descent=0
                while out_bis && (f2<old_f || descent==0) && !stop
                    descent+=1;saved=copy(solution);old_f=f2
                    direction,dy=_model_nurbs_root_direction(jacobian,gh,ff,gradient,inverse,direction,constraints)
                    if dy < -1e-32
                        if gradient
                            lambda=min(lambda2/sqrt(-dy),1.0)
                        else
                            lambda=1.0;lambda2=0.5/sqrt(_model_nurbs_root_dot(direction,direction))
                        end
                        solution=solution.+direction.*lambda
                        solution,delta,constraints,out_bis,new=_model_nurbs_root_bounds(solution,saved,lower,upper,tolerance)
                        gradient_save=copy(gh)
                        if new
                            values=_model_nurbs_root_values!(function_value,solution)
                            values===nothing && return function_value
                            ff,jacobian,gh,f2,gnr=values
                        end
                        lambda2=gnr
                    else
                        stop=true
                    end
                    while f2/minimum>0.005 && f2<old_f && !stop
                        descent+=1
                        if f2<old_f && dy<0.0
                            saved=copy(solution);old_f=f2
                            solution=solution.+direction.*lambda
                            solution,delta,constraints,out_bis,new=_model_nurbs_root_bounds(solution,saved,lower,upper,tolerance)
                        end
                        gradient_save=copy(gh)
                        if new
                            values=_model_nurbs_root_values!(function_value,solution)
                            values===nothing && return function_value
                            ff,jacobian,gh,f2,gnr=values
                        end
                        lambda2=gnr;dy=_model_nurbs_root_dot(gh,direction)
                        stop=dy>=0.0 || descent>=10 || out_bis
                    end
                    stop=dy>=0.0 || descent>=10
                end
                if (f2/minimum>0.005 && f2>=old_f) || f2>=minimum
                    descent+=1
                    good,delta=_model_nurbs_root_minimize2(function_value,saved,delta,old_f,f2,gradient_save,gh,tolerance)
                    if !good
                        solution=copy(saved);out=false
                    else
                        solution=saved.+delta
                        solution,delta,constraints,out,new=_model_nurbs_root_bounds(solution,saved,lower,upper,tolerance)
                        if new
                            values=_model_nurbs_root_values!(function_value,solution)
                            values===nothing && return function_value
                            ff,jacobian,gh,f2,gnr=values
                        end
                    end
                    dy=_model_nurbs_root_dot(gh,direction)
                end
            end
        end
        push!(save,f2)
        verify=gradient || (count>1 ? save[count]<1e-4*save[count-1] : f2<1e-6*save[1])
        if verify && all(abs.(previous.-solution).<=tolerance)
            minimum<f2 && (solution=copy(saved))
            _model_nurbs_root_value!(function_value,solution);return function_value
        end
        if f2-minimum<=tol_function
            if count>5
                if f2>=0.95*save[count-4]
                    gradient ? (return function_value) : (gradient=true)
                else
                    gradient=false
                end
            else
                gradient=false
            end
            gnr>0.9*old_gradient && f2>0.5*minimum && (gradient=true)
            if !gradient && !verify && all(abs.(previous.-solution).<=tolerance)
                _model_nurbs_root_value!(function_value,solution);return function_value
            end
        elseif !gradient
            gradient=true;solution=copy(previous)
            values=_model_nurbs_root_values!(function_value,solution)
            values===nothing && return function_value
            ff,jacobian,gh,f2,gnr=values
        else
            return function_value
        end
    end
    return function_value
end
