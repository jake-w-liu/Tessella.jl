# Operation-owned OCCT-style polynomial cache for a copied nonperiodic tensor
# patch. This evaluator is private to projection; public CAD evaluation retains
# its existing de Boor path and the mesh producer's numeric policy is unchanged.
mutable struct _ModelNURBSProjectionCache
    patch::_ModelNURBSSurface
    spans::NTuple{2,Int}
    centers::NTuple{2,Float64}
    half_spans::NTuple{2,Float64}
    swapped::Bool
    global_rational::Bool
    coefficients::Vector{Vector{Float64}}
end

function _ModelNURBSProjectionCache(patch::_ModelNURBSSurface)
    weights=patch.surface.weights
    first_weight=weights[1,1]
    rational=any(weight->weight!=first_weight,weights)
    return _ModelNURBSProjectionCache(patch,(0,0),(0.0,0.0),(0.0,0.0),
        false,rational,Vector{Float64}[])
end

function _model_nurbs_cache_bohm(poles::Vector{Vector{Float64}},degree::Int,
                                 knots,parameter::Float64)
    out=copy.(poles)
    twice_degree_minus_i=2*degree+1
    for i in 1:degree
        twice_degree_minus_i-=1
        upper_knot=twice_degree_minus_i
        for j in degree:-1:i
            upper_knot-=1
            denominator=knots[upper_knot+1]-knots[j]
            coefficient=iszero(denominator) ? 0.0 : 1.0/denominator
            for component in eachindex(out[j])
                out[j+1][component]=(out[j+1][component]-out[j][component])*coefficient
            end
        end
    end
    for i in 1:degree
        coefficient=parameter-knots[i]
        for j in i:-1:1,component in eachindex(out[j])
            out[j][component]=out[j][component]+coefficient*out[j+1][component]
        end
    end
    coefficient=Float64(degree)
    remaining=degree
    for i in 1:degree
        for component in eachindex(out[i+1])
            out[i+1][component]*=coefficient
        end
        remaining-=1
        coefficient*=remaining
    end
    return out
end

function _model_nurbs_cache_polynomial(coefficients,parameter::Float64,derivatives::Int)
    width=length(coefficients[1])
    out=Vector{Float64}[copy(coefficients[end])]
    for _ in 1:derivatives
        push!(out,zeros(width))
    end
    for row in length(coefficients)-1:-1:1
        coefficient=coefficients[row]
        # PLib's order-two specialization evaluates the second derivative
        # before the first derivative and value; it does not rescale afterward.
        for component in 1:width
            derivatives>=2 && (out[3][component]=
                out[3][component]*parameter+out[2][component]*2.0)
            derivatives>=1 && (out[2][component]=
                out[2][component]*parameter+out[1][component])
            out[1][component]=out[1][component]*parameter+coefficient[component]
        end
    end
    return out
end

function _model_nurbs_cache_prepare!(cache::_ModelNURBSProjectionCache,u::Float64,v::Float64)
    surface=cache.patch.surface
    pu,pv=surface.degree_u,surface.degree_v
    U,V=surface.knots_u,surface.knots_v
    spans=(clamp(searchsortedlast(U,u),pu+1,length(U)-pu-1),
           clamp(searchsortedlast(V,v),pv+1,length(V)-pv-1))
    spans==cache.spans && return cache
    centers=(U[spans[1]]+0.5*(U[spans[1]+1]-U[spans[1]]),
             V[spans[2]]+0.5*(V[spans[2]+1]-V[spans[2]]))
    half_spans=(0.5*(U[spans[1]+1]-U[spans[1]]),
                0.5*(V[spans[2]+1]-V[spans[2]]))
    starts=(spans[1]-pu,spans[2]-pv)
    swapped=pu>pv
    min_axis=swapped ? 2 : 1
    max_axis=3-min_axis
    min_degree,max_degree=swapped ? (pv,pu) : (pu,pv)
    knots=(U,V)
    min_knots=knots[min_axis][spans[min_axis]-min_degree+1:spans[min_axis]+min_degree]
    max_knots=knots[max_axis][spans[max_axis]-max_degree+1:spans[max_axis]+max_degree]
    first_weight=surface.weights[starts...]
    weight_epsilon=nextfloat(first_weight)-first_weight
    local_rational=false
    for j in 0:pv,i in 0:pu
        if abs(surface.weights[starts[1]+i,starts[2]+j]-first_weight)>weight_epsilon
            local_rational=true
            break
        end
    end
    width=local_rational ? 4 : 3
    poles=Vector{Float64}[]
    for a in 0:min_degree
        row=Float64[]
        for b in 0:max_degree
            i,j=swapped ? (starts[1]+b,starts[2]+a) : (starts[1]+a,starts[2]+b)
            point=surface.controls[i,j]
            if local_rational
                weight=surface.weights[i,j]
                append!(row,(point[1]*weight,point[2]*weight,point[3]*weight,weight))
            else
                append!(row,point)
            end
        end
        push!(poles,row)
    end
    min_poles=_model_nurbs_cache_bohm(poles,min_degree,min_knots,centers[min_axis])
    transformed=Vector{Float64}[]
    for row in min_poles
        pieces=Vector{Float64}[row[b*width+1:(b+1)*width] for b in 0:max_degree]
        max_poles=_model_nurbs_cache_bohm(pieces,max_degree,max_knots,centers[max_axis])
        push!(transformed,reduce(vcat,max_poles))
    end
    coefficients=Vector{Float64}[]
    max_factor=1.0
    for b in 0:max_degree
        row=Float64[]
        min_factor=1.0
        for a in 0:min_degree
            factor=min_factor*max_factor
            for component in 1:width
                push!(row,transformed[a+1][b*width+component]*factor)
            end
            if cache.global_rational && !local_rational
                push!(row,a==0 && b==0 ? 1.0 : 0.0)
            end
            min_factor*=half_spans[min_axis]/(a+1)
        end
        push!(coefficients,row)
        max_factor*=half_spans[max_axis]/(b+1)
    end
    cache.spans=spans
    cache.centers=centers
    cache.half_spans=half_spans
    cache.swapped=swapped
    cache.coefficients=coefficients
    return cache
end

function _model_nurbs_cache_evaluate(cache::_ModelNURBSProjectionCache,
                                      u::Float64,v::Float64,order::Int)
    0<=order<=2 || throw(ArgumentError("NURBS projection cache: derivative order must be 0, 1 or 2"))
    _model_nurbs_cache_prepare!(cache,u,v)
    inverse=(1.0/cache.half_spans[1],1.0/cache.half_spans[2])
    parameter=order==0 ?
        ((u-cache.centers[1])/cache.half_spans[1],(v-cache.centers[2])/cache.half_spans[2]) :
        ((u-cache.centers[1])*inverse[1],(v-cache.centers[2])*inverse[2])
    min_parameter,max_parameter=cache.swapped ? reverse(parameter) : parameter
    width=cache.global_rational ? 4 : 3
    max_values=_model_nurbs_cache_polynomial(cache.coefficients,max_parameter,order)
    min_degree=cache.swapped ? cache.patch.surface.degree_v : cache.patch.surface.degree_u
    homogeneous=Dict{NTuple{2,Int},Vector{Float64}}()
    for b in 0:order
        pieces=Vector{Float64}[max_values[b+1][a*width+1:(a+1)*width] for a in 0:min_degree]
        min_values=_model_nurbs_cache_polynomial(pieces,min_parameter,order-b)
        for a in 0:order-b
            homogeneous[a,b]=min_values[a+1]
        end
    end
    if order==0
        point=homogeneous[0,0]
        return cache.global_rational ?
            (point[1]/point[4],point[2]/point[4],point[3]/point[4]) :
            (point[1],point[2],point[3])
    end
    result=Dict{NTuple{2,Int},NTuple{3,Float64}}()
    if cache.global_rational
        inverse_weight=1.0/homogeneous[0,0][4]
        # The homogeneous storage has min-axis derivatives within each
        # max-axis row. BSplSLib::RationalDerivative visits those rows first;
        # this preserves its subtraction order for mixed derivatives.
        zero_h=(0.0,0.0,0.0,0.0)
        for b in 0:order,a in 0:order
            h=get(homogeneous,(a,b),zero_h)
            value=Float64[h[1],h[2],h[3]]
            for pp in 0:b-1,qq in 0:a
                coefficient=binomial(b,pp)*binomial(a,qq)*get(homogeneous,(a-qq,b-pp),zero_h)[4]
                prior=result[qq,pp]
                for component in 1:3
                    value[component]-=coefficient*prior[component]
                end
            end
            for qq in 0:a-1
                coefficient=binomial(a,qq)*homogeneous[a-qq,0][4]
                prior=result[qq,b]
                for component in 1:3
                    value[component]-=coefficient*prior[component]
                end
            end
            result[a,b]=(value[1]*inverse_weight,value[2]*inverse_weight,value[3]*inverse_weight)
        end
    else
        for (key,value) in homogeneous
            result[key]=(value[1],value[2],value[3])
        end
    end
    u_slot=cache.swapped ? (0,1) : (1,0)
    v_slot=cache.swapped ? (1,0) : (0,1)
    point=result[0,0]
    du=ntuple(i->result[u_slot][i]*inverse[1],3)
    dv=ntuple(i->result[v_slot][i]*inverse[2],3)
    order==1 && return (point,du,dv)
    uu_slot=cache.swapped ? (0,2) : (2,0)
    vv_slot=cache.swapped ? (2,0) : (0,2)
    duu=ntuple(i->result[uu_slot][i]*(inverse[1]*inverse[1]),3)
    dvv=ntuple(i->result[vv_slot][i]*(inverse[2]*inverse[2]),3)
    duv=ntuple(i->result[1,1][i]*(inverse[1]*inverse[2]),3)
    return (point,du,dv,duu,dvv,duv)
end
