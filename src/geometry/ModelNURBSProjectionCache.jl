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
    # Work space of one evaluation: the value and derivative rows in the
    # longer-degree parameter, the homogeneous derivatives `(a,b)` and the
    # Cartesian ones. An evaluation used to build each of them as new lists
    # and dictionaries, 5.5 kB at a time.
    rows::Vector{Float64}
    homogeneous::Vector{Float64}
    cartesian::Vector{NTuple{3,Float64}}
end

function _ModelNURBSProjectionCache(patch::_ModelNURBSSurface)
    weights=patch.surface.weights
    first_weight=weights[1,1]
    rational=any(weight->weight!=first_weight,weights)
    return _ModelNURBSProjectionCache(patch,(0,0),(0.0,0.0),(0.0,0.0),
        false,rational,Vector{Float64}[],Float64[],zeros(36),
        fill((0.0,0.0,0.0),9))
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

# Value and derivatives up to `order` at `(u,v)`, from the polynomial of the
# current span. The arithmetic is `_model_nurbs_cache_polynomial`'s, operation
# for operation (PLib's Horner scheme, which updates the second derivative
# before the first and the value), in fixed work space instead of new lists:
# the results are the same to the last bit.
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
    coefficients=cache.coefficients
    count=length(coefficients)
    len=length(coefficients[1])
    min_degree=cache.swapped ? cache.patch.surface.degree_v : cache.patch.surface.degree_u
    rows=cache.rows
    length(rows)<3len && resize!(rows,3len)
    # Rows `k*len+1 : (k+1)*len` hold derivative `k` in the longer-degree
    # parameter, each a list of `min_degree+1` coefficients of `width` numbers.
    @inbounds begin
        highest=coefficients[count]
        for c in 1:len
            rows[c]=highest[c]
        end
        for c in len+1:(order+1)*len
            rows[c]=0.0
        end
        for row in count-1:-1:1
            coefficient=coefficients[row]
            for c in 1:len
                order>=2 && (rows[2len+c]=rows[2len+c]*max_parameter+rows[len+c]*2.0)
                order>=1 && (rows[len+c]=rows[len+c]*max_parameter+rows[c])
                rows[c]=rows[c]*max_parameter+coefficient[c]
            end
        end
        # Homogeneous derivative `(a,b)`, `a` in the shorter-degree parameter:
        # four numbers at `((b*3+a)*4+1 : ...)`. Only `a+b <= order` is set.
        homogeneous=cache.homogeneous
        for b in 0:order
            derivatives=order-b
            base=b*len
            for component in 1:width
                value=rows[base+min_degree*width+component];first=0.0;second=0.0
                for a in min_degree-1:-1:0
                    derivatives>=2 && (second=second*min_parameter+first*2.0)
                    derivatives>=1 && (first=first*min_parameter+value)
                    value=value*min_parameter+rows[base+a*width+component]
                end
                homogeneous[(b*3)*4+component]=value
                derivatives>=1 && (homogeneous[(b*3+1)*4+component]=first)
                derivatives>=2 && (homogeneous[(b*3+2)*4+component]=second)
            end
        end
        if order==0
            return cache.global_rational ?
                (homogeneous[1]/homogeneous[4],homogeneous[2]/homogeneous[4],
                 homogeneous[3]/homogeneous[4]) :
                (homogeneous[1],homogeneous[2],homogeneous[3])
        end
        # Cartesian derivative `(a,b)` at `b*3+a+1`.
        cartesian=cache.cartesian
        if cache.global_rational
            inverse_weight=1.0/homogeneous[4]
            # The homogeneous storage has min-axis derivatives within each
            # max-axis row. BSplSLib::RationalDerivative visits those rows first;
            # this preserves its subtraction order for mixed derivatives. A
            # derivative of total order above `order` reads as zero.
            for b in 0:order,a in 0:order
                present=a+b<=order
                at=(b*3+a)*4
                x=present ? homogeneous[at+1] : 0.0
                y=present ? homogeneous[at+2] : 0.0
                z=present ? homogeneous[at+3] : 0.0
                for pp in 0:b-1,qq in 0:a
                    weight=(a-qq)+(b-pp)<=order ?
                        homogeneous[((b-pp)*3+(a-qq))*4+4] : 0.0
                    coefficient=binomial(b,pp)*binomial(a,qq)*weight
                    prior=cartesian[pp*3+qq+1]
                    x-=coefficient*prior[1];y-=coefficient*prior[2];z-=coefficient*prior[3]
                end
                for qq in 0:a-1
                    coefficient=binomial(a,qq)*homogeneous[(a-qq)*4+4]
                    prior=cartesian[b*3+qq+1]
                    x-=coefficient*prior[1];y-=coefficient*prior[2];z-=coefficient*prior[3]
                end
                cartesian[b*3+a+1]=(x*inverse_weight,y*inverse_weight,z*inverse_weight)
            end
        else
            for b in 0:order,a in 0:order-b
                at=(b*3+a)*4
                cartesian[b*3+a+1]=(homogeneous[at+1],homogeneous[at+2],homogeneous[at+3])
            end
        end
        # `(1,0)` is the derivative in the shorter-degree parameter.
        first_u,first_v=cache.swapped ? (cartesian[4],cartesian[2]) :
                                        (cartesian[2],cartesian[4])
        point=cartesian[1]
        du=(first_u[1]*inverse[1],first_u[2]*inverse[1],first_u[3]*inverse[1])
        dv=(first_v[1]*inverse[2],first_v[2]*inverse[2],first_v[3]*inverse[2])
        order==1 && return (point,du,dv)
        second_u,second_v=cache.swapped ? (cartesian[7],cartesian[3]) :
                                          (cartesian[3],cartesian[7])
        mixed=cartesian[5]
        uu=inverse[1]*inverse[1];vv=inverse[2]*inverse[2];uv=inverse[1]*inverse[2]
        duu=(second_u[1]*uu,second_u[2]*uu,second_u[3]*uu)
        dvv=(second_v[1]*vv,second_v[2]*vv,second_v[3]*vv)
        duv=(mixed[1]*uv,mixed[2]*uv,mixed[3]*uv)
        return (point,du,dv,duu,dvv,duv)
    end
end
