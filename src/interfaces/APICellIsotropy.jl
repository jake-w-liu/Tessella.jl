# Private adaptive native minIsotropy provider for public higher-order queries.
# It uses ideal-gradient Bernstein products and rational lower bounds, unlike
# sampled minSICN, sampled minSIGE, or determinant coefficient minSJ.
struct _NativeBezierICNDomain
    determinant::Vector{Float64}
    matrix::Matrix{Float64}
    minL::Float64
    maxL::Float64
    minB::Float64
    maxB::Float64
end

const _NATIVE_CELL_ISOTROPY_RULES=Dict{Tuple{Int,Symbol,Int,Int,Int},_NativeBezierJacobianRule}()

function _native_bezier_raised_space(space,factor)
    if space.family===:pyr
        return _native_bezier_space(:pyr,factor*space.degree,factor*space.nij,factor*space.nk)
    end
    return _native_bezier_space(space.family,factor*space.degree)
end

function _native_bezier_product_weight(space,a,b,c=nothing)
    degree=space.degree;dimension=space.family===:tri || space.family===:qua ? 2 : 3
    simplex_dimension=space.family===:tri ? 2 : space.family===:tet ? 3 : 0
    factor=c===nothing ? 2 : 3
    remaining_a=degree;remaining_b=degree;remaining_c=degree;remaining_total=factor*degree
    numerator=1.0;denominator=1.0
    for axis in 1:dimension
        order=space.family===:pyr ? (axis==3 ? space.nk : space.nij) : degree
        na=axis<=simplex_dimension ? remaining_a : order
        nb=axis<=simplex_dimension ? remaining_b : order
        nc=axis<=simplex_dimension ? remaining_c : order
        nt=axis<=simplex_dimension ? remaining_total : factor*order
        # nChoosek returns C int; the source multiplies these before converting
        # the product to double. Keep that association and integer width.
        weight=_native_bezier_choose(na,a[axis])*_native_bezier_choose(nb,b[axis])
        if c!==nothing
            weight*=_native_bezier_choose(nc,c[axis])
        end
        total=a[axis]+b[axis]+(c===nothing ? 0 : c[axis])
        numerator*=Float64(weight);denominator*=Float64(_native_bezier_choose(nt,total))
        remaining_a-=a[axis];remaining_b-=b[axis]
        c===nothing || (remaining_c-=c[axis])
        remaining_total-=total
    end
    return numerator,denominator
end

function _native_bezier_square(space,coefficients)
    output_space=_native_bezier_raised_space(space,2)
    indices=Dict(exponent=>index for (index,exponent) in enumerate(output_space.exponents))
    result=zeros(Float64,length(output_space.exponents))
    for i in eachindex(coefficients),j in i:length(coefficients)
        a=space.exponents[i];b=space.exponents[j]
        numerator,denominator=_native_bezier_product_weight(space,a,b)
        i<j && (numerator*=2.0)
        target=indices[(a[1]+b[1],a[2]+b[2],a[3]+b[3])]
        result[target]+=(numerator/denominator)*coefficients[i]*coefficients[j]
    end
    return result
end

function _native_bezier_cube(space,coefficients)
    output_space=_native_bezier_raised_space(space,3)
    indices=Dict(exponent=>index for (index,exponent) in enumerate(output_space.exponents))
    result=zeros(Float64,length(output_space.exponents))
    for i in eachindex(coefficients),j in i:length(coefficients),k in j:length(coefficients)
        a=space.exponents[i];b=space.exponents[j];c=space.exponents[k]
        numerator,denominator=_native_bezier_product_weight(space,a,b,c)
        if i<j && j<k
            numerator*=6.0
        elseif i<j || j<k
            numerator*=3.0
        end
        target=indices[(a[1]+b[1]+c[1],a[2]+b[2]+c[2],a[3]+b[3]+c[3])]
        result[target]+=(numerator/denominator)*coefficients[i]*coefficients[j]*coefficients[k]
    end
    return result
end

function _native_bezier_rational_lower(numerator,denominator)
    length(numerator)==length(denominator) || throw(ArgumentError("native rational Bernstein sizes differ"))
    upper=floatmax(Float64);lower=-upper
    for position in eachindex(numerator)
        n=numerator[position];d=denominator[position]
        if d==0.0
            n<0.0 && return -floatmax(Float64)
        elseif d>0.0
            ratio=n/d;ratio<upper && (upper=ratio)
        else
            ratio=n/d;ratio>lower && (lower=ratio)
        end
    end
    return lower>upper ? -floatmax(Float64) : upper
end

@inline function _native_bezier_icn_ready(domain,minL,maxL)
    tolerance=1.e-3+(1.e-2-1.e-3)*_native_bezier_max(domain.minB,0.0)
    return minL-domain.minB<tolerance
end

function _native_bezier_icn_domain(detspace,gradspace,det,matrix,dimension)
    minL=floatmax(Float64);maxL=-minL
    for (detrow,gradrow) in zip(detspace.corners,gradspace.corners)
        norm=0.0
        for column in axes(matrix,2)
            norm+=matrix[gradrow,column]*matrix[gradrow,column]
        end
        value=dimension==2 ? 2.0*det[detrow]/norm :
            3.0*_gm_pow(det[detrow],2.0/3.0)/norm
        value<minL && (minL=value);value>maxL && (maxL=value)
    end
    lower=0.0
    tolerance=1.e-3+(1.e-2-1.e-3)*_native_bezier_max(lower,0.0)
    if !(minL-lower<tolerance) && all(value->!(value<0.0),det)
        denominator=if dimension==2
            values=zeros(Float64,length(det))
            for column in axes(matrix,2)
                squared=_native_bezier_square(gradspace,@view(matrix[:,column]))
                for position in eachindex(values)
                    values[position]+=squared[position]
                end
            end
            values
        else
            prox=zeros(Float64,size(matrix,1))
            for row in axes(matrix,1)
                for column in axes(matrix,2)
                    prox[row]+=matrix[row,column]*matrix[row,column]
                end
                prox[row]=sqrt(prox[row])
            end
            _native_bezier_cube(gradspace,prox)
        end
        fraction=_native_bezier_rational_lower(det,denominator)
        lower=dimension==2 ? 2.0*fraction : 3.0*_gm_pow(fraction*fraction,1.0/3.0)
    end
    return _NativeBezierICNDomain(det,matrix,minL,maxL,lower,0.0)
end

function _native_bezier_icn_extrema(detspace,gradspace,det,matrix,dimension)
    first=_native_bezier_icn_domain(detspace,gradspace,det,matrix,dimension)
    domains=_NativeBezierICNDomain[first];minL=first.minL;maxL=first.maxL
    _native_bezier_heap_make!(domains,false)
    iteration=0
    while !_native_bezier_icn_ready(domains[1],minL,maxL) && iteration+1<1000
        parent=domains[1];last=pop!(domains)
        isempty(domains) || _native_bezier_heap_adjust!(domains,1,length(domains),last,false)
        subd=_native_bezier_split(detspace,parent.determinant)
        subm=Vector{Matrix{Float64}}(undef,length(subd))
        for child in eachindex(subm)
            subm[child]=Matrix{Float64}(undef,size(parent.matrix)...)
        end
        for column in axes(parent.matrix,2)
            children=_native_bezier_split(gradspace,collect(@view(parent.matrix[:,column])))
            for child in eachindex(children),row in eachindex(children[child])
                subm[child][row,column]=children[child][row]
            end
        end
        for child in eachindex(subd)
            domain=_native_bezier_icn_domain(detspace,gradspace,subd[child],subm[child],dimension)
            minL=_native_bezier_min(minL,domain.minL);maxL=_native_bezier_max(maxL,domain.maxL)
            push!(domains,domain);_native_bezier_heap_push!(domains,domain,false)
        end
        iteration+=1
    end
    minB=domains[1].minB;minL=domains[1].minL
    for position in 2:length(domains)
        domain=domains[position]
        minB=_native_bezier_min(minB,domain.minB);minL=_native_bezier_min(minL,domain.minL)
    end
    factor=0.5*(minB+minL)
    return factor*minL+(1.0-factor)*minB
end

function _native_bezier_matrix_coefficients(space,samples)
    if space.family in (:tri,:tet)
        coefficients=zeros(Float64,size(samples)...)
        _native_bezier_product!(coefficients,space.inverse,samples)
        return coefficients
    end
    result=Matrix{Float64}(undef,size(samples)...)
    for column in axes(samples,2)
        coefficients=_native_bezier_coefficients(space,collect(@view(samples[:,column])))
        for row in eachindex(coefficients)
            result[row,column]=coefficients[row]
        end
    end
    return result
end

function _native_cell_isotropy_spaces(spec)
    p=spec.order;family=spec.family
    if family===:pyr
        return _native_bezier_space(family,p-1,p,p-1),
            _native_bezier_space(family,3p-3,3p,3p-3)
    end
    gradient_degree=family in (:tri,:tet) ? p-1 : p
    determinant_degree=p*spec.dim-(family===:tri ? 2 : family===:tet ? 3 : 0)
    return _native_bezier_space(family,gradient_degree),
        _native_bezier_space(family,determinant_degree)
end

function _native_cell_isotropy_rule_build(msh,space,caller)
    spec=msh_spec(msh);family=spec.family
    data=family===:pyr ? _native_bezier_pyramid_reference(msh,caller) : _native_cell_polynomial(msh,caller)
    gradients=ntuple(_->Matrix{Float64}(undef,length(space.points),spec.nnodes),3)
    for (sample,point) in enumerate(space.points)
        _,shape=family===:pyr ? _native_bezier_pyramid_shape(data,point) :
            _native_cell_query_shape(msh,point,data,caller)
        for node in 1:spec.nnodes,axis in 1:3
            gradients[axis][sample,node]=shape[node][axis]
        end
    end
    _native_cell_ideal_gradients!(gradients,family)
    primary=msh_type(family,1)
    primary_data=family===:pyr ? _native_bezier_pyramid_reference(primary,caller) : _native_cell_polynomial(primary,caller)
    _,primary_shape=family===:pyr ? _native_bezier_pyramid_shape(primary_data,(0.0,0.0,0.0)) :
        _native_cell_query_shape(primary,(0.0,0.0,0.0),primary_data,caller)
    primtuple=ntuple(_->Matrix{Float64}(undef,1,length(primary_shape)),3)
    for node in eachindex(primary_shape),axis in 1:3
        primtuple[axis][1,node]=primary_shape[node][axis]
    end
    _native_cell_ideal_gradients!(primtuple,family)
    primary_gradients=Matrix{Float64}(undef,3,length(primary_shape))
    for node in eachindex(primary_shape),axis in 1:3
        primary_gradients[axis,node]=primtuple[axis][1,node]
    end
    return _NativeBezierJacobianRule(space,gradients,primary_gradients,spec.dim)
end

function _native_cell_isotropy_rule(msh,space,caller)
    key=(Int(msh),space.family,space.degree,space.nij,space.nk)
    lock(_NATIVE_BEZIER_LOCK)
    try
        haskey(_NATIVE_CELL_ISOTROPY_RULES,key) && return _NATIVE_CELL_ISOTROPY_RULES[key]
        rule=_native_cell_isotropy_rule_build(msh,space,caller)
        _NATIVE_CELL_ISOTROPY_RULES[key]=rule
        return rule
    finally
        unlock(_NATIVE_BEZIER_LOCK)
    end
end

function _native_cell_quality_isotropy(mesh,cells,cell,msh,caller)
    spec=msh_spec(msh);_native_cell_numeric_contract(spec,caller)
    rule=_native_bezier_jacobian_rule(msh,caller)
    samples=_native_bezier_sample_jacobian(mesh,cells,cell,rule,false)
    coefficients=_native_bezier_coefficients(rule.space,samples)
    lower,upper=_native_bezier_extrema(rule.space,coefficients)
    ((lower<=0.0 && upper>=0.0) || upper<0.0) && return 0.0
    gradspace,detspace=_native_cell_isotropy_spaces(spec)
    detrule=_native_cell_isotropy_rule(msh,detspace,caller)
    detsamples=_native_bezier_sample_jacobian(mesh,cells,cell,detrule,false)
    if detsamples[1]<0.0
        for position in eachindex(detsamples)
            detsamples[position]*=-1.0
        end
    end
    det=_native_bezier_coefficients(detspace,detsamples)
    gradrule=_native_cell_isotropy_rule(msh,gradspace,caller)
    coordinates=Matrix{Float64}(undef,size(cells,1),3)
    for node in axes(coordinates,1)
        point=_mixed_node(mesh,cells,cell,node)
        for axis in 1:3
            coordinates[node,axis]=point[axis]
        end
    end
    allgrad=Matrix{Float64}(undef,length(gradspace.points),3*spec.dim)
    for axis in 1:spec.dim
        derivative=zeros(Float64,length(gradspace.points),3)
        _native_bezier_product!(derivative,gradrule.gradients[axis],coordinates)
        for component in 1:3,row in axes(allgrad,1)
            allgrad[row,3*(axis-1)+component]=derivative[row,component]
        end
    end
    matrix=_native_bezier_matrix_coefficients(gradspace,allgrad)
    return _native_bezier_icn_extrema(detspace,gradspace,det,matrix,spec.dim)
end
