# Private public-query Bernstein spaces. Certified GeometryMaps/FEM providers
# have their own numerical and enclosure contracts and do not use these helpers.
using ..GmshLibm: _gm_pow

struct _NativeBezierSpace
    family::Symbol
    degree::Int
    nij::Int
    nk::Int
    exponents::Vector{NTuple{3,Int}}
    points::Vector{NTuple{3,Float64}}
    corners::Vector{Int}
    inverse::Union{Nothing,Matrix{Float64}}
end

struct _NativeBezierJacobianRule
    space::_NativeBezierSpace
    gradients::NTuple{3,Matrix{Float64}}
    primary_gradients::Matrix{Float64}
    dimension::Int
end

const _NATIVE_BEZIER_SPACES=Dict{Tuple{Symbol,Int,Int,Int},_NativeBezierSpace}()
const _NATIVE_BEZIER_JACOBIANS=Dict{Int,_NativeBezierJacobianRule}()
const _NATIVE_BEZIER_LOCK=ReentrantLock()

# The native Windows runtime's pow dispatches finite integer exponents to
# binary repeated squaring. Bernstein exponents are nonnegative integers;
# applying a general real-exponent libm routine changes the interpolation
# matrix by ulps which its inverse can amplify substantially.
function _native_bezier_integer_power(value::Float64,exponent::Int)
    result=isodd(exponent) ? value : 1.0
    exponent>>=1
    while exponent>0
        value*=value
        isodd(exponent) && (result*=value)
        exponent>>=1
    end
    return result
end

function _native_bezier_choose(n,k)
    (n<k || k<0) && throw(ArgumentError("invalid native binomial coefficient"))
    k>n÷2 && (k=n-k)
    k==1 && return Int32(n)
    k==0 && return Int32(1)
    value=Int32(1)
    for index in 1:k
        value=div(value*Int32(n),Int32(index))
        n-=1
    end
    return value
end

function _native_bezier_space_build(family,degree,nij,nk)
    exponents=NTuple{3,Int}[]
    points=NTuple{3,Float64}[]
    if family===:tri
        for j in 0:degree,i in 0:degree-j
            push!(exponents,(i,j,0))
        end
    elseif family===:tet
        for k in 0:degree,j in 0:degree-k,i in 0:degree-j-k
            push!(exponents,(i,j,k))
        end
    elseif family===:qua
        for j in 0:degree,i in 0:degree
            push!(exponents,(i,j,0))
        end
    elseif family===:hex
        for k in 0:degree,j in 0:degree,i in 0:degree
            push!(exponents,(i,j,k))
        end
    elseif family===:pyr
        for k in nk:-1:0,j in 0:nij,i in 0:nij
            push!(exponents,(i,j,k))
        end
    else
        throw(ArgumentError("unsupported native Bernstein family $family"))
    end
    if family===:pyr
        divisor=max(nij,nk)
        for (i,j,k) in exponents
            z=Float64(nk-k)/divisor
            width=1.0-z;scale=2.0/nij*width
            push!(points,(-width+Float64(i)*scale,-width+Float64(j)*scale,z))
        end
    elseif degree==0
        push!(points,(0.0,0.0,0.0))
    else
        scale=family in (:tri,:tet) ? 1.0/degree : 2.0/degree
        for exponent in exponents
            if family in (:tri,:tet)
                push!(points,(Float64(exponent[1])*scale,Float64(exponent[2])*scale,
                    Float64(exponent[3])*scale))
            else
                push!(points,(Float64(exponent[1])*scale-1.0,Float64(exponent[2])*scale-1.0,
                    family===:qua ? 0.0 : Float64(exponent[3])*scale-1.0))
            end
        end
    end
    count=length(exponents);n=degree+1
    corners=if degree==0 && family!==:pyr
        Int[1]
    elseif family===:tri
        Int[1,degree+1,count]
    elseif family===:tet
        Int[1,degree+1,(degree+2)*(degree+1)÷2,count]
    elseif family===:qua
        Int[1,n,count,count-degree]
    elseif family===:hex
        Int[1,n,n*n,n*degree+1,n*n*degree+1,n*n*degree+degree+1,count,
            (degree+2)*(degree+1)*degree+1]
    else
        m=nij+1
        nk==0 ? Int[1,nij+1,m*m,m*nij+1] :
            Int[1,nij+1,m*m,m*nij+1,m*m*nk+1,m*m*nk+nij+1,count,m*m*nk+m*nij+1]
    end
    inverse=nothing
    if family in (:tri,:tet)
        matrix=Matrix{Float64}(undef,count,count)
        simplex_dimension=family===:tri ? 2 : 3
        scale=degree==0 ? 0.0 : 1.0/degree
        for coefficient in 1:count,sample in 1:count
            value=1.0;complement=1.0;remaining=degree
            for axis in 1:simplex_dimension
                point=Float64(exponents[sample][axis])*scale
                power=exponents[coefficient][axis]
                value*=Float64(_native_bezier_choose(remaining,power))*_native_bezier_integer_power(point,power)
                complement-=point;remaining-=power
            end
            value*=_native_bezier_integer_power(complement,remaining)
            matrix[sample,coefficient]=value
        end
        inverse=_native_cell_inverse(matrix)
        all(isfinite,inverse) || throw(ArgumentError("native Bernstein inverse is not finite"))
    end
    return _NativeBezierSpace(family,degree,nij,nk,exponents,points,corners,inverse)
end

function _native_bezier_space(family,degree,nij=0,nk=0)
    key=(family,degree,nij,nk)
    lock(_NATIVE_BEZIER_LOCK)
    try
        haskey(_NATIVE_BEZIER_SPACES,key) && return _NATIVE_BEZIER_SPACES[key]
        space=_native_bezier_space_build(family,degree,nij,nk)
        _NATIVE_BEZIER_SPACES[key]=space
        return space
    finally
        unlock(_NATIVE_BEZIER_LOCK)
    end
end

function _native_bezier_jacobian_build(msh,caller)
    spec=msh_spec(msh);family=spec.family;p=spec.order
    _native_cell_numeric_contract(spec,caller)
    family in (:tri,:qua,:tet,:hex,:pyr) || throw(ArgumentError(
        "$caller: native Bernstein determinant space is unavailable for $family"))
    degree=family===:tri ? 2p-2 : family===:qua ? 2p-1 : family===:tet || family===:pyr ? 3p-3 : 3p-1
    space=family===:pyr ? _native_bezier_space(family,degree,degree+2,degree) :
        _native_bezier_space(family,degree)
    gradients=ntuple(_->Matrix{Float64}(undef,length(space.points),spec.nnodes),3)
    data=family===:pyr ? _native_bezier_pyramid_reference(msh,caller) :
        _native_cell_polynomial(msh,caller)
    for (sample,point) in enumerate(space.points)
        _,shape_gradients=family===:pyr ? _native_bezier_pyramid_shape(data,point) :
            _native_cell_query_shape(msh,point,data,caller)
        for node in 1:spec.nnodes,axis in 1:3
            gradients[axis][sample,node]=shape_gradients[node][axis]
        end
    end
    primary=msh_type(family,1)
    primary_data=family===:pyr ? _native_bezier_pyramid_reference(primary,caller) :
        _native_cell_polynomial(primary,caller)
    _,primary_shape=family===:pyr ? _native_bezier_pyramid_shape(primary_data,(0.0,0.0,0.0)) :
        _native_cell_query_shape(primary,(0.0,0.0,0.0),primary_data,caller)
    primary_gradients=Matrix{Float64}(undef,3,length(primary_shape))
    for node in eachindex(primary_shape),axis in 1:3
        primary_gradients[axis,node]=primary_shape[node][axis]
    end
    return _NativeBezierJacobianRule(space,gradients,primary_gradients,spec.dim)
end

function _native_bezier_jacobian_rule(msh,caller)
    lock(_NATIVE_BEZIER_LOCK)
    try
        haskey(_NATIVE_BEZIER_JACOBIANS,msh) && return _NATIVE_BEZIER_JACOBIANS[msh]
        rule=_native_bezier_jacobian_build(msh,caller)
        _NATIVE_BEZIER_JACOBIANS[msh]=rule
        return rule
    finally
        unlock(_NATIVE_BEZIER_LOCK)
    end
end

@inline _native_bezier_det(a,b,c)=a[1]*(b[2]*c[3]-c[2]*b[3])-
    b[1]*(a[2]*c[3]-c[2]*a[3])+c[1]*(a[2]*b[3]-b[2]*a[3])

# Positive fullMatrix::mult products are distinct from the negative Schur
# updates used by the interpolation inverse. Eigen dispatches a one-column
# result to GEMV, whose independent column blocks accumulate into each row.
function _native_bezier_gemv!(target,left,right)
    rows,depth=size(left)
    block_columns=depth<128 ? depth : stride(left,2)*8<32000 ? 16 : 4
    @inbounds for first in 1:block_columns:depth
        last=min(depth,first+block_columns-1)
        for row in 1:rows
            value=0.0
            for k in first:last
                value+=left[row,k]*right[k,1]
            end
            target[row,1]+=value
        end
    end
    return target
end

function _native_bezier_gemm_block!(target,left,right)
    rows,columns=size(target);depth=size(left,2)
    tail=(rows÷4)*4+1;full_columns=(columns÷4)*4;peeled=(depth÷8)*8
    @inbounds for column in 1:columns,row in 1:rows
        value=0.0
        paired=rows%4>=2 && tail<=row<tail+2 && column<=full_columns
        if paired
            even=0.0;odd=0.0
            for k in 1:2:peeled
                even+=left[row,k]*right[k,column]
                odd+=left[row,k+1]*right[k+1,column]
            end
            value=even+odd
        end
        for k in (paired ? peeled+1 : 1):depth
            value+=left[row,k]*right[k,column]
        end
        target[row,column]+=value
    end
    return target
end

function _native_bezier_product!(target,left,right)
    isempty(target) && return target
    rows,columns=size(target);depth=size(left,2)
    columns==1 && return _native_bezier_gemv!(target,left,right)
    kc,mc,nc=_native_cell_blocking(depth,rows,columns,1,_native_cell_cache_sizes())
    for row in 1:mc:rows,k in 1:kc:depth,column in 1:nc:columns
        last_row=min(rows,row+mc-1);last_k=min(depth,k+kc-1)
        last_column=min(columns,column+nc-1)
        _native_bezier_gemm_block!(@view(target[row:last_row,column:last_column]),
            @view(left[row:last_row,k:last_k]),@view(right[k:last_k,column:last_column]))
    end
    return target
end

function _native_bezier_sample_jacobian(mesh,cells,cell,rule,scaled)
    nodes=size(cells,1);samples=length(rule.space.points)
    coordinates=Matrix{Float64}(undef,nodes,3)
    for node in 1:nodes
        point=_mixed_node(mesh,cells,cell,node)
        for axis in 1:3
            coordinates[node,axis]=point[axis]
        end
    end
    derivatives=ntuple(_->zeros(Float64,samples,3),3)
    for axis in 1:rule.dimension
        _native_bezier_product!(derivatives[axis],rule.gradients[axis],coordinates)
    end
    primary=zeros(Float64,3,3)
    for node in axes(rule.primary_gradients,2),axis in 1:3,coordinate in 1:3
        primary[axis,coordinate]+=rule.primary_gradients[axis,node]*coordinates[node,coordinate]
    end
    normal=(0.0,0.0,0.0);scale=1.0
    if rule.dimension==2
        a=(primary[1,1],primary[1,2],primary[1,3]);b=(primary[2,1],primary[2,2],primary[2,3])
        normal=(a[2]*b[3]-a[3]*b[2],-a[1]*b[3]+a[3]*b[1],a[1]*b[2]-a[2]*b[1])
        magnitude=sqrt(normal[1]*normal[1]+normal[2]*normal[2]+normal[3]*normal[3])
        reciprocal=1.0/magnitude
        normal=(normal[1]*reciprocal,normal[2]*reciprocal,normal[3]*reciprocal)
        if scaled
            magnitude==0.0 && return zeros(Float64,samples)
            factor=1.0/magnitude
            normal=(normal[1]*factor,normal[2]*factor,normal[3]*factor)
        end
    elseif scaled
        a=(primary[1,1],primary[1,2],primary[1,3]);b=(primary[2,1],primary[2,2],primary[2,3])
        c=(primary[3,1],primary[3,2],primary[3,3])
        scale=1.0/abs(_native_bezier_det(a,b,c))
    end
    result=Vector{Float64}(undef,samples)
    for sample in 1:samples
        a=(derivatives[1][sample,1],derivatives[1][sample,2],derivatives[1][sample,3])
        b=(derivatives[2][sample,1],derivatives[2][sample,2],derivatives[2][sample,3])
        c=rule.dimension==2 ? normal :
            (derivatives[3][sample,1],derivatives[3][sample,2],derivatives[3][sample,3])
        result[sample]=_native_bezier_det(a,b,c)
    end
    if rule.dimension==3 && scaled
        for sample in eachindex(result)
            result[sample]*=scale
        end
    end
    return result
end

# Ainsworth--Sanchez conversion used by native tensor spaces. The recurrence
# updates divided differences and Bernstein rows in descending index order.
function _native_bezier_convert_line!(source,degree,start,stride,target)
    count=degree+1
    x=Float64[Float64(i)/degree for i in 0:degree]
    degree==0 && (x[1]=0.0)
    f=Float64[source[start+index*stride] for index in 0:degree]
    weights=zeros(Float64,count);weights[1]=1.0
    coefficients=zeros(Float64,count);coefficients[1]=f[1]
    for stage in 1:degree
        for index in degree:-1:stage
            f[index+1]=(f[index+1]-f[index])/(x[index+1]-x[index-stage+1])
        end
        for index in stage:-1:1
            factor=Float64(index)/stage
            weights[index+1]=factor*weights[index]*(1.0-x[stage])-
                (1.0-factor)*weights[index+1]*x[stage]
            coefficients[index+1]=factor*coefficients[index]+(1.0-factor)*coefficients[index+1]+
                f[stage+1]*weights[index+1]
        end
        weights[1]=-weights[1]*x[stage]
        coefficients[1]=coefficients[1]+f[stage+1]*weights[1]
    end
    for index in 0:degree
        target[start+index*stride]=coefficients[index+1]
    end
    return target
end

function _native_bezier_coefficients(space,samples)
    family=space.family;degree=space.degree;n=degree+1
    if family in (:tri,:tet)
        coefficients=zeros(Float64,length(samples),1)
        _native_bezier_product!(coefficients,space.inverse,reshape(samples,:,1))
        return vec(coefficients)
    end
    coefficients=similar(samples)
    if family===:qua
        for i in 0:n-1
            _native_bezier_convert_line!(samples,degree,i+1,n,coefficients)
        end
        for j in 0:n-1
            _native_bezier_convert_line!(coefficients,degree,j*n+1,1,coefficients)
        end
    elseif family===:hex || family===:pyr
        ni=family===:pyr ? space.nij+1 : n
        nz=family===:pyr ? space.nk+1 : n
        for ij in 0:ni*ni-1
            _native_bezier_convert_line!(samples,nz-1,ij+1,ni*ni,coefficients)
        end
        for i in 0:ni-1,k in 0:nz-1
            _native_bezier_convert_line!(coefficients,ni-1,i+k*ni*ni+1,ni,coefficients)
        end
        for jk in 0:ni*nz-1
            _native_bezier_convert_line!(coefficients,ni-1,jk*ni+1,1,coefficients)
        end
    else
        throw(ArgumentError("unsupported native Bernstein coefficient family $family"))
    end
    return coefficients
end

@inline _native_bezier_ij(i,j,n)=i+j*n-j*(j-1)÷2+1
@inline _native_bezier_ijk(i,j,k,n)=(n+2)*(n+1)*n÷6-
    (n-k+2)*(n-k+1)*(n-k)÷6+j*(n-k)-j*(j-1)÷2+i+1

function _native_bezier_triangle_split(coefficients,n)
    a=copy(coefficients)
    for stage in 1:n-1,j in 0:n-stage-1,i in n-1-j:-1:stage
        at=_native_bezier_ij(i,j,n);previous=_native_bezier_ij(i-1,j,n)
        a[at]=0.5*(a[previous]+a[at])
    end
    for stage in 1:n-1,j in n-1:-1:stage,i in 0:n-j-1
        at=_native_bezier_ij(i,j,n);previous=_native_bezier_ij(i,j-1,n)
        a[at]=0.5*(a[previous]+a[at])
    end
    b=copy(a)
    for stage in 1:n-1,j in 0:n-stage-1,i in 0:n-stage-j-1
        at=_native_bezier_ij(i,j,n)
        b[at]=b[_native_bezier_ij(i+1,j,n)]+b[_native_bezier_ij(i,j+1,n)]-b[at]
    end
    c=copy(b)
    for stage in 1:n-1,j in 0:n-stage-1,i in n-1-j:-1:stage
        at=_native_bezier_ij(i,j,n)
        c[at]=c[_native_bezier_ij(i-1,j,n)]+c[_native_bezier_ij(i-1,j+1,n)]-c[at]
    end
    d=copy(b)
    for stage in 1:n-1,j in n-1:-1:stage,i in 0:n-j-1
        at=_native_bezier_ij(i,j,n)
        d[at]=d[_native_bezier_ij(i,j-1,n)]+d[_native_bezier_ij(i+1,j-1,n)]-d[at]
    end
    return [a,b,c,d]
end

function _native_bezier_tet_transform!(coefficients,n,mode)
    for stage in 1:n-1
        if mode===:u || mode===:edge03at1
            for k in 0:n-stage-1,j in 0:n-stage-k-1,i in n-1-j-k:-1:stage
                at=_native_bezier_ijk(i,j,k,n);a=_native_bezier_ijk(i-1,j,k,n)
                coefficients[at]=mode===:u ? 0.5*(coefficients[a]+coefficients[at]) :
                    coefficients[a]+coefficients[_native_bezier_ijk(i-1,j,k+1,n)]-coefficients[at]
            end
        elseif mode===:v || mode===:edge03at2
            for k in 0:n-stage-1,j in n-1-k:-1:stage,i in 0:n-j-k-1
                at=_native_bezier_ijk(i,j,k,n);a=_native_bezier_ijk(i,j-1,k,n)
                coefficients[at]=mode===:v ? 0.5*(coefficients[a]+coefficients[at]) :
                    coefficients[a]+coefficients[_native_bezier_ijk(i,j-1,k+1,n)]-coefficients[at]
            end
        elseif mode===:w || mode===:edge12at3
            for k in n-1:-1:stage,j in 0:n-k-1,i in 0:n-j-k-1
                at=_native_bezier_ijk(i,j,k,n)
                coefficients[at]=mode===:w ? 0.5*(coefficients[_native_bezier_ijk(i,j,k-1,n)]+coefficients[at]) :
                    coefficients[_native_bezier_ijk(i+1,j,k-1,n)]+
                    coefficients[_native_bezier_ijk(i,j+1,k-1,n)]-coefficients[at]
            end
        elseif mode===:edge12at0
            for k in 0:n-stage-1,j in 0:n-stage-k-1,i in 0:n-stage-j-k-1
                at=_native_bezier_ijk(i,j,k,n)
                coefficients[at]=coefficients[_native_bezier_ijk(i+1,j,k,n)]+
                    coefficients[_native_bezier_ijk(i,j+1,k,n)]-coefficients[at]
            end
        else
            throw(ArgumentError("unknown native Tet subdivision mode"))
        end
    end
    return coefficients
end

function _native_bezier_tet_split(coefficients,n)
    a=copy(coefficients)
    _native_bezier_tet_transform!(a,n,:u)
    _native_bezier_tet_transform!(a,n,:v)
    _native_bezier_tet_transform!(a,n,:w)
    b=copy(a);_native_bezier_tet_transform!(b,n,:edge12at0)
    c=copy(b);_native_bezier_tet_transform!(c,n,:edge03at1)
    d=copy(b);_native_bezier_tet_transform!(d,n,:edge03at2)
    e=copy(d);_native_bezier_tet_transform!(e,n,:edge03at1)
    f=copy(c);_native_bezier_tet_transform!(f,n,:edge12at3)
    g=copy(d);_native_bezier_tet_transform!(g,n,:edge12at3)
    h=copy(e);_native_bezier_tet_transform!(h,n,:edge12at0)
    return [a,b,c,d,e,f,g,h]
end

function _native_bezier_casteljau!(workspace,n,start,stride)
    for stage in 1:n-1,index in stage:2:2n-stage-1
        at=start+index*stride
        workspace[at]=0.5*(workspace[at-stride]+workspace[at+stride])
    end
    return workspace
end

function _native_bezier_tensor_split(space,coefficients)
    dimension=space.family===:qua ? 2 : 3
    # Literal Pyramid subdivision uses degrees as widths and keeps copied tails.
    ni=space.family===:pyr ? space.nij : space.degree+1
    nz=space.family===:pyr ? space.nk : space.degree+1
    dimension==2 && (nz=1)
    wi=2ni-1;wz=dimension==2 ? 1 : 2nz-1
    workspace=zeros(Float64,wi*wi*wz)
    for i in 0:ni-1,j in 0:ni-1,k in 0:nz-1
        workspace[2i+2j*wi+2k*wi*wi+1]=coefficients[i+j*ni+k*ni*ni+1]
    end
    if dimension==3
        for i in 0:2:wi-1,j in 0:2:wi-1
            _native_bezier_casteljau!(workspace,nz,i+j*wi+1,wi*wi)
        end
    end
    for i in 0:2:wi-1,k in 0:wz-1
        _native_bezier_casteljau!(workspace,ni,i+k*wi*wi+1,wi)
    end
    for j in 0:wi-1,k in 0:wz-1
        _native_bezier_casteljau!(workspace,ni,j*wi+k*wi*wi+1,1)
    end
    children=Vector{Vector{Float64}}(undef,dimension==2 ? 4 : 8)
    for child in eachindex(children)
        output=copy(coefficients)
        firsti=((child-1)&1)*(ni-1)
        firstj=(((child-1)>>1)&1)*(ni-1)
        firstk=dimension==2 ? 0 : (((child-1)>>2)&1)*(nz-1)
        for i in 0:ni-1,j in 0:ni-1,k in 0:nz-1
            output[i+j*ni+k*ni*ni+1]=workspace[firsti+i+(firstj+j)*wi+(firstk+k)*wi*wi+1]
        end
        children[child]=output
    end
    return children
end

function _native_bezier_split(space,coefficients)
    space.family===:tri && return _native_bezier_triangle_split(coefficients,space.degree+1)
    space.family===:tet && return _native_bezier_tet_split(coefficients,space.degree+1)
    return _native_bezier_tensor_split(space,coefficients)
end

struct _NativeBezierDetDomain
    coefficients::Vector{Float64}
    minL::Float64
    maxL::Float64
    minB::Float64
    maxB::Float64
end

# std::min/std::max retain their first operand on ties or unordered comparisons.
# Julia's min/max propagate a later NaN and impose a signed-zero ordering.
@inline _native_bezier_min(a::Float64,b::Float64)=b<a ? b : a
@inline _native_bezier_max(a::Float64,b::Float64)=a<b ? b : a

function _native_bezier_min_element(values)
    result=values[firstindex(values)]
    for position in firstindex(values)+1:lastindex(values)
        values[position]<result && (result=values[position])
    end
    return result
end

function _native_bezier_domain(space,coefficients)
    minL=coefficients[space.corners[1]];maxL=minL
    for position in 2:length(space.corners)
        value=coefficients[space.corners[position]]
        minL=_native_bezier_min(minL,value);maxL=_native_bezier_max(maxL,value)
    end
    minB=coefficients[1];maxB=minB
    for position in 2:length(coefficients)
        value=coefficients[position]
        minB=_native_bezier_min(minB,value);maxB=_native_bezier_max(maxB,value)
    end
    return _NativeBezierDetDomain(coefficients,minL,maxL,minB,maxB)
end

@inline function _native_bezier_ready(domain,minL,maxL)
    tolerance=_native_bezier_max(abs(minL),abs(maxL))*1.e-3
    return (minL<=0.0 || domain.minB>0.0) && (maxL>=0.0 || domain.maxB<0.0) &&
        minL-domain.minB<tolerance && domain.maxB-maxL<tolerance
end

@inline _native_bezier_heap_less(a,b,maximum)=maximum ? a.maxB<b.maxB : a.minB>b.minB

function _native_bezier_heap_push!(domains,value,maximum,position=length(domains),top=1)
    while position>top
        parent=(position-2)÷2+1
        _native_bezier_heap_less(domains[parent],value,maximum) || break
        domains[position]=domains[parent];position=parent
    end
    domains[position]=value
    return domains
end

function _native_bezier_heap_adjust!(domains,position,count,value,maximum)
    top=position;child=position
    while child-1<(count-1)÷2
        child=2child+1
        _native_bezier_heap_less(domains[child],domains[child-1],maximum) && (child-=1)
        domains[position]=domains[child];position=child
    end
    if iseven(count) && child-1==(count-2)÷2
        child=2child+1
        domains[position]=domains[child-1];position=child-1
    end
    return _native_bezier_heap_push!(domains,value,maximum,position,top)
end

function _native_bezier_heap_make!(domains,maximum)
    count=length(domains)
    for parent in count÷2:-1:1
        _native_bezier_heap_adjust!(domains,parent,count,domains[parent],maximum)
    end
    return domains
end

function _native_bezier_extrema(space,coefficients)
    first=_native_bezier_domain(space,coefficients)
    domains=_NativeBezierDetDomain[first];minL=first.minL;maxL=first.maxL
    for maximum in (false,true)
        _native_bezier_heap_make!(domains,maximum)
        iteration=0
        while !_native_bezier_ready(domains[1],minL,maxL) && iteration+1<1000
            parent=domains[1];last=pop!(domains)
            isempty(domains) || _native_bezier_heap_adjust!(domains,1,length(domains),last,maximum)
            for child in _native_bezier_split(space,parent.coefficients)
                domain=_native_bezier_domain(space,child)
                minL=_native_bezier_min(minL,domain.minL);maxL=_native_bezier_max(maxL,domain.maxL)
                push!(domains,domain)
                _native_bezier_heap_push!(domains,domain,maximum)
            end
            iteration+=1
        end
    end
    minimum=floatmax(Float64);maximum=-minimum
    for domain in domains
        minimum=_native_bezier_min(minimum,domain.minB)
        maximum=_native_bezier_max(maximum,domain.maxB)
    end
    return minimum,maximum
end

function _native_cell_quality_bezier(mesh,cells,cell,msh,name,caller)
    rule=_native_bezier_jacobian_rule(msh,caller)
    samples=_native_bezier_sample_jacobian(mesh,cells,cell,rule,name=="minSJ")
    coefficients=_native_bezier_coefficients(rule.space,samples)
    name=="minSJ" && return _native_bezier_min_element(coefficients)
    lower,upper=_native_bezier_extrema(rule.space,coefficients)
    return name=="minDetJac" ? lower : upper
end
