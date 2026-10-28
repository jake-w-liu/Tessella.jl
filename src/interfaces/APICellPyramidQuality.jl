# Native Bergot interpolation for public N-cell qualities. The standalone
# Pyramid FEM interpolant retains its own stable dot/inverse/fallback policy.
struct _NativeBezierPyramidReference
    order::Int
    modes::Vector{NTuple{3,Int}}
    coefficients::Matrix{Float64}
end

const _NATIVE_BEZIER_PYRAMIDS=Dict{Int,_NativeBezierPyramidReference}()

function _native_bezier_legendre(order,value)
    values=Vector{Float64}(undef,order+1);values[1]=1.0
    for i in 0:order-1
        result=(2.0*i+1.0)*value*values[i+1]
        i>0 && (result-=Float64(i)*values[i])
        values[i+2]=result/Float64(i+1)
    end
    return values
end

function _native_bezier_legendre_derivatives(order,value,values)
    derivatives=Vector{Float64}(undef,order+1)
    if value==1.0 || value==-1.0
        for k in 0:order
            derivatives[k+1]=0.5*k*(k+1)
        end
        if value==-1.0 && order>=2
            for k in 2:2:order
                derivatives[k+1]=-derivatives[k+1]
            end
        end
        return derivatives
    end
    derivatives[1]=0.0;denominator=1.0-value*value
    for i in 1:order
        first=-value*i;second=Float64(i)
        derivatives[i+1]=(first*values[i+1]+second*values[i])/denominator
    end
    return derivatives
end

function _native_bezier_jacobi(order,alpha,beta,value)
    total=alpha+beta;difference=alpha*alpha-beta*beta
    values=Vector{Float64}(undef,order+1);values[1]=1.0
    if order>=1
        values[2]=0.5*(2.0*(alpha+1.0)+(total+2.0)*(value-1.0))
    end
    for i in 1:order-1
        ii=Float64(i);twice=2.0*ii
        a1=2.0*(ii+1.0)*(ii+total+1.0)*(twice+total)
        a2=(twice+total+1.0)*difference
        a3=1.0
        for k in 0:2
            a3*=twice+total+Float64(k)
        end
        a4=2.0*(ii+alpha)*(ii+beta)*(twice+total+2.0)
        values[i+2]=((a2+a3*value)*values[i+1]-a4*values[i])/a1
    end
    return values
end

function _native_bezier_jacobi_derivatives(order,alpha,beta,value,values)
    total=alpha+beta;derivatives=Vector{Float64}(undef,order+1)
    if value==1.0 || value==-1.0
        coefficient=Int(value==1.0 ? alpha : beta)
        last=max(order,1)+coefficient
        factorials=Vector{Float64}(undef,last+1);factorials[1]=1.0
        for i in 1:last
            factorials[i+1]=i*factorials[i]
        end
        derivatives[1]=0.0
        for k in 1:order
            derivatives[k+1]=0.5*(k+total+1.0)*factorials[k+coefficient+1]/
                (factorials[coefficient+2]*factorials[k])
        end
        if value==-1.0 && order>=2
            for k in 2:2:order
                derivatives[k+1]=-derivatives[k+1]
            end
        end
        return derivatives
    end
    derivatives[1]=0.0
    order>=1 && (derivatives[2]=0.5*(total+2.0))
    for i in 2:order
        ii=Float64(i);a=2.0*ii+total
        denominator=a*(1.0-value*value)
        first=ii*(alpha-beta-a*value);second=2.0*(ii+alpha)*(ii+beta)
        derivatives[i+1]=(first*values[i+1]+second*values[i])/denominator
    end
    return derivatives
end

function _native_bezier_bergot_modes(order,modes,point,gradient)
    u,v,w=point;width=1.0-w
    reduced_u=w==1.0 ? 0.0 : u/width;reduced_v=w==1.0 ? 0.0 : v/width
    uf=_native_bezier_legendre(order,reduced_u);vf=_native_bezier_legendre(order,reduced_v)
    ug=gradient ? _native_bezier_legendre_derivatives(order,reduced_u,uf) : Float64[]
    vg=gradient ? _native_bezier_legendre_derivatives(order,reduced_v,vf) : Float64[]
    wf=Vector{Vector{Float64}}(undef,order+1);wg=similar(wf)
    transformed=2.0*w-1.0
    for maximum in 0:order
        wf[maximum+1]=_native_bezier_jacobi(order-maximum,2.0*maximum+2.0,0.0,transformed)
        if gradient
            wg[maximum+1]=_native_bezier_jacobi_derivatives(order-maximum,2.0*maximum+2.0,
                0.0,transformed,wf[maximum+1])
        end
    end
    values=Vector{Float64}(undef,length(modes));derivatives=zeros(Float64,3,length(modes))
    for (index,(i,j,k)) in enumerate(modes)
        maximum=max(i,j);f=wf[maximum+1][k+1]
        values[index]=uf[i+1]*vf[j+1]*f*_native_bezier_integer_power(width,maximum)
        gradient || continue
        g=wg[maximum+1][k+1]
        if maximum==0
            derivatives[3,index]=2.0*g
        elseif maximum==1
            if i==0
                derivatives[2,index]=f;derivatives[3,index]=2.0*v*g
            elseif j==0
                derivatives[1,index]=f;derivatives[3,index]=2.0*u*g
            else
                derivatives[1,index]=reduced_v*f;derivatives[2,index]=reduced_u*f
                derivatives[3,index]=reduced_u*reduced_v*f+2.0*reduced_u*v*g
            end
        else
            power2=_native_bezier_integer_power(width,maximum-2);power1=power2*width
            derivatives[1,index]=ug[i+1]*vf[j+1]*f*power1
            derivatives[2,index]=uf[i+1]*vg[j+1]*f*power1
            derivatives[3,index]=f*power2*(u*ug[i+1]*vf[j+1]+v*uf[i+1]*vg[j+1])+
                uf[i+1]*vf[j+1]*power1*(2.0*width*g-maximum*f)
        end
    end
    return values,derivatives
end

# Bergot (i,j) mode validity: the complete pyramid keeps every pair while the
# incomplete (serendipity) space drops pairs with i+j>order except the two
# "shell" modes with i+j==order+1 and min(i,j)==1 — upstream BergotBasis::
# validIJ. Order-two serendipity therefore keeps 13 modes for Pyr13.
_native_bergot_valid_ij(i,j,order,serendipity)=
    serendipity ? (i+j<=order || (i+j==order+1 && (i==1 || j==1))) :
        (i<=order && j<=order)

function _native_bezier_pyramid_modes(order,serendipity)
    modes=NTuple{3,Int}[]
    for i in 0:order,j in 0:order
        _native_bergot_valid_ij(i,j,order,serendipity) || continue
        for k in 0:order-max(i,j)
            push!(modes,(i,j,k))
        end
    end
    return modes
end

function _native_bezier_pyramid_reference_build(msh,caller)
    spec=msh_spec(msh);_native_cell_numeric_contract(spec,caller)
    modes=_native_bezier_pyramid_modes(spec.order,spec.serendipity)
    length(modes)==spec.nnodes || throw(ArgumentError("$caller: native Pyramid mode count differs"))
    points=_native_cell_reference_nodes(spec)
    matrix=Matrix{Float64}(undef,spec.nnodes,spec.nnodes)
    for column in 1:spec.nnodes
        values,_=_native_bezier_bergot_modes(spec.order,modes,
            (points[1,column],points[2,column],points[3,column]),false)
        for row in 1:spec.nnodes
            matrix[row,column]=values[row]
        end
    end
    return _NativeBezierPyramidReference(spec.order,modes,_native_cell_inverse(matrix))
end

function _native_bezier_pyramid_reference(msh,caller)
    lock(_NATIVE_BEZIER_LOCK)
    try
        haskey(_NATIVE_BEZIER_PYRAMIDS,msh) && return _NATIVE_BEZIER_PYRAMIDS[msh]
        data=_native_bezier_pyramid_reference_build(msh,caller)
        _NATIVE_BEZIER_PYRAMIDS[msh]=data
        return data
    finally
        unlock(_NATIVE_BEZIER_LOCK)
    end
end

function _native_bezier_pyramid_shape(data,point)
    modes,mode_gradients=_native_bezier_bergot_modes(data.order,data.modes,point,true)
    count=length(modes);values=zeros(Float64,count);gradients=zeros(Float64,3,count)
    for row in 1:count,column in 1:count
        coefficient=data.coefficients[row,column]
        values[row]+=coefficient*modes[column]
        for axis in 1:3
            gradients[axis,row]+=coefficient*mode_gradients[axis,column]
        end
    end
    return values,[(gradients[1,row],gradients[2,row],gradients[3,row]) for row in 1:count]
end

# Reduced-coordinate mode derivatives for the collapsed pyramid map. With
# u=(1-w)r, v=(1-w)s every nonzero term of G_r carries a (1-w) factor, so the
# Jacobian columns A=G_r/(1-w), B=G_s/(1-w), C=G_w are finite polynomials in
# (r,s,w). Adding r*A+s*B to the third column recovers the physical
# w-derivative, so det(A,B,C) is exactly the physical Jacobian determinant —
# a polynomial of tensor degree at most (3p-1,3p-1,3p-3) for order p.
# Bernstein bounds on the [-1,1]^2 x [0,1] box therefore cover the entire
# pyramid domain, including every limiting apex direction; unlike the
# physical (u,v,w) sampling used by the quality queries, no sub-range or
# non-conservative convention is involved.
function _native_bezier_bergot_reduced(order,modes,point)
    r,s,w=point;width=1.0-w
    uf=_native_bezier_legendre(order,r);vf=_native_bezier_legendre(order,s)
    ug=_native_bezier_legendre_derivatives(order,r,uf)
    vg=_native_bezier_legendre_derivatives(order,s,vf)
    transformed=2.0*w-1.0
    wf=Vector{Vector{Float64}}(undef,order+1);wg=similar(wf)
    for maximum in 0:order
        wf[maximum+1]=_native_bezier_jacobi(order-maximum,2.0*maximum+2.0,0.0,
            transformed)
        wg[maximum+1]=_native_bezier_jacobi_derivatives(order-maximum,
            2.0*maximum+2.0,0.0,transformed,wf[maximum+1])
    end
    count=length(modes)
    amodes=zeros(Float64,count);bmodes=zeros(Float64,count)
    cmodes=zeros(Float64,count)
    for (index,(i,j,k)) in enumerate(modes)
        maximum=max(i,j);f=wf[maximum+1][k+1];g=wg[maximum+1][k+1]
        if maximum==0
            # M = J_k(2w-1): r and s derivatives vanish exactly.
            cmodes[index]=uf[i+1]*vf[j+1]*2.0*g
            continue
        end
        power=_native_bezier_integer_power(width,maximum-1)
        amodes[index]=ug[i+1]*vf[j+1]*f*power
        bmodes[index]=uf[i+1]*vg[j+1]*f*power
        cmodes[index]=uf[i+1]*vf[j+1]*power*(2.0*width*g-maximum*f)
    end
    return amodes,bmodes,cmodes
end

function _native_bezier_pyramid_reduced_shape(data,point)
    amodes,bmodes,cmodes=_native_bezier_bergot_reduced(data.order,data.modes,point)
    count=length(amodes)
    a=zeros(Float64,count);b=zeros(Float64,count);c=zeros(Float64,count)
    for row in 1:count,column in 1:count
        coefficient=data.coefficients[row,column]
        a[row]+=coefficient*amodes[column]
        b[row]+=coefficient*bmodes[column]
        c[row]+=coefficient*cmodes[column]
    end
    return a,b,c
end
