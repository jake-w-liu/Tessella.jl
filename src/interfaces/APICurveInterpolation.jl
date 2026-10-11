
# Operation grouping follows the dynamic, column-major Eigen inverse path used
# by Gmsh's fullMatrix<double>. The admitted matrices fit its unblocked LU.
function _dim01_native_line_inverse(matrix::Matrix{Float64})
    panel_width=4
    n=size(matrix,1)
    size(matrix,2)==n && 1<=n<=11 || throw(ArgumentError("invalid small inverse matrix"))
    lu=copy(matrix)
    result=zeros(Float64,n,n)
    for i in 1:n
        result[i,i]=1.0
    end
    for k in 1:n
        pivot=k
        biggest=abs(lu[k,k])
        for row in k+1:n
            score=abs(lu[row,k])
            if score>biggest
                pivot=row
                biggest=score
            end
        end
        biggest!=0.0 || throw(ArgumentError("singular small inverse matrix"))
        if pivot!=k
            for column in 1:n
                lu[k,column],lu[pivot,column]=lu[pivot,column],lu[k,column]
                result[k,column],result[pivot,column]=result[pivot,column],result[k,column]
            end
        end
        diagonal=lu[k,k]
        for row in k+1:n
            lu[row,k]/=diagonal
        end
        for column in k+1:n,row in k+1:n
            lu[row,column]-=lu[row,k]*lu[k,column]
        end
    end
    for start in 1:panel_width:n
        stop=min(n,start+panel_width-1)
        for i in start:stop,column in 1:n
            result[i,column]*=1.0
            value=result[i,column]
            for row in i+1:stop
                result[row,column]-=value*lu[row,i]
            end
        end
        for column in 1:n,row in stop+1:n
            value=0.0
            for k in start:stop
                value+=lu[row,k]*result[k,column]
            end
            result[row,column]+=(-1.0)*value
        end
    end
    for stop in n:-panel_width:1
        start=max(1,stop-panel_width+1)
        for i in stop:-1:start
            reciprocal=1.0/lu[i,i]
            for column in 1:n
                result[i,column]*=reciprocal
                value=result[i,column]
                for row in start:i-1
                    result[row,column]-=value*lu[row,i]
                end
            end
        end
        for column in 1:n,row in 1:start-1
            value=0.0
            for k in start:stop
                value+=lu[row,k]*result[k,column]
            end
            result[row,column]+=(-1.0)*value
        end
    end
    return result
end


# Operation-local numeric layer only. The caller supplies the independently
# verified native floating inverse implementation. No primary constants,
# global coefficient cache or ideal polynomial _dim01_native_line_coefficients are substituted.

function _dim01_native_line_order(order)
    order isa Integer && !(order isa Bool) && 1<=order<=10 ||
        throw(ArgumentError("native Line interpolation order must be 1 through 10"))
    return Int(order)
end

# Numeric.h::_dim01_native_line_power literal operation grouping for all admitted exponents.
@inline function _dim01_native_line_power(a::Float64,n::Int)
    n==0 && return 1.0
    n==1 && return a
    n==2 && return a*a
    n==3 && return a*a*a
    if n==4
        a2=a*a;return a2*a2
    elseif n==5
        a2=a*a;return a2*a2*a
    elseif n==6
        a3=a*a*a;return a3*a3
    elseif n==7
        a3=a*a*a;return a3*a3*a
    elseif n==8
        a2=a*a;a4=a2*a2;return a4*a4
    elseif n==9
        a3=a*a*a;return a3*a3*a3
    elseif n==10
        a2=a*a;a4=a2*a2;return a4*a4*a2
    end
    throw(ArgumentError("native Line interpolation _dim01_native_line_exponent is outside zero through ten"))
end

function _dim01_native_line_reference_nodes(order)
    p=_dim01_native_line_order(order);nodes=Vector{Float64}(undef,p+1)
    scale=2.0/p
    nodes[1]=0.0*scale-1.0
    nodes[2]=Float64(p)*scale-1.0
    for slot in 3:p+1
        nodes[slot]=Float64(slot-2)*scale-1.0
    end
    return nodes
end

@inline _dim01_native_line_exponent(p::Int,slot::Int)=slot==1 ? 0 : slot==2 ? p : slot-2

function _dim01_native_line_vandermonde(order)
    p=_dim01_native_line_order(order);nodes=_dim01_native_line_reference_nodes(p)
    result=Matrix{Float64}(undef,p+1,p+1)
    # polynomialBasis.cpp builds V(monomial slot,node slot), then stores V^-1
    # as _dim01_native_line_coefficients(basis slot,monomial slot).
    for row in 1:p+1,column in 1:p+1
        dd=1.0
        dd*=_dim01_native_line_power(nodes[column],_dim01_native_line_exponent(p,row))
        result[row,column]=dd
    end
    return result
end

function _dim01_native_line_coefficients(order)
    p=_dim01_native_line_order(order)
    result=_dim01_native_line_inverse(_dim01_native_line_vandermonde(p))
    result isa Matrix{Float64} && size(result)==(p+1,p+1) && all(isfinite,result) ||
        throw(ArgumentError("native Line floating inverse returned an invalid coefficient matrix"))
    return result
end

function _dim01_native_line_weights(order,u::Float64,coefficient::Matrix{Float64})
    p=_dim01_native_line_order(order)
    isfinite(u) && size(coefficient)==(p+1,p+1) ||
        throw(ArgumentError("native Line shape interpolation input is invalid"))
    terms=Vector{Float64}(undef,p+1)
    for slot in 1:p+1
        value=1.0
        value*=_dim01_native_line_power(u,_dim01_native_line_exponent(p,slot))
        terms[slot]=value
    end
    result=Vector{Float64}(undef,p+1)
    for row in 1:p+1
        value=0.0
        for column in 1:p+1
            value+=coefficient[row,column]*terms[column]
        end
        result[row]=value
    end
    return result
end


# polynomialBasis.cpp::df: multiply coefficient, power and exponent in that
# order, then add each monomial contribution sequentially.
function _dim01_native_line_gradients(order,u::Float64,coefficient::Matrix{Float64})
    p=_dim01_native_line_order(order)
    isfinite(u) && size(coefficient)==(p+1,p+1) ||
        throw(ArgumentError("native Line gradient interpolation input is invalid"))
    powers=Vector{Float64}(undef,p)
    for column in 2:p+1
        exponent=_dim01_native_line_exponent(p,column)
        powers[column-1]=_dim01_native_line_power(u,exponent-1)
    end
    result=Vector{Float64}(undef,p+1)
    for row in 1:p+1
        value=0.0
        for column in 2:p+1
            exponent=_dim01_native_line_exponent(p,column)
            value+=(coefficient[row,column]*powers[column-1])*exponent
        end
        result[row]=value
    end
    return result
end
