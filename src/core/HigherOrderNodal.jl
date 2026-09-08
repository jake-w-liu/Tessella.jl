# Higher-order nodal bases. This file is included inside MeshFunctionSpaces.

using LinearAlgebra: I, inv, opnorm

const _NODAL_COEFFICIENT_PRECISION=256
const _NODAL_ROUNDING_FACTOR=256eps(Float64)

struct _PolynomialInterpolant
    exponents::Vector{NTuple{3,Int}}
    coefficients::Matrix{Float64}
    high_precision_coefficients::Matrix{BigFloat}
    nodes::Matrix{Float64}
end

struct _PyramidInterpolant
    modes::Vector{NTuple{3,Int}}
    coefficients::Matrix{Float64}
    extended_coefficients::Matrix{BigFloat}
    nodes::Matrix{Float64}
end

const _POLYNOMIAL_INTERPOLANTS=Dict{Int,_PolynomialInterpolant}()
const _PYRAMID_INTERPOLANTS=Dict{Int,_PyramidInterpolant}()
const _NODAL_INTERPOLANT_LOCK=ReentrantLock()

function _reference_root_index(p::Int,x::Float64)
    p==0 && return -1
    for index in 0:p
        x==index/p && return index
    end
    return -1
end

function _exact_reference_coordinate(p::Int,x::Float64)
    root=_reference_root_index(p,x)
    return root<0 ? _Exact(x) : _Exact(root)/_Exact(p)
end

function _cardinal_pair_exact(p::Int,index::Int,x::Float64,
                              caller::AbstractString,point::Int,node::Int)
    value=_Exact(1)
    derivative=_Exact(0)
    exact_x=_exact_reference_coordinate(p,x)
    for other in 0:p
        other==index && continue
        factor=(_Exact(p)*exact_x-other)/(index-other)
        factor_derivative=_Exact(p)/_Exact(index-other)
        derivative=derivative*factor+value*factor_derivative
        value*=factor
    end
    return _exact_basis_result(value,caller,point,node,0),
           _exact_basis_result(derivative,caller,point,node,1)
end

function _cardinal_pair(p::Int,index::Int,x::Float64,
                        caller::AbstractString,point::Int,node::Int)
    value=1.0
    derivative=0.0
    derivative_scale=0.0
    scaled=p*x
    root=_reference_root_index(p,x)
    unstable_factor=false
    for other in 0:p
        other==index && continue
        numerator=scaled-other
        factor_scale=max(1.0,abs(scaled),abs(other))
        unstable_factor|=isfinite(numerator) && root!=other &&
                         abs(numerator)<=_NODAL_ROUNDING_FACTOR*factor_scale
        factor=numerator/(index-other)
        factor_derivative=p/(index-other)
        derivative_scale=derivative_scale*abs(factor)+
                         abs(value*factor_derivative)
        derivative=derivative*factor+value*factor_derivative
        value*=factor
    end
    value_ok=isfinite(value) && value!=0.0
    if value==0.0
        value_ok=root>=0 && root!=index
    end
    derivative_ok=isfinite(derivative) &&
        (derivative_scale==0.0 ||
         abs(derivative)>_NODAL_ROUNDING_FACTOR*derivative_scale)
    (value_ok && derivative_ok && !unstable_factor) &&
        return value,derivative
    return _cardinal_pair_exact(p,index,x,caller,point,node)
end

function _falling_pair_exact(p::Int,degree::Int,x::Float64,
                             caller::AbstractString,point::Int,node::Int)
    value=_Exact(1)
    derivative=_Exact(0)
    exact_x=_exact_reference_coordinate(p,x)
    for offset in 0:degree-1
        factor=(_Exact(p)*exact_x-offset)/(offset+1)
        factor_derivative=_Exact(p)/_Exact(offset+1)
        derivative=derivative*factor+value*factor_derivative
        value*=factor
    end
    return _exact_basis_result(value,caller,point,node,0),
           _exact_basis_result(derivative,caller,point,node,1)
end

function _falling_pair(p::Int,degree::Int,x::Float64,
                       caller::AbstractString,point::Int,node::Int)
    degree==0 && return 1.0,0.0
    value=1.0
    derivative=0.0
    derivative_scale=0.0
    scaled=p*x
    root=_reference_root_index(p,x)
    unstable_factor=false
    for offset in 0:degree-1
        numerator=scaled-offset
        factor_scale=max(1.0,abs(scaled),abs(offset))
        unstable_factor|=isfinite(numerator) && root!=offset &&
                         abs(numerator)<=_NODAL_ROUNDING_FACTOR*factor_scale
        factor=numerator/(offset+1)
        factor_derivative=p/(offset+1)
        derivative_scale=derivative_scale*abs(factor)+
                         abs(value*factor_derivative)
        derivative=derivative*factor+value*factor_derivative
        value*=factor
    end
    value_ok=isfinite(value) && value!=0.0
    if value==0.0
        value_ok=0<=root<degree
    end
    derivative_ok=isfinite(derivative) &&
        (derivative_scale==0.0 ||
         abs(derivative)>_NODAL_ROUNDING_FACTOR*derivative_scale)
    (value_ok && derivative_ok && !unstable_factor) &&
        return value,derivative
    return _falling_pair_exact(p,degree,x,caller,point,node)
end

function _cardinal_table(p::Int,x::Float64,derivative_scale::Float64,
                         caller::AbstractString,point::Int)
    values=Vector{Float64}(undef,p+1)
    derivatives=Vector{Float64}(undef,p+1)
    for index in 0:p
        value,derivative=
            _cardinal_pair(p,index,x,caller,point,index+1)
        values[index+1]=value
        derivatives[index+1]=derivative_scale*derivative
    end
    return values,derivatives
end

function _falling_table(p::Int,x::Float64,caller::AbstractString,point::Int)
    values=Vector{Float64}(undef,p+1)
    derivatives=Vector{Float64}(undef,p+1)
    for degree in 0:p
        values[degree+1],derivatives[degree+1]=
            _falling_pair(p,degree,x,caller,point,degree+1)
    end
    return values,derivatives
end

function _internal_cardinal_pair(p::Int,index::Int,x::Float64,
                                 caller::AbstractString,point::Int,node::Int)
    value=1.0
    derivative=0.0
    derivative_scale=0.0
    scaled=p*x
    root=_reference_root_index(p,x)
    unstable_factor=false
    for other in 1:p-1
        other==index && continue
        numerator=scaled-other
        factor_scale=max(1.0,abs(scaled),abs(other))
        unstable_factor|=isfinite(numerator) && root!=other &&
                         abs(numerator)<=_NODAL_ROUNDING_FACTOR*factor_scale
        factor=numerator/(index-other)
        factor_derivative=p/(index-other)
        derivative_scale=derivative_scale*abs(factor)+
                         abs(value*factor_derivative)
        derivative=derivative*factor+value*factor_derivative
        value*=factor
    end
    value_ok=isfinite(value) && value!=0.0
    if value==0.0
        value_ok=1<=root<p && root!=index
    end
    derivative_ok=isfinite(derivative) &&
        (derivative_scale==0.0 ||
         abs(derivative)>_NODAL_ROUNDING_FACTOR*derivative_scale)
    (value_ok && derivative_ok && !unstable_factor) &&
        return value,derivative

    exact_value=_Exact(1)
    exact_derivative=_Exact(0)
    exact_x=_exact_reference_coordinate(p,x)
    for other in 1:p-1
        other==index && continue
        factor=(_Exact(p)*exact_x-other)/(index-other)
        factor_derivative=_Exact(p)/_Exact(index-other)
        exact_derivative=exact_derivative*factor+
                         exact_value*factor_derivative
        exact_value*=factor
    end
    return _exact_basis_result(exact_value,caller,point,node,0),
           _exact_basis_result(exact_derivative,caller,point,node,1)
end

function _stable_nodal_sum(terms,caller::AbstractString,point::Int,node::Int,
                           component::Int)
    value=0.0
    scale=0.0
    for term in terms
        value+=term
        scale+=abs(term)
    end
    if isfinite(value) && isfinite(scale) &&
       (scale==0.0 || abs(value)>_NODAL_ROUNDING_FACTOR*scale)
        return value
    end
    exact=sum(_Exact(term) for term in terms;init=_Exact(0))
    return _exact_basis_result(exact,caller,point,node,component)
end

@inline function _simplex_alphas(q::NTuple{D,Int},p::Int) where {D}
    return (p-sum(q),q...)
end

function _simplex_value(values,alphas,caller,point,node)
    factors=ntuple(index->values[index][alphas[index]+1],length(alphas))
    return _basis_product(1.0,factors,caller,point,node,0)
end

function _simplex_gradient_component(values,derivatives,alphas,gradients,
                                     axis::Int,caller,point,node)
    terms=ntuple(length(alphas)) do active
        factors=ntuple(length(alphas)) do index
            index==active ? derivatives[index][alphas[index]+1] :
                            values[index][alphas[index]+1]
        end
        _basis_product(
            gradients[active][axis],factors,caller,point,node,axis)
    end
    return _stable_nodal_sum(terms,caller,point,node,axis)
end

function _write_complete_line!(result,coordinates,point_count,p,gradient,
                               caller)
    layout=_line_monomials(p)
    cursor=0
    for point in 1:point_count
        base=3point-2
        values,derivatives=_cardinal_table(
            p,(coordinates[base]+1)/2,0.5,caller,point)
        for (node,q) in enumerate(layout)
            if gradient
                cursor+=1;result[cursor]=derivatives[q[1]+1]
                cursor+=1;result[cursor]=0.0
                cursor+=1;result[cursor]=0.0
            else
                cursor+=1;result[cursor]=values[q[1]+1]
            end
        end
    end
    return result
end

function _write_complete_simplex!(result,coordinates,point_count,p,gradient,
                                  layout,dimension::Int,caller)
    gradients=dimension==2 ?
        ((-1.0,-1.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0)) :
        ((-1.0,-1.0,-1.0),(1.0,0.0,0.0),
         (0.0,1.0,0.0),(0.0,0.0,1.0))
    cursor=0
    for point in 1:point_count
        base=3point-2
        lambdas=dimension==2 ?
            (_origin_barycentric(
                 (coordinates[base],coordinates[base+1]),caller,point),
             coordinates[base],coordinates[base+1]) :
            (_origin_barycentric(
                 (coordinates[base],coordinates[base+1],coordinates[base+2]),
                 caller,point),coordinates[base],coordinates[base+1],
             coordinates[base+2])
        tables=map(lambda->_falling_table(p,lambda,caller,point),lambdas)
        values=map(first,tables)
        derivatives=map(last,tables)
        for (node,q) in enumerate(layout)
            alphas=_simplex_alphas(q,p)
            if gradient
                for axis in 1:3
                    cursor+=1
                    result[cursor]=_simplex_gradient_component(
                        values,derivatives,alphas,gradients,axis,
                        caller,point,node)
                end
            else
                cursor+=1
                result[cursor]=_simplex_value(
                    values,alphas,caller,point,node)
            end
        end
    end
    return result
end

function _write_complete_tensor!(result,coordinates,point_count,p,gradient,
                                 layout,::Val{D},caller) where {D}
    cursor=0
    for point in 1:point_count
        base=3point-2
        tables=ntuple(Val(D)) do axis
            _cardinal_table(
                p,(coordinates[base+axis-1]+1)/2,0.5,caller,point)
        end
        for (node,q) in enumerate(layout)
            factors=ntuple(
                axis->tables[axis][1][q[axis]+1],Val(D))
            if gradient
                for axis in 1:3
                    cursor+=1
                    if axis>D
                        result[cursor]=0.0
                    else
                        derivative_factors=ntuple(Val(D)) do other
                            other==axis ? tables[other][2][q[other]+1] :
                                          factors[other]
                        end
                        result[cursor]=_basis_product(
                            1.0,derivative_factors,caller,point,node,axis)
                    end
                end
            else
                cursor+=1
                result[cursor]=_basis_product(
                    1.0,factors,caller,point,node,0)
            end
        end
    end
    return result
end

function _write_complete_prism!(result,coordinates,point_count,p,gradient,
                                caller)
    layout=_pri_monomials(p)
    triangle_gradients=((-1.0,-1.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0))
    cursor=0
    for point in 1:point_count
        base=3point-2
        lambdas=(_origin_barycentric(
                     (coordinates[base],coordinates[base+1]),caller,point),
                 coordinates[base],coordinates[base+1])
        triangle_tables=map(
            lambda->_falling_table(p,lambda,caller,point),lambdas)
        triangle_values=map(first,triangle_tables)
        triangle_derivatives=map(last,triangle_tables)
        line_values,line_derivatives=_cardinal_table(
            p,(coordinates[base+2]+1)/2,0.5,caller,point)
        for (node,q) in enumerate(layout)
            alphas=_simplex_alphas((q[1],q[2]),p)
            triangle_value=_simplex_value(
                triangle_values,alphas,caller,point,node)
            line_value=line_values[q[3]+1]
            if gradient
                for axis in 1:2
                    cursor+=1
                    triangle_gradient=_simplex_gradient_component(
                        triangle_values,triangle_derivatives,alphas,
                        triangle_gradients,axis,caller,point,node)
                    result[cursor]=_basis_product(
                        1.0,(triangle_gradient,line_value),
                        caller,point,node,axis)
                end
                cursor+=1
                result[cursor]=_basis_product(
                    1.0,(triangle_value,line_derivatives[q[3]+1]),
                    caller,point,node,3)
            else
                cursor+=1
                result[cursor]=_basis_product(
                    1.0,(triangle_value,line_value),caller,point,node,0)
            end
        end
    end
    return result
end

# Gmsh 4.15.2 cannot reliably construct its incomplete order-3--9 Pyramid
# bases. Those catalog types contain only vertices and edge-interior nodes, so
# the native extension corrects the rational linear vertex basis with exact
# order-p edge traces.
const _PYRAMID_EDGE_PARAMETERS=(
    (0.5,0.0,0.0,0.5),
    (0.0,0.5,0.0,0.5),
    (0.0,0.0,1.0,0.0),
    (0.0,0.5,0.0,0.5),
    (0.0,0.0,1.0,0.0),
    (-0.5,0.0,0.0,0.5),
    (0.0,0.0,1.0,0.0),
    (0.0,0.0,1.0,0.0),
)

function _write_serendipity_pyramid!(result,coordinates,point_count,p,
                                     gradient,caller)
    node_count=5+8(p-1)
    cursor=0
    values=Vector{Float64}(undef,node_count)
    gradients=Matrix{Float64}(undef,3,node_count)
    for point in 1:point_count
        base=3point-2
        u,v,w=coordinates[base],coordinates[base+1],coordinates[base+2]
        primary_values=_first_order_values(
            Val(:pyr),u,v,w,caller,point)
        primary_gradients=_first_order_gradients(
            Val(:pyr),u,v,w,caller,point)
        for node in 1:5
            values[node]=primary_values[node]
            for axis in 1:3
                gradients[axis,node]=primary_gradients[node][axis]
            end
        end
        edge_node=5
        for (edge,(first_zero,second_zero)) in enumerate(_PYR_EDGES)
            first=first_zero+1
            second=second_zero+1
            du,dv,dw,offset=_PYRAMID_EDGE_PARAMETERS[edge]
            parameter=offset+du*u+dv*v+dw*w
            parameter_gradient=(du,dv,dw)
            for internal in 1:p-1
                edge_node+=1
                polynomial,polynomial_derivative=_internal_cardinal_pair(
                    p,internal,parameter,caller,point,edge_node)
                target=internal/p
                normalization=target*(1-target)
                edge_value=_basis_quotient(
                    1.0,(primary_values[first],primary_values[second],
                         polynomial),(normalization,),
                    caller,point,edge_node,0)
                values[edge_node]=edge_value
                for axis in 1:3
                    first_term=_basis_quotient(
                        1.0,(primary_gradients[first][axis],
                             primary_values[second],polynomial),
                        (normalization,),caller,point,edge_node,axis)
                    second_term=_basis_quotient(
                        1.0,(primary_values[first],
                             primary_gradients[second][axis],polynomial),
                        (normalization,),caller,point,edge_node,axis)
                    parameter_term=_basis_quotient(
                        1.0,(primary_values[first],primary_values[second],
                             polynomial_derivative,
                             parameter_gradient[axis]),
                        (normalization,),caller,point,edge_node,axis)
                    edge_gradient=_stable_nodal_sum(
                        (first_term,second_term,parameter_term),
                        caller,point,edge_node,axis)
                    gradients[axis,edge_node]=edge_gradient
                    gradients[axis,first]=_stable_nodal_sum(
                        (gradients[axis,first],-(1-target)*edge_gradient),
                        caller,point,first,axis)
                    gradients[axis,second]=_stable_nodal_sum(
                        (gradients[axis,second],-target*edge_gradient),
                        caller,point,second,axis)
                end
                values[first]=_stable_nodal_sum(
                    (values[first],-(1-target)*edge_value),
                    caller,point,first,0)
                values[second]=_stable_nodal_sum(
                    (values[second],-target*edge_value),
                    caller,point,second,0)
            end
        end
        edge_node==node_count || error(
            "MeshFunctionSpaces: internal pyramid edge-node count mismatch")
        if gradient
            for node in 1:node_count,axis in 1:3
                cursor+=1;result[cursor]=gradients[axis,node]
            end
        else
            for node in 1:node_count
                cursor+=1;result[cursor]=values[node]
            end
        end
    end
    return result
end

function _serendipity_exponents(family::Symbol,p::Int)
    if family===:tri
        return [(q[1],q[2],0) for q in _tri_monomials(p,true)]
    elseif family===:tet
        return collect(_tet_monomials(p,true))
    elseif family===:qua
        result=NTuple{3,Int}[(0,0,0),(1,0,0),(1,1,0),(0,1,0)]
        for degree in 2:p
            append!(result,((degree,0,0),(degree,1,0),
                            (1,degree,0),(0,degree,0)))
        end
        return result
    elseif family===:hex
        result=NTuple{3,Int}[
            (0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1)]
        templates=((2,0,0),(2,0,1),(2,1,1),(2,1,0),
                   (0,2,0),(0,2,1),(1,2,1),(1,2,0),
                   (0,0,2),(0,1,2),(1,1,2),(1,0,2))
        for degree in 2:p,template in templates
            push!(result,ntuple(
                axis->template[axis]==2 ? degree : template[axis],3))
        end
        return result
    elseif family===:pri
        result=NTuple{3,Int}[
            (0,0,0),(1,0,0),(0,1,0),
            (0,0,1),(1,0,1),(0,1,1)]
        templates=((2,0,0),(2,0,1),(0,2,0),(0,2,1),
                   (0,0,2),(1,0,2),(0,1,2))
        for degree in 2:p,template in templates
            push!(result,ntuple(
                axis->template[axis]==2 ? degree : template[axis],3))
        end
        for first in 1:p-1
            second=p-first
            append!(result,((first,second,0),(first,second,1)))
        end
        return result
    end
    error("MeshFunctionSpaces: unsupported internal serendipity family $family")
end

@inline function _monomial_value(x,y,z,exponent)
    return x^exponent[1]*y^exponent[2]*z^exponent[3]
end

@inline function _monomial_derivative(x,y,z,exponent,axis::Int)
    degree=exponent[axis]
    degree==0 && return zero(x)
    first=axis==1 ? degree*x^(degree-1) : x^exponent[1]
    second=axis==2 ? degree*y^(degree-1) : y^exponent[2]
    third=axis==3 ? degree*z^(degree-1) : z^exponent[3]
    return first*second*third
end

function _build_polynomial_interpolant(element_type::Int)
    spec=msh_spec(element_type)
    exponents=_serendipity_exponents(spec.family,spec.order)
    length(exponents)==spec.nnodes || error(
        "MeshFunctionSpaces: internal serendipity mode-count mismatch for " *
        "element type $element_type")
    nodes=lagrange_nodes(element_type)
    return setprecision(BigFloat,_NODAL_COEFFICIENT_PRECISION) do
        count=length(exponents)
        vandermonde=Matrix{BigFloat}(undef,count,count)
        for column in 1:count,row in 1:count
            vandermonde[row,column]=_monomial_value(
                BigFloat(nodes[1,column]),BigFloat(nodes[2,column]),
                BigFloat(nodes[3,column]),exponents[row])
        end
        coefficients=inv(vandermonde)
        residual=max(opnorm(coefficients*vandermonde-I,Inf),
                     opnorm(vandermonde*coefficients-I,Inf))
        residual<=big"1e-50" || error(
            "MeshFunctionSpaces: high-precision interpolation inversion " *
            "failed for element type $element_type (residual=$residual)")
        float_coefficients=Float64.(coefficients)
        all(isfinite,float_coefficients) || error(
            "MeshFunctionSpaces: interpolation coefficients overflowed for " *
            "element type $element_type")
        return _PolynomialInterpolant(
            exponents,float_coefficients,coefficients,nodes)
    end
end

function _matching_nodal_point(nodes::Matrix{Float64},u::Float64,v::Float64,
                               w::Float64)
    for node in axes(nodes,2)
        nodes[1,node]==u && nodes[2,node]==v && nodes[3,node]==w && return node
    end
    return 0
end

function _polynomial_interpolant(element_type::Int)
    return lock(_NODAL_INTERPOLANT_LOCK) do
        get!(_POLYNOMIAL_INTERPOLANTS,element_type) do
            _build_polynomial_interpolant(element_type)
        end
    end
end

function _checked_big_nodal(value::BigFloat,caller::AbstractString,point::Int,
                            node::Int,component::Int)
    result=Float64(value)
    description=component==0 ?
        "evaluation point $point node $node basis value" :
        "evaluation point $point node $node gradient component $component"
    isfinite(result) || throw(ArgumentError(
        "$caller: $description is not Float64-representable"))
    result==0.0 && !iszero(value) && throw(ArgumentError(
        "$caller: $description is nonzero but below Float64 resolution"))
    return result
end

function _interpolant_dot(coefficients::Matrix{Float64},row::Int,terms)
    value=0.0
    scale=0.0
    for column in eachindex(terms)
        term=coefficients[row,column]*terms[column]
        value+=term
        scale+=abs(term)
    end
    if isfinite(value) && isfinite(scale) &&
       (scale==0.0 || abs(value)>
        _NODAL_ROUNDING_FACTOR*length(terms)*scale)
        return value
    end
    return nothing
end

function _big_interpolant_dot(coefficients::Matrix{BigFloat},row::Int,terms,
                              caller::AbstractString,point::Int,node::Int,
                              component::Int)
    value=sum(eachindex(terms);init=BigFloat(0)) do column
        coefficients[row,column]*terms[column]
    end
    return _checked_big_nodal(value,caller,point,node,component)
end

function _write_polynomial_interpolant!(result,coordinates,point_count,
                                        element_type,gradient,caller)
    data=_polynomial_interpolant(element_type)
    count=length(data.exponents)
    terms=Vector{Float64}(undef,count)
    derivative_terms=Matrix{Float64}(undef,3,count)
    fallback_nodes=Int[]
    fallback_indices=Int[]
    cursor=0
    for point in 1:point_count
        base=3point-2
        u,v,w=coordinates[base],coordinates[base+1],coordinates[base+2]
        if !gradient
            matching=_matching_nodal_point(data.nodes,u,v,w)
            if matching!=0
                for node in 1:count
                    cursor+=1
                    result[cursor]=node==matching ? 1.0 : 0.0
                end
                continue
            end
        end
        if gradient
            for axis in 1:3,(index,exponent) in enumerate(data.exponents)
                derivative_terms[axis,index]=
                    _monomial_derivative(u,v,w,exponent,axis)
            end
            point_start=cursor
            empty!(fallback_indices)
            for node in 1:count,axis in 1:3
                local_index=3node-3+axis
                cursor+=1
                value=_interpolant_dot(
                    data.coefficients,node,@view(derivative_terms[axis,:]))
                if value===nothing
                    push!(fallback_indices,local_index)
                    result[cursor]=0.0
                else
                    result[cursor]=value
                end
            end
            isempty(fallback_indices) || setprecision(
                    BigFloat,_NODAL_COEFFICIENT_PRECISION) do
                ub,vb,wb=BigFloat(u),BigFloat(v),BigFloat(w)
                big_terms=Matrix{BigFloat}(undef,3,count)
                for axis in 1:3,(index,exponent) in enumerate(data.exponents)
                    big_terms[axis,index]=
                        _monomial_derivative(ub,vb,wb,exponent,axis)
                end
                for local_index in fallback_indices
                    node=(local_index-1)÷3+1
                    axis=(local_index-1)%3+1
                    result[point_start+local_index]=_big_interpolant_dot(
                        data.high_precision_coefficients,node,
                        @view(big_terms[axis,:]),caller,point,node,axis)
                end
            end
        else
            for (index,exponent) in enumerate(data.exponents)
                terms[index]=_monomial_value(u,v,w,exponent)
            end
            point_start=cursor
            empty!(fallback_nodes)
            for node in 1:count
                cursor+=1
                value=_interpolant_dot(data.coefficients,node,terms)
                if value===nothing
                    push!(fallback_nodes,node)
                    result[cursor]=0.0
                else
                    result[cursor]=value
                end
            end
            isempty(fallback_nodes) || setprecision(
                    BigFloat,_NODAL_COEFFICIENT_PRECISION) do
                ub,vb,wb=BigFloat(u),BigFloat(v),BigFloat(w)
                big_terms=BigFloat[
                    _monomial_value(ub,vb,wb,exponent)
                    for exponent in data.exponents]
                for node in fallback_nodes
                    result[point_start+node]=_big_interpolant_dot(
                        data.high_precision_coefficients,node,big_terms,
                        caller,point,node,0)
                end
            end
        end
    end
    return result
end

function _orthogonal_values(order::Int,x)
    values=Vector{typeof(x)}(undef,order+1)
    derivatives=Vector{typeof(x)}(undef,order+1)
    values[1]=one(x);derivatives[1]=zero(x)
    order==0 && return values,derivatives
    values[2]=x;derivatives[2]=one(x)
    for degree in 2:order
        scale=one(x)/degree
        values[degree+1]=scale*((2degree-1)*x*values[degree]-
                               (degree-1)*values[degree-1])
        derivatives[degree+1]=scale*((2degree-1)*
            (values[degree]+x*derivatives[degree])-
            (degree-1)*derivatives[degree-1])
    end
    return values,derivatives
end

function _jacobi_values(order::Int,alpha::Int,beta::Int,x)
    values=Vector{typeof(x)}(undef,order+1)
    derivatives=Vector{typeof(x)}(undef,order+1)
    values[1]=one(x);derivatives[1]=zero(x)
    order==0 && return values,derivatives
    values[2]=((alpha-beta)+(alpha+beta+2)*x)/2
    derivatives[2]=(alpha+beta+2)*one(x)/2
    for degree in 2:order
        twice=2degree+alpha+beta
        denominator=2degree*(degree+alpha+beta)*(twice-2)
        first_scale=(twice-1)
        linear=(twice*(twice-2)*x+alpha^2-beta^2)
        previous_scale=2*(degree+alpha-1)*(degree+beta-1)*twice
        values[degree+1]=(
            first_scale*linear*values[degree]-
            previous_scale*values[degree-1])/denominator
        derivatives[degree+1]=(
            first_scale*(twice*(twice-2)*values[degree]+
                         linear*derivatives[degree])-
            previous_scale*derivatives[degree-1])/denominator
    end
    return values,derivatives
end

@inline _pyramid_valid_indices(p::Int,serendipity::Bool,i::Int,j::Int)=
    !serendipity || i+j<=p || (i+j==p+1 && (i==1 || j==1))

function _pyramid_modes(p::Int,serendipity::Bool)
    modes=NTuple{3,Int}[]
    for i in 0:p,j in 0:p
        _pyramid_valid_indices(p,serendipity,i,j) || continue
        for k in 0:p-max(i,j)
            push!(modes,(i,j,k))
        end
    end
    return modes
end

function _pyramid_mode_data(p::Int,serendipity::Bool,u,v,w,
                            modes::Vector{NTuple{3,Int}},gradient::Bool)
    values=Vector{typeof(u)}(undef,length(modes))
    gradients=gradient ? Matrix{typeof(u)}(undef,3,length(modes)) :
                         Matrix{typeof(u)}(undef,0,0)
    q=one(w)-w
    reduced_u=iszero(q) ? zero(u) : u/q
    reduced_v=iszero(q) ? zero(v) : v/q
    u_values,u_derivatives=_orthogonal_values(p,reduced_u)
    v_values,v_derivatives=_orthogonal_values(p,reduced_v)
    jacobi=Vector{Tuple{Vector{typeof(u)},Vector{typeof(u)}}}(undef,p+1)
    transformed_w=2w-one(w)
    for m in 0:p
        jacobi[m+1]=_jacobi_values(p-m,2m+2,0,transformed_w)
    end
    for (index,(i,j,k)) in enumerate(modes)
        m=max(i,j)
        w_values,w_derivatives=jacobi[m+1]
        polynomial=w_values[k+1]
        values[index]=u_values[i+1]*v_values[j+1]*polynomial*q^m
        gradient || continue
        if m==0
            gradients[1,index]=zero(u)
            gradients[2,index]=zero(u)
            gradients[3,index]=2w_derivatives[k+1]
        elseif m==1
            if i==0
                gradients[1,index]=zero(u)
                gradients[2,index]=polynomial
                gradients[3,index]=2v*w_derivatives[k+1]
            elseif j==0
                gradients[1,index]=polynomial
                gradients[2,index]=zero(u)
                gradients[3,index]=2u*w_derivatives[k+1]
            else
                gradients[1,index]=reduced_v*polynomial
                gradients[2,index]=reduced_u*polynomial
                gradients[3,index]=reduced_u*reduced_v*polynomial+
                    2reduced_u*v*w_derivatives[k+1]
            end
        else
            power_m2=q^(m-2)
            power_m1=power_m2*q
            gradients[1,index]=u_derivatives[i+1]*v_values[j+1]*
                               polynomial*power_m1
            gradients[2,index]=u_values[i+1]*v_derivatives[j+1]*
                               polynomial*power_m1
            gradients[3,index]=polynomial*power_m2*(
                u*u_derivatives[i+1]*v_values[j+1]+
                v*u_values[i+1]*v_derivatives[j+1])+
                u_values[i+1]*v_values[j+1]*power_m1*(
                    2q*w_derivatives[k+1]-m*polynomial)
        end
    end
    return values,gradients
end

function _build_pyramid_interpolant(element_type::Int)
    spec=msh_spec(element_type)
    modes=_pyramid_modes(spec.order,spec.serendipity)
    length(modes)==spec.nnodes || error(
        "MeshFunctionSpaces: internal pyramid mode-count mismatch for " *
        "element type $element_type")
    nodes=lagrange_nodes(element_type)
    count=length(modes)
    vandermonde=Matrix{Float64}(undef,count,count)
    for column in 1:count
        values,_=_pyramid_mode_data(
            spec.order,spec.serendipity,nodes[1,column],nodes[2,column],
            nodes[3,column],modes,false)
        vandermonde[:,column]=values
    end
    coefficients=inv(vandermonde)
    all(isfinite,coefficients) || error(
        "MeshFunctionSpaces: pyramid interpolation coefficients overflowed " *
        "for element type $element_type")
    residual=max(opnorm(coefficients*vandermonde-I,Inf),
                 opnorm(vandermonde*coefficients-I,Inf))
    residual<=2e-8 || error(
        "MeshFunctionSpaces: pyramid interpolation inversion failed for " *
        "element type $element_type (residual=$residual)")
    extended_coefficients=setprecision(
        BigFloat,_NODAL_COEFFICIENT_PRECISION) do
        BigFloat.(coefficients)
    end
    return _PyramidInterpolant(
        modes,coefficients,extended_coefficients,nodes)
end

function _pyramid_interpolant(element_type::Int)
    return lock(_NODAL_INTERPOLANT_LOCK) do
        get!(_PYRAMID_INTERPOLANTS,element_type) do
            _build_pyramid_interpolant(element_type)
        end
    end
end

function _write_pyramid_interpolant!(result,coordinates,point_count,
                                     element_type,gradient,caller)
    spec=msh_spec(element_type)
    data=_pyramid_interpolant(element_type)
    count=length(data.modes)
    fallback_nodes=Int[]
    fallback_indices=Int[]
    cursor=0
    for point in 1:point_count
        base=3point-2
        u,v,w=coordinates[base],coordinates[base+1],coordinates[base+2]
        if !gradient
            matching=_matching_nodal_point(data.nodes,u,v,w)
            if matching!=0
                for node in 1:count
                    cursor+=1
                    result[cursor]=node==matching ? 1.0 : 0.0
                end
                continue
            end
        end
        values,gradients=_pyramid_mode_data(
            spec.order,spec.serendipity,u,v,w,data.modes,gradient)
        if gradient
            point_start=cursor
            empty!(fallback_indices)
            for node in 1:count,axis in 1:3
                terms=@view gradients[axis,:]
                local_index=3node-3+axis
                cursor+=1
                value=_interpolant_dot(data.coefficients,node,terms)
                if value===nothing
                    push!(fallback_indices,local_index)
                    result[cursor]=0.0
                else
                    result[cursor]=value
                end
            end
            isempty(fallback_indices) || setprecision(
                    BigFloat,_NODAL_COEFFICIENT_PRECISION) do
                _,big_gradients=_pyramid_mode_data(
                    spec.order,spec.serendipity,BigFloat(u),BigFloat(v),
                    BigFloat(w),data.modes,true)
                for local_index in fallback_indices
                    node=(local_index-1)÷3+1
                    axis=(local_index-1)%3+1
                    result[point_start+local_index]=_big_interpolant_dot(
                        data.extended_coefficients,node,
                        @view(big_gradients[axis,:]),caller,point,node,axis)
                end
            end
        else
            point_start=cursor
            empty!(fallback_nodes)
            for node in 1:count
                cursor+=1
                value=_interpolant_dot(data.coefficients,node,values)
                if value===nothing
                    push!(fallback_nodes,node)
                    result[cursor]=0.0
                else
                    result[cursor]=value
                end
            end
            isempty(fallback_nodes) || setprecision(
                    BigFloat,_NODAL_COEFFICIENT_PRECISION) do
                big_values,_=_pyramid_mode_data(
                    spec.order,spec.serendipity,BigFloat(u),BigFloat(v),
                    BigFloat(w),data.modes,false)
                for node in fallback_nodes
                    result[point_start+node]=_big_interpolant_dot(
                        data.extended_coefficients,node,big_values,
                        caller,point,node,0)
                end
            end
        end
    end
    return result
end

function _write_higher_order_nodal!(result,coordinates,point_count::Int,
                                    basis_type::Int,gradient::Bool,
                                    caller::AbstractString)
    spec=msh_spec(basis_type)
    isempty(result) && return result
    if spec.family===:pnt || spec.order==1
        return _write_nodal_basis!(
            result,coordinates,point_count,spec.family,gradient,caller)
    elseif spec.order==0
        fill!(result,0.0)
        gradient || fill!(result,1.0)
        return result
    elseif spec.family===:pyr && spec.serendipity && spec.order>2
        return _write_serendipity_pyramid!(
            result,coordinates,point_count,spec.order,gradient,caller)
    elseif spec.family===:pyr
        return _write_pyramid_interpolant!(
            result,coordinates,point_count,basis_type,gradient,caller)
    elseif spec.serendipity
        return _write_polynomial_interpolant!(
            result,coordinates,point_count,basis_type,gradient,caller)
    elseif spec.family===:lin
        return _write_complete_line!(
            result,coordinates,point_count,spec.order,gradient,caller)
    elseif spec.family===:tri
        return _write_complete_simplex!(
            result,coordinates,point_count,spec.order,gradient,
            _tri_monomials(spec.order),2,caller)
    elseif spec.family===:tet
        return _write_complete_simplex!(
            result,coordinates,point_count,spec.order,gradient,
            _tet_monomials(spec.order),3,caller)
    elseif spec.family===:qua
        return _write_complete_tensor!(
            result,coordinates,point_count,spec.order,gradient,
            _qua_monomials(spec.order),Val(2),caller)
    elseif spec.family===:hex
        return _write_complete_tensor!(
            result,coordinates,point_count,spec.order,gradient,
            _hex_monomials(spec.order),Val(3),caller)
    elseif spec.family===:pri
        return _write_complete_prism!(
            result,coordinates,point_count,spec.order,gradient,caller)
    end
    error("MeshFunctionSpaces: unsupported internal nodal family $(spec.family)")
end

function _nodal_bubble_count(basis_type::Int)
    spec=msh_spec(basis_type)
    p=spec.order
    spec.family===:pnt && return 0
    p==0 && return 1
    spec.serendipity && return 0
    spec.family===:lin && return p-1
    spec.family===:tri && return (p-1)*(p-2)÷2
    spec.family===:qua && return (p-1)^2
    spec.family===:tet && return (p-1)*(p-2)*(p-3)÷6
    spec.family===:pri && return (p-1)^2*(p-2)÷2
    spec.family===:hex && return (p-1)^3
    spec.family===:pyr && return (p-2)*(p-1)*(2p-3)÷6
    error("MeshFunctionSpaces: unsupported internal nodal family $(spec.family)")
end
