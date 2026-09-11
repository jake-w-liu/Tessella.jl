"""
    MeshFunctionSpaces

Reference finite-element function spaces and global degree-of-freedom keys.
Actual- and explicit-order nodal functions cover every fixed Point, Line,
Triangle, Quadrangle, Tetrahedron, Hexahedron, Prism, and Pyramid type.
Hierarchical H1 and lowest-order H(curl) functions, cached orientations, and
populated keys cover the finalized linear-simplex [`Mesh`](@ref). The
implementation follows Gmsh 4.15.2's reference coordinates, output layout,
lexicographic orientation indices, nodal keys, and edge keys. Trihedron bases
fail explicitly.
"""
module MeshFunctionSpaces

using ..MeshTypes: Mesh
using ..Elements: lagrange_nodes, msh_spec, _hex_monomials,
                  _line_monomials, _msh_type_exact, _pri_monomials,
                  _PYR_EDGES, _qua_monomials, _tet_monomials,
                  _tri_monomials
using ..MeshEntityTopology: MeshEdgeTopology, _mesh_edges,
                            _simplex_edge_patterns
using ..MeshPointLocation: mesh_element_block, mesh_element_record
using ..MeshReferenceGeometry: _Exact, _checked_element_tag,
                               _checked_element_type,
                               _checked_local_coordinates, _exact_to_float,
                               _stable_midpoint

export mesh_basis_functions, mesh_basis_orientation,
       mesh_basis_orientations, mesh_key_dimension, mesh_keys,
       mesh_keys_for_element, mesh_keys_information,
       mesh_number_of_keys, mesh_number_of_orientations

struct _FunctionSpace
    kind::Symbol
    components::Int
    hierarchical::Bool
    key_dimension::Int
    key_order::Int
end

const _LAGRANGE=_FunctionSpace(:lagrange,1,false,0,-1)
const _GRAD_LAGRANGE=_FunctionSpace(:grad_lagrange,3,false,0,-1)
const _H1_1=_FunctionSpace(:lagrange,1,true,0,1)
const _GRAD_H1_1=_FunctionSpace(:grad_lagrange,3,true,0,1)
const _HCURL_0=_FunctionSpace(:hcurl,3,true,1,0)
const _CURL_HCURL_0=_FunctionSpace(:curl_hcurl,3,true,1,0)
const _FIXED_NODAL_FAMILIES=(:pnt,:lin,:tri,:qua,:tet,:hex,:pri,:pyr)
const _QUADRANGLE_SIGNS=((-1.0,-1.0),(1.0,-1.0),
                         (1.0,1.0),(-1.0,1.0))
const _HEXAHEDRON_SIGNS=((-1.0,-1.0,-1.0),(1.0,-1.0,-1.0),
                         (1.0,1.0,-1.0),(-1.0,1.0,-1.0),
                         (-1.0,-1.0,1.0),(1.0,-1.0,1.0),
                         (1.0,1.0,1.0),(-1.0,1.0,1.0))

function _explicit_nodal_space(name::String,prefix::String,kind::Symbol,
                               caller::AbstractString)
    startswith(name,prefix) || return nothing
    suffix=SubString(name,nextind(name,lastindex(prefix)))
    isempty(suffix) && return nothing
    order=0
    digits=0
    for character in suffix
        '0'<=character<='9' || throw(ArgumentError(
            "$caller: unsupported function_space_type $(repr(name)); " *
            "explicit Lagrange orders must use decimal digits"))
        digits+=1
        digits<=2 || throw(ArgumentError(
            "$caller: explicit interpolation order is outside the supported " *
            "Gmsh 4.15.2 range 0:10"))
        order=10order+(Int(character)-Int('0'))
        order<=10 || throw(ArgumentError(
            "$caller: interpolation order $order is outside the supported " *
            "Gmsh 4.15.2 range 0:10"))
    end
    components=kind===:lagrange ? 1 : 3
    return _FunctionSpace(kind,components,false,0,order)
end

function _function_space(value,caller::AbstractString)
    value isa AbstractString || throw(ArgumentError(
        "$caller: function_space_type must be a string"))
    name=String(value)
    occursin('\0',name) && throw(ArgumentError(
        "$caller: function_space_type must not contain NUL"))
    (name=="Lagrange" || name=="IsoParametric") && return _LAGRANGE
    (name=="GradLagrange" || name=="GradIsoParametric") &&
        return _GRAD_LAGRANGE
    explicit=_explicit_nodal_space(
        name,"GradLagrange",:grad_lagrange,caller)
    explicit===nothing || return explicit
    explicit=_explicit_nodal_space(name,"Lagrange",:lagrange,caller)
    explicit===nothing || return explicit
    name=="H1Legendre1" && return _H1_1
    name=="GradH1Legendre1" && return _GRAD_H1_1
    name=="HcurlLegendre0" && return _HCURL_0
    name=="CurlHcurlLegendre0" && return _CURL_HCURL_0
    throw(ArgumentError(
        "$caller: unsupported function_space_type $(repr(name)); supported " *
        "nodal spaces are Lagrange, IsoParametric, Lagrange0 through " *
        "Lagrange10, GradLagrange, GradIsoParametric, and GradLagrange0 " *
        "through GradLagrange10; supported hierarchical spaces are " *
        "H1Legendre1, GradH1Legendre1, HcurlLegendre0, and " *
        "CurlHcurlLegendre0"))
end

function _basis_element_contract(element_type_value,space::_FunctionSpace,
                                 caller::AbstractString)
    element_type=_checked_element_type(element_type_value,caller)
    if space.hierarchical
        element_type in (1,2,4) || throw(ArgumentError(
            "$caller: element type $element_type is not a supported linear " *
            "segment, triangle, or tetrahedron for hierarchical spaces"))
        return element_type,msh_spec(element_type).family,
               _vertex_count(element_type),element_type
    end

    spec=msh_spec(element_type)
    spec.family===:trih && throw(ArgumentError(
        "$caller: Trihedron element type $element_type has no nodal basis " *
        "in Gmsh 4.15.2"))
    spec.family in _FIXED_NODAL_FAMILIES || throw(ArgumentError(
        "$caller: element type $element_type belongs to unsupported " *
        "$(spec.family) family"))
    basis_type=if space.key_order==-1
        element_type
    elseif spec.family===:pnt
        15
    else
        _msh_type_exact(spec.family,space.key_order,false,caller)
    end
    basis_spec=msh_spec(basis_type)
    return element_type,spec.family,basis_spec.nnodes,basis_type
end

@inline function _vertex_count(element_type::Int)
    element_type==1 && return 2
    element_type==2 && return 3
    return 4
end

@inline function _edge_count(element_type::Int)
    element_type==1 && return 1
    element_type==2 && return 3
    return 6
end

@inline _orientation_count(element_type::Int)=
    element_type==1 ? 2 : element_type==2 ? 6 : 24

function _checked_orientation_sequence(values,norientations::Int,
                                       hierarchical::Bool,
                                       caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: wanted_orientations must be a vector or tuple of integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    if isempty(values)
        return hierarchical ? collect(0:norientations-1) : Int[0]
    end
    result=Vector{Int}(undef,length(values))
    seen=Set{Int}()
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: wanted_orientations[$index] must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: wanted_orientations[$index] must not be Bool"))
        orientation=try
            Int(value)
        catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: wanted_orientations[$index] exceeds Int bounds"))
        end
        0<=orientation<norientations || throw(ArgumentError(
            "$caller: orientation $orientation is outside " *
            "0:$(norientations-1)"))
        orientation in seen && throw(ArgumentError(
            "$caller: duplicate wanted orientation $orientation"))
        push!(seen,orientation)
        result[index]=orientation
    end
    return result
end

function _checked_result_length(factors::Int...;caller::AbstractString)
    result=1
    try
        for factor in factors
            result=Base.checked_mul(result,factor)
        end
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "$caller: requested basis result exceeds the platform array range"))
    end
    return result
end

function _origin_barycentric(values,caller::AbstractString,point::Int)
    fast=1.0
    scale=1.0
    for value in values
        fast-=value
        scale=max(scale,abs(value))
    end
    if isfinite(fast) && abs(fast)>64eps(Float64)*scale
        return fast
    end
    exact=_Exact(1)
    for value in values
        exact-=_Exact(value)
    end
    result=_exact_to_float(
        exact,caller,"evaluation point $point origin barycentric coordinate")
    result==0.0 && exact!=0 && throw(ArgumentError(
        "$caller: evaluation point $point origin barycentric coordinate is " *
        "nonzero but below Float64 resolution"))
    return result
end

function _exact_basis_result(value::_Exact,caller::AbstractString,
                             point::Int,node::Int,component::Int)
    description=component==0 ?
        "evaluation point $point node $node basis value" :
        "evaluation point $point node $node gradient component $component"
    result=_exact_to_float(value,caller,description)
    result==0.0 && value!=0 && throw(ArgumentError(
        "$caller: $description is nonzero but below Float64 resolution"))
    return result
end

function _basis_product(scale::Float64,factors::Tuple,
                        caller::AbstractString,point::Int,node::Int,
                        component::Int)
    value=scale
    zero=scale==0.0
    for factor in factors
        zero|=factor==0.0
        value*=factor
    end
    isfinite(value) && (value!=0.0 || zero) && return value
    exact=_Exact(scale)
    for factor in factors
        exact*=_Exact(factor)
    end
    return _exact_basis_result(exact,caller,point,node,component)
end

function _basis_quotient(scale::Float64,numerators::Tuple,
                         denominators::Tuple,caller::AbstractString,
                         point::Int,node::Int,component::Int)
    value=scale
    zero=scale==0.0
    for factor in numerators
        zero|=factor==0.0
        value*=factor
    end
    for factor in denominators
        factor!=0.0 || error(
            "MeshFunctionSpaces: internal zero nodal-basis denominator")
        value/=factor
    end
    isfinite(value) && (value!=0.0 || zero) && return value
    exact=_Exact(scale)
    for factor in numerators
        exact*=_Exact(factor)
    end
    for factor in denominators
        exact/=_Exact(factor)
    end
    return _exact_basis_result(exact,caller,point,node,component)
end

@inline function _pyramid_factor(q::Float64,w::Float64,sign::Float64,
                                 coordinate::Float64)
    value=q+sign*coordinate
    scale=max(1.0,abs(w),abs(coordinate))
    return value,isfinite(value) && abs(value)>64eps(Float64)*scale
end

function _pyramid_basis_value(sign_u::Float64,sign_v::Float64,
                              u::Float64,v::Float64,w::Float64,q::Float64,
                              caller::AbstractString,point::Int,node::Int)
    first,first_stable=_pyramid_factor(q,w,sign_u,u)
    second,second_stable=_pyramid_factor(q,w,sign_v,v)
    if first_stable && second_stable
        return _basis_quotient(
            0.25,(first,second),(q,),caller,point,node,0)
    end
    exact_q=_Exact(1)-_Exact(w)
    exact_first=exact_q+_Exact(sign_u)*_Exact(u)
    exact_second=exact_q+_Exact(sign_v)*_Exact(v)
    return _exact_basis_result(
        exact_first*exact_second/(4exact_q),caller,point,node,0)
end

function _pyramid_horizontal_gradient(sign::Float64,factor_sign::Float64,
                                      coordinate::Float64,w::Float64,
                                      q::Float64,caller::AbstractString,
                                      point::Int,node::Int,component::Int)
    factor,stable=_pyramid_factor(q,w,factor_sign,coordinate)
    stable && return _basis_quotient(
        0.25sign,(factor,),(q,),caller,point,node,component)
    exact_q=_Exact(1)-_Exact(w)
    exact_factor=exact_q+_Exact(factor_sign)*_Exact(coordinate)
    return _exact_basis_result(
        _Exact(sign)*exact_factor/(4exact_q),
        caller,point,node,component)
end

function _pyramid_vertical_gradient(sign_u::Float64,sign_v::Float64,
                                    u::Float64,v::Float64,w::Float64,
                                    q::Float64,
                                    caller::AbstractString,point::Int,
                                    node::Int)
    first=sign_u*sign_v*u*v
    second=q*q
    numerator=first-second
    denominator=4second
    scale=abs(first)+abs(second)
    if isfinite(numerator) && isfinite(denominator) && denominator!=0.0 &&
       isfinite(scale) &&
       (scale==0.0 || abs(numerator)>64eps(Float64)*scale)
        value=numerator/denominator
        isfinite(value) && value!=0.0 && return value
    end
    exact_q=_Exact(1)-_Exact(w)
    exact=(_Exact(sign_u*sign_v)*_Exact(u)*_Exact(v)-exact_q^2)/
          (4exact_q^2)
    return _exact_basis_result(exact,caller,point,node,3)
end

@inline _first_order_values(::Val{:pnt},u,v,w,caller,point)=(1.0,)
@inline _first_order_gradients(::Val{:pnt},u,v,w,caller,point)=
    ((0.0,0.0,0.0),)

@inline function _first_order_values(::Val{:lin},u,v,w,caller,point)
    return ((1.0-u)/2,(1.0+u)/2)
end

@inline _first_order_gradients(::Val{:lin},u,v,w,caller,point)=
    ((-0.5,0.0,0.0),(0.5,0.0,0.0))

@inline function _first_order_values(::Val{:tri},u,v,w,caller,point)
    return (_origin_barycentric((u,v),caller,point),u,v)
end

@inline _first_order_gradients(::Val{:tri},u,v,w,caller,point)=
    ((-1.0,-1.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0))

@inline function _first_order_values(::Val{:tet},u,v,w,caller,point)
    return (_origin_barycentric((u,v,w),caller,point),u,v,w)
end

@inline _first_order_gradients(::Val{:tet},u,v,w,caller,point)=
    ((-1.0,-1.0,-1.0),(1.0,0.0,0.0),
     (0.0,1.0,0.0),(0.0,0.0,1.0))

function _first_order_values(::Val{:qua},u,v,w,caller,point)
    return ntuple(4) do node
        sign_u,sign_v=_QUADRANGLE_SIGNS[node]
        _basis_product(
            0.25,(1.0+sign_u*u,1.0+sign_v*v),
            caller,point,node,0)
    end
end

function _first_order_gradients(::Val{:qua},u,v,w,caller,point)
    return ntuple(4) do node
        sign_u,sign_v=_QUADRANGLE_SIGNS[node]
        (_basis_product(
             0.25sign_u,(1.0+sign_v*v,),caller,point,node,1),
         _basis_product(
             0.25sign_v,(1.0+sign_u*u,),caller,point,node,2),
         0.0)
    end
end

function _first_order_values(::Val{:hex},u,v,w,caller,point)
    return ntuple(8) do node
        sign_u,sign_v,sign_w=_HEXAHEDRON_SIGNS[node]
        _basis_product(
            0.125,(1.0+sign_u*u,1.0+sign_v*v,1.0+sign_w*w),
            caller,point,node,0)
    end
end

function _first_order_gradients(::Val{:hex},u,v,w,caller,point)
    return ntuple(8) do node
        sign_u,sign_v,sign_w=_HEXAHEDRON_SIGNS[node]
        first=1.0+sign_u*u
        second=1.0+sign_v*v
        third=1.0+sign_w*w
        (_basis_product(
             0.125sign_u,(second,third),caller,point,node,1),
         _basis_product(
             0.125sign_v,(first,third),caller,point,node,2),
         _basis_product(
             0.125sign_w,(first,second),caller,point,node,3))
    end
end

function _first_order_values(::Val{:pri},u,v,w,caller,point)
    triangle=(_origin_barycentric((u,v),caller,point),u,v)
    return ntuple(6) do node
        triangle_node=(node-1)%3+1
        sign_w=node<=3 ? -1.0 : 1.0
        _basis_product(
            0.5,(triangle[triangle_node],1.0+sign_w*w),
            caller,point,node,0)
    end
end

function _first_order_gradients(::Val{:pri},u,v,w,caller,point)
    triangle=(_origin_barycentric((u,v),caller,point),u,v)
    triangle_gradients=((-1.0,-1.0),(1.0,0.0),(0.0,1.0))
    return ntuple(6) do node
        triangle_node=(node-1)%3+1
        sign_w=node<=3 ? -1.0 : 1.0
        factor=1.0+sign_w*w
        gradient=triangle_gradients[triangle_node]
        (_basis_product(
             0.5,(gradient[1],factor),caller,point,node,1),
         _basis_product(
             0.5,(gradient[2],factor),caller,point,node,2),
         _basis_product(
             0.5sign_w,(triangle[triangle_node],),caller,point,node,3))
    end
end

function _first_order_values(::Val{:pyr},u,v,w,caller,point)
    q=1.0-w
    q==0.0 && return (0.0,0.0,0.0,0.0,1.0)
    base=ntuple(4) do node
        sign_u,sign_v=_QUADRANGLE_SIGNS[node]
        _pyramid_basis_value(
            sign_u,sign_v,u,v,w,q,caller,point,node)
    end
    return (base...,w)
end

function _first_order_gradients(::Val{:pyr},u,v,w,caller,point)
    q=1.0-w
    if q==0.0
        return ((-0.25,-0.25,-0.25),(0.25,-0.25,-0.25),
                (0.25,0.25,-0.25),(-0.25,0.25,-0.25),
                (0.0,0.0,1.0))
    end
    base=ntuple(4) do node
        sign_u,sign_v=_QUADRANGLE_SIGNS[node]
        (_pyramid_horizontal_gradient(
             sign_u,sign_v,v,w,q,caller,point,node,1),
         _pyramid_horizontal_gradient(
             sign_v,sign_u,u,w,q,caller,point,node,2),
         _pyramid_vertical_gradient(
             sign_u,sign_v,u,v,w,q,caller,point,node))
    end
    return (base...,(0.0,0.0,1.0))
end

function _write_nodal_family!(result,coordinates,point_count::Int,
                              family::Val,gradient::Bool,
                              caller::AbstractString)
    cursor=0
    if gradient
        @inbounds for point in 1:point_count
            offset=3point-2
            gradients=_first_order_gradients(
                family,coordinates[offset],coordinates[offset+1],
                coordinates[offset+2],caller,point)
            for value in gradients,component in 1:3
                cursor+=1
                result[cursor]=value[component]
            end
        end
    else
        @inbounds for point in 1:point_count
            offset=3point-2
            values=_first_order_values(
                family,coordinates[offset],coordinates[offset+1],
                coordinates[offset+2],caller,point)
            for value in values
                cursor+=1
                result[cursor]=value
            end
        end
    end
    cursor==length(result) || error(
        "MeshFunctionSpaces: internal nodal-basis result length mismatch")
    return result
end

function _write_nodal_basis!(result,coordinates,point_count::Int,
                             family::Symbol,gradient::Bool,
                             caller::AbstractString)
    family===:pnt && return _write_nodal_family!(
        result,coordinates,point_count,Val(:pnt),gradient,caller)
    family===:lin && return _write_nodal_family!(
        result,coordinates,point_count,Val(:lin),gradient,caller)
    family===:tri && return _write_nodal_family!(
        result,coordinates,point_count,Val(:tri),gradient,caller)
    family===:qua && return _write_nodal_family!(
        result,coordinates,point_count,Val(:qua),gradient,caller)
    family===:tet && return _write_nodal_family!(
        result,coordinates,point_count,Val(:tet),gradient,caller)
    family===:hex && return _write_nodal_family!(
        result,coordinates,point_count,Val(:hex),gradient,caller)
    family===:pri && return _write_nodal_family!(
        result,coordinates,point_count,Val(:pri),gradient,caller)
    family===:pyr && return _write_nodal_family!(
        result,coordinates,point_count,Val(:pyr),gradient,caller)
    error("MeshFunctionSpaces: unsupported internal nodal family $family")
end

include("HigherOrderNodal.jl")

function _barycentric_coordinates(element_type::Int,u::Float64,v::Float64,
                                  w::Float64,caller::AbstractString,point::Int)
    if element_type==1
        return ((1.0-u)/2,(1.0+u)/2)
    elseif element_type==2
        return (_origin_barycentric((u,v),caller,point),u,v)
    end
    return (_origin_barycentric((u,v,w),caller,point),u,v,w)
end

@inline function _reference_gradients(element_type::Int)
    element_type==1 && return ((-0.5,0.0,0.0),(0.5,0.0,0.0))
    element_type==2 && return (
        (-1.0,-1.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0))
    return ((-1.0,-1.0,-1.0),(1.0,0.0,0.0),
            (0.0,1.0,0.0),(0.0,0.0,1.0))
end

@inline function _canonical_edge_pairs(element_type::Int)
    element_type==1 && return ((1,2),)
    element_type==2 && return ((1,2),(2,3),(1,3))
    return ((1,2),(2,3),(1,3),(1,4),(3,4),(2,4))
end

function _orientation_permutation(vertex_count::Int,orientation::Int)
    remaining=collect(1:vertex_count)
    permutation=Vector{Int}(undef,vertex_count)
    rank=orientation
    @inbounds for position in 1:vertex_count
        block=factorial(vertex_count-position)
        selection=rank÷block+1
        rank%=block
        permutation[position]=remaining[selection]
        deleteat!(remaining,selection)
    end
    return permutation
end

function _exact_barycentric_coordinate(element_type::Int,vertex::Int,
                                       u::Float64,v::Float64,w::Float64)
    if element_type==1
        return vertex==1 ? (_Exact(1)-_Exact(u))/2 :
                           (_Exact(1)+_Exact(u))/2
    elseif element_type==2
        vertex==1 && return _Exact(1)-_Exact(u)-_Exact(v)
        return vertex==2 ? _Exact(u) : _Exact(v)
    end
    vertex==1 && return _Exact(1)-_Exact(u)-_Exact(v)-_Exact(w)
    vertex==2 && return _Exact(u)
    vertex==3 && return _Exact(v)
    return _Exact(w)
end

function _scaled_whitney_component(element_type::Int,
                                   first_vertex::Int,lambda_i::Float64,
                                   gradient_i::Float64,
                                   second_vertex::Int,lambda_j::Float64,
                                   gradient_j::Float64,
                                   u::Float64,v::Float64,w::Float64,
                                   caller::AbstractString,point::Int,
                                   edge::Int,component::Int)
    first_term=lambda_i*gradient_j
    second_term=lambda_j*gradient_i
    difference=first_term-second_term
    value=2.0*difference
    scale=max(abs(first_term),abs(second_term))
    if isfinite(value) && (scale==0.0 ||
                           abs(difference)>64eps(Float64)*scale)
        return value
    end
    exact_i=_exact_barycentric_coordinate(
        element_type,first_vertex,u,v,w)
    exact_j=_exact_barycentric_coordinate(
        element_type,second_vertex,u,v,w)
    exact=2*(exact_i*_Exact(gradient_j)-exact_j*_Exact(gradient_i))
    result=_exact_to_float(
        exact,caller,"evaluation point $point edge $edge component $component")
    result==0.0 && exact!=0 && throw(ArgumentError(
        "$caller: evaluation point $point edge $edge component $component is " *
        "nonzero but below Float64 resolution"))
    return result
end

@inline function _curl_component(gradient_i,gradient_j,component::Int)
    component==1 && return 4.0*(gradient_i[2]*gradient_j[3]-
                               gradient_i[3]*gradient_j[2])
    component==2 && return 4.0*(gradient_i[3]*gradient_j[1]-
                               gradient_i[1]*gradient_j[3])
    return 4.0*(gradient_i[1]*gradient_j[2]-
                gradient_i[2]*gradient_j[1])
end

"""Return supported Gmsh-shaped reference basis values and orientation count."""
function mesh_basis_functions(element_type_value,local_coord,
                              function_space_type,wanted_orientations=Int32[];
                              caller::AbstractString="mesh_basis_functions")
    element_type=_checked_element_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    element_type,_,nodal_count,basis_type=
        _basis_element_contract(element_type,space,caller)
    coordinates,point_count=_checked_local_coordinates(local_coord,caller)
    total_orientations=space.hierarchical ?
        _orientation_count(element_type) : 1
    orientations=_checked_orientation_sequence(
        wanted_orientations,total_orientations,space.hierarchical,caller)
    function_count=space.hierarchical ?
        (space.key_dimension==0 ? _vertex_count(element_type) :
                                  _edge_count(element_type)) : nodal_count
    result_length=_checked_result_length(
        length(orientations),point_count,function_count,space.components;
        caller=caller)
    result=Vector{Float64}(undef,result_length)
    if !space.hierarchical
        _write_higher_order_nodal!(
            result,coordinates,point_count,basis_type,
            space.kind===:grad_lagrange,caller)
        return Int32(space.components),result,Int32(1)
    end
    gradients=_reference_gradients(element_type)
    pairs=_canonical_edge_pairs(element_type)
    cursor=0
    @inbounds for orientation in orientations
        permutation=space.hierarchical ?
            _orientation_permutation(_vertex_count(element_type),orientation) :
            Int[]
        for point in 1:point_count
            base=3point-2
            lambda=_barycentric_coordinates(
                element_type,coordinates[base],coordinates[base+1],
                coordinates[base+2],caller,point)
            if space.kind==:lagrange
                for value in lambda
                    cursor+=1
                    result[cursor]=value
                end
            elseif space.kind==:grad_lagrange
                for gradient in gradients,component in 1:3
                    cursor+=1
                    result[cursor]=gradient[component]
                end
            else
                for (edge,(first_vertex,second_vertex)) in enumerate(pairs)
                    sign=permutation[first_vertex]<permutation[second_vertex] ?
                        1.0 : -1.0
                    for component in 1:3
                        cursor+=1
                        result[cursor]=sign*(space.kind==:hcurl ?
                            _scaled_whitney_component(
                                element_type,first_vertex,
                                lambda[first_vertex],
                                gradients[first_vertex][component],
                                second_vertex,
                                lambda[second_vertex],
                                gradients[second_vertex][component],
                                coordinates[base],coordinates[base+1],
                                coordinates[base+2],
                                caller,point,edge,component) :
                            _curl_component(
                                gradients[first_vertex],
                                gradients[second_vertex],component))
                    end
                end
            end
        end
    end
    return Int32(space.components),result,Int32(total_orientations)
end

"""Return the number of reference orientations for a supported function space."""
function mesh_number_of_orientations(element_type_value,function_space_type;
                                     caller::AbstractString=
                                         "mesh_number_of_orientations")
    element_type=_checked_element_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    element_type,_,_,_=_basis_element_contract(element_type,space,caller)
    return Int32(space.hierarchical ? _orientation_count(element_type) : 1)
end

function _orientation_rank(node_tags)
    count=length(node_tags)
    rank=0
    @inbounds for first in 1:count-1
        smaller=0
        for second in first+1:count
            smaller+=node_tags[second]<node_tags[first]
        end
        rank+=smaller*factorial(count-first)
    end
    return Int32(rank)
end

function _cell_orientation(cells::AbstractMatrix{Int32},cell::Int,
                           space::_FunctionSpace)
    space.hierarchical || return Int32(0)
    return _orientation_rank(@view cells[:,cell])
end

function _checked_orientation_range(element_range::UnitRange{Int},
                                    element_count::Int,
                                    caller::AbstractString)
    first_element=first(element_range)
    last_element=last(element_range)
    (first_element>=1 && last_element<=element_count &&
     first_element<=last_element+1) || throw(ArgumentError(
        "$caller: element range $element_range is outside 1:$element_count"))
    return first_element,last_element
end

function _orientations_in_range!(result::Vector{Int32},
                                 cells::AbstractMatrix{Int32},
                                 space::_FunctionSpace,
                                 first_element::Int,last_element::Int)
    slot=0
    @inbounds for cell in first_element:last_element
        slot+=1
        result[slot]=_cell_orientation(cells,cell,space)
    end
    return nothing
end

"""Return one orientation index per cached element of a supported type."""
function mesh_basis_orientations(mesh::Mesh,element_type_value,
                                 function_space_type;
                                 caller::AbstractString=
                                     "mesh_basis_orientations")
    element_type=_checked_element_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    element_type,_,_,_=_basis_element_contract(element_type,space,caller)
    block=mesh_element_block(mesh,element_type)
    block===nothing && return Int32[]
    _,cells=block
    result=Vector{Int32}(undef,size(cells,2))
    _orientations_in_range!(result,cells,space,1,size(cells,2))
    return result
end

"""Return one orientation index per cached element in a 1-based block range.

This is the contiguous-block partition backing the `task`/`num_tasks` session
query contract: callers select `begin=(task*count)÷num_tasks` through
`end=((task+1)*count)÷num_tasks`. Validation and orientation contracts match
[`mesh_basis_orientations`](@ref); only the selected positions are evaluated,
so an empty range returns an empty vector without touching the mesh.
"""
function mesh_basis_orientations(mesh::Mesh,element_type_value,
                                 function_space_type,
                                 element_range::UnitRange{Int};
                                 caller::AbstractString=
                                     "mesh_basis_orientations")
    element_type=_checked_element_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    element_type,_,_,_=_basis_element_contract(element_type,space,caller)
    block=mesh_element_block(mesh,element_type)
    block===nothing && return Int32[]
    _,cells=block
    element_count=size(cells,2)
    first_element,last_element=
        _checked_orientation_range(element_range,element_count,caller)
    selected=max(last_element-first_element+1,0)
    result=Vector{Int32}(undef,selected)
    _orientations_in_range!(result,cells,space,first_element,last_element)
    return result
end

"""Return the orientation index for one dense cached element tag."""
function mesh_basis_orientation(mesh::Mesh,element_tag_value,
                                function_space_type;
                                caller::AbstractString=
                                    "mesh_basis_orientation")
    tag=_checked_element_tag(mesh,element_tag_value,caller)
    space=_function_space(function_space_type,caller)
    record=mesh_element_record(mesh,tag)
    _basis_element_contract(record.element_type,space,caller)
    return space.hierarchical ? _orientation_rank(record.node_tags) : Int32(0)
end

"""Return the node (0) or edge (1) dimension that owns keys for a space."""
function mesh_key_dimension(element_type_value,function_space_type;
                            caller::AbstractString="mesh_key_dimension")
    element_type=_checked_element_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    _basis_element_contract(element_type,space,caller)
    return space.key_dimension
end

"""Return the number of keys per supported reference element."""
function mesh_number_of_keys(element_type_value,function_space_type;
                             caller::AbstractString="mesh_number_of_keys")
    element_type=_checked_element_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    element_type,_,nodal_count,_=
        _basis_element_contract(element_type,space,caller)
    return Int32(space.hierarchical ?
                 (space.key_dimension==0 ? _vertex_count(element_type) :
                                           _edge_count(element_type)) :
                 nodal_count)
end

function _checked_bool(value,caller::AbstractString,name::AbstractString)
    value isa Bool || throw(ArgumentError("$caller: $name must be Bool"))
    return value
end

function _coordinate_result_length(key_count::Int,return_coord::Bool,
                                   caller::AbstractString)
    return_coord || return 0
    return _checked_result_length(3,key_count;caller=caller)
end

function _node_keys(mesh::Mesh,cells::AbstractMatrix{Int32},
                    return_coord::Bool,caller::AbstractString)
    key_count=length(cells)
    type_keys=fill(Int32(0),key_count)
    entity_keys=Vector{UInt64}(undef,key_count)
    coordinates=Vector{Float64}(
        undef,_coordinate_result_length(key_count,return_coord,caller))
    cursor=0
    @inbounds for cell in axes(cells,2),local_node in axes(cells,1)
        cursor+=1
        node=Int(cells[local_node,cell])
        entity_keys[cursor]=UInt64(node)
        if return_coord
            offset=3cursor-3
            coordinates[offset+1]=mesh.coords[1,node]
            coordinates[offset+2]=mesh.coords[2,node]
            coordinates[offset+3]=mesh.coords[3,node]
        end
    end
    return type_keys,entity_keys,coordinates
end

function _edge_keys(mesh::Mesh,cells::AbstractMatrix{Int32},
                    element_type::Int,topology::Union{Nothing,MeshEdgeTopology},
                    return_coord::Bool,caller::AbstractString)
    topology===nothing && throw(ArgumentError(
        "$caller: edge keys require a populated mesh edge catalog"))
    patterns=_simplex_edge_patterns(element_type)
    key_count=_checked_result_length(
        size(cells,2),length(patterns);caller=caller)
    node_pairs=Vector{UInt64}(undef,2key_count)
    cursor=0
    @inbounds for cell in axes(cells,2),pattern in patterns
        cursor+=1
        node_pairs[2cursor-1]=UInt64(cells[pattern[1],cell])
        node_pairs[2cursor]=UInt64(cells[pattern[2],cell])
    end
    entity_keys,_=_mesh_edges(topology,mesh,node_pairs,caller)
    type_keys=fill(Int32(1),key_count)
    coordinates=Vector{Float64}(
        undef,_coordinate_result_length(key_count,return_coord,caller))
    if return_coord
        @inbounds for key in 1:key_count
            first_node=Int(node_pairs[2key-1])
            second_node=Int(node_pairs[2key])
            offset=3key-3
            for component in 1:3
                coordinates[offset+component]=_stable_midpoint(
                    mesh.coords[component,first_node],
                    mesh.coords[component,second_node])
            end
        end
    end
    return type_keys,entity_keys,coordinates
end

function _keys_for_cells(mesh::Mesh,cells::AbstractMatrix{Int32},
                         element_type::Int,space::_FunctionSpace,
                         nodal_count::Int,
                         topology::Union{Nothing,MeshEdgeTopology},
                         return_coord::Bool,caller::AbstractString)
    if !space.hierarchical && size(cells,1)!=nodal_count
        throw(ArgumentError(
            "$caller: element type $element_type stores $(size(cells,1)) " *
            "nodal keys, but the requested basis requires $nodal_count; " *
            "the cached mesh owns only its stored interpolation nodes"))
    end
    return space.key_dimension==0 ?
        _node_keys(mesh,cells,return_coord,caller) :
        _edge_keys(
            mesh,cells,element_type,topology,return_coord,caller)
end

"""Return detached keys for all cached elements of one supported type."""
function mesh_keys(mesh::Mesh,element_type_value,function_space_type,
                   topology::Union{Nothing,MeshEdgeTopology}=nothing;
                   return_coord=true,caller::AbstractString="mesh_keys")
    element_type=_checked_element_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    element_type,_,nodal_count,_=
        _basis_element_contract(element_type,space,caller)
    coordinates_requested=_checked_bool(
        return_coord,caller,"return_coord")
    block=mesh_element_block(mesh,element_type)
    block===nothing && return Int32[],UInt64[],Float64[]
    _,cells=block
    return _keys_for_cells(
        mesh,cells,element_type,space,nodal_count,topology,
        coordinates_requested,caller)
end

"""Return detached keys for one dense cached element tag."""
function mesh_keys_for_element(
    mesh::Mesh,element_tag_value,function_space_type,
    topology::Union{Nothing,MeshEdgeTopology}=nothing;
    return_coord=true,caller::AbstractString="mesh_keys_for_element")
    tag=_checked_element_tag(mesh,element_tag_value,caller)
    space=_function_space(function_space_type,caller)
    coordinates_requested=_checked_bool(
        return_coord,caller,"return_coord")
    record=mesh_element_record(mesh,tag)
    element_type,_,nodal_count,_=
        _basis_element_contract(record.element_type,space,caller)
    cells=reshape(Int32.(record.node_tags),length(record.node_tags),1)
    return _keys_for_cells(
        mesh,cells,element_type,space,nodal_count,topology,
        coordinates_requested,caller)
end

function _checked_type_keys(values,expected::Int,caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: type_keys must be a vector or tuple of integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    result=Vector{Int32}(undef,length(values))
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: type_keys[$index] must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: type_keys[$index] must not be Bool"))
        converted=try Int(value) catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: type_keys[$index] exceeds Int bounds"))
        end
        converted==expected || throw(ArgumentError(
            "$caller: type_keys[$index] must be $expected for this function space"))
        result[index]=Int32(converted)
    end
    return result
end

function _checked_entity_keys(values,caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: entity_keys must be a vector or tuple of positive integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    result=Vector{UInt64}(undef,length(values))
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: entity_keys[$index] must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: entity_keys[$index] must not be Bool"))
        converted=try UInt64(value) catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: entity_keys[$index] must fit UInt64"))
        end
        converted>0 || throw(ArgumentError(
            "$caller: entity_keys[$index] must be positive"))
        result[index]=converted
    end
    return result
end

"""Return `(entity dimension, polynomial order)` metadata for complete key groups."""
function mesh_keys_information(type_keys,entity_keys,element_type_value,
                               function_space_type;
                               caller::AbstractString="mesh_keys_information")
    element_type=_checked_element_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    element_type,_,nodal_count,basis_type=
        _basis_element_contract(element_type,space,caller)
    expected_type=space.key_dimension
    types=_checked_type_keys(type_keys,expected_type,caller)
    entities=_checked_entity_keys(entity_keys,caller)
    length(types)==length(entities) || throw(ArgumentError(
        "$caller: type_keys and entity_keys must have equal lengths"))
    keys_per_element=space.hierarchical ?
        (space.key_dimension==0 ? _vertex_count(element_type) :
                                  _edge_count(element_type)) : nodal_count
    length(types)%keys_per_element==0 || throw(ArgumentError(
        "$caller: key count $(length(types)) must be divisible by " *
        "$keys_per_element for element type $element_type"))
    if space.hierarchical
        return fill(
            (Int32(space.key_dimension),Int32(space.key_order)),length(types))
    end
    bubble_count=_nodal_bubble_count(basis_type)
    nonbubble_count=nodal_count-bubble_count
    dimension=Int32(msh_spec(basis_type).dim)
    result=Vector{Tuple{Int32,Int32}}(undef,length(types))
    order=Int32(space.key_order)
    for group_start in 1:nodal_count:length(result)
        fill!(@view(result[group_start:group_start+nonbubble_count-1]),
              (Int32(0),order))
        if bubble_count>0
            first_bubble=group_start+nonbubble_count
            last_bubble=group_start+nodal_count-1
            fill!(@view(result[first_bubble:last_bubble]),(dimension,order))
        end
    end
    return result
end

end # module MeshFunctionSpaces
