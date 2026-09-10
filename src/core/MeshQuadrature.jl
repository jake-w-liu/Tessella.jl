"""
    MeshQuadrature

Session-independent reference quadrature for every fixed-node Gmsh family with
defined rules: Point, Line, Triangle, Quadrangle, Tetrahedron, Hexahedron,
Prism, and Pyramid. Results use Gmsh's flattened `(u,v,w)` layout. Low-degree
`Gauss` rules preserve Gmsh 4.15.2's economical tables and tensor-rule
transitions, while bounded `CompositeGauss` rules are generated natively from
Gauss--Legendre points, Duffy maps, and the pyramid's Gauss--Jacobi rule.
"""
module MeshQuadrature

using LinearAlgebra: SymTridiagonal, eigen
using ..Elements: msh_spec

export mesh_integration_points

const _MAX_GAUSS_LEGENDRE_POINTS=128
const _MAX_QUADRATURE_POINTS=1_000_000
const _MAX_TRIANGLE_ECONOMICAL_ORDER=20
const _MAX_TETRAHEDRON_ECONOMICAL_ORDER=21

function _quadrature_element_spec(value,caller::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$caller: element_type must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: element_type must not be Bool"))
    element_type=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError ||
         err isa MethodError) || rethrow()
        throw(ArgumentError(
            "$caller: element_type exceeds the platform Int range"))
    end
    spec=try
        msh_spec(element_type)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        throw(ArgumentError(
            "$caller: element_type $element_type is not an ordinary " *
            "fixed-node Gmsh element"))
    end
    spec.family===:trih && throw(ArgumentError(
        "$caller: Trihedron element type $element_type has no integration " *
        "rules in Gmsh 4.15.2"))
    spec.family in (:pnt,:lin,:tri,:qua,:tet,:hex,:pri,:pyr) ||
        throw(ArgumentError(
            "$caller: element type $element_type belongs to unsupported " *
            "$(spec.family) family"))
    return spec
end

function _integration_rule(value,caller::AbstractString)
    value isa AbstractString || throw(ArgumentError(
        "$caller: integration_type must be a string"))
    name=String(value)
    bytes=codeunits(name)
    any(==(0x00),bytes) && throw(ArgumentError(
        "$caller: integration_type must not contain NUL"))

    composite,prefix_length=if startswith(name,"CompositeGauss")
        true,ncodeunits("CompositeGauss")
    elseif startswith(name,"Gauss")
        false,ncodeunits("Gauss")
    else
        throw(ArgumentError(
            "$caller: unknown integration_type $(repr(name)); expected " *
            "GaussN or CompositeGaussN"))
    end

    order=0
    for index in (prefix_length+1):length(bytes)
        byte=bytes[index]
        0x30<=byte<=0x39 || throw(ArgumentError(
            "$caller: integration_type $(repr(name)) must end in an " *
            "optional nonnegative ASCII decimal order"))
        digit=Int(byte-0x30)
        order<=(typemax(Int)-digit)÷10 || throw(ArgumentError(
            "$caller: integration order exceeds the platform Int range"))
        order=10order+digit
    end
    return composite,order
end

@inline function _legendre_pair(degree::Int,x::Float64)
    degree==0 && return 1.0,0.0
    previous=1.0
    current=x
    for index in 2:degree
        following=((2index-1)*x*current-(index-1)*previous)/index
        previous,current=current,following
    end
    return current,previous
end

function _gauss_legendre(point_count::Int)
    1<=point_count<=_MAX_GAUSS_LEGENDRE_POINTS || error(
        "MeshQuadrature: internal Gauss--Legendre point count is out of bounds")
    point_count==1 && return Float64[0.0],Float64[2.0]

    points=Vector{Float64}(undef,point_count)
    weights=similar(points)
    half=point_count÷2
    for index in 1:half
        root=cospi((index-0.25)/(point_count+0.5))
        converged=false
        for _ in 1:64
            polynomial,previous=_legendre_pair(point_count,root)
            derivative=point_count*(root*polynomial-previous)/(root^2-1)
            next_root=root-polynomial/derivative
            isfinite(next_root) || error(
                "MeshQuadrature: Gauss--Legendre iteration became nonfinite")
            if abs(next_root-root)<=4eps(Float64)*max(1.0,abs(next_root))
                root=next_root
                converged=true
                break
            end
            root=next_root
        end
        converged || error(
            "MeshQuadrature: Gauss--Legendre iteration did not converge")
        polynomial,previous=_legendre_pair(point_count,root)
        derivative=point_count*(root*polynomial-previous)/(root^2-1)
        weight=2/((1-root^2)*derivative^2)
        (isfinite(weight) && weight>0) || error(
            "MeshQuadrature: Gauss--Legendre weight is not finite and positive")
        points[index]=-root
        points[point_count+1-index]=root
        weights[index]=weight
        weights[point_count+1-index]=weight
    end
    if isodd(point_count)
        middle=half+1
        polynomial,previous=_legendre_pair(point_count,0.0)
        derivative=-point_count*previous
        weight=2/derivative^2
        (polynomial==0.0 && isfinite(weight) && weight>0) || error(
            "MeshQuadrature: invalid central Gauss--Legendre point")
        points[middle]=0.0
        weights[middle]=weight
    end
    return points,weights
end

function _axis_count(family::Symbol,order::Int,caller::AbstractString)
    order<=2*_MAX_GAUSS_LEGENDRE_POINTS-1 || throw(ArgumentError(
        "$caller: integration order $order exceeds the bounded " *
        "Gauss--Legendre limit"))
    point_count=family===:lin ? order÷2+1 :
                family===:tri ? (order+3)÷2 : (order+4)÷2
    point_count<=_MAX_GAUSS_LEGENDRE_POINTS || throw(ArgumentError(
        "$caller: integration order $order requires $point_count " *
        "Gauss--Legendre points per axis; the limit is " *
        "$_MAX_GAUSS_LEGENDRE_POINTS"))
    total=family===:lin ? point_count :
          family===:tri ? point_count^2 : point_count^3
    total<=_MAX_QUADRATURE_POINTS || throw(ArgumentError(
        "$caller: integration order $order requires $total quadrature " *
        "points; the limit is $_MAX_QUADRATURE_POINTS"))
    return point_count,total
end

function _tensor_axis_count(family::Symbol,order::Int,caller::AbstractString)
    order<=2*_MAX_GAUSS_LEGENDRE_POINTS-1 || throw(ArgumentError(
        "$caller: integration order $order exceeds the bounded " *
        "one-dimensional quadrature limit"))
    point_count=family===:pri ? (order+3)÷2 : order÷2+1
    point_count<=_MAX_GAUSS_LEGENDRE_POINTS || throw(ArgumentError(
        "$caller: integration order $order requires $point_count " *
        "quadrature points per axis; the limit is " *
        "$_MAX_GAUSS_LEGENDRE_POINTS"))
    dimension=family===:qua ? 2 : 3
    total=point_count
    for _ in 2:dimension
        total=Base.checked_mul(total,point_count)
    end
    total<=_MAX_QUADRATURE_POINTS || throw(ArgumentError(
        "$caller: integration order $order requires $total quadrature " *
        "points; the limit is $_MAX_QUADRATURE_POINTS"))
    return point_count,total
end

function _flatten_rule(points,weights)
    point_count=length(points)
    point_count==length(weights) || error(
        "MeshQuadrature: inconsistent economical quadrature table")
    coordinates=Vector{Float64}(undef,3point_count)
    result_weights=Vector{Float64}(undef,point_count)
    for index in 1:point_count
        point=points[index]
        offset=3(index-1)
        coordinates[offset+1]=point[1]
        coordinates[offset+2]=point[2]
        coordinates[offset+3]=point[3]
        result_weights[index]=weights[index]
    end
    return coordinates,result_weights
end

# Authority: the economical Solin tables used by pinned Gmsh 4.15.2.
include("MeshQuadratureSimplex.jl")

# Authority: the low-order economical tables in Gmsh 4.15.2's
# GaussQuadratureQuad.cpp and GaussQuadratureHex.cpp. Higher orders use the
# canonical tensor constructors below, as Gmsh does.
function _quadrangle_gauss(order::Int)
    order==1 && return _flatten_rule(
        ((0.816496580928,0.0,0.0),
         (-0.408248290464,0.840896415255,0.0),
         (-0.408248290464,-0.840896415255,0.0)),
        (1.3333333333333,1.3333333333333,1.3333333333333))
    order==2 || error(
        "MeshQuadrature: internal economical quadrangle order is invalid")
    return _flatten_rule(
        ((0.0,0.0,0.0),
         (0.0,0.9660917830792959,0.0),
         (0.0,-0.9660917830792959,0.0),
         (0.7745966692414834,0.7745966692414834,0.0),
         (0.7745966692414834,-0.7745966692414834,0.0),
         (-0.7745966692414834,0.7745966692414834,0.0),
         (-0.7745966692414834,-0.7745966692414834,0.0)),
        (1.1428571428571428,
         0.31746031746031744,0.31746031746031744,
         0.5555555555555556,0.5555555555555556,
         0.5555555555555556,0.5555555555555556))
end

function _hexahedron_gauss_one()
    return _flatten_rule(
        ((0.40824826,0.70710678,-0.57735027),
         (0.40824826,-0.70710678,-0.57735027),
         (-0.40824826,0.70710678,0.57735027),
         (-0.40824826,-0.70710678,0.57735027),
         (-0.81649658,0.0,-0.57735027),
         (0.81649658,0.0,0.57735027)),
        (1.3333333333,1.3333333333,1.3333333333,
         1.3333333333,1.3333333333,1.3333333333))
end

function _line_rule(order::Int,caller::AbstractString)
    point_count,_=_axis_count(:lin,order,caller)
    points,weights=_gauss_legendre(point_count)
    coordinates=zeros(Float64,3point_count)
    for index in 1:point_count
        coordinates[3index-2]=points[index]
    end
    return coordinates,weights
end

function _triangle_composite(order::Int,caller::AbstractString)
    axis_count,total=_axis_count(:tri,order,caller)
    points,axis_weights=_gauss_legendre(axis_count)
    coordinates=Vector{Float64}(undef,3total)
    weights=Vector{Float64}(undef,total)
    index=0
    for first in 1:axis_count,second in 1:axis_count
        u=0.5*(1.0+points[first])
        remainder=1.0-u
        v=0.5*(1.0+points[second])*remainder
        jacobian=0.25*remainder
        index+=1
        offset=3(index-1)
        coordinates[offset+1]=u
        coordinates[offset+2]=v
        coordinates[offset+3]=0.0
        weights[index]=jacobian*axis_weights[first]*axis_weights[second]
    end
    return coordinates,weights
end

function _tetrahedron_composite(order::Int,caller::AbstractString)
    axis_count,total=_axis_count(:tet,order,caller)
    points,axis_weights=_gauss_legendre(axis_count)
    coordinates=Vector{Float64}(undef,3total)
    weights=Vector{Float64}(undef,total)
    index=0
    for first in 1:axis_count,second in 1:axis_count,third in 1:axis_count
        u=0.5*(1.0+points[first])
        first_remainder=1.0-u
        v=0.5*(1.0+points[second])*first_remainder
        second_remainder=1.0-u-v
        w=0.5*(1.0+points[third])*second_remainder
        jacobian=0.125*first_remainder*second_remainder
        index+=1
        offset=3(index-1)
        coordinates[offset+1]=u
        coordinates[offset+2]=v
        coordinates[offset+3]=w
        weights[index]=jacobian*axis_weights[first]*axis_weights[second]*
                       axis_weights[third]
    end
    return coordinates,weights
end

function _triangle_gauss(order::Int,caller::AbstractString)
    order<=_MAX_TRIANGLE_ECONOMICAL_ORDER &&
        return _triangle_economical_rule(order)
    return _triangle_composite(order,caller)
end

function _tetrahedron_gauss(order::Int,caller::AbstractString)
    order<=9 && return _tetrahedron_economical_rule(order)
    order<=_MAX_TETRAHEDRON_ECONOMICAL_ORDER &&
        return _tetrahedron_lattice_rule(order)
    return _tetrahedron_composite(order,caller)
end

function _cartesian_tensor_rule(family::Symbol,order::Int,
                                caller::AbstractString)
    axis_count,total=_tensor_axis_count(family,order,caller)
    points,axis_weights=_gauss_legendre(axis_count)
    coordinates=Vector{Float64}(undef,3total)
    weights=Vector{Float64}(undef,total)
    index=0
    if family===:qua
        for first in 1:axis_count,second in 1:axis_count
            index+=1
            offset=3(index-1)
            coordinates[offset+1]=points[first]
            coordinates[offset+2]=points[second]
            coordinates[offset+3]=0.0
            weights[index]=axis_weights[first]*axis_weights[second]
        end
    else
        for first in 1:axis_count,second in 1:axis_count,
            third in 1:axis_count
            index+=1
            offset=3(index-1)
            coordinates[offset+1]=points[first]
            coordinates[offset+2]=points[second]
            coordinates[offset+3]=points[third]
            weights[index]=axis_weights[first]*axis_weights[second]*
                           axis_weights[third]
        end
    end
    return coordinates,weights
end

function _prism_rule(order::Int,composite::Bool,caller::AbstractString)
    axis_count,composite_total=_tensor_axis_count(:pri,order,caller)
    triangle_coordinates,triangle_weights=composite ?
        _triangle_composite(order,caller) : _triangle_gauss(order,caller)
    total=Base.checked_mul(length(triangle_weights),axis_count)
    total<=_MAX_QUADRATURE_POINTS || throw(ArgumentError(
        "$caller: integration order $order requires $total quadrature " *
        "points; the limit is $_MAX_QUADRATURE_POINTS"))
    composite && total!=composite_total && error(
        "MeshQuadrature: inconsistent composite prism point count")
    line_points,line_weights=_gauss_legendre(axis_count)
    coordinates=Vector{Float64}(undef,3total)
    weights=Vector{Float64}(undef,total)
    index=0
    for triangle_point in eachindex(triangle_weights),
        line_point in 1:axis_count
        index+=1
        source_offset=3(triangle_point-1)
        target_offset=3(index-1)
        coordinates[target_offset+1]=triangle_coordinates[source_offset+1]
        coordinates[target_offset+2]=triangle_coordinates[source_offset+2]
        coordinates[target_offset+3]=line_points[line_point]
        weights[index]=triangle_weights[triangle_point]*
                       line_weights[line_point]
    end
    return coordinates,weights
end

function _gauss_jacobi_20(point_count::Int)
    # Golub--Welsch for weight (1-x)^2 on [-1,1]. The simplified recurrence
    # coefficients avoid gamma functions and remain well-scaled at the public
    # point-count bound.
    diagonal=Vector{Float64}(undef,point_count)
    for degree in 0:(point_count-1)
        diagonal[degree+1]=-1/((degree+1)*(degree+2))
    end
    off_diagonal=Vector{Float64}(undef,max(0,point_count-1))
    for degree in 1:(point_count-1)
        off_diagonal[degree]=degree*(degree+2)/(
            (degree+1)*sqrt((2degree+1)*(2degree+3)))
    end
    decomposition=eigen(SymTridiagonal(diagonal,off_diagonal))
    points=decomposition.values
    weights=Vector{Float64}(undef,point_count)
    for index in 1:point_count
        weights[index]=(8/3)*abs2(decomposition.vectors[1,index])
        (-1<points[index]<1 && isfinite(weights[index]) && weights[index]>0) ||
            error("MeshQuadrature: Gauss--Jacobi eigensolve produced an " *
                  "invalid point or weight")
    end
    return points,weights
end

function _pyramid_rule(order::Int,caller::AbstractString)
    axis_count,total=_tensor_axis_count(:pyr,order,caller)
    line_points,line_weights=_gauss_legendre(axis_count)
    jacobi_points,jacobi_weights=_gauss_jacobi_20(axis_count)
    coordinates=Vector{Float64}(undef,3total)
    weights=Vector{Float64}(undef,total)
    index=0
    for vertical in 1:axis_count,first in 1:axis_count,
        second in 1:axis_count
        scale=0.5*(1-jacobi_points[vertical])
        index+=1
        offset=3(index-1)
        coordinates[offset+1]=scale*line_points[first]
        coordinates[offset+2]=scale*line_points[second]
        coordinates[offset+3]=0.5*(1+jacobi_points[vertical])
        weights[index]=0.125*line_weights[first]*line_weights[second]*
                       jacobi_weights[vertical]
    end
    return coordinates,weights
end

"""
    mesh_integration_points(element_type, integration_type)

Return detached `(local_coordinates, weights)` for any fixed-node family with a
Gmsh 4.15.2 integration rule. `integration_type` is `GaussN` or
`CompositeGaussN`, with an omitted `N` meaning zero. Triangle `Gauss` rules use
economical tables through order 20, while Tetrahedron rules use economical
tables through order 21; higher orders and `CompositeGaussN` use bounded tensor
constructions. Prism rules combine the corresponding Triangle and Line rules.
Trihedra have no rule in the pinned Gmsh release.
"""
function mesh_integration_points(element_type,integration_type;
                                 caller::AbstractString=
                                     "mesh_integration_points")
    spec=_quadrature_element_spec(element_type,caller)
    composite,order=_integration_rule(integration_type,caller)
    family=spec.family
    family===:pnt && return Float64[0.0,0.0,0.0],Float64[1.0]
    family===:lin && return _line_rule(order,caller)
    if family in (:tri,:tet)
        return family===:tri ?
               (composite ? _triangle_composite(order,caller) :
                            _triangle_gauss(order,caller)) :
               (composite ? _tetrahedron_composite(order,caller) :
                            _tetrahedron_gauss(order,caller))
    elseif family===:qua
        !composite && order in (1,2) && return _quadrangle_gauss(order)
        return _cartesian_tensor_rule(:qua,order,caller)
    elseif family===:hex
        !composite && order==1 && return _hexahedron_gauss_one()
        return _cartesian_tensor_rule(:hex,order,caller)
    elseif family===:pri
        return _prism_rule(order,composite,caller)
    else
        return _pyramid_rule(order,caller)
    end
end

end # module MeshQuadrature
