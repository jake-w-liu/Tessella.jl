"""
    MeshQuadrature

Session-independent reference quadrature for fixed-node Point, Line, Triangle,
and Tetrahedron elements. Results use Gmsh's flattened `(u,v,w)` layout. Low-
degree `Gauss` rules preserve Gmsh 4.15.2's economical simplex tables, while
bounded `CompositeGauss` rules are generated natively from Gauss--Legendre
points and Duffy maps.
"""
module MeshQuadrature

using ..Elements: msh_spec

export mesh_integration_points

const _MAX_GAUSS_LEGENDRE_POINTS=128
const _MAX_QUADRATURE_POINTS=1_000_000
const _MAX_ECONOMICAL_SIMPLEX_ORDER=5

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
    spec.family in (:pnt,:lin,:tri,:tet) || throw(ArgumentError(
        "$caller: element type $element_type belongs to unsupported " *
        "$(spec.family) family; supported fixed-node families are Point, " *
        "Line, Triangle, and Tetrahedron"))
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
    return composite,order,name
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
function _triangle_gauss(order::Int)
    order<=1 && return _flatten_rule(
        ((0.333333333333333,0.333333333333333,0.0),),
        (0.500000000000000,))
    order==2 && return _flatten_rule(
        ((0.166666666666667,0.166666666666667,0.0),
         (0.166666666666667,0.666666666666667,0.0),
         (0.666666666666667,0.166666666666667,0.0)),
        (0.166666666666667,0.166666666666667,0.166666666666667))
    order==3 && return _flatten_rule(
        ((0.333333333333333,0.333333333333333,0.0),
         (0.200000000000000,0.200000000000000,0.0),
         (0.200000000000000,0.600000000000000,0.0),
         (0.600000000000000,0.200000000000000,0.0)),
        (-0.281250000000000,0.260416666666667,
         0.260416666666667,0.260416666666667))
    order==4 && return _flatten_rule(
        ((0.445948490915965,0.445948490915965,0.0),
         (0.445948490915965,0.108103018168070,0.0),
         (0.108103018168070,0.445948490915965,0.0),
         (0.091576213509771,0.091576213509771,0.0),
         (0.091576213509771,0.816847572980459,0.0),
         (0.816847572980459,0.091576213509771,0.0)),
        (0.111690794839005,0.111690794839005,0.111690794839005,
         0.054975871827661,0.054975871827661,0.054975871827661))
    return _flatten_rule(
        ((0.333333333333333,0.333333333333333,0.0),
         (0.470142064105115,0.470142064105115,0.0),
         (0.470142064105115,0.059715871789770,0.0),
         (0.059715871789770,0.470142064105115,0.0),
         (0.101286507323456,0.101286507323456,0.0),
         (0.101286507323456,0.797426985353087,0.0),
         (0.797426985353087,0.101286507323456,0.0)),
        (0.112500000000000,
         0.066197076394253,0.066197076394253,0.066197076394253,
         0.062969590272414,0.062969590272414,0.062969590272414))
end

function _tetrahedron_gauss(order::Int)
    order<=1 && return _flatten_rule(
        ((0.25,0.25,0.25),),(0.166666666666667,))
    order==2 && return _flatten_rule(
        ((0.138196601125,0.138196601125,0.138196601125),
         (0.585410196625,0.138196601125,0.138196601125),
         (0.138196601125,0.585410196625,0.138196601125),
         (0.138196601125,0.138196601125,0.585410196625)),
        (0.0416666666666667,0.0416666666666667,
         0.0416666666666667,0.0416666666666667))
    order==3 && return _flatten_rule(
        ((0.25,0.25,0.25),
         (0.166666666667,0.166666666667,0.166666666667),
         (0.166666666667,0.166666666667,0.500000000000),
         (0.166666666667,0.500000000000,0.166666666667),
         (0.500000000000,0.166666666667,0.166666666667)),
        (-0.133333333333333,0.075000000000000,0.075000000000000,
         0.075000000000000,0.075000000000000))
    order==4 && return _flatten_rule(
        ((0.2500000000000,0.2500000000000,0.2500000000000),
         (0.0714285714286,0.0714285714286,0.0714285714286),
         (0.0714285714286,0.0714285714286,0.7857142857140),
         (0.0714285714286,0.7857142857140,0.0714285714286),
         (0.7857142857140,0.0714285714286,0.0714285714286),
         (0.3994035761670,0.3994035761670,0.1005964238330),
         (0.3994035761670,0.1005964238330,0.3994035761670),
         (0.1005964238330,0.3994035761670,0.3994035761670),
         (0.3994035761670,0.1005964238330,0.1005964238330),
         (0.1005964238330,0.3994035761670,0.1005964238330),
         (0.1005964238330,0.1005964238330,0.3994035761670)),
        (-0.0131555555555,
         0.0076222222222,0.0076222222222,0.0076222222222,
         0.0076222222222,
         0.0248888888888,0.0248888888888,0.0248888888888,
         0.0248888888888,0.0248888888888,0.0248888888888))
    return _flatten_rule(
        ((0.0927352503109,0.0927352503109,0.0927352503109),
         (0.7217942490670,0.0927352503109,0.0927352503109),
         (0.0927352503109,0.7217942490670,0.0927352503109),
         (0.0927352503109,0.0927352503109,0.7217942490670),
         (0.3108859192630,0.3108859192630,0.3108859192630),
         (0.0673422422101,0.3108859192630,0.3108859192630),
         (0.3108859192630,0.0673422422101,0.3108859192630),
         (0.3108859192630,0.3108859192630,0.0673422422101),
         (0.4544962958740,0.4544962958740,0.0455037041256),
         (0.4544962958740,0.0455037041256,0.4544962958740),
         (0.0455037041256,0.4544962958740,0.4544962958740),
         (0.4544962958740,0.0455037041256,0.0455037041256),
         (0.0455037041256,0.4544962958740,0.0455037041256),
         (0.0455037041256,0.0455037041256,0.4544962958740)),
        (0.01224884051940,0.01224884051940,
         0.01224884051940,0.01224884051940,
         0.01878132095300,0.01878132095300,
         0.01878132095300,0.01878132095300,
         0.00709100346285,0.00709100346285,0.00709100346285,
         0.00709100346285,0.00709100346285,0.00709100346285))
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

"""
    mesh_integration_points(element_type, integration_type)

Return detached `(local_coordinates, weights)` for a fixed-node Point, Line,
Triangle, or Tetrahedron type. `integration_type` is `GaussN` or
`CompositeGaussN`, with an omitted `N` meaning zero. Economical triangle and
tetrahedron `Gauss` tables cover orders zero through five; use
`CompositeGaussN` for a higher bounded order.
"""
function mesh_integration_points(element_type,integration_type;
                                 caller::AbstractString=
                                     "mesh_integration_points")
    spec=_quadrature_element_spec(element_type,caller)
    composite,order,name=_integration_rule(integration_type,caller)
    family=spec.family
    family===:pnt && return Float64[0.0,0.0,0.0],Float64[1.0]
    family===:lin && return _line_rule(order,caller)
    if !composite
        order<=_MAX_ECONOMICAL_SIMPLEX_ORDER || throw(ArgumentError(
            "$caller: economical $(repr(name)) is not implemented for " *
            "$(family===:tri ? "Triangle" : "Tetrahedron") elements; " *
            "use CompositeGauss$order"))
        return family===:tri ? _triangle_gauss(order) :
               _tetrahedron_gauss(order)
    end
    return family===:tri ? _triangle_composite(order,caller) :
           _tetrahedron_composite(order,caller)
end

end # module MeshQuadrature
