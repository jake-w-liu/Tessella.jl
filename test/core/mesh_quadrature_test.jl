using Test
using Tessella
using Tessella.MeshQuadrature
using Tessella.Elements: MSH_CATALOG, msh_spec

function _line_moment(power::Int)
    return isodd(power) ? 0.0 : 2.0/(power+1)
end

function _simplex_moment(exponents::Tuple)
    numerator=prod(factorial(big(exponent)) for exponent in exponents)
    denominator=factorial(big(sum(exponents)+length(exponents)))
    return Float64(numerator//denominator)
end

function _cartesian_moment(exponents::Tuple)
    return prod(_line_moment(exponent) for exponent in exponents)
end

function _pyramid_moment(exponents::NTuple{3,Int})
    first,second,vertical=exponents
    (isodd(first) || isodd(second)) && return 0.0
    numerator=4factorial(big(vertical))*factorial(big(first+second+2))
    denominator=(first+1)*(second+1)*
                factorial(big(first+second+vertical+3))
    return Float64(numerator//denominator)
end

function _quadrature_moment(coordinates,weights,exponents::Tuple)
    dimension=length(exponents)
    return sum(eachindex(weights);init=0.0) do index
        offset=3index-2
        value=weights[index]
        for axis in 1:dimension
            value*=coordinates[offset+axis-1]^exponents[axis]
        end
        value
    end
end

function _expect_quadrature_argument_error(f::Function)
    try
        f()
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return nothing
    end
    error("quadrature request unexpectedly succeeded")
end

function _check_rule_moments(element_type::Int,name::String,order::Int;
                             atol::Float64)
    coordinates,weights=mesh_integration_points(element_type,name)
    family=msh_spec(element_type).family
    if family===:lin
        for power in 0:order
            @test isapprox(
                _quadrature_moment(coordinates,weights,(power,)),
                _line_moment(power);atol=atol,rtol=atol)
        end
    elseif family===:tri
        for first in 0:order,second in 0:(order-first)
            exponents=(first,second)
            @test isapprox(
                _quadrature_moment(coordinates,weights,exponents),
                _simplex_moment(exponents);atol=atol,rtol=atol)
        end
    elseif family===:qua
        for first in 0:order,second in 0:(order-first)
            exponents=(first,second)
            @test isapprox(
                _quadrature_moment(coordinates,weights,exponents),
                _cartesian_moment(exponents);atol=atol,rtol=atol)
        end
    elseif family===:tet
        for first in 0:order,second in 0:(order-first),
            third in 0:(order-first-second)
            exponents=(first,second,third)
            @test isapprox(
                _quadrature_moment(coordinates,weights,exponents),
                _simplex_moment(exponents);atol=atol,rtol=atol)
        end
    elseif family===:hex
        for first in 0:order,second in 0:(order-first),
            third in 0:(order-first-second)
            exponents=(first,second,third)
            @test isapprox(
                _quadrature_moment(coordinates,weights,exponents),
                _cartesian_moment(exponents);atol=atol,rtol=atol)
        end
    elseif family===:pri
        for first in 0:order,second in 0:(order-first),
            third in 0:(order-first-second)
            exponents=(first,second,third)
            expected=_simplex_moment((first,second))*
                     _line_moment(third)
            @test isapprox(
                _quadrature_moment(coordinates,weights,exponents),
                expected;atol=atol,rtol=atol)
        end
    else
        @assert family===:pyr
        for first in 0:order,second in 0:(order-first),
            third in 0:(order-first-second)
            exponents=(first,second,third)
            @test isapprox(
                _quadrature_moment(coordinates,weights,exponents),
                _pyramid_moment(exponents);atol=atol,rtol=atol)
        end
    end
    return coordinates,weights
end

@testset "Gmsh-shaped reference quadrature" begin
    @test mesh_integration_points(15,"Gauss5")==
          (Float64[0,0,0],Float64[1])
    @test mesh_integration_points(1,"Gauss0")==
          (Float64[0,0,0],Float64[2])
    @test mesh_integration_points(2,"Gauss2")==(
        Float64[
            0.166666666666667,0.166666666666667,0,
            0.166666666666667,0.666666666666667,0,
            0.666666666666667,0.166666666666667,0],
        fill(0.166666666666667,3))
    @test mesh_integration_points(4,"Gauss3")==(
        Float64[
            0.25,0.25,0.25,
            0.166666666667,0.166666666667,0.166666666667,
            0.166666666667,0.166666666667,0.5,
            0.166666666667,0.5,0.166666666667,
            0.5,0.166666666667,0.166666666667],
        Float64[-0.133333333333333,0.075,0.075,0.075,0.075])
    @test mesh_integration_points(3,"Gauss1")==(
        Float64[
            0.816496580928,0,0,
            -0.408248290464,0.840896415255,0,
            -0.408248290464,-0.840896415255,0],
        fill(1.3333333333333,3))
    @test mesh_integration_points(3,"Gauss2")==(
        Float64[
            0,0,0,
            0,0.9660917830792959,0,
            0,-0.9660917830792959,0,
            0.7745966692414834,0.7745966692414834,0,
            0.7745966692414834,-0.7745966692414834,0,
            -0.7745966692414834,0.7745966692414834,0,
            -0.7745966692414834,-0.7745966692414834,0],
        Float64[1.1428571428571428,
                0.31746031746031744,0.31746031746031744,
                0.5555555555555556,0.5555555555555556,
                0.5555555555555556,0.5555555555555556])
    @test mesh_integration_points(5,"Gauss1")==(
        Float64[
            0.40824826,0.70710678,-0.57735027,
            0.40824826,-0.70710678,-0.57735027,
            -0.40824826,0.70710678,0.57735027,
            -0.40824826,-0.70710678,0.57735027,
            -0.81649658,0,-0.57735027,
            0.81649658,0,0.57735027],
        fill(1.3333333333,6))
    pyramid_zero=mesh_integration_points(7,"Gauss0")
    @test pyramid_zero[1]==Float64[0,0,0.25]
    @test isapprox(only(pyramid_zero[2]),4/3;atol=4eps(Float64),rtol=0)

    economical_counts=Dict(
        1=>[1,1,2,2,3,3],
        2=>[1,1,3,4,6,7],
        3=>[1,3,7,4,9,9],
        4=>[1,1,4,5,11,14],
        5=>[1,6,8,8,27,27],
        6=>[1,2,6,12,18,28],
        7=>[1,1,8,8,27,27])
    for element_type in (1,2,3,4,5,6,7),order in 0:5
        moment_tolerance=element_type==5 && order==1 ? 3e-10 : 8e-13
        # Gmsh's pinned economical quadrangle Gauss2 table is anisotropic:
        # it integrates constants and linears, but not the y^2 moment. Exact
        # table coverage above preserves that externally observable contract.
        checked_order=element_type==3 && order==2 ? 1 : order
        coordinates,weights=_check_rule_moments(
            element_type,"Gauss$order",checked_order;atol=moment_tolerance)
        @test length(weights)==economical_counts[element_type][order+1]
        @test length(coordinates)==3length(weights)
        @test all(isfinite,coordinates)
        @test all(isfinite,weights)
    end

    triangle_counts=(1,1,3,4,6,7,12,13,16,19,25,27,33,37,42,
                     48,52,61,70,73,79)
    tetrahedron_counts=(1,1,4,5,11,14,24,31,43,53,126,126,210,
                        210,330,330,495,495,715,715,1001,1001)
    for order in 6:20
        # The pinned decimal Gauss20 table has a 1.15e-10 analytic moment
        # residual; the exact-array differential separately guards every entry.
        triangle_tolerance=order==20 ? 2e-10 : 8e-13
        triangle=_check_rule_moments(
            2,"Gauss$order",order;atol=triangle_tolerance)
        prism=_check_rule_moments(
            6,"Gauss$order",order;atol=2triangle_tolerance)
        @test length(triangle[2])==triangle_counts[order+1]
        @test length(prism[2])==
              triangle_counts[order+1]*((order+3)÷2)
        @test all(isfinite,triangle[1])
        @test all(isfinite,triangle[2])
        @test all(isfinite,prism[1])
        @test all(isfinite,prism[2])
    end
    for order in 6:21
        # Gmsh's published lattice decimals accumulate up to an 8.8e-11
        # residual; exact ordering and values are covered by the differential.
        tolerance=order>=10 ? 2e-10 : 8e-13
        tetrahedron=_check_rule_moments(
            4,"Gauss$order",order;atol=tolerance)
        @test length(tetrahedron[2])==tetrahedron_counts[order+1]
        @test all(isfinite,tetrahedron[1])
        @test all(isfinite,tetrahedron[2])
    end
    for even_order in 10:2:20
        @test mesh_integration_points(4,"Gauss$even_order")==
              mesh_integration_points(4,"Gauss$(even_order+1)")
    end
    @test mesh_integration_points(2,"Gauss21")==
          mesh_integration_points(2,"CompositeGauss21")
    @test mesh_integration_points(4,"Gauss22")==
          mesh_integration_points(4,"CompositeGauss22")
    @test mesh_integration_points(6,"Gauss21")==
          mesh_integration_points(6,"CompositeGauss21")

    for order in (0,1,2,3,4,5,6,12,20)
        for element_type in (1,2,3,4,5,6,7)
            coordinates,weights=_check_rule_moments(
                element_type,"CompositeGauss$order",order;atol=8e-13)
            @test length(coordinates)==3length(weights)
            @test all(isfinite,coordinates)
            @test all(>(0),weights)
            if element_type==1
                @test all(value->-1<=value<=1,coordinates[1:3:end])
                @test all(iszero,coordinates[2:3:end])
                @test all(iszero,coordinates[3:3:end])
            elseif element_type==2
                u=coordinates[1:3:end]
                v=coordinates[2:3:end]
                @test all(u .>= 0)
                @test all(v .>= 0)
                @test all(index->u[index]+v[index]<=1,eachindex(u,v))
                @test all(iszero,coordinates[3:3:end])
            elseif element_type==4
                u=coordinates[1:3:end]
                v=coordinates[2:3:end]
                w=coordinates[3:3:end]
                @test all(u .>= 0)
                @test all(v .>= 0)
                @test all(w .>= 0)
                @test all(
                    index->u[index]+v[index]+w[index]<=1,
                    eachindex(u,v,w))
            elseif element_type==3
                @test all(value->-1<=value<=1,coordinates[1:3:end])
                @test all(value->-1<=value<=1,coordinates[2:3:end])
                @test all(iszero,coordinates[3:3:end])
            elseif element_type==5
                @test all(value->-1<=value<=1,coordinates)
            elseif element_type==6
                u=coordinates[1:3:end]
                v=coordinates[2:3:end]
                @test all(u .>= 0)
                @test all(v .>= 0)
                @test all(index->u[index]+v[index]<=1,eachindex(u,v))
                @test all(value->-1<=value<=1,coordinates[3:3:end])
            else
                u=coordinates[1:3:end]
                v=coordinates[2:3:end]
                w=coordinates[3:3:end]
                @test all(value->0<=value<=1,w)
                @test all(index->abs(u[index])<=1-w[index],eachindex(u,w))
                @test all(index->abs(v[index])<=1-w[index],eachindex(v,w))
            end
        end
    end

    @test mesh_integration_points(2,"Gauss")==
          mesh_integration_points(2,"Gauss0")
    @test mesh_integration_points(4,"CompositeGauss")==
          mesh_integration_points(4,"CompositeGauss0")
    @test mesh_integration_points(2,"Gauss0005")==
          mesh_integration_points(2,"Gauss5")
    substring=SubString("xxCompositeGauss05",3)
    @test mesh_integration_points(4,substring)==
          mesh_integration_points(4,"CompositeGauss5")

    family_reference=Dict(
        family=>mesh_integration_points(type,"CompositeGauss4")
        for (family,type) in ((:pnt,15),(:lin,1),(:tri,2),(:qua,3),
                              (:tet,4),(:hex,5),(:pri,6),(:pyr,7)))
    for (element_type,spec) in MSH_CATALOG
        spec.family in keys(family_reference) || continue
        @test mesh_integration_points(element_type,"CompositeGauss4")==
              family_reference[spec.family]
    end

    detached=mesh_integration_points(2,"Gauss5")
    original=deepcopy(detached)
    detached[1][1]=99
    detached[2][1]=99
    @test mesh_integration_points(2,"Gauss5")==original

    lattice_detached=mesh_integration_points(4,"Gauss21")
    lattice_original=deepcopy(lattice_detached)
    lattice_detached[1][1]=99
    lattice_detached[2][1]=99
    @test mesh_integration_points(4,"Gauss21")==lattice_original

    largest_line=mesh_integration_points(1,"CompositeGauss255")
    @test length(largest_line[2])==128
    @test issorted(largest_line[1][1:3:end])
    @test largest_line[1][1:3:end]==-reverse(largest_line[1][1:3:end])
    @test largest_line[2]==reverse(largest_line[2])
    @test isapprox(sum(largest_line[2]),2.0;atol=8eps(Float64),rtol=0)
    @test isapprox(
        _quadrature_moment(largest_line...,(254,)),2/255;
        atol=2e-15,rtol=2e-13)
    @test isapprox(
        _quadrature_moment(largest_line...,(255,)),0.0;
        atol=5e-18,rtol=0)

    largest_triangle=mesh_integration_points(2,"CompositeGauss254")
    @test length(largest_triangle[2])==128^2
    @test isapprox(sum(largest_triangle[2]),0.5;atol=8e-15,rtol=0)
    @test isapprox(
        _quadrature_moment(largest_triangle...,(254,0)),
        _simplex_moment((254,0));atol=2e-17,rtol=3e-13)

    let largest_tetrahedron=
            mesh_integration_points(4,"CompositeGauss197")
        @test length(largest_tetrahedron[2])==1_000_000
        @test isapprox(
            sum(largest_tetrahedron[2]),1/6;atol=8e-15,rtol=0)
    end
    let largest_quadrangle=
            mesh_integration_points(3,"CompositeGauss255")
        @test length(largest_quadrangle[2])==128^2
        @test isapprox(
            sum(largest_quadrangle[2]),4.0;atol=8e-14,rtol=0)
        @test isapprox(
            _quadrature_moment(largest_quadrangle...,(254,0)),4/255;
            atol=2e-15,rtol=3e-13)
    end
    let largest_hexahedron=
            mesh_integration_points(5,"CompositeGauss199")
        @test length(largest_hexahedron[2])==1_000_000
        @test isapprox(
            sum(largest_hexahedron[2]),8.0;atol=8e-13,rtol=0)
        @test isapprox(
            _quadrature_moment(largest_hexahedron...,(198,0,0)),8/199;
            atol=2e-15,rtol=3e-13)
    end
    let largest_prism=mesh_integration_points(6,"CompositeGauss198")
        @test length(largest_prism[2])==1_000_000
        @test isapprox(sum(largest_prism[2]),1.0;atol=8e-14,rtol=0)
        @test isapprox(
            _quadrature_moment(largest_prism...,(198,0,0)),
            2*_simplex_moment((198,0));atol=2e-17,rtol=3e-13)
    end
    let largest_pyramid=mesh_integration_points(7,"CompositeGauss199")
        @test length(largest_pyramid[2])==1_000_000
        @test isapprox(
            sum(largest_pyramid[2]),4/3;atol=8e-14,rtol=0)
        @test isapprox(
            _quadrature_moment(largest_pyramid...,(0,0,199)),
            _pyramid_moment((0,0,199));atol=2e-17,rtol=3e-13)
    end

    mesh_integration_points(4,"CompositeGauss40")
    @test @allocated(mesh_integration_points(
        4,"CompositeGauss40"))<=400_000

    for value in (true,1.0,"1",typemax(UInt128))
        @test_throws ArgumentError mesh_integration_points(value,"Gauss2")
    end
    for element_type in (999,34,140)
        @test_throws ArgumentError mesh_integration_points(
            element_type,"Gauss2")
    end
    for value in (
        :Gauss2,nothing,"","gauss2"," Gauss2","Gauss 2",
        "Gauss+2","Gauss-1","Gauss٢","CompositeGauss2x",
        "Gauss2\0","Gauss999999999999999999999999999999999999")
        @test_throws ArgumentError mesh_integration_points(2,value)
    end
    for (element_type,name) in (
        (2,"Gauss255"),(4,"Gauss198"),(6,"Gauss199"),
        (1,"CompositeGauss256"),(2,"CompositeGauss255"),
        (3,"CompositeGauss256"),(4,"CompositeGauss198"),
        (5,"CompositeGauss200"),(6,"CompositeGauss199"),
        (7,"CompositeGauss200"))
        @test_throws ArgumentError mesh_integration_points(element_type,name)
    end
    oversized_callbacks=Function[
        ()->mesh_integration_points(2,"Gauss255"),
        ()->mesh_integration_points(4,"Gauss198"),
        ()->mesh_integration_points(6,"Gauss199"),
        ()->mesh_integration_points(3,"CompositeGauss256"),
        ()->mesh_integration_points(5,"CompositeGauss200"),
        ()->mesh_integration_points(6,"CompositeGauss199"),
        ()->mesh_integration_points(7,"CompositeGauss200"),
    ]
    foreach(_expect_quadrature_argument_error,oversized_callbacks)
    for callback in oversized_callbacks
        @test @allocated(_expect_quadrature_argument_error(callback))<=50_000
    end
    error_message=try
        mesh_integration_points(2,"Gauss255")
        ""
    catch err
        sprint(showerror,err)
    end
    @test occursin("requires 129 Gauss--Legendre points per axis",error_message)
    @test occursin("the limit is 128",error_message)
    trihedron_message=try
        mesh_integration_points(140,"Gauss2")
        ""
    catch err
        sprint(showerror,err)
    end
    @test occursin("no integration rules in Gmsh 4.15.2",trihedron_message)
    @test isempty(Test.detect_ambiguities(
        Tessella.MeshQuadrature;recursive=true))
end

# Mutant analysis:
# - Reusing the line point-count formula on triangles/tetrahedra is rejected by
#   the requested-degree analytic moment checks for CompositeGauss rules.
# - Omitting either Duffy Jacobian factor is rejected by the constant moments
#   (reference area 1/2 and volume 1/6).
# - Flattening axis-major coordinates instead of point-major triples is rejected
#   by asymmetric mixed moments and the exact economical arrays.
# - Returning shared table storage is rejected by the detached-result mutation.
# - Parsing signs, whitespace, Unicode digits, or overflowing decimal suffixes is
#   rejected by the malformed-rule cases before any quadrature allocation.
# - Reversing tensor loop order is rejected by exact non-simplex Gmsh arrays and
#   asymmetric Cartesian, prism, and pyramid moments.
# - Omitting the pyramid's `(1-w)^2` Gauss--Jacobi measure or one Duffy scale is
#   rejected by its constant, mixed, and vertical analytic moments.
# - Reusing the ordinary line point count in the prism is rejected by Gauss1 and
#   every odd-order composite and high-order Gauss point count.
# - Truncating economical simplex tables at order five is rejected by analytic
#   moments and exact counts through Triangle/Prism order 20 and Tetrahedron
#   order 21; off-by-one tensor transitions are rejected by equality checks.
# - Dispatching by interpolation order is rejected by catalog-wide comparisons
#   across complete, serendipity, and order-zero family types.
# - Allocating before the family-specific total-point preflight is rejected by
#   the first over-limit quadrangle, hexahedron, prism, and pyramid requests.
