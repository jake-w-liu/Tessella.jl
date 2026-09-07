using Test
using Tessella
using Tessella.MeshQuadrature
using Tessella.Elements: MSH_CATALOG

function _line_moment(power::Int)
    return isodd(power) ? 0.0 : 2.0/(power+1)
end

function _simplex_moment(exponents::Tuple)
    numerator=prod(factorial(big(exponent)) for exponent in exponents)
    denominator=factorial(big(sum(exponents)+length(exponents)))
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

function _check_rule_moments(element_type::Int,name::String,order::Int;
                             atol::Float64)
    coordinates,weights=mesh_integration_points(element_type,name)
    if element_type==1
        for power in 0:order
            @test isapprox(
                _quadrature_moment(coordinates,weights,(power,)),
                _line_moment(power);atol=atol,rtol=atol)
        end
    elseif element_type==2
        for first in 0:order,second in 0:(order-first)
            exponents=(first,second)
            @test isapprox(
                _quadrature_moment(coordinates,weights,exponents),
                _simplex_moment(exponents);atol=atol,rtol=atol)
        end
    else
        for first in 0:order,second in 0:(order-first),
            third in 0:(order-first-second)
            exponents=(first,second,third)
            @test isapprox(
                _quadrature_moment(coordinates,weights,exponents),
                _simplex_moment(exponents);atol=atol,rtol=atol)
        end
    end
    return coordinates,weights
end

@testset "Gmsh-shaped simplex reference quadrature" begin
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

    economical_counts=Dict(
        1=>[1,1,2,2,3,3],
        2=>[1,1,3,4,6,7],
        4=>[1,1,4,5,11,14])
    for element_type in (1,2,4),order in 0:5
        coordinates,weights=_check_rule_moments(
            element_type,"Gauss$order",order;atol=8e-13)
        @test length(weights)==economical_counts[element_type][order+1]
        @test length(coordinates)==3length(weights)
        @test all(isfinite,coordinates)
        @test all(isfinite,weights)
    end

    for order in (0,1,2,3,4,5,6,12,20)
        for element_type in (1,2,4)
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
            else
                u=coordinates[1:3:end]
                v=coordinates[2:3:end]
                w=coordinates[3:3:end]
                @test all(u .>= 0)
                @test all(v .>= 0)
                @test all(w .>= 0)
                @test all(
                    index->u[index]+v[index]+w[index]<=1,
                    eachindex(u,v,w))
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
        for (family,type) in ((:pnt,15),(:lin,1),(:tri,2),(:tet,4)))
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

    mesh_integration_points(4,"CompositeGauss40")
    @test @allocated(mesh_integration_points(
        4,"CompositeGauss40"))<=400_000

    for value in (true,1.0,"1",typemax(UInt128))
        @test_throws ArgumentError mesh_integration_points(value,"Gauss2")
    end
    for element_type in (999,34,3,5,6,7,140)
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
        (2,"Gauss6"),(4,"Gauss6"),
        (1,"CompositeGauss256"),(2,"CompositeGauss255"),
        (4,"CompositeGauss198"))
        @test_throws ArgumentError mesh_integration_points(element_type,name)
    end
    error_message=try
        mesh_integration_points(4,"Gauss6")
        ""
    catch err
        sprint(showerror,err)
    end
    @test occursin("use CompositeGauss6",error_message)
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
