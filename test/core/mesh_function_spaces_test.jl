using Test
using Tessella
using Tessella.MeshFunctionSpaces
using Tessella.Elements: MSH_CATALOG, msh_spec
using Tessella.MeshEntityTopology: _mesh_edge_topology_for_cells

function _function_space_fixture()
    return Mesh(
        Float64[0 2 0 0;
                0 0 3 0;
                0 0 0 4];
        segs=reshape(Int32[2,1],2,1),
        tris=reshape(Int32[3,1,2],3,1),
        tets=reshape(Int32[4,2,1,3],4,1))
end

function _lexicographic_permutations(count)
    result=Vector{Vector{Int32}}()
    function visit(prefix,remaining)
        if isempty(remaining)
            push!(result,Int32.(prefix))
            return
        end
        for index in eachindex(remaining)
            visit([prefix;remaining[index]],
                  remaining[[i for i in eachindex(remaining) if i!=index]])
        end
    end
    visit(Int[],collect(1:count))
    return result
end

function _orientation_fixture(element_type)
    vertex_count=element_type==1 ? 2 : element_type==2 ? 3 : 4
    coordinates=Float64[0 1 0 0;
                        0 0 1 0;
                        0 0 0 1][:,1:vertex_count]
    cells=hcat(_lexicographic_permutations(vertex_count)...)
    element_type==1 && return Mesh(coordinates;segs=cells)
    element_type==2 && return Mesh(coordinates;tris=cells)
    return Mesh(coordinates;tets=cells)
end

function _first_order_reference_points(element_type)
    _,dimension,_,_,coordinates,primary_nodes=
        Tessella.API.mesh.get_element_properties(element_type)
    result=zeros(Float64,3primary_nodes)
    if dimension>0
        for node in 1:primary_nodes,axis in 1:dimension
            result[3node-3+axis]=coordinates[dimension*(node-1)+axis]
        end
    end
    return result,Int(primary_nodes)
end

function _exact_pyramid_reference(point,gradient::Bool)
    u,v,w=Rational{BigInt}.(point)
    q=1-w
    q==0 && return gradient ? Rational{BigInt}[
        -1//4,-1//4,-1//4, 1//4,-1//4,-1//4,
         1//4,1//4,-1//4, -1//4,1//4,-1//4, 0,0,1] :
        Rational{BigInt}[0,0,0,0,1]
    signs=((-1,-1),(1,-1),(1,1),(-1,1))
    exact=Rational{BigInt}[]
    if gradient
        for (sign_u,sign_v) in signs
            push!(exact,
                  sign_u*(q+sign_v*v)/(4q),
                  sign_v*(q+sign_u*u)/(4q),
                  (sign_u*sign_v*u*v-q^2)/(4q^2))
        end
        append!(exact,(0,0,1))
    else
        for (sign_u,sign_v) in signs
            push!(exact,(q+sign_u*u)*(q+sign_v*v)/(4q))
        end
        push!(exact,w)
    end
    return map(exact) do value
        setprecision(BigFloat,256) do
            Float64(BigFloat(value))
        end
    end
end

@testset "first-order linear-simplex reference function spaces" begin
    points=Float64[0.2,0.3,0.1, -0.4,0.8,-0.2]
    expected_lagrange=Dict(
        1=>[0.4,0.6, 0.7,0.3],
        2=>[0.5,0.2,0.3, 0.6,-0.4,0.8],
        4=>[0.4,0.2,0.3,0.1, 0.8,-0.4,0.8,-0.2])
    expected_gradients=Dict(
        1=>[-0.5,0,0, 0.5,0,0],
        2=>[-1,-1,0, 1,0,0, 0,1,0],
        4=>[-1,-1,-1, 1,0,0, 0,1,0, 0,0,1])
    for element_type in (1,2,4)
        components,values,orientations=mesh_basis_functions(
            element_type,points,"Lagrange")
        @test components==1
        @test orientations==1
        @test isapprox(
            values,expected_lagrange[element_type];atol=0,rtol=eps(Float64))
        @test mesh_basis_functions(
            element_type,points,"IsoParametric")[2]==values
        @test mesh_basis_functions(
            element_type,points,"Lagrange1",Int32[0])[2]==values

        gradient=mesh_basis_functions(
            element_type,points,"GradLagrange")[2]
        @test gradient==repeat(expected_gradients[element_type],2)
        @test mesh_basis_functions(
            element_type,points,"GradIsoParametric")[2]==gradient
        @test mesh_basis_functions(
            element_type,points,"GradLagrange1")[2]==gradient

        hierarchical=mesh_basis_functions(
            element_type,points,"H1Legendre1")
        @test hierarchical[1]==1
        @test hierarchical[3]==factorial(element_type==1 ? 2 :
                                          element_type==2 ? 3 : 4)
        @test hierarchical[2]==repeat(values,Int(hierarchical[3]))
        hierarchical_gradient=mesh_basis_functions(
            element_type,points,"GradH1Legendre1",Int32[0])
        @test hierarchical_gradient==
            (Int32(3),gradient,Int32(hierarchical[3]))
    end

    triangle_hcurl=mesh_basis_functions(
        2,Float64[0.2,0.3,0.1],"HcurlLegendre0",Int32[0])
    @test triangle_hcurl==(
        Int32(3),
        Float64[1.4,0.4,0, -0.6,0.4,0, 0.6,1.6,0],
        Int32(6))
    @test mesh_basis_functions(
        2,Float64[0.2,0.3,0.1],"CurlHcurlLegendre0",Int32[0])==
        (Int32(3),Float64[0,0,4, 0,0,4, 0,0,-4],Int32(6))

    all_triangle=mesh_basis_functions(
        2,Float64[0.2,0.3,0.1],"HcurlLegendre0")[2]
    selected_triangle=mesh_basis_functions(
        2,Float64[0.2,0.3,0.1],"HcurlLegendre0",Int32[5,0,3])[2]
    @test selected_triangle==vcat(
        all_triangle[46:54],all_triangle[1:9],all_triangle[28:36])
    @test mesh_number_of_orientations(1,"HcurlLegendre0")==2
    @test mesh_number_of_orientations(2,"H1Legendre1")==6
    @test mesh_number_of_orientations(4,"CurlHcurlLegendre0")==24
    @test mesh_number_of_orientations(4,"Lagrange")==1

    # The exact fallback preserves a unit residual that ordinary left-to-right
    # arithmetic rounds away at this scale.
    extreme=mesh_basis_functions(
        2,Float64[floatmax(Float64),-floatmax(Float64),0],"Lagrange")[2]
    @test extreme[1]==1.0
    @test extreme[2]==floatmax(Float64)
    @test extreme[3]==-floatmax(Float64)
    @test mesh_basis_functions(
        1,Float64[floatmax(Float64),0,0],
        "HcurlLegendre0",Int32[0])[2]==Float64[1,0,0]

    for malformed in (
        Float64[0],Float64[0,0],Float64[0,0,0,0],
        Float64[NaN,0,0],Float64[Inf,0,0])
        @test_throws ArgumentError mesh_basis_functions(
            2,malformed,"Lagrange")
    end
    for invalid in (
        ()->mesh_basis_functions(true,[0,0,0],"Lagrange"),
        ()->mesh_basis_functions(10,[0,0,0],"Lagrange"),
        ()->mesh_basis_functions(999,[0,0,0],"Lagrange"),
        ()->mesh_basis_functions(2,[0,0,0],:Lagrange),
        ()->mesh_basis_functions(2,[0,0,0],"Lagrange2"),
        ()->mesh_basis_functions(2,[0,0,0],"HcurlLegendre1"),
        ()->mesh_basis_functions(2,[0,0,0],"Lagrange",Int32[1]),
        ()->mesh_basis_functions(2,[0,0,0],"HcurlLegendre0",Int32[-1]),
        ()->mesh_basis_functions(2,[0,0,0],"HcurlLegendre0",Int32[6]),
        ()->mesh_basis_functions(2,[0,0,0],"HcurlLegendre0",Int32[1,1]),
        ()->mesh_basis_functions(2,[0,0,0],"HcurlLegendre0",[true]),
    )
        @test_throws ArgumentError invalid()
    end

    @test isempty(Docs.undocumented_names(
        Tessella.MeshFunctionSpaces;private=false))
    @test isempty(Test.detect_ambiguities(
        Tessella.MeshFunctionSpaces;recursive=true))
end

@testset "fixed-family order-one nodal reference functions" begin
    representatives=(15,1,2,3,4,5,6,7)
    interior=Dict(
        15=>Float64[0,0,0],1=>Float64[0.17,0,0],
        2=>Float64[0.17,0.23,0],3=>Float64[0.17,-0.23,0],
        4=>Float64[0.17,0.23,0.11],5=>Float64[0.17,-0.23,0.11],
        6=>Float64[0.17,0.23,-0.11],7=>Float64[0.1,-0.2,0.3])

    for element_type in representatives
        reference_points,node_count=
            _first_order_reference_points(element_type)
        components,values,orientations=mesh_basis_functions(
            element_type,reference_points,"Lagrange1")
        @test components==1
        @test orientations==1
        @test length(values)==node_count^2
        for point in 1:node_count,node in 1:node_count
            @test isapprox(
                values[(point-1)*node_count+node],point==node ? 1.0 : 0.0;
                atol=8e-15,rtol=0)
        end

        point=interior[element_type]
        nodal=mesh_basis_functions(
            element_type,point,"Lagrange1")[2]
        gradients=mesh_basis_functions(
            element_type,point,"GradLagrange1")[2]
        @test isapprox(sum(nodal),1.0;atol=8e-15,rtol=0)
        @test length(gradients)==3node_count
        for axis in 1:3
            @test isapprox(
                sum(gradients[axis:3:end]),0.0;atol=8e-15,rtol=0)
        end
        @test all(isfinite,nodal)
        @test all(isfinite,gradients)

        if element_type!=15
            step=1e-6
            tolerance=element_type==7 ? 4e-8 : 2e-9
            for axis in 1:3
                plus=copy(point);minus=copy(point)
                plus[axis]+=step;minus[axis]-=step
                plus_values=mesh_basis_functions(
                    element_type,plus,"Lagrange1")[2]
                minus_values=mesh_basis_functions(
                    element_type,minus,"Lagrange1")[2]
                for node in 1:node_count
                    finite_difference=(plus_values[node]-minus_values[node])/
                                      (2step)
                    @test isapprox(
                        gradients[3node-3+axis],finite_difference;
                        atol=tolerance,rtol=tolerance)
                end
            end
        end
    end

    for element_type in representatives,
        (plain,explicit) in (("Lagrange","Lagrange1"),
                             ("IsoParametric","Lagrange1"),
                             ("GradLagrange","GradLagrange1"),
                             ("GradIsoParametric","GradLagrange1"))
        point=interior[element_type]
        @test mesh_basis_functions(element_type,point,plain)==
              mesh_basis_functions(element_type,point,explicit)
    end

    catalog_points=Float64[
        0.13,-0.21,0.37,
        -0.4,0.7,-0.2,
        0.0,0.0,1.0]
    family_representative=Dict(
        msh_spec(element_type).family=>element_type
        for element_type in representatives)
    family_reference=Dict(
        family=>mesh_basis_functions(
            element_type,catalog_points,"Lagrange1")
        for (family,element_type) in family_representative)
    family_gradient_reference=Dict(
        family=>mesh_basis_functions(
            element_type,catalog_points,"GradLagrange1")
        for (family,element_type) in family_representative)
    family_counts=Dict(
        family=>Int(length(reference[2])÷(length(catalog_points)÷3))
        for (family,reference) in family_reference)
    for (element_type,spec) in MSH_CATALOG
        spec.family===:trih && continue
        expected=family_reference[spec.family]
        @test mesh_basis_functions(
            element_type,catalog_points,"Lagrange1")==expected
        @test mesh_basis_functions(
            element_type,catalog_points,"GradLagrange1")==
            family_gradient_reference[spec.family]
        @test mesh_number_of_orientations(
            element_type,"Lagrange1")==1
        @test mesh_number_of_keys(
            element_type,"Lagrange1")==family_counts[spec.family]
        @test mesh_key_dimension(element_type,"Lagrange1")==0
        key_count=family_counts[spec.family]
        @test mesh_keys_information(
            zeros(Int32,key_count),UInt64.(1:key_count),
            element_type,"Lagrange1")==
            fill((Int32(0),Int32(1)),key_count)
    end

    @test mesh_basis_functions(7,[0,0,1],"Lagrange1")[2]==
          Float64[0,0,0,0,1]
    apex_gradient=Float64[
        -0.25,-0.25,-0.25,
         0.25,-0.25,-0.25,
         0.25,0.25,-0.25,
        -0.25,0.25,-0.25,
         0,0,1]
    @test mesh_basis_functions(7,[0,0,1],"GradLagrange1")[2]==
          apex_gradient
    @test mesh_basis_functions(
        7,[1e-8,-1e-8,1],"GradLagrange1")[2]==apex_gradient
    for height in (prevfloat(1.0),nextfloat(1.0))
        q=1-height
        @test mesh_basis_functions(
            7,[0,0,height],"Lagrange1")[2]==
            Float64[q/4,q/4,q/4,q/4,height]
        @test mesh_basis_functions(
            7,[0,0,height],"GradLagrange1")[2]==apex_gradient
    end
    maximum=floatmax(Float64)
    for extreme_pyramid in (
            Float64[maximum,maximum,maximum],
            Float64[-maximum,maximum,maximum],
            Float64[maximum,-maximum,-maximum],
            Float64[-maximum,-maximum,-maximum]),
        (space,gradient) in (("Lagrange1",false),
                             ("GradLagrange1",true))
        expected=_exact_pyramid_reference(extreme_pyramid,gradient)
        values=mesh_basis_functions(7,extreme_pyramid,space)[2]
        @test all(isfinite,expected)
        @test values==expected
        @test all(isfinite,values)
    end

    empty_basis=mesh_basis_functions(14,Float64[],"Lagrange1")
    @test empty_basis==(Int32(1),Float64[],Int32(1))
    detached=mesh_basis_functions(12,catalog_points,"Lagrange1")
    original=copy(detached[2])
    detached[2][1]=99
    @test mesh_basis_functions(
        12,catalog_points,"Lagrange1")[2]==original

    for invalid in (
        ()->mesh_basis_functions(10,[0,0,0],"Lagrange"),
        ()->mesh_basis_functions(12,[0,0,0],"GradIsoParametric"),
        ()->mesh_basis_functions(3,[0,0,0],"H1Legendre1"),
        ()->mesh_basis_functions(5,[0,0,0],"HcurlLegendre0"),
        ()->mesh_basis_functions(140,[0,0,0],"Lagrange1"),
        ()->mesh_basis_functions(
            5,[floatmax(Float64),floatmax(Float64),floatmax(Float64)],
            "Lagrange1"),
    )
        @test_throws ArgumentError invalid()
    end
end

@testset "linear-simplex orientations and finite-element keys" begin
    mesh=_function_space_fixture()
    @test mesh_basis_orientations(mesh,1,"Lagrange")==Int32[0]
    @test mesh_basis_orientations(mesh,1,"HcurlLegendre0")==Int32[1]
    @test mesh_basis_orientations(mesh,2,"HcurlLegendre0")==Int32[4]
    @test mesh_basis_orientations(mesh,4,"H1Legendre1")==Int32[20]
    @test mesh_basis_orientation(mesh,1,"H1Legendre1")==1
    @test mesh_basis_orientation(mesh,2,"HcurlLegendre0")==4
    @test mesh_basis_orientation(mesh,3,"CurlHcurlLegendre0")==20
    @test mesh_basis_orientation(mesh,3,"GradLagrange")==0
    for element_type in (1,2,4)
        permutation_mesh=_orientation_fixture(element_type)
        expected=Int32.(0:factorial(size(
            element_type==1 ? permutation_mesh.segs :
            element_type==2 ? permutation_mesh.tris : permutation_mesh.tets,1))-1)
        @test mesh_basis_orientations(
            permutation_mesh,element_type,"HcurlLegendre0")==expected
        @test mesh_basis_orientations(
            permutation_mesh,element_type,"H1Legendre1")==expected
    end

    @test mesh_number_of_keys(1,"Lagrange")==2
    @test mesh_number_of_keys(2,"H1Legendre1")==3
    @test mesh_number_of_keys(4,"HcurlLegendre0")==6
    @test mesh_key_dimension(4,"GradLagrange")==0
    @test mesh_key_dimension(4,"CurlHcurlLegendre0")==1

    nodal=mesh_keys(mesh,2,"Lagrange")
    @test nodal[1]==Int32[0,0,0]
    @test nodal[2]==UInt64[3,1,2]
    @test nodal[3]==[0,3,0, 0,0,0, 2,0,0]
    nodal[1][1]=99;nodal[2][1]=99;nodal[3][1]=99
    @test mesh_keys(mesh,2,"Lagrange")[2]==UInt64[3,1,2]
    @test mesh_keys(mesh,2,"Lagrange";return_coord=false)[3]==Float64[]

    topology=_mesh_edge_topology_for_cells(mesh,nothing,mesh.segs,1)
    topology=_mesh_edge_topology_for_cells(mesh,topology,mesh.tris,2)
    triangle=mesh_keys(mesh,2,"HcurlLegendre0",topology)
    @test triangle[1]==Int32[1,1,1]
    @test triangle[2]==UInt64[2,1,3]
    @test triangle[3]==[0,1.5,0, 1,0,0, 1,1.5,0]
    topology=_mesh_edge_topology_for_cells(mesh,topology,mesh.tets,4)
    @test _mesh_edge_topology_for_cells(mesh,topology,mesh.tets,4)===topology
    tetrahedron=mesh_keys_for_element(
        mesh,3,"CurlHcurlLegendre0",topology)
    @test tetrahedron[2]==UInt64[4,1,5,6,2,3]
    @test tetrahedron[3]==
        [1,0,2, 1,0,0, 0,0,2, 0,1.5,2, 0,1.5,0, 1,1.5,0]

    @test mesh_keys_information(
        Int32[0,0,0],UInt64[3,1,2],2,"Lagrange")==
        Tuple{Int32,Int32}[(0,-1),(0,-1),(0,-1)]
    @test mesh_keys_information(
        Int32[0,0,0],UInt64[3,1,2],2,"H1Legendre1")==
        Tuple{Int32,Int32}[(0,1),(0,1),(0,1)]
    @test mesh_keys_information(
        Int32[1,1,1],UInt64[2,1,3],2,"HcurlLegendre0")==
        Tuple{Int32,Int32}[(1,0),(1,0),(1,0)]
    @test mesh_keys_information(
        Int32[],UInt64[],2,"HcurlLegendre0")==Tuple{Int32,Int32}[]

    for invalid in (
        ()->mesh_keys(mesh,2,"HcurlLegendre0"),
        ()->mesh_keys(mesh,2,"Lagrange";return_coord=1),
        ()->mesh_keys_for_element(mesh,0,"Lagrange"),
        ()->mesh_keys_for_element(mesh,4,"Lagrange"),
        ()->mesh_keys_information(Int32[0],UInt64[1],2,"Lagrange"),
        ()->mesh_keys_information(Int32[0,0,0],UInt64[1,2],2,"Lagrange"),
        ()->mesh_keys_information(Int32[1,1,0],UInt64[1,2,3],2,
                                  "HcurlLegendre0"),
        ()->mesh_keys_information(Int32[1,1,1],UInt64[1,0,3],2,
                                  "HcurlLegendre0"),
    )
        @test_throws ArgumentError invalid()
    end
end
