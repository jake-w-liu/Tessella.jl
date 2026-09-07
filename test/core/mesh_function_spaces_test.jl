using Test
using Tessella
using Tessella.MeshFunctionSpaces
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
        ()->mesh_basis_functions(3,[0,0,0],"Lagrange"),
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
