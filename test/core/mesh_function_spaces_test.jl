using Test
using Tessella
using Tessella.MeshFunctionSpaces
using Tessella.Elements: MSH_CATALOG, lagrange_nodes, msh_spec, msh_type
using Tessella.MeshEntityTopology: _mesh_edge_topology_for_cells

@noinline _h1_quad_full_allocated(points)=
    @allocated mesh_basis_functions(3,points,"H1Legendre1")
@noinline _h1_quad_selected_allocated(points)=
    @allocated mesh_basis_functions(
        3,points,"H1Legendre1",Int32[0])

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

function _exact_cardinal_value(order,index,x)
    exact_x=Rational{BigInt}(x)
    return Float64(prod(
        (order*exact_x-other)/(index-other)
        for other in 0:order if other!=index;
        init=Rational{BigInt}(1)))
end

function _exact_cardinal_derivative(order,index,x)
    exact_x=Rational{BigInt}(x)
    value=Rational{BigInt}(1)
    derivative=Rational{BigInt}(0)
    for other in 0:order
        other==index && continue
        factor=(order*exact_x-other)/(index-other)
        factor_derivative=order//(index-other)
        derivative=derivative*factor+value*factor_derivative
        value*=factor
    end
    return Float64(derivative)
end

function _exact_falling_value(order,degree,x)
    exact_x=Rational{BigInt}(x)
    return Float64(prod(
        (order*exact_x-offset)/(offset+1)
        for offset in 0:degree-1;
        init=Rational{BigInt}(1)))
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
        ()->mesh_basis_functions(999,[0,0,0],"Lagrange"),
        ()->mesh_basis_functions(2,[0,0,0],:Lagrange),
        ()->mesh_basis_functions(2,[0,0,0],"Lagrange11"),
        ()->mesh_basis_functions(2,[0,0,0],"Lagrange-1"),
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

@testset "order-one H1 on non-simplex reference families" begin
    points=Float64[0.25,-0.5,0.0, 0.1,0.2,-0.3]
    expected=Dict(3=>(4,24),5=>(8,40320),6=>(6,720),15=>(1,1),
                  16=>(4,24),17=>(8,40320),18=>(6,720))
    for (element_type,(vertex_count,orientation_count)) in expected
        lag=mesh_basis_functions(
            element_type,points,"Lagrange1",Int32[0])[2]
        @test length(lag)==2vertex_count
        h1=mesh_basis_functions(element_type,points,"H1Legendre1")
        @test h1[1]==1
        @test h1[3]==orientation_count
        # Every orientation repeats the vertex-function block, exactly as
        # Gmsh 4.15.2 orders its orientation-major H1 output.
        @test h1[2]==repeat(lag,orientation_count)
        grad=mesh_basis_functions(
            element_type,points,"GradLagrange1",Int32[0])[2]
        # Point GradH1Legendre1 is the verified zero gradient even though the
        # pinned release answers [1,0,0] there; the differential pins that
        # Gmsh divergence exactly.
        @test mesh_basis_functions(
            element_type,points,"GradH1Legendre1",Int32[0])==
            (Int32(3),grad,Int32(orientation_count))
        selected=mesh_basis_functions(
            element_type,points,"H1Legendre1",
            orientation_count==1 ? Int32[0] : Int32[orientation_count-1,0])[2]
        @test selected==vcat(
            orientation_count==1 ? Float64[] :
                h1[2][end-2vertex_count+1:end],h1[2][1:2vertex_count])
        @test mesh_number_of_orientations(
            element_type,"H1Legendre1")==orientation_count
        @test mesh_number_of_orientations(
            element_type,"GradH1Legendre1")==orientation_count
        @test mesh_number_of_keys(
            element_type,"H1Legendre1")==vertex_count
        @test mesh_key_dimension(
            element_type,"GradH1Legendre1")==0
        @test mesh_keys_information(
            zeros(Int32,vertex_count),UInt64.(1:vertex_count),
            element_type,"H1Legendre1")==
            fill((Int32(0),element_type==15 ? Int32(0) : Int32(1)),
                 vertex_count)
    end

    # Vertex interpolation is exact on the quadrangle reference: each vertex
    # block selects its own hat function.
    quad_points,quad_count=_first_order_reference_points(3)
    @test quad_count==4
    quad_h1=mesh_basis_functions(
        3,quad_points,"H1Legendre1",Int32[0])[2]
    for vertex in 1:4
        block=quad_h1[4vertex-3:4vertex]
        @test block==Float64[j==vertex for j in 1:4]
    end

    # Hexahedron orientation metadata uses 8! even though the pinned
    # release's getNumberOfOrientations metadata query diverges there.
    @test mesh_basis_functions(
        5,Float64[0.25,-0.5,0.5],"H1Legendre1")[3]==40320
    empty_hex=mesh_basis_functions(5,Float64[],"H1Legendre1")
    @test empty_hex==(Int32(1),Float64[],Int32(40320))
    detached=mesh_basis_functions(3,points,"H1Legendre1")
    original=copy(detached[2])
    detached[2][1]=99
    @test mesh_basis_functions(3,points,"H1Legendre1")[2]==original
    repeated=mesh_basis_functions(3,points,"H1Legendre1")[2]
    @test repeated==original

    # Slice allocation scales with the selected orientations, rejecting a
    # full-compute-then-slice implementation: one of 24 blocks costs less
    # than the 24-block payload, and the full block stays near its payload.
    full_payload=24*2*4*8
    selected_bytes=_h1_quad_selected_allocated(points)
    @test selected_bytes<full_payload
    @test _h1_quad_full_allocated(points)<full_payload+4096

    for invalid in (
        ()->mesh_basis_functions(7,[0,0,0],"H1Legendre1"),
        ()->mesh_basis_functions(7,[0,0,0],"GradH1Legendre1"),
        ()->mesh_basis_functions(19,[0,0,0],"H1Legendre1"),
        ()->mesh_basis_functions(140,[0,0,0],"H1Legendre1"),
        ()->mesh_basis_functions(3,[0,0,0],"HcurlLegendre0"),
        ()->mesh_basis_functions(5,[0,0,0],"HcurlLegendre0"),
        ()->mesh_basis_functions(6,[0,0,0],"CurlHcurlLegendre0"),
        ()->mesh_basis_functions(15,[0,0,0],"HcurlLegendre0"),
        ()->mesh_basis_functions(3,[0,0,0],"H1Legendre2"),
        ()->mesh_basis_functions(3,[0,0,0],"H1Legendre1",Int32[24]),
        ()->mesh_basis_functions(5,[0,0,0],"H1Legendre1",Int32[-1]),
        ()->mesh_basis_functions(6,[0,0,0],"H1Legendre1",Int32[0,0]),
        ()->mesh_basis_functions(
            3,[0,0,0],"H1Legendre1",Int32[0,1,2,3,4,5,6,7,8,9,
                                            10,11,12,13,14,15,16,17,18,19,
                                            20,21,22,23,24]),
    )
        @test_throws ArgumentError invalid()
    end
end

@testset "higher-order nodal reference functions" begin
    line=mesh_basis_functions(8,[0.2,0,0],"Lagrange")
    @test line[1:2:end]==(Int32(1),Int32(1))
    @test isapprox(line[2],Float64[-0.08,0.12,0.96];atol=eps(),rtol=0)
    @test mesh_basis_functions(1,[0.2,0,0],"Lagrange2")==line
    line_gradient=mesh_basis_functions(8,[0.2,0,0],"GradLagrange")
    @test line_gradient[1:2:end]==(Int32(3),Int32(1))
    @test isapprox(
        line_gradient[2],Float64[-0.3,0,0, 0.7,0,0, -0.4,0,0];
        atol=eps(),rtol=0)

    triangle=mesh_basis_functions(9,[0.2,0.3,0],"IsoParametric")
    @test isapprox(
        triangle[2],Float64[0,-0.12,-0.12,0.4,0.24,0.6];
        atol=2eps(Float64),rtol=0)
    @test mesh_number_of_keys(2,"Lagrange2")==6
    @test mesh_keys_information(
        zeros(Int32,10),UInt64.(1:10),21,"Lagrange")==
        vcat(fill((Int32(0),Int32(-1)),9),
             [(Int32(2),Int32(-1))])

    # Multiplication can round p*x to an integer even when the Float64 x is
    # not the corresponding rational node. Such points require the exact path
    # instead of being mistaken for Kronecker zeros.
    near_third=nextfloat(1/3)
    line_values=mesh_basis_functions(
        1,[2near_third-1,0,0],"Lagrange3")[2]
    line_layout=(0,3,1,2)
    exact_line=Float64[
        _exact_cardinal_value(3,index,near_third) for index in line_layout]
    @test line_values[[1,2,4]]==exact_line[[1,2,4]]
    @test all(value->!iszero(value),line_values[[1,2,4]])
    line_gradients=mesh_basis_functions(
        1,[2near_third-1,0,0],"GradLagrange3")[2]
    exact_line_gradients=Float64[
        _exact_cardinal_derivative(3,index,near_third)/2
        for index in line_layout]
    @test line_gradients[[1,4,10]]==exact_line_gradients[[1,2,4]]

    triangle_point=Float64[near_third,0.2,0]
    triangle_nodes=lagrange_nodes(21)
    triangle_node=findfirst(
        node->triangle_nodes[1,node]==2/3 &&
              triangle_nodes[2,node]==0.0,axes(triangle_nodes,2))
    triangle_node===nothing && error("missing order-three triangle node")
    triangle_value=mesh_basis_functions(
        21,triangle_point,"Lagrange")[2][triangle_node]
    origin=1-near_third-triangle_point[2]
    expected_triangle=_exact_falling_value(3,1,origin)*
                      _exact_falling_value(3,2,near_third)
    @test triangle_value==expected_triangle
    @test !iszero(triangle_value)

    incomplete_pyramid=mesh_basis_functions(
        125,[2near_third-1,-1,0],"Lagrange")[2]
    expected_edge=_exact_cardinal_value(3,2,near_third)
    @test isapprox(
        incomplete_pyramid[7],expected_edge;atol=0,rtol=8eps(Float64))
    @test !iszero(incomplete_pyramid[7])
end

@testset "catalog-wide actual- and explicit-order nodal contracts" begin
    interior=Dict(
        :pnt=>Float64[0,0,0],
        :lin=>Float64[0.13,0,0],
        :tri=>Float64[0.2,0.3,0],
        :qua=>Float64[0.13,-0.21,0],
        :tet=>Float64[0.1,0.2,0.15],
        :hex=>Float64[0.13,-0.21,0.17],
        :pri=>Float64[0.2,0.3,-0.17],
        :pyr=>Float64[0.05,-0.08,0.3])

    # This is the first type-105 request in this process: concurrent calls
    # exercise single-owner construction of the cached interpolation data.
    concurrent_tasks=[Threads.@spawn mesh_basis_functions(
        105,interior[:hex],"GradLagrange") for _ in 1:8]
    concurrent_results=fetch.(concurrent_tasks)
    @test all(result->result==concurrent_results[1],concurrent_results)

    for (element_type,spec) in sort!(collect(MSH_CATALOG);by=first)
        spec.family===:trih && continue
        point=interior[spec.family]
        values=mesh_basis_functions(
            element_type,point,"Lagrange")
        gradients=mesh_basis_functions(
            element_type,point,"GradIsoParametric")
        @test values[1:2:end]==(Int32(1),Int32(1))
        @test gradients[1:2:end]==(Int32(3),Int32(1))
        @test length(values[2])==spec.nnodes
        @test length(gradients[2])==3spec.nnodes
        @test all(isfinite,values[2])
        @test all(isfinite,gradients[2])
        @test isapprox(sum(values[2]),1.0;atol=2e-8,rtol=2e-8)
        for axis in 1:3
            @test isapprox(
                sum(gradients[2][axis:3:end]),0.0;atol=2e-7,rtol=0)
        end

        nodes=lagrange_nodes(element_type)
        if spec.order>=1
            @test isapprox(nodes*values[2],point;atol=2e-8,rtol=2e-8)
            expected_jacobian=zeros(3,3)
            for axis in 1:spec.dim
                expected_jacobian[axis,axis]=1
            end
            nodal_gradients=reshape(gradients[2],3,spec.nnodes)
            @test isapprox(
                nodes*transpose(nodal_gradients),expected_jacobian;
                atol=2e-7,rtol=2e-7)
        end
        nodal=mesh_basis_functions(
            element_type,collect(vec(nodes)),"IsoParametric")[2]
        nodal_error=0.0
        for point_index in 1:spec.nnodes,node in 1:spec.nnodes
            expected=point_index==node ? 1.0 : 0.0
            nodal_error=max(
                nodal_error,
                abs(nodal[(point_index-1)*spec.nnodes+node]-expected))
        end
        @test nodal_error<=2e-7
        @test mesh_number_of_keys(element_type,"Lagrange")==spec.nnodes
        @test mesh_number_of_orientations(element_type,"Lagrange")==1
    end

    for (element_type,family,order) in (
        (8,:lin,10),(20,:tri,10),(16,:qua,10),(137,:tet,10),
        (17,:hex,9),(18,:pri,9),(19,:pyr,9))
        complete_type=msh_type(family,order)
        point=interior[family]
        @test mesh_basis_functions(
            element_type,point,"Lagrange$order")==
            mesh_basis_functions(complete_type,point,"Lagrange")
        @test mesh_basis_functions(
            element_type,point,"GradLagrange$order")==
            mesh_basis_functions(complete_type,point,"GradLagrange")
    end
    @test mesh_basis_functions(15,[0,0,0],"Lagrange10")==
          (Int32(1),Float64[1],Int32(1))
    @test mesh_basis_functions(5,interior[:hex],"Lagrange0")==
          (Int32(1),Float64[1],Int32(1))
    @test mesh_basis_functions(5,interior[:hex],"GradLagrange0")==
          (Int32(3),Float64[0,0,0],Int32(1))
    @test mesh_keys_information(
        Int32[0],UInt64[1],88,"Lagrange")==
        [(Int32(3),Int32(-1))]

    # Gmsh 4.15.2 cannot construct its order-3--9 prism bases. The complete
    # prism is independently the tensor product of its triangle and line
    # bases, with local-node correspondence established from public nodes.
    prism_point=Float64[0.217,0.193,-0.271]
    for order in 3:9
        prism_type=msh_type(:pri,order)
        triangle_type=msh_type(:tri,order)
        line_type=msh_type(:lin,order)
        prism_nodes=lagrange_nodes(prism_type)
        triangle_nodes=lagrange_nodes(triangle_type)
        line_nodes=lagrange_nodes(line_type)
        triangle_indices=Dict(
            (triangle_nodes[1,node],triangle_nodes[2,node])=>node
            for node in axes(triangle_nodes,2))
        line_indices=Dict(
            line_nodes[1,node]=>node for node in axes(line_nodes,2))
        triangle_values=mesh_basis_functions(
            triangle_type,[prism_point[1],prism_point[2],0],"Lagrange")[2]
        line_values=mesh_basis_functions(
            line_type,[prism_point[3],0,0],"Lagrange")[2]
        expected=Float64[
            triangle_values[triangle_indices[
                (prism_nodes[1,node],prism_nodes[2,node])]]*
            line_values[line_indices[prism_nodes[3,node]]]
            for node in axes(prism_nodes,2)]
        @test mesh_basis_functions(
            prism_type,prism_point,"Lagrange")[2]==expected
    end

    # The incomplete prism polynomial space contains u^p*w. Its nodal values
    # must reproduce that mode and its analytic gradient at an interior point.
    incomplete_prism_point=Float64[0.41,0.23,0.37]
    for order in 3:9
        element_type=msh_type(:pri,order;serendipity=true)
        nodes=lagrange_nodes(element_type)
        values=mesh_basis_functions(
            element_type,incomplete_prism_point,"Lagrange")[2]
        gradients=mesh_basis_functions(
            element_type,incomplete_prism_point,"GradLagrange")[2]
        nodal_mode=vec(nodes[1,:].^order.*nodes[3,:])
        u,w=incomplete_prism_point[1],incomplete_prism_point[3]
        @test isapprox(sum(values.*nodal_mode),u^order*w;
                       atol=2e-12,rtol=2e-12)
        @test isapprox(sum(gradients[1:3:end].*nodal_mode),
                       order*u^(order-1)*w;atol=2e-11,rtol=2e-11)
        @test isapprox(sum(gradients[2:3:end].*nodal_mode),0.0;
                       atol=2e-11,rtol=0)
        @test isapprox(sum(gradients[3:3:end].*nodal_mode),u^order;
                       atol=2e-11,rtol=2e-11)
    end

    # The native incomplete-pyramid extension restricts to the exact order-p
    # line basis on every edge. This checks non-node points, so the exact-node
    # fast path cannot mask an incorrect edge construction.
    pyramid_edges=((0,1),(0,3),(0,4),(1,2),
                   (1,4),(2,3),(2,4),(3,4))
    edge_parameter=0.371
    for order in 3:9
        element_type=msh_type(:pyr,order;serendipity=true)
        nodes=lagrange_nodes(element_type)
        line_values=mesh_basis_functions(
            1,[2edge_parameter-1,0,0],"Lagrange$order")[2]
        for (edge,(first_zero,second_zero)) in enumerate(pyramid_edges)
            point=collect(
                (1-edge_parameter).*nodes[:,first_zero+1].+
                edge_parameter.*nodes[:,second_zero+1])
            values=mesh_basis_functions(
                element_type,point,"Lagrange")[2]
            expected=zeros(length(values))
            expected[first_zero+1]=line_values[1]
            expected[second_zero+1]=line_values[2]
            first_internal=6+(edge-1)*(order-1)
            expected[first_internal:first_internal+order-2]=line_values[3:end]
            @test isapprox(values,expected;atol=4e-13,rtol=4e-13)
        end
    end

    # The collapsed pyramid coordinate is singular at the apex. Values and
    # analytic gradients must approach the explicitly defined apex limit
    # continuously instead of inheriting the large roundoff seen in Gmsh's
    # near-apex Float64 path.
    for (element_type,spec) in MSH_CATALOG
        spec.family===:pyr && spec.order>=1 || continue
        apex_values=mesh_basis_functions(
            element_type,[0,0,1],"Lagrange")[2]
        near_values=mesh_basis_functions(
            element_type,[0,0,prevfloat(1.0)],"Lagrange")[2]
        apex_gradients=mesh_basis_functions(
            element_type,[0,0,1],"GradLagrange")[2]
        near_gradients=mesh_basis_functions(
            element_type,[0,0,prevfloat(1.0)],"GradLagrange")[2]
        @test isapprox(near_values,apex_values;atol=2e-12,rtol=2e-12)
        @test isapprox(
            near_gradients,apex_gradients;atol=2e-11,rtol=2e-11)
    end

    for element_type in (66,46,51,75,98,110,124,56,61,83,105,117,131)
        spec=msh_spec(element_type)
        point=interior[spec.family]
        gradients=mesh_basis_functions(
            element_type,point,"GradLagrange")[2]
        step=1e-6
        maximum_error=0.0
        for axis in 1:3
            plus=copy(point);minus=copy(point)
            plus[axis]+=step;minus[axis]-=step
            plus_values=mesh_basis_functions(
                element_type,plus,"Lagrange")[2]
            minus_values=mesh_basis_functions(
                element_type,minus,"Lagrange")[2]
            for node in 1:spec.nnodes
                finite_difference=(plus_values[node]-minus_values[node])/(2step)
                maximum_error=max(
                    maximum_error,
                    abs(gradients[3node-3+axis]-finite_difference))
            end
        end
        # The degree-nine/ten incomplete polynomial spaces have large
        # monomial coefficients; a centered Float64 difference amplifies the
        # two value-rounding errors by 1/(2h).
        @test maximum_error<=2e-6
    end

    detached=mesh_basis_functions(105,interior[:hex],"Lagrange")
    reference=copy(detached[2])
    detached[2][1]=99
    @test mesh_basis_functions(105,interior[:hex],"Lagrange")[2]==reference

    for invalid in (
        ()->mesh_basis_functions(2,[0,0,0],"Lagrange000"),
        ()->mesh_basis_functions(2,[0,0,0],"Lagrange1x"),
        ()->mesh_basis_functions(2,[0,0,0],"GradLagrange-1"),
        ()->mesh_basis_functions(98,[0,0,0],"Lagrange10"),
        ()->mesh_basis_functions(124,[0,0,0],"GradLagrange10"),
        ()->mesh_basis_functions(140,[0,0,0],"Lagrange2"),
    )
        @test_throws ArgumentError invalid()
    end

    # Mutants rejected above:
    # - dispatching by the input type instead of the explicit order fails the
    #   complete-type equivalence checks from serendipity inputs;
    # - permuting any higher-order node fails the catalog-wide Kronecker check;
    # - omitting tensor, simplex, prism, or pyramid chain-rule terms fails the
    #   coordinate-gradient reproduction and held-out finite differences;
    # - changing the native unsupported-prism or incomplete-pyramid spaces
    #   fails tensor-product, polynomial-reproduction, or edge-trace oracles;
    # - removing the collapsed-coordinate limits fails the apex-continuity
    #   checks for complete and incomplete pyramids;
    # - an unlocked or externally mutable coefficient cache fails concurrent
    #   first construction or the detached-result check.
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
        ()->mesh_basis_functions(98,[0,0,0],"Lagrange10"),
        ()->mesh_basis_functions(130,[0,0,0],"GradLagrange10"),
        ()->mesh_basis_functions(7,[0,0,0],"H1Legendre1"),
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
    @test mesh_basis_orientations(mesh,3,"H1Legendre1")==Int32[]
    @test mesh_basis_orientations(mesh,5,"H1Legendre1")==Int32[]
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
    @test mesh_number_of_keys(3,"H1Legendre1")==4
    @test mesh_number_of_keys(5,"H1Legendre1")==8
    @test mesh_number_of_keys(6,"H1Legendre1")==6
    @test mesh_number_of_orientations(3,"H1Legendre1")==24
    @test mesh_number_of_orientations(5,"H1Legendre1")==40320
    @test mesh_number_of_orientations(6,"H1Legendre1")==720
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
    @test mesh_keys(mesh,5,"H1Legendre1")==
        (Int32[],UInt64[],Float64[])
    @test mesh_keys_information(
        Int32[0,0,0,0],UInt64[1,2,3,4],3,"H1Legendre1")==
        Tuple{Int32,Int32}[(0,1),(0,1),(0,1),(0,1)]
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
        ()->mesh_keys(mesh,2,"Lagrange2"),
        ()->mesh_keys(mesh,2,"Lagrange";return_coord=1),
        ()->mesh_keys_for_element(mesh,0,"Lagrange"),
        ()->mesh_keys_for_element(mesh,4,"Lagrange"),
        ()->mesh_keys_for_element(mesh,2,"GradLagrange2"),
        ()->mesh_keys_information(Int32[0],UInt64[1],2,"Lagrange"),
        ()->mesh_keys_information(Int32[0,0,0],UInt64[1,2],2,"Lagrange"),
        ()->mesh_keys_information(
            Int32[0,0,0],UInt64[1,2,3],3,"H1Legendre1"),
        ()->mesh_keys_information(Int32[1,1,0],UInt64[1,2,3],2,
                                  "HcurlLegendre0"),
        ()->mesh_keys_information(Int32[1,1,1],UInt64[1,0,3],2,
                                  "HcurlLegendre0"),
    )
        @test_throws ArgumentError invalid()
    end
end
