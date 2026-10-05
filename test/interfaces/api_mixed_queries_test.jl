using Test
using LinearAlgebra
using Tessella
using Tessella.Elements: MixedMesh, ElementBlock, lagrange_nodes

const _MIXED_QUERY_API=Tessella.API
const _MIXED_QUERY_A=Float64[2.0 0.25 0.5;0.5 3.0 0.25;0.25 0.5 4.0]
const _MIXED_QUERY_ORIGIN=Float64[10,-5,3]

function _mixed_query_fixture(msh;warped=false)
    coordinates=_MIXED_QUERY_A*lagrange_nodes(msh) .+ _MIXED_QUERY_ORIGIN
    if warped
        coordinates[3,end]+=0.7
    end
    return MixedMesh(coordinates,[ElementBlock(msh,
        reshape(Int32.(1:size(coordinates,2)),size(coordinates,2),1))])
end

function _mixed_query_install!(mesh)
    lock(_MIXED_QUERY_API.STATE_LOCK) do
        _MIXED_QUERY_API._replace_mesh_cache_locked!(
            _MIXED_QUERY_API._copy_mesh(mesh))
    end
end

@testset "native mixed reference queries through API" begin
    _MIXED_QUERY_API.finalize()
    try
        _MIXED_QUERY_API.initialize()
        for (msh,uvw) in ((1,(-0.25,0.0,0.0)),(2,(0.2,0.3,0.0)),
                         (3,(-0.2,0.3,0.0)),(4,(0.1,0.2,0.3)),
                         (5,(-0.2,0.3,-0.4)),(6,(0.2,0.3,-0.4)),
                         (7,(-0.2,0.3,0.4)),(15,(0.0,0.0,0.0)))
            fixture=_mixed_query_fixture(msh)
            _mixed_query_install!(fixture)
            physical=_MIXED_QUERY_A*collect(uvw)+_MIXED_QUERY_ORIGIN
            result=_MIXED_QUERY_API.mesh.get_element_by_coordinates(physical...,-1,true)
            @test result[1:3]==(UInt64(1),Int32(msh),UInt64.(1:size(fixture.coords,2)))
            @test collect(result[4:6])≈collect(uvw) atol=2e-14
            @test collect(_MIXED_QUERY_API.mesh.get_local_coordinates_in_element(
                1,physical...))≈collect(uvw) atol=2e-14
            @test _MIXED_QUERY_API.mesh.get_elements_by_coordinates(
                physical...,-1,true)==UInt64[1]
            locator=_MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]
            @test locator!==nothing
            _MIXED_QUERY_API.mesh.get_elements_by_coordinates(physical...,-1,true)
            @test _MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]===locator
            result[3][1]=999
            @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(
                physical...,-1,true)[3][1]==1

            jacobians,determinants,coordinates=_MIXED_QUERY_API.mesh.get_jacobian(1,collect(uvw))
            @test coordinates≈physical atol=2e-14
            @test length(jacobians)==9 && length(determinants)==1
            dimension=Tessella.Elements.msh_dimension(msh)
            if dimension>0
                @test reshape(jacobians,3,3)[:,1:dimension]≈_MIXED_QUERY_A[:,1:dimension] atol=2e-14
            else
                @test reshape(jacobians,3,3)==Matrix{Float64}(I,3,3)
            end
            @test determinants[1]>0
            @test _MIXED_QUERY_API.mesh.get_jacobians(msh,collect(uvw))==
                (jacobians,determinants,coordinates)
            @test _MIXED_QUERY_API.mesh.get_jacobians(msh,collect(uvw),-1,0,2)==
                (Float64[],Float64[],Float64[])
            @test _MIXED_QUERY_API.mesh.get_jacobians(msh,collect(uvw),-1,1,2)==
                (jacobians,determinants,coordinates)
            @test _MIXED_QUERY_API.mesh.get_basis_functions_orientation(msh,"Lagrange")==Int32[0]
            @test _MIXED_QUERY_API.mesh.get_basis_functions_orientation_for_element(1,"Lagrange")==0
            keys=_MIXED_QUERY_API.mesh.get_keys(msh,"Lagrange")
            @test keys==(zeros(Int32,size(fixture.coords,2)),
                         UInt64.(1:size(fixture.coords,2)),vec(fixture.coords))
            @test _MIXED_QUERY_API.mesh.get_keys_for_element(1,"Lagrange")==keys
            @test isempty(_MIXED_QUERY_API.mesh.get_keys(msh,"Lagrange",-1,false)[3])
            if msh!=7
                @test _MIXED_QUERY_API.mesh.get_basis_functions_orientation(msh,"H1Legendre2")==Int32[0]
                hierarchical=_MIXED_QUERY_API.mesh.get_keys(msh,"H1Legendre2")
                @test length(hierarchical[1])==_MIXED_QUERY_API.mesh.get_number_of_keys(msh,"H1Legendre2")
                @test _MIXED_QUERY_API.mesh.get_keys_for_element(1,"H1Legendre2")==hierarchical
            end
            if msh in (1,2,4)
                @test all(isfinite,_MIXED_QUERY_API.mesh.get_element_qualities([1],"volume"))
            elseif msh!=15
                @test_throws ArgumentError _MIXED_QUERY_API.mesh.get_element_qualities([1],"volume")
            end
            if msh!=15
                @test _MIXED_QUERY_API.mesh.get_element_qualities([1],"minEdge")[1]>0
                @test _MIXED_QUERY_API.mesh.get_element_qualities([1],"maxEdge")[1]>=
                    _MIXED_QUERY_API.mesh.get_element_qualities([1],"minEdge")[1]
            end
        end

        # Non-affine quadrangle, hexahedron, prism and pyramid maps are inverted
        # as their native families, rather than through a simplex decomposition.
        for (msh,uvw) in ((3,(-0.3,0.4,0.0)),(5,(-0.3,0.4,0.2)),
                         (6,(0.2,0.3,0.4)),(7,(-0.2,0.1,0.4)))
            _mixed_query_install!(_mixed_query_fixture(msh;warped=true))
            _,_,physical=_MIXED_QUERY_API.mesh.get_jacobian(1,collect(uvw))
            recovered=collect(_MIXED_QUERY_API.mesh.get_local_coordinates_in_element(1,physical...))
            @test recovered[1:2]≈collect(uvw)[1:2] atol=2e-13
            @test recovered[3]≈uvw[3] atol=2e-13
            @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(physical...,-1,true)[2]==msh
        end

        _mixed_query_install!(_mixed_query_fixture(7))
        tip=_MIXED_QUERY_API.LAST_MESH[].coords[:,5]
        @test _MIXED_QUERY_API.mesh.get_local_coordinates_in_element(1,tip...)==(0.0,0.0,1.0)
        @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(tip...,-1,true)[2]==7
        @test_throws ArgumentError _MIXED_QUERY_API.mesh.get_keys(7,"H1Legendre1")
        @test_throws ArgumentError _MIXED_QUERY_API.mesh.get_jacobian(0,[0,0,0])
        @test_throws ArgumentError _MIXED_QUERY_API.mesh.get_jacobian(1,[NaN,0,0])
        @test_throws ArgumentError _MIXED_QUERY_API.mesh.get_element_by_coordinates(Inf,0,0)
        @test_throws ArgumentError _MIXED_QUERY_API.mesh.get_element_by_coordinates(1000,1000,1000,-1,true)

        for (msh,uvw) in ((8,(-0.25,0.0,0.0)),(9,(0.2,0.3,0.0)),
                         (10,(-0.2,0.3,0.0)),(11,(0.1,0.2,0.3)),
                         (12,(-0.2,0.3,-0.4)),(13,(0.2,0.3,-0.4)),
                         (14,(-0.2,0.3,0.4)))
            _mixed_query_install!(_mixed_query_fixture(msh))
            physical=_MIXED_QUERY_A*collect(uvw)+_MIXED_QUERY_ORIGIN
            @test collect(_MIXED_QUERY_API.mesh.get_local_coordinates_in_element(1,physical...))≈
                collect(uvw) atol=2e-13
            @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(physical...,-1,true)[2]==msh
            frame,determinants,coordinates=_MIXED_QUERY_API.mesh.get_jacobian(1,collect(uvw))
            dim=Tessella.Elements.msh_dimension(msh)
            @test reshape(frame,3,3)[:,1:dim]≈_MIXED_QUERY_A[:,1:dim] atol=2e-13
            @test coordinates≈physical atol=2e-13
            @test determinants[1]>0
            @test length(_MIXED_QUERY_API.mesh.get_keys(msh,"Lagrange")[1])==
                size(lagrange_nodes(msh),2)
            primary_edge_counts=Dict(8=>1,9=>3,10=>4,11=>6,12=>12,13=>9,14=>8)
            first_midpoint=Dict(8=>3,9=>4,10=>5,11=>5,12=>9,13=>7,14=>6)
            all_edge_nodes=_MIXED_QUERY_API.mesh.get_element_edge_nodes(msh)
            @test length(all_edge_nodes)==3primary_edge_counts[msh]
            @test all_edge_nodes[1:3]==UInt64[1,2,first_midpoint[msh]]
            @test length(_MIXED_QUERY_API.mesh.get_element_edge_nodes(msh,-1,true))==
                2primary_edge_counts[msh]
            if msh==11
                @test _MIXED_QUERY_API.mesh.get_element_face_nodes(msh,3)[1:6]==
                    UInt64[1,3,2,7,6,5]
            elseif msh==12
                @test _MIXED_QUERY_API.mesh.get_element_face_nodes(msh,4)[1:9]==
                    UInt64[1,4,3,2,10,14,12,9,21]
            end
            if msh!=14
                @test _MIXED_QUERY_API.mesh.get_basis_functions_orientation_for_element(1,"H1Legendre2")==0
                @test length(_MIXED_QUERY_API.mesh.get_keys(msh,"H1Legendre2")[1])==
                    _MIXED_QUERY_API.mesh.get_number_of_keys(msh,"H1Legendre2")
            else
                @test _MIXED_QUERY_API.mesh.get_barycenters(msh,-1,false,true)≈
                    _MIXED_QUERY_ORIGIN+_MIXED_QUERY_A*[0,0,1/5] atol=2e-13
                @test _MIXED_QUERY_API.mesh.get_barycenters(msh,-1,false,false)≈
                    _MIXED_QUERY_ORIGIN+_MIXED_QUERY_A*[0,0,3/14] atol=2e-13
                @test _MIXED_QUERY_API.mesh.get_barycenters(msh,-1,true,true)≈
                    5 .* _MIXED_QUERY_ORIGIN+_MIXED_QUERY_A*[0,0,1] atol=2e-13
            end
        end

        # A curved quadratic line extends past every nodal y coordinate.
        # Its true map must still be considered by the locator.
        curved=MixedMesh(Float64[-1 1 0;0 1 1;0 0 0],
                         [ElementBlock(8,reshape(Int32[1,2,3],3,1))])
        _mixed_query_install!(curved)
        @test collect(_MIXED_QUERY_API.mesh.get_local_coordinates_in_element(1,0.5,1.125,0.0))≈
            [0.5,0.0,0.0] atol=2e-13
        @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(0.5,1.125,0.0,1,true)[2]==8
        @test_throws ArgumentError _MIXED_QUERY_API.mesh.get_element_qualities([1],"volume")

        # Exact curved triangular map (u,v) -> (u,v,u*v), including its
        # independently derived tangent columns and area Jacobian.
        reference_triangle=lagrange_nodes(9)
        curved_triangle=copy(reference_triangle)
        curved_triangle[3,:].=reference_triangle[1,:].*reference_triangle[2,:]
        _mixed_query_install!(MixedMesh(curved_triangle,
            [ElementBlock(9,reshape(Int32.(1:6),6,1))]))
        @test collect(_MIXED_QUERY_API.mesh.get_local_coordinates_in_element(1,0.2,0.3,0.06))≈
            [0.2,0.3,0.0] atol=2e-13
        @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(0.2,0.3,0.06,2,true)[2]==9
        frame,determinants,coordinates=_MIXED_QUERY_API.mesh.get_jacobian(1,[0.2,0.3,0.0])
        @test reshape(frame,3,3)[:,1:2]≈Float64[1 0;0 1;0.3 0.2] atol=2e-13
        @test determinants[1]≈sqrt(1.13) atol=2e-13
        @test coordinates≈[0.2,0.3,0.06] atol=2e-13

        reflected_hex=_mixed_query_fixture(12)
        reflected_hex.coords[1,:].*=-1
        _mixed_query_install!(reflected_hex)
        @test _MIXED_QUERY_API.mesh.get_jacobian(1,[0.1,0.2,0.3])[2][1]<0

        # Nearly parallel quad tangents need the exact Gram-system fallback.
        thin=MixedMesh(Float64[0 1 2 1;0 0 1e-10 1e-10;0 0 0 0],
                       [ElementBlock(3,reshape(Int32[1,2,3,4],4,1))])
        _mixed_query_install!(thin)
        @test collect(_MIXED_QUERY_API.mesh.get_local_coordinates_in_element(1,0.9,4e-11,0.0))≈
            [0.0,-0.2,0.0] atol=2e-13
        @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(0.9,4e-11,0.0,2,true)[2]==3

        # A disjoint curved quadratic hex has no real inverse for the first
        # cell's center. Its conservative candidate bounds must not abort the
        # successful lookup in the other cell.
        reference=lagrange_nodes(12)
        disjoint=hcat(reference,copy(reference))
        disjoint[1,28:54].=100 .+ reference[1,:] .+ 0.3 .* reference[1,:].^2
        _mixed_query_install!(MixedMesh(disjoint,[ElementBlock(12,
            hcat(Int32.(1:27),Int32.(28:54)))]))
        @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(0,0,0,3,true)[1]==1
        @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(100,0,0,3,true)[1]==2
        jacobians,determinants,physical=_MIXED_QUERY_API.mesh.get_jacobians(12,[0,0,0])
        @test determinants≈[1.0,1.0] atol=2e-13
        @test physical≈[0.0,0.0,0.0,100.0,0.0,0.0] atol=2e-13
        @test jacobians[1:9]≈jacobians[10:18] atol=2e-13
        @test _MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]!==nothing
        current_model=_MIXED_QUERY_API.model.get_current()
        _MIXED_QUERY_API.model.add("mixed_query_resource_check")
        @test _MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]===nothing
        _MIXED_QUERY_API.model.set_current(current_model)
        @test _MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]===nothing
        @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(0,0,0,3,true)[1]==1
        _MIXED_QUERY_API.model.remove()
        @test _MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]===nothing
        _MIXED_QUERY_API.mesh.clear()
        @test _MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]===nothing
        _MIXED_QUERY_API.model.add_discrete_entity(3,7)
        reference=lagrange_nodes(14)
        _MIXED_QUERY_API.mesh.add_nodes(3,7,collect(1:14),vec(reference))
        _MIXED_QUERY_API.mesh.add_elements_by_type(7,14,[17],collect(1:14))
        @test _MIXED_QUERY_API.mesh.get_barycenters(14,-1,false,true)≈[0,0,1/5] atol=2e-13
        @test _MIXED_QUERY_API.mesh.get_barycenters(14,-1,false,false)≈[0,0,3/14] atol=2e-13
        @test _MIXED_QUERY_API.mesh.get_barycenters(14,-1,true,true)≈[0,0,1] atol=2e-13
        @test _MIXED_QUERY_API.mesh.get_barycenters(14,-1,true,false)≈[0,0,3] atol=2e-13
    finally
        _MIXED_QUERY_API.finalize()
    end
end


@testset "unclassified mixed Point location keeps geometric adapters" begin
    _MIXED_QUERY_API.finalize()
    try
        _MIXED_QUERY_API.initialize()
        mesh=MixedMesh(Float64[0 1 0;0 0 1;0 0 0],
            [ElementBlock(15,reshape(Int32[1],1,1)),
             ElementBlock(1,reshape(Int32[1,2],2,1)),
             ElementBlock(2,reshape(Int32[1,2,3],3,1))])
        _mixed_query_install!(mesh)
        published=_MIXED_QUERY_API.LAST_MESH[]
        @test _MIXED_QUERY_API.LAST_MESH_CLASS[]===nothing
        @test _MIXED_QUERY_API.mesh.get_elements_by_coordinates(0.,0.,0.,-1,true)==UInt64[3,2,1]
        locator=_MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]
        @test _MIXED_QUERY_API.mesh.get_elements_by_coordinates(0.,0.,0.,-1,true)==UInt64[3,2,1]
        @test _MIXED_QUERY_API.LAST_MIXED_MESH_LOCATOR[]===locator
        for (dimension,element) in ((0,1),(1,2),(2,3))
            @test _MIXED_QUERY_API.mesh.get_elements_by_coordinates(0.,0.,0.,dimension,true)==UInt64[element]
        end
        @test _MIXED_QUERY_API.mesh.get_element_by_coordinates(0.,0.,0.,0,true)==
            (UInt64(1),Int32(15),UInt64[1],0.,0.,0.)
        jac=(Float64[1,0,0,0,1,0,0,0,1],Float64[1],Float64[0,0,0])
        @test _MIXED_QUERY_API.mesh.get_jacobian(1,[0.,0.,0.])==jac
        @test _MIXED_QUERY_API.mesh.get_jacobians(15,[0.,0.,0.])==jac
        @test _MIXED_QUERY_API.mesh.get_keys(15,"Lagrange")==
            (Int32[0],UInt64[1],Float64[0,0,0])
        @test _MIXED_QUERY_API.LAST_MESH[]===published
        @test _MIXED_QUERY_API.LAST_MESH_CLASS[]===nothing
        @test_throws ArgumentError _MIXED_QUERY_API.mesh.get_elements_by_coordinates(4.,4.,4.,-1,true)
    finally
        _MIXED_QUERY_API.finalize()
    end
end
