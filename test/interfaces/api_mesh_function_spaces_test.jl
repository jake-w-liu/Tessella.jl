using Test
using Tessella
using Tessella.MeshTypes: mesh_crc

const _MESH_FUNCTION_API=Tessella.API

function _mesh_function_api_fixture()
    return Mesh(
        Float64[0 2 0 0;
                0 0 3 0;
                0 0 0 4];
        segs=reshape(Int32[2,1],2,1),
        tris=reshape(Int32[3,1,2],3,1),
        tets=reshape(Int32[4,2,1,3],4,1))
end

function _install_mesh_function_fixture!(mesh)
    lock(_MESH_FUNCTION_API.STATE_LOCK) do
        _MESH_FUNCTION_API._replace_mesh_cache_locked!(
            _MESH_FUNCTION_API._copy_mesh(mesh))
    end
    return nothing
end

@testset "function spaces and keys through API" begin
    _MESH_FUNCTION_API.finalize()
    @test _MESH_FUNCTION_API.mesh.get_number_of_orientations(
        4,"HcurlLegendre0")==24
    @test _MESH_FUNCTION_API.mesh.get_number_of_keys(
        4,"HcurlLegendre0")==6
    @test _MESH_FUNCTION_API.mesh.get_number_of_orientations(
        14,"Lagrange1")==1
    @test _MESH_FUNCTION_API.mesh.get_number_of_keys(
        12,"GradLagrange1")==8
    @test _MESH_FUNCTION_API.mesh.get_basis_functions(
        2,[0.2,0.3,0.1],"Lagrange")==
        (Int32(1),[0.5,0.2,0.3],Int32(1))
    quadrangle=_MESH_FUNCTION_API.mesh.get_basis_functions(
        3,[0.1,0.2,0.3],"Lagrange")
    @test quadrangle[1]==1
    @test quadrangle[3]==1
    @test isapprox(
        quadrangle[2],[0.18,0.22,0.33,0.27];atol=eps(Float64),rtol=0)
    @test _MESH_FUNCTION_API.mesh.get_basis_functions(
        14,[0,0,1],"Lagrange1")==
        (Int32(1),Float64[0,0,0,0,1],Int32(1))
    quadratic=_MESH_FUNCTION_API.mesh.get_basis_functions(
        10,[0.1,0.2,0.0],"Lagrange")
    @test quadratic[1:2:end]==(Int32(1),Int32(1))
    @test length(quadratic[2])==9
    @test isapprox(sum(quadratic[2]),1.0;atol=8eps(),rtol=0)
    @test _MESH_FUNCTION_API.mesh.get_basis_functions(
        3,[0.1,0.2,0.0],"Lagrange2")==quadratic
    @test _MESH_FUNCTION_API.mesh.get_number_of_keys(20,"Lagrange")==9
    @test _MESH_FUNCTION_API.mesh.get_number_of_orientations(
        3,"H1Legendre1")==24
    @test _MESH_FUNCTION_API.mesh.get_number_of_orientations(
        5,"H1Legendre1")==40320
    @test _MESH_FUNCTION_API.mesh.get_number_of_orientations(
        6,"GradH1Legendre1")==720
    @test _MESH_FUNCTION_API.mesh.get_number_of_orientations(
        15,"H1Legendre1")==1
    @test _MESH_FUNCTION_API.mesh.get_number_of_keys(
        5,"H1Legendre1")==8
    @test _MESH_FUNCTION_API.mesh.get_number_of_keys(
        6,"GradH1Legendre1")==6
    nonsimplex_h1=_MESH_FUNCTION_API.mesh.get_basis_functions(
        3,[0.1,0.2,0.3],"H1Legendre1",Int32[0])
    @test nonsimplex_h1[1]==1
    @test nonsimplex_h1[3]==24
    @test nonsimplex_h1[2]==quadrangle[2]
    @test _MESH_FUNCTION_API.mesh.get_basis_functions(
        3,[0.1,0.2,0.0],"GradH1Legendre1",Int32[0])[2]==
        _MESH_FUNCTION_API.mesh.get_basis_functions(
            3,[0.1,0.2,0.0],"GradLagrange1",Int32[0])[2]
    @test_throws ArgumentError _MESH_FUNCTION_API.mesh.get_basis_functions(
        7,[0.1,0.2,0.0],"H1Legendre1")
    @test_throws ArgumentError _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
        2,"Lagrange")
    @test_throws ArgumentError _MESH_FUNCTION_API.mesh.get_keys(
        2,"Lagrange")

    try
        _MESH_FUNCTION_API.initialize()
        @test_throws ArgumentError _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            2,"Lagrange")
        @test_throws ArgumentError _MESH_FUNCTION_API.mesh.get_keys(
            2,"Lagrange")

        fixture=_mesh_function_api_fixture()
        baseline=mesh_crc(fixture)
        _install_mesh_function_fixture!(fixture)

        @test _MESH_FUNCTION_API.mesh.get_basis_functions(
            12,[0.1,0.2,0.3],"GradLagrange1")[1:2:end]==
            (Int32(3),Int32(1))
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            12,"Lagrange1")==Int32[]
        @test _MESH_FUNCTION_API.mesh.get_keys(
            12,"Lagrange1")==
            (Int32[],UInt64[],Float64[])
        # Non-simplex H1 spaces are reference-only: absent type blocks stay
        # empty and never synthesize cached topology.
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            3,"H1Legendre1")==Int32[]
        @test _MESH_FUNCTION_API.mesh.get_keys(
            5,"H1Legendre1")==
            (Int32[],UInt64[],Float64[])
        @test _MESH_FUNCTION_API.mesh.get_keys_information(
            zeros(Int32,4),UInt64.(1:4),3,"H1Legendre1")==
            fill((Int32(0),Int32(1)),4)
        @test mesh_crc(_MESH_FUNCTION_API.mesh.get())==baseline
        @test _MESH_FUNCTION_API.mesh.get_keys_information(
            zeros(Int32,5),UInt64.(1:5),14,"Lagrange1")==
            fill((Int32(0),Int32(1)),5)
        @test _MESH_FUNCTION_API.mesh.get_all_edges()==
            (UInt64[],UInt64[])
        @test mesh_crc(_MESH_FUNCTION_API.mesh.get())==baseline

        nodal=_MESH_FUNCTION_API.mesh.get_keys(2,"Lagrange")
        @test nodal==(
            Int32[0,0,0],UInt64[3,1,2],
            Float64[0,3,0, 0,0,0, 2,0,0])
        @test _MESH_FUNCTION_API.mesh.get_all_edges()==
            (UInt64[],UInt64[])
        @test _MESH_FUNCTION_API.mesh.get_keys(
            2,"GradLagrange1",-1,false)==
            (Int32[0,0,0],UInt64[3,1,2],Float64[])
        nodal[1][1]=99;nodal[2][1]=99;nodal[3][1]=99
        @test _MESH_FUNCTION_API.mesh.get_keys_for_element(
            2,"H1Legendre1")[2]==UInt64[3,1,2]

        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            1,"HcurlLegendre0")==Int32[1]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            2,"H1Legendre1")==Int32[4]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            4,"CurlHcurlLegendre0")==Int32[20]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation_for_element(
                1,"HcurlLegendre0")==1
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation_for_element(
                2,"GradLagrange")==0
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation_for_element(
                3,"HcurlLegendre0")==20

        segment=_MESH_FUNCTION_API.mesh.get_keys(
            1,"HcurlLegendre0")
        @test segment==(
            Int32[1],UInt64[1],Float64[1,0,0])
        @test _MESH_FUNCTION_API.mesh.get_all_edges()==
            (UInt64[1],UInt64[2,1])
        triangle=_MESH_FUNCTION_API.mesh.get_keys(
            2,"CurlHcurlLegendre0")
        @test triangle==(
            Int32[1,1,1],UInt64[2,1,3],
            Float64[0,1.5,0, 1,0,0, 1,1.5,0])
        tetrahedron=_MESH_FUNCTION_API.mesh.get_keys(
            4,"HcurlLegendre0")
        @test tetrahedron==(
            Int32[1,1,1,1,1,1],UInt64[4,1,5,6,2,3],
            Float64[1,0,2, 1,0,0, 0,0,2,
                    0,1.5,2, 0,1.5,0, 1,1.5,0])
        @test _MESH_FUNCTION_API.mesh.get_all_edges()==(
            UInt64[1,2,3,4,5,6],
            UInt64[2,1, 3,1, 2,3, 4,2, 1,4, 3,4])
        @test mesh_crc(_MESH_FUNCTION_API.mesh.get())==baseline

        @test _MESH_FUNCTION_API.mesh.get_keys_information(
            triangle[1],triangle[2],2,"HcurlLegendre0")==
            Tuple{Int32,Int32}[(1,0),(1,0),(1,0)]
        @test _MESH_FUNCTION_API.mesh.get_keys_information(
            Int32[0,0,0],UInt64[3,1,2],2,"Lagrange")==
            Tuple{Int32,Int32}[(0,-1),(0,-1),(0,-1)]

        before=_MESH_FUNCTION_API.mesh.get_all_edges()
        for invalid in (
            ()->_MESH_FUNCTION_API.mesh.get_keys(
                2,"Lagrange2"),
            ()->_MESH_FUNCTION_API.mesh.get_keys_for_element(
                2,"GradLagrange2"),
            ()->_MESH_FUNCTION_API.mesh.get_keys(
                4,"HcurlLegendre0",-1,1),
            ()->_MESH_FUNCTION_API.mesh.get_keys(
                4,"HcurlLegendre0",0),
            ()->_MESH_FUNCTION_API.mesh.get_keys_for_element(
                3,"HcurlLegendre0",1),
            ()->_MESH_FUNCTION_API.mesh.get_keys_for_element(
                0,"HcurlLegendre0"),
            ()->_MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
                4,"HcurlLegendre0",0),
            ()->_MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
                4,"HcurlLegendre0",-1,-1,1),
            ()->_MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
                4,"HcurlLegendre0",-1,0,0),
            ()->_MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
                4,"HcurlLegendre0",-1,true,1),
            ()->_MESH_FUNCTION_API.mesh.get_keys_information(
                Int32[1],UInt64[1],4,"HcurlLegendre0"),
            ()->_MESH_FUNCTION_API.mesh.get_basis_functions(
                10,[0,0,0],"Lagrange11"),
            ()->_MESH_FUNCTION_API.mesh.get_basis_functions(
                7,[0,0,0],"H1Legendre1"),
            ()->_MESH_FUNCTION_API.mesh.get_basis_functions(
                140,[0,0,0],"Lagrange1"),
        )
            @test_throws ArgumentError invalid()
            @test _MESH_FUNCTION_API.mesh.get_all_edges()==before
            @test mesh_crc(_MESH_FUNCTION_API.mesh.get())==baseline
        end

        # A single-element query creates only that element's edges. Calling the
        # triangle afterward reuses the same global identifiers.
        _install_mesh_function_fixture!(fixture)
        single=_MESH_FUNCTION_API.mesh.get_keys_for_element(
            3,"HcurlLegendre0")
        @test single[2]==UInt64[1,2,3,4,5,6]
        @test _MESH_FUNCTION_API.mesh.get_all_edges()==(
            UInt64[1,2,3,4,5,6],
            UInt64[4,2, 2,1, 1,4, 3,4, 3,1, 3,2])
        completed_catalog=_MESH_FUNCTION_API.LAST_MESH_EDGES[]
        @test _MESH_FUNCTION_API.mesh.get_keys_for_element(
            3,"HcurlLegendre0",false)[2]==UInt64[1,2,3,4,5,6]
        @test _MESH_FUNCTION_API.LAST_MESH_EDGES[]===completed_catalog
        @test _MESH_FUNCTION_API.mesh.get_keys(
            2,"HcurlLegendre0")[2]==UInt64[5,2,6]

        # Existing caller-owned edge identifiers remain canonical DOF keys.
        _install_mesh_function_fixture!(fixture)
        _MESH_FUNCTION_API.mesh.add_edges(UInt64[99],UInt64[4,2])
        @test _MESH_FUNCTION_API.mesh.get_keys_for_element(
            3,"HcurlLegendre0",false)==(
                Int32[1,1,1,1,1,1],UInt64[99,2,3,4,5,6],Float64[])
        @test _MESH_FUNCTION_API.mesh.get_edges(UInt64[4,2])==
            (UInt64[99],Int32[-1])

        # Asking for an absent type is empty and does not synthesize topology.
        triangle_only=Mesh(fixture.coords;tris=fixture.tris)
        _install_mesh_function_fixture!(triangle_only)
        @test _MESH_FUNCTION_API.mesh.get_keys(
            4,"HcurlLegendre0")==
            (Int32[],UInt64[],Float64[])
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            4,"HcurlLegendre0")==Int32[]
        @test _MESH_FUNCTION_API.mesh.get_all_edges()==
            (UInt64[],UInt64[])

        # Task partitioning selects contiguous Gmsh blocks. The single cached
        # tetrahedron falls outside task 0 of 2 (empty, not an error) and
        # inside task 1 of 2, whose 0-based block is [0,1).
        _install_mesh_function_fixture!(fixture)
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            4,"CurlHcurlLegendre0",-1,0,1)==Int32[20]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            4,"CurlHcurlLegendre0",-1,0,2)==Int32[]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            4,"CurlHcurlLegendre0",-1,1,2)==Int32[20]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            4,"CurlHcurlLegendre0",-1,2,2)==Int32[]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            12,"Lagrange1",-1,1,2)==Int32[]

        # Five ascending segments give five zero hierarchical ranks; the
        # uneven 2/3 split pins the contiguous formula and the union
        # reproduces the complete query.
        chain_coordinates=zeros(3,6)
        chain_coordinates[1,:].=0:5
        chain_segments=Matrix{Int32}(undef,2,5)
        for segment in 1:5
            chain_segments[1,segment]=Int32(segment)
            chain_segments[2,segment]=Int32(segment+1)
        end
        _install_mesh_function_fixture!(
            Mesh(chain_coordinates;segs=chain_segments))
        chain_full=_MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            1,"HcurlLegendre0")
        @test chain_full==Int32[0,0,0,0,0]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            1,"HcurlLegendre0",-1,0,2)==Int32[0,0]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            1,"HcurlLegendre0",-1,1,2)==Int32[0,0,0]
        @test _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
            1,"HcurlLegendre0",-1,1,3)==chain_full[2:3]
        unioned=Int32[]
        for task in 0:2
            append!(unioned,
                    _MESH_FUNCTION_API.mesh.get_basis_functions_orientation(
                        1,"HcurlLegendre0",-1,task,3))
        end
        @test unioned==chain_full
        @test all(0 .<= chain_full .<= 1)
        _install_mesh_function_fixture!(fixture)

        basis=_MESH_FUNCTION_API.mesh.get_basis_functions(
            4,[0.2,0.3,0.1],"HcurlLegendre0",Int32[20])
        @test basis[1]==3
        @test basis[3]==24
        @test length(basis[2])==18
        basis[2][1]=99
        @test _MESH_FUNCTION_API.mesh.get_basis_functions(
            4,[0.2,0.3,0.1],"HcurlLegendre0",Int32[20])[2][1]!=99

        @test isempty(Docs.undocumented_names(
            Tessella.API.mesh;private=false))
        @test isempty(Test.detect_ambiguities(
            Tessella.API.mesh;recursive=true))
    finally
        _MESH_FUNCTION_API.finalize()
    end
end
