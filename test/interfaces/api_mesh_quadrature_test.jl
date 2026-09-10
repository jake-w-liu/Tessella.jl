using Test
using Tessella
using Tessella.MeshTypes: mesh_crc

const _QUADRATURE_API=Tessella.API

function _quadrature_api_fixture()
    return Mesh(
        Float64[0 1 0 0;
                0 0 1 0;
                0 0 0 1];
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1))
end

@testset "reference quadrature through API" begin
    _QUADRATURE_API.finalize()
    expected=_QUADRATURE_API.mesh.get_integration_points(2,"Gauss2")
    @test length(expected[1])==9
    @test length(expected[2])==3
    @test isapprox(sum(expected[2]),0.5;atol=2e-15,rtol=0)
    quadrangle=_QUADRATURE_API.mesh.get_integration_points(3,"Gauss2")
    @test length(quadrangle[2])==7
    @test isapprox(sum(quadrangle[2]),4.0;atol=2e-14,rtol=0)

    try
        _QUADRATURE_API.initialize()
        fixture=_quadrature_api_fixture()
        lock(_QUADRATURE_API.STATE_LOCK) do
            _QUADRATURE_API._replace_mesh_cache_locked!(
                _QUADRATURE_API._copy_mesh(fixture))
        end
        baseline=mesh_crc(_QUADRATURE_API.mesh.get())
        tetrahedron=_QUADRATURE_API.mesh.get_integration_points(
            11,"CompositeGauss8")
        @test length(tetrahedron[2])==216
        @test isapprox(sum(tetrahedron[2]),1/6;atol=2e-15,rtol=0)
        tetrahedron[1][1]=99
        tetrahedron[2][1]=99
        @test _QUADRATURE_API.mesh.get_integration_points(
            11,"CompositeGauss8")[1][1]!=99
        pyramid=_QUADRATURE_API.mesh.get_integration_points(
            14,"CompositeGauss8")
        @test length(pyramid[2])==125
        @test isapprox(sum(pyramid[2]),4/3;atol=2e-14,rtol=0)
        @test mesh_crc(_QUADRATURE_API.mesh.get())==baseline
    finally
        _QUADRATURE_API.finalize()
    end

    @test _QUADRATURE_API.mesh.get_integration_points(15,"Gauss99")==
          (Float64[0,0,0],Float64[1])
    high_triangle=_QUADRATURE_API.mesh.get_integration_points(2,"Gauss20")
    @test length(high_triangle[2])==79
    @test isapprox(sum(high_triangle[2]),0.5;atol=5e-13,rtol=0)
    @test_throws ArgumentError _QUADRATURE_API.mesh.get_integration_points(
        2,"Gauss255")
    message=try
        _QUADRATURE_API.mesh.get_integration_points(140,"Gauss2")
        ""
    catch err
        sprint(showerror,err)
    end
    @test occursin("API.mesh.get_integration_points",message)
    @test occursin("no integration rules",message)
end
