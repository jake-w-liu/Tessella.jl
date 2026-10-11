@testset "explicit190 plane UV chart" begin
    surface=Tessella.import_iges(joinpath(@__DIR__,"..","fixtures","brep_plane190_uv_rotated_box.igs");fill=false)
    @test Tessella.validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==8 && Tessella.MeshTypes.ntris(surface)==12
    volume=sum(LinearAlgebra.dot(surface.coords[:,cell[1]],
        LinearAlgebra.cross(surface.coords[:,cell[2]],surface.coords[:,cell[3]]))/6
        for cell in eachcol(surface.tris);init=0.0)
    @test volume≈1. atol=1e-12
end
