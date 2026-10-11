@testset "IGES143 planar Face boundaries" begin
    surface=Tessella.import_iges(joinpath(@__DIR__,"..","fixtures","brep_bounded143_through_hole.igs");fill=false)
    @test Tessella.validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==16 && Tessella.MeshTypes.ntris(surface)==32
    signed_volume=sum(LinearAlgebra.dot(surface.coords[:,cell[1]],
        LinearAlgebra.cross(surface.coords[:,cell[2]],surface.coords[:,cell[3]]))/6
        for cell in eachcol(surface.tris);init=0.0)
    @test signed_volume≈.75 atol=1e-12
end
