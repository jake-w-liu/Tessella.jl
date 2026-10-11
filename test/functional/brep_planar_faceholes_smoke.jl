@testset "planar Face holes import" begin
    surface=Tessella.import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_planar_through_hole.step");fill=false)
    @test Tessella.validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==16 && Tessella.MeshTypes.ntris(surface)==32
    signed_volume=sum(LinearAlgebra.dot(surface.coords[:,cell[1]],
        LinearAlgebra.cross(surface.coords[:,cell[2]],surface.coords[:,cell[3]]))/6
        for cell in eachcol(surface.tris);init=0.0)
    @test signed_volume≈.75 atol=1e-12
    # Both material caps retain the open square passage.
    cap_areas=zeros(2)
    for (index,height) in enumerate((0.,1.)),cell in eachcol(surface.tris)
        all(row->surface.coords[3,row]==height,cell) || continue
        cap_areas[index]+=LinearAlgebra.norm(LinearAlgebra.cross(
            surface.coords[:,cell[2]]-surface.coords[:,cell[1]],
            surface.coords[:,cell[3]]-surface.coords[:,cell[1]]))/2
    end
    @test all(area->isapprox(area,.75;atol=1e-12),cap_areas)
end
