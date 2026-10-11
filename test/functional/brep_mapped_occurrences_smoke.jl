@testset "STEP mapped occurrences retain real placements" begin
    mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_mapped_boxes.step");fill=false)
    @test mesh.entity_data.entities[(3,1)].bbox==(3.,0.,0.,4.,1.,1.) &&
        mesh.entity_data.entities[(3,2)].bbox==(5.,0.,0.,6.,1.,1.)
    @test validate(mesh).ok && size(mesh.coords,2)==16
end
