module BRepMixedRelationshipSmokeChecks
using Test, Tessella
let mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_relationship_mixed_roots.step");fill=false)
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok
    @test [mesh.entity_data.entities[(3,owner)].bbox for owner in 1:5]==
        [(10.,0.,0.,11.,1.,1.),(2.,0.,0.,3.,1.,1.),(5.,0.,0.,6.,2.,1.),
         (10.,3.,5.,11.,4.,6.),(20.,0.,0.,21.,1.,1.)]
    @test size(mesh.coords)==(3,40) && allunique(mesh.entity_data.external_node_tags) &&
        occursin("MAPPED_ITEM #3050",mesh.physical_names[(3,4)])
end
end
