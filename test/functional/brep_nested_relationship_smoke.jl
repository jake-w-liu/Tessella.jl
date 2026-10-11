module BRepNestedRelationshipSmokeChecks
using Test, Tessella
let mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_relationship_nested_mm.step");fill=false)
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok
    @test [mesh.entity_data.entities[(3,owner)].bbox for owner in 1:2]==
        [(10.,5.,0.,11.,6.,1.),(10.,7.,0.,11.,8.,1.)]
    @test size(mesh.coords)==(3,16) && allunique(mesh.entity_data.external_node_tags) &&
        occursin("NAUO #2020 nested/PRODUCT #2000",mesh.physical_names[(3,1)])
end
end
