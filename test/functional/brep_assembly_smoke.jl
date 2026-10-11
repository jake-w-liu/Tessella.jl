module BRepAssemblySmokeChecks
using Test, Tessella
let mesh=import_iges(joinpath(@__DIR__,"..","fixtures","brep_two_blocks.iges");fill=false)
    @test mesh isa Tessella.MixedMesh
    @test validate(mesh).ok
    @test Set(keys(mesh.entity_data.entities))==Set([(2,1),(3,1),(2,2),(3,2)]) &&
          all(mesh.entity_data.entities[(3,owner)].boundaries==[owner] for owner in 1:2)
    classified=mesh.entity_data.node_entities
    first_nodes=Set(mesh.entity_data.external_node_tags[row] for row in eachindex(classified) if classified[row]==(2,1))
    second_nodes=Set(mesh.entity_data.external_node_tags[row] for row in eachindex(classified) if classified[row]==(2,2))
    @test allunique(mesh.entity_data.external_node_tags) && isempty(intersect(first_nodes,second_nodes)) &&
          length(first_nodes)==8 && length(second_nodes)==8 && all(owner->owner[1]==2,classified)
end
end
