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
# Overlapping bounding boxes admit genuinely disjoint solids; actual contact
# still rejects. Full geometric-certificate coverage lives in
# test/geometry/brep_assembly_functional_test.jl.
let source=read(joinpath(@__DIR__,"..","fixtures","brep_disjoint_overlap.step"),String)
    mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_disjoint_overlap.step");fill=false)
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok &&
          Set(keys(mesh.entity_data.entities))==Set([(2,1),(3,1),(2,2),(3,2)])
    mktempdir() do directory
        path=joinpath(directory,"coincident.step")
        write(path,replace(source,"FACETED_BREP('box',#140)"=>"FACETED_BREP('box',#40)"))
        @test_throws ArgumentError import_step(path;fill=false)
    end
end
# Conforming coincident facets weld into one shared Surface bound by both
# Regions; the interface cell emits once and rim nodes deduplicate.
let mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_contact_twin.step");fill=false)
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok
    data=mesh.entity_data
    @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)]) &&
          sort(data.entities[(3,1)].boundaries)==Int32[1,3] &&
          sort(data.entities[(3,2)].boundaries)==Int32[2,3] &&
          allunique(eachcol(mesh.coords))
    counts=Dict{Int32,Int}()
    for (bi,block) in enumerate(mesh.blocks)
        for owner in data.block_entities[bi]
            counts[owner]=get(counts,owner,0)+1
        end
    end
    @test counts==Dict(Int32(1)=>3,Int32(2)=>3,Int32(3)=>1)
end
end
