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
# Coincident vertex-only contact is conforming: admitted unwelded — both
# solids keep their own coincident boundary node (OCC import semantics).
let mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_contact_vertex.step");fill=false)
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok &&
          Set(keys(mesh.entity_data.entities))==Set([(2,1),(3,1),(2,2),(3,2)]) &&
          !allunique(eachcol(mesh.coords))
end
# A solid strictly enclosed by another is a nested Region: under fill the
# enclosed boundary welds into a shared Surface bound by both Regions and
# the joint fill tiles each region exactly.
let mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_nested.step");fill=true)
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok
    data=mesh.entity_data
    @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)]) &&
          sort(data.entities[(3,1)].boundaries)==Int32[1,3] &&
          sort(data.entities[(3,2)].boundaries)==Int32[2,3] &&
          allunique(eachcol(mesh.coords))
end
# Transversally crossing boundaries fragment into three Regions welded at
# two shared interface Surfaces — A only, B only and the intersection.
let mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_crossing_boxes.step");fill=false)
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok
    data=mesh.entity_data
    @test Set(keys(data.entities))==Set([(2,t) for t in 1:5] ∪ [(3,t) for t in 1:3]) &&
          sort(data.entities[(3,3)].boundaries)==Int32[3,4,5] &&
          allunique(eachcol(mesh.coords))
end
# Coplanar partial-face overlap imprints: box B moved onto A's x=1 face
# welds the exact overlap patch into a shared Surface bound by both Regions.
let source=read(joinpath(@__DIR__,"..","fixtures","brep_two_boxes.step"),String)
    partial=foldl(replace,["(3.,0.,0.)"=>"(1.,0.5,0.)","(4.,0.,0.)"=>"(2.,0.5,0.)",
                "(4.,2.,0.)"=>"(2.,1.5,0.)","(3.,2.,0.)"=>"(1.,1.5,0.)",
                "(3.,0.,1.)"=>"(1.,0.5,1.)","(4.,0.,1.)"=>"(2.,0.5,1.)",
                "(4.,2.,1.)"=>"(2.,1.5,1.)","(3.,2.,1.)"=>"(1.,1.5,1.)"];init=source)
    mktempdir() do directory
        path=joinpath(directory,"imprint.step");write(path,partial)
        mesh=import_step(path;fill=false)
        @test mesh isa Tessella.MixedMesh && validate(mesh).ok
        data=mesh.entity_data
        @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)]) &&
              data.entities[(2,3)].bbox==(1.0,0.5,0.0,1.0,1.0,1.0) &&
              allunique(eachcol(mesh.coords))
    end
end
end
