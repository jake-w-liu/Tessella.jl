module StepUnitsConversionChecks
using Test,Tessella,LinearAlgebra
const E=Tessella.Elements
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function material_volumes(mesh,kind)
    volumes=Dict{Int32,Float64}()
    for (block,owners) in zip(mesh.blocks,mesh.entity_data.block_entities)
        block.msh==kind || continue
        for (cell,owner) in zip(eachcol(block.nodes),owners)
            a,b,c=(mesh.coords[:,row] for row in cell[1:3])
            volume=kind==2 ? dot(a,cross(b,c))/6 :
                det(hcat(b-a,c-a,mesh.coords[:,cell[4]]-a))/6
            volumes[owner]=get(volumes,owner,0.)+volume
        end
    end
    return volumes
end

function actual_records(mesh)
    data=mesh.entity_data
    nodes=Dict(uid=>(data.node_entities[row],Tuple(mesh.coords[:,row]))
        for (row,uid) in enumerate(data.external_node_tags))
    cells=Dict(uid=>(Int(block.msh),data.block_entities[index][column],
            Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (index,block) in enumerate(mesh.blocks)
        for (column,uid) in enumerate(data.external_element_tags[index]))
    return nodes,cells
end

function append_entities(source,extra)
    tail=findlast("ENDSEC;",source)
    return source[1:first(tail)-1]*extra*source[first(tail):end]
end

@testset "ordinary roots convert declared SI and conversion-based length units" begin
    cases=(("brep_occ_box.step",(0.,0.,0.,1.,1.,1.)),
        ("brep_unit_metre.step",(0.,0.,0.,1000.,1000.,1000.)),
        ("brep_unit_inch.step",(0.,0.,0.,25.4,25.4,25.4)),
        ("brep_unit_yard.step",(0.,0.,0.,914.4,914.4,914.4)),
        ("brep_unit_angular_prefix.step",(0.,0.,0.,1.,1.,1.)),
        ("brep_unit_uncertainty_alias.step",(0.,0.,0.,1.,1.,1.)))
    for (file,bounds) in cases
        mesh=import_step(joinpath(FIXTURES,file);fill=false)
        @test mesh isa E.Mesh && validate(mesh).ok
        @test vec(minimum(mesh.coords,dims=2))≈collect(bounds[1:3]) &&
            vec(maximum(mesh.coords,dims=2))≈collect(bounds[4:6])
    end
end

@testset "one actual metre fill retains positive oriented tets" begin
    filled=import_step(joinpath(FIXTURES,"brep_unit_metre.step"))
    @test filled isa E.Mesh && validate(filled).ok
    @test size(filled.tets,2)>0 && size(filled.coords,2)>8
    @test all(det(hcat((filled.coords[:,cell[i]]-filled.coords[:,cell[1]] for i in 2:4)...))>0
        for cell in eachcol(filled.tets))
    volume=sum(det(hcat((filled.coords[:,cell[i]]-filled.coords[:,cell[1]] for i in 2:4)...))/6
        for cell in eachcol(filled.tets))
    @test isapprox(volume,1e9;atol=1e-4)
end

@testset "mapped occurrences convert the shared representation context units" begin
    cases=(("brep_mapped_boxes.step",((3.,0.,0.,4.,1.,1.),(5.,0.,0.,6.,1.,1.))),
        ("brep_mapped_metre.step",((3000.,0.,0.,4000.,1000.,1000.),(5000.,0.,0.,6000.,1000.,1000.))),
        ("brep_mapped_inch.step",((76.2,0.,0.,101.6,25.4,25.4),(127.,0.,0.,152.4,25.4,25.4))))
    for (file,bounds) in cases
        mesh=import_step(joinpath(FIXTURES,file);fill=false)
        @test mesh isa E.MixedMesh && validate(mesh).ok
        @test all(collect(mesh.entity_data.entities[(3,owner)].bbox)≈
            collect(bounds[owner]) for owner in 1:2)
        @test allunique(mesh.entity_data.external_node_tags) &&
            allunique(reduce(vcat,mesh.entity_data.external_element_tags))
    end
end

@testset "conversion-based relationship units reproduce the declared SI result" begin
    metres=import_step(joinpath(FIXTURES,"brep_relationship_metres.step");fill=false)
    converted=import_step(joinpath(FIXTURES,"brep_relationship_metres_conversion.step");fill=false)
    @test converted isa E.MixedMesh && validate(converted).ok
    @test converted.coords==metres.coords
    @test converted.physical_names==metres.physical_names
    @test [converted.entity_data.entities[(3,owner)].bbox for owner in 1:2]==
        [metres.entity_data.entities[(3,owner)].bbox for owner in 1:2]
    @test all(isapprox(value,1e9;atol=1e-4) for value in values(material_volumes(converted,2)))
    filled=import_step(joinpath(FIXTURES,"brep_relationship_metres_conversion.step"))
    @test validate(filled).ok &&
        all(isapprox(value,1e9;atol=1e-4) for value in values(material_volumes(filled,4)))
    mktempdir() do directory
        for binary in (false,true)
            path=joinpath(directory,"converted-$binary.msh")
            E.write_mixed_msh(path,filled;version=4.1,binary)
            restored=E.read_mixed_msh(path)
            @test validate(restored).ok && restored.physical_names==filled.physical_names
            @test actual_records(restored)==actual_records(filled)
        end
    end
end

@testset "conflicting member contexts and malformed conversions reject atomically" begin
    for file in ("brep_unit_bad_measure.step","brep_unit_bad_category.step",
        "brep_unit_negative.step","brep_unit_cycle.step")
        @test_throws ArgumentError import_step(joinpath(FIXTURES,file);fill=false)
    end
    source=read(joinpath(FIXTURES,"brep_occ_box.step"),String)
    conflict=append_entities(source,
        "#800 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#15),#801);\n"*
        "#801 = ( GEOMETRIC_REPRESENTATION_CONTEXT(3)\n"*
        "GLOBAL_UNIT_ASSIGNED_CONTEXT((#802,#347,#348))\n"*
        "REPRESENTATION_CONTEXT('other','other') );\n"*
        "#802 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(\$,.METRE.) );\n")
    bad_alias=append_entities(replace(source,
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#346"=>
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#901"),
        "#901 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(\$,.METRE.) );\n")
    mktempdir() do directory
        path=joinpath(directory,"invalid.step")
        for text in (conflict,bad_alias)
            write(path,text)
            @test_throws ArgumentError import_step(path;fill=false)
            @test read(path,String)==text
        end
    end
end
end
