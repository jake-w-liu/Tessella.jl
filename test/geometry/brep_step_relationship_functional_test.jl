module StepRelationshipAssemblyChecks
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

@testset "actual direct STEP product usages instantiate separate owned solids" begin
    for (file,offsets) in (("brep_relationship_boxes.step",(3.,5.)),
            ("brep_relationship_shifted_boxes.step",(2.,4.)))
        mesh=import_step(joinpath(FIXTURES,file);fill=false)
        @test mesh isa E.MixedMesh && validate(mesh).ok
        @test [mesh.entity_data.entities[(3,owner)].bbox for owner in 1:2]==
            [(offset,0.,0.,offset+1,1.,1.) for offset in offsets]
        @test size(mesh.coords)==(3,16) && sum(size(block.nodes,2) for block in mesh.blocks)==24
        @test all(isapprox(material_volumes(mesh,2)[Int32(owner)],1.;atol=2e-13) for owner in 1:2)
        @test allunique(mesh.entity_data.external_node_tags) &&
            allunique(reduce(vcat,mesh.entity_data.external_element_tags))
        for owner in 1:2
            actual=Set(Tuple(mesh.coords[:,row]) for row in axes(mesh.coords,2)
                if mesh.entity_data.node_entities[row]==(2,Int32(owner)))
            expected=Set((x+offsets[owner],y,z) for x in (0.,1.),y in (0.,1.),z in (0.,1.))
            @test actual==expected
        end
        @test occursin("PRODUCT #1000 Assembly/NAUO #1020 left/PRODUCT #7",mesh.physical_names[(3,1)]) &&
            occursin("NAUO #1030 right/PRODUCT #7",mesh.physical_names[(3,2)])
    end
end

@testset "one actual assembly fill retains positive volume and MSH4 occurrence identities" begin
    filled=import_step(joinpath(FIXTURES,"brep_relationship_boxes.step"))
    @test validate(filled).ok && Set(keys(material_volumes(filled,4)))==Set(Int32[1,2])
    @test all(isapprox(volume,1.;atol=2e-12) for volume in values(material_volumes(filled,4)))
    @test all(det(hcat((filled.coords[:,cell[i]]-filled.coords[:,cell[1]] for i in 2:4)...))>0
        for block in filled.blocks if block.msh==4 for cell in eachcol(block.nodes))
    mktempdir() do directory
        for binary in (false,true)
            path=joinpath(directory,"placed-$binary.msh")
            E.write_mixed_msh(path,filled;version=4.1,binary)
            restored=E.read_mixed_msh(path)
            @test validate(restored).ok && restored.physical_names==filled.physical_names
            @test actual_records(restored)==actual_records(filled)
        end
    end
    # The separate mapped-item producer is still reached without relationship entities.
    mapped=import_step(joinpath(FIXTURES,"brep_mapped_boxes.step");fill=false)
    @test validate(mapped).ok && mapped.entity_data.entities[(3,2)].bbox==(5.,0.,0.,6.,1.,1.)
end

@testset "incomplete or unsupported assembly graphs reject before source-solid fallback" begin
    source=read(joinpath(FIXTURES,"brep_relationship_boxes.step"),String)
    root_context="#1045 = ( GEOMETRIC_REPRESENTATION_CONTEXT(3)"
    invalid=[
        replace(source,"NEXT_ASSEMBLY_USAGE_OCCURRENCE('1','left','',#1002,#5,\$)"=>
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE('1','left','',#1002,#1002,\$)"),
        replace(source,"NEXT_ASSEMBLY_USAGE_OCCURRENCE('2','right','',#1002,#5,\$)"=>
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE('2','right','',#5,#1002,\$)"),
        replace(source,"REPRESENTATION_RELATIONSHIP('','',#10,#1005)"=>
            "REPRESENTATION_RELATIONSHIP('','',#1005,#10)"),
        replace(source,"ITEM_DEFINED_TRANSFORMATION('','',#11,#1010)"=>
            "ITEM_DEFINED_TRANSFORMATION('','',#11,#9999)"),
        append_entities(replace(source,"ITEM_DEFINED_TRANSFORMATION('','',#11,#1010)"=>
            "ITEM_DEFINED_TRANSFORMATION('','',#11,#1100)"),
            "#1100 = AXIS2_PLACEMENT_3D('',#1008,#13,#14);\n"),
        replace(source,"SHAPE_REPRESENTATION('Assembly',(#1006,#1010,#1013),#1045)"=>
            "SHAPE_REPRESENTATION('Assembly',(#1006,#1010,#1013),#345)"),
        replace(source,root_context=>"#1045 = ( GEOMETRIC_REPRESENTATION_CONTEXT(2)"),
        replace(source,"SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(.CENTI.,.METRE.)"),
        append_entities(replace(source,"((#346,#347,#348)) REPRESENTATION_CONTEXT('Assembly Context'"=>
            "((#1100,#347,#348)) REPRESENTATION_CONTEXT('Assembly Context'"),
            "#1100 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.CENTI.,.METRE.) );\n"),
        replace(source,"CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#1023,#1021)"=>
            "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#1023,#9999)"),
        append_entities(source,"#1100 = ITEM_DEFINED_TRANSFORMATION('','',#11,#1010);\n"),
        append_entities(source,"#1100 = MANIFOLD_SOLID_BREP('unplaced',#16);\n"),
        append_entities(source,"#1100 = SHAPE_REPRESENTATION('unselected',(#1010),#1045);\n"),
        append_entities(source,"#1100 = SHAPE_DEFINITION_REPRESENTATION(#4,#10);\n"),
        replace(source,"CARTESIAN_POINT('',(6.,0.,0.))"=>"CARTESIAN_POINT('',(4.,0.,0.))")]
    mktempdir() do directory
        path=joinpath(directory,"invalid.step")
        for text in invalid
            write(path,text)
            @test_throws ArgumentError import_step(path;fill=false)
            @test read(path,String)==text
        end
    end
end
end
