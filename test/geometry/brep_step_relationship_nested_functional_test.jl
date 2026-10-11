module StepRelationshipNestedChecks
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

function corners(mesh,owner)
    return Set(Tuple(mesh.coords[:,row]) for row in axes(mesh.coords,2)
        if mesh.entity_data.node_entities[row]==(2,Int32(owner)))
end
expected(lo,width)=Set((lo[1]+x,lo[2]+y,lo[3]+z) for x in (0.,width),y in (0.,width),z in (0.,width))

function append_entities(source,extra)
    tail=findlast("ENDSEC;",source)
    return source[1:first(tail)-1]*extra*source[first(tail):end]
end

@testset "actual nested product usages expand literal root-to-leaf paths" begin
    cases=(("brep_relationship_nested_mm.step",((10.,5.,0.),(10.,7.,0.)),1.),
        ("brep_relationship_nested_mixed_units.step",((9990.,2029.,0.),(9990.,2058.,0.)),1.))
    for (file,origins,width) in cases
        mesh=import_step(joinpath(FIXTURES,file);fill=false)
        @test mesh isa E.MixedMesh && validate(mesh).ok
        @test [mesh.entity_data.entities[(3,owner)].bbox for owner in 1:2]==
            [(origin...,origin.+1...) for origin in origins]
        @test size(mesh.coords)==(3,16) && sum(size(block.nodes,2) for block in mesh.blocks)==24
        @test all(isapprox(material_volumes(mesh,2)[Int32(owner)],width;atol=2e-10) for owner in 1:2)
        @test allunique(mesh.entity_data.external_node_tags) &&
            allunique(reduce(vcat,mesh.entity_data.external_element_tags))
        for owner in 1:2
            @test corners(mesh,owner)==expected(origins[owner],1.)
            @test mesh.entity_data.entities[(3,owner)].boundaries==[Int32(owner)]
        end
        @test occursin("PRODUCT #1000 Assembly/NAUO #2020 nested/PRODUCT #2000/"*
            "NAUO #1020 left/PRODUCT #7",mesh.physical_names[(3,1)])
        @test occursin("NAUO #1030 right/PRODUCT #7",mesh.physical_names[(3,2)])
        @test all(ncodeunits(mesh.physical_names[(3,owner)])<=
            E.MSH_PHYSICAL_NAME_MAX_BYTES-ncodeunits(" boundary") for owner in 1:2)
    end
end

@testset "shared intermediate and leaf definitions expand through every usage path" begin
    source=read(joinpath(FIXTURES,"brep_relationship_nested_mm.step"),String)
    diamond=replace(append_entities(source,
        "#2040 = NEXT_ASSEMBLY_USAGE_OCCURRENCE('4','again','',#1002,#2002,\$);\n"*
        "#2041 = PRODUCT_DEFINITION_SHAPE('Placement','again occurrence',#2040);\n"*
        "#2042 = ITEM_DEFINED_TRANSFORMATION('','',#2007,#2100);\n"*
        "#2043 = ( REPRESENTATION_RELATIONSHIP('','',#2005,#1005)\n"*
        "  REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#2042)\n"*
        "  SHAPE_REPRESENTATION_RELATIONSHIP() );\n"*
        "#2044 = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#2043,#2041);\n"*
        "#2100 = AXIS2_PLACEMENT_3D('',#2101,#13,#14);\n"*
        "#2101 = CARTESIAN_POINT('',(20.,2.,0.));\n"),
        "SHAPE_REPRESENTATION('Assembly',(#1006,#2010),#1045)"=>
            "SHAPE_REPRESENTATION('Assembly',(#1006,#2010,#2100),#1045)")
    mktempdir() do directory
        path=joinpath(directory,"diamond.step");write(path,diamond)
        mesh=import_step(path;fill=false)
        @test mesh isa E.MixedMesh && validate(mesh).ok
        @test Set(mesh.entity_data.entities[(3,owner)].bbox for owner in 1:4)==
            Set((x,y,0.,x+1,y+1,1.) for (x,y) in ((10.,5.),(10.,7.),(23.,2.),(25.,2.)))
        @test size(mesh.coords)==(3,32) &&
            sum(size(block.nodes,2) for block in mesh.blocks)==48
        @test all(isapprox(value,1.;atol=2e-13) for value in values(material_volumes(mesh,2)))
        @test allunique(mesh.entity_data.external_node_tags) &&
            allunique(reduce(vcat,mesh.entity_data.external_element_tags))
        @test count(name->occursin("NAUO #2020 nested",name),
            (mesh.physical_names[(3,owner)] for owner in 1:4))==2
        @test count(name->occursin("NAUO #2040 again",name),
            (mesh.physical_names[(3,owner)] for owner in 1:4))==2
        @test all(name->occursin("PRODUCT #2000",name) &&
            occursin("PRODUCT #7",name) && occursin("MANIFOLD_SOLID_BREP #15",name),
            (mesh.physical_names[(3,owner)] for owner in 1:4))
    end
end

@testset "one actual nested fill retains positive volumes and MSH4 path identities" begin
    filled=import_step(joinpath(FIXTURES,"brep_relationship_nested_mm.step"))
    @test validate(filled).ok && Set(keys(material_volumes(filled,4)))==Set(Int32[1,2])
    @test all(isapprox(volume,1.;atol=2e-12) for volume in values(material_volumes(filled,4)))
    @test all(det(hcat((filled.coords[:,cell[i]]-filled.coords[:,cell[1]] for i in 2:4)...))>0
        for block in filled.blocks if block.msh==4 for cell in eachcol(block.nodes))
    mktempdir() do directory
        for binary in (false,true)
            path=joinpath(directory,"nested-$binary.msh")
            E.write_mixed_msh(path,filled;version=4.1,binary)
            restored=E.read_mixed_msh(path)
            @test validate(restored).ok && restored.physical_names==filled.physical_names
            @test actual_records(restored)==actual_records(filled)
        end
    end
end

@testset "malformed nested graphs, cycles and orphan bindings reject explicitly" begin
    source=read(joinpath(FIXTURES,"brep_relationship_nested_mm.step"),String)
    cycle=replace(source,
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('1','left','',#2002,#5,\$)"=>
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('1','left','',#5,#2002,\$)")
    self_root=replace(source,
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','nested','',#1002,#2002,\$)"=>
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','nested','',#2002,#2002,\$)")
    orphan_pd=append_entities(source,"#1100 = PRODUCT_DEFINITION('design','',#6,#9);\n")
    nonmember=replace(source,"ITEM_DEFINED_TRANSFORMATION('','',#2007,#2010)"=>
        "ITEM_DEFINED_TRANSFORMATION('','',#2007,#1010)")
    reversed=replace(source,"REPRESENTATION_RELATIONSHIP('','',#2005,#1005)"=>
        "REPRESENTATION_RELATIONSHIP('','',#1005,#2005)")
    mixed_parent=replace(source,"SHAPE_REPRESENTATION('Subassembly',(#2007,#1010,#1013),#2045)"=>
        "SHAPE_REPRESENTATION('Subassembly',(#2007,#1010,#1013,#15),#2045)")
    empty_leaf=replace(source,"ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#15),#345)"=>
        "ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11),#345)")
    dropped_cdsr=replace(source,"#1024 = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#1023,#1021);\n"=>"")
    unreachable_usage=replace(append_entities(source,
        "#2040 = NEXT_ASSEMBLY_USAGE_OCCURRENCE('4','again','',#5,#5,\$);\n"),
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','nested','',#1002,#2002,\$)"=>
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','nested','',#1002,#5,\$)")
    invalid=[cycle,self_root,orphan_pd,nonmember,reversed,mixed_parent,empty_leaf,
        dropped_cdsr,unreachable_usage]
    mktempdir() do directory
        path=joinpath(directory,"invalid.step")
        for text in invalid
            write(path,text)
            @test_throws ArgumentError import_step(path;fill=false)
            @test read(path,String)==text
        end
    end
end

@testset "ordinary direct assemblies and mapped occurrences are unchanged" begin
    direct=import_step(joinpath(FIXTURES,"brep_relationship_boxes.step");fill=false)
    @test validate(direct).ok && size(direct.coords)==(3,16)
    @test [direct.entity_data.entities[(3,owner)].bbox for owner in 1:2]==
        [(3.,0.,0.,4.,1.,1.),(5.,0.,0.,6.,1.,1.)]
    @test occursin("PRODUCT #1000 Assembly/NAUO #1020 left/PRODUCT #7",
        direct.physical_names[(3,1)])
    mapped=import_step(joinpath(FIXTURES,"brep_mapped_boxes.step");fill=false)
    @test validate(mapped).ok && mapped.entity_data.entities[(3,2)].bbox==(5.,0.,0.,6.,1.,1.)
end
end
