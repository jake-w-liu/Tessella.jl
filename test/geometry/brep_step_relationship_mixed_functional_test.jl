module StepRelationshipMixedChecks
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

@testset "mixed roots expand literal solids usages and mapped items" begin
    mesh=import_step(joinpath(FIXTURES,"brep_relationship_mixed_roots.step");fill=false)
    @test mesh isa E.MixedMesh && validate(mesh).ok
    @test [mesh.entity_data.entities[(3,owner)].bbox for owner in 1:5]==
        [(10.,0.,0.,11.,1.,1.),(2.,0.,0.,3.,1.,1.),(5.,0.,0.,6.,2.,1.),
         (10.,3.,5.,11.,4.,6.),(20.,0.,0.,21.,1.,1.)]
    @test size(mesh.coords)==(3,40) &&
        sum(size(block.nodes,2) for block in mesh.blocks)==60
    @test [material_volumes(mesh,2)[Int32(owner)] for owner in 1:5]≈[1.,1.,2.,1.,1.]
    @test allunique(mesh.entity_data.external_node_tags) &&
        allunique(reduce(vcat,mesh.entity_data.external_element_tags))
    names=[mesh.physical_names[(3,owner)] for owner in 1:5]
    @test names[1]=="STEP PRODUCT #2000 Assembly/FACETED_BREP #241"
    @test names[2]=="STEP PRODUCT #2000 Assembly/NAUO #2020 part/"*
        "PRODUCT #1000 Part/FACETED_BREP #41"
    @test names[3]=="STEP PRODUCT #2000 Assembly/NAUO #2020 part/"*
        "PRODUCT #1000 Part/FACETED_BREP #141"
    @test names[4]=="STEP PRODUCT #2000 Assembly/NAUO #3030 mapped/"*
        "PRODUCT #3000 MappedPart/MAPPED_ITEM #3050/FACETED_BREP #341"
    @test names[5]=="STEP PRODUCT #5000 Loose/FACETED_BREP #441"
    @test all(name->ncodeunits(name)<=
        E.MSH_PHYSICAL_NAME_MAX_BYTES-ncodeunits(" boundary"),names)
end

@testset "shared representation maps emit separate literal occurrences" begin
    source=read(joinpath(FIXTURES,"brep_relationship_mixed_roots.step"),String)
    shared=replace(append_entities(source,
        "#3061 = AXIS2_PLACEMENT_3D('',#3065,#13,#14);\n"*
        "#3065 = CARTESIAN_POINT('',(18.,0.,0.));\n"*
        "#3500 = MAPPED_ITEM('second',#3055,#3061);\n"),
        "SHAPE_REPRESENTATION('MappedPart',(#3006,#3050,#3060),#3045)"=>
            "SHAPE_REPRESENTATION('MappedPart',(#3006,#3050,#3060,#3500,#3061),#3045)")
    mktempdir() do directory
        path=joinpath(directory,"shared.step");write(path,shared)
        mesh=import_step(path;fill=false)
        @test mesh isa E.MixedMesh && validate(mesh).ok
        @test Set(mesh.entity_data.entities[(3,owner)].bbox for owner in 1:6)==
            Set([(10.,0.,0.,11.,1.,1.),(2.,0.,0.,3.,1.,1.),(5.,0.,0.,6.,2.,1.),
                (10.,3.,5.,11.,4.,6.),(18.,3.,0.,19.,4.,1.),(20.,0.,0.,21.,1.,1.)])
        @test [material_volumes(mesh,2)[Int32(owner)] for owner in 1:6]≈
            [1.,1.,2.,1.,1.,1.]
        @test count(name->occursin("MAPPED_ITEM #3050",name),
            (mesh.physical_names[(3,owner)] for owner in 1:6))==1
        @test count(name->occursin("MAPPED_ITEM #3500",name),
            (mesh.physical_names[(3,owner)] for owner in 1:6))==1
    end
end

@testset "Cartesian mapped targets inside product representations" begin
    source=read(joinpath(FIXTURES,"brep_relationship_mixed_roots.step"),String)
    cartesian=replace(replace(source,
        "#3060 = AXIS2_PLACEMENT_3D('',#3062,#3063,#3064);"=>
        "#3060 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',#3063,#3064,#3062,1.,#13);"),
        "SHAPE_REPRESENTATION('MappedPart',(#3006,#3050,#3060),#3045)"=>
            "SHAPE_REPRESENTATION('MappedPart',(#3006,#3050),#3045)")
    mktempdir() do directory
        path=joinpath(directory,"cartesian.step");write(path,cartesian)
        mesh=import_step(path;fill=false)
        @test mesh isa E.MixedMesh && validate(mesh).ok
        @test mesh.entity_data.entities[(3,4)].bbox==(14.,8.,0.,15.,9.,1.)
        @test material_volumes(mesh,2)[Int32(4)]≈1.
    end
end

@testset "one actual mixed fill retains positive volumes and MSH4 path identities" begin
    filled=import_step(joinpath(FIXTURES,"brep_relationship_mixed_roots.step"))
    @test validate(filled).ok && Set(keys(material_volumes(filled,4)))==Set(Int32[1,2,3,4,5])
    @test [material_volumes(filled,4)[Int32(owner)] for owner in 1:5]≈[1.,1.,2.,1.,1.]
    @test all(det(hcat((filled.coords[:,cell[i]]-filled.coords[:,cell[1]] for i in 2:4)...))>0
        for block in filled.blocks if block.msh==4 for cell in eachcol(block.nodes))
    mktempdir() do directory
        for binary in (false,true)
            path=joinpath(directory,"mixed-$binary.msh")
            E.write_mixed_msh(path,filled;version=4.1,binary)
            restored=E.read_mixed_msh(path)
            @test validate(restored).ok && restored.physical_names==filled.physical_names
            @test actual_records(restored)==actual_records(filled)
        end
    end
end

@testset "orphan and malformed mixed graphs reject atomically" begin
    source=read(joinpath(FIXTURES,"brep_relationship_mixed_roots.step"),String)
    bad_map=replace(source,"#3050 = MAPPED_ITEM('instance',#3055,#3060);"=>
        "#3050 = MAPPED_ITEM('instance',#3006,#3060);")
    nonmember=replace(source,
        "SHAPE_REPRESENTATION('MappedPart',(#3006,#3050,#3060),#3045)"=>
        "SHAPE_REPRESENTATION('MappedPart',(#3006,#3050),#3045)")
    empty_source=replace(source,
        "SHAPE_REPRESENTATION('source',(#341),#3059)"=>
        "SHAPE_REPRESENTATION('source',(#3056),#3059)")
    orphan_item=append_entities(source,"#3500 = MAPPED_ITEM('x',#3055,#3060);\n")
    orphan_map=append_entities(source,"#3501 = REPRESENTATION_MAP(#3056,#3057);\n")
    orphan_rep=append_entities(source,"#3502 = SHAPE_REPRESENTATION('x',(#341),#3059);\n")
    unplaced=append_entities(source,"#3503 = FACETED_BREP('extra',#36);\n")
    leafless=replace(source,"SHAPE_REPRESENTATION('Part',(#1006,#41,#141),#345)"=>
        "SHAPE_REPRESENTATION('Part',(#1006),#345)")
    overlap=replace(replace(source,"20.,"=>"10.,"),"21.,"=>"11.,")
    cycle=replace(append_entities(source,
        "#3500 = MAPPED_ITEM('self',#3501,#3056);\n"*
        "#3501 = REPRESENTATION_MAP(#3056,#3057);\n"),
        "SHAPE_REPRESENTATION('source',(#341),#3059)"=>
            "SHAPE_REPRESENTATION('source',(#341,#3500,#3056),#3059)")
    bad_item=replace(append_entities(source,"#3502 = DIRECTION('',(1.,1.,0.));\n"),
        "SHAPE_REPRESENTATION('source',(#341),#3059)"=>
            "SHAPE_REPRESENTATION('source',(#341,#3502),#3059)")
    invalid=[bad_map,nonmember,empty_source,orphan_item,orphan_map,orphan_rep,
        unplaced,leafless,overlap,cycle,bad_item]
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
