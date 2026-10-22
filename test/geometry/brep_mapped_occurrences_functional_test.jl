module StepMappedOccurrenceChecks
using Test, Tessella, LinearAlgebra
const E=Tessella.Elements
const FIXTURE=joinpath(@__DIR__,"..","fixtures","brep_mapped_boxes.step")

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

@testset "actual STEP mapped boxes retain occurrence placement" begin
    mesh=import_step(FIXTURE;fill=false)
    @test mesh isa E.MixedMesh && Tessella.MeshTypes.validate(mesh).ok
    @test [mesh.entity_data.entities[(3,i)].bbox for i in 1:2]==
        [(3.,0.,0.,4.,1.,1.),(5.,0.,0.,6.,1.,1.)]
    @test all(isapprox(material_volumes(mesh,2)[Int32(i)],1.;atol=1e-12) for i in 1:2)
    @test size(mesh.coords)==(3,16) && sum(size(block.nodes,2) for block in mesh.blocks)==24
    @test occursin("MAPPED_ITEM #906",mesh.physical_names[(3,1)]) &&
        occursin("MAPPED_ITEM #907",mesh.physical_names[(3,2)])
    @test allunique(mesh.entity_data.external_node_tags) &&
        allunique(vcat(mesh.entity_data.external_element_tags...))
    @test all(mesh.entity_data.node_entities[row]==(2,Int32(i))
        for i in 1:2 for row in 8i-7:8i)
    filled=import_step(FIXTURE)
    @test Tessella.MeshTypes.validate(filled).ok
    @test all(isapprox(value,1.;atol=1e-12) for value in values(material_volumes(filled,4)))
    @test all(det(hcat((filled.coords[:,cell[i]]-filled.coords[:,cell[1]] for i in 2:4)...))>0
        for block in filled.blocks if block.msh==4 for cell in eachcol(block.nodes))
    mktempdir() do dir
        for binary in (false,true)
            path=joinpath(dir,"occurrences.msh")
            E.write_mixed_msh(path,filled;version=4.1,binary)
            restored=E.read_mixed_msh(path)
            @test restored.coords==filled.coords && restored.physical_names==filled.physical_names
            @test restored.entity_data.node_entities==filled.entity_data.node_entities &&
                restored.entity_data.external_node_tags==filled.entity_data.external_node_tags
        end
    end
end

function variant(source,extra="";replacements=[])
    for (before,after) in replacements
        source=replace(source,before=>after)
    end
    isempty(extra) && return source
    tail=findlast("ENDSEC;",source)
    return source[1:first(tail)-1]*extra*source[first(tail):end]
end

@testset "map origin and nested placements compose actual source rows" begin
    source=read(FIXTURE,String)
    mktempdir() do dir
        path=joinpath(dir,"mapped.step")
        shifted=variant(source,"#910 = CARTESIAN_POINT('',(1.,0.,0.));\n"*
            "#911 = AXIS2_PLACEMENT_3D('',#910,#13,#14);\n";
            replacements=[("REPRESENTATION_MAP(#11,#10)","REPRESENTATION_MAP(#911,#10)")])
        write(path,shifted)
        mesh=import_step(path;fill=false)
        @test [mesh.entity_data.entities[(3,i)].bbox for i in 1:2]==
            [(2.,0.,0.,3.,1.,1.),(4.,0.,0.,5.,1.,1.)]
        @test all(isapprox(value,1.;atol=1e-12) for value in values(material_volumes(mesh,2)))
        nested=variant(source,"#910 = CARTESIAN_POINT('',(10.,2.,0.));\n"*
            "#911 = AXIS2_PLACEMENT_3D('',#910,#13,#14);\n"*
            "#912 = REPRESENTATION_MAP(#11,#908);\n"*
            "#913 = MAPPED_ITEM('parent',#912,#911);\n"*
            "#914 = SHAPE_REPRESENTATION('parent',(#913),#345);\n";
            replacements=[("SHAPE_DEFINITION_REPRESENTATION(#4,#908)","SHAPE_DEFINITION_REPRESENTATION(#4,#914)")])
        write(path,nested);mesh=import_step(path;fill=false)
        @test [mesh.entity_data.entities[(3,i)].bbox for i in 1:2]==
            [(13.,2.,0.,14.,3.,1.),(15.,2.,0.,16.,3.,1.)]
        @test occursin("MAPPED_ITEM #913/MAPPED_ITEM #907",mesh.physical_names[(3,2)])
        single=replace(source,"(#906,#907)"=>"(#906)")
        single=replace(single,"#907 = MAPPED_ITEM('rotated',#900,#905);"=>"")
        write(path,single);mesh=import_step(path;fill=false)
        @test mesh isa Mesh && Tessella.MeshTypes.bounding_box(mesh)==((3.,0.,0.),(4.,1.,1.))
        @test validate(mesh).ok
        default_x=replace(source,"AXIS2_PLACEMENT_3D('',#902,#13,#904)"=>
            "AXIS2_PLACEMENT_3D('',#902,#14,\u0024)")
        write(path,default_x);mesh=import_step(path;fill=false)
        @test mesh.entity_data.entities[(3,2)].bbox==(6.,0.,0.,7.,1.,1.)
    end
end

@testset "unsupported occurrence topology fails explicitly" begin
    source=read(FIXTURE,String)
    cases=[
        replace(source,"REPRESENTATION_MAP(#11,#10)"=>"REPRESENTATION_MAP(#11,#908)"),
        replace(source,"#906 = MAPPED_ITEM('translated',#900,#903)"=>"#906 = MAPPED_ITEM('translated',#9999,#903)"),
        replace(source,"AXIS2_PLACEMENT_3D('',#902,#13,#904)"=>"AXIS2_PLACEMENT_3D('',#902,#13,#13)"),
        # Moving the rotated occurrence's origin to (4,0,0) places its box
        # exactly on the translated occurrence — duplicate solids reject.
        # (Placing it at (3,0,0) is legal: the rotated box lands on x in
        # [2,3] and shares the x=3 face as a conforming welded contact.)
        replace(source,"CARTESIAN_POINT('',(6.,0.,0.))"=>"CARTESIAN_POINT('',(4.,0.,0.))"),
        replace(source,"SHAPE_REPRESENTATION('two occurrences',(#906,#907),#345)"=>"SHAPE_REPRESENTATION('two occurrences',(#906,#907),#9999)"),
        variant(source,"#910 = MAPPED_ITEM('unselected',#900,#903);\n"),
        variant(source,"#910 = SHAPE_DEFINITION_REPRESENTATION(#4,#10);\n"),
        replace(source,"#905 = AXIS2_PLACEMENT_3D('',#902,#13,#904)"=>"#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('',\u0024,\u0024,#902,2.,\u0024)"),
        variant(source,"#910 = ITEM_DEFINED_TRANSFORMATION('','',#11,#903);\n")]
    mktempdir() do dir
        path=joinpath(dir,"bad.step")
        for candidate in cases
            write(path,candidate)
            @test_throws ArgumentError import_step(path;fill=false)
        end
        @test read(FIXTURE,String)==source
    end
end

@testset "touching mapped occurrences weld a shared interface" begin
    source=read(FIXTURE,String)
    # The rotated occurrence's origin at (3,0,0) puts its box on x in [2,3]
    # — a coincident quad-face contact with the translated box on x=3, not
    # a volume overlap. The assembly certificate welds the shared footprint
    # into one Surface bound by both Regions, like any other conforming
    # contact.
    touching=replace(source,"CARTESIAN_POINT('',(6.,0.,0.))"=>
        "CARTESIAN_POINT('',(3.,0.,0.))")
    mktempdir() do dir
        path=joinpath(dir,"touching.step")
        write(path,touching)
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh && Tessella.MeshTypes.validate(mesh).ok
            @test Set(keys(mesh.entity_data.entities))==
                Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
            @test sort(mesh.entity_data.entities[(3,1)].boundaries)==Int32[1,3]
            @test sort(mesh.entity_data.entities[(3,2)].boundaries)==Int32[2,3]
            @test occursin("shared boundary",mesh.physical_names[(2,3)])
            @test allunique(eachcol(mesh.coords))
        end
    end
end
end
