module StepCartesianRigidChecks
using Test,Tessella,LinearAlgebra
const E=Tessella.Elements
const FIXTURE=joinpath(@__DIR__,"..","fixtures","brep_mapped_boxes.step")
const TARGET="#905 = AXIS2_PLACEMENT_3D('',#902,#13,#904);"
function occurrence_variant(source,target;extra="",shift=false)
    text=replace(source,TARGET=>target)
    if shift
        text=replace(text,"REPRESENTATION_MAP(#11,#10)"=>"REPRESENTATION_MAP(#912,#10)")
        extra*="#911 = CARTESIAN_POINT('',(1.,0.,0.));\n#912 = AXIS2_PLACEMENT_3D('',#911,#13,#14);\n"
    end
    position=findlast("ENDSEC;",text)
    return text[1:first(position)-1]*extra*text[first(position):end]
end
function signed_volumes(mesh,kind)
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
corners=Tuple{Float64,Float64,Float64}[(x,y,z) for x in (0.,1.),y in (0.,1.),z in (0.,1.)]
function expected_corners(case)
    first_box=Set((x+(case==:shifted_origin ? 2. : 3.),y,z) for (x,y,z) in corners)
    second_box=case==:default ? Set((6.0 + x,y,z) for (x,y,z) in corners) :
        case==:arbitrary ? Set((6.0 + 0.8*x - 0.6*y,0.6*x + 0.8*y,z) for (x,y,z) in corners) :
        case==:doubled ? Set((6.0 + 2.0*x,2.0*y,2.0*z) for (x,y,z) in corners) :
        case==:mirrored ? Set((6.0 - x,-y,-z) for (x,y,z) in corners) :
        Set((6.0 - y,x,z) for (x,y,z) in corners)
    return first_box,second_box
end
@testset "proper eight-field Cartesian mapped targets preserve actual rigid occurrences" begin
    source=read(FIXTURE,String)
    quarter="#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',#904,#910,#902,1.,#13);"
    axis2="#910 = DIRECTION('',(-1.,0.,0.));\n"
    cases=[(:default,occurrence_variant(source,"#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',\$,\$,\$,#902,\$,\$);")),
        (:quarter,occurrence_variant(source,quarter;extra=axis2)),
        (:inferred,occurrence_variant(source,replace(quarter,",1.,#13)"=>",1.,\$)");extra=axis2)),
        (:shifted_origin,occurrence_variant(source,quarter;extra=axis2,shift=true)),
        (:arbitrary,occurrence_variant(source,"#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',#910,#911,#902,1.,#13);";
            extra="#910 = DIRECTION('',(0.8,0.6,0.));\n#911 = DIRECTION('',(-0.6,0.8,0.));\n")),
        (:doubled,occurrence_variant(source,"#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',\$,\$,\$,#902,2.,\$);")),
        (:mirrored,occurrence_variant(source,"#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',\$,\$,\$,#902,-1.,\$);"))]
    mktempdir() do directory
        path=joinpath(directory,"operators.step")
        for (case,text) in cases
            write(path,text)
            mesh=import_step(path;fill=false)
            @test mesh isa E.MixedMesh && validate(mesh).ok && size(mesh.coords)==(3,16)
            expected=expected_corners(case)
            for owner in 1:2
                actual=[mesh.coords[:,row] for row in axes(mesh.coords,2)
                    if mesh.entity_data.node_entities[row]==(2,Int32(owner))]
                @test length(actual)==8 && all(any(norm(point-collect(witness))<=2e-14 for point in actual)
                    for witness in expected[owner])
            end
            measured=signed_volumes(mesh,2)
            @test isapprox(measured[Int32(1)],1.;atol=2e-13) &&
                isapprox(measured[Int32(2)],case==:doubled ? 8.0 : 1.0;atol=2e-13)
            @test allunique(mesh.entity_data.external_node_tags) &&
                allunique(reduce(vcat,mesh.entity_data.external_element_tags))
            @test occursin("MAPPED_ITEM #907",mesh.physical_names[(3,2)]) &&
                sum(size(block.nodes,2) for block in mesh.blocks)==24
        end
        # One actual fill, using a quarter-turn with analytically unchanged volume.
        write(path,only(text for (case,text) in cases if case==:quarter))
        filled=import_step(path)
        @test validate(filled).ok && Set(keys(signed_volumes(filled,4)))==Set(Int32[1,2])
        @test all(isapprox(volume,1.;atol=2e-12) for volume in values(signed_volumes(filled,4)))
        @test all(det(hcat((filled.coords[:,cell[i]]-filled.coords[:,cell[1]] for i in 2:4)...))>0
            for block in filled.blocks if block.msh==4 for cell in eachcol(block.nodes))
        for binary in (false,true)
            output=joinpath(directory,"operator-$binary.msh")
            E.write_mixed_msh(output,filled;version=4.1,binary)
            restored=E.read_mixed_msh(output)
            @test validate(restored).ok && restored.coords==filled.coords
            println("CARTESIAN_MSH4 binary=",binary," source_families=",
                [block.msh for block in filled.blocks]," restored_families=",
                [block.msh for block in restored.blocks]);flush(stdout)
            @test actual_records(restored)==actual_records(filled)
        end
        invalid=[
            occurrence_variant(source,"#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('', \$,\$,#902,2.,\$);"),
            occurrence_variant(source,replace(quarter,",1.,#13)"=>",0.,#13)");extra=axis2),
            occurrence_variant(source,replace(quarter,"#904,#910"=>"#13,#910");extra=axis2),
            occurrence_variant(source,replace(quarter,"#902,1."=>"#9999,1.");extra=axis2),
            occurrence_variant(source,replace(quarter,"('','',''"=>"(42,'',''");extra=axis2)]
        for text in invalid
            write(path,text)
            @test_throws ArgumentError import_step(path;fill=false)
        end
    end
end
end
