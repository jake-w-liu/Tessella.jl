module StepRelationshipSIChecks
using Test,Tessella,LinearAlgebra
const E=Tessella.Elements
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function volumes(mesh,kind)
    result=Dict{Int32,Float64}()
    for (block,owners) in zip(mesh.blocks,mesh.entity_data.block_entities)
        block.msh==kind || continue
        for (cell,owner) in zip(eachcol(block.nodes),owners)
            a,b,c=(mesh.coords[:,row] for row in cell[1:3])
            value=kind==2 ? dot(a,cross(b,c))/6 :
                det(hcat(b-a,c-a,mesh.coords[:,cell[4]]-a))/6
            result[owner]=get(result,owner,0.)+value
        end
    end
    return result
end

function records(mesh)
    data=mesh.entity_data
    nodes=Dict(uid=>(data.node_entities[row],Tuple(mesh.coords[:,row]))
        for (row,uid) in enumerate(data.external_node_tags))
    cells=Dict(uid=>(Int(block.msh),data.block_entities[index][column],
            Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (index,block) in enumerate(mesh.blocks)
        for (column,uid) in enumerate(data.external_element_tags[index]))
    return nodes,cells
end

function owner_corners(mesh,owner)
    return Set(Tuple(mesh.coords[:,row]) for row in axes(mesh.coords,2)
        if mesh.entity_data.node_entities[row]==(2,Int32(owner)))
end
expected_corners(lo,width)=Set((lo+x,y,z) for x in (0.,width),y in (0.,width),z in (0.,width))

@testset "actual SI declarations scale source rows and the two frame origins independently" begin
    cases=(("brep_relationship_centimetres.step",10.,(30.,50.)),
        ("brep_relationship_metres.step",1000.,(3000.,5000.)),
        ("brep_relationship_mixed_units.step",1.,(30.,59.)),
        ("brep_relationship_shifted_mixed_units.step",1.,(29.,58.)))
    for (file,width,offsets) in cases
        mesh=import_step(joinpath(FIXTURES,file);fill=false)
        @test mesh isa E.MixedMesh && validate(mesh).ok
        @test [mesh.entity_data.entities[(3,owner)].bbox for owner in 1:2]==
            [(offset,0.,0.,offset+width,width,width) for offset in offsets]
        @test size(mesh.coords)==(3,16) && sum(size(block.nodes,2) for block in mesh.blocks)==24
        actual=volumes(mesh,2)
        @test Set(keys(actual))==Set(Int32[1,2]) && all(isapprox(value,width^3;rtol=2e-12) for value in values(actual))
        @test allunique(mesh.entity_data.external_node_tags) && allunique(reduce(vcat,mesh.entity_data.external_element_tags))
        for owner in 1:2
            @test owner_corners(mesh,owner)==expected_corners(offsets[owner],width)
        end
        @test occursin("NAUO #1020 left",mesh.physical_names[(3,1)]) &&
            occursin("NAUO #1030 right",mesh.physical_names[(3,2)])
    end
    # Actual unchanged millimetre control, without repeating the old volume fill.
    millimetres=import_step(joinpath(FIXTURES,"brep_relationship_boxes.step");fill=false)
    @test validate(millimetres).ok && size(millimetres.coords)==(3,16)
    @test owner_corners(millimetres,1)==expected_corners(3.,1.)
    @test owner_corners(millimetres,2)==expected_corners(5.,1.)
end

@testset "one centimetre fill retains positive physical volumes and exact MSH4 records" begin
    filled=import_step(joinpath(FIXTURES,"brep_relationship_centimetres.step"))
    @test validate(filled).ok && Set(keys(volumes(filled,4)))==Set(Int32[1,2])
    @test all(isapprox(value,1000.;rtol=2e-12) for value in values(volumes(filled,4)))
    @test all(det(hcat((filled.coords[:,cell[i]]-filled.coords[:,cell[1]] for i in 2:4)...))>0
        for block in filled.blocks if block.msh==4 for cell in eachcol(block.nodes))
    mktempdir() do directory
        for binary in (false,true)
            path=joinpath(directory,"si-$binary.msh")
            E.write_mixed_msh(path,filled;version=4.1,binary)
            restored=E.read_mixed_msh(path)
            @test validate(restored).ok && restored.physical_names==filled.physical_names
            @test records(restored)==records(filled)
        end
    end
end

@testset "unwitnessed SI declarations and nonfinite physical transport reject" begin
    source=read(joinpath(FIXTURES,"brep_relationship_boxes.step"),String)
    cases=(replace(source,"SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(.QUETTA.,.METRE.)"),
        replace(source,"SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(.MILLI.,.RADIAN.)"),
        replace(source,"PLANE_ANGLE_UNIT()"=>"LENGTH_UNIT()"),
        replace(source,"((#346,#347,#348))"=>"((#9999,#347,#348))"),
        replace(replace(source,"SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(.EXA.,.METRE.)"),
            "CARTESIAN_POINT('',(3.,0.,0.))"=>"CARTESIAN_POINT('',(1.E300,0.,0.))"))
    mktempdir() do directory
        path=joinpath(directory,"invalid.step")
        for text in cases
            write(path,text)
            @test_throws ArgumentError import_step(path;fill=false)
            @test read(path,String)==text
        end
    end
end
end
