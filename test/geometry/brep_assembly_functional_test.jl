module BRepAssemblyFunctionalChecks
using Test, Tessella, LinearAlgebra
const E=Tessella.Elements
const MT=Tessella.MeshTypes
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function volumes(mesh)
    result=Dict{Int,Float64}()
    positive=true
    for (bi,block) in enumerate(mesh.blocks)
        block.msh==4 || continue
        for (column,cell) in enumerate(eachcol(block.nodes))
            a,b,c,d=(mesh.coords[:,row] for row in cell)
            volume=det(hcat(b-a,c-a,d-a))/6
            positive &= volume>0
            owner=Int(mesh.entity_data.block_entities[bi][column])
            result[owner]=get(result,owner,0.0)+volume
        end
    end
    return result,positive
end

function surface_volumes(mesh)
    result=Dict{Int,Float64}()
    for (bi,block) in enumerate(mesh.blocks)
        block.msh==2 || continue
        for (column,cell) in enumerate(eachcol(block.nodes))
            a,b,c=(mesh.coords[:,row] for row in cell)
            owner=Int(mesh.entity_data.block_entities[bi][column])
            result[owner]=get(result,owner,0.0)+dot(a,cross(b,c))/6
        end
    end
    return result
end

function boundary_and_classification(mesh)
    actual=Dict{Int,Set{Tuple}}();incidences=Dict{Int,Dict{Tuple,Int}}()
    node_sets=Dict{Int,Set{Int32}}()
    metadata=mesh.entity_data
    valid=true
    for (bi,block) in enumerate(mesh.blocks)
        for (column,cell) in enumerate(eachcol(block.nodes))
            owner=Int(metadata.block_entities[bi][column])
            union!(get!(node_sets,owner,Set{Int32}()),cell)
            valid &= all(row->metadata.node_entities[row][2]==owner,cell)
            if block.msh==2
                push!(get!(actual,owner,Set{Tuple}()),Tuple(sort(cell)))
            elseif block.msh==4
                counts=get!(incidences,owner,Dict{Tuple,Int}())
                for local_face in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
                    face=Tuple(sort(cell[collect(local_face)]))
                    counts[face]=get(counts,face,0)+1
                end
            end
        end
    end
    for (owner,counts) in incidences
        valid &= get(actual,owner,Set{Tuple}())==Set(face for (face,count) in counts if count==1)
        valid &= all(count->count in (1,2),values(counts))
    end
    for (owner,faces) in actual
        boundary=Set(row for face in faces for row in face)
        valid &= all(row->metadata.node_entities[row]==(row in boundary ? 2 : 3,Int32(owner)),node_sets[owner])
    end
    return valid,isempty(intersect(node_sets[1],node_sets[2]))
end

function signature(mesh)
    data=mesh.entity_data
    return Dict(data.external_element_tags[bi][column]=>
        (block.msh,data.block_entities[bi][column],
         Tuple(data.external_node_tags[row] for row in cell))
        for (bi,block) in enumerate(mesh.blocks)
        for (column,cell) in enumerate(eachcol(block.nodes)))
end

function verify_import(importer,file,expected,identities)
    println("Assembly fixture: ",file);flush(stdout)
    path=joinpath(FIXTURES,file)
    surface=importer(path;fill=false)
    @test surface isa E.MixedMesh
    @test MT.validate(surface).ok
    @test all(isapprox(surface_volumes(surface)[owner],volume;atol=1e-12)
              for (owner,volume) in expected)
    @test all(block.msh==2 for block in surface.blocks)
    mesh=importer(path)
    @test mesh isa E.MixedMesh
    @test MT.validate(mesh).ok
    measured,positive=volumes(mesh)
    @test positive
    @test all(isapprox(measured[owner],volume;atol=1e-12) for (owner,volume) in expected)
    @test Set(keys(mesh.entity_data.entities))==Set([(2,1),(3,1),(2,2),(3,2)])
    @test all(mesh.entity_data.entities[(3,owner)].boundaries==[owner] for owner in 1:2)
    @test all(occursin(identity,mesh.physical_names[(3,owner)]) for (owner,identity) in enumerate(identities))
    @test all(mesh.entity_data.entities[(dim,owner)].physical_tags==[owner] for dim in 2:3,owner in 1:2)
    @test allunique(mesh.entity_data.external_node_tags)
    @test allunique(vcat(mesh.entity_data.external_element_tags...))
    conforms,separate=boundary_and_classification(mesh)
    @test conforms
    @test separate
    @test signature(importer(path))==signature(mesh)
    mktempdir() do directory
        output=joinpath(directory,"assembly.msh")
        E.write_mixed_msh(output,mesh;version=4.1)
        restored=E.read_mixed_msh(output)
        @test MT.validate(restored).ok
        @test signature(restored)==signature(mesh)
        before=Dict(uid=>mesh.entity_data.node_entities[row] for (row,uid) in enumerate(mesh.entity_data.external_node_tags))
        after=Dict(uid=>restored.entity_data.node_entities[row] for (row,uid) in enumerate(restored.entity_data.external_node_tags))
        @test before==after
        @test restored.physical_names==mesh.physical_names
        restored_volumes,_=volumes(restored)
        @test all(isapprox(restored_volumes[owner],volume;atol=1e-12) for (owner,volume) in expected)
    end
    return mesh
end

@testset "multiple explicit CAD solid roots" begin
    box_mesh=verify_import(import_step,"brep_two_boxes.step",Dict(1=>1.0,2=>2.0),["#41","#141"])
    block_mesh=verify_import(import_iges,"brep_two_blocks.iges",Dict(1=>1.0,2=>2.0),["#1","#3"])
    verify_import(import_iges,"brep_two_msbo_tets.igs",Dict(1=>1/6,2=>1/3),["#59","#119"])
    @test block_mesh.entity_data.entities[(3,2)].bbox==(3.,0.,0.,5.,1.,1.)
    @test box_mesh.entity_data.entities[(3,2)].bbox==(3.,0.,0.,4.,2.,1.)
end

function simple_volume(mesh)
    sum(det(hcat(mesh.coords[:,cell[2]]-mesh.coords[:,cell[1]],
                 mesh.coords[:,cell[3]]-mesh.coords[:,cell[1]],
                 mesh.coords[:,cell[4]]-mesh.coords[:,cell[1]]))/6
        for cell in eachcol(mesh.tets);init=0.0)
end

function replace_parameter(source,pointer,payload)
    join((length(line)>=73 && line[73]=='P' && tryparse(Int,strip(line[65:72]))==pointer ?
          rpad(payload,64)*line[65:end] : line for line in split(source,'\n')),'\n')
end

@testset "single-solid compatibility" begin
    for (importer,file,nodes,faces,volume) in (
        (import_step,"brep_poly_loop_box.step",8,12,1.0),
        (import_iges,"unit_box.iges",8,12,2.0),
        (import_iges,"brep_occ_box.igs",8,12,1.0),
        (import_iges,"brep_iges_msbo_tet.igs",4,4,1/6))
        surface=importer(joinpath(FIXTURES,file);fill=false)
        @test surface isa MT.Mesh
        @test MT.nnodes(surface)==nodes && MT.ntris(surface)==faces
        mesh=importer(joinpath(FIXTURES,file))
        @test mesh isa MT.Mesh && MT.validate(mesh).ok
        @test simple_volume(mesh)≈volume atol=1e-12
    end
end

function shifted_box(source,lo,hi)
    lox,loy,loz=lo;hix,hiy,hiz=hi
    corners=((lox,loy,loz),(hix,loy,loz),(hix,hiy,loz),(lox,hiy,loz),
             (lox,loy,hiz),(hix,loy,hiz),(hix,hiy,hiz),(lox,hiy,hiz))
    olds=((0.6,0.6,0.6),(1.2,0.6,0.6),(1.2,1.2,0.6),(0.6,1.2,0.6),
          (0.6,0.6,1.2),(1.2,0.6,1.2),(1.2,1.2,1.2),(0.6,1.2,1.2))
    for (i,(p,o)) in enumerate(zip(corners,olds))
        source=replace(source,"#$(100+i) = CARTESIAN_POINT('',($(o[1]),$(o[2]),$(o[3])));"=>
                               "#$(100+i) = CARTESIAN_POINT('',($(p[1]),$(p[2]),$(p[3])));")
    end
    return source
end

@testset "overlapping boxes admit only genuinely disjoint solids" begin
    fixture=joinpath(FIXTURES,"brep_disjoint_overlap.step")
    expected=Dict(1=>1/6,2=>0.216)
    mesh=verify_import(import_step,"brep_disjoint_overlap.step",expected,["#41","#141"])
    @test mesh.entity_data.entities[(3,2)].bbox==(0.6,0.6,0.6,1.2,1.2,1.2)
    source=read(fixture,String)
    mktempdir() do directory
        path=joinpath(directory,"mutant.step")
        # Tet vertex (1,0,0) lands inside the box's x=1 face: vertex contact.
        write(path,shifted_box(source,(1.,-0.4,-0.4),(2.,0.4,0.4)))
        @test_throws ArgumentError import_step(path;fill=false)
        # Box corner enters the tet interior: shells interpenetrate.
        write(path,shifted_box(source,(0.2,0.2,0.2),(0.8,0.8,0.8)))
        @test_throws ArgumentError import_step(path;fill=false)
        # Box sits entirely inside the tet: enclosure without facet contact.
        write(path,shifted_box(source,(0.05,0.05,0.05),(0.15,0.15,0.15)))
        @test_throws ArgumentError import_step(path;fill=false)
        # Box slides so one face overlaps part of the tet's z=0 face.
        write(path,shifted_box(source,(-0.5,-0.5,-0.5),(0.5,0.5,0.)))
        @test_throws ArgumentError import_step(path;fill=false)
        # Box edge runs along a tet edge segment: edge-only contact.
        write(path,shifted_box(source,(0.3,-0.2,-0.2),(0.7,0.,0.)))
        @test_throws ArgumentError import_step(path;fill=false)
        # A thin box piercing the slant face without covering a vertex:
        # non-coincident penetration remains rejected.
        write(path,shifted_box(source,(0.4,-0.1,0.4),(0.6,0.1,0.6)))
        @test_throws ArgumentError import_step(path;fill=false)
    end
    # A second tet sharing the slant face verbatim is a coincident-facet
    # contact and still rejects atomically.
    box_start=findfirst("#101 = CARTESIAN_POINT",source)[1]
    box_stop=findfirst("#141 = FACETED_BREP('box',#140);",source)[end]
    box_section=source[box_start:box_stop]
    mirror=join(("#201 = CARTESIAN_POINT('',(1.,0.,0.));",
                "#202 = CARTESIAN_POINT('',(0.,1.,0.));",
                "#203 = CARTESIAN_POINT('',(0.,0.,1.));",
                "#204 = CARTESIAN_POINT('',(1.,1.,1.));",
                "#210 = POLY_LOOP('',(#201,#203,#202));",
                "#211 = POLY_LOOP('',(#201,#202,#204));",
                "#212 = POLY_LOOP('',(#202,#203,#204));",
                "#213 = POLY_LOOP('',(#203,#201,#204));",
                "#220 = FACE_BOUND('',#210,.T.);",
                "#221 = FACE_BOUND('',#211,.T.);",
                "#222 = FACE_BOUND('',#212,.T.);",
                "#223 = FACE_BOUND('',#213,.T.);",
                "#230 = FACE('',(#220));",
                "#231 = FACE('',(#221));",
                "#232 = FACE('',(#222));",
                "#233 = FACE('',(#223));",
                "#240 = CLOSED_SHELL('',(#230,#231,#232,#233));",
                "#241 = FACETED_BREP('mirror',#240);"),'\n')
    mktempdir() do directory
        path=joinpath(directory,"twin.step")
        write(path,replace(source,box_section=>mirror))
        @test_throws ArgumentError import_step(path;fill=false)
    end
end

@testset "assembly errors leave existing import products intact" begin
    good=import_step(joinpath(FIXTURES,"brep_two_boxes.step");fill=false)
    before=E.mixed_crc(good)
    source=read(joinpath(FIXTURES,"brep_two_boxes.step"),String)
    iges=read(joinpath(FIXTURES,"brep_two_blocks.iges"),String)
    msbo=read(joinpath(FIXTURES,"brep_two_msbo_tets.igs"),String)
    mktempdir() do directory
        path=joinpath(directory,"bad.step")
        write(path,replace(source,"FACETED_BREP('wide',#140)"=>"FACETED_BREP('wide',#999)"))
        @test_throws ArgumentError import_step(path)
        write(path,replace(source,"#140 = CLOSED_SHELL('',(#130,#131,#132,#133,#134,#135))"=>
                                 "#140 = CLOSED_SHELL('',(#130,#131,#132,#133,#134))"))
        @test_throws ArgumentError import_step(path;fill=false)
        write(path,replace(source,"ENDSEC;\nEND-ISO-10303-21;"=>"#200 = CSG_SOLID('unsupported',#41);\nENDSEC;\nEND-ISO-10303-21;"))
        @test_throws ArgumentError import_step(path;fill=false)
        path=joinpath(directory,"bad.iges")
        write(path,replace_parameter(iges,5,"124,1,0,0,1,0,1,0,0,0,0,1,0;"))
        @test_throws ArgumentError import_iges(path;fill=false) # touching boxes
        write(path,replace_parameter(iges,5,"124,1,0,0,0.5,0,1,0,0,0,0,1,0;"))
        @test_throws ArgumentError import_iges(path;fill=false) # overlap
        write(path,replace_parameter(iges,3,"150,-2.,1.,1.,0.,0.,0.;"))
        @test_throws ArgumentError import_iges(path)
        write(path,replace_parameter(msbo,59,"514,57,1,0;"))
        @test_throws ArgumentError import_iges(path;fill=false)
        # Minimal old flat layout has no source identities for two regions.
        write(path,rpad("150,1,1,1,0,0,0;",72)*"P      1\n"*
                   rpad("150,1,1,1,3,0,0;",72)*"P      2\n")
        @test_throws ArgumentError import_iges(path;fill=false)
        write(path,rpad("150,1,1,1,0,0,0;",64)*lpad(1,8)*"P      1\n"*
                   rpad("150,1,1,1,3,0,0;",64)*lpad(1,8)*"P      2\n")
        @test_throws ArgumentError import_iges(path;fill=false)
    end
    @test E.mixed_crc(good)==before
    @test_throws ArgumentError import_step(joinpath(FIXTURES,"unsupported.step"))
    @test_throws ArgumentError import_iges(joinpath(FIXTURES,"unsupported.iges"))
end
end
