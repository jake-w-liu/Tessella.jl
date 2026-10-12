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

function replace_parameter(source,pointer,payload)
    join((length(line)>=73 && line[73]=='P' && tryparse(Int,strip(line[65:72]))==pointer ?
          rpad(payload,64)*line[65:end] : line for line in split(source,'\n')),'\n')
end

@testset "multiple explicit CAD solid roots" begin
    box_mesh=verify_import(import_step,"brep_two_boxes.step",Dict(1=>1.0,2=>2.0),["#41","#141"])
    block_mesh=verify_import(import_iges,"brep_two_blocks.iges",Dict(1=>1.0,2=>2.0),["#1","#3"])
    verify_import(import_iges,"brep_two_msbo_tets.igs",Dict(1=>1/6,2=>1/3),["#59","#119"])
    @test block_mesh.entity_data.entities[(3,2)].bbox==(3.,0.,0.,5.,1.,1.)
    @test box_mesh.entity_data.entities[(3,2)].bbox==(3.,0.,0.,4.,2.,1.)
    # A solid-root type-124 applies to the triangulated product, so any
    # finite nonsingular affine is admitted: uniform scale tiles volumes by
    # s^3, shear preserves them, a mirror flips winding but keeps outward
    # normals, and a singular matrix rejects atomically.
    blocks=read(joinpath(FIXTURES,"brep_two_blocks.iges"),String)
    mktempdir() do directory
        path=joinpath(directory,"affine.iges")
        write(path,replace_parameter(blocks,5,"124,2,0,0,3,0,2,0,0,0,0,2,0;"))
        mesh=import_iges(path;fill=false)
        @test mesh isa E.MixedMesh && MT.validate(mesh).ok &&
            mesh.entity_data.entities[(3,2)].bbox==(3.,0.,0.,7.,2.,2.)
        measured,positive=volumes(import_iges(path))
        @test positive && isapprox(measured[2],16.0;atol=1e-12)
        write(path,replace_parameter(blocks,5,"124,-1,0,0,3,0,1,0,0,0,0,1,0;"))
        mesh=import_iges(path;fill=false)
        @test mesh isa E.MixedMesh && MT.validate(mesh).ok
        write(path,replace_parameter(blocks,5,"124,1,.5,0,3,0,1,0,0,0,0,1,0;"))
        mesh=import_iges(path;fill=false)
        @test mesh isa E.MixedMesh && MT.validate(mesh).ok &&
            mesh.entity_data.entities[(3,2)].bbox==(3.,0.,0.,5.5,1.,1.)
        write(path,replace_parameter(blocks,5,"124,0,0,0,3,0,0,0,0,0,0,0,0;"))
        @test_throws ArgumentError import_iges(path;fill=false)
    end
end

function simple_volume(mesh)
    sum(det(hcat(mesh.coords[:,cell[2]]-mesh.coords[:,cell[1]],
                 mesh.coords[:,cell[3]]-mesh.coords[:,cell[1]],
                 mesh.coords[:,cell[4]]-mesh.coords[:,cell[1]]))/6
        for cell in eachcol(mesh.tets);init=0.0)
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
        # Tet vertex (1,0,0) lands inside the box's x=1 face: an isolated
        # vertex-in-face T-junction. The touched facet retessellates around
        # the vertex so the contact becomes conforming; the solids keep
        # separate Surfaces and Regions over the zero-area touch.
        write(path,shifted_box(source,(1.,-0.4,-0.4),(2.,0.4,0.4)))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            @test Set(keys(m.entity_data.entities))==Set(
                [(2,1),(3,1),(2,2),(3,2)])
            fill || continue
            measured,positive=volumes(m)
            @test positive
            @test isapprox(measured[1],1/6;atol=1e-12) &&
                  isapprox(measured[2],0.64;atol=1e-12)
        end
        # Box corner enters the tet interior: the shells cross transversally
        # and fragment into tet-only, box-only and intersection Regions. The
        # intersection is the 0.4-edge tetrahedron at the inner box corner.
        write(path,shifted_box(source,(0.2,0.2,0.2),(0.8,0.8,0.8)))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            @test Set(keys(m.entity_data.entities))==Set(
                [(2,t) for t in 1:5] ∪ [(3,t) for t in 1:3])
            fill || continue
            measured,positive=volumes(m)
            @test positive
            @test isapprox(measured[1],0.156;atol=1e-12) &&
                  isapprox(measured[2],0.216-0.4^3/6;atol=1e-12) &&
                  isapprox(measured[3],0.4^3/6;atol=1e-12)
        end
        # Box sits entirely inside the tet: enclosure without facet contact
        # is a nested assembly — under fill the enclosed boundary welds into
        # a shared Surface and both Regions tile the volume exactly.
        write(path,shifted_box(source,(0.05,0.05,0.05),(0.15,0.15,0.15)))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            @test Set(keys(m.entity_data.entities))==Set(
                fill ? [(2,1),(3,1),(2,2),(3,2),(2,3)] : [(2,1),(3,1),(2,2),(3,2)])
            @test allunique(eachcol(m.coords))
        end
        # Box slides under the tet so its top face partially overlaps the
        # tet's z=0 face. The faces' normals oppose, so the coplanar imprint
        # welds the exact overlap square [0,0.5]x[0,0.5] into a shared
        # Surface while the rest of the box face stays its own boundary.
        write(path,shifted_box(source,(-0.5,-0.5,-0.5),(0.5,0.5,0.)))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            data=m.entity_data
            @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
            @test data.entities[(2,3)].bbox==(0.0,0.0,0.0,0.5,0.5,0.0)
            @test sort(data.entities[(3,1)].boundaries)==Int32[1,3]
            @test sort(data.entities[(3,2)].boundaries)==Int32[2,3]
            @test allunique(eachcol(m.coords))
            fill || continue
            measured,positive=volumes(m)
            @test positive
            @test isapprox(measured[1],1/6;atol=1e-12) &&
                  isapprox(measured[2],0.5;atol=1e-12)
        end
        # Box edge runs along a segment of the tet edge (0,0,0)-(1,0,0):
        # a one-dimensional T-junction. The overlay splits the tet edge at
        # x=0.3 and x=0.7 so the contact becomes conforming; the zero-area
        # overlap produces no shared Surface and both solids keep their own
        # coincident boundary nodes.
        write(path,shifted_box(source,(0.3,-0.2,-0.2),(0.7,0.,0.)))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            @test Set(keys(m.entity_data.entities))==Set([(2,1),(3,1),(2,2),(3,2)])
            fill || @test size(m.coords,2)==14
            @test !allunique(eachcol(m.coords))
            fill || continue
            measured,positive=volumes(m)
            @test positive
            @test isapprox(measured[1],1/6;atol=1e-12) &&
                  isapprox(measured[2],0.4*0.2*0.2;atol=1e-12)
        end
        # A thin box piercing the slant face without covering a vertex
        # crosses it transversally: the fragments admit as a three-Region
        # product whose pieces tile both source volumes exactly.
        write(path,shifted_box(source,(0.4,-0.1,0.4),(0.6,0.1,0.6)))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            @test Set(keys(m.entity_data.entities))==Set(
                [(2,t) for t in 1:5] ∪ [(3,t) for t in 1:3])
            fill || continue
            measured,positive=volumes(m)
            @test positive
            @test isapprox(measured[1],1/6-0.007/6;atol=1e-12) &&
                  isapprox(measured[2],0.008-0.007/6;atol=1e-12) &&
                  isapprox(measured[3],0.007/6;atol=1e-12)
        end
        # A box corner coinciding with the tet vertex (1,0,0) is conforming
        # vertex-only contact: admitted unwelded — each solid keeps its own
        # coincident boundary node and no shared entity is created.
        write(path,shifted_box(source,(1.,-1.,-1.),(2.,0.,0.)))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            @test Set(keys(m.entity_data.entities))==Set([(2,1),(3,1),(2,2),(3,2)])
            fill || @test size(m.coords,2)==12
            @test !allunique(eachcol(m.coords))
        end
        # A box edge coinciding with the tet edge (0,0,0)-(1,0,0) is
        # conforming edge-only contact: admitted unwelded likewise.
        write(path,shifted_box(source,(0.,-1.,-1.),(1.,0.,0.)))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            @test Set(keys(m.entity_data.entities))==Set([(2,1),(3,1),(2,2),(3,2)])
            fill || @test size(m.coords,2)==12
            @test !allunique(eachcol(m.coords))
        end
    end
    # A second tet sharing the slant face verbatim is a conforming coincident
    # interface: the assembly welds it into one shared Surface entity bound by
    # both Regions, emits the interface cell once, and deduplicates the rim
    # nodes by exact position.
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
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            data=mesh.entity_data
            @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
            @test sort(data.entities[(3,1)].boundaries)==Int32[1,3]
            @test sort(data.entities[(3,2)].boundaries)==Int32[2,3]
            # One shared-interface cell under tag 3; three rim nodes weld.
            surface_cells=Dict{Int32,Vector{NTuple{3,Int}}}()
            for (bi,block) in enumerate(mesh.blocks)
                block.msh==2 || continue
                for (column,cell) in enumerate(eachcol(block.nodes))
                    owner=data.block_entities[bi][column]
                    push!(get!(surface_cells,owner,NTuple{3,Int}[]),
                          Tuple(sort(Int.(cell))))
                end
            end
            @test length(surface_cells[3])==1
            @test length(surface_cells[1])==3 && length(surface_cells[2])==3
            @test length(surface_cells[1])+length(surface_cells[2])+
                  length(surface_cells[3])==7
            @test allunique(eachcol(mesh.coords))
            shared_rows=Set(surface_cells[3][1])
            @test all(data.node_entities[row]==(2,Int32(3)) for row in shared_rows)
            measured,positive=volumes(mesh)
            @test positive
            fill && @test isapprox(measured[1],1/6;atol=1e-12) &&
                          isapprox(measured[2],1/3;atol=1e-12)
        end
    end
    # Two boxes sharing a quad face weld only when both sides triangulate the
    # face along the same diagonal; a mismatched diagonal is non-conforming
    # contact and rejects. The shared face carries two coincident facet
    # pairs inside one plane patch.
    source=read(joinpath(FIXTURES,"brep_two_boxes.step"),String)
    box_section=source[findfirst("#101 = CARTESIAN_POINT",source)[1]:findfirst(
        "#141 = FACETED_BREP('wide',#140);",source)[end]]
    side_box=x1loop->join(("#101 = CARTESIAN_POINT('',(1.,0.,0.));",
          "#102 = CARTESIAN_POINT('',(2.,0.,0.));",
          "#103 = CARTESIAN_POINT('',(2.,1.,0.));",
          "#104 = CARTESIAN_POINT('',(1.,1.,0.));",
          "#105 = CARTESIAN_POINT('',(1.,0.,1.));",
          "#106 = CARTESIAN_POINT('',(2.,0.,1.));",
          "#107 = CARTESIAN_POINT('',(2.,1.,1.));",
          "#108 = CARTESIAN_POINT('',(1.,1.,1.));",
          "#110 = POLY_LOOP('',(#101,#104,#103,#102));",
          "#111 = POLY_LOOP('',(#105,#106,#107,#108));",
          "#112 = POLY_LOOP('',(#101,#102,#106,#105));",
          "#113 = POLY_LOOP('',(#104,#108,#107,#103));",
          "#114 = POLY_LOOP('',($x1loop));",
          "#115 = POLY_LOOP('',(#102,#103,#107,#106));",
          "#120 = FACE_BOUND('',#110,.T.);",
          "#121 = FACE_BOUND('',#111,.T.);",
          "#122 = FACE_BOUND('',#112,.T.);",
          "#123 = FACE_BOUND('',#113,.T.);",
          "#124 = FACE_BOUND('',#114,.T.);",
          "#125 = FACE_BOUND('',#115,.T.);",
          "#130 = FACE('',(#120));",
          "#131 = FACE('',(#121));",
          "#132 = FACE('',(#122));",
          "#133 = FACE('',(#123));",
          "#134 = FACE('',(#124));",
          "#135 = FACE('',(#125));",
          "#140 = CLOSED_SHELL('',(#130,#131,#132,#133,#134,#135));",
          "#141 = FACETED_BREP('side',#140);"),'\n')
    mktempdir() do directory
        path=joinpath(directory,"quad.step")
        # Same-diagonal shared quad: admitted with a two-cell interface.
        write(path,replace(source,box_section=>side_box("#105,#108,#104,#101")))
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            data=mesh.entity_data
            @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
            @test sort(data.entities[(3,1)].boundaries)==Int32[1,3]
            @test sort(data.entities[(3,2)].boundaries)==Int32[2,3]
            counts=Dict{Int32,Int}()
            for (bi,block) in enumerate(mesh.blocks)
                block.msh==2 || continue
                for owner in data.block_entities[bi]
                    counts[owner]=get(counts,owner,0)+1
                end
            end
            @test counts==Dict(Int32(1)=>10,Int32(2)=>10,Int32(3)=>2)
            @test allunique(eachcol(mesh.coords))
        end
        # Mismatched diagonal: the two sides' shared-face triangulations
        # cross. The coplanar imprint retessellates both faces over the
        # common arrangement — both diagonals split the quad into four
        # triangles around its center — and welds them as one shared
        # Surface bounded by both Regions.
        write(path,replace(source,box_section=>side_box("#101,#105,#108,#104")))
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            data=mesh.entity_data
            @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
            @test data.entities[(2,3)].bbox==(1.0,0.0,0.0,1.0,1.0,1.0)
            @test sort(data.entities[(3,1)].boundaries)==Int32[1,3]
            @test sort(data.entities[(3,2)].boundaries)==Int32[2,3]
            counts=Dict{Int32,Int}()
            for (bi,block) in enumerate(mesh.blocks)
                block.msh==2 || continue
                for owner in data.block_entities[bi]
                    counts[owner]=get(counts,owner,0)+1
                end
            end
            @test counts==Dict(Int32(1)=>10,Int32(2)=>10,Int32(3)=>4)
            @test allunique(eachcol(mesh.coords))
        end
    end
end

@testset "nested enclosed solids weld the enclosed boundary" begin
    # Outer box [0,4]^3 encloses inner box [1,3]^3. The inner Region is a
    # separate material region inside the container's volume.
    mesh=import_step(joinpath(FIXTURES,"brep_nested.step");fill=false)
    @test mesh isa E.MixedMesh
    @test MT.validate(mesh).ok
    data=mesh.entity_data
    @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2)])
    @test data.entities[(3,1)].boundaries==Int32[1]
    @test data.entities[(3,2)].boundaries==Int32[2]
    @test size(mesh.coords,2)==16
    mesh=import_step(joinpath(FIXTURES,"brep_nested.step");fill=true)
    @test mesh isa E.MixedMesh
    @test MT.validate(mesh).ok
    data=mesh.entity_data
    @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
    @test sort(data.entities[(3,1)].boundaries)==Int32[1,3]
    @test sort(data.entities[(3,2)].boundaries)==Int32[2,3]
    counts=Dict{Int32,Int}();tets=Dict{Int32,Int}()
    for (bi,block) in enumerate(mesh.blocks)
        target=block.msh==2 ? counts : block.msh==4 ? tets : continue
        for owner in data.block_entities[bi]
            target[owner]=get(target,owner,0)+1
        end
    end
    # The enclosed boundary emits once as shared Surface 3; the inner
    # Region's own Surface 2 carries no cells.
    @test counts==Dict(Int32(1)=>12,Int32(3)=>12)
    @test sort!(collect(keys(tets)))==Int32[1,2]
    @test allunique(eachcol(mesh.coords))
    # Both Regions share the enclosed interface nodes.
    @test all(e->e==(2,3),
        (data.node_entities[row] for row in axes(mesh.coords,2)
         if all(c->1.0<=c<=3.0,mesh.coords[:,row])))
    # A solid merely touching the container's face coplanarly is still
    # enclosed: its x=4 face is an aligned coincident wall of the merged
    # material — the region classifier welds the flush contact the same way,
    # so Regions tile the union exactly ({A}=64-8, {A,S}=8). One poking
    # through is a transversal crossing and admits as a fragment. Shifting
    # the inner solid's x coordinates — largest first so a produced
    # coordinate never collides with a pending replacement — rebuilds each
    # mutant.
    source=read(joinpath(FIXTURES,"brep_nested.step"),String)
    inner=["(1.,1.,1.)","(3.,1.,1.)","(3.,3.,1.)","(1.,3.,1.)",
           "(1.,1.,3.)","(3.,1.,3.)","(3.,3.,3.)","(1.,3.,3.)"]
    coords(p)=parse.(Float64,split(chop(p;head=1,tail=1),','))
    fmt(v)=string(Int(v))*"."
    shift(dx)=begin
        s=source
        for p in sort(inner;by=p->-coords(p)[1])
            v=coords(p);v[1]+=dx
            s=replace(s,p=>"("*join(fmt.(v),",")*")")
        end
        s
    end
    mktempdir() do directory
        path=joinpath(directory,"shift.step")
        write(path,shift(1.0))   # inner face flush with the container's face
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            @test Set(keys(mesh.entity_data.entities))==Set(
                [(2,1),(3,1),(2,2),(3,2),(2,3)])
            fill || continue
            measured,positive=volumes(mesh)
            @test positive
            @test isapprox(measured[1],56.0;atol=1e-12) &&
                  isapprox(measured[2],8.0;atol=1e-12)
        end
        # The inner solid pokes through the container's face: a transversal
        # crossing fragments the pair into outside-only, inside-only and
        # intersection Regions tiling both volumes.
        write(path,shift(2.0))
        for fill in (false,true)
            m=import_step(path;fill)
            @test m isa E.MixedMesh && MT.validate(m).ok
            @test Set(keys(m.entity_data.entities))==Set(
                [(2,t) for t in 1:5] ∪ [(3,t) for t in 1:3])
            fill || continue
            measured,positive=volumes(m)
            @test positive
            @test isapprox(measured[1],60.0;atol=1e-12) &&
                  isapprox(measured[2],4.0;atol=1e-12) &&
                  isapprox(measured[3],4.0;atol=1e-12)
        end
    end
end

@testset "transversal crossings fragment into region products" begin
    # Unit box A and box B shifted by (0.5,0.5,0.5) interpenetrate: the two
    # boundaries cross along a closed polyline and the assembly admits as
    # three Regions — A only, B only, A∩B — welded at shared interface
    # Surfaces. Volumes 7/8 + 7/8 + 1/8 tile the union exactly.
    file="brep_crossing_boxes.step"
    mesh=import_step(joinpath(FIXTURES,file);fill=false)
    @test mesh isa E.MixedMesh
    @test MT.validate(mesh).ok
    data=mesh.entity_data
    @test Set(keys(data.entities))==Set(
        [(2,t) for t in 1:5] ∪ [(3,t) for t in 1:3])
    @test sort(data.entities[(3,1)].boundaries)==Int32[1,4]
    @test sort(data.entities[(3,2)].boundaries)==Int32[2,5]
    @test sort(data.entities[(3,3)].boundaries)==Int32[3,4,5]
    @test size(mesh.coords,2)==22
    @test allunique(eachcol(mesh.coords))
    @test occursin("∩",mesh.physical_names[(3,3)])
    counts=Dict{Int32,Int}()
    for (bi,block) in enumerate(mesh.blocks)
        block.msh==2 || continue
        for owner in data.block_entities[bi]
            counts[owner]=get(counts,owner,0)+1
        end
    end
    # Own-region surfaces carry the outside pieces; the two shared
    # interface Surfaces carry the coincident cuts once each.
    @test counts==Dict(Int32(1)=>18,Int32(2)=>18,Int32(4)=>6,Int32(5)=>6)
    mesh=import_step(joinpath(FIXTURES,file);fill=true)
    @test mesh isa E.MixedMesh
    @test MT.validate(mesh).ok
    data=mesh.entity_data
    @test Set(keys(data.entities))==Set(
        [(2,t) for t in 1:5] ∪ [(3,t) for t in 1:3])
    measured,positive=volumes(mesh)
    @test positive
    @test isapprox(measured[1],0.875;atol=1e-12) &&
          isapprox(measured[2],0.875;atol=1e-12) &&
          isapprox(measured[3],0.125;atol=1e-12)
    @test allunique(eachcol(mesh.coords))
    # Same-pair mixed contact: the replacement folds B's y faces onto A's
    # (aligned coplanar overlap — an interior wall of the merged material
    # the region classifier dedupes) while its z faces keep crossing A's
    # transversally. The fragment resolves both into three Regions:
    # {A},{B} and {A∩B}=[0.5,1]x[0,1]x[0.5,1]=0.25.
    source=read(joinpath(FIXTURES,file),String)
    mktempdir() do directory
        path=joinpath(directory,"mutant.step")
        write(path,replace(replace(source,",0.5,"=>",0.,"),",1.5,"=>",1.,"))
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            @test sort(Int.(t for (d,t) in keys(mesh.entity_data.entities)
                              if d==3))==Int32[1,2,3]
            fill || continue
            measured,positive=volumes(mesh)
            @test positive
            @test isapprox(measured[1],0.75;atol=1e-12) &&
                  isapprox(measured[2],0.75;atol=1e-12) &&
                  isapprox(measured[3],0.25;atol=1e-12)
        end
        # B touching A's x=1 face over the y,z in [0.5,1] patch is a
        # partial coplanar face overlap: the imprint welds exactly the
        # overlap patch as a shared Surface bounded by both Regions.
        write(path,replace(replace(source,"(0.5,"=>"(1.,"),"(1.5,"=>"(2.,"))
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            data=mesh.entity_data
            @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
            @test data.entities[(2,3)].bbox==(1.0,0.5,0.5,1.0,1.0,1.0)
            @test sort(data.entities[(3,1)].boundaries)==Int32[1,3]
            @test sort(data.entities[(3,2)].boundaries)==Int32[2,3]
            @test allunique(eachcol(mesh.coords))
        end
    end
end

@testset "coplanar face overlaps imprint a shared surface" begin
    # Box A is [0,1]^3; box B starts at [3,4]x[0,2]x[0,1] disjoint. Moving B
    # to [1,2]x[0.5,1.5]x[0,1] puts its x=1 face coplanar with A's x=1 face
    # over the y in [0.5,1] patch with opposed outward normals. The coplanar
    # imprint retessellates both faces over one exact planar arrangement and
    # welds the overlap patch into a shared Surface bounded by both Regions;
    # the off-plane pieces stay each solid's own boundary.
    source=read(joinpath(FIXTURES,"brep_two_boxes.step"),String)
    move_box(reps)=foldl(replace,reps;init=source)
    partial=move_box(["(3.,0.,0.)"=>"(1.,0.5,0.)","(4.,0.,0.)"=>"(2.,0.5,0.)",
                      "(4.,2.,0.)"=>"(2.,1.5,0.)","(3.,2.,0.)"=>"(1.,1.5,0.)",
                      "(3.,0.,1.)"=>"(1.,0.5,1.)","(4.,0.,1.)"=>"(2.,0.5,1.)",
                      "(4.,2.,1.)"=>"(2.,1.5,1.)","(3.,2.,1.)"=>"(1.,1.5,1.)"])
    mktempdir() do directory
        path=joinpath(directory,"imprint.step")
        write(path,partial)
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            data=mesh.entity_data
            @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
            @test data.entities[(2,3)].bbox==(1.0,0.5,0.0,1.0,1.0,1.0)
            @test sort(data.entities[(3,1)].boundaries)==Int32[1,3]
            @test sort(data.entities[(3,2)].boundaries)==Int32[2,3]
            @test allunique(eachcol(mesh.coords))
            fill || continue
            measured,positive=volumes(mesh)
            @test positive
            @test isapprox(measured[1],1.0;atol=1e-12) &&
                  isapprox(measured[2],1.0;atol=1e-12)
        end
        # The same overlap with same-side interiors — B straddles A's volume
        # so the coincident y,z faces share aligned outward normals — is
        # positive-volume interpenetration. The region classifier resolves
        # it into {A}=[0,0.5]x[0,1]^2, {B}=[1,1.5]x[0,1]^2 and
        # {A∩B}=[0.5,1]x[0,1]^2 — three halves tiling the union exactly.
        write(path,move_box(["(3.,0.,0.)"=>"(0.5,0.,0.)","(4.,0.,0.)"=>"(1.5,0.,0.)",
                            "(4.,2.,0.)"=>"(1.5,1.,0.)","(3.,2.,0.)"=>"(0.5,1.,0.)",
                            "(3.,0.,1.)"=>"(0.5,0.,1.)","(4.,0.,1.)"=>"(1.5,0.,1.)",
                            "(4.,2.,1.)"=>"(1.5,1.,1.)","(3.,2.,1.)"=>"(0.5,1.,1.)"]))
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            fill || continue
            measured,positive=volumes(mesh)
            @test positive
            @test isapprox(measured[1],0.5;atol=1e-12) &&
                  isapprox(measured[2],0.5;atol=1e-12) &&
                  isapprox(measured[3],0.5;atol=1e-12)
        end
    end
end

@testset "off-plane T-junctions imprint isolated contact features" begin
    # The tet of brep_disjoint_overlap is defined by four poly-loop vertices
    # (#1-#4), so repositioning them builds arbitrary contact shapes against
    # the fixed box [0.6,1.2]^3. An edge resting on a face interior is a
    # chord contact: the touched face retessellates with the chord as a
    # constrained sub-edge, the contact becomes conforming and stays
    # unwelded (zero area), and both solids keep their own entities.
    fixture=joinpath(FIXTURES,"brep_disjoint_overlap.step")
    source=read(fixture,String)
    function move_tet(p1,p2,p3,p4)
        olds=("#1 = CARTESIAN_POINT('',(0.,0.,0.));",
              "#2 = CARTESIAN_POINT('',(1.,0.,0.));",
              "#3 = CARTESIAN_POINT('',(0.,1.,0.));",
              "#4 = CARTESIAN_POINT('',(0.,0.,1.));")
        out=source
        for (i,(o,p)) in enumerate(zip(olds,(p1,p2,p3,p4)))
            out=replace(out,o=>"#$i = CARTESIAN_POINT('',($(p[1]),$(p[2]),$(p[3])));")
        end
        return out
    end
    mktempdir() do directory
        path=joinpath(directory,"junction.step")
        # Edge (#3,#4) lies inside the box's z=1.2 top face; the blade body
        # sits above the box. Volume det/6 = 0.007875.
        write(path,move_tet((0.7,0.8,1.9),(1.0,0.95,1.5),
            (0.7,0.8,1.2),(1.05,0.75,1.2)))
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            data=mesh.entity_data
            @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2)])
            # Unwelded contact: each solid keeps its own coincident node.
            @test !allunique(eachcol(mesh.coords))
            fill || continue
            measured,positive=volumes(mesh)
            @test positive
            @test isapprox(measured[1],0.007875;atol=1e-12) &&
                  isapprox(measured[2],0.216;atol=1e-12)
        end
        # Vertex #4 lands inside the box's z=1.2 face with the tet above: a
        # vertex-in-face contact imprinted as a Steiner vertex.
        write(path,move_tet((0.9,0.9,1.9),(1.4,0.9,1.6),
            (0.9,1.4,1.5),(0.9,0.9,1.2)))
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            @test Set(keys(mesh.entity_data.entities))==Set(
                [(2,1),(3,1),(2,2),(3,2)])
            @test !allunique(eachcol(mesh.coords))
            fill || continue
            measured,positive=volumes(mesh)
            @test positive
            @test isapprox(measured[1],0.029166666666666653;atol=1e-12) &&
                  isapprox(measured[2],0.216;atol=1e-12)
        end
        # The same blade edge plus a tip (#2) poking through the face into
        # the box mixes a transversal crossing with the chord imprint on the
        # same pair — and contact vertex #4 sits on the box face's diagonal
        # edge in decimal coordinates (1.05+0.75=1.8), a designed incidence
        # Float64 places ~1e-16 off the edge. The classifier resolves both
        # spatially into {tet}, {box} and {tet∩box} Regions tiling the union.
        write(path,move_tet((0.7,0.8,1.9),(0.9,0.9,1.0),
            (0.7,0.8,1.2),(1.05,0.75,1.2)))
        for fill in (false,true)
            mesh=import_step(path;fill)
            @test mesh isa E.MixedMesh
            @test MT.validate(mesh).ok
            @test Set(keys(mesh.entity_data.entities))==Set(
                [(2,1),(2,2),(2,3),(2,4),(2,5),(3,1),(3,2),(3,3)])
            fill || continue
            measured,positive=volumes(mesh)
            @test positive
            @test isapprox(measured[1],0.004083333333333333;atol=1e-12) &&
                  isapprox(measured[2],0.21483333333333332;atol=1e-12) &&
                  isapprox(measured[3],0.0011666666666666663;atol=1e-12)
        end
    end
end

# Emit a minimal FACETED_BREP STEP document for polygon-soup solids: each
# face is a closed polygon given by vertex coordinates with an outward
# winding (matching the hand-written brep_*.step fixtures).
function brep_step_document(solids)
    lines=["ISO-10303-21;","HEADER;",
        "FILE_NAME('generated','2026-10-27T00:00:00',('Tessella'),('Tessella')," *
        "'Tessella','Tessella','');",
        "FILE_DESCRIPTION(('Generated FACETED_BREPs'),'2;1');",
        "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));",
        "ENDSEC;","DATA;"]
    next=1
    ids(name)=begin id=next;next+=1;"#$id" end
    for (name,faces) in solids
        faceids=String[]
        for face in faces
            pts=[ids("p") for p in face]
            for (pid,p) in zip(pts,face)
                push!(lines,"$pid = CARTESIAN_POINT('',($(p[1]),$(p[2]),$(p[3])));")
            end
            loop=ids("l");push!(lines,"$loop = POLY_LOOP('',($(join(pts,","))));")
            bound=ids("b");push!(lines,"$bound = FACE_BOUND('',$loop,.T.);")
            f=ids("f");push!(lines,"$f = FACE('',($bound));");push!(faceids,f)
        end
        shell=ids("s")
        push!(lines,"$shell = CLOSED_SHELL('',($(join(faceids,","))));")
        brep=ids("r");push!(lines,"$brep = FACETED_BREP('$name',$shell);")
    end
    push!(lines,"ENDSEC;","END-ISO-10303-21;")
    join(lines,'\n')
end

# Outward-wound quad faces of box [lo,hi] (same winding as brep_nested.step).
function brep_step_box(lo,hi)
    x1,y1,z1=lo;x2,y2,z2=hi
    [[(x1,y1,z1),(x1,y2,z1),(x2,y2,z1),(x2,y1,z1)],
     [(x1,y1,z2),(x2,y1,z2),(x2,y2,z2),(x1,y2,z2)],
     [(x1,y1,z1),(x2,y1,z1),(x2,y1,z2),(x1,y1,z2)],
     [(x1,y2,z1),(x1,y2,z2),(x2,y2,z2),(x2,y2,z1)],
     [(x1,y1,z1),(x1,y1,z2),(x1,y2,z2),(x1,y2,z1)],
     [(x2,y1,z1),(x2,y2,z1),(x2,y2,z2),(x2,y1,z2)]]
end

# Outward-wound triangle faces of tet v1..v4 (same pattern as
# brep_disjoint_overlap.step).
brep_step_tet(v1,v2,v3,v4)=
    [[v1,v3,v2],[v1,v2,v4],[v1,v4,v3],[v2,v3,v4]]

@testset "mixed contacts decompose into arbitrary regions" begin
    # A welded pair plus a second transversal crossing on a different pair:
    # box A welded to B shares the x=4 face while C pokes through B's x=6
    # wall. The region classifier resolves welds and crossings spatially —
    # Regions {A},{B},{C},{B∩C} tile the union exactly.
    cases=Tuple{String,Vector{Tuple{String,Vector{Vector{NTuple{3,Float64}}}}},Dict{Int,Float64}}[]
    push!(cases,("weld plus crossing",[
        ("A",brep_step_box((0.,0.,0.),(4.,4.,4.))),
        ("B",brep_step_box((4.,0.,0.),(6.,4.,4.))),
        ("C",brep_step_box((5.,1.,1.),(7.,3.,3.)))],
        Dict(1=>64.0,2=>28.0,3=>4.0,4=>4.0)))
    # An enclosed box inside A while D crosses A's wall: the fragment's
    # enclosed pieces regroup into the nested Region {A,S}, the crossing
    # cuts {A∩D} — welds, crossings and enclosure compose spatially.
    push!(cases,("crossing plus enclosure",[
        ("A",brep_step_box((0.,0.,0.),(4.,4.,4.))),
        ("S",brep_step_box((1.,1.,1.),(2.,2.,2.))),
        ("D",brep_step_box((3.5,0.5,0.5),(6.,3.5,3.5)))],
        Dict(1=>58.5,2=>18.0,3=>1.0,4=>4.5)))
    # A crossing pair plus a spike whose apex rests inside A's top face on
    # another pair: imprint and crossing resolve in the same decomposition.
    # Pyramid S: apex (1,1,4), base z=2.6 triangle with area 0.32 —
    # volume 0.32*1.4/3.
    push!(cases,("crossing plus imprint",[
        ("A",brep_step_box((0.,0.,0.),(4.,4.,4.))),
        ("S",brep_step_tet((1.,1.,4.),(0.6,0.6,2.6),(1.4,0.6,2.6),(1.,1.4,2.6))),
        ("D",brep_step_box((3.5,0.5,0.5),(6.,3.5,3.5)))],
        Dict(1=>59.35066666666666,2=>18.0,3=>0.14933333333333335,4=>4.5)))
    # An enclosed pyramid touching A's top face at a single apex point: the
    # {A} region boundary self-touches there, so the region surface carries
    # two coincident apex rows and the fill certifies the pinched link.
    push!(cases,("one-point pinch",[
        ("A",brep_step_box((0.,0.,0.),(4.,4.,4.))),
        ("S",brep_step_tet((1.,1.,4.),(0.6,0.6,2.6),(1.4,0.6,2.6),(1.,1.4,2.6)))],
        Dict(1=>63.85066666666666,2=>0.14933333333333335)))
    # A bipyramid enclosed in A touching the top and bottom faces at one
    # apex each — two independent pinch vertices on the same region.
    bipy_faces=let a=(1.,1.,4.),r=((0.5,0.8,2.),(1.5,0.8,2.),(1.,1.6,2.)),b=(1.,1.,0.)
        [[a,r[1],r[2]],[a,r[2],r[3]],[a,r[3],r[1]],
         [b,r[2],r[1]],[b,r[3],r[2]],[b,r[1],r[3]]]
    end
    push!(cases,("two-point pinch",[
        ("A",brep_step_box((0.,0.,0.),(4.,4.,4.))),
        ("B",bipy_faces)],
        Dict(1=>63.46666666666666,2=>0.5333333333333333)))
    mktempdir() do directory
        path=joinpath(directory,"mixed.step")
        for (name,solids,expected) in cases
            println("mixed assembly: ",name);flush(stdout)
            write(path,brep_step_document(solids))
            for fill in (false,true)
                mesh=import_step(path;fill)
                @test mesh isa E.MixedMesh
                @test MT.validate(mesh).ok
                fill || continue
                measured,positive=volumes(mesh)
                @test positive
                @test all(isapprox(measured[t],v;atol=1e-9)
                          for (t,v) in expected) &&
                      length(measured)==length(expected)
            end
        end
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
        # Sliding block 2 to x in [1,3] puts its x=1 face coincident with
        # block 1's x=1 face over the full quad; the differently triangulated
        # sides are a coplanar overlap, imprinted and welded into shared
        # Surface 3 bounded by both Regions.
        write(path,replace_parameter(iges,5,"124,1,0,0,1,0,1,0,0,0,0,1,0;"))
        for fill in (false,true)
            mesh=import_iges(path;fill)
            @test mesh isa E.MixedMesh && MT.validate(mesh).ok
            data=mesh.entity_data
            @test Set(keys(data.entities))==Set([(2,1),(3,1),(2,2),(3,2),(2,3)])
            @test data.entities[(2,3)].bbox==(1.0,0.0,0.0,1.0,1.0,1.0)
            @test sort(data.entities[(3,1)].boundaries)==Int32[1,3]
            @test sort(data.entities[(3,2)].boundaries)==Int32[2,3]
        end
        # Sliding block 2 to x in [0.5,2.5] interpenetrates block 1 over
        # [0.5,1]: a genuine volumetric Boolean the fragment decomposes into
        # three Regions {1},{2},{1,2} with volumes 0.5,1.5,0.5.
        write(path,replace_parameter(iges,5,"124,1,0,0,0.5,0,1,0,0,0,0,1,0;"))
        for fill in (false,true)
            mesh=import_iges(path;fill)
            @test mesh isa E.MixedMesh && MT.validate(mesh).ok
            fill || continue
            measured,positive=volumes(mesh)
            @test positive
            @test isapprox(measured[1],0.5;atol=1e-12) &&
                  isapprox(measured[2],1.5;atol=1e-12) &&
                  isapprox(measured[3],0.5;atol=1e-12)
        end
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
