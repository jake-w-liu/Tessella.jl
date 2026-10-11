module BRepPlane190UVFunctionalChecks
using Test,Tessella,LinearAlgebra
const M=Tessella.MeshTypes
const E=Tessella.Elements
const FIXTURES=joinpath(@__DIR__,"..","fixtures")
const R=[.8 -.36 .48; .6 .48 -.64; 0. .8 .6]
const T=[2.5,-1.5,4.]

function shell_volume(mesh)
    sum(dot(mesh.coords[:,cell[1]],cross(mesh.coords[:,cell[2]],mesh.coords[:,cell[3]]))/6
        for cell in eachcol(mesh.tris);init=0.0)
end

function shell_topology(mesh)
    edges=Dict{Tuple,Int}()
    for (a,b,c) in eachcol(mesh.tris),edge in ((a,b),(b,c),(c,a))
        key=Tuple(sort(collect(edge)))
        edges[key]=get(edges,key,0)+1
    end
    return all(==(2),values(edges)),M.nnodes(mesh)-length(edges)+M.ntris(mesh)
end

function cap_area(mesh,height)
    area=0.0
    for cell in eachcol(mesh.tris)
        all(row->mesh.coords[3,row]==height,cell) || continue
        a,b,c=(mesh.coords[:,row] for row in cell)
        area+=norm(cross(b-a,c-a))/2
    end
    return area
end

function material_volume(mesh)
    values=[det(hcat(mesh.coords[:,cell[2]]-mesh.coords[:,cell[1]],
                     mesh.coords[:,cell[3]]-mesh.coords[:,cell[1]],
                     mesh.coords[:,cell[4]]-mesh.coords[:,cell[1]]))/6 for cell in eachcol(mesh.tets)]
    return sum(values),!isempty(values) && all(>(0),values)
end

function contains(mesh,point)
    for cell in eachcol(mesh.tets)
        a,b,c,d=(mesh.coords[:,row] for row in cell)
        weights=hcat(b-a,c-a,d-a)\(point-a)
        all(>=(-1e-12),weights) && sum(weights)<=1+1e-12 && return true
    end
    return false
end

function boundary_conformity(mesh)
    faces=Dict{Tuple,Int}()
    for cell in eachcol(mesh.tets),local_face in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        key=Tuple(sort(cell[collect(local_face)]))
        faces[key]=get(faces,key,0)+1
    end
    exterior=[face for (face,count) in faces if count==1]
    positive=all(count->count in (1,2),values(faces)) && !isempty(exterior)
    caps=zeros(2)
    for face in exterior
        xyz=mesh.coords[:,collect(face)]
        a,b,c=eachcol(xyz)
        if all(==(0.),xyz[3,:]) || all(==(1.),xyz[3,:])
            center=(a+b+c)/3
            positive &= !(0.25<center[1]<0.75 && 0.25<center[2]<0.75)
            caps[xyz[3,1]==0. ? 1 : 2]+=norm(cross(b-a,c-a))/2
        else
            outer=any(axis->all(==(0.),xyz[axis,:]) || all(==(1.),xyz[axis,:]),1:2)
            inner=any(axis->(all(==(.25),xyz[axis,:]) || all(==(.75),xyz[axis,:])) &&
                              all(value->.25<=value<=.75,xyz[3-axis,:]),1:2)
            positive &= outer || inner
        end
    end
    return positive && all(area->isapprox(area,.75;atol=1e-12),caps)
end


# Change literal P-record bodies while preserving real IGES DE/P identities.
function replace_record(source,de,payload)
    lines=String.(split(source,'\n'))
    indices=findall(line->length(line)>=80 && line[73]=='P' &&
                          parse(Int,strip(line[65:72]))==de,lines)
    isempty(indices) && error("missing fixture P record $de")
    length(payload)<=64*length(indices) || error("replacement exceeds fixture P allocation")
    for (index,row) in enumerate(indices)
        first=64*(index-1)+1;last=min(64*index,length(payload))
        fragment=first<=last ? payload[first:last] : ""
        lines[row]=rpad(fragment,64)*lines[row][65:end]
    end
    return join(lines,'\n')
end

function change_directory(source,de,field,value;second=false)
    lines=String.(split(source,'\n'));sequence=de+(second ? 1 : 0)
    rows=findall(line->length(line)>=80 && line[73]=='D' && parse(Int,strip(line[74:80]))==sequence,lines)
    length(rows)==1 || error("missing fixture directory row $sequence")
    row=only(rows);first=8*(field-1)+1;last=8*field
    lines[row]=lines[row][1:first-1]*lpad(string(value),8)*lines[row][last+1:end]
    join(lines,'\n')
end

function box_boundary_areas(mesh;volume=false)
    cells=Tuple[]
    if volume
        incidence=Dict{Tuple,Int}()
        for cell in eachcol(mesh.tets),face in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
            key=Tuple(sort(cell[collect(face)]));incidence[key]=get(incidence,key,0)+1
        end
        append!(cells,[face for (face,count) in incidence if count==1])
    else
        append!(cells,Tuple.(eachcol(mesh.tris)))
    end
    areas=zeros(3,2);valid=!isempty(cells)
    for cell in cells
        xyz=mesh.coords[:,collect(cell)];local_xyz=transpose(R)*(xyz.-T)
        selected=false
        for axis in 1:3,(side,height) in enumerate((0.,1.))
            if all(value->isapprox(value,height;atol=1e-12),local_xyz[axis,:])
                a,b,c=eachcol(xyz)
                areas[axis,side]+=norm(cross(b-a,c-a))/2;selected=true
            end
        end
        valid &= selected
    end
    return valid,areas
end

@testset "declared190 UV chart maps a rotated translated box" begin
    println("Plane190 UV-only rotated box");flush(stdout)
    path=joinpath(FIXTURES,"brep_plane190_uv_rotated_box.igs")
    surface=import_iges(path;fill=false)
    @test surface isa M.Mesh && validate(surface).ok
    @test M.nnodes(surface)==8 && M.ntris(surface)==12
    @test shell_topology(surface)==(true,2)
    @test shell_volume(surface)≈1. atol=1e-12
    corners=[T+R*[x,y,z] for x in (0.,1.),y in (0.,1.),z in (0.,1.)]
    @test all(corner->count(row->isapprox(surface.coords[:,row],corner;atol=1e-12),axes(surface.coords,2))==1,corners)
    valid,areas=box_boundary_areas(surface)
    @test valid && all(area->isapprox(area,1.;atol=1e-12),areas)
    mesh=import_iges(path)
    total,positive=material_volume(mesh)
    @test validate(mesh).ok && positive
    @test total≈1. atol=1e-12
    valid,areas=box_boundary_areas(mesh;volume=true)
    @test valid && all(area->isapprox(area,1.;atol=1e-12),areas)
    @test contains(mesh,T+R*[.5,.5,.5])
    @test mesh_crc(import_iges(path;fill=false))==mesh_crc(surface)
end

@testset "declared190 UV chart preserves planar holes" begin
    println("Plane190 UV-only through-hole");flush(stdout)
    path=joinpath(FIXTURES,"brep_plane190_uv_through_hole.igs")
    surface=import_iges(path;fill=false)
    @test validate(surface).ok && M.nnodes(surface)==16 && M.ntris(surface)==32
    @test shell_topology(surface)==(true,0)
    @test shell_volume(surface)≈.75 atol=1e-12
    @test cap_area(surface,0.)≈.75 atol=1e-12
    @test cap_area(surface,1.)≈.75 atol=1e-12
    prior=import_iges(joinpath(FIXTURES,"brep_planar_through_hole_trimmed.igs");fill=false)
    @test mesh_crc(surface)==mesh_crc(prior)
    filled=import_iges(path)
    total,positive=material_volume(filled)
    @test validate(filled).ok && positive
    @test total≈.75 atol=1e-12
    @test boundary_conformity(filled)
    @test all(height->!contains(filled,[.5,.5,height]),(.1,.5,.9))
    @test contains(filled,[.1,.1,.5])
end

@testset "chart context crosses the real assembly import route" begin
    surface=import_iges(joinpath(FIXTURES,"brep_plane190_uv_assembly.igs");fill=false)
    @test surface isa E.MixedMesh && validate(surface).ok
    @test count(entity->entity.dim==3,values(surface.entity_data.entities))==2
    owner_volumes=Dict{Int,Float64}()
    for (bi,block) in enumerate(surface.blocks)
        block.msh==2 || continue
        for (column,cell) in enumerate(eachcol(block.nodes))
            a,b,c=(surface.coords[:,row] for row in cell)
            owner=Int(surface.entity_data.block_entities[bi][column])
            owner_volumes[owner]=get(owner_volumes,owner,0.)+dot(a,cross(b,c))/6
        end
    end
    @test length(owner_volumes)==2 && all(volume->isapprox(volume,1.;atol=1e-12),values(owner_volumes))
    old=import_iges(joinpath(FIXTURES,"brep_two_blocks.iges");fill=false)
    @test old isa E.MixedMesh && validate(old).ok
    @test count(entity->entity.dim==3,values(old.entity_data.entities))==2
end

@testset "190 UV chart guards and 3D model priority are atomic" begin
    path=joinpath(FIXTURES,"brep_plane190_uv_rotated_box.igs")
    source=read(path,String);valid=import_iges(path;fill=false);before=mesh_crc(valid)
    cases=[
        change_directory(source,7,5,0;second=true), # unparameterized Form0
        replace_record(source,7,"190,1,3,0;"), # missing reference direction
        replace_record(source,5,"123,-1.44,1.92,-1.8;"), # parallel reference/normal
        replace_record(source,3,"123,0,0,0;"), # zero normal
        replace_record(source,9,"110,0,0,1,0,-1,0;"), # nonzero UV z
        replace_record(source,19,"142,0,29,17,0,1;"), # wrong referenced surface
        change_directory(source,7,7,999), # surface transform unsupported
        change_directory(source,1,7,999), # origin transform unsupported
        change_directory(source,9,7,999), # UV curve transform unsupported
        change_directory(source,9,1,100), # nonlinear entity type
        change_directory(source,9,5,1;second=true), # unbounded line form
        replace_record(source,17,"102,4,17,11,13,15;"), # composite cycle
        replace_record(source,9,"110,0,0,0,0,-.5,0;"), # unclosed literal UV wire
        replace_record(source,5,"123,0,0,0;") # zero reference direction
    ]
    mktempdir() do directory
        bad=joinpath(directory,"chart.igs")
        for text in cases
            write(bad,text)
            @test_throws ArgumentError import_iges(bad;fill=false)
        end
        # A declared190 chart normalizes direction ratios and projects the
        # small normal component out of an approximately perpendicular refdir.
        near=replace_record(source,5,"123,3.999999999999856,3.000000000000192,-0.00000000000018;")
        write(bad,near)
        near_surface=import_iges(bad;fill=false)
        @test shell_volume(near_surface)≈1. atol=1e-12
        dual_source=read(joinpath(FIXTURES,"brep_curve142_dual_space_box.igs"),String)
        dual=import_iges(joinpath(FIXTURES,"brep_curve142_dual_space_box.igs");fill=false)
        # Unsupported UV coordinates do not supersede a valid declared3D curve.
        write(bad,replace_record(dual_source,11,"110,0,0,1,0,-1,1;"))
        @test mesh_crc(import_iges(bad;fill=false))==mesh_crc(dual)
        # The earlier UV-only rejection expands only for a certified190 chart.
        write(bad,replace_record(dual_source,29,"142,0,7,27,0,1;"))
        @test mesh_crc(import_iges(bad;fill=false))==mesh_crc(dual)
        implicit=read(joinpath(FIXTURES,"brep_planar_through_hole_trimmed.igs"),String)
        write(bad,replace_record(implicit,13,"142,0,1,11,0,1;"))
        @test_throws ArgumentError import_iges(bad;fill=false)
    end
    @test mesh_crc(valid)==before
end
end
