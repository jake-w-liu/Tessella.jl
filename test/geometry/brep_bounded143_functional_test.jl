module BRepBounded143FunctionalChecks
using Test,Tessella,LinearAlgebra
const M=Tessella.MeshTypes
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

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

function replace_inner_wire(source,vertices)
    for (index,de) in enumerate((13,15,17,19))
        a=vertices[index];b=vertices[mod1(index+1,4)]
        payload=join(("110",string.(a)...,string.(b)...),',')*";"
        source=replace_record(source,de,payload)
    end
    return source
end

@testset "IGES143 declared planar boundaries" begin
    println("IGES143 through-hole surface and material");flush(stdout)
    path=joinpath(FIXTURES,"brep_bounded143_through_hole.igs")
    surface=import_iges(path;fill=false)
    @test surface isa M.Mesh && validate(surface).ok
    @test M.nnodes(surface)==16 && M.ntris(surface)==32
    @test shell_topology(surface)==(true,0)
    @test shell_volume(surface)≈.75 atol=1e-12
    @test cap_area(surface,0.)≈.75 atol=1e-12
    @test cap_area(surface,1.)≈.75 atol=1e-12
    prior=import_iges(joinpath(FIXTURES,"brep_planar_through_hole_trimmed.igs");fill=false)
    @test surface.coords==prior.coords
    @test Set(Tuple(sort(cell)) for cell in eachcol(surface.tris))==
          Set(Tuple(sort(cell)) for cell in eachcol(prior.tris))
    mesh=import_iges(path)
    volume,positive=material_volume(mesh)
    @test validate(mesh).ok && positive
    @test volume≈.75 atol=1e-12
    @test boundary_conformity(mesh)
    @test all(height->!contains(mesh,[.5,.5,height]),(.1,.5,.9))
    @test contains(mesh,[.1,.1,.5])
    @test mesh_crc(import_iges(path;fill=false))==mesh_crc(surface)
end

@testset "IGES143 one-boundary Faces retain the closed box" begin
    println("IGES143 one-boundary box");flush(stdout)
    path=joinpath(FIXTURES,"brep_bounded143_single_wire_box.igs")
    surface=import_iges(path;fill=false)
    @test surface isa M.Mesh && validate(surface).ok
    @test M.nnodes(surface)==8 && M.ntris(surface)==12
    @test shell_topology(surface)==(true,2)
    @test shell_volume(surface)≈1. atol=1e-12
    mesh=import_iges(path)
    volume,positive=material_volume(mesh)
    @test validate(mesh).ok && positive
    @test volume≈1. atol=1e-12
    old_surface=import_iges(joinpath(FIXTURES,"brep_occ_box.igs");fill=false)
    @test M.nnodes(old_surface)==8 && M.ntris(old_surface)==12
    @test shell_volume(old_surface)≈1. atol=1e-12
end

@testset "IGES143 malformed boundary admission is atomic" begin
    path=joinpath(FIXTURES,"brep_bounded143_through_hole.igs")
    source=read(path,String)
    valid=import_iges(path;fill=false);before=mesh_crc(valid)
    inside=[(.25,.25,0.),(.75,.25,0.),(.75,.75,0.),(.25,.75,0.)]
    cases=[
        replace_record(source,23,"143,0,1,3,21,11;"), # truncated boundary list
        replace_record(source,23,"143,0,1,2,999,11;"), # absent boundary identity
        replace_record(source,23,"143,0,1,2,21,21;"), # ambiguous repeated wire
        replace_record(source,23,"143,2,1,2,21,11;"), # invalid surface representation
        replace_record(source,23,"514,0,1,2,21,11;"), # DE/P type mismatch
        replace_record(source,21,"141,2,1,1,4,13,1,0,15,1,0,17,1,0,19,1,0;"),
        replace_record(source,21,"141,0,1,25,4,13,1,0,15,1,0,17,1,0,19,1,0;"), # other base surface
        replace_record(source,21,"141,0,1,1,4,999,1,0,15,1,0,17,1,0,19,1,0;"), # missing model curve
        replace_inner_wire(source,[(x+2,y,z) for (x,y,z) in inside]), # no containing outer wire
        replace_inner_wire(source,inside[[1,3,2,4]]), # crossed inner wire
        replace_inner_wire(source,reverse(inside)), # incorrect inner winding
        replace_inner_wire(source,[(x-.25,y,z) for (x,y,z) in inside]), # touching outer wire
        replace_inner_wire(source,[(x,y,.1) for (x,y,z) in inside]) # nonplanar Face domain
    ]
    mktempdir() do directory
        bad=joinpath(directory,"bad.igs")
        for text in cases
            write(bad,text)
            @test_throws ArgumentError import_iges(bad;fill=false)
        end
        # Representation1 can retain an explicit model-space boundary.
        alternate=replace_record(source,23,"143,1,1,2,21,11;")
        write(bad,alternate)
        @test mesh_crc(import_iges(bad;fill=false))==before
    end
    @test mesh_crc(valid)==before
end

@testset "IGES142 model curves remain distinct from UV curves" begin
    dual_path=joinpath(FIXTURES,"brep_curve142_dual_space_box.igs")
    model_path=joinpath(FIXTURES,"brep_curve142_model_space_box.igs")
    dual=import_iges(dual_path;fill=false)
    model=import_iges(model_path;fill=false)
    @test validate(dual).ok && validate(model).ok
    @test mesh_crc(dual)==mesh_crc(model)
    @test M.nnodes(dual)==8 && M.ntris(dual)==12
    @test shell_topology(dual)==(true,2)
    @test shell_volume(dual)≈1. atol=1e-12
    filled=import_iges(dual_path)
    volume,positive=material_volume(filled)
    @test validate(filled).ok && positive
    @test volume≈1. atol=1e-12
    source=read(dual_path,String);before=mesh_crc(dual)
    mktempdir() do directory
        path=joinpath(directory,"boundary.igs")
        write(path,replace_record(source,29,"142,0,7,27,0,1;"))
        @test mesh_crc(import_iges(path;fill=false))==before
        write(path,replace_record(source,29,"142,0,7,27,999,2;"))
        @test_throws ArgumentError import_iges(path;fill=false)
        write(path,replace_record(source,29,"142,0,7,27,25,1;"))
        @test mesh_crc(import_iges(path;fill=false))==before
    end
    @test mesh_crc(dual)==before
end

@testset "historical native torus dual-space boundary compatibility" begin
    surface=import_iges(joinpath(FIXTURES,"brep_occ_torus.igs");fill=false)
    @test surface isa M.Mesh && validate(surface).ok
    @test shell_volume(surface)≈π^2 rtol=.02
end
end
