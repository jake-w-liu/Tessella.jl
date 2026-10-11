module BRepPlane190TransformFunctionalChecks
using Test,Tessella,LinearAlgebra
const M=Tessella.MeshTypes
const FIXTURES=joinpath(@__DIR__,"..","fixtures")
const R=[.8 -.36 .48; .6 .48 -.64; 0. .8 .6]
const HOLE_R=[0. -1. 0.; 1. 0. 0.; 0. 0. 1.]
const T=[2.5,-1.5,4.]

surface_volume(mesh)=sum(dot(mesh.coords[:,cell[1]],
    cross(mesh.coords[:,cell[2]],mesh.coords[:,cell[3]]))/6
    for cell in eachcol(mesh.tris);init=0.)

function material_volume(mesh)
    volumes=[det(hcat(mesh.coords[:,cell[2]]-mesh.coords[:,cell[1]],
                     mesh.coords[:,cell[3]]-mesh.coords[:,cell[1]],
                     mesh.coords[:,cell[4]]-mesh.coords[:,cell[1]]))/6
             for cell in eachcol(mesh.tets)]
    sum(volumes),!isempty(volumes) && all(>(0),volumes)
end

function shell_topology(mesh)
    incidence=Dict{Tuple,Int}()
    for (a,b,c) in eachcol(mesh.tris),edge in ((a,b),(b,c),(c,a))
        # Bind the actual local edge, independently of emission order.
        key=Tuple(sort(collect(edge)))
        incidence[key]=get(incidence,key,0)+1
    end
    all(==(2),values(incidence)),M.nnodes(mesh)-length(incidence)+M.ntris(mesh)
end

function mapped_corners(mesh;hole=false,rotation=R)
    points=[[x,y,z] for z in (0.,1.) for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
    hole && append!(points,[[x,y,z] for z in (0.,1.) for (x,y) in ((.25,.25),(.75,.25),(.75,.75),(.25,.75))])
    length(points)==M.nnodes(mesh) && all(point->count(row->
        isapprox(mesh.coords[:,row],T+rotation*point;atol=1e-12),axes(mesh.coords,2))==1,points)
end

function boundary_certificate(mesh;volume=false,hole=false,rotation=R)
    cells=Tuple[];valid=true
    if volume
        incidence=Dict{Tuple,Int}()
        for cell in eachcol(mesh.tets),face in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
            key=Tuple(sort(cell[collect(face)]))
            incidence[key]=get(incidence,key,0)+1
        end
        valid=all(count->count in (1,2),values(incidence))
        append!(cells,[face for (face,count) in incidence if count==1])
    else
        append!(cells,Tuple.(eachcol(mesh.tris)))
    end
    caps=zeros(3,2);valid &= !isempty(cells)
    for cell in cells
        xyz=mesh.coords[:,collect(cell)];local_xyz=transpose(rotation)*(xyz.-T)
        a,b,c=eachcol(xyz);area=norm(cross(b-a,c-a))/2
        selected=false
        for axis in 1:3,(side,height) in enumerate((0.,1.))
            if all(value->isapprox(value,height;atol=1e-12),local_xyz[axis,:])
                caps[axis,side]+=area;selected=true
                if hole && axis==3
                    center=sum(local_xyz;dims=2)/3
                    valid &= !(0.25<center[1]<0.75 && 0.25<center[2]<0.75)
                end
            end
        end
        if hole && !selected
            selected=any(axis->any(height->all(value->isapprox(value,height;atol=1e-12),local_xyz[axis,:]),(.25,.75)) &&
                all(value->.25-1e-12<=value<=.75+1e-12,local_xyz[3-axis,:]),1:2)
        end
        valid &= selected
    end
    expected=ones(3,2);hole && (expected[3,:].=.75)
    valid && isapprox(caps,expected;atol=1e-12)
end

function contains(mesh,local_point;rotation=R)
    point=T+rotation*local_point
    for cell in eachcol(mesh.tets)
        a,b,c,d=(mesh.coords[:,row] for row in cell)
        weights=hcat(b-a,c-a,d-a)\(point-a)
        all(>=(-1e-12),weights) && sum(weights)<=1+1e-12 && return true
    end
    false
end

function replace_record(source,de,payload)
    lines=String.(split(source,'\n'))
    indices=findall(line->length(line)>=80 && line[73]=='P' && parse(Int,strip(line[65:72]))==de,lines)
    isempty(indices) && error("missing fixture P record $de")
    length(payload)<=64*length(indices) || error("replacement exceeds fixture P allocation")
    for (index,row) in enumerate(indices)
        first=64*(index-1)+1;last=min(64*index,length(payload))
        fragment=first<=last ? payload[first:last] : ""
        lines[row]=rpad(fragment,64)*lines[row][65:end]
    end
    join(lines,'\n')
end

function change_directory(source,de,field,value;second=false)
    lines=String.(split(source,'\n'));sequence=de+(second ? 1 : 0)
    row=only(findall(line->length(line)>=80 && line[73]=='D' && parse(Int,strip(line[74:80]))==sequence,lines))
    first=8*(field-1)+1;last=8*field
    lines[row]=lines[row][1:first-1]*lpad(string(value),8)*lines[row][last+1:end]
    join(lines,'\n')
end

@testset "Plane190 defining matrix maps the actual box" begin
    println("Plane190 DE124 box");flush(stdout)
    path=joinpath(FIXTURES,"brep_plane190_transform_box.igs")
    surface=import_iges(path;fill=false)
    @test validate(surface).ok && M.nnodes(surface)==8 && M.ntris(surface)==12
    @test mapped_corners(surface)
    @test shell_topology(surface)==(true,2)
    @test surface_volume(surface)≈1. atol=1e-12
    @test boundary_certificate(surface)
    filled=import_iges(path);total,positive=material_volume(filled)
    @test validate(filled).ok && positive
    @test total≈1. atol=1e-12
    @test boundary_certificate(filled;volume=true)
end

@testset "Plane190 defining matrix preserves a real hole" begin
    println("Plane190 DE124 through-hole");flush(stdout)
    path=joinpath(FIXTURES,"brep_plane190_transform_exact_hole.igs")
    surface=import_iges(path;fill=false)
    @test validate(surface).ok && M.nnodes(surface)==16 && M.ntris(surface)==32
    @test mapped_corners(surface;hole=true,rotation=HOLE_R)
    @test shell_topology(surface)==(true,0)
    @test surface_volume(surface)≈.75 atol=1e-12
    @test boundary_certificate(surface;hole=true,rotation=HOLE_R)
    filled=import_iges(path);total,positive=material_volume(filled)
    @test validate(filled).ok && positive
    @test total≈.75 atol=1e-12
    @test boundary_certificate(filled;volume=true,hole=true,rotation=HOLE_R)
    @test all(height->!contains(filled,[.5,.5,height];rotation=HOLE_R),(.1,.5,.9)) &&
          contains(filled,[.1,.1,.5];rotation=HOLE_R)
end

@testset "arbitrary defining rotation maps the hole surface" begin
    surface=import_iges(joinpath(FIXTURES,"brep_plane190_transform_hole.igs");fill=false)
    @test validate(surface).ok && mapped_corners(surface;hole=true) && shell_topology(surface)==(true,0)
    @test isapprox(surface_volume(surface),.75;atol=1e-12) && boundary_certificate(surface;hole=true)
end

@testset "transformed plane keeps real model curves and old chart behavior" begin
    path=joinpath(FIXTURES,"brep_plane190_transform_model_box.igs")
    surface=import_iges(path;fill=false);before=mesh_crc(surface)
    @test validate(surface).ok && mapped_corners(surface) && isapprox(surface_volume(surface),1.;atol=1e-12)
    mktempdir() do directory
        control=joinpath(directory,"model.igs")
        write(control,replace_record(read(path,String),9,"110,0,0,1,0,-1,1;"))
        model=import_iges(control;fill=false)
        @test mapped_corners(model) && boundary_certificate(model) && isapprox(surface_volume(model),1.;atol=1e-12)
    end
    @test mesh_crc(surface)==before
    old=import_iges(joinpath(FIXTURES,"brep_plane190_uv_rotated_box.igs");fill=false)
    @test validate(old).ok && mapped_corners(old) && isapprox(surface_volume(old),1.;atol=1e-12)
end

@testset "unsupported defining matrices and dependency transforms reject atomically" begin
    path=joinpath(FIXTURES,"brep_plane190_transform_box.igs")
    source=read(path,String);matrix=133
    surface=import_iges(path;fill=false);before=mesh_crc(surface)
    cases=[
        change_directory(source,7,7,999),
        replace_record(source,matrix,"123,1,0,0,0,0,1,0,0,0,0,1,0;"),
        replace_record(source,matrix,"124,2,0,0,0,0,1,0,0,0,0,1,0;"),
        replace_record(source,matrix,"124,-1,0,0,0,0,1,0,0,0,0,1,0;"),
        change_directory(source,matrix,5,1;second=true),
        change_directory(source,matrix,5,10;second=true),
        change_directory(source,matrix,7,999),
        change_directory(source,matrix,7,matrix),
        replace_record(source,matrix,"124,1,0,0;"),
        change_directory(source,1,7,matrix),
        change_directory(source,9,7,matrix),
        change_directory(source,7,5,0;second=true)
    ]
    mktempdir() do directory
        bad=joinpath(directory,"matrix.igs")
        for text in cases
            write(bad,text)
            @test_throws ArgumentError import_iges(bad;fill=false)
        end
    end
    @test mesh_crc(surface)==before
end
end
