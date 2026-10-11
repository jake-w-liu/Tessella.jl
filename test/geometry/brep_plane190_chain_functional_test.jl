module BRepPlane190ChainFunctionalChecks
using Test,Tessella,LinearAlgebra
const M=Tessella.MeshTypes
const FIXTURES=joinpath(@__DIR__,"..","fixtures")
const CHILD_R=[0. -1. 0.;1. 0. 0.;0. 0. 1.]
const CHILD_T=[.5,-.25,1.]
const PARENT_R=[.8 -.36 .48;.6 .48 -.64;0. .8 .6]
const PARENT_T=[2.,-1.25,3.]

surface_volume(mesh)=sum(dot(mesh.coords[:,cell[1]],cross(mesh.coords[:,cell[2]],mesh.coords[:,cell[3]]))/6
    for cell in eachcol(mesh.tris);init=0.)

function mapped_corners(mesh;third=false)
    extra=third ? [.125,-.5,.25] : zeros(3)
    M.nnodes(mesh)==8 && all(point->count(row->isapprox(mesh.coords[:,row],
        PARENT_R*(CHILD_R*point+CHILD_T)+PARENT_T+extra;atol=1e-12),axes(mesh.coords,2))==1,
        ([x,y,z] for x in (0.,1.),y in (0.,1.),z in (0.,1.)))
end

function shell_topology(mesh)
    edges=Dict{Tuple,Int}()
    for (a,b,c) in eachcol(mesh.tris),edge in ((a,b),(b,c),(c,a))
        key=Tuple(sort(collect(edge)));edges[key]=get(edges,key,0)+1
    end
    all(==(2),values(edges)),M.nnodes(mesh)-length(edges)+M.ntris(mesh)
end

function boundary_areas(mesh)
    areas=zeros(3,2);valid=M.ntris(mesh)>0
    for cell in eachcol(mesh.tris)
        xyz=mesh.coords[:,cell]
        local_xyz=transpose(CHILD_R)*(transpose(PARENT_R)*(xyz.-PARENT_T).-CHILD_T)
        a,b,c=eachcol(xyz);area=norm(cross(b-a,c-a))/2;selected=false
        for axis in 1:3,(side,height) in enumerate((0.,1.))
            if all(value->isapprox(value,height;atol=1e-12),local_xyz[axis,:])
                areas[axis,side]+=area;selected=true
            end
        end
        valid &= selected
    end
    valid && all(area->isapprox(area,1.;atol=1e-12),areas)
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

@testset "literal defining matrices compose in declared order" begin
    println("Plane190 ordered two-matrix chart");flush(stdout)
    mesh=import_iges(joinpath(FIXTURES,"brep_plane190_chain_box.igs");fill=false)
    @test validate(mesh).ok && M.nnodes(mesh)==8 && M.ntris(mesh)==12
    @test mapped_corners(mesh)
    @test shell_topology(mesh)==(true,2)
    @test surface_volume(mesh)≈1. atol=1e-12
    @test boundary_areas(mesh)
    collapsed=import_iges(joinpath(FIXTURES,"brep_plane190_chain_collapsed_box.igs");fill=false)
    @test validate(collapsed).ok && mapped_corners(collapsed)
    @test isapprox(surface_volume(collapsed),1.;atol=1e-12) && boundary_areas(collapsed)
    three=import_iges(joinpath(FIXTURES,"brep_plane190_chain_three_box.igs");fill=false)
    @test validate(three).ok && mapped_corners(three;third=true) && isapprox(surface_volume(three),1.;atol=1e-12)
end

@testset "chain transport retains real model curves and single-matrix behavior" begin
    path=joinpath(FIXTURES,"brep_plane190_chain_model_box.igs")
    model=import_iges(path;fill=false)
    @test validate(model).ok && mapped_corners(model) && isapprox(surface_volume(model),1.;atol=1e-12)
    mktempdir() do directory
        control=joinpath(directory,"model.igs")
        write(control,replace_record(read(path,String),9,"110,0,0,1,0,-1,1;"))
        @test mesh_crc(import_iges(control;fill=false))==mesh_crc(model)
    end
    old=import_iges(joinpath(FIXTURES,"brep_plane190_transform_box.igs");fill=false)
    @test validate(old).ok && isapprox(surface_volume(old),1.;atol=1e-12)
end

@testset "invalid matrix chains fail atomically" begin
    path=joinpath(FIXTURES,"brep_plane190_chain_box.igs")
    source=read(path,String);original=import_iges(path;fill=false);before=mesh_crc(original)
    cases=[
        change_directory(source,133,7,999), # missing parent
        change_directory(source,133,7,133), # self-cycle
        change_directory(source,135,7,133), # two-node cycle
        replace_record(source,135,"123,1,0,0,0,0,1,0,0,0,0,1,0;"), # P/DE type mismatch
        change_directory(source,135,5,1;second=true), # unsupported matrix Form1
        replace_record(source,135,"124,-1,0,0,0,0,1,0,0,0,0,1,0;"), # improper reflection
        replace_record(source,135,"124,2,0,0,0,0,1,0,0,0,0,1,0;"), # nonrigid scale
        replace_record(source,135,"124,1,0,0;"), # truncated parent
        replace_record(replace_record(source,133,"124,1,0,0,1e308,0,1,0,0,0,0,1,0;"),
                       135,"124,1,0,0,1e308,0,1,0,0,0,0,1,0;"), # nonfinite composition
        replace_record(replace_record(source,133,"124,1.0000004,0,0,0,0,1,0,0,0,0,1,0;"),
                       135,"124,1.0000004,0,0,0,0,1,0,0,0,0,1,0;"), # accumulated nonrigidity
        change_directory(source,5,7,135), # other frame dependency remains unsupported
        change_directory(source,133,7,-135) # invalid parent identity
    ]
    mktempdir() do directory
        bad=joinpath(directory,"chain.igs")
        for text in cases
            write(bad,text)
            @test_throws ArgumentError import_iges(bad;fill=false)
        end
    end
    @test mesh_crc(original)==before
end
end
