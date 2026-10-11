module BRepPlanarFaceHoleFunctionalChecks
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

@testset "planar BRep Face outer and inner wires" begin
    surfaces=M.Mesh[]
    for (importer,file) in ((import_step,"brep_planar_through_hole.step"),
                            (import_step,"brep_occ_planar_through_hole.step"),
                            (import_iges,"brep_planar_through_hole_msbo.igs"),
                            (import_iges,"brep_planar_through_hole_trimmed.igs"))
        println("Holed Face fixture: ",file);flush(stdout)
        surface=importer(joinpath(FIXTURES,file);fill=false)
        push!(surfaces,surface)
        @test surface isa M.Mesh && validate(surface).ok
        @test M.nnodes(surface)==16 && M.ntris(surface)==32
        closed,euler=shell_topology(surface)
        @test closed && euler==0
        @test shell_volume(surface)≈.75 atol=1e-12
        @test cap_area(surface,0.)≈.75 atol=1e-12
        @test cap_area(surface,1.)≈.75 atol=1e-12
        mesh=importer(joinpath(FIXTURES,file))
        volume,positive=material_volume(mesh)
        @test validate(mesh).ok && positive
        @test volume≈.75 atol=1e-12
        @test boundary_conformity(mesh)
        @test all(height->!contains(mesh,[.5,.5,height]),(.1,.5,.9))
        @test contains(mesh,[.1,.1,.5])
    end
    @test all(surface->surface.coords==first(surfaces).coords,surfaces)
    # Native CAD root traversal can select another valid planar diagonal.
    @test all(surface->Set(Tuple(sort(cell)) for cell in eachcol(surface.tris))==
                       Set(Tuple(sort(cell)) for cell in eachcol(first(surfaces).tris)),surfaces[[1,3,4]])
end

@testset "invalid Face wires reject without changing a valid import" begin
    source=read(joinpath(FIXTURES,"brep_planar_through_hole.step"),String)
    surface=import_step(joinpath(FIXTURES,"brep_planar_through_hole.step");fill=false)
    before=mesh_crc(surface)
    mutations=(
        "POLY_LOOP('',(#5,#6,#7,#8))"=>"POLY_LOOP('',(#5,#7,#6,#8))", # crossed inner loop
        "POLY_LOOP('',(#5,#6,#7,#8))"=>"POLY_LOOP('',(#8,#7,#6,#5))", # wrong winding
        "#5 = CARTESIAN_POINT('',(0.25,0.25,0.0))"=>"#5 = CARTESIAN_POINT('',(2.,0.25,0.0))", # outside
        "#5 = CARTESIAN_POINT('',(0.25,0.25,0.0))"=>"#5 = CARTESIAN_POINT('',(0.,0.25,0.0))", # touching
        "#5 = CARTESIAN_POINT('',(0.25,0.25,0.0))"=>"#5 = CARTESIAN_POINT('',(0.25,0.25,0.1))", # off plane
        "#300 = FACE('',(#201,#200))"=>"#300 = FACE('',(#201,#200,#201))", # duplicate wire
        "#201 = FACE_BOUND('',#101,.T.)"=>"#201 = FACE_OUTER_BOUND('',#101,.T.)")
    mktempdir() do directory
        path=joinpath(directory,"bad.step")
        for mutation in mutations
            write(path,replace(source,mutation))
            @test_throws ArgumentError import_step(path;fill=false)
        end
        unmarked=replace(source,"FACE_OUTER_BOUND"=>"FACE_BOUND")
        # Outer-role inference must still reject ambiguous or invalid domains.
        for mutation in (mutations[1],mutations[3],mutations[6])
            write(path,replace(unmarked,mutation))
            @test_throws ArgumentError import_step(path;fill=false)
        end
    end
    @test mesh_crc(surface)==before
end

@testset "unholed torus Face compatibility" begin
    surface=import_iges(joinpath(FIXTURES,"brep_occ_torus.igs");fill=false)
    @test surface isa M.Mesh && validate(surface).ok
    @test shell_volume(surface)≈π^2 rtol=.02
end
end
