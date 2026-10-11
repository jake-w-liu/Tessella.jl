module PrismaticFillFunctionalChecks
using Test,Tessella,LinearAlgebra
const M=Tessella.MeshTypes
const P=Tessella.Mesh3D
const FIXTURES=joinpath(@__DIR__,"..","fixtures")
const R=[.8 -.36 .48;.6 .48 -.64;0. .8 .6]
const HOLE_R=[0. -1. 0.;1. 0. 0.;0. 0. 1.]
const T=[2.5,-1.5,4.]
function volumes(mesh)
    [det(hcat(mesh.coords[:,cell[2]]-mesh.coords[:,cell[1]],
              mesh.coords[:,cell[3]]-mesh.coords[:,cell[1]],
              mesh.coords[:,cell[4]]-mesh.coords[:,cell[1]]))/6 for cell in eachcol(mesh.tets)]
end
support(cell)=Tuple(sort(collect(cell)))
function boundary_support(mesh)
    incidence=Dict{Tuple,Int}()
    for cell in eachcol(mesh.tets),face in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        key=support(cell[collect(face)]);incidence[key]=get(incidence,key,0)+1
    end
    all(count->count in (1,2),values(incidence)) || return Set{Tuple}()
    Set(key for (key,count) in incidence if count==1)
end
function contains(mesh,point)
    for cell in eachcol(mesh.tets)
        a,b,c,d=(mesh.coords[:,row] for row in cell)
        weights=hcat(b-a,c-a,d-a)\(point-a)
        all(>=(-1e-12),weights) && sum(weights)<=1+1e-12 && return true
    end
    false
end
function triangular_prism(;cyclic=false)
    coordinates=[0. 1. 0. 0. 1. 0.;0. 0. 1. 0. 0. 1.;0. 0. 0. 1. 1. 1.]
    triangles=NTuple{3,Int32}[(1,3,2),(4,5,6)]
    for (a,b) in ((1,2),(2,3),(3,1))
        before,after=cyclic ? (a,b) : minmax(a,b)
        # Only the literal supports select the diagonal; input winding is
        # independent of the global staircase vertex order.
        if before==a
            push!(triangles,(Int32(a),Int32(b),Int32(a+3)))
            push!(triangles,(Int32(b),Int32(b+3),Int32(a+3)))
        else
            push!(triangles,(Int32(a),Int32(b),Int32(b+3)))
            push!(triangles,(Int32(a),Int32(b+3),Int32(a+3)))
        end
    end
    connectivity=Matrix{Int32}(undef,3,length(triangles))
    for (column,cell) in enumerate(triangles);connectivity[:,column]=collect(cell);end
    Mesh(coordinates;tris=connectivity)
end

@testset "retained rotated Plane190 hole has certified fast volume fill" begin
    path=joinpath(FIXTURES,"brep_plane190_transform_hole.igs")
    surface=import_iges(path;fill=false);before=copy(surface.coords);facets=copy(surface.tris)
    @test M.nnodes(surface)==16 && M.ntris(surface)==32 && validate(surface).ok
    direct=P._try_prismatic_surface_fill(surface)
    @test direct!==nothing
    direct===nothing && error("retained rotated prism was not recognized")
    @test M.nnodes(direct)==16 && M.ntets(direct)==24 && validate(direct).ok
    @test isequal(direct.coords,before) && boundary_support(direct)==Set(support(cell) for cell in eachcol(facets))
    signed=volumes(direct)
    @test all(>(0),signed) && isapprox(sum(signed),.75;atol=1e-12)
    @test !contains(direct,T+R*[.5,.5,.5]) && contains(direct,T+R*[.1,.1,.5])
    direct.coords[1,1]=99.
    @test isequal(surface.coords,before) && surface.tris==facets
    # Exercise the unchanged public CAD volume front door, including the
    # failed ordinary/fan producers before this bounded fast path.
    filled=import_iges(path)
    @test M.nnodes(filled)==16 && M.ntets(filled)==24 && validate(filled).ok
    @test isequal(filled.coords,before) && boundary_support(filled)==Set(support(cell) for cell in eachcol(facets))
    @test all(>(0),volumes(filled)) && isapprox(sum(volumes(filled)),.75;atol=1e-12)
    @test !contains(filled,T+R*[.5,.5,.5]) && contains(filled,T+R*[.1,.1,.5])
end

@testset "translated and outer-shape prism controls" begin
    surface=import_iges(joinpath(FIXTURES,"brep_plane190_transform_hole.igs");fill=false)
    for (name,control,expected) in (
        ("translated",Mesh(surface.coords.+[8.,-3.,2.];tris=surface.tris),.75),
        ("outer box",import_iges(joinpath(FIXTURES,"brep_plane190_transform_box.igs");fill=false),1.),
        ("triangle",triangular_prism(),.5))
        @testset "$name" begin
            candidate=P._try_prismatic_surface_fill(control)
            @test candidate!==nothing
            candidate===nothing && error("$name prism was not recognized")
            @test validate(candidate).ok && isequal(candidate.coords,control.coords)
            @test boundary_support(candidate)==Set(support(cell) for cell in eachcol(control.tris)) &&
                  all(>(0),volumes(candidate)) && isapprox(sum(volumes(candidate)),expected;atol=1e-12)
        end
    end
end

@testset "different cap triangulations use the certified public fallback" begin
    path=joinpath(FIXTURES,"brep_plane190_transform_exact_hole.igs")
    surface=import_iges(path;fill=false);before=copy(surface.coords);facets=copy(surface.tris)
    # This retained quarter-turn fixture has different literal triangle
    # supports on its caps. Recognition must preserve the existing fallback.
    @test P._try_prismatic_surface_fill(surface)===nothing
    filled=import_iges(path)
    @test validate(filled).ok && M.ntets(filled)>0
    @test all(>(0),volumes(filled)) && isapprox(sum(volumes(filled)),.75;atol=1e-12)
    @test first(P._certify_surface_fill(surface,filled))
    @test !contains(filled,T+HOLE_R*[.5,.5,.5]) && contains(filled,T+HOLE_R*[.1,.1,.5])
    @test isequal(surface.coords,before) && surface.tris==facets
end

@testset "unsupported and cyclic recognition retain the fallback" begin
    cyclic=triangular_prism(;cyclic=true);before=copy(cyclic.coords);facets=copy(cyclic.tris)
    @test P._try_prismatic_surface_fill(cyclic)===nothing
    @test isequal(cyclic.coords,before) && cyclic.tris==facets
    skew=triangular_prism();skew.coords[1,4]+=.2
    @test P._try_prismatic_surface_fill(skew)===nothing
    source=triangular_prism()
    excessive=Mesh(hcat(source.coords,zeros(3,123));tris=source.tris)
    @test P._try_prismatic_surface_fill(excessive)===nothing
end
end
