module CurvedTrimChecks
using Test,Tessella,LinearAlgebra
using Tessella.BRep: parse_step_entities
using Tessella.MeshTypes: node,nnodes,ntris,ntets,tet_volume,validate
using Tessella.Model: GeoModel,add_cylinder!,model_entities,model_closest_point
const E=Tessella.Elements
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function _signed_surface_volume(mesh)
    volume=0.0
    for k in 1:size(mesh.tris,2)
        a,b,c=(mesh.coords[:,mesh.tris[i,k]] for i in 1:3)
        volume+=dot(a,cross(b,c))/6
    end
    return volume
end

_tet_volumes(mesh)=[det(hcat((mesh.coords[:,cell[i]]-mesh.coords[:,cell[1]]
                    for i in 2:4)...))/6 for cell in eachcol(mesh.tets)]

# Every undirected surface edge must be bounded by exactly two triangles,
# once in each direction — the shared-rim welds included.
function _edge_incidence(mesh)
    counts=Dict{NTuple{2,Int32},Int}()
    for t in eachcol(mesh.tris),(a,b) in
            ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
        key=a<b ? (a,b) : (b,a)
        counts[key]=get(counts,key,0)+1
    end
    return counts
end

_polygon_area(r,n)=0.5*n*r*r*sin(2π/n)

function _with_source(body,fixture,replacements...)
    source=read(joinpath(FIXTURES,fixture),String)
    for (old,new) in replacements
        occursin(old,source) || throw(ArgumentError("fixture rewrite anchor missing: $old"))
        source=replace(source,old=>new)
    end
    return mktemp() do path,io
        write(io,source);close(io);body(path)
    end
end

@testset "drilled-plate entity audit" begin
    ents=parse_step_entities(read(joinpath(FIXTURES,"brep_occ_drilled_plate.step"),String))
    kinds=Dict{String,Int}()
    for e in values(ents); kinds[e.kind]=get(kinds,e.kind,0)+1; end
    @test kinds["MANIFOLD_SOLID_BREP"]==1
    @test kinds["CLOSED_SHELL"]==1
    @test kinds["ADVANCED_FACE"]==7
    @test kinds["CYLINDRICAL_SURFACE"]==1
    @test kinds["PLANE"]==6
    @test kinds["CIRCLE"]==2
    @test kinds["LINE"]==13
    @test kinds["EDGE_CURVE"]==15
    @test kinds["FACE_OUTER_BOUND"]==7
    @test kinds["FACE_BOUND"]==2
end

@testset "drilled plate imports as a welded tube + trimmed-cap shell" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_drilled_plate.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tets,2)==0
    @test nnodes(surface)==56 && ntris(surface)==112
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] &&
          vec(maximum(surface.coords,dims=2))≈[4.,4.,1.]
    # The two rim circles contribute 48 nodes; all lie on the wall R=0.5
    # about the (2,2) hole axis and on the cap planes z∈{0,1}.
    rims=[surface.coords[:,i] for i in axes(surface.coords,2) if
          isapprox(hypot(surface.coords[1,i]-2.0,surface.coords[2,i]-2.0),0.5;atol=1e-12)]
    @test length(rims)==48
    @test all(p->p[3]==0.0 || p[3]==1.0,rims)
    @test all(values(_edge_incidence(surface)).==(2))
    # Polygonal hole area is exact for the 24-sided rim grid.
    @test _signed_surface_volume(surface)≈16.0-_polygon_area(0.5,24) atol=1e-12

    filled=import_step(joinpath(FIXTURES,"brep_occ_drilled_plate.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test sum(volumes)≈16.0-_polygon_area(0.5,24) atol=1e-12
end

@testset "drilled-plate wall nodes lie on the analytic cylinder surface" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_drilled_plate.step");fill=false)
    model=GeoModel()
    add_cylinder!(model,2.,2.,0.,0.,0.,1.,0.5)
    surfaces=[tag for (dim,tag) in model_entities(model) if dim==2]
    @test length(surfaces)>=3
    # Identify the lateral Cylinder surface: it is the unique face containing
    # a mid-height point of the hole wall.
    probe=[2.5,2.0,0.5]
    dmin,tag=findmin([begin
        closest,_=model_closest_point(model,2,t,probe)
        norm(closest.-probe)
    end for t in surfaces])
    @test dmin<1e-9
    # Every imported rim-grid node projects onto that analytic surface with
    # zero residual — the curved trim is welded at exactly the wall circle.
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        isapprox(hypot(p[1]-2.0,p[2]-2.0),0.5;atol=1e-12) || continue
        closest,_=model_closest_point(model,2,tag,Vector(p))
        @test norm(closest.-p)<1e-9
    end
end

@testset "hollow cylinder imports with two welded tube walls" begin
    ents=parse_step_entities(read(joinpath(FIXTURES,"brep_occ_hollow_cylinder.step"),String))
    @test count(e->e.kind=="CYLINDRICAL_SURFACE",values(ents))==2
    @test count(e->e.kind=="CIRCLE",values(ents))==4
    @test count(e->e.kind=="ADVANCED_FACE",values(ents))==4

    surface=import_step(joinpath(FIXTURES,"brep_occ_hollow_cylinder.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test nnodes(surface)==96 && ntris(surface)==192
    @test all(values(_edge_incidence(surface)).==(2))
    # Annulus area between the R=1 and R=0.5 polygons times length 2.
    expected=(_polygon_area(1.0,24)-_polygon_area(0.5,24))*2.0
    @test _signed_surface_volume(surface)≈expected atol=1e-12
    inner=[surface.coords[:,i] for i in axes(surface.coords,2) if
           isapprox(hypot(surface.coords[1,i],surface.coords[2,i]),0.5;atol=1e-12)]
    outer=[surface.coords[:,i] for i in axes(surface.coords,2) if
           isapprox(hypot(surface.coords[1,i],surface.coords[2,i]),1.0;atol=1e-12)]
    @test length(inner)==48 && length(outer)==48

    filled=import_step(joinpath(FIXTURES,"brep_occ_hollow_cylinder.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test sum(volumes)≈expected atol=1e-12
end

@testset "curved-trim certificate rejects malformed or unmatched trims" begin
    # A hole circle whose radius disagrees with the wall is malformed.
    _with_source("brep_occ_drilled_plate.step",
                 "CYLINDRICAL_SURFACE('',#185,0.5)"=>
                 "CYLINDRICAL_SURFACE('',#185,0.6)") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("radius does not match",sprint(showerror,err))
    end
    # Dropping the top cap's hole wire leaves the rim unclaimed — the tube
    # wall can no longer close and the shell rejects explicitly through the
    # trimmed-sector certificate's boundary-wire diagnostic.
    _with_source("brep_occ_drilled_plate.step",
                 "#127 = ADVANCED_FACE('',(#118,#121),#126,.T.)"=>
                 "#127 = ADVANCED_FACE('',(#118),#126,.T.)") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("multiple boundary wires",sprint(showerror,err)) ||
              occursin("full rim",sprint(showerror,err))
    end
end

end # module
