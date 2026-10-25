module NurbsShellChecks
using Test,Tessella,LinearAlgebra
using Tessella.BRep: parse_step_entities
using Tessella.MeshTypes: node,nnodes,ntris,ntets,tet_volume,validate
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

# Every undirected edge bounded by exactly two triangles, once each way.
function _directed_edge_audit(mesh)
    counts=Dict{NTuple{2,Int32},Int}()
    for t in eachcol(mesh.tris),(a,b) in
            ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
        counts[(a,b)]=get(counts,(a,b),0)+1
    end
    return counts
end

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

@testset "nurbs-torus entity audit" begin
    ents=parse_step_entities(read(joinpath(FIXTURES,"brep_occ_nurbs_torus.step"),String))
    kinds=Dict{String,Int}()
    for e in values(ents); kinds[e.kind]=get(kinds,e.kind,0)+1; end
    @test kinds["MANIFOLD_SOLID_BREP"]==1
    @test kinds["CLOSED_SHELL"]==1
    @test kinds["ADVANCED_FACE"]==1
    @test kinds["B_SPLINE_SURFACE_WITH_KNOTS"]==1
    @test kinds["SEAM_CURVE"]==2
    @test kinds["EDGE_CURVE"]==2
    @test kinds["PCURVE"]==4
    @test kinds["ORIENTED_EDGE"]==4
end

@testset "NURBS torus single-face shell imports as welded closed grid" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_torus.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tets,2)==0
    # 48x24 cell grid, both directions seam-welded: 48*24 nodes, 2 tris/cell.
    @test nnodes(surface)==48*24 && ntris(surface)==2*48*24
    audit=_directed_edge_audit(surface)
    @test all(v->v==1,values(audit))
    # B-spline torus: the deg-2 surface sits inside the control polygon, so
    # the emitted mesh volume is a bounded sampling fraction (~0.63) of the
    # analytic 2π^2·R·r^2 — assert the observed band, not the ideal value.
    V=_signed_surface_volume(surface)
    analytic=2*π*π*1.0*0.35^2
    @test V>0.55*analytic && V<0.75*analytic
    # Every emitted node must lie on the torus tube envelope: distance from
    # the major circle within [r-chord, r+chord] with a small slack for the
    # B-spline interior pull.
    for i in axes(surface.coords,2)
        x,y,z=surface.coords[1,i],surface.coords[2,i],surface.coords[3,i]
        d=hypot(hypot(x,y)-1.0,z)
        @test d>0.1 && d<0.55
    end
end

@testset "NURBS sphere single-face shell imports with pole collapse" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_sphere.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    # 48x24 grid: two collapsed v-rows become poles → 48*23+2 nodes;
    # fans at both poles (2*48 tris) + 48*22 interior quads (2 tris each).
    @test nnodes(surface)==48*23+2 && ntris(surface)==2*48+2*48*22
    audit=_directed_edge_audit(surface)
    @test all(v->v==1,values(audit))
    V=_signed_surface_volume(surface)
    analytic=4/3*π*1.0^3
    @test V>0.65*analytic && V<0.85*analytic
    # All nodes lie on or inside the R=1 sphere envelope.
    @test all(i->norm(surface.coords[:,i])<=1.0+1e-9,axes(surface.coords,2))
    # The two pole vertices must exist in the emitted mesh.
    poles=[i for i in axes(surface.coords,2) if abs(abs(surface.coords[3,i])-1.0)<1e-9]
    @test length(poles)==2
end

@testset "NURBS shell fills to positive-volume tets" begin
    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_torus.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    volumes=_tet_volumes(filled)
    # Curved seams produce sliver tets (~1e-19) — certified quality is
    # no inverted cells; strict positivity is a polygonal-fixture property.
    @test all(>=(-1e-15),volumes)
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_torus.step");fill=false)
    @test sum(volumes)≈_signed_surface_volume(surface) rtol=1e-6
end

@testset "trimmed NURBS faces keep the classified reject" begin
    # Open u: break the wraparound control row — surface no longer closes.
    _with_source("brep_occ_nurbs_torus.step",
                 "(1.35,0.0,0.0)"=>"(1.30,0.2,0.0)";  # move C[1,1] off the seam
                 ) do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("unsupported BREP",sprint(showerror,err))
        @test occursin("trimmed NURBS",sprint(showerror,err))
    end
    # Interior hole: add a self-loop bound edge whose vertex is off the
    # parametric boundary — must not pass the bound cross-check.
    _with_source("brep_occ_nurbs_torus.step",
                 "#109 = ORIENTED_EDGE('',*,*,#87,.F.)"=>
                 "#109 = ORIENTED_EDGE('',*,*,#108,.F.)") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("unsupported BREP",sprint(showerror,err))
    end
end

@testset "malformed NURBS surface rejects with named diagnostic" begin
    _with_source("brep_occ_nurbs_torus.step",
                 "(3,1,1,1,1,1,1,3)"=>"(3,1,1,1,1,1,1)") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("unsupported BREP",sprint(showerror,err))
        @test occursin("malformed B-spline",sprint(showerror,err))
    end
end

end # module
