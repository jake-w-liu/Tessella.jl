#!/usr/bin/env julia
# P6: Tessella holed unit square vs analytic area 0.75 and vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: ntris, node

const GEO=joinpath(@__DIR__,"t4_hole.geo")
result=execute_geo(GEO; mesh_dim=2)
mesh=result.mesh
mesh===nothing && error("Tessella t4 hole produced no mesh")
ntris(mesh)>0 || error("Tessella t4 hole has no triangles")
area=sum(begin
    a=node(mesh,mesh.tris[1,t]); b=node(mesh,mesh.tris[2,t]); c=node(mesh,mesh.tris[3,t])
    abs((b[1]-a[1])*(c[2]-a[2])-(c[1]-a[1])*(b[2]-a[2]))/2
end for t in 1:ntris(mesh))
abs(area-0.75)<=1e-12 || error("Tessella t4 hole area $area != 0.75")

function find_gmsh_api()
    explicit=get(ENV,"GMSH_JULIA_API","")
    !isempty(explicit) && isfile(explicit) && return explicit
    executable=Sys.which("gmsh")
    executable===nothing && error("gmsh is not on PATH")
    prefix=dirname(dirname(realpath(executable)))
    for path in (joinpath(prefix,"lib","gmsh.jl"),
                 "/opt/homebrew/opt/gmsh/lib/gmsh.jl")
        isfile(path) && return path
    end
    error("could not locate gmsh.jl")
end
include(find_gmsh_api())
gmsh.initialize(["gmsh","-v","0"])
try
    gmsh.open(GEO)
    gmsh.model.mesh.generate(2)
    types, tags, _ = gmsh.model.mesh.getElements(2)
    isempty(types) && error("Gmsh t4 hole has no surface elements")
    gmsh_tris=sum(length, tags; init=0)
    gmsh_tris>0 || error("Gmsh t4 hole has no surface elements")
    println("GMSH_PARITY_T4_HOLE_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_area=0.75 gmsh_tris=$gmsh_tris tessella_tris=$(ntris(mesh))")
finally
    gmsh.finalize()
end
