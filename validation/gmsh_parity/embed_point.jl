#!/usr/bin/env julia
# P6: Tessella Point-In-Surface embed vs analytic area 1 and vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: ntris, nnodes, node

const GEO=joinpath(@__DIR__,"embed_point.geo")
result=execute_geo(GEO; mesh_dim=2)
mesh=result.mesh
mesh===nothing && error("Tessella embed produced no mesh")
any(i->hypot(mesh.coords[1,i]-0.5,mesh.coords[2,i]-0.5,mesh.coords[3,i])<=1e-12,
    1:nnodes(mesh)) || error("Tessella embed missing interior node (0.5,0.5)")
area=sum(begin
    a=node(mesh,mesh.tris[1,t]); b=node(mesh,mesh.tris[2,t]); c=node(mesh,mesh.tris[3,t])
    abs((b[1]-a[1])*(c[2]-a[2])-(c[1]-a[1])*(b[2]-a[2]))/2
end for t in 1:ntris(mesh))
abs(area-1)<=1e-12 || error("Tessella embed area $area != 1")

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
    _, coord, _ = gmsh.model.mesh.getNodes()
    found=false
    for i in 1:3:length(coord)
        hypot(coord[i]-0.5, coord[i+1]-0.5, coord[i+2])<=1e-9 && (found=true; break)
    end
    found || error("Gmsh embed missing interior node (0.5,0.5)")
    types, tags, _ = gmsh.model.mesh.getElements(2)
    gmsh_tris=sum(length, tags; init=0)
    gmsh_tris>0 || error("Gmsh embed has no surface elements")
    println("GMSH_PARITY_EMBED_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_area=1 gmsh_tris=$gmsh_tris tessella_tris=$(ntris(mesh)) tessella_nodes=$(nnodes(mesh))")
finally
    gmsh.finalize()
end
