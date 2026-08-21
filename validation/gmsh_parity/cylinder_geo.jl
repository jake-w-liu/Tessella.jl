#!/usr/bin/env julia
# P6: Tessella OCC Cylinder volume vs the inscribed 24-gon prism, vs Gmsh tet count.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: ntets, tet_volume, node

const GEO=joinpath(@__DIR__,"cylinder.geo")
result=execute_geo(GEO; mesh_dim=3)
mesh=result.mesh
mesh===nothing && error("Tessella cylinder produced no mesh")
ntets(mesh)>0 || error("Tessella cylinder has no tets")
V=sum(tet_volume(node(mesh,mesh.tets[1,t]),node(mesh,mesh.tets[2,t]),
                 node(mesh,mesh.tets[3,t]),node(mesh,mesh.tets[4,t]))
      for t in 1:ntets(mesh))
prism=0.5*24*sin(2π/24)*2
abs(V-prism)<=1e-12 || error("Tessella cylinder volume $V != prism $prism")

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
    gmsh.option.setNumber("Mesh.MeshSizeMax", 1.0)
    gmsh.model.mesh.generate(3)
    types, tags, _ = gmsh.model.mesh.getElements(3)
    isempty(types) && error("Gmsh cylinder has no volume elements")
    gmsh_tets=sum(length, tags; init=0)
    gmsh_tets>0 || error("Gmsh cylinder has no volume elements")
    println("GMSH_PARITY_CYLINDER_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_prism=$prism gmsh_tets=$gmsh_tets tessella_tets=$(ntets(mesh))")
finally
    gmsh.finalize()
end
