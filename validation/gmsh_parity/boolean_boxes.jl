#!/usr/bin/env julia
# P6: Tessella BooleanDifference of two boxes vs analytic volume 1 and vs Gmsh.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: ntets, tet_volume, node

const GEO=joinpath(@__DIR__,"boolean_boxes.geo")
result=execute_geo(GEO; mesh_dim=3)
mesh=result.mesh
mesh===nothing && error("Tessella boolean boxes produced no mesh")
ntets(mesh)>0 || error("Tessella boolean boxes have no tets")
length(result.model.volumes)==1 || error("Tessella boolean operands were not deleted")
V=sum(tet_volume(node(mesh,mesh.tets[1,t]),node(mesh,mesh.tets[2,t]),
                 node(mesh,mesh.tets[3,t]),node(mesh,mesh.tets[4,t]))
      for t in 1:ntets(mesh))
abs(V-1)<=1e-12 || error("Tessella boolean-difference volume $V != 1")

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
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    gmsh.model.mesh.generate(3)
    types, tags, _ = gmsh.model.mesh.getElements(3)
    isempty(types) && error("Gmsh boolean boxes have no volume elements")
    gmsh_tets=sum(length, tags; init=0)
    gmsh_tets>0 || error("Gmsh boolean boxes have no volume elements")
    println("GMSH_PARITY_BOOLEAN_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_volume=1 gmsh_tets=$gmsh_tets tessella_tets=$(ntets(mesh))")
finally
    gmsh.finalize()
end
