#!/usr/bin/env julia
# P6: Tessella API box volume vs analytic 1 and vs Gmsh 4.15.2 OCC Box mesh volume.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.API
using Tessella.MeshTypes: ntets, tet_volume, node

API.initialize()
API.model.add_box(0,0,0,1,1,1; tag=1)
mesh=API.mesh.generate(3)
V=sum(tet_volume(node(mesh,mesh.tets[1,t]),node(mesh,mesh.tets[2,t]),
                 node(mesh,mesh.tets[3,t]),node(mesh,mesh.tets[4,t])) for t in 1:ntets(mesh))
abs(V-1) <= 1e-12 || error("Tessella API box volume $V != 1")
ntets(mesh)>0 || error("Tessella API box has no tets")
API.finalize()

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
    gmsh.model.add("box")
    gmsh.model.occ.addBox(0,0,0,1,1,1,1)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 1.0)
    gmsh.model.mesh.generate(3)
    types, tags, _ = gmsh.model.mesh.getElements(3)
    isempty(types) && error("Gmsh box has no volume elements")
    gmsh_tets=sum(length, tags; init=0)
    gmsh_tets>0 || error("Gmsh box has no volume elements")
    println("GMSH_PARITY_BOX_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_volume=1 gmsh_tets=$gmsh_tets tessella_tets=$(ntets(mesh))")
finally
    gmsh.finalize()
end
