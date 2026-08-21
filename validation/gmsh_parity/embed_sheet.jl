#!/usr/bin/env julia
# P6: Tessella Surface-In-Volume open sheet vs analytic volume 1 and vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: ntets, nnodes, tet_volume, node
using Tessella.Mesh3D: mesh_covers_triangle3

const GEO=joinpath(@__DIR__,"embed_sheet.geo")
result=execute_geo(GEO; mesh_dim=3)
mesh=result.mesh
mesh===nothing && error("Tessella embed-sheet produced no mesh")
ntets(mesh)>0 || error("Tessella embed-sheet has no tets")
mesh_covers_triangle3(mesh,(0.2,0.2,0.5),(0.8,0.2,0.5),(0.5,0.8,0.5)) ||
    error("Tessella embed-sheet is not a union of tet faces")
V=sum(tet_volume(node(mesh,mesh.tets[1,t]),node(mesh,mesh.tets[2,t]),
                 node(mesh,mesh.tets[3,t]),node(mesh,mesh.tets[4,t]))
      for t in 1:ntets(mesh))
abs(V-1)<=1e-12 || error("Tessella embed-sheet volume $V != 1")

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
    isempty(types) && error("Gmsh embed-sheet has no volume elements")
    gmsh_tets=sum(length, tags; init=0)
    gmsh_tets>0 || error("Gmsh embed-sheet has no volume elements")
    _, coord, _ = gmsh.model.mesh.getNodes()
    function has_node(x,y,z)
        for i in 1:3:length(coord)
            hypot(coord[i]-x, coord[i+1]-y, coord[i+2]-z)<=1e-8 && return true
        end
        return false
    end
    has_node(0.2,0.2,0.5) || error("Gmsh embed-sheet missing (0.2,0.2,0.5)")
    has_node(0.8,0.2,0.5) || error("Gmsh embed-sheet missing (0.8,0.2,0.5)")
    has_node(0.5,0.8,0.5) || error("Gmsh embed-sheet missing (0.5,0.8,0.5)")
    println("GMSH_PARITY_EMBED_SHEET_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_volume=1 gmsh_tets=$gmsh_tets tessella_tets=$(ntets(mesh)) tessella_nodes=$(nnodes(mesh))")
finally
    gmsh.finalize()
end
