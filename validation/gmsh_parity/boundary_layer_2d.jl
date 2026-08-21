#!/usr/bin/env julia
# P6: Tessella 2-D boundary-layer quads vs analytic strip area and vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: triangle_area

const GEO=joinpath(@__DIR__,"boundary_layer_2d.geo")
const HW=0.04
const RA=1.3
const NL=3
H=HW*(RA^NL-1)/(RA-1)
coords=Float64[0 1; 0 0; 0 0]
segs=Int32[1; 2;;]
edge=Mesh(coords; segs=segs)
mesh=mesh_boundary_layer_2d(edge; hwall=HW, ratio=RA, nlayers=NL)
quads=only(b for b in mesh.blocks if b.msh==3)
tessella_quads=size(quads.nodes,2)
let area=0.0
    for cell in axes(quads.nodes,2)
        p1=ntuple(d->mesh.coords[d,quads.nodes[1,cell]],3)
        p2=ntuple(d->mesh.coords[d,quads.nodes[2,cell]],3)
        p3=ntuple(d->mesh.coords[d,quads.nodes[3,cell]],3)
        p4=ntuple(d->mesh.coords[d,quads.nodes[4,cell]],3)
        area+=triangle_area(p1,p2,p4)+triangle_area(p4,p2,p3)
    end
    abs(area-H)<=1e-12 || error("Tessella 2-D BL area $area != $H")
end
tessella_quads==NL || error("Tessella 2-D BL quad count $tessella_quads != $NL")

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
    isempty(types) && error("Gmsh 2-D BL has no surface elements")
    gmsh_quads=0
    gmsh_tris=0
    for (t,tg) in zip(types, tags)
        t==3 && (gmsh_quads+=length(tg))
        t==2 && (gmsh_tris+=length(tg))
    end
    gmsh_quads>0 || error("Gmsh 2-D BL produced no quadrangles")
    gmsh_tris>0 || error("Gmsh 2-D BL produced no remaining triangles")
    println("GMSH_PARITY_BL2D_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_area=$H gmsh_quads=$gmsh_quads gmsh_tris=$gmsh_tris tessella_quads=$tessella_quads")
finally
    gmsh.finalize()
end
