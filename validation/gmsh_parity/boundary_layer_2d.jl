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
function quadrangle_area(m)
    block=only(b for b in m.blocks if b.msh==3)
    area=0.0
    for cell in axes(block.nodes,2)
        p1=ntuple(d->m.coords[d,block.nodes[1,cell]],3)
        p2=ntuple(d->m.coords[d,block.nodes[2,cell]],3)
        p3=ntuple(d->m.coords[d,block.nodes[3,cell]],3)
        p4=ntuple(d->m.coords[d,block.nodes[4,cell]],3)
        area+=triangle_area(p1,p2,p4)+triangle_area(p4,p2,p3)
    end
    return area
end
tessella_area=quadrangle_area(mesh)
abs(tessella_area-H)<=1e-12 || error(
    "Tessella 2-D BL area $tessella_area != $H")
tessella_quads==NL || error("Tessella 2-D BL quad count $tessella_quads != $NL")

# Independent rigid-motion oracle for the general-plane path. A rotated and
# translated copy must preserve all topology and transform every node of the
# historical x/y result.
ex=(inv(sqrt(2.0)),inv(sqrt(2.0)),0.0)
ey=(-inv(sqrt(6.0)),inv(sqrt(6.0)),2inv(sqrt(6.0)))
normal=(ex[2]*ey[3]-ex[3]*ey[2],
        ex[3]*ey[1]-ex[1]*ey[3],
        ex[1]*ey[2]-ex[2]*ey[1])
origin=(1.25,-2.5,3.75)
rigid(point)=(origin[1]+point[1]*ex[1]+point[2]*ey[1],
              origin[2]+point[1]*ex[2]+point[2]*ey[2],
              origin[3]+point[1]*ex[3]+point[2]*ey[3])
tilted_endpoints=(rigid((0.0,0.0,0.0)),rigid((1.0,0.0,0.0)))
tilted_coords=hcat(collect.(tilted_endpoints)...)
tilted_edge=Mesh(tilted_coords;segs=segs)
tilted=mesh_boundary_layer_2d(
    tilted_edge;hwall=HW,ratio=RA,nlayers=NL,plane_normal=normal)
length(tilted.blocks)==length(mesh.blocks) || error(
    "general-plane boundary layer changed block count")
all(first.msh==second.msh && first.nodes==second.nodes
    for (first,second) in zip(tilted.blocks,mesh.blocks)) || error(
        "general-plane boundary layer changed topology")
function rigid_node_error(source,destination,transform)
    maximum_error=0.0
    for index in axes(source.coords,2)
        expected=transform((source.coords[1,index],source.coords[2,index],
                            source.coords[3,index]))
        for coordinate in 1:3
            maximum_error=max(
                maximum_error,
                abs(destination.coords[coordinate,index]-expected[coordinate]))
        end
    end
    return maximum_error
end
tilted_max_error=rigid_node_error(mesh,tilted,rigid)
tilted_max_error<=2eps(Float64) || error(
    "general-plane rigid-motion node error $tilted_max_error")
abs(quadrangle_area(tilted)-H)<=128eps(Float64) || error(
    "general-plane boundary-layer area changed")
vertical=Mesh(Float64[1 1;2 2;3 5];segs=segs)
vertical_mesh=mesh_boundary_layer_2d(
    vertical;hwall=0.125,ratio=1.5,nlayers=2,plane_normal=(1,0,0))
general_plane_crc=Tessella.Elements.mixed_crc(vertical_mesh).sha
general_plane_crc==
    "49e70bbffb8fd2121a526c35a5ca51ab87ea19790bd715d52a147a21535b3e19" ||
    error("general-plane boundary-layer CRC changed: $general_plane_crc")

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
    println("GMSH_PARITY_BL2D_OK gmsh=$(gmsh.GMSH_API_VERSION) "*
            "tessella_area=$H gmsh_quads=$gmsh_quads gmsh_tris=$gmsh_tris "*
            "tessella_quads=$tessella_quads general_plane_max_error=$tilted_max_error "*
            "general_plane_crc=$general_plane_crc")
finally
    gmsh.finalize()
end
