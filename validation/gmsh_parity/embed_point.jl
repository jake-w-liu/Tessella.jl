#!/usr/bin/env julia
# P6: Tessella Point-In-Surface embed vs analytic area 1 and vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: ntris, nnodes, node
using Tessella.Elements: mixed_crc, read_mixed_msh, write_mixed_msh

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
projected=model_to_mixed(result.model,mesh,1)
validate(projected).ok || error("Tessella embedded-point projection is invalid")
point_block=only(findall(block->block.msh==15,projected.blocks))
Int32(5) in projected.entity_data.block_entities[point_block] ||
    error("Tessella embedded-point projection has no Point[5] element")
projected.entity_data.entities[(2,1)].embedded_curves==Int32[] ||
    error("Tessella encoded a point as an embedded curve")
projected_crc=mixed_crc(projected)
projected_crc.sha=="222619f8e92298ab72ece09cae6dd9f8300781c9d8a7dc587fc5b3de50219fc3" ||
    error("Tessella embedded-point projection CRC changed: $(projected_crc.sha)")

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
    gmsh.model.mesh.getEmbedded(2,1)==Tuple{Int32,Int32}[(0,5)] ||
        error("Gmsh source model did not retain Point[5] In Surface[1]")

    mktempdir() do directory
        for binary in (false,true)
            projected_path=joinpath(directory,"tessella-point-$binary.msh")
            write_mixed_msh(
                projected_path,projected;version=4.1,binary=binary)
            reread=read_mixed_msh(projected_path)
            mixed_crc(reread)==projected_crc ||
                error("Tessella embedded-point MSH4 $binary CRC mismatch")
            gmsh.clear()
            gmsh.open(projected_path)
            gmsh.model.mesh.getEmbedded(2,1)==Tuple{Int32,Int32}[] ||
                error("Gmsh reconstructed a nonexistent point embedding")
            _,point_elements,_=gmsh.model.mesh.getElements(0,5)
            sum(length,point_elements;init=0)==1 ||
                error("Gmsh did not retain the projected Point[5] element")
        end
    end
    println("GMSH_PARITY_EMBED_OK gmsh=$(gmsh.GMSH_API_VERSION) " *
            "tessella_area=1 gmsh_tris=$gmsh_tris tessella_tris=$(ntris(mesh)) " *
            "tessella_nodes=$(nnodes(mesh)) projected_crc=$(projected_crc.sha) " *
            "msh4_ascii_binary=ok")
finally
    gmsh.finalize()
end
