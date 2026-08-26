#!/usr/bin/env julia
# P6: Tessella Line-In-Surface embed vs analytic area 1 and vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: ntris, nnodes, node
using Tessella.Elements: mixed_crc, read_mixed_msh, write_mixed_msh

const GEO=joinpath(@__DIR__,"embed_line.geo")
result=execute_geo(GEO; mesh_dim=2)
mesh=result.mesh
mesh===nothing && error("Tessella embed-line produced no mesh")
Tessella.Model._mesh_covers_segment(mesh,(0.25,0.5,0.0),(0.75,0.5,0.0)) ||
    error("Tessella embed-line is not a chain of mesh edges")
area=sum(begin
    a=node(mesh,mesh.tris[1,t]); b=node(mesh,mesh.tris[2,t]); c=node(mesh,mesh.tris[3,t])
    abs((b[1]-a[1])*(c[2]-a[2])-(c[1]-a[1])*(b[2]-a[2]))/2
end for t in 1:ntris(mesh))
abs(area-1)<=1e-12 || error("Tessella embed-line area $area != 1")
projected=model_to_mixed(result.model,mesh,1)
validate(projected).ok || error("Tessella embedded-line projection is invalid")
projected.entity_data.entities[(2,1)].embedded_curves==Int32[5] ||
    error("Tessella projection did not attach Curve[5] to Surface[1]")
line_block=only(findall(block->block.msh==1,projected.blocks))
Int32(5) in projected.entity_data.block_entities[line_block] ||
    error("Tessella embedded-line projection has no Curve[5] elements")
projected_crc=mixed_crc(projected)
projected_crc.sha=="0655fe3edb4344be584d2fe12b8d57637f65090524542d2bc1b515e148d55ea5" ||
    error("Tessella embedded-line projection CRC changed: $(projected_crc.sha)")

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
    function has_node(x,y)
        for i in 1:3:length(coord)
            hypot(coord[i]-x, coord[i+1]-y)<=1e-9 && return true
        end
        return false
    end
    has_node(0.25,0.5) || error("Gmsh embed-line missing node (0.25,0.5)")
    has_node(0.75,0.5) || error("Gmsh embed-line missing node (0.75,0.5)")
    types, tags, _ = gmsh.model.mesh.getElements(2)
    gmsh_tris=sum(length, tags; init=0)
    gmsh_tris>0 || error("Gmsh embed-line has no surface elements")
    gmsh.model.mesh.getEmbedded(2,1)==Tuple{Int32,Int32}[(1,5)] ||
        error("Gmsh source model did not retain Curve[5] In Surface[1]")

    mktempdir() do directory
        gmsh.option.setNumber("Mesh.SaveAll",1)
        gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        for binary in (false,true)
            gmsh.option.setNumber("Mesh.Binary",binary ? 1 : 0)
            native_path=joinpath(directory,"gmsh-line-$binary.msh")
            gmsh.write(native_path)
            native=read_mixed_msh(native_path)
            native.entity_data.entities[(2,1)].embedded_curves==Int32[5] ||
                error("Tessella did not decode Gmsh's embedded-curve $binary record")

            projected_path=joinpath(directory,"tessella-line-$binary.msh")
            write_mixed_msh(
                projected_path,projected;version=4.1,binary=binary)
            reread=read_mixed_msh(projected_path)
            mixed_crc(reread)==projected_crc ||
                error("Tessella embedded-line MSH4 $binary CRC mismatch")
            gmsh.clear()
            gmsh.open(projected_path)
            expected_embedding=binary ? Tuple{Int32,Int32}[] :
                                        Tuple{Int32,Int32}[(1,5)]
            gmsh.model.mesh.getEmbedded(2,1)==expected_embedding ||
                error("Gmsh embedded-curve $binary reopen result differs")
            boundary=gmsh.model.getBoundary([(2,1)],false,true,false)
            boundary==Tuple{Int32,Int32}[(1,1),(1,2),(1,3),(1,4)] ||
                error("Gmsh treated Curve[5] as a surface boundary")
            _,line_elements,_=gmsh.model.mesh.getElements(1,5)
            sum(length,line_elements;init=0)>0 ||
                error("Gmsh did not retain projected Curve[5] elements")

            # Reload the source geometry before Gmsh writes the next native
            # variant; opening an MSH replaces the active model.
            gmsh.clear()
            gmsh.open(GEO)
            gmsh.model.mesh.generate(2)
            gmsh.option.setNumber("Mesh.SaveAll",1)
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        end
    end
    println("GMSH_PARITY_EMBED_LINE_OK gmsh=$(gmsh.GMSH_API_VERSION) " *
            "tessella_area=1 gmsh_tris=$gmsh_tris tessella_tris=$(ntris(mesh)) " *
            "tessella_nodes=$(nnodes(mesh)) projected_crc=$(projected_crc.sha) " *
            "gmsh_ascii_relation=retained gmsh_binary_relation=upstream_gap " *
            "tessella_ascii_binary=lossless")
finally
    gmsh.finalize()
end
