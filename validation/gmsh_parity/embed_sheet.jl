#!/usr/bin/env julia
# P6: Tessella Surface-In-Volume open sheet vs analytic volume 1 and vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: ntets, nnodes, tet_volume, node
using Tessella.Mesh3D: mesh_covers_triangle3
using Tessella.Elements: mixed_crc, read_mixed_msh, write_mixed_msh

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
projected=model_to_mixed(result.model,mesh,3,1)
validate(projected).ok || error("Tessella embedded-sheet projection is invalid")
surface_block=only(findall(block->block.msh==2,projected.blocks))
Int32(101) in projected.entity_data.block_entities[surface_block] ||
    error("Tessella embedded-sheet projection has no Surface[101] elements")
haskey(projected.entity_data.entities,(3,1)) ||
    error("Tessella embedded-sheet projection has no Volume[1] entity")
projected_crc=mixed_crc(projected)
projected_crc.sha=="2ccb9e42be0322ab810c601a99527b241c7e02d7269aa8c51b306f61347b6027" ||
    error("Tessella embedded-sheet projection CRC changed: $(projected_crc.sha)")

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
    gmsh.model.mesh.getEmbedded(3,1)==Tuple{Int32,Int32}[(2,101)] ||
        error("Gmsh source model did not retain Surface[101] In Volume[1]")

    mktempdir() do directory
        gmsh.option.setNumber("Mesh.SaveAll",1)
        gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        msh2_crcs=Set{String}()
        for binary in (false,true)
            gmsh.option.setNumber("Mesh.Binary",binary ? 1 : 0)
            native_path=joinpath(directory,"gmsh-sheet-$binary.msh")
            gmsh.write(native_path)
            native=read_mixed_msh(native_path)
            haskey(native.entity_data.entities,(2,101)) ||
                error("Tessella did not read Gmsh's Surface[101] entity")
            haskey(native.entity_data.entities,(3,1)) ||
                error("Tessella did not read Gmsh's Volume[1] entity")

            for version in (2.2,4.1)
                projected_path=joinpath(
                    directory,"tessella-sheet-$version-$binary.msh")
                write_mixed_msh(
                    projected_path,projected;version=version,binary=binary)
                reread=read_mixed_msh(projected_path)
                if version==4.1
                    mixed_crc(reread)==projected_crc ||
                        error("Tessella embedded-sheet MSH4 $binary CRC mismatch")
                else
                    reread.entity_data===nothing ||
                        error("Tessella MSH2 unexpectedly retained MSH4 entities")
                    sheet_block=only(findall(block->block.msh==2,reread.blocks))
                    Int32(101) in reread.elementary_entities[sheet_block] ||
                        error("Tessella MSH2 lost Surface[101] ownership")
                    push!(msh2_crcs,mixed_crc(reread).sha)
                end
                gmsh.clear()
                gmsh.open(projected_path)
                isempty(gmsh.model.mesh.getEmbedded(3,1)) ||
                    error("Gmsh reconstructed a nonexistent volume-embedding record")
                _,sheet_elements,_=gmsh.model.mesh.getElements(2,101)
                sum(length,sheet_elements;init=0)>0 ||
                    error("Gmsh did not retain projected Surface[101] elements")
                _,volume_elements,_=gmsh.model.mesh.getElements(3,1)
                sum(length,volume_elements;init=0)>0 ||
                    error("Gmsh did not retain projected Volume[1] elements")
            end

            gmsh.clear()
            gmsh.open(GEO)
            gmsh.option.setNumber("Mesh.MeshSizeMax",0.5)
            gmsh.model.mesh.generate(3)
            gmsh.option.setNumber("Mesh.SaveAll",1)
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        end
        length(msh2_crcs)==1 || error(
            "Tessella embedded-sheet MSH2 modes produced different CRCs: $msh2_crcs")
        only(msh2_crcs)==
            "8141bdc49a658a1770944e32527e2971b598f56ea2b0a4bae27e128bee0da640" ||
            error("Tessella embedded-sheet MSH2 CRC changed: $(only(msh2_crcs))")
    end
    println("GMSH_PARITY_EMBED_SHEET_OK gmsh=$(gmsh.GMSH_API_VERSION) " *
            "tessella_volume=1 gmsh_tets=$gmsh_tets tessella_tets=$(ntets(mesh)) " *
            "tessella_nodes=$(nnodes(mesh)) projected_crc=$(projected_crc.sha) " *
            "msh2_msh4_ascii_binary=ok serialized_volume_relation=absent")
finally
    gmsh.finalize()
end
