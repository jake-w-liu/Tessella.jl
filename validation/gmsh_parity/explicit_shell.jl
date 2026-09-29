#!/usr/bin/env julia
# P6: explicit planar Surface Loop/Volume lifecycle vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: nnodes, ntets, node, tet_volume
using Tessella.Elements: mixed_crc, read_mixed_msh, write_mixed_msh

const GEO=joinpath(@__DIR__,"explicit_shell.geo")
const EXPECTED_BOUNDARIES=Int32[1,2,3,4,5,6]

result=execute_geo(GEO;mesh_dim=3)
result.model.surface_loops[1]==collect(1:6) ||
    error("Tessella explicit shell did not retain Surface Loop[1]")
result.model.volumes[1]==[1] ||
    error("Tessella explicit shell did not retain Volume[1]")
mesh=result.mesh
mesh===nothing && error("Tessella explicit shell produced no mesh")
ntets(mesh)>0 || error("Tessella explicit shell has no tetrahedra")
tessella_volume=sum(tet_volume(
    node(mesh,mesh.tets[1,cell]),node(mesh,mesh.tets[2,cell]),
    node(mesh,mesh.tets[3,cell]),node(mesh,mesh.tets[4,cell]))
    for cell in 1:ntets(mesh))
abs(tessella_volume-1)<=1e-12 ||
    error("Tessella explicit shell volume $tessella_volume != 1")
projected=model_to_mixed(result.model,mesh,3,1)
validate(projected).ok || error("Tessella explicit-shell projection is invalid")
projected.entity_data.entities[(3,1)].boundaries==EXPECTED_BOUNDARIES ||
    error("Tessella explicit-shell projection lost volume boundaries")
projected_crc=mixed_crc(projected)
projected_crc.sha==
    "9bce88e319c67236317df64b876739a62f80982ed86eca028bd1e7bda022bcb6" ||
    error("Tessella explicit-shell projection CRC changed: $(projected_crc.sha)")

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
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "explicit-shell differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    gmsh_boundary=gmsh.model.getBoundary([(3,1)],false,true,false)
    gmsh_boundary==[(Int32(2),surface) for surface in EXPECTED_BOUNDARIES] ||
        error("Gmsh explicit Volume[1] boundary mismatch: $gmsh_boundary")
    gmsh.model.mesh.generate(3)
    types,element_tags,element_nodes=gmsh.model.mesh.getElements(3,1)
    types==Int32[4] || error("Gmsh explicit shell emitted element types $types")
    length(element_tags)==1 && !isempty(only(element_tags)) ||
        error("Gmsh explicit shell has no volume elements")
    node_tags,node_coordinates,_=gmsh.model.mesh.getNodes()
    coordinates=Dict{UInt64,NTuple{3,Float64}}()
    for (index,tag) in pairs(node_tags)
        offset=3*index-2
        coordinates[UInt64(tag)]=(node_coordinates[offset],
                                  node_coordinates[offset+1],
                                  node_coordinates[offset+2])
    end
    gmsh_volume=0.0
    connectivity=only(element_nodes)
    length(connectivity)%4==0 || error("Gmsh tetrahedron connectivity is malformed")
    for offset in 1:4:length(connectivity)
        points=ntuple(slot->coordinates[UInt64(connectivity[offset+slot-1])],4)
        gmsh_volume+=tet_volume(points...)
    end
    abs(gmsh_volume-1)<=1e-10 ||
        error("Gmsh explicit shell volume $gmsh_volume != 1")

    mktempdir() do directory
        gmsh.option.setNumber("Mesh.SaveAll",1)
        gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        msh2_crcs=Set{String}()
        for binary in (false,true)
            gmsh.option.setNumber("Mesh.Binary",binary ? 1 : 0)
            native_path=joinpath(directory,"gmsh-explicit-$binary.msh")
            gmsh.write(native_path)
            native=read_mixed_msh(native_path)
            native.entity_data.entities[(3,1)].boundaries==EXPECTED_BOUNDARIES ||
                error("Tessella did not read Gmsh's explicit volume boundaries")

            for version in (2.2,4.1)
                projected_path=joinpath(
                    directory,"tessella-explicit-$version-$binary.msh")
                write_mixed_msh(
                    projected_path,projected;version=version,binary=binary)
                reread=read_mixed_msh(projected_path)
                if version==4.1
                    mixed_crc(reread)==projected_crc ||
                        error("Tessella explicit-shell MSH4 $binary CRC mismatch")
                    reread.entity_data.entities[(3,1)].boundaries==
                        EXPECTED_BOUNDARIES ||
                        error("Tessella MSH4 lost explicit volume boundaries")
                else
                    reread.entity_data===nothing ||
                        error("Tessella MSH2 unexpectedly retained MSH4 entities")
                    push!(msh2_crcs,mixed_crc(reread).sha)
                end

                gmsh.clear()
                gmsh.open(projected_path)
                reopened_boundary=gmsh.model.getBoundary(
                    [(3,1)],false,true,false)
                expected_reopened=version==4.1 ?
                    [(Int32(2),surface) for surface in EXPECTED_BOUNDARIES] :
                    Tuple{Int32,Int32}[]
                reopened_boundary==expected_reopened || error(
                    "Gmsh boundary mismatch for MSH$version binary=$binary: " *
                    "$reopened_boundary")
                for (dim,tags) in ((0,1:8),(1,1:12),(2,1:6),(3,1:1)),tag in tags
                    _,cells,_=gmsh.model.mesh.getElements(dim,tag)
                    sum(length,cells;init=0)>0 || error(
                        "Gmsh lost projected entity ($dim,$tag) elements for " *
                        "MSH$version binary=$binary")
                end
            end

            gmsh.clear()
            gmsh.open(GEO)
            gmsh.model.mesh.generate(3)
            gmsh.option.setNumber("Mesh.SaveAll",1)
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        end
        length(msh2_crcs)==1 || error(
            "Tessella explicit-shell MSH2 modes produced different CRCs: " *
            "$msh2_crcs")
        only(msh2_crcs)==
            "2bbf8be73f7f4b3204327f03b1247334e1a141483874b87238bcd783749e7c23" || error(
            "Tessella explicit-shell MSH2 CRC changed: $(only(msh2_crcs))")
    end
    println("GMSH_PARITY_EXPLICIT_SHELL_OK gmsh=$(gmsh.GMSH_API_VERSION) " *
            "tessella_volume=$tessella_volume gmsh_volume=$gmsh_volume " *
            "tessella_tets=$(ntets(mesh)) tessella_nodes=$(nnodes(mesh)) " *
            "gmsh_tets=$(length(only(element_tags))) " *
            "projected_crc=$(projected_crc.sha) msh2_msh4_ascii_binary=ok")
finally
    gmsh.finalize()
end
