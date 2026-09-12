#!/usr/bin/env julia
# P2: MSH 4.1 binary size_t width vs Gmsh 4.15.2.
#
# Gmsh writes binary MSH with an 8-byte size_t on 64-bit builds and its
# reader rejects the 4-byte data-size word outright. Tessella decodes both
# widths on input and can emit the narrower encoding for Tessella-only
# serialization via write_mixed_msh(...; binary=true, gmsh_compatible=false,
# size_t_bytes=4). This probe checks that a narrow write round-trips
# losslessly through Tessella's reader and that pinned Gmsh still rejects it
# — the exact boundary the gmsh_compatible gate encodes.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
const Elements = Tessella.Elements

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

function mesh_signature(mesh)
    return (mesh.coords,
            [(b.msh,copy(b.nodes)) for b in mesh.blocks],
            mesh.entity_data===nothing ? nothing :
                (mesh.entity_data.external_node_tags,
                 mesh.entity_data.external_element_tags),
            [(l.dim,l.slave_entity,l.master_entity,
              l.slave_nodes,l.master_nodes) for l in mesh.periodic_links])
end

mktempdir() do dir
    gmsh.initialize(["gmsh","-v","0"])
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "size_t width parity requires Gmsh API 4.15.2, " *
        "got $(gmsh.GMSH_API_VERSION)")
    mesh=try
        gmsh.model.add("width")
        gmsh.model.occ.addBox(0,0,0,1,1,1)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMax",0.4)
        gmsh.model.mesh.generate(3)
        gmsh.option.setNumber("Mesh.Binary",1)
        source=joinpath(dir,"gmsh-wide.msh")
        gmsh.write(source)
        read=Elements.read_mixed_msh(source)
        Elements.validate(read).ok ||
            error("Gmsh wide binary mesh failed structural validation")

        # Narrow write: Tessella-only serialization, lossless on re-read.
        narrow=joinpath(dir,"tessella-narrow.msh")
        Elements.write_mixed_msh(narrow,read;binary=true,
                                 gmsh_compatible=false,size_t_bytes=4)
        open(narrow,"r") do io
            readline(io)=="\$MeshFormat" ||
                error("narrow output is missing \$MeshFormat")
            split(readline(io))==["4.1","1","4"] ||
                error("narrow output did not declare a 4-byte data size")
        end
        again=Elements.read_mixed_msh(narrow)
        Elements.validate(again).ok ||
            error("narrow round trip failed structural validation")
        mesh_signature(again)==mesh_signature(read) ||
            error("narrow round trip changed the mesh")

        # Pinned Gmsh must reject the narrow file: its 64-bit reader requires
        # sizeof(size_t)=8, which is why the writer gates this behind
        # gmsh_compatible=false.
        rejected=false
        try
            gmsh.clear()
            gmsh.open(narrow)
        catch
            rejected=true
        end
        rejected || error("Gmsh 4.15.2 accepted a 4-byte size_t file")

        wide=joinpath(dir,"tessella-wide.msh")
        Elements.write_mixed_msh(wide,read;binary=true)
        read
    finally
        gmsh.finalize()
    end

    # A failed open poisons Gmsh's model state, so verify the default-width
    # file loads identically in a fresh session.
    gmsh.initialize(["gmsh","-v","0"])
    try
        gmsh.open(joinpath(dir,"tessella-wide.msh"))
        _,etags,_=gmsh.model.mesh.getElements()
        gmsh_total=sum(length.(etags))
        tessella_total=sum(Elements._block_ncells(b) for b in mesh.blocks)
        gmsh_total==tessella_total || error(
            "wide binary output element count changed under Gmsh: " *
            "$gmsh_total != $tessella_total")
    finally
        gmsh.finalize()
    end
end
println("msh_size_t_width: narrow round trip lossless; Gmsh 4.15.2 rejects " *
        "the 4-byte width and reads the 8-byte default")
