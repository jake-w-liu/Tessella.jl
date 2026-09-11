#!/usr/bin/env julia
# P2/P5: $PartitionedEntities/$GhostElements structural metadata vs Gmsh 4.15.2.
#
# Gmsh partitions a meshed box with ghost cells and writes MSH 4.1 in ASCII
# and binary form. Tessella must parse the partition sections structurally —
# partition count, partitioned entity parent links and memberships, ghost
# entity pairs, and ghost element ownership — re-emit them into ASCII and
# binary MSH 4.1 without loss, and stay readable by Gmsh itself.

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

entity_equal(a,b)=all(getfield(a,f)==getfield(b,f)
                     for f in fieldnames(Elements.MixedEntity))
ghost_equal(a,b)=all(getfield(a,f)==getfield(b,f)
                    for f in fieldnames(Elements.MixedGhostElement))

function partition_equal(a,b)
    (a===nothing)!=(b===nothing) && return false
    a===nothing && return true
    a.num_partitions==b.num_partitions || return false
    a.ghost_entities==b.ghost_entities || return false
    keys(a.entities)==keys(b.entities) || return false
    all(entity_equal(a.entities[k],b.entities[k]) for k in keys(a.entities)) ||
        return false
    length(a.ghost_elements)==length(b.ghost_elements) || return false
    return all(ghost_equal(x,y) for (x,y) in
               zip(a.ghost_elements,b.ghost_elements))
end

function gmsh_summary(path)
    gmsh.clear()
    gmsh.open(path)
    entities=length(gmsh.model.getEntities())
    partitions=gmsh.model.getNumberOfPartitions()
    nodes=length(gmsh.model.mesh.getNodes()[1])
    _,etags,_=gmsh.model.mesh.getElements()
    elements=sum(length.(etags))
    return (entities,partitions,nodes,elements)
end

mktempdir() do dir
    gmsh.initialize(["gmsh","-v","0"])
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "partition metadata parity requires Gmsh API 4.15.2, " *
        "got $(gmsh.GMSH_API_VERSION)")
    try
        gmsh.model.add("partitioned")
        gmsh.model.occ.addBox(0,0,0,1,1,1)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMax",0.35)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.partition(3)
        paths=Dict{Symbol,String}()
        for (binary,key) in ((false,:ascii),(true,:binary))
            gmsh.option.setNumber("Mesh.Binary",binary ? 1 : 0)
            path=joinpath(dir,"partitioned-$key.msh")
            gmsh.write(path)
            paths[key]=path
        end
        for (key,path) in paths
            reference=gmsh_summary(path)
            mesh=Elements.read_mixed_msh(path)
            Elements.validate(mesh).ok ||
                error("$key: read mesh failed structural validation")
            partition=mesh.partition_data
            partition===nothing &&
                error("$key: partition metadata not parsed")
            partition.num_partitions==3 ||
                error("$key: partition count $(partition.num_partitions) != 3")
            isempty(partition.entities) &&
                error("$key: no partitioned entities parsed")
            for (key2,entity) in partition.entities
                key2==(entity.dim,Int(entity.tag)) ||
                    error("$key: partitioned entity key mismatch $key2")
                entity.parent[1]>=0 ||
                    error("$key: partitioned entity $key2 lost its parent")
            end
            for (version,binary,tag) in
                ((4.1,false,"v4"),(4.1,true,"v4b"))
                out=joinpath(dir,"out-$key-$tag.msh")
                Elements.write_mixed_msh(out,mesh;version=version,binary=binary)
                again=Elements.read_mixed_msh(out)
                Elements.validate(again).ok ||
                    error("$key/$tag: re-emitted mesh failed validation")
                partition_equal(again.partition_data,partition) ||
                    error("$key/$tag: partition metadata changed on round trip")
                summary=gmsh_summary(out)
                summary==reference || error(
                    "$key/$tag: Gmsh reads $summary != $reference")
            end
        end
    finally
        gmsh.finalize()
    end
end
println("msh_partition_metadata: all probes matched across ascii/binary " *
        "sources and v4.1 ascii/binary outputs")
