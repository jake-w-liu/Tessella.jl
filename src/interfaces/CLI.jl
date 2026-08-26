"""
    CLI

Command-line façade: `tessella file.geo -2|-3` executes the bounded `.geo`
subset and writes `file.msh`. Supported periodic and embedded surface entities
and embedded volume entities, including nested sheet constraints and explicit
planar volume shells, are written with classified MSH4 metadata. Unknown flags
and OCC-only files are blockers.
"""
module CLI

using ..GeoExec: execute_geo
using ..IO: write_msh
using ..Model: model_periodic_constraints, model_to_mixed
using ..Elements: write_mixed_msh

export main

const _MAX_ARGUMENTS=10_000
const _MAX_ARGUMENT_BYTES=1_000_000

function _copy_arguments(args)
    tokens=String[];sizehint!(tokens,length(args))
    total=0
    for raw in args
        total=try
            Base.checked_add(total,ncodeunits(raw))
        catch err
            err isa InterruptException && rethrow()
            err isa OverflowError || rethrow()
            throw(ArgumentError("tessella: command-line argument size overflow"))
        end
        total<=_MAX_ARGUMENT_BYTES || throw(ArgumentError(
            "tessella: command-line arguments exceed $_MAX_ARGUMENT_BYTES bytes"))
        push!(tokens,String(raw))
    end
    return tokens
end

function _path_argument(value::AbstractString,what::AbstractString)
    path=String(value)
    isempty(path) && throw(ArgumentError("tessella: $what path must not be empty"))
    occursin('\0',path) && throw(ArgumentError("tessella: $what path contains a NUL byte"))
    return path
end

function _default_output(path::String)
    stem,extension=splitext(path)
    return isempty(extension) ? path*".msh" : stem*".msh"
end

function _same_file_target(a::String,b::String)
    normpath(abspath(a))==normpath(abspath(b)) && return true
    return ispath(a) && ispath(b) && samefile(a,b)
end

"""
    main(args) -> String

Execute one `.geo` input with optional `-2` or `-3` meshing and `-o output`.
Parsing-only mode returns the input path. Meshing writes an atomic MSH 4.1 file
and returns its path. Native periodic/embedded surface meshes and embedded volume
meshes include their supported classified entity, cell, and node records, including
nested point/curve constraints on an embedded sheet and explicit volume boundaries.
Periodic volume output is
blocked; periodic surface output requires exactly one selected surface containing
every slave curve. Duplicate/conflicting flags, multiple inputs, ignored output
arguments, and any output that aliases the input are rejected.
"""
function main(args::AbstractVector{<:AbstractString})
    length(args)<=_MAX_ARGUMENTS || throw(ArgumentError(
        "tessella: more than $_MAX_ARGUMENTS command-line arguments"))
    isempty(args) && throw(ArgumentError("tessella: missing input file"))
    tokens=_copy_arguments(args)
    path=nothing;dim=0;outfile=nothing
    dimension_seen=false;output_seen=false
    i=1
    while i<=length(tokens)
        a=tokens[i]
        if a in ("-2","-3")
            dimension_seen && throw(ArgumentError(
                "tessella: duplicate or conflicting dimension flag $a"))
            dim=a=="-2" ? 2 : 3
            dimension_seen=true
        elseif a=="-o"
            output_seen && throw(ArgumentError("tessella: duplicate -o flag"))
            i+=1;i<=length(tokens) || throw(ArgumentError("tessella: -o needs a path"))
            startswith(tokens[i],'-') && throw(ArgumentError(
                "tessella: -o needs a path, got flag $(tokens[i])"))
            outfile=_path_argument(tokens[i],"output")
            output_seen=true
        elseif startswith(a,"-")
            throw(ArgumentError("tessella: unsupported flag $a"))
        else
            path===nothing || throw(ArgumentError("tessella: multiple input files"))
            path=_path_argument(a,"input")
        end
        i+=1
    end
    path===nothing && throw(ArgumentError("tessella: missing input file"))
    dim==0 && outfile!==nothing && throw(ArgumentError(
        "tessella: -o requires -2 or -3"))
    destination=dim==0 ? nothing : (outfile===nothing ? _default_output(path) : outfile)
    destination!==nothing && _same_file_target(path,destination) && throw(ArgumentError(
        "tessella: output path must not alias the input file"))

    result=execute_geo(path;mesh_dim=dim)
    if dim>0
        result.mesh===nothing && throw(ErrorException("tessella: no mesh produced"))
        constraints=model_periodic_constraints(result.model)
        if !isempty(constraints) && dim!=2
            throw(ArgumentError(
                "tessella: periodic metadata projection is limited to surface meshes"))
        end
        target=nothing
        embedded=any(
            pair->pair.first[1]==dim && !isempty(pair.second),
            pairs(result.model.embeds))
        targets=dim==2 ? result.model.surfaces : result.model.volumes
        if length(targets)==1
            target=only(keys(targets))
            embedded=!isempty(get(
                result.model.embeds,(dim,target),NTuple{2,Int}[]))
        end
        explicit_shell=dim==3 && target!==nothing &&
                       !isempty(result.model.volumes[target])
        if isempty(constraints) && !embedded && !explicit_shell
            write_msh(destination,result.mesh;version=4.1)
        else
            entity_name=dim==2 ? "surface" : "volume"
            target===nothing && throw(ArgumentError(
                "tessella: classified $entity_name projection requires exactly " *
                "one selected $entity_name"))
            projected=model_to_mixed(result.model,result.mesh,dim,target)
            if dim==2
                projected_slaves=Set(Int(link.slave_entity)
                    for link in projected.periodic_links if link.dim==1)
                missing_slaves=sort!(Int[constraint.slave_entity
                    for constraint in constraints
                    if !(Int(constraint.slave_entity) in projected_slaves)])
                isempty(missing_slaves) || throw(ArgumentError(
                    "tessella: selected Surface[$target] does not contain periodic " *
                    "slave Curve tags $missing_slaves"))
            end
            write_mixed_msh(destination,projected;version=4.1)
        end
        return destination
    end
    return path
end

end # module
