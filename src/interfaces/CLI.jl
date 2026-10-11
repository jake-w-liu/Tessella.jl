"""
    CLI

Command-line façade: `tessella file.geo -1|-2|-3` executes the bounded `.geo`
subset and writes `file.msh`. Supported periodic and embedded surface entities
and embedded volume entities, including nested sheet constraints and periodic
planar boundaries of explicit volume shells, are written with classified MSH4
metadata. Unknown flags and OCC-only files are blockers.
"""
module CLI

using ..GeoExec: execute_geo, geo_entity_mesh
using ..IO: write_msh
using ..Model: model_periodic_constraints, model_to_mixed
using ..Elements: write_mixed_msh, MixedMesh, simplex_to_mixed
using ..GmshLibm: _with_win_cli_precision

export main

const _HELP="""
Usage: tessella file.geo [-1|-2|-3] [-o output.msh] [-format msh2|msh4] [-bin]
       [-clscale value] [-clmin value] [-clmax value] [-order int]
       tessella -help | -version

Mesh size bounds apply before the positive size factor. CLI numeric options
have priority over .geo assignments. Curve orders 1-10 are supported;
Surface/Volume generation supports the orders documented in README.md.
"""

const _NUMERIC_FLAGS=Dict("-clscale"=>"Mesh.MeshSizeFactor",
    "-clmin"=>"Mesh.MeshSizeMin", "-clmax"=>"Mesh.MeshSizeMax",
    "-order"=>"Mesh.ElementOrder")

function _numeric_argument(flag::String,raw::String)
    value=tryparse(Float64,raw)
    value!==nothing && isfinite(value) || throw(ArgumentError(
        "tessella: $flag needs a finite number"))
    if flag=="-order"
        isinteger(value) && 1<=value<=10 || throw(ArgumentError(
            "tessella: -order needs an integer from 1 to 10"))
    elseif flag=="-clmin"
        value>=0 || throw(ArgumentError("tessella: -clmin must be nonnegative"))
    else
        value>0 || throw(ArgumentError("tessella: $flag must be positive"))
    end
    return value
end

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

Execute one `.geo` input with optional `-1`, `-2` or `-3` meshing and `-o output`.
`-format msh2|msh4` selects MSH 2.2 or 4.1 (the default), and `-bin` selects
binary output. Output flags require a meshing dimension.
`-clscale`, `-clmin`, `-clmax` and `-order` set mesh sizing and element order;
these values retain priority over script assignments, including option aliases.
`-help`/`--help` and `-version`/`--version` print and return their text without
an input file. Informational flags must be used alone.
Accepted syntax and blockers follow [`Tessella.GeoExec.execute_geo`](@ref).
Parsing-only mode returns the input path. Meshing writes an atomic MSH file
and returns its path. Native periodic/embedded surface meshes and embedded volume
meshes include their supported classified entity, cell, and node records, including
nested point/curve constraints on an embedded sheet and explicit volume boundaries.
Periodic `-2` output retains periodic-curve graphs with reused masters and acyclic
chains. Periodic `-3` output supports planar boundary-surface relations on one
explicit volume shell and retains their induced point/curve relations. Stored
periodic volume relations are mesh-inert and unserialized, matching Gmsh 4.15.2.
Independent periodic curves in volumes remain blocked.
One-dimensional output preserves Point/Curve ownership, physical groups,
periodic node maps and the configured Curve element order.
Duplicate or conflicting flags, multiple inputs, ignored output arguments, and any
output that aliases the input are rejected.
On Windows x86_64, power expressions select the standalone Gmsh executable's
64-bit x87 arithmetic precision in a task-local execution context. Each native
power operation restores the complete caller control word before returning.
Direct API and `execute_geo` calls retain caller precision.
"""
function main(args::AbstractVector{<:AbstractString})
    return _with_win_cli_precision() do
        _main(args)
    end
end

function _main(args::AbstractVector{<:AbstractString})
    length(args)<=_MAX_ARGUMENTS || throw(ArgumentError(
        "tessella: more than $_MAX_ARGUMENTS command-line arguments"))
    isempty(args) && throw(ArgumentError("tessella: missing input file"))
    tokens=_copy_arguments(args)
    information=filter(a->a in ("-help","--help","-version","--version"),tokens)
    if !isempty(information)
        length(tokens)==1 || throw(ArgumentError(
            "tessella: informational flags must be used alone"))
        message=only(tokens) in ("-help","--help") ? _HELP :
            string(Base.pkgversion(parentmodule(@__MODULE__)))
        println(message)
        return message
    end
    numeric_options=Dict{String,Float64}()
    path=nothing;dim=0;outfile=nothing;version=4.1;binary=false
    dimension_seen=false;output_seen=false
    format_seen=false;binary_seen=false
    i=1
    while i<=length(tokens)
        a=tokens[i]
        if a in ("-1","-2","-3")
            dimension_seen && throw(ArgumentError(
                "tessella: duplicate or conflicting dimension flag $a"))
            dim=a=="-1" ? 1 : a=="-2" ? 2 : 3
            dimension_seen=true
        elseif a=="-o"
            output_seen && throw(ArgumentError("tessella: duplicate -o flag"))
            i+=1;i<=length(tokens) || throw(ArgumentError("tessella: -o needs a path"))
            startswith(tokens[i],'-') && throw(ArgumentError(
                "tessella: -o needs a path, got flag $(tokens[i])"))
            outfile=_path_argument(tokens[i],"output")
            output_seen=true
        elseif a=="-format"
            format_seen && throw(ArgumentError("tessella: duplicate -format flag"))
            i+=1;i<=length(tokens) || throw(ArgumentError("tessella: -format needs msh2 or msh4"))
            format=tokens[i]
            format in ("msh2","msh4") || throw(ArgumentError(
                "tessella: unsupported output format $format; supported: msh2, msh4"))
            version=format=="msh2" ? 2.2 : 4.1
            format_seen=true
        elseif haskey(_NUMERIC_FLAGS,a)
            option=_NUMERIC_FLAGS[a]
            haskey(numeric_options,option) && throw(ArgumentError(
                "tessella: duplicate $a flag"))
            i+=1;i<=length(tokens) || throw(ArgumentError("tessella: $a needs a value"))
            numeric_options[option]=_numeric_argument(a,tokens[i])
        elseif a=="-bin"
            binary_seen && throw(ArgumentError("tessella: duplicate -bin flag"))
            binary=true;binary_seen=true
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
        "tessella: -o requires -1, -2 or -3"))
    dim==0 && (format_seen || binary_seen) && throw(ArgumentError(
        "tessella: output format flags require -1, -2 or -3"))
    destination=dim==0 ? nothing : (outfile===nothing ? _default_output(path) : outfile)
    destination!==nothing && _same_file_target(path,destination) && throw(ArgumentError(
        "tessella: output path must not alias the input file"))

    get(numeric_options,"Mesh.MeshSizeMin",0.0)<=
        get(numeric_options,"Mesh.MeshSizeMax",1e22) || throw(ArgumentError(
            "tessella: minimum mesh size exceeds maximum"))
    result=execute_geo(path;mesh_dim=dim,_forced_numeric_options=numeric_options)
    if dim>0
        result.mesh===nothing && throw(ErrorException("tessella: no mesh produced"))
        if dim==1
            authority=result.model.meshing.mesh_authority
            authority!==nothing && authority.mesh!==nothing || throw(ErrorException(
                "tessella: one-dimensional generation omitted its classified mesh"))
            write_mixed_msh(destination,authority.mesh;version,binary)
            return destination
        end
        constraints=model_periodic_constraints(result.model)
        target=nothing
        # The parsed dimension is settled before closure capture. Capturing
        # mutable dim would box every use in this execution function.
        mesh_dimension=dim
        embedded=any(
            pair->pair.first[1]==mesh_dimension && !isempty(pair.second),
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
            if result.mesh isa MixedMesh
                write_mixed_msh(destination,result.mesh;version,binary)
            elseif binary
                write_mixed_msh(destination,simplex_to_mixed(result.mesh);version,binary)
            else
                write_msh(destination,result.mesh;version)
            end
        else
            entity_name=dim==2 ? "surface" : "volume"
            target===nothing && throw(ArgumentError(
                "tessella: classified $entity_name projection requires exactly " *
                "one selected $entity_name"))
            # `model_to_mixed` projects a single entity's own mesh — the
            # merged execution product carries every entity's cells, so the
            # target's part is the correct input.
            projected=model_to_mixed(result.model,
                geo_entity_mesh(result,dim,target),dim,target)
            projected_relations=Set(
                (link.dim,Int(link.slave_entity),Int(link.master_entity))
                for link in projected.periodic_links)
            # Stored volume relations are mesh-inert and never serialize into
            # $Periodic, matching Gmsh 4.15.2, so they are not projected.
            missing_relations=sort!([
                (constraint.dim,Int(constraint.slave_entity),
                 Int(constraint.master_entity))
                for constraint in constraints
                if constraint.dim!=3 &&
                   !((constraint.dim,Int(constraint.slave_entity),
                      Int(constraint.master_entity)) in projected_relations)])
            isempty(missing_relations) || throw(ArgumentError(
                "tessella: selected $entity_name projection omits periodic " *
                "relations $missing_relations"))
            write_mixed_msh(destination,projected;version,binary)
        end
        return destination
    end
    return path
end

end # module
