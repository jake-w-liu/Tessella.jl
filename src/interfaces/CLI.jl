"""
    CLI

Command-line façade: `tessella file.geo -1|-2|-3` executes the bounded `.geo`
subset and writes `file.msh`. Supported periodic and embedded surface entities
and embedded volume entities, including nested sheet constraints and periodic
planar boundaries of explicit volume shells, are written with classified MSH4
metadata. Unknown flags and OCC-only files are blockers.

The session path (`-refine`, `-reclassify`, `-reparam`, `-part`, `-save`,
`-open`, `-merge`, `-new`, and `.msh`/`.pos` inputs) runs on the `API` session:
`API.open`/`API.merge` import the files, `API.mesh.generate`/`refine`/
`classify_surfaces`/`create_geometry`/`partition` perform the batch action, and
`API.write` exports the current model mesh. `-convert` rewrites each listed
file into its `_new` counterpart (binary MSH 4.1 for meshes, parsed `.pos`
text for views), matching the upstream conversion contract. `-part` assigns
partition labels like `API.mesh.partition`; `-part_topo`/`-part_ghosts` opt
into partition topology/ghost production where the discrete producer supports
it and reject loudly elsewhere. `.geo` `-convert` rejects (no unrolled writer).
`-reclassify`/`-reparam` route to the real classify/create-geometry pipeline —
BRep-imported models can hit bounded export paths that reject rather than write
wrong data.
"""
module CLI

using ..GeoExec: execute_geo, geo_entity_mesh
using ..IO: write_msh
using ..Model: model_periodic_constraints, model_to_mixed
using ..Elements: write_mixed_msh, read_mixed_msh, MixedMesh, simplex_to_mixed, validate
using ..GmshLibm: _with_win_cli_precision
import ..IO
import ..Post
import ..API

export main

const _HELP="""
Usage: tessella file.geo [-1|-2|-3] [-o output.msh] [-format msh2|msh4] [-bin]
       [-clscale value] [-clmin value] [-clmax value] [-order int]
       [-tol value] [-clcurv n] [-epslc1d value] [-edgelmin value]
       [-clextend n] [-smooth n] [-algo del2d|front2d|del3d] [-quad]
       [-optimize] [-optimize_netgen] [-optimize_threshold value]
       [-save_all] [-nopopup] [-v int] [-nt int] [-setnumber name value]
       [-setstring name value] [-setlist name n v1 ... vn] [-string str]
       [-option file.geo]
       tessella file.msh -check
       tessella file.(geo|msh|pos) [-refine] [-part n] [-save]
       [-reclassify deg] [-reparam deg] [-open] [-merge] [-new]
       [-part_topo|-part_no_topo] [-part_ghosts|-part_no_ghosts]
       tessella -convert file.(msh|pos) [file2 ...]
       tessella -help | -version | -info | -help_options

Mesh size bounds apply before the positive size factor. CLI numeric options
have priority over .geo assignments. Curve orders 1-10 are supported;
Surface/Volume generation supports the orders documented in README.md.
-check validates MSH structure and duplicate cells without generation or output.
-setnumber/-setstring/-setlist and -string/-option act as .geo parser input
before the input file; on non-.geo inputs only `-setnumber` option names
(Family.Member) apply, through the session option table.
GUI-only, display, solver-daemon and unimplemented upstream flags reject.
"""

# Flags that set one numeric mesh/session option. The value is stored through
# the `.geo` option setter on the standalone path (coercion and aliases
# included) and through `API.option` on the session path.
const _NUMERIC_FLAGS=Dict(
    "-clscale"=>("Mesh.MeshSizeFactor",:positive),
    "-clmin"=>("Mesh.MeshSizeMin",:nonnegative),
    "-clmax"=>("Mesh.MeshSizeMax",:positive),
    "-order"=>("Mesh.ElementOrder",:order),
    "-degree"=>("Mesh.ElementOrder",:order),
    "-tol"=>("Geometry.Tolerance",:positive),
    "-clcurv"=>("Mesh.MeshSizeFromCurvature",:nonnegative),
    "-epslc1d"=>("Mesh.LcIntegrationPrecision",:positive),
    "-edgelmin"=>("Mesh.ToleranceEdgeLength",:positive),
    "-clextend"=>("Mesh.MeshSizeExtendFromBoundary",:finite),
    "-smooth"=>("Mesh.Smoothing",:nonnegint),
    "-optimize_threshold"=>("Mesh.OptimizeThreshold",:finite),
    "-v"=>("General.Verbosity",:finite),
    "-verbose"=>("General.Verbosity",:finite),
    "-nt"=>("General.NumThreads",:nonnegint))

# Valueless flags that store `1` into the given numeric option.
const _BOOLEAN_FLAGS=Dict(
    "-optimize"=>"Mesh.Optimize",
    "-optimize_netgen"=>"Mesh.OptimizeNetgen",
    "-save_all"=>"Mesh.SaveAll",
    "-saveall"=>"Mesh.SaveAll",
    "-nopopup"=>"General.NoPopup",
    "-renumber"=>"Mesh.Renumber",
    "-term"=>"General.Terminal")

# `-algo` names we can back with a real implemented algorithm.
const _ALGORITHMS=Dict(
    "del2d"=>("Mesh.Algorithm",5.0),"delaunay"=>("Mesh.Algorithm",5.0),
    "tri"=>("Mesh.Algorithm",5.0),
    "front2d"=>("Mesh.Algorithm",6.0),"frontal"=>("Mesh.Algorithm",6.0),
    "del3d"=>("Mesh.Algorithm3D",1.0),"delaunay3d"=>("Mesh.Algorithm3D",1.0),
    "gmsh3d"=>("Mesh.Algorithm3D",1.0))

# Upstream flags that are explicitly not implemented in Tessella. Each gets a
# loud rejection naming the missing functionality instead of falling into the
# generic "unsupported flag" error (which stays for genuinely unknown words).
const _REJECTED_FLAGS=Dict{String,String}(
    "-0"=>"model .geo output (unrolled geometry) is not implemented",
    "-4"=>"adaptive meshing (batch mode 4) is not implemented",
    "-barycentric_refine"=>"barycentric refinement is not implemented",
    "-parametric"=>"parametric node output control (Mesh.SaveParametric) is not implemented",
    "-save_parametric"=>"parametric node output control (Mesh.SaveParametric) is not implemented",
    "-save_topology"=>"topology-only output (Mesh.SaveTopology) is not implemented",
    "-save_entities"=>"topology-only output (Mesh.SaveTopology) is not implemented",
    "-switch_tags"=>"element tag switching is not implemented",
    "-preserve_numbering_msh2"=>"MSH2 element-number preservation is not implemented",
    "-preserveNumberingMsh2"=>"MSH2 element-number preservation is not implemented",
    "-part_weight"=>"weighted graph partitioning is not implemented",
    "-partWeight"=>"weighted graph partitioning is not implemented",
    "-part_split"=>"per-partition file output is not implemented",
    "-part_split_files"=>"per-partition file output is not implemented",
    "-oneFilePerPart"=>"per-partition file output is not implemented",
    "-part_topo_pro"=>"partition topology .pro output is not implemented",
    "-savePartTopology"=>"partition topology .pro output is not implemented",
    "-part_physicals"=>"partition physical-group toggling is not implemented",
    "-part_no_physicals"=>"partition physical-group toggling is not implemented",
    "-optimize_ho"=>"high-order mesh optimization is not implemented",
    "-hoOptimize"=>"high-order mesh optimization is not implemented",
    "-ho_min"=>"high-order optimization is not implemented",
    "-ho_max"=>"high-order optimization is not implemented",
    "-ho_nlayers"=>"high-order optimization is not implemented",
    "-aniso_max"=>"the BAMG backend is not implemented",
    "-anisoMax"=>"the BAMG backend is not implemented",
    "-smooth_ratio"=>"the BAMG backend is not implemented",
    "-smoothRatio"=>"the BAMG backend is not implemented",
    "-swapangle"=>"boundary-edge swap optimization is not implemented",
    "-rand"=>"random mesh perturbation is not implemented",
    "-scale"=>"geometrical scaling is not implemented",
    "-meshscale"=>"mesh scaling is not implemented",
    "-clcurviso"=>"isotropic curvature sizing is not implemented",
    "-claniso"=>"anisotropic curvature sizing is not implemented",
    "-ignore_periodicity"=>"periodic-constraint ignoring is not implemented",
    "-ignorePeriodicity"=>"periodic-constraint ignoring is not implemented",
    "-ignore_periocity"=>"periodic-constraint ignoring is not implemented",
    "-bgm"=>"background mesh file import is not implemented",
    "-match"=>"geometry/mesh matching is not implemented",
    "-numsubedges"=>"high-order element display subdivision is GUI-only",
    "-debugSurface"=>"surface debugging is GUI-only",
    "-statreport"=>"mesh statistics reports are not implemented",
    "-append_statreport"=>"mesh statistics reports are not implemented",
    "-run"=>"ONELAB solver clients are not implemented",
    "-onelab"=>"ONELAB solver clients are not implemented",
    "-socket"=>"ONELAB solver clients are not implemented",
    "-listen"=>"solver socket listeners are not implemented",
    "-minterpreter"=>"Octave interpreters are not implemented",
    "-pyinterpreter"=>"Python interpreters are not implemented",
    "-cpu"=>"CPU-time reporting is not implemented",
    "-log"=>"message log files are not implemented",
    "-watch"=>"file watching is not implemented",
    "-bg"=>"background images are GUI-only",
    "-n"=>"view visibility options are GUI-only",
    "-nodb"=>"double buffering is GUI-only",
    "-fontsize"=>"GUI font options are GUI-only",
    "-deltafontsize"=>"GUI font options are GUI-only",
    "-theme"=>"GUI themes are GUI-only",
    "-scheme"=>"GUI themes are GUI-only",
    "-display"=>"the display server is GUI-only",
    "-camera"=>"camera mode is GUI-only",
    "-stereo"=>"stereo rendering is GUI-only",
    "-gamepad"=>"gamepad input is GUI-only",
    "-a"=>"initial GUI modules are GUI-only",
    "-g"=>"initial GUI modules are GUI-only",
    "-m"=>"initial GUI modules are GUI-only",
    "-s"=>"initial GUI modules are GUI-only",
    "-p"=>"initial GUI modules are GUI-only",
    "-windows"=>"GUI window management is GUI-only",
    "-tiles"=>"GUI window tiling is GUI-only",
    "-dual"=>"dual-mesh display is GUI-only",
    "-voronoi"=>"Voronoi display is GUI-only",
    "-noview"=>"view visibility is GUI-only",
    "-nomesh"=>"mesh visibility is GUI-only",
    "-link"=>"view linking is GUI-only",
    "-smoothview"=>"view smoothing is GUI-only",
    "-combine"=>"view combination is GUI-only")

# Flags that upstream removed before 4.15.2 or that never existed there —
# rejected with an explicit "not an upstream option" message.
const _OBSOLETE_FLAGS=("-newreg","-newvol","-newsur","-post_pro","-clear",
    "-liststring","-listnumber","-min","-min0","-min1","-min2","-max",
    "-view","-size_field","-quadlayout","-autoex","-openFile")

function _numeric_argument(flag::String,raw::String,kind::Symbol)
    value=tryparse(Float64,raw)
    value!==nothing && isfinite(value) || throw(ArgumentError(
        "tessella: $flag needs a finite number"))
    if kind==:order
        isinteger(value) && 1<=value<=10 || throw(ArgumentError(
            "tessella: $flag needs an integer from 1 to 10"))
    elseif kind==:positive
        value>0 || throw(ArgumentError("tessella: $flag must be positive"))
    elseif kind==:nonnegative
        value>=0 || throw(ArgumentError("tessella: $flag must be nonnegative"))
    elseif kind==:nonnegint
        isinteger(value) && value>=0 || throw(ArgumentError(
            "tessella: $flag needs a nonnegative integer"))
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

function _help_options()
    io=IOBuffer()
    for (family,table) in sort!(collect(pairs(IO._GEO_NUMBER_OPTION_TABLE));by=first)
        println(io,"$family number options:")
        for (member,default) in sort!(collect(pairs(table));by=first)
            println(io,"  $family.$member = $default")
        end
    end
    for (family,table) in sort!(collect(pairs(IO._GEO_STRING_OPTION_TABLE));by=first)
        println(io,"$family string options:")
        for (member,default) in sort!(collect(pairs(table));by=first)
            println(io,"  $family.$member = \"$default\"")
        end
    end
    return String(take!(io))
end

function _info()
    return "tessella $(Base.pkgversion(parentmodule(@__MODULE__))) " *
        "(Tessella.jl, Julia $VERSION, $(Sys.MACHINE))"
end

# The injection flags produce literal `.geo` statements that run before the
# input file — upstream evaluates `-string`, `-option`, `-setnumber`,
# `-setstring` and `-setlist` on the parser before opening files, so prepending
# the equivalent statements reproduces the same ordering (script assignments
# in the file itself still take precedence, like upstream's option writes).
function _escape_geo_string(value::String,flag::String)
    any(c->c in ('"','\n','\r','\0'),value) && throw(ArgumentError(
        "tessella: $flag value cannot contain quotes or line breaks"))
    return "\""*value*"\""
end

function _option_name(token::String)
    # `-setnumber`/`-setstring` distinguish `Family.Member` option names from
    # parser variables only on non-.geo inputs (on .geo inputs both are just
    # prepended statements and the executor validates them).
    occursin('.',token) && return true
    return false
end

function _inject_geo(path::String,preludes)
    # Materialize the input with the injected statements prepended, in the
    # input's own directory so relative `Merge`/`Include` and file outputs
    # keep resolving against the same location. `(name, value)` option
    # preludes become ordinary option assignments — the `.geo` option setter
    # applies its usual coercion and alias mirroring.
    directory=dirname(abspath(path))
    temporary=tempname(directory)*".geo"
    open(temporary,"w") do io
        for prelude in preludes
            prelude isa String && (println(io,prelude);continue)
            println(io,"$(prelude[1]) = $(prelude[2]);")
        end
        open(path) do input
            write(io,input)
        end
    end
    return temporary
end

function _convert_file(path::String)
    isfile(path) || throw(ArgumentError(
        "tessella: -convert input file $path does not exist"))
    extension=lowercase(splitext(path)[2])
    destination=path*"_new"
    _same_file_target(path,destination) && throw(ArgumentError(
        "tessella: -convert output path must not alias the input file"))
    if extension==".msh"
        mesh=read_mixed_msh(path)
        write_mixed_msh(destination,mesh;version=4.1,binary=true)
        return destination
    elseif extension==".pos"
        views=Post.read_pos(path)
        isempty(views) && throw(ArgumentError(
            "tessella: -convert found no views in $path"))
        Post.write_pos(destination,views)
        return destination
    elseif extension==".geo"
        # Upstream writes an unrolled `.geo` for `.geo` inputs — Tessella has
        # no unrolled-geometry writer, so this is an explicit rejection.
        throw(ArgumentError(
            "tessella: -convert of .geo inputs needs the unrolled .geo " *
            "exporter, which is not implemented"))
    end
    throw(ArgumentError(
        "tessella: -convert supports .msh and .pos inputs, got $path"))
end

"""
    main(args) -> String

Execute the bounded Gmsh-style command line. A single `.geo` input runs the
`.geo` executor with optional `-1`, `-2` or `-3` meshing and `-o output`.
`-format`/`-f`/`msh2|msh4` or `-vmsh` select MSH 2.2 or 4.1 (the default), and
`-bin` selects binary output. Output flags require an output-producing action.
`-clscale`, `-clmin`, `-clmax`, `-order`/`-degree`, `-tol`, `-clcurv`,
`-epslc1d`, `-edgelmin`, `-clextend`, `-smooth`, `-optimize_threshold`,
`-v`/`-verbose` and `-nt` set real numeric options; they retain priority over
script assignments, including option aliases. `-optimize`, `-optimize_netgen`,
`-save_all`/`-saveall`, `-renumber`, `-nopopup`, `-term` and `-quad` set the
matching boolean/algorithm options; `-algo name` selects an implemented
algorithm (`del2d`, `front2d`, `del3d` and their upstream aliases).
`-setnumber`, `-setstring`, `-setlist`/`-setlistofnumbers`, `-string` and
`-option` inject parser input before a `.geo` file (on non-`.geo` inputs only
`-setnumber` option names apply, through `API.option`).
`-refine`, `-reclassify deg`, `-reparam deg`, `-save`, `-part n`,
`-part_topo`/`-part_no_topo`, `-part_ghosts`/`-part_no_ghosts`, `-open`,
`-merge`, `-new`, and `.msh`/`.pos` inputs run on the `API` session; `-part n`
partitions the current mesh after meshing (or the loaded mesh) exactly like
`API.mesh.partition` — `-part_topo`/`-part_ghosts` request partition
topology/ghost production and reject on unsupported inputs.
`-convert files...` rewrites `.msh`/`.pos` files to `file_new` (binary MSH 4.1
or parsed view text) and must come last.
`file.msh -check` validates the mesh without writing. `-`, `-parse_and_exit`,
`-help`/`--help`, `-version`/`--version`, `-info`/`--info`, `-help_options`
and `-pid` are informational; the first four must be used alone.
Duplicate or conflicting flags, unmarked multiple inputs, ignored output
arguments, and any output that aliases its input are rejected. Upstream
GUI-only, solver-daemon and unimplemented flags are rejected by name.
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
    information=filter(a->a in ("-help","--help","-version","--version",
        "-info","--info","-help_options"),tokens)
    if !isempty(information)
        length(tokens)==1 || throw(ArgumentError(
            "tessella: informational flags must be used alone"))
        flag=only(tokens)
        message=flag in ("-help","--help") ? _HELP :
            flag in ("-info","--info") ? _info() :
            flag=="-help_options" ? _help_options() :
            string(Base.pkgversion(parentmodule(@__MODULE__)))
        println(message)
        return message
    end
    numeric_options=Dict{String,Float64}()
    # Pre-file parser input, in flag order: literal `.geo` statements and
    # `(name, value)` numeric option assignments. On `.geo` inputs both are
    # prepended to the executed source; on non-.geo inputs only the option
    # assignments apply, through `API.option`.
    preludes=Union{String,Tuple{String,Float64}}[]
    # File-list markers: `-open`/`-merge` are persistent modes for following
    # files; `-new` creates a fresh model before the next file, which then
    # merges (a fresh model merge is an open). Bare positional files after
    # the first still require an explicit marker.
    files=Tuple{String,Symbol}[]
    marker=:none;markers_seen=false
    dim=0;outfile=nothing;version=4.1;binary=false
    dimension_seen=false;output_seen=false
    format_seen=false;binary_seen=false;check_seen=false
    refine_seen=false;save_seen=false;parse_seen=false
    reclassify=nothing;reparam=nothing
    part=nothing;part_topo=nothing;part_ghosts=nothing
    algo_seen=false;quad_seen=false
    i=1
    while i<=length(tokens)
        a=tokens[i]
        if a=="-convert"
            isempty(files) || throw(ArgumentError(
                "tessella: -convert cannot follow an input file"))
            i+1<=length(tokens) || throw(ArgumentError(
                "tessella: -convert needs at least one file"))
            outputs=[_convert_file(_path_argument(token,"-convert"))
                     for token in tokens[i+1:end]]
            return last(outputs)
        elseif a in ("-","-parse_and_exit")
            parse_seen && throw(ArgumentError(
                "tessella: duplicate $a flag"))
            parse_seen=true
        elseif a in ("-1","-2","-3")
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
        elseif a in ("-format","-f")
            format_seen && throw(ArgumentError(
                "tessella: duplicate output format flag"))
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: $a needs msh2 or msh4"))
            format=tokens[i]
            format in ("msh2","msh4") || throw(ArgumentError(
                "tessella: unsupported output format $format; supported: msh2, msh4"))
            version=format=="msh2" ? 2.2 : 4.1
            format_seen=true
        elseif a=="-vmsh"
            format_seen && throw(ArgumentError(
                "tessella: duplicate output format flag"))
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: -vmsh needs a file version"))
            requested=tryparse(Float64,tokens[i])
            requested!==nothing && requested in (2.2,4.1) || throw(ArgumentError(
                "tessella: -vmsh supports MSH versions 2.2 and 4.1"))
            version=requested
            format_seen=true
        elseif haskey(_NUMERIC_FLAGS,a)
            option,kind=_NUMERIC_FLAGS[a]
            haskey(numeric_options,option) && throw(ArgumentError(
                "tessella: duplicate or conflicting option flag $a"))
            i+=1;i<=length(tokens) || throw(ArgumentError("tessella: $a needs a value"))
            numeric_options[option]=_numeric_argument(a,tokens[i],kind)
        elseif haskey(_BOOLEAN_FLAGS,a)
            option=_BOOLEAN_FLAGS[a]
            haskey(numeric_options,option) && throw(ArgumentError(
                "tessella: duplicate or conflicting option flag $a"))
            numeric_options[option]=1.0
        elseif a=="-quad"
            quad_seen && throw(ArgumentError("tessella: duplicate -quad flag"))
            numeric_options["Mesh.RecombineAll"]=1.0
            numeric_options["Mesh.RecombinationAlgorithm"]=2.0
            quad_seen=true
        elseif a=="-algo"
            algo_seen && throw(ArgumentError("tessella: duplicate -algo flag"))
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: -algo needs an algorithm name"))
            name=tokens[i]
            haskey(_ALGORITHMS,name) || throw(ArgumentError(
                "tessella: unsupported mesh algorithm \"$name\"; supported: " *
                "del2d, front2d, del3d"))
            option,value=_ALGORITHMS[name]
            haskey(numeric_options,option) && throw(ArgumentError(
                "tessella: -algo conflicts with an earlier $option option"))
            numeric_options[option]=value
            algo_seen=true
        elseif a=="-check"
            check_seen && throw(ArgumentError("tessella: duplicate -check flag"))
            check_seen=true
        elseif a=="-bin"
            binary_seen && throw(ArgumentError("tessella: duplicate -bin flag"))
            binary=true;binary_seen=true
        elseif a=="-pid"
            println(getpid())
        elseif a=="-refine"
            refine_seen && throw(ArgumentError("tessella: duplicate -refine flag"))
            refine_seen=true
        elseif a=="-save"
            save_seen && throw(ArgumentError("tessella: duplicate -save flag"))
            save_seen=true
        elseif a=="-reclassify" || a=="-reparam"
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: $a needs an angle in degrees"))
            angle=tryparse(Float64,tokens[i])
            angle!==nothing && isfinite(angle) || throw(ArgumentError(
                "tessella: $a needs a finite angle in degrees"))
            a=="-reclassify" && reclassify!==nothing && throw(ArgumentError(
                "tessella: duplicate -reclassify flag"))
            a=="-reparam" && reparam!==nothing && throw(ArgumentError(
                "tessella: duplicate -reparam flag"))
            a=="-reclassify" ? (reclassify=angle) : (reparam=angle)
        elseif a=="-part"
            part!==nothing && throw(ArgumentError("tessella: duplicate -part flag"))
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: -part needs a partition count"))
            count=tryparse(Float64,tokens[i])
            count!==nothing && isinteger(count) && count>=1 || throw(ArgumentError(
                "tessella: -part needs a positive integer partition count"))
            part=Int(count)
        elseif a=="-part_topo"
            part_topo!==nothing && throw(ArgumentError(
                "tessella: duplicate or conflicting partition topology flag"))
            part_topo=true
        elseif a=="-part_no_topo"
            part_topo!==nothing && throw(ArgumentError(
                "tessella: duplicate or conflicting partition topology flag"))
            part_topo=false
        elseif a=="-part_ghosts"
            part_ghosts!==nothing && throw(ArgumentError(
                "tessella: duplicate or conflicting partition ghost flag"))
            part_ghosts=true
        elseif a=="-part_no_ghosts"
            part_ghosts!==nothing && throw(ArgumentError(
                "tessella: duplicate or conflicting partition ghost flag"))
            part_ghosts=false
        elseif a in ("-new","-open","-merge")
            marker=Symbol(a[2:end]);markers_seen=true
        elseif a=="-setnumber"
            i+=2;i<=length(tokens) || throw(ArgumentError(
                "tessella: -setnumber needs a name and a value"))
            name=tokens[i-1]
            (isempty(name) || startswith(name,'-')) && throw(ArgumentError(
                "tessella: -setnumber needs a name"))
            value=tryparse(Float64,tokens[i])
            value!==nothing && isfinite(value) || throw(ArgumentError(
                "tessella: -setnumber needs a finite value"))
            _option_name(name) ? push!(preludes,(name,value)) :
                push!(preludes,"$name = $(tokens[i]);")
        elseif a=="-setstring"
            i+=2;i<=length(tokens) || throw(ArgumentError(
                "tessella: -setstring needs a name and a value"))
            name=tokens[i-1]
            (isempty(name) || startswith(name,'-')) && throw(ArgumentError(
                "tessella: -setstring needs a name"))
            startswith(tokens[i],'-') && throw(ArgumentError(
                "tessella: -setstring needs a string value"))
            # Option names and parser variables are both statements for the
            # `.geo` parser; on non-.geo inputs the session path rejects them.
            push!(preludes,"$name = $(_escape_geo_string(tokens[i],a));")
        elseif a in ("-setlist","-setlistofnumbers")
            i+=2;i<=length(tokens) || throw(ArgumentError(
                "tessella: $a needs a name, a count and values"))
            name=tokens[i-1]
            (isempty(name) || startswith(name,'-')) && throw(ArgumentError(
                "tessella: $a needs a name"))
            count=tryparse(Float64,tokens[i])
            count!==nothing && isinteger(count) && count>=0 || throw(
                ArgumentError("tessella: $a needs a nonnegative count"))
            needed=Int(count)
            i+needed<=length(tokens) || throw(ArgumentError(
                "tessella: $a expected $needed values"))
            values=tokens[i+1:i+needed];i+=needed
            all(token->tryparse(Float64,token)!==nothing,values) || throw(
                ArgumentError("tessella: $a values must be finite numbers"))
            push!(preludes,"$name = {$(join(values,","))};")
        elseif a=="-string"
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: -string needs a statement"))
            push!(preludes,tokens[i])
        elseif a=="-option"
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: -option needs an option file"))
            option_file=_path_argument(tokens[i],"option")
            isfile(option_file) || throw(ArgumentError(
                "tessella: -option file $option_file does not exist"))
            push!(preludes,
                "Include $(_escape_geo_string(abspath(option_file),a));")
        else
            if startswith(a,"-")
                haskey(_REJECTED_FLAGS,a) && throw(ArgumentError(
                    "tessella: $a is not supported: $(_REJECTED_FLAGS[a])"))
                a in _OBSOLETE_FLAGS && throw(ArgumentError(
                    "tessella: $a is not an upstream Gmsh 4.15.2 option"))
                throw(ArgumentError("tessella: unsupported flag $a"))
            end
            if !isempty(files) && !markers_seen
                throw(ArgumentError("tessella: multiple input files"))
            end
            file=_path_argument(a,"input")
            push!(files,(file,marker))
            # `-new` applies to its next file only; afterwards files merge.
            marker = marker==:new ? :merge : marker
        end
        i+=1
    end
    isempty(files) && throw(ArgumentError("tessella: missing input file"))
    # Upstream batch actions are mutually exclusive — the strict CLI rejects
    # combinations instead of keeping only the last one.
    batch_actions=(dimension_seen ? 1 : 0)+(refine_seen ? 1 : 0)+
        (save_seen ? 1 : 0)+(reclassify===nothing ? 0 : 1)+
        (reparam===nothing ? 0 : 1)
    batch_actions>1 && throw(ArgumentError(
        "tessella: conflicting batch action flags"))
    inputs=[file for (file,_) in files]
    extension=lowercase(splitext(first(inputs))[2])
    if check_seen
        (dimension_seen || output_seen || format_seen || binary_seen ||
         !isempty(numeric_options) || !isempty(preludes) ||
         refine_seen || save_seen || parse_seen ||
         reclassify!==nothing || reparam!==nothing || part!==nothing ||
         length(files)>1) && throw(ArgumentError(
            "tessella: -check cannot be combined with generation, sizing, " *
            "session or output options"))
        length(files)==1 || throw(ArgumentError(
            "tessella: -check accepts a single input"))
        extension==".msh" || throw(ArgumentError(
            "tessella: -check requires an MSH input"))
        mesh=read_mixed_msh(first(inputs))
        diagnostic=validate(mesh)
        diagnostic.ok || throw(ArgumentError(
            "tessella: inconsistent mesh: "*join(diagnostic.messages,"; ")))
        return first(inputs)
    end
    parse_seen && (dimension_seen || save_seen || refine_seen ||
        reclassify!==nothing || reparam!==nothing || part!==nothing ||
        output_seen || format_seen || binary_seen) && throw(ArgumentError(
            "tessella: parse-and-exit cannot be combined with meshing or " *
            "output options"))
    # `.msh`/`.pos` inputs, multiple files, and post-mesh actions run on the
    # API session. Everything else — including extensionless `.geo` sources —
    # keeps the standalone `.geo` execution path.
    session=extension in (".msh",".pos") || length(files)>1 || refine_seen ||
        save_seen || reclassify!==nothing || reparam!==nothing ||
        part!==nothing
    if !session
        dim==0 && outfile!==nothing && throw(ArgumentError(
            "tessella: -o requires -1, -2 or -3"))
        dim==0 && (format_seen || binary_seen) && throw(ArgumentError(
            "tessella: output format flags require -1, -2 or -3"))
        destination=dim==0 ? nothing :
            (outfile===nothing ? _default_output(first(inputs)) : outfile)
        destination!==nothing &&
            _same_file_target(first(inputs),destination) && throw(ArgumentError(
                "tessella: output path must not alias the input file"))
        get(numeric_options,"Mesh.MeshSizeMin",0.0)<=
            get(numeric_options,"Mesh.MeshSizeMax",1e22) || throw(ArgumentError(
                "tessella: minimum mesh size exceeds maximum"))
        source=first(inputs)
        temporary=nothing
        result=try
            if !isempty(preludes)
                temporary=_inject_geo(source,preludes)
                source=temporary
            end
            execute_geo(source;mesh_dim=dim,_forced_numeric_options=numeric_options)
        finally
            temporary===nothing || rm(temporary;force=true)
        end
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
        return first(inputs)
    end
    return _session_main(files,preludes,numeric_options;
        dim,outfile,version,binary,output_flags=format_seen || binary_seen,
        refine=refine_seen,save=save_seen,
        reclassify,reparam,part,part_topo,part_ghosts)
end

# Session path: files route through `API.open`/`API.merge`/`API.model.add`,
# option writes through `API.option`, meshing through `API.mesh.generate`,
# post-mesh actions through `API.mesh.refine`/`classify_surfaces`/
# `create_geometry`/`partition`, and output through `API.write`.
function _session_main(files,preludes,numeric_options;
                       dim,outfile,version,binary,output_flags,refine,save,
                       reclassify,reparam,part,part_topo,part_ghosts)
    need_output=dim>0 || refine || save || reclassify!==nothing ||
        reparam!==nothing || part!==nothing
    !need_output && outfile!==nothing && throw(ArgumentError(
        "tessella: -o requires a meshing, conversion or save action"))
    !need_output && output_flags && throw(ArgumentError(
        "tessella: output format flags require a meshing, conversion or save action"))
    # Literal `.geo` statements need a leading `.geo` file; `(name, value)`
    # option preludes become session options when the session table knows the
    # name and `.geo` assignments otherwise (`General.NoPopup` and friends).
    first_is_geo=lowercase(splitext(files[1][1])[2])==".geo"
    option_preludes=Tuple{String,Float64}[]
    lines=String[]
    for prelude in preludes
        if prelude isa Tuple
            haskey(API.DEFAULT_OPTIONS,prelude[1]) ?
                push!(option_preludes,prelude) :
                push!(lines,"$(prelude[1]) = $(prelude[2]);")
        else
            push!(lines,prelude)
        end
    end
    !first_is_geo && !isempty(lines) && throw(ArgumentError(
        "tessella: -setstring/-setlist/-string/-option inject .geo parser " *
        "input and require a .geo input file"))
    # Session-table options apply through `API.option`; options that live only
    # in the `.geo` table (`General.NoPopup`, `General.Terminal`) are injected
    # into a leading `.geo` file and are inert on non-.geo inputs.
    session_numeric=Dict{String,Float64}()
    geo_numeric_preludes=Union{String,Tuple{String,Float64}}[]
    for (name,value) in numeric_options
        haskey(API.DEFAULT_OPTIONS,name) ? (session_numeric[name]=value) :
            push!(geo_numeric_preludes,"$name = $value;")
    end
    destination=outfile===nothing ? _default_output(files[1][1]) : outfile
    extension=lowercase(splitext(files[1][1])[2])
    if extension!=".geo" && need_output
        outfile===nothing && throw(ArgumentError(
            "tessella: non-.geo inputs require -o for output"))
    end
    if need_output
        for (file,) in files
            _same_file_target(file,destination) && throw(ArgumentError(
                "tessella: output path must not alias the input file"))
        end
    end
    get(numeric_options,"Mesh.MeshSizeMin",0.0)<=
        get(numeric_options,"Mesh.MeshSizeMax",1e22) || throw(ArgumentError(
            "tessella: minimum mesh size exceeds maximum"))
    API.initialize()
    temporary=nothing
    try
        for (name,value) in option_preludes
            API.option(name,value)
        end
        for (index,(file,file_marker)) in enumerate(files)
            source=file
            injected=vcat(lines,geo_numeric_preludes)
            if index==1 && first_is_geo && !isempty(injected)
                temporary=_inject_geo(file,injected)
                source=temporary
            end
            # `-new` creates a fresh current model before its file merges —
            # merging into an empty model is an open. `-open` opens; the
            # default and `-merge` merge.
            file_marker==:new && API.model.add(index<length(files) ?
                splitext(basename(files[index+1][1]))[1] : "newmodel")
            if index==1 || file_marker==:open
                API.open(source)
            else
                API.merge(source)
            end
        end
        for (name,value) in session_numeric
            API.option(name,value)
        end
        API.option("Mesh.MshFileVersion",version)
        API.option("Mesh.Binary",binary ? 1.0 : 0.0)
        part_topo===nothing ||
            API.option("Mesh.PartitionCreateTopology",part_topo ? 1.0 : 0.0)
        part_ghosts===nothing ||
            API.option("Mesh.PartitionCreateGhostCells",part_ghosts ? 1.0 : 0.0)
        dim>0 && API.mesh.generate(dim)
        refine && API.mesh.refine()
        reclassify===nothing ||
            API.mesh.classify_surfaces(reclassify*pi/180,true,false,pi)
        if reparam!==nothing
            API.mesh.classify_surfaces(reparam*pi/180,true,true,pi)
            API.mesh.create_geometry()
        end
        part===nothing || API.mesh.partition(part)
        need_output && API.write(destination)
    finally
        temporary===nothing || rm(temporary;force=true)
        API.finalize()
    end
    return need_output ? destination : files[1][1]
end

end # module
