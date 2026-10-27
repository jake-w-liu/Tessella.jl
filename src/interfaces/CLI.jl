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
it and reject loudly elsewhere. `-part_split` (aliases `-part_split_files`,
`-oneFilePerPart`) writes one `<stem>_<p>.msh` file per partition instead of
the combined output when `-part` is active, carrying each partition's own
entities' cells plus its ghost cells as regular cells on its ghost entity —
the upstream `Mesh.PartitionSplitMeshFiles` layout. `.geo` `-convert` rejects
(no unrolled writer).
`-reclassify`/`-reparam` route to the real classify/create-geometry pipeline —
BRep-imported models can hit bounded export paths that reject rather than write
wrong data.
"""
module CLI

using ..GeoExec: execute_geo, geo_entity_mesh
using ..IO: write_msh, _msh_rewrite_save_options
using ..Model: model_periodic_constraints, model_to_mixed, model_entities
using ..Elements: write_mixed_msh, read_mixed_msh, MixedMesh, MixedEntity,
                  _scaled_mesh_copy, write_partitioned_topology,
                  MixedEntityData, MixedPartitionData, MixedGhostElement, MixedPeriodicLink,
                  ElementBlock, SpecialElementBlock, MixedElementBlock, ElementRef,
                  msh_dimension, msh_spec, simplex_to_mixed, barycentric_refine, validate
using ..MeshTypes: Mesh
using ..GeoUnrolled: write_geo_unrolled
using ..GmshLibm: _with_win_cli_precision
using Printf: @sprintf
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
       [-optimize_ho|-hoOptimize] [-ho_min v] [-ho_max v] [-ho_nlayers v]
       [-aniso_max value] [-rand value] [-smooth_ratio value] [-clcurviso]
       [-ignore_periodicity] [-bgm file.pos]
       [-save_parametric|-parametric] [-save_topology|-save_entities]
       [-save_all] [-nopopup] [-v int] [-nt int] [-setnumber name value]
       [-setstring name value] [-setlist name n v1 ... vn] [-string str]
       [-option file.geo] [-log file] [-cpu] [-statreport file]
       [-append_statreport file] [-switch_tags] [-preserve_numbering_msh2]
       tessella file.msh -check
       tessella file.(geo|msh|pos) [-refine] [-barycentric_refine] [-part n] [-save] [-0]
       [-reclassify deg] [-reparam deg] [-open] [-merge] [-new]
       [-part_topo|-part_no_topo] [-part_ghosts|-part_no_ghosts]
       [-part_split|-part_split_files|-oneFilePerPart]
       tessella -convert file.(msh|pos) [file2 ...]
       tessella -help | -version | -info | -help_options

Mesh size bounds apply before the positive size factor. CLI numeric options
have priority over .geo assignments. Curve orders 1-10 are supported;
Surface/Volume generation supports the orders documented in README.md.
-check validates MSH structure and duplicate cells without generation or output.
-part_split (aliases -part_split_files and -oneFilePerPart) replaces combined
output with one <stem>_<partition>.msh file per partition when -part is active;
without -part it is inert, like upstream Mesh.PartitionSplitMeshFiles.
-setnumber/-setstring/-setlist and -string/-option act as .geo parser input
before the input file; on non-.geo inputs only `-setnumber` option names
(Family.Member) apply, through the session option table.
-save_parametric/-parametric and -save_topology/-save_entities honor
Mesh.SaveParametric and Mesh.SaveTopology on MSH 2.2 and MSH 4.1 output (parametric node records and
topology-only \$Entities, matching Gmsh 4.15.2). -optimize_ho/-hoOptimize and
the -ho_min/-ho_max/-ho_nlayers options drive the session high-order
relocation optimizer. -barycentric_refine subdivides the current mesh's
model-dimension simplices around their barycenter (triangles on 2-D models,
tetrahedra on 3-D models), matching Gmsh 4.15.2. -0 writes the executed
model's unrolled .geo script to <stem>.geo_unrolled (or -o), matching Gmsh
4.15.2 — periodic constraints and non-built-in entity kinds are dropped.
-aniso_max stores Mesh.AnisoMax (caps the 1-D background-mesh metric);
-rand/-smooth_ratio store Mesh.RandomFactor/Mesh.SmoothRatio; -clcurviso
stores Mesh.MeshSizeFromCurvatureIsotropic; -ignore_periodicity stores
Mesh.IgnorePeriodicity. -bgm merges a .pos view file and uses its last view
as the PostView background-mesh field.
-log file mirrors Msg-channel output into the file; -cpu appends the
(From start: Wall, CPU, Mem) resource suffix to Info/Error messages.
-statreport writes (truncating) and -append_statreport appends the upstream
2-D statistics header/row after 2-D generation. -switch_tags swaps MSH 2.2
elementary/physical tags on input; -preserve_numbering_msh2 emits the stored
element numbering on MSH 2.2 output.
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
    # `Mesh.AnisoMax` is consumed by `metric_aniso_cap` on the 1-D BGM metric.
    "-aniso_max"=>("Mesh.AnisoMax",:finite),
    "-anisoMax"=>("Mesh.AnisoMax",:finite),
    # `Mesh.RandomFactor`/`Mesh.SmoothRatio` are option stores upstream: the
    # only documented consumer is the absent BAMG backend, so the flag's whole
    # observable contract is the stored value.
    "-rand"=>("Mesh.RandomFactor",:finite),
    "-smooth_ratio"=>("Mesh.SmoothRatio",:finite),
    "-smoothRatio"=>("Mesh.SmoothRatio",:finite),
    "-epslc1d"=>("Mesh.LcIntegrationPrecision",:positive),
    "-edgelmin"=>("Mesh.ToleranceEdgeLength",:positive),
    "-clextend"=>("Mesh.MeshSizeExtendFromBoundary",:finite),
    "-smooth"=>("Mesh.Smoothing",:nonnegint),
    "-optimize_threshold"=>("Mesh.OptimizeThreshold",:finite),
    "-ho_min"=>("Mesh.HighOrderThresholdMin",:finite),
    "-ho_max"=>("Mesh.HighOrderThresholdMax",:finite),
    "-ho_nlayers"=>("Mesh.HighOrderNumLayers",:nonnegint),
    "-v"=>("General.Verbosity",:finite),
    "-verbose"=>("General.Verbosity",:finite),
    "-nt"=>("General.NumThreads",:nonnegint),
    # `-scale`/-meshscale store `Mesh.ScalingFactor` — applied to node
    # coordinates at export by every mesh writer upstream
    # (`CreateFile.cpp` passes it to writeMSH/writeSTL/writeVTK/…).
    "-scale"=>("Mesh.ScalingFactor",:finite),
    "-meshscale"=>("Mesh.ScalingFactor",:finite))

# Valueless flags that store `1` into the given numeric option.
const _BOOLEAN_FLAGS=Dict(
    "-optimize"=>"Mesh.Optimize",
    "-optimize_netgen"=>"Mesh.OptimizeNetgen",
    "-save_all"=>"Mesh.SaveAll",
    "-saveall"=>"Mesh.SaveAll",
    "-save_parametric"=>"Mesh.SaveParametric",
    "-parametric"=>"Mesh.SaveParametric",
    "-save_topology"=>"Mesh.SaveTopology",
    # Upstream `-save_entities` aliases `-save_topology` (same CTX store).
    "-save_entities"=>"Mesh.SaveTopology",
    "-optimize_ho"=>"Mesh.HighOrderOptimize",
    "-hoOptimize"=>"Mesh.HighOrderOptimize",
    "-nopopup"=>"General.NoPopup",
    "-renumber"=>"Mesh.Renumber",
    # `Mesh.IgnorePeriodicity` gates periodic-section handling on MSH 2.x
    # reads; the flag stores 1 (the documented and compiled-in default).
    "-ignore_periodicity"=>"Mesh.IgnorePeriodicity",
    "-ignorePeriodicity"=>"Mesh.IgnorePeriodicity",
    # `MeshSizeFromCurvatureIsotropic` forces isotropic curvature sizing —
    # the only variant the 1-D curvature estimator computes, so the stored
    # option is consistent with the produced grading.
    "-clcurviso"=>"Mesh.MeshSizeFromCurvatureIsotropic",
    "-term"=>"General.Terminal",
    # `-savePartTopology` writes `<name>_topology.pro` beside a partitioned
    # MSH output (upstream `partitionSaveTopologyFile`).
    "-part_topo_pro"=>"Mesh.PartitionTopologyFile",
    "-savePartTopology"=>"Mesh.PartitionTopologyFile",
    # `Mesh.SwitchElementTags` swaps the elementary and physical tag fields
    # while an MSH 2.2 file is read (`Merge`, `API.open`/`API.merge`) —
    # upstream's `_storeElementInEntity` `reg`/`physical` exchange.
    "-switch_tags"=>"Mesh.SwitchElementTags",
    # `Mesh.PreserveNumberingMsh2` emits the stored element numbering on
    # MSH 2.2 output instead of sequential order positions.
    "-preserve_numbering_msh2"=>"Mesh.PreserveNumberingMsh2",
    "-preserveNumberingMsh2"=>"Mesh.PreserveNumberingMsh2")

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
    "-4"=>"adaptive meshing (batch mode 4) is not implemented",
    # Upstream `-part_weight` switches to partition algorithm 3 and feeds
    # per-element-family weights into its graph partitioner; no Tessella
    # partition producer consumes `Mesh.Partition*Weight`, so storing the
    # values would silently do nothing.
    "-part_weight"=>"weighted graph partitioning is not implemented",
    "-partWeight"=>"weighted graph partitioning is not implemented",

    "-part_physicals"=>"partition physical-group toggling is not implemented",
    "-part_no_physicals"=>"partition physical-group toggling is not implemented",
    # `Mesh.AllowSwapAngle` is consumed upstream only by the BDS
    # discrete-surface boundary-edge swap optimizer, which Tessella's
    # surface pipelines do not expose.
    "-swapangle"=>"boundary-edge swap optimization is not implemented",
    "-match"=>"geometry/mesh matching is not implemented",
    "-numsubedges"=>"high-order element display subdivision is GUI-only",
    "-debugSurface"=>"surface mesh debugging is not implemented",
    "-run"=>"ONELAB solver clients are not implemented",
    "-onelab"=>"ONELAB solver clients are not implemented",
    "-socket"=>"ONELAB solver clients are not implemented",
    "-listen"=>"solver socket listeners are not implemented",
    "-minterpreter"=>"Octave interpreters are not implemented",
    "-pyinterpreter"=>"Python interpreters are not implemented",
    "-watch"=>"file watching is GUI-only",
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
# rejected with an explicit "not an upstream option" message. `-claniso` and
# `-ignore_periocity` never existed: the upstream parser answers them with
# "Skipping unknown option" (`-ignore_periocity` survives only as a typo in
# the `-help` listing; the accepted spelling is `-ignore_periodicity`).
const _OBSOLETE_FLAGS=("-newreg","-newvol","-newsur","-post_pro","-clear",
    "-liststring","-listnumber","-min","-min0","-min1","-min2","-max",
    "-view","-size_field","-quadlayout","-autoex","-openFile",
    "-claniso","-ignore_periocity")

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
        # Upstream converts `.geo` inputs by writing the unrolled model
        # script — the same artifact `-0` produces.
        result=execute_geo(path;mesh_dim=0)
        write_geo_unrolled(destination,result.model)
        return destination
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
`-aniso_max`/`-anisoMax`, `-rand`, `-smooth_ratio`/`-smoothRatio`,
`-v`/`-verbose` and `-nt` set real numeric options; they retain priority over
script assignments, including option aliases. `-optimize`, `-optimize_netgen`,
`-save_all`/`-saveall`, `-renumber`, `-nopopup`, `-term`, `-clcurviso`,
`-ignore_periodicity`/`-ignorePeriodicity` and `-quad` set the matching
boolean/algorithm options; `-algo name` selects an implemented algorithm
(`del2d`, `front2d`, `del3d` and their upstream aliases). `-bgm file.pos`
loads a `.pos` view file before the input runs and installs its last view as
the `PostView` background-mesh field, on both the `.geo` and session paths.
`-setnumber`, `-setstring`, `-setlist`/`-setlistofnumbers`, `-string` and
`-option` inject parser input before a `.geo` file (on non-`.geo` inputs only
`-setnumber` option names apply, through `API.option`).
`-refine`, `-reclassify deg`, `-reparam deg`, `-save`, `-part n`,
`-part_topo`/`-part_no_topo`, `-part_ghosts`/`-part_no_ghosts`, `-open`,
`-merge`, `-new`, and `.msh`/`.pos` inputs run on the `API` session; `-part n`
partitions the current mesh after meshing (or the loaded mesh) exactly like
`API.mesh.partition` — `-part_topo`/`-part_ghosts` request partition
topology/ghost production and reject on unsupported inputs. `-part_split`
(aliases `-part_split_files`, `-oneFilePerPart`) swaps the combined output for
one MSH file per partition named `<output stem>_<p>.msh`: each carries the
partition's own cells (entities whose `partitions` list contains `p`) plus the
ghost cells it receives, emitted as regular cells on its ghost entity — the
upstream `Mesh.PartitionSplitMeshFiles` layout. Without `-part`, or when the
partitioned result holds a single partition, the flag is inert and the normal
output is written.
`-convert files...` rewrites `.msh`/`.pos` files to `file_new` (binary MSH 4.1
or parsed view text) and must come last.
`file.msh -check` validates the mesh without writing. `-`, `-parse_and_exit`,
`-help`/`--help`, `-version`/`--version`, `-info`/`--info`, `-help_options`
and `-pid` are informational; the first four must be used alone.
`-log file` mirrors every Msg-channel message into `file` (upstream's
`Msg::SetLogFileName`), and `-cpu` appends upstream's
`(From start: Wall %gs, CPU %gs, Mem %gMb)` resource suffix to Info and Error
messages (`Msg::SetInfoCpu`/`SetInfoMem`). `-statreport file` truncates `file`
and writes the 2-D statistics header plus one row after 2-D generation;
`-append_statreport file` appends the row only — both reproduce upstream's
`PrintMesh2dStatistics` layout. `-switch_tags` stores `Mesh.SwitchElementTags`
to swap elementary/physical fields on MSH 2.2 input, and
`-preserve_numbering_msh2`/`-preserveNumberingMsh2` stores
`Mesh.PreserveNumberingMsh2` to emit the stored element numbering on MSH 2.2
output.
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
    # `-cpu` resource reporting measures from the start of the batch run;
    # `-log`'s file lives for the whole `main` call. Both are process-global
    # like upstream's `Msg` state, so the `finally` restores them for callers
    # (API sessions, repeated `main` invocations) that never saw the flags.
    IO._MSG_START_NS[]=time_ns()
    try
        return _run(tokens)
    finally
        IO._msg_log_close()
        IO._MSG_INFO_RESOURCES[]=false
    end
end

function _run(tokens::Vector{String})
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
    refine_seen=false;save_seen=false;parse_seen=false;barycentric_seen=false
    unroll_seen=false
    reclassify=nothing;reparam=nothing
    part=nothing;part_topo=nothing;part_ghosts=nothing;part_split=nothing
    algo_seen=false;quad_seen=false
    # `-statreport`/`-append_statreport` — `(path, :write|:append)`; the last
    # occurrence wins like upstream's `createAppendMeshStatReport` store.
    statreport=nothing
    # `-bgm` files, in argument order: `(path, index of the file's last view)`.
    # The index is a 0-based position in the merged view list — only the
    # views these flags merge precede it, since the injected `Merge`s run
    # before every other prelude and before the script body.
    bgm_specs=Tuple{String,Int}[]
    bgm_view_base=0;bgm_seen=false
    i=1
    while i<=length(tokens)
        a=tokens[i]
        if a=="-convert"
            isempty(files) || throw(ArgumentError(
                "tessella: -convert cannot follow an input file"))
            i+1<=length(tokens) || throw(ArgumentError(
                "tessella: -convert needs at least one file"))
            # An explicit loop avoids the comprehension's element-type
            # inference — `collect` would query `return_types` on
            # `_convert_file`, forcing the `execute_geo` graph to infer.
            destination=nothing
            for token in tokens[i+1:end]
                destination=_convert_file(_path_argument(token,"-convert"))
            end
            return something(destination)
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
        elseif a=="-cpu"
            # `Msg::SetInfoCpu(true)` + `Msg::SetInfoMem(true)` — Info and
            # StatusBar messages carry the resource suffix. Upstream accepts
            # repeats; the flag is idempotent.
            IO._MSG_INFO_RESOURCES[]=true
        elseif a=="-log"
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: -log needs a file path"))
            # `Msg::SetLogFileName` closes a previously opened log file.
            IO._msg_log_close()
            IO._msg_log_open(_path_argument(tokens[i],"-log"))
        elseif a in ("-statreport","-append_statreport")
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: $a needs a file path"))
            statreport=(_path_argument(tokens[i],a),
                a=="-statreport" ? :write : :append)
        elseif a=="-refine"
            refine_seen && throw(ArgumentError("tessella: duplicate -refine flag"))
            refine_seen=true
        elseif a=="-barycentric_refine"
            barycentric_seen && throw(ArgumentError(
                "tessella: duplicate -barycentric_refine flag"))
            barycentric_seen=true
        elseif a=="-0"
            unroll_seen && throw(ArgumentError("tessella: duplicate -0 flag"))
            unroll_seen=true
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
        elseif a in ("-part_split","-part_split_files","-oneFilePerPart")
            part_split!==nothing && throw(ArgumentError(
                "tessella: duplicate or conflicting partition-split flag"))
            part_split=true
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
        elseif a=="-bgm"
            bgm_seen=true
            i+=1;i<=length(tokens) || throw(ArgumentError(
                "tessella: -bgm needs a background-mesh file"))
            bgm_file=_path_argument(tokens[i],"background mesh")
            lowercase(splitext(bgm_file)[2])==".pos" || throw(ArgumentError(
                "tessella: -bgm supports .pos background-mesh files only"))
            if !isfile(bgm_file)
                # Upstream reports `Invalid background mesh` and meshes on —
                # the flag's file is optional input, not a hard dependency.
                println(stderr,"tessella: -bgm file $bgm_file does not " *
                    "exist; continuing without a background mesh")
            else
                views=try
                    Post.read_pos(bgm_file)
                catch err
                    err isa InterruptException && rethrow()
                    throw(ArgumentError(
                        "tessella: -bgm could not parse $bgm_file"))
                end
                # A view-less file loads but selects nothing — upstream's
                # "Invalid background mesh (no view)" arm. The `Merge` still
                # runs on the standalone path.
                isempty(views) ||
                    push!(bgm_specs,
                        (abspath(bgm_file),bgm_view_base+length(views)-1))
                bgm_view_base+=length(views)
            end
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
        (reparam===nothing ? 0 : 1)+(barycentric_seen ? 1 : 0)+
        (unroll_seen ? 1 : 0)
    batch_actions>1 && throw(ArgumentError(
        "tessella: conflicting batch action flags"))
    inputs=[file for (file,_) in files]
    extension=lowercase(splitext(first(inputs))[2])
    if check_seen
        (dimension_seen || output_seen || format_seen || binary_seen ||
         !isempty(numeric_options) || !isempty(preludes) || bgm_seen ||
         refine_seen || barycentric_seen || unroll_seen || save_seen ||
         parse_seen ||
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
        barycentric_seen || unroll_seen ||
        reclassify!==nothing || reparam!==nothing || part!==nothing ||
        output_seen || format_seen || binary_seen) && throw(ArgumentError(
            "tessella: parse-and-exit cannot be combined with meshing or " *
            "output options"))
    # `.msh`/`.pos` inputs, multiple files, and post-mesh actions run on the
    # API session. High-order optimization lives in the session meshing plan,
    # so `-optimize_ho`/`-hoOptimize` routes there as well. Everything else —
    # including extensionless `.geo` sources — keeps the standalone `.geo`
    # execution path.
    session=extension in (".msh",".pos") || length(files)>1 || refine_seen ||
        barycentric_seen ||
        save_seen || reclassify!==nothing || reparam!==nothing ||
        part!==nothing ||
        get(numeric_options,"Mesh.HighOrderOptimize",0.0)!=0.0
    if !session
        if unroll_seen
            # `-0` writes the model's unrolled `.geo` script — the geometry
            # executes without meshing, and the default destination is the
            # input stem with a `.geo_unrolled` suffix (upstream convention).
            destination=outfile===nothing ?
                splitext(first(inputs))[1]*".geo_unrolled" : outfile
            _same_file_target(first(inputs),destination) && throw(ArgumentError(
                "tessella: -0 output path must not alias the input file"))
            result=execute_geo(first(inputs);mesh_dim=0)
            write_geo_unrolled(destination,result.model)
            return destination
        end
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
        # `-bgm` files load as `.pos` views before the input script runs —
        # upstream reads the background-mesh file at startup — and each file's
        # last view becomes a `PostView` background field. The injected
        # `Merge`/`Background Mesh View` statements run ahead of every other
        # prelude, so the recorded view index is exact.
        injected=preludes
        if !isempty(bgm_specs)
            statements=String[]
            for (bgm_path,bgm_index) in bgm_specs
                push!(statements,"Merge $(_escape_geo_string(
                    replace(bgm_path,'\\'=>'/'),"-bgm"));")
                push!(statements,"Background Mesh View[$bgm_index];")
            end
            injected=vcat(statements,preludes)
        end
        # The `.geo` numeric-option store is internal to the execution;
        # `Mesh.SaveParametric`/`Mesh.SaveTopology` values assigned inside the
        # script or by flags reach the writer through this recorder.
        save_options=Dict{String,Float64}()
        save_recorder=(family,member,value)->begin
            family=="Mesh" &&
                member in ("SaveParametric","SaveTopology","ScalingFactor",
                           "PartitionTopologyFile","PreserveNumberingMsh2") &&
                (save_options[member]=Float64(value))
            return nothing
        end
        exec_t0=time_ns()
        result=try
            if !isempty(injected)
                temporary=_inject_geo(source,injected)
                source=temporary
            end
            execute_geo(source;mesh_dim=dim,
                _forced_numeric_options=numeric_options,
                _initial_numeric_option_callback=save_recorder)
        finally
            temporary===nothing || rm(temporary;force=true)
        end
        save_parametric=get(save_options,"SaveParametric",nothing)
        save_topology=get(save_options,"SaveTopology",nothing)
        # `Mesh.PreserveNumberingMsh2` (`-preserve_numbering_msh2`/in-script
        # option) emits the stored element numbering on MSH 2.2 output.
        preserve_numbering=!iszero(Float64(get(save_options,
            "PreserveNumberingMsh2",0.0)))
        statreport!==nothing && _geo_stat_report(
            statreport[1],statreport[2],result,dim,first(inputs),
            Float64(time_ns()-exec_t0)/1.e9,
            Dict{String,Float64}("General.Verbosity"=>
                get(numeric_options,"General.Verbosity",5.0),
                "General.Terminal"=>
                get(numeric_options,"General.Terminal",1.0)))
        # `Mesh.ScalingFactor` (`-scale`/in-script option) multiplies the
        # serialized vertex coordinates — same hook upstream's writers take.
        scaling=Float64(get(save_options,"ScalingFactor",1.0))
        # `Mesh.PartitionTopologyFile` (`-savePartTopology`/in-script option)
        # writes `<name>_topology.pro` beside a partitioned MSH output.
        topology=!iszero(Float64(get(save_options,"PartitionTopologyFile",0.0)))
        emit_topology=(written)->begin
            topology && written isa MixedMesh &&
                written.partition_data!==nothing &&
                Int(written.partition_data.num_partitions)>0 &&
                write_partitioned_topology(
                    joinpath(dirname(destination),
                        splitext(basename(destination))[1]*"_topology.pro"),
                    written)
            return nothing
        end
        if dim>0
            result.mesh===nothing && throw(ErrorException("tessella: no mesh produced"))
            if dim==1
                authority=result.model.meshing.mesh_authority
                authority!==nothing && authority.mesh!==nothing || throw(ErrorException(
                    "tessella: one-dimensional generation omitted its classified mesh"))
                write_msh(destination,_scaled_mesh_copy(authority.mesh,scaling);
                          version,binary,
                          save_parametric,save_topology,preserve_numbering)
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
                # `-save_parametric`/`-save_topology` need classified
                # ownership/parameters — prefer the meshing authority's
                # native `MixedMesh` when one was produced.
                if save_parametric!==nothing || save_topology!==nothing
                    classified=_cli_classified_mesh(result)
                    if classified isa MixedMesh
                        classified_out=_scaled_mesh_copy(classified,scaling)
                        write_msh(destination,classified_out;version,binary,
                                  save_parametric,save_topology,
                                  preserve_numbering)
                        emit_topology(classified_out)
                    else
                        (version==2.2 &&
                            (save_parametric==1.0 || save_topology==1.0)) &&
                            throw(ArgumentError(
                                "tessella: -save_parametric/-save_topology " *
                                "requires a classified mesh (the run produced " *
                                "no entity metadata)"))
                        binary ?
                            write_mixed_msh(destination,
                                simplex_to_mixed(_scaled_mesh_copy(result.mesh,scaling));
                                version,binary,preserve_numbering) :
                            write_msh(destination,_scaled_mesh_copy(result.mesh,scaling);
                                      version)
                    end
                elseif result.mesh isa MixedMesh
                    mixed_out=_scaled_mesh_copy(result.mesh,scaling)
                    write_mixed_msh(destination,mixed_out;version,binary,
                        preserve_numbering)
                    emit_topology(mixed_out)
                elseif binary
                    write_mixed_msh(destination,
                        simplex_to_mixed(_scaled_mesh_copy(result.mesh,scaling));
                        version,binary,preserve_numbering)
                else
                    write_msh(destination,_scaled_mesh_copy(result.mesh,scaling);version)
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
                projected_out=_scaled_mesh_copy(projected,scaling)
                write_msh(destination,projected_out;
                          version,binary,
                          save_parametric,save_topology,preserve_numbering)
                emit_topology(projected_out)
            end
            return destination
        end
        return first(inputs)
    end
    return _session_main(files,preludes,numeric_options;
        dim,outfile,version,binary,output_flags=format_seen || binary_seen,
        refine=refine_seen,barycentric=barycentric_seen,save=save_seen,
        statreport,
        unroll=unroll_seen,bgm=[path for (path,_) in bgm_specs],
        reclassify,reparam,part,part_topo,part_ghosts,part_split)
end

# Classified `MixedMesh` preferred for `-save_parametric`/`-save_topology`
# output: the meshing authority's native record retains entity ownership and
# parametric node coordinates the merged `result.mesh` may not carry.
function _cli_classified_mesh(result)
    authority=result.model.meshing.mesh_authority
    authority!==nothing && authority.mesh isa MixedMesh && return authority.mesh
    result.mesh isa MixedMesh && return result.mesh
    return result.mesh
end

# ---------------------------------------------------------------------------
# `-statreport` / `-append_statreport`
#
# Upstream's `PrintMesh2dStatistics` runs at the end of the 2-D generation
# pass: `-statreport` truncates the report file and writes the header (header
# only when the model is empty), `-append_statreport` appends one row with no
# header. Tessella emits the same row layout, filled only with values the
# produced mesh honestly carries: `#faces`/`#fail` count the model's surfaces
# and the ones the run left unmeshed, `#t` and the `Q` aggregates evaluate
# upstream's `qmTriangle::gamma` (`2·r_in/r_circ`) per triangle, `#Q>90`
# counts `gamma > 0.9`, `#e` is the surface triangulation's unique-edge
# count, and `CPU` is the generation wall time. Upstream never assigns the
# per-face edge-length and efficiency-index fields (the oracle emits whatever
# the uninitialized memory held — typically zeros, occasionally garbage), so
# `tau`/`#Egood`/`#Egood/#e` stay 0 as "not computed" markers.
# ---------------------------------------------------------------------------

@inline _stat_sub(p,q)=(p[1]-q[1],p[2]-q[2],p[3]-q[3])
@inline _stat_norm(v)=hypot(v[1],v[2],v[3])
@inline _stat_cross(a,b)=(a[2]*b[3]-a[3]*b[2],
                         a[3]*b[1]-a[1]*b[3],
                         a[1]*b[2]-a[2]*b[1])
@inline _stat_edge(i,j)=i<j ? (Int32(i),Int32(j)) : (Int32(j),Int32(i))

function _stat_unit(v)
    magnitude=_stat_norm(v)
    magnitude>0 || return nothing
    return (v[1]/magnitude,v[2]/magnitude,v[3]/magnitude)
end

# `qmTriangle::gamma` — `2 * (2·sinA·sinB·sinC / (sinA+sinB+sinC))` on the
# corner sines, i.e. `2·r_in/r_circ`.
function _stat_tri_gamma(pa,pb,pc)
    a=_stat_unit(_stat_sub(pc,pb));b=_stat_unit(_stat_sub(pa,pc))
    c=_stat_unit(_stat_sub(pb,pa))
    (a===nothing || b===nothing || c===nothing) && return 0.0
    sina=_stat_norm(_stat_cross(b,c));sinb=_stat_norm(_stat_cross(c,a))
    sinc=_stat_norm(_stat_cross(a,b))
    sina==0.0 && sinb==0.0 && sinc==0.0 && return 0.0
    return 2*(2*sina*sinb*sinc/(sina+sinb+sinc))
end

# Triangle index triples in `mesh`, keyed by owner entity tag (0 when the
# mesh carries no entity classification).
function _stat_surface_tris(mesh)
    tris=Dict{Int32,Vector{NTuple{3,Int32}}}()
    if mesh isa Mesh
        for j in axes(mesh.tris,2)
            push!(get!(()->NTuple{3,Int32}[],tris,Int32(0)),
                (Int32(mesh.tris[1,j]),Int32(mesh.tris[2,j]),
                 Int32(mesh.tris[3,j])))
        end
        return mesh.coords,tris
    elseif mesh isa MixedMesh
        data=mesh.entity_data
        for (bi,block) in enumerate(mesh.blocks)
            block isa ElementBlock || continue
            spec=msh_spec(block.msh)
            (spec.dim==2 && spec.family===:tri) || continue
            owners=data===nothing ? nothing : data.block_entities[bi]
            for j in axes(block.nodes,2)
                push!(get!(()->NTuple{3,Int32}[],tris,
                        owners===nothing ? Int32(0) : Int32(owners[j])),
                    (Int32(block.nodes[1,j]),Int32(block.nodes[2,j]),
                     Int32(block.nodes[3,j])))
            end
        end
        return mesh.coords,tris
    end
    return nothing
end

# Append gamma qualities and unique edges for every triangle in `mesh` and
# return the per-owner triangle counts keyed like `_stat_surface_tris`.
# `allowed` restricts accumulation to the listed owner tags (`nothing` keeps
# every cell); counts still report every owner so callers can count failures.
function _stat_accumulate_part!(qualities::Vector{Float64},
                                edges::Set{NTuple{2,Int32}},
                                mesh;allowed=nothing)
    collected=_stat_surface_tris(mesh)
    collected===nothing && return Dict{Int32,Int}()
    coords,tris=collected
    counts=Dict{Int32,Int}()
    for (owner,cells) in tris
        counts[owner]=get(counts,owner,0)+length(cells)
        (allowed===nothing || owner in allowed) || continue
        for (i1,i2,i3) in cells
            push!(qualities,_stat_tri_gamma(
                coords[:,Int(i1)],coords[:,Int(i2)],coords[:,Int(i3)]))
            push!(edges,_stat_edge(i1,i2),_stat_edge(i2,i3),
                  _stat_edge(i3,i1))
        end
    end
    return counts
end

const _STAT_HEADER="2D stats\tname\t\t#faces\t\t#fail\t\t#t\t\t" *
    "Qavg\t\tQbest\t\tQworst\t\t#Q>90\t\t#Q>90/#t\t#e\t\ttau\t\t" *
    "#Egood\t\t#Egood/#e\tCPU"

# Write (or append) the single 2-D statistics row. `face_counts` maps model
# surface tags to their produced triangle counts; `unattributed` holds counts
# for cells whose owner could not be resolved — they contribute to the global
# `#t`/`Q`/`#e` aggregates and count as meshed faces only when every surface
# would otherwise report zero.
function _emit_stat_report(path::String,mode::Symbol,name::String,
                           faces::Vector{Int},face_counts::Dict{Int,Int},
                           unattributed::Int,qualities::Vector{Float64},
                           n_edges::Int,wall::Float64,model_empty::Bool,
                           options::Dict{String,Float64})
    io=try
        open(path,mode===:append ? "a" : "w")
    catch err
        err isa InterruptException && rethrow()
        # Upstream reports the failure on the Msg channel and keeps going.
        API._native_detached_diagnostic_emit(:error,
            "Could not open file '$path'",options)
        return nothing
    end
    newline=IO._MSG_LOG_NEWLINE
    try
        if mode===:write
            print(io,_STAT_HEADER,newline)
            model_empty && return nothing
        end
        nt=length(qualities)
        # Upstream divides the accumulated stats by the produced counts —
        # a zero-element pass prints NaN, so the honest emission does too.
        avg=nt>0 ? sum(qualities)/nt : NaN
        best=nt>0 ? maximum(qualities) : 0.0
        worst=nt>0 ? minimum(qualities) : 1.0
        n_good=count(>(0.9),qualities)
        nfail=count(tag->get(face_counts,tag,0)==0,faces)
        # Unattributed cells land on no face; when they are the only product
        # the honest statement is that the face(s) meshed.
        unattributed>0 && nfail==length(faces) && (nfail=0)
        n_e=n_edges
        # The upstream efficiency-index and good-length fields are never
        # populated (uninitialized memory); emit 0 as the uncomputed marker.
        print(io,@sprintf("\t%16s\t%d\t\t%d\t\t",name,length(faces),nfail))
        print(io,@sprintf("%d\t\t%8.7f\t%8.7f\t%8.7f\t%d\t\t%8.7f\t",
            nt,avg,best,worst,n_good,
            nt>0 ? n_good/nt : NaN))
        print(io,@sprintf("%d\t\t%8.7f\t%d\t\t%8.7f\t%8.1f",
            n_e,0.0,0,n_e>0 ? 0.0 : NaN,wall),newline)
        # `Msg::Info("Efficiency index for surface mesh tau=%g ", ...)`.
        API._native_detached_diagnostic_emit(:info,
            @sprintf("Efficiency index for surface mesh tau=%g ",0.0),
            options)
        return nothing
    finally
        close(io)
    end
end

# `.geo` path: per-face products come from `result.actual_mesh_parts` via
# `geo_entity_mesh`, which throws for a face the run left unmeshed.
function _geo_stat_report(path::String,mode::Symbol,result,dim::Int,
                          input::String,wall::Float64,
                          options::Dict{String,Float64})
    faces=sort!(Int[tag for tag in keys(result.model.surfaces)])
    # The report is produced only when the 2-D pass ran: `-2`/`-3` requests
    # it, and an in-script `Mesh 2` produces dim-2 actual parts.
    produced=dim>=2 || any(part->part[1]==2,result.actual_mesh_parts)
    produced || return nothing
    model=result.model
    model_empty=isempty(model.points) && isempty(model.curves) &&
        isempty(model.surfaces) && isempty(model.volumes) &&
        isempty(model.discrete)
    qualities=Float64[];edges=Set{NTuple{2,Int32}}()
    face_counts=Dict{Int,Int}()
    for tag in faces
        part=try
            geo_entity_mesh(result,2,tag)
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        # Every cell of an entity part belongs to `tag` regardless of the
        # part's internal owner keys (dense `Mesh` parts carry none).
        counts=part===nothing ? Dict{Int32,Int}() :
            _stat_accumulate_part!(qualities,edges,part)
        face_counts[tag]=sum(values(counts);init=0)
    end
    name=isempty(model.name) ? splitext(basename(input))[1] : model.name
    _emit_stat_report(path,mode,name,faces,face_counts,0,
        qualities,length(edges),wall,model_empty,options)
    return nothing
end

# Session path: the generated mesh's own entity classification gives the
# per-face attribution; the model's dim-2 entity list gives `#faces`. Discrete
# entities are skipped like upstream's `firstFace()` CAD-only iteration.
function _session_stat_report(path::String,mode::Symbol,wall::Float64,
                              options::Dict{String,Float64})
    mesh,faces,model_empty=lock(API.STATE_LOCK) do
        current=API._model_locked()
        entities=model_entities(current)
        face_tags=sort!(Int[Int(tag) for (d,tag) in entities
            if d==2 && !haskey(current.discrete,(d,tag))])
        (API.LAST_MESH[],face_tags,isempty(entities))
    end
    qualities=Float64[];edges=Set{NTuple{2,Int32}}()
    # Upstream aggregates only the non-discrete faces' statistics — cells
    # owned by discrete or non-surface entities do not enter the row.
    # Owner 0 (unclassified cells) has no upstream analog and stays honest.
    allowed=union!(Set{Int32}(Int32.(faces)),Int32[0])
    counts=mesh===nothing ? Dict{Int32,Int}() :
        _stat_accumulate_part!(qualities,edges,mesh;allowed)
    face_counts=Dict{Int,Int}();unattributed=0
    for tag in faces
        face_counts[tag]=get(counts,Int32(tag),0)
    end
    unattributed=get(counts,Int32(0),0)
    _emit_stat_report(path,mode,API.MODEL_NAME[],faces,face_counts,
        unattributed,qualities,length(edges),wall,model_empty,options)
    return nothing
end

# `-barycentric_refine` refines the session's current mesh in place — the
# Gmsh 4.15.2 action subdivides only linear simplices owned by entities at
# the model dimension (triangles on 2-D models, tetrahedra on 3-D models);
# every other cell passes through. The transform runs on the stored mesh
# product and the stale index caches are invalidated so `API.write`
# serializes the refined mesh.
function _session_barycentric_refine!()
    lock(API.STATE_LOCK) do
        cached=API.LAST_MESH[]
        # Upstream refines whatever mesh is loaded; an unmeshed model
        # (e.g. a `.geo` opened without generation) is a no-op.
        cached===nothing && return nothing
        refined=(cached isa MixedMesh || cached isa Mesh) ?
                barycentric_refine(cached) : throw(ArgumentError(
                    "tessella: -barycentric_refine cannot refine $(typeof(cached))"))
        refined===cached && return nothing
        # Rebuild the session classification from the refined mesh's own
        # entity metadata — `API.write` rejects an unclassified cache.
        class=nothing
        if refined isa MixedMesh && refined.entity_data!==nothing
            class=API._higher_classification_from_data(
                refined,"tessella: -barycentric_refine")
            class=API._mixed_rebind_class(class,refined;
                node_parametric=refined.entity_data.node_parametric)
        end
        API._replace_mesh_cache_locked!(refined,class)
        API.LAST_MESH_LOCATOR[]=nothing
        API.LAST_MESH_EDGES[]=nothing
        API.LAST_MESH_FACES[]=nothing
        isdefined(API,:LAST_MESH_HIGH_ORDER) &&
            (API.LAST_MESH_HIGH_ORDER[]=nothing)
        return nothing
    end
end

# Session path: files route through `API.open`/`API.merge`/`API.model.add`,
# option writes through `API.option`, meshing through `API.mesh.generate`,
# post-mesh actions through `API.mesh.refine`/`classify_surfaces`/
# `create_geometry`/`partition`, and output through `API.write`.
function _session_main(files,preludes,numeric_options;
                       dim,outfile,version,binary,output_flags,refine,save,
                       reclassify,reparam,part,part_topo,part_ghosts,part_split,
                       barycentric=false,unroll=false,bgm=String[],
                       statreport=nothing)
    need_output=dim>0 || refine || barycentric || unroll || save ||
        reclassify!==nothing ||
        reparam!==nothing || part!==nothing
    !need_output && outfile!==nothing && throw(ArgumentError(
        "tessella: -o requires a meshing, conversion or save action"))
    !need_output && output_flags && throw(ArgumentError(
        "tessella: output format flags require a meshing, conversion or save action"))
    # `Mesh.SaveParametric`/`Mesh.SaveTopology`/`Mesh.HighOrderNumLayers` live
    # in the `.geo` option table but outside `API.DEFAULT_OPTIONS`; the save
    # pair is applied to the serialized file after `API.write` and all three
    # are recorded in the session option table directly — so `-setnumber`
    # preludes behave like session options rather than `.geo` injections.
    save_option_values=Dict{String,Float64}()
    session_save_names=("Mesh.SaveParametric","Mesh.SaveTopology",
                        "Mesh.HighOrderNumLayers","Mesh.PreserveNumberingMsh2")
    # Options consumed while files are opened/merged must be in `API.OPTIONS`
    # before the file loop — `Mesh.SwitchElementTags` swaps MSH 2.2
    # elementary/physical fields during import.
    session_preopen_names=("Mesh.SwitchElementTags",)
    preopen_values=Dict{String,Float64}()
    # Literal `.geo` statements need a leading `.geo` file; `(name, value)`
    # option preludes become session options when the session table knows the
    # name and `.geo` assignments otherwise (`General.NoPopup` and friends).
    first_is_geo=lowercase(splitext(files[1][1])[2])==".geo"
    option_preludes=Tuple{String,Float64}[]
    lines=String[]
    for prelude in preludes
        if prelude isa Tuple
            if prelude[1] in session_preopen_names
                preopen_values[prelude[1]]=prelude[2]
            elseif prelude[1] in session_save_names
                save_option_values[prelude[1]]=prelude[2]
            elseif haskey(API.DEFAULT_OPTIONS,prelude[1])
                push!(option_preludes,prelude)
            else
                push!(lines,"$(prelude[1]) = $(prelude[2]);")
            end
        else
            push!(lines,prelude)
        end
    end
    !first_is_geo && !isempty(lines) && throw(ArgumentError(
        "tessella: -setstring/-setlist/-string/-option inject .geo parser " *
        "input and require a .geo input file"))
    for (name,value) in numeric_options
        name in session_preopen_names && (preopen_values[name]=value;continue)
        name in session_save_names &&
            (save_option_values[name]=value)
    end
    save_parametric=get(save_option_values,"Mesh.SaveParametric",nothing)
    save_topology=get(save_option_values,"Mesh.SaveTopology",nothing)
    save_flags_set=save_parametric!==nothing || save_topology!==nothing
    # Session-table options apply through `API.option`; options that live only
    # in the `.geo` table (`General.NoPopup`, `General.Terminal`) are injected
    # into a leading `.geo` file and are inert on non-.geo inputs.
    session_numeric=Dict{String,Float64}()
    geo_numeric_preludes=Union{String,Tuple{String,Float64}}[]
    for (name,value) in numeric_options
        name in session_preopen_names && continue
        name in session_save_names && continue
        haskey(API.DEFAULT_OPTIONS,name) ? (session_numeric[name]=value) :
            push!(geo_numeric_preludes,"$name = $value;")
    end
    # `-0` always defaults to `<stem>.geo_unrolled`, even on mesh inputs.
    destination=outfile===nothing ?
        (unroll ? splitext(files[1][1])[1]*".geo_unrolled" :
         _default_output(files[1][1])) : outfile
    extension=lowercase(splitext(files[1][1])[2])
    if extension!=".geo" && need_output && !unroll
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
        # Options consumed during open/merge (e.g. `Mesh.SwitchElementTags`
        # steering the MSH 2.2 reader) land before the first file is read.
        for (name,value) in preopen_values
            API.OPTIONS[name]=value
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
        # `-bgm`: merge the `.pos` views and point a `PostView` background
        # field at each file's last view — the same semantics the standalone
        # path's `Merge` + `Background Mesh View[i]` prelude expresses.
        for bgm_path in bgm
            API.merge(bgm_path)
            count=length(API.SESSION_VIEWS)
            count>0 || continue
            tag=API.mesh.field.add("PostView")
            API.mesh.field.set_number(tag,"ViewIndex",count-1)
            API.mesh.field.set_as_background_mesh(tag)
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
        # `.geo`-table options the session table does not declare are stored
        # directly — the same place `API.option` would have put them.
        for (name,value) in save_option_values
            API.OPTIONS[name]=value
        end
        gen_t0=statreport!==nothing && dim>=2 ? time_ns() : UInt64(0)
        dim>0 && API.mesh.generate(dim)
        if statreport!==nothing && dim>=2
            _session_stat_report(statreport[1],statreport[2],
                Float64(time_ns()-gen_t0)/1.e9,copy(API.OPTIONS))
        end
        refine && API.mesh.refine()
        barycentric && _session_barycentric_refine!()
        reclassify===nothing ||
            API.mesh.classify_surfaces(reclassify*pi/180,true,false,pi)
        if reparam!==nothing
            API.mesh.classify_surfaces(reparam*pi/180,true,true,pi)
            API.mesh.create_geometry()
        end
        part===nothing || API.mesh.partition(part)
        if unroll
            # `-0` on a mesh input serializes the loaded model's discrete
            # point records — the only entities the unrolled format carries.
            mesh=lock(API.STATE_LOCK) do
                API.LAST_MESH[]
            end
            mesh isa MixedMesh || throw(ArgumentError(
                "tessella: -0 on a mesh input requires entity metadata"))
            write_geo_unrolled(destination,mesh)
            return destination
        end
        if need_output
            # `-part_split` only takes effect on an active partition request —
            # upstream's Mesh.PartitionSplitMeshFiles is inert otherwise.
            out_extension=lowercase(splitext(destination)[2])
            if part_split===true && part!==nothing
                _write_partition_files(destination;version,binary,
                    save_parametric,save_topology)
            elseif save_flags_set &&
                   out_extension in (".msh",".msh2",".msh22",".msh4",".msh41")
                write_version=out_extension in (".msh2",".msh22") ? 2.2 :
                              out_extension in (".msh4",".msh41") ? 4.1 :
                              version
                API.option("Mesh.MshFileVersion",write_version)
                parent=dirname(destination)
                isempty(parent) && (parent=".")
                tmp,tmpio=mktemp(parent;cleanup=false)
                close(tmpio)
                tmpmsh=tmp*".msh"
                try
                    API.write(tmpmsh)
                    source=lock(API.STATE_LOCK) do
                        API.LAST_MESH[]
                    end
                    _msh_rewrite_save_options(tmpmsh;
                        version=write_version,binary=binary,
                        save_parametric=save_parametric,
                        save_topology=save_topology,
                        source=source)
                    mv(tmpmsh,destination;force=true)
                finally
                    isfile(tmpmsh) && rm(tmpmsh;force=true)
                end
            else
                API.write(destination)
            end
        end
    finally
        temporary===nothing || rm(temporary;force=true)
        API.finalize()
    end
    return need_output ? destination : files[1][1]
end

# ---------------------------------------------------------------------------
# Per-partition output: -part_split / -part_split_files / -oneFilePerPart
#
# Upstream `Mesh.PartitionSplitMeshFiles` replaces the combined write with one
# `<stem>_<p>.msh` file per partition — the stem derives from `-o` or the
# input name and the extension is always `.msh`. Each file keeps the full
# entity and partition declarations but only carries the cells owned by that
# partition — entities whose `partitions` list contains it — plus the ghost
# cells it receives, emitted as regular cells on its ghost entity. Its
# `$GhostElements` section lists those copies with their owner partitions.
# ---------------------------------------------------------------------------

_partition_ncells(block::ElementBlock)=size(block.nodes,2)
_partition_ncells(block::SpecialElementBlock)=length(block.tags)
_partition_dim(block::ElementBlock)=msh_dimension(block.msh)
_partition_dim(block::SpecialElementBlock)=
    Elements.MSH_SPECIAL_TYPES[block.msh].dim
_partition_nodes(block::ElementBlock,cell::Int)=block.nodes[:,cell]
_partition_nodes(block::SpecialElementBlock,cell::Int)=
    view(block.connectivity,
         Int(block.offsets[cell]):Int(block.offsets[cell+1])-1)

# Flat element index (1:N across blocks, the `$GhostElements` convention) to
# its (block, cell) pair.
function _partition_flat_cell(starts::Vector{Int},index::Int)
    for block in 1:length(starts)-1
        index<starts[block+1] && return (block,index-starts[block]+1)
    end
    return (length(starts)-1,index-starts[end-1]+1)
end

# Extract partition `partition` of a partitioned `MixedMesh` as a standalone
# mesh view: the partition's own cells plus its received ghost cells placed as
# regular cells on its ghost entity. Node data sections are dropped — their
# dense element/node payloads cannot be re-indexed by a subset extraction —
# while ancillary sections, periodic links and physical names carry over.
function _partition_view(mesh::MixedMesh,partition::Int;
                         partition_metadata::Bool=true)
    data=mesh.entity_data;pdata=mesh.partition_data
    caller="tessella: -part_split"
    data===nothing && throw(ArgumentError(
        "$caller: the partitioned mesh carries no entity classification"))
    pdata===nothing && throw(ArgumentError(
        "$caller: the mesh carries no partition entities"))
    count=Int(pdata.num_partitions)
    1<=partition<=count || throw(ArgumentError(
        "$caller: partition $partition is outside 1:$count"))
    # Ghost entities live at the mesh's top dimension; their declared entity
    # rows carry the `(dim, tag)` placement upstream keeps in ghost records.
    topdim=maximum((_partition_dim(block) for block in mesh.blocks);init=0)
    ghost_owner=Dict{Tuple{Int,Int},Int}()    # (dim,tag) -> receiving partition
    ghost_key=Dict{Int,Tuple{Int,Int}}()      # partition -> (dim,tag)
    for (tag,part) in pdata.ghost_entities
        dim=topdim
        for d in 0:3
            (haskey(data.entities,(d,Int(tag))) ||
             haskey(pdata.entities,(d,Int(tag)))) && (dim=d;break)
        end
        ghost_owner[(dim,Int(tag))]=Int(part)
        ghost_key[Int(part)]=(dim,Int(tag))
    end
    starts=Vector{Int}(undef,length(mesh.blocks)+1);starts[1]=1
    for (bi,block) in pairs(mesh.blocks)
        starts[bi+1]=starts[bi]+_partition_ncells(block)
    end
    nel=starts[end]-1
    selected=[falses(_partition_ncells(block)) for block in mesh.blocks]
    for (bi,block) in enumerate(mesh.blocks)
        entities=data.block_entities[bi];dim=_partition_dim(block)
        for j in 1:_partition_ncells(block)
            key=(dim,Int(entities[j]))
            record=get(pdata.entities,key,nothing)
            if record===nothing
                # Cells parked on this partition's ghost entity are already
                # placed ghost copies — re-splitting a split file keeps them.
                get(ghost_owner,key,0)==partition && (selected[bi][j]=true)
            else
                Int32(partition) in record.partitions && (selected[bi][j]=true)
            end
        end
    end
    ghost_cells=Tuple{Int,Int}[]
    for record in pdata.ghost_elements
        Int32(partition) in record.ghost_partitions || continue
        source=Int(record.element)
        1<=source<=nel || throw(ArgumentError(
            "$caller: ghost element index $source is outside 1:$nel"))
        (bi,j)=_partition_flat_cell(starts,source)
        key=(_partition_dim(mesh.blocks[bi]),Int(data.block_entities[bi][j]))
        if get(ghost_owner,key,0)==partition || selected[bi][j]
            continue   # placed ghost copy or own cell — already counted
        end
        push!(ghost_cells,(bi,j))
    end
    # Retained nodes cover every selected and duplicated cell.
    used=falses(size(mesh.coords,2))
    for (bi,block) in enumerate(mesh.blocks),j in eachindex(selected[bi])
        selected[bi][j] || continue
        for node in _partition_nodes(block,j);used[Int(node)]=true;end
    end
    for (bi,j) in ghost_cells
        for node in _partition_nodes(mesh.blocks[bi],j);used[Int(node)]=true;end
    end
    rows=findall(used)
    inverse=zeros(Int32,size(mesh.coords,2))
    inverse[rows].=Int32.(eachindex(rows))
    # New block/cell positions for link remapping and ghost-record indices:
    # `ghost_target` maps each source flat index to its flat index in the view.
    block_index=Dict{Int,Int}();cell_index=Dict{Tuple{Int,Int},Int32}()
    ghost_target=Dict{Int,Int}();flat=0
    for (bi,block) in enumerate(mesh.blocks)
        any(selected[bi]) || continue
        block_index[bi]=length(block_index)+1;position=0
        for j in eachindex(selected[bi])
            selected[bi][j] || continue
            position+=1;flat+=1
            cell_index[(bi,j)]=Int32(position)
            ghost_target[starts[bi]+j-1]=flat
        end
    end
    blocks=MixedElementBlock[]
    owners=Vector{Int32}[];element_tags=Vector{UInt64}[]
    elementary=Vector{Int32}[]
    function remap_ref(ref)
        ref==ElementRef() && return ElementRef()
        target=(Int(ref.block),Int(ref.cell))
        haskey(cell_index,target) || return ElementRef()
        return ElementRef(block_index[target[1]],cell_index[target])
    end
    for (bi,block) in enumerate(mesh.blocks)
        columns=findall(selected[bi]);isempty(columns) && continue
        if block isa ElementBlock
            push!(blocks,ElementBlock(block.msh,
                inverse[block.nodes[:,columns]],block.tags[columns]))
        else
            cells=[inverse[collect(_partition_nodes(block,j))] for j in columns]
            parents=ElementRef[remap_ref(block.parent_refs[j]) for j in columns]
            domains=Matrix{ElementRef}(undef,2,length(columns))
            for side in 1:2,(cj,j) in enumerate(columns)
                domains[side,cj]=remap_ref(block.domain_refs[side,j])
            end
            push!(blocks,SpecialElementBlock(block.msh,cells,block.tags[columns];
                parent_refs=parents,domain_refs=domains))
        end
        push!(owners,data.block_entities[bi][columns])
        push!(element_tags,data.external_element_tags[bi][columns])
        mesh.elementary_entities===nothing ||
            push!(elementary,mesh.elementary_entities[bi][columns])
    end
    # The view's entity table holds the source's global entities — upstream
    # keeps the original `$Entities` section — plus the ghost entity the
    # partition's ghost copies land on (Tessella declares ghost entities
    # there; a Gmsh file that omitted it gets a synthesized record). MSH 2.x
    # output carries no `$PartitionedEntities`, so partition entities must be
    # merged in for cell-entity resolution.
    entities=copy(data.entities)
    partition_metadata || merge!(entities,pdata.entities)
    if !isempty(ghost_cells)
        key=get(ghost_key,partition,nothing)
        key===nothing && throw(ArgumentError(
            "$caller: partition $partition receives ghost cells but has no " *
            "declared ghost entity"))
        gdim,gtag=key
        order=Int[];grouped=Dict{Int,Vector{Tuple{Int,Int}}}()
        for (bi,j) in ghost_cells
            msh=mesh.blocks[bi].msh
            msh in order || push!(order,msh)
            push!(get!(()->Tuple{Int,Int}[],grouped,msh),(bi,j))
        end
        for msh in order
            cells=grouped[msh]
            if mesh.blocks[cells[1][1]] isa ElementBlock
                conn=hcat((inverse[_partition_nodes(mesh.blocks[bi],j)]
                           for (bi,j) in cells)...)
                push!(blocks,ElementBlock(msh,conn))
            else
                conn=[inverse[collect(_partition_nodes(mesh.blocks[bi],j))]
                      for (bi,j) in cells]
                push!(blocks,SpecialElementBlock(msh,conn))
            end
            push!(owners,fill(Int32(gtag),length(cells)))
            push!(element_tags,
                [data.external_element_tags[bi][j] for (bi,j) in cells])
            mesh.elementary_entities===nothing ||
                push!(elementary,fill(Int32(gtag),length(cells)))
            for (bi,j) in cells
                flat+=1;ghost_target[starts[bi]+j-1]=flat
            end
        end
        # Declare the ghost entity when the entity table lacks it, computing
        # its bounds from the placed ghost cells' coordinates.
        haskey(entities,key) || begin
            nodes=Int[]
            for (bi,j) in ghost_cells
                append!(nodes,Int.(_partition_nodes(mesh.blocks[bi],j)))
            end
            bounds=isempty(nodes) ? ntuple(_->0.0,6) : ntuple(
                axis->axis<=3 ?
                    minimum(mesh.coords[axis,row] for row in nodes) :
                    maximum(mesh.coords[axis-3,row] for row in nodes),6)
            entities[key]=MixedEntity(gdim,gtag,bounds)
        end
    end
    metadata=MixedEntityData(entities;
        node_entities=data.node_entities[rows],
        node_parametric=data.node_parametric[rows],
        external_node_tags=data.external_node_tags[rows],
        block_entities=owners,external_element_tags=element_tags,
        _native_mesh_tags=data.native_mesh_tags)
    # Upstream writes one record per ghost copy in the file: the owner
    # partition plus a single-destination list naming this file's partition.
    records=MixedGhostElement[]
    for record in pdata.ghost_elements
        Int32(partition) in record.ghost_partitions || continue
        push!(records,MixedGhostElement(
            ghost_target[Int(record.element)],record.partition,
            Int32[partition]))
    end
    partition2=partition_metadata ?
        MixedPartitionData(pdata.num_partitions;entities=pdata.entities,
            ghost_entities=pdata.ghost_entities,ghost_elements=records) :
        nothing
    links=MixedPeriodicLink[]
    for link in mesh.periodic_links
        all(node->used[Int(node)],link.slave_nodes) &&
        all(node->used[Int(node)],link.master_nodes) || continue
        push!(links,MixedPeriodicLink(
            link.dim,link.slave_entity,link.master_entity,
            inverse[link.slave_nodes],inverse[link.master_nodes];
            affine=link.affine))
    end
    # MSH 2.x has no partition metadata; carry the entity tags as elementary
    # tags so ghost cells still land on the ghost entity there.
    elementary_out=mesh.elementary_entities===nothing ?
        (partition_metadata ? nothing : copy.(owners)) : elementary
    return MixedMesh(mesh.coords[:,rows],blocks;
        entity_data=metadata,physical_names=mesh.physical_names,
        elementary_entities=elementary_out,
        periodic_links=links,ancillary_sections=mesh.ancillary_sections,
        partition_data=partition2)
end

# Map each exported cell to its assigned partition label. `element_partitions`
# indexes the cached mesh's flat element order; exported cells translate
# through their public element tags so reordering and `Mesh.SaveAll` filtering
# keep the mapping exact.
function _partition_label_assignment(mesh::MixedMesh,partition,public)
    data=mesh.entity_data
    data===nothing && return copy(partition.element_partitions)
    uid_index=Dict{UInt64,Int}()
    for index in eachindex(partition.element_partitions)
        uid_index[API._cache_element_tag(public,index)]=index
    end
    labels=Vector{Int32}(undef,
        sum(_partition_ncells,mesh.blocks;init=0))
    flat=0
    for (bi,block) in enumerate(mesh.blocks)
        for uid in data.external_element_tags[bi]
            flat+=1
            index=get(uid_index,uid,0)
            index>0 || throw(ArgumentError(
                "tessella: -part_split found no partition label for " *
                "element $uid"))
            labels[flat]=partition.element_partitions[index]
        end
    end
    return labels
end

# Label-mode `-part` carries no partition entities: each split file simply
# holds the elements assigned to the partition on their original entities.
# `owner_remap` folds a produced partition-child owner back to its declared
# parent so the file never references entities it does not carry.
function _partition_label_view(mesh::MixedMesh,labels::Vector{Int32},
                               partition::Int;v2::Bool=false,owner_remap=nothing)
    data=mesh.entity_data
    remap_owner(entity)=owner_remap===nothing ? entity :
        Int32(get(owner_remap,entity,(entity[1],entity[2]))[2])
    selected=[falses(_partition_ncells(block)) for block in mesh.blocks]
    flat=0
    for (bi,block) in enumerate(mesh.blocks),j in 1:_partition_ncells(block)
        flat+=1;Int32(partition)==labels[flat] && (selected[bi][j]=true)
    end
    used=falses(size(mesh.coords,2))
    for (bi,block) in enumerate(mesh.blocks),j in eachindex(selected[bi])
        selected[bi][j] || continue
        for node in _partition_nodes(block,j);used[Int(node)]=true;end
    end
    rows=findall(used)
    inverse=zeros(Int32,size(mesh.coords,2))
    inverse[rows].=Int32.(eachindex(rows))
    blocks=MixedElementBlock[]
    owners=Vector{Int32}[];element_tags=Vector{UInt64}[]
    elementary=Vector{Int32}[]
    for (bi,block) in enumerate(mesh.blocks)
        columns=findall(selected[bi]);isempty(columns) && continue
        if block isa ElementBlock
            push!(blocks,ElementBlock(block.msh,
                inverse[block.nodes[:,columns]],block.tags[columns]))
        else
            cells=[inverse[collect(_partition_nodes(block,j))] for j in columns]
            push!(blocks,SpecialElementBlock(block.msh,cells,
                block.tags[columns]))
        end
        data===nothing || begin
            dim=_partition_dim(block)
            push!(owners,[remap_owner((dim,Int(owner)))
                for owner in data.block_entities[bi][columns]])
            push!(element_tags,data.external_element_tags[bi][columns])
        end
        mesh.elementary_entities===nothing || begin
            dim=_partition_dim(block)
            push!(elementary,[remap_owner((dim,Int(owner)))
                for owner in mesh.elementary_entities[bi][columns]])
        end
    end
    node_owners=data===nothing ? nothing : data.node_entities[rows]
    if owner_remap!==nothing && node_owners!==nothing
        node_owners=[(Int(owner[1]),remap_owner((Int(owner[1]),Int(owner[2]))))
            for owner in node_owners]
    end
    metadata=data===nothing ? nothing : MixedEntityData(data.entities;
        node_entities=node_owners,
        node_parametric=data.node_parametric[rows],
        external_node_tags=data.external_node_tags[rows],
        block_entities=owners,external_element_tags=element_tags,
        _native_mesh_tags=data.native_mesh_tags)
    links=MixedPeriodicLink[]
    for link in mesh.periodic_links
        all(node->used[Int(node)],link.slave_nodes) &&
        all(node->used[Int(node)],link.master_nodes) || continue
        push!(links,MixedPeriodicLink(
            link.dim,link.slave_entity,link.master_entity,
            inverse[link.slave_nodes],inverse[link.master_nodes];
            affine=link.affine))
    end
    elementary_out=mesh.elementary_entities===nothing ?
        (v2 && data!==nothing ? copy.(owners) : nothing) : elementary
    return MixedMesh(mesh.coords[:,rows],blocks;
        entity_data=metadata,physical_names=mesh.physical_names,
        elementary_entities=elementary_out,
        periodic_links=links,ancillary_sections=mesh.ancillary_sections)
end

# A produced partition whose entity layer is only `createNewEntities`' bare
# ownership children (no interface topology, no ghost layer) splits exactly
# like the label record did: every cell folds to its declared parent entity
# and each per-partition file carries no partition metadata.
function _partition_bare_labels(mesh::MixedMesh,pdata)
    data=mesh.entity_data
    data===nothing && return nothing
    isempty(pdata.ghost_entities) && isempty(pdata.ghost_elements) ||
        return nothing
    remap=Dict{Tuple{Int,Int},Tuple{Int,Int32}}()
    for (key,entity) in pdata.entities
        isempty(entity.boundaries) || return nothing
        entity.parent[1]==key[1] || return nothing
        length(entity.partitions)==1 || return nothing
        remap[key]=(key[1],entity.parent[2])
    end
    labels=Vector{Int32}(undef,sum(_partition_ncells,mesh.blocks;init=0))
    flat=0
    for (bi,block) in enumerate(mesh.blocks)
        owners=data.block_entities[bi];dim=_partition_dim(block)
        for j in 1:_partition_ncells(block)
            flat+=1
            record=get(pdata.entities,(dim,Int(owners[j])),nothing)
            record===nothing && return nothing
            labels[flat]=record.partitions[1]
        end
    end
    return remap,labels
end

# `-part_split`: write `stem_1.msh`…`stem_N.msh` instead of `destination`. The
# combined model write is staged into a throwaway file so the split works on
# exactly what `API.write` would serialize — MSH 4.1 keeps partition and ghost
# metadata through that read-back regardless of the requested output version.
function _write_partition_files(destination::String;version,binary,
                                save_parametric=nothing,save_topology=nothing)
    parent=dirname(abspath(destination))
    isdir(parent) || throw(ArgumentError(
        "tessella: -part_split output directory does not exist: $parent"))
    temporary=tempname(parent)*".msh"
    mesh=try
        API.option("Mesh.MshFileVersion",4.1)
        API.write(temporary)
        read_mixed_msh(temporary)
    finally
        API.option("Mesh.MshFileVersion",version)
        rm(temporary;force=true)
    end
    stem,_=splitext(destination)
    pdata=mesh.partition_data
    bare=pdata===nothing ? nothing : _partition_bare_labels(mesh,pdata)
    if pdata===nothing || bare!==nothing
        npart,labels=if pdata===nothing
            partition,public=lock(API.STATE_LOCK) do
                record=API.LAST_MESH_PARTITION[]
                (record===nothing || record.mesh!==API.LAST_MESH[]) &&
                    return (nothing,nothing)
                return (record,API._mesh_public_tags(record.mesh))
            end
            partition===nothing && throw(ArgumentError(
                "tessella: -part_split found no serializable partition state"))
            partition.num_partitions,
            _partition_label_assignment(mesh,partition,public)
        else
            Int(pdata.num_partitions),bare[2]
        end
        npart>1 || return API.write(destination)
        for part in 1:Int(npart)
            view=_partition_label_view(mesh,labels,part;v2=version!=4.1,
                owner_remap=bare===nothing ? nothing : bare[1])
            write_msh(stem*"_$part.msh",view;version,binary,
                      save_parametric,save_topology)
        end
        return destination
    end
    pdata.num_partitions>1 || return API.write(destination)
    for part in 1:Int(pdata.num_partitions)
        view=_partition_view(mesh,part;partition_metadata=version==4.1)
        write_msh(stem*"_$part.msh",view;version,binary,
                  save_parametric,save_topology)
    end
    return destination
end

end # module
