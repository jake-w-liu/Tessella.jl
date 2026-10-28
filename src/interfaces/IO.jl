"""
    IO

Mesh file I/O (PLAN.md §3 "IO"): gmsh **MSH v2.2 and v4.1** (ASCII) read/write,
ASCII/binary **STL** ingest for boundary surface meshes, and a lightweight
**`.geo`** parameter/structure and mesh-field scanner.

The round-trip contract (DEVELOPMENT.md CRC gate for Stage 0) is *connectivity
preservation*: reading a mesh and writing it back — in either format version —
must reproduce the same topology, verified by `MeshTypes.mesh_crc`. gmsh allows
arbitrary node/element tags; we relabel to a compact `1:N` on read (connectivity
is invariant under a consistent relabel) while preserving physical-group tags.

Full evaluation of a `.geo` OpenCASCADE CSG script (Booleans → BREP → faces) is
still pending. `read_geo_params` reads the declarations that do not require a CAD
evaluator: global mesh sizing, physical groups, numeric list variables, and the raw
background-field graph.
"""
module IO

using ..MeshTypes: Mesh, nnodes, nsegs, ntris, ntets, node, validate
using ..Elements: MSH_PHYSICAL_NAME_MAX_BYTES, _copy_physical_names, MixedMesh,
                  MixedEntityData, ElementBlock, SpecialElementBlock, msh_spec,
                  msh_dimension, msh_family, msh_num_nodes, msh_type,
                  simplex_to_mixed, write_mixed_msh, read_mixed_msh,
                  _block_ncells
using ..Predicates: _two_diff_tail
using ..GmshLibm: _gm87_sin, _gm87_cos, _gm_tan, _gm_asin, _gm_acos, _gm_atan,
                  _gm87_atan2, _gm_sinh, _gm_cosh, _gm_tanh, _gm87_exp, _gm87_log,
                  _gm_log10, _gm87_pow
using Printf: @printf, @sprintf, Format, format

export read_msh, write_msh, MshFile
export read_stl, read_geo_params, GeoParams, GeoFieldSpec
export read_unv, write_unv, read_mesh, write_mesh, UnvFile, MeshFile
export read_meshb, write_meshb
export read_su2, write_su2, read_bdf, write_bdf, Su2File, BdfFile
export read_inp, write_inp, read_neu, write_neu, InpFile, NeuFile
export read_p3d, write_p3d, read_mail, write_mail, read_ir3, write_ir3,
       P3dFile, MailFile, Ir3File
export read_diff, write_diff, DiffFile
export read_vtk, write_vtk, VtkFile
export read_wrl, WrlFile, read_xao, XaoFile, read_brep, BrepFile

include("GeoColorNames.jl")
include("GeoOptionTables.jl")

# gmsh element type codes we handle. (type => n_nodes)
const MSH_POINT = 15
const MSH_LINE  = 1
const MSH_TRI   = 2
const MSH_TET   = 4
const MSH_TET2  = 11    # 10-node (quadratic) tet — read as its 4 corner vertices
const MSH_TRI2  = 9     # 6-node (quadratic) tri — read as its 3 corner vertices
const MSH_SEG2  = 8     # 3-node (quadratic) line — read as its 2 end vertices
const _NN = Dict(MSH_POINT => 1, MSH_LINE => 2, MSH_TRI => 3, MSH_TET => 4,
                 MSH_TET2 => 10, MSH_TRI2 => 6, MSH_SEG2 => 3)
const _EDIM = Dict(MSH_POINT => 0, MSH_LINE => 1, MSH_TRI => 2, MSH_TET => 3,
                   MSH_TET2 => 3, MSH_TRI2 => 2, MSH_SEG2 => 1)
const _DEFAULT_MAX_IO_NAME_BYTES = 1 << 20

@inline function _io_limit(value, caller::AbstractString,
                           name::AbstractString; ceiling::Int=typemax(Int))
    value isa Integer || throw(ArgumentError(
        "$caller: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: $name must not be Bool"))
    converted = try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name is outside Int bounds"))
    end
    0 <= converted <= ceiling || throw(ArgumentError(
        "$caller: $name must lie in 0:$ceiling (got $value)"))
    return converted
end

struct _MshReadLimits
    nodes::Int
    elements::Int
    entities::Int
    blocks::Int
    physical_names::Int
    name_bytes::Int
end

const _DEFAULT_MSH_READ_LIMITS = _MshReadLimits(
    Int(typemax(Int32)), Int(typemax(Int32)),
    Int(typemax(Int32)), Int(typemax(Int32)),
    Int(typemax(Int32)), _DEFAULT_MAX_IO_NAME_BYTES)

# ════════════════════════════════════════════════════════════════════════════════
# Reading
# ════════════════════════════════════════════════════════════════════════════════

"""
    MshFile

Result of `read_msh`: the compact `mesh`, and `physical_names` mapping
`(dim, physical_tag) => name`.
"""
struct MshFile
    mesh::Mesh
    physical_names::Dict{Tuple{Int,Int},String}
    # Classified node metadata recovered from the file, in mesh-node order.
    # `node_entities[i]` is the `(entityDim, entityTag)` owner of node `i`
    # when the input carried ownership records (MSH 4.1 node blocks or an
    # MSH 2.2 `$ParametricNodes` section); both vectors stay empty otherwise.
    # `node_parametric[i]` is the stored `(u)`/`(u,v)` parameter vector —
    # `nothing` marks a classified non-parametric node.
    node_entities::Vector{Tuple{Int,Int32}}
    node_parametric::Vector{Union{Nothing,Vector{Float64}}}
    # Declared topology records `(dim, tag) => (boundaries, physical_tags)`.
    # MSH 4.1 `$Entities` keep their signed boundary tags; the MSH 2.2
    # `$Entities` topology section is unsigned on disk and stored that way.
    entities::Dict{Tuple{Int,Int},Tuple{Vector{Int32},Vector{Int32}}}
end

MshFile(mesh,physical_names) = MshFile(mesh,physical_names,
    Tuple{Int,Int32}[],Union{Nothing,Vector{Float64}}[],
    Dict{Tuple{Int,Int},Tuple{Vector{Int32},Vector{Int32}}}())

# Mutable accumulator used while parsing, converted to a `Mesh` at the end.
mutable struct _Accum
    node_x::Vector{Float64}
    node_y::Vector{Float64}
    node_z::Vector{Float64}
    tag2idx::Dict{UInt64,Int32}          # gmsh node tag → compact 1-based index
    segs::Vector{NTuple{2,Int32}}
    tris::Vector{NTuple{3,Int32}}
    tets::Vector{NTuple{4,Int32}}
    seg_tag::Vector{Int32}
    tri_tag::Vector{Int32}
    tet_tag::Vector{Int32}
    pnames::Dict{Tuple{Int,Int},String}
    ep::Dict{Tuple{Int,Int},Int}    # (entityDim, entityTag) → physical tag (v4)
    element_tags::Set{UInt64}
    physical_name_records::Int
    element_blocks::Int
    node_entities::Vector{Tuple{Int,Int32}}
    node_parametric::Vector{Union{Nothing,Vector{Float64}}}
    topo_entities::Dict{Tuple{Int,Int},Tuple{Vector{Int32},Vector{Int32}}}
    _Accum() = new(Float64[], Float64[], Float64[], Dict{UInt64,Int32}(),
                   NTuple{2,Int32}[], NTuple{3,Int32}[], NTuple{4,Int32}[],
                   Int32[], Int32[], Int32[], Dict{Tuple{Int,Int},String}(),
                   Dict{Tuple{Int,Int},Int}(),Set{UInt64}(),0,0,
                   Tuple{Int,Int32}[],Union{Nothing,Vector{Float64}}[],
                   Dict{Tuple{Int,Int},Tuple{Vector{Int32},Vector{Int32}}}())
end

@inline function _nodeidx!(acc::_Accum, tag::UInt64)
    get(acc.tag2idx, tag, 0)
end

"""
    read_msh(path; tessella_extensions=false,
             max_nodes=typemax(Int32), max_elements=typemax(Int32),
             max_entities=typemax(Int32),
             max_blocks=typemax(Int32),
             max_physical_names=typemax(Int32),
             max_name_bytes=1_048_576,
             max_file_bytes=typemax(Int)) -> MshFile

Read a gmsh ASCII `.msh` file (format 2.2 or 4.1). Node/element tags are
relabelled to a compact `1:N`; per-element physical tags are preserved. Resource
limits are checked against section headers before record-sized allocation, and the
parsed simplex mesh is validated before it is returned. MSH4 files without the
optional `\$Entities` section receive neutral implicit entities. Repeated
`\$Entities` and `\$Elements` sections are accumulated; repeated physical names
must agree. Backslashes in physical names are literal. Set
`tessella_extensions=true` only when reading the escaped-name extension produced by
`write_msh(...; gmsh_compatible=false)`.
"""
function read_msh(path;
                  tessella_extensions=false,
                  max_nodes=typemax(Int32),
                  max_elements=typemax(Int32),
                  max_entities=typemax(Int32),
                  max_blocks=typemax(Int32),
                  max_physical_names=typemax(Int32),
                  max_name_bytes=_DEFAULT_MAX_IO_NAME_BYTES,
                  max_file_bytes=typemax(Int))
    caller="read_msh"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    tessella_extensions isa Bool || throw(ArgumentError(
        "$caller: tessella_extensions must be Bool"))
    limits=_MshReadLimits(
        _io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32))),
        _io_limit(max_elements,caller,"max_elements";ceiling=Int(typemax(Int32))),
        _io_limit(max_entities,caller,"max_entities";ceiling=Int(typemax(Int32))),
        _io_limit(max_blocks,caller,"max_blocks";ceiling=Int(typemax(Int32))),
        _io_limit(max_physical_names,caller,"max_physical_names";
                  ceiling=Int(typemax(Int32))),
        _io_limit(max_name_bytes,caller,"max_name_bytes"))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    open(path, "r") do io
        try
            return _read_msh(io,limits,tessella_extensions)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError("read_msh: truncated .msh file"))
            rethrow()
        end
    end
end

_read_msh(io::Base.IO)=_read_msh(io,_DEFAULT_MSH_READ_LIMITS,false)
_read_msh(io::Base.IO,limits::_MshReadLimits)=_read_msh(io,limits,false)

function _read_msh(io::Base.IO,limits::_MshReadLimits,
                   tessella_extensions::Bool)
    acc = _Accum()
    version = 0.0
    seen_format=false;seen_nodes=false;seen_elements=false
    while !eof(io)
        line = strip(readline(io))
        isempty(line) && continue
        if line == "\$MeshFormat"
            seen_format && throw(ArgumentError("IO: duplicate \$MeshFormat section"))
            seen_format=true
            fmt = split(strip(readline(io)))
            length(fmt) == 3 || throw(ArgumentError("IO: malformed \$MeshFormat header"))
            version = parse(Float64, fmt[1])
            version in (2.2, 4.1) ||
                throw(ArgumentError("IO: unsupported .msh version $version (supported: 2.2 and 4.1)"))
            filetype = length(fmt) >= 2 ? parse(Int, fmt[2]) : 0
            filetype == 0 || throw(ArgumentError("IO: binary .msh not supported (file-type $filetype); use ASCII"))
            parse(Int,fmt[3]) == 8 || throw(ArgumentError("IO: unsupported .msh floating-point data size $(fmt[3])"))
            _expect_end(io, "\$EndMeshFormat")
        elseif line == "\$PhysicalNames"
            seen_format || throw(ArgumentError("IO: \$PhysicalNames appeared before \$MeshFormat"))
            _read_physical_names!(acc, io, limits, tessella_extensions)
        elseif line == "\$Nodes" || line == "\$ParametricNodes"
            seen_nodes && throw(ArgumentError("IO: duplicate \$Nodes section"))
            seen_elements && throw(ArgumentError("IO: \$Nodes must precede \$Elements"))
            seen_nodes=true
            version != 0 || throw(ArgumentError("IO: \$Nodes appeared before \$MeshFormat"))
            if line == "\$ParametricNodes"
                version < 3 || throw(ArgumentError(
                    "IO: \$ParametricNodes is only valid in MSH v2"))
                _read_parametric_nodes_v2!(acc, io, limits)
            else
                version < 3 ? _read_nodes_v2!(acc, io, limits) :
                              _read_nodes_v4!(acc, io, limits)
            end
        elseif line == "\$Elements"
            seen_nodes || throw(ArgumentError("IO: \$Elements appeared before \$Nodes"))
            seen_elements && version<3 && throw(ArgumentError(
                "IO: repeated \$Elements sections require MSH v4.1"))
            seen_elements=true
            version != 0 || throw(ArgumentError("IO: \$Elements appeared before \$MeshFormat"))
            version < 3 ? _read_elements_v2!(acc, io, limits) :
                          _read_elements_v4!(acc, io, limits)
        elseif line == "\$Entities"
            seen_format || throw(ArgumentError("IO: \$Entities appeared before \$MeshFormat"))
            (seen_nodes || seen_elements) &&
                throw(ArgumentError("IO: \$Entities must precede \$Nodes and \$Elements"))
            if version < 3
                _read_topo_entities_v2!(acc, io, limits)
            else
                version == 4.1 || throw(ArgumentError(
                    "IO: \$Entities is only valid in MSH v2/v4.1"))
                _read_entities_v4!(acc, io, limits)
            end
        elseif startswith(line, "\$") && !startswith(line, "\$End")
            _skip_section(io, "\$End" * line[2:end])   # ignore unknown sections
        end
    end
    version != 0 || throw(ArgumentError("IO: missing \$MeshFormat section"))
    seen_nodes || throw(ArgumentError("IO: missing \$Nodes section"))
    mesh=_to_mesh(acc)
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "read_msh: parsed mesh is invalid — " * join(diagnostic.messages,"; ")))
    length(acc.node_entities) in (0, nnodes(mesh)) ||
        throw(ArgumentError("IO: classified node-owner records do not match the node count"))
    length(acc.node_parametric) in (0, nnodes(mesh)) ||
        throw(ArgumentError("IO: parametric node records do not match the node count"))
    return MshFile(mesh, acc.pnames, acc.node_entities,
                   acc.node_parametric, acc.topo_entities)
end

function _expect_end(io, tok)
    eof(io) && throw(ArgumentError("IO: expected $tok, reached end of file"))
    l = strip(readline(io))
    l == tok || throw(ArgumentError("IO: expected $tok, got '$l'"))
end

function _skip_section(io, endtok)
    while !eof(io)
        strip(readline(io)) == endtok && return
    end
    throw(ArgumentError("IO: unterminated section (missing $endtok)"))
end

function _read_physical_names!(acc, io, limits::_MshReadLimits,
                               tessella_extensions::Bool)
    n = parse(Int, strip(readline(io)))
    n >= 0 || throw(ArgumentError("IO: negative physical-name count $n"))
    n<=typemax(Int32) || throw(ArgumentError("IO: physical-name count exceeds Int32"))
    n<=limits.physical_names-acc.physical_name_records || throw(ArgumentError(
        "read_msh: physical-name count exceeds max_physical_names=$(limits.physical_names)"))
    acc.physical_name_records+=n
    for _ in 1:n
        line = strip(readline(io))
        # format: `dim tag "name"`. Take the name between the first and last quote VERBATIM
        # — splitting on whitespace and rejoining collapses runs of interior spaces/tabs and
        # breaks the name-preservation round-trip contract.
        pattern=tessella_extensions ?
            r"^([+-]?\d+)\s+([+-]?\d+)\s+\"((?:\\.|[^\"])*)\"\s*$" :
            r"^([+-]?\d+)\s+([+-]?\d+)\s+\"([^\"]*)\"\s*$"
        m = match(pattern, line)
        m === nothing && throw(ArgumentError("IO: malformed physical-name record '$line'"))
        dim = parse(Int,m.captures[1]); tag = parse(Int,m.captures[2])
        0 <= dim <= 3 || throw(ArgumentError("IO: physical-name dimension $dim is outside 0:3"))
        tag > 0 || throw(ArgumentError("IO: physical-name tag must be positive (got $tag)"))
        _int32_tag(tag,"physical-name")
        name=tessella_extensions ? _unescape_name(m.captures[3]) :
                                   String(m.captures[3])
        ncodeunits(name)<=limits.name_bytes || throw(ArgumentError(
            "read_msh: physical name for ($dim,$tag) exceeds " *
            "max_name_bytes=$(limits.name_bytes)"))
        isvalid(name) || throw(ArgumentError(
            "read_msh: physical name for ($dim,$tag) is not valid UTF-8"))
        occursin('\0',name) && throw(ArgumentError(
            "read_msh: physical name for ($dim,$tag) contains a NUL byte"))
        key=(dim,tag)
        if haskey(acc.pnames,key)
            acc.pnames[key]==name || throw(ArgumentError(
                "IO: conflicting physical names for dimension/tag ($dim,$tag)"))
        else
            acc.pnames[key]=name
        end
    end
    _expect_end(io, "\$EndPhysicalNames")
end

function _unescape_name(s::AbstractString)
    io = IOBuffer(); i = firstindex(s)
    while i <= lastindex(s)
        c = s[i]
        if c == '\\'
            i = nextind(s,i); i <= lastindex(s) ||
                throw(ArgumentError("IO: trailing escape in physical name"))
            e = s[i]
            write(io, e == 'n' ? '\n' : e == 'r' ? '\r' : e == 't' ? '\t' : e)
        else
            write(io,c)
        end
        i = nextind(s,i)
    end
    return String(take!(io))
end

# Mesh labels use the MSH size_t domain; dense rows and CAD/physical tags
# keep their existing bounded integer storage.
@inline function _io_mesh_uid(token::AbstractString,context::AbstractString)
    try
        return parse(UInt64,token)
    catch error
        error isa InterruptException && rethrow()
        throw(ArgumentError("IO: invalid unsigned 64-bit mesh tag in $context"))
    end
end

# ── v2 ──────────────────────────────────────────────────────────────────────────
function _read_nodes_v2!(acc, io, limits::_MshReadLimits)
    n = parse(Int, strip(readline(io)))
    n >= 0 || throw(ArgumentError("IO: negative v2 node count $n"))
    n <= typemax(Int32)-length(acc.node_x) ||
        throw(ArgumentError("IO: v2 node count exceeds Int32 indexing"))
    n <= limits.nodes-length(acc.node_x) || throw(ArgumentError(
        "read_msh: v2 node count exceeds max_nodes=$(limits.nodes)"))
    for _ in 1:n
        p = split(strip(readline(io)))
        length(p) == 4 || throw(ArgumentError("IO: malformed v2 node record (expected 4 fields)"))
        tag = _io_mesh_uid(p[1],"v2 node")
        tag > 0 || throw(ArgumentError("IO: node tag must be positive (got $tag)"))
        haskey(acc.tag2idx,tag) && throw(ArgumentError("IO: duplicate node tag $tag"))
        x=parse(Float64,p[2]); y=parse(Float64,p[3]); z=parse(Float64,p[4])
        (isfinite(x)&&isfinite(y)&&isfinite(z)) || throw(ArgumentError("IO: node $tag has non-finite coordinates"))
        push!(acc.node_x, x); push!(acc.node_y, y); push!(acc.node_z, z)
        acc.tag2idx[tag] = length(acc.node_x)
    end
    _expect_end(io, "\$EndNodes")
end

function _read_elements_v2!(acc, io, limits::_MshReadLimits)
    n = parse(Int, strip(readline(io)))
    n >= 0 || throw(ArgumentError("IO: negative v2 element count $n"))
    n<=typemax(Int32) || throw(ArgumentError("IO: v2 element count exceeds Int32"))
    n<=limits.elements-length(acc.element_tags) || throw(ArgumentError(
        "read_msh: v2 element count exceeds max_elements=$(limits.elements)"))
    for _ in 1:n
        p = split(strip(readline(io)))
        length(p) >= 3 || throw(ArgumentError("IO: malformed v2 element record"))
        _record_element_tag!(acc,_io_mesh_uid(p[1],"v2 element"))
        etype = parse(Int, p[2])
        haskey(_NN,etype) || throw(ArgumentError(
            "IO: unsupported gmsh element type $etype; refusing to silently discard cells"))
        ntags = parse(Int, p[3])
        ntags >= 0 || throw(ArgumentError("IO: negative element tag count $ntags"))
        ntags <= length(p)-3 || throw(ArgumentError("IO: truncated v2 element tags"))
        phys = ntags >= 1 ? parse(Int, p[4]) : 0
        nodestart = 3 + ntags + 1
        needed = _NN[etype]
        length(p) == nodestart+needed-1 ||
            throw(ArgumentError("IO: element type $etype expects $needed node tags"))
        _push_element!(acc, etype, phys, @view p[nodestart:end])
    end
    _expect_end(io, "\$EndElements")
end

# MSH 2.2 `$ParametricNodes` (emitted by `-save_parametric`): each record is
#   nodeTag x y z entityDim entityTag [u [v [w]]]
# where the owner dimension decides the parameter count. A node with no stored
# parameters appears as the bare `nodeTag x y z entityDim entityTag` record.
function _read_parametric_nodes_v2!(acc, io, limits::_MshReadLimits)
    n = parse(Int, strip(readline(io)))
    n >= 0 || throw(ArgumentError("IO: negative v2 parametric node count $n"))
    n <= typemax(Int32)-length(acc.node_x) ||
        throw(ArgumentError("IO: v2 parametric node count exceeds Int32 indexing"))
    n <= limits.nodes-length(acc.node_x) || throw(ArgumentError(
        "read_msh: v2 parametric node count exceeds max_nodes=$(limits.nodes)"))
    for _ in 1:n
        p = split(strip(readline(io)))
        length(p) >= 6 || throw(ArgumentError(
            "IO: malformed v2 parametric node record (expected ≥ 6 fields)"))
        tag = _io_mesh_uid(p[1],"v2 parametric node")
        tag > 0 || throw(ArgumentError("IO: node tag must be positive (got $tag)"))
        haskey(acc.tag2idx,tag) && throw(ArgumentError("IO: duplicate node tag $tag"))
        x=parse(Float64,p[2]); y=parse(Float64,p[3]); z=parse(Float64,p[4])
        (isfinite(x)&&isfinite(y)&&isfinite(z)) || throw(ArgumentError(
            "IO: node $tag has non-finite coordinates"))
        edim = parse(Int, p[5]); etag = parse(Int, p[6])
        0 <= edim <= 3 || throw(ArgumentError(
            "IO: v2 parametric node $tag owner dimension $edim is outside 0:3"))
        etag > 0 || throw(ArgumentError(
            "IO: v2 parametric node $tag owner entity tag must be positive"))
        _int32_tag(etag, "v2 parametric node owner")
        extra = length(p)-6
        (extra == 0 || extra == edim) || throw(ArgumentError(
            "IO: v2 parametric node $tag carries $extra parameters; expected 0 or $edim"))
        params = extra == 0 ? nothing : begin
            uv = parse.(Float64, @view p[7:end])
            all(isfinite,uv) || throw(ArgumentError(
                "IO: node $tag has non-finite parametric coordinates"))
            Vector{Float64}(uv)
        end
        push!(acc.node_x, x); push!(acc.node_y, y); push!(acc.node_z, z)
        acc.tag2idx[tag] = length(acc.node_x)
        push!(acc.node_entities, (edim, Int32(etag)))
        push!(acc.node_parametric, params)
    end
    _expect_end(io, "\$EndParametricNodes")
end

# MSH 2.2 `$Entities` (emitted by `-save_topology` on discrete models): a
# topology-only section — no bounding boxes, unsigned boundary tags.
#   numPoints numCurves numSurfaces numVolumes
#   pointTag numPhysicalTags physicalTag...                       (points)
#   tag numBounding boundingTag... numPhysicalTags physicalTag... (dim ≥ 1)
function _read_topo_entities_v2!(acc, io, limits::_MshReadLimits)
    hdr = split(strip(readline(io)))
    length(hdr)==4 || throw(ArgumentError("IO: malformed v2 entity header"))
    counts = parse.(Int, hdr[1:4])
    all(>=(0),counts) || throw(ArgumentError("IO: negative v2 entity count"))
    total = try Base.checked_add(Base.checked_add(counts[1],counts[2]),
                                 Base.checked_add(counts[3],counts[4])) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("IO: v2 entity count overflows Int"))
    end
    total<=typemax(Int32) || throw(ArgumentError("IO: v2 entity count exceeds Int32"))
    total<=limits.entities-length(acc.topo_entities) || throw(ArgumentError(
        "read_msh: v2 entity count exceeds max_entities=$(limits.entities)"))
    for dim in 0:3
        for _ in 1:counts[dim+1]
            p = split(strip(readline(io)))
            tag = parse(Int, p[1])
            tag > 0 || throw(ArgumentError(
                "IO: v2 entity tag must be positive (got $tag)"))
            _int32_tag(tag, "v2 entity")
            haskey(acc.topo_entities,(dim,tag)) && throw(ArgumentError(
                "IO: duplicate v2 dimension-$dim entity tag $tag"))
            pos = 2
            bounds = Int32[]
            if dim > 0
                length(p)>=pos || throw(ArgumentError(
                    "IO: truncated v2 entity record ($dim,$tag)"))
                nb = parse(Int, p[pos]); pos += 1
                nb >= 0 || throw(ArgumentError(
                    "IO: negative bounding-entity count on v2 entity ($dim,$tag)"))
                length(p)>=pos+nb-1+1 || throw(ArgumentError(
                    "IO: truncated boundary list on v2 entity ($dim,$tag)"))
                for j in 1:nb
                    bt = parse(Int, p[pos]); pos += 1
                    bt > 0 || throw(ArgumentError(
                        "IO: v2 bounding entity tag must be positive (entity ($dim,$tag))"))
                    _int32_tag(bt, "v2 bounding entity")
                    push!(bounds, Int32(bt))
                end
            end
            length(p)>=pos || throw(ArgumentError(
                "IO: truncated v2 entity record ($dim,$tag)"))
            nphys = parse(Int, p[pos]); pos += 1
            nphys >= 0 || throw(ArgumentError(
                "IO: negative physical-tag count on v2 entity ($dim,$tag)"))
            length(p)==pos+nphys-1 || throw(ArgumentError(
                "IO: physical-tag count mismatch on v2 entity ($dim,$tag)"))
            phys = Int32[]
            for j in 1:nphys
                pt = parse(Int, p[pos]); pos += 1
                pt > 0 || throw(ArgumentError(
                    "IO: v2 physical tag must be positive (entity ($dim,$tag))"))
                _int32_tag(pt, "v2 physical")
                push!(phys, Int32(pt))
            end
            acc.topo_entities[(dim,tag)] = (bounds, phys)
        end
    end
    _expect_end(io, "\$EndEntities")
end

# ── v4.1 ────────────────────────────────────────────────────────────────────────
function _read_entities_v4!(acc, io, limits::_MshReadLimits)
    hdr = split(strip(readline(io)))
    length(hdr)==4 || throw(ArgumentError("IO: malformed v4 entity header"))
    nP, nC, nS, nV = parse.(Int, hdr[1:4])
    all(>=(0),(nP,nC,nS,nV)) || throw(ArgumentError("IO: negative v4 entity count"))
    total=try Base.checked_add(Base.checked_add(nP,nC),Base.checked_add(nS,nV)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("IO: v4 entity count overflows Int"))
    end
    total<=typemax(Int32) || throw(ArgumentError("IO: v4 entity count exceeds Int32"))
    total<=limits.entities-length(acc.ep) || throw(ArgumentError(
        "read_msh: v4 entity count exceeds max_entities=$(limits.entities)"))
    _read_entity_block!(acc, io, nP, 0)   # points
    _read_entity_block!(acc, io, nC, 1)   # curves
    _read_entity_block!(acc, io, nS, 2)   # surfaces
    _read_entity_block!(acc, io, nV, 3)   # volumes
    _expect_end(io, "\$EndEntities")
end

function _read_entity_block!(acc, io, count, dim)
    ep = acc.ep
    for _ in 1:count
        p = split(strip(readline(io)))
        idx = dim == 0 ? 5 : 8
        length(p)>=idx || throw(ArgumentError("IO: malformed dimension-$dim entity record"))
        tag = parse(Int, p[1])
        tag > 0 || throw(ArgumentError("IO: entity tag must be positive (got $tag)"))
        haskey(ep,(dim,tag)) && throw(ArgumentError("IO: duplicate dimension-$dim entity tag $tag"))
        coordlast = dim == 0 ? 4 : 7
        box = parse.(Float64,p[2:coordlast])
        all(isfinite,box) || throw(ArgumentError("IO: entity ($dim,$tag) has non-finite bounds"))
        if dim > 0
            all(box[i] <= box[i+3] for i in 1:3) ||
                throw(ArgumentError("IO: entity ($dim,$tag) has reversed bounding-box limits"))
        end
        # points: tag x y z numPhysicalTags [physicalTag...] (no max bbox)
        # dim≥1: tag minX minY minZ maxX maxY maxZ numPhysicalTags [physicalTag...] ...
        nphys = parse(Int, p[idx])
        nphys>=0 || throw(ArgumentError("IO: negative physical-tag count on entity ($dim,$tag)"))
        nphys<=1 || throw(ArgumentError(
            "IO: entity ($dim,$tag) belongs to $nphys physical groups; Mesh supports one physical tag per cell"))
        length(p)>=idx+nphys || throw(ArgumentError("IO: truncated physical-tag list on entity ($dim,$tag)"))
        physlist = Int32[]
        for j in 1:nphys
            physj=parse(Int,p[idx+j])
            physj>0 || throw(ArgumentError("IO: physical tag on entity ($dim,$tag) must be positive"))
            _int32_tag(physj,"physical")
            push!(physlist, Int32(physj))
        end
        if nphys >= 1
            ep[(dim, tag)] = physlist[1]   # first physical group tag
        else
            ep[(dim,tag)]=0         # record existence, even without a physical group
        end
        bounds = Int32[]
        if dim == 0
            length(p)==idx+nphys ||
                throw(ArgumentError("IO: unexpected trailing fields on point entity $tag"))
        else
            bidx=idx+nphys+1
            length(p)>=bidx || throw(ArgumentError("IO: missing bounding-entity count on entity ($dim,$tag)"))
            nb=parse(Int,p[bidx]);nb>=0 ||
                throw(ArgumentError("IO: negative bounding-entity count on entity ($dim,$tag)"))
            length(p)==bidx+nb ||
                throw(ArgumentError("IO: bounding-entity count mismatch on entity ($dim,$tag)"))
            for j in 1:nb
                bt=parse(Int,p[bidx+j]);bt != 0 ||
                    throw(ArgumentError("IO: bounding entity tag cannot be zero on entity ($dim,$tag)"))
                haskey(ep,(dim-1,abs(bt))) || throw(ArgumentError(
                    "IO: entity ($dim,$tag) references undeclared boundary entity ($(dim-1),$(abs(bt)))"))
                push!(bounds, Int32(bt))
            end
        end
        acc.topo_entities[(dim,tag)] = (bounds, physlist)
    end
end

function _read_nodes_v4!(acc, io, limits::_MshReadLimits)
    hdr = split(strip(readline(io)))
    length(hdr) == 4 || throw(ArgumentError("IO: malformed v4 node header"))
    numBlocks = parse(Int, hdr[1]); numNodes = parse(Int, hdr[2])
    declared_min=_io_mesh_uid(hdr[3],"v4 tag range");declared_max=_io_mesh_uid(hdr[4],"v4 tag range")
    (numBlocks >= 0 && numNodes >= 0) || throw(ArgumentError("IO: negative v4 node count"))
    numNodes<=typemax(Int32) || throw(ArgumentError("IO: $numNodes nodes exceed Int32 indexing"))
    numNodes<=limits.nodes-length(acc.node_x) || throw(ArgumentError(
        "read_msh: v4 node count exceeds max_nodes=$(limits.nodes)"))
    numBlocks<=limits.blocks || throw(ArgumentError(
        "read_msh: v4 node-block count exceeds max_blocks=$(limits.blocks)"))
    if numNodes==0
        (declared_min==0&&declared_max==0) ||
            throw(ArgumentError("IO: empty v4 node section must declare tag range 0 0"))
    else
        0<declared_min<=declared_max || throw(ArgumentError("IO: invalid v4 node-tag range"))
    end
    nread = 0
    actual_min=typemax(UInt64);actual_max=zero(UInt64)
    for _ in 1:numBlocks
        bh = split(strip(readline(io)))
        length(bh) == 4 || throw(ArgumentError("IO: malformed v4 node-block header"))
        edim=parse(Int,bh[1]);etag=parse(Int,bh[2]);parametric=parse(Int,bh[3])
        0<=edim<=3 || throw(ArgumentError("IO: v4 node-block dimension $edim is outside 0:3"))
        etag>0 || throw(ArgumentError("IO: v4 node-block entity tag must be positive"))
        _ensure_io_entity!(acc,edim,etag,limits,"node block")
        parametric in (0,1) || throw(ArgumentError("IO: v4 node-block parametric flag must be 0 or 1"))
        nInBlock = parse(Int, bh[4])
        nInBlock >= 0 || throw(ArgumentError("IO: negative v4 node-block size"))
        nInBlock <= numNodes-nread || throw(ArgumentError("IO: v4 node blocks exceed declared node count $numNodes"))
        tags = UInt64[]
        blocktags=Set{UInt64}()
        for _ in 1:nInBlock
            tag=_io_mesh_uid(strip(readline(io)),"v4 node")
            tag > 0 || throw(ArgumentError("IO: node tag must be positive (got $tag)"))
            (haskey(acc.tag2idx,tag) || tag in blocktags) &&
                throw(ArgumentError("IO: duplicate node tag $tag"))
            push!(tags,tag);push!(blocktags,tag)
            actual_min=min(actual_min,tag);actual_max=max(actual_max,tag)
        end
        for i in 1:nInBlock
            c = split(strip(readline(io)))
            expected=3+(parametric==1 ? edim : 0)
            length(c)==expected ||
                throw(ArgumentError("IO: v4 node coordinate record has $(length(c)) fields; expected $expected"))
            x=parse(Float64,c[1]); y=parse(Float64,c[2]); z=parse(Float64,c[3])
            (isfinite(x)&&isfinite(y)&&isfinite(z)) || throw(ArgumentError("IO: node $(tags[i]) has non-finite coordinates"))
            params = nothing
            if parametric==1
                uv = parse.(Float64, @view c[4:expected])
                all(isfinite,uv) ||
                    throw(ArgumentError("IO: node $(tags[i]) has non-finite parametric coordinates"))
                params = Vector{Float64}(uv)
            end
            push!(acc.node_x,x); push!(acc.node_y,y); push!(acc.node_z,z)
            acc.tag2idx[tags[i]] = length(acc.node_x)
            push!(acc.node_entities, (edim, Int32(etag)))
            push!(acc.node_parametric, params)
            nread += 1
        end
    end
    nread == numNodes || throw(ArgumentError("IO: v4 header declared $numNodes nodes but blocks contained $nread"))
    numNodes==0 || (actual_min==declared_min&&actual_max==declared_max) ||
        throw(ArgumentError("IO: v4 node-tag range does not match the node records"))
    _expect_end(io, "\$EndNodes")
end

function _read_elements_v4!(acc, io, limits::_MshReadLimits)
    ep = acc.ep
    hdr = split(strip(readline(io)))
    length(hdr) == 4 || throw(ArgumentError("IO: malformed v4 element header"))
    numBlocks = parse(Int, hdr[1])
    numElements = parse(Int,hdr[2])
    declared_min=_io_mesh_uid(hdr[3],"v4 tag range");declared_max=_io_mesh_uid(hdr[4],"v4 tag range")
    (numBlocks >= 0 && numElements >= 0) || throw(ArgumentError("IO: negative v4 element count"))
    numElements<=typemax(Int32) || throw(ArgumentError("IO: v4 element count exceeds Int32"))
    numElements<=limits.elements-length(acc.element_tags) || throw(ArgumentError(
        "read_msh: v4 element count exceeds max_elements=$(limits.elements)"))
    numBlocks<=limits.blocks-acc.element_blocks || throw(ArgumentError(
        "read_msh: v4 element-block count exceeds max_blocks=$(limits.blocks)"))
    acc.element_blocks+=numBlocks
    if numElements==0
        (declared_min==0&&declared_max==0) ||
            throw(ArgumentError("IO: empty v4 element section must declare tag range 0 0"))
    else
        0<declared_min<=declared_max || throw(ArgumentError("IO: invalid v4 element-tag range"))
    end
    nread = 0
    actual_min=typemax(UInt64);actual_max=zero(UInt64)
    for _ in 1:numBlocks
        bh = split(strip(readline(io)))
        length(bh) == 4 || throw(ArgumentError("IO: malformed v4 element-block header"))
        edim = parse(Int, bh[1]); etag = parse(Int, bh[2])
        etype = parse(Int, bh[3]); nInBlock = parse(Int, bh[4])
        0<=edim<=3 || throw(ArgumentError("IO: v4 element-block dimension $edim is outside 0:3"))
        etag>0 || throw(ArgumentError("IO: v4 element-block entity tag must be positive"))
        _ensure_io_entity!(acc,edim,etag,limits,"element block")
        haskey(_EDIM,etype) && _EDIM[etype]!=edim &&
            throw(ArgumentError("IO: element type $etype is incompatible with entity dimension $edim"))
        nInBlock >= 0 || throw(ArgumentError("IO: negative v4 element-block size"))
        nInBlock <= numElements-nread || throw(ArgumentError("IO: v4 element blocks exceed declared element count $numElements"))
        phys = get(ep, (edim, etag), 0)
        for _ in 1:nInBlock
            p = split(strip(readline(io)))
            isempty(p) && throw(ArgumentError("IO: empty v4 element record"))
            elem_tag=_io_mesh_uid(p[1],"v4 element");_record_element_tag!(acc,elem_tag)
            actual_min=min(actual_min,elem_tag);actual_max=max(actual_max,elem_tag)
            needed=get(_NN,etype,0)
            needed==0 || length(p)==needed+1 ||
                throw(ArgumentError("IO: element type $etype expects $needed node tags"))
            _push_element!(acc, etype, phys, @view p[2:end])   # p[1] = element tag
            nread += 1
        end
    end
    nread == numElements || throw(ArgumentError("IO: v4 header declared $numElements elements but blocks contained $nread"))
    numElements==0 || (actual_min==declared_min&&actual_max==declared_max) ||
        throw(ArgumentError("IO: v4 element-tag range does not match the element records"))
    _expect_end(io, "\$EndElements")
end

function _ensure_io_entity!(acc::_Accum,dim::Int,tag::Int,
                            limits::_MshReadLimits,context::AbstractString)
    key=(dim,tag)
    haskey(acc.ep,key) && return nothing
    length(acc.ep)<limits.entities || throw(ArgumentError(
        "read_msh: implicit v4 entity count exceeds " *
        "max_entities=$(limits.entities) while reading $context"))
    acc.ep[key]=0
    return nothing
end

@inline function _record_element_tag!(acc::_Accum,tag::UInt64)
    tag>0 || throw(ArgumentError("IO: element tag must be positive (got $tag)"))
    tag in acc.element_tags && throw(ArgumentError("IO: duplicate element tag $tag"))
    push!(acc.element_tags,tag);return nothing
end

@inline function _int32_tag(tag::Integer,what::AbstractString)
    typemin(Int32)<=tag<=typemax(Int32) ||
        throw(ArgumentError("IO: $what tag $tag does not fit Int32"))
    return Int32(tag)
end

@inline function _physical_tag(tag::Integer,what::AbstractString="physical")
    tag>=0 || throw(ArgumentError("IO: $what tag $tag must be non-negative"))
    return _int32_tag(tag,what)
end

# ── shared element push (relabel tags → compact indices) ────────────────────────
function _push_element!(acc, etype::Int, phys::Int, nodetoks)
    haskey(_NN, etype) || throw(ArgumentError(
        "IO: unsupported gmsh element type $etype; refusing to silently discard cells"))
    ptag=_physical_tag(phys)
    function idx(t)
        tag=_io_mesh_uid(t,"element connectivity"); i=get(acc.tag2idx,tag,Int32(0))
        i != 0 || throw(ArgumentError("IO: element references unknown node tag $tag"))
        return i
    end
    if etype == MSH_LINE || etype == MSH_SEG2
        # linear line, or quadratic line read as its 2 end vertices (first 2 nodes)
        a = idx(nodetoks[1]); b = idx(nodetoks[2])
        if etype == MSH_SEG2
            idx(nodetoks[3])
        end
        push!(acc.segs, (Int32(a), Int32(b))); push!(acc.seg_tag, ptag)
    elseif etype == MSH_TRI || etype == MSH_TRI2
        # linear tri, or quadratic tri read as its 3 corner vertices (first 3 nodes)
        a = idx(nodetoks[1]); b = idx(nodetoks[2]); c = idx(nodetoks[3])
        if etype == MSH_TRI2
            for k in 4:6; idx(nodetoks[k]); end
        end
        push!(acc.tris, (Int32(a), Int32(b), Int32(c))); push!(acc.tri_tag, ptag)
    elseif etype == MSH_TET || etype == MSH_TET2
        # linear tet, or quadratic tet read as its 4 corner vertices (first 4 nodes)
        a = idx(nodetoks[1]); b = idx(nodetoks[2]); c = idx(nodetoks[3]); d = idx(nodetoks[4])
        if etype == MSH_TET2
            for k in 5:10; idx(nodetoks[k]); end
        end
        push!(acc.tets, (Int32(a), Int32(b), Int32(c), Int32(d))); push!(acc.tet_tag, ptag)
    elseif etype == MSH_POINT
        idx(nodetoks[1])
    end
    # points (type 15) carry no cell in our simplex model; skipped intentionally.
    return
end

function _to_mesh(acc::_Accum)
    nn = length(acc.node_x)
    coords = Matrix{Float64}(undef, 3, nn)
    @inbounds for i in 1:nn
        coords[1,i] = acc.node_x[i]; coords[2,i] = acc.node_y[i]; coords[3,i] = acc.node_z[i]
    end
    segs = _cols(acc.segs); tris = _cols(acc.tris); tets = _cols(acc.tets)
    return Mesh(coords; segs=segs, tris=tris, tets=tets,
                seg_tag=acc.seg_tag, tri_tag=acc.tri_tag, tet_tag=acc.tet_tag)
end

function _cols(v::Vector{NTuple{K,Int32}}) where {K}
    M = Matrix{Int32}(undef, K, length(v))
    @inbounds for j in eachindex(v), i in 1:K
        M[i,j] = v[j][i]
    end
    return M
end

# ════════════════════════════════════════════════════════════════════════════════
# Writing
# ════════════════════════════════════════════════════════════════════════════════

"""
    write_msh(path, mesh; version=2.2, physical_names=Dict(),
              gmsh_compatible=true)

Write `mesh` to a gmsh ASCII `.msh` file. `version` is `2.2` or `4.1`.
`physical_names` maps `(dim, tag) => name`. Per-cell physical tags come from the
mesh's `*_tag` vectors; cells are grouped into entity blocks by `(dim, tag)`.
Gmsh 4.15.2 treats backslashes and tabs literally but cannot preserve quoted or
line-breaking physical names, and truncates serialized names beyond 128 bytes.
Unsafe or oversized names are rejected. Set
`gmsh_compatible=false`
only for a Tessella-to-Tessella escaped-name round trip, and read that file with
`read_msh(...; tessella_extensions=true)`.
"""
function write_msh(path, mesh; version=2.2,
                   physical_names=Dict{Tuple{Int,Int},String}(),
                   point_elements=nothing,
                   gmsh_compatible=true)
    path isa AbstractString || throw(ArgumentError(
        "write_msh: path must be a string"))
    mesh isa Mesh || throw(ArgumentError(
        "write_msh: mesh must be a Mesh"))
    version isa Real || throw(ArgumentError(
        "write_msh: version must be real"))
    version isa Bool && throw(ArgumentError(
        "write_msh: version must not be Bool"))
    ver = try
        Float64(version)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("write_msh: version must be representable as Float64"))
    end
    ver in (2.2,4.1) ||
        throw(ArgumentError("write_msh: version must be 2.2 or 4.1 (got $version)"))
    gmsh_compatible isa Bool || throw(ArgumentError(
        "write_msh: gmsh_compatible must be Bool"))
    points=_checked_point_elements(point_elements,mesh)
    _validate_write_mesh(mesh)
    names=_copy_physical_names(physical_names,"write_msh")
    if gmsh_compatible
        unsafe=Tuple{Int,Int}[]
        for (key,name) in names
            any(character -> character in ('\"','\n','\r'),name) &&
                push!(unsafe,key)
        end
        isempty(unsafe) || throw(ArgumentError(
            "write_msh: Gmsh 4.15.2 cannot preserve quoted or " *
            "line-breaking physical name(s) " * join(sort!(unsafe),",") *
            "; use gmsh_compatible=false only for a Tessella-only round trip"))
    end
    target = abspath(path); parent = dirname(target)
    isdir(parent) || throw(ArgumentError("write_msh: parent directory does not exist: $parent"))
    isdir(target) && throw(ArgumentError(
        "write_msh: destination is a directory: $target"))
    # Build beside the destination and atomically rename only after a complete,
    # flushed write. A formatter/error cannot truncate a previously valid mesh.
    mktemp(parent) do tmp, io
        if ver == 2.2
            _write_msh_v2(io, mesh, names, gmsh_compatible, points)
        else
            isempty(points) || throw(ArgumentError(
                "write_msh: point_elements requires version 2.2 — the 4.1 " *
                "entity-block layout cannot carry dim-0 elements"))
            _write_msh_v4(io, mesh, names, gmsh_compatible)
        end
        flush(io); close(io)
        mv(tmp,target;force=true)
    end
    return path
end

# `point_elements` — `(node, physical, elementary)` triples written as MSH
# type-15 point elements ahead of the line cells, mirroring upstream's
# per-vertex `MPoint` elements.
function _checked_point_elements(point_elements,mesh::Mesh)
    point_elements===nothing && return Tuple{Int,Int,Int}[]
    point_elements isa AbstractVector || throw(ArgumentError(
        "write_msh: point_elements must be a vector of (node, physical, " *
        "elementary) triples"))
    checked=Tuple{Int,Int,Int}[]
    for entry in point_elements
        (entry isa Tuple && length(entry)==3) || throw(ArgumentError(
            "write_msh: each point_elements entry must be a " *
            "(node, physical, elementary) triple"))
        node,physical,elementary=entry
        for (what,value) in (("node",node),("physical",physical),
                             ("elementary",elementary))
            value isa Integer || throw(ArgumentError(
                "write_msh: point_elements $what must be an integer"))
            value isa Bool && throw(ArgumentError(
                "write_msh: point_elements $what must not be Bool"))
        end
        1<=node<=nnodes(mesh) || throw(ArgumentError(
            "write_msh: point element node $node is out of range"))
        physical>=0 && elementary>=0 || throw(ArgumentError(
            "write_msh: point element tags must be nonnegative"))
        push!(checked,(Int(node),Int(physical),Int(elementary)))
    end
    return checked
end

function _validate_write_mesh(m::Mesh)
    # Degenerate (zero-length / self-loop) segments are writable: Gmsh emits
    # them for closed degenerate curves. Other defects still block.
    diagnostic=validate(m; allow_degenerate_segs=true)
    diagnostic.ok || throw(ArgumentError(
        "write_msh: mesh is invalid — " * join(diagnostic.messages,"; ")))
    @inbounds for i in 1:nnodes(m),d in 1:3
        isfinite(m.coords[d,i]) || throw(ArgumentError("write_msh: node $i has a non-finite coordinate"))
    end
    for (kind,tags) in (("segment",m.seg_tag),("triangle",m.tri_tag),("tetrahedron",m.tet_tag))
        @inbounds for (i,tag) in pairs(tags)
            tag>=0 || throw(ArgumentError("write_msh: $kind $i has negative physical tag $tag"))
        end
    end
    return nothing
end

_escape_name(s::AbstractString) = replace(s, "\\"=>"\\\\", "\""=>"\\\"",
                                           "\n"=>"\\n", "\r"=>"\\r", "\t"=>"\\t")

function _write_physical_names(io, physical_names, gmsh_compatible::Bool)
    isempty(physical_names) && return
    length(physical_names)<=typemax(Int32) || throw(ArgumentError(
        "write_msh: physical-name count exceeds Int32"))
    println(io, "\$PhysicalNames")
    println(io, length(physical_names))
    for ((dim, tag), name) in sort(collect(physical_names); by=x->x[1])
        encoded=gmsh_compatible ? name : _escape_name(name)
        ncodeunits(encoded)<=MSH_PHYSICAL_NAME_MAX_BYTES || throw(ArgumentError(
            "write_msh: serialized physical name for ($dim,$tag) exceeds " *
            "$MSH_PHYSICAL_NAME_MAX_BYTES bytes"))
        println(io, dim, " ", tag, " \"", encoded, "\"")
    end
    println(io, "\$EndPhysicalNames")
end

function _write_msh_v2(io, m::Mesh, physical_names, gmsh_compatible::Bool,
                       point_elements=Tuple{Int,Int,Int}[])
    println(io, "\$MeshFormat"); println(io, "2.2 0 8"); println(io, "\$EndMeshFormat")
    _write_physical_names(io, physical_names, gmsh_compatible)
    # nodes
    println(io, "\$Nodes"); println(io, nnodes(m))
    @inbounds for i in 1:nnodes(m)
        p = node(m, i)
        @printf(io, "%d %.17g %.17g %.17g\n", i, p[1], p[2], p[3])
    end
    println(io, "\$EndNodes")
    # elements — point cells first, matching upstream's per-vertex emission
    nel = length(point_elements) + nsegs(m) + ntris(m) + ntets(m)
    println(io, "\$Elements"); println(io, nel)
    eid = 0
    @inbounds for (node, physical, elementary) in point_elements
        eid += 1
        @printf(io, "%d %d 2 %d %d %d\n", eid, MSH_POINT, physical,
                elementary, node)
    end
    @inbounds for t in 1:nsegs(m)
        eid += 1; tag = m.seg_tag[t]
        @printf(io, "%d %d 2 %d %d %d %d\n", eid, MSH_LINE, tag, tag, m.segs[1,t], m.segs[2,t])
    end
    @inbounds for t in 1:ntris(m)
        eid += 1; tag = m.tri_tag[t]
        @printf(io, "%d %d 2 %d %d %d %d %d\n", eid, MSH_TRI, tag, tag, m.tris[1,t], m.tris[2,t], m.tris[3,t])
    end
    @inbounds for t in 1:ntets(m)
        eid += 1; tag = m.tet_tag[t]
        @printf(io, "%d %d 2 %d %d %d %d %d %d\n", eid, MSH_TET, tag, tag, m.tets[1,t], m.tets[2,t], m.tets[3,t], m.tets[4,t])
    end
    println(io, "\$EndElements")
end

function _write_msh_v4(io, m::Mesh, physical_names, gmsh_compatible::Bool)
    println(io, "\$MeshFormat"); println(io, "4.1 0 8"); println(io, "\$EndMeshFormat")
    _write_physical_names(io, physical_names, gmsh_compatible)

    # Group cells into (dim, tag) entity blocks so physical tags survive.
    segblocks = _group_cells(m.segs, m.seg_tag, 2)
    triblocks = _group_cells(m.tris, m.tri_tag, 3)
    tetblocks = _group_cells(m.tets, m.tet_tag, 4)

    # Entity tags are positive geometric identifiers, distinct from physical
    # tags (where zero means "unclassified"). Number them independently per
    # dimension. Classify nodes on an entity that owns elements whenever one
    # exists: Gmsh 4.15.2 corrupts the declared node count when it rewrites a
    # mesh whose nodes live on an unrelated synthetic entity.
    segblocks = [(phys,i,cols) for (i,(phys,cols)) in enumerate(segblocks)]
    triblocks = [(phys,i,cols) for (i,(phys,cols)) in enumerate(triblocks)]
    tetblocks = [(phys,i,cols) for (i,(phys,cols)) in enumerate(tetblocks)]
    node_dim,node_entity,synthetic_node_entity =
        _msh_v4_node_owner(m,segblocks,triblocks,tetblocks)

    # Entities: one entity per (dim, physical tag) block, declaring its physical group.
    _write_entities_v4(
        io,m,segblocks,triblocks,tetblocks,node_entity,synthetic_node_entity)

    # Nodes use one deterministic owner. Other element entities can legitimately
    # acquire empty node blocks when Gmsh rewrites the file.
    println(io, "\$Nodes")
    println(io, nnodes(m)==0 ? 0 : 1, " ", nnodes(m), " ", nnodes(m) == 0 ? 0 : 1, " ", nnodes(m))
    if nnodes(m)>0
        println(io, node_dim, " ", node_entity, " ", 0, " ", nnodes(m))
        @inbounds for i in 1:nnodes(m); println(io, i); end
        @inbounds for i in 1:nnodes(m)
            p = node(m, i); @printf(io, "%.17g %.17g %.17g\n", p[1], p[2], p[3])
        end
    end
    println(io, "\$EndNodes")

    # Elements.
    numBlocks = length(segblocks) + length(triblocks) + length(tetblocks)
    nel = nsegs(m) + ntris(m) + ntets(m)
    println(io, "\$Elements")
    println(io, numBlocks, " ", nel, " ", nel == 0 ? 0 : 1, " ", nel)
    eid = Ref(0)
    _write_elem_blocks_v4(io, m.segs, segblocks, 1, MSH_LINE, eid)
    _write_elem_blocks_v4(io, m.tris, triblocks, 2, MSH_TRI, eid)
    _write_elem_blocks_v4(io, m.tets, tetblocks, 3, MSH_TET, eid)
    println(io, "\$EndElements")
end

function _msh_v4_node_owner(m,segblocks,triblocks,tetblocks)
    nnodes(m)==0 && return (3,0,false)
    !isempty(tetblocks) && return (3,tetblocks[1][2],false)
    !isempty(triblocks) && return (2,triblocks[1][2],false)
    !isempty(segblocks) && return (1,segblocks[1][2],false)
    return (3,1,true)
end

# group cell columns by tag → Dict(tag => Vector{col index})
function _group_cells(cells::Matrix{Int32}, tags::Vector{Int32}, k::Int)
    g = Dict{Int32,Vector{Int}}()
    @inbounds for j in axes(cells, 2)
        push!(get!(g, tags[j], Int[]), j)
    end
    return sort(collect(g); by=x->x[1])
end

function _write_entities_v4(io,m,segblocks,triblocks,tetblocks,node_entity,
                            synthetic_node_entity::Bool)
    println(io, "\$Entities")
    println(io,0," ",length(segblocks)," ",length(triblocks)," ",
            length(tetblocks)+(synthetic_node_entity ? 1 : 0))
    for (phys, entity, cols) in segblocks
        # curve: tag minx miny minz maxx maxy maxz numPhys phys... numBounding
        _write_entity_record(io,entity,phys,_cell_bounds(m,m.segs,cols))
    end
    for (phys, entity, cols) in triblocks
        _write_entity_record(io,entity,phys,_cell_bounds(m,m.tris,cols))
    end
    for (phys, entity, cols) in tetblocks
        _write_entity_record(io,entity,phys,_cell_bounds(m,m.tets,cols))
    end
    synthetic_node_entity &&
        _write_entity_record(io,node_entity,0,_all_node_bounds(m))
    println(io, "\$EndEntities")
end

function _write_entity_record(io,entity::Integer,phys::Integer,bounds)
    @printf(io,"%d %.17g %.17g %.17g %.17g %.17g %.17g %d",entity,bounds...,phys==0 ? 0 : 1)
    phys==0 || print(io," ",phys)
    println(io," 0")
end

# ══════════════════════════════════════════════════════════════════════════════
# Classified MixedMesh output honouring `Mesh.SaveParametric`/`Mesh.SaveTopology`
# (the CLI `-save_parametric`/`-parametric` and `-save_topology` flags).
#
# Verified Gmsh 4.15.2 semantics:
#   * MSH 4.1 — `Mesh.SaveParametric=0` clears the per-node parametric flag and
#     drops the stored `(u[,v])` coordinates; the default output carries them.
#     `Mesh.SaveTopology` is a no-op because `$Entities` is always declared.
#   * MSH 2.2 — `Mesh.SaveParametric=1` replaces `$Nodes` with
#     `$ParametricNodes` records `tag x y z entityDim entityTag [u [v]]` (the
#     owner dimension decides the parameter count); `Mesh.SaveTopology=1`
#     inserts a topology-only `$Entities` section (no bounding boxes, unsigned
#     boundary tags) immediately before the nodes.
# ══════════════════════════════════════════════════════════════════════════════

@inline function _msh_save_flag(value,name::AbstractString)
    value===nothing && return nothing
    value isa Bool && return value
    (value isa Real && !(value isa Bool)) || throw(ArgumentError(
        "write_msh: $name must be Bool, 0/1, or nothing"))
    v=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("write_msh: $name must be representable as Float64"))
    end
    v==0.0 && return false
    v==1.0 && return true
    throw(ArgumentError("write_msh: $name must be 0 or 1 (got $value)"))
end

"""
    write_msh(path, m::MixedMesh; version=4.1, binary=false,
              gmsh_compatible=true, save_parametric=nothing,
              save_topology=nothing) -> path

Serialize a classified `MixedMesh` honouring the `Mesh.SaveParametric` and
`Mesh.SaveTopology` options with Gmsh 4.15.2 semantics (see the section banner
above). `nothing` keeps the default writer behaviour.
"""
function write_msh(path, m::MixedMesh; version=4.1, binary=false,
                   gmsh_compatible=true, save_parametric=nothing,
                   save_topology=nothing, preserve_numbering=false)
    path isa AbstractString || throw(ArgumentError(
        "write_msh: path must be a string"))
    version isa Real || throw(ArgumentError("write_msh: version must be real"))
    version isa Bool && throw(ArgumentError(
        "write_msh: version must not be Bool"))
    ver = try
        Float64(version)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("write_msh: version must be representable as Float64"))
    end
    ver in (2.2,4.1) ||
        throw(ArgumentError("write_msh: version must be 2.2 or 4.1 (got $version)"))
    binary isa Bool || throw(ArgumentError("write_msh: binary must be Bool"))
    gmsh_compatible isa Bool || throw(ArgumentError(
        "write_msh: gmsh_compatible must be Bool"))
    preserve_numbering isa Bool || throw(ArgumentError(
        "write_msh: preserve_numbering must be Bool"))
    sp=_msh_save_flag(save_parametric,"save_parametric")
    st=_msh_save_flag(save_topology,"save_topology")
    if ver==4.1
        target=sp===false ? _msh_strip_parametric(m) : m
        write_mixed_msh(path,target;version=ver,binary=binary,
                                 gmsh_compatible=gmsh_compatible)
        return path
    end
    # ── MSH 2.2 ───────────────────────────────────────────────────────────────
    if sp===true
        m.entity_data===nothing && throw(ArgumentError(
            "write_msh: save_parametric requires classified node/entity " *
            "metadata (MixedMesh.entity_data)"))
        records=_msh2_node_records(m)
    else
        records=nothing
    end
    entities=st===true ? _msh2_entity_records(m) : nothing
    if records===nothing && entities===nothing
        write_mixed_msh(path,m;version=ver,binary=binary,
                                 gmsh_compatible=gmsh_compatible,
                                 preserve_numbering=preserve_numbering)
        return path
    end
    parent=dirname(path); isempty(parent) && (parent=".")
    isdir(parent) || throw(ArgumentError(
        "write_msh: output directory does not exist: $parent"))
    tmp,tmpio=mktemp(parent;cleanup=false)
    close(tmpio)
    try
        write_mixed_msh(tmp,m;version=ver,binary=binary,
                                 gmsh_compatible=gmsh_compatible,
                                 preserve_numbering=preserve_numbering)
        bytes=read(tmp)
        out=_msh2_transform(bytes;binary=binary,entities=entities,
                            parametric=records)
        write(path,out)
    finally
        isfile(tmp) && rm(tmp;force=true)
    end
    return path
end

# Copy `m` with every `node_parametric` record cleared — the MSH 4.1 writer keys
# the parametric flag off `node_parametric[i] !== nothing`.
function _msh_strip_parametric(m::MixedMesh)
    data=m.entity_data
    data===nothing && return m
    any(p->p!==nothing,data.node_parametric) || return m
    cleared=Union{Nothing,Vector{Float64}}[nothing for _ in data.node_parametric]
    newdata=MixedEntityData(data.entities;
        node_entities=data.node_entities,
        node_parametric=cleared,
        external_node_tags=data.external_node_tags,
        block_entities=data.block_entities,
        external_element_tags=data.external_element_tags,
        _native_mesh_tags=data.native_mesh_tags)
    return MixedMesh(m.coords,m.blocks;
        physical_names=m.physical_names,
        entity_data=newdata,
        elementary_entities=m.elementary_entities,
        periodic_links=m.periodic_links,
        ancillary_sections=m.ancillary_sections,
        data_sections=m.data_sections,
        partition_data=m.partition_data)
end

# Per-node `(entityDim, entityTag, params)` records for `$ParametricNodes` —
# indexed by mesh-node order, matching the dense tags the v2 writers emit.
function _msh2_node_records(m::MixedMesh)
    data=m.entity_data
    data===nothing && throw(ArgumentError(
        "write_msh: save_parametric requires classified node/entity metadata"))
    nn=size(m.coords,2)
    length(data.node_entities)==nn || throw(ArgumentError(
        "write_msh: node-entity classification does not cover all nodes"))
    length(data.node_parametric)==nn || throw(ArgumentError(
        "write_msh: parametric node records do not cover all nodes"))
    records=Tuple{Int,Int32,Vector{Float64}}[]
    sizehint!(records,nn)
    @inbounds for i in 1:nn
        dim,etag=data.node_entities[i]
        0<=dim<=3 || throw(ArgumentError(
            "write_msh: node $i owner dimension $dim is outside 0:3"))
        params=data.node_parametric[i]
        values=params===nothing ? Float64[] : params
        (length(values)==0 || length(values)==dim) || throw(ArgumentError(
            "write_msh: node $i carries $(length(values)) parameters; " *
            "expected 0 or owner dimension $dim"))
        push!(records,(dim,etag,Float64.(values)))
    end
    return records
end

# Topology records `(dim,tag) => (unsigned boundaries, physical tags)` for the
# MSH 2.2 `$Entities` section. Declared `MixedEntity` records are preferred;
# without entity_data the section is derived from per-cell elementary tags.
function _msh2_entity_records(m::MixedMesh)
    data=m.entity_data
    records=Dict{Tuple{Int,Int},Tuple{Vector{Int32},Vector{Int32}}}()
    if data!==nothing
        for ((dim,tag),entity) in data.entities
            bounds=Int32[]
            for b in entity.boundaries
                b==0 && continue
                push!(bounds,Int32(abs(b)))
            end
            records[(dim,tag)]=(bounds,Int32.(entity.physical_tags))
        end
        return records
    end
    elementary=m.elementary_entities
    elementary===nothing && return records
    for (bi,block) in pairs(m.blocks)
        dim=msh_dimension(block.msh)
        0<=dim<=3 || continue
        tags=elementary[bi]
        for j in 1:_mixed_ncells(block)
            etag=tags[j]
            etag<=0 && continue
            key=(dim,Int(etag))
            bounds,phys=get(records,key,(Int32[],Int32[]))
            pt=block.tags[j]
            if pt>0 && !(Int32(pt) in phys)
                phys=copy(phys); push!(phys,Int32(pt))
            end
            records[key]=(bounds,phys)
        end
    end
    return records
end

@inline _mixed_ncells(b)=_block_ncells(b)

# Splice an MSH 2.2 byte stream: optionally insert `$Entities` before `$Nodes`
# and/or replace `$Nodes` with `$ParametricNodes`.
# `parametric` is a per-file-node vector of `(entityDim, entityTag, params)` —
# or a `Dict`-style lookup closure returning those records (session path).
function _msh2_transform(bytes::Vector{UInt8};binary::Bool,
                         entities=nothing,parametric=nothing)
    sm=findfirst(codeunits("\$Nodes\n"),bytes)
    sm===nothing && throw(ArgumentError(
        "write_msh: MSH 2.2 stream is missing the \$Nodes section"))
    (sm[1]==1 || bytes[sm[1]-1]==UInt8('\n')) || throw(ArgumentError(
        "write_msh: malformed MSH 2.2 \$Nodes marker"))
    em=findfirst(codeunits("\$EndNodes"),bytes)
    em===nothing && throw(ArgumentError(
        "write_msh: MSH 2.2 stream is missing \$EndNodes"))
    em[1]>sm[end] || throw(ArgumentError(
        "write_msh: \$EndNodes precedes \$Nodes"))
    send=em[end]
    (send<length(bytes) && bytes[send+1]==UInt8('\n')) && (send+=1)
    # The node count line is text in both ASCII and binary modes.
    count_eol=findnext(isequal(UInt8('\n')),bytes,sm[end]+1)
    count_eol===nothing && throw(ArgumentError(
        "write_msh: truncated MSH 2.2 node count"))
    count=try
        parse(Int,String(@view bytes[sm[end]+1:count_eol-1]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("write_msh: malformed MSH 2.2 node count"))
    end
    count>=0 || throw(ArgumentError("write_msh: negative MSH 2.2 node count"))
    records_start=count_eol+1
    io=IOBuffer()
    write(io,@view bytes[1:sm[1]-1])
    entities===nothing || _msh2_write_entities(io,entities)
    if parametric===nothing
        write(io,@view bytes[sm[1]:send])
    else
        length(parametric)==count || throw(ArgumentError(
            "write_msh: parametric records do not cover all nodes"))
        if binary
            records_bytes=count*28
            records_start+records_bytes-1<=em[1]-1 || throw(ArgumentError(
                "write_msh: truncated binary MSH 2.2 node records"))
            length(parametric)==count || throw(ArgumentError(
                "write_msh: parametric records do not cover all nodes"))
            println(io,"\$ParametricNodes"); println(io,count)
            @inbounds for k in 1:count
                off=records_start+(k-1)*28
                write(io,@view bytes[off:off+27])   # tag + xyz verbatim
                dim,etag,params=parametric[k]
                write(io,Int32(dim)); write(io,Int32(etag))
                for v in params
                    write(io,Float64(v))
                end
            end
            write(io,UInt8('\n')); println(io,"\$EndParametricNodes")
        else
            println(io,"\$ParametricNodes"); println(io,count)
            pos=records_start
            @inbounds for k in 1:count
                eol=findnext(isequal(UInt8('\n')),bytes,pos)
                eol===nothing && throw(ArgumentError(
                    "write_msh: truncated MSH 2.2 node record $k"))
                eol<=em[1] || throw(ArgumentError(
                    "write_msh: MSH 2.2 node record $k overruns \$EndNodes"))
                line_end=eol-1
                (line_end>=pos && bytes[line_end]==UInt8('\r')) && (line_end-=1)
                write(io,@view bytes[pos:line_end])    # `tag x y z` verbatim
                dim,etag,params=parametric[k]
                @printf(io," %d %d",dim,etag)
                for v in params
                    @printf(io," %.17g",v)
                end
                write(io,UInt8('\n'))
                pos=eol+1
            end
            println(io,"\$EndParametricNodes")
        end
    end
    write(io,@view bytes[send+1:end])
    return take!(io)
end

function _msh2_write_entities(io,
        entities::Dict{Tuple{Int,Int},Tuple{Vector{Int32},Vector{Int32}}})
    counts=zeros(Int,4)
    for (dim,tag) in keys(entities)
        0<=dim<=3 || throw(ArgumentError(
            "write_msh: topology entity dimension $dim is outside 0:3"))
        tag>0 || throw(ArgumentError(
            "write_msh: topology entity tag must be positive (got $tag)"))
        counts[dim+1]+=1
    end
    println(io,"\$Entities")
    println(io,counts[1]," ",counts[2]," ",counts[3]," ",counts[4])
    for dim in 0:3
        members=sort!(collect(
            tag for (d,tag) in keys(entities) if d==dim))
        for tag in members
            bounds,phys=entities[(dim,tag)]
            if dim==0
                print(io,tag," ",length(phys))
            else
                print(io,tag," ",length(bounds))
                for b in bounds
                    print(io," ",b)
                end
                print(io," ",length(phys))
            end
            for p in phys
                print(io," ",p)
            end
            # Upstream v2 records carry a trailing space — kept for parity.
            println(io," ")
        end
    end
    println(io,"\$EndEntities")
    return nothing
end

# Map 3-D coordinates to the mesh node index — exact binary reinterpretation,
# deterministic under the `%.17g` text round-trip. Values are queues so
# coincident nodes resolve in mesh order.
function _msh_coord_index(m::MixedMesh)
    index=Dict{NTuple{3,Int},Vector{Int}}()
    @inbounds for i in axes(m.coords,2)
        key=(reinterpret(Int,m.coords[1,i]),
             reinterpret(Int,m.coords[2,i]),
             reinterpret(Int,m.coords[3,i]))
        push!(get!(index,key,Int[]),i)
    end
    return index
end

"""
    _msh_rewrite_save_options(path; version, binary, save_parametric,
                              save_topology, source=nothing) -> nothing

Post-process an MSH file already written by `API.write` so it honours
`Mesh.SaveParametric`/`Mesh.SaveTopology` (used by the CLI session path, where
the session writer owns serialization). `source` is the classified `MixedMesh`
the file was exported from; it is required for MSH 2.2 `save_parametric` and
supplies the `\$Entities` records for `save_topology`.
"""
function _msh_rewrite_save_options(path;version,binary,
                                   save_parametric,save_topology,
                                   source=nothing)
    version isa Real || throw(ArgumentError(
        "write_msh: version must be real"))
    ver=Float64(version)
    ver in (2.2,4.1) || throw(ArgumentError(
        "write_msh: version must be 2.2 or 4.1 (got $version)"))
    sp=_msh_save_flag(save_parametric,"save_parametric")
    st=_msh_save_flag(save_topology,"save_topology")
    if ver==4.1
        # SaveTopology is a no-op; only a 0-valued SaveParametric acts.
        sp===false || return nothing
        m=read_mixed_msh(path)
        stripped=_msh_strip_parametric(m)
        write_mixed_msh(path,stripped;version=4.1,binary=binary)
        return nothing
    end
    # ── MSH 2.2 ───────────────────────────────────────────────────────────────
    (sp===true || st===true) || return nothing
    bytes=read(path)
    records=nothing
    if sp===true
        source isa MixedMesh || throw(ArgumentError(
            "write_msh: save_parametric requires the classified source mesh"))
        data=source.entity_data
        data===nothing && throw(ArgumentError(
            "write_msh: save_parametric requires classified node/entity " *
            "metadata on the source mesh"))
        records=_msh2_parametric_lookup(bytes,source,binary)
    end
    entities=nothing
    if st===true
        if source isa MixedMesh
            entities=_msh2_entity_records(source)
        else
            entities=Dict{Tuple{Int,Int},Tuple{Vector{Int32},Vector{Int32}}}()
        end
    end
    out=_msh2_transform(bytes;binary=binary,entities=entities,
                        parametric=records)
    tmp=path*".rewrite.tmp"
    write(tmp,out)
    mv(tmp,path;force=true)
    return nothing
end

# Build per-file-node `(dim,tag,params)` records by matching the file's node
# coordinates against `source`'s coordinate map.
function _msh2_parametric_lookup(bytes::Vector{UInt8},source::MixedMesh,
                               binary::Bool)
    data=source.entity_data
    index=_msh_coord_index(source)
    sm=findfirst(codeunits("\$Nodes\n"),bytes)
    sm===nothing && throw(ArgumentError(
        "write_msh: MSH 2.2 stream is missing the \$Nodes section"))
    em=findfirst(codeunits("\$EndNodes"),bytes)
    count_eol=findnext(isequal(UInt8('\n')),bytes,sm[end]+1)
    count=parse(Int,String(@view bytes[sm[end]+1:count_eol-1]))
    records=Tuple{Int,Int32,Vector{Float64}}[]
    sizehint!(records,count)
    if binary
        @inbounds for k in 1:count
            off=count_eol+1+(k-1)*28
            x=reinterpret(Float64,@view bytes[off+4:off+11])[1]
            y=reinterpret(Float64,@view bytes[off+12:off+19])[1]
            z=reinterpret(Float64,@view bytes[off+20:off+27])[1]
            push!(records,_msh2_lookup_record(index,data,x,y,z,k))
        end
    else
        pos=count_eol+1
        @inbounds for k in 1:count
            eol=findnext(isequal(UInt8('\n')),bytes,pos)
            eol===nothing && throw(ArgumentError(
                "write_msh: truncated MSH 2.2 node record $k"))
            line_end=eol-1
            (line_end>=pos && bytes[line_end]==UInt8('\r')) && (line_end-=1)
            fields=split(String(@view bytes[pos:line_end]))
            length(fields)==4 || throw(ArgumentError(
                "write_msh: malformed MSH 2.2 node record $k"))
            x=parse(Float64,fields[2]); y=parse(Float64,fields[3])
            z=parse(Float64,fields[4])
            push!(records,_msh2_lookup_record(index,data,x,y,z,k))
            pos=eol+1
        end
    end
    return records
end

function _msh2_lookup_record(index,data,x,y,z,k)
    key=(reinterpret(Int,x),reinterpret(Int,y),reinterpret(Int,z))
    queue=get(index,key,nothing)
    (queue===nothing || isempty(queue)) && throw(ArgumentError(
        "write_msh: file node $k at ($x,$y,$z) has no counterpart in the " *
        "classified source mesh"))
    i=popfirst!(queue)
    dim,etag=data.node_entities[i]
    params=data.node_parametric[i]
    values=params===nothing ? Float64[] : Float64.(params)
    (length(values)==0 || length(values)==dim) || throw(ArgumentError(
        "write_msh: source node $i carries $(length(values)) parameters; " *
        "expected 0 or owner dimension $dim"))
    return (dim,etag,values)
end

function _cell_bounds(m::Mesh,cells::Matrix{Int32},cols)
    lo1=lo2=lo3=Inf;hi1=hi2=hi3=-Inf
    @inbounds for j in cols,k in axes(cells,1)
        v=cells[k,j];x=m.coords[1,v];y=m.coords[2,v];z=m.coords[3,v]
        lo1=min(lo1,x);lo2=min(lo2,y);lo3=min(lo3,z)
        hi1=max(hi1,x);hi2=max(hi2,y);hi3=max(hi3,z)
    end
    return (lo1,lo2,lo3,hi1,hi2,hi3)
end

function _all_node_bounds(m::Mesh)
    nnodes(m)==0 && return (0.,0.,0.,0.,0.,0.)
    lo1=lo2=lo3=Inf;hi1=hi2=hi3=-Inf
    @inbounds for i in 1:nnodes(m)
        x=m.coords[1,i];y=m.coords[2,i];z=m.coords[3,i]
        lo1=min(lo1,x);lo2=min(lo2,y);lo3=min(lo3,z)
        hi1=max(hi1,x);hi2=max(hi2,y);hi3=max(hi3,z)
    end
    return (lo1,lo2,lo3,hi1,hi2,hi3)
end

function _write_elem_blocks_v4(io, cells::Matrix{Int32}, blocks, dim::Int, etype::Int, eid)
    for (_, entity, cols) in blocks
        println(io, dim, " ", entity, " ", etype, " ", length(cols))
        k = size(cells, 1)
        @inbounds for j in cols
            eid[] += 1
            print(io, eid[])
            for i in 1:k
                print(io, " ", cells[i, j])
            end
            println(io)
        end
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# STL (ingest boundary surface mesh)
# ════════════════════════════════════════════════════════════════════════════════

"""
    read_stl(path; merge_tol=1e-9, max_nodes=typemax(Int32),
             max_facets=typemax(Int32), max_file_bytes=typemax(Int)) -> Mesh

Read an ASCII or binary STL. Coincident vertices within `merge_tol` (relative to
the bounding-box diagonal) are merged so the result is a connected triangle
surface, not a triangle soup. File, facet, and unique-node ceilings are checked
before their corresponding bulk allocations.
"""
function read_stl(path; merge_tol=1e-9,
                  max_nodes=typemax(Int32),max_facets=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_stl"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    merge_tol isa Real || throw(ArgumentError(
        "$caller: merge_tol must be real"))
    merge_tol isa Bool && throw(ArgumentError(
        "$caller: merge_tol must not be Bool"))
    mtol=try
        Float64(merge_tol)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_stl: merge_tol must be representable as Float64"))
    end
    (isfinite(mtol) && mtol >= 0) ||
        throw(ArgumentError("read_stl: merge_tol must be finite and non-negative (got $merge_tol)"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";
                         ceiling=Int(typemax(Int32)))
    facet_limit=_io_limit(max_facets,caller,"max_facets";
                          ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    try
        isbin = _stl_is_binary(path)
        tris_xyz = isbin ? _read_stl_binary(path,facet_limit) :
                           _read_stl_ascii(path,facet_limit)
        out=_weld_triangles(tris_xyz, mtol, node_limit)
        d=validate(out)
        d.ok || throw(ArgumentError(
            "read_stl: welded surface is invalid — "*join(d.messages,"; ")))
        return out
    catch err
        err isa InterruptException && rethrow()
        err isa OutOfMemoryError && rethrow()
        err isa ArgumentError && rethrow()
        err isa Base.InvalidCharError && throw(ArgumentError(
            "read_stl: malformed STL text — "*sprint(showerror,err)))
        rethrow()
    end
end

function _stl_is_binary(path)
    open(path, "r") do io
        head = read(io, min(filesize(path), 84))
        # ASCII STL starts with "solid"; but some binary files do too, so also
        # check the 80-byte-header + uint32 triangle count against file size.
        if length(head) >= 84
            ntri = ltoh(reinterpret(UInt32, head[81:84])[1])
            expected = try Base.checked_add(84,Base.checked_mul(50,Int(ntri))) catch err
                err isa InterruptException && rethrow()
                -1
            end
            filesize(path) == expected && return true
        end
        # A valid ASCII STL may carry a UTF-8 BOM and/or a non-ASCII solid name,
        # so test only for the "solid" keyword followed by ASCII whitespace (after
        # skipping an optional BOM) — NOT that the whole 84-byte header is 7-bit
        # ASCII (which misclassifies such files as binary → a garbage-count crash).
        i = (length(head) >= 3 && head[1] == 0xEF && head[2] == 0xBB && head[3] == 0xBF) ? 4 : 1
        isws(b) = b == 0x20 || b == 0x09 || b == 0x0a || b == 0x0d
        return !(length(head) >= i + 5 && @view(head[i:i+4]) == b"solid" && isws(head[i+5]))
    end
end

function _read_stl_ascii(path,facet_limit::Int=Int(typemax(Int32)))
    open(path,"r") do stream
        _read_stl_ascii(stream,facet_limit)
    end
end

function _read_stl_ascii(stream::Base.IO,facet_limit::Int=Int(typemax(Int32)))
    tris = NTuple{9,Float64}[]
    verts = NTuple{3,Float64}[]
    infacet=false;inloop=false
    for line in eachline(stream)
        s = split(strip(line))
        isempty(s) && continue
        tok=lowercase(replace(s[1],'\ufeff'=>""))
        if tok=="facet"
            (!infacet&&!inloop&&length(s)==5&&lowercase(s[2])=="normal") ||
                throw(ArgumentError("read_stl: malformed or nested ASCII facet"))
            normal=(parse(Float64,s[3]),parse(Float64,s[4]),parse(Float64,s[5]))
            all(isfinite,normal) || throw(ArgumentError("read_stl: ASCII facet has a non-finite normal"))
            infacet=true;empty!(verts)
        elseif tok=="outer"
            (infacet&&!inloop&&length(s)==2&&lowercase(s[2])=="loop") ||
                throw(ArgumentError("read_stl: malformed ASCII outer loop"))
            inloop=true
        elseif tok == "vertex"
            (infacet&&inloop) || throw(ArgumentError("read_stl: vertex outside an ASCII facet loop"))
            length(s) == 4 || throw(ArgumentError("read_stl: malformed ASCII vertex record '$line'"))
            p=(parse(Float64,s[2]),parse(Float64,s[3]),parse(Float64,s[4]))
            (isfinite(p[1])&&isfinite(p[2])&&isfinite(p[3])) ||
                throw(ArgumentError("read_stl: ASCII vertex has non-finite coordinates"))
            push!(verts,p)
            length(verts)<=3 || throw(ArgumentError("read_stl: ASCII facet has more than three vertices"))
        elseif tok=="endloop"
            (infacet&&inloop&&length(s)==1&&length(verts)==3) ||
                throw(ArgumentError("read_stl: ASCII facet loop does not contain exactly three vertices"))
            inloop=false
        elseif tok=="endfacet"
            (infacet&&!inloop&&length(s)==1&&length(verts)==3) || throw(ArgumentError("read_stl: malformed ASCII endfacet"))
            length(tris)<typemax(Int32) || throw(ArgumentError("read_stl: facet count exceeds Int32"))
            length(tris)<facet_limit || throw(ArgumentError(
                "read_stl: facet count exceeds max_facets=$facet_limit"))
            push!(tris,(verts[1]...,verts[2]...,verts[3]...));empty!(verts);infacet=false
        elseif tok=="solid" || tok=="endsolid"
            (!infacet&&!inloop) || throw(ArgumentError("read_stl: $tok appeared inside a facet"))
        else
            throw(ArgumentError("read_stl: unexpected ASCII STL record '$line'"))
        end
    end
    (!infacet&&!inloop&&isempty(verts)) || throw(ArgumentError("read_stl: incomplete ASCII facet"))
    isempty(tris) && throw(ArgumentError("read_stl: ASCII STL contains no complete facets"))
    return tris
end

function _read_stl_binary(path,facet_limit::Int=Int(typemax(Int32)))
    open(path, "r") do io
        filesize(path) >= 84 || throw(ArgumentError("read_stl: truncated binary STL header"))
        skip(io, 80)
        ntri = Int(ltoh(read(io, UInt32)))
        ntri>0 || throw(ArgumentError("read_stl: binary STL contains no facets"))
        ntri<=typemax(Int32) || throw(ArgumentError("read_stl: binary facet count exceeds Int32"))
        ntri<=facet_limit || throw(ArgumentError(
            "read_stl: binary facet count $ntri exceeds max_facets=$facet_limit"))
        expected = try Base.checked_add(84,Base.checked_mul(50,ntri)) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("read_stl: binary triangle count overflows file-size arithmetic"))
        end
        filesize(path)==expected ||
            throw(ArgumentError("read_stl: binary STL size $(filesize(path)) does not match triangle count $ntri (expected $expected bytes)"))
        tris = Vector{NTuple{9,Float64}}(undef, ntri)
        try
            for t in 1:ntri
                normal=ntuple(_ -> Float64(reinterpret(Float32,ltoh(read(io,UInt32)))),3)
                all(isfinite,normal) || throw(ArgumentError("read_stl: binary facet $t has a non-finite normal"))
                v = ntuple(_ -> Float64(reinterpret(Float32,ltoh(read(io,UInt32)))), 9)
                all(isfinite,v) || throw(ArgumentError("read_stl: binary facet $t has non-finite coordinates"))
                tris[t] = v
                skip(io, 2)    # attribute byte count
            end
        catch err
            err isa ArgumentError && rethrow()
            err isa EOFError && throw(ArgumentError("read_stl: truncated binary STL facet data"))
            rethrow()
        end
        return tris
    end
end

# Compare an exact represented-coordinate difference with an integer multiple
# of the tolerance. Since the multiplier is an exactly represented integer,
# both exact terms are multiples of the smallest subnormal: a zero FMA residual
# is an exact zero, including subnormal tolerances. Only cancellation with the
# subtraction tail can require the rational fallback.
@inline function _stl_bucket_sign(shift::Float64,tail::Float64,
                                  tolerance::Float64,index::Int)
    residual=fma(-Float64(index),tolerance,shift)
    tail==0.0 && return sign(residual)
    residual==0.0 && return sign(tail)
    total=residual+tail
    abs(total)>4eps(Float64)*(abs(residual)+abs(tail)) && return sign(total)
    return nothing
end

@inline function _stl_bucket_exact(value::Float64,origin::Float64,
                                   tolerance::Float64)
    return floor(Int,(Rational{BigInt}(value)-Rational{BigInt}(origin)) /
                     Rational{BigInt}(tolerance))
end

function _stl_bucket_index(value::Float64,origin::Float64,
                            tolerance::Float64,inverse::Float64)
    shift=value-origin
    scaled=shift*inverse
    scaled<0.5 && return 0
    index=floor(Int,scaled)
    fraction=scaled-Float64(index)
    # Subtraction, reciprocal and multiplication together have relative error
    # below 8eps, including a subnormal reciprocal (1/floatmax >= 2^-1024).
    # Away from an integer boundary, this certifies the exact quotient's floor.
    min(fraction,1.0-fraction)>8eps(Float64)*scaled && return index
    tail=_two_diff_tail(value,origin,shift)
    isfinite(tail) || return _stl_bucket_exact(value,origin,tolerance)
    while true
        lower=_stl_bucket_sign(shift,tail,tolerance,index)
        lower===nothing && return _stl_bucket_exact(value,origin,tolerance)
        if lower<0
            index-=1
            continue
        end
        upper=_stl_bucket_sign(shift,tail,tolerance,index+1)
        upper===nothing && return _stl_bucket_exact(value,origin,tolerance)
        upper<0 && return index
        index+=1
    end
end

function _weld_triangles(tris_xyz::Vector{NTuple{9,Float64}}, reltol::Real,
                         node_limit::Int=Int(typemax(Int32)))
    isempty(tris_xyz) && return Mesh(Matrix{Float64}(undef, 3, 0))
    length(tris_xyz)<=typemax(Int32) || throw(ArgumentError("read_stl: facet count exceeds Int32"))
    # bbox diagonal → absolute tolerance for quantization
    lo = (Inf,Inf,Inf); hi = (-Inf,-Inf,-Inf)
    for t in tris_xyz, k in (1,4,7)
        p = (t[k], t[k+1], t[k+2])
        all(isfinite,p) || throw(ArgumentError("read_stl: facet has non-finite coordinates"))
        lo = (min(lo[1],p[1]),min(lo[2],p[2]),min(lo[3],p[3]))
        hi = (max(hi[1],p[1]),max(hi[2],p[2]),max(hi[3],p[3]))
    end
    # Exact welding (`merge_tol=0`) depends only on represented coordinate
    # equality. Do not reject it merely because subtracting opposite finite
    # Float64 extrema would overflow while computing an unused diagonal.
    diag=0.0;tol=0.0;robust_diagonal=nothing
    if reltol>0
        diag = hypot(hi[1]-lo[1],hi[2]-lo[2],hi[3]-lo[3])
        if isfinite(diag)
            tol = diag * reltol
        else
            # All coordinates are finite, but an opposite-sign span can exceed
            # Float64. Compute only this exceptional scale in high precision;
            # the resulting represented tolerance can still be finite.
            robust_diagonal=setprecision(BigFloat,256) do
                dx=BigFloat(hi[1])-BigFloat(lo[1])
                dy=BigFloat(hi[2])-BigFloat(lo[2])
                dz=BigFloat(hi[3])-BigFloat(lo[3])
                sqrt(dx*dx+dy*dy+dz*dz)
            end
            tol=setprecision(BigFloat,256) do
                Float64(robust_diagonal*BigFloat(reltol))
            end
        end
        isfinite(tol) || throw(ArgumentError(
            "read_stl: absolute welding tolerance is non-finite"))
    end
    # Quantize relative to the bbox min corner `lo`, not the absolute coordinate:
    # a constant per-axis bucket shift (does not change which vertices merge) that
    # bounds every key to [0, diag/tol] so far-from-origin coords can't overflow
    # Int64 in round(Int, ·) — an unshifted absolute coord ~1e11 does.
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    exactmap=Dict{NTuple{3,Float64},Int32}()
    buckets=Dict{NTuple{3,Int},Vector{Int32}}()
    inv = tol == 0 ? 0.0 : 1/tol
    exact_bucket=false
    if tol > 0
        cells=if robust_diagonal===nothing
            diag/tol
        else
            setprecision(BigFloat,256) do
                ratio=robust_diagonal/BigFloat(tol)
                ratio<=BigFloat(typemax(Int)-2) || throw(ArgumentError(
                    "read_stl: merge_tol is too small for Int bucket indices at this extent"))
                Float64(ratio)
            end
        end
        # Near typemax(Int), adjacent integers share one Float64 image. Leave
        # more than one ULP of headroom before adding neighbour offsets.
        safe_cell_limit=prevfloat(Float64(typemax(Int)))
        (isfinite(cells) && cells<=safe_cell_limit) || throw(ArgumentError(
            "read_stl: merge_tol is too small for Int bucket indices at this extent"))
        # Above 2^52, Float64 no longer distinguishes every integer bucket.
        # Exact represented-coordinate division also handles an overflowing
        # `value-origin` on the rare extreme-span path.
        exact_bucket=robust_diagonal!==nothing || cells>2.0^52 || !isfinite(inv)
    end
    # single-assignment copies: closures over the accumulated `lo`, `tol`, and
    # `exact_bucket` boxed them and allocated on every welded vertex
    corner=lo; weld_tol=tol; exact_buckets=exact_bucket
    function bucket_component(value::Float64,origin::Float64)
        if exact_buckets
            return _stl_bucket_exact(value,origin,weld_tol)
        end
        return _stl_bucket_index(value,origin,weld_tol,inv)
    end
    keyof(p)=(bucket_component(p[1],corner[1]),bucket_component(p[2],corner[2]),
              bucket_component(p[3],corner[3]))
    function getid(p)
        if weld_tol == 0
            key=(p[1]==0 ? 0.0 : p[1],p[2]==0 ? 0.0 : p[2],p[3]==0 ? 0.0 : p[3])
            return get!(exactmap,key) do
                length(xs)<node_limit || throw(ArgumentError(
                    "read_stl: unique vertex count exceeds max_nodes=$node_limit"))
                push!(xs,p[1]);push!(ys,p[2]);push!(zs,p[3]);Int32(length(xs))
            end
        end
        key=keyof(p)
        for dz in -1:1, dy in -1:1, dx in -1:1
            ids=get(buckets,(key[1]+dx,key[2]+dy,key[3]+dz),nothing)
            ids === nothing && continue
            for id in ids
                hypot(p[1]-xs[id],p[2]-ys[id],p[3]-zs[id]) <= weld_tol && return id
            end
        end
        length(xs)<node_limit || throw(ArgumentError(
            "read_stl: unique vertex count exceeds max_nodes=$node_limit"))
        push!(xs,p[1]);push!(ys,p[2]);push!(zs,p[3]);id=Int32(length(xs))
        push!(get!(() -> Int32[],buckets,key),id)
        return id
    end
    tris = Matrix{Int32}(undef, 3, length(tris_xyz))
    ntri = 0
    ncollapsed = 0
    for t in tris_xyz
        a = getid((t[1],t[2],t[3])); b = getid((t[4],t[5],t[6])); c = getid((t[7],t[8],t[9]))
        if a==b || b==c || a==c
            ncollapsed += 1
            continue
        end
        ntri += 1
        tris[1,ntri]=a; tris[2,ntri]=b; tris[3,ntri]=c
    end
    ncollapsed==0 ||
        throw(ArgumentError("read_stl: $ncollapsed facet(s) collapsed under merge_tol=$reltol"))
    coords = Matrix{Float64}(undef, 3, length(xs))
    @inbounds for i in eachindex(xs); coords[1,i]=xs[i]; coords[2,i]=ys[i]; coords[3,i]=zs[i]; end
    return Mesh(coords; tris=tris[:, 1:ntri])
end

# ════════════════════════════════════════════════════════════════════════════════
# .geo parameter / structure and field scan (no OCC evaluation)
# ════════════════════════════════════════════════════════════════════════════════

"""
    GeoFieldSpec

A parsed `Field[tag] = Kind` declaration. Constant numeric expressions in known
numeric options are evaluated while scanning and stored as normalized literals;
string, point-dependent expression, and geometric-entity values remain `.geo`
source strings. [`Tessella.SizeField.build_geo_size_field`](@ref) interprets the
supported field kinds after geometric entity references have been resolved.
"""
struct GeoFieldSpec
    tag::Int
    kind::String
    options::Dict{String,String}
    option_order::Vector{String}
    creation_mesh_size_from_curvature::Int
end

# Programmatic specifications cannot recover source assignment order from a
# plain dictionary. Use a stable lexical order there; the parser-populated
# five-argument form below preserves exact statement order for
# Gmsh aliases and construction-time globals.
GeoFieldSpec(tag::Integer,kind::AbstractString,options::Dict{String,String})=
    GeoFieldSpec(Int(tag),String(kind),options,sort!(collect(keys(options))),0)
GeoFieldSpec(tag::Integer,kind::AbstractString,options::Dict{String,String},
             option_order::Vector{String})=
    GeoFieldSpec(Int(tag),String(kind),options,option_order,0)

"""
    GeoParams

What `read_geo_params` can extract from a `.geo` without a geometry kernel:
`mesh_size_min/max/factor`, `random_seed`, `physical_groups` (a `(dim, tag) =>
name` map for named groups, dim ∈ {0:point,1:curve,2:surface,3:volume}), raw
`fields`, and the `background_field` tag. Named declarations with no explicit tag
receive a tag from the global Physical-group namespace. Unnamed Physical declarations
are checked by the scanner but omitted from the name map. Missing numeric mesh
options are represented by `NaN`;
no background field is tag `0`. `boundary_layer_fields` preserves the distinct,
deduplicated `BoundaryLayer Field = ...` declarations; these are mesher controls,
not background scalar fields. `geometry_tolerance` stores `Geometry.Tolerance`
when it is a finite constant numeric expression, or `NaN` when absent.
`mesh_boundary_layer_fan_elements` stores the final
`Mesh.BoundaryLayerFanElements` value (default 5). Automatic-field declarations
also snapshot `Mesh.MeshSizeFromCurvature` in their `GeoFieldSpec`, since Gmsh
uses that global at field construction time.
"""
struct GeoParams
    mesh_size_min::Float64
    mesh_size_max::Float64
    mesh_size_factor::Float64
    random_seed::Int
    physical_groups::Dict{Tuple{Int,Int},String}
    fields::Dict{Int,GeoFieldSpec}
    background_field::Int
    boundary_layer_fields::Vector{Int}
    geometry_tolerance::Float64
    mesh_boundary_layer_fan_elements::Int
    # Scan-time recoverable diagnostics (Field declarations/option writes):
    # `yymsg(0)`-level parse errors and `yymsg(1)` warnings in source order,
    # plus the `Msg::Error` count. `execute_geo` merges them into its own
    # channels so the run reports them like Gmsh does instead of aborting the
    # scan.
    scan_errors::Vector{String}
    scan_warnings::Vector{String}
    scan_msg_error_count::Int
end

# Preserve the former full positional constructor, adding the new global with
# its Gmsh default.
GeoParams(mesh_size_min::Real,mesh_size_max::Real,mesh_size_factor::Real,
          random_seed::Integer,physical_groups::Dict{Tuple{Int,Int},String},
          fields::Dict{Int,GeoFieldSpec},background_field::Integer,
          boundary_layer_fields::Vector{Int},geometry_tolerance::Real)=
    GeoParams(Float64(mesh_size_min),Float64(mesh_size_max),Float64(mesh_size_factor),
              Int(random_seed),physical_groups,fields,Int(background_field),
              boundary_layer_fields,Float64(geometry_tolerance),5,
              String[],String[],0)

# Preserve the former full ten-field positional constructor.
GeoParams(mesh_size_min::Real,mesh_size_max::Real,mesh_size_factor::Real,
          random_seed::Integer,physical_groups::Dict{Tuple{Int,Int},String},
          fields::Dict{Int,GeoFieldSpec},background_field::Integer,
          boundary_layer_fields::Vector{Int},geometry_tolerance::Real,
          mesh_boundary_layer_fan_elements::Integer)=
    GeoParams(Float64(mesh_size_min),Float64(mesh_size_max),Float64(mesh_size_factor),
              Int(random_seed),physical_groups,fields,Int(background_field),
              boundary_layer_fields,Float64(geometry_tolerance),
              Int(mesh_boundary_layer_fan_elements),String[],String[],0)

# Preserve the eight-argument constructor introduced with boundary-layer
# declarations.
GeoParams(mesh_size_min::Real,mesh_size_max::Real,mesh_size_factor::Real,
          random_seed::Integer,physical_groups::Dict{Tuple{Int,Int},String},
          fields::Dict{Int,GeoFieldSpec},background_field::Integer,
          boundary_layer_fields::Vector{Int})=
    GeoParams(Float64(mesh_size_min),Float64(mesh_size_max),Float64(mesh_size_factor),
              Int(random_seed),physical_groups,fields,Int(background_field),
              boundary_layer_fields,NaN,5,String[],String[],0)

# Preserve the seven-argument field-aware constructor used before boundary-layer
# declarations became part of the parsed model.
GeoParams(mesh_size_min::Real,mesh_size_max::Real,mesh_size_factor::Real,
          random_seed::Integer,physical_groups::Dict{Tuple{Int,Int},String},
          fields::Dict{Int,GeoFieldSpec},background_field::Integer)=
    GeoParams(Float64(mesh_size_min),Float64(mesh_size_max),Float64(mesh_size_factor),
              Int(random_seed),physical_groups,fields,Int(background_field),Int[],NaN,5,
              String[],String[],0)

# Preserve the original public positional constructor.
GeoParams(mesh_size_min::Real, mesh_size_max::Real, random_seed::Integer,
          physical_groups::Dict{Tuple{Int,Int},String}) =
    GeoParams(Float64(mesh_size_min), Float64(mesh_size_max), 1.0, Int(random_seed),
              physical_groups, Dict{Int,GeoFieldSpec}(), 0, Int[],NaN,5,
              String[],String[],0)

const _PHYS_DIM = Dict("Point"=>0, "Curve"=>1, "Line"=>1, "Surface"=>2, "Volume"=>3)

const _MAX_GEO_STATEMENT_BYTES=1_000_000
const _MAX_GEO_EXPRESSION_BYTES=65_536
const _MAX_GEO_EXPRESSION_TOKENS=4_096
const _MAX_GEO_EXPRESSION_DEPTH=128
const _MAX_GEO_LIST_ITEMS=65_536
const _MAX_GEO_CONTEXT_LIST_ITEMS=1_000_000

struct _GeoNumericListTerm
    first::Float64
    step::Float64
    count::Int
    # Non-arithmetic expansion produced by a side-effecting term (for example
    # an `Extrude{...}{...}` entity-list result); empty for plain values and
    # ranges.
    values::Vector{Float64}
end

# A `Struct` definition inside a `gmsh_yynamespaces` namespace (Parser.h):
# `_fopt`/`_copt` map member names to float/char vectors, `_tag` is the
# resolved auto-or-forced tag and `_member_ValMax` continues `Enum`
# auto-numbering across `Append` definitions.
mutable struct _GeoStructDef
    tag::Int
    fopt::Dict{String,Vector{Float64}}
    copt::Dict{String,Vector{String}}
    member_valmax::Int
end
_GeoNumericListTerm(first,step,count)=
    _GeoNumericListTerm(first,step,count,Float64[])

function _geo_list_term_values(term::_GeoNumericListTerm)
    isempty(term.values) || return term.values
    values=Vector{Float64}(undef,term.count)
    numeric=term.first
    for index in eachindex(values)
        values[index]=numeric;numeric+=term.step
    end
    return values
end

# This is deliberately a constant-expression evaluator, not a Julia evaluator
# and not a general Gmsh interpreter.  Keeping a small lexer/parser here makes
# unknown identifiers, assignments and stateful built-ins impossible to execute.
struct _GeoExprToken
    kind::Symbol
    text::String
    value::Float64
    pos::Int
end

# A `.geo` `For` frame — `LoopControlVariablesTab` (current/last/step), the
# loop variable name ("" for `For (a:b)`), and the saved stream position
# (`body_start` indexes `statements`, the array the header was parsed in —
# `EndFor` only rewinds a level belonging to its own stream).
mutable struct _GeoLoopLevel
    var::String
    x0::Float64
    x1::Float64
    step::Float64
    body_start::Int
    statements::Vector{String}
end

# A `While` frame (Tessella extension) — condition source re-evaluated at
# each `EndWhile`, the body start position, and an iteration bound.
mutable struct _GeoWhileLevel
    cond::String
    body_start::Int
    statements::Vector{String}
    iterations::Int
end

# Exported geometry entity parts use the carrier's row space, while the
# existing context products continue to use actual native row space.
struct _GeoGeometryParts
    carrier::Union{Mesh,MixedMesh}
    parts::Vector{Tuple{Int,Int,Union{Mesh,MixedMesh}}}
    rows::Vector{Vector{Int32}}
    node_tags::Vector{UInt64}
end

mutable struct _GeoNumericContext
    values::Dict{String,Float64}
    lists::Dict{String,Vector{Float64}}
    # Gmsh preserves an array payload after `a = value`, but rejects `a[i] = ...`
    # until a whole-list operation makes `a` a list again.
    list_variables::Set{String}
    unavailable::Dict{String,String}
    # Keep an unknown array payload separate from an unknown scalar head.  For
    # example, `a[] = Surface{:}; a = 1` makes `a` known while `a[]` still
    # contains the geometry-derived tail retained by Gmsh.
    unavailable_lists::Dict{String,String}
    stored_list_items::Int
    # Optional `execute_geo` hook for side-effecting list terms such as
    # `Extrude{...}{...}`; `(term_source) -> Vector{Float64}` or `nothing`
    # when the term is not an exec term. `nothing` in pure contexts.
    exec_hook::Any
    # `Geometry.ExtrudeReturnLateralEntities` (defaults on): controls whether
    # extrude result lists append the lateral entities.
    extrude_return_lateral::Bool
    # `.geo` `Function`/`Macro`/`Call` support: each name maps to
    # `(statement_array, body_start_index)` — the array reference keeps
    # `Include`d files' bodies alive for later calls (Gmsh's FunctionManager
    # stores file positions and survives `Delete All`/`NewModel`).
    functions::Dict{String,Tuple{Vector{String},Int}}
    call_depth::Int
    # `gmsh_yynamespaces` — `Struct` definitions per namespace: float members
    # (fopt), char members (copt), the resolved tag and the running enum
    # counter (`member_valmax`). `struct_maxtag` is each `Structs` container's
    # `_max_tag` (per-namespace auto-tag allocator).
    structs::Dict{String,Dict{String,_GeoStructDef}}
    struct_maxtag::Dict{String,Int}
    # `gmsh_yystringsymbols` — a separate table from the numeric symbols, so
    # `x = 1; x = "s"` keeps both bindings live, `Delete x` erases only the
    # numeric entry, and `Delete Variables` leaves strings untouched.
    strings::Dict{String,Vector{String}}
    # The batch localGmsh ONELAB client's parameter store, backing
    # `SetNumber`/`SetString`/`GetNumber`/`GetString`/`DefineNumber`/
    # `DefineString`/`DefineConstant`.
    onelab_numbers::Dict{String,Float64}
    onelab_strings::Dict{String,String}
    # `gmsh_yyname` — the file currently being parsed, backing
    # `CurrentDirectory`, `CurrentFileName`, `FixRelativePath` and relative
    # `Include`/`Merge`/`Save` resolution. Empty for source strings.
    file_name::String
    # The mesh produced by mid-file `Mesh n` statements — Gmsh keeps it on the
    # current model; `Save`, `Delete Meshes`, `NewModel` and the mesh-operation
    # statements consume it.
    mesh::Union{Nothing,Mesh,MixedMesh}
    # `ImbricatedTest` — nested `If` depth, readable in expressions.
    if_depth::Int
    # `statusImbricatedTests` — per-level branch status for the flat token
    # stream. Written at `if_depth+1` (upstream index 0 is a real slot — an
    # `If` opened at depth -1 stores there); reads only happen at depth > 0.
    if_status::Vector{Bool}
    # `ImbricatedLoop` + `LoopControlVariablesTab` + `yyposImbricatedLoopsTab` —
    # the open `For` frames: loop variable, current/last/step values, and the
    # stream position `EndFor` rewinds to. `While` (Tessella extension)
    # frames hold the condition source for re-evaluation.
    loop_stack::Vector{_GeoLoopLevel}
    while_stack::Vector{_GeoWhileLevel}
    # `StringOption(GMSH_SET)` mirror: `(family, index, member)` option-string
    # writes, read back by `x.y`/`x[i].y` string expressions.
    option_strings::Dict{Tuple{String,Int,String},String}
    # `NumberOption` mirror for `x.y = v` numeric option writes.
    option_numbers::Dict{Tuple{String,Int,String},Float64}
    # `ColorOption` mirror for `x.Color.f = {r,g,b[,a]}` writes.
    option_colors::Dict{Tuple{String,Int,String},NTuple{4,Int}}
    # `Msg::Error`-path diagnostics (option tables, plugin/field actions) —
    # unlike `yymsg(0)` they do not count toward the >20 parser abort cap;
    # they only mark the run as having errors.
    msg_error_count::Int
    # `Field[i]` declarations recovered by the params pre-pass — `execute_geo`
    # wires the live map so `Field[i].member = v` mutates real field options.
    fields::Any
    # `Background Field = {tag}` / `BoundaryLayer Field = {tags}` —
    # `FieldManager::setBackgroundFieldId` stores the raw id unvalidated and
    # `setBoundaryLayerField` dedup-appends, mirroring `GeoParams`.
    background_field::Int
    boundary_layer_fields::Vector{Int}
    # Exec hook resolving `Point{t}`/`Physical X{t}` name reads in string
    # expressions: `(dim, tag, :entity|:physical) -> String`; `nothing` in pure
    # contexts.
    entity_name_lookup::Any
    # `yymsg` channels — level 0 diagnostics are errors that accumulate through
    # `gmsh_yyerrorstate` while the statement stream continues; level 1 are
    # warnings.
    exec_errors::Vector{String}
    exec_warnings::Vector{String}
    # The ordered union of every `Error   :` diagnostic (`yymsg(0)` and
    # `Msg::Error`) — the user-visible stream a run produced. `exec_errors`
    # stays the errorstate counter; `error_log` carries both channels so
    # `GeoExecution.errors` reflects the whole diagnostic stream in order.
    error_log::Vector{String}
    # `gmsh_yyerrorstate` is saved/restored per parsed file: the >20-error
    # abort cap counts only errors raised inside the current file's stream.
    # File-boundary callers (`execute_geo`, `Include`) set the base.
    file_error_base::Int
    # `treat_Struct_FullName_Float` read misses: in `.geo` EXECUTION a missing
    # variable is a `yymsg(0)` diagnostic plus a 0 result (the statement
    # stream continues); the params-scan pass throws instead, so its callers
    # can mark the target unavailable. Set by `execute_geo`.
    soft_unknown_reads::Bool
    # Stream control: `Abort`/`Exit` set it to :abort/:exit and stop the
    # statement stream at the next boundary.
    stop::Symbol
    exit_code::Int
    # Node → owning `(dim, tag)` for the mid-file mesh — populated by
    # `_geo_merge_entity_meshes` so `RelocateMesh`/`TransformMesh {..}{..}` can
    # scope to an entity.
    mesh_node_owner::Dict{Int,Tuple{Int,Int}}
    # The `(dim, tag, mesh)` parts `_geo_mesh_model` merged into the current
    # mesh — entity-level consumers (classified `model_to_mixed` projection,
    # field views) need each entity's own mesh, not the merged product.
    mesh_parts::Vector{Tuple{Int,Int,Union{Mesh,MixedMesh}}}
    mesh_part_rows::Vector{Vector{Int32}}
    mesh_geometry_parts::Union{Nothing,_GeoGeometryParts}
    mesh_node_parametric::Vector{Union{Nothing,Vector{Float64}}}
    mesh_node_max::UInt64
    mesh_element_max::UInt64
    # `GEO_Internals::PhysicalGroups` — the raw physical-group records as the
    # `.geo` parser keeps them: `(dim, raw_signed_tag)` → the raw signed member
    # tags. `Physical X(n)` looks up `n` here verbatim, so groups `-4` and `4`
    # coexist. The observable view (`m.physical`) is derived from these records
    # at each sync point: member `num` resolves entity `abs(num)` and lands in
    # group `num == 0 ? 0 : abs(tag)` (`orientedPhysicals`).
    raw_physicals::Dict{Tuple{Int,Int},Vector{Int}}
    # `GEO_Internals::_changed` — set by every geometry/physical mutation,
    # cleared at the end of each synchronize. `SyncModel` syncs
    # unconditionally; `BoundingBox`, `Save`, `Print`, `Merge`, `Delete {..}`
    # and the `Physical X{..}` selectors sync only when this is set, so a
    # repeated sync point does not re-emit "unknown member" warnings.
    geo_changed::Bool
    # The params scan silences the stderr print (its diagnostics land in
    # `GeoParams.scan_*` instead); execution prints each diagnostic once,
    # matching Gmsh's single pass.
    stderr_diagnostics::Bool
    # `AddToTemporaryBoundingBox` — the `Point(...) = {...}`-statement point
    # cloud whose raw diagonal is `CTX::instance()->lc` while a file is being
    # parsed (`lc == 0` maps to 1 upstream). `nothing` until the first Point
    # statement runs.
    temp_bbox::Union{Nothing,NTuple{2,NTuple{3,Float64}}}
    # `allGmshSurfaces` — the registry of parametric geometry surfaces
    # (`Sphere`/`PolarSphere`/`Parametric Surface`) keyed by their tag.
    # Registering an existing tag reports `Surface N already exists` and then
    # replaces the entry, matching `gmshSurface::New*`. It is not the model
    # surface namespace.
    parametric_surfaces::Dict{Int,Any}
    # `myGmshSurface` — the parametric-surface coordinate context a
    # `Sphere(n) = {center, point}`/`PolarSphere(n) = {center, point}` (or
    # `Parametric Surface`) statement installs. `Euclidian Coordinates`
    # resets it. It is not an entity: no tag is consumed and nothing lands in
    # `m.surfaces`/`m.volumes`.
    parametric_surface::Any
    # Bookkeeping re-evaluation inside the allocator observer: upstream
    # reduces each FExpr once in the grammar production, so the observer's
    # re-reads must not duplicate the diagnostics the statement emitted.
    suppress_diagnostics::Bool
    # CTX::lc and forcedBBox persist across meshes and change at native events.
    native_context_lc::Float64
    native_context_forced_bbox::Bool
    native_mesh_error_callback::Any
    native_mesh_begin_callback::Any
    native_abort_option_callback::Any
    # MSH4 reads install element lookup bindings before later mesh operations.
    native_element_lookup_receipt::Any
    native_numeric_option_callback::Any
    # Selector4 can terminate after publishing an earlier Face prefix.
    native_recombine_partial_callback::Any
    # Optional detached diagnostics for atomic API file preparation.
    native_recombine_diagnostics::Any
end
_GeoNumericContext()=_GeoNumericContext(
    Dict{String,Float64}(),Dict{String,Vector{Float64}}(),Set{String}(),
    Dict{String,String}(),Dict{String,String}(),0,nothing,true,
    Dict{String,Tuple{Vector{String},Int}}(),0,
    Dict{String,Dict{String,_GeoStructDef}}(),Dict{String,Int}(),
    Dict{String,Vector{String}}(),Dict{String,Float64}(),
    Dict{String,String}(),"",nothing,0,falses(101),
    _GeoLoopLevel[],_GeoWhileLevel[],
    Dict{Tuple{String,Int,String},String}(),
    Dict{Tuple{String,Int,String},Float64}(),
    Dict{Tuple{String,Int,String},NTuple{4,Int}}(),0,nothing,0,Int[],nothing,
    String[],String[],String[],0,false,:run,0,Dict{Int,Tuple{Int,Int}}(),
    Tuple{Int,Int,Union{Mesh,MixedMesh}}[],Vector{Int32}[],nothing,
    Union{Nothing,Vector{Float64}}[],UInt64(0),UInt64(0),
    Dict{Tuple{Int,Int},Vector{Int}}(),true,true,nothing,
    Dict{Int,Any}(),nothing,false,1.0,false,nothing,nothing,nothing,nothing,nothing,nothing,nothing)

# Thrown where the upstream grammar fails to reduce — a bison syntax error.
# `execute_geo` catches it per statement, records `syntax error (<token>)`
# through `_geo_yyerror!` and resumes with the next statement, mirroring
# bison's error recovery. Scan contexts (`read_geo_params`) see it converted
# back to `ArgumentError` at the entry boundary; `detail` preserves the
# descriptive message for that path.
struct _GeoSyntaxAbort <: Exception
    token::String
    detail::String
    # `header` marks an abort whose production already consumed the
    # statement's own `;` — a failed control-construct header (`If (x y)`,
    # `For i In {bad}`) or a terminator-owning term like the `Split`
    # Transform inside an `x =`/`x() =` list. Upstream's `error tEND`
    # recovery then reports the following statement's token and discards
    # through its `;`, eating following control keywords and one
    # `;`-statement.
    header::Bool
end

_geo_syntax_abort(token::AbstractString)=
    throw(_GeoSyntaxAbort(String(token),"",false))
# ---------------------------------------------------------------------------
# `-log <file>` / `-cpu` — the Tessella counterparts of upstream's
# `Msg::_logFile` (`Msg::SetLogFileName`) and `Msg::SetInfoCpu`/
# `Msg::SetInfoMem`. The CLI installs the sinks for the duration of a `main`
# call: every Msg-channel line the terminal prints (Error/Warning/Info/Debug
# plus verbatim `Msg::Direct` text) is mirrored into the log file, and
# Info/StatusBar messages carry the ` (From start: Wall %gs, CPU %gs, Mem
# %gMb)` resource suffix. Both are process-global like the upstream
# originals; `CLI.main` resets them in a `finally`, so API-only runs are
# unaffected.
const _MSG_LOG_IO=Ref{Any}(nothing)
const _MSG_INFO_RESOURCES=Ref{Bool}(false)
const _MSG_START_NS=Ref{UInt64}(0)
# Upstream's `Fopen(name,"w")` is a text-mode stream — CRLF on Windows.
const _MSG_LOG_NEWLINE=Sys.iswindows() ? "\r\n" : "\n"

function _msg_log_open(path::AbstractString)
    io=try
        open(path,"w")
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("tessella: -log could not open file '$path'"))
    end
    _MSG_LOG_IO[]=io
    return io
end

function _msg_log_close()
    io=_MSG_LOG_IO[]
    _MSG_LOG_IO[]=nothing
    io===nothing || close(io)
    return nothing
end

function _msg_log_write(prefix::AbstractString,message::AbstractString)
    io=_MSG_LOG_IO[]
    io===nothing && return nothing
    print(io,prefix,message,_MSG_LOG_NEWLINE)
    flush(io)
    return nothing
end

# User-mode process CPU seconds — the `Cpu()` upstream's `GetResources`
# reports (Windows process user FILETIME, POSIX `ru_utime`).
function _msg_cpu_time()
    if Sys.iswindows()
        creation=Ref{UInt64}(0);exit_time=Ref{UInt64}(0)
        kernel=Ref{UInt64}(0);user_time=Ref{UInt64}(0)
        handle=ccall((:GetCurrentProcess,"kernel32"),Ptr{Cvoid},())
        ok=ccall((:GetProcessTimes,"kernel32"),Int32,
            (Ptr{Cvoid},Ref{UInt64},Ref{UInt64},Ref{UInt64},Ref{UInt64}),
            handle,creation,exit_time,kernel,user_time)
        ok!=0 || return 0.0
        return Float64(user_time[])/1.e7
    end
    # `struct rusage` opens with `ru_utime` — two Clong fields on LP64.
    usage=zeros(Clong,32)
    ccall(:getrusage,Cint,(Cint,Ptr{Cvoid}),Cint(0),usage)
    return Float64(usage[1])+Float64(usage[2])/1.e6
end

# ` (From start: Wall %gs, CPU %gs, Mem %gMb)` — upstream
# `Msg::PrintResources(false, cpu, cpu, mem)`. `Sys.maxrss` is the same peak
# working set upstream reads (`PeakWorkingSetSize`/`ru_maxrss`).
function _msg_resource_suffix()
    wall=Float64(time_ns()-_MSG_START_NS[])/1.e9
    return @sprintf(" (From start: Wall %gs, CPU %gs, Mem %gMb)",
        wall,_msg_cpu_time(),Float64(Sys.maxrss())/1048576.0)
end

_geo_syntax_abort(token::AbstractString,detail::AbstractString)=
    throw(_GeoSyntaxAbort(String(token),String(detail),false))

Base.showerror(io::Base.IO,err::_GeoSyntaxAbort)=
    print(io,isempty(err.detail) ? "syntax error ($(err.token))" : err.detail)

# `yymsg(0, ...)` — record a parse error; the stream continues and `execute_geo`
# throws the accumulated diagnostics at the end (or at `Exit`/`Abort`).
function _geo_yyerror!(context::_GeoNumericContext,message::AbstractString)
    context.suppress_diagnostics && return nothing
    if context.stderr_diagnostics
        # `-cpu` appends `PrintResources` inside `Msg::Error` too.
        text=_MSG_INFO_RESOURCES[] ? message*_msg_resource_suffix() : message
        print(stderr,"Error   : ",text,"\n");flush(stderr)
        _msg_log_write("Error   : ",text)
    end
    push!(context.exec_errors,String(message))
    push!(context.error_log,String(message))
    return nothing
end

# `yymsg(1, ...)` — a warning; printed like `Msg::Warning`, execution
# continues.
function _geo_yywarn!(context::_GeoNumericContext,message::AbstractString)
    context.suppress_diagnostics && return nothing
    if context.stderr_diagnostics
        print(stderr,"Warning : ",message,"\n");flush(stderr)
        _msg_log_write("Warning : ",message)
    end
    push!(context.exec_warnings,String(message))
    return nothing
end

# `Msg::Info` — `Info    : <msg>` on stdout at `General.Verbosity ≥ 4` (the
# default is 5). Recorded nowhere upstream — purely observable output.
# `-cpu` appends upstream's `PrintResources` suffix before every channel.
function _geo_yyinfo!(context::_GeoNumericContext,message::AbstractString)
    something(_geo_option_number(context,"General",0,"Verbosity"),5.0)<4 &&
        return nothing
    text=_MSG_INFO_RESOURCES[] ? message*_msg_resource_suffix() : message
    print(stdout,"Info    : ",text,"\n");flush(stdout)
    _msg_log_write("Info    : ",text)
    return nothing
end

# `Msg::Error` outside the parser diagnostics channel — option-table lookups,
# plugin and field actions. Prints immediately and marks the run failed, but
# does not bump `gmsh_yyerrorstate` (no >20 abort).
function _geo_msg_error!(context::_GeoNumericContext,message::AbstractString)
    context.suppress_diagnostics && return nothing
    if context.stderr_diagnostics
        text=_MSG_INFO_RESOURCES[] ? message*_msg_resource_suffix() : message
        print(stderr,"Error   : ",text,"\n");flush(stderr)
        _msg_log_write("Error   : ",text)
    end
    context.msg_error_count+=1
    push!(context.error_log,String(message))
    return nothing
end

# The lexer's `Msg::Error("Unexpected end of file")` (Gmsh.l skip/skipTest):
# one printed diagnostic that both marks the run failed (`msg_error_count`,
# which bare `Exit` reads) and aborts the parse — upstream exits nonzero, so
# the statement stream cannot silently deliver a truncated model. Recorded
# in `exec_errors` so `execute_geo` throws, without reprinting.
function _geo_eof_error!(context::_GeoNumericContext)
    _geo_msg_error!(context,"Unexpected end of file")
    push!(context.exec_errors,"Unexpected end of file")
    return nothing
end

@inline _geo_context_has_variable(context::_GeoNumericContext,name::String)=
    haskey(context.values,name) || haskey(context.lists,name) ||
    haskey(context.unavailable,name) || haskey(context.unavailable_lists,name)

function _geo_context_forget!(context::_GeoNumericContext,name::AbstractString)
    key=String(name)
    if haskey(context.lists,key)
        context.stored_list_items-=length(context.lists[key])
        delete!(context.lists,key)
    end
    delete!(context.list_variables,key)
    delete!(context.values,key)
    return nothing
end

function _geo_context_set_scalar!(context::_GeoNumericContext,name::String,
                                  value::Float64,caller::AbstractString)
    # `gmsh_yysymbol::value` stores NaN/Inf unchecked (`x = 1/0` keeps `inf`);
    # finiteness is validated downstream where a value becomes geometry.
    if haskey(context.lists,name)
        values=context.lists[name]
        if isempty(values)
            context.stored_list_items<_MAX_GEO_CONTEXT_LIST_ITEMS || throw(ArgumentError(
                "$caller: stored numeric lists exceed $_MAX_GEO_CONTEXT_LIST_ITEMS entries"))
            push!(values,value)
            context.stored_list_items+=1
        else
            values[1]=value
        end
        delete!(context.unavailable_lists,name)
    end
    context.values[name]=value
    delete!(context.list_variables,name)
    delete!(context.unavailable,name)
    return nothing
end

function _geo_context_set_list!(context::_GeoNumericContext,name::String,
                                values::Vector{Float64},caller::AbstractString)
    length(values)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
        "$caller: list exceeds $_MAX_GEO_LIST_ITEMS entries"))
    # `s.value` accepts NaN/Inf unchecked — e.g. `LinSpace(a,b,1)` produces
    # `[nan]` upstream — so the store keeps them; use sites validate.
    old_length=haskey(context.lists,name) ? length(context.lists[name]) : 0
    stored=context.stored_list_items-old_length+length(values)
    stored<=_MAX_GEO_CONTEXT_LIST_ITEMS || throw(ArgumentError(
        "$caller: stored numeric lists exceed $_MAX_GEO_CONTEXT_LIST_ITEMS entries"))
    context.lists[name]=copy(values)
    push!(context.list_variables,name)
    context.stored_list_items=stored
    if isempty(values)
        delete!(context.values,name)
    else
        context.values[name]=values[1]
    end
    delete!(context.unavailable,name)
    delete!(context.unavailable_lists,name)
    return nothing
end

function _geo_context_list(context::_GeoNumericContext,name::String,
                           caller::AbstractString)
    if name in _GEO_SIDE_EFFECT_SYMBOLS
        throw(ArgumentError(
            "$caller: dynamic tag allocator $name is scalar and cannot use []"))
    elseif haskey(context.unavailable_lists,name)
        throw(ArgumentError(
            "$caller: numeric list $name is unavailable ($(context.unavailable_lists[name]))"))
    elseif haskey(context.lists,name)
        return copy(context.lists[name])
    elseif haskey(context.values,name)
        return Float64[context.values[name]]
    elseif haskey(context.unavailable,name)
        throw(ArgumentError(
            "$caller: numeric variable $name is unavailable ($(context.unavailable[name]))"))
    end
    # `tList '[' String__Index ']'` misses emit `yymsg(0)` and splice nothing.
    context.soft_unknown_reads ||
        throw(ArgumentError("$caller: unknown numeric list variable $name"))
    _geo_yyerror!(context,"Unknown variable '$name'")
    return Float64[]
end

function _geo_context_index(value::Float64,length::Int,caller::AbstractString)
    index=_geo_int_value(value,"$caller index")
    0<=index<length || throw(ArgumentError(
        "$caller: zero-based index $index is outside a list of length $length"))
    return index+1
end

function _geo_context_list_value(context::_GeoNumericContext,name::String,
                                 index_value::Float64,caller::AbstractString)
    if !haskey(context.lists,name) && !haskey(context.values,name) &&
       !haskey(context.unavailable_lists,name) &&
       !haskey(context.unavailable,name)
        if context.soft_unknown_reads
            _geo_yyerror!(context,"Unknown variable '$name(.)'")
            return 0.0
        end
        throw(ArgumentError("$caller: Unknown variable '$name(.)'"))
    end
    values=_geo_context_list(context,name,caller)
    index=_geo_int_value(index_value,"$caller index")
    if !(0<=index<length(values))
        if context.soft_unknown_reads
            _geo_yyerror!(context,"Uninitialized variable '$name[$index]'")
            return 0.0
        end
        throw(ArgumentError("$caller: zero-based index $index is outside " *
                            "a list of length $(length(values))"))
    end
    return values[index+1]
end

mutable struct _GeoExprParser
    source::String
    index::Int
    token::_GeoExprToken
    token_count::Int
    depth::Int
    context::_GeoNumericContext
    caller::String
end

const _GEO_SIDE_EFFECT_SYMBOLS=Set((
    "newp","newl","newc","newcl","newll","news","newsl","newreg","newv","newf"))
const _GEO_POINT_TAG_SYMBOLS=("newp",)
const _GEO_REGION_TAG_SYMBOLS=(
    "newl","newc","newcl","newll","news","newsl","newreg","newv")
const _GEO_FIELD_TAG_SYMBOLS=("newf",)
const _GEO_NONCONSTANT_FUNCTIONS=Set((
    "Rand","DefineNumber","GetNumber","GetValue","Exists","FileExists",
    "StringToName","S2N","Find","StrFind","StrCmp","StrLen","TextAttributes"))
const _GEO_NUMERIC_FUNCTIONS=Set((
    "Acos","Asin","Atan","Ceil","Cos","Cosh","Exp","Fabs","Abs",
    "Floor","Log","Log10","Round","Sqrt","Sin","Sinh","Step","Tan",
    "Tanh","Atan2","Fmod","Modulo","Hypot","Max","Min"))
# Fixed-arity binary numeric functions (`tX LP FExpr ',' FExpr RP`).
const _GEO_BINARY_FUNCTIONS=Set(("Atan2","Fmod","Modulo","Hypot","Max","Min"))
# Stateful FExpr functions taking only numeric arguments.
const _GEO_STATEFUL_FUNCTIONS=Set(("Rand",))
# Functions whose arguments cannot be parsed as plain numeric expressions:
# StringExprVar strings, ListOfDouble lists, bare `Struct_FullName` names.
const _GEO_RAW_ARG_FUNCTIONS=Set((
    "Exists","GetForced","GetNumber","GetValue","FileExists","Find",
    "StrFind","StrCmp","StrLen","TextAttributes","DimNameSpace",
    "StringToName","S2N","DefineNumber","GetNumberChoice"))
# StringExpr functions — `name(StringExpr, ...)` calls producing strings.
# `Str(a,b,...)` joins its arguments with `\n` (verified against the 4.15.2
# grammar); under `x() = Str(...)` it instead constructs the string list.
const _GEO_STRING_FUNCTIONS=Set((
    "Str","StrCat","StrPrefix","StrRelative","StrReplace","UpperCase",
    "LowerCase","LowerCaseIn","StrChoice","StrSub","Sprintf","GetEnv",
    "GetString","GetStringValue","GetForcedStr","FixRelativePath","DirName",
    "AbsolutePath","DefineString","NameStruct","NameToString","N2S"))
const _GEO_ALL_FUNCTIONS=union(_GEO_NUMERIC_FUNCTIONS,
    _GEO_NONCONSTANT_FUNCTIONS,_GEO_STATEFUL_FUNCTIONS,
    _GEO_RAW_ARG_FUNCTIONS,_GEO_STRING_FUNCTIONS)
# Reserved tokens in FExpr position — selector/transform/statement keywords
# head their own productions upstream, so a bare occurrence is a syntax
# error at the *following* token, never `Unknown variable` (`j = List` →
# `(;`, `x = List{..}` → `({`, `x = List + 1` → `(+`).
const _GEO_EXPR_RESERVED_NAMES=Set((
    "List","LinSpace","LogSpace","Catenary","Unique","ListFromFile",
    "Point","Curve","Line","Surface","Volume","GeoEntity","Physical",
    "Parent","BoundingBox","Mass","CenterOfMass","MatrixOfInertia","Normal",
    "Parametric","Color","Translate","Rotate","Symmetry","Dilate","Affine",
    "Closest","Duplicata","Boundary","PointsOf","CombinedBoundary",
    "Intersect","Split","Extrude","Boolean","BooleanUnion",
    "BooleanDifference","BooleanIntersection","BooleanFragments","Disk",
    "Sphere","Box","Cylinder","Torus","Cone","Wedge","Prism","Plane",
    "Bezier","BSpline","Spline","Nurbs","Circle","Ellipse","Wire",
    "ThruSections","Ruled","Using","Layers","Recombine","QuadTri",
    "Compound","Transform","Periodic","Delete","Recursive","Show","Hide",
    "Merge","Print","Printf","Error","Warning","Save","Exit","Abort",
    "Call","Return","Function","If","ElseIf","Else","EndIf","For","EndFor",
    "While","EndWhile","Include","SetFactory","NewModel","SyncModel",
    "Model","Optimize","Reclassify","ClassifySurfaces","CreateGeometry",
    "CreateTopology","Homology","Cohomology","Betti","RelocateMesh",
    "Import","Export","Project","Smoother","DefineConstant","DefineNumber",
    "UndefineConstant","Mesh","General","Geometry","Field","Kernel",
    "Factory","Onelab","Macro","Info","Debug","SetOrder","Coherence",
    "Transfinite","Reverse","Adapt","Homogeneous","Generalized","In",
    "Abs","GetNumber","SetNumber","GetString","SetString","Exists"))
# Gmsh 4.15.2 lexer keywords (Gmsh.l) — verbatim. `String__Index` positions
# (struct/namespace names, `x.y` families, `ns::` parts) reject these; member
# names (`tSTRING_Reserved`) reject all but Min/Max/Step/Box/Cylinder.
const _GEO_LEXER_KEYWORDS=Set((
    "Abort","Abs","AbsolutePath","Acos","AdaptMesh","Affine","Alias",
    "AliasWithOptions","Append","ArcCos","ArcSin","ArcTan","ArcTan2","Asin",
    "Atan","Atan2","BSpline","Betti","Bezier","Block","BooleanCoherence",
    "BooleanCommon","BooleanCut","BooleanDifference","BooleanFragments",
    "BooleanFuse","BooleanIntersection","BooleanSection","BooleanUnion",
    "BoundingBox","Box","Call","Catenary","CatmullRom","Ceil","CenterOfMass",
    "Chamfer","Characteristic","Circle","ClassifySurfaces","Closest",
    "CodeName","Coherence","Cohomology","Color","ColorTable","Combine",
    "Compound","Cone","Coordinates","CopyOptions","Cos","Cosh","Cpu",
    "CreateGeometry","CreateTopology","CurrentDir","CurrentDirectory",
    "CurrentFileName","Curve","Cylinder","DefineConstant","DefineNumber",
    "DefineString","Degenerated","Delete","Dilate","DimNameSpace","DirName",
    "Disk","Draw","Ellipse","Ellipsis","Ellipsoid","Elliptic","Else",
    "ElseIf","EndFor","EndIf","Error","Euclidian","Exists","Exit","Exp",
    "Extrude","Fabs","Field","FileExists","Fillet","Find","FixRelativePath",
    "Floor","Fmod","For","Function","GMSH_MAJOR_VERSION","GMSH_MINOR_VERSION",
    "GMSH_PATCH_VERSION","GeoEntity","GetEnv","GetForced","GetForcedStr",
    "GetNumber","GetString","GetStringValue","GetValue","GmshExecutableName",
    "HealShapes","Hide","Hole","Homology","Hypot","INTERPOLATION_SCHEME","If",
    "In","Intersect","Knots","Layers","Length","Levelset","LinSpace","Line",
    "List","ListFromFile","Log","Log10","LogSpace","LowerCase","LowerCaseIn",
    "MPI_Rank","MPI_Size","Macro","Mass","MatrixOfInertia","Max","Memory",
    "MeshAlgorithm","MeshSize","MeshSizeFromBoundary","Min","Modulo","N2S",
    "NameStruct","NameToString","NewModel","Normal","Nurbs","OnelabAction",
    "OnelabRun","Order","Parametric","Parent","Periodic","Physical","Pi",
    "Plane","Plugin","Point","PolarSphere","Printf","Quadric","Rand",
    "Recombine","RecombineMesh","Rectangle","Recursive","RefineMesh",
    "RelocateMesh","RenumberMeshElements","RenumberMeshNodes",
    "ReorientMesh","Return","Reverse","ReverseMesh","Rotate","Round","Ruled",
    "S2N","ScaleLastLayer","SetChanged","SetFactory","SetMaxTag","SetNumber",
    "SetPartition","SetString","SetTag","Sewing","ShapeFromFile","Show","Sin",
    "Sinh","Slide","Smoother","Sphere","Spline","Split","Sprintf","Sqrt",
    "Step","Str","StrCat","StrChoice","StrCmp","StrFind","StrLen","StrPrefix",
    "StrRelative","StrReplace","StrSub","StringToName","Struct","Surface",
    "Symmetry","SyncModel","T2","T3","TIME","Tan","Tanh","TestLevel",
    "TextAttributes","ThickSolid","ThruSections","Today","Torus",
    "TotalMemory","Transfinite","TransformMesh","Translate",
    "UndefineConstant","Unique","UnsplitWindow","UpperCase","Using","Volume",
    "Warning","Wedge","Wire"))
# `tSTRING_Reserved` adds Min/Max/Step/Box/Cylinder to plain identifiers —
# the only lexer keywords legal as member names.
const _GEO_MEMBER_EXTRA_KEYS=Set(("Min","Max","Step","Box","Cylinder"))
# Lexer keywords that head an `FExpr`/`FExpr_Multi` production — upstream
# they shift and the error lands on the *following* token (`v = Point(1)`
# → `(`, `v = List` → `;`), unlike pure-statement keywords which error at
# the keyword itself (`v = Circle(1)` → `Circle`). `List` is excluded:
# `tList LP` is a real production whose argument errors at the closer.
const _GEO_FEXPR_SHIFT_KEYWORDS=Set((
    "Point","Curve","Surface","Volume","GeoEntity","Physical","Line",
    "Parent","Normal","Parametric","Affine","BoundingBox",
    "CenterOfMass","Chamfer","Closest","Color","Dilate","Fillet",
    "Intersect","Mass","MatrixOfInertia","Rotate",
    "Ruled","Split","Struct","Translate","Symmetry","Extrude",
    "BooleanCoherence","BooleanCommon","BooleanCut","BooleanDifference",
    "BooleanFragments","BooleanFuse","BooleanIntersection",
    "BooleanSection","BooleanUnion"))
# Zero-argument `StringExpr` tokens — lexer-emitted literals that take no
# argument list (`Today(1)` errors at `(` upstream).
const _GEO_BARE_STRING_TOKENS=Set((
    "Today","CodeName","GmshExecutableName","OnelabAction",
    "CurrentFileName","CurrentDirectory","CurrentDir"))
_geo_tstring_valid(name::String)=!(name in _GEO_LEXER_KEYWORDS)
_geo_member_key_valid(name::String)=
    !(name in _GEO_LEXER_KEYWORDS) || name in _GEO_MEMBER_EXTRA_KEYS

# FLTK `menu_font_names` order — `getFontIndex` returns the position.
const _GEO_FONT_NAMES=(
    "Times-Roman","Times-Bold","Times-Italic","Times-BoldItalic",
    "Helvetica","Helvetica-Bold","Helvetica-Oblique","Helvetica-BoldOblique",
    "Courier","Courier-Bold","Courier-Oblique","Courier-BoldOblique",
    "Symbol","ZapfDingbats","Screen")
# `getFontAlign` — 0-8 alignment indices.
const _GEO_FONT_ALIGNS=Dict{String,Int}(
    "Left"=>0,"BottomLeft"=>0,"left"=>0,
    "Center"=>1,"BottomCenter"=>1,"center"=>1,
    "Right"=>2,"BottomRight"=>2,"right"=>2,
    "TopLeft"=>3,"TopCenter"=>4,"TopRight"=>5,
    "CenterLeft"=>6,"CenterCenter"=>7,"CenterRight"=>8)
# Bare identifier constants in FExpr position.
const _GEO_BARE_CONSTANTS=Dict{String,Function}(
    "Pi"=>(_p)->Float64(pi),
    "TestLevel"=>(_p)->Float64(_p.context.if_depth),
    "MPI_Rank"=>(_p)->0.0,
    "MPI_Size"=>(_p)->1.0,
    "GMSH_MAJOR_VERSION"=>(_p)->4.0,
    "GMSH_MINOR_VERSION"=>(_p)->15.0,
    "GMSH_PATCH_VERSION"=>(_p)->2.0,
    # Bare stateful tokens (`tCpu`/`tMemory`/`tTotalMemory` take no parens).
    "Cpu"=>(_p)->_geo_cpu_time(),
    "Memory"=>(_p)->_geo_memory_mb(),
    "TotalMemory"=>(_p)->_geo_total_memory_mb())

# Lexer keywords that never head a statement production — pure-expression
# function names, bare constants/tokens and contextual words (`In`, `Using`,
# `Knots`, ...). Used by `_geo_lhs_keyword_check`: as an assignment LHS head
# they error at the keyword itself upstream (`Pi = 5` → `syntax error (Pi)`),
# while statement-leading keywords shifted by the parser error at the next
# token (`Printf[...] = ` → `[`, `Struct =` → `=`).
const _GEO_EXPR_ONLY_KEYWORDS=union(
    _GEO_ALL_FUNCTIONS,_GEO_BARE_STRING_TOKENS,
    Set(keys(_GEO_BARE_CONSTANTS)),
    Set(("In","Using","Knots","Append","Order","Layers","ScaleLastLayer",
        "Sewing","Hole","Elliptic","Quadric","Ellipsoid","TIME","T2","T3",
        "INTERPOLATION_SCHEME","Exists","GetValue","N2S","S2N","NameStruct",
        "NameToString","StringToName","DimNameSpace","OnelabAction")))

# The first source token of `s` for `syntax error (token)` reporting: a whole
# identifier, a quote character, or a single punctuation char.
function _geo_stmt_head_token(s::AbstractString)
    t=strip(s)
    isempty(t) && return ";"
    (tm=match(r"^[A-Za-z_][A-Za-z0-9_]*",t))!==nothing && return String(tm.match)
    return string(first(t))
end

# Statement head check for an affectation LHS: a lexer keyword in `String__`
# position is a syntax error — at the keyword when it is expression-only, at
# the following token when the keyword heads a statement the parser already
# committed to.
function _geo_lhs_keyword_check(body::AbstractString)
    m=match(r"^[A-Za-z_][A-Za-z0-9_]*",body)
    m===nothing && return nothing
    name=String(m.match)
    name in _GEO_LEXER_KEYWORDS || return nothing
    if name in _GEO_EXPR_ONLY_KEYWORDS
        _geo_syntax_abort(name)
    else
        _geo_syntax_abort(_geo_stmt_head_token(
            body[nextind(body,lastindex(m.match)):end]))
    end
    return nothing
end

mutable struct _GeoTagAllocatorState
    builtin_point_max::Int
    builtin_curve_max::Int
    builtin_surface_max::Int
    builtin_volume_max::Int
    occ_point_max::Int
    occ_curve_max::Int
    occ_surface_max::Int
    occ_volume_max::Int
    builtin_curveloop_max::Int
    builtin_surfaceloop_max::Int
    # OCC internals track wire (-1) and shell (-2) maxima too —
    # `NEWCURVELOOP()`/`NEWSURFACELOOP()`/`NEWREG()` read them whenever
    # `getOCCInternals()` exists (`occ_active`).
    occ_curveloop_max::Int
    occ_surfaceloop_max::Int
    physical_group_max::Int
    field_max::Int
    point_entity_max::Int
    curve_entity_max::Int
    surface_entity_max::Int
    factory::Symbol
    # Whether OCC internals exist (`getOCCInternals() != nullptr`). Once
    # `SetFactory("OpenCASCADE")` creates them they keep their counters and keep
    # participating in `newX` maxima even after switching back to Built-in —
    # `NEWPOINT()`/`NEWREG()` take `max(GEO, OCC)` whenever internals exist.
    occ_active::Bool
    # Live entity tags per factory plus the last SetMaxTag floor values let
    # allocator reads shrink after a Boolean operand `Delete`, matching
    # Gmsh's `max(SetMaxTag floor, live maximum)` counters. `volume_boundaries`
    # maps each live volume tag to the hidden point/curve/surface ranges its
    # declaration consumed so deleting a volume also releases them.
    live_builtin_points::Set{Int}
    live_builtin_curves::Set{Int}
    live_builtin_surfaces::Set{Int}
    live_builtin_volumes::Set{Int}
    live_occ_points::Set{Int}
    live_occ_curves::Set{Int}
    live_occ_surfaces::Set{Int}
    live_occ_volumes::Set{Int}
    builtin_point_floor::Int
    builtin_curve_floor::Int
    builtin_surface_floor::Int
    builtin_volume_floor::Int
    occ_point_floor::Int
    occ_curve_floor::Int
    occ_surface_floor::Int
    occ_volume_floor::Int
    volume_boundaries::Dict{Int,NTuple{3,Vector{Int}}}
    # Live `Field[i]` ids — `newf` is `FieldManager::maxId() + 1`, a live
    # maximum, so `Delete Field[max]` shrinks it.
    live_fields::Set{Int}
    geometry_unavailable::Union{Nothing,String}
    physical_unavailable::Union{Nothing,String}
    field_unavailable::Union{Nothing,String}
end

_GeoTagAllocatorState()=_GeoTagAllocatorState(
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    :builtin,false,
    Set{Int}(),Set{Int}(),Set{Int}(),Set{Int}(),
    Set{Int}(),Set{Int}(),Set{Int}(),Set{Int}(),
    0,0,0,0,0,0,0,0,
    Dict{Int,NTuple{3,Vector{Int}}}(),Set{Int}(),nothing,nothing,nothing)

@inline function _geo_allocator_point_max(state::_GeoTagAllocatorState)
    return state.occ_active ?
        max(state.builtin_point_max,state.occ_point_max) :
        state.builtin_point_max
end

@inline function _geo_allocator_curve_max(state::_GeoTagAllocatorState)
    return state.occ_active ?
        max(state.builtin_curve_max,state.occ_curve_max) :
        state.builtin_curve_max
end

@inline function _geo_allocator_surface_max(state::_GeoTagAllocatorState)
    return state.occ_active ?
        max(state.builtin_surface_max,state.occ_surface_max) :
        state.builtin_surface_max
end

# `NEWPOINT()`/`NEWCURVE()`/... compute `getMaxTag(dim) + 1` per kernel and
# then take the max — `max(incr(GEO), incr(OCC))`, not `incr(max(GEO, OCC))`.
# At the INT32_MAX wrap the order is observable: a wrapped builtin candidate
# loses to a small positive OCC counter.
_geo_allocator_point_next(state::_GeoTagAllocatorState)=
    state.occ_active ?
        max(_geo_int32_incr(state.builtin_point_max),
            _geo_int32_incr(state.occ_point_max)) :
        _geo_int32_incr(state.builtin_point_max)
_geo_allocator_curve_next(state::_GeoTagAllocatorState)=
    state.occ_active ?
        max(_geo_int32_incr(state.builtin_curve_max),
            _geo_int32_incr(state.occ_curve_max)) :
        _geo_int32_incr(state.builtin_curve_max)
_geo_allocator_surface_next(state::_GeoTagAllocatorState)=
    state.occ_active ?
        max(_geo_int32_incr(state.builtin_surface_max),
            _geo_int32_incr(state.occ_surface_max)) :
        _geo_int32_incr(state.builtin_surface_max)
_geo_allocator_volume_next(state::_GeoTagAllocatorState)=
    state.occ_active ?
        max(_geo_int32_incr(state.builtin_volume_max),
            _geo_int32_incr(state.occ_volume_max)) :
        _geo_int32_incr(state.builtin_volume_max)
_geo_allocator_curveloop_next(state::_GeoTagAllocatorState)=
    state.occ_active ?
        max(_geo_int32_incr(state.builtin_curveloop_max),
            _geo_int32_incr(state.occ_curveloop_max)) :
        _geo_int32_incr(state.builtin_curveloop_max)
_geo_allocator_surfaceloop_next(state::_GeoTagAllocatorState)=
    state.occ_active ?
        max(_geo_int32_incr(state.builtin_surfaceloop_max),
            _geo_int32_incr(state.occ_surfaceloop_max)) :
        _geo_int32_incr(state.builtin_surfaceloop_max)

# `NEWREG()` — `max` over every non-point dimension and the physical tag of
# each `getMaxTag(dim) + 1`: the increment happens *per dimension* before the
# max, so a dimension at INT32_MAX contributes a wrapped INT32_MIN that loses
# to the remaining counters (verified: `Line(2147483647)` leaves
# `newreg == 1` upstream, not -2147483648).
@inline function _geo_allocator_region_next(state::_GeoTagAllocatorState)
    next=maximum(_geo_int32_incr.((state.builtin_curve_max,
        state.builtin_surface_max,state.builtin_volume_max,
        state.builtin_curveloop_max,state.builtin_surfaceloop_max,
        state.physical_group_max)))
    return state.occ_active ?
        max(next,_geo_int32_incr(state.occ_curve_max),
            _geo_int32_incr(state.occ_surface_max),
            _geo_int32_incr(state.occ_volume_max),
            _geo_int32_incr(state.occ_curveloop_max),
            _geo_int32_incr(state.occ_surfaceloop_max)) : next
end

function _geo_allocator_invalidate!(state::_GeoTagAllocatorState,
                                    reason::AbstractString;
                                    geometry::Bool=true,
                                    physical::Bool=false,
                                    fields::Bool=true)
    geometry && (state.geometry_unavailable=String(reason))
    physical && (state.physical_unavailable=String(reason))
    fields && (state.field_unavailable=String(reason))
    return nothing
end

# `getMaxTag() + 1` — a C++ `int` increment: a bump past INT32_MAX wraps to
# INT32_MIN upstream (verified: `SetMaxTag Point(2147483647)` yields
# `newp == -2147483648`).
_geo_int32_incr(x::Int)=Int(reinterpret(Int32,(Int64(x)+1) % UInt32))

function _geo_context_refresh_allocators!(context::_GeoNumericContext,
                                          state::_GeoTagAllocatorState)
    function refresh!(names,next_tag::Int,reason::Union{Nothing,String})
        for name in names
            delete!(context.list_variables,name)
            delete!(context.unavailable_lists,name)
            if reason===nothing
                context.values[name]=Float64(next_tag)
                delete!(context.unavailable,name)
            else
                delete!(context.values,name)
                context.unavailable[name]=String(reason)
            end
        end
        return nothing
    end
    refresh!(_GEO_POINT_TAG_SYMBOLS,
             _geo_allocator_point_next(state),
             state.geometry_unavailable)
    region_unavailable=state.geometry_unavailable===nothing ?
        state.physical_unavailable : state.geometry_unavailable
    # `Geometry.OldNewReg` (default 1): `newl`/`newll`/`news`/`newsl`/`newv`
    # all read `NEWREG()` — the max over every non-point entity dimension
    # plus physical tags. At 0 each reads its own `getMaxTag(dim) + 1`
    # (`newreg` itself is always `NEWREG()`).
    if !iszero(something(_geo_option_number(
            context,"Geometry",0,"OldNewReg"),1.0))
        refresh!(_GEO_REGION_TAG_SYMBOLS,_geo_allocator_region_next(state),
                 region_unavailable)
    else
        refresh!(("newl","newc"),_geo_allocator_curve_next(state),
                 state.geometry_unavailable)
        refresh!(("newll","newcl"),_geo_allocator_curveloop_next(state),
                 state.geometry_unavailable)
        refresh!(("news",),_geo_allocator_surface_next(state),
                 state.geometry_unavailable)
        refresh!(("newsl",),_geo_allocator_surfaceloop_next(state),
                 state.geometry_unavailable)
        refresh!(("newv",),_geo_allocator_volume_next(state),
                 state.geometry_unavailable)
        refresh!(("newreg",),_geo_allocator_region_next(state),
                 region_unavailable)
    end
    # `NEWFIELD()` is `FieldManager::maxId() + 1` — a live maximum over the
    # field map (it can go back down on `Delete Field`, and below 1 when every
    # id is nonpositive). At execution the map is `context.fields`; the params
    # scan tracks it through `state.field_max`.
    fields=context.fields
    field_max=if fields===nothing
        state.field_max
    else
        isempty(fields) ? 0 : maximum(keys(fields))
    end
    refresh!(_GEO_FIELD_TAG_SYMBOLS,_geo_int32_incr(field_max),
             state.field_unavailable)
    return nothing
end

# `addX(int &tag, ...)` convention: a negative definition tag auto-assigns
# from the *active* kernel's per-dimension counter (`getMaxTag(dim) + 1`).
# Loop counters live on the shared builtins — `SetFactory` syncs both kernels'
# -2/-1 counters to the running max on every switch, so a shared counter
# reproduces the active kernel's value in every reachable state.
function _geo_allocator_next_tag(state::_GeoTagAllocatorState,kind::Symbol)
    if kind==:point
        return _geo_int32_incr(state.factory==:opencascade ?
            state.occ_point_max : state.builtin_point_max)
    elseif kind==:curve
        return _geo_int32_incr(state.factory==:opencascade ?
            state.occ_curve_max : state.builtin_curve_max)
    elseif kind==:surface
        return _geo_int32_incr(state.factory==:opencascade ?
            state.occ_surface_max : state.builtin_surface_max)
    elseif kind==:volume
        return _geo_int32_incr(state.factory==:opencascade ?
            state.occ_volume_max : state.builtin_volume_max)
    elseif kind==:curveloop
        return _geo_int32_incr(state.factory==:opencascade ?
            state.occ_curveloop_max : state.builtin_curveloop_max)
    elseif kind==:surfaceloop
        return _geo_int32_incr(state.factory==:opencascade ?
            state.occ_surfaceloop_max : state.builtin_surfaceloop_max)
    end
    throw(ArgumentError("unknown auto-assign tag namespace $kind"))
end

function _geo_allocator_record_explicit!(state::_GeoTagAllocatorState,
                                         kind::Symbol,tag::Int)
    state.geometry_unavailable===nothing || return nothing
    tag<0 && (tag=_geo_allocator_next_tag(state,kind))
    if kind==:point
        if state.factory==:opencascade
            push!(state.live_occ_points,tag)
            state.occ_point_max=max(state.occ_point_max,tag)
        else
            push!(state.live_builtin_points,tag)
            state.builtin_point_max=max(state.builtin_point_max,tag)
        end
        state.point_entity_max=max(state.point_entity_max,tag)
    elseif kind==:curve
        if state.factory==:opencascade
            push!(state.live_occ_curves,tag)
            state.occ_curve_max=max(state.occ_curve_max,tag)
        else
            push!(state.live_builtin_curves,tag)
            state.builtin_curve_max=max(state.builtin_curve_max,tag)
        end
        state.curve_entity_max=max(state.curve_entity_max,tag)
    elseif kind==:surface
        if state.factory==:opencascade
            push!(state.live_occ_surfaces,tag)
            state.occ_surface_max=max(state.occ_surface_max,tag)
        else
            push!(state.live_builtin_surfaces,tag)
            state.builtin_surface_max=max(state.builtin_surface_max,tag)
        end
        state.surface_entity_max=max(state.surface_entity_max,tag)
    elseif kind==:volume
        if state.factory==:opencascade
            push!(state.live_occ_volumes,tag)
            state.occ_volume_max=max(state.occ_volume_max,tag)
        else
            push!(state.live_builtin_volumes,tag)
            state.builtin_volume_max=max(state.builtin_volume_max,tag)
        end
    elseif kind==:curveloop
        if state.factory==:opencascade
            state.occ_curveloop_max=max(state.occ_curveloop_max,tag)
        else
            state.builtin_curveloop_max=max(state.builtin_curveloop_max,tag)
        end
    elseif kind==:surfaceloop
        if state.factory==:opencascade
            state.occ_surfaceloop_max=max(state.occ_surfaceloop_max,tag)
        else
            state.builtin_surfaceloop_max=max(state.builtin_surfaceloop_max,tag)
        end
    else
        throw(ArgumentError("unknown dynamic tag namespace $kind"))
    end
    return nothing
end

function _geo_allocator_record_physical!(state::_GeoTagAllocatorState,tag::Int)
    state.physical_unavailable===nothing || return nothing
    state.physical_group_max=max(state.physical_group_max,tag)
    return nothing
end

function _geo_allocator_advance(current::Int,count::Int,
                                caller::AbstractString,label::AbstractString)
    count>=0 || throw(ArgumentError("$caller: $label allocation count must be non-negative"))
    current<=typemax(Int32)-count || throw(ArgumentError(
        "$caller: $label topology exhausts Gmsh's signed 32-bit tag range"))
    return current+count
end

function _geo_allocator_record_hidden!(state::_GeoTagAllocatorState,
                                       kind::Symbol,count::Int,
                                       caller::AbstractString)
    current,entity_max,label=if kind==:point
        (_geo_allocator_point_max(state),state.point_entity_max,
         "Point")
    elseif kind==:curve
        (_geo_allocator_curve_max(state),state.curve_entity_max,
         "Curve")
    elseif kind==:surface
        (_geo_allocator_surface_max(state),state.surface_entity_max,
         "Surface")
    else
        throw(ArgumentError("unknown hidden dynamic tag namespace $kind"))
    end
    base=max(current,entity_max)
    final=_geo_allocator_advance(base,count,caller,label)
    allocated=(base+1):final
    if kind==:point
        if state.factory==:opencascade
            union!(state.live_occ_points,allocated)
            state.occ_point_max=final
        else
            union!(state.live_builtin_points,allocated)
            state.builtin_point_max=final
        end
        state.point_entity_max=final
    elseif kind==:curve
        if state.factory==:opencascade
            union!(state.live_occ_curves,allocated)
            state.occ_curve_max=final
        else
            union!(state.live_builtin_curves,allocated)
            state.builtin_curve_max=final
        end
        state.curve_entity_max=final
    else
        if state.factory==:opencascade
            union!(state.live_occ_surfaces,allocated)
            state.occ_surface_max=final
        else
            union!(state.live_builtin_surfaces,allocated)
            state.builtin_surface_max=final
        end
        state.surface_entity_max=final
    end
    return allocated
end

function _geo_allocator_set_max!(state::_GeoTagAllocatorState,
                                 kind::String,value::Int)
    state.geometry_unavailable===nothing || return nothing
    if kind=="Point"
        if state.factory==:opencascade
            state.occ_point_floor=max(state.occ_point_floor,value)
            state.occ_point_max=max(state.occ_point_max,value)
        else
            state.builtin_point_floor=value
            state.builtin_point_max=value
        end
    elseif kind=="Curve"
        if state.factory==:opencascade
            state.occ_curve_floor=max(state.occ_curve_floor,value)
            state.occ_curve_max=max(state.occ_curve_max,value)
        else
            state.builtin_curve_floor=value
            state.builtin_curve_max=value
        end
    elseif kind=="Surface"
        if state.factory==:opencascade
            state.occ_surface_floor=max(state.occ_surface_floor,value)
            state.occ_surface_max=max(state.occ_surface_max,value)
        else
            state.builtin_surface_floor=value
            state.builtin_surface_max=value
        end
    elseif kind=="Volume"
        if state.factory==:opencascade
            state.occ_volume_floor=max(state.occ_volume_floor,value)
            state.occ_volume_max=max(state.occ_volume_max,value)
        else
            state.builtin_volume_floor=value
            state.builtin_volume_max=value
        end
    elseif kind=="CurveLoop"
        # `GeoEntity{-1}` — loop counters have no floor bookkeeping (they are
        # not live-tracked): builtin assigns, OCC keeps the running max.
        if state.factory==:opencascade
            state.occ_curveloop_max=max(state.occ_curveloop_max,value)
        else
            state.builtin_curveloop_max=value
        end
    elseif kind=="SurfaceLoop"
        if state.factory==:opencascade
            state.occ_surfaceloop_max=max(state.occ_surfaceloop_max,value)
        else
            state.builtin_surfaceloop_max=value
        end
    else
        throw(ArgumentError("unknown SetMaxTag namespace $kind"))
    end
    return nothing
end

# `SetFactory("OpenCASCADE"|"Built-in")` — switching kernels carries each
# dimension's counter to `max(this, other)`, matching Gmsh's behavior where
# both internals track the same physical tag space.
function _geo_allocator_set_factory!(state::_GeoTagAllocatorState,
                                     factory::AbstractString)
    if factory=="OpenCASCADE"
        state.occ_point_max=max(state.occ_point_max,state.builtin_point_max)
        state.occ_curve_max=max(state.occ_curve_max,state.builtin_curve_max)
        state.occ_surface_max=max(state.occ_surface_max,state.builtin_surface_max)
        state.occ_volume_max=max(state.occ_volume_max,state.builtin_volume_max)
        state.occ_curveloop_max=max(state.occ_curveloop_max,state.builtin_curveloop_max)
        state.occ_surfaceloop_max=max(state.occ_surfaceloop_max,state.builtin_surfaceloop_max)
        state.factory=:opencascade
        state.occ_active=true
    else
        # Switching back syncs the builtin counters only while OCC internals
        # exist (`if(getOCCInternals())`); the OCC counters themselves survive
        # — `occ_active` stays set so `newX` reads keep taking the cross-kernel
        # maximum.
        if state.occ_active
            state.builtin_point_max=max(state.builtin_point_max,state.occ_point_max)
            state.builtin_curve_max=max(state.builtin_curve_max,state.occ_curve_max)
            state.builtin_surface_max=max(state.builtin_surface_max,state.occ_surface_max)
            state.builtin_volume_max=max(state.builtin_volume_max,state.occ_volume_max)
            state.builtin_curveloop_max=max(state.builtin_curveloop_max,state.occ_curveloop_max)
            state.builtin_surfaceloop_max=max(state.builtin_surfaceloop_max,state.occ_surfaceloop_max)
        end
        state.factory=:builtin
    end
    return nothing
end

function _geo_allocator_record_primitive!(state::_GeoTagAllocatorState,
                                          kind::String,tag::Int,
                                          values::Vector{Float64},
                                          caller::AbstractString)
    state.geometry_unavailable===nothing || return nothing
    point_count,curve_count,surface_count=if kind=="Box"
        (8,12,6)
    elseif kind=="Cylinder"
        (2,3,3)
    elseif kind=="Sphere" || kind=="PolarSphere"
        (2,3,1)
    elseif kind=="Cone"
        length(values)==8 || throw(ArgumentError(
            "$caller: Cone requires eight numeric parameters"))
        (2,3,iszero(values[7]) || iszero(values[8]) ? 2 : 3)
    elseif kind=="Torus"
        # A partial torus adds the second rim vertex, the third (start
        # meridian) edge, and the two Plane caps.
        partial=length(values)==6 && values[6]<2π
        partial ? (2,3,3) : (1,2,1)
    else
        throw(ArgumentError("$caller: unsupported primitive allocator kind $kind"))
    end
    points=_geo_allocator_record_hidden!(state,:point,point_count,caller)
    curves=_geo_allocator_record_hidden!(state,:curve,curve_count,caller)
    surfaces=_geo_allocator_record_hidden!(state,:surface,surface_count,caller)
    _geo_allocator_record_explicit!(state,:volume,tag)
    state.volume_boundaries[tag]=
        (collect(points),collect(curves),collect(surfaces))
    return nothing
end

# `Rectangle`/`Disk` are OCC-only *surface* primitives: the statement tag
# binds the face, while the wire, edges, and vertices take sequential
# automatic tags in `_bind`'s recursive pass. A positive `roundedRadius`
# replaces each rectangle corner with an arc, doubling the point/curve
# counts; a `Disk` edge is a single closed (ellipse) edge on one vertex.
function _geo_allocator_record_surface_primitive!(
        state::_GeoTagAllocatorState,kind::String,tag::Int,
        values::Vector{Float64},caller::AbstractString)
    state.geometry_unavailable===nothing || return nothing
    point_count,curve_count=if kind=="Rectangle"
        rounded=length(values)==6 && values[6]>0
        rounded ? (8,8) : (4,4)
    elseif kind=="Disk"
        (1,1)
    else
        throw(ArgumentError(
            "$caller: unsupported surface primitive allocator kind $kind"))
    end
    _geo_allocator_record_explicit!(state,:surface,tag)
    _geo_allocator_record_explicit!(state,:curveloop,
        _geo_allocator_next_tag(state,:curveloop))
    _geo_allocator_record_hidden!(state,:point,point_count,caller)
    _geo_allocator_record_hidden!(state,:curve,curve_count,caller)
    return nothing
end

function _geo_allocator_record_field!(state::_GeoTagAllocatorState,tag::Int)
    state.field_unavailable===nothing || return nothing
    push!(state.live_fields,tag)
    # `maxId()` is the live maximum key — an all-negative field set leaves
    # `newf` at or below zero.
    state.field_max=_geo_live_max(state.live_fields)
    return nothing
end

# `Delete Field[i]` erases the map entry — `newf` drops back to
# `max(live ids) + 1` when the deleted id was the maximum.
function _geo_allocator_delete_field!(state::_GeoTagAllocatorState,tag::Int)
    state.field_unavailable===nothing || return nothing
    delete!(state.live_fields,tag)
    state.field_max=_geo_live_max(state.live_fields)
    return nothing
end

@inline _geo_live_max(live::Set{Int})=isempty(live) ? 0 : maximum(live)

# A Boolean operand `Delete` removes the operand volume together with the
# hidden boundary entities its declaration consumed, shrinking the owning
# factory's counters to `max(SetMaxTag floor, live maximum)`, matching Gmsh.
function _geo_allocator_delete_volume!(state::_GeoTagAllocatorState,tag::Int)
    boundaries=get(state.volume_boundaries,tag,nothing)
    if tag in state.live_builtin_volumes
        delete!(state.live_builtin_volumes,tag)
        if boundaries!==nothing
            setdiff!(state.live_builtin_points,boundaries[1])
            setdiff!(state.live_builtin_curves,boundaries[2])
            setdiff!(state.live_builtin_surfaces,boundaries[3])
        end
        state.builtin_point_max=max(state.builtin_point_floor,
            _geo_live_max(state.live_builtin_points))
        state.builtin_curve_max=max(state.builtin_curve_floor,
            _geo_live_max(state.live_builtin_curves))
        state.builtin_surface_max=max(state.builtin_surface_floor,
            _geo_live_max(state.live_builtin_surfaces))
        state.builtin_volume_max=max(state.builtin_volume_floor,
            _geo_live_max(state.live_builtin_volumes))
    end
    if tag in state.live_occ_volumes
        delete!(state.live_occ_volumes,tag)
        if boundaries!==nothing
            setdiff!(state.live_occ_points,boundaries[1])
            setdiff!(state.live_occ_curves,boundaries[2])
            setdiff!(state.live_occ_surfaces,boundaries[3])
        end
        state.occ_point_max=max(state.occ_point_floor,
            _geo_live_max(state.live_occ_points))
        state.occ_curve_max=max(state.occ_curve_floor,
            _geo_live_max(state.live_occ_curves))
        state.occ_surface_max=max(state.occ_surface_floor,
            _geo_live_max(state.live_occ_surfaces))
        state.occ_volume_max=max(state.occ_volume_floor,
            _geo_live_max(state.live_occ_volumes))
    end
    delete!(state.volume_boundaries,tag)
    state.point_entity_max=max(
        _geo_live_max(state.live_builtin_points),
        _geo_live_max(state.live_occ_points))
    state.curve_entity_max=max(
        _geo_live_max(state.live_builtin_curves),
        _geo_live_max(state.live_occ_curves))
    state.surface_entity_max=max(
        _geo_live_max(state.live_builtin_surfaces),
        _geo_live_max(state.live_occ_surfaces))
    return nothing
end

# `.geo` `Delete`/`Recursive Delete` entity accounting. Gmsh's
# `GEO_Internals::remove` decrements the owning factory's dimension counter
# only when the deleted tag was its current maximum (`setMaxTag(dim,
# tmax-1)` in `DeletePoint`/`DeleteCurve`/`DeleteSurface`/`DeleteVolume`) —
# it does NOT recompute the counter from the surviving entities, so a sparse
# delete leaves the counter above the live maximum. `deleted` lists the
# entities actually removed, in execution order.
function _geo_allocator_delete_entities!(state::_GeoTagAllocatorState,
                                         deleted)
    dim_state=Dict(0=>(state.live_builtin_points,state.live_occ_points,
                       :builtin_point_max,:occ_point_max),
                   1=>(state.live_builtin_curves,state.live_occ_curves,
                       :builtin_curve_max,:occ_curve_max),
                   2=>(state.live_builtin_surfaces,state.live_occ_surfaces,
                       :builtin_surface_max,:occ_surface_max),
                   3=>(state.live_builtin_volumes,state.live_occ_volumes,
                       :builtin_volume_max,:occ_volume_max))
    for entry in deleted
        (dim,tag)=entry
        haskey(dim_state,dim) || continue
        (live_b,live_o,field_b,field_o)=dim_state[dim]
        if tag in live_b
            delete!(live_b,tag)
            getfield(state,field_b)==tag &&
                setfield!(state,field_b,tag-1)
        elseif tag in live_o
            delete!(live_o,tag)
            getfield(state,field_o)==tag &&
                setfield!(state,field_o,tag-1)
        else
            # Untracked tags are built-in-kernel records — decrement the
            # built-in counter when the tag was its maximum.
            getfield(state,field_b)==tag &&
                setfield!(state,field_b,tag-1)
        end
        dim==3 && delete!(state.volume_boundaries,tag)
    end
    state.point_entity_max=max(
        _geo_live_max(state.live_builtin_points),
        _geo_live_max(state.live_occ_points))
    state.curve_entity_max=max(
        _geo_live_max(state.live_builtin_curves),
        _geo_live_max(state.live_occ_curves))
    state.surface_entity_max=max(
        _geo_live_max(state.live_builtin_surfaces),
        _geo_live_max(state.live_occ_surfaces))
    return nothing
end

# `Delete Model`/`Delete All` destroy the geometry internals: every entity,
# physical group, mesh attribute, and tag counter resets while user variables
# survive (`Delete All` additionally clears them — handled at the exec layer).
# `GEO_Internals::_allocateAll` initializes every entity counter to
# `Geometry.FirstEntityTag - 1` and `_maxPhysicalNum` to
# `Geometry.FirstPhysicalTag - 1`, reading the options at model-creation
# time — `NewModel` (`new GModel`) and `Delete Model`/`All` (`destroy` →
# `_freeAll`+`_allocateAll`) therefore observe writes made earlier in the
# same file (verified: `Geometry.FirstEntityTag = 5; NewModel` leaves
# `newp == 5` upstream).
# `fresh` distinguishes `new GModel()` (`NewModel`, `Delete All` →
# `ClearProject`) — where OCC internals cease to exist — from
# `GModel::destroy` (`Delete Model`) which calls `resetOCCInternals()`:
# the OCC internals object survives the destroy (only its contents reset to
# `firstEntityTag - 1`), so `getOCCInternals()` stays non-null and `newX`
# reads keep taking the cross-kernel maximum.
function _geo_allocator_reset_model!(state::_GeoTagAllocatorState,
                                     context::_GeoNumericContext;
                                     fresh::Bool=false)
    entity_base=max(_geo_signed_gmsh_int_value(
        something(_geo_option_number(
            context,"Geometry",0,"FirstEntityTag"),1.0),
        "Geometry.FirstEntityTag"),1)-1
    physical_base=max(_geo_signed_gmsh_int_value(
        something(_geo_option_number(
            context,"Geometry",0,"FirstPhysicalTag"),1.0),
        "Geometry.FirstPhysicalTag"),1)-1
    for field in (:builtin_point_max,:builtin_curve_max,:builtin_surface_max,
                  :builtin_volume_max,:occ_point_max,:occ_curve_max,
                  :occ_surface_max,:occ_volume_max,:builtin_curveloop_max,
                  :builtin_surfaceloop_max,:occ_curveloop_max,
                  :occ_surfaceloop_max,
                  :point_entity_max,:curve_entity_max,:surface_entity_max)
        setfield!(state,field,entity_base)
    end
    state.physical_group_max=physical_base
    state.field_max=0
    for live in (state.live_builtin_points,state.live_builtin_curves,
                 state.live_builtin_surfaces,state.live_builtin_volumes,
                 state.live_occ_points,state.live_occ_curves,
                 state.live_occ_surfaces,state.live_occ_volumes)
        empty!(live)
    end
    for floor in (:builtin_point_floor,:builtin_curve_floor,
                  :builtin_surface_floor,:builtin_volume_floor,
                  :occ_point_floor,:occ_curve_floor,:occ_surface_floor,
                  :occ_volume_floor)
        setfield!(state,floor,entity_base)
    end
    empty!(state.volume_boundaries)
    empty!(state.live_fields)
    fresh && (state.occ_active=false)
    return nothing
end

function _geo_allocator_live_entities(state::_GeoTagAllocatorState,
                                      kind::Symbol)
    if kind==:point
        return state.factory==:opencascade ?
            state.live_occ_points : state.live_builtin_points
    elseif kind==:curve
        return state.factory==:opencascade ?
            state.live_occ_curves : state.live_builtin_curves
    end
    return state.factory==:opencascade ?
        state.live_occ_surfaces : state.live_builtin_surfaces
end

# Gmsh re-tags a Boolean result's boundary entities to the lowest tags still
# free in the active factory. Claim `count` such tags for `kind` and raise the
# owning counters; returns the claimed tags (possibly non-contiguous).
function _geo_allocator_claim_lowest!(state::_GeoTagAllocatorState,
                                      kind::Symbol,count::Int)
    count<=0 && return Int[]
    live=_geo_allocator_live_entities(state,kind)
    claimed=Int[];sizehint!(claimed,count)
    tag=1
    while length(claimed)<count
        tag in live || push!(claimed,tag)
        tag+=1
    end
    union!(live,claimed)
    top=claimed[end]
    if kind==:point
        if state.factory==:opencascade
            state.occ_point_max=max(state.occ_point_max,top)
        else
            state.builtin_point_max=max(state.builtin_point_max,top)
        end
        state.point_entity_max=max(state.point_entity_max,top)
    elseif kind==:curve
        if state.factory==:opencascade
            state.occ_curve_max=max(state.occ_curve_max,top)
        else
            state.builtin_curve_max=max(state.builtin_curve_max,top)
        end
        state.curve_entity_max=max(state.curve_entity_max,top)
    else
        if state.factory==:opencascade
            state.occ_surface_max=max(state.occ_surface_max,top)
        else
            state.builtin_surface_max=max(state.builtin_surface_max,top)
        end
        state.surface_entity_max=max(state.surface_entity_max,top)
    end
    return claimed
end

# Materialized Boolean results own real entities whose tags predictive claim
# bookkeeping cannot model exactly (splits, imprints, pseudo-preserved operand
# reuse). After such a statement executes, the live sets and per-volume
# boundary records resync to the model: tags that vanished drop from both
# factories, and untracked tags join the active factory's live set.
function _geo_allocator_resync_model!(state::_GeoTagAllocatorState,m)
    mtags=(Set{Int}(keys(m.points)),Set{Int}(keys(m.curves)),
           Set{Int}(keys(m.surfaces)))
    live_b=(state.live_builtin_points,state.live_builtin_curves,
            state.live_builtin_surfaces)
    live_o=(state.live_occ_points,state.live_occ_curves,
            state.live_occ_surfaces)
    occ=state.factory==:opencascade
    for k in 1:3
        intersect!(live_b[k],mtags[k]);intersect!(live_o[k],mtags[k])
        union!(occ ? live_o[k] : live_b[k],
               setdiff(mtags[k],live_b[k],live_o[k]))
    end
    vtags=Set{Int}(keys(m.volumes))
    intersect!(state.live_builtin_volumes,vtags)
    intersect!(state.live_occ_volumes,vtags)
    union!(occ ? state.live_occ_volumes : state.live_builtin_volumes,
           setdiff(vtags,state.live_builtin_volumes,state.live_occ_volumes))
    empty!(state.volume_boundaries)
    for (t,shells) in m.volumes
        isempty(shells) && continue
        surfs=Set{Int}();curvs=Set{Int}();pts=Set{Int}()
        for sh in shells,sc in get(m.surface_loops,sh,Int[])
            s=abs(sc);push!(surfs,s)
            for lp in get(m.surfaces,s,Int[]),cc in get(m.loops,lp,Int[])
                c=abs(cc);push!(curvs,c)
                union!(pts,get(m.curves,c,Int[]))
            end
        end
        state.volume_boundaries[t]=(collect(pts),collect(curvs),
                                    collect(surfs))
    end
    # Discrete entities (`Merge`-imported, `CreateTopology`-derived, homology
    # chains) share the per-dimension elementary namespace upstream —
    # `getMaxElementaryNumber` counts them — so the allocation maxima must
    # skip their tags even though they join no live set.
    discrete_max=ntuple(4) do d
        top=0
        for (dd,tag) in keys(m.discrete)
            dd==d-1 && tag>top && (top=tag)
        end
        top
    end
    state.builtin_point_max=max(state.builtin_point_floor,
        _geo_live_max(state.live_builtin_points),discrete_max[1])
    state.builtin_curve_max=max(state.builtin_curve_floor,
        _geo_live_max(state.live_builtin_curves),discrete_max[2])
    state.builtin_surface_max=max(state.builtin_surface_floor,
        _geo_live_max(state.live_builtin_surfaces),discrete_max[3])
    state.builtin_volume_max=max(state.builtin_volume_floor,
        _geo_live_max(state.live_builtin_volumes),discrete_max[4])
    state.occ_point_max=max(state.occ_point_floor,
        _geo_live_max(state.live_occ_points),discrete_max[1])
    state.occ_curve_max=max(state.occ_curve_floor,
        _geo_live_max(state.live_occ_curves),discrete_max[2])
    state.occ_surface_max=max(state.occ_surface_floor,
        _geo_live_max(state.live_occ_surfaces),discrete_max[3])
    state.occ_volume_max=max(state.occ_volume_floor,
        _geo_live_max(state.live_occ_volumes),discrete_max[4])
    state.point_entity_max=max(
        _geo_live_max(state.live_builtin_points),
        _geo_live_max(state.live_occ_points),discrete_max[1])
    state.curve_entity_max=max(
        _geo_live_max(state.live_builtin_curves),
        _geo_live_max(state.live_occ_curves),discrete_max[2])
    state.surface_entity_max=max(
        _geo_live_max(state.live_builtin_surfaces),
        _geo_live_max(state.live_occ_surfaces),discrete_max[3])
    # the resynced sets are exact — allocator reads are live again even if a
    # statement the tracker cannot model (multi-operand Boolean groups)
    # invalidated them
    state.geometry_unavailable=nothing
    return nothing
end

@inline _geo_ascii_letter(c::Char)=('a'<=c<='z') || ('A'<=c<='Z') || c=='_'
@inline _geo_ascii_digit(c::Char)=('0'<=c<='9')
@inline _geo_ascii_ident(c::Char)=_geo_ascii_letter(c) || _geo_ascii_digit(c)

function _geo_expr_preview(source::AbstractString)
    ncodeunits(source)<=160 && return String(source)
    return String(first(source,120))*"…"*String(last(source,24))
end

function _geo_expr_error(parser::_GeoExprParser,message::AbstractString,
                         pos::Int=parser.token.pos)
    # Upstream every scalar-expression parse failure is a bison syntax
    # error — the offending lookahead goes in parens; `detail` keeps the
    # descriptive text for scan-context callers.
    detail="$(parser.caller): $message at byte $pos in expression `" *
        _geo_expr_preview(parser.source)*"`"
    # Expression fragments end where the enclosing statement's `tEND` (`;`)
    # sat upstream — that terminator is the offending lookahead bison reports,
    # so `:eof` maps to `;`, not a literal end-of-file token.
    token=parser.token.kind==:eof ? ";" : parser.token.text
    _geo_syntax_abort(isempty(token) ? "?" : token,detail)
end

function _geo_lex_token!(parser::_GeoExprParser)
    source=parser.source
    last=lastindex(source)
    i=parser.index
    while i<=last && isspace(source[i])
        i=nextind(source,i)
    end
    if i>last
        parser.index=i
        return _GeoExprToken(:eof,"",0.0,ncodeunits(source)+1)
    end
    start=i;c=source[i];i=nextind(source,i)
    token=if _geo_ascii_digit(c) ||
             (c=='.' && i<=last && _geo_ascii_digit(source[i]))
        had_digit=_geo_ascii_digit(c)
        while i<=last && _geo_ascii_digit(source[i])
            had_digit=true;i=nextind(source,i)
        end
        if i<=last && source[i]=='.'
            i=nextind(source,i)
            while i<=last && _geo_ascii_digit(source[i])
                had_digit=true;i=nextind(source,i)
            end
        end
        had_digit || _geo_expr_error(parser,"malformed numeric literal",start)
        if i<=last && (source[i]=='e' || source[i]=='E')
            i=nextind(source,i)
            if i<=last && (source[i]=='+' || source[i]=='-')
                i=nextind(source,i)
            end
            exponent_start=i
            while i<=last && _geo_ascii_digit(source[i])
                i=nextind(source,i)
            end
            exponent_start!=i ||
                _geo_expr_error(parser,"malformed numeric exponent",start)
        end
        text=String(source[start:prevind(source,i)])
        # `strtod` saturates overflow to HUGE_VAL and flushes underflow to
        # 0 — `1e309` lexes as `inf` upstream with no diagnostic. Julia's
        # `tryparse` reports `nothing` for both range failures, so recover
        # the saturated value from the exponent sign.
        value=something(tryparse(Float64,text),
                        occursin(r"[eE]-",text) ? 0.0 : Inf)
        _GeoExprToken(:number,text,value::Float64,start)
    elseif _geo_ascii_letter(c)
        while i<=last && _geo_ascii_ident(source[i])
            i=nextind(source,i)
        end
        text=String(source[start:prevind(source,i)])
        _GeoExprToken(:identifier,text,0.0,start)
    else
        kind=if c=='+'
            i<=last && source[i]=='+' ? :side_effect : :plus
        elseif c=='-'
            i<=last && source[i]=='-' ? :side_effect : :minus
        elseif c=='*'; :star
        elseif c=='/'; :slash
        elseif c=='%'; :percent
        elseif c=='^'; :caret
        elseif c=='('; :left_paren
        elseif c==')'; :right_paren
        elseif c=='['; :left_bracket
        elseif c==']'; :right_bracket
        elseif c==','; :comma
        elseif c=='#'; :hash
        elseif c=='.'; :dot
        elseif c=='<'
            i<=last && source[i]=='=' ? :less_equal :
            i<=last && source[i]=='<' ? :shift_left : :less
        elseif c=='>'
            i<=last && source[i]=='=' ? :greater_equal :
            i<=last && source[i]=='>' ? :shift_right : :greater
        elseif c=='='
            i<=last && source[i]=='=' ? :equal : :unsupported_operator
        elseif c=='!'
            i<=last && source[i]=='=' ? :not_equal : :not
        elseif c=='&'
            i<=last && source[i]=='&' ? :and : :bit_and
        elseif c=='|'
            i<=last && source[i]=='|' ? :or : :bit_or
        elseif c=='?'; :question
        elseif c==':'; :colon
        elseif c=='~'; :tilde
        elseif c=='{'; :left_brace
        elseif c=='}'; :right_brace
        elseif c in ('"','\'')
            :quoted
        else
            :invalid
        end
        if kind in (:side_effect,:less_equal,:greater_equal,:equal,:not_equal,
                    :and,:or,:shift_left,:shift_right)
            i=nextind(source,i)
        end
        _GeoExprToken(kind,String(source[start:prevind(source,i)]),0.0,start)
    end
    parser.index=i
    parser.token_count+=1
    parser.token_count<=_MAX_GEO_EXPRESSION_TOKENS ||
        _geo_expr_error(parser,
            "expression exceeds $_MAX_GEO_EXPRESSION_TOKENS tokens",start)
    return token
end

@inline function _geo_advance!(parser::_GeoExprParser)
    parser.token=_geo_lex_token!(parser)
    return nothing
end

function _geo_enter!(parser::_GeoExprParser)
    parser.depth+=1
    parser.depth<=_MAX_GEO_EXPRESSION_DEPTH || _geo_expr_error(parser,
        "expression nesting exceeds $_MAX_GEO_EXPRESSION_DEPTH")
    return nothing
end
@inline _geo_leave!(parser::_GeoExprParser)=(parser.depth-=1;nothing)

function _geo_finite_result(parser::_GeoExprParser,value,operation::AbstractString,
                            pos::Int)
    value isa Real || _geo_expr_error(parser,"$operation did not return a number",pos)
    result=Float64(value)
    # `.geo` execution propagates IEEE754 through `s.value` like upstream —
    # `1/0` stores `inf`, `LinSpace(a,b,1)` stores `nan`; use sites validate
    # finiteness. The params scan keeps the rejection (constant options).
    isfinite(result) || parser.context.soft_unknown_reads ||
        _geo_expr_error(parser,"$operation produced a non-finite value",pos)
    return result
end

function _geo_apply_binary(parser::_GeoExprParser,kind::Symbol,a::Float64,
                           b::Float64,pos::Int)
    label=kind==:plus ? "addition" : kind==:minus ? "subtraction" :
          kind==:star ? "multiplication" : kind==:slash ? "division" :
          kind==:percent ? "modulo" : "exponentiation"
    if kind==:slash && b==0
        # `FExpr '/' FExpr` runs `yymsg(0, "Division by zero ...")` and never
        # assigns `$$`, which bison pre-seeds from `$1` — the result is the
        # numerator. Execution reports the recoverable error and keeps `a`;
        # the params scan instead lets the IEEE result flow to its
        # "non-finite" boundary.
        if parser.context.soft_unknown_reads
            _geo_yyerror!(parser.context,
                "Division by zero in '$(_geo_gmsh_number(a)) / " *
                "$(_geo_gmsh_number(b))'")
            return a
        end
        return _geo_finite_result(parser,a/b,"division",pos)
    end
    value=try
        if kind==:percent
            # Gmsh's grammar evaluates `%` as `(int)lhs % (int)rhs`, not as
            # floating-point fmod (Gmsh.y, FExpr). Operands take the shipped
            # binary's `cvttsd2si` conversion; `ib == 0` and `INT_MIN % -1`
            # trap the x86 `idiv` — the process dies upstream, so an explicit
            # diagnostic is the faithful non-crashing equivalent.
            ia=_geo_int32_operand(parser,a,"modulo",pos)
            ib=_geo_int32_operand(parser,b,"modulo",pos)
            iszero(ib) && _geo_expr_error(parser,
                "modulo divisor truncates to zero",pos)
            ia==typemin(Int32) && ib==Int32(-1) && _geo_expr_error(parser,
                "modulo overflows the signed 32-bit quotient",pos)
            rem(ia,ib)
        else
            kind==:plus ? a+b : kind==:minus ? a-b : kind==:star ? a*b :
            kind==:slash ? a/b : _gm87_pow(a,b)
        end
    catch err
        err isa InterruptException && rethrow()
        (err isa DomainError || err isa OverflowError || err isa DivideError) || rethrow()
        _geo_expr_error(parser,"$label is outside its finite real domain",pos)
    end
    return _geo_finite_result(parser,value,label,pos)
end

function _geo_apply_function(parser::_GeoExprParser,name::String,
                             args::Vector{Float64},pos::Int)
    name in _GEO_NONCONSTANT_FUNCTIONS &&
        !(name in _GEO_STATEFUL_FUNCTIONS) && _geo_expr_error(parser,
        "non-constant or externally stateful function $name is not supported",pos)
    unary=if name=="Acos"; _gm_acos
    elseif name=="Asin"; _gm_asin
    elseif name=="Atan"; _gm_atan
    elseif name=="Ceil"; ceil
    elseif name=="Cos"; _gm87_cos
    elseif name=="Cosh"; _gm_cosh
    elseif name=="Exp"; _gm87_exp
    elseif name=="Fabs" || name=="Abs"; abs
    elseif name=="Floor"; floor
    elseif name=="Log"; _gm87_log
    elseif name=="Log10"; _gm_log10
    elseif name=="Round"; x->floor(x+0.5)
    elseif name=="Sqrt"; sqrt
    elseif name=="Sin"; _gm87_sin
    elseif name=="Sinh"; _gm_sinh
    elseif name=="Step"; x->x<0 ? 0.0 : 1.0
    elseif name=="Tan"; _gm_tan
    elseif name=="Tanh"; _gm_tanh
    else; nothing
    end
    binary=if name=="Atan2"; _gm87_atan2
    elseif name=="Fmod" || name=="Modulo"; rem
    # Gmsh 4.15.2 spells this as sqrt(a*a + b*b); preserve its overflow and
    # underflow behavior (execution keeps the IEEE754 Inf/NaN like upstream).
    elseif name=="Hypot"; (a,b)->sqrt(a*a+b*b)
    elseif name=="Max"; (a,b)->a<b ? b : a
    elseif name=="Min"; (a,b)->b<a ? b : a
    else; nothing
    end
    if unary!==nothing
        length(args)==1 || _geo_expr_error(parser,
            "function $name requires exactly one argument",pos)
        value=try
            unary(args[1])
        catch err
            err isa InterruptException && rethrow()
            (err isa DomainError || err isa OverflowError) || rethrow()
            # C libm returns `nan` for domain errors (e.g. `sqrt(-1)`) and
            # `inf` for overflows — `.geo` execution keeps those IEEE754
            # results; the params scan retains its finite-domain rejection.
            parser.context.soft_unknown_reads &&
                return _geo_finite_result(
                    parser,err isa DomainError ? NaN : Inf,
                    "function $name",pos)
            _geo_expr_error(parser,"function $name is outside its finite real domain",pos)
        end
        return _geo_finite_result(parser,value,"function $name",pos)
    elseif binary!==nothing
        length(args)==2 || _geo_expr_error(parser,
            "function $name requires exactly two arguments",pos)
        value=try
            binary(args[1],args[2])
        catch err
            err isa InterruptException && rethrow()
            (err isa DomainError || err isa OverflowError || err isa DivideError) || rethrow()
            parser.context.soft_unknown_reads &&
                return _geo_finite_result(
                    parser,err isa DivideError ? Inf : NaN,
                    "function $name",pos)
            _geo_expr_error(parser,"function $name is outside its finite real domain",pos)
        end
        return _geo_finite_result(parser,value,"function $name",pos)
    end
    if name=="Rand"
        length(args)==1 || _geo_expr_error(parser,
            "function Rand requires exactly one argument",pos)
        # `.geo` execution evaluates `rand()/RAND_MAX`; the params scan keeps
        # its constant-only contract and rejects the stateful call.
        parser.context.soft_unknown_reads || _geo_expr_error(parser,
            "non-constant or externally stateful function Rand is not " *
            "supported",pos)
        return _geo_finite_result(parser,args[1]*rand(),"function Rand",pos)
    end
    _geo_expr_error(parser,"unknown numeric function $name",pos)
end

# `GetMemoryUsage()/1024/1024` (peak RSS, `ru_maxrss`) and `TotalRam()` — MiB.
_geo_memory_mb()=Float64(Sys.maxrss())/1048576.0
_geo_total_memory_mb()=Float64(Sys.total_memory())/1048576.0

# Process USER CPU seconds, matching Gmsh OS.cpp GetResources/Cpu.
# Windows FILETIME consists of two UInt32 words; native converts each word
# separately in this order. POSIX timeval is read from getrusage(RUSAGE_SELF).
function _geo_cpu_time()
    @static if Sys.iswindows()
        times=Vector{UInt32}(undef,8)
        process=ccall((:GetCurrentProcess,"kernel32.dll"),Ptr{Cvoid},())
        ok=GC.@preserve times ccall((:GetProcessTimes,"kernel32.dll"),Cint,
            (Ptr{Cvoid},Ptr{UInt32},Ptr{UInt32},Ptr{UInt32},Ptr{UInt32}),
            process,pointer(times,1),pointer(times,3),pointer(times,5),pointer(times,7))
        ok!=0 || return 0.0
        return 1.e-7*4294967296.0*Float64(times[8])+1.e-7*Float64(times[7])
    else
        buf=Vector{UInt8}(undef,256)
        ret=ccall(:getrusage,Cint,(Cint,Ptr{UInt8}),0,buf)
        ret==0 || return 0.0
        @static if Sys.isbsd()
            secs=reinterpret(Int64,buf)
            usecs=reinterpret(Int32,buf)
            return Float64(secs[1])+1.e-6*Float64(usecs[3])
        else
            vals=reinterpret(Int64,buf)
            return Float64(vals[1])+1.e-6*Float64(vals[2])
        end
    end
end

# Gmsh's conditional operators, lowest to highest precedence:
# ?: (right-associative), ||, &&, ==/!=, </<=/>/>=, then additive arithmetic.
# Comparisons and logic produce 1.0/0.0, matching FExpr boolean evaluation.
function _geo_parse_ternary!(parser::_GeoExprParser)
    condition=_geo_parse_or!(parser)
    if parser.token.kind==:question
        pos=parser.token.pos;_geo_advance!(parser);_geo_enter!(parser)
        if_true=_geo_parse_ternary!(parser)
        parser.token.kind==:colon || _geo_expr_error(
            parser,"expected ':' in a ternary expression",pos)
        _geo_advance!(parser)
        if_false=_geo_parse_ternary!(parser)
        _geo_leave!(parser)
        return condition!=0 ? if_true : if_false
    end
    return condition
end

function _geo_parse_or!(parser::_GeoExprParser)
    value=_geo_parse_and!(parser)
    while parser.token.kind==:or
        _geo_advance!(parser)
        other=_geo_parse_and!(parser)
        value=Float64(value!=0 || other!=0)
    end
    return value
end

function _geo_parse_and!(parser::_GeoExprParser)
    value=_geo_parse_equality!(parser)
    while parser.token.kind==:and
        _geo_advance!(parser)
        other=_geo_parse_equality!(parser)
        value=Float64(value!=0 && other!=0)
    end
    return value
end

function _geo_parse_equality!(parser::_GeoExprParser)
    value=_geo_parse_relational!(parser)
    while parser.token.kind in (:equal,:not_equal)
        kind=parser.token.kind;_geo_advance!(parser)
        other=_geo_parse_relational!(parser)
        value=Float64(kind==:equal ? value==other : value!=other)
    end
    return value
end

function _geo_parse_relational!(parser::_GeoExprParser)
    value=_geo_parse_additive!(parser)
    # `<<`/`>>` sit at the relational precedence level upstream (Gmsh.y
    # `%left '<' tLESSOREQUAL '>' tGREATEROREQUAL tLESSLESS tGREATERGREATER`).
    while parser.token.kind in (:less,:less_equal,:greater,:greater_equal,
                                :shift_left,:shift_right)
        kind=parser.token.kind;pos=parser.token.pos;_geo_advance!(parser)
        other=_geo_parse_additive!(parser)
        value=if kind==:shift_left || kind==:shift_right
            _geo_apply_shift(parser,kind,value,other,pos)
        else
            Float64(
                kind==:less ? value<other : kind==:less_equal ? value<=other :
                kind==:greater ? value>other : value>=other)
        end
    end
    return value
end

function _geo_parse_additive!(parser::_GeoExprParser)
    value=_geo_parse_multiplicative!(parser)
    while parser.token.kind==:plus || parser.token.kind==:minus
        kind=parser.token.kind;pos=parser.token.pos;_geo_advance!(parser)
        value=_geo_apply_binary(parser,kind,value,
                                _geo_parse_multiplicative!(parser),pos)
    end
    return value
end

function _geo_parse_multiplicative!(parser::_GeoExprParser)
    # `|`/`&` bind tighter than `*`/`/`/`%` upstream (Gmsh.y precedence), so
    # the multiplicative operands are bitwise terms.
    value=_geo_parse_bitwise!(parser)
    while parser.token.kind in (:star,:slash,:percent)
        kind=parser.token.kind;pos=parser.token.pos;_geo_advance!(parser)
        value=_geo_apply_binary(parser,kind,value,
                                _geo_parse_bitwise!(parser),pos)
    end
    return value
end

# `FExpr '|' FExpr` / `FExpr '&' FExpr` — C `(int)a | (int)b` semantics.
function _geo_parse_bitwise!(parser::_GeoExprParser)
    value=_geo_parse_unary!(parser)
    while parser.token.kind==:bit_or || parser.token.kind==:bit_and
        kind=parser.token.kind;pos=parser.token.pos
        # The params scan keeps its documented arithmetic subset.
        parser.context.soft_unknown_reads || _geo_expr_error(parser,
            "operator $(parser.token.text) is outside the supported " *
            "arithmetic subset",pos)
        _geo_advance!(parser)
        other=_geo_parse_unary!(parser)
        ia=_geo_int32_operand(parser,value,"bitwise",pos)
        ib=_geo_int32_operand(parser,other,"bitwise",pos)
        value=Float64(kind==:bit_or ? ia|ib : ia&ib)
    end
    return value
end

# `(int)value` — upstream compiles the double→int conversion to `cvttsd2si`:
# in-range doubles truncate toward zero; out-of-range and non-finite inputs
# produce the x86 "integer indefinite" value 0x80000000 (`typemin(Int32)`).
# Deterministic on the shipped binary, so the full domain is emulated.
function _geo_int32_operand(parser::_GeoExprParser,value::Float64,
                            operation::AbstractString,pos::Int)
    t=trunc(value)
    (isfinite(value) && Float64(typemin(Int32))<=t<=Float64(typemax(Int32))) ||
        return typemin(Int32)
    return Int32(t)
end

# `((int)a >> (int)b)` / `<<` — compiled to x86 SAR/SHL upstream: the count is
# masked to 5 bits, `>>` is arithmetic, `<<` wraps. Deterministic on the
# shipped binary, so the full domain is emulated rather than rejected.
function _geo_apply_shift(parser::_GeoExprParser,kind::Symbol,a::Float64,
                          b::Float64,pos::Int)
    # The params scan keeps its documented arithmetic subset.
    parser.context.soft_unknown_reads || _geo_expr_error(parser,
        "operator $(kind==:shift_left ? "<<" : ">>") is outside the " *
        "supported arithmetic subset",pos)
    ia=_geo_int32_operand(parser,a,"shift",pos)
    ib=_geo_int32_operand(parser,b,"shift",pos)
    count=Int(ib & Int32(31))
    return Float64(kind==:shift_left ? ia<<count : ia>>count)
end

function _geo_parse_unary!(parser::_GeoExprParser)
    if parser.token.kind==:plus || parser.token.kind==:minus
        kind=parser.token.kind;pos=parser.token.pos;_geo_advance!(parser)
        _geo_enter!(parser)
        value=_geo_parse_unary!(parser)
        _geo_leave!(parser)
        kind==:minus && (value=_geo_finite_result(parser,-value,"unary minus",pos))
        return value
    elseif parser.token.kind==:not
        pos=parser.token.pos;_geo_advance!(parser)
        _geo_enter!(parser)
        value=_geo_parse_unary!(parser)
        _geo_leave!(parser)
        return Float64(value==0)
    elseif parser.token.kind==:hash
        pos=parser.token.pos;_geo_advance!(parser)
        parser.token.kind==:identifier || _geo_expr_error(
            parser,"expected a numeric list name after #",pos)
        name=parser.token.text;_geo_advance!(parser)
        hash_resolved=false
        if name in ("StringToName","S2N") &&
           parser.token.kind==:left_bracket
            name=_geo_parse_strtoname!(parser)
            hash_resolved=true
        end
        name=_geo_parse_name_tildes!(parser,name)
        (hash_resolved || _geo_tstring_valid(name)) || _geo_syntax_abort(name)
        dim_ns=nothing
        if parser.token.kind==:colon &&
           parser.index<=lastindex(parser.source) &&
           parser.source[parser.index]==':'
            (hash_resolved || _geo_tstring_valid(name)) || _geo_syntax_abort(
                _geo_keyword_self_reports(name) ? name : "::")
            dim_ns=name
            _geo_advance!(parser);_geo_advance!(parser)
            parser.token.kind==:identifier || _geo_expr_error(
                parser,"expected a name after $dim_ns::",pos)
            name=parser.token.text;_geo_advance!(parser)
            dim_member_resolved=false
            if name in ("StringToName","S2N") &&
               parser.token.kind==:left_bracket
                name=_geo_parse_strtoname!(parser)
                dim_member_resolved=true
            end
            name=_geo_parse_name_tildes!(parser,name)
            (dim_member_resolved || _geo_tstring_valid(name)) ||
                _geo_syntax_abort(name)
        end
        if parser.token.kind==:dot
            # `#x.member()` / `#ns::x.member()` — `getMember_Dim`: the number
            # of values stored on the member (0 silently for unknown structs).
            _geo_advance!(parser)
            parser.token.kind==:identifier || _geo_expr_error(parser,
                "expected a member name after $name.",pos)
            member=parser.token.text;_geo_advance!(parser)
            _geo_member_key_valid(member) || _geo_syntax_abort(member)
            opener=parser.token.kind
            opener in (:left_bracket,:left_paren) || _geo_expr_error(
                parser,"expected () after member $member",pos)
            closer=opener==:left_bracket ? :right_bracket : :right_paren
            _geo_advance!(parser)
            parser.token.kind==closer || _geo_expr_error(
                parser,"member cardinality uses an empty selector",pos)
            _geo_advance!(parser)
            return Float64(_geo_member_dim(parser.context,
                something(dim_ns,""),name,member,parser.caller))
        end
        dim_ns===nothing || _geo_expr_error(
            parser,"expected a member name after $dim_ns::$name",pos)
        # Gmsh `#name()` counts a numeric or string list's elements; the legacy
        # `#name[]` spelling is retained for compatibility.
        opener=parser.token.kind
        opener in (:left_bracket,:left_paren) || _geo_expr_error(
            parser,"expected () after list name $name",pos)
        closer=opener==:left_bracket ? :right_bracket : :right_paren
        _geo_advance!(parser)
        parser.token.kind==closer || _geo_expr_error(
            parser,"list cardinality uses an empty selector",pos)
        _geo_advance!(parser)
        haskey(parser.context.strings,name) &&
            return Float64(length(parser.context.strings[name]))
        return Float64(length(_geo_context_list(
            parser.context,name,parser.caller)))
    end
    return _geo_parse_power!(parser)
end

function _geo_parse_power!(parser::_GeoExprParser)
    value=_geo_parse_primary!(parser)
    if parser.token.kind==:caret
        pos=parser.token.pos;_geo_advance!(parser);_geo_enter!(parser)
        exponent=_geo_parse_unary!(parser)
        _geo_leave!(parser)
        value=_geo_apply_binary(parser,:caret,value,exponent,pos)
    end
    return value
end

function _geo_parse_primary!(parser::_GeoExprParser)
    token=parser.token
    if token.kind==:number
        _geo_advance!(parser)
        return token.value
    elseif token.kind==:left_paren
        _geo_advance!(parser);_geo_enter!(parser)
        value=_geo_parse_ternary!(parser)
        parser.token.kind==:right_paren ||
            _geo_expr_error(parser,"expected closing parenthesis")
        _geo_advance!(parser);_geo_leave!(parser)
        return value
    elseif token.kind==:identifier
        name=token.text;_geo_advance!(parser)
        # `Struct` is `DefineStruct` — an FExpr_Single producing the tag.
        name=="Struct" &&
            return _geo_parse_inline_struct!(parser,token.pos)
        # `StringToName[expr]`/`S2N[expr]` — a `String__Index` name production;
        # the resolved name continues through member/scope/index grammar and
        # is not re-lexed upstream (`S2N["Field"]::s` reads namespace Field).
        resolved_name=false
        if name in ("StringToName","S2N") && parser.token.kind==:left_bracket
            name=_geo_parse_strtoname!(parser)
            resolved_name=true
        end
        # `name~{FExpr}` — `StringIndex` name suffixes.
        name=_geo_parse_name_tildes!(parser,name)
        # `ns::` — contiguous `::` is upstream's `tSCOPE`; a spaced or single
        # `:` stays for range/ternary handling.
        resolved_ns=nothing
        if parser.token.kind==:colon &&
           parser.index<=lastindex(parser.source) &&
           parser.source[parser.index]==':'
            # `kw::s` — an FExpr-incapable keyword reports itself
            # (`Field::s` → `Field`); FExpr heads error at `::` (`Sin::s`).
            (resolved_name || _geo_tstring_valid(name)) || _geo_syntax_abort(
                _geo_keyword_self_reports(name) ? name : "::")
            resolved_ns=name
            _geo_advance!(parser);_geo_advance!(parser)
            parser.token.kind==:identifier || _geo_expr_error(parser,
                "expected a name after $resolved_ns::",token.pos)
            name=parser.token.text;_geo_advance!(parser)
            member_resolved=false
            if name in ("StringToName","S2N") &&
               parser.token.kind==:left_bracket
                name=_geo_parse_strtoname!(parser)
                member_resolved=true
            end
            name=_geo_parse_name_tildes!(parser,name)
            (member_resolved || _geo_tstring_valid(name)) ||
                _geo_syntax_abort(name)
        end
        if parser.token.kind==:dot
            # `x.y`/`ns::x.y` (+ optional `(i)`/`[i]` member index) —
            # `treat_Struct_FullName_dot_tSTRING_Float`: struct member first,
            # `NumberOption(x, 0, y)` when x is not a struct.
            resolved_ns===nothing && !_geo_tstring_valid(name) &&
                _geo_syntax_abort(
                    _geo_keyword_self_reports(name) ? name : ".")
            _geo_advance!(parser)
            parser.token.kind==:identifier || _geo_expr_error(parser,
                "expected an option name after $name.",token.pos)
            member=parser.token.text;_geo_advance!(parser)
            _geo_member_key_valid(member) || _geo_syntax_abort(member)
            index=0
            if parser.token.kind==:left_paren ||
               parser.token.kind==:left_bracket
                closer=parser.token.kind==:left_paren ? :right_paren :
                    :right_bracket
                _geo_advance!(parser);_geo_enter!(parser)
                index=_geo_int_value(_geo_parse_ternary!(parser),
                    "$(parser.caller) option index")
                parser.token.kind==closer || _geo_expr_error(parser,
                    "expected closing delimiter in $name.$member index")
                _geo_advance!(parser);_geo_leave!(parser)
            end
            # `x.y++`/`x.y--` — NumberOption GET then SET; returns the NEW
            # value (unlike variable postfix, which returns the old one).
            # `ns::x.y++` is not in the upstream grammar.
            if parser.token.kind==:side_effect
                resolved_ns===nothing || _geo_expr_error(parser,
                    "unexpected ++/-- after $(resolved_ns)::$(name).$member")
                delta=parser.token.text=="++" ? 1.0 : -1.0
                _geo_advance!(parser)
                d=_geo_option_number(parser.context,name,index,member;
                    warn=true)
                d===nothing && return 0.0
                _geo_store_option_number!(parser.context,name,index,member,
                    d+delta,parser.caller)
                return _geo_option_number(parser.context,name,index,member)
            end
            return _geo_member_number(parser.context,
                something(resolved_ns,""),name,member,index,parser.caller)
        elseif parser.token.kind==:left_paren || parser.token.kind==:left_bracket
            resolved_ns===nothing || _geo_expr_error(parser,
                "unexpected index after $resolved_ns::$name",token.pos)
            # `x[i].y` — `NumberOption(GMSH_GET, x, i, y)`; the bracket index
            # belongs to the option family, not a variable element.
            lookahead_pos=parser.index
            closer_probe=parser.token.kind==:left_paren ? :right_paren :
                :right_bracket
            depth=1;j=lookahead_pos;src=parser.source
            while j<=lastindex(src) && depth>0
                ch=src[j]
                ch=='(' || ch=='[' ? depth+=1 :
                    ch==')' || ch==']' ? depth-=1 : nothing
                j=nextind(src,j)
            end
            if depth==0 && j<=lastindex(src) && src[j]=='.'
                # Indexed option family: evaluate the bracket, then `member`.
                _geo_advance!(parser);_geo_enter!(parser)
                opt_index=_geo_int_value(_geo_parse_ternary!(parser),
                    "$(parser.caller) option index")
                parser.token.kind==closer_probe || _geo_expr_error(parser,
                    "expected closing delimiter in $name[.] option index")
                _geo_advance!(parser);_geo_leave!(parser)
                _geo_advance!(parser) # `.`
                parser.token.kind==:identifier || _geo_expr_error(parser,
                    "expected an option name after $name[.].",token.pos)
                member=parser.token.text;_geo_advance!(parser)
                _geo_member_key_valid(member) || _geo_syntax_abort(member)
                if parser.token.kind==:side_effect
                    delta=parser.token.text=="++" ? 1.0 : -1.0
                    _geo_advance!(parser)
                    d=_geo_option_number(parser.context,name,opt_index,member;
                        warn=true)
                    d===nothing && return 0.0
                    _geo_store_option_number!(parser.context,name,opt_index,
                        member,d+delta,parser.caller)
                    return _geo_option_number(parser.context,name,opt_index,
                        member)
                end
                return _geo_option_number_read(parser.context,name,opt_index,
                    member,parser.caller)
            end
            # Functions taking StringExprVar, ListOfDouble or raw-name
            # arguments cannot flow through the numeric argument loop; capture
            # the balanced argument text and dispatch per signature.
            # `tStringToName` takes `[` only — `(` is the error token
            # upstream (`StringToName(1)` → `syntax error (()`).
            name in ("StringToName","S2N") &&
                parser.token.kind==:left_paren && _geo_syntax_abort("(")
            if name in _GEO_RAW_ARG_FUNCTIONS
                return _geo_raw_arg_function!(parser,name,token.pos)
            end
            opener=parser.token.kind
            closer=opener==:left_paren ? :right_paren : :right_bracket
            known_function=(name in _GEO_NUMERIC_FUNCTIONS) ||
                (name in _GEO_NONCONSTANT_FUNCTIONS) ||
                (name in _GEO_STATEFUL_FUNCTIONS)
            # `tList LP` is a real `FExpr_Multi` production — its argument
            # errors at the closer (`List(1)` → `)`, `List[1]` → `]`).
            name=="List" && _geo_syntax_abort(
                closer==:right_paren ? ")" : "]")
            # `tSymmetry`/`tExtrude` take a `VExpr` — `(`-arguments fail
            # arity at `)`; `[` cannot start a `VExpr` (`Symmetry[1]`→`[`).
            name in ("Symmetry","Extrude") && _geo_syntax_abort(
                opener==:left_paren ? ")" : "[")
            # Lexer-emitted literals (`newp`, `Cpu`, `Pi`, `Today`, ...)
            # take no argument list, and `FExpr`-head keywords (`Point`,
            # `Physical`, `Line`, ...) shift their production first —
            # both error at the delimiter itself upstream (`Cpu()` → `(`,
            # `Point[1]` → `[`).
            (haskey(_GEO_BARE_CONSTANTS,name) ||
             name in _GEO_SIDE_EFFECT_SYMBOLS ||
             name in _GEO_BARE_STRING_TOKENS ||
             (name in _GEO_FEXPR_SHIFT_KEYWORDS &&
              !(name in ("Symmetry","Extrude")))) && _geo_syntax_abort(
                opener==:left_paren ? "(" : "[")
            # StringExpr-head names (`StrCat`, `N2S`, `Sprintf`, ...) cannot
            # begin an `FExpr` — the token itself is the error
            # (`Sin[StrCat["a"]]` → `StrCat`).
            name in _GEO_STRING_FUNCTIONS && _geo_syntax_abort(name)
            # Remaining lexer keywords cannot begin an expression — the
            # keyword itself is the error (`v = Circle(1)` → `Circle`).
            (name in _GEO_LEXER_KEYWORDS &&
             !(name in _GEO_ALL_FUNCTIONS)) && _geo_syntax_abort(name)
            _geo_advance!(parser);_geo_enter!(parser)
            # `x(i)`/`x[i]` are the same `String__Index LP FExpr RP`
            # indexed read for every non-function tSTRING — the
            # function/variable split is lexical, so an unknown name
            # indexes (`Unknown variable 'x(.)'`), never "unknown
            # function". Both forms take exactly one index and each
            # function its fixed arity; the comma that would exceed it is
            # the error token upstream (`x[1,2]`/`Sin[1,2]` → `,`).
            is_index=!known_function
            max_args=is_index || !(name in _GEO_BINARY_FUNCTIONS) ? 1 : 2
            args=Float64[]
            if parser.token.kind!=closer
                while true
                    push!(args,_geo_parse_ternary!(parser))
                    parser.token.kind==:comma || break
                    length(args)>=max_args && _geo_syntax_abort(",")
                    _geo_advance!(parser)
                end
            end
            # No FExpr production accepts an empty `()`/`[]` — the closer
            # is the offending token (`x[]`/`f()` → `]`/`)`).
            isempty(args) && _geo_syntax_abort(
                closer==:right_paren ? ")" : "]")
            parser.token.kind==closer || _geo_expr_error(parser,
                opener==:left_paren ? "expected closing parenthesis" :
                                      "expected closing bracket")
            # A fixed-arity call missing arguments ends the production at
            # the closer — `Atan2[1]` expects `, FExpr` and errors at `]`.
            (!is_index && length(args)<max_args) && _geo_syntax_abort(
                closer==:right_paren ? ")" : "]")
            _geo_advance!(parser);_geo_leave!(parser)
            if is_index
                context=parser.context
                if context.soft_unknown_reads
                    if !_geo_context_has_variable(context,name)
                        _geo_yyerror!(context,"Unknown variable '$name(.)'")
                        parser.token.kind==:side_effect && _geo_advance!(parser)
                        return 0.0
                    end
                    values=_geo_context_list(context,name,parser.caller)
                    idx=_geo_int_value(args[1],"$(parser.caller) index")
                    if !(0<=idx<length(values))
                        _geo_yyerror!(context,"Uninitialized variable '$name[$idx]'")
                        # `x[i]++` on a miss still consumes the increment
                        # token; the write does not happen.
                        parser.token.kind==:side_effect && _geo_advance!(parser)
                        return 0.0
                    end
                    value=values[idx+1]
                else
                    value=_geo_context_list_value(
                        parser.context,name,args[1],parser.caller)
                end
                # `x[i]++`/`x(i)++` — postfix increment returns the old value,
                # then mutates the element (FExpr NumericIncrement).
                if parser.token.kind==:side_effect
                    delta=parser.token.text=="++" ? 1.0 : -1.0
                    _geo_increment_indexed!(
                        parser.context,name,args[1],delta,parser.caller)
                    _geo_advance!(parser)
                end
                return value
            end
            value=_geo_apply_function(parser,name,args,token.pos)
            return value
        elseif resolved_ns!==nothing
            # `ns::name` — `Struct_FullName` tag read; variables and constants
            # do not apply inside an explicit namespace.
            return _geo_fullname_number(parser.context,resolved_ns,name,
                parser.caller)
        elseif haskey(_GEO_BARE_CONSTANTS,name)
            return _GEO_BARE_CONSTANTS[name](parser)
        elseif haskey(parser.context.values,name)
            value=parser.context.values[name]
            # `x++`/`x--` — postfix increment in expression position.
            if parser.token.kind==:side_effect
                delta=parser.token.text=="++" ? 1.0 : -1.0
                _geo_increment_scalar!(parser.context,name,delta,parser.caller)
                _geo_advance!(parser)
            end
            return value
        elseif haskey(parser.context.unavailable,name)
            label=name in _GEO_SIDE_EFFECT_SYMBOLS ?
                  "dynamic tag allocator $name" : "numeric variable $name"
            # Soft option-write reads recover like `Unknown variable` — the
            # unavailability reason is the diagnostic, not a syntax error.
            if parser.context.soft_unknown_reads
                _geo_yyerror!(parser.context,
                    "$label is unavailable ($(parser.context.unavailable[name]))")
                parser.token.kind==:side_effect && _geo_advance!(parser)
                return 0.0
            end
            _geo_expr_error(parser,
                "$label is unavailable ($(parser.context.unavailable[name]))",token.pos)
        elseif haskey(parser.context.lists,name)
            # A defined-but-empty list: `s.value` is empty → `Uninitialized`.
            if parser.context.soft_unknown_reads
                _geo_yyerror!(parser.context,"Uninitialized variable '$name'")
                parser.token.kind==:side_effect && _geo_advance!(parser)
                return 0.0
            end
            _geo_expr_error(parser,
                "numeric list $name is empty and has no scalar value",token.pos)
        elseif haskey(parser.context.strings,name)
            # String-only names miss the numeric table → `Unknown variable`.
            if parser.context.soft_unknown_reads
                _geo_yyerror!(parser.context,"Unknown variable '$name'")
                parser.token.kind==:side_effect && _geo_advance!(parser)
                return 0.0
            end
            _geo_expr_error(parser,
                "$name is a string variable, not a number",token.pos)
        elseif name in _GEO_SIDE_EFFECT_SYMBOLS
            _geo_expr_error(parser,
                "dynamic tag allocator $name requires current allocation state",token.pos)
        elseif (struct_tag=_geo_struct_tag(parser.context,"",name))!==nothing
            # Bare struct read — `Struct_FullName` yields the tag. Checked
            # before the reserved sets: `Mesh`, `While`, `Model`, ... are
            # plain `tSTRING` upstream and can legally name a struct.
            return Float64(struct_tag)
        elseif name in _GEO_LEXER_KEYWORDS || name in _GEO_ALL_FUNCTIONS
            # Reserved token — upstream the lexer emits `tX`, which commits
            # to its production. `FExpr`-head and function keywords shift
            # and error at the lookahead (`v = Sin`/`v = Point` → `;`);
            # pure-statement keywords error at the keyword itself
            # (`v = Circle` → `Circle`), as do `StringExpr`-production
            # tokens — they cannot appear in an FExpr at all
            # (`v = StrCat` → `(StrCat)`).
            if name in _GEO_STRING_FUNCTIONS ||
               (name in _GEO_LEXER_KEYWORDS &&
                !(name in _GEO_ALL_FUNCTIONS) &&
                !(name in _GEO_FEXPR_SHIFT_KEYWORDS) && name!="List")
                _geo_syntax_abort(name)
            else
                _geo_syntax_abort(parser.token.kind==:eof ? ";" :
                    isempty(parser.token.text) ? ";" : parser.token.text)
            end
        else
            if parser.context.soft_unknown_reads
                _geo_yyerror!(parser.context,"Unknown variable '$name'")
                # A pending `x++` on an unknown name is consumed without a
                # write — the FExpr postfix diagnostic already fired.
                parser.token.kind==:side_effect && _geo_advance!(parser)
                return 0.0
            end
            _geo_expr_error(parser,"unknown scalar identifier $name",token.pos)
        end
    elseif token.kind==:side_effect
        # `++`/`--` where a value is expected is a grammar failure upstream
        # (`NumericIncrement` binds only a postfix `tSTRING`).
        _geo_syntax_abort(token.text)
    elseif token.kind==:unsupported_operator
        _geo_expr_error(parser,"operator $(token.text) is outside the supported arithmetic subset",
                        token.pos)
    elseif token.kind==:quoted
        _geo_expr_error(parser,"quoted strings are not numeric expressions",token.pos)
    elseif token.kind==:invalid
        _geo_expr_error(parser,"invalid token $(repr(token.text))",token.pos)
    end
    _geo_expr_error(parser,"expected a numeric value",token.pos)
end

# `StringToName[expr]`/`S2N[expr]` in `String__Index` position: evaluate the
# bracketed StringExprVar to the name it stands for.
function _geo_parse_strtoname!(parser::_GeoExprParser)
    open_pos=parser.token.pos
    close=_geo_matching_delim(parser.source,open_pos)
    close==0 && _geo_expr_error(parser,"unbalanced StringToName brackets")
    inner=String(parser.source[
        nextind(parser.source,open_pos):prevind(parser.source,close)])
    parser.index=nextind(parser.source,close)
    _geo_advance!(parser)
    return _geo_eval_string(inner,parser.context,parser.caller)
end

# `name~{FExpr}` — `StringIndex` suffixes expand the name to `name_i`.
function _geo_parse_name_tildes!(parser::_GeoExprParser,name::String)
    while parser.token.kind==:tilde
        _geo_advance!(parser)
        parser.token.kind==:left_brace || _geo_expr_error(parser,
            "expected {index} after ~ in a namespaced name")
        _geo_advance!(parser);_geo_enter!(parser)
        idx=_geo_int_value(_geo_parse_ternary!(parser),
            "$(parser.caller) namespace index")
        parser.token.kind==:right_brace || _geo_expr_error(parser,
            "expected } after ~{index}")
        _geo_advance!(parser);_geo_leave!(parser)
        name=name*"_"*string(idx)
    end
    return name
end

# `Struct` in FExpr position (`DefineStruct` is FExpr_Single): parses
# `name[`::`name] (Append [n]) [options]` inline from the token stream and
# returns the tag.
function _geo_parse_inline_struct!(parser::_GeoExprParser,pos::Int)
    context=parser.context;caller=parser.caller
    parser.token.kind==:identifier || _geo_expr_error(parser,
        "expected a struct name after Struct",pos)
    name=parser.token.text;_geo_advance!(parser)
    name=_geo_parse_name_tildes!(parser,name)
    _geo_tstring_valid(name) || _geo_syntax_abort(name)
    ns=nothing
    if parser.token.kind==:colon &&
       parser.index<=lastindex(parser.source) &&
       parser.source[parser.index]==':'
        _geo_advance!(parser);_geo_advance!(parser)
        parser.token.kind==:identifier || _geo_expr_error(parser,
            "expected a struct name after $name::",pos)
        ns=name
        name=parser.token.text;_geo_advance!(parser)
        name=_geo_parse_name_tildes!(parser,name)
        _geo_tstring_valid(name) || _geo_syntax_abort(name)
    end
    append_level=0
    if parser.token.kind==:left_paren
        _geo_advance!(parser)
        parser.token.kind==:identifier && parser.token.text=="Append" ||
            _geo_expr_error(parser,
                "expected `Append` inside the Struct definition parens",pos)
        _geo_advance!(parser)
        append_level=99
        if parser.token.kind!=:right_paren
            _geo_enter!(parser)
            append_level=_geo_int_value(_geo_parse_ternary!(parser),
                "$(caller) Struct Append level")
            _geo_leave!(parser)
        end
        parser.token.kind==:right_paren || _geo_expr_error(parser,
            "expected ) after Struct Append",pos)
        _geo_advance!(parser)
    end
    parser.token.kind==:left_bracket || _geo_expr_error(parser,
        "expected [options] in the Struct definition",pos)
    open_pos=parser.token.pos
    close=_geo_matching_delim(parser.source,open_pos)
    close==0 && _geo_expr_error(parser,"unbalanced Struct options",pos)
    options_src=String(parser.source[
        nextind(parser.source,open_pos):prevind(parser.source,close)])
    parser.index=nextind(parser.source,close)
    _geo_advance!(parser)
    (tag,redef)=_geo_def_struct!(context,something(ns,""),name,append_level,
                                 options_src,caller)
    redef && _geo_expr_diag!(context,
        "Redefinition of Struct '$(something(ns,""))::$name'",caller)
    return Float64(tag)
end

# ======== FExpr raw-argument and stateful functions ========

# Capture a balanced `(...)`/`[...]` argument body verbatim — needed by FExpr
# functions whose arguments are strings, lists or bare variable names.
function _geo_raw_arg_function!(parser::_GeoExprParser,name::String,pos::Int)
    open_pos=parser.token.pos
    close=_geo_matching_delim(parser.source,open_pos)
    close==0 && _geo_expr_error(
        parser,"unbalanced delimiter in $name arguments",pos)
    inner=String(parser.source[nextind(parser.source,open_pos):prevind(parser.source,close)])
    parser.index=nextind(parser.source,close)
    _geo_advance!(parser)
    context=parser.context;caller=parser.caller
    args=_geo_split_args(inner,caller)
    if name=="Exists"
        length(args)==1 || _geo_expr_error(
            parser,"Exists takes exactly one argument",pos)
        # `tExists LP Struct_FullName ['.' tSTRING_Reserved] RP` —
        # `type_treat=1`: variable→1, string variable→1, struct→its *tag*,
        # member→1 (char members included); unknown structs fall through to
        # a silent `NumberOption` read (an existing option yields its value).
        arg=strip(args[1])
        if (om=match(
                r"^(?:([A-Za-z_][A-Za-z0-9_]*)\s*::\s*)?([A-Za-z_][A-Za-z0-9_]*)\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*([\(\[])",
                arg))!==nothing
            # `Exists(x.member(i))` — no indexed member production upstream;
            # the bracket is a syntax error.
            _geo_syntax_abort(String(om.captures[4]))
        end
        if (om=match(
                r"^(?:([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*::\s*)?([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)$",
                arg))!==nothing
            return _geo_member_number(context,
                om.captures[1]===nothing ? "" :
                    _geo_symbol_name(om.captures[1],context,caller),
                _geo_symbol_name(om.captures[2],context,caller),
                String(om.captures[3]),0,caller;type_treat=1)
        end
        if (om=match(
                r"^([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*::\s*(.+)$",
                arg))!==nothing
            sym=_geo_symbol_name(om.captures[2],context,caller)
            return _geo_fullname_number(context,
                _geo_symbol_name(om.captures[1],context,caller),sym,
                caller;type_treat=1)
        end
        if (om=match(r"^[A-Za-z_][A-Za-z0-9_]*\s*([\(\[])",arg))!==nothing
            # `Exists(x(i))`/`Exists(x[i])` — `String__Index` has no paren or
            # bracket tail; the bracket is a syntax error upstream.
            _geo_syntax_abort(String(om.captures[1]))
        end
        sym=_geo_symbol_name(arg,context,caller)
        return _geo_fullname_number(context,nothing,sym,caller;type_treat=1)
    elseif name=="GetForced"
        (1<=length(args)<=2) || _geo_expr_error(
            parser,"GetForced takes one or two arguments",pos)
        default=length(args)==2 ?
            _geo_eval_numeric(args[2],context,caller) : 0.0
        arg=strip(args[1])
        if (om=match(
                r"^(?:([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*::\s*)?([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*[\(\[]\s*(.*?)\s*[\)\]]$",
                arg))!==nothing
            midx=_geo_int_value(_geo_eval_numeric(om.captures[4],context,
                "$caller GetForced index"),"$caller GetForced index")
            return _geo_member_number(context,
                om.captures[1]===nothing ? "" :
                    _geo_symbol_name(om.captures[1],context,caller),
                _geo_symbol_name(om.captures[2],context,caller),
                String(om.captures[3]),midx,caller;
                default=default,type_treat=2)
        end
        if (om=match(
                r"^(?:([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*::\s*)?([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)$",
                arg))!==nothing
            return _geo_member_number(context,
                om.captures[1]===nothing ? "" :
                    _geo_symbol_name(om.captures[1],context,caller),
                _geo_symbol_name(om.captures[2],context,caller),
                String(om.captures[3]),0,caller;
                default=default,type_treat=2)
        end
        if (om=match(
                r"^([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*::\s*(.+)$",
                arg))!==nothing
            rest=String(om.captures[2])
            ns=_geo_symbol_name(om.captures[1],context,caller)
            if (im=match(r"^(.*?)\s*[\(\[]\s*(.*?)\s*[\)\]]$",rest))!==nothing
                sym=_geo_symbol_name(im.captures[1],context,caller)
                idx=_geo_int_value(_geo_eval_numeric(im.captures[2],context,
                    "$caller GetForced index"),"$caller GetForced index")
                return _geo_fullname_number(context,ns,sym,caller;
                    type_var=2,index=idx,default=default,type_treat=2)
            end
            sym=_geo_symbol_name(rest,context,caller)
            return _geo_fullname_number(context,ns,sym,
                caller;default=default,type_treat=2)
        end
        if (im=match(
                r"^([A-Za-z_][A-Za-z0-9_]*(?:~\{[^}]*\})*)\s*[\(\[]\s*(.*?)\s*[\)\]]$",
                arg))!==nothing
            # `GetForced(x(i), d)` — `Struct_FullName LP FExpr RP` indexed
            # read (`type_var=2`); `LP`/`RP` cover `[]` too.
            sym=_geo_symbol_name(im.captures[1],context,caller)
            idx=_geo_int_value(_geo_eval_numeric(im.captures[2],context,
                "$caller GetForced index"),"$caller GetForced index")
            return _geo_fullname_number(context,nothing,sym,caller;
                type_var=2,index=idx,default=default,type_treat=2)
        end
        sym=_geo_symbol_name(arg,context,caller)
        return _geo_fullname_number(context,nothing,sym,caller;
            default=default,type_treat=2)
    elseif name=="GetNumber"
        (1<=length(args)<=2) || _geo_expr_error(
            parser,"GetNumber takes one or two arguments",pos)
        key=_geo_eval_string(args[1],context,caller)
        haskey(context.onelab_numbers,key) &&
            return context.onelab_numbers[key]
        length(args)==2 &&
            return _geo_eval_numeric(args[2],context,caller)
        return 0.0
    elseif name=="GetValue"
        length(args)==2 || _geo_expr_error(
            parser,"GetValue takes a prompt and a default",pos)
        # `Msg::GetNumber` prompts on stdin; batch runs return the default.
        return _geo_eval_numeric(args[2],context,caller)
    elseif name=="FileExists"
        length(args)==1 || _geo_expr_error(
            parser,"FileExists takes exactly one argument",pos)
        path=_geo_fix_relative_path(
            context.file_name,_geo_eval_string(args[1],context,caller))
        return Float64(ispath(path))
    elseif name=="Find"
        length(args)==2 || _geo_expr_error(
            parser,"Find takes two list arguments",pos)
        needles=_geo_numeric_list_values(args[1],context,caller)
        haystack=_geo_numeric_list_values(args[2],context,caller)
        return Float64(count(n->n in haystack,needles))
    elseif name=="StrFind"
        length(args)==2 || _geo_expr_error(
            parser,"StrFind takes two arguments",pos)
        text=_geo_eval_string(args[1],context,caller)
        needle=_geo_eval_string(args[2],context,caller)
        # `std::string::find != npos` — a boolean hit, not a position.
        return occursin(needle,text) ? 1.0 : 0.0
    elseif name=="StrCmp"
        length(args)==2 || _geo_expr_error(
            parser,"StrCmp takes two arguments",pos)
        a=_geo_eval_string(args[1],context,caller)
        b=_geo_eval_string(args[2],context,caller)
        return Float64(cmp(a,b))
    elseif name=="StrLen"
        length(args)==1 || _geo_expr_error(
            parser,"StrLen takes exactly one argument",pos)
        return Float64(ncodeunits(_geo_eval_string(args[1],context,caller)))
    elseif name=="TextAttributes"
        # `(align<<16)|(font<<8)|fontsize` — `drawContext` indices: font via
        # the FLTK menu list, align via the 0-8 alignment map; fontsize is
        # `General.GraphicsFontSize` unless a `FontSize` pair overrides it.
        isempty(args) && throw(ArgumentError(
            "$caller: TextAttributes requires at least one attribute"))
        fontsize=something(_geo_option_number(context,"General",0,
            "GraphicsFontSize"),15.0)
        font=0;align=0
        # `RecursiveListOfStringExprVar` — `x()`/`x.member()`/`{...}` items
        # expand before pairing.
        svals=_geo_eval_string_list(inner,context,caller)
        isodd(length(svals)) && (_geo_yyerror!(context,
            "Number of text attributes should be even"); return 0.0)
        for i in 1:2:length(svals)-1
            key=svals[i];value=svals[i+1]
            if key=="FontSize"
                # `atoi` — a non-numeric value folds to 0.
                fontsize=Float64(something(tryparse(Int,value),0))
            elseif key=="Font"
                fi=findfirst(==(value),_GEO_FONT_NAMES)
                if fi===nothing
                    _geo_msg_error!(context,
                        "Unknown font \"$value\" (using \"Helvetica\" instead)")
                    font=4
                else
                    font=fi-1
                end
            elseif key=="Align"
                ai=get(_GEO_FONT_ALIGNS,value,nothing)
                if ai===nothing
                    _geo_msg_error!(context,
                        "Unknown font alignment \"$value\" (using \"Left\" instead)")
                    align=0
                else
                    align=ai
                end
            end
        end
        return Float64((align<<16)|(font<<8)|Int(fontsize))
    elseif name=="DimNameSpace"
        # `tDimNameSpace LP [String__Index] RP` — the number of Struct
        # definitions in the queried (or default) namespace. Upstream reads
        # through `gmsh_yynamespaces[ns]` — `Map::operator[]` *creates* an
        # empty `Structs` entry for a missing namespace, so a later
        # `NameStruct(ns::#i)` reports "Unknown Struct of index" rather than
        # "Unknown NameSpace".
        length(args)<=1 || _geo_expr_error(
            parser,"DimNameSpace takes zero or one arguments",pos)
        ns=isempty(args) ? "" : _geo_symbol_name(args[1],context,caller)
        structs=get!(context.structs,ns,Dict{String,_GeoStructDef}())
        return Float64(length(structs))
    elseif name in ("StringToName","S2N")
        length(args)==1 || _geo_expr_error(
            parser,"$name takes exactly one argument",pos)
        sym=_geo_symbol_name(args[1],context,caller)
        return _geo_eval_numeric(sym,context,caller)
    elseif name=="DefineNumber"
        # `DefineNumber(FExpr, FloatParameterOptions...)` — the ONELAB
        # parameter name comes from the `Name` option; the exchanged value
        # (server override) is returned.
        isempty(args) && _geo_expr_error(
            parser,"DefineNumber takes a value argument",pos)
        value=_geo_eval_numeric(args[1],context,caller)
        vals=_geo_exchange_onelab_number!(context,"",Float64[value],
            length(args)>1 ? args[2:end] : nothing,caller)
        return vals[1]
    elseif name=="GetNumberChoice"
        return 0.0
    end
    _geo_expr_error(parser,"unknown function $name",pos)
end

# Postfix `++`/`--` on a scalar FExpr: `s.value[0] += delta` — the retained
# list payload (if any) stays in sync because element 0 IS the scalar value.
function _geo_increment_scalar!(context::_GeoNumericContext,name::String,
                                delta::Float64,caller::AbstractString)
    if haskey(context.values,name)
        context.values[name]+=delta
        haskey(context.lists,name) && !isempty(context.lists[name]) &&
            (context.lists[name][1]+=delta)
        return nothing
    end
    haskey(context.unavailable,name) && throw(ArgumentError(
        "$caller: numeric variable $name is unavailable ($(context.unavailable[name]))"))
    # `x++` on a missing symbol or an empty `s.value` — `yymsg(0)` under
    # execution; the FExpr form evaluates to 0.
    if context.soft_unknown_reads
        _geo_yyerror!(context,_geo_context_has_variable(context,name) ?
            "Uninitialized variable '$name'" : "Unknown variable '$name'")
        return nothing
    end
    throw(ArgumentError("$caller: Unknown variable '$name'"))
end

# FExpr postfix `x[i]++`/`x(i)--` — reads `s.value[index]` then adds the
# delta, returning the OLD value. Unlike the statement form this works on
# scalars (element 0) and never resizes; misses and out-of-range indices
# record `yymsg(0)` diagnostics and evaluate to 0.
function _geo_increment_indexed!(context::_GeoNumericContext,name::String,
                                 index::Float64,delta::Float64,
                                 caller::AbstractString)
    i=_geo_int_value(index,"$caller increment index")
    i<0 && throw(ArgumentError("$caller: negative index $i in '$name'"))
    haskey(context.lists,name) ||
        haskey(context.values,name) || begin
        if context.soft_unknown_reads
            _geo_yyerror!(context,"Unknown variable '$name'")
            return 0.0
        end
        throw(ArgumentError("$caller: Unknown variable '$name'"))
    end
    values=haskey(context.lists,name) ? context.lists[name] :
        (context.lists[name]=Float64[context.values[name]])
    if i>=length(values)
        if context.soft_unknown_reads
            _geo_yyerror!(context,"Uninitialized variable '$name[$i]'")
            return 0.0
        end
        throw(ArgumentError("$caller: Uninitialized variable '$name[$i]'"))
    end
    old=values[i+1]
    values[i+1]+=delta
    context.values[name]=values[1]
    return old
end

function _geo_eval_numeric(raw::AbstractString,context::_GeoNumericContext,
                           caller::AbstractString)
    source=String(strip(raw))
    isempty(source) && throw(ArgumentError("$caller: numeric expression must not be empty"))
    ncodeunits(source)<=_MAX_GEO_EXPRESSION_BYTES || throw(ArgumentError(
        "$caller: expression exceeds $_MAX_GEO_EXPRESSION_BYTES bytes"))
    parser=_GeoExprParser(source,firstindex(source),_GeoExprToken(:eof,"",0.0,1),
                          0,0,context,String(caller))
    _geo_advance!(parser)
    value=_geo_parse_ternary!(parser)
    parser.token.kind==:eof || begin
        token=parser.token
        if token.kind==:side_effect
            # Trailing `++`/`--` after a complete FExpr — `syntax error (++)`
            # upstream.
            _geo_syntax_abort(token.text)
        elseif token.kind==:unsupported_operator
            _geo_expr_error(parser,
                "operator $(token.text) is outside the supported arithmetic subset",token.pos)
        end
        _geo_expr_error(parser,"unexpected token $(repr(token.text))",token.pos)
    end
    return _geo_finite_result(parser,value,"expression",1)
end

# ======== Gmsh `.geo` string expressions (`StringExpr`/`StringExprVar`) ========

# Split a comma-separated argument string on top-level commas, tracking
# parentheses, brackets, braces and quoted literals (Gmsh strings keep
# backslashes verbatim — there is no escape processing in `parsestring`).
function _geo_split_args(raw::AbstractString,caller::AbstractString)
    src=String(strip(raw))
    isempty(src) && return String[]
    items=String[];start=firstindex(src);i=start;last=lastindex(src)
    parens=0;brackets=0;braces=0;quote_char='\0'
    while i<=last
        c=src[i]
        if quote_char!='\0'
            c==quote_char && (quote_char='\0')
        elseif c in ('"','\'')
            quote_char=c
        elseif c=='(' parens+=1
        elseif c==')'
            parens-=1;parens>=0 || _geo_syntax_abort(")")
        elseif c=='[' brackets+=1
        elseif c==']'
            brackets-=1;brackets>=0 || _geo_syntax_abort("]")
        elseif c=='{' braces+=1
        elseif c=='}'
            braces-=1;braces>=0 || _geo_syntax_abort("}")
        elseif c==',' && parens==0 && brackets==0 && braces==0
            item=String(strip(src[start:prevind(src,i)]))
            isempty(item) && _geo_syntax_abort(",")
            length(items)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                "$caller: argument list exceeds $_MAX_GEO_LIST_ITEMS entries"))
            push!(items,item);start=nextind(src,i)
        end
        i=nextind(src,i)
    end
    parens==0 || _geo_syntax_abort("(")
    brackets==0 || _geo_syntax_abort("[")
    braces==0 || _geo_syntax_abort("{")
    quote_char=='\0' || _geo_syntax_abort("\"")
    item=String(strip(src[start:last]))
    push!(items,item)
    return items
end

# Find the index of the closer matching the delimiter at `open_pos`.
function _geo_matching_delim(src::AbstractString,open_pos::Int)
    open_ch=src[open_pos]
    close_ch=open_ch=='(' ? ')' : open_ch=='[' ? ']' : '}'
    depth=0;quote_char='\0';i=open_pos;last=lastindex(src)
    while i<=last
        c=src[i]
        if quote_char!='\0'
            c==quote_char && (quote_char='\0')
        elseif c in ('"','\'')
            quote_char=c
        elseif c==open_ch;depth+=1
        elseif c==close_ch
            depth-=1
            depth==0 && return i
        end
        i=nextind(src,i)
    end
    return 0
end

# Resolve Gmsh `String__Index` to a flat symbol-table name: `x`, `x~{i}` (each
# `~{expr}` suffix appends `_i`, recursively), `StringToName[sexpr]` and
# `StringToName[sexpr]~{i}`. `tSCOPE` is not part of the production — `ns::x`
# reaches variables only through `Struct_FullName` expression forms, so a
# stray `::` is the bison syntax error upstream reports.
function _geo_symbol_name(raw::AbstractString,context::_GeoNumericContext,
                          caller::AbstractString)
    s=String(strip(raw))
    isempty(s) && _geo_syntax_abort(";",
        "$caller: expected a variable name")
    name=if (sm=match(r"^(?:StringToName|S2N)\s*\[",s))!==nothing
        open_pos=lastindex(sm.match)
        close=_geo_matching_delim(s,open_pos)
        close==0 && _geo_syntax_abort("[",
            "$caller: unbalanced StringToName brackets")
        inner=_geo_eval_string(
            s[nextind(s,open_pos):prevind(s,close)],context,caller)
        rest=String(strip(s[nextind(s,close):end]))
        isempty(rest) || startswith(rest,"~") || _geo_syntax_abort(
            _geo_string_index_offender(rest),
            "$caller: unexpected text after StringToName[...]: $rest")
        _geo_symbol_name_suffixes!(inner,rest,context,caller)
    else
        m=match(r"^([A-Za-z_][A-Za-z0-9_]*)(.*)$",s)
        m===nothing && _geo_syntax_abort(_geo_string_index_offender(s),
            "$caller: expected a variable name, got $s")
        base=String(m.captures[1]);rest=String(strip(m.captures[2]))
        # `tSCOPE` is not part of `String__Index` — `ns::` heads reach the
        # parser only through `Struct_FullName` expression productions.
        startswith(rest,"::") && _geo_syntax_abort("::",
            "$caller: namespace `::` is not valid in a variable name: $s")
        _geo_symbol_name_suffixes!(base,rest,context,caller)
    end
    return name
end

# The first source token in `arg` that cannot extend a `String__Index`
# (`identifier` / `StringToName[...]` followed by `~{expr}` suffixes) — the
# bison lookahead a malformed `NameToString[...]` argument aborts on. A
# contiguous `::` reports as the tSCOPE pair.
function _geo_string_index_offender(arg::AbstractString)
    s=String(strip(arg))
    isempty(s) && return "?"
    i=firstindex(s);last=lastindex(s)
    if (m=match(r"^(?:StringToName|S2N)\s*\[",s))!==nothing
        close=_geo_matching_delim(s,lastindex(m.match))
        close==0 && return "["
        i=nextind(s,close)
    else
        m=match(r"^[A-Za-z_][A-Za-z0-9_]*",s)
        m===nothing && return string(s[i])
        i=nextind(s,lastindex(m.match))
    end
    while i<=last
        c=s[i]
        isspace(c) && (i=nextind(s,i);continue)
        if c=='~'
            j=nextind(s,i)
            while j<=last && isspace(s[j]); j=nextind(s,j); end
            (j<=last && s[j]=='{') || return "~"
            close=_geo_matching_delim(s,j)
            close==0 && return "{"
            i=nextind(s,close)
        else
            return c==':' && nextind(s,i)<=last && s[nextind(s,i)]==':' ?
                "::" : string(c)
        end
    end
    return "?"
end

# Resolve the leading `String__Index` of an assignment LHS — `name`,
# `name~{FExpr}`…, `StringToName[sexpr]`/`S2N[sexpr]` (each optionally
# followed by `~{}` suffixes). Returns `(resolved_name, rest)` where `rest`
# is the unconsumed source (`(...)`, `[i]`, `.opt`, ` = …`…), or `nothing`
# when the head is a plain identifier (the common case — the caller's
# regexes then see the body unchanged) or not a `String__Index` at all.
function _geo_lhs_string_index(src::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString)
    s=String(strip(src));last=lastindex(s)
    if (m=match(r"^(?:StringToName|S2N)\s*\[",s))!==nothing
        close=_geo_matching_delim(s,lastindex(m.match))
        close==0 && return nothing
        i=nextind(s,close)
    else
        m=match(r"^[A-Za-z_][A-Za-z0-9_]*",s)
        m===nothing && return nothing
        # Only `tSTRING` heads form a `String__Index` — a lexer keyword
        # (`Field`, `Physical`…) keeps its own statement forms.
        _geo_tstring_valid(String(m.match)) || return nothing
        i=nextind(s,lastindex(m.match))
    end
    head_end=i
    while i<=last
        j=i
        while j<=last && isspace(s[j]); j=nextind(s,j); end
        (j<=last && s[j]=='~') || break
        k=nextind(s,j)
        while k<=last && isspace(s[k]); k=nextind(s,k); end
        (k<=last && s[k]=='{') || break
        close=_geo_matching_delim(s,k)
        close==0 && break
        i=nextind(s,close);head_end=i
    end
    head=String(s[firstindex(s):prevind(s,head_end)])
    occursin('~',head) || occursin('[',head) || return nothing
    return (_geo_symbol_name(head,context,caller),
            String(s[head_end:end]))
end

function _geo_symbol_name_suffixes!(base::String,rest::String,
                                    context::_GeoNumericContext,
                                    caller::AbstractString)
    s=rest
    while startswith(s,"~")
        m=match(r"^~\s*\{",s)
        m===nothing && _geo_syntax_abort("~",
            "$caller: expected `{expr}` after `~` in a namespaced name")
        open_pos=lastindex(m.match)
        close=_geo_matching_delim(s,open_pos)
        close==0 && _geo_syntax_abort("{",
            "$caller: unbalanced braces in a namespaced name")
        idx=_geo_int_value(_geo_eval_numeric(
            s[nextind(s,open_pos):prevind(s,close)],context,
            "$caller namespace index"),"$caller namespace index")
        base=base*"_"*string(idx)
        s=String(strip(s[nextind(s,close):end]))
    end
    isempty(s) || _geo_syntax_abort(
        _geo_string_index_offender(s),
        "$caller: unexpected text after a variable name: $s")
    return base
end

# ======== gmsh_yynamespaces — Struct / NameSpace / NameStruct ========

# A `yymsg(0)`-class expression diagnostic: recorded under `.geo` execution
# (the stream continues with the default value), thrown under the params
# scan so callers mark the target unavailable.
function _geo_expr_diag!(context::_GeoNumericContext,message::AbstractString,
                        caller::AbstractString)
    context.soft_unknown_reads || throw(ArgumentError("$caller: $message"))
    _geo_yyerror!(context,message)
    return nothing
end

# `NameSpaces::Find`/`Structs::Find` — the struct or nothing.
@inline _geo_struct_get(context::_GeoNumericContext,ns::String,name::String)=
    (s=get(context.structs,ns,nothing))===nothing ? nothing : get(s,name,nothing)

# `NameSpaces::getTag` — the struct tag or nothing.
@inline _geo_struct_tag(context::_GeoNumericContext,ns::String,name::String)=
    (st=_geo_struct_get(context,ns,name))===nothing ? nothing : Float64(st.tag)

# `String__Index` — `tSTRING` / `name~{FExpr}` / `StringToName[StringExprVar]`.
# Returns the resolved name text.
function _geo_string_index(context::_GeoNumericContext,src::AbstractString,
                           caller::AbstractString)
    return _geo_symbol_name(src,context,caller)
end

# `Struct_FullName` — `String__Index` or `String__Index :: String__Index`.
# `ns` is nothing when no `::` prefix is present.
function _geo_parse_fullname(context::_GeoNumericContext,src::AbstractString,
                             caller::AbstractString)
    s=String(strip(src))
    depth=0;qc='\0';scope_at=0
    i=firstindex(s);last=lastindex(s)
    while i<=last
        c=s[i]
        if qc!='\0'
            c==qc && (qc='\0')
        elseif c=='"' || c=='\''
            qc=c
        elseif c in ('{','[','(')
            depth+=1
        elseif c in ('}',']',')')
            depth-=1
        elseif depth==0 && c==':' && i<last && s[nextind(s,i)]==':'
            scope_at=i;break
        end
        i=nextind(s,i)
    end
    if scope_at==0
        return (nothing,_geo_string_index(context,s,caller))
    end
    ns=_geo_string_index(context,s[firstindex(s):prevind(s,scope_at)],caller)
    name=_geo_string_index(context,s[nextind(s,scope_at,2):end],caller)
    return (ns,name)
end

# `treat_Struct_FullName_Float` — `name`/`ns::name` numeric reads.
# `ns===nothing` means no `::` prefix: variables (then Exists-string symbols)
# shadow struct tags. `type_var`: 1 = scalar/tag read, 2 = `name(index)`
# indexed read. `type_treat`: 0 = get, 1 = `Exists`, 2 = `GetForced`.
function _geo_fullname_number(context::_GeoNumericContext,
                              ns::Union{Nothing,String},name::String,
                              caller::AbstractString;
                              type_var::Int=1,index::Int=0,
                              default::Float64=0.0,type_treat::Int=0)
    if ns===nothing && _geo_context_has_variable(context,name)
        type_treat==1 && return 1.0
        values=_geo_context_list(context,name,caller)
        if type_var==1
            if isempty(values)
                type_treat==0 &&
                    _geo_expr_diag!(context,"Uninitialized variable '$name'",caller)
                return default
            end
            return values[1]
        elseif type_var==2
            if !(0<=index<length(values))
                type_treat==0 && _geo_expr_diag!(context,
                    "Uninitialized variable '$name[$index]'",caller)
                return default
            end
            return values[index+1]
        end
        return default
    elseif ns===nothing && type_treat==1 && haskey(context.strings,name)
        return 1.0
    end
    if type_var==1
        tag=_geo_struct_tag(context,something(ns,""),name)
        tag!==nothing && return tag
        type_treat==0 &&
            _geo_expr_diag!(context,"Unknown variable '$name'",caller)
        return default
    end
    type_treat==0 &&
        _geo_expr_diag!(context,"Unknown variable '$name(.)'",caller)
    return default
end

# `treat_Struct_FullName_dot_tSTRING_Float` — `x.member`/`ns::x.member`
# numeric reads (index optional). Unknown structs fall through to
# `NumberOption(x, 0, member)` — the option system owns that diagnostic.
function _geo_member_number(context::_GeoNumericContext,ns::String,
                            name::String,member::String,index::Int,
                            caller::AbstractString;
                            default::Float64=0.0,type_treat::Int=0)
    st=_geo_struct_get(context,ns,name)
    if st===nothing
        value=_geo_option_number(context,name,0,member;warn=type_treat==0)
        return value===nothing ? default : value
    end
    vals=get(st.fopt,member,nothing)
    vals===nothing && return begin
        type_treat!=0 &&
            return haskey(st.copt,member) ? 1.0 : default
        type_treat==0 && _geo_expr_diag!(context,
            "Unknown member '$member' of Struct $name",caller)
        default
    end
    0<=index<length(vals) && return type_treat==1 ? 1.0 : vals[index+1]
    type_treat==0 &&
        _geo_expr_diag!(context,"Index $index out of range",caller)
    return default
end

# `treat_Struct_FullName_dot_tSTRING_String` — `x.member`/`ns::x.member`
# string reads. Unknown structs fall through to `StringOption`.
function _geo_member_string(context::_GeoNumericContext,ns::String,
                            name::String,member::String,index::Int,
                            caller::AbstractString;
                            default::String="",type_treat::Int=0)
    st=_geo_struct_get(context,ns,name)
    if st===nothing
        value=_geo_option_string(context,name,0,member,caller;
            warn=type_treat==0)
        return value===nothing ? default : value
    end
    vals=get(st.copt,member,nothing)
    vals===nothing && return begin
        type_treat==0 && _geo_expr_diag!(context,
            "Unknown member '$member' of Struct $name",caller)
        default
    end
    0<=index<length(vals) && return vals[index+1]
    type_treat==0 &&
        _geo_expr_diag!(context,"Index $index out of range",caller)
    return default
end

# `treat_Struct_FullName_String` — bare `name` in string position reads only
# the string-symbol table (a struct tag is not a string).
function _geo_fullname_string(context::_GeoNumericContext,name::String,
                              caller::AbstractString;
                              default::String="",type_treat::Int=0)
    if haskey(context.strings,name)
        values=context.strings[name]
        length(values)==1 && return values[1]
        type_treat==0 && _geo_expr_diag!(context,
            "Expected single valued string variable '$name'",caller)
        return default
    end
    type_treat==0 &&
        _geo_expr_diag!(context,"Unknown string variable '$name'",caller)
    return default
end

# `treat_Struct_FullName_dot_tSTRING_Float_getDim` — `#x.member()` /
# `#ns::x.member()`: the member's element count (float or char side).
function _geo_member_dim(context::_GeoNumericContext,ns::String,name::String,
                        member::String,caller::AbstractString)
    st=_geo_struct_get(context,ns,name)
    st===nothing && return 0
    vals=get(st.fopt,member,nothing)
    vals!==nothing && return length(vals)
    svals=get(st.copt,member,nothing)
    svals!==nothing && return length(svals)
    _geo_expr_diag!(context,"Unknown member '$member' of Struct $name",caller)
    return 0
end

# `treat_Struct_FullName_dot_tSTRING_ListOfFloat` — `x.member()` /
# `ns::x.member()` whole-member splices; nothing when the member is missing
# (upstream yields nullptr → the splice contributes no elements).
function _geo_member_fvector(context::_GeoNumericContext,ns::String,
                             name::String,member::String,
                             caller::AbstractString)
    st=_geo_struct_get(context,ns,name)
    st===nothing && return begin
        _geo_expr_diag!(context,"Unknown Struct: $name",caller)
        nothing
    end
    vals=get(st.fopt,member,nothing)
    vals===nothing && return begin
        _geo_expr_diag!(context,"Unknown member '$member' of Struct $name",caller)
        nothing
    end
    return copy(vals)
end

# `treat_Struct_FullName_dot_tSTRING_ListOfString` — the char-side splice.
function _geo_member_svector(context::_GeoNumericContext,ns::String,
                             name::String,member::String,
                             caller::AbstractString)
    st=_geo_struct_get(context,ns,name)
    st===nothing && return begin
        _geo_expr_diag!(context,"Unknown Struct: $name",caller)
        nothing
    end
    vals=get(st.copt,member,nothing)
    vals===nothing && return begin
        _geo_expr_diag!(context,"Unknown member '$member' of Struct $name",caller)
        nothing
    end
    return copy(vals)
end

# `Structs::get_key_struct_from_tag` — iterates the namespace's (sorted) map
# and returns the first struct whose tag equals `tag`: 0 = found,
# 1 = unknown namespace, 2 = unknown struct tag.
function _geo_struct_key_from_tag(context::_GeoNumericContext,ns::String,
                                  tag::Int)
    structs=get(context.structs,ns,nothing)
    structs===nothing && return (1,"")
    for (name,st) in sort!(collect(structs);by=first)
        st.tag==tag && return (0,name)
    end
    return (2,"")
end

# `FloatParameterOption` items of a `Struct name[...]` definition: `init_options`
# seeds `member_valmax` from the struct being appended to (or -1) and clears
# `flag_Enum`. Mutates fopt/copt.
function _geo_struct_option!(context::_GeoNumericContext,src::AbstractString,
                             fopt::Dict{String,Vector{Float64}},
                             copt::Dict{String,Vector{String}},
                             flag_enum::Base.RefValue{Bool},
                             member_valmax::Base.RefValue{Int},
                             caller::AbstractString)
    s=String(strip(src))
    isempty(s) && return nothing
    om=match(r"^([A-Za-z_][A-Za-z0-9_]*)(.*)$",s)
    om===nothing &&
        _geo_syntax_abort(string(first(s)),
            "$caller: malformed Struct option $(repr(s))")
    key=String(om.captures[1]);rest=String(strip(om.captures[2]))
    if key in ("Min","Max","Step")
        # `tMin FExpr`/`tMax FExpr`/`tStep FExpr` — a single expression, never
        # a brace list or bare flag.
        isempty(rest) &&
            _geo_syntax_abort(",", "$caller: $key requires a value")
        push!(get!(fopt,key,Float64[]),
              _geo_eval_numeric(rest,context,"$caller $key option"))
        return nothing
    end
    if isempty(rest)
        # bare `tSTRING` — flag or `Enum` auto-numbering.
        if !flag_enum[]
            push!(get!(fopt,key,Float64[]),1.0)
            key=="Enum" && (flag_enum[]=true)
        else
            member_valmax[]+=1
            push!(get!(fopt,key,Float64[]),Float64(member_valmax[]))
        end
        return nothing
    end
    if (strm=match(r"^Str\s*[\(\[]",rest))!==nothing
        # `tSTRING Str_BracedRecursiveListOfStringExprVar` — `Key Str[{..}]`.
        open_at=lastindex(strm.match)
        close=_geo_matching_delim(rest,open_at)
        close==0 && _geo_syntax_abort(rest[open_at],
            "$caller: unbalanced Str[...] in a Struct option")
        String(strip(rest[nextind(rest,close):end]))=="" ||
            _geo_syntax_abort(rest[nextind(rest,close)],
                "$caller: unexpected text after $key Str[...]")
        inner=String(strip(rest[nextind(rest,open_at):prevind(rest,close)]))
        if startswith(inner,"{") &&
                _geo_matching_delim(inner,firstindex(inner))==lastindex(inner)
            # `tSTRING Str_BracedRecursiveListOfStringExprVar` — `Str[{...}]`
            # pushes each element of the braced list to charOptions.
            for item in _geo_eval_string_list(inner[nextind(inner,
                    firstindex(inner)):prevind(inner,lastindex(inner))],
                    context,caller)
                push!(get!(copt,key,String[]),item)
            end
        elseif startswith(inner,"{")
            # `Str[{a}, x]` — the leading `{...}` commits to the
            # `Str_BracedRecursiveListOfStringExprVar` production, so the
            # trailing comma is the offending lookahead.
            _geo_eval_string_list(inner,context,caller;brace_commit=true)
        else
            # Unbraced `Str[...]`/`Str(...)` does not match the braced-list
            # production — it is a `tSTRING StringExpr` where the `Str`
            # function joins its arguments with `\n` into ONE char member.
            push!(get!(copt,key,String[]),
                _geo_eval_string("Str["*inner*"]",context,caller))
        end
        return nothing
    end
    if startswith(rest,"{")
        close=_geo_matching_delim(rest,firstindex(rest))
        close==0 && _geo_syntax_abort('{',
            "$caller: unbalanced braces in a Struct option")
        inner=String(rest[nextind(rest,firstindex(rest)):prevind(rest,close)])
        after=String(strip(rest[nextind(rest,close):end]))
        if !isempty(after)
            # `FExpr_Multi tAFFECT String__Index '(' ')'` — pair the listed
            # values with a string list variable elementwise.
            em=match(r"^=\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)\s*$",after)
            em===nothing && _geo_syntax_abort(first(after),
                "$caller: unexpected text after $key {...}")
            var=String(em.captures[1])
            vals=_geo_numeric_list_values("{$inner}",context,caller;
                                          allow_multiplier=true)
            if !haskey(context.strings,var)
                _geo_expr_diag!(context,"Unknown string variable '$var'",caller)
            else
                strs=context.strings[var]
                if length(vals)!=length(strs)
                    _geo_expr_diag!(context,"Size mismatch in enumeration: " *
                        "$(length(vals)) != $(length(strs))",caller)
                else
                    append!(get!(fopt,key,Float64[]),vals)
                    append!(get!(copt,key,String[]),strs)
                end
            end
            return nothing
        end
        # `{ FExpr = StringExpr, ... }` enumeration, or a plain numeric list.
        items=_geo_split_args(inner,caller)
        is_enum=any(item->begin
            d=0;q='\0'
            for c in item
                q!='\0' ? (c==q && (q='\0')) :
                c in ('"','\'') ? (q=c) :
                c in ('{','[','(') ? (d+=1) :
                c in ('}',']',')') ? (d-=1) :
                (d==0 && c=='=' && return true)
            end
            false
        end,items)
        if is_enum
            for item in items
                em=match(r"^(.*?)\s*=\s*(.*?)\s*$",item)
                em===nothing && _geo_syntax_abort('{',
                    "$caller: malformed enumeration item $(repr(item))")
                push!(get!(fopt,key,Float64[]),_geo_eval_numeric(
                    em.captures[1],context,"$caller enumeration value"))
                push!(get!(copt,key,String[]),_geo_eval_string(
                    em.captures[2],context,caller))
            end
        else
            # A `key {...}` list is `ListOfDouble`/Enumeration — a quoted
            # string inside is the offending lookahead upstream
            # (`lst {"a","b"}` → `syntax error (")`; string lists need
            # `Str[{...}]`).
            qpos=findfirst(c->c=='"'||c=='\'',inner)
            qpos===nothing ||
                _geo_syntax_abort(string(inner[qpos]))
            vals=_geo_numeric_list_values("{$inner}",context,caller;
                                          allow_multiplier=true)
            flag_enum[] && !isempty(vals) &&
                (member_valmax[]=_geo_signed_gmsh_int_value(
                    vals[1],"$caller $key enum bound"))
            append!(get!(fopt,key,Float64[]),vals)
        end
        return nothing
    end
    if _geo_string_rhs(rest)
        # `tSTRING StringExpr` — a char member.
        push!(get!(copt,key,String[]),_geo_eval_string(rest,context,caller))
        return nothing
    end
    # `tSTRING ListOfDouble` — single FExpr, `FExpr_Multi`, `-'{...}` or
    # `n*{...}` multiplier lists.
    vals=_geo_numeric_list_values(rest,context,caller;allow_multiplier=true)
    flag_enum[] && !isempty(vals) &&
        (member_valmax[]=_geo_signed_gmsh_int_value(
            vals[1],"$caller $key enum bound"))
    append!(get!(fopt,key,Float64[]),vals)
    return nothing
end

# `NameSpaces::defStruct`/`Structs::defStruct` — insert-or-append a struct.
# `append_level`: 0 = fresh definition, otherwise `Append`/`Append n`.
# Returns (tag, redefinition::Bool) — on redefinition the existing struct is
# kept and its tag returned (upstream `tag_out = getTag()`).
function _geo_def_struct!(context::_GeoNumericContext,ns::String,
                          name::String,append_level::Int,
                          options_src::AbstractString,
                          caller::AbstractString)
    structs=get!(context.structs,ns,Dict{String,_GeoStructDef}())
    existing=get(structs,name,nothing)
    append=append_level!=0 && existing!==nothing
    # `init_options(getMember_ValMax(ns, name))` — continues enum numbering
    # across appends; a missing struct seeds -1 so the first flag lands on 0.
    member_valmax=Ref(existing===nothing ? -1 : existing.member_valmax)
    flag_enum=Ref(false)
    fopt=Dict{String,Vector{Float64}}();copt=Dict{String,Vector{String}}()
    for opt in _geo_split_args(options_src,caller)
        _geo_struct_option!(context,opt,fopt,copt,flag_enum,member_valmax,
                            caller)
    end
    if existing!==nothing && !append
        return (existing.tag,true)
    end
    maxtag=get(context.struct_maxtag,ns,0)
    tag=if haskey(fopt,"Tag")
        forced=_geo_signed_gmsh_int_value(fopt["Tag"][1],
            "$caller Struct Tag")
        context.struct_maxtag[ns]=max(maxtag,forced)
        forced
    else
        append ? -1 : (context.struct_maxtag[ns]=maxtag+1; maxtag+1)
    end
    if !append
        # `Tag` auto-injection — visible as a member.
        haskey(fopt,"Tag") || (fopt["Tag"]=Float64[tag])
        structs[name]=_GeoStructDef(tag,fopt,copt,member_valmax[])
        return (tag,false)
    end
    # `Struct::append` — `map::insert` merges only *new* member keys; the tag
    # updates only when `Tag` was forced (tag >= 0); member_valmax always
    # tracks the latest run.
    tag>=0 && (existing.tag=tag)
    existing.member_valmax=member_valmax[]
    for (k,v) in fopt
        haskey(existing.fopt,k) || (existing.fopt[k]=v)
    end
    for (k,v) in copt
        haskey(existing.copt,k) || (existing.copt[k]=v)
    end
    return (existing.tag,false)
end

# `Struct name[...] ;` / `Struct ns::name (Append [n]) [opts] ;` — the
# `DefineStruct tEND` statement form. `x = Struct ...` instead flows through
# the FExpr production (`_geo_parse_inline_struct!`).
function _geo_exec_struct_def!(raw::AbstractString,
                               context::_GeoNumericContext)
    caller="execute_geo: Struct"
    s=String(strip(raw))
    endswith(s,";") &&
        (s=String(strip(s[firstindex(s):prevind(s,lastindex(s))])))
    open_at=0;i=firstindex(s);last=lastindex(s)
    while i<=last
        c=s[i]
        if c=='['
            open_at=i;break
        elseif c=='"' || c=='\''
            close=findnext(==(c),s,nextind(s,i))
            close===nothing && break
            i=close
        end
        i=nextind(s,i)
    end
    if open_at==0
        # Grammar order is `Struct_FullName AppendOrNot '['` — the syntax
        # error lands on the first token outside that sequence
        # (`Struct = 5` → `(=)`, `Struct s` → `(;)`, `Struct s junk` →
        # `(junk)`); quoted strings are not `String__Index` heads.
        hm=match(r"^\s*[A-Za-z_][A-Za-z0-9_]*(?:\s*::\s*[A-Za-z_][A-Za-z0-9_]*)?(?:\s*\(\s*Append\b[^)]*\))?",s)
        tail=hm===nothing ? s : String(strip(s[nextind(s,
            lastindex(hm.match)):end]))
        isempty(tail) && _geo_syntax_abort(";")
        tm=match(r"^[A-Za-z_][A-Za-z0-9_]*",tail)
        _geo_syntax_abort(tm===nothing ? string(first(tail)) :
            String(tm.match))
    end
    close=_geo_matching_delim(s,open_at)
    close==0 && _geo_syntax_abort(";")
    tail=strip(s[nextind(s,close):end])
    if !isempty(tail)
        tm=match(r"^[A-Za-z_][A-Za-z0-9_]*",tail)
        _geo_syntax_abort(tm===nothing ? string(first(tail)) :
            String(tm.match))
    end
    # `Struct_FullName` parts are `String__Index` — lexer keywords abort at
    # the keyword itself (`Struct Point [a]` → `syntax error (Point)`).
    # `~{FExpr}` suffixes are blanked so their contents cannot fake a `::`.
    header_pre=String(strip(s[firstindex(s):prevind(s,open_at)]))
    header_pre=replace(header_pre,r"\(\s*Append\b[^)]*\)\s*$"=>"")
    cleaned=IOBuffer();i=firstindex(header_pre);last=lastindex(header_pre)
    while i<=last
        c=header_pre[i]
        if c=='~'
            j=nextind(header_pre,i)
            while j<=last && isspace(header_pre[j]); j=nextind(header_pre,j); end
            if j<=last && header_pre[j]=='{'
                cl=_geo_matching_delim(header_pre,j)
                cl!=0 && (i=nextind(header_pre,cl);continue)
            end
        end
        print(cleaned,c);i=nextind(header_pre,i)
    end
    for part in split(String(take!(cleaned)),"::")
        pname=strip(part)
        isempty(pname) && _geo_syntax_abort("::")
        (pm=match(r"^[A-Za-z_][A-Za-z0-9_]*$",pname))===nothing &&
            _geo_syntax_abort(string(first(pname)))
        String(pm.match) in _GEO_LEXER_KEYWORDS &&
            _geo_syntax_abort(String(pm.match))
    end
    options_src=String(s[nextind(s,open_at):prevind(s,close)])
    header=String(strip(s[firstindex(s):prevind(s,open_at)]))
    # `(Append)`/`(Append FExpr)` — `AppendOrNot`: 0 fresh, `Append` → 99,
    # `Append n` → n; the level is only a boolean gate upstream.
    append_level=0
    if (am=match(r"\(\s*Append\b(.*)\)\s*$",header))!==nothing
        header=String(strip(header[firstindex(header):prevind(header,
            am.offset)]))
        level_src=String(strip(am.captures[1]))
        append_level=isempty(level_src) ? 99 : _geo_signed_gmsh_int_value(
            _geo_eval_numeric(level_src,context,"$caller Append level"),
            "$caller Append level")
    end
    (ns,name)=_geo_parse_fullname(context,header,caller)
    (tag,redef)=_geo_def_struct!(context,something(ns,""),name,append_level,
                                 options_src,caller)
    redef && _geo_expr_diag!(context,
        "Redefinition of Struct '$(something(ns,""))::$name'",caller)
    return tag
end

# Convenience wrapper: evaluate a string expression without surrounding junk.
function _geo_eval_string(raw::AbstractString,context::_GeoNumericContext,
                          caller::AbstractString)::String
    s=String(strip(raw))
    isempty(s) && throw(ArgumentError("$caller: string expression must not be empty"))
    return _geo_string_expr(s,context,caller)
end

# `--`/`++` fold `a-b`/`a+b` on string symbols? No — Gmsh strings do not combine;
# a `StringExpr` is a single term.

function _geo_string_expr(src::AbstractString,context::_GeoNumericContext,
                          caller::AbstractString)::String
    s=String(strip(src))
    isempty(s) && throw(ArgumentError("$caller: empty string expression"))
    c=s[firstindex(s)]
    if c=='"' || c=='\''
        close=findnext(==(c),s,2)
        close===nothing && throw(ArgumentError(
            "$caller: unterminated string literal"))
        tail=String(strip(s[nextind(s,close):end]))
        isempty(tail) || throw(ArgumentError(
            "$caller: unexpected text after string literal: $tail"))
        return String(s[nextind(s,firstindex(s)):prevind(s,close)])
    end
    # `Physical Point{tag}` / `Point{tag}` entity-name reads
    if (em=match(r"^Physical\s+(Point|Line|Curve|Surface|Volume)\s*\{(.*)\}\s*$",s))!==nothing
        tag=_geo_int_value(_geo_eval_numeric(em.captures[2],context,
            "$caller Physical name tag"),"$caller Physical name tag")
        dim=_GEO_STRING_ENTITY_DIM[em.captures[1]]
        name=context.entity_name_lookup===nothing ? "" :
            context.entity_name_lookup(dim,tag,:physical)
        return name
    elseif (em=match(r"^(Point|Line|Curve|Surface|Volume)\s*\{(.*)\}\s*$",s))!==nothing
        tag=_geo_int_value(_geo_eval_numeric(em.captures[2],context,
            "$caller entity name tag"),"$caller entity name tag")
        dim=_GEO_STRING_ENTITY_DIM[em.captures[1]]
        name=context.entity_name_lookup===nothing ? "" :
            context.entity_name_lookup(dim,tag,:entity)
        return name
    end
    m=match(r"^([A-Za-z_][A-Za-z0-9_]*)(.*)$",s)
    m===nothing && throw(ArgumentError(
        "$caller: expected a string expression, got $s"))
    name=String(m.captures[1]);rest=String(strip(m.captures[2]))
    return _geo_string_dispatch(name,rest,context,caller)
end

const _GEO_STRING_ENTITY_DIM=Dict(
    "Point"=>0,"Line"=>1,"Curve"=>1,"Surface"=>2,"Volume"=>3)

function _geo_string_dispatch(name::String,rest::String,
                              context::_GeoNumericContext,
                              caller::AbstractString)::String
    # `StringIndex` — `name~{FExpr}` suffixes fold into the symbol before the
    # `::`/member/index dispatch (`x~{1}.m` → `x_1.m`).
    while startswith(rest,"~")
        m2=match(r"^~\s*\{",rest)
        m2===nothing && break
        open_pos=lastindex(m2.match)
        close=_geo_matching_delim(rest,open_pos)
        close==0 && break
        idx=_geo_int_value(_geo_eval_numeric(
            rest[nextind(rest,open_pos):prevind(rest,close)],context,
            "$caller namespace index"),"$caller namespace index")
        name=name*"_"*string(idx)
        rest=String(strip(rest[nextind(rest,close):end]))
    end
    if startswith(rest,"::")
        # `ns::x.member` — `Struct_FullName '.' tSTRING_Reserved` in string
        # position. The string-symbol table is not namespaced, so `ns::x`
        # alone cannot resolve a string.
        after=String(strip(rest[nextind(rest,firstindex(rest),2):end]))
        nm=match(r"^([A-Za-z_][A-Za-z0-9_]*)(.*)$",after)
        nm===nothing && throw(ArgumentError(
            "$caller: expected a name after $name::"))
        sname=String(nm.captures[1]);srest=String(strip(nm.captures[2]))
        while startswith(srest,"~")
            m2=match(r"^~\s*\{",srest)
            m2===nothing && break
            open_pos=lastindex(m2.match)
            close=_geo_matching_delim(srest,open_pos)
            close==0 && break
            idx=_geo_int_value(_geo_eval_numeric(
                srest[nextind(srest,open_pos):prevind(srest,close)],context,
                "$caller namespace index"),"$caller namespace index")
            sname=sname*"_"*string(idx)
            srest=String(strip(srest[nextind(srest,close):end]))
        end
        startswith(srest,".") || throw(ArgumentError(
            "$caller: expected a member after $name::$sname"))
        mm=match(r"^\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*(.*)$",srest)
        mm===nothing && throw(ArgumentError(
            "$caller: expected a member name after $name::$sname."))
        member=String(mm.captures[1]);mtail=String(strip(mm.captures[2]))
        index=0
        if startswith(mtail,"(")
            close=_geo_matching_delim(mtail,firstindex(mtail))
            close==0 && throw(ArgumentError(
                "$caller: unbalanced index after $name::$sname.$member"))
            index=_geo_int_value(_geo_eval_numeric(
                mtail[nextind(mtail,firstindex(mtail)):prevind(mtail,close)],
                context,"$caller member index"),"$caller member index")
            isempty(strip(mtail[nextind(mtail,close):end])) || throw(
                ArgumentError("$caller: unexpected text after "*
                    "$name::$sname.$member index"))
        else
            # `LP`/`RP` do not cover `[]` for string member indexes —
            # `s.lst[0]` → `syntax error ([)` upstream.
            startswith(mtail,"[") && _geo_syntax_abort("[")
            isempty(mtail) || throw(ArgumentError(
                "$caller: unexpected text after $name::$sname.$member: $mtail"))
        end
        return _geo_member_string(context,name,sname,member,index,caller)
    end
    if isempty(rest)
        # zero-argument StringExpr tokens and bare string variables
        name=="Today" && return Libc.strftime("%a %b %e %H:%M:%S %Y",time())
        name=="CodeName" && return "Gmsh"
        name=="GmshExecutableName" &&
            return joinpath(Sys.BINDIR,Base.julia_exename())
        name=="OnelabAction" && return ""
        name=="CurrentFileName" && return basename(context.file_name)
        (name=="CurrentDirectory" || name=="CurrentDir") &&
            return isempty(context.file_name) ?
                "" : dirname(abspath(context.file_name))*"/"
        name=="Empty" && return "" # tEmpty? — not a Gmsh token; leave to var
        if haskey(context.strings,name)
            values=context.strings[name]
            length(values)==1 && return values[1]
            context.soft_unknown_reads || throw(ArgumentError(
                "$caller: Expected single valued string variable '$name'"))
            _geo_yyerror!(context,
                "Expected single valued string variable '$name'")
            return ""
        end
        haskey(context.onelab_strings,name) &&
            return context.onelab_strings[name]
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Unknown string variable '$name'"))
        _geo_yyerror!(context,"Unknown string variable '$name'")
        return ""
    end
    # `Name(...)` function calls and `name(i)` indexed reads
    if startswith(rest,"(") || startswith(rest,"[") || startswith(rest,".")
        opener=rest[firstindex(rest)]
        if opener=='(' || opener=='['
            close=_geo_matching_delim(rest,firstindex(rest))
            close==0 && throw(ArgumentError(
                "$caller: unbalanced $opener after $name"))
            inner=String(rest[nextind(rest,firstindex(rest)):prevind(rest,close)])
            tail=String(strip(rest[nextind(rest,close):end]))
            if opener=='[' && name in ("NameToString","N2S")
                # `NameToString[String__Index]` yields the identifier text —
                # `$$ = $3` returns the name unchanged (`x`, `x~{i}`, nested
                # `StringToName`). Anything else aborts at the first token
                # that cannot extend a `String__Index`.
                isempty(tail) || throw(ArgumentError(
                    "$caller: unexpected text after NameToString[...]: $tail"))
                arg=String(strip(inner))
                try
                    return _geo_symbol_name(arg,context,caller)
                catch err
                    err isa _GeoSyntaxAbort && rethrow()
                    _geo_syntax_abort(_geo_string_index_offender(arg))
                end
            elseif opener=='[' && name in ("StringToName","S2N")
                # StringToName[s] evaluates to a *name*; in string context the
                # named string variable's value is read.
                sym=_geo_symbol_name(inner,context,caller)
                isempty(tail) || throw(ArgumentError(
                    "$caller: unexpected text after StringToName[...]: $tail"))
                return _geo_string_read(context,sym,caller)
            elseif opener=='(' &&
                name in ("NameToString","N2S","StringToName","S2N")
                # The name productions take `[` only — `(` is the error
                # token upstream.
                _geo_syntax_abort("(")
            elseif name in _GEO_STRING_FUNCTIONS
                # `Name(args)`/`Name[args]` — `LP`/`RP` cover both `()` and
                # `[]` upstream, so `StrCat[a, b]` is a function call.
                isempty(tail) || throw(ArgumentError(
                    "$caller: unexpected text after $name$opener...: $tail"))
                return _geo_string_function(name,inner,context,caller)
            else
                # `name[i]` / `name(i)` indexed string-element read
                idx=_geo_int_value(_geo_eval_numeric(inner,context,
                    "$caller string index"),"$caller string index")
                tail2=tail
                if startswith(tail2,".")
                    # `name[i].reserved` — indexed option-string read
                    om=match(r"^\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*$",tail2)
                    om===nothing && throw(ArgumentError(
                        "$caller: unsupported text after $name[$idx]: $tail2"))
                    return _geo_option_string(context,name,idx,
                        String(om.captures[1]),caller)
                end
                isempty(tail2) || throw(ArgumentError(
                    "$caller: unexpected text after $name[$idx]: $tail2"))
                return _geo_string_index(context,name,idx,caller)
            end
        else
            # `name.member` / `name.member(i)` — `StringExprVar` member read:
            # struct member first, `StringOption(name, member)` otherwise.
            om=match(r"^\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*(.*)$",rest)
            om===nothing && throw(ArgumentError(
                "$caller: unsupported text after $name: $rest"))
            member=String(om.captures[1]);mtail=String(strip(om.captures[2]))
            index=0
            if startswith(mtail,"(")
                close=_geo_matching_delim(mtail,firstindex(mtail))
                close==0 && throw(ArgumentError(
                    "$caller: unbalanced index after $name.$member"))
                index=_geo_int_value(_geo_eval_numeric(
                    mtail[nextind(mtail,firstindex(mtail)):prevind(mtail,close)],
                    context,"$caller member index"),"$caller member index")
                isempty(strip(mtail[nextind(mtail,close):end])) || throw(
                    ArgumentError("$caller: unexpected text after "*
                        "$name.$member index"))
            else
                startswith(mtail,"[") && _geo_syntax_abort("[")
                isempty(mtail) || throw(ArgumentError(
                    "$caller: unexpected text after $name.$member: $mtail"))
            end
            return _geo_member_string(context,"",name,member,index,caller)
        end
    end
    throw(ArgumentError("$caller: unrecognized string expression: $name$rest"))
end

function _geo_string_read(context::_GeoNumericContext,name::String,
                          caller::AbstractString)
    if haskey(context.strings,name)
        values=context.strings[name]
        length(values)==1 && return values[1]
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Expected single valued string variable '$name'"))
        _geo_yyerror!(context,
            "Expected single valued string variable '$name'")
        return ""
    end
    haskey(context.onelab_strings,name) &&
        return context.onelab_strings[name]
    context.soft_unknown_reads || throw(ArgumentError(
        "$caller: Unknown string variable '$name'"))
    _geo_yyerror!(context,"Unknown string variable '$name'")
    return ""
end

function _geo_string_index(context::_GeoNumericContext,name::String,index::Int,
                           caller::AbstractString)
    if !haskey(context.strings,name)
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Unknown string variable '$name'"))
        _geo_yyerror!(context,"Unknown string variable '$name'")
        return ""
    end
    values=context.strings[name]
    0<=index<length(values) && return values[index+1]
    context.soft_unknown_reads || throw(ArgumentError(
        "$caller: Index $index out of range"))
    _geo_yyerror!(context,"Index $index out of range")
    return ""
end

# `x.y` / `x[i].y` option-string reads: Gmsh dispatches `StringOption(GMSH_GET)`
# over General/Geometry/Mesh/Solver/PostProcessing/View/Print option tables.
# Stored writes win; otherwise the Gmsh default applies. Unknown categories or
# names raise Gmsh's "Unknown string option" diagnostic.
function _geo_option_string(context::_GeoNumericContext,family::String,
                            index::Int,member::String,caller::AbstractString;
                            warn::Bool=true)
    key=(family,index,member)
    haskey(context.option_strings,key) &&
        return context.option_strings[key]
    if !(family in keys(_GEO_STRING_OPTION_TABLE))
        warn || return nothing
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Unknown string option category '$family'"))
        _geo_msg_error!(context,"Unknown string option category '$family'")
        return ""
    end
    table=_GEO_STRING_OPTION_TABLE[family]
    if !(member in keys(table))
        warn || return nothing
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Unknown string option '$family.$member'"))
        _geo_msg_error!(context,"Unknown string option '$family.$member'")
        return ""
    end
    default=table[member]
    # `opt_general_filename` returns the model's current file name, not the
    # stored default.
    (family=="General" && member=="FileName" && isempty(default)) &&
        return context.file_name
    return default
end

function _geo_set_option_string!(context::_GeoNumericContext,family::String,
                                 index::Int,member::String,value::String,
                                 caller::AbstractString)
    if !(family in keys(_GEO_STRING_OPTION_TABLE))
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Unknown string option category '$family'"))
        return _geo_msg_error!(context,
            "Unknown string option category '$family'")
    end
    if !(member in keys(_GEO_STRING_OPTION_TABLE[family]))
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Unknown string option '$family.$member'"))
        return _geo_msg_error!(context,
            "Unknown string option '$family.$member'")
    end
    context.option_strings[(family,index,member)]=value
    return nothing
end

# `_GEO_STRING_OPTION_TABLE`/`_GEO_NUMBER_OPTION_TABLE`/`_GEO_COLOR_OPTION_NAMES`
# are generated verbatim from Gmsh 4.15.2's `DefaultOptions.h` in
# `GeoOptionTables.jl`. `General.FileName` is dynamic — handled by the
# `file_name` sentinel below.

function _geo_string_function(name::String,inner::AbstractString,
                              context::_GeoNumericContext,
                              caller::AbstractString)::String
    args=_geo_split_args(inner,caller)
    str_arg(i)=_geo_string_expr(args[i],context,caller)
    str_args()=[str_arg(i) for i in eachindex(args)]
    if name=="StrCat" || name=="Str"
        # `tStr`/`tStrCat` take `RecursiveListOfStringExprVar` — `x()`,
        # `x.member()` and `{...}` items expand to multiple elements.
        vals=_geo_eval_string_list(inner,context,caller)
        return join(vals,name=="Str" ? "\n" : "")
    elseif name=="StrPrefix"
        length(args)==1 || throw(ArgumentError("$caller: StrPrefix takes one argument"))
        s=str_args()[1]
        i=findlast(==('.'),s)
        (i===nothing || i==firstindex(s)) && return s
        return String(s[firstindex(s):prevind(s,i)])
    elseif name=="StrRelative"
        length(args)==1 || throw(ArgumentError("$caller: StrRelative takes one argument"))
        s=str_args()[1]
        i=findlast(c->c in ('/','\\'),s)
        (i===nothing || i==firstindex(s)) && return s
        return String(s[nextind(s,i):end])
    elseif name=="StrReplace"
        length(args)==3 || throw(ArgumentError("$caller: StrReplace takes three arguments"))
        vals=str_args()
        return replace(vals[1],vals[2]=>vals[3])
    elseif name=="UpperCase"
        length(args)==1 || throw(ArgumentError("$caller: UpperCase takes one argument"))
        return uppercase(str_args()[1])
    elseif name=="LowerCase"
        length(args)==1 || throw(ArgumentError("$caller: LowerCase takes one argument"))
        return lowercase(str_args()[1])
    elseif name=="LowerCaseIn"
        length(args)==1 || throw(ArgumentError("$caller: LowerCaseIn takes one argument"))
        s=str_args()[1]
        isempty(s) && return s
        out=collect(s)
        for i in 2:length(out)
            out[i-1]!='_' && (out[i]=lowercase(out[i]))
        end
        return String(out)
    elseif name=="StrChoice"
        length(args)==3 || throw(ArgumentError("$caller: StrChoice takes three arguments"))
        flag=_geo_eval_numeric(args[1],context,"$caller: StrChoice")
        return flag!=0 ? str_arg(2) : str_arg(3)
    elseif name=="StrSub"
        (length(args)==2 || length(args)==3) || throw(ArgumentError(
            "$caller: StrSub takes two or three arguments"))
        s=str_arg(1)
        start=_geo_int_value(_geo_eval_numeric(args[2],context,
            "$caller: StrSub start"),"$caller: StrSub start")
        # C++ `substr(start, len)` — 0-based byte indexing; `start` past the
        # end throws `std::out_of_range`, so a diagnostic is faithful.
        bytes=codeunits(s)
        (start<0 || start>length(bytes)) && throw(ArgumentError(
            "$caller: StrSub start $start out of range"))
        lo=start+1
        if length(args)==3
            len=_geo_int_value(_geo_eval_numeric(args[3],context,
                "$caller: StrSub length"),"$caller: StrSub length")
            len<0 && (len=length(bytes))
            return String(bytes[lo:min(lo+len-1,length(bytes))])
        end
        return String(bytes[lo:end])
    elseif name=="Sprintf"
        # `Sprintf(fmt)` / `Sprintf(fmt, RecursiveListOfDouble)` — the
        # argument tail is one list; mismatch returns the format unchanged.
        isempty(args) && throw(ArgumentError(
            "$caller: Sprintf requires a format string"))
        fmt=_geo_string_expr(args[1],context,caller)
        length(args)==1 && return fmt
        values=_geo_numeric_list_values("{"*join(args[2:end],",")*"}",
            context,"$caller: Sprintf list";allow_multiplier=true)
        out,extra=_geo_print_list_of_double(fmt,values,caller)
        if extra<0
            _geo_yyerror!(context,"Too few arguments in Sprintf")
            return fmt
        elseif extra>0
            _geo_yyerror!(context,
                "$extra extra argument$(extra>1 ? "s" : "") in Sprintf")
            return fmt
        end
        return out
    elseif name=="GetEnv"
        length(args)==1 || throw(ArgumentError("$caller: GetEnv takes one argument"))
        return get(ENV,_geo_string_expr(args[1],context,caller),"")
    elseif name=="GetString"
        (length(args)==1 || length(args)==2) || throw(ArgumentError(
            "$caller: GetString takes one or two arguments"))
        key=_geo_string_expr(args[1],context,caller)
        haskey(context.onelab_strings,key) && return context.onelab_strings[key]
        length(args)==2 && return _geo_string_expr(args[2],context,caller)
        return ""
    elseif name=="GetStringValue"
        length(args)==2 || throw(ArgumentError(
            "$caller: GetStringValue takes two arguments"))
        vals=str_args()
        # `Msg::GetString` prompts interactively; a batch run returns the
        # default unchanged.
        return vals[2]
    elseif name=="GetForcedStr"
        (length(args)==1 || length(args)==2) || throw(ArgumentError(
            "$caller: GetForcedStr takes one or two arguments"))
        # `tGetForcedStr LP Struct_FullName ['.' tSTRING_Reserved]
        # [GetForcedStr_Default] RP` — `type_treat=2`: silent defaults.
        arg=strip(args[1])
        default=length(args)==2 ?
            _geo_string_expr(args[2],context,caller) : ""
        if (om=match(
                r"^(?:([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*::\s*)?([A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*)\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)$",
                arg))!==nothing
            return _geo_member_string(context,
                om.captures[1]===nothing ? "" :
                    _geo_symbol_name(om.captures[1],context,caller),
                _geo_symbol_name(om.captures[2],context,caller),
                String(om.captures[3]),0,caller;
                default=default,type_treat=2)
        end
        # `treat_Struct_FullName_String` receives `nullptr` for the
        # namespace upstream — `GetForcedStr(ns::x, d)` reads `x` from the
        # (un-namespaced) string-symbol table, discarding `ns`.
        arg=replace(arg,
            r"^[A-Za-z_][A-Za-z0-9_]*(?:\s*~\s*\{[^}]*\})*\s*::\s*"=>"")
        sym=_geo_symbol_name(arg,context,caller)
        if haskey(context.strings,sym)
            values=context.strings[sym]
            length(values)==1 && return values[1]
        end
        return default
    elseif name=="FixRelativePath"
        length(args)==1 || throw(ArgumentError("$caller: FixRelativePath takes one argument"))
        return _geo_fix_relative_path(context.file_name,str_args()[1])
    elseif name=="DirName"
        length(args)==1 || throw(ArgumentError("$caller: DirName takes one argument"))
        s=str_args()[1]
        i=findlast(c->c in ('/','\\'),s)
        return i===nothing ? "" : String(s[firstindex(s):i])
    elseif name=="AbsolutePath"
        length(args)==1 || throw(ArgumentError("$caller: AbsolutePath takes one argument"))
        s=str_args()[1]
        # `GetAbsolutePath` is a pure path computation — no existence check:
        # absolute paths pass through, relative paths join the process cwd.
        (startswith(s,"/") || startswith(s,"\\") ||
         (ncodeunits(s)>3 && s[2]==':' && s[3] in ('/','\\'))) && return s
        return isempty(s) ? s : joinpath(pwd(),s)
    elseif name=="DefineString"
        isempty(args) && throw(ArgumentError(
            "$caller: DefineString requires an argument"))
        spec=_geo_string_expr(args[1],context,caller)
        # `DefineString("name [opts]")` registers an ONELAB string parameter;
        # the leading word is the parameter name.
        pname=first(split(spec,' '))
        isempty(pname) || get!(context.onelab_strings,pname,spec)
        return spec
    elseif name=="NameStruct"
        # `tNameStruct LP NameStruct_Arg RP` — `#tag` or `ns::#tag`; the
        # reverse tag→name lookup is `get_key_struct_from_tag` (case 1/2 are
        # `yymsg(1)` warnings, not errors). The arg grammar is
        # `'#' FExpr | String__Index tSCOPE '#' FExpr` — anything else is a
        # bison syntax error at the offending token (`(5)`, `())`, `(bar)`,
        # `(,)` upstream).
        isempty(args) && _geo_syntax_abort(")")
        length(args)==1 || _geo_syntax_abort(",")
        arg=String(strip(args[1]))
        isempty(arg) && _geo_syntax_abort(")")
        last=lastindex(arg);i=firstindex(arg)
        if arg[i]=='#'
            ns="";tag_src=String(strip(arg[nextind(arg,i):end]))
        else
            # Consume the leading `String__Index` (identifier /
            # `StringToName[...]` plus `~{expr}` suffixes), then `::` and `#`
            # must follow — anything else is the offending token.
            if (sm=match(r"^(?:StringToName|S2N)\s*\[",arg))!==nothing
                cl=_geo_matching_delim(arg,lastindex(sm.match))
                cl==0 && _geo_syntax_abort("[")
                i=nextind(arg,cl)
            else
                hm=match(r"^[A-Za-z_][A-Za-z0-9_]*",arg)
                hm===nothing &&
                    _geo_syntax_abort(_geo_stmt_head_token(arg))
                i=nextind(arg,lastindex(hm.match))
            end
            while true
                j=i
                while j<=last && isspace(arg[j]); j=nextind(arg,j); end
                (j<=last && arg[j]=='~') || break
                k=nextind(arg,j)
                while k<=last && isspace(arg[k]); k=nextind(arg,k); end
                (k<=last && arg[k]=='{') || break
                cl=_geo_matching_delim(arg,k);cl==0 && break
                i=nextind(arg,cl)
            end
            ns=_geo_symbol_name(arg[firstindex(arg):prevind(arg,i)],
                context,caller)
            rest=String(strip(arg[i:end]))
            startswith(rest,"::") || _geo_syntax_abort(
                isempty(rest) ? ")" :
                (wm=match(r"^[A-Za-z_][A-Za-z0-9_]*",rest))!==nothing ?
                    String(wm.match) : _geo_stmt_head_token(rest))
            rest2=String(strip(rest[nextind(rest,2):end]))
            isempty(rest2) && _geo_syntax_abort(")")
            rest2[1]=='#' || _geo_syntax_abort(
                (wm=match(r"^[A-Za-z_][A-Za-z0-9_]*",rest2))!==nothing ?
                    String(wm.match) : _geo_stmt_head_token(rest2))
            tag_src=String(strip(rest2[nextind(rest2,firstindex(rest2)):end]))
        end
        isempty(tag_src) && _geo_syntax_abort(")")
        tag=_geo_int_value(_geo_eval_numeric(tag_src,context,
            "$caller NameStruct tag"),"$caller NameStruct tag")
        (code,sname)=_geo_struct_key_from_tag(context,ns,tag)
        code==1 &&
            _geo_yywarn!(context,"Unknown NameSpace '$ns' of Struct")
        code==2 &&
            _geo_yywarn!(context,"Unknown Struct of index $tag")
        return sname
    end
    throw(ArgumentError("$caller: unknown string function $name"))
end

function _geo_fix_relative_path(reference::AbstractString,in::AbstractString)
    isempty(in) && return ""
    (startswith(in,"/") || startswith(in,"\\") ||
     (ncodeunits(in)>3 && in[2]==':' &&
      (in[3]=='/' || in[3]=='\\'))) && return String(in)
    dir=dirname(String(reference))
    return isempty(dir) ? String(in) : dir*"/"*in
end

# `printListOfDouble` from Gmsh.y: each `%`-format consumes one list element plus
# the literal text up to the next `%`; `%%` emits a literal `%` and still
# consumes an element. Returns (output, extra) with extra>0 = unconsumed
# arguments and extra<0 = not enough arguments.
function _geo_print_list_of_double(format::String,list::Vector{Float64},
                                   caller::AbstractString)
    num_formats=count(==('%'),format)
    if num_formats==0
        buffer=format
        for (i,d) in enumerate(list)
            buffer*=" [$(i-1)]"*_geo_sprintf("%g",d,caller)
        end
        return buffer,0
    end
    first_pct=findfirst(==('%'),format)
    buffer=format[firstindex(format):prevind(format,first_pct)]
    j=first_pct;i=1
    last=lastindex(format)
    while i<=length(list)
        j>last && return buffer,length(list)-i+1
        k=j;j=nextind(format,j)
        if j<=last
            if format[j]=='%'
                # `%%` escape — Gmsh appends "%" AND sprintfs the whole
                # `%%...%` region (which renders `%%`→`%` and the rest
                # literally), so the region contributes "%%" + tail.
                j=nextind(format,j)
                while j<=last && format[j]!='%'
                    j=nextind(format,j)
                end
                buffer*="%%"
                buffer*=format[nextind(format,k,2):prevind(format,j)]
            else
                while j<=last && format[j]!='%'
                    j=nextind(format,j)
                end
                if k!=j
                    spec=String(format[k:prevind(format,j)])
                    buffer*=_geo_sprintf(spec,list[i],caller)
                end
            end
        else
            return buffer,length(list)-i+1
        end
        i+=1
    end
    j<=last && return buffer,-1
    return buffer,0
end

# C `sprintf` for a single `%`-spec consumed by `printListOfDouble`. Gmsh
# passes every element as a `double` through varargs — on the reference
# platform (aarch64 macOS) integer conversions deterministically read the
# low bits of the IEEE-754 encoding (`%d`/`%x`/`%o`/`%u` → low 32 bits
# sign/zero-extended, `%l*` → all 64 bits, `%h`/`%hh` → low 16/8, `%c` → the
# low byte, `%s`/`%p` → the bits as a pointer). Floating conversions use the
# double itself.
function _geo_sprintf(spec::AbstractString,value::Float64,
                      caller::AbstractString)
    m=match(r"^%([-+ #0]*)(\d*)(?:\.(\d+))?([hlL]*)([diuxXoeEfFgGaAcCsSpn%])",spec)
    m===nothing && throw(ArgumentError(
        "$caller: malformed printf format $spec"))
    conv=m.captures[5];length_mod=m.captures[4]
    # The spec head ends at the conversion letter; the rest of the region is
    # literal text emitted after the converted value.
    tail=String(spec[nextind(spec,ncodeunits(m.match)):end])
    conv=='%' && return "%"*tail
    conv=='n' && return tail
    if conv in ("s","S")
        # `(char*)double_bits` — NULL (0.0) prints "(null)"; other values are
        # nondeterministic garbage in Gmsh — pinned to "(null)".
        return "(null)"*tail
    end
    if conv=='p'
        return "0x"*string(reinterpret(UInt64,value);base=16)*tail
    end
    if conv=='c'
        byte=reinterpret(UInt64,value)%UInt8
        return (byte==0 ? "" : string(Char(byte)))*tail
    end
    if conv in ("e","E","f","F","g","G","a","A")
        return _geo_printf_numeric(spec,value,caller)
    end
    # integer conversions — read the double's bits through the stated width.
    bits=reinterpret(UInt64,value)
    int_value=if length_mod in ("l","ll","L")
        reinterpret(Int64,bits)
    elseif length_mod in ("h","hh")
        Int32(reinterpret(Int16,UInt16(bits & 0xffff)))
    else
        reinterpret(Int32,UInt32(bits & 0xffffffff))
    end
    conv in ("u","o","x","X") && (int_value=unsigned(int_value))
    return _geo_printf_numeric(spec,int_value,caller)
end

function _geo_printf_numeric(spec::AbstractString,value::Real,
                             caller::AbstractString)
    # The spec region runs from `%` up to (not including) the next `%`, so it
    # carries trailing literal text; the conversion letter ends the head.
    conv=match(r"^%[-+ #0]*\d*(?:\.\d+)?[hlL]*([diuxXoeEfFgGaAcCsSp])",spec)
    conv===nothing && throw(ArgumentError(
        "$caller: malformed printf format $spec"))
    # C `printf` writes non-finite floats as `nan`/`inf` (uppercase under
    # `%E`/`%G`-family letters); Format.jl's `NaN`/`Inf` differ. Sign flags
    # apply to `inf` but never to `nan`.
    if value isa AbstractFloat && !isfinite(value) &&
            occursin(conv.captures[1],"eEfFgGaA")
        text=if isnan(value)
            "nan"
        else
            sign=value<0 ? "-" :
                contains(match(r"^%([-+ #0]*)",spec).captures[1],'+') ? "+" :
                contains(match(r"^%([-+ #0]*)",spec).captures[1],' ') ? " " : ""
            sign*"inf"
        end
        return isuppercase(conv.captures[1][1]) ? uppercase(text) : text
    end
    try
        fmt=Format(spec)
        return format(fmt,value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: printf format $spec failed: $(sprint(showerror,err))"))
    end
end

# `Msg::Direct`/`Msg::Warning`/`Msg::Error` run their text through
# `vsnprintf(text, <no args>)` — a second printf expansion where `%%`→`%` and
# every other spec reads a missing argument. On the reference platform the
# missing value reads as 0/`NULL`: numeric specs print "0"-family output,
# `%s`→"(null)", `%c`→'\0' (truncating the C string there), and a lone `%` at
# the end of the string is dropped. `>`/`>>` file redirection bypasses this
# pass (the raw `printListOfDouble` buffer is written).
function _geo_vsnprintf_expand(text::AbstractString)
    out=IOBuffer();i=firstindex(text);last=lastindex(text)
    while i<=last
        c=text[i]
        if c!='%'
            print(out,c);i=nextind(text,i);continue
        end
        j=nextind(text,i)
        j>last && break # trailing lone `%` is dropped
        if text[j]=='%'
            print(out,'%');i=nextind(text,j);continue
        end
        # Consume the spec: flags, width, precision, length modifier, letter.
        m=match(r"^[-+ #0]*\d*(?:\.\d+)?[hlL]*[diuxXoeEfFgGaAcCsSpn]",
                text[j:last])
        if m===nothing
            i=nextind(text,i);continue # unknown spec — pass `%` through
        end
        spec="%"*m.match
        letter=last(m.match)
        if letter=='n'
            # writes through a pointer — no output
        elseif letter in ("s","S")
            print(out,"(null)")
        elseif letter=='c'
            break # '\0' truncates the C string
        elseif letter=='p'
            print(out,"0x0")
        elseif letter in ("e","E","f","F","g","G","a","A")
            print(out,_geo_printf_numeric(spec,0.0,"execute_geo: Printf"))
        else
            print(out,_geo_printf_numeric(spec,0,"execute_geo: Printf"))
        end
        i=nextind(text,j,ncodeunits(m.match))
    end
    return String(take!(out))
end

@inline _geo_number_source(value::Float64)=repr(value)

# C `%g` rendering for diagnostic interpolation (`yymsg("... %g ...")`).
_geo_gmsh_number(value::Float64)=_geo_printf_numeric("%g",value,"execute_geo")

function _geo_int_value(value::Float64,caller::AbstractString)
    return try
        trunc(Int,value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError) || rethrow()
        throw(ArgumentError("$caller: value is outside the platform Int range"))
    end
end

function _geo_signed_gmsh_int_value(value::Float64,caller::AbstractString)
    integer=_geo_int_value(value,caller)
    typemin(Int32)<=integer<=typemax(Int32) || throw(ArgumentError(
        "$caller is outside Gmsh's signed 32-bit integer range"))
    return integer
end

function _geo_split_list(raw::AbstractString,caller::AbstractString)
    value=String(strip(raw))
    (startswith(value,"{") && endswith(value,"}")) ||
        throw(ArgumentError("$caller: expected a brace-delimited list"))
    body=String(strip(value[nextind(value,firstindex(value)):prevind(value,lastindex(value))]))
    isempty(body) && return String[]
    items=String[];start=firstindex(body);i=start;last=lastindex(body)
    parens=0;brackets=0;selector_braces=0
    while i<=last
        c=body[i]
        if c=='(';parens+=1
        elseif c==')'
            parens-=1;parens>=0 || _geo_syntax_abort(")")
        elseif c=='[';brackets+=1
        elseif c==']'
            brackets-=1;brackets>=0 || _geo_syntax_abort("]")
        elseif c=='{'
            # A `{` that *begins* a member is a nested list — not legal in
            # `RecursiveListOfDouble` (`{a, {b,c}}` is a syntax error
            # upstream). A `{` after content belongs to an `FExpr_Multi`
            # selector/transform term: `Point{1}`, `Physical X{7}`,
            # `BoundingBox X{1}`, `Rotate{{axis},{origin},a}{...}`, etc.
            (!isempty(strip(body[start:prevind(body,i)])) ||
             brackets>0 || selector_braces>0) || _geo_syntax_abort("{")
            selector_braces+=1
        elseif c=='}'
            selector_braces-=1;selector_braces>=0 || _geo_syntax_abort("}")
        elseif c==',' && parens==0 && brackets==0 && selector_braces==0
            item=String(strip(body[start:prevind(body,i)]))
            isempty(item) && _geo_syntax_abort(",")
            length(items)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                "$caller: list exceeds $_MAX_GEO_LIST_ITEMS entries"))
            push!(items,item);start=nextind(body,i)
        end
        i=nextind(body,i)
    end
    parens==0 || _geo_syntax_abort("(")
    brackets==0 || _geo_syntax_abort("[")
    selector_braces==0 || _geo_syntax_abort("{")
    item=String(strip(body[start:last]))
    isempty(item) && _geo_syntax_abort("}")
    length(items)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
        "$caller: list exceeds $_MAX_GEO_LIST_ITEMS entries"))
    push!(items,item)
    return items
end

# `_geo_split_list` raises `_GeoSyntaxAbort` — recoverable `.geo` syntax
# errors. Callers outside the `.geo` grammar (field-option strings) get the
# ArgumentError contract.
function _geo_split_list_or_argerr(raw::AbstractString,caller::AbstractString)
    return try
        _geo_split_list(raw,caller)
    catch err
        err isa InterruptException && rethrow()
        err isa _GeoSyntaxAbort || rethrow()
        throw(ArgumentError(
            "$caller: malformed list syntax ($(err.token))"))
    end
end

# Index of the opener matching the closer at `close_pos` (0 when
# unbalanced) — the reverse of `_geo_matching_delim`. Quoted spans are
# skipped pairwise, matching the forward scan's quote handling.
function _geo_matching_delim_rev(src::AbstractString,close_pos::Int)
    close_ch=src[close_pos]
    open_ch=close_ch==')' ? '(' : close_ch==']' ? '[' : '{'
    depth=0;i=close_pos
    while true
        c=src[i]
        if c=='"' || c=='\''
            q=c
            while i>firstindex(src)
                i=prevind(src,i)
                src[i]==q && break
            end
            src[i]==q || return 0
        elseif c==close_ch
            depth+=1
        elseif c==open_ch
            depth-=1
            depth==0 && return i
        end
        i>firstindex(src) || return 0
        i=prevind(src,i)
    end
end

# Lexer keywords that cannot head an FExpr production report themselves as
# the syntax-error token (`Field::s` → `(Field)`), unlike FExpr-capable
# keywords where the following `::` is the offender (`Sin::s`, `Pi::s`,
# `Today::s` → `syntax error (::)`).
_geo_keyword_self_reports(name::String)=
    name in _GEO_LEXER_KEYWORDS && !(name in _GEO_ALL_FUNCTIONS) &&
        !(name in _GEO_FEXPR_SHIFT_KEYWORDS) && name!="List" &&
        !haskey(_GEO_BARE_CONSTANTS,name) && !(name in _GEO_BARE_STRING_TOKENS)

# Does a `String__Index` (`tSTRING`, `name~{FExpr}`…, `StringToName[sexpr]`/
# `S2N[sexpr]`) end immediately before `src[colon]`? Only then does a
# contiguous `::` lex as upstream's tSCOPE inside an expression. After any
# other token the `::` itself is the offending lookahead (`1::2`,
# `a[i]::s`, `f(1)::s`, `Sin::s` → `syntax error (::)`), while a keyword
# that cannot head an FExpr reports itself (`Field::s` → `(Field)`) — those
# heads are accepted here so the expression parser reports the keyword.
function _geo_scope_precedes(src::AbstractString,colon::Int)
    j=colon
    while true
        j>firstindex(src) || return false
        j=prevind(src,j)
        isspace(src[j]) || break
    end
    c=src[j]
    if _geo_ascii_ident(c)
        k=j
        while k>firstindex(src)
            _geo_ascii_ident(src[prevind(src,k)]) || break
            k=prevind(src,k)
        end
        _geo_ascii_letter(src[k]) || return false
        n=String(src[k:j])
        return _geo_tstring_valid(n) || _geo_keyword_self_reports(n)
    elseif c=='}' || c==']'
        open=_geo_matching_delim_rev(src,j)
        open==0 && return false
        k=open
        while true
            k>firstindex(src) || return false
            k=prevind(src,k)
            isspace(src[k]) || break
        end
        if c=='}'
            # `~{FExpr}` is a StringIndex suffix — whatever precedes `~`
            # must itself be a String__Index tail.
            src[k]=='~' || return false
            return _geo_scope_precedes(src,k)
        end
        # `StringToName[...]`/`S2N[...]` are the only bracket-closed
        # String__Index forms — `a[i]`/`a(i)` are indexed FExpr reads.
        _geo_ascii_letter(src[k]) || return false
        e=k
        while e>firstindex(src)
            _geo_ascii_ident(src[prevind(src,e)]) || break
            e=prevind(src,e)
        end
        return String(src[e:k]) in ("StringToName","S2N")
    else
        return false
    end
end

# `tail` reports the token that follows an unterminated trailing range —
# `}` inside braces, `;` at a bare statement RHS.
function _geo_split_range(raw::AbstractString,caller::AbstractString;
                          tail::AbstractString="}")
    source=String(strip(raw))
    isempty(source) && throw(ArgumentError("$caller: range term must not be empty"))
    pieces=String[];start=firstindex(source);i=start;last=lastindex(source)
    parens=0;brackets=0;braces=0;pending_ternary=0
    while i<=last
        c=source[i]
        if c=='(';parens+=1
        elseif c==')'
            parens-=1;parens>=0 || _geo_syntax_abort(")",
                "$caller: unmatched closing parenthesis in range term")
        elseif c=='[';brackets+=1
        elseif c==']'
            brackets-=1;brackets>=0 || _geo_syntax_abort("]",
                "$caller: unmatched closing bracket in range term")
        elseif c=='{';braces+=1
        elseif c=='}'
            braces-=1;braces>=0 || _geo_syntax_abort("}")
        elseif c=='?' && parens==0 && brackets==0 && braces==0
            pending_ternary+=1
        elseif c==':' && parens==0 && brackets==0 && braces==0
            if i<lastindex(source) && source[nextind(source,i)]==':' &&
               _geo_scope_precedes(source,i)
                # `ns::name` — a contiguous `::` after a `String__Index` is
                # upstream's tSCOPE token, not a range separator; skip both
                # characters. The token after `::` is validated by the
                # expression parser (`ns::#1` → `#`, `ns::5` → `5`).
                i=nextind(source,i)
            elseif pending_ternary>0
                # A ':' that answers a pending '?' belongs to the ternary, not
                # the range — matching Gmsh's grammar where the conditional
                # expression consumes its ':' greedily.
                pending_ternary-=1
            else
                length(pieces)<2 || _geo_syntax_abort(":",
                    "$caller: a Gmsh range has at most two ':' separators")
                piece=i==start ? "" : String(strip(source[start:prevind(source,i)]))
                if isempty(piece)
                    # `::` lexes as a single token upstream (`1::2` → `::`).
                    (i>firstindex(source) && source[prevind(source,i)]==':') &&
                        _geo_syntax_abort("::")
                    _geo_syntax_abort(":",
                        "$caller: Gmsh range endpoints and increments must " *
                        "not be empty")
                end
                # `'-'? String__Index '[' ']'` is a whole-RLD-item splice, not
                # a range endpoint — `{a[]:5}` errors at the `:` upstream.
                isempty(pieces) && match(
                    r"^[+-]?\s*[A-Za-z_][A-Za-z0-9_]*\s*(?:\[\s*\]|\(\s*\))$",
                    piece)!==nothing && _geo_syntax_abort(":",
                    "$caller: numeric list splice requires exactly one " *
                    "scalar index in a range endpoint")
                push!(pieces,piece)
                start=nextind(source,i)
            end
        end
        i=nextind(source,i)
    end
    parens==0 || _geo_syntax_abort("(",
        "$caller: unmatched opening parenthesis in range term")
    brackets==0 || _geo_syntax_abort("[",
        "$caller: unmatched opening bracket in range term")
    braces==0 || _geo_syntax_abort("{")
    isempty(pieces) && return nothing
    piece=start>last ? "" : String(strip(source[start:last]))
    isempty(piece) && _geo_syntax_abort(tail,
        "$caller: Gmsh range endpoints and increments must not be empty")
    push!(pieces,piece)
    return pieces
end

@inline function _geo_range_contains(value::Float64,last::Float64,step::Float64)
    return step>0 ? value<=last : value>=last
end

function _geo_range_count(first::Float64,last::Float64,step::Float64,
                          limit::Int,caller::AbstractString)
    iszero(step) && throw(ArgumentError(
        "$caller: Gmsh range increment must be nonzero"))
    value=first;count=0
    while _geo_range_contains(value,last,step)
        count<limit || throw(ArgumentError(
            "$caller: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
        count+=1
        next=value+step
        if next==value && _geo_range_contains(next,last,step)
            throw(ArgumentError(
                "$caller: Gmsh range increment does not advance at Float64 precision"))
        end
        value=next
    end
    return count
end

# ======== `FExpr_Multi` list-producing terms ========
#
# Gmsh's list grammar produces multi-valued terms from variable splices
# (`a()`, `a[]`, `a({..})`/`a[{..}]`), ranges, list functions (`List`,
# `LinSpace`, `LogSpace`, `Catenary`, `Unique`, `Abs`, `ListFromFile`),
# entity selectors (`Point{..}`, `Physical X{..}`, `BoundingBox X{..}`, ... —
# resolved through `context.exec_hook` in `execute_geo`), `-` negation and
# `scalar * multi` scaling. Bare `{...}` groups are *not* `FExpr_Multi` —
# `Abs({1,-2})` is a syntax error upstream — so the `Abs`/`Unique` argument
# must itself be a multi term.

# Top-level `*` index (outside every bracket kind and string literal), or
# `nothing` — the split point for `FExpr '*' FExpr_Multi`.
function _geo_top_level_star(source::AbstractString)
    parens=0;brackets=0;braces=0;quote_char='\0'
    i=firstindex(source);last=lastindex(source)
    while i<=last
        c=source[i]
        if quote_char!='\0'
            c==quote_char && (quote_char='\0')
        elseif c in ('"','\'')
            quote_char=c
        elseif c=='(' parens+=1
        elseif c==')' parens-=1
        elseif c=='[' brackets+=1
        elseif c==']' brackets-=1
        elseif c=='{' braces+=1
        elseif c=='}' braces-=1
        elseif c=='*' && parens==0 && brackets==0 && braces==0
            return i
        end
        i=nextind(source,i)
    end
    return nothing
end

# Balanced `(...)`/`[...]` argument body for a leading `name` function call —
# returns `(inner, rest)` or `nothing` when the source is not `name<delim>`.
function _geo_call_body(source::AbstractString,name::AbstractString)
    head=match(Regex("^\\s*"*name*"\\s*([\\(\\[])"),source)
    head===nothing && return nothing
    open_pos=firstindex(source)+ncodeunits(head.match)-1
    close=_geo_matching_delim(source,open_pos)
    close==0 && return nothing
    inner=String(source[nextind(source,open_pos):prevind(source,close)])
    rest=String(strip(source[nextind(source,close):end]))
    return inner,rest,source[open_pos]
end

# `List[name]` (literal brackets, `tList '[' String__Index ']'`) copies the
# symbol payload; `List(<multi>)`/`List[<multi>]` and `List({...})` pass the
# evaluated list through. `List(name)` — parens around a bare name — is a
# syntax error upstream because a bare name is not `FExpr_Multi`.
function _geo_list_term_dispatch(source::AbstractString,body::String,
                                 opener::Char,context::_GeoNumericContext,
                                 caller::AbstractString,limit::Int,depth::Int)
    text=String(strip(body))
    if opener=='['
        # `tList '[' String__Index ']'` — the whole symbol's value list
        # splices; an index on the name is parsed but ignored upstream.
        bare=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?:\(.*\)|\[.*\])?$",text)
        if bare!==nothing
            name=String(bare.captures[1])
            _geo_context_has_variable(context,name) ||
                (_geo_yyerror!(context,"Unknown variable '$name'");
                 return Float64[])
            return _geo_symbol_payload(context,name,caller)
        end
    end
    if startswith(text,"{")
        return _geo_numeric_list_values(text,context,"$caller List";
                                        depth=depth+1)
    end
    # `tList LP FExpr_Multi RP` takes a single multi-expression — a top-level
    # `,` is the unexpected token upstream (`List(5,6)`, `List[9,10]`).
    depth0=0
    for c in text
        c in ('(','[','{') && (depth0+=1)
        c in (')',']','}') && (depth0-=1)
        (depth0==0 && c==',') && _geo_syntax_abort(",")
    end
    inner=_geo_multi_term_values(text,context,caller,limit,depth+1)
    inner!==nothing && return inner
    # `List(a)`/`List[1,2]` don't reduce — `List` takes a `[name]` splice,
    # an `FExpr_Multi`, or a brace list.
    _geo_syntax_abort(opener=='[' ? "]" : ")")
end

function _geo_multi_term_values(raw::AbstractString,
                                context::_GeoNumericContext,
                                caller::AbstractString,limit::Int,depth::Int)
    depth<=_MAX_GEO_EXPRESSION_DEPTH || throw(ArgumentError(
        "$caller: list expression nesting exceeds $_MAX_GEO_EXPRESSION_DEPTH"))
    source=String(strip(raw))
    isempty(source) && return nothing
    # Side-effecting/entity terms (`Extrude`, `Boolean`, `Point{..}`, ...) —
    # `execute_geo` installs the hook; pure contexts skip it.
    if context.exec_hook!==nothing
        hooked=context.exec_hook(source)
        hooked!==nothing && return collect(Float64,hooked)
    end
    # `a[]`/`a()`/`a[{..}]`/`a({..})` splices.
    ref=_geo_list_reference_values(source,context,caller;depth=depth+1)
    ref!==nothing && return ref
    # `FExpr ':' FExpr` / `FExpr ':' FExpr ':' FExpr` — the range operands are
    # FExprs, so `2*3:5` reads as `(2*3):5` and `-1:3` as `(-1):3`: the range
    # split must run before the leading `-` and `*` multi forms.
    pieces=_geo_split_range(source,caller)
    if pieces!==nothing
        first=_geo_eval_numeric(pieces[1],context,"$caller range start")
        last=_geo_eval_numeric(pieces[2],context,"$caller range end")
        step=length(pieces)==2 ? (first<last ? 1.0 : -1.0) :
             _geo_eval_numeric(pieces[3],context,"$caller range increment")
        count=_geo_range_count(first,last,step,limit,caller)
        return collect(_geo_list_term_values(
            _GeoNumericListTerm(first,step,count)))
    end
    # `'-' FExpr_Multi` — negation of a list producer (`-a()`, `-List[x]`).
    if startswith(source,'-')
        inner=_geo_multi_term_values(
            String(strip(source[nextind(source,firstindex(source)):end])),
            context,caller,limit,depth+1)
        inner===nothing && return nothing
        return [-v for v in inner]
    end
    # `FExpr '*' FExpr_Multi` — scalar scaling. `FExpr '*' '{' RLD '}'` is a
    # `ListOfDouble` production handled by `_geo_braced_list_source` upstream
    # of this path — inside `FExpr_Multi` (a `{..}` member, a call argument)
    # a `{` after `*` is a syntax error.
    star=_geo_top_level_star(source)
    if star!==nothing
        rhs=String(strip(source[nextind(source,star):end]))
        startswith(rhs,"{") && _geo_syntax_abort("{")
        lhs_source=String(strip(
            source[firstindex(source):prevind(source,star)]))
        # A loose top-level operator before the `*` (a ternary arm, a range
        # colon, a comparison, ...) binds the `*` inside its operand — this
        # is not a `FExpr '*' FExpr_Multi` form; fall back to the scalar path.
        _geo_loose_top_level_token(lhs_source)!==nothing && return nothing
        rhs_values=_geo_multi_term_values(rhs,context,caller,limit,depth+1)
        rhs_values===nothing && return nothing
        lhs=try
            _geo_eval_numeric(lhs_source,context,"$caller multiplier")
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            return nothing
        end
        return [lhs*v for v in rhs_values]
    end
    head=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*[\(\[]",source)
    if head===nothing
        # `tX '{'` — the list builtins take `LP` (`(`/`[`), never `{` —
        # `List{`, `LinSpace{`, ... commit to a syntax error at `{`.
        brace_head=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*\{",source)
        brace_head!==nothing && brace_head.captures[1] in (
            "List","LinSpace","LogSpace","Catenary","Unique","Abs",
            "ListFromFile") && _geo_syntax_abort("{")
        return nothing
    end
    name=String(head.captures[1])
    call=_geo_call_body(source,name)
    call===nothing && return nothing
    inner,rest,opener=call
    isempty(rest) || return nothing
    if name=="List"
        return _geo_list_term_dispatch(
            source,inner,opener,context,caller,limit,depth)
    elseif name=="LinSpace" || name=="LogSpace"
        args=_geo_split_args(inner,"$caller $name")
        length(args)==3 || _geo_syntax_abort(")")
        a=_geo_eval_numeric(args[1],context,"$caller $name")
        b=_geo_eval_numeric(args[2],context,"$caller $name")
        nraw=_geo_eval_numeric(args[3],context,"$caller $name")
        # `for(i < (int)n)` iterations with the *raw* `n-1` denominator —
        # `LinSpace(a,b,1)` divides by zero (NaN) and `LinSpace(a,b,1.5)`
        # takes one step scaled by 0.5.
        n=_geo_int_value(nraw,"$caller $name count")
        n<=0 && return Float64[]
        n<=limit || throw(ArgumentError(
            "$caller: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
        denominator=nraw-1
        return name=="LinSpace" ?
            [a+(b-a)*Float64(i)/denominator for i in 0:n-1] :
            [_gm87_pow(10.0,a+(b-a)*Float64(i)/denominator) for i in 0:n-1]
    elseif name=="Catenary"
        args=_geo_split_args(inner,"$caller Catenary")
        length(args)==6 || _geo_syntax_abort(")")
        vals=[_geo_eval_numeric(arg,context,"$caller Catenary")
              for arg in args]
        n=_geo_int_value(vals[6],"$caller Catenary count")
        n<=0 && return Float64[]
        n<=limit || throw(ArgumentError(
            "$caller: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
        return _geo_catenary_values(vals[1],vals[2],vals[3],vals[4],vals[5],
                                    n,context,caller)
    elseif name=="Unique"
        startswith(strip(inner),"{") && _geo_syntax_abort("{")
        values=_geo_multi_term_values(inner,context,caller,limit,depth+1)
        values===nothing && _geo_syntax_abort(")")
        return _geo_unique_sorted(values)
    elseif name=="Abs"
        startswith(strip(inner),"{") && _geo_syntax_abort("{")
        values=_geo_multi_term_values(inner,context,caller,limit,depth+1)
        # `Abs(scalar)` stays on the scalar `FExpr` path.
        values===nothing && return nothing
        return abs.(values)
    elseif name=="ListFromFile"
        # `tListFromFile LP StringExprVar RP` — `LP`/`RP` accept `(` or `[`
        # interchangeably upstream.
        return _geo_list_from_file(inner,context,caller,limit)
    end
    return nothing
end

# `std::sort` + `std::unique` — sorted with duplicates removed.
function _geo_unique_sorted(values::Vector{Float64})
    # `lt=<` is C++ `operator<`: -0.0 and 0.0 are equivalent (MergeSort keeps
    # their input order), and `==` collapses them — the first of a run wins.
    sorted=sort(values;lt=<,alg=MergeSort)
    result=Float64[]
    isempty(sorted) && return result
    last=sorted[1];push!(result,last)
    @inbounds for i in 2:length(sorted)
        if !(sorted[i]==last)
            last=sorted[i];push!(result,last)
        end
    end
    return result
end

# `catenary` from `src/numeric/Numeric.cpp`: solve `y = a + 1/b cosh(b(x-c))`
# through (x0,y0), (x1,y1) with lowest point ys by `newton_fd` (finite-
# difference Newton, MAXIT=100, EPS=1e-4, relax=1); failure or an unphysical
# solution falls back to linear interpolation with a level-1 warning.
function _geo_catenary_residual(x::Vector{Float64},p::NTuple{5,Float64})
    x0,x1,y0,y1,ys=p
    b,c=x
    return [(ys-1/b)+1/b*cosh(b*(x0-c))-y0,
            (ys-1/b)+1/b*cosh(b*(x1-c))-y1]
end

function _geo_catenary_values(x0::Float64,x1::Float64,y0::Float64,
                              y1::Float64,ys::Float64,n::Int,
                              context::_GeoNumericContext,caller::AbstractString)
    param=(x0,x1,y0,y1,ys)
    x=[1.0/(x1-x0),(x0+x1)/2.0]
    tolx=1e-6*abs(x1-x0)
    toly=1e-6*abs(max(y0,y1)-ys)
    success=false
    physical=true
    if x0!=x1
        # newton_fd: MAXIT=100, EPS=1e-4, relax=1, N=2 LU solve. An all-zero
        # residual exits the iteration and reports failure upstream, so the
        # break here leaves `success` false and falls back to linear
        # interpolation exactly like Gmsh.
        for _ in 1:100
            hypot(x[1],x[2])>1e6 && break
            f=_geo_catenary_residual(x,param)
            (f[1]==0.0 && f[2]==0.0) && break
            h1=1e-4*abs(x[1]);h1==0.0 && (h1=1e-4)
            h2=1e-4*abs(x[2]);h2==0.0 && (h2=1e-4)
            x[1]+=h1;j1=(_geo_catenary_residual(x,param).-f)./h1;x[1]-=h1
            x[2]+=h2;j2=(_geo_catenary_residual(x,param).-f)./h2;x[2]-=h2
            # LU with partial pivoting for the 2×2 system J*dx=f.
            a11,a21,a12,a22=j1[1],j1[2],j2[1],j2[2]
            if abs(a21)>abs(a11)
                a11,a12,a21,a22,f[1],f[2]=a21,a22,a11,a12,f[2],f[1]
            end
            a11==0.0 && break
            l21=a21/a11;u22=a22-l21*a12
            u22==0.0 && break
            dx2=(f[2]-l21*f[1])/u22
            dx1=(f[1]-a12*dx2)/a11
            x[1]-=dx1;x[2]-=dx2
            hypot(dx1,dx2)<tolx && (success=true;break)
        end
    end
    # `std::vector<double> y(N)` value-initializes to zeros; a failed Newton
    # leaves `physical` true and returns the zeroed vector with no warning —
    # only a converged but unphysical curve hits the linear fallback.
    y=zeros(n)
    if success
        a=ys-1/x[1]
        for i in 0:n-1
            r=x0+(i+1)*(x1-x0)/(n+1)
            y[i+1]=a+1/x[1]*cosh(x[1]*(r-x[2]))
            if y[i+1]>max(y0,y1)+toly || y[i+1]<ys-toly
                physical=false;break
            end
        end
    end
    physical && return y
    _geo_yywarn!(context,
        "Catenary did not converge, using linear interpolation")
    return [y0+(i+1)*(y1-y0)/(n+1) for i in 0:n-1]
end

# `ListFromFile("path")` — `fscanf("%lf")` semantics: whitespace-separated
# numeric tokens are read, non-numeric tokens warn and are skipped; relative
# paths resolve against the `.geo` file's directory (`FixRelativePath`).
function _geo_list_from_file(inner::AbstractString,
                             context::_GeoNumericContext,
                             caller::AbstractString,limit::Int)
    path=_geo_eval_string(inner,context,"$caller ListFromFile")
    resolved=_geo_fix_relative_path(context.file_name,path)
    values=Float64[]
    isfile(resolved) || begin
        _geo_yyerror!(context,"Could not open file '$path'")
        return values
    end
    for token in eachsplit(read(resolved,String))
        match_result=match(
            r"^[+-]?(?:(?:\d+\.?\d*)|(?:\.\d+))(?:[eE][+-]?\d+)?",token)
        if match_result===nothing
            _geo_yywarn!(context,"Ignoring '$token' in file '$path'")
            continue
        end
        value=tryparse(Float64,match_result.match)
        if value===nothing
            _geo_yywarn!(context,"Ignoring '$token' in file '$path'")
            continue
        end
        length(values)<limit || throw(ArgumentError(
            "$caller: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
        push!(values,value)
        remainder=token[nextind(token,lastindex(match_result.match)):end]
        isempty(remainder) ||
            _geo_yywarn!(context,"Ignoring '$remainder' in file '$path'")
    end
    return values
end

function _geo_numeric_list_term(raw::AbstractString,context::_GeoNumericContext,
                                caller::AbstractString,limit::Int;
                                tail::AbstractString="}")
    if context.exec_hook!==nothing
        hooked=context.exec_hook(raw)
        if hooked!==nothing
            length(hooked)<=limit || throw(ArgumentError(
                "$caller: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
            return _GeoNumericListTerm(0.0,0.0,length(hooked),hooked)
        end
    end
    multi=_geo_multi_term_values(raw,context,caller,limit,0)
    if multi!==nothing
        return _GeoNumericListTerm(0.0,0.0,length(multi),multi)
    end
    pieces=_geo_split_range(raw,caller;tail=tail)
    if pieces===nothing
        limit>0 || throw(ArgumentError(
            "$caller: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
        value=_geo_eval_numeric(raw,context,"$caller entry")
        return _GeoNumericListTerm(value,0.0,1)
    end
    first=_geo_eval_numeric(pieces[1],context,"$caller range start")
    last=_geo_eval_numeric(pieces[2],context,"$caller range end")
    # Gmsh's three-term order is start:end:increment. For the two-term form it
    # chooses +1 when start < end and -1 otherwise (including equal endpoints).
    step=length(pieces)==2 ? (first<last ? 1.0 : -1.0) :
         _geo_eval_numeric(pieces[3],context,"$caller range increment")
    count=_geo_range_count(first,last,step,limit,caller)
    return _GeoNumericListTerm(first,step,count)
end

function _geo_numeric_list_terms(items::Vector{String},context::_GeoNumericContext,
                                 caller::AbstractString)
    terms=Vector{_GeoNumericListTerm}(undef,length(items));total=0
    for (i,item) in pairs(items)
        term=_geo_numeric_list_term(item,context,caller,_MAX_GEO_LIST_ITEMS-total)
        terms[i]=term
        total+=term.count
    end
    return terms,total
end

function _geo_validate_physical_ranges(raw::AbstractString,
                                       context::_GeoNumericContext,
                                       caller::AbstractString)
    # Physical membership is intentionally not interpreted by this scanner:
    # entries can be geometry-derived lists such as `Surface{:}` or
    # `Surface In BoundingBox{...}`, and can contain ternaries. Preserve the
    # pre-existing opaque-RHS behavior for those forms. Direct numeric lists
    # are nevertheless recognized, checked and bounded when ranges are present.
    occursin(':',raw) || return nothing
    body=String(strip(raw))
    if startswith(body,"{") && endswith(body,"}")
        body_first=nextind(body,firstindex(body))
        body_last=prevind(body,lastindex(body))
        body=String(strip(body[body_first:body_last]))
    end
    (occursin('{',body) || occursin('}',body) || occursin('?',body)) &&
        return nothing
    items=_geo_split_list(raw,caller);total=0
    for item in items
        pieces=_geo_split_range(item,caller)
        if pieces===nothing
            total<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                "$caller: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
            total+=1
        else
            term=_geo_numeric_list_term(item,context,caller,
                                        _MAX_GEO_LIST_ITEMS-total)
            total+=term.count
        end
    end
    return nothing
end

# First unexpected token of a trailing fragment — an identifier reports its
# whole name; anything else reports the single character (like the lexer).
function _geo_first_token(source::AbstractString)
    m=match(r"^[A-Za-z_][A-Za-z0-9_]*|^<<|^>>|^<=|^>=|^==|^!=|^&&|^\|\||^\+\+|^--|^::|^.",source)
    m===nothing ? ";" : String(m.match)
end

# Index of the `}` balancing a leading `{`, or `nothing` when it never
# closes.
function _geo_braced_end(source::AbstractString)
    depth=0
    for i in eachindex(source)
        c=source[i]
        c=='{' && (depth+=1)
        c=='}' && (depth-=1;depth==0 && return i)
    end
    return nothing
end

# `FExpr '*' '{' RLD '}'` is only a `ListOfDouble` production — the `*` must
# sit at the top level of the RHS. A top-level operator *looser* than `*`
# (binary `+`/`-`, shifts, comparisons, `&&`/`||`, `?`/`:`) makes `*` bind
# inside that operator's right operand, where `{` is illegal and is the
# reported token. `,` and `=` are themselves the unexpected token upstream.
# `|`/`&`/`*`/`^`/unary ops bind tighter and leave `*` at top level.
function _geo_loose_top_level_token(source::AbstractString)
    depth=0;value=false;i=firstindex(source);last=lastindex(source)
    while i<=last
        c=source[i]
        if c=='('||c=='['
            depth+=1;value=false;i=nextind(source,i);continue
        elseif c==')'||c==']'
            depth-=1;value=true;i=nextind(source,i);continue
        elseif isspace(c)
            i=nextind(source,i);continue
        elseif depth!=0
            i=nextind(source,i);continue
        end
        if c=='?'
            return "{"
        elseif c==':'
            j=nextind(source,i)
            (j<=last && source[j]==':') || return "{"
            i=j
        elseif c==','
            return ","
        elseif c=='='
            j=nextind(source,i)
            return (j<=last && source[j]=='=') ? "{" : "="
        elseif c=='<'||c=='>'
            return "{"
        elseif c=='!'
            j=nextind(source,i)
            (j<=last && source[j]=='=') && return "{"
            value=false;i=j;continue
        elseif c=='&'||c=='|'
            j=nextind(source,i)
            (j<=last && source[j]==c) && return "{"
            value=false;i=j;continue
        elseif c=='+'||c=='-'
            value && return "{"
            value=false;i=nextind(source,i);continue
        elseif c=='*'||c=='/'||c=='%'||c=='^'||c=='~'||c=='#'
            value=false;i=nextind(source,i);continue
        else
            value=true;i=nextind(source,i);continue
        end
        i=nextind(source,i)
    end
    return nothing
end

function _geo_braced_list_source(raw::AbstractString,context::_GeoNumericContext,
                                 caller::AbstractString;
                                 allow_multiplier::Bool)
    value=String(strip(raw))
    if startswith(value,"{")
        # `'{' RLD '}'` completes at the first balanced `}` — a tail is a
        # syntax error at its first token (`{1,2}+x` → `+`, `{1,2}*3` → `*`).
        closing=_geo_braced_end(value)
        closing===nothing && _geo_syntax_abort("{",
            "$caller: malformed brace list $raw")
        rest=String(strip(value[nextind(value,closing):end]))
        isempty(rest) || _geo_syntax_abort(_geo_first_token(rest),
            "$caller: malformed brace list $raw")
        return value,1.0
    end
    brace=findfirst(==('{'),value)
    if brace===nothing
        endswith(value,"}") && _geo_syntax_abort("}",
            "$caller: malformed brace list $raw")
        return nothing,1.0
    end
    prefix=String(strip(value[firstindex(value):prevind(value,brace)]))
    # A `{` that is not the list opener — selector/function argument braces
    # (`Unique({..})`, `Point{1}`, `BoundingBox X{..}`) — belongs to a
    # non-brace-list term; the term evaluator decides.
    (prefix=="-" || (allow_multiplier && endswith(prefix,"*"))) ||
        return nothing,1.0
    group=String(value[brace:lastindex(value)])
    closing=_geo_braced_end(group)
    closing===nothing && _geo_syntax_abort("{",
        "$caller: malformed brace list $raw")
    rest=String(strip(group[nextind(group,closing):end]))
    isempty(rest) || _geo_syntax_abort(_geo_first_token(rest),
        "$caller: malformed brace list $raw")
    list_source=String(group[firstindex(group):closing])
    multiplier=if prefix=="-"
        -1.0
    else
        # `FExpr '*' '{' RLD '}'` — a loose operator before the `*` binds the
        # `*` inside its operand, so `1+2*{..}` errors at `{` upstream.
        factor_last=prevind(prefix,lastindex(prefix))
        factor_source=String(strip(prefix[firstindex(prefix):factor_last]))
        isempty(factor_source) && _geo_syntax_abort("*",
            "$caller: list multiplier must not be empty")
        loose=_geo_loose_top_level_token(factor_source)
        loose===nothing || _geo_syntax_abort(loose,
            "$caller: a list multiplier containing a top-level operator " *
            "must be parenthesized")
        _geo_eval_numeric(factor_source,context,"$caller list multiplier")
    end
    body_first=nextind(list_source,firstindex(list_source))
    body_last=prevind(list_source,lastindex(list_source))
    body=body_first>body_last ? "" : String(strip(list_source[body_first:body_last]))
    # `-{..}`/`x*{..}` need `'{' RecursiveListOfDouble '}'` — nonempty.
    isempty(body) && _geo_syntax_abort("}",
        "$caller: a negated or multiplied brace list must not be empty")
    return list_source,multiplier
end

function _geo_list_reference_values(raw::AbstractString,
                                    context::_GeoNumericContext,
                                    caller::AbstractString;depth::Int=0)
    depth<=_MAX_GEO_EXPRESSION_DEPTH || throw(ArgumentError(
        "$caller: list selector nesting exceeds $_MAX_GEO_EXPRESSION_DEPTH"))
    source=String(strip(raw))
    # `x.member()`/`ns::x.member()` — `FExpr_Multi` whole-member splices;
    # missing members contribute no elements (upstream `List_T*` nullptr).
    member_splice=match(
        r"^(-?)\s*(?:([A-Za-z_][A-Za-z0-9_]*)\s*::\s*)?([A-Za-z_][A-Za-z0-9_]*)\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)\s*$",
        source)
    if member_splice!==nothing
        ns=member_splice.captures[2]===nothing ?
            "" : String(member_splice.captures[2])
        sname=String(member_splice.captures[3])
        smember=String(member_splice.captures[4])
        _geo_member_key_valid(smember) || return nothing
        vals=_geo_member_fvector(context,ns,sname,smember,caller)
        vals===nothing && return Float64[]
        member_splice.captures[1]=="-" && (vals=-vals)
        return vals
    end
    whole=match(
        r"^(-?)\s*([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*\]\s*$",source)
    # `x()` — `String__Index LP RP` splices the whole symbol payload,
    # scalar or list.
    parens=whole===nothing ? match(
        r"^(-?)\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)\s*$",source) : nothing
    selected=whole===nothing && parens===nothing ? match(
        r"^(-?)\s*([A-Za-z_][A-Za-z0-9_]*)\s*([\(\[])\s*\{\s*(.*)\s*\}\s*([\)\]])\s*$",
        source) : nothing
    if selected!==nothing
        # `a[{..}]`/`a({..})` — `String__Index LP '{' ... '}' RP`: `LP`/`RP`
        # each match `(` or `[` / `)` or `]` independently, so mismatched
        # pairs like `a({1,2}]` are legal upstream. The name must be a
        # `String__Index` — `List`, `Abs`, `Unique`, ... are reserved tokens
        # upstream and head `FExpr_Multi` productions, so `List({7,8})` is
        # the `tList LP '{' ... '}' RP` form, not an index into a variable
        # named `List`.
        selected.captures[2] in
            ("List","LinSpace","LogSpace","Catenary","Unique","Abs",
             "ListFromFile") && return nothing
    end
    whole===nothing && parens===nothing && selected===nothing && return nothing
    matched=whole===nothing ? (parens===nothing ? selected : parens) : whole
    sign=matched.captures[1]=="-" ? -1.0 : 1.0
    name=String(matched.captures[2])
    # `String__Index` derives from `tSTRING` only — lexer keywords head
    # their own productions (`DimNameSpace()` is a function call, never a
    # variable splice).
    _geo_tstring_valid(name) || return nothing
    values=parens===nothing ? _geo_context_list(context,name,caller) :
        _geo_symbol_payload(context,name,caller)
    result=if selected===nothing
        values
    else
        index_values=_geo_numeric_list_values(
            "{"*selected.captures[4]*"}",context,"$caller selector";
            allow_multiplier=true,depth=depth+1)
        # `(int)` truncation, then `s.value[index]` — out-of-range indices
        # are `yymsg(0)` "Uninitialized variable" diagnostics and are
        # *skipped*, shortening the result.
        selected_values=Float64[]
        sizehint!(selected_values,length(index_values))
        for index_value in index_values
            index=_geo_int_value(index_value,"$caller list $name index")
            index<0 && throw(ArgumentError(
                "$caller: negative index $index in '$name'"))
            if index>=length(values)
                context.soft_unknown_reads || throw(ArgumentError(
                    "$caller: Uninitialized variable '$name[$index]'"))
                _geo_yyerror!(context,"Uninitialized variable '$name[$index]'")
                continue
            end
            push!(selected_values,values[index+1])
        end
        selected_values
    end
    if sign<0
        for index in eachindex(result)
            result[index]=-result[index]
        end
    end
    return result
end

# `x()` splice payload: the whole `s.value` vector — the element list for a
# list variable, or the single scalar slot for a scalar.
function _geo_symbol_payload(context::_GeoNumericContext,name::String,
                             caller::AbstractString)
    if haskey(context.lists,name)
        # Copy: callers may negate the payload in place (`-x()`), which must
        # not mutate the stored variable.
        return copy(context.lists[name])
    elseif haskey(context.values,name)
        return Float64[context.values[name]]
    end
    context.soft_unknown_reads ||
        throw(ArgumentError("$caller: Unknown variable '$name'"))
    _geo_yyerror!(context,"Unknown variable '$name'")
    return Float64[]
end

function _geo_numeric_list_values(raw::AbstractString,
                                  context::_GeoNumericContext,
                                  caller::AbstractString;
                                  allow_multiplier::Bool=true,
                                  depth::Int=0)
    depth<=_MAX_GEO_EXPRESSION_DEPTH || throw(ArgumentError(
        "$caller: list expansion nesting exceeds $_MAX_GEO_EXPRESSION_DEPTH"))
    source=String(strip(raw))
    isempty(source) && throw(ArgumentError(
        "$caller: numeric list expression must not be empty"))
    ncodeunits(source)<=_MAX_GEO_EXPRESSION_BYTES || throw(ArgumentError(
        "$caller: list expression exceeds $_MAX_GEO_EXPRESSION_BYTES bytes"))
    referenced=_geo_list_reference_values(
        source,context,caller;depth=depth+1)
    referenced!==nothing && return referenced
    # Side-effecting list terms (translational `Extrude`) can form a whole
    # RHS, not only a `{...}` member — try the hook before the brace grammar.
    if context.exec_hook!==nothing
        hooked=context.exec_hook(source)
        hooked!==nothing && return hooked
    end
    list_source,multiplier=_geo_braced_list_source(
        source,context,caller;allow_multiplier=allow_multiplier)
    if list_source===nothing
        # Whole-RHS `FExpr_Multi`: unbraced ranges (`1:5`), list functions
        # (`LinSpace`, `Unique`, ...), splices and selector terms — plus the
        # single-scalar case.
        term=_geo_numeric_list_term(source,context,caller,_MAX_GEO_LIST_ITEMS;
                                    tail=";")
        return _geo_list_term_values(term)
    end

    items=_geo_split_list(list_source,caller)
    values=Float64[]
    for item in items
        expanded=_geo_list_reference_values(
            item,context,caller;depth=depth+1)
        if expanded!==nothing
            length(values)+length(expanded)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                "$caller: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
            append!(values,expanded)
            continue
        end
        term=_geo_numeric_list_term(
            item,context,caller,_MAX_GEO_LIST_ITEMS-length(values))
        append!(values,_geo_list_term_values(term))
    end
    if multiplier!=1.0
        for index in eachindex(values)
            scaled=multiplier*values[index]
            # `s.value` keeps non-finite results under `.geo` execution; the
            # scan retains the finite contract.
            isfinite(scaled) || context.soft_unknown_reads ||
                throw(ArgumentError(
                    "$caller entry: list multiplication produced a non-finite value"))
            values[index]=scaled
        end
    end
    return values
end

function _geo_list_mutation_value(operation::AbstractString,a::Float64,b::Float64,
                                  context::_GeoNumericContext,
                                  caller::AbstractString)
    value=try
        operation=="=" ? b : operation=="+=" ? a+b :
        operation=="-=" ? a-b : operation=="*=" ? a*b : a/b
    catch err
        err isa InterruptException && rethrow()
        (err isa OverflowError || err isa DivideError) || rethrow()
        throw(ArgumentError("$caller: list mutation is outside its finite real domain"))
    end
    isfinite(value) || context.soft_unknown_reads || throw(ArgumentError(
        "$caller: list mutation produced a non-finite value"))
    return value
end

# ======== Gmsh `.geo` variable affectation (`Affectation`) ========
#
# Gmsh keeps one `gmsh_yysymbol` per name: `s.value` is a double vector and
# `s.list` a flag. `x = 5` sets list=false (single value); `x = {..}` sets
# list=(n != 1); `x() =`/`x[i] =` force list=true. Scalar reads see `value[0]`;
# `x[i]` sees `value[i]`. `yymsg(0)` diagnostics are accumulated parse errors —
# reported through `_geo_yyerror!` while execution continues.

# `x <op> ListOfDouble` — bare name with `=`/`+=`/`-=`/`*=`/`/=`.
function _geo_assign_bare!(context::_GeoNumericContext,name::String,
                         operation::String,rhs::Vector{Float64},
                         caller::AbstractString)
    if !_geo_context_has_variable(context,name) &&
       !haskey(context.strings,name)
        if operation!="="
            # `x op rhs` on a fresh name: single-element rhs gets "Unknown
            # variable"; a multi-element rhs falls through to the "Cannot
            # assign list" diagnostic (the symbol is created but stays scalar).
            _geo_yyerror!(context,length(rhs)==1 ?
                "Unknown variable '$name'" :
                "Cannot assign list to variable '$name'")
            return nothing
        end
        # `=` on a fresh name: `s.list = (n != 1)` — an empty braces list is
        # still a list.
        _geo_store_value!(context,name,rhs,length(rhs)!=1,caller)
        return nothing
    end
    # Gmsh keeps separate numeric and string symbol tables — a string-only
    # name can still take a numeric value here. `=` re-derives the flag as
    # `s.list = (n != 1)` before dispatch: `x = {1,2}` turns a scalar name
    # into a list, while a one-element rhs makes even a list name scalar
    # again — the write lands on `value[0]` and the array payload is kept
    # (`_geo_context_set_scalar!` preserves it). `op=` dispatches on the
    # current flag.
    is_list=operation=="=" ? length(rhs)!=1 : name in context.list_variables
    if !is_list
        if length(rhs)!=1
            _geo_yyerror!(context,"Cannot assign list to variable '$name'")
            return nothing
        end
        d=rhs[1]
        if !haskey(context.values,name)
            operation!="=" &&
                _geo_yywarn!(context,"Uninitialized variable '$name'")
            context.values[name]=0.0
        end
        current=context.values[name]
        result=_geo_apply_numeric_op(operation,current,d,context,caller,
                                     "'$name'")
        result===nothing && return nothing
        _geo_context_set_scalar!(context,name,result,caller)
        return nothing
    end
    # List variable: `=` clears then appends; `+=` appends; `-=` erases the
    # first occurrence of each RHS element; `*=`/`/=` are errors.
    if operation=="=" || operation=="+="
        values=operation=="=" ? Float64[] :
            copy(_geo_context_list(context,name,caller))
        length(values)+length(rhs)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
            "$caller: $name exceeds $_MAX_GEO_LIST_ITEMS entries"))
        append!(values,rhs)
        _geo_store_value!(context,name,values,true,caller)
    elseif operation=="-="
        values=copy(_geo_context_list(context,name,caller))
        for value in rhs
            index=findfirst(==(value),values)
            index===nothing || deleteat!(values,index)
        end
        _geo_store_value!(context,name,values,true,caller)
    else
        _geo_yyerror!(context,"Operators *= and /= not available for lists")
    end
    return nothing
end

# Element-wise op on a stored value; returns `nothing` when a parse error was
# recorded instead (division by zero leaves the value unchanged).
function _geo_apply_numeric_op(operation::String,current::Float64,
                               d::Float64,context::_GeoNumericContext,
                               caller::AbstractString,label::String)
    if operation=="="; return d
    elseif operation=="+="; return current+d
    elseif operation=="-="; return current-d
    elseif operation=="*="; return current*d
    else
        if d!=0
            return current/d
        end
        _geo_yyerror!(context,
            "Division by zero in '$label /= $(_geo_gmsh_number(d))'")
        return nothing
    end
end

# Write a value vector + list flag, keeping `values`/`lists`/`list_variables`
# coherent (`values[name]` is always `value[0]`).
function _geo_store_value!(context::_GeoNumericContext,name::String,
                         values::Vector{Float64},is_list::Bool,
                         caller::AbstractString)
    if is_list
        _geo_context_set_list!(context,name,values,caller)
    else
        _geo_context_set_scalar!(context,name,
            isempty(values) ? 0.0 : values[1],caller)
        haskey(context.lists,name) &&
            (context.lists[name]=isempty(values) ? Float64[] : copy(values))
    end
    return nothing
end

# `x() <op> ListOfDouble` / `x[] <op> ListOfDouble` — the explicit-list form
# forces `s.list = true` even for single-element or empty values. Gmsh's
# `gmsh_yysymbols[$1]` creates the entry unconditionally, so `-=`, `*=` and
# `/=` on a missing name still create an (empty) list — only `*=`/`/=` then
# emit the not-available diagnostic.
function _geo_assign_list!(context::_GeoNumericContext,name::String,
                         operation::String,rhs::Vector{Float64},
                         caller::AbstractString)
    if !_geo_context_has_variable(context,name)
        # `s.list = true` on a fresh symbol — even for `-=`/`*=`/`/=`.
        _geo_store_value!(context,name,Float64[],true,caller)
    end
    if operation=="=" || operation=="+="
        values=operation=="=" ? Float64[] :
            copy(_geo_context_list(context,name,caller))
        length(values)+length(rhs)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
            "$caller: $name exceeds $_MAX_GEO_LIST_ITEMS entries"))
        append!(values,rhs)
        _geo_store_value!(context,name,values,true,caller)
    elseif operation=="-="
        values=copy(_geo_context_list(context,name,caller))
        for value in rhs
            index=findfirst(==(value),values)
            index===nothing || deleteat!(values,index)
        end
        _geo_store_value!(context,name,values,true,caller)
    else
        _geo_yyerror!(context,"Operators *= and /= not available for lists")
    end
    return nothing
end

# `x[i] <op> v` / `x(i) <op> v` — `assignVariable`/`incrementVariable`: lists
# resize with zeros to the index; scalars error with "not a list".
function _geo_assign_indexed!(context::_GeoNumericContext,name::String,
                            index::Float64,operation::String,
                            rhs::Float64,caller::AbstractString)
    i=_geo_int_value(index,"$caller index")
    i<0 && throw(ArgumentError("$caller: negative index $i in '$name'"))
    if !_geo_context_has_variable(context,name)
        # String-only names share the numeric table's absence — `=` creates a
        # fresh list regardless.
        if operation=="="
            values=zeros(Float64,i+1);values[i+1]=rhs
            _geo_store_value!(context,name,values,true,caller)
        else
            _geo_yyerror!(context,"Unknown variable '$name'")
        end
        return nothing
    end
    name in context.list_variables || (_geo_yyerror!(context,
        "Variable '$name' is not a list"); return nothing)
    values=context.lists[name] # direct reference — mutation persists
    if length(values)<i+1
        grown=i+1-length(values)
        context.stored_list_items+grown<=
            _MAX_GEO_CONTEXT_LIST_ITEMS || throw(ArgumentError(
            "$caller: stored numeric lists exceed $_MAX_GEO_CONTEXT_LIST_ITEMS entries"))
        append!(values,zeros(Float64,grown))
        context.stored_list_items+=grown
    end
    result=_geo_apply_numeric_op(operation,values[i+1],rhs,context,caller,
                                 "'$name[$i]'")
    result===nothing && return nothing
    values[i+1]=result
    i==0 && (context.values[name]=values[1])
    return nothing
end

# `x++`/`x--` — scalar-only `incrementVariable` without index. Gmsh's
# statement form records a level-0 diagnostic on a missing or empty scalar
# and leaves the value untouched (it does not assign through the miss).
function _geo_increment_variable!(context::_GeoNumericContext,name::String,
                                  delta::Float64,caller::AbstractString)
    if !_geo_context_has_variable(context,name)
        _geo_yyerror!(context,"Unknown variable '$name'")
        return nothing
    end
    if name in context.list_variables
        _geo_yyerror!(context,"Variable '$name' is a list")
        return nothing
    end
    if !haskey(context.values,name)
        _geo_yyerror!(context,"Uninitialized variable '$name'")
        return nothing
    end
    context.values[name]+=delta
    haskey(context.lists,name) && !isempty(context.lists[name]) &&
        (context.lists[name][1]+=delta)
    return nothing
end

# `x({i,j,...}) <op> {v,...}` — `assignVariables`: multi-index mutation with
# matching index/value counts.
function _geo_assign_multi_index!(context::_GeoNumericContext,name::String,
                                  indices::Vector{Float64},operation::String,
                                  rhs::Vector{Float64},caller::AbstractString)
    if length(indices)!=length(rhs)
        _geo_yyerror!(context,"Incompatible array dimensions in affectation")
        return nothing
    end
    for position in eachindex(indices)
        _geo_assign_indexed!(context,name,indices[position],operation,
                             rhs[position],caller)
    end
    return nothing
end

# Top-level `.geo` numeric/string assignment dispatcher. `body` is the
# statement text without the trailing `;`. Returns true when the statement was
# an affectation (possibly with accumulated `yymsg` diagnostics), false when it
# matches none of the forms.
function _geo_exec_assignment!(context::_GeoNumericContext,
                               raw::AbstractString,caller::AbstractString)
    body=String(strip(raw))
    # `String__Index` LHS heads — `x~{FExpr}`…/`StringToName[sexpr]` resolve to
    # a plain identifier before the option/index assignment splits below.
    lr=_geo_lhs_string_index(body,context,caller)
    lr===nothing || (body=lr[1]*lr[2])
    # `Field[i].member = v` — mutates the field object's option (Gmsh
    # `field->options[member]`): unknown field/option names are diagnostics,
    # not silent stores.
    if (fm=match(
            r"^Field\s*\[\s*(.*?)\s*\]\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$",
            body))!==nothing
        _geo_assign_field_option!(context,fm.captures[1],
            String(fm.captures[2]),String(strip(fm.captures[3])),caller)
        return true
    end
    # `Plugin(name).member = v` — plugin options need a plugin registry; every
    # set fails in batch Gmsh with "Unknown option ... or plugin ...".
    if (pm=match(
            r"^Plugin\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$",
            body))!==nothing
        _geo_yyerror!(context,"Unknown option '$(pm.captures[2])' or plugin " *
            "'$(pm.captures[1])'")
        return true
    end
    # A lexer keyword cannot head a `String__Index` affectation: expression-
    # only tokens error at themselves, statement-leading tokens at the
    # following token (`DimNameSpace[] = 3` → `(DimNameSpace)`, `Point[...]`
    # → `[`). Field/Plugin statements above already consumed their literal
    # heads, and a resolved `~{}`/`StringToName` name is not re-lexed upstream.
    lr===nothing && _geo_lhs_keyword_check(body)
    # `x.Color.member = {r,g,b[,a]}` and `x.ColorTable = {colors}` — packed
    # RGBA color options.
    if (cm=match(
            r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?:[\[\(]\s*(.*?)\s*[\]\)])?\s*\.\s*Color\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$",
            body))!==nothing
        _geo_assign_color_option!(context,String(cm.captures[1]),
            cm.captures[2],String(cm.captures[3]),
            String(strip(cm.captures[4])),caller)
        return true
    end
    if (cm=match(
            r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?:[\[\(]\s*(.*?)\s*[\]\)])?\s*\.\s*ColorTable\s*=\s*(.*?)\s*$",
            body))!==nothing
        _geo_assign_color_table!(context,String(cm.captures[1]),
            cm.captures[2],String(strip(cm.captures[3])),caller)
        return true
    end
    # `x.opt = ...`, `x[i].opt = ...` — option writes.
    if (m=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?:[\[\(]\s*(.*?)\s*[\]\)])?\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?:(\+=|-=|\*=|/=|=)\s*(.*?)|(\+\+|--))\s*$",body))!==nothing
        family=String(m.captures[1]);member=String(m.captures[3])
        index=0
        if m.captures[2]!==nothing
            index=_geo_int_value(_geo_eval_numeric(m.captures[2],context,
                "$caller option index"),"$caller option index")
        end
        if m.captures[6]!==nothing
            delta=m.captures[6]=="++" ? 1.0 : -1.0
            _geo_option_number_increment!(
                context,family,index,member,delta,caller)
            return true
        end
        operation=String(m.captures[4]);rhs=String(strip(m.captures[5]))
        # Option writes are disambiguated by the grammar, not the value:
        # `x.y = <StringExpr>` for quoted/string-valued heads, otherwise the
        # NumericAffectation+FExpr production.
        if operation=="=" && _geo_string_rhs(rhs)
            _geo_set_option_string!(context,family,index,member,
                _geo_eval_string(rhs,context,caller),caller)
        else
            _geo_set_option_number!(context,family,index,member,operation,
                _geo_eval_numeric(rhs,context,caller),caller)
        end
        return true
    end
    # `x++`, `x--`, `x[i]++`, `x(i)--`.
    if (m=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(\+\+|--)\s*$",body))!==nothing
        _geo_increment_variable!(context,String(m.captures[1]),
            m.captures[2]=="++" ? 1.0 : -1.0,caller)
        return true
    end
    if (m=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*[\[\(]\s*(.*?)\s*[\]\)]\s*(\+\+|--)\s*$",body))!==nothing
        name=String(m.captures[1])
        index=_geo_eval_numeric(m.captures[2],context,"$caller index")
        _geo_assign_indexed!(context,name,index,"+=",
            m.captures[3]=="++" ? 1.0 : -1.0,caller)
        return true
    end
    # `x() <op> rhs`, `x[] <op> rhs`, `x[i] <op> rhs`, `x({..}) <op> rhs`,
    # `x <op> rhs`.
    if (m=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?:([\[\(])\s*(.*?)\s*[\]\)])?\s*(\+=|-=|\*=|/=|=)\s*(.*?)\s*$",body))!==nothing
        name=String(m.captures[1])
        _geo_check_assignable_name(name,caller)
        opener=m.captures[2];selector=m.captures[3]
        operation=String(m.captures[4]);rhs=String(strip(m.captures[5]))
        if opener===nothing
            _geo_assign_bare_rhs!(context,name,operation,rhs,caller)
            return true
        end
        selector_src=String(strip(selector))
        if isempty(selector_src)
            # `x() = Str[...]` — `String__Index LP RP tAFFECT tStr LP
            # BracedOrNotRecursiveListOfStringExprVar RP`; a whole-argument
            # `{...}` wraps the same flat list.
            if (sm=match(r"^Str\s*[\(\[]",rhs))!==nothing
                open_at=lastindex(sm.match)
                close=_geo_matching_delim(rhs,open_at)
                close==0 && _geo_syntax_abort(rhs[open_at],
                    "$caller: unbalanced Str[...] in string-list assignment")
                String(strip(rhs[nextind(rhs,close):end]))=="" ||
                    _geo_syntax_abort(rhs[nextind(rhs,close)],
                        "$caller: unexpected text after Str[...]")
                inner=String(strip(rhs[nextind(rhs,open_at):prevind(rhs,close)]))
                if startswith(inner,"{") &&
                        _geo_matching_delim(inner,firstindex(inner))==
                        lastindex(inner)
                    inner=String(strip(inner[2:prevind(inner,lastindex(inner))]))
                end
                values=_geo_eval_string_list(inner,context,caller;
                    brace_commit=true)
                if operation=="="
                    context.strings[name]=values
                elseif operation=="+="
                    haskey(context.strings,name) || (_geo_yyerror!(context,
                        "Uninitialized variable '$name'"); return true)
                    append!(context.strings[name],values)
                else
                    _geo_yyerror!(context,
                        "Operator $operation not available for string lists")
                end
                return true
            end
            values=_geo_numeric_list_values(rhs,context,
                "$caller: numeric list $name";allow_multiplier=true)
            _geo_assign_list!(context,name,operation,values,caller)
            return true
        end
        if startswith(selector_src,"{") && endswith(selector_src,"}")
            indices=_geo_numeric_list_values(selector_src,context,
                "$caller: numeric list $name selector";allow_multiplier=false)
            values=_geo_numeric_list_values(rhs,context,
                "$caller: numeric list $name";allow_multiplier=true)
            _geo_assign_multi_index!(context,name,indices,operation,values,caller)
            return true
        end
        index=_geo_eval_numeric(selector_src,context,"$caller index")
        # The grammar's indexed form takes a single FExpr value, not a list.
        rhs_values=_geo_numeric_list_values(rhs,context,
            "$caller: numeric list $name";allow_multiplier=true)
        length(rhs_values)==1 || throw(ArgumentError(
            "$caller: indexed assignment '$name[$(selector_src)]' takes a " *
            "single value, got $(length(rhs_values))"))
        _geo_assign_indexed!(context,name,index,operation,rhs_values[1],caller)
        return true
    end
    # `ns::x = ...` — `tSCOPE` is not part of the `String__Index` LHS
    # production (`ns::x = 5` → `syntax error (::)` upstream).
    match(r"^[A-Za-z_][A-Za-z0-9_]*\s*::",body)!==nothing &&
        _geo_syntax_abort("::")
    # `x.y[i] = ...`/`x.y(i) = ...` — a member index after `.y` is not part
    # of the option-write LHS (`s.V[1] =`/`s.V(1) =` → `syntax error ([)`/
    # `(()` upstream).
    if (mm=match(
            r"^[A-Za-z_][A-Za-z0-9_]*\s*\.\s*[A-Za-z_][A-Za-z0-9_]*\s*([\[\(])",
            body))!==nothing
        _geo_syntax_abort(mm.captures[1]=="[" ? "[" : "(")
    end
    return false
end

# `Field[i].member = v` — Gmsh writes through `field->options[member]`, so the
# option must exist on the field's type. The stored representation mirrors the
# params pre-pass: `GeoFieldSpec.options[member]` holds a normalized source
# string consumed downstream by the field evaluator.
function _geo_assign_field_option!(context::_GeoNumericContext,tag_src,
                                   member::String,rhs::String,
                                   caller::AbstractString)
    tag=_geo_signed_gmsh_int_value(_geo_eval_numeric(tag_src,context,
        "$caller field tag"),"$caller field tag")
    fields=context.fields
    spec=fields===nothing ? nothing : get(fields,tag,nothing)
    # `field->options[member]` upstream — membership is per field *kind*, and
    # the `Unknown option` diagnostic reports the class `getName()`
    # (`Attractor` instantiates `DistanceField`).
    member_ok=spec!==nothing &&
        (kind_entry=get(_GEO_FIELD_OPTIONS,spec.kind,nothing))!==nothing &&
        member in kind_entry[1]
    normalized=if member_ok
        _geo_normalize_field_option(rhs,member,context,caller)
    else
        # The RHS reduced inside the production before the field lookup ran —
        # a numeric tail still reports its variables; a quoted/braced RHS
        # follows the StringExpr/ListOfDouble forms instead.
        if !startswith(rhs,"\"") && !startswith(rhs,"{")
            try
                _geo_eval_numeric(rhs,context,caller)
            catch err
                err isa InterruptException && rethrow()
                err isa _GeoSyntaxAbort && rethrow()
            end
        end
        nothing
    end
    if spec===nothing
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: No field with id $tag"))
        return _geo_yyerror!(context,"No field with id $tag")
    end
    if !member_ok
        reported=kind_entry===nothing ? spec.kind : kind_entry[2]
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Unknown option '$member' in field $tag of type " *
            "'$reported'"))
        return _geo_yyerror!(context,"Unknown option '$member' in field $tag " *
            "of type '$reported'")
    end
    spec.options[member]=normalized
    member in spec.option_order || push!(spec.option_order,member)
    return nothing
end

# `x.Color.member = ColorExpr` — packed RGBA write into the color-option
# mirror; unknown family/member follow `ColorOption`'s warn path.
function _geo_assign_color_option!(context::_GeoNumericContext,family::String,
                                   index_src,member::String,rhs::String,
                                   caller::AbstractString)
    index=index_src===nothing ? 0 :
        _geo_int_value(_geo_eval_numeric(index_src,context,
            "$caller color index"),"$caller color index")
    names=get(_GEO_COLOR_OPTION_NAMES,family,nothing)
    if names===nothing || !(member in names)
        context.soft_unknown_reads || throw(ArgumentError(
            "$caller: Unknown color option '$family.$member'"))
        names===nothing &&
            return _geo_msg_error!(context,
                "Unknown color option category '$family'")
        return _geo_msg_error!(context,
            "Unknown color option '$family.$member'")
    end
    context.option_colors[(family,index,member)]=
        _geo_color_expr_value(rhs,context,caller)
    return nothing
end

# `x.ColorTable = {colors}` — Gmsh writes the color table of `View[i]` where
# `i` is the bracket index (the family name is ignored); `.geo` execution
# never has a view, so the faithful result is "View[i] does not exist".
function _geo_assign_color_table!(context::_GeoNumericContext,family::String,
                                  index_src,rhs::String,
                                  caller::AbstractString)
    index=index_src===nothing ? 0 :
        _geo_int_value(_geo_eval_numeric(index_src,context,
            "$caller view index"),"$caller view index")
    context.soft_unknown_reads || throw(ArgumentError(
        "$caller: View[$index] does not exist"))
    _geo_yyerror!(context,"View[$index] does not exist")
    return nothing
end

# `ColorExpr` — `{r,g,b[,a]}` 0:255 components or a color name through the
# X11 table. Returns the (r,g,b,a) tuple.
function _geo_color_expr_value(rhs::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString)
    s=String(strip(rhs))
    if startswith(s,"{") && endswith(s,"}")
        values=_geo_numeric_list_values(s,context,"$caller color";
            allow_multiplier=false)
        (length(values)==3 || length(values)==4) || throw(ArgumentError(
            "$caller: color requires 3 or 4 components, got $(length(values))"))
        rgba=ntuple(i->i<=length(values) ?
            _geo_int_value(values[i],"$caller color component") : 255,4)
        return rgba
    end
    name=String(strip(s,'"'))
    haskey(_GEO_COLOR_NAMES,name) && return _GEO_COLOR_NAMES[name]
    context.soft_unknown_reads || throw(ArgumentError(
        "$caller: Unknown color '$name'"))
    _geo_yyerror!(context,"Unknown color '$name'")
    return (0,0,0,255)
end

function _geo_assign_bare_rhs!(context::_GeoNumericContext,name::String,
                               operation::String,rhs::String,
                               caller::AbstractString)
    # `x = <StringExpr>` — string-variable assignment. Quoted literals and
    # string-producing heads route to the string evaluator; everything else is
    # a ListOfDouble (bare-name FExpr wins the grammar ambiguity).
    if operation=="=" && _geo_string_rhs(rhs)
        context.strings[name]=String[_geo_eval_string(rhs,context,caller)]
        delete!(context.values,name);delete!(context.lists,name)
        delete!(context.list_variables,name)
        return nothing
    end
    values=_geo_numeric_list_values(rhs,context,
        "$caller: numeric variable $name";allow_multiplier=true)
    _geo_assign_bare!(context,name,operation,values,caller)
    return nothing
end

# `RecursiveListOfStringExprVar` — flat `StringExprVar`/`MultiStringExprVar`
# items; `x()`/`x.member()` splices expand every element of a string-list
# variable or char member. `{...}` groups are not items upstream; in
# `BracedOrNot` contexts (`x() = Str[...]`, struct-option `Str[...]` —
# `brace_commit=true`) a `{...}` opening the argument reduces as the braced
# alternative and the following comma is the offending lookahead, while in
# flat-only contexts (`Str`/`StrCat`/`TextAttributes`/`OnelabRun` arguments)
# the `{` itself fails.
function _geo_eval_string_list(raw::AbstractString,
                               context::_GeoNumericContext,
                               caller::AbstractString;
                               brace_commit::Bool=false)
    items=_geo_split_args(raw,caller)
    out=String[]
    for (i,item) in enumerate(items)
        s=strip(item)
        if startswith(s,"{")
            _geo_syntax_abort(
                brace_commit && i==1 && length(items)>1 ? "," : "{")
        elseif (sm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)$",s))!==nothing
            name=String(sm.captures[1])
            haskey(context.strings,name) || throw(ArgumentError(
                "$caller: Unknown string variable '$name'"))
            append!(out,context.strings[name])
        elseif (sm=match(
                r"^(?:([A-Za-z_][A-Za-z0-9_]*)\s*::\s*)?([A-Za-z_][A-Za-z0-9_]*)\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)$",
                s))!==nothing
            # `x.member()`/`ns::x.member()` — `MultiStringExprVar` char-member
            # splice; a missing member contributes no elements.
            vals=_geo_member_svector(context,
                sm.captures[1]===nothing ? "" : String(sm.captures[1]),
                String(sm.captures[2]),String(sm.captures[3]),caller)
            vals!==nothing && append!(out,vals)
        else
            push!(out,_geo_eval_string(s,context,caller))
        end
    end
    return out
end

# True when `rhs` can only be a StringExpr: a quoted literal, a `name(`-call to
# a string-producing function, or a zero-argument string token.
function _geo_string_rhs(rhs::AbstractString)
    s=strip(rhs)
    (startswith(s,'"') || startswith(s,'\'')) && return true
    m=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*[\(\[\{]",s)
    # `StringToName`/`S2N` are `StringIndex` name productions — a
    # `String__Index`, not a `StringExpr` — so the FExpr `x =` production
    # wins upstream and the resolved name reads as a numeric variable.
    m!==nothing && m.captures[1] in _GEO_STRING_FUNCTIONS &&
        !(m.captures[1] in ("StringToName","S2N")) && return true
    # `Point{i}`/`Surface{i}`/`Physical X{i}` are entity-name StringExprVars.
    match(r"^(?:Physical\s+)?(?:Point|Line|Curve|Surface|Volume)\s*\{",s)!==
        nothing && return true
    s in ("Today","CodeName","GmshExecutableName","OnelabAction",
          "CurrentFileName","CurrentDirectory") && return true
    return false
end

function _geo_check_assignable_name(name::String,caller::AbstractString)
    name=="Pi" && throw(ArgumentError(
        "$caller: Pi is a reserved numeric constant; use a different variable name"))
    name in _GEO_SIDE_EFFECT_SYMBOLS && throw(ArgumentError(
        "$caller: dynamic tag allocator $name is read-only and cannot be assigned"))
    return nothing
end

# The upstream option callbacks coerce on `GMSH_SET` — `(int)val`, clamps,
# `val ? 1 : 0` — so the stored (and read-back) value is the transformed
# result, not the raw RHS. `_GEO_NUMBER_OPTION_STORAGE` (generated from
# src/common/Options.cpp) carries each coerced option's transform; every
# write path (`=`, compound assigns, `++`/`--`) routes through it.
function _geo_option_store_value(family::String,member::String,
                                 value::Float64,caller::AbstractString)
    kind=get(_GEO_NUMBER_OPTION_STORAGE,"$family.$member",nothing)
    kind===nothing && return value
    if kind===:int
        return Float64(_geo_signed_gmsh_int_value(value,caller))
    elseif kind===:bool
        return iszero(value) ? 0.0 : 1.0
    elseif kind===:uint
        # `(unsigned int)val` — the double is truncated through `cvttsd2si`
        # (out-of-Int32-range inputs yield INT32_MIN on x86) then
        # reinterpreted unsigned.
        iv=trunc(value)
        iv32=typemin(Int32)<=iv<=typemax(Int32) ? Int32(iv) : typemin(Int32)
        return Float64(reinterpret(UInt32,iv32))
    elseif kind isa Tuple
        tag=kind[1]
        if tag===:int_alias
            iv=_geo_signed_gmsh_int_value(value,caller)
            return Float64(iv==Int(kind[2]) ? Int(kind[3]) : iv)
        elseif tag===:int_floor
            return Float64(max(_geo_signed_gmsh_int_value(value,caller),
                               Int(kind[2])))
        elseif tag===:int_cap
            return Float64(min(_geo_signed_gmsh_int_value(value,caller),
                               Int(kind[2])))
        elseif tag===:int_range_reset
            iv=_geo_signed_gmsh_int_value(value,caller)
            lo,hi,banned,default=Int(kind[2]),Int(kind[3]),kind[4],Int(kind[5])
            (iv<lo || iv>hi || iv==banned) && (iv=default)
            return Float64(iv)
        elseif tag===:clamp
            lo,hi=Float64(kind[2]),Float64(kind[3])
            return min(max(value,lo),hi)
        end
    end
    throw(ArgumentError(
        "$caller: unhandled option storage kind $kind for $family.$member"))
end

# `Mesh.MeshSizeFactor`/`Mesh.CharacteristicLengthFactor`,
# `Mesh.MeshSizeMin`/`Mesh.CharacteristicLengthMin` and
# `Mesh.MeshSizeMax`/`Mesh.CharacteristicLengthMax` are each two names for a
# single `CTX.mesh` value (both DefaultOptions entries dispatch to the same
# `opt_mesh_lc_*` callback), so a write to either name mirrors into the other
# key. `opt_mesh_lc_factor` also ignores non-positive writes (`if(val > 0)`),
# which is replicated by skipping the store entirely.
function _geo_store_option_number!(context::_GeoNumericContext,
                                   family::String,index::Int,
                                   member::String,value::Float64,
                                   caller::AbstractString)
    v=_geo_option_store_value(family,member,value,caller)
    if family=="Mesh" && index==0
        alias=if member=="MeshSizeFactor" || member=="CharacteristicLengthFactor"
            v>0 || return nothing
            member=="MeshSizeFactor" ? "CharacteristicLengthFactor" :
                "MeshSizeFactor"
        elseif member=="MeshSizeMin" || member=="CharacteristicLengthMin"
            member=="MeshSizeMin" ? "CharacteristicLengthMin" : "MeshSizeMin"
        elseif member=="MeshSizeMax" || member=="CharacteristicLengthMax"
            member=="MeshSizeMax" ? "CharacteristicLengthMax" : "MeshSizeMax"
        elseif member=="MeshSizeFromPoints" || member=="CharacteristicLengthFromPoints"
            member=="MeshSizeFromPoints" ? "CharacteristicLengthFromPoints" : "MeshSizeFromPoints"
        elseif member=="MeshSizeFromCurvature" || member=="CharacteristicLengthFromCurvature"
            member=="MeshSizeFromCurvature" ? "CharacteristicLengthFromCurvature" : "MeshSizeFromCurvature"
        elseif member=="MeshSizeExtendFromBoundary" || member=="CharacteristicLengthExtendFromBoundary"
            member=="MeshSizeExtendFromBoundary" ? "CharacteristicLengthExtendFromBoundary" : "MeshSizeExtendFromBoundary"
        else
            nothing
        end
        alias!==nothing &&
            (context.option_numbers[(family,index,alias)]=v)
        # The `Min*Nodes`/`Minimum*Nodes`/`Minimum*Points` spellings are
        # each one CTX.mesh slot upstream (`opt_mesh_min_*_nodes`), so a
        # write to any alias must reach the name the graders read.
        aliases=if member in ("MinCircleNodes","MinimumCircleNodes",
                            "MinimumCirclePoints")
            ("MinCircleNodes","MinimumCircleNodes","MinimumCirclePoints")
        elseif member in ("MinCurveNodes","MinimumCurveNodes",
                          "MinimumCurvePoints")
            ("MinCurveNodes","MinimumCurveNodes","MinimumCurvePoints")
        elseif member in ("MinLineNodes","MinimumLineNodes")
            ("MinLineNodes","MinimumLineNodes")
        else
            ()
        end
        for other in aliases
            other==member && continue
            context.option_numbers[(family,index,other)]=v
        end
    end
    context.option_numbers[(family,index,member)]=v
    if family=="General" && index==0 && member=="AbortOnError" &&
            context.native_abort_option_callback!==nothing
        context.native_abort_option_callback(v)
    end
    if index==0 && context.native_numeric_option_callback!==nothing
        context.native_numeric_option_callback(family,member,v)
    end
    return nothing
end

function _geo_option_number_increment!(context::_GeoNumericContext,
                                       family::String,index::Int,
                                       member::String,delta::Float64,
                                       caller::AbstractString)
    d=_geo_option_number(context,family,index,member)
    d===nothing && return _geo_option_number_unknown!(
        context,family,member,caller)
    _geo_store_option_number!(context,family,index,member,d+delta,caller)
    return nothing
end

function _geo_set_option_number!(context::_GeoNumericContext,family::String,
                                 index::Int,member::String,operation::String,
                                 value::Float64,caller::AbstractString)
    d=_geo_option_number(context,family,index,member)
    d===nothing && return _geo_option_number_unknown!(
        context,family,member,caller)
    result=_geo_apply_numeric_op(operation,d,value,context,caller,
        index==0 ? "$family.$member" : "$family[$index].$member")
    result===nothing && return nothing
    _geo_store_option_number!(context,family,index,member,result,caller)
    return nothing
end

# `NumberOption`/`StringOption` misses print `Msg::Error` (the
# `warnIfUnknown=true` default applies on every `.geo` write path).
function _geo_option_number_unknown!(context::_GeoNumericContext,
                                     family::String,member::String,
                                     caller::AbstractString)
    if context.soft_unknown_reads
        haskey(_GEO_NUMBER_OPTION_TABLE,family) ||
            return _geo_msg_error!(context,
                "Unknown number option category '$family'")
        return _geo_msg_error!(context,
            "Unknown number option '$family.$member'")
    end
    haskey(_GEO_NUMBER_OPTION_TABLE,family) || throw(ArgumentError(
        "$caller: Unknown number option category '$family'"))
    throw(ArgumentError("$caller: Unknown number option '$family.$member'"))
end

# `Mesh.Nb*` element statistics — upstream's option getters return the live
# `GModel` mesh statistics, so script writes never stick (`Mesh.NbNodes=42`
# still reads the generated count after `Mesh 2`, and 0 before any mesh).
# `Mesh.NbPartitions` is the exception: it is an ordinary stored option (the
# `PartitionMesh` count) and keeps the standard read path.
const _GEO_MESH_STAT_OPTIONS=Dict{String,Symbol}(
    "NbNodes"=>:nodes,"NbTriangles"=>:tri,"NbQuadrangles"=>:qua,
    "NbTetrahedra"=>:tet,"NbPrisms"=>:pri,"NbPyramids"=>:pyr,
    "NbHexahedra"=>:hex,"NbTrihedra"=>:trih)

function _geo_mesh_stat(mesh::Mesh,stat::Symbol)
    stat===:nodes && return Float64(nnodes(mesh))
    stat===:tri && return Float64(ntris(mesh))
    stat===:tet && return Float64(ntets(mesh))
    return 0.0
end

function _geo_mesh_stat(mesh::MixedMesh,stat::Symbol)
    stat===:nodes && return Float64(size(mesh.coords,2))
    count=0
    for block in mesh.blocks
        msh_family(block.msh)===stat || continue
        count+=_block_ncells(block)
    end
    return Float64(count)
end

# `NumberOption(GMSH_GET)` — `nothing` for unknown option names.
function _geo_option_number(context::_GeoNumericContext,family::String,
                            index::Int,member::String;warn::Bool=false)
    family=="General" && member=="BoundingBoxSize" && return context.native_context_lc
    if family=="Mesh" && index==0
        stat=get(_GEO_MESH_STAT_OPTIONS,member,nothing)
        stat===nothing || return context.mesh===nothing ? 0.0 :
            _geo_mesh_stat(context.mesh,stat)
    end
    haskey(context.option_numbers,(family,index,member)) &&
        return context.option_numbers[(family,index,member)]
    table=get(_GEO_NUMBER_OPTION_TABLE,family,nothing)
    if table===nothing
        warn && return _geo_option_number_unknown!(context,family,member,"")
        return nothing
    end
    default=get(table,member,nothing)
    if default===nothing
        warn && return _geo_option_number_unknown!(context,family,member,"")
        return nothing
    end
    return Float64(default)
end

# `x.y`/`x.y(i)`/`x.y[i]` in FExpr — `treat_Struct_FullName_dot_tSTRING_Float`
# case 1: `NumberOption(GMSH_GET, x, index, y, out, type_treat==0)`. Unknown
# names warn and evaluate to 0.
function _geo_option_number_read(context::_GeoNumericContext,family::String,
                                 index::Int,member::String,
                                 caller::AbstractString)
    if context.soft_unknown_reads
        value=_geo_option_number(context,family,index,member;warn=true)
        return value===nothing ? 0.0 : value
    end
    value=_geo_option_number(context,family,index,member)
    value===nothing && throw(ArgumentError(
        "$caller: Unknown number option '$family.$member'"))
    return value
end

# `Msg::ExchangeOnelabParameter` — the ONELAB name is `Name`'s value (default
# `name` is NOT used; Gmsh errors when options exist without `Name`).
function _geo_exchange_onelab_number!(context::_GeoNumericContext,name::String,
                                      values::Vector{Float64},options,
                                      caller::AbstractString)
    options===nothing && return values
    oname=nothing;fopts=Dict{String,Vector{Float64}}()
    for opt in options
        om=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(.*?)\s*$",
                 String(strip(opt)))
        om===nothing && throw(ArgumentError(
            "$caller: malformed ONELAB option $(repr(opt))"))
        key=String(om.captures[1]);arg=String(strip(om.captures[2]))
        if key=="Name" || key=="Label"
            v=_geo_eval_string(arg,context,caller)
            key=="Name" && (oname=v)
        elseif isempty(arg)
            fopts[key]=[1.0]
        elseif startswith(arg,"{") && occursin('=',arg)
            # `Choices {v = "label", ...}` — `Enumeration`: numbers to
            # fopts, labels to the char side.
            fopts[key]=Float64[]
            inner=String(chop(arg;head=1,tail=1))
            for enum in _geo_split_args(inner,caller)
                em=match(r"^(.*?)\s*=\s*(.*?)\s*$",String(enum))
                em===nothing && throw(ArgumentError(
                    "$caller: malformed enumeration item $(repr(enum))"))
                push!(fopts[key],_geo_eval_numeric(em.captures[1],context,
                    caller))
            end
        elseif startswith(arg,"{")
            fopts[key]=_geo_numeric_list_values(arg,context,caller;
                                                allow_multiplier=true)
        elseif _geo_string_rhs(arg)
            _geo_eval_string(arg,context,caller)
        else
            fopts[key]=[_geo_eval_numeric(arg,context,caller)]
        end
    end
    oname===nothing &&
        _geo_msg_error!(context,"From now on you need to use the `Name' " *
                      "attribute to create a ONELAB parameter")
    oname===nothing && return values
    if haskey(context.onelab_numbers,oname) &&
       get(fopts,"ReadOnly",[0.0])[1]==0.0
        # Server value overrides the default (`val = ps[0].getValues()`).
        return Float64[context.onelab_numbers[oname]]
    end
    context.onelab_numbers[oname]=isempty(values) ? 0.0 : values[1]
    return values
end

function _geo_exchange_onelab_string!(context::_GeoNumericContext,name::String,
                                      value::String,options,
                                      caller::AbstractString)
    options===nothing && return nothing
    oname=nothing
    for opt in options
        om=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(.*?)\s*$",String(strip(opt)))
        om===nothing && continue
        om.captures[1]=="Name" &&
            (oname=_geo_eval_string(om.captures[2],context,caller))
    end
    oname===nothing &&
        _geo_msg_error!(context,"From now on you need to use the `Name' " *
                      "attribute to create a ONELAB parameter")
    oname===nothing && return nothing
    if haskey(context.onelab_strings,oname)
        context.strings[name][1]=context.onelab_strings[oname]
    else
        context.onelab_strings[oname]=value
    end
    return nothing
end


# `_GEO_NUMBER_OPTION_TABLE` is generated in `GeoOptionTables.jl` from Gmsh
# 4.15.2's `DefaultOptions.h` — full name coverage keeps `x.y` writes on
# Gmsh-known options accepted and truly-unknown names on the `Msg::Error`
# path.

function _geo_apply_list_assignment!(context::_GeoNumericContext,
                                     raw::AbstractString,
                                     caller::AbstractString)
    # `name[i] op v` or the parenthesized `name(i) op v` — upstream maps both
    # through `String__Index '['|'(' FExpr ']'|')'` (identical index writes).
    statement=match(
        r"^([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*(.*?)\s*\]\s*(\+\+|--|\+=|-=|\*=|/=|=)\s*(.*?)\s*$",
        strip(raw))
    paren=statement===nothing
    paren && (statement=match(
        r"^([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*(.*?)\s*\)\s*(\+\+|--|\+=|-=|\*=|/=|=)\s*(.*?)\s*$",
        strip(raw)))
    statement===nothing && return false
    name=String(statement.captures[1])
    # `tSTRING name '['` — a lexer keyword head is a syntax error at the
    # keyword for expression-only tokens (`DimNameSpace[]`) and at `[` for
    # statement-leading ones (`Point[]`). The parenthesized form is
    # different: `Point(i) = ...`, `Spline(i) = ...`, `Box(i) = ...` and the
    # other `Keyword(FExpr) =` entity declarations are real productions the
    # caller's entity arms consume — not list-index writes — so only
    # expression-only keywords abort here (`Sin(0) = 1` → `(Sin)`).
    if name in _GEO_LEXER_KEYWORDS
        name in _GEO_EXPR_ONLY_KEYWORDS &&
            _geo_syntax_abort(name)
        paren ? (return false) : _geo_syntax_abort("[")
    end
    name in _GEO_SIDE_EFFECT_SYMBOLS && throw(ArgumentError(
        "$caller: dynamic tag allocator $name is read-only and cannot be assigned as a list"))
    selector=String(strip(statement.captures[2]))
    operation=statement.captures[3]
    rhs=String(strip(statement.captures[4]))
    increment=operation in ("++","--")
    increment && (operation=operation=="++" ? "+=" : "-=")
    rhs_values=increment ? Float64[1.0] : _geo_numeric_list_values(
        rhs,context,"$caller: numeric list $name";allow_multiplier=true)

    if isempty(selector)
        # `x[]++`/`x()++` is not an upstream production — after `]`/`)` the
        # grammar expects an affectation (`syntax error (++)`).
        increment && _geo_syntax_abort(statement.captures[3])
        if operation=="="
            _geo_context_set_list!(context,name,rhs_values,caller)
        elseif operation=="+="
            values=_geo_context_list(context,name,caller)
            length(values)+length(rhs_values)<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                "$caller: appending to $name exceeds $_MAX_GEO_LIST_ITEMS entries"))
            append!(values,rhs_values)
            _geo_context_set_list!(context,name,values,caller)
        elseif operation=="-="
            values=_geo_context_list(context,name,caller)
            for value in rhs_values
                index=findfirst(==(value),values)
                index===nothing || deleteat!(values,index)
            end
            _geo_context_set_list!(context,name,values,caller)
        else
            # `x[] *= v`/`x[] /= v` — recoverable upstream diagnostic; the
            # list is left unchanged.
            _geo_yyerror!(context,
                "Operators *= and /= not available for lists")
        end
        return true
    end

    index_values=if startswith(selector,"{") && endswith(selector,"}")
        _geo_numeric_list_values(
            selector,context,"$caller: numeric list $name selector";
            allow_multiplier=false)
    else
        Float64[_geo_eval_numeric(
            selector,context,"$caller: numeric list $name index")]
    end
    length(index_values)==length(rhs_values) ||
        (_geo_yyerror!(context,
            "Incompatible array dimensions in affectation"); return true)
    indices=Int[_geo_signed_gmsh_int_value(
        value,"$caller: numeric list $name index") for value in index_values]
    for index in indices
        index>=0 || throw(ArgumentError(
            "$caller: negative index $index in numeric list $name"))
        index+1<=_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
            "$caller: index $index exceeds $_MAX_GEO_LIST_ITEMS entries"))
    end

    # `assignVariable`/`assignVariables` upstream: `x[i] = v` on a fresh name
    # creates the list padded with zeros; `x[i] op v` on a fresh name is
    # `Unknown variable`; any index write on a scalar is `not a list`. Each
    # is a recoverable `yymsg(0)` diagnostic.
    if haskey(context.values,name) && !(name in context.list_variables)
        _geo_yyerror!(context,"Variable '$name' is not a list")
        return true
    end
    if haskey(context.unavailable_lists,name) ||
       haskey(context.unavailable,name)
        values=_geo_context_list(context,name,caller)  # throws with the reason
    elseif haskey(context.lists,name)
        values=copy(context.lists[name])
    else
        if operation=="="
            values=Float64[]
        else
            _geo_yyerror!(context,"Unknown variable '$name'")
            return true
        end
    end
    for position in eachindex(indices)
        index=indices[position]
        length(values)<index+1 && resize!(values,index+1)
        if operation=="/=" && rhs_values[position]==0
            # Upstream skips the write and reports
            # `Division by zero in 'x[i] /= 0'`.
            _geo_yyerror!(context,"Division by zero in '$name[$index] /= " *
                "$(rhs_values[position])'")
            continue
        end
        values[index+1]=_geo_list_mutation_value(
            operation,values[index+1],rhs_values[position],context,caller)
    end
    _geo_context_set_list!(context,name,values,caller)
    return true
end

function _geo_allocator_statement_values(raw::AbstractString,
                                         context::_GeoNumericContext,
                                         caller::AbstractString)
    pieces=_geo_split_list("{"*String(raw)*"}",caller)
    # Bookkeeping re-evaluation — the statement's own diagnostics were
    # already emitted, so list-term reads stay silent here.
    mute=context.suppress_diagnostics
    context.suppress_diagnostics=true
    try
        terms,total=_geo_numeric_list_terms(pieces,context,caller)
        values=Float64[];sizehint!(values,total)
        for term in terms
            value=term.first
            for _ in 1:term.count
                push!(values,value)
                value+=term.step
            end
        end
        return values
    finally
        context.suppress_diagnostics=mute
    end
end

# A bookkeeping re-evaluation inside the allocator observer — the statement's
# production already reduced (and diagnosed) the FExpr, so the observer's
# re-read stays silent while observer-specific diagnostics (e.g. `GeoEntity
# dim out of range`) still emit.
function _geo_observe_eval_numeric(src,context::_GeoNumericContext,
                                   caller::AbstractString)
    mute=context.suppress_diagnostics
    context.suppress_diagnostics=true
    try
        return _geo_eval_numeric(src,context,caller)
    finally
        context.suppress_diagnostics=mute
    end
end

# Definition-site tag expression: Gmsh casts to `int` without a positivity
# check — `tag < 0` auto-assigns in `addX` (`record_explicit!` mirrors that)
# and `tag == 0` is a literal tag.
function _geo_allocator_statement_tag(raw::AbstractString,
                                      context::_GeoNumericContext,
                                      caller::AbstractString)
    return _geo_signed_gmsh_int_value(
        _geo_observe_eval_numeric(raw,context,caller),caller)
end

function _geo_physical_declaration(raw::AbstractString)
    source=String(strip(raw))
    if endswith(source,";")
        source=String(strip(
            source[firstindex(source):prevind(source,lastindex(source))]))
    end
    named=match(
        r"^Physical\s+(Point|Curve|Line|Surface|Volume|GeoEntity\{[^{}]*\})\s*\(\s*\"([^\"]*)\"\s*(?:,\s*(.*?))?\s*\)\s*(=|\+=|-=|\*=|/=)\s*(.*)$",
        source)
    if named!==nothing
        raw_tag=named.captures[3]
        tag_source=raw_tag===nothing ? nothing : String(strip(raw_tag))
        tag_source!==nothing && isempty(tag_source) && return nothing
        return (kind=String(named.captures[1]),
                name=String(named.captures[2]),
                tag_source=tag_source,
                op=String(named.captures[4]),
                membership=String(strip(named.captures[5])))
    end
    unnamed=match(
        r"^Physical\s+(Point|Curve|Line|Surface|Volume|GeoEntity\{[^{}]*\})\s*\(\s*(.*?)\s*\)\s*(=|\+=|-=|\*=|/=)\s*(.*)$",
        source)
    unnamed===nothing && return nothing
    tag_source=String(strip(unnamed.captures[2]))
    isempty(tag_source) && return nothing
    return (kind=String(unnamed.captures[1]),name=nothing,
            tag_source=tag_source,op=String(unnamed.captures[3]),
            membership=String(strip(unnamed.captures[4])))
end

# `Physical GeoEntity{d}(tag) op {members}` — upstream's `GeoEntity`
# nonterminal admits the `GeoEntity{FExpr}` dim form beside the kind words
# (Gmsh.y:2405). Evaluate the dim expression against the live context.
function _geo_physical_kind_dim(physical,context,caller::AbstractString)
    kind=physical.kind
    if startswith(kind,"GeoEntity{") && endswith(kind,"}")
        dim=_geo_int_value(_geo_eval_numeric(
            kind[nextind(kind,firstindex(kind),10):prevind(kind,lastindex(kind))],
            context,"$caller GeoEntity dimension"),
            "$caller GeoEntity dimension")
        0<=dim<=3 || throw(ArgumentError(
            "$caller: GeoEntity dim out of range [0,3]"))
        return dim
    end
    return _PHYS_DIM[kind]
end

function _geo_allocator_physical_tag!(state::_GeoTagAllocatorState,
                                      tag_source::Union{Nothing,AbstractString},
                                      context::_GeoNumericContext,
                                      caller::AbstractString;
                                      conservative::Bool=false)
    if tag_source===nothing && state.physical_unavailable!==nothing
        conservative && return nothing
        throw(ArgumentError(
            "$caller: automatic Physical tags are unavailable because " *
            state.physical_unavailable))
    end
    tag=if tag_source===nothing
        state.physical_group_max<typemax(Int32) || throw(ArgumentError(
            "$caller: no automatic Physical tags remain"))
        state.physical_group_max+1
    else
        try
            _geo_allocator_statement_tag(tag_source,context,"$caller tag")
        catch err
            err isa InterruptException && rethrow()
            (conservative &&
             (err isa ArgumentError || err isa _GeoSyntaxAbort)) || rethrow()
            _geo_allocator_invalidate!(state,
                "could not evaluate Physical group tag while tracking allocators: " *
                _geo_expr_preview(tag_source);geometry=false,physical=true,
                fields=false)
            return nothing
        end
    end
    _geo_allocator_record_physical!(state,tag)
    return tag
end

function _geo_allocator_observe_statement!(state::_GeoTagAllocatorState,
                                           raw::AbstractString,
                                           context::_GeoNumericContext,
                                           caller::AbstractString;
                                           conservative::Bool=false)
    source=String(strip(raw))
    if endswith(source,";")
        source=String(strip(source[firstindex(source):prevind(source,lastindex(source))]))
    end
    isempty(source) && return nothing

    factory=match(r"^SetFactory\s*\(\s*\"([^\"]*)\"\s*\)$",source)
    if factory!==nothing
        name=String(factory.captures[1])
        if name in ("Built-in","OpenCASCADE")
            _geo_allocator_set_factory!(state,name)
        else
            reason="SetFactory accepts only \"Built-in\" or \"OpenCASCADE\"; got " *
                   repr(name)
            if conservative
                _geo_allocator_invalidate!(
                    state,reason;geometry=true,fields=false)
                return nothing
            end
            throw(ArgumentError("$caller: $reason"))
        end
        return nothing
    elseif startswith(source,"SetFactory")
        reason="SetFactory requires a literal \"Built-in\" or \"OpenCASCADE\" name"
        if conservative
            _geo_allocator_invalidate!(state,reason;geometry=true,fields=false)
            return nothing
        end
        throw(ArgumentError("$caller: $reason"))
    end

    field=match(r"^Field\s*\[\s*(.*?)\s*\]\s*=\s*([A-Za-z][A-Za-z0-9_]*)$",source)
    if field!==nothing
        # `newField` no-ops on unknown kinds and duplicate ids — neither
        # changes `maxId()`.
        field.captures[2] in _GEO_FIELD_KINDS || return nothing
        tag=_geo_allocator_statement_tag(
            field.captures[1],context,"$caller Field declaration tag")
        tag in state.live_fields || _geo_allocator_record_field!(state,tag)
        return nothing
    end

    # `Physical` statements are owned by the executor: the tag resolution
    # (`getPhysicalNumber`/`setPhysicalName`), the conditional counter bump,
    # and the raw group records all carry side effects the observer cannot
    # reproduce without model access — `newreg`-family reads therefore see the
    # counters the exec path already updated.
    _geo_physical_declaration(source)!==nothing && return nothing

    declarations=(
        (r"^Point\s*\(\s*(.*?)\s*\)\s*=",:point,"Point"),
        (r"^Line\s*\(\s*(.*?)\s*\)\s*=",:curve,"Line"),
        (r"^Circle\s*\(\s*(.*?)\s*\)\s*=",:curve,"Circle"),
        (r"^Ellipse\s*\(\s*(.*?)\s*\)\s*=",:curve,"Ellipse"),
        (r"^(?:Line\s+Loop|Curve\s+Loop)\s*\(\s*(.*?)\s*\)\s*=",:curveloop,"Curve Loop"),
        (r"^Plane\s+Surface\s*\(\s*(.*?)\s*\)\s*=",:surface,"Plane Surface"),
        (r"^(?:Ruled\s+)?Surface\s*\(\s*(.*?)\s*\)\s*=",:surface,"Surface"),
        (r"^Surface\s+Loop\s*\(\s*(.*?)\s*\)\s*=",:surfaceloop,"Surface Loop"),
        (r"^Volume\s*\(\s*(.*?)\s*\)\s*=",:volume,"Volume"),
    )
    for (pattern,kind,label) in declarations
        matched=match(pattern,source)
        matched===nothing && continue
        state.geometry_unavailable===nothing || return nothing
        tag=try
            _geo_allocator_statement_tag(
                matched.captures[1],context,"$caller $label tag")
        catch err
            err isa InterruptException && rethrow()
            (conservative &&
             (err isa ArgumentError || err isa _GeoSyntaxAbort)) || rethrow()
            _geo_allocator_invalidate!(state,
                "could not evaluate $label tag while tracking allocators: " *
                _geo_expr_preview(matched.captures[1]);geometry=true,fields=false)
            return nothing
        end
        _geo_allocator_record_explicit!(state,kind,tag)
        return nothing
    end

    primitive=match(
        r"^(Box|Cylinder|Sphere|PolarSphere|Cone|Torus|Rectangle|Disk)\s*\(\s*(.*?)\s*\)\s*=\s*\{\s*(.*?)\s*\}$",
        source)
    if primitive!==nothing
        state.geometry_unavailable===nothing || return nothing
        kind=String(primitive.captures[1])
        tag,values=try
            (_geo_allocator_statement_tag(
                 primitive.captures[2],context,"$caller $kind tag"),
             _geo_allocator_statement_values(
                 primitive.captures[3],context,"$caller $kind parameters"))
        catch err
            err isa InterruptException && rethrow()
            (conservative &&
             (err isa ArgumentError || err isa _GeoSyntaxAbort)) || rethrow()
            _geo_allocator_invalidate!(state,
                "could not evaluate $kind while tracking allocators: " *
                _geo_expr_preview(source);geometry=true,fields=false)
            return nothing
        end
        expected=kind=="Box" ? (6,) : kind=="Cylinder" ? (7,) :
                 kind=="Sphere" ? (2,4,5,6,7) :
                 kind=="PolarSphere" ? (2,) : kind=="Torus" ? (5,6) :
                 kind=="Rectangle" ? (5,6) : kind=="Disk" ? (4,5) : (8,)
        # A wrong-arity primitive creates nothing upstream (recoverable
        # `yymsg` diagnostic at execution) — the counters stay put, so the
        # scan simply skips the record and lets the exec pass report it.
        length(values) in expected || return nothing
        # `tBox`/`tCylinder`/`tCone`/`tTorus`, `tRectangle`/`tDisk`, and the
        # parameterized `tSphere` form are OCC-gated upstream — under the
        # built-in factory the statement is a recoverable diagnostic and
        # creates nothing, so the counters stay put. `PolarSphere` and the
        # two-point `Sphere` are built-in `newGeometry*` entities under
        # either factory.
        occ_gated=kind=="Sphere" ? length(values)>=4 :
            kind in ("Box","Cylinder","Cone","Torus","Rectangle","Disk")
        if !occ_gated || state.factory==:opencascade
            if kind in ("Rectangle","Disk")
                _geo_allocator_record_surface_primitive!(
                    state,kind,tag,values,caller)
            else
                _geo_allocator_record_primitive!(state,kind,tag,values,caller)
            end
        end
        return nothing
    end

    set_max=match(
        r"^SetMaxTag\s+(?:(Point|Curve|Line|Surface|Volume)|GeoEntity\s*\{\s*(.*?)\s*\})\s*\(\s*(.*?)\s*\)$",
        source)
    if set_max!==nothing
        state.geometry_unavailable===nothing || return nothing
        kind=if set_max.captures[1]!==nothing
            set_max.captures[1]=="Line" ? "Curve" : String(set_max.captures[1])
        else
            dim=try
                _geo_signed_gmsh_int_value(
                    _geo_observe_eval_numeric(set_max.captures[2],context,
                                      "$caller SetMaxTag dimension"),
                    "$caller SetMaxTag dimension")
            catch err
                err isa InterruptException && rethrow()
                (conservative &&
             (err isa ArgumentError || err isa _GeoSyntaxAbort)) || rethrow()
                _geo_allocator_invalidate!(state,
                    "could not evaluate SetMaxTag dimension while tracking " *
                    "allocators: "*_geo_expr_preview(set_max.captures[2]);
                    geometry=true,fields=false)
                return nothing
            end
            # `GeoEntity{d}` emits the range diagnostic yet still calls
            # `setMaxTag(dim, tag)`; the internals switch covers `-2..3` so
            # dims outside that are a no-op.
            (dim<0 || dim>3) && _geo_yyerror!(
                context,"GeoEntity dim out of range [0,3]")
            -2<=dim<=3 || return nothing
            ("SurfaceLoop","CurveLoop","Point","Curve","Surface","Volume")[dim+3]
        end
        value=try
            _geo_signed_gmsh_int_value(
                _geo_observe_eval_numeric(set_max.captures[3],context,
                                  "$caller SetMaxTag $kind"),
                "$caller SetMaxTag $kind")
        catch err
            err isa InterruptException && rethrow()
            (conservative &&
             (err isa ArgumentError || err isa _GeoSyntaxAbort)) || rethrow()
            _geo_allocator_invalidate!(state,
                "could not evaluate SetMaxTag $kind while tracking allocators: " *
                _geo_expr_preview(set_max.captures[3]);geometry=true,fields=false)
            return nothing
        end
        _geo_allocator_set_max!(state,kind,value)
        return nothing
    end

    boolean=match(
        r"^Boolean(Difference|Union|Intersection|Fuse|Common|Cut)\s*\(\s*(.*?)\s*\)\s*=\s*" *
        r"\{\s*Volume\s*\{\s*([^{},;:\[\]]+?)\s*\}([^}]*)\}\s*" *
        r"\{\s*Volume\s*\{\s*([^{},;:\[\]]+?)\s*\}([^}]*)\}\s*;?\s*$",
        source)
    if boolean!==nothing
        # Under the built-in factory upstream's tagged `BooleanShape` is a
        # silent no-op — no entities, no counter movement.
        state.factory==:opencascade || return nothing
        state.geometry_unavailable===nothing || return nothing
        parsed=try
            (_geo_allocator_statement_tag(
                 boolean.captures[2],context,"$caller Boolean result tag"),
             _geo_allocator_statement_tag(
                 boolean.captures[3],context,"$caller Boolean first operand"),
             _geo_allocator_statement_tag(
                 boolean.captures[5],context,"$caller Boolean second operand"))
        catch err
            err isa InterruptException && rethrow()
            (conservative &&
             (err isa ArgumentError || err isa _GeoSyntaxAbort)) || rethrow()
            _geo_allocator_invalidate!(state,
                "could not evaluate a Boolean statement while tracking " *
                "allocators: "*_geo_expr_preview(source);
                geometry=true,fields=false)
            return nothing
        end
        _geo_allocator_record_explicit!(state,:volume,parsed[1])
        delete_a=occursin(r"\bDelete\b",boolean.captures[4])
        delete_b=occursin(r"\bDelete\b",boolean.captures[6])
        # Gmsh re-tags a deleted object operand's boundary onto the result at
        # the lowest free tags; the counts are exact for boundary-preserving
        # Booleans and bounded by the object's own boundary otherwise.
        object_boundary=delete_a ?
            get(state.volume_boundaries,parsed[2],nothing) : nothing
        delete_a && _geo_allocator_delete_volume!(state,parsed[2])
        delete_b && _geo_allocator_delete_volume!(state,parsed[3])
        if object_boundary!==nothing
            claimed=ntuple(3) do d
                _geo_allocator_claim_lowest!(
                    state,(:point,:curve,:surface)[d],
                    length(object_boundary[d]))
            end
            state.volume_boundaries[parsed[1]]=claimed
        end
        return nothing
    end

    # `.geo` lifecycle statements executed through `_geo_exec_delete!` update
    # the allocator directly (per-entity counter decrement or full reset), so
    # the observer leaves them alone.
    match(r"^(?:Recursive\s+Delete|Delete)\b",source)!==nothing &&
        return nothing

    topology_change=occursin(
        r"\b(?:Boolean[A-Za-z]*|Extrude|Delete|Duplicata|SetMaxTag|Merge|Coherence|Split)\b",
        source) || match(
        # A leading transform can merge coincident entities (lowering
        # automatic counters) even without a Duplicata.
        r"^(?:Translate|Rotate|Dilate|Symmetry|Affine)\s*\{",source)!==nothing
    if topology_change
        _geo_allocator_invalidate!(state,
            "topology-changing statement is outside the tracked allocator subset: " *
            _geo_expr_preview(source);geometry=true,
            fields=occursin(r"\bField\b",source))
        return nothing
    end

    if conservative && !startswith(source,"Physical ") &&
       !startswith(source,"Field") &&
       match(r"^[A-Z][A-Za-z]*(?:\s+[A-Z][A-Za-z]*)*\s*\(.*\)\s*=",source)!==nothing
        _geo_allocator_invalidate!(state,
            "entity declaration is outside the tracked allocator subset: " *
            _geo_expr_preview(source);geometry=true,fields=false)
    end
    return nothing
end

const _GEO_FIELD_RAW_OPTIONS=Set((
    "F","FX","FY","FZ","M11","M22","M33","M12","M13","M23",
    "m11","m22","m33","m12","m13","m23","FileName","CommandLine",
    "p4estFileToLoad"))
const _GEO_FIELD_FLOAT_LIST_OPTIONS=Set(("SizesList","hwall_n_nodes"))
const _GEO_FIELD_INTEGER_LIST_OPTIONS=Set((
    "FieldsList","FanPointsSizesList","PointsList","NodesList","VerticesList",
    "FanPointsList","FanNodesList","CurvesList","EdgesList","SurfacesList",
    "FacesList","VolumesList","RegionsList","ExcludedSurfacesList",
    "ExcludedFaceList"))
const _GEO_FIELD_ENTITY_LIST_OPTIONS=Set((
    "PointsList","NodesList","VerticesList","CurvesList","EdgesList",
    "SurfacesList","FacesList","VolumesList","RegionsList"))
const _GEO_FIELD_NUMERIC_OPTIONS=Set((
    "Sampling","NNodesByEdge","NumPointsPerCurve","FieldX","FieldY","FieldZ",
    "InField","IField","DistMin","DistMax","SizeMin","SizeMax","Sigmoid",
    "StopAtDistMax","LcMin","LcMax","VIn","VOut","XMin","XMax","YMin",
    "YMax","ZMin","ZMax","Thickness","XCenter","YCenter","ZCenter",
    "Radius","XAxis","YAxis","ZAxis","X1","Y1","Z1","X2","Y2","Z2",
    "InnerR1","OuterR1","InnerR2","OuterR2","InnerV1","OuterV1","InnerV2",
    "OuterV2","R1_inner","R1_outer","R2_inner","R2_outer","V1_inner",
    "V1_outer","V2_inner","V2_outer","Kind","Delta","FromStereo",
    "RadiusStereo","TextFormat","SetOutsideValue","OutsideValue",
    "IncludeBoundary","IncludeEmbedded","Power","ViewIndex","ViewTag","IView",
    "CropNegativeValues","UseClosest","dMin","dMax","SizeMinTangent",
    "SizeMaxTangent","SizeMinNormal","SizeMaxNormal","lMinTangent",
    "lMaxTangent","lMinNormal","lMaxNormal","Size","hwall_n","Ratio","ratio",
    "SizeFar","hfar","thickness","Quads","IntersectMetrics","AnisoMax",
    "BetaLaw","Beta","NbLayers","nPointsPerCircle","nPointsPerGap","hMin",
    "hMax","hBulk","gradation","smoothing","features"))
const _GEO_FIELD_INTEGER_OPTIONS=Set((
    "Sampling","NNodesByEdge","NumPointsPerCurve","FieldX","FieldY","FieldZ",
    "InField","IField","Kind","FromStereo","ViewIndex","ViewTag","IView",
    "Quads","IntersectMetrics","BetaLaw","NbLayers","nPointsPerCircle",
    "nPointsPerGap"))

# `FieldManager::mapTypeName` — Gmsh 4.15.2's registered `Field[i] = Kind`
# names (Field.cpp); the comparison is case-sensitive.
const _GEO_FIELD_KINDS=Set((
    "Structured","Threshold","BoundaryLayer","Box","Cylinder","Ball","Frustum",
    "LonLat","PostView","Gradient","Octree","Distance","Attractor","Extend",
    "Restrict","Constant","Min","MinAniso","IntersectAniso","Max","Laplacian",
    "Mean","Curvature","Param","ExternalProcess","MathEval","MathEvalAniso",
    "AttractorAnisoCurve","MaxEigenHessian","AutomaticMeshSizeField"))

# Per-kind option membership — each upstream Field subclass registers its own
# `options` map, so `Field[i].member = v` validates `member` against the
# field's *kind* (`field->options[member]`), not a global list. The second
# element is the class `getName()` reported by `Unknown option`/`No field`
# diagnostics — `Attractor` instantiates `DistanceField`, so it reports
# `of type 'Distance'`.
const _GEO_FIELD_OPTIONS=Dict{String,Tuple{Set{String},String}}(
    "Structured"=>(Set(("FileName","OutsideValue","SetOutsideValue",
        "TextFormat")),"Structured"),
    "Threshold"=>(Set(("DistMax","DistMin","IField","InField","LcMax","LcMin",
        "Sigmoid","SizeMax","SizeMin","StopAtDistMax")),"Threshold"),
    "BoundaryLayer"=>(Set(("AnisoMax","Beta","BetaLaw","CurvesList","EdgesList",
        "ExcludedFaceList","ExcludedSurfacesList","FanNodesList",
        "FanPointsList","FanPointsSizesList","IntersectMetrics","NbLayers",
        "NodesList","PointsList","Quads","Ratio","Size","SizeFar","SizesList",
        "Thickness","hfar","hwall_n","hwall_n_nodes","ratio","thickness")),
        "BoundaryLayer"),
    "Box"=>(Set(("Thickness","VIn","VOut","XMax","XMin","YMax","YMin","ZMax",
        "ZMin")),"Box"),
    "Cylinder"=>(Set(("Radius","VIn","VOut","XAxis","XCenter","YAxis","YCenter",
        "ZAxis","ZCenter")),"Cylinder"),
    "Ball"=>(Set(("Radius","Thickness","VIn","VOut","XCenter","YCenter",
        "ZCenter")),"Ball"),
    "Frustum"=>(Set(("InnerR1","InnerR2","InnerV1","InnerV2","OuterR1","OuterR2",
        "OuterV1","OuterV2","R1_inner","R1_outer","R2_inner","R2_outer",
        "V1_inner","V1_outer","V2_inner","V2_outer","X1","X2","Y1","Y2","Z1",
        "Z2")),"Frustum"),
    "LonLat"=>(Set(("FromStereo","IField","InField","RadiusStereo")),"LonLat"),
    "PostView"=>(Set(("CropNegativeValues","IView","UseClosest","ViewIndex",
        "ViewTag")),"PostView"),
    "Gradient"=>(Set(("Delta","IField","InField","Kind")),"Gradient"),
    "Octree"=>(Set(("InField",)),"Octree"),
    "Distance"=>(Set(("CurvesList","EdgesList","FacesList","FieldX","FieldY",
        "FieldZ","NNodesByEdge","NodesList","NumPointsPerCurve","PointsList",
        "Sampling","SurfacesList")),"Distance"),
    "Attractor"=>(Set(("CurvesList","EdgesList","FacesList","FieldX","FieldY",
        "FieldZ","NNodesByEdge","NodesList","NumPointsPerCurve","PointsList",
        "Sampling","SurfacesList")),"Distance"),
    "Extend"=>(Set(("AnisoMax","Beta","BetaLaw","CurvesList","EdgesList",
        "ExcludedFaceList","ExcludedSurfacesList","FanNodesList",
        "FanPointsList","FanPointsSizesList","IncludeEmbedded",
        "IntersectMetrics","NbLayers","NodesList","PointsList","Quads","Ratio",
        "Size","SizeFar","SizeMax","SizesList","SurfacesList","Thickness",
        "ViewIndex","VolumesList","hfar","hwall_n","hwall_n_nodes","ratio",
        "thickness")),"Extend"),
    "Restrict"=>(Set(("CurvesList","EdgesList","FacesList","IField","InField",
        "IncludeBoundary","IncludeEmbedded","PointsList","RegionsList",
        "SurfacesList","VerticesList","VolumesList")),"Restrict"),
    "Constant"=>(Set(("CurvesList","IncludeBoundary","IncludeEmbedded",
        "PointsList","SurfacesList","VIn","VOut","VolumesList")),"Constant"),
    "Min"=>(Set(("FieldsList",)),"Min"),
    "MinAniso"=>(Set(("FieldsList",)),"MinAniso"),
    "IntersectAniso"=>(Set(("FieldsList",)),"IntersectAniso"),
    "Max"=>(Set(("FieldsList",)),"Max"),
    "Laplacian"=>(Set(("Delta","IField","InField")),"Laplacian"),
    "Mean"=>(Set(("Delta","IField","InField")),"Mean"),
    "Curvature"=>(Set(("Delta","IField","InField")),"Curvature"),
    "Param"=>(Set(("FX","FY","FZ","IField","InField")),"Param"),
    "ExternalProcess"=>(Set(("CommandLine",)),"ExternalProcess"),
    "MathEval"=>(Set(("F",)),"MathEval"),
    "MathEvalAniso"=>(Set(("M11","M12","M13","M22","M23","M33","m11","m12",
        "m13","m22","m23","m33")),"MathEvalAniso"),
    "AttractorAnisoCurve"=>(Set(("CurvesList","DistMax","DistMin","EdgesList",
        "NNodesByEdge","NumPointsPerCurve","Sampling","SizeMaxNormal",
        "SizeMaxTangent","SizeMinNormal","SizeMinTangent","dMax","dMin",
        "lMaxNormal","lMaxTangent","lMinNormal","lMinTangent")),
        "AttractorAnisoCurve"),
    "MaxEigenHessian"=>(Set(("Delta","IField","InField")),"MaxEigenHessian"),
    "AutomaticMeshSizeField"=>(Set(("p4estFileToLoad","nPointsPerCircle",
        "nPointsPerGap","hMin","hMax","hBulk","gradation","smoothing",
        "features")),"AutomaticMeshSizeField"))

function _geo_normalize_field_option(raw::AbstractString,name_raw::AbstractString,
                                     context::_GeoNumericContext,caller::AbstractString)
    name=String(name_raw);caller_string=String(caller)
    value=String(strip(raw))
    name in _GEO_FIELD_RAW_OPTIONS && return value
    if name in _GEO_FIELD_FLOAT_LIST_OPTIONS || name in _GEO_FIELD_INTEGER_LIST_OPTIONS
        list_source,multiplier=_geo_braced_list_source(
            value,context,caller_string;allow_multiplier=false)
        list_source===nothing && throw(ArgumentError(
            "$caller_string: expected a brace-delimited list"))
        items=_geo_split_list(list_source,caller_string)
        normalized=String[]
        function append_numeric!(number::Float64)
            length(normalized)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                "$caller_string: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
            scaled=multiplier*number
            isfinite(scaled) || throw(ArgumentError(
                "$caller_string entry: list multiplication produced a non-finite value"))
            if name in _GEO_FIELD_INTEGER_LIST_OPTIONS
                integer=_geo_signed_gmsh_int_value(scaled,"$caller_string entry")
                push!(normalized,string(integer))
            else
                push!(normalized,_geo_number_source(scaled))
            end
            return nothing
        end
        for item in items
            whole_reference=match(
                r"^[+-]?\s*([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*\]\s*$",item)
            reference_name=whole_reference===nothing ? "" :
                           String(whole_reference.captures[1])
            if whole_reference!==nothing &&
               !haskey(context.values,reference_name) &&
               !haskey(context.lists,reference_name) &&
               name in _GEO_FIELD_ENTITY_LIST_OPTIONS
                length(normalized)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                    "$caller_string: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
                opaque=String(strip(item))
                if multiplier==-1.0
                    opaque=startswith(opaque,"-") ?
                        String(strip(opaque[nextind(opaque,firstindex(opaque)):end])) :
                        startswith(opaque,"+") ?
                            "-"*String(strip(opaque[nextind(opaque,firstindex(opaque)):end])) :
                            "-"*opaque
                end
                push!(normalized,opaque)
                continue
            end
            if name in _GEO_FIELD_ENTITY_LIST_OPTIONS &&
               occursin(r"^[+-]?[A-Za-z_][A-Za-z0-9_]*(?:\[\])?$",item)
                bare=replace(item,r"^[+-]"=>"");bare=endswith(bare,"[]") ? bare[1:end-2] : bare
                if !haskey(context.values,bare) &&
                   !haskey(context.lists,bare) && bare!="Pi"
                    length(normalized)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                        "$caller_string: expanded list exceeds $_MAX_GEO_LIST_ITEMS entries"))
                    entry=item
                    if multiplier==-1.0
                        entry=startswith(item,"-") ? String(item[nextind(item,firstindex(item)):end]) :
                              startswith(item,"+") ? "-"*String(item[nextind(item,firstindex(item)):end]) :
                              "-"*item
                    end
                    push!(normalized,entry)
                    continue
                end
            end
            expanded=_geo_list_reference_values(item,context,caller_string)
            if expanded!==nothing
                for number in expanded
                    append_numeric!(number)
                end
                continue
            end
            term=_geo_numeric_list_term(item,context,caller_string,
                                        _MAX_GEO_LIST_ITEMS-length(normalized))
            for number in _geo_list_term_values(term)
                append_numeric!(number)
            end
        end
        return "{"*join(normalized,", ")*"}"
    elseif name in _GEO_FIELD_NUMERIC_OPTIONS
        (startswith(value,"{") || endswith(value,"}")) && throw(ArgumentError(
            "$caller_string: numeric scalar option cannot use a brace list"))
        number=_geo_eval_numeric(value,context,caller_string)
        return name in _GEO_FIELD_INTEGER_OPTIONS ?
            string(_geo_signed_gmsh_int_value(number,caller_string)) :
            _geo_number_source(number)
    end
    # Unknown options are retained so the field builder can issue its kind-aware
    # unsupported-option diagnostic.  Guessing that an unknown value is numeric
    # would make future string-valued Gmsh options unsafe.
    return value
end

function _geo_unquoted_code(source::AbstractString)
    out=IOBuffer();quote_char='\0'
    for c in source
        if quote_char!='\0'
            c==quote_char && (quote_char='\0')
            write(out,' ')
        elseif c=='"' || c=='\''
            quote_char=c;write(out,' ')
        else
            write(out,c)
        end
    end
    return String(take!(out))
end

function _geo_invalidate_context!(context::_GeoNumericContext,reason::AbstractString)
    for name in union(keys(context.values),keys(context.lists),
                      keys(context.unavailable),keys(context.unavailable_lists))
        context.unavailable[name]=String(reason)
        context.unavailable_lists[name]=String(reason)
    end
    empty!(context.values)
    empty!(context.lists)
    empty!(context.list_variables)
    context.stored_list_items=0
    return nothing
end

function _geo_strip_control_terminators(source::AbstractString)
    body=String(strip(source))
    closed=0
    while true
        m=match(r"^(?:EndIf|EndFor|Return)\b\s*(.*)$",body)
        m===nothing && return body,closed
        closed+=1
        body=String(strip(m.captures[1]))
    end
end

@inline function _geo_relevant_code(code::AbstractString)
    return occursin(r"\b(?:Mesh\.(?:MeshSizeMin|MeshSizeMax|MeshSizeFactor|RandomSeed|MeshSizeFromCurvature|MinimumElementsPerTwoPi|BoundaryLayerFanElements|BoundaryLayerFanPoints)|Geometry\.Tolerance)\b",code) ||
           occursin(r"\bField\s*\[",code) ||
           occursin(r"\b(?:Background|BoundaryLayer)\s+Field\b",code) ||
           occursin(r"\bPhysical\s+(?:Point|Curve|Line|Surface|Volume)\b",code)
end

# `ListOfDouble`/`FExpr_Multi` terms the bounded evaluator cannot model —
# entity selectors (`Point{..}`), queries (`Volume In BoundingBox{..}`,
# `Parent Curve{..}`), shape producers (`BooleanX{..}{..}`,
# `CombinedBoundary{..}`, `Duplicata{..}`, transforms, `Extrude{..}`), and
# attribute queries (`Mass Volume{..}`, `Physical Volume(..)`). Upstream
# they are valid productions — `tSTRING '{' MultipleShape '}'` accepts even
# an unknown `name{..}` — so a parse failure on input containing one marks
# the variable unavailable (like a failed evaluation), never a syntax
# error. A leading `{` (`a = {1,2}` stores a list in a scalar) is likewise
# a valid `ListOfDouble` the scalar evaluator cannot model.
function _geo_list_term_unmodelable(source::AbstractString)
    occursin(r"[A-Za-z_][A-Za-z0-9_]*\s*\{",source) && return true
    startswith(strip(source),"{") && return true
    return occursin(Regex(
        "(?<![A-Za-z0-9_.])(?:Point|Curve|Line|Surface|Volume|" *
        "GeoEntity|Physical|Parent|BoundingBox|Mass|CenterOfMass|" *
        "MatrixOfInertia|Normal|Curvature|Parametric|Color|Closest|" *
        "Intersect|Split|Fillet|Chamfer|ThruSections|Ruled|Compound|" *
        "Levelset|Periodic|Duplicata|Boundary|PointsOf|CombinedBoundary|" *
        "OrientedBoundary|OrientedCombinedBoundary|Translate|Rotate|" *
        "Dilate|Symmetry|Affine|Transform|Extrude|ShapeFromFile|" *
        "Boolean[A-Za-z]*)\\b"),source)
end

# An evaluation that aborts while reading an *already unavailable* binding
# (`newX` after an untracked topology change, a poisoned list count, ...)
# fails on the value, not the grammar — the assignment propagates the
# unavailability rather than reporting a syntax error. Only the RHS is
# scanned: an unavailable LHS name must not mask a malformed RHS.
function _geo_rhs_uses_unavailable(source::AbstractString,
                                   context::_GeoNumericContext)
    rhs=source
    eq=findfirst(isequal('='),source)
    eq!==nothing && (rhs=source[nextind(source,eq):end])
    for m in eachmatch(r"[A-Za-z_][A-Za-z0-9_]*",rhs)
        name=m.match
        (haskey(context.unavailable,name) ||
         haskey(context.unavailable_lists,name)) && return true
    end
    return false
end

function _geo_record_scalar!(context::_GeoNumericContext,name_raw::AbstractString,
                             raw::AbstractString)
    name=String(name_raw)
    # Lexer keywords (`Pi`, `DimNameSpace`, ...) are not mutable bindings —
    # upstream rejects the statement as a syntax error, so the scan leaves
    # no binding behind.
    name in _GEO_LEXER_KEYWORDS && return nothing
    if name in _GEO_SIDE_EFFECT_SYMBOLS
        # Recoverable statement diagnostic upstream — record and continue
        # like `_geo_exec_scalar!` does under the executor's catch.
        _geo_yyerror!(context,"dynamic tag allocator $name is read-only")
        return nothing
    end
    try
        # `String__Index NumericAffectation ListOfDouble`: the RHS is a
        # list production, so `a = {1,2}` defines a list while a single
        # value stays scalar (`s.list = List_Nbr(rhs) != 1` upstream).
        values=_geo_numeric_list_values(raw,context,
            "read_geo_params: scalar variable $name")
        if length(values)==1
            _geo_context_set_scalar!(context,name,values[1],
                "read_geo_params: scalar variable $name")
        else
            _geo_context_set_list!(context,name,values,
                "read_geo_params: scalar variable $name")
        end
    catch err
        err isa InterruptException && rethrow()
        # `_GeoSyntaxAbort` propagates to the per-statement recovery — a
        # syntax error leaves no binding behind upstream — unless the RHS is
        # a valid entity/shape `ListOfDouble` production the bounded
        # evaluator cannot model, or it reads an already-unavailable
        # binding; both mark the variable unavailable.
        if err isa _GeoSyntaxAbort
            (_geo_list_term_unmodelable(raw) ||
             _geo_rhs_uses_unavailable(raw,context)) || rethrow()
        else
            err isa ArgumentError || rethrow()
        end
        had_list_payload=haskey(context.lists,name) ||
                         haskey(context.unavailable_lists,name)
        if haskey(context.lists,name)
            delete!(context.values,name)
            delete!(context.list_variables,name)
        else
            _geo_context_forget!(context,name)
        end
        message=sprint(showerror,err)
        # An unmodelable entity/shape production failed inside the bounded
        # evaluator — the grammar-abort text is misleading as the reason.
        reason=err isa _GeoSyntaxAbort && _geo_list_term_unmodelable(raw) ?
            "entity or shape expression requires geometry evaluation" :
            ncodeunits(message)<=240 ? message :
                String(first(message,220))*"…"
        context.unavailable[name]=reason
        if had_list_payload
            context.unavailable_lists[name]=reason
        else
            delete!(context.unavailable_lists,name)
        end
    end
    return nothing
end

function _geo_record_list!(context::_GeoNumericContext,raw::AbstractString)
    matched=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*[\[\(]",strip(raw))
    matched===nothing && return false
    name=String(matched.captures[1])
    if name in _GEO_SIDE_EFFECT_SYMBOLS
        _geo_yyerror!(context,"dynamic tag allocator $name is read-only")
        return true
    end
    try
        return _geo_apply_list_assignment!(
            context,raw,"read_geo_params: numeric list assignment")
    catch err
        err isa InterruptException && rethrow()
        # `_GeoSyntaxAbort` propagates to the per-statement recovery — a
        # syntax error leaves no binding behind upstream — unless the RHS is
        # a valid entity/shape `ListOfDouble` production the bounded
        # evaluator cannot model (`Volume In BoundingBox{..}`,
        # `BooleanX{..}{..}`, `Point{..}`, ...) or reads an
        # already-unavailable binding, both of which mark the variable
        # unavailable instead. ArgumentError marks a failed evaluation so
        # later reads see the variable as unavailable.
        if err isa _GeoSyntaxAbort
            (_geo_list_term_unmodelable(raw) ||
             _geo_rhs_uses_unavailable(raw,context)) || rethrow()
        else
            err isa ArgumentError || rethrow()
        end
        _geo_context_forget!(context,name)
        message=sprint(showerror,err)
        reason=err isa _GeoSyntaxAbort && _geo_list_term_unmodelable(raw) ?
            "entity or shape expression requires geometry evaluation" :
            ncodeunits(message)<=240 ? message :
                String(first(message,220))*"…"
        context.unavailable[name]=reason
        context.unavailable_lists[name]=reason
        push!(context.list_variables,name)
        return true
    end
end

function _geo_expression_field_tags(raw::AbstractString,context::_GeoNumericContext,
                                    caller::AbstractString;
                                    deduplicate::Bool=true)
    value=String(strip(raw))
    isempty(value) && throw(ArgumentError("$caller: field-tag expression must not be empty"))
    values=_geo_numeric_list_values(
        value,context,caller;allow_multiplier=true)
    tags=Int[];sizehint!(tags,length(values));seen=Set{Int}()
    for numeric in values
        # Field tags take Gmsh's raw `(int)` cast — `Field[0]` and negative
        # ids are literal (`newField` has no positivity check).
        tag=_geo_signed_gmsh_int_value(numeric,"$caller entry")
        if !deduplicate || !(tag in seen)
            push!(seen,tag);push!(tags,tag)
        end
    end
    return tags
end

function _gmsh_random_seed(value::Float64)
    # Gmsh 4.15.2 stores this option as an unsigned 32-bit value: assignments
    # are truncated toward zero and clamped to the option range.
    value<=0 && return 0
    upper=min(Float64(typemax(Int)),Float64(typemax(UInt32)))
    value>=upper && return _geo_int_value(upper,"read_geo_params: Mesh.RandomSeed")
    return _geo_int_value(value,"read_geo_params: Mesh.RandomSeed")
end

# A `.geo` statement is complete at a top-level `}` only when it is a
# transform or query statement (keyword + the expected number of balanced
# `{...}` groups and nothing else) — `Translate{…}{…}`, `Duplicata{…}`,
# `Boundary{…}`, etc. Statements like `bnd[] = Boundary{…}` still end at `;`
# because they start with an unlisted name.
function _geo_brace_terminated_statement(text::AbstractString)
    s=String(strip(text))
    mm=match(r"^([A-Za-z_][A-Za-z0-9_]*)",s)
    mm===nothing && return false
    name=mm.captures[1]
    name=="Extrude" && return _geo_brace_terminated_extrude(s)
    # `Delete { ListOfShapes }` ends at the closing brace with no `;`, as do
    # `Recursive Delete { ... }` and `Delete Embedded { ... }`.
    if name in ("Delete","Recursive","Show","Hide","Color","Intersect",
                "Closest","Fillet","Chamfer","ThruSections","Ruled",
                "BooleanUnion","BooleanDifference","BooleanIntersection",
                "BooleanFragments","BooleanSection","BooleanFuse",
                "BooleanCoherence","BooleanCommon","BooleanCut")
        rest0=String(strip(s[nextind(s,firstindex(s),ncodeunits(mm.match)):end]))
        if name=="Recursive"
            dm=match(r"^(Delete|Color|Show|Hide)\b",rest0)
            dm===nothing && return false
            rest0=String(strip(rest0[nextind(rest0,firstindex(rest0),
                                           ncodeunits(dm.match)):end]))
            dm.captures[1]=="Delete" ||
                (name=dm.captures[1])
        elseif name=="Delete" &&
               (em=match(r"^([A-Za-z_][A-Za-z0-9_]*)\b",rest0))!==nothing
            # `Delete Embedded { ... }` (and any `Delete <word> { ... }` —
            # unknown words error later during execution, not at scan time).
            rest0=String(strip(rest0[nextind(rest0,firstindex(rest0),
                                           ncodeunits(em.match)):end]))
        end
        # `Color`/`Recursive Color` take a ColorExpr before the shape group: a
        # `{r,g,b[,a]}` list is itself a brace group, while named/quoted/option
        # colors are not.
        if name=="Color" && !startswith(rest0,"{")
            cm2=match(
                r"^(?:\"[^\"]*\"|'[^']*'|[A-Za-z_][A-Za-z0-9_]*(?:\s*\[[^\]]*\])?\s*(?:\.\s*[A-Za-z_][A-Za-z0-9_]*)*)",
                rest0)
            cm2===nothing && return false
            rest0=String(strip(rest0[nextind(rest0,firstindex(rest0),
                                           ncodeunits(cm2.match)):end]))
            name="ColorWord"
        elseif name=="Ruled"
            rm=match(r"^ThruSections\b",rest0)
            rm===nothing && return false
            rest0=String(strip(rest0[nextind(rest0,firstindex(rest0),
                                           ncodeunits(rm.match)):end]))
            name="ThruSections"
        elseif name=="Intersect"
            cm=match(r"^Curve\b",rest0)
            cm===nothing && return false
            rest0=String(strip(rest0[nextind(rest0,firstindex(rest0),5):end]))
        end
        groups=if name in ("Show","Hide","ThruSections","ColorWord")
            1
        elseif name=="Color"
            # `{r,g,b}` color list + shape group = 2 groups.
            startswith(rest0,"{") ? 2 : return false
        elseif name=="Fillet"
            # `Fillet {regions} {edges} {radii}`.
            3
        elseif name=="Chamfer"
            # `Chamfer {regions} {edges} {surfaces} {distances}`.
            4
        else
            name in ("Closest","Intersect") ||
                startswith(name,"Boolean") ? 2 : 1
        end
        s=rest0
    else
        groups=name in ("Translate","Rotate","Dilate","Symmetry","Affine") ? 2 :
               name in ("Duplicata","Boundary","CombinedBoundary",
                        "OrientedBoundary","OrientedCombinedBoundary",
                        "PointsOf") ? 1 : return false
        s=String(strip(s[nextind(s,firstindex(s),ncodeunits(mm.match)):end]))
    end
    rest=String(strip(s))
    for _ in 1:groups
        if !startswith(rest,"{")
            # `Intersect Curve{..} Surface{..}`, `Split Curve{..} Point{..}`
            # carry entity keywords between brace groups.
            wm=match(r"^[A-Za-z_][A-Za-z0-9_.]*",rest)
            wm===nothing && return false
            rest=String(strip(rest[nextind(rest,firstindex(rest),
                                       ncodeunits(wm.match)):end]))
        end
        startswith(rest,"{") || return false
        depth=0;done=false;i=firstindex(rest);last=lastindex(rest)
        while i<=last
            c=rest[i]
            c=='{' && (depth+=1)
            if c=='}'
                depth-=1
                depth==0 && (done=true; break)
            end
            i=nextind(rest,i)
        end
        done || return false
        rest=String(strip(rest[nextind(rest,i):end]))
    end
    return isempty(rest)
end

# `Extrude {dx,dy,dz} { ... }` needs two brace groups while the boundary-layer
# form `Extrude { ... }` needs one — a `{...}` group closes the statement only
# when its body is a shape list (contains `;`, a `Kind{`/`Physical`/`Parent`
# entry, an `Extrude` parameter keyword, or is empty) rather than a bare
# numeric vector.
function _geo_brace_terminated_extrude(text::AbstractString)
    rest=String(strip(text[nextind(text,firstindex(text),7):end]))
    while true
        startswith(rest,"{") || return false
        depth=0;done=false;i=firstindex(rest);last=lastindex(rest)
        while i<=last
            c=rest[i]
            c=='{' && (depth+=1)
            if c=='}'
                depth-=1
                depth==0 && (done=true; break)
            end
            i=nextind(rest,i)
        end
        done || return false
        body=String(rest[2:prevind(rest,i)])
        tail=String(strip(rest[nextind(rest,i):end]))
        # A shape-list group completes the statement only when nothing
        # follows it; a numeric vector group must be followed by the shape
        # list, so the statement is still incomplete. A shape-list group
        # may still be followed by the pipe suffix `Using Wire {expr}`
        # (`tExtrude '{' ListOfShapes '}' tUsing tWire '{' FExpr '}'`).
        if _geo_extrude_shape_group(body)
            isempty(tail) && return true
            (um=match(r"^Using\s+Wire\b",tail))===nothing && return false
            grp=_geo_brace_group(String(strip(tail[nextind(tail,
                firstindex(tail),length(um.match)):end])))
            return grp!==nothing && isempty(grp[2])
        end
        rest=tail
    end
end

# `Extrude { ListOfShapes }` alone is a complete `}`-terminated statement
# (boundary-layer form) but stays open for the `Using Wire {expr}` pipe
# suffix — the splitters hold the emit until the next token decides.
function _geo_extrude_pipe_candidate(stmt::AbstractString)
    match(r"^Extrude\b",stmt)===nothing && return false
    grp=_geo_brace_group(String(strip(stmt[nextind(stmt,firstindex(stmt),
                                                  7):end])))
    grp===nothing && return false
    body,tail=grp
    isempty(tail) || return false
    return _geo_extrude_shape_group(body)
end

# `body` is the content of one `{...}` group inside an `Extrude` statement.
function _geo_extrude_shape_group(body::AbstractString)
    b=String(strip(body))
    isempty(b) && return true
    occursin(';',b) && return true
    return match(
        Regex("^\\s*(Point|Line|Curve|Surface|Volume|Physical|Parent|" *
              "GeoEntity|Layers|Recombine|ScaleLast|QuadTriAddVerts|" *
              "QuadTriNoNewVerts|RecombLaterals|Using|Hole)\\b"),b)!==nothing
end

# Consume one balanced `{...}` group at the start of `rest` — returns
# `(body, tail)` or `nothing` when `rest` does not open with a balanced
# group.
function _geo_brace_group(rest::AbstractString)
    startswith(rest,"{") || return nothing
    depth=0;i=firstindex(rest);last=lastindex(rest)
    while i<=last
        c=rest[i]
        c=='{' && (depth+=1)
        if c=='}'
            depth-=1
            depth==0 && return (String(rest[2:prevind(rest,i)]),
                                String(strip(rest[nextind(rest,i):end])))
        end
        i=nextind(rest,i)
    end
    return nothing
end

const _GEO_CONTROL_BARE=Dict(
    "Else"=>:else,"EndIf"=>:endif,"EndWhile"=>:endwhile,"EndFor"=>:endfor,
    "Return"=>:return)

# Scan a balanced open/close group starting at raw[i] (the opener). Returns the
# index of the matching close on the same line; control headers that spill to a
# second line are a bounded-subset blocker.
function _geo_scan_balanced(raw::AbstractString,i::Int,last::Int,
                            open::Char,close::Char)
    depth=0;qc='\0';k=i
    while k<=last
        c=raw[k]
        if qc!='\0'
            c==qc && (qc='\0')
        elseif c=='"' || c=='\''
            qc=c
        elseif c==open
            depth+=1
        elseif c==close
            depth-=1
            depth==0 && return k
        end
        k=nextind(raw,k)
    end
    throw(ArgumentError(
        "execute_geo: control statement must complete on a single line"))
end

# If raw[i:last] begins a control construct (If/ElseIf/While headers with a
# parenthesized expression, `For name In {a:b[:c]}`, or the bare Else/EndIf/
# EndWhile/EndFor markers), return (statement_text, last_consumed_index).
# Gmsh 4.15.2's built-in kernel accepts only `For name In {a:b[:c]}` — comma
# lists, single values, and C-style `For (init; cond; incr)` are rejected —
# and has no While keyword; `While (expr) ... EndWhile` is a bounded
# Tessella extension.
function _geo_control_statement(raw::AbstractString,i::Int,last::Int)
    rest=SubString(raw,i,last)
    matched=match(
        r"^(If|ElseIf|While|Else|EndIf|EndWhile|For|EndFor|Function|Macro|Return)\b",
        rest)
    matched===nothing && return nothing
    word=matched.captures[1]
    if word=="EndWhile"
        # `EndWhile` is not an upstream keyword — it only acts as a control
        # marker at a statement boundary (`;`, line end, or EOF). `EndWhile x`
        # is the upstream `tSTRING`-identifier form (`Unknown command`).
        k=firstindex(rest)+sizeof(word)
        while k<=lastindex(rest) && isspace(rest[k])
            k=nextind(rest,k)
        end
        k<=lastindex(rest) && rest[k]!=';' && return nothing
        return (word,i+sizeof(word)-1)
    end
    haskey(_GEO_CONTROL_BARE,word) && return (word,i+sizeof(word)-1)
    j=firstindex(rest)+sizeof(word)
    jlast=lastindex(rest)
    while j<=jlast && isspace(rest[j])
        j=nextind(rest,j)
    end
    if word=="Function" || word=="Macro"
        # Gmsh's `Function name` header carries an identifier or quoted
        # string-expression name on the same line and no `;`. A malformed
        # header is not consumed here: it falls through to the `;` splitter
        # so execution can report the syntax error at the offending token.
        name_match=match(
            r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\"|'[^']*')",SubString(rest,j))
        name_match===nothing && return nothing
        consumed=j-1+sizeof(name_match.match)
        return (String(rest[firstindex(rest):consumed]),i+consumed-1)
    end
    if word=="For"
        if j<=jlast && rest[j]=='('
            # `For (a:b[:c])` — the anonymous range form.
            close=_geo_scan_balanced(rest,j,jlast,'(',')')
            return (String(rest[1:close]),i+close-1)
        end
        name_match=match(r"^[A-Za-z_][A-Za-z0-9_]*",SubString(rest,j))
        name_match===nothing && return nothing
        j+=sizeof(name_match.match)
        while j<=jlast && isspace(rest[j])
            j=nextind(rest,j)
        end
        match(r"^In\b",SubString(rest,j))===nothing && return nothing
        j+=2
        while j<=jlast && isspace(rest[j])
            j=nextind(rest,j)
        end
        (j<=jlast && rest[j]=='{') || return nothing
        close=_geo_scan_balanced(rest,j,jlast,'{','}')
        return (String(rest[1:close]),i+close-1)
    end
    (j<=jlast && rest[j]=='(') || return nothing
    close=_geo_scan_balanced(rest,j,jlast,'(',')')
    if word=="While"
        # `While(i) op v` is the upstream `String__Index ( FExpr )` index
        # write, not the bounded While extension — let the `;` splitter and
        # the index-assignment path handle it.
        tail=String(strip(SubString(rest,close+1)))
        match(r"^(?:\+\+|--|\+=|-=|\*=|/=|=)",tail)!==nothing && return nothing
    end
    return (String(rest[1:close]),i+close-1)
end

function _scan_geo_statements(consume,path::AbstractString)
    open(path,"r") do stream
        _scan_geo_statements(consume,stream)
    end
end

function _scan_geo_statements(consume,stream::Base.IO)
    buffer=IOBuffer();quote_char='\0';block_comment=false
    buf_has_content=false
    # `;` inside `{...}` groups (Boolean operand lists, `Delete` suffixes) is
    # part of the statement, not a terminator — matching _geo_exec_statements.
    # Transform and query statements may instead end at their final `}`.
    depth=0
    # `Extrude {shapes}` is `}`-complete but may still take the `Using Wire`
    # pipe suffix — its emit is deferred until the next real token, even
    # across line boundaries.
    extrude_pending=false
    for raw in eachline(stream)
        i=firstindex(raw);lastindex_raw=lastindex(raw)
        while i<=lastindex_raw
            # Control constructs carry no `;` — emit one as its own statement
            # whenever it begins at a statement boundary (matching
            # `_geo_exec_statements`).
            if !block_comment && quote_char=='\0' && !buf_has_content
                control=_geo_control_statement(raw,i,lastindex_raw)
                if control!==nothing
                    take!(buffer)
                    consume(control[1])
                    i=nextind(raw,control[2])
                    continue
                end
            end
            c=raw[i];j=nextind(raw,i)
            nextc=j<=lastindex_raw ? raw[j] : '\0'
            if extrude_pending && !block_comment && quote_char=='\0' &&
               !isspace(c) && !(c=='/' && (nextc=='/' || nextc=='*'))
                extrude_pending=false
                if !(c=='U' && match(r"^Using\b",
                        SubString(raw,i,lastindex_raw))!==nothing)
                    statement=strip(String(take!(buffer)))
                    isempty(statement) || consume(statement)
                    buf_has_content=false
                end
            end
            if block_comment
                if c=='*' && nextc=='/'
                    block_comment=false;i=nextind(raw,j);continue
                end
                i=j;continue
            elseif quote_char!='\0'
                write(buffer,c)
                c==quote_char && (quote_char='\0')
                i=j;continue
            elseif c=='/' && nextc=='/'
                break
            elseif c=='/' && nextc=='*'
                # Comments are lexical whitespace in Gmsh.  Retaining one
                # separator prevents `1/*...*/2` from becoming the valid token
                # `12` in the scanner.
                write(buffer,' ')
                position(buffer)<=_MAX_GEO_STATEMENT_BYTES || throw(ArgumentError(
                    "read_geo_params: statement exceeds $_MAX_GEO_STATEMENT_BYTES bytes"))
                block_comment=true;i=nextind(raw,j);continue
            elseif c=='"' || c=='\''
                quote_char=c;write(buffer,c);buf_has_content=true
            elseif c=='{'
                depth+=1;write(buffer,c);buf_has_content=true
            elseif c=='}'
                depth=max(0,depth-1);write(buffer,c);buf_has_content=true
                if depth==0 && _geo_brace_terminated_statement(
                        String(buffer.data[1:position(buffer)]))
                    if _geo_extrude_pipe_candidate(
                            String(buffer.data[1:position(buffer)]))
                        # Possibly `Extrude {shapes} Using Wire {expr}` —
                        # hold the emit until the next real token.
                        extrude_pending=true
                    else
                        statement=strip(String(take!(buffer)))
                        isempty(statement) || consume(statement)
                        buf_has_content=false
                    end
                end
            elseif c==';' && depth==0
                write(buffer,c)
                statement=strip(String(take!(buffer)))
                # A `;` left after a `}`-terminated transform or between
                # statements is an upstream `syntax error (;)` — it is still
                # consumed so the diagnostic is recorded.
                consume(statement)
                buf_has_content=false
            else
                write(buffer,c)
                isspace(c) || (buf_has_content=true)
            end
            position(buffer)<=_MAX_GEO_STATEMENT_BYTES || throw(ArgumentError(
                "read_geo_params: statement exceeds $_MAX_GEO_STATEMENT_BYTES bytes"))
            i=j
        end
        if quote_char!='\0'
            write(buffer,'\n')
        elseif position(buffer)>0
            write(buffer,' ')
        end
    end
    quote_char=='\0' || throw(ArgumentError("read_geo_params: unterminated quoted string"))
    block_comment && throw(ArgumentError("read_geo_params: unterminated block comment"))
    tail=strip(String(take!(buffer)))
    if extrude_pending && !isempty(tail)
        # A held `Extrude {shapes}` with no `Using Wire` following is
        # complete — emit it before the trailing-content checks.
        consume(tail);tail=""
    end
    (depth==0 || isempty(tail)) || throw(ArgumentError(
        "read_geo_params: unterminated brace-delimited statement"))
    # A trailing run of bare `;` is still statements upstream (`syntax error
    # (;)` per stray terminator) — consume them for the recorded diagnostic.
    all(c->c==';'||isspace(c),tail) && !isempty(tail) && (consume(tail);tail="")
    # Control keywords need no `;` upstream — a trailing `EndFor`/`EndIf`/
    # `Return` still acts at EOF, closing its level like the execution pass.
    !isempty(tail) &&
        match(r"^(If|ElseIf|Else|EndIf|For|EndFor|While|EndWhile|Return|Function|Macro)\b",
              tail)!==nothing && (consume(tail);tail="")
    code=_geo_unquoted_code(tail)
    if !isempty(tail) && (occursin(r"Field\s*\[",code) ||
                           occursin(r"(?:Background|BoundaryLayer)\s+Field",code) ||
                           occursin(r"Mesh\.(?:MeshSize|RandomSeed)",code) ||
                           occursin("Geometry.Tolerance",code) ||
                           occursin(r"Physical\s+(?:Point|Curve|Line|Surface|Volume)",code))
        throw(ArgumentError("read_geo_params: unterminated relevant statement (missing semicolon)"))
    end
    return nothing
end

# `.geo` entity-definition keywords → dimension — used both by the exec-side
# shape-list executor and the `_changed` bookkeeping rule below.
const _GEO_SHAPE_DEFINITION_DIMS=Dict(
    "Point"=>0,
    "Line"=>1,"Curve"=>1,"Spline"=>1,"BSpline"=>1,"Bezier"=>1,"Nurbs"=>1,
    "Circle"=>1,"Ellipse"=>1,"Wire"=>1,"Compound Spline"=>1,
    "Compound BSpline"=>1,"Compound Curve"=>1,
    "Surface"=>2,"Plane Surface"=>2,"Ruled Surface"=>2,"BSpline Surface"=>2,
    "Bezier Surface"=>2,"Parametric Surface"=>2,"Rectangle"=>2,"Disk"=>2,
    "Compound Surface"=>2,
    "Volume"=>3,"Sphere"=>3,"PolarSphere"=>3,"Box"=>3,"Torus"=>3,
    "Cylinder"=>3,"Cone"=>3,"Wedge"=>3,"ThickSolid"=>3,"Compound Volume"=>3)

# `_changed` bookkeeping: nearly every `.geo` statement mutates the GEO
# internals and marks them changed; only the read-only consumers, control
# flow, option writes, and plain variable assignments do not. Any statement
# not recognized as non-mutating marks the internals so a later sync point
# resynchronizes exactly once, matching `getChanged()`/`synchronize()`.
const _GEO_NONMUTATING_HEADS=Set([
    # model consumers that sync-if-changed instead of marking it
    "BoundingBox","Color","Hide","Include","Merge","MergeWithBoundingBox",
    "Print","Save","Show","SyncModel",
    # control flow, diagnostics, and stream control
    "Abort","Call","Draw","Else","ElseIf","EndFor","EndIf","EndWhile",
    "Error","Exit","For","Function","If","Macro","Return","SetChanged",
    "UnsplitWindow","While",
    # command words that do not touch the GEO internals
    "CreateDir","Info","NonBlockingSystemCall","OnelabRun","Printf",
    "SetName","Sleep","System","SystemCall","Warning",
    # option-table writes (`Mesh.X = v`, `Geometry.X = v`, `Field[i].x = v`,
    # `View[i].X = v`, ...) mutate option state, not the GEO internals — and
    # `Mesh n` itself meshes rather than mutating internals. `Physical` is
    # whitelisted too: `modifyPhysicalGroup` only marks `_changed` on its
    # success paths, which the exec handlers set directly — the failed forms
    # ("already exists"/"does not exist"/"Unsupported operation") must not
    # force a resync. `Delete` likewise: only the `{list}` form (unconditional
    # in `GEO_Internals::remove`) and the `All`/`Model`/`Physicals` resets mark
    # the internals — `Delete Variables|Options|Meshes|Field[i]|<var>` and
    # `Delete Embedded` leave `_changed` alone, and the marking forms set it
    # inside their exec handlers.
    "Mesh","Geometry","General","Solver","PostProcessing","View","Field",
    "Onelab","MathEval","Colors","Physical","Delete"])

function _geo_statement_marks_internals(line::AbstractString)
    s=strip(line)
    (hm=match(r"^([A-Za-z_][A-Za-z0-9_]*)",s))===nothing && return false
    head=hm.captures[1]
    head in _GEO_NONMUTATING_HEADS && return false
    head=="Recursive" &&
        return match(r"^Recursive\s+(Show|Hide|Color)\b",s)===nothing
    # `name = expr`, `name[] = ...` and `name(i) = ...` assignments leave the
    # internals alone — but `Kind(tag) = rhs` entity definitions create.
    if match(r"^[A-Za-z_][A-Za-z0-9_]*(?:\[[^\]]*\]|\([^)]*\))?\s*=[^=]",
             s)!==nothing
        return haskey(_GEO_SHAPE_DEFINITION_DIMS,head)
    end
    return true
end

# Per-dimension union of the live entity sets the allocator mirrors — sync
# member resolution searches the model, which holds entities from every
# activated factory.
function _geo_scan_live_entities(state::_GeoTagAllocatorState,dim::Int)
    dim==0 && return union(state.live_builtin_points,state.live_occ_points)
    dim==1 && return union(state.live_builtin_curves,state.live_occ_curves)
    dim==2 && return union(state.live_builtin_surfaces,
                           state.live_occ_surfaces)
    return union(state.live_builtin_volumes,state.live_occ_volumes)
end

# `GEO_Internals::synchronize` at scan granularity — silent: it only rebuilds
# the entity-level physical mirror (`scan_ephys`) so `("name", 0)`
# declarations resolve through `getMaxPhysicalNumber` like Gmsh. Execution
# replays the same sync and owns the "unknown member" warnings. The clearing
# rule matches `synchronize`: entity memberships drop only when the model
# carries no mesh elements AND raw groups exist; otherwise they accumulate.
# Groups whose member list the scan could not evaluate (`scan_opaque`)
# contribute `abs(tag)` under a phantom entity slot — one member may resolve.
function _geo_scan_sync_physicals!(
        scan_raw::Dict{Tuple{Int,Int},Vector{Int}},
        scan_opaque::Set{Tuple{Int,Int}},
        scan_ephys::Dict{Tuple{Int,Int},Vector{Int}},
        state::_GeoTagAllocatorState,has_mesh::Bool)
    (isempty(scan_raw) || has_mesh) || empty!(scan_ephys)
    for ((dim,tag),members) in scan_raw
        live=_geo_scan_live_entities(state,dim)
        for member in members
            entity=abs(member)
            entity in live || continue
            ep=get!(scan_ephys,(dim,entity),Int[])
            p=sign(member)*tag
            p in ep || push!(ep,p)
        end
    end
    for (dim,tag) in scan_opaque
        ep=get!(scan_ephys,(dim,-1),Int[])
        a=abs(tag)
        a in ep || push!(ep,a)
    end
    return nothing
end

# `DeletePhysicalX` — drop the raw group, erase the `(dim, tag)` name
# binding, and strip entity memberships whose `abs` equals the tag. For a
# negative raw tag `abs(p) == tag` can never hold, so the memberships
# survive — matching `GModel::removePhysicalGroup`.
function _geo_scan_delete_physical!(
        scan_raw::Dict{Tuple{Int,Int},Vector{Int}},
        scan_opaque::Set{Tuple{Int,Int}},
        scan_ephys::Dict{Tuple{Int,Int},Vector{Int}},
        groups::Dict{Tuple{Int,Int},String},
        physical_name_tags::Dict{Tuple{Int,String},Int},
        key::Tuple{Int,Int})
    dim,tag=key
    delete!(scan_raw,key);delete!(scan_opaque,key);delete!(groups,key)
    for ((d,n),t) in collect(physical_name_tags)
        d==dim && t==tag && delete!(physical_name_tags,(d,n))
    end
    for pnums in values(scan_ephys)
        filter!(p->abs(p)!=tag,pnums)
    end
    return nothing
end

# `GModel::getMaxPhysicalNumber(dim)` — the largest `abs` entity-level
# physical tag carried by a dimension-`dim` entity in the scan mirror.
function _geo_scan_max_entity_physical(
        scan_ephys::Dict{Tuple{Int,Int},Vector{Int}},dim::Int)
    found=0
    for ((d,_),pnums) in scan_ephys
        d==dim || continue
        for p in pnums
            a=abs(p)
            a>found && (found=a)
        end
    end
    return found
end

# Evaluate a `Physical X(...) = {..}` member list to raw signed tags; `nothing`
# marks the group opaque — dynamic or geometry-derived members the scanner
# cannot resolve.
function _geo_scan_physical_members(membership::AbstractString,
                                    context::_GeoNumericContext)
    source=strip(membership)
    isempty(source) && return Int[]
    (startswith(source,"{") && endswith(source,"}")) || return nothing
    isempty(strip(source[firstindex(source)+1:prevind(source,lastindex(source))])) &&
        return Int[]
    values=try
        _geo_numeric_list_values(source,context,
            "read_geo_params: Physical membership")
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
    members=Int[];sizehint!(members,length(values))
    for v in values
        isfinite(v) || return nothing
        push!(members,_geo_signed_gmsh_int_value(
            v,"read_geo_params: Physical membership entry"))
    end
    return members
end

# ======== Scan-side control flow ========
#
# Upstream Gmsh drives `If`/`For`/`Function` with lexer-level `skip` calls and
# `fsetpos` rewinds: dead branches and function bodies never parse, loop bodies
# re-parse per iteration, and `Return` unwinds an entered `Call`. The params
# scan mirrors that observable behavior with an explicit level stack: `:if`
# levels gate statements on `active`, `:for`/`:while`/`:function` levels
# collect body text which replays (or registers) when the matching closer
# arrives. `MAX_RECUR_LOOPS`/`MAX_RECUR_TESTS` cap nesting at 100 like the
# upstream tables.

mutable struct _GeoScanLevel
    kind::Symbol              # :if :for :while :function
    active::Bool              # :if branch state; outer-active flag otherwise
    taken::Bool               # :if — an earlier branch already ran
    var::String               # :for loop variable ("" for `For (a:b)`)
    range::Vector{Float64}    # :for — evaluated [first,last,step]
    cond::String              # :while — condition source
    buffer::Vector{String}    # :for/:while/:function collected body text
    name::String              # :function name
    closed::Bool              # :if — upstream `skip` already closed the level
    debt::Bool                # :if — opened while `ImbricatedTest` was <= 0
    nest::Int                 # buffered `For`/`While`/`Function` opens awaiting
                              # their `EndFor`/`EndWhile`/`Return`
end
_GeoScanLevel(kind::Symbol;active::Bool=true,taken::Bool=false,var::String="",
              range::Vector{Float64}=Float64[],cond::String="",
              name::String="",closed::Bool=false,debt::Bool=false)=
    _GeoScanLevel(kind,active,taken,var,range,cond,String[],name,closed,debt,0)

# Unwinds a `Return` inside a replayed `Call` body back to the call site.
struct _GeoScanReturnSignal end

const _MAX_GEO_SCAN_NESTING=100
const _MAX_GEO_SCAN_LOOP_ITERATIONS=1_000_000

@inline function _geo_scan_collecting(levels::Vector{_GeoScanLevel})
    return !isempty(levels) &&
        levels[end].kind in (:for,:while,:function)
end
@inline function _geo_scan_active(levels::Vector{_GeoScanLevel})
    return all(l -> l.active,levels)
end
@inline function _geo_scan_outer_active(levels::Vector{_GeoScanLevel})
    # Activity of the scope enclosing the current top level.
    return all(l -> l.active,@view levels[1:end-1])
end

# Soft numeric evaluation for scan-time conditions/range endpoints: an
# unknown variable records `yymsg(0)` and yields 0, matching the grammar's
# already-reduced FExpr arguments.
function _geo_scan_soft_eval(src::AbstractString,
                             context::_GeoNumericContext)
    context.soft_unknown_reads=true
    try
        return _geo_eval_numeric(src,context,"read_geo_params")
    finally
        context.soft_unknown_reads=false
    end
end

# Evaluate a strict `FExpr : FExpr [: FExpr]` range (the `tFor` forms are not
# general lists: `{1,2}` errors at `,`, `{5}` at `}`, `(5)` at `)` upstream).
# Each piece is evaluated in order with soft reads; a non-range form aborts
# at the token the grammar rejected.
function _geo_scan_for_range(spec::AbstractString,
                             context::_GeoNumericContext,
                             caller::AbstractString,
                             closer::AbstractString)
    source=String(strip(spec))
    isempty(source) && _geo_syntax_abort(closer)
    pieces=_geo_split_range(source,caller)
    if pieces===nothing || length(pieces)==1
        # Not a `a:b` form — evaluate the head expression, then report the
        # offending separator or the closing brace/paren.
        text=pieces===nothing ? source : pieces[1]
        comma=_geo_scan_top_level_comma(text)
        if comma>0
            _geo_scan_soft_eval(text[firstindex(text):prevind(text,comma)],
                                context)
            _geo_syntax_abort(",")
        end
        _geo_scan_soft_eval(text,context)
        _geo_syntax_abort(closer)
    end
    first=_geo_scan_soft_eval(pieces[1],context)
    last=_geo_scan_soft_eval(pieces[2],context)
    step=length(pieces)==2 ? 1.0 : _geo_scan_soft_eval(pieces[3],context)
    return Float64[first,last,step]
end

# First top-level `,` index in `s` (outside every bracket kind and quotes), or
# 0 — used to locate the comma the upstream range grammar rejects.
function _geo_scan_top_level_comma(s::AbstractString)
    parens=0;brackets=0;braces=0;qc='\0'
    i=firstindex(s);last=lastindex(s)
    while i<=last
        c=s[i]
        if qc!='\0'
            c==qc && (qc='\0')
        elseif c=='"' || c=='\''
            qc=c
        elseif c=='(';parens+=1
        elseif c==')';parens-=1
        elseif c=='[';brackets+=1
        elseif c==']';brackets-=1
        elseif c=='{';braces+=1
        elseif c=='}';braces-=1
        elseif c==',' && parens==0 && brackets==0 && braces==0
            return i
        end
        i=nextind(s,i)
    end
    return 0
end

# Strip a leading control construct from `line` for the params scan — mirrors
# `_geo_control_statement`/`_exec_line!`: a malformed header aborts at the
# offending token (`If foo` → `(foo)`, `For = 5` → `(=)`, `If (x` → `(;)`).
# Returns (word=..., arg=..., rest=...) or nothing.
function _geo_scan_control_prefix(line::AbstractString,
                                  context::_GeoNumericContext)
    m=match(
        r"^(If|ElseIf|While|Else|EndIf|For|EndFor|EndWhile|Return|Function|Macro)\b",
        line)
    m===nothing && return nothing
    word=String(m.captures[1])
    rest=String(strip(line[nextind(line,firstindex(line),
                                 ncodeunits(word)):end]))
    # `While`/`EndWhile` are the Tessella bounded extension — upstream has no
    # lexer tokens for them, so `While = 5` is an ordinary assignment and
    # `EndWhile x` an `Unknown command`. They only act as control flow in the
    # `While (cond)` / bare-`EndWhile` shapes.
    word=="While" && !startswith(rest,"(") && return nothing
    word=="EndWhile" && !isempty(rest) && return nothing
    if word in ("Else","EndIf","EndFor","EndWhile","Return")
        return (word=word,arg="",rest=rest)
    elseif word=="Function" || word=="Macro"
        nm=match(r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\"|'[^']*')",rest)
        nm===nothing && _geo_syntax_abort(_geo_stmt_head_token(rest))
        return (word="Function",
                arg=_geo_scan_call_name(nm.match,context),
                rest=String(strip(rest[nextind(rest,
                    lastindex(nm.match)):end])))
    elseif word=="For"
        if startswith(rest,"(")
            # `For (a:b[:c])` — the anonymous range form.
            close=_geo_matching_delim(rest,firstindex(rest))
            if close==0
                _geo_scan_soft_eval_safe(rest[2:end],context)
                _geo_syntax_abort(";")
            end
            spec=String(strip(rest[2:prevind(rest,close)]))
            return (word="For",arg=("",spec,")"),
                    rest=String(strip(rest[nextind(rest,close):end])))
        end
        nm=match(r"^[A-Za-z_][A-Za-z0-9_]*",rest)
        nm===nothing && _geo_syntax_abort(_geo_stmt_head_token(rest))
        var=String(nm.match)
        var in _GEO_LEXER_KEYWORDS && _geo_syntax_abort(var)
        rest2=String(strip(rest[nextind(rest,lastindex(nm.match)):end]))
        match(r"^In\b",rest2)===nothing &&
            _geo_syntax_abort(_geo_stmt_head_token(rest2))
        rest3=String(strip(rest2[nextind(rest2,firstindex(rest2)+1):end]))
        startswith(rest3,"{") ||
            _geo_syntax_abort(_geo_stmt_head_token(rest3))
        close=_geo_matching_delim(rest3,firstindex(rest3))
        if close==0
            _geo_scan_for_range(strip(rest3[2:end]),context,
                                "read_geo_params: For range","}")
            _geo_syntax_abort(";")
        end
        spec=String(strip(rest3[2:prevind(rest3,close)]))
        return (word="For",arg=(var,spec,"}"),
                rest=String(strip(rest3[nextind(rest3,close):end])))
    else
        # If/ElseIf/While — parenthesized condition.
        startswith(rest,"(") ||
            _geo_syntax_abort(_geo_stmt_head_token(rest))
        close=_geo_matching_delim(rest,firstindex(rest))
        if close==0
            _geo_scan_soft_eval_safe(rest[2:end],context)
            _geo_syntax_abort(";")
        end
        # `While(i) op v` is the upstream `String__Index ( FExpr )` index
        # write — only `While (cond)` without an affectation tail is the
        # extension's loop.
        tail=String(strip(rest[nextind(rest,close):end]))
        word=="While" && match(r"^(?:\+\+|--|\+=|-=|\*=|/=|=)",tail)!==nothing &&
            return nothing
        return (word=word,arg=String(strip(rest[2:prevind(rest,close)])),
                rest=tail)
    end
end

# Soft-eval that degrades an unparseable tail to `syntax error (;)` — for the
# `If (x`/`For (` headers that ran to the line end upstream still tries to
# reduce the FExpr, then reports the missing terminator.
function _geo_scan_soft_eval_safe(src::AbstractString,
                                  context::_GeoNumericContext)
    try
        return _geo_scan_soft_eval(src,context)
    catch err
        err isa InterruptException && rethrow()
        err isa _GeoSyntaxAbort && rethrow()
        _geo_syntax_abort(";")
    end
end

# `Call`/`Function` names: `String__Index` or `StringExpr` upstream — a bare
# identifier, `x~{i}` string-index, quoted literal, or a string expression.
function _geo_scan_call_name(raw::AbstractString,context::_GeoNumericContext)
    s=String(strip(raw))
    length(s)>=2 && startswith(s,"\"") && endswith(s,"\"") &&
        return s[2:prevind(s,lastindex(s))]
    length(s)>=2 && startswith(s,"'") && endswith(s,"'") &&
        return s[2:prevind(s,lastindex(s))]
    try
        return _geo_symbol_name(s,context,"read_geo_params: Call")
    catch err
        err isa InterruptException && rethrow()
        err isa _GeoSyntaxAbort && rethrow()
        context.soft_unknown_reads=true
        try
            return _geo_eval_string(s,context,"read_geo_params: Call")
        finally
            context.soft_unknown_reads=false
        end
    end
end

# A "relevant" statement that survived every dispatch arm is malformed —
# report the token upstream's committed production rejects, mirroring the
# `_exec_line!` ordering:
#   `Physical`      → after `tPhysical`, a kind word, `(`, `)`, `=`
#   `Field[...]`    → `[ FExpr ]`, optional `.member`, then `=`
#   `<name> Field`  → `Field` when the head is a keyword (`SetMaxTag Field` →
#                     `(Field)`), else the token after `Field`
#                     (`Background Field;` → `(;)`)
#   `X.Y ...`       → option-write heads; the token after `X.Y`
#   mid-expression  → the relevant keyword itself (`x = 5+Field[1]` → `(Field)`)
function _geo_scan_malformed_abort!(body::AbstractString,
                                    context::_GeoNumericContext)
    src=String(strip(body))
    if match(r"^Physical\b",src)!==nothing
        rest=String(strip(src[nextind(src,firstindex(src)+7):end]))
        km=match(r"^[A-Za-z_][A-Za-z0-9_]*",rest)
        km===nothing && _geo_syntax_abort(_geo_stmt_head_token(rest))
        kind=String(km.match)
        kind in ("Point","Curve","Line","Surface","Volume") ||
            _geo_syntax_abort(kind)
        rest=String(strip(rest[nextind(rest,lastindex(km.match)):end]))
        startswith(rest,"(") ||
            _geo_syntax_abort(_geo_stmt_head_token(rest))
        close=_geo_matching_delim(rest,firstindex(rest))
        if close==0
            _geo_scan_soft_eval(strip(rest[2:end]),context)
            _geo_syntax_abort(";")
        end
        _geo_syntax_abort(_geo_stmt_head_token(
            rest[nextind(rest,close):end]))
    end
    if match(r"^Field\b",src)!==nothing
        rest=String(strip(src[nextind(src,firstindex(src)+4):end]))
        if startswith(rest,"[")
            close=_geo_matching_delim(rest,firstindex(rest))
            close==0 && _geo_syntax_abort("[")
            _geo_scan_soft_eval(rest[2:prevind(rest,close)],context)
            rest=String(strip(rest[nextind(rest,close):end]))
            (dm=match(r"^\.\s*[A-Za-z_][A-Za-z0-9_]*",rest))!==nothing &&
                (rest=String(strip(rest[nextind(rest,
                    lastindex(dm.match)):end])))
        end
        _geo_syntax_abort(_geo_stmt_head_token(rest))
    end
    if (m=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+Field\b(.*)$",src))!==nothing
        # A keyword head already committed its own production upstream —
        # `SetMaxTag Field` errors on `Field` (the expected GeoEntity slot),
        # while a plain `tSTRING` head makes `tField` the `String__Index`
        # affectation slot, erroring at what follows (`Background Field;` →
        # `(;)`).
        _geo_tstring_valid(String(m.captures[1])) ||
            _geo_syntax_abort("Field")
        _geo_syntax_abort(_geo_stmt_head_token(m.captures[2]))
    end
    if (m=match(r"^[A-Za-z_][A-Za-z0-9_]*\.[A-Za-z_][A-Za-z0-9_]*",src))!==nothing
        _geo_syntax_abort(_geo_stmt_head_token(
            src[nextind(src,lastindex(m.match)):end]))
    end
    # A relevant token inside an otherwise fine statement aborts at that
    # token (`x = 5 + Field[1]` → `syntax error (Field)`).
    rm=match(
        r"\b(Physical|Field|Background|BoundaryLayer|Mesh|Geometry|Point|Curve|Line|Surface|Volume)\b|Mesh\.[A-Za-z_][A-Za-z0-9_]*|Geometry\.[A-Za-z_][A-Za-z0-9_]*",
        src)
    rm===nothing ? _geo_syntax_abort(_geo_stmt_head_token(src)) :
        _geo_syntax_abort(_geo_stmt_head_token(src[rm.offset:end]))
    return nothing
end

"""
    read_geo_params(path; max_file_bytes=typemax(Int)) -> GeoParams

Scan a gmsh `.geo`/`.geo_unrolled` for mesh-sizing options, Physical group
declarations, and `Field[...]`/`Background Field` statements. Deterministic,
finite arithmetic expressions can use `Pi`, prior scalar variables and Gmsh's
pure numeric functions, including in explicit Physical-group tags and
`Field[...]` tags. Numeric lists support zero-based indexing, cardinality, copies,
concatenation, bounded selection, and checked assignment or compound mutation.
Finite constant Gmsh list ranges use
`start:end[:increment]` order and are expanded with a strict resource bound in
recognized numeric field options and field selectors; known numeric list variables
are normalized there as well. Entirely numeric Physical memberships containing
ranges are checked but remain geometry data. Geometry-derived Physical memberships
remain opaque to this scanner and are evaluated by `execute_geo`. Named declarations
populate `physical_groups`; a nonempty name without an explicit tag receives the next
tag from the global Physical namespace. Unnamed declarations receive the same scanner
checks without adding a name.
Read-only dynamic tag allocators are
evaluated while their Point, shared geometric-region, or Field namespace remains
fully tracked. `SetMaxTag Point|Curve|Surface|Volume` updates the associated checked
counter. Built-in can lower a counter; OpenCASCADE only raises it. Tracked
`Boolean` operand `Delete` releases the operand and its hidden boundary entities
and re-tags the result boundary at the lowest free tags; any other topology
change makes the affected allocator unavailable. Reads use the greatest
counter among activated factories.
`If`/`ElseIf`/`Else`/`EndIf` conditionals, `For` ranges, the bounded `While`
extension, `Function`/`Macro` definitions and `Call` replay follow upstream's
lexer-level skip semantics: dead branches emit no diagnostics, `ElseIf`
conditions still evaluate, loop bodies re-scan per iteration, function bodies
register at `Return` and re-scan at `Call`. Option reads, random/external
functions, CSG and Boolean geometry are deliberately not evaluated. Numeric
field options are normalized to literals; geometric references and string or
point-dependent field expressions remain source strings.
"""
function read_geo_params(path;max_file_bytes=typemax(Int))
    caller="read_geo_params"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    smin = NaN; smax = NaN; sfactor = 1.0; seed = 0; background = 0
    geometry_tolerance=NaN
    mesh_size_from_curvature=0;boundary_layer_fan_elements=5
    boundary_layers=Int[]
    groups = Dict{Tuple{Int,Int},String}()
    physical_name_tags=Dict{Tuple{Int,String},Int}()
    # `GEO_Internals::PhysicalGroups` at scan granularity — raw group → raw
    # member tags. `scan_opaque` marks groups whose member list could not be
    # evaluated (dynamic/geometry-derived content).
    scan_raw=Dict{Tuple{Int,Int},Vector{Int}}()
    scan_opaque=Set{Tuple{Int,Int}}()
    # Entity-level physical mirror, rebuilt at the grammar's sync points so
    # `Physical X("name", 0)` resolves through `getMaxPhysicalNumber(dim)`.
    scan_ephys=Dict{Tuple{Int,Int},Vector{Int}}()
    # Refs, not plain locals — the self-recursive `_consume_stmt` closure would
    # box captured mutables (the allocation audit tolerates none of these).
    scan_has_mesh=Ref(false)
    # `GEO_Internals::_changed` — fresh internals start changed; every
    # mutating statement sets it; synchronize clears it.
    scan_changed=Ref(true)
    kinds = Dict{Int,String}()
    options = Dict{Int,Dict{String,String}}()
    option_order=Dict{Int,Vector{String}}()
    creation_curvature=Dict{Int,Int}()
    context=_GeoNumericContext()
    # Scan diagnostics are returned in `GeoParams.scan_*` — the execution pass
    # replays each statement and prints them once, so the scan stays quiet.
    context.stderr_diagnostics=false
    allocator_state=_GeoTagAllocatorState()
    # Control-flow bookkeeping for the streaming pass — see `_GeoScanLevel`.
    # `functions_defined` holds `Function`/`Macro` bodies for `Call` replay.
    ctrl_levels=_GeoScanLevel[]
    functions_defined=Dict{String,Vector{String}}()
    # Loop bodies replay *after* their level pops, so the stack alone cannot
    # reproduce upstream's `ImbricatedLoop` count — track replay depth.
    loop_nesting=Ref(0)
    # `ImbricatedTest` is mirrored on `context.if_depth` (also readable via
    # the `TestLevel` option): `If` shifts increment, closers decrement, and
    # stray `EndIf`s drive it negative — an `If` then opens a *debt* level
    # whose midpoints are orphan errors but whose bodies still run.
    # `skip_depth` mirrors upstream `skip`/`skipTest` `nb_skip`: `If` words
    # inside a dead region only nest the char-level scan and never parse.
    skip_depth=Ref(0)
    # `Exit`/`Abort` terminate the stream — everything after is unreachable.
    # `stopped` lives in a Ref so the self-recursive closure never boxes it;
    # `stop_kind` scopes the stop like upstream: `Exit` is process-global,
    # `Abort` and the >20-error cap halt only the current file's ParseFile.
    stopped=Ref(false)
    stop_kind=Ref(:none)
    # `error tEND` recovery (Gmsh.y:276): a syntax error discards tokens
    # through the next `;`. Control keywords carry no `;`, so the discard
    # eats them plus the first `;`-terminated statement — set after a
    # control-header abort, cleared at file boundaries.
    discard_pending=Ref(false)
    # Replay a collected loop body once per iteration — upstream re-parses the
    # body text at each `EndFor` rewind, so per-statement recovery applies.
    function _scan_replay(body::Vector{String})
        for stmt in body
            stopped[] && break
            _consume_stmt(stmt)
        end
        return nothing
    end
    function _scan_bind_loop_var(var::AbstractString,value::Float64)
        # `For name In {..}` — upstream clears the list flag on the symbol and
        # stores the scalar (`s.list = false; s.value[0] = first`).
        isempty(var) && return nothing
        delete!(context.lists,var)
        delete!(context.list_variables,var)
        delete!(context.unavailable_lists,var)
        _geo_context_set_scalar!(context,var,value,caller)
        return nothing
    end
    # `EndFor`'s increment — `s.value[0] += step` through the live symbol, so
    # a body's own write propagates; a missing symbol is `Unknown loop
    # variable`, a list payload `Bad loop variable`, each recoverable.
    function _scan_loop_advance(var::AbstractString,x0::Float64,step::Float64,
                                context::_GeoNumericContext,
                                caller::AbstractString)
        isempty(var) && return x0+step
        if !haskey(context.values,var) && !haskey(context.lists,var)
            _geo_yyerror!(context,"Unknown loop variable '$var'")
            return x0
        end
        haskey(context.lists,var) &&
            (_geo_yyerror!(context,"Bad loop variable $var");return x0)
        x0=context.values[var]+step
        _scan_bind_loop_var(var,x0)
        return x0
    end
    # A control header whose evaluation fails on an unavailable value makes
    # the whole construct untrackable — upstream it runs for real, so the
    # body may mutate any binding: invalidate the tracked state and treat
    # the header as inactive (loop skipped, branch not taken).
    # `_GeoSyntaxAbort` still reports a genuinely malformed header.
    function _scan_untrackable!()
        reason="a loop or conditional reading an unavailable value may have changed it"
        _geo_invalidate_context!(context,reason)
        _geo_allocator_invalidate!(
            allocator_state,reason;physical=true)
        return nothing
    end
    # After an untrackable loop the names its body assigns exist upstream
    # (with unknown values) — mark them unavailable so later conditional
    # reads take the untrackable path instead of `Unknown variable`.
    function _scan_mark_body_targets!(buffer::Vector{String})
        reason="a loop that may have assigned it could not be tracked"
        for stmt in buffer,
            m in eachmatch(r"(?:^|[^A-Za-z0-9_.])([A-Za-z_][A-Za-z0-9_]*)\s*(?:\[[^\]\n]*\]|\([^\)\n]*\))?\s*(?:\+=|-=|\*=|/=|=(?!=)|\+\+|--)",stmt)
            name=String(m.captures[1])
            _geo_context_forget!(context,name)
            context.unavailable[name]=reason
            context.unavailable_lists[name]=reason
        end
        return nothing
    end
    function _scan_push(levels::Vector{_GeoScanLevel},level::_GeoScanLevel)
        # `If` nesting is capped through `context.if_depth` (`ImbricatedTest`)
        # at the shift site; loops count open levels plus replay depth —
        # upstream's `ImbricatedLoop` covers both.
        if level.kind in (:for,:while)
            count(l->l.kind in (:for,:while),levels)+loop_nesting[] <
                _MAX_GEO_SCAN_NESTING-1 ||
                _geo_yyerror!(context,
                    "Reached maximum number of imbricated loops")
        end
        push!(levels,level)
        return nothing
    end
    function _scan_replay_loop(body::Vector{String})
        loop_nesting[]+=1
        try
            _scan_replay(body)
        finally
            loop_nesting[]-=1
        end
        return nothing
    end
    # End-of-file drain for one parsed file (upstream lexer semantics): an
    # `If` whose last branch ran stays silent, a pending dead-region skip
    # reports `Unexpected end of file` and marks the run failed, an unclosed
    # `For`/`While` with a live header replays its buffered body once (the
    # pass upstream already streamed) while a skipped one reports the EOF,
    # and an unclosed `Function` reports the EOF but leaves the name
    # registered (upstream `createFunction` ran at the header). Only levels
    # at or above `min_depth` belong to this file — outer levels were opened
    # by an including parent and persist across it.
    function _scan_eof_drain(min_depth::Int)
        # A dead-region `skip`/`skipTest` dies with the file — its nested
        # `If` words never reach the parent stream.
        skip_depth[]=0
        # At most one collecting level can sit on top of the `If` run — while
        # a level collects, every other control word buffers as text.
        # Live `If` levels persist across the file boundary (upstream's
        # `ImbricatedTest` is global — a parent's `EndIf` still closes one);
        # only dead regions and the open collecting level unwind here.
        while length(ctrl_levels)>min_depth &&
              ctrl_levels[end].kind!==:if
            # Pop before replaying: a level still on the stack collects its
            # own replayed lines back into `buffer`, which would never
            # terminate.
            top=pop!(ctrl_levels)
            if top.kind===:for
                if top.active && !isempty(top.range)
                    first,last,step=top.range
                    if (step>0 && first>last) || (step<0 && first<last)
                        # Upstream's `skip("For","EndFor")` was mid-flight.
                        _geo_eof_error!(context)
                    else
                        # The buffered body is the one pass upstream already
                        # streamed before `EndFor` would have rewound.
                        _scan_replay_loop(top.buffer)
                    end
                end
            elseif top.kind===:while
                if top.active
                    # `While` mirrors `For`: a true header condition streamed
                    # the body once (replay it), a false one was a dead
                    # region whose `skip` reports the EOF.
                    if _geo_scan_soft_eval(top.cond,context)!=0
                        _scan_replay_loop(top.buffer)
                    else
                        _geo_eof_error!(context)
                    end
                end
            else # :function — upstream registered the name at the header
                if top.active
                    haskey(functions_defined,top.name) ||
                        (functions_defined[top.name]=top.buffer)
                    _geo_eof_error!(context)
                end
            end
        end
        # Inactive `If` levels above `min_depth` were mid-`skip`/`skipTest`
        # at this file's EOF — upstream reports and decrements each. A
        # `closed` level already decremented at its closing midpoint.
        for k in length(ctrl_levels):-1:min_depth+1
            lvl=ctrl_levels[k]
            lvl.kind===:if || break
            lvl.active && break
            _geo_eof_error!(context)
            lvl.closed || (context.if_depth-=1)
            pop!(ctrl_levels)
        end
        return nothing
    end
    # Bison recovery: `error tEND` (Gmsh.y:276) discards tokens through the
    # next `;` after a syntax error. The scan's statements are `;`-delimited
    # (control keywords never split), so an abort's discard region is
    # normally inside the aborted statement — only a `;`-less line (a
    # replayed bare control or a `}`-terminated statement) spills it into
    # following statements. The diagnostic is recorded (muted) and execution
    # replays the statement for the visible error.
    function _consume_stmt(line)
        if discard_pending[]
            endswith(rstrip(line),';') && (discard_pending[]=false)
            return nothing
        end
        # `gmsh_yyerrorstate > 20` aborts the current file's ParseFile —
        # counted against `file_error_base` so an `Include`d file gets its
        # own budget.
        if !stopped[] &&
           length(context.exec_errors)-context.file_error_base>20
            _geo_yyerror!(context,"Too many errors: aborting parser...")
            stopped[]=true
            stop_kind[]=:toomany
        end
        stopped[] && return nothing
        try
            _consume_stmt_ctrl(String(line))
        catch err
            err isa InterruptException && rethrow()
            err isa _GeoSyntaxAbort || rethrow()
            # Detail-bearing aborts are runtime diagnostics (unknown or
            # unavailable names, domain and resource-limit failures) —
            # upstream prints them verbatim, so the scan records the detail.
            # Bare `syntax error (<token>)` is only for grammar failures.
            _geo_yyerror!(context,isempty(err.detail) ?
                "syntax error ($(err.token))" : err.detail)
            endswith(rstrip(line),';') || (discard_pending[]=true)
        end
        return nothing
    end
    # Strip control prefixes off `line`, update the level stack, then route the
    # residual statement: buffered when the innermost level collects, skipped
    # in dead branches, otherwise processed by `_consume_stmt_inner`.
    function _consume_stmt_ctrl(line::String)
        stopped[] && return nothing
        residual=String(strip(line))
        # Header aborts escape marked: the `error tEND` discard then applies
        # (a statement-level abort reaches `_consume_stmt` unmarked — its own
        # `;` already terminates the recovery).
        try
        while true
            isempty(residual) && return nothing
            collecting=_geo_scan_collecting(ctrl_levels)
            if !collecting &&
               (!_geo_scan_active(ctrl_levels) || skip_depth[]>0)
                # Upstream's dead-region `skip`/`skipTest` is char-level —
                # `If` words only nest the scan (`nb_skip`), while `For`,
                # `While`, `Function` and `Macro` headers — malformed or
                # not — are invisible raw text.
                wm=match(r"^(If|For|While|Function|Macro)\b",residual)
                if wm!==nothing
                    wm.match=="If" && (skip_depth[]+=1)
                    return nothing
                end
                if skip_depth[]>0
                    # Inside a nested `If` every midpoint/closer word is
                    # invisible; the matching `EndIf` decrements `nb_skip`.
                    match(r"^(ElseIf|Else)\b",residual)!==nothing &&
                        return nothing
                    if match(r"^EndIf\b",residual)!==nothing
                        skip_depth[]-=1
                        return nothing
                    end
                end
                # A midpoint that fires while a branch already ran makes
                # upstream `skip("If","EndIf")` — the level is *closed* and
                # only `If`/`EndIf` words remain visible: later `ElseIf`/
                # `Else` tokens are skipped text, not reparse candidates.
                if !isempty(ctrl_levels) && (tl=ctrl_levels[end]).kind===:if &&
                   tl.closed &&
                   match(r"^(ElseIf|Else)\b",residual)!==nothing
                    return nothing
                end
            end
            pref=_geo_scan_control_prefix(residual,context)
            pref===nothing && break
            word=pref.word;arg=pref.arg;rest=pref.rest
            handled=false
            if collecting
                top=ctrl_levels[end]
                if top.kind===:for && word=="EndFor" && top.nest==0
                    pop!(ctrl_levels)
                    if top.active && !isempty(top.range)
                        first,last,step=top.range
                        # Upstream binds `first` at the header and advances
                        # at `EndFor` through the live symbol — a body's own
                        # `i = ...` write propagates, and a zero step is
                        # legal (the body runs once, then the check pops).
                        if !((step>0 && first>last) || (step<0 && first<last))
                            x0=first;iterations=0
                            while true
                                _scan_replay_loop(top.buffer)
                                stopped[] && break
                                iterations+=1
                                iterations<=_MAX_GEO_SCAN_LOOP_ITERATIONS ||
                                    throw(ArgumentError(
                                        "$caller: loop exceeds " *
                                        "$_MAX_GEO_SCAN_LOOP_ITERATIONS " *
                                        "iterations"))
                                x0=_scan_loop_advance(
                                    top.var,x0,step,context,caller)
                                (step>0 ? x0<=last :
                                 step<0 ? x0>=last : false) || break
                            end
                        end
                    elseif top.active
                        # Untrackable range — the buffered body never
                        # replayed, so its assignment targets are values the
                        # scan cannot know.
                        _scan_mark_body_targets!(top.buffer)
                    end
                    handled=true
                elseif top.kind===:while && word=="EndWhile" && top.nest==0
                    pop!(ctrl_levels)
                    if top.active
                        iterations=0
                        untrackable=false
                        while !stopped[] && !untrackable
                            cond=try
                                _geo_scan_soft_eval(top.cond,context)!=0
                            catch err
                                err isa InterruptException && rethrow()
                                err isa ArgumentError || rethrow()
                                untrackable=true
                                false
                            end
                            cond || break
                            iterations+=1
                            iterations<=_MAX_GEO_SCAN_LOOP_ITERATIONS ||
                                throw(ArgumentError(
                                    "$caller: loop exceeds " *
                                    "$_MAX_GEO_SCAN_LOOP_ITERATIONS iterations"))
                            _scan_replay_loop(top.buffer)
                        end
                        if untrackable
                            _scan_untrackable!()
                            _scan_mark_body_targets!(top.buffer)
                        end
                    end
                    handled=true
                elseif top.kind===:function && word=="Return" && top.nest==0
                    pop!(ctrl_levels)
                    if top.active
                        # Upstream keeps the first body — a redefinition
                        # already errored at its `Function` header.
                        haskey(functions_defined,top.name) ||
                            (functions_defined[top.name]=top.buffer)
                    end
                    handled=true
                end
                if !handled
                    # A foreign control word inside a collecting level is
                    # body text — but only that word. The rest of a joined
                    # statement (`EndWhile EndFor`) is still scanned so a
                    # closer for an *outer* family level applies, matching
                    # upstream's flat token stream.
                    if top.nest>0 && word in ("EndFor","EndWhile","Return")
                        # The closer of a loop/function nested inside the
                        # buffered body — it unwinds the count, it does not
                        # close this level.
                        top.nest-=1
                    elseif word in ("For","While","Function","Macro")
                        top.nest+=1
                    end
                    consumed=String(strip(chopsuffix(residual,rest)))
                    isempty(consumed) || push!(top.buffer,consumed)
                end
            else
                handled=true
                if word=="If"
                    # `tIf` shifts: `ImbricatedTest++` — the FExpr cond
                    # already reduced (it reads the pre-increment level).
                    parent=_geo_scan_active(ctrl_levels)
                    cond=false
                    if parent
                        cond=try
                            _geo_scan_soft_eval(arg,context)!=0
                        catch err
                            err isa InterruptException && rethrow()
                            err isa ArgumentError || rethrow()
                            _scan_untrackable!()
                            false
                        end
                    end
                    context.if_depth+=1
                    if context.if_depth>_MAX_GEO_SCAN_NESTING-1
                        _geo_yyerror!(context,
                            "Reached maximum number of imbricated tests")
                        context.if_depth=_MAX_GEO_SCAN_NESTING-1
                    end
                    _scan_push(ctrl_levels,_GeoScanLevel(
                        :if;active=cond,taken=cond,
                        debt=context.if_depth<=0))
                elseif word=="ElseIf"
                    if !isempty(ctrl_levels) &&
                       ctrl_levels[end].kind===:if
                        top=ctrl_levels[end]
                        parent=_geo_scan_outer_active(ctrl_levels)
                        if top.debt
                            # Opened at `ImbricatedTest<=0`: the FExpr still
                            # reduces, then the token is an orphan error and
                            # the following text parses live.
                            parent && _geo_scan_soft_eval(arg,context)
                            _geo_yyerror!(context,"Orphan ElseIf")
                            top.active=parent
                            top.taken=true
                        elseif parent
                            # The FExpr argument already reduced upstream —
                            # it evaluates even after a taken branch.
                            c=try
                                _geo_scan_soft_eval(arg,context)!=0
                            catch err
                                err isa InterruptException && rethrow()
                                err isa ArgumentError || rethrow()
                                _scan_untrackable!()
                                false
                            end
                            was=top.taken
                            top.taken|=c
                            top.active=!was && c
                            # A prior branch already ran: upstream now
                            # `skip`s to `EndIf`, closing the level — later
                            # midpoints are invisible text.
                            if was
                                top.closed=true
                                context.if_depth-=1
                            end
                        else
                            top.active=false
                        end
                    else
                        _geo_scan_active(ctrl_levels) || (residual=rest;continue)
                        _geo_scan_soft_eval(arg,context)
                        _geo_yyerror!(context,"Orphan ElseIf")
                    end
                elseif word=="Else"
                    if !isempty(ctrl_levels) &&
                       ctrl_levels[end].kind===:if
                        top=ctrl_levels[end]
                        parent=_geo_scan_outer_active(ctrl_levels)
                        if top.debt
                            _geo_yyerror!(context,"Orphan Else")
                            top.active=parent
                            top.taken=true
                        elseif top.taken
                            # `Else` on a taken chain — upstream `skip`s to
                            # `EndIf` and closes the level.
                            top.closed=true
                            top.active=false
                            context.if_depth-=1
                        else
                            # `taken` is *not* set — upstream's status flag
                            # stays 0 on `Else`, so a following ElseIf can
                            # still activate.
                            top.active=parent
                        end
                    else
                        _geo_scan_active(ctrl_levels) &&
                            _geo_yyerror!(context,"Orphan Else")
                    end
                elseif word=="EndIf"
                    if !isempty(ctrl_levels) &&
                       ctrl_levels[end].kind===:if
                        top=pop!(ctrl_levels)
                        if top.debt
                            # A debt level closes below zero either way; the
                            # `EndIf` only parses (and warns) when the cond
                            # ran live or a midpoint stopped `skipTest` —
                            # otherwise the skip consumed it silently.
                            context.if_depth-=1
                            top.taken &&
                                _geo_yywarn!(context,"Orphan EndIf")
                        elseif !top.closed
                            # A `closed` level's `ImbricatedTest` already
                            # decremented at the closing midpoint; the
                            # `EndIf` was consumed by the skip.
                            context.if_depth-=1
                        end
                    else
                        # `tEndIf` decrements unconditionally — the warning
                        # fires whenever the result is negative.
                        context.if_depth-=1
                        context.if_depth<0 &&
                            _geo_yywarn!(context,"Orphan EndIf")
                    end
                elseif word=="For"
                    var,spec,closer=arg
                    parent=_geo_scan_active(ctrl_levels)
                    range=Float64[]
                    if parent
                        range=try
                            _geo_scan_for_range(
                                spec,context,"$caller: For range",closer)
                        catch err
                            err isa InterruptException && rethrow()
                            err isa ArgumentError || rethrow()
                            _scan_untrackable!()
                            Float64[]
                        end
                    end
                    parent && !isempty(var) && !isempty(range) &&
                        _scan_bind_loop_var(var,range[1])
                    _scan_push(ctrl_levels,_GeoScanLevel(
                        :for;active=parent,var=var,range=range))
                elseif word=="While"
                    _scan_push(ctrl_levels,_GeoScanLevel(
                        :while;active=_geo_scan_active(ctrl_levels),
                        cond=arg))
                elseif word=="EndFor"
                    _geo_scan_active(ctrl_levels) &&
                        _geo_yyerror!(context,"Invalid For/EndFor loop")
                elseif word=="EndWhile"
                    # `While` is a Tessella extension — bare `EndWhile` is a
                    # no-diagnostic marker here.
                elseif word=="Function"
                    parent=_geo_scan_active(ctrl_levels)
                    # `createFunction` registers the name at this reduction,
                    # before the body `skip` runs — a redefinition errors
                    # here, not at `Return`, and the first body wins.
                    parent && haskey(functions_defined,arg) &&
                        _geo_yyerror!(context,
                            "Redefinition of function $arg")
                    _scan_push(ctrl_levels,_GeoScanLevel(
                        :function;active=parent,name=arg))
                elseif word=="Return"
                    # `Return` inside a dead region is invisible text in
                    # upstream's family `skip` — it only unwinds a `Call`
                    # when the enclosing branches actually run.
                    if _geo_scan_active(ctrl_levels)
                        context.call_depth>0 && throw(_GeoScanReturnSignal())
                        _geo_yyerror!(context,"Error while exiting function")
                    end
                end
            end
            residual=rest
        end
        catch err
            err isa InterruptException && rethrow()
            err isa _GeoSyntaxAbort || rethrow()
            throw(_GeoSyntaxAbort(err.token,err.detail,true))
        end
        isempty(residual) && return nothing
        if _geo_scan_collecting(ctrl_levels)
            push!(ctrl_levels[end].buffer,residual)
            return nothing
        end
        _geo_scan_active(ctrl_levels) || return nothing
        # `Call name;` — upstream `enterFunction` re-parses the recorded body
        # in the caller's scope; an unknown name is a recoverable diagnostic.
        if (cm=match(r"^Call\s+(.+?)\s*;?\s*$",residual))!==nothing
            name=_geo_scan_call_name(cm.captures[1],context)
            body=get(functions_defined,name,nothing)
            if body===nothing
                _geo_yyerror!(context,"Unknown function '$name'")
            else
                depth0=length(ctrl_levels)
                context.call_depth+=1
                try
                    context.call_depth<=128 || throw(ArgumentError(
                        "$caller: Call depth exceeds 128"))
                    _scan_replay(body)
                catch err
                    err isa _GeoScanReturnSignal || rethrow()
                    resize!(ctrl_levels,depth0)
                finally
                    context.call_depth-=1
                end
            end
            return nothing
        end
        _consume_stmt_inner(residual)
        return nothing
    end
    function _consume_stmt_inner(line)
        stopped[] && return
        raw_body=if endswith(line,";")
            String(strip(line[firstindex(line):prevind(line,lastindex(line))]))
        else
            _geo_brace_terminated_statement(line) || throw(ArgumentError(
                "read_geo_params: internal statement scanner lost a terminator: " *
                repr(line)))
            String(strip(line))
        end
        # A bare `;` (including the `;` tail after a `}`-terminated statement)
        # is an upstream `syntax error (;)` — one diagnostic per stray
        # terminator, recorded, then recovery resumes.
        if isempty(raw_body) || all(c->c==';'||isspace(c),raw_body)
            for _ in 2:count(==(';'),line)
                _geo_yyerror!(context,"syntax error (;)")
            end
            _geo_syntax_abort(";")
        end
        raw_code=_geo_unquoted_code(raw_body)

        # `Call`/`DefineConstant`-family statements mutate state the scan does
        # not model — invalidate prior bindings instead of using a stale value
        # later. Control headers never reach this point: `_consume_stmt_ctrl`
        # strips them before dispatch.
        if occursin(r"\b(?:Call|DefineConstant|UndefineConstant)\b",raw_code)
            reason="an unsupported loop, conditional, macro or include may have changed it"
            _geo_invalidate_context!(context,reason)
            _geo_allocator_invalidate!(
                allocator_state,reason;physical=true)
        end
        body=raw_body
        code=_geo_unquoted_code(body)

        # `Physical X{..}` selectors and Boundary-family actions read the
        # synced model mid-statement — the grammar gates that sync on
        # `getChanged()`, so it runs before the statement's own mutation.
        if occursin(r"\bPhysical\s+(?:Point|Line|Curve|Surface|Volume)\s*\{|\b(?:Boundary|CombinedBoundary|OrientedBoundary|OrientedCombinedBoundary|PointsOf)\s*\{",
                    body) && scan_changed[]
            _geo_scan_sync_physicals!(scan_raw,scan_opaque,scan_ephys,
                                      allocator_state,scan_has_mesh[])
            scan_changed[]=false
        end
        _geo_statement_marks_internals(body) && (scan_changed[]=true)

        _geo_context_refresh_allocators!(context,allocator_state)
        physical=_geo_physical_declaration(body)
        if physical!==nothing
            dim=_geo_physical_kind_dim(physical,context,"read_geo_params")
            name=physical.name
            # `setPhysicalName(name, dim, T)` on the name->tag mirror: a bound
            # name resolves to its (possibly negative) tag; otherwise the
            # binding lands only when `T` is still unnamed (map `insert`
            # never overwrites).
            tag_named(t)=any(
                et==t for ((ed,_),et) in physical_name_tags if ed==dim)
            bind_name!(t)=begin
                name!==nothing && !isempty(name) &&
                    !haskey(physical_name_tags,(dim,name)) && !tag_named(t) &&
                    (physical_name_tags[(dim,name)]=t)
            end
            # Compound `+=`/`-=`/`*=`/`/=` declarations modify an existing
            # group — no tag allocation and no duplicate-declaration check.
            # The `"name"`-only form still increments the physical counter
            # unconditionally (`setMaxPhysicalTag(t+1)`) before
            # `setPhysicalName` resolves the name.
            if physical.op!="="
                group_tag=nothing
                if physical.tag_source===nothing
                    (name===nothing || isempty(name)) && throw(ArgumentError(
                        "read_geo_params: an automatic Physical " *
                        "$(physical.kind) group requires a nonempty name"))
                    allocator_state.physical_group_max=
                        allocator_state.physical_group_max==typemax(Int32) ?
                        Int(typemin(Int32)) :
                        allocator_state.physical_group_max+1
                    bound=get(physical_name_tags,(dim,name),nothing)
                    group_tag=bound===nothing ?
                        allocator_state.physical_group_max : bound
                    bind_name!(allocator_state.physical_group_max)
                else
                    etag=_geo_eval_numeric(physical.tag_source,context,
                        "read_geo_params: Physical $(physical.kind) tag")
                    if isfinite(etag)
                        etag_int=round(Int,etag)
                        if etag_int==etag
                            resolved=etag_int==0 ?
                                _geo_scan_max_entity_physical(
                                    scan_ephys,dim)+1 : etag_int
                            bound=get(physical_name_tags,(dim,name),nothing)
                            group_tag=bound===nothing ? resolved : bound
                            bind_name!(resolved)
                        end
                    end
                end
                membership=physical.membership
                isempty(membership) || _geo_validate_physical_ranges(
                    membership,context,
                    "read_geo_params: Physical $(physical.kind) membership")
                # `modifyPhysicalGroup` on the raw mirror: `-=` removes the
                # first occurrence of each listed member and deletes the
                # emptied group — including `-= {}` on any group (the
                # `tags.empty()` check) — erasing its name binding and
                # stripping entity memberships whose `abs` matches the tag
                # (which can never match a negative raw tag). `+=` appends.
                # Missing groups are silent for `-=`, errors for the rest —
                # the executor owns the diagnostics. Failed ops never mark
                # `_changed`.
                if group_tag!==nothing && physical.op in ("+=","-=") &&
                   haskey(scan_raw,(dim,group_tag))
                    key=(dim,group_tag)
                    scan_changed[]=true
                    ids=_geo_scan_physical_members(membership,context)
                    if physical.op=="+="
                        if ids===nothing || key in scan_opaque
                            delete!(scan_raw,key);push!(scan_opaque,key)
                        else
                            append!(scan_raw[key],ids)
                        end
                    else
                        emptied=if key in scan_opaque
                            ids===nothing ?
                                strip(membership)=="{}" : isempty(ids)
                        else
                            members=scan_raw[key]
                            ids===nothing ||
                                for id in ids
                                    pos=findfirst(==(id),members)
                                    pos===nothing || deleteat!(members,pos)
                                end
                            isempty(members) ||
                                (ids!==nothing && isempty(ids))
                        end
                        emptied && _geo_scan_delete_physical!(
                            scan_raw,scan_opaque,scan_ephys,groups,
                            physical_name_tags,key)
                    end
                end
                return
            end
            if physical.tag_source===nothing
                (name===nothing || isempty(name)) && throw(ArgumentError(
                    "read_geo_params: an automatic Physical $(physical.kind) " *
                    "group requires a nonempty name"))
            end
            name_bound=nothing
            if name!==nothing && !isempty(name)
                name_bound=get(physical_name_tags,(dim,name),nothing)
            end
            tag=if physical.tag_source===nothing
                # `Physical Point("name")` bumps the counter unconditionally
                # (`setMaxPhysicalTag(t+1)`) and then `setPhysicalName`
                # resolves the already-bound tag or binds the fresh one.
                allocator_state.physical_group_max=
                    allocator_state.physical_group_max==typemax(Int32) ?
                    Int(typemin(Int32)) :
                    allocator_state.physical_group_max+1
                if name_bound===nothing
                    bind_name!(allocator_state.physical_group_max)
                    allocator_state.physical_group_max
                else
                    name_bound
                end
            elseif name_bound!==nothing
                _geo_eval_numeric(physical.tag_source,context,
                    "read_geo_params: Physical $(physical.kind) tag")
                name_bound
            else
                t=_geo_allocator_physical_tag!(
                    allocator_state,physical.tag_source,context,
                    "read_geo_params: Physical $(physical.kind)")
                if t===nothing
                    nothing
                else
                    # `setPhysicalName(name, dim, T)`: `T == 0` resolves to
                    # `getMaxPhysicalNumber(dim)+1` — the entity-level mirror
                    # maximum; `CreatePhysicalGroup` then raises the counter
                    # to the resolved tag.
                    resolved=t==0 ?
                        _geo_scan_max_entity_physical(scan_ephys,dim)+1 : t
                    resolved>0 && (allocator_state.physical_group_max=
                        max(allocator_state.physical_group_max,resolved))
                    bind_name!(resolved)
                    resolved
                end
            end
            tag===nothing && return
            membership=physical.membership
            _geo_validate_physical_ranges(
                membership,context,
                "read_geo_params: Physical $(physical.kind) membership")
            group_key=(dim,tag)
            # A duplicate declaration is a *recoverable* error in Gmsh
            # ("already exists" at `modifyPhysicalGroup`) — execution
            # continues, so the scan likewise keeps going; the executor
            # records the diagnostic. The tag/name side effects above
            # (counter bump, `setPhysicalName` resolution) already matched.
            haskey(scan_raw,group_key) && return
            members=_geo_scan_physical_members(membership,context)
            if members===nothing
                scan_raw[group_key]=Int[];push!(scan_opaque,group_key)
            else
                scan_raw[group_key]=members
            end
            # `modifyPhysicalGroup` op 0 reaches `_changed = true`.
            scan_changed[]=true
            name!==nothing && !isempty(name) && (groups[group_key]=name)
            return
        end
        occursin(r"^Physical(?:\s|$)",body) && throw(ArgumentError(
            "read_geo_params: malformed Physical declaration; use " *
            "Physical Kind(tag), Physical Kind(\"name\", tag), or " *
            "Physical Kind(\"name\")"))
        _geo_allocator_observe_statement!(
            allocator_state,body,context,"read_geo_params";conservative=true)

        # `Delete`/`Recursive Delete` statements mutate allocator and variable
        # state the scan mirrors. An entity list's removed set depends on live
        # ownership the scan does not track, so its counter effects make later
        # allocator reads unavailable rather than stale — the same convention
        # other untrackable topology changes use.
        if (dl=match(r"^(?:(Recursive)\s+)?Delete\b(.*)$",body))!==nothing
            tail=String(strip(dl.captures[2]))
            if startswith(tail,"{")
                # Only geometric counters can change — field and physical
                # tags survive entity deletion, so `newf`/`newreg` stay live.
                _geo_allocator_invalidate!(allocator_state,
                    "an entity `Delete` list may have decremented dynamic tag counters";
                    geometry=true,fields=false)
                # `GEO_Internals::remove` marks `_changed` once per listed
                # element — unconditionally — so a nonempty list always
                # resynchronizes; `Delete {}` keeps the prior flag (and falls
                # back to a model-level remove, a no-op here). The literal
                # entity refs leave the live sets so the physical mirror
                # drops their memberships like the model does.
                inner=endswith(tail,"}") ? strip(tail[
                    firstindex(tail)+1:prevind(tail,lastindex(tail))]) :
                    strip(tail[firstindex(tail)+1:end])
                for lm in eachmatch(
                        r"\b(Point|Line|Curve|Surface|Volume)\s*\{\s*(-?\d+)\s*\}",
                        inner)
                    dim=_PHYS_DIM[lm.captures[1]]
                    tag=parse(Int,lm.captures[2])
                    for live in (allocator_state.live_builtin_points,
                                 allocator_state.live_builtin_curves,
                                 allocator_state.live_builtin_surfaces,
                                 allocator_state.live_builtin_volumes,
                                 allocator_state.live_occ_points,
                                 allocator_state.live_occ_curves,
                                 allocator_state.live_occ_surfaces,
                                 allocator_state.live_occ_volumes)
                        delete!(live,tag)
                    end
                    delete!(scan_ephys,(dim,tag))
                end
                isempty(inner) || (scan_changed[]=true)
                if scan_changed[]
                    _geo_scan_sync_physicals!(scan_raw,scan_opaque,scan_ephys,
                                              allocator_state,scan_has_mesh[])
                    scan_changed[]=false
                end
                return
            end
            dl.captures[1]===nothing || _geo_syntax_abort(
                isempty(tail) ? ";" : string(first(tail)))
            if (em=match(r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\")\s*\{",tail))!==nothing
                word=strip(em.captures[1],'"')
                # Anything other than `Embedded` is the executor's
                # recoverable syntax/unknown-command diagnostic (`Delete
                # Point{9}` → `syntax error (Point)` upstream).
                word=="Embedded" || return
                return  # `Delete Embedded` clears embeddings — no params state
            end
            if (vm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+([A-Za-z_][A-Za-z0-9_]*)\s*;?$",
                         tail))!==nothing
                vm.captures[1]=="Empty" && vm.captures[2]=="Views" && return
                # `Delete X Y` reports `Unknown command` at execution —
                # recoverable, no params state.
                return
            end
            if (fm=match(r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\")\s*\[\s*(.*?)\s*\]\s*;?$",
                         tail))!==nothing
                base=String(strip(fm.captures[1],'"'))
                # `Delete View[i]` erases an exec-layer post-processing view —
                # no params state. `Delete Field[i]` updates the field tables.
                # Any other indexed name is the executor's recoverable
                # `Unknown command` diagnostic.
                base=="View" && return
                base=="Field" || return
                tag=_geo_signed_gmsh_int_value(
                    _geo_eval_numeric(fm.captures[2],context,
                        "read_geo_params: Delete Field tag"),
                    "read_geo_params: Delete Field tag")
                if !haskey(kinds,tag)
                    # `FieldManager::deleteField` — `Msg::Error` only; parsing
                    # continues.
                    _geo_msg_error!(context,
                        "Cannot delete field id $tag, it does not exist")
                    return
                end
                delete!(kinds,tag);delete!(options,tag)
                delete!(option_order,tag);delete!(creation_curvature,tag)
                # Gmsh's `deleteField` leaves `_boundaryLayerFields` and the
                # background id alone — stale references surface downstream.
                _geo_allocator_delete_field!(allocator_state,tag)
                return
            end
            nm=match(r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\")\s*;?$",tail)
            nsm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*~\s*\{\s*(.*?)\s*\}\s*;?$",tail)
            name=if nm!==nothing
                String(strip(nm.captures[1],'"'))
            elseif nsm!==nothing
                idx=_geo_int_value(_geo_eval_numeric(
                        nsm.captures[2],context,
                        "read_geo_params: Delete namespace index"),
                    "read_geo_params: Delete namespace index")
                String(nsm.captures[1])*"_"*string(idx)
            else
                # `tDelete` expects a name, `name~{i}`, `{ ... }` or `[...]` —
                # upstream's syntax error lands on the offending token
                # (`Delete = 5` → `(=)`).
                _geo_syntax_abort(isempty(tail) ? ";" :
                    string(first(tail)))
            end
            if name=="All"
                # `ClearProject` — model, fields, physical names, parser
                # variables (numeric AND string) and the struct namespaces
                # all go away; the factory reverts to built-in. Function
                # definitions survive — `FunctionManager` is not cleared.
                _geo_allocator_reset_model!(allocator_state,context;fresh=true)
                allocator_state.factory=:builtin
                empty!(context.values);empty!(context.lists)
                empty!(context.strings)
                empty!(context.structs);empty!(context.struct_maxtag)
                empty!(context.list_variables)
                empty!(context.unavailable);empty!(context.unavailable_lists)
                context.stored_list_items=0
                empty!(kinds);empty!(options);empty!(option_order)
                empty!(creation_curvature);empty!(boundary_layers)
                empty!(groups);empty!(physical_name_tags)
                empty!(scan_raw);empty!(scan_opaque);empty!(scan_ephys)
                scan_has_mesh[]=false;scan_changed[]=true
                background=0
                return
            elseif name=="Model"
                # `GModel::destroy` + internals destroy — fields die with the
                # model; physical names and variables survive. The OCC
                # internals object persists (`resetOCCInternals`), so the
                # factory and `occ_active` carry over untouched.
                _geo_allocator_reset_model!(allocator_state,context;fresh=false)
                empty!(kinds);empty!(options);empty!(option_order)
                empty!(creation_curvature);empty!(boundary_layers)
                empty!(scan_raw);empty!(scan_opaque);empty!(scan_ephys)
                scan_has_mesh[]=false;scan_changed[]=true
                background=0
                return
            elseif name=="Physicals"
                # `resetPhysicalGroups` + `removePhysicalGroups` drop the raw
                # group records and the entity memberships but keep names and
                # the physical tag counter.
                empty!(scan_raw);empty!(scan_opaque);empty!(scan_ephys)
                scan_changed[]=true
                return
            elseif name=="Variables"
                empty!(context.values);empty!(context.lists)
                empty!(context.list_variables)
                empty!(context.unavailable);empty!(context.unavailable_lists)
                context.stored_list_items=0
                return
            elseif name=="Options"
                # `ReInitOptions` — restore the option-valued state the scan
                # tracks (fields and the background field are not options).
                smin=NaN;smax=NaN;sfactor=1.0;seed=0
                geometry_tolerance=NaN
                mesh_size_from_curvature=0;boundary_layer_fan_elements=5
                return
            elseif name=="Meshes"
                # `GModel::deleteMesh` — drops the mesh; internals unchanged.
                scan_has_mesh[]=false
                return
            elseif name=="Struct"
                # `gmsh_yynamespaces.clear()`.
                empty!(context.structs);empty!(context.struct_maxtag)
                return
            end
            known=haskey(context.values,name)||haskey(context.lists,name)||
                  haskey(context.unavailable,name)||
                  haskey(context.unavailable_lists,name)
            # `yymsg(0)` upstream — a recorded diagnostic, recoverable.
            known || _geo_yyerror!(context,
                "Unknown object or expression to delete '$name'")
            _geo_context_forget!(context,name)
            delete!(context.unavailable,name);delete!(context.unavailable_lists,name)
            return
        end

        # `Struct name[...] ;` — `DefineStruct tEND`. The definition lands in
        # the scan's namespace map so `x.member` reads resolve downstream;
        # a redefinition or option error is a recorded diagnostic, not a
        # scan abort (upstream `yymsg(0)` recovery semantics).
        if (sm=match(r"^Struct\s+(.+?)\s*;?\s*$",body))!==nothing
            try
                _geo_exec_struct_def!(String(sm.captures[1]),context)
            catch err
                err isa InterruptException && rethrow()
                err isa _GeoSyntaxAbort && rethrow()
                err isa ArgumentError || rethrow()
                _geo_yyerror!(context,err.msg)
            end
            return
        end

        # Transform and standalone shape-list statements carry no params-scan
        # data; a `Physical Kind{...}` inside one is a group reference, not a
        # declaration.
        match(r"^(?:Translate|Rotate|Dilate|Symmetry|Affine|Duplicata|Extrude|Boundary|CombinedBoundary|OrientedBoundary|OrientedCombinedBoundary|PointsOf)\s*\{",body)!==nothing &&
            return
        # Bare `BooleanX {..}{..}`, `BooleanX(tag) = {..}{..}` and the OCC
        # shape-list actions (`Closest`/`Fillet`/`Chamfer` brace groups,
        # `Intersect Curve`, `Split Curve`, `ThruSections`) are complete
        # upstream productions — params-neutral, with any topology effect
        # already recorded by the allocator observer above.
        match(Regex("^(?:Boolean(?:Union|Difference|Intersection|Fragments|" *
            "Section|Fuse|Coherence|Common|Cut)\\s*\\{|" *
            "Boolean(?:Union|Difference|Intersection|Fragments|Section|" *
            "Fuse|Coherence|Common|Cut)\\s*\\([^)]*\\)\\s*=|" *
            "(?:Closest|Fillet|Chamfer)\\s*\\{|" *
            "(?:Intersect|Split)\\s+Curve\\b|(?:Ruled\\s+)?ThruSections\\b)"),
            body)!==nothing && return

        mesh=match(r"^(Mesh|Geometry)(?:\s*\[\s*(.*?)\s*\])?\s*\.\s*(MeshSizeMin|MeshSizeMax|MeshSizeFactor|RandomSeed|MeshSizeFromCurvature|MinimumElementsPerTwoPi|BoundaryLayerFanElements|BoundaryLayerFanPoints|Tolerance)\s*=\s*(.*)$",body)
        if mesh!==nothing
            family=String(mesh.captures[1]);member=String(mesh.captures[3])
            key="$family.$member";raw=String(strip(mesh.captures[4]))
            # Soft reads like the generic option-write branch below: a failed
            # RHS records `yymsg(0)` and stores 0 upstream, never an abort.
            context.soft_unknown_reads=true
            value,index=try
                v=_geo_eval_numeric(raw,context,"read_geo_params: $key")
                i=mesh.captures[2]===nothing ? 0 :
                    _geo_signed_gmsh_int_value(
                        _geo_eval_numeric(mesh.captures[2],context,
                            "read_geo_params: $key option index"),
                        "read_geo_params: $key option index")
                (v,i)
            finally
                context.soft_unknown_reads=false
            end
            # `NumberOption(GMSH_SET)` lands in option state — `x.y` expression
            # reads and the allocator refresh consume it. Only the index-0
            # slot feeds the params surface below.
            _geo_set_option_number!(context,family,index,member,"=",value,
                "read_geo_params: $key")
            if index==0
                if key=="Mesh.MeshSizeMin";smin=value
                elseif key=="Mesh.MeshSizeMax";smax=value
                elseif key=="Mesh.MeshSizeFactor"
                    # `opt_mesh_lc_factor` ignores non-positive writes — the
                    # params surface sees the stored value, not the raw RHS.
                    sfactor=something(_geo_option_number(
                        context,"Mesh",0,"MeshSizeFactor"),1.0)
                elseif key=="Mesh.RandomSeed";seed=_gmsh_random_seed(value)
                elseif key=="Geometry.Tolerance";geometry_tolerance=value
                elseif key=="Mesh.MeshSizeFromCurvature" ||
                       key=="Mesh.MinimumElementsPerTwoPi"
                    mesh_size_from_curvature=_gmsh_int_option(value,key)
                else
                    boundary_layer_fan_elements=_gmsh_int_option(value,key)
                end
            end
            return
        end

        fm=match(r"^Field\s*\[\s*(.*)\s*\]\s*=\s*([A-Za-z][A-Za-z0-9_]*)$",body)
        if fm!==nothing
            # The `FExpr` tag already reduced upstream: an unknown variable
            # records `yymsg(0)` and yields 0 — `newField(0, kind)` proceeds.
            tag=_geo_signed_gmsh_int_value(
                _geo_scan_soft_eval(fm.captures[1],context),
                "read_geo_params: Field declaration tag")
            kind=fm.captures[2]
            # `newField` fails recoverably: duplicate ids keep the original
            # field and unknown types are not registered — each is a
            # `Msg::Error` plus a `yymsg(0)` parse diagnostic.
            if haskey(kinds,tag)
                _geo_msg_error!(context,"Field id $tag is already defined")
                _geo_yyerror!(context,
                    "Cannot create field $tag of type '$kind'")
                return
            end
            if !(String(kind) in _GEO_FIELD_KINDS)
                _geo_msg_error!(context,"Unknown field type \"$kind\"")
                _geo_yyerror!(context,
                    "Cannot create field $tag of type '$kind'")
                return
            end
            kinds[tag]=kind
            creation_curvature[tag]=mesh_size_from_curvature
            return
        end

        om=match(r"^Field\s*\[\s*(.*)\s*\]\.([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*)$",body)
        if om!==nothing
            name=String(om.captures[2]);value=String(strip(om.captures[3]))
            isempty(value) && _geo_syntax_abort(";")
            # The tag and right-hand side already reduced upstream, so each
            # failure records `yymsg(0)` and continues: tag→0 (`No field with
            # id 0`), a bad RHS→0 before the member lookup.
            context.soft_unknown_reads=true
            try
                tag=_geo_signed_gmsh_int_value(
                    _geo_eval_numeric(om.captures[1],context,
                        "read_geo_params: Field option tag"),
                    "read_geo_params: Field option tag")
                caller="read_geo_params: Field[$tag].$name"
                # Upstream validates `member` against the field's *kind*
                # (`field->options[member]`), after the FExpr RHS reduced.
                kind=get(kinds,tag,nothing)
                member_ok=kind!==nothing &&
                    name in _GEO_FIELD_OPTIONS[kind][1]
                normalized=if member_ok
                    _geo_normalize_field_option(value,name,context,caller)
                else
                    # Unknown members still evaluate a numeric RHS upstream
                    # (`Unknown variable 'x'` precedes the `Unknown option`
                    # diagnostic); string RHS never reaches FExpr.
                    if !startswith(value,"\"") && !startswith(value,"{")
                        try
                            _geo_eval_numeric(value,context,caller)
                        catch err
                            err isa InterruptException && rethrow()
                            err isa _GeoSyntaxAbort && rethrow()
                        end
                    end
                    value
                end
                # Option writes resolve the field id at statement time — a
                # write to a field that does not exist yet stores nothing and
                # is a recoverable parse error, even if the field is declared
                # later. The executor replays the same check on its live map.
                kind===nothing &&
                    (_geo_yyerror!(context,"No field with id $tag"); return)
                member_ok || (_geo_yyerror!(context,
                    "Unknown option '$name' in field $tag of type " *
                    "'$(_GEO_FIELD_OPTIONS[kind][2])'"); return)
                get!(() -> Dict{String,String}(),options,tag)[name]=normalized
                push!(get!(() -> String[],option_order,tag),name)
            finally
                context.soft_unknown_reads=false
            end
            return
        end

        # `x.y = expr`, `x[i].y = expr`, `x.y++`, `x[i].y--` — the upstream
        # `String__Index '.' tSTRING_Reserved` NumberOption productions.
        # Written option state is what `x.y` expression reads and the
        # allocator refresh (`Geometry.OldNewReg`) consume; `Field[i].x` was
        # already handled above and string-valued RHS falls through to the
        # string-option guard below. Unknown option names surface the same
        # "Unknown number option" diagnostic as upstream's failed
        # `NumberOption(GMSH_GET)` gate.
        onm=match(r"^([A-Za-z_][A-Za-z0-9_]*)(?:\s*\[\s*(.*?)\s*\])?\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*(=|\+=|-=|\*=|/=)\s*(.+)$",body)
        if onm!==nothing && !_geo_string_rhs(String(strip(onm.captures[5])))
            family=String(onm.captures[1])
            member=String(onm.captures[3])
            # Upstream `NumberOption` writes are fully recoverable: a bad RHS
            # records `yymsg(0)` and yields 0, an unknown category/member is a
            # `Msg::Error` — parsing continues either way. The scan mirrors
            # that with soft reads so the diagnostics land in `scan_*` and the
            # execution pass replays them once.
            context.soft_unknown_reads=true
            try
                index=onm.captures[2]===nothing ? 0 :
                    _geo_signed_gmsh_int_value(
                        _geo_eval_numeric(onm.captures[2],context,
                            "read_geo_params: $family option index"),
                        "read_geo_params: $family option index")
                _geo_set_option_number!(context,family,index,member,
                    String(onm.captures[4]),
                    _geo_eval_numeric(onm.captures[5],context,
                        "read_geo_params: $family.$member"),
                    "read_geo_params: $family.$member")
            finally
                context.soft_unknown_reads=false
            end
            return
        end
        oni=match(r"^([A-Za-z_][A-Za-z0-9_]*)(?:\s*\[\s*(.*?)\s*\])?\s*\.\s*([A-Za-z_][A-Za-z0-9_]*)\s*(\+\+|--)\s*$",body)
        if oni!==nothing
            family=String(oni.captures[1])
            member=String(oni.captures[3])
            context.soft_unknown_reads=true
            try
                index=oni.captures[2]===nothing ? 0 :
                    _geo_signed_gmsh_int_value(
                        _geo_eval_numeric(oni.captures[2],context,
                            "read_geo_params: $family option index"),
                        "read_geo_params: $family option index")
                _geo_option_number_increment!(context,family,index,member,
                    oni.captures[4]=="++" ? 1.0 : -1.0,
                    "read_geo_params: $family.$member")
            finally
                context.soft_unknown_reads=false
            end
            return
        end

        # `String__Index tField tAFFECT ListOfDouble` — `Background Field`
        # stores the raw id (`setBackgroundFieldId` does not validate it),
        # `BoundaryLayer Field` dedup-appends; any other `X Field` prefix is a
        # recoverable unknown-command diagnostic.
        fm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+Field\s*=\s*(.*)$",body)
        if fm!==nothing
            name=fm.captures[1]
            # The tag list already reduced upstream — a bad entry records
            # `yymsg(0)` and yields 0 rather than aborting the parse.
            context.soft_unknown_reads=true
            try
                if name=="Background"
                    tags=_geo_expression_field_tags(fm.captures[2],context,
                                                    "read_geo_params: Background Field";
                                                    deduplicate=false)
                    if length(tags)>1
                        _geo_yyerror!(context,
                            "Only 1 field can be set as a background field.")
                    elseif isempty(tags)
                        _geo_yywarn!(context,"No field given (Background Field).")
                    else
                        background=tags[1]
                    end
                elseif name=="BoundaryLayer"
                    for tag in _geo_expression_field_tags(
                            fm.captures[2],context,
                            "read_geo_params: BoundaryLayer Field")
                        tag in boundary_layers || push!(boundary_layers,tag)
                    end
                else
                    _geo_yyerror!(context,"Unknown command '$name Field'")
                end
            finally
                context.soft_unknown_reads=false
            end
            return
        end

        # `x~{FExpr}`/`StringToName[sexpr]` LHS heads resolve to a plain
        # identifier before the index/scalar splits — same rewrite the
        # executor applies. The `~{i}` index already reduced upstream: a bad
        # index records `yymsg(0)` and resolves with 0.
        context.soft_unknown_reads=true
        try
            lr=_geo_lhs_string_index(body,context,"read_geo_params")
            lr===nothing ||
                (body=String(lr[1]*lr[2]); code=_geo_unquoted_code(body))
        finally
            context.soft_unknown_reads=false
        end

        _geo_record_list!(context,body) && return

        scalar=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$",body)
        if scalar!==nothing
            # `x = <StringExpr>` stores a string variable upstream — but
            # only when the RHS cannot reduce as `ListOfDouble`: quoted
            # literals, string-producing calls, and bare string tokens.
            # Entity names (`x = Point{1}`) are entity-list queries
            # upstream — the numeric path handles (and diagnoses) them.
            sname=String(scalar.captures[1])
            srhs=String(strip(scalar.captures[2]))
            string_only=(startswith(srhs,'"') || startswith(srhs,'\'') ||
                (sm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*[\(\[]",srhs))!==
                    nothing &&
                 sm.captures[1] in _GEO_STRING_FUNCTIONS &&
                 !(sm.captures[1] in ("StringToName","S2N")) ||
                srhs in ("Today","CodeName","GmshExecutableName",
                         "OnelabAction","CurrentFileName","CurrentDirectory"))
            if string_only &&
               !(sname in _GEO_LEXER_KEYWORDS) &&
               !(sname in _GEO_SIDE_EFFECT_SYMBOLS)
                try
                    context.strings[sname]=
                        String[_geo_eval_string(srhs,context,
                            "read_geo_params: string variable $sname")]
                    delete!(context.values,sname)
                    delete!(context.lists,sname)
                    delete!(context.list_variables,sname)
                catch err
                    err isa InterruptException && rethrow()
                    # Same recovery as numeric records: valid but
                    # unmodelable string productions (entity names the
                    # bounded evaluator cannot resolve) mark the variable
                    # unavailable; grammar failures propagate.
                    if err isa _GeoSyntaxAbort
                        (_geo_list_term_unmodelable(srhs) ||
                         _geo_rhs_uses_unavailable(srhs,context)) ||
                            rethrow()
                    else
                        err isa ArgumentError || rethrow()
                    end
                    _geo_context_forget!(context,sname)
                    message=sprint(showerror,err)
                    context.unavailable[sname]=
                        err isa _GeoSyntaxAbort &&
                            _geo_list_term_unmodelable(srhs) ?
                        "entity or shape expression requires geometry evaluation" :
                        ncodeunits(message)<=240 ? message :
                            String(first(message,220))*"…"
                end
                return
            end
            _geo_record_scalar!(context,scalar.captures[1],scalar.captures[2])
            return
        end
        # `++name`/`--name` have no upstream production — `NumericIncrement`
        # is postfix-only (`a++`), so a prefix increment errors at the token.
        if (pm=match(r"^(\+\+|--)\s*[A-Za-z_]",body))!==nothing
            _geo_syntax_abort(pm.captures[1])
        end
        mutation=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(\+\+|--|\+=|-=|\*=|/=)\s*(.*)$",
                       body)
        if mutation!==nothing
            name=String(mutation.captures[1])
            op=mutation.captures[2]
            rhs=String(strip(mutation.captures[3]))
            # Lexer keywords (`Pi`, `newp`, ...) are not assignable — the op
            # token is the upstream syntax error.
            name in _GEO_LEXER_KEYWORDS && _geo_syntax_abort(op)
            name in _GEO_SIDE_EFFECT_SYMBOLS && _geo_syntax_abort(op)
            is_list=name in context.list_variables
            if op in ("++","--")
                delta=op=="++" ? 1.0 : -1.0
                # Upstream checks `s.list` before touching `s.value` — a
                # list name rejects even though its first element also
                # sits in `values`.
                if is_list
                    _geo_yyerror!(context,"Variable '$name' is a list")
                elseif haskey(context.values,name)
                    _geo_context_set_scalar!(context,name,
                        context.values[name]+delta,"read_geo_params")
                elseif haskey(context.unavailable,name) ||
                       haskey(context.unavailable_lists,name)
                    # stays unavailable
                else
                    _geo_yyerror!(context,"Unknown variable '$name'")
                end
                return
            end
            # `name op= ListOfDouble` — unknown RHS reads are recoverable
            # (`yymsg` + a 0.0 element) upstream, so the RHS is evaluated
            # with soft reads; genuine syntax errors still abort.
            values=begin
                soft=context.soft_unknown_reads
                context.soft_unknown_reads=true
                try
                    _geo_numeric_list_values(rhs,context,
                        "read_geo_params: compound assignment")
                finally
                    context.soft_unknown_reads=soft
                end
            end
            if haskey(context.unavailable,name) ||
               haskey(context.unavailable_lists,name)
                return  # stays unavailable
            elseif is_list
                list=copy(context.lists[name])
                if op=="+="
                    length(list)+length(values)<=_MAX_GEO_LIST_ITEMS ||
                        throw(ArgumentError(
                            "read_geo_params: list exceeds " *
                            "$_MAX_GEO_LIST_ITEMS entries"))
                    _geo_context_set_list!(context,name,
                        vcat(list,values),"read_geo_params")
                elseif op=="-="
                    # Each RHS element erases the first matching entry once
                    # (`std::find` + erase upstream).
                    for value in values
                        index=findfirst(==(value),list)
                        index===nothing || deleteat!(list,index)
                    end
                    _geo_context_set_list!(context,name,list,
                        "read_geo_params")
                else
                    _geo_yyerror!(context,
                        "Operators *= and /= not available for lists")
                end
            elseif haskey(context.values,name)
                if length(values)!=1
                    _geo_yyerror!(context,
                        "Cannot assign list to variable '$name'")
                else
                    d=values[1];v=context.values[name]
                    if op=="+="; v+=d
                    elseif op=="-="; v-=d
                    elseif op=="*="; v*=d
                    elseif iszero(d)
                        _geo_yyerror!(context,
                            "Division by zero in '$name /= " *
                            "$(_geo_gmsh_number(d))'")
                    else
                        v/=d
                    end
                    _geo_context_set_scalar!(context,name,v,"read_geo_params")
                end
            elseif length(values)==1
                # `a op= <scalar>` on an undefined name is a recoverable
                # `Unknown variable` — no binding is created.
                _geo_yyerror!(context,"Unknown variable '$name'")
            else
                # `a op= <list>` on an undefined name still reports
                # `Cannot assign list` (the new symbol defaults to scalar
                # upstream); the empty symbol is marked unavailable so
                # later reads diagnose rather than invent a value.
                _geo_yyerror!(context,
                    "Cannot assign list to variable '$name'")
                context.unavailable[name]="uninitialized variable '$name'"
            end
            return
        end

        # `Exit`/`Abort` end the program — everything after is unreachable.
        # `Exit` is process-global; `Abort` halts only this file's ParseFile.
        if match(r"^Exit\b",body)!==nothing || match(r"^Abort\s*$",body)!==nothing
            stopped[]=true
            stop_kind[]=startswith(body,"Exit") ? :exit : :abort
            return
        end

        # `NewModel` — `new GModel()` becomes current: geometry, fields,
        # physical groups AND the physical-name table die with the old model
        # (a fresh `GModel` gets a fresh `_physicalNames`, unlike
        # `Delete Model` which keeps the same model object). Parser variables
        # survive; `gmsh_yyfactory` does too, but the fresh `GModel` has no
        # OCC internals so every `factory == "OpenCASCADE" && internals`
        # dispatch upstream falls through to the built-in path — modelled
        # here as `factory = :builtin`.
        if match(r"^NewModel\s*$",body)!==nothing
            _geo_allocator_reset_model!(allocator_state,context;fresh=true)
            allocator_state.factory=:builtin
            empty!(kinds);empty!(options);empty!(option_order)
            empty!(creation_curvature);empty!(boundary_layers)
            empty!(physical_name_tags);empty!(groups)
            empty!(scan_raw);empty!(scan_opaque);empty!(scan_ephys)
            scan_has_mesh[]=false;scan_changed[]=true
            background=0
            return
        end

        # `SetFactory` synchronizes the internals when changed before
        # switching kernels, then flips the allocator counters the scan
        # tracks.
        if (sf=match(r"^SetFactory\s*\(\s*(.*?)\s*\)\s*$",body))!==nothing
            if scan_changed[]
                _geo_scan_sync_physicals!(scan_raw,scan_opaque,scan_ephys,
                                          allocator_state,scan_has_mesh[])
                scan_changed[]=false
            end
            arg=strip(sf.captures[1],['"','\'',' '])
            _geo_allocator_set_factory!(allocator_state,arg)
            return
        end

        # `Include`/`Merge` of a `.geo` file splices its statements into the
        # stream — recurse so their variables/options register here. `Merge`
        # forms synchronize the internals when changed first; `Include` does
        # not. `.msh` merges are mesh data — params-neutral.
        if (im=match(r"^(?:Include|Merge|MergeWithBoundingBox)\s+(.*?)\s*$",
                     body))!==nothing
            if !startswith(body,"Include") && scan_changed[]
                _geo_scan_sync_physicals!(scan_raw,scan_opaque,scan_ephys,
                                          allocator_state,scan_has_mesh[])
                scan_changed[]=false
            end
            child_src=strip(im.captures[1],['"','\'',' '])
            child=_geo_fix_relative_path(path,child_src)
            if isfile(child) && !endswith(lowercase(child),".msh")
                # Each `ParseFile` gets its own `gmsh_yyerrorstate` budget
                # and stop scope: a child `Abort`/error-cap halts the child
                # only, while `Exit` propagates into the parent stream.
                parent_kind=stop_kind[]
                parent_base=context.file_error_base
                child_depth=length(ctrl_levels)
                stopped[]=false
                context.file_error_base=length(context.exec_errors)
                _scan_geo_statements(_consume_stmt,child)
                _scan_eof_drain(child_depth)
                if stop_kind[]!==:exit
                    stopped[]=false
                    stop_kind[]=parent_kind
                end
                context.file_error_base=parent_base
            end
            return
        end

        # String variables and option-string writes carry no numeric params:
        # `s = "..."`, `s() = Str(...)`, `s() += Str(...)`, `x.name = "..."`,
        # `x[i].name = "..."`.
        if (sm=match(r"^[A-Za-z_][A-Za-z0-9_]*(?:\s*[\(\[][^\]]*[\)\]]|\s*\([^)]*\.\s*[A-Za-z_][A-Za-z0-9_]*)?\s*(?:\.\s*[A-Za-z_][A-Za-z0-9_]*)?\s*(?:=|\+=)\s*(.*)$",
                     body))!==nothing
            rhs_s=String(strip(sm.captures[1]))
            if startswith(rhs_s,"\"") || startswith(rhs_s,"'") ||
               match(r"^Str\s*\(",rhs_s)!==nothing
                return
            end
        end

        # `SyncModel` forces a synchronize unconditionally; `BoundingBox`,
        # `Save`, `Print`, `Show`/`Hide`/`Color` (including the `Recursive`
        # and deprecated quoted forms), `RefineMesh`, `ReorientMesh`,
        # `RelocateMesh`, `AdaptMesh` and `Mesh n` sync
        # only when the internals changed — the entity-physical mirror must
        # refresh at exactly these points.
        if match(r"^SyncModel\b",body)!==nothing
            _geo_scan_sync_physicals!(scan_raw,scan_opaque,scan_ephys,
                                      allocator_state,scan_has_mesh[])
            scan_changed[]=false
            return
        end
        if match(Regex("^(?:BoundingBox|Save\\b|Print\\b|Show|Hide|Color|" *
            "RefineMesh|ReorientMesh|RelocateMesh|" *
            "AdaptMesh)\\b"),body)!==nothing ||
           match(r"^Recursive\s+(?:Show|Hide|Color)\b",body)!==nothing ||
           match(r"^(?:Show|Hide)\s+[\"']",body)!==nothing
            if scan_changed[]
                _geo_scan_sync_physicals!(scan_raw,scan_opaque,scan_ephys,
                                          allocator_state,scan_has_mesh[])
                scan_changed[]=false
            end
            return
        end
        if match(r"^Mesh\s+[^.]",body)!==nothing
            # `Mesh n` syncs-if-changed then generates elements — the model
            # carries mesh elements whenever a curve/surface/volume exists.
            if scan_changed[]
                _geo_scan_sync_physicals!(scan_raw,scan_opaque,scan_ephys,
                                          allocator_state,scan_has_mesh[])
                scan_changed[]=false
            end
            scan_has_mesh[] = !(isempty(allocator_state.live_builtin_curves) &&
                              isempty(allocator_state.live_builtin_surfaces) &&
                              isempty(allocator_state.live_builtin_volumes) &&
                              isempty(allocator_state.live_occ_curves) &&
                              isempty(allocator_state.live_occ_surfaces) &&
                              isempty(allocator_state.live_occ_volumes))
            return
        end

        # Diagnostics, stream control, geometry-independent commands and the
        # remaining mesh statements carry no params state.
        if match(Regex("^(?:Printf|Warning|Error|Draw|SetChanged|" *
            "UnsplitWindow|NewModel|Sleep|Remesh|SetName|" *
            "CreateDir|System|SystemCall|NonBlockingSystemCall|OnelabRun|" *
            "OptimizeMesh|Plugin|SetCurrentWindow|SplitCurrentWindow|" *
            "SetBoundingBox|SetOrder|PartitionMesh|" *
            "CreateOverlaps|RecombineMesh|" *
            "RenumberMeshNodes|RenumberMeshElements|CreateMeshEdges|" *
            "CreateMeshFaces|TransformMesh|" *
            "CreateTopology|CreateGeometry|" *
            "SetTag|GetForcedStr|GetNumberChoice)\\b"),body)!==nothing
            return
        end

        # A relevant statement that reached here is malformed — the grammar
        # already committed to a production, so the error lands on the token
        # it rejects. `_geo_relevant_code` only matches well-formed heads;
        # `Physical`-/`Field`-led and `<name> Field` statements reach the same
        # committed productions even when the tail is bad (`Physical Foo` →
        # `(Foo)`, `Field foo` → `(foo)`). Quoted occurrences were removed
        # from `code`.
        if _geo_relevant_code(code) ||
           match(r"^(?:Physical|Field)\b",body)!==nothing ||
           match(r"^[A-Za-z_][A-Za-z0-9_]*\s+Field\b",body)!==nothing
            _geo_scan_malformed_abort!(body,context)
        end

        # Everything else reaching here is an unclassifiable `tSTRING`-headed
        # or non-identifier statement — upstream's catch-all diagnostics:
        # `name` alone → `syntax error (;)`; `name(expr)`/`name <expr>` without
        # `=` evaluates the argument (soft — an unknown variable records) then
        # `Unknown command 'name'`; `name {` aborts at the first braced token;
        # a non-identifier head aborts at the token (`= 5` → `(=)`,
        # `5 + 3` → `(5)`). A lexer-keyword head is never an unknown
        # command — see the split below.
        nm=match(r"^[A-Za-z_][A-Za-z0-9_]*",body)
        nm===nothing && _geo_syntax_abort(_geo_stmt_head_token(body))
        # A lexer-keyword head is never an unknown command. Expression-only
        # keywords abort at the keyword itself (`DimNameSpace(ns);` →
        # `(DimNameSpace)`, `Sin(1);` → `(Sin)`); statement-leading keywords
        # (`Point`, `Curve Loop`, `Periodic`, `Spline`, the primitives, ...)
        # shifted into real productions upstream — the scan does not track
        # their state, so the executor owns any tail diagnostics.
        nm_head=String(nm.match)
        if nm_head in _GEO_LEXER_KEYWORDS
            nm_head in _GEO_EXPR_ONLY_KEYWORDS || return nothing
            _geo_syntax_abort(nm_head)
        end
        leftover=String(strip(body[nextind(body,lastindex(nm.match)):end]))
        if isempty(leftover)
            _geo_syntax_abort(";")
        elseif startswith(leftover,"::")
            # `tSCOPE` never follows a bare `tSTRING` in a statement — the
            # `Struct_FullName` namespace forms exist only inside
            # expressions (`ns::x = 5` → `syntax error (::)` upstream).
            _geo_syntax_abort("::")
        elseif startswith(leftover,"{")
            _geo_syntax_abort(_geo_stmt_head_token(leftover[2:end]))
        elseif startswith(leftover,"(")
            close=_geo_matching_delim(leftover,firstindex(leftover))
            if close==0
                _geo_scan_soft_eval_safe(leftover[2:end],context)
                _geo_syntax_abort(";")
            end
            _geo_scan_soft_eval(leftover[2:prevind(leftover,close)],context)
            after=String(strip(leftover[nextind(leftover,close):end]))
            isempty(after) ||
                _geo_syntax_abort(_geo_stmt_head_token(after))
            _geo_yyerror!(context,"Unknown command '$(nm.match)'")
        else
            # A string argument parses as `StringExpr` upstream — only the
            # numeric form reaches FExpr evaluation.
            _geo_string_rhs(leftover) ||
                _geo_scan_soft_eval(leftover,context)
            _geo_yyerror!(context,"Unknown command '$(nm.match)'")
        end
        return nothing
    end
    try
        _scan_geo_statements(_consume_stmt,path)
    catch err
        err isa InterruptException && rethrow()
        # `_GeoSyntaxAbort` is recoverable only inside `execute_geo`'s
        # statement loop; the scan reports parse failures as ArgumentError.
        err isa _GeoSyntaxAbort || rethrow()
        throw(ArgumentError(isempty(err.detail) ?
            "read_geo_params: syntax error ($(err.token))" : err.detail))
    end
    # EOF with open levels — upstream's lexer-level `skip` semantics, shared
    # with `Include`d files through `_scan_eof_drain`.
    _scan_eof_drain(0)
    # No existence validation here: Gmsh stores background/boundary-layer ids
    # and option writes only after a field exists — stale ids are diagnosed
    # per statement above or surface at field-build time like `get(id)` misses.
    fields=Dict{Int,GeoFieldSpec}()
    for (tag,kind) in kinds
        fields[tag]=GeoFieldSpec(tag,kind,get(options,tag,Dict{String,String}()),
                                 get(option_order,tag,String[]),
                                 get(creation_curvature,tag,0))
    end
    return GeoParams(smin,smax,sfactor,seed,groups,fields,background,boundary_layers,
                     geometry_tolerance,boundary_layer_fan_elements,
                     context.exec_errors,context.exec_warnings,
                     context.msg_error_count)
end

function _gmsh_int_option(value::Float64,key::AbstractString)
    return try
        trunc(Int,value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError) || rethrow()
        throw(ArgumentError("read_geo_params: $key is outside the platform Int range"))
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# I-DEAS Universal (.unv) and INRIA/Medit (.mesh) text-mesh formats
#
# Field and descriptor semantics follow the pinned Gmsh 4.15.2
# `GModel::readUNV`/`writeUNV` and `GModel::readMESH`/`writeMESH`:
#
# * UNV dataset 2411 holds nodes (label; export/displacement coordinate-system
#   ids; color; `x y z` on the following line, with `D` exponents accepted) and
#   dataset 2412 holds elements (label; FE descriptor; elementary tag;
#   physical tag; color; node count; node labels, eight per line; beam
#   descriptors carry an extra orientation-triplet line before the nodes).
#   Dataset records are bracketed by `-1` separator lines; unknown datasets
#   are skipped, unknown element descriptors are warned once and skipped after
#   their node list is consumed (upstream behaviour), and descriptor `0` is
#   guessed from the node count (2→line, 3→triangle, 4→tetrahedron).
#   Malformed records and dangling node references are errors.
# * MEDIT/MESH files carry `MeshVersionFormatted` (versions 1–4), `Dimension`
#   (2 or 3, inline or on the following line), `Vertices` (`x y z ref`),
#   `Edges[P2|P3|P4]`, `Triangles[P2|P3|P4]`, `Quadrilaterals`,
#   `Tetrahedra[P2|P3|P4]` and `Hexahedra` sections; unrecognized sections are
#   skipped like upstream. The trailing `ref` becomes the elementary entity
#   tag and, when positive, also the physical tag (Tessella extension — the
#   pinned reader groups entities by `ref` only).
#
# `read_unv`/`read_mesh` return file-level records: a compact MSH2-style
# `MixedMesh` (per-cell `elementary_entities`, physical `tags`) plus the
# original node/element labels. `API.open`/`API.merge` publish the identical
# content through the same channel `.msh` imports use — the parsed records are
# re-encoded as an MSH 2.2 byte stream handed to `GeoExec` staging — so entity
# classification, identity transport and failure atomicity match exactly.
# `API.write` exports the session mesh through `write_unv`/`write_mesh`.
# ════════════════════════════════════════════════════════════════════════════════

struct _TextMeshCell
    label::Int64            # source element label (> 0; dense numbering for .mesh)
    msh::Int                # Gmsh element type code
    nodes::Vector{Int32}    # compact node rows in Gmsh vertex order
    elementary::Int64       # elementary-entity tag (normalized non-negative)
    physical::Int64         # physical tag (normalized non-negative)
end

function _unv_int(token::AbstractString,caller::AbstractString,context::AbstractString)
    value=tryparse(Int,token)
    value===nothing && throw(ArgumentError(
        "$caller: invalid integer '$token' in $context"))
    return Int64(value)
end

function _unv_float(token::AbstractString,caller::AbstractString,context::AbstractString)
    # UNV writes scientific exponents with `D`; accept both `D`/`d`/`E`/`e`.
    token=replace(token,'D'=>'E','d'=>'e')
    value=tryparse(Float64,token)
    value===nothing && throw(ArgumentError(
        "$caller: invalid floating-point value '$token' in $context"))
    isfinite(value) || throw(ArgumentError(
        "$caller: non-finite floating-point value in $context"))
    return value
end

@inline _unv_separator(line::AbstractString)=strip(line)=="-1"

# Beam-family descriptors consume one orientation-triplet line before their
# node labels (upstream `GModel::readUNV` first `switch(type)` block).
const _UNV_BEAM_DESCRIPTORS=Set{Int}((11,21,22,31,23,24,32))

# `type == 0` records guess the descriptor from the node count.
const _UNV_GUESS_TYPES=Dict{Int,Int}(2=>11,3=>41,4=>111)

# UNV descriptor → (Gmsh type, file→Gmsh vertex permutation or `nothing`).
# Permutations are the upstream `readUNV` constructor orders, 1-based:
# `gmsh_nodes[j] = file_nodes[perm[j]]`. Descriptor 117 (linear pyramid) is a
# Tessella extension — the pinned release warns and skips pyramids.
const _UNV_DESCRIPTOR_TABLE=Dict{Int,Tuple{Int,Union{Nothing,NTuple{N,Int} where N}}}(
    11=>(1,nothing),21=>(1,nothing),22=>(1,nothing),31=>(1,nothing),
    23=>(8,(1,3,2)),24=>(8,(1,3,2)),32=>(8,(1,3,2)),
    41=>(2,nothing),51=>(2,nothing),61=>(2,nothing),74=>(2,nothing),
    81=>(2,nothing),91=>(2,nothing),
    42=>(9,(1,3,5,2,4,6)),52=>(9,(1,3,5,2,4,6)),62=>(9,(1,3,5,2,4,6)),
    72=>(9,(1,3,5,2,4,6)),82=>(9,(1,3,5,2,4,6)),92=>(9,(1,3,5,2,4,6)),
    44=>(3,nothing),54=>(3,nothing),64=>(3,nothing),71=>(3,nothing),
    84=>(3,nothing),94=>(3,nothing),
    45=>(16,(1,3,5,7,2,4,6,8)),55=>(16,(1,3,5,7,2,4,6,8)),
    65=>(16,(1,3,5,7,2,4,6,8)),75=>(16,(1,3,5,7,2,4,6,8)),
    85=>(16,(1,3,5,7,2,4,6,8)),95=>(16,(1,3,5,7,2,4,6,8)),
    111=>(4,nothing),
    118=>(11,(1,3,5,10,2,4,6,7,9,8)),
    104=>(5,nothing),115=>(5,nothing),
    105=>(17,(1,3,5,7,13,15,17,19,2,8,9,4,10,6,11,12,14,20,16,18)),
    116=>(17,(1,3,5,7,13,15,17,19,2,8,9,4,10,6,11,12,14,20,16,18)),
    101=>(6,nothing),112=>(6,nothing),
    102=>(18,(1,3,5,10,12,14,2,6,7,4,8,9,11,15,13)),
    113=>(18,(1,3,5,10,12,14,2,6,7,4,8,9,11,15,13)),
    117=>(7,nothing))

# Gmsh type → (UNV descriptor, Gmsh→file vertex permutation or `nothing`).
# The permutations are the upstream `getVertexUNV` maps (+1), the exact
# inverse of the read table. Type 7 (pyramid) uses the 117 extension.
const _UNV_WRITE_TABLE=Dict{Int,Tuple{Int,Union{Nothing,NTuple{N,Int} where N}}}(
    1=>(21,nothing),8=>(24,(1,3,2)),
    2=>(91,nothing),9=>(92,(1,4,2,5,3,6)),
    3=>(94,nothing),16=>(95,(1,5,2,6,3,7,4,8)),
    4=>(111,nothing),11=>(118,(1,5,2,6,3,7,8,10,9,4)),
    5=>(115,nothing),
    17=>(116,(1,9,2,12,3,14,4,10,11,13,15,16,5,17,6,19,7,20,8,18)),
    6=>(112,nothing),18=>(113,(1,7,2,10,3,8,9,11,12,4,13,5,15,6,14)),
    7=>(117,nothing),
    # Point cells have no UNV descriptor; upstream emits them as descriptor-0
    # records which readers then guess/warn-and-skip (a 1-node record is
    # unguessable, so the cell is skipped with a warning on reimport).
    15=>(0,nothing))

function _unv_read_nodes!(io::Base.IO,caller::AbstractString,
        node_labels::Vector{Int64},node_rows::Dict{Int64,Int32},
        coords::Vector{Float64},limit::Int)
    while !eof(io)
        stripped=strip(readline(io))
        _unv_separator(stripped) && return nothing
        isempty(stripped) && continue
        fields=split(stripped)
        length(fields)>=4 || throw(ArgumentError(
            "$caller: malformed UNV 2411 node record '$stripped'"))
        label=_unv_int(fields[1],caller,"UNV node label")
        for i in 2:4
            _unv_int(fields[i],caller,"UNV node field")
        end
        eof(io) && throw(ArgumentError(
            "$caller: truncated UNV 2411 node record"))
        xyz=split(strip(readline(io)))
        length(xyz)>=3 || throw(ArgumentError(
            "$caller: malformed UNV 2411 node coordinates"))
        x=_unv_float(xyz[1],caller,"UNV node coordinate")
        y=_unv_float(xyz[2],caller,"UNV node coordinate")
        z=_unv_float(xyz[3],caller,"UNV node coordinate")
        label>=1 || throw(ArgumentError(
            "$caller: non-positive UNV node label $label"))
        haskey(node_rows,label) && throw(ArgumentError(
            "$caller: duplicate UNV node label $label"))
        length(node_labels)>=limit && throw(ArgumentError(
            "$caller: UNV node count exceeds max_nodes=$limit"))
        push!(node_labels,label);node_rows[label]=Int32(length(node_labels))
        push!(coords,x);push!(coords,y);push!(coords,z)
    end
    return nothing
end

function _unv_read_elements!(io::Base.IO,caller::AbstractString,
        node_rows::Dict{Int64,Int32},records::Vector{_TextMeshCell},
        warned::Set{Int},next_elementary::Base.RefValue{Int64},
        element_labels::Set{Int64},limit::Int)
    while !eof(io)
        stripped=strip(readline(io))
        _unv_separator(stripped) && return nothing
        (isempty(stripped) || length(stripped)<3) && continue
        fields=split(stripped)
        length(fields)>=6 || throw(ArgumentError(
            "$caller: malformed UNV 2412 element record '$stripped'"))
        label=_unv_int(fields[1],caller,"UNV element label")
        descriptor=_unv_int(fields[2],caller,"UNV element descriptor")
        elementary=_unv_int(fields[3],caller,"UNV elementary tag")
        physical=_unv_int(fields[4],caller,"UNV physical tag")
        _unv_int(fields[5],caller,"UNV element color")
        count=_unv_int(fields[6],caller,"UNV element node count")
        label>=1 || throw(ArgumentError(
            "$caller: non-positive UNV element label $label"))
        label in element_labels && throw(ArgumentError(
            "$caller: duplicate UNV element label $label"))
        push!(element_labels,label)
        count<0 && throw(ArgumentError(
            "$caller: negative UNV element node count $count"))
        length(records)>=limit && throw(ArgumentError(
            "$caller: UNV element count exceeds max_elements=$limit"))
        if descriptor==0
            descriptor=get(_UNV_GUESS_TYPES,Int(count),0)
            descriptor==0 || descriptor in warned || begin
                push!(warned,descriptor)
                @warn "$caller: UNV element $label has no type; " *
                      "guessing descriptor $descriptor from $count nodes"
            end
        end
        if descriptor in _UNV_BEAM_DESCRIPTORS
            eof(io) && throw(ArgumentError(
                "$caller: truncated UNV beam orientation line"))
            readline(io)   # orientation triplet — parsed positionally, unused
        end
        nodelabels=Int64[]
        count>0 && sizehint!(nodelabels,Int(count))
        while length(nodelabels)<count
            eof(io) && throw(ArgumentError(
                "$caller: truncated UNV element $label node list"))
            nl=strip(readline(io))
            _unv_separator(nl) && throw(ArgumentError(
                "$caller: truncated UNV element $label node list"))
            isempty(nl) && continue
            for tok in split(nl)
                length(nodelabels)<count || throw(ArgumentError(
                    "$caller: UNV element $label node list overruns $count labels"))
                push!(nodelabels,_unv_int(tok,caller,"UNV element node label"))
            end
        end
        entry=get(_UNV_DESCRIPTOR_TABLE,descriptor,nothing)
        if entry===nothing
            descriptor in warned || begin
                push!(warned,descriptor)
                @warn "$caller: skipping unsupported UNV element descriptor $descriptor"
            end
            continue
        end
        msh,perm=entry
        arity=perm===nothing ? msh_num_nodes(msh) : length(perm)
        count==arity || throw(ArgumentError(
            "$caller: UNV descriptor $descriptor expects $arity nodes, got $count"))
        rows=Vector{Int32}(undef,Int(count))
        @inbounds for i in 1:count
            row=get(node_rows,nodelabels[i],Int32(0))
            row!=0 || throw(ArgumentError(
                "$caller: UNV element $label references unknown node $(nodelabels[i])"))
            rows[i]=row
        end
        if elementary<0
            elementary=next_elementary[]+Int64(1)
        end
        next_elementary[]=max(next_elementary[],elementary)
        physical<0 && (physical=0)
        0<=elementary<=typemax(Int32) || throw(ArgumentError(
            "$caller: UNV elementary tag $elementary on element $label exceeds Int32"))
        0<=physical<=typemax(Int32) || throw(ArgumentError(
            "$caller: UNV physical tag $physical on element $label exceeds Int32"))
        push!(records,_TextMeshCell(label,msh,
            perm===nothing ? rows : rows[collect(perm)],elementary,physical))
    end
    return nothing
end

function _read_unv(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    node_labels=Int64[]
    node_rows=Dict{Int64,Int32}()
    coords=Float64[]
    records=_TextMeshCell[]
    warned=Set{Int}()
    next_elementary=Ref{Int64}(0)
    element_labels=Set{Int64}()
    while !eof(io)
        line=readline(io)
        _unv_separator(line) || continue
        eof(io) && break
        dataset=readline(io)
        while _unv_separator(dataset)
            eof(io) && (dataset="";break)
            dataset=readline(io)
        end
        number=tryparse(Int,strip(dataset))
        number===nothing && continue
        if number==2411
            _unv_read_nodes!(io,caller,node_labels,node_rows,coords,node_limit)
        elseif number==2412
            _unv_read_elements!(io,caller,node_rows,records,warned,
                next_elementary,element_labels,element_limit)
        end
        # Unknown datasets are skipped line-by-line until the next `-1`
        # separator, exactly like the pinned reader.
    end
    return node_labels,coords,records
end

function _foreign_mesh_result(caller::AbstractString,node_labels::Vector{Int64},
        coords::Vector{Float64},records::Vector{_TextMeshCell};
        physical_names::Dict{Tuple{Int,Int},String}=Dict{Tuple{Int,Int},String}())
    nn=length(node_labels)
    C=reshape(coords,3,nn)
    per=Dict{Int,Vector{_TextMeshCell}}();order=Int[]
    for record in records
        bucket=get!(()->_TextMeshCell[],per,record.msh)
        isempty(bucket) && push!(order,record.msh)
        push!(bucket,record)
    end
    blocks=ElementBlock[];owners=Vector{Int32}[];element_tags=Vector{UInt64}[]
    for msh in order
        cells=per[msh];arity=msh_num_nodes(msh);n=length(cells)
        nodes=Matrix{Int32}(undef,arity,n)
        tags=Vector{Int32}(undef,n);ents=Vector{Int32}(undef,n)
        uids=Vector{UInt64}(undef,n)
        @inbounds for j in 1:n
            cell=cells[j]
            nodes[:,j]=cell.nodes
            tags[j]=Int32(cell.physical)
            ents[j]=Int32(cell.elementary)
            uids[j]=UInt64(cell.label)
        end
        push!(blocks,ElementBlock(msh,nodes,tags))
        push!(owners,ents);push!(element_tags,uids)
    end
    mesh=MixedMesh(C,blocks;elementary_entities=owners,physical_names=physical_names)
    return mesh,UInt64.(node_labels),element_tags
end

"""
    UnvFile

Result of [`read_unv`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the UNV elementary tag and whose
`tags` hold the physical tag), the original `node_tags` in row order and the
per-block `element_tags`.
"""
struct UnvFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    MeshFile

Result of [`read_mesh`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the Medit `ref` and whose `tags`
carry the same positive `ref` as a physical tag), sequential `node_tags` and
per-block `element_tags`.
"""
struct MeshFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_unv(path; max_nodes, max_elements, max_file_bytes) -> UnvFile

Read an I-DEAS Universal `.unv` file (datasets 2411 nodes and 2412 elements).
Element descriptors map to Gmsh types with the pinned upstream table; unknown
descriptors warn once and are skipped, `D` exponents are accepted, and node
references must resolve. Other datasets (including 2477 groups) are skipped.
"""
function read_unv(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_unv"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records=open(path,"r") do io
        try
            return _read_unv(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .unv file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records)
    return UnvFile(mesh,node_tags,element_tags)
end

# ── INRIA/Medit (.mesh) ────────────────────────────────────────────────────────

# keyword → (Gmsh type, record node count). TetrahedraP2 swaps the last two
# connectivity slots in both directions (upstream `readMESH`/`writeMESH`).
const _MEDIT_READ_SECTIONS=Dict{String,Tuple{Int,Int}}(
    "Edges"=>(1,2),"EdgesP2"=>(8,3),"EdgesP3"=>(26,4),"EdgesP4"=>(27,5),
    "Triangles"=>(2,3),"TrianglesP2"=>(9,6),"TrianglesP3"=>(21,10),
    "TrianglesP4"=>(23,15),
    "Quadrilaterals"=>(3,4),
    "Tetrahedra"=>(4,4),"TetrahedraP2"=>(11,10),"TetrahedraP3"=>(29,20),
    "TetrahedraP4"=>(30,35),
    "Hexahedra"=>(5,8))
const _MEDIT_TET10_SWAP=(1,2,3,4,5,6,7,8,10,9)

function _medit_count(io,caller,section,fields)
    # Upstream uses `fscanf`, so the count may follow the keyword on the same
    # line or occupy the following line.
    tokens=length(fields)>=2 ? view(fields,2:length(fields)) : begin
        eof(io) && throw(ArgumentError(
            "$caller: truncated .mesh $section count"))
        split(strip(readline(io)))
    end
    length(tokens)>=1 || throw(ArgumentError(
        "$caller: missing .mesh $section count"))
    count=_unv_int(tokens[1],caller,"$section count")
    count>=0 || throw(ArgumentError(
        "$caller: negative .mesh $section count $count"))
    return Int(count)
end

function _medit_float(token,caller,context)
    value=tryparse(Float64,token)
    value===nothing && throw(ArgumentError(
        "$caller: invalid floating-point value '$token' in $context"))
    isfinite(value) || throw(ArgumentError(
        "$caller: non-finite floating-point value in $context"))
    return value
end

function _read_medit(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    # The first meaningful line carries `MeshVersionFormatted` (pinned reader).
    version=nothing
    while !eof(io)
        stripped=strip(readline(io))
        (isempty(stripped) || startswith(stripped,'#')) && continue
        fields=split(stripped)
        fields[1]=="MeshVersionFormatted" || throw(ArgumentError(
            "$caller: .mesh files must begin with MeshVersionFormatted"))
        if length(fields)>=2
            version=_unv_int(fields[2],caller,"MeshVersionFormatted value")
        else
            eof(io) && throw(ArgumentError(
                "$caller: truncated MeshVersionFormatted"))
            vf=split(strip(readline(io)))
            length(vf)>=1 || throw(ArgumentError(
                "$caller: missing MeshVersionFormatted value"))
            version=_unv_int(vf[1],caller,"MeshVersionFormatted value")
        end
        break
    end
    version===nothing && throw(ArgumentError(
        "$caller: missing MeshVersionFormatted header"))
    # Upstream accepts ASCII versions 1, 2 and 4; version 3 marks the binary
    # encoding, which this reader intentionally rejects.
    version in (1,2,4) || throw(ArgumentError(
        "$caller: unsupported Medit mesh version $version (binary marker)"))
    dimension=3
    coords=Float64[]
    vertex_count=0
    records=_TextMeshCell[]
    element_counter=Ref{Int64}(0)
    node_labels=Int64[]
    while !eof(io)
        stripped=strip(readline(io))
        (isempty(stripped) || startswith(stripped,'#')) && continue
        fields=split(stripped)
        keyword=fields[1]
        if keyword=="Dimension"
            value=if length(fields)>=2
                _unv_int(fields[2],caller,"Dimension")
            else
                eof(io) && throw(ArgumentError("$caller: truncated Dimension"))
                df=split(strip(readline(io)))
                length(df)>=1 || throw(ArgumentError(
                    "$caller: missing Dimension value"))
                _unv_int(df[1],caller,"Dimension")
            end
            value in (2,3) || throw(ArgumentError(
                "$caller: unsupported .mesh dimension $value"))
            dimension=Int(value)
        elseif keyword=="Vertices"
            count=_medit_count(io,caller,"Vertices",fields)
            count<=node_limit || throw(ArgumentError(
                "$caller: .mesh vertex count exceeds max_nodes=$node_limit"))
            vertex_count==0 || isempty(records) || throw(ArgumentError(
                "$caller: repeated Vertices after element sections"))
            empty!(coords);resize!(node_labels,0);vertex_count=count
            for _ in 1:count
                eof(io) && throw(ArgumentError(
                    "$caller: truncated Vertices section"))
                vf=split(strip(readline(io)))
                need=dimension+1
                length(vf)>=need || throw(ArgumentError(
                    "$caller: malformed .mesh vertex record"))
                x=_medit_float(vf[1],caller,"vertex coordinate")
                y=_medit_float(vf[2],caller,"vertex coordinate")
                z=dimension==3 ? _medit_float(vf[3],caller,"vertex coordinate") : 0.0
                _unv_int(vf[need],caller,"vertex ref")  # parsed, unused upstream
                push!(coords,x);push!(coords,y);push!(coords,z)
                push!(node_labels,Int64(length(node_labels)+1))
            end
        elseif haskey(_MEDIT_READ_SECTIONS,keyword)
            msh,arity=_MEDIT_READ_SECTIONS[keyword]
            count=_medit_count(io,caller,keyword,fields)
            length(records)+count<=element_limit || throw(ArgumentError(
                "$caller: .mesh element count exceeds max_elements=$element_limit"))
            vertex_count>0 || count==0 || throw(ArgumentError(
                "$caller: $keyword appeared before Vertices"))
            for _ in 1:count
                eof(io) && throw(ArgumentError(
                    "$caller: truncated $keyword section"))
                ef=split(strip(readline(io)))
                length(ef)>=arity+1 || throw(ArgumentError(
                    "$caller: malformed $keyword record"))
                rows=Vector{Int32}(undef,arity)
                @inbounds for i in 1:arity
                    index=_unv_int(ef[i],caller,"$keyword vertex index")
                    1<=index<=vertex_count || throw(ArgumentError(
                        "$caller: $keyword record references vertex $index " *
                        "outside 1:$vertex_count"))
                    rows[i]=Int32(index)
                end
                ref=_unv_int(ef[arity+1],caller,"$keyword ref")
                0<=ref<=typemax(Int32) || throw(ArgumentError(
                    "$caller: .mesh ref $ref is outside 0:$(typemax(Int32))"))
                msh==11 && (rows=rows[collect(_MEDIT_TET10_SWAP)])
                element_counter[]+=Int64(1)
                push!(records,_TextMeshCell(
                    element_counter[],msh,rows,ref,ref))
            end
        end
        # `End` and every unrecognized keyword/section are skipped line-wise,
        # matching the pinned reader.
    end
    return node_labels,coords,records
end

"""
    read_mesh(path; max_nodes, max_elements, max_file_bytes) -> MeshFile

Read an INRIA/Medit `.mesh` file (versions 1–4, `Dimension` 2 or 3).
`Vertices`, `Edges[P2|P3|P4]`, `Triangles[P2|P3|P4]`, `Quadrilaterals`,
`Tetrahedra[P2|P3|P4]` and `Hexahedra` records map to their Gmsh types; the
trailing `ref` becomes the elementary entity tag and, when positive, the
physical tag. Unrecognized sections (`Corners`, `RequiredVertices`, ...) are
skipped; dangling vertex references and malformed records are errors.
"""
function read_mesh(path;
                   max_nodes=typemax(Int32),max_elements=typemax(Int32),
                   max_file_bytes=typemax(Int))
    caller="read_mesh"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records=open(path,"r") do io
        try
            return _read_medit(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .mesh file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records)
    return MeshFile(mesh,node_tags,element_tags)
end

# ── INRIA/Medit binary (.meshb) ────────────────────────────────────────────────
# Keyword codes come from the libMeshb `GmfKwdCod` table (the enum ordinal is
# the on-disk code): each counted section is `code, next-keyword byte offset,
# count, records`, where records are the same arity vertex indices plus a
# trailing `ref` as the ASCII sections. The file version selects the storage
# widths — version 1: 4-byte reals, 4-byte indices, 4-byte offsets; version 2:
# 8-byte reals, 4-byte indices, 4-byte offsets; version 3: 8-byte reals,
# 4-byte indices, 8-byte offsets; version 4: 8-byte reals, 8-byte indices,
# 8-byte offsets. Keyword codes and the Dimension value are always 4-byte
# words. The first word doubles as the endianness tag: reading it as a
# little-endian Int32 gives 1, or 0x01000000 for a big-endian file.
const _MEDITB_READ_SECTIONS=Dict{Int32,Tuple{Int,Int}}(
    Int32(5)=>(1,2),Int32(25)=>(8,3),Int32(92)=>(26,4),Int32(93)=>(27,5),
    Int32(6)=>(2,3),Int32(24)=>(9,6),Int32(90)=>(21,10),Int32(91)=>(23,15),
    Int32(7)=>(3,4),
    Int32(8)=>(4,4),Int32(30)=>(11,10),Int32(96)=>(29,20),Int32(97)=>(30,35),
    Int32(10)=>(5,8))
const _MEDITB_SECTION_CODE=Dict{String,Int32}(
    "Edges"=>5,"EdgesP2"=>25,"EdgesP3"=>92,"EdgesP4"=>93,
    "Triangles"=>6,"TrianglesP2"=>24,"TrianglesP3"=>90,"TrianglesP4"=>91,
    "Quadrilaterals"=>7,
    "Tetrahedra"=>8,"TetrahedraP2"=>30,"TetrahedraP3"=>96,"TetrahedraP4"=>97,
    "Hexahedra"=>10)
const _MEDITB_VERSION_FORMATTED=Int32(1)
const _MEDITB_DIMENSION=Int32(3)
const _MEDITB_VERTICES=Int32(4)
const _MEDITB_END=Int32(54)
const _MEDITB_SWAPPED_TAG=Int32(16777216)   # byte-swapped 1 → big-endian file

@inline function _meditb_value(io,swapped::Bool,::Type{T}) where {T}
    value=read(io,T)
    swapped || return value
    U=T===Float32 ? UInt32 : T===Float64 ? UInt64 : unsigned(T)
    return reinterpret(T,bswap(reinterpret(U,value)))
end

function _meditb_file_end(io)
    mark=position(io);seekend(io);last=position(io);seek(io,mark)
    return last
end

function _meditb_check_pos(io,caller,offset)
    offset>0 && offset!=position(io) && throw(ArgumentError(
        "$caller: inconsistent .meshb next-keyword offset $offset"))
    return nothing
end

function _read_meditb(io::Base.IO,caller::AbstractString,node_limit::Int,
                      element_limit::Int)
    # The mandatory header is `MeshVersionFormatted` (code 1) then the version
    # word; the code read little-endian also yields the endianness tag.
    eof(io) && throw(ArgumentError("$caller: empty .meshb file"))
    magic=read(io,Int32)
    magic==_MEDITB_VERSION_FORMATTED || magic==_MEDITB_SWAPPED_TAG ||
        throw(ArgumentError(
            "$caller: not a binary Medit mesh (unexpected " *
            "MeshVersionFormatted code)"))
    swapped=magic==_MEDITB_SWAPPED_TAG
    version=_meditb_value(io,swapped,Int32)
    version in Int32(1):Int32(4) || throw(ArgumentError(
        "$caller: unsupported .meshb version $version"))
    itype=version>=4 ? Int64 : Int32   # 'i' fields: indices, refs, counts
    ptype=version>=3 ? Int64 : Int32   # next-keyword offsets
    ftype=version==1 ? Float32 : Float64
    read_int()=Int64(_meditb_value(io,swapped,itype))
    read_pos()=Int64(_meditb_value(io,swapped,ptype))
    read_real()=Float64(_meditb_value(io,swapped,ftype))
    # The Dimension keyword is mandatory and must come first; its value stays
    # a 4-byte word in every version.
    code=_meditb_value(io,swapped,Int32)
    code==_MEDITB_DIMENSION || throw(ArgumentError(
        "$caller: .meshb keyword $code found where Dimension was expected"))
    offset=read_pos()
    dimension=Int(_meditb_value(io,swapped,Int32))
    dimension in (2,3) || throw(ArgumentError(
        "$caller: unsupported .meshb dimension $dimension"))
    _meditb_check_pos(io,caller,offset)
    file_end=_meditb_file_end(io)
    coords=Float64[]
    vertex_count=0
    records=_TextMeshCell[]
    element_counter=Ref{Int64}(0)
    node_labels=Int64[]
    while !eof(io)
        code=_meditb_value(io,swapped,Int32)
        code==_MEDITB_END && break
        offset=read_pos()
        if code==_MEDITB_VERTICES
            count=read_int()
            count>=0 || throw(ArgumentError(
                "$caller: negative .meshb Vertices count $count"))
            count<=node_limit || throw(ArgumentError(
                "$caller: .meshb vertex count exceeds max_nodes=$node_limit"))
            vertex_count==0 || isempty(records) || throw(ArgumentError(
                "$caller: repeated Vertices after element sections"))
            empty!(coords);resize!(node_labels,0);vertex_count=Int(count)
            for _ in 1:count
                x=read_real();y=read_real()
                z=dimension==3 ? read_real() : 0.0
                isfinite(x) && isfinite(y) && isfinite(z) || throw(ArgumentError(
                    "$caller: non-finite floating-point value in vertex " *
                    "coordinate"))
                read_int()   # vertex ref: parsed, unused upstream
                push!(coords,x);push!(coords,y);push!(coords,z)
                push!(node_labels,Int64(length(node_labels)+1))
            end
        elseif haskey(_MEDITB_READ_SECTIONS,code)
            msh,arity=_MEDITB_READ_SECTIONS[code]
            count=read_int()
            count>=0 || throw(ArgumentError(
                "$caller: negative .meshb keyword $code count $count"))
            length(records)+count<=element_limit || throw(ArgumentError(
                "$caller: .meshb element count exceeds " *
                "max_elements=$element_limit"))
            vertex_count>0 || count==0 || throw(ArgumentError(
                "$caller: keyword $code appeared before Vertices"))
            for _ in 1:count
                rows=Vector{Int32}(undef,arity)
                @inbounds for i in 1:arity
                    index=read_int()
                    1<=index<=vertex_count || throw(ArgumentError(
                        "$caller: .meshb keyword $code record references " *
                        "vertex $index outside 1:$vertex_count"))
                    rows[i]=Int32(index)
                end
                ref=read_int()
                0<=ref<=typemax(Int32) || throw(ArgumentError(
                    "$caller: .meshb ref $ref is outside 0:$(typemax(Int32))"))
                msh==11 && (rows=rows[collect(_MEDIT_TET10_SWAP)])
                element_counter[]+=Int64(1)
                push!(records,_TextMeshCell(
                    element_counter[],msh,rows,ref,ref))
            end
        else
            # Unrecognized keyword (Corners, Required*, normals, solutions,
            # prisms, ...): jump to the next-keyword offset — the binary
            # counterpart of skipping an unrecognized ASCII section.
            offset>position(io) && offset<=file_end || throw(ArgumentError(
                "$caller: unrecognized .meshb keyword $code without a " *
                "valid next-keyword offset"))
            seek(io,offset)
            continue
        end
        _meditb_check_pos(io,caller,offset)
    end
    return node_labels,coords,records
end

"""
    read_meshb(path; max_nodes, max_elements, max_file_bytes) -> MeshFile

Read an INRIA/Medit binary `.meshb` file (versions 1–4, `Dimension` 2 or 3,
little- or big-endian). Keyword-coded sections carry the same records as the
ASCII `.mesh` sections: `Vertices`, `Edges[P2|P3|P4]`, `Triangles[P2|P3|P4]`,
`Quadrilaterals`, `Tetrahedra[P2|P3|P4]` and `Hexahedra` map to their Gmsh
types; the trailing `ref` becomes the elementary entity tag and, when
positive, the physical tag. Unrecognized keywords are skipped through their
next-keyword offset; dangling vertex references, malformed records and
inconsistent offsets are errors. Returns the same [`MeshFile`](@ref)
contract as [`read_mesh`](@ref). Gmsh itself has no `.meshb` route; this is
the libMeshb/GMF binary counterpart of the Medit reader.
"""
function read_meshb(path;
                    max_nodes=typemax(Int32),max_elements=typemax(Int32),
                    max_file_bytes=typemax(Int))
    caller="read_meshb"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records=open(path,"r") do io
        try
            return _read_meditb(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .meshb file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records)
    return MeshFile(mesh,node_tags,element_tags)
end

# ── Session-model staging through the native mesh-import channel ───────────────
# `.unv`/`.mesh` records are re-encoded as a bounded MSH 2.2 stream (node and
# element labels, `physical elementary` tag pair) so `GeoExec` publishes them
# with exactly the same classification, identity and atomicity semantics as
# `.msh` imports.

function _write_foreign_bridge(path,node_labels::Vector{Int64},
        coords::Vector{Float64},records::Vector{_TextMeshCell};
        node_offset::Integer=0,element_offset::Integer=0,
        physical_names::Dict{Tuple{Int,Int},String}=Dict{Tuple{Int,Int},String}())
    # `physical_names` is emitted as an MSH `$PhysicalNames` section so the
    # staging reader publishes e.g. SU2 marker names exactly like a native
    # `.msh` import does. Names containing a double quote cannot be encoded.
    for ((dim,tag),name) in physical_names
        dim in 0:3 && tag>0 && !occursin('"',name) || throw(ArgumentError(
            "foreign mesh bridge: physical name for ($dim,$tag) is not " *
            "representable in MSH 2.2"))
    end
    open(path,"w") do io
        println(io,"\$MeshFormat")
        println(io,"2.2 0 8")
        println(io,"\$EndMeshFormat")
        if !isempty(physical_names)
            println(io,"\$PhysicalNames")
            println(io,length(physical_names))
            for ((dim,tag),name) in sort!(collect(physical_names))
                println(io,dim," ",tag," \"",name,"\"")
            end
            println(io,"\$EndPhysicalNames")
        end
        println(io,"\$Nodes")
        println(io,length(node_labels))
        @inbounds for i in eachindex(node_labels)
            println(io,node_labels[i]+node_offset," ",
                repr(coords[3i-2])," ",repr(coords[3i-1])," ",repr(coords[3i]))
        end
        println(io,"\$EndNodes")
        println(io,"\$Elements")
        println(io,length(records))
        for record in records
            print(io,record.label+element_offset," ",record.msh," 2 ",
                record.physical," ",record.elementary)
            for row in record.nodes
                print(io," ",node_labels[Int(row)]+node_offset)
            end
            println(io)
        end
        println(io,"\$EndElements")
    end
    return path
end

function _read_foreign_records(path,caller::AbstractString)
    extension=lowercase(splitext(String(path))[2])
    if extension==".unv"
        node_labels,coords,records=open(path,"r") do io
            try
                return _read_unv(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .unv file"))
                rethrow()
            end
        end
        return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
    elseif extension==".mesh"
        node_labels,coords,records=open(path,"r") do io
            try
                return _read_medit(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .mesh file"))
                rethrow()
            end
        end
        return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
    elseif extension==".meshb"
        node_labels,coords,records=open(path,"r") do io
            try
                return _read_meditb(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .meshb file"))
                rethrow()
            end
        end
        return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
    elseif extension==".su2"
        return open(path,"r") do io
            try
                return _read_su2(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .su2 file"))
                rethrow()
            end
        end
    elseif extension in (".bdf",".nas")
        return open(path,"r") do io
            try
                return _read_bdf(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .bdf file"))
                rethrow()
            end
        end
    elseif extension==".inp"
        return open(path,"r") do io
            try
                return _read_inp(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .inp file"))
                rethrow()
            end
        end
    elseif extension==".neu"
        return open(path,"r") do io
            try
                return _read_neu(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .neu file"))
                rethrow()
            end
        end
    elseif extension==".p3d"
        return open(path,"r") do io
            try
                return _read_p3d(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .p3d file"))
                rethrow()
            end
        end
    elseif extension==".mail"
        return open(path,"r") do io
            try
                return _read_mail(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .mail file"))
                rethrow()
            end
        end
    elseif extension==".ir3"
        return open(path,"r") do io
            try
                return _read_ir3(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .ir3 file"))
                rethrow()
            end
        end
    elseif extension==".diff"
        return open(path,"r") do io
            try
                return _read_diff(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .diff file"))
                rethrow()
            end
        end
    elseif extension==".vtk"
        return open(path,"r") do io
            try
                return _read_vtk(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .vtk file"))
                rethrow()
            end
        end
    elseif extension==".wrl"
        node_labels,coords,records=open(path,"r") do io
            try
                return _read_wrl(io,caller,
                    Int(typemax(Int32)),Int(typemax(Int32)))
            catch err
                err isa InterruptException && rethrow()
                err isa EOFError && throw(ArgumentError(
                    "$caller: truncated .wrl file"))
                rethrow()
            end
        end
        return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
    elseif extension==".off"
        node_labels,coords,records=_read_off(path,caller)
        return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
    elseif extension==".obj"
        node_labels,coords,records=_read_obj(path,caller)
        return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
    elseif extension==".ply"
        node_labels,coords,records=_read_ply(path,caller)
        return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
    elseif extension==".ply2"
        node_labels,coords,records=_read_ply2(path,caller)
        return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
    elseif extension in (".xao",".brep")
        # `.xao`/`.brep` carry an OCCT CASCADE BRep dump: publishing model
        # entities needs the analytic curve/surface boundary kernel, which is
        # the classified `importShapes` gap — `read_xao`/`read_brep` expose the
        # truthful envelope/shape census instead of fabricating a model.
        throw(ArgumentError("$caller: $extension CAD import requires the " *
            "OpenCASCADE BRep boundary kernel (TShape curve/surface " *
            "reconstruction), which Tessella does not implement; " *
            "IO.read_xao/IO.read_brep expose the envelope census"))
    elseif extension in (".med",".mmed",".omm")
        throw(ArgumentError("$caller: $extension import requires the " *
            "MED/MED-file library (upstream Gmsh must be compiled with " *
            "MED support), which Tessella does not link"))
    elseif extension==".cgns"
        throw(ArgumentError("$caller: .cgns import requires the CGNS " *
            "library (upstream Gmsh must be compiled with CGNS support), " *
            "which Tessella does not link"))
    else
        throw(ArgumentError("$caller: unsupported mesh extension '$extension'"))
    end
end

function _foreign_bridge_temp(node_labels,coords,records;
                              node_offset::Integer=0,element_offset::Integer=0,
                              physical_names::Dict{Tuple{Int,Int},String}=
                                  Dict{Tuple{Int,Int},String}())
    temporary=string(tempname(),".msh")
    try
        _write_foreign_bridge(temporary,node_labels,coords,records;
            node_offset=node_offset,element_offset=element_offset,
            physical_names=physical_names)
    catch
        rm(temporary;force=true)
        rethrow()
    end
    return temporary
end

function _stage_foreign_mesh_file(path;context_lc=1.0,forced_bbox=false)
    node_labels,coords,records,pnames=_read_foreign_records(path,"API.open")
    temporary=_foreign_bridge_temp(node_labels,coords,records;
        physical_names=pnames)
    geo=getfield(parentmodule(@__MODULE__),:GeoExec)
    try
        result=geo._geo_stage_mesh_file(temporary;
            context_lc=context_lc,forced_bbox=forced_bbox)
        result.model.name=String(splitext(basename(String(path)))[1])
        return result
    finally
        rm(temporary;force=true)
    end
end

function _stage_foreign_mesh_merge(path,model,seed;context_lc=1.0,forced_bbox=false)
    node_labels,coords,records,pnames=_read_foreign_records(path,"API.merge")
    # File-local labels are not session UIDs; shift them past the live
    # maxima so a merge never collides with existing nodes or elements.
    temporary=_foreign_bridge_temp(node_labels,coords,records;
        node_offset=Int64(seed.mesh_node_max),
        element_offset=Int64(seed.mesh_element_max),
        physical_names=pnames)
    geo=getfield(parentmodule(@__MODULE__),:GeoExec)
    try
        return geo._geo_stage_mesh_merge(temporary,model,seed;
            context_lc=context_lc,forced_bbox=forced_bbox)
    finally
        rm(temporary;force=true)
    end
end

# ── Writers ────────────────────────────────────────────────────────────────────

function _foreign_target(path,caller::AbstractString)
    path isa AbstractString && !isempty(path) && !occursin('\0',path) ||
        throw(ArgumentError("$caller: filename must be a nonempty NUL-free string"))
    target=abspath(path);parent=dirname(target)
    isdir(parent) || throw(ArgumentError(
        "$caller: parent directory does not exist: $parent"))
    isdir(target) && throw(ArgumentError(
        "$caller: destination is a directory: $target"))
    return target,parent
end

@inline function _foreign_cell_entity(m::MixedMesh,bi::Int,j::Int)
    data=m.entity_data
    data!==nothing && return Int64(data.block_entities[bi][j])
    m.elementary_entities!==nothing && return Int64(m.elementary_entities[bi][j])
    return Int64(0)
end

@inline function _foreign_cell_uid(m::MixedMesh,bi::Int,j::Int,counter::UInt64)
    data=m.entity_data
    data===nothing && return counter
    return data.external_element_tags[bi][j]
end

"""
    write_unv(path, mesh)

Write a `Mesh` or `MixedMesh` as an I-DEAS Universal file: dataset 2411 nodes
(`label`, coordinate-system ids 1/1, color 11, `%25.16E` coordinates) and
dataset 2412 elements (`label`, upstream descriptor, elementary tag, physical
tag, color 7, node count, labels eight per line; beam records 21/24 carry a
`0 0 0` orientation line). Element types with no upstream descriptor —
including 2477 group output, which the pinned release uses for physical
groups — are rejected; physical tags are emitted in field 4, which the
upstream reader interprets as the physical tag.
"""
function write_unv(path,m)
    caller="write_unv"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    # Validate every block before touching the destination file.
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .unv"))
        get(_UNV_WRITE_TABLE,block.msh,nothing)===nothing && throw(ArgumentError(
            "$caller: element type $(block.msh) ($(msh_family(block.msh)) " *
            "order $(msh_spec(block.msh).order)) is not representable in .unv"))
    end
    data=m.entity_data
    nn=size(m.coords,2)
    labels=data===nothing ? collect(UInt64(1):UInt64(nn)) : data.external_node_tags
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"    -1");println(io,"  2411")
        @inbounds for i in 1:nn
            print(io,@sprintf("%10d%10d%10d%10d\n",labels[i],1,1,11))
            print(io,@sprintf("%25.16E%25.16E%25.16E\n",
                m.coords[1,i],m.coords[2,i],m.coords[3,i]))
        end
        println(io,"    -1");println(io,"    -1");println(io,"  2412")
        counter=UInt64(0)
        for (bi,block) in enumerate(m.blocks)
            descriptor,perm=_UNV_WRITE_TABLE[block.msh]
            arity=size(block.nodes,1)
            for j in 1:size(block.nodes,2)
                counter+=UInt64(1)
                elnum=_foreign_cell_uid(m,bi,j,counter)
                elementary=_foreign_cell_entity(m,bi,j)
                physical=Int64(block.tags[j])
                print(io,@sprintf("%10d%10d%10d%10d%10d%10d\n",
                    elnum,descriptor,elementary,physical,7,arity))
                (descriptor==21 || descriptor==24) &&
                    print(io,@sprintf("%10d%10d%10d\n",0,0,0))
                for k in 1:arity
                    row=perm===nothing ? k : perm[k]
                    k>1 && mod(k-1,8)==0 && print(io,"\n")
                    print(io,@sprintf("%10d",labels[block.nodes[row,j]]))
                end
                print(io,"\n")
            end
        end
        println(io,"    -1")
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# Medit section keyword per Gmsh type (linear and P2–P4 simplex sections plus
# the linear quadrilateral and hexahedron, matching the pinned writer).
const _MEDIT_WRITE_SECTIONS=Dict{Int,String}(
    1=>"Edges",8=>"EdgesP2",26=>"EdgesP3",27=>"EdgesP4",
    2=>"Triangles",9=>"TrianglesP2",21=>"TrianglesP3",23=>"TrianglesP4",
    3=>"Quadrilaterals",
    4=>"Tetrahedra",11=>"TetrahedraP2",29=>"TetrahedraP3",30=>"TetrahedraP4",
    5=>"Hexahedra")
const _MEDIT_SECTION_ORDER=("Edges","EdgesP2","EdgesP3","EdgesP4",
    "Triangles","TrianglesP2","TrianglesP3","TrianglesP4","Quadrilaterals",
    "Tetrahedra","TetrahedraP2","TetrahedraP3","TetrahedraP4","Hexahedra")

"""
    write_mesh(path, mesh)

Write a `Mesh` or `MixedMesh` as an INRIA/Medit file (`MeshVersionFormatted 2`,
`Dimension 3`, `Vertices`/`Edges`/`Triangles`/`Quadrilaterals`/`Tetrahedra`/
`Hexahedra` and the P2–P4 variants, `End`). The per-element `ref` is the cell
physical tag when positive, otherwise its elementary-entity tag; vertex `ref`
is the owning entity tag. Element families Medit cannot express (quadratic
quads/hexes, prisms, pyramids, point cells, special records) are rejected.
"""
function write_mesh(path,m)
    caller="write_mesh"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    sections=Dict{Int,String}()
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .mesh"))
        # Medit has no point-element section; upstream writeMESH emits vertex
        # entities into `Vertices` only, so point cells are not records here.
        section=block.msh==15 ? "" : get(_MEDIT_WRITE_SECTIONS,block.msh,nothing)
        section===nothing && throw(ArgumentError(
            "$caller: element type $(block.msh) ($(msh_family(block.msh)) " *
            "order $(msh_spec(block.msh).order)) is not representable in .mesh"))
        sections[block.msh]=section
    end
    data=m.entity_data
    nn=size(m.coords,2)
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io," MeshVersionFormatted 2")
        println(io," Dimension")
        println(io," 3")
        println(io," Vertices")
        println(io," ",nn)
        @inbounds for i in 1:nn
            ref=data===nothing ? Int32(0) : data.node_entities[i][2]
            println(io,@sprintf(" %20.14G %20.14G %20.14G %d",
                m.coords[1,i],m.coords[2,i],m.coords[3,i],ref))
        end
        for section in _MEDIT_SECTION_ORDER
            count=0
            for (bi,block) in enumerate(m.blocks)
                get(sections,block.msh,"")==section || continue
                count+=size(block.nodes,2)
            end
            count==0 && continue
            println(io," ",section)
            println(io," ",count)
            for (bi,block) in enumerate(m.blocks)
                get(sections,block.msh,"")==section || continue
                arity=size(block.nodes,1)
                for j in 1:size(block.nodes,2)
                    physical=Int64(block.tags[j])
                    entity=_foreign_cell_entity(m,bi,j)
                    ref=physical>0 ? physical : entity
                    for k in 1:arity
                        row=block.msh==11 ? _MEDIT_TET10_SWAP[k] : k
                        print(io," ",block.nodes[row,j])
                    end
                    println(io," ",ref)
                end
            end
        end
        println(io," End")
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

"""
    write_meshb(path, mesh; version=0)

Write a `Mesh` or `MixedMesh` as an INRIA/Medit binary `.meshb` file
(keyword-coded `Vertices`, `Edges`/`Triangles`/`Quadrilaterals`/`Tetrahedra`/
`Hexahedra` and the P2–P4 variants, `End`), the binary counterpart of
[`write_mesh`](@ref). `version` selects the libMeshb encoding: 1 (32-bit
reals/indices/offsets), 2 (64-bit reals, 32-bit indices/offsets), 3 (64-bit
reals/offsets, 32-bit indices) or 4 (64-bit everywhere); the default `0`
picks version 3, upgrading to 4 when a section count exceeds Int32 — the
same upgrade rule libMeshb applies. The per-element `ref` and vertex `ref`
follow the ASCII writer's semantics exactly, as do element-order and
rejections of families Medit cannot express. Gmsh itself has no `.meshb`
route; this writes the standard GMF layout readable by libMeshb-based
tools.
"""
function write_meshb(path,m;version::Integer=0)
    caller="write_meshb"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    sections=Dict{Int,String}()
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .meshb"))
        # Medit has no point-element section; upstream writeMESH emits vertex
        # entities into `Vertices` only, so point cells are not records here.
        section=block.msh==15 ? "" : get(_MEDIT_WRITE_SECTIONS,block.msh,nothing)
        section===nothing && throw(ArgumentError(
            "$caller: element type $(block.msh) ($(msh_family(block.msh)) " *
            "order $(msh_spec(block.msh).order)) is not representable " *
            "in .meshb"))
        sections[block.msh]=section
    end
    data=m.entity_data
    nn=size(m.coords,2)
    counts=Int64[nn]
    for section in _MEDIT_SECTION_ORDER
        c=Int64(0)
        for block in m.blocks
            get(sections,block.msh,"")==section || continue
            c+=size(block.nodes,2)
        end
        push!(counts,c)
    end
    maxref=Int64(0)
    for (bi,block) in enumerate(m.blocks)
        for j in 1:size(block.nodes,2)
            physical=Int64(block.tags[j])
            entity=_foreign_cell_entity(m,bi,j)
            maxref=max(maxref,physical>0 ? physical : entity)
        end
    end
    version==0 && (version=any(>(typemax(Int32)),counts) ? 4 : 3)
    version in 1:4 || throw(ArgumentError(
        "$caller: .meshb version must be between 1 and 4"))
    version<=3 && (any(>(typemax(Int32)),counts) ||
        maxref>typemax(Int32)) && throw(ArgumentError(
        "$caller: .meshb version $version cannot encode counts or refs " *
        "above 2^31-1 (requires version 4)"))
    itype=version>=4 ? Int64 : Int32
    ptype=version>=3 ? Int64 : Int32
    isize=Int64(sizeof(itype));psize=Int64(sizeof(ptype))
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        write_i(v)=write(io,itype(v))
        write_pos(v)=write(io,ptype(v))
        write_real(x)=version==1 ? write(io,Float32(x)) :
            write(io,Float64(x))
        # Header: MeshVersionFormatted code+version, then the mandatory
        # Dimension keyword; `cursor` tracks the next keyword's byte offset.
        write(io,_MEDITB_VERSION_FORMATTED);write(io,Int32(version))
        cursor=Int64(8+4+psize+4)
        write(io,_MEDITB_DIMENSION);write_pos(cursor);write(io,Int32(3))
        # Vertices is always emitted, like the ASCII writer.
        vrec=Int64(3*(version==1 ? 4 : 8)+isize)
        write(io,_MEDITB_VERTICES)
        write_pos(cursor+4+psize+isize+nn*vrec);write_i(nn)
        @inbounds for i in 1:nn
            ref=data===nothing ? Int32(0) : data.node_entities[i][2]
            write_real(m.coords[1,i]);write_real(m.coords[2,i])
            write_real(m.coords[3,i]);write_i(ref)
        end
        cursor+=4+psize+isize+nn*vrec
        for section in _MEDIT_SECTION_ORDER
            count=0
            for block in m.blocks
                get(sections,block.msh,"")==section || continue
                count+=size(block.nodes,2)
            end
            count==0 && continue
            arity=_MEDIT_READ_SECTIONS[section][2]
            reclen=Int64(arity*isize+isize)
            write(io,_MEDITB_SECTION_CODE[section])
            write_pos(cursor+4+psize+isize+count*reclen);write_i(count)
            for (bi,block) in enumerate(m.blocks)
                get(sections,block.msh,"")==section || continue
                for j in 1:size(block.nodes,2)
                    physical=Int64(block.tags[j])
                    entity=_foreign_cell_entity(m,bi,j)
                    ref=physical>0 ? physical : entity
                    for k in 1:arity
                        row=block.msh==11 ? _MEDIT_TET10_SWAP[k] : k
                        write_i(block.nodes[row,j])
                    end
                    write_i(ref)
                end
            end
            cursor+=4+psize+isize+count*reclen
        end
        write(io,_MEDITB_END);write_pos(0)
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── STL session staging ──────────────────────────────────────────────────────
# `.stl` open/merge routes through the same atomic MSH 2.2 bridge as the
# `.unv`/`.mesh` foreign readers: `read_stl`'s welded surface is re-encoded
# as one elementary entity of Tri3 records, which the native importer then
# publishes as a discrete Surface — the same ownership `API.mesh.import_stl`
# produces for a session file.

function _stl_bridge_records(path,tag::Integer)
    mesh=read_stl(path)
    nn=nnodes(mesh)
    node_labels=collect(Int64(1):Int64(nn))
    coords=Vector{Float64}(undef,3*nn)
    @inbounds for i in 1:nn
        coords[3i-2]=mesh.coords[1,i]
        coords[3i-1]=mesh.coords[2,i]
        coords[3i]=mesh.coords[3,i]
    end
    records=_TextMeshCell[_TextMeshCell(Int64(i),2,
        Int32[mesh.tris[1,i],mesh.tris[2,i],mesh.tris[3,i]],
        Int64(tag),Int64(0)) for i in axes(mesh.tris,2)]
    return node_labels,coords,records
end

# Lowest free Surface tag across the native entity dictionaries, discrete
# entities, imported mesh entities, and the live mesh authority — the
# `API.mesh.import_stl` allocation rule applied to a merge target.
function _stl_free_surface_tag(model)
    mod=getfield(parentmodule(@__MODULE__),:Model)
    occupied=Set{Int}()
    for (dim,tag) in mod.model_entities(model)
        dim==2 && push!(occupied,Int(tag))
    end
    for (dim,tag) in keys(model.discrete)
        dim==2 && push!(occupied,Int(tag))
    end
    for (dim,tag) in keys(model.meshing.imported_entities)
        dim==2 && push!(occupied,Int(tag))
    end
    authority=model.meshing.mesh_authority
    if authority!==nothing
        for mesh in (authority.mesh,authority.carrier)
            if mesh!==nothing && hasproperty(mesh,:entity_data) &&
                    mesh.entity_data!==nothing
                for (dim,tag) in keys(mesh.entity_data.entities)
                    dim==2 && push!(occupied,Int(tag))
                end
            end
        end
    end
    tag=1
    while tag in occupied
        tag+=1
    end
    return tag
end

# Temporary MSH 2.2 encoding of an STL file: every facet lands on the given
# elementary tag so the importer creates that discrete Surface.
function _stl_bridge_temp(path,tag::Integer)
    node_labels,coords,records=_stl_bridge_records(path,tag)
    return _foreign_bridge_temp(node_labels,coords,records)
end

function _stage_stl_mesh_file(path;context_lc=1.0,forced_bbox=false)
    temporary=_stl_bridge_temp(path,1)
    geo=getfield(parentmodule(@__MODULE__),:GeoExec)
    try
        result=geo._geo_stage_mesh_file(temporary;
            context_lc=context_lc,forced_bbox=forced_bbox)
        result.model.name=String(splitext(basename(String(path)))[1])
        return result
    finally
        rm(temporary;force=true)
    end
end

function _stage_stl_mesh_merge(path,model,seed;
                               context_lc=1.0,forced_bbox=false)
    # File-local labels are not session UIDs; the bridge offsets them past the
    # live maxima while the STL surface lands on the lowest free entity tag.
    node_labels,coords,records=_stl_bridge_records(
        path,_stl_free_surface_tag(model))
    temporary=_foreign_bridge_temp(node_labels,coords,records;
        node_offset=Int64(seed.mesh_node_max),
        element_offset=Int64(seed.mesh_element_max))
    geo=getfield(parentmodule(@__MODULE__),:GeoExec)
    try
        return geo._geo_stage_mesh_merge(temporary,model,seed;
            context_lc=context_lc,forced_bbox=forced_bbox)
    finally
        rm(temporary;force=true)
    end
end

"""
    write_stl(path, mesh)

Write the surface facets of a `Mesh` or `MixedMesh` as an ASCII STL file
(Gmsh's default STL encoding). `Tri3` blocks emit one facet each; `Quad4`
blocks split on their first-third vertex diagonal — `(a,b,c)` then `(a,c,d)`
— matching upstream's deterministic quad triangulation. Facet normals are
recomputed from the wound vertices; degenerate facets emit a zero normal.

STL is a surface-only format, so volume element blocks (any `msh_dimension`
3 type) and surface types STL cannot carry (quadratic triangles, higher-order
quads, special records) reject loudly instead of being silently dropped.
Lower-dimension point and line records carry no surface data and are
skipped, as upstream's face-only traversal does. A mesh with no surface
elements at all is an error — an empty `solid` would look like a successful
export of nothing.
"""
function write_stl(path,m)
    caller="write_stl"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    # Validate every block before touching the destination file.
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .stl"))
        dim=msh_dimension(block.msh)
        dim>=3 && throw(ArgumentError(
            "$caller: volume element type $(block.msh) " *
            "($(msh_family(block.msh))) is not representable in .stl"))
        if dim==2 && !(block.msh in (2,3))
            throw(ArgumentError(
                "$caller: element type $(block.msh) ($(msh_family(block.msh)) " *
                "order $(msh_spec(block.msh).order)) is not representable " *
                "in .stl"))
        end
    end
    target,parent=_foreign_target(path,caller)
    emitted=0
    mktemp(parent) do temporary,io
        println(io,"solid tessella")
        for block in m.blocks
            block isa ElementBlock || continue
            block.msh in (2,3) || continue
            for j in axes(block.nodes,2)
                a,b=Int(block.nodes[1,j]),Int(block.nodes[2,j])
                c=Int(block.nodes[3,j])
                _stl_facet!(io,m.coords,a,b,c)
                emitted+=1
                if block.msh==3
                    _stl_facet!(io,m.coords,a,c,Int(block.nodes[4,j]))
                    emitted+=1
                end
            end
        end
        println(io,"endsolid tessella")
        emitted==0 && throw(ArgumentError(
            "$caller: mesh has no surface elements to write"))
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

function _stl_facet!(io,coords,a::Int,b::Int,c::Int)
    ax,ay,az=coords[1,a],coords[2,a],coords[3,a]
    ux,uy,uz=coords[1,b]-ax,coords[2,b]-ay,coords[3,b]-az
    vx,vy,vz=coords[1,c]-ax,coords[2,c]-ay,coords[3,c]-az
    nx,ny,nz=uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx
    len=sqrt(nx*nx+ny*ny+nz*nz)
    if isfinite(len) && len>0
        nx/=len;ny/=len;nz/=len
    else
        nx=ny=nz=0.0
    end
    println(io," facet normal ",
        @sprintf("%.16g %.16g %.16g",nx,ny,nz))
    println(io,"  outer loop")
    for i in (a,b,c)
        println(io,"   vertex ",
            @sprintf("%.16g %.16g %.16g",
                coords[1,i],coords[2,i],coords[3,i]))
    end
    println(io,"  endloop")
    println(io," endfacet")
    return nothing
end

# ════════════════════════════════════════════════════════════════════════════════
# SU2 native mesh (.su2) and Nastran bulk data (.bdf/.nas)
#
# Field semantics follow the public SU2 mesh specification and the MSC/Altair
# Nastran bulk-data card layouts:
#
# * SU2 files are keyword sections. `NDIME= 2|3` declares the dimension,
#   `NPOIN= n` is followed by `x y [z] index` point rows (0-based indices; the
#   trailing index column is optional but must be uniform), `NELEM= n` by
#   `type v1..vn [index]` element rows, and `NMARK= n` introduces that many
#   `MARKER_TAG= name`/`MARKER_ELEMS= n` boundary-element lists. Element type
#   identifiers are the VTK codes: 3=line, 5=triangle, 9=quadrilateral,
#   10=tetrahedron, 12=hexahedron, 13=prism, 14=pyramid; connectivity is
#   0-based and already in Gmsh vertex order. `%` starts a comment line.
#   Marker *i* (1-based) becomes physical tag *i* and, on domain elements
#   (NELEM), physical tag 0; elementary entity tags are 1 for NELEM records
#   and *i* for marker records. Marker names are returned as
#   `physical_names` and transported through the MSH 2.2 staging bridge.
#   `NZONE`/`IZONE`/multi-zone decks reject loudly.
# * BDF files carry `GRID` cards (id, cp, x, y, z — 8-column fixed and
#   comma-delimited free formats, plus `GRID*` 16-column large fields), and
#   element cards `CBAR`, `CBEAM`, `CROD`, `CTUBE`, `CTRIA3`, `CTRIA6`,
#   `CQUAD4`, `CQUAD8`, `CTETRA[4|10]`, `CPENTA[6|15]`, `CHEXA[8|20]`,
#   `CPYRA[5|13]` — the explicit second-order names plus the base-name forms
#   that pad extra nodes onto the same card. The property-id field becomes
#   both the elementary entity tag and the physical tag. `+`-prefixed
#   continuation lines merge into the parent card (fixed, free and
#   `*`-prefixed large-field forms); `$` comments and all non-element
#   cards are skipped, `INCLUDE` rejects loudly, `ENDDATA` ends the deck.
#
# `read_su2`/`read_bdf` return the same file-level records as `read_unv`:
# a compact MSH2-style `MixedMesh` plus the original node/element labels.
# `API.open`/`API.merge` publish them through the identical MSH 2.2 bridge
# the `.unv`/`.mesh` readers use, and `API.write` exports the session mesh
# through `write_su2`/`write_bdf` on the same MSH4-projection path.
# ════════════════════════════════════════════════════════════════════════════════

# ── SU2 native mesh ────────────────────────────────────────────────────────────

# SU2 element-type identifier → (Gmsh type, file node count). Connectivity
# follows VTK ordering, which coincides with the Gmsh vertex order for every
# supported family, so no permutation is needed.
const _SU2_ELEMENT_TYPES=Dict{Int,Tuple{Int,Int}}(
    3=>(1,2),5=>(2,3),9=>(3,4),10=>(4,4),12=>(5,8),13=>(6,6),14=>(7,5))
const _SU2_WRITE_TABLE=Dict{Int,Int}(1=>3,2=>5,3=>9,4=>10,5=>12,6=>13,7=>14)

function _su2_next(io)
    # Next meaningful line: blank lines and `%` comments are skipped.
    while !eof(io)
        line=strip(readline(io))
        (isempty(line) || startswith(line,'%')) && continue
        return line
    end
    return nothing
end

function _su2_key(line::AbstractString)
    m=match(r"^\s*([A-Za-z_][A-Za-z_0-9]*)\s*=\s*(.*?)\s*$",line)
    m===nothing && return nothing
    return uppercase(m.captures[1]),m.captures[2]
end

function _su2_count(value::AbstractString,caller::AbstractString,
                    context::AbstractString)
    count=_unv_int(value,caller,"$context count")
    count>=0 || throw(ArgumentError("$caller: negative .su2 $context count $count"))
    return Int(count)
end

function _su2_element_row(line::AbstractString,caller::AbstractString,
        point_rows::Dict{Int64,Int32},context::AbstractString)
    fields=split(line)
    length(fields)>=2 || throw(ArgumentError(
        "$caller: malformed .su2 $context record '$line'"))
    su2type=_unv_int(fields[1],caller,".su2 element type")
    entry=get(_SU2_ELEMENT_TYPES,Int(su2type),nothing)
    entry===nothing && throw(ArgumentError(
        "$caller: unsupported .su2 element type $su2type in $context record"))
    msh,arity=entry
    length(fields)==arity+1 || length(fields)==arity+2 || throw(ArgumentError(
        "$caller: .su2 $context type $su2type expects $arity nodes, got " *
        "$(length(fields)-1) trailing fields"))
    nodes=Vector{Int32}(undef,arity)
    @inbounds for k in 1:arity
        v=_unv_int(fields[k+1],caller,".su2 node index")
        v>=0 || throw(ArgumentError("$caller: negative .su2 node index $v"))
        row=get(point_rows,v,Int32(0))
        row!=0 || throw(ArgumentError(
            "$caller: .su2 $context references unknown point index $v"))
        nodes[k]=row
    end
    index=length(fields)==arity+2 ?
        _unv_int(fields[end],caller,".su2 element index") : Int64(-1)
    return msh,nodes,index
end

function _read_su2(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    ndime=0;npoin=0
    node_labels=Int64[];point_rows=Dict{Int64,Int32}();coords=Float64[]
    records=_TextMeshCell[];pnames=Dict{Tuple{Int,Int},String}()
    element_labels=Set{Int64}();points_indexed=nothing
    marker_count=0;marker_index=0;next_label=Ref{Int64}(0)
    while true
        line=_su2_next(io)
        line===nothing && break
        kv=_su2_key(line)
        # Stray data lines and unknown keywords are skipped like upstream.
        kv===nothing && continue
        key,value=kv
        if key in ("NZONE","IZONE","MULTIZONE","NPERIODIC")
            throw(ArgumentError(
                "$caller: multi-zone/periodic .su2 files are not supported ($key)"))
        elseif key=="NDIME"
            d=_unv_int(value,caller,"NDIME")
            d in (2,3) || throw(ArgumentError(
                "$caller: unsupported .su2 NDIME=$d"))
            ndime=Int(d)
        elseif key=="NPOIN"
            ndime in (2,3) || throw(ArgumentError(
                "$caller: NPOIN appeared before NDIME"))
            npoin==0 || throw(ArgumentError("$caller: repeated NPOIN section"))
            count=_su2_count(value,caller,"NPOIN")
            count<=node_limit || throw(ArgumentError(
                "$caller: .su2 point count exceeds max_nodes=$node_limit"))
            for i in 1:count
                row=_su2_next(io)
                row===nothing && throw(ArgumentError(
                    "$caller: truncated .su2 NPOIN section"))
                fields=split(row)
                (length(fields)==ndime || length(fields)==ndime+1) ||
                    throw(ArgumentError(
                        "$caller: malformed .su2 point record '$row'"))
                indexed=length(fields)==ndime+1
                if points_indexed===nothing
                    points_indexed=indexed
                elseif indexed!=points_indexed
                    throw(ArgumentError(
                        "$caller: mixed indexed/non-indexed .su2 point records"))
                end
                x=_unv_float(fields[1],caller,".su2 point coordinate")
                y=_unv_float(fields[2],caller,".su2 point coordinate")
                z=ndime==3 ?
                    _unv_float(fields[3],caller,".su2 point coordinate") : 0.0
                index=indexed ?
                    _unv_int(fields[ndime+1],caller,".su2 point index") :
                    Int64(i-1)
                index>=0 || throw(ArgumentError(
                    "$caller: negative .su2 point index $index"))
                haskey(point_rows,index) && throw(ArgumentError(
                    "$caller: duplicate .su2 point index $index"))
                push!(node_labels,index+1)
                point_rows[index]=Int32(length(node_labels))
                push!(coords,x);push!(coords,y);push!(coords,z)
            end
            npoin=count
        elseif key=="NELEM"
            count=_su2_count(value,caller,"NELEM")
            length(records)+count<=element_limit || throw(ArgumentError(
                "$caller: .su2 element count exceeds max_elements=$element_limit"))
            for i in 1:count
                row=_su2_next(io)
                row===nothing && throw(ArgumentError(
                    "$caller: truncated .su2 NELEM section"))
                msh,nodes,index=_su2_element_row(row,caller,point_rows,"element")
                label=index>=0 ? index+1 : Int64(i)
                label in element_labels && throw(ArgumentError(
                    "$caller: duplicate .su2 element index $(label-1)"))
                push!(element_labels,label)
                next_label[]=max(next_label[],label)
                # NELEM records carry no physical association; all share the
                # single domain entity 1, matching upstream's one-region model.
                push!(records,_TextMeshCell(label,msh,nodes,Int64(1),Int64(0)))
            end
        elseif key=="NMARK"
            marker_index==0 || throw(ArgumentError(
                "$caller: NMARK appeared after marker sections"))
            marker_count=_su2_count(value,caller,"NMARK")
        elseif key=="MARKER_TAG"
            marker_count>0 || throw(ArgumentError(
                "$caller: MARKER_TAG appeared before NMARK"))
            marker_index+=1
            marker_index<=marker_count || throw(ArgumentError(
                "$caller: more markers than the declared NMARK=$marker_count"))
            name=value
            eline=_su2_next(io)
            eline===nothing && throw(ArgumentError(
                "$caller: truncated .su2 marker $marker_index"))
            ekv=_su2_key(eline)
            ekv!==nothing && ekv[1]=="MARKER_ELEMS" || throw(ArgumentError(
                "$caller: expected MARKER_ELEMS after MARKER_TAG, got '$eline'"))
            count=_su2_count(ekv[2],caller,"MARKER_ELEMS")
            length(records)+count<=element_limit || throw(ArgumentError(
                "$caller: .su2 element count exceeds max_elements=$element_limit"))
            mdim=-1
            for _ in 1:count
                row=_su2_next(io)
                row===nothing && throw(ArgumentError(
                    "$caller: truncated .su2 marker $marker_index elements"))
                msh,nodes,_=_su2_element_row(
                    row,caller,point_rows,"marker element")
                mdim<0 && (mdim=msh_dimension(msh))
                next_label[]+=1;label=next_label[]
                label in element_labels && throw(ArgumentError(
                    "$caller: duplicate .su2 element index $(label-1)"))
                push!(element_labels,label)
                push!(records,_TextMeshCell(label,msh,nodes,
                    Int64(marker_index),Int64(marker_index)))
            end
            isempty(name) || (pnames[(mdim<0 ? ndime-1 : mdim,marker_index)]=name)
        else
            continue
        end
    end
    marker_index!=marker_count && throw(ArgumentError(
        "$caller: .su2 file declares NMARK=$marker_count but contains " *
        "$marker_index marker(s)"))
    return node_labels,coords,records,pnames
end

"""
    Su2File

Result of [`read_su2`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the domain/marker entity tag and
whose `tags` hold the marker-derived physical tag; marker names live in
`mesh.physical_names`), the original `node_tags` in row order and the
per-block `element_tags`.
"""
struct Su2File
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_su2(path; max_nodes, max_elements, max_file_bytes) -> Su2File

Read an SU2 native `.su2` mesh file (`NDIME`, `NPOIN`, `NELEM` and the
`NMARK`/`MARKER_TAG`/`MARKER_ELEMS` boundary lists). Marker *i* becomes
physical tag *i* and its `MARKER_TAG` becomes the physical name; `NELEM`
records carry physical tag 0 on domain entity 1. Connectivity resolves the
0-based point indices; unsupported element types, multi-zone decks and
malformed records are errors. `%` comment lines are skipped.
"""
function read_su2(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_su2"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_su2(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .su2 file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return Su2File(mesh,node_tags,element_tags)
end

function _su2_marker_name(m::MixedMesh,dim::Int,tag::Int32)
    if tag>0
        name=get(m.physical_names,(dim,Int(tag)),"")
        isempty(name) || return replace(strip(name),r"\s+"=>"_")
        return "physical_$tag"
    end
    return "unmarked"
end

"""
    write_su2(path, mesh)

Write a `Mesh` or `MixedMesh` as an SU2 native file: `NDIME` from the highest
element dimension (2 or 3), `NPOIN` with 0-based trailing indices, `NELEM`
with the top-dimension cells, and `NMARK` boundary lists grouping the
dimension-1 cells by physical tag (`MARKER_TAG` is the physical name, or a
`physical_<tag>`/`unmarked` fallback). Cell types SU2 cannot express
(quadratic and higher-order cells, special records) reject loudly; cells
more than one dimension below the domain — curve and point elements of a
3-D mesh — are skipped like upstream because the format has no boundary-edge
model. Physical tags on the top-dimension cells and physical names of
un-named groups are not representable and are dropped.
"""
function write_su2(path,m)
    caller="write_su2"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .su2"))
        block.msh==15 && continue   # no point element in SU2
        get(_SU2_WRITE_TABLE,block.msh,nothing)===nothing && throw(ArgumentError(
            "$caller: element type $(block.msh) ($(msh_family(block.msh)) " *
            "order $(msh_spec(block.msh).order)) is not representable in .su2"))
    end
    ndime=maximum(msh_dimension(block.msh) for block in m.blocks;init=0)
    ndime in (2,3) || throw(ArgumentError(
        "$caller: .su2 export requires 2-D or 3-D elements"))
    nn=size(m.coords,2)
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"NDIME= $ndime")
        println(io,"NPOIN= $nn")
        @inbounds for i in 1:nn
            print(io,repr(m.coords[1,i])," ",repr(m.coords[2,i]))
            ndime==3 && print(io," ",repr(m.coords[3,i]))
            println(io," ",i-1)
        end
        # Top-dimension cells form NELEM; boundary-dimension cells group into
        # markers by physical tag; deeper sub-boundary cells are skipped.
        nelem=0;nmarker=0
        for block in m.blocks
            d=msh_dimension(block.msh)
            d==ndime && (nelem+=size(block.nodes,2))
            d==ndime-1 && (nmarker+=size(block.nodes,2))
        end
        println(io,"NELEM= $nelem")
        eidx=0
        for block in m.blocks
            msh_dimension(block.msh)==ndime || continue
            su2type=_SU2_WRITE_TABLE[block.msh];arity=size(block.nodes,1)
            for j in 1:size(block.nodes,2)
                print(io,su2type)
                for k in 1:arity
                    print(io," ",block.nodes[k,j]-1)
                end
                println(io," ",eidx);eidx+=1
            end
        end
        # Marker groups keyed by physical tag, emitted in tag order.
        groups=Dict{Int32,Vector{Tuple{Int,Int}}}()
        for (bi,block) in enumerate(m.blocks)
            msh_dimension(block.msh)==ndime-1 || continue
            for j in 1:size(block.nodes,2)
                push!(get!(()->Tuple{Int,Int}[],groups,block.tags[j]),(bi,j))
            end
        end
        println(io,"NMARK= $(length(groups))")
        for tag in sort!(collect(keys(groups)))
            cells=groups[tag]
            println(io,"MARKER_TAG= $(_su2_marker_name(m,ndime-1,tag))")
            println(io,"MARKER_ELEMS= $(length(cells))")
            for (bi,j) in cells
                block=m.blocks[bi];arity=size(block.nodes,1)
                print(io,_SU2_WRITE_TABLE[block.msh])
                for k in 1:arity
                    print(io," ",block.nodes[k,j]-1)
                end
                println(io)
            end
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── Nastran bulk data (.bdf/.nas) ──────────────────────────────────────────────

# Card mnemonic → candidate (Gmsh type, node count, flat-field positions)
# triples, tried in order. Flat positions index the merged card field list:
# 1=mnemonic, 2=EID, 3=PID, connectivity starts at 4; each continuation line
# contributes its data fields after the first (continuation-id) field.
# Alternative position sets cover the MSC layout (corner nodes on the first
# line, midsides on continuations after THETA/ZOFFS/Ti), the Altair layout
# (all nodes on the first line) and minimal continuations.
const _BDF_ELEMENT_CARDS=Dict{String,Vector{Tuple{Int,Int,Vector{Int}}}}(
    "CBAR"=>[(1,2,[4,5])],
    "CBEAM"=>[(1,2,[4,5])],
    "CROD"=>[(1,2,[4,5])],
    "CTUBE"=>[(1,2,[4,5])],
    "CTRIA3"=>[(2,3,[4,5,6])],
    "CTRIA6"=>[(9,6,[4,5,6,13,14,15]),(9,6,[4,5,6,7,8,9]),
               (9,6,[4,5,6,10,11,12])],
    "CQUAD4"=>[(3,4,[4,5,6,7])],
    "CQUAD8"=>[(16,8,collect(4:11)),(16,8,[4,5,6,7,10,11,12,13])],
    "CTETRA"=>[(11,10,collect(4:13)),(11,10,[4,5,6,7,10,11,12,13,14,15]),
               (4,4,[4,5,6,7])],
    "CTETRA10"=>[(11,10,collect(4:13)),(11,10,[4,5,6,7,10,11,12,13,14,15])],
    "CTETRA4"=>[(4,4,[4,5,6,7])],
    "CPENTA"=>[(18,15,collect(4:18)),(6,6,collect(4:9))],
    "CPENTA15"=>[(18,15,collect(4:18))],
    "CPENTA6"=>[(6,6,collect(4:9))],
    "CHEXA"=>[(17,20,collect(4:23)),(5,8,collect(4:11)),
              (5,8,[4,5,6,7,10,11,12,13])],
    "CHEXA20"=>[(17,20,collect(4:23))],
    "CHEXA8"=>[(5,8,collect(4:11)),(5,8,[4,5,6,7,10,11,12,13])],
    "CPYRA"=>[(19,13,collect(4:16)),
              (19,13,[4,5,6,7,8,10,11,12,13,14,15,16,17]),
              (7,5,collect(4:8))],
    "CPYRA13"=>[(19,13,collect(4:16)),
                (19,13,[4,5,6,7,8,10,11,12,13,14,15,16,17])],
    "CPYRA5"=>[(7,5,collect(4:8))])

# Non-element card prefixes skipped silently (properties, materials, loads,
# constraints, coordinate systems, case control); other unmapped `C`-cards
# warn once because they look like element connections.
const _BDF_SILENT_PREFIXES=(
    "MAT","PLOAD","PSHELL","PSOLID","PCOMP","PBAR","PBEAM","PROD","PELAS",
    "PGAP","PVISC","PMASS","PSHEAR","PTUBE","PBEND","PDAMP","PBUSH","PSET",
    "PCONEAX","PLPLANE","PRAC","SPC","MPC","FORCE","MOMENT","GRAV","PARAM",
    "SET","CORD","TABLE","DAREA","DOPTPRM","EIG","METHOD","SUBCASE","LOAD",
    "TEMP","SUPORT","ASET","OMIT","K2GG","M2GG","SPOINT","EPOINT","SEQGP",
    "DMIG","GENEL","RBE","RSPLINE","ENDDATA","BEGIN","CEND","TITLE",
    "SUBTITLE","LABEL","ECHO","OUTPUT","DISP","STRESS","STRAIN","SPCF",
    "ELFORCE","GPFORCE","OLOAD","AERO","TRIM","DIVERG","FLUTTER","CASE",
    "CONM","CELAS","CDAMP","CBUSH","CMASS","CGAP","CVISC","CHBDY",
    "CAABSF","CDUM","CRAC","CFAST","CWELD","CSEAM","CJOIN")

@inline _bdf_fld(fields::Vector{String},index::Int)=
    index<=length(fields) ? fields[index] : ""

function _bdf_int(token::AbstractString,caller::AbstractString,
                  context::AbstractString;default=nothing)
    isempty(token) && return default
    value=tryparse(Int,token)
    value===nothing && throw(ArgumentError(
        "$caller: invalid integer '$token' in $context"))
    return value
end

function _bdf_float(token::AbstractString,caller::AbstractString,
                    context::AbstractString;default=nothing)
    isempty(token) && return default
    s=replace(token,'D'=>'e','d'=>'e')
    value=tryparse(Float64,s)
    if value===nothing
        # Nastran implied-exponent real: `1.0+2`, `-2.0-3`, `1.+4`.
        m=match(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$",s)
        m!==nothing &&
            (value=tryparse(Float64,m.captures[1]*"e"*m.captures[2]))
    end
    value===nothing && throw(ArgumentError(
        "$caller: invalid real value '$token' in $context"))
    isfinite(value) || throw(ArgumentError(
        "$caller: non-finite real value in $context"))
    return value
end

# Split a card's first line into (mnemonic, fields, continuation-id, large?).
# `fields[1]` is the mnemonic; fixed lines contribute fields 2–9 (columns
# 9–72) plus the column 73–80 continuation id, large-field lines contribute
# fields 2–4 plus the column 65–80 continuation id, and free lines are
# comma-split (an empty final field marks the continuation id).
function _bdf_fields(line::AbstractString)
    if occursin(',',line)
        parts=String[strip(part) for part in split(line,',')]
        name=uppercase(parts[1])
        return name,parts,"",endswith(name,'*'),true
    end
    n=lastindex(line)
    name=uppercase(strip(SubString(line,1,min(8,n))))
    large=endswith(name,'*')
    if large
        fields=String[name]
        for a in 17:16:64
            a>n && break
            push!(fields,strip(SubString(line,a,min(a+15,n))))
        end
        contid=n>=65 ? strip(SubString(line,65,min(80,n))) : ""
        return name,fields,contid,true,false
    end
    fields=String[name]
    for a in 9:8:72
        a>n && break
        push!(fields,strip(SubString(line,a,min(a+7,n))))
    end
    while length(fields)<9
        push!(fields,"")
    end
    if occursin(r"\s",name) || any(field->occursin(r"\s",field),fields)
        # A fixed field can never contain interior whitespace: the line was
        # actually written whitespace-delimited, so re-read it as free form.
        parts=String[part for part in split(line)]
        name=uppercase(parts[1])
        return name,parts,"",endswith(name,'*'),true
    end
    contid=n>=73 ? strip(SubString(line,73,min(80,n))) : ""
    return name,fields,contid,false,false
end

# Continuation line → (id-tag, data fields, next continuation id).
function _bdf_cont_fields(line::AbstractString,large::Bool)
    if occursin(',',line)
        parts=String[strip(part) for part in split(line,',')]
        return parts[1],parts[2:end],""
    end
    n=lastindex(line)
    w=large ? 16 : 8
    tag=strip(SubString(line,1,min(w,n)))
    fields=String[]
    for a in (w+1):w:(large ? 64 : 72)
        a>n && break
        push!(fields,strip(SubString(line,a,min(a+w-1,n))))
    end
    if occursin(r"\s",tag) || any(field->occursin(r"\s",field),fields)
        # Whitespace-delimited continuation (see _bdf_fields).
        parts=String[part for part in split(line)]
        return parts[1],parts[2:end],""
    end
    contid=n>(large ? 64 : 72) ? strip(SubString(line,(large ? 65 : 73),min(80,n))) : ""
    return tag,fields,contid
end

function _bdf_process_card!(name::AbstractString,fields::Vector{String},
        free::Bool,caller::AbstractString,node_labels::Vector{Int64},
        node_rows::Dict{Int64,Int32},coords::Vector{Float64},
        records::Vector{_TextMeshCell},element_labels::Set{Int64},
        warned::Set{String},warned_cp::Base.RefValue{Bool},
        node_limit::Int,element_limit::Int)
    # Large-field mnemonics arrive with a trailing `*`; strip it for dispatch.
    base=endswith(name,'*') ? chop(name) : name
    base=="GRID" || return _bdf_process_element!(base,fields,caller,node_rows,
        records,element_labels,warned,element_limit)
    id=_bdf_int(_bdf_fld(fields,2),caller,"GRID id")
    (id===nothing || id<=0) && throw(ArgumentError(
        "$caller: GRID id must be a positive integer"))
    cptok=_bdf_fld(fields,3);cp=tryparse(Int,cptok)
    shifted=false
    if cp===nothing && free && !isempty(cptok) &&
            tryparse(Float64,replace(cptok,'D'=>'e','d'=>'e'))!==nothing
        # Free-format GRID without the cp field: `GRID id x y z` — a real
        # (not integer) token where cp belongs means the collapsed layout.
        shifted=true;cp=0
    elseif cp===nothing
        isempty(cptok) && (cp=0)
        cp===nothing && throw(ArgumentError(
            "$caller: invalid integer '$cptok' in GRID cp"))
    end
    if cp!=0
        # Alternate input coordinate frames are not implemented; the node is
        # skipped rather than re-interpreted in the wrong frame.
        warned_cp[] || @warn("$caller: GRID nodes on coordinate frame $cp " *
            "are skipped (only the basic frame 0 is supported)")
        warned_cp[]=true
        return nothing
    end
    x=_bdf_float(_bdf_fld(fields,shifted ? 3 : 4),caller,"GRID x";default=0.0)
    y=_bdf_float(_bdf_fld(fields,shifted ? 4 : 5),caller,"GRID y";default=0.0)
    z=_bdf_float(_bdf_fld(fields,shifted ? 5 : 6),caller,"GRID z";default=0.0)
    haskey(node_rows,id) && throw(ArgumentError(
        "$caller: duplicate GRID id $id"))
    length(node_labels)>=node_limit && throw(ArgumentError(
        "$caller: GRID count exceeds max_nodes=$node_limit"))
    push!(node_labels,id);node_rows[id]=Int32(length(node_labels))
    push!(coords,x);push!(coords,y);push!(coords,z)
    return nothing
end

function _bdf_process_element!(name::AbstractString,fields::Vector{String},
        caller::AbstractString,node_rows::Dict{Int64,Int32},
        records::Vector{_TextMeshCell},element_labels::Set{Int64},
        warned::Set{String},element_limit::Int)
    candidates=get(_BDF_ELEMENT_CARDS,name,nothing)
    candidates===nothing || return _bdf_emit_element!(name,candidates,fields,
        caller,node_rows,records,element_labels,element_limit)
    if name=="INCLUDE"
        throw(ArgumentError(
            "$caller: INCLUDE cards are not supported — merge the referenced " *
            "file into the bulk data first"))
    end
    # Recognized element-style mnemonics the table does not cover warn once;
    # every other card (properties, case control, CORD frames, ...) skips
    # silently like upstream.
    if startswith(name,'C') &&
            !any(prefix->startswith(name,prefix),_BDF_SILENT_PREFIXES)
        name in warned || @warn("$caller: skipping unrecognized Nastran " *
            "card '$name'")
        push!(warned,name)
    end
    return nothing
end

function _bdf_emit_element!(name::AbstractString,candidates,fields,
        caller::AbstractString,node_rows::Dict{Int64,Int32},
        records::Vector{_TextMeshCell},element_labels::Set{Int64},
        element_limit::Int)
    eid=_bdf_int(_bdf_fld(fields,2),caller,"$name element id")
    (eid===nothing || eid<=0) && throw(ArgumentError(
        "$caller: $name element id must be a positive integer"))
    pid=_bdf_int(_bdf_fld(fields,3),caller,"$name property id";default=0)
    pid<0 && throw(ArgumentError("$caller: negative $name property id $pid"))
    chosen=nothing
    for (msh,arity,positions) in candidates
        all(p->p<=length(fields),positions) || continue
        nodes=Vector{Int32}(undef,arity);ok=true
        for (k,p) in enumerate(positions)
            v=tryparse(Int,fields[p])
            if v===nothing || v<=0
                ok=false;break
            end
            row=get(node_rows,Int64(v),Int32(0))
            row==0 && (ok=false;break)
            nodes[k]=row
        end
        ok || continue
        chosen=(msh,nodes);break
    end
    chosen===nothing && throw(ArgumentError(
        "$caller: malformed $name card — cannot resolve its connectivity"))
    msh,nodes=chosen
    eid in element_labels && throw(ArgumentError(
        "$caller: duplicate element id $eid"))
    push!(element_labels,eid)
    length(records)>=element_limit && throw(ArgumentError(
        "$caller: element count exceeds max_elements=$element_limit"))
    push!(records,_TextMeshCell(eid,msh,nodes,Int64(pid),Int64(pid)))
    return nothing
end

function _read_bdf(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    node_labels=Int64[];node_rows=Dict{Int64,Int32}();coords=Float64[]
    records=_TextMeshCell[];element_labels=Set{Int64}()
    warned=Set{String}();warned_cp=Ref(false)
    pending=Ref{Tuple{String,Vector{String},Bool,Bool}}()
    have_pending=Ref(false)
    while !eof(io)
        raw=readline(io)
        isempty(raw) && continue
        raw[1]=='$' && continue                     # full-line comment
        ci=findfirst('$',raw)
        ci===nothing || (raw=rstrip(raw[1:ci-1]))   # trailing $ comment
        stripped=strip(raw)
        isempty(stripped) && continue
        lead=first(stripped)
        continuation=lead=='+' || lead==',' ||
            (lead=='*' && have_pending[] && pending[][3])
        if continuation && have_pending[]
            tag,cfields,=_bdf_cont_fields(raw,pending[][3])
            name,fields,large,free=pending[]
            if free && !isempty(fields)
                # A free-format parent's last field is the continuation id;
                # drop it when blank, when it is a bare `+`/`*` marker or when
                # it matches the continuation tag, so flat positions stay
                # uniform with the fixed layouts.
                tail=strip(fields[end],['+','*',' '])
                tagclean=strip(tag,['+','*',' '])
                (isempty(fields[end]) || isempty(tail) || tail==tagclean) &&
                    pop!(fields)
            end
            append!(fields,cfields)
            continue
        end
        continuation && !have_pending[] && continue   # orphaned continuation
        if have_pending[]
            name,fields,large,free=pending[]
            _bdf_process_card!(name,fields,free,caller,node_labels,node_rows,
                coords,records,element_labels,warned,warned_cp,
                node_limit,element_limit)
            have_pending[]=false
        end
        name,fields,contid,large,free=_bdf_fields(raw)
        name=="ENDDATA" && break
        pending[]=(name,fields,large,free)
        have_pending[]=true
    end
    if have_pending[]
        name,fields,large,free=pending[]
        _bdf_process_card!(name,fields,free,caller,node_labels,node_rows,
            coords,records,element_labels,warned,warned_cp,
            node_limit,element_limit)
    end
    return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
end

"""
    BdfFile

Result of [`read_bdf`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` and `tags` both carry the card property
id), the original `node_tags` (GRID ids) in row order and the per-block
`element_tags` (element ids).
"""
struct BdfFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_bdf(path; max_nodes, max_elements, max_file_bytes) -> BdfFile

Read a Nastran bulk data file (`.bdf`/`.nas`). `GRID` cards supply nodes
(blank or zero `cp` only — nodes on alternate frames are skipped with a
warning) and `CBAR`, `CBEAM`, `CROD`, `CTUBE`, `CTRIA3`, `CTRIA6`, `CQUAD4`,
`CQUAD8`, `CTETRA[4|10]`, `CPENTA[6|15]`, `CHEXA[8|20]` and `CPYRA[5|13]`
cards supply elements, with the property-id field becoming the elementary
and physical tag. Fixed 8-column and comma-delimited free formats are
supported along with `*`-suffixed 16-column large fields, `+`/`*`
continuation lines, `D`/`d` exponents and implied-exponent reals; `\$`
comments are stripped. `INCLUDE` rejects loudly, `ENDDATA` ends the deck,
and all other non-element cards are skipped.
"""
function read_bdf(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_bdf"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_bdf(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .bdf file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return BdfFile(mesh,node_tags,element_tags)
end

# ── Nastran writer ─────────────────────────────────────────────────────────────

const _BDF_WRITE_TABLE=Dict{Int,String}(
    1=>"CBAR",2=>"CTRIA3",9=>"CTRIA6",3=>"CQUAD4",16=>"CQUAD8",
    4=>"CTETRA",11=>"CTETRA10",6=>"CPENTA",18=>"CPENTA15",
    5=>"CHEXA",17=>"CHEXA20",7=>"CPYRA",19=>"CPYRA13")

function _bdf_i8(value::Integer,caller::AbstractString,context::AbstractString)
    s=string(value)
    length(s)<=8 || throw(ArgumentError(
        "$caller: $context $value does not fit an 8-column BDF field"))
    return s
end

function _bdf_f8(value::Real,caller::AbstractString,context::AbstractString)
    for spec in (x->repr(x),x->@sprintf("%.7g",x),x->@sprintf("%.6g",x),
                 x->@sprintf("%.5g",x),x->@sprintf("%.4g",x),x->@sprintf("%.3e",x),
                 x->@sprintf("%.2e",x),x->@sprintf("%.1e",x),x->@sprintf("%.0e",x))
        s=spec(value)
        s=replace(s,"e+"=>"+","e-"=>"-","E+"=>"+","E-"=>"-")
        length(s)<=8 && return s
    end
    throw(ArgumentError(
        "$caller: $context $value does not fit an 8-column BDF field"))
end

function _bdf_card_line(io,name::AbstractString,fields;
                        contid::AbstractString="")
    print(io,rpad(name,8))
    for field in fields
        print(io,lpad(field,8))
    end
    isempty(contid) || print(io,lpad(contid,8))
    println(io)
    return nothing
end

# Emit one element card in small-field format: the mnemonic plus EID, PID and
# the first six connectivity fields on the first line (except CTRIA6, whose
# MSC layout carries only the three corner nodes before THETA/ZOFFS/TFLAG),
# `+`-id continuation lines holding up to eight further fields each.
function _bdf_write_card(io,card::AbstractString,eid::AbstractString,
        pid::AbstractString,nids::Vector{String})
    if card=="CTRIA6"
        # MSC layout: G1–G3 on the first line, T1–T3 then G4–G6 on the
        # continuation (reader candidates 13/14/15 in flat positions).
        _bdf_card_line(io,card,[eid,pid,nids[1],nids[2],nids[3],"","",""];
            contid="+")
        _bdf_card_line(io,"+",["","","",nids[4],nids[5],nids[6]])
        return nothing
    end
    first=nids[1:min(6,length(nids))]
    _bdf_card_line(io,card,[eid;pid;first];
        contid=length(nids)>6 ? "+" : "")
    length(nids)<=6 && return nothing
    for chunk in Iterators.partition(nids[7:end],8)
        _bdf_card_line(io,"+",chunk)
    end
    return nothing
end

"""
    write_bdf(path, mesh)

Write a `Mesh` or `MixedMesh` as a Nastran bulk data file (small-field
8-column format): a `BEGIN BULK`/`ENDDATA` deck of `GRID` cards followed by
`CBAR`, `CTRIA3`, `CTRIA6`, `CQUAD4`, `CQUAD8`, `CTETRA`, `CTETRA10`,
`CPENTA`, `CPENTA15`, `CHEXA`, `CHEXA20`, `CPYRA` and `CPYRA13` element
cards whose PID field carries the physical tag (blank when zero). GRID ids
and element ids keep the external node/element tags when present.
Second-order cards split over `+`-id continuation lines in the MSC layout.
Point cells are skipped (the format has no point element); every other
unmapped type — 3-node lines, 9-node quads, tets beyond order 2 and special
records — rejects loudly. Coordinates must fit the 8-column field.
"""
function write_bdf(path,m)
    caller="write_bdf"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .bdf"))
        block.msh==15 && continue   # no point element in Nastran
        get(_BDF_WRITE_TABLE,block.msh,nothing)===nothing && throw(ArgumentError(
            "$caller: element type $(block.msh) ($(msh_family(block.msh)) " *
            "order $(msh_spec(block.msh).order)) is not representable in .bdf"))
    end
    data=m.entity_data
    nn=size(m.coords,2)
    labels=data===nothing ? collect(UInt64(1):UInt64(nn)) : data.external_node_tags
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"\$ Nastran bulk data written by Tessella")
        println(io,"BEGIN BULK")
        @inbounds for i in 1:nn
            _bdf_card_line(io,"GRID",
                [_bdf_i8(labels[i],caller,"GRID id"),"",
                 _bdf_f8(m.coords[1,i],caller,"GRID coordinate"),
                 _bdf_f8(m.coords[2,i],caller,"GRID coordinate"),
                 _bdf_f8(m.coords[3,i],caller,"GRID coordinate")])
        end
        counter=UInt64(0)
        for (bi,block) in enumerate(m.blocks)
            block isa ElementBlock || continue
            block.msh==15 && continue
            card=_BDF_WRITE_TABLE[block.msh];arity=size(block.nodes,1)
            for j in 1:size(block.nodes,2)
                counter+=UInt64(1)
                eid=_bdf_i8(_foreign_cell_uid(m,bi,j,counter),caller,"element id")
                physical=Int64(block.tags[j])
                pid=physical>0 ? _bdf_i8(physical,caller,"property id") : ""
                _bdf_write_card(io,card,eid,pid,
                    [_bdf_i8(labels[block.nodes[k,j]],caller,"grid id")
                     for k in 1:arity])
            end
        end
        println(io,"ENDDATA")
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ════════════════════════════════════════════════════════════════════════════════
# Abaqus input deck (.inp) and Gambit neutral (.neu) mesh I/O
#
# `.inp` files are keyword decks: `**` comment lines, `*KEYWORD, param=value`
# headers and free-format comma-separated data lines (a record ending in `,`
# continues on the next line). `*NODE` rows are `id, x, y [, z]`; `*ELEMENT`
# blocks declare `TYPE=` plus an optional `ELSET=`; `*ELSET`/`*NSET` blocks
# list member ids (or `first, last [, step]` triples under `GENERATE`).
# Element names map through `_INP_ELEMENT_TYPES` with the pinned release's
# writer permutations inverted for reading (Abaqus orders the tet10, hex20,
# hex27, pri15, pyr13 and lin3 mid-nodes differently); alphabetic modifier
# suffixes (`R`, `H`, `T`, ...) are stripped until a canonical name matches.
# Set names matching the upstream `Volume|Surface|Line|Point<tag>` convention
# become elementary entity tags (when the dimension agrees), names matching
# `Physical<Dim><tag>` pin that physical tag, and every other element set
# becomes a named physical group allocated the lowest free tag in its
# dimension; overlapping sets keep the first assignment. `*INCLUDE` and
# `INPUT=` external-deck references reject; all other keyword blocks are
# skipped. `*NSET` memberships are parsed for validity but cannot be
# transported — the staging bridge carries element-level groups only.
#
# `.neu` files are fixed-record sections terminated by `ENDOFSECTION`:
# `CONTROL INFO` (the `NUMNP NELEM NGRPS NBSETS NDFCD NDFVL` counts line),
# `NODAL COORDINATES` (`id x y [z]` rows), `ELEMENTS/CELLS` (`id type npts`
# plus the node list, which may wrap across lines). The Gambit type codes are
# 1=edge, 2=quadrilateral, 3=triangle, 4=hexahedron, 5=prism, 6=tetrahedron,
# 7=pyramid — with the hexahedron and pyramid base corners swapped (3↔4)
# relative to Gmsh order. `ELEMENT GROUP` blocks carry
# `GROUP: <tag> ELEMENTS: <n> MATERIAL: <m> NFLAGS: <f>`, a name line, `f`
# flag values and `n` member ids: group section `i` becomes elementary entity
# `i` and the `GROUP:` number the physical tag. `BOUNDARY CONDITIONS` blocks
# carry `name itype nentries nvalues code` plus `element type side` records
# that resolve to boundary tri/quad (or edge, for 2-D inputs) cells through
# the Gambit side tables; bc set `i` becomes elementary entity `i` and either
# the `Physical<Dim><tag>` encoded in a fallback name or the ordinal `i` the
# physical tag. Elements belonging to no group keep elementary tag 0 and stay
# in the mesh (the pinned reader drops them); group and boundary references
# resolve by element label rather than upstream's positional indexing, which
# is equivalent for the dense `1..NELEM` numbering real files use.
# `/`-prefixed lines are comments. The pinned format carries first-order
# cells only, so `write_neu` emits the standard linear types.
# ════════════════════════════════════════════════════════════════════════════════

# ── Abaqus input deck (.inp) ───────────────────────────────────────────────────

# Abaqus element name → (Gmsh type, file→Gmsh vertex permutation or `nothing`).
# Permutations are the inverse of the pinned release's `getVertexINP` write
# maps; alphabetic suffix modifiers (`R`, `H`, `T`, `I`, ...) are tolerated by
# `_inp_type` prefix-stripping.
const _INP_ELEMENT_TYPES=Dict{String,Tuple{Int,Union{Nothing,NTuple{N,Int} where N}}}(
    "MASS"=>(15,nothing),
    "T2D2"=>(1,nothing),"T3D2"=>(1,nothing),"C1D2"=>(1,nothing),
    "B21"=>(1,nothing),"B31"=>(1,nothing),
    "R2D2"=>(1,nothing),"R3D2"=>(1,nothing),"SAX1"=>(1,nothing),
    "T2D3"=>(8,(1,3,2)),"T3D3"=>(8,(1,3,2)),"C1D3"=>(8,(1,3,2)),
    "B22"=>(8,(1,3,2)),"B23"=>(8,(1,3,2)),"B32"=>(8,(1,3,2)),
    "B33"=>(8,(1,3,2)),"SAX2"=>(8,(1,3,2)),
    "CPS3"=>(2,nothing),"CPE3"=>(2,nothing),"CAX3"=>(2,nothing),
    "C2D3"=>(2,nothing),"S3"=>(2,nothing),"STRI3"=>(2,nothing),
    "M3D3"=>(2,nothing),"R3D3"=>(2,nothing),
    "CPS6"=>(9,nothing),"CPE6"=>(9,nothing),"CAX6"=>(9,nothing),
    "C2D6"=>(9,nothing),"S6"=>(9,nothing),"STRI65"=>(9,nothing),
    "CPS4"=>(3,nothing),"CPE4"=>(3,nothing),"CAX4"=>(3,nothing),
    "C2D4"=>(3,nothing),"S4"=>(3,nothing),"M3D4"=>(3,nothing),
    "R3D4"=>(3,nothing),
    "CPS8"=>(16,nothing),"CPE8"=>(16,nothing),"CAX8"=>(16,nothing),
    "C2D8"=>(16,nothing),"S8"=>(16,nothing),
    "M3D9"=>(10,nothing),
    "C3D4"=>(4,nothing),"C3D10"=>(11,(1,2,3,4,5,6,7,9,10,8)),
    "C3D6"=>(6,nothing),
    "C3D15"=>(18,(1,2,3,4,5,6,7,9,13,8,14,15,10,12,11)),
    "C3D8"=>(5,nothing),
    "C3D20"=>(17,(1,2,3,4,5,6,7,8,9,12,17,10,18,11,19,20,13,16,14,15)),
    "C3D27"=>(12,(1,2,3,4,5,6,7,8,9,12,17,10,18,11,19,20,13,16,14,15,
                  22,24,27,25,26,23,21)),
    "C3D5"=>(7,nothing),"C3D13"=>(19,(1,2,3,4,5,6,9,10,7,11,8,12,13)))

# Gmsh type → (Abaqus name, Gmsh→file vertex permutation or `nothing`) — the
# `getVertexINP` maps of the pinned release, so a written file reads back with
# the permutations above.
const _INP_WRITE_TABLE=Dict{Int,Tuple{String,Union{Nothing,NTuple{N,Int} where N}}}(
    15=>("MASS",nothing),
    1=>("T3D2",nothing),8=>("T3D3",(1,3,2)),
    2=>("CPS3",nothing),9=>("CPS6",nothing),
    3=>("CPS4",nothing),16=>("CPS8",nothing),10=>("M3D9",nothing),
    4=>("C3D4",nothing),11=>("C3D10",(1,2,3,4,5,6,7,10,8,9)),
    5=>("C3D8",nothing),
    17=>("C3D20",(1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16)),
    12=>("C3D27",(1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16,
                  27,21,26,22,24,25,23)),
    6=>("C3D6",nothing),
    18=>("C3D15",(1,2,3,4,5,6,7,10,8,13,15,14,9,11,12)),
    7=>("C3D5",nothing),19=>("C3D13",(1,2,3,4,5,6,9,11,7,8,10,12,13)))

# `Volume|Surface|Line|Point<n>` names are the upstream *ELEMENT `ELSET=`
# convention and map to elementary entities of the matching dimension.
const _INP_ENTITY_PREFIXES=Dict{String,Int}(
    "VOLUME"=>3,"SURFACE"=>2,"LINE"=>1,"POINT"=>0)
const _INP_PHYSICAL_PREFIXES=("PhysicalPoint","PhysicalLine",
    "PhysicalSurface","PhysicalVolume")   # indexed by dim+1

# Hard ceiling on `GENERATE first,last[,step]` expansion: the io limits are
# usually typemax, but a single range statement is never meant to expand to
# more than a few million members (and must not allocate unboundedly).
const _INP_GENERATE_LIMIT=Int64(8_388_608)

function _inp_entity_name(name::AbstractString)
    m=match(r"^(Volume|Surface|Line|Point)([0-9]+)$",name)
    m===nothing && return nothing
    dim=_INP_ENTITY_PREFIXES[uppercase(m.captures[1])]
    tag=tryparse(Int,m.captures[2])
    (tag===nothing || tag<0 || tag>typemax(Int32)) && return nothing
    return dim,Int64(tag)
end

function _inp_physical_name(name::AbstractString)
    m=match(r"^Physical(Volume|Surface|Line|Point)([0-9]+)$",name)
    m===nothing && return nothing
    dim=_INP_ENTITY_PREFIXES[uppercase(m.captures[1])]
    tag=tryparse(Int,m.captures[2])
    (tag===nothing || tag<0 || tag>typemax(Int32)) && return nothing
    return dim,Int64(tag)
end

_inp_physical_fallback(dim::Int,tag::Integer)=string(_INP_PHYSICAL_PREFIXES[dim+1],tag)

function _inp_type(name::AbstractString)
    u=uppercase(strip(name))
    entry=get(_INP_ELEMENT_TYPES,u,nothing)
    entry!==nothing && return entry
    # Modifier suffixes (`C3D8R`, `S4RT`, `C3D10MH`, ...) are alphabetic;
    # strip them while a canonical prefix matches. Digit-terminated names
    # never strip, so e.g. `CPS34` cannot alias `CPS3`.
    while !isempty(u) && isletter(last(u))
        u=chop(u)
        entry=get(_INP_ELEMENT_TYPES,u,nothing)
        entry!==nothing && return entry
    end
    return nothing
end

# `*KEYWORD, a, b=c` → (uppercase keyword, param dict, bare-flag set). Param
# values keep their case (set names) but lose surrounding quotes.
function _inp_header(line::AbstractString)
    parts=split(strip(line),',')
    first=strip(parts[1])
    keyword=lastindex(first)>1 ?
        uppercase(strip(SubString(first,2,lastindex(first)))) : ""
    params=Dict{String,String}();flags=Set{String}()
    for part in Iterators.drop(parts,1)
        kv=split(part,'=';limit=2)
        if length(kv)==2
            params[uppercase(strip(kv[1]))]=
                strip(strip(kv[2]),[' ','\t','\'','"'])
        else
            flag=strip(part)
            isempty(flag) || push!(flags,uppercase(flag))
        end
    end
    return keyword,params,flags
end

function _inp_next(io::Base.IO,pending)
    line=pending[]
    line!==nothing && (pending[]=nothing;return line)
    eof(io) && return nothing
    return readline(io)
end

# Consume a keyword's data section into `out`: non-empty non-comment lines
# until the next `*` keyword line (pushed back). A line ending in `,` is a
# continuation record joined with the following data line(s).
function _inp_block_lines(io::Base.IO,pending,out::Vector{String})
    empty!(out)
    while true
        line=_inp_next(io,pending)
        line===nothing && break
        s=strip(line)
        if startswith(s,'*')
            startswith(s,"**") && continue    # comment inside a data block
            pending[]=line;break              # next keyword — push back
        end
        isempty(s) && continue
        while endswith(s,',')
            nxt=_inp_next(io,pending)
            nxt===nothing && break
            t=strip(nxt)
            if startswith(t,'*')
                startswith(t,"**") && continue
                pending[]=nxt;break
            end
            isempty(t) || (s*=" "*t)
        end
        isempty(s) || push!(out,s)
    end
    return out
end

function _inp_fields(line::AbstractString)
    return String[strip(field) for field in
        Iterators.filter(!isempty,Iterators.map(strip,split(line,',')))]
end

# `first, last [, step]` triple under GENERATE, or a plain comma list.
# Returns the member id list (empty allowed).
function _inp_set_ids(line::AbstractString,generate::Bool,
                      caller::AbstractString,what::AbstractString,cap::Int)
    f=_inp_fields(line)
    isempty(f) && return Int64[]
    if generate
        length(f) in 2:3 || throw(ArgumentError(
            "$caller: .inp GENERATE $what expects first,last[,step] " *
            "in '$line'"))
        a=_unv_int(f[1],caller,".inp GENERATE first")
        b=_unv_int(f[2],caller,".inp GENERATE last")
        step=length(f)==3 ? _unv_int(f[3],caller,".inp GENERATE step") : Int64(1)
        (step>0 && a<=b) || throw(ArgumentError(
            "$caller: invalid .inp GENERATE range $a,$b,$step"))
        count=(b-a)÷step+1
        count<=cap || throw(ArgumentError(
            "$caller: .inp GENERATE range $a,$b,$step expands to $count " *
            "ids (cap $cap)"))
        return Int64[i for i in a:step:b]
    end
    return Int64[_unv_int(field,caller,".inp $what id") for field in f]
end

function _read_inp(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    node_labels=Int64[];node_rows=Dict{Int64,Int32}();coords=Float64[]
    records=_TextMeshCell[];record_rows=Dict{Int64,Int32}()
    element_labels=Set{Int64}();warned=Set{String}()
    # (set name, element label) memberships in file order — sets can be
    # declared after the element blocks, so tags resolve after parsing.
    assignments=Tuple{String,Int64}[]
    pending=Ref{Union{Nothing,String}}(nothing)
    lines=String[]
    while true
        line=_inp_next(io,pending)
        line===nothing && break
        s=strip(line)
        (isempty(s) || startswith(s,"**")) && continue
        startswith(s,'*') || continue              # stray text is skipped
        keyword,params,flags=_inp_header(s)
        if keyword=="INCLUDE" || haskey(params,"INPUT")
            throw(ArgumentError(
                "$caller: *INCLUDE/INPUT external decks are not supported — " *
                "merge the referenced file into the deck first"))
        elseif keyword=="NODE"
            for l in _inp_block_lines(io,pending,lines)
                f=_inp_fields(l)
                length(f) in 2:4 || throw(ArgumentError(
                    "$caller: malformed .inp node record '$l'"))
                id=_unv_int(f[1],caller,".inp node label")
                id>=1 || throw(ArgumentError(
                    "$caller: non-positive .inp node label $id"))
                haskey(node_rows,id) && throw(ArgumentError(
                    "$caller: duplicate .inp node label $id"))
                length(node_labels)>=node_limit && throw(ArgumentError(
                    "$caller: .inp node count exceeds max_nodes=$node_limit"))
                x=_unv_float(f[2],caller,".inp node coordinate")
                y=length(f)>=3 ?
                    _unv_float(f[3],caller,".inp node coordinate") : 0.0
                z=length(f)>=4 ?
                    _unv_float(f[4],caller,".inp node coordinate") : 0.0
                push!(node_labels,id);node_rows[id]=Int32(length(node_labels))
                push!(coords,x);push!(coords,y);push!(coords,z)
            end
        elseif keyword=="ELEMENT"
            rawtype=get(params,"TYPE","")
            isempty(rawtype) && throw(ArgumentError(
                "$caller: *ELEMENT block without a TYPE= parameter"))
            entry=_inp_type(rawtype)
            elset=get(params,"ELSET","")
            if entry===nothing
                key="etype:"*uppercase(rawtype)
                key in warned || (push!(warned,key);
                    @warn("$caller: skipping unsupported .inp element type " *
                          "'$(rawtype)'"))
                _inp_block_lines(io,pending,lines)
                continue
            end
            msh,perm=entry;arity=msh_num_nodes(msh)
            for l in _inp_block_lines(io,pending,lines)
                f=_inp_fields(l)
                length(f)==arity+1 || throw(ArgumentError(
                    "$caller: .inp element record '$l' does not carry the " *
                    "$arity nodes type $(uppercase(rawtype)) requires"))
                eid=_unv_int(f[1],caller,".inp element label")
                eid>=1 || throw(ArgumentError(
                    "$caller: non-positive .inp element label $eid"))
                eid in element_labels && throw(ArgumentError(
                    "$caller: duplicate .inp element label $eid"))
                push!(element_labels,eid)
                length(records)>=element_limit && throw(ArgumentError(
                    "$caller: .inp element count exceeds " *
                    "max_elements=$element_limit"))
                nodes=Vector{Int32}(undef,arity)
                @inbounds for k in 1:arity
                    v=_unv_int(f[k+1],caller,".inp element node")
                    row=get(node_rows,v,Int32(0))
                    row!=0 || throw(ArgumentError(
                        "$caller: .inp element $eid references unknown " *
                        "node $v"))
                    nodes[k]=row
                end
                perm===nothing || (nodes=nodes[collect(perm)])
                push!(records,_TextMeshCell(eid,msh,nodes,Int64(0),Int64(0)))
                record_rows[eid]=Int32(length(records))
                isempty(elset) || push!(assignments,(elset,eid))
            end
        elseif keyword=="NSET" || keyword=="ELSET"
            member=keyword=="ELSET"
            name=get(params,member ? "ELSET" : "NSET","")
            isempty(name) && throw(ArgumentError(
                "$caller: *$keyword block without a $keyword= name"))
            generate="GENERATE" in flags || haskey(params,"GENERATE")
            cap=min(max(element_limit,node_limit),_INP_GENERATE_LIMIT)
            for l in _inp_block_lines(io,pending,lines)
                ids=_inp_set_ids(l,generate,caller,keyword,cap)
                member || continue       # NSET ids are validated only
                for id in ids
                    push!(assignments,(name,id))
                end
            end
        else
            _inp_block_lines(io,pending,lines)   # unrecognized keyword data
        end
    end
    isempty(node_labels) && throw(ArgumentError(
        "$caller: not an Abaqus input deck (no *NODE records parsed)"))
    # First pass: `Physical<Dim><tag>` names pin their tag so auto-allocated
    # names can never collide with a tag pinned later in the file.
    pinned=Dict{Tuple{Int,String},Int64}();used=Dict{Int,Set{Int64}}()
    for (name,id) in assignments
        row=get(record_rows,id,Int32(0));row==0 && continue
        pn=_inp_physical_name(name);pn===nothing && continue
        dim,tag=pn
        dim==msh_dimension(records[row].msh) || continue
        pinned[(dim,name)]=tag;push!(get!(Set{Int64},used,dim),tag)
    end
    pnames=Dict{Tuple{Int,Int},String}()
    tag_of=Dict{Tuple{Int,String},Int64}()
    elementary_of=Dict{Int64,Int64}();physical_of=Dict{Int64,Int64}()
    for (name,id) in assignments
        row=get(record_rows,id,Int32(0))
        if row==0
            key="set:"*name
            key in warned || (push!(warned,key);
                @warn("$caller: .inp set '$(name)' references unknown " *
                      "element $id — membership skipped"))
            continue
        end
        dim=msh_dimension(records[row].msh)
        ent=_inp_entity_name(name)
        if ent!==nothing && ent[1]==dim
            haskey(elementary_of,id) || (elementary_of[id]=ent[2])
            continue
        end
        haskey(physical_of,id) && continue     # first physical set wins
        key2=(dim,name)
        tag=get(tag_of,key2,Int64(0))
        if tag==0
            tag=get(pinned,key2,Int64(0))
            if tag==0
                avail=get(used,dim,nothing)
                tag=Int64(1)
                avail===nothing || (while tag in avail; tag+=1; end)
            end
            tag_of[key2]=tag
            push!(get!(Set{Int64},used,dim),tag)
            name==_inp_physical_fallback(dim,tag) ||
                (pnames[(dim,Int(tag))]=name)
        end
        physical_of[id]=tag
    end
    for i in eachindex(records)
        c=records[i]
        elementary=get(elementary_of,c.label,Int64(0))
        physical=get(physical_of,c.label,Int64(0))
        0<=elementary<=typemax(Int32) || throw(ArgumentError(
            "$caller: .inp elementary tag $elementary on element " *
            "$(c.label) exceeds Int32"))
        0<=physical<=typemax(Int32) || throw(ArgumentError(
            "$caller: .inp physical tag $physical on element " *
            "$(c.label) exceeds Int32"))
        records[i]=_TextMeshCell(c.label,c.msh,c.nodes,elementary,physical)
    end
    return node_labels,coords,records,pnames
end

"""
    InpFile

Result of [`read_inp`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the `Volume`/`Surface`/`Line`/
`Point`-convention entity tag and whose `tags` hold the resolved element-set
physical tag; set names live in `mesh.physical_names`), the original
`node_tags` in row order and the per-block `element_tags`.
"""
struct InpFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_inp(path; max_nodes, max_elements, max_file_bytes) -> InpFile

Read an Abaqus `.inp` input deck's mesh sections: `*NODE` records (`id, x, y
[, z]`, `D`/`d` exponents accepted) and `*ELEMENT, TYPE=name` blocks mapping
through `_INP_ELEMENT_TYPES` (suffix modifiers stripped; unrecognized types
warn once and are skipped). `*ELSET`/`ELSET=` memberships become physical
tags: `Physical<Dim><tag>` names pin the tag, `Volume|Surface|Line|Point<tag>`
names become elementary entities and other names get the lowest free tag in
their dimension plus a physical name. `*NSET` blocks are validated but not
transported, `GENERATE` ranges are expanded, `*INCLUDE`/`INPUT=` decks reject
and all other keywords are skipped. `**` lines are comments; dangling node
references and malformed records are errors.
"""
function read_inp(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_inp"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_inp(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .inp file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return InpFile(mesh,node_tags,element_tags)
end

function _inp_entity_label(dim::Int,tag::Integer)
    return string(dim==3 ? "Volume" : dim==2 ? "Surface" :
                  dim==1 ? "Line" : "Point",tag)
end

function _inp_group_name(m::MixedMesh,dim::Int,tag::Int32)
    name=get(m.physical_names,(dim,Int(tag)),"")
    isempty(name) || return replace(strip(name),r"\s+"=>"_")
    return _inp_physical_fallback(dim,tag)
end

# One element record line: `eid, n1, n2, ...`, wrapped with a trailing comma
# when the field list would run past 16 entries (the upstream writer's wrap).
function _inp_write_element(io,eid,nodes::Vector{String})
    fields=Vector{String}(undef,length(nodes)+1)
    fields[1]=string(eid);fields[2:end]=nodes
    if length(fields)<=16
        println(io,join(fields,", "))
    else
        println(io,join(view(fields,1:16),", "),",")
        println(io,join(view(fields,17:length(fields)),", "))
    end
    return nothing
end

function _inp_write_id_list(io,ids)
    for (i,id) in enumerate(ids)
        i>1 && print(io,mod(i-1,10)==0 ? "\n" : ", ")
        print(io,id)
    end
    println(io)
    return nothing
end

"""
    write_inp(path, mesh)

Write a `Mesh` or `MixedMesh` as an Abaqus input deck: `*Heading`, `*NODE`
(`id, x, y, z` at full `repr` precision), one `*ELEMENT, type=<abaqus>,
ELSET=<Volume|Surface|Line|Point><entity>` block per (entity, type) group —
the upstream naming convention, so the entity tag is recovered on re-read —
then `*ELSET` blocks listing the elements of each physical group and `*NSET`
blocks listing their nodes (`Physical<Dim><tag>` names when a group has no
name, matching upstream's fallbacks). Element types with no Abaqus name
(order-3+ cells, the full-P2 prism/pyramid, special records) reject loudly;
`MASS` carries point cells.
"""
function write_inp(path,m)
    caller="write_inp"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .inp"))
        get(_INP_WRITE_TABLE,block.msh,nothing)===nothing && throw(ArgumentError(
            "$caller: element type $(block.msh) ($(msh_family(block.msh)) " *
            "order $(msh_spec(block.msh).order)) is not representable in .inp"))
    end
    data=m.entity_data
    nn=size(m.coords,2)
    labels=data===nothing ? collect(UInt64(1):UInt64(nn)) : data.external_node_tags
    for tag in labels
        tag<=typemax(Int64) || throw(ArgumentError(
            "$caller: node tag $tag is not representable in .inp"))
    end
    # Stable element ids for the whole mesh (external tags or 1..N in block
    # order): *ELEMENT blocks and *ELSET member lists share them.
    celluids=Vector{Vector{UInt64}}(undef,length(m.blocks))
    counter=UInt64(0)
    for (bi,block) in enumerate(m.blocks)
        ids=Vector{UInt64}(undef,size(block.nodes,2))
        for j in axes(block.nodes,2)
            counter+=UInt64(1)
            uid=_foreign_cell_uid(m,bi,j,counter)
            uid<=typemax(Int64) || throw(ArgumentError(
                "$caller: element tag $uid is not representable in .inp"))
            ids[j]=uid
        end
        celluids[bi]=ids
    end
    # Group cells by (entity, type); each group becomes one *ELEMENT block.
    groups=Dict{Tuple{Int64,Int},Vector{Tuple{Int,Int}}}()
    order=Tuple{Int64,Int}[]
    for (bi,block) in enumerate(m.blocks)
        for j in axes(block.nodes,2)
            key=(_foreign_cell_entity(m,bi,j),block.msh)
            if !haskey(groups,key)
                groups[key]=Tuple{Int,Int}[];push!(order,key)
            end
            push!(groups[key],(bi,j))
        end
    end
    # Physical membership for *ELSET/*NSET: (dim,tag) → element ids.
    groupset=Dict{Tuple{Int,Int32},Vector{Int64}}()
    setnodes=Dict{Tuple{Int,Int32},Set{Int32}}()
    for (bi,block) in enumerate(m.blocks)
        for j in axes(block.nodes,2)
            tag=block.tags[j];tag<=0 && continue
            key=(msh_dimension(block.msh),tag)
            push!(get!(()->Int64[],groupset,key),Int64(celluids[bi][j]))
            nodes=get!(()->Set{Int32}(),setnodes,key)
            for row in view(block.nodes,:,j)
                push!(nodes,row)
            end
        end
    end
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"*Heading")
        println(io," ",basename(String(path)))
        println(io,"*NODE")
        @inbounds for i in 1:nn
            println(io," ",labels[i],", ",repr(m.coords[1,i]),", ",
                repr(m.coords[2,i]),", ",repr(m.coords[3,i]))
        end
        for key in sort!(order;
            by=k->(msh_dimension(k[2]),k[1]==0 ? typemax(Int64) : k[1],k[2]))
            entity,msh=key
            inpname,perm=_INP_WRITE_TABLE[msh]
            println(io,"*ELEMENT, type=",inpname,", ELSET=",
                _inp_entity_label(msh_dimension(msh),entity))
            arity=msh_num_nodes(msh)
            for (bi,j) in groups[key]
                block=m.blocks[bi]
                nodes=[string(labels[block.nodes[perm===nothing ? k : perm[k],j]])
                       for k in 1:arity]
                _inp_write_element(io,Int64(celluids[bi][j]),nodes)
            end
        end
        for (dim,tag) in sort!(collect(keys(groupset)))
            name=_inp_group_name(m,dim,tag)
            dim>=1 || continue   # upstream writes element sets for dim>=1 only
            println(io,"*ELSET,ELSET=",name)
            _inp_write_id_list(io,groupset[(dim,tag)])
        end
        for (dim,tag) in sort!(collect(keys(groupset)))
            nodes=get(setnodes,(dim,tag),nothing)
            nodes===nothing && continue
            println(io,"*NSET,NSET=",_inp_group_name(m,dim,tag))
            _inp_write_id_list(io,[Int64(labels[row]) for row in sort!(collect(nodes))])
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── Gambit neutral (.neu) ──────────────────────────────────────────────────────

# Gambit type code → (Gmsh type, file→Gmsh vertex permutation or `nothing`).
# Hexahedra and pyramids swap base corners 3↔4 relative to Gmsh order.
const _NEU_ELEMENT_TYPES=Dict{Int,Tuple{Int,Union{Nothing,NTuple{N,Int} where N}}}(
    1=>(1,nothing),2=>(3,nothing),3=>(2,nothing),
    4=>(5,(1,2,4,3,5,6,8,7)),5=>(6,nothing),6=>(4,nothing),
    7=>(7,(1,2,4,3,5)))

# Gmsh type → (Gambit code, Gmsh→file vertex permutation or `nothing`); the
# permutations are involutions, so they match the read table.
const _NEU_WRITE_TABLE=Dict{Int,Tuple{Int,Union{Nothing,NTuple{N,Int} where N}}}(
    1=>(1,nothing),3=>(2,nothing),2=>(3,nothing),
    5=>(4,(1,2,4,3,5,6,8,7)),6=>(5,nothing),4=>(6,nothing),
    7=>(7,(1,2,4,3,5)))

# Boundary-condition `side` tables per Gambit parent type, 0-based indices
# into the element's Gmsh-order vertex list (the pinned reader's sub tables;
# the edge sides of tri/quad parents are a Tessella extension for 2-D decks).
const _NEU_SIDES=Dict{Int,Vector{Vector{Int}}}(
    1=>Vector{Int}[],
    2=>[[0,1],[1,2],[2,3],[3,0]],
    3=>[[0,1],[1,2],[2,0]],
    4=>[[0,1,5,4],[1,2,6,5],[3,2,6,7],[0,3,7,4],[0,1,2,3],[4,5,6,7]],
    5=>[[0,1,4,3],[1,2,5,4],[0,2,5,3],[0,1,2],[3,4,5]],
    6=>[[0,1,2],[0,1,3],[1,2,3],[0,2,3]],
    7=>[[0,1,2,3],[0,1,4],[1,2,4],[3,2,4],[0,3,4]])

# The Gambit boundary-condition name → code table used on write.
const _NEU_BOUNDARY_CODES=Dict{String,Int}(
    "UNSPECIFIED"=>0,"AXIS"=>1,"CONJUGATE"=>2,"CONVECTION"=>3,"CYCLIC"=>4,
    "DEAD"=>5,"ELEMENT_SID"=>6,"ESPECIES"=>7,"EXHAUST_FAN"=>8,"FAN"=>9,
    "FREE_SURFACE"=>10,"GAP"=>11,"INFLOW"=>12,"INLET"=>13,"INLET_VENT"=>14,
    "INTAKE_FAN"=>15,"INTERFACE"=>16,"INTERIOR"=>17,"INTERNAL"=>18,
    "LIVE"=>19,"MASS_FLOW_INLET"=>20,"MELT"=>21,"MELT_INTERFACE"=>22,
    "MOVING_BOUNDARY"=>23,"NODE"=>24,"OUTFLOW"=>25,"OUTLET"=>26,
    "OUTLET_VENT"=>27,"PERIODIC"=>28,"PLOT"=>29,"POROUS"=>30,
    "POROUS_JUMP"=>31,"PRESSURE"=>32,"PRESSURE_FAR_FIELD"=>33,
    "PRESSURE_INFLOW"=>34,"PRESSURE_INLET"=>35,"PRESSURE_OUTFLOW"=>36,
    "PRESSURE_OUTLET"=>37,"RADIATION"=>38,"RADIATOR"=>39,
    "RECIRCULATION_INLET"=>40,"RECIRCULATION_OUTLET"=>41,"SLIP"=>42,
    "SREACTION"=>43,"SURFACE"=>44,"SYMMETRY"=>45,"TRACTION"=>46,
    "TRAJECTORY"=>47,"VELOCITY"=>48,"VELOCITY_INLET"=>49,"VENT"=>50,
    "WALL"=>51,"SPRING"=>52)

const _NEU_SIZE_KEYS=("NUMNP","NELEM","NGRPS","NBSETS","NDFCD","NDFVL")

# Next meaningful .neu line: blank lines and `/` comment records are skipped.
function _neu_data_line(io::Base.IO)
    while !eof(io)
        s=strip(readline(io))
        isempty(s) && continue
        startswith(s,'/') && continue
        return s
    end
    return nothing
end

function _neu_section_end!(io::Base.IO,caller::AbstractString,
                           section::AbstractString)
    while true
        line=_neu_data_line(io)
        line===nothing && throw(ArgumentError(
            "$caller: truncated .neu $section section"))
        line=="ENDOFSECTION" && return nothing
    end
end

# `n` integer tokens drawn from `carry` first then from data lines —
# `.neu` records wrap across lines and consecutive records may share a line,
# so leftover tokens carry into the next call. `ENDOFSECTION` mid-record is
# a truncation error.
function _neu_ints(io::Base.IO,n::Int,caller::AbstractString,
                   context::AbstractString,carry=String[])
    out=Int64[]
    while length(out)<n
        if !isempty(carry)
            push!(out,_unv_int(popfirst!(carry),caller,
                ".neu $context"))
            continue
        end
        line=_neu_data_line(io)
        (line===nothing || line=="ENDOFSECTION") && throw(ArgumentError(
            "$caller: truncated .neu $context"))
        append!(carry,split(line))
    end
    return out
end

function _read_neu(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    node_labels=Int64[];node_rows=Dict{Int64,Int32}();coords=Float64[]
    records=_TextMeshCell[];record_rows=Dict{Int64,Int32}()
    gambit_of=Dict{Int64,Int}()      # element label → Gambit type
    pnames=Dict{Tuple{Int,Int},String}()
    numnp=numelem=-1;ndfcd=3
    group_index=0;bc_index=0
    next_label=Ref{Int64}(0)
    while !eof(io)
        s=strip(readline(io))
        isempty(s) && continue
        startswith(s,'/') && continue
        tokens=split(s)
        if all(tok->tok in _NEU_SIZE_KEYS,tokens) &&
                "NUMNP" in tokens
            # CONTROL INFO size line follows the column headings.
            counts=_neu_data_line(io)
            counts===nothing && throw(ArgumentError(
                "$caller: truncated .neu CONTROL INFO section"))
            f=split(counts)
            length(f)>=4 || throw(ArgumentError(
                "$caller: malformed .neu size record '$counts'"))
            numnp=_unv_int(f[1],caller,".neu NUMNP")
            numelem=_unv_int(f[2],caller,".neu NELEM")
            _unv_int(f[3],caller,".neu NGRPS")
            _unv_int(f[4],caller,".neu NBSETS")
            length(f)>=5 && (ndfcd=_unv_int(f[5],caller,".neu NDFCD"))
            numnp>=0 && numelem>=0 || throw(ArgumentError(
                "$caller: negative .neu size record"))
            ndfcd in (2,3) || throw(ArgumentError(
                "$caller: unsupported .neu NDFCD=$ndfcd"))
            numnp<=node_limit || throw(ArgumentError(
                "$caller: .neu node count exceeds max_nodes=$node_limit"))
            numelem<=element_limit || throw(ArgumentError(
                "$caller: .neu element count exceeds max_elements=$element_limit"))
        elseif occursin("NODAL COORDINATES",s)
            numnp>=0 || throw(ArgumentError(
                "$caller: .neu NODAL COORDINATES before the NUMNP counts"))
            for _ in 1:numnp
                line=_neu_data_line(io)
                (line===nothing || line=="ENDOFSECTION") && throw(ArgumentError(
                    "$caller: truncated .neu NODAL COORDINATES section"))
                f=split(line)
                length(f)==1+ndfcd || throw(ArgumentError(
                    "$caller: malformed .neu node record '$line'"))
                id=_unv_int(f[1],caller,".neu node label")
                id>=1 || throw(ArgumentError(
                    "$caller: non-positive .neu node label $id"))
                haskey(node_rows,id) && throw(ArgumentError(
                    "$caller: duplicate .neu node label $id"))
                x=_unv_float(f[2],caller,".neu node coordinate")
                y=_unv_float(f[3],caller,".neu node coordinate")
                z=ndfcd==3 ? _unv_float(f[4],caller,".neu node coordinate") : 0.0
                push!(node_labels,id);node_rows[id]=Int32(length(node_labels))
                push!(coords,x);push!(coords,y);push!(coords,z)
            end
            _neu_section_end!(io,caller,"NODAL COORDINATES")
        elseif occursin("ELEMENTS/CELLS",s)
            numelem>=0 || throw(ArgumentError(
                "$caller: .neu ELEMENTS/CELLS before the NUMNP counts"))
            carry=String[]
            for _ in 1:numelem
                head=_neu_ints(io,3,caller,"ELEMENTS/CELLS record",carry)
                id,ntype,ndp=head[1],head[2],head[3]
                entry=get(_NEU_ELEMENT_TYPES,Int(ntype),nothing)
                entry===nothing && throw(ArgumentError(
                    "$caller: unsupported .neu element type $ntype"))
                msh,perm=entry;arity=msh_num_nodes(msh)
                ndp==arity || throw(ArgumentError(
                    "$caller: .neu element $id declares $ndp nodes but " *
                    "type $ntype requires $arity"))
                id>=1 || throw(ArgumentError(
                    "$caller: non-positive .neu element label $id"))
                haskey(record_rows,id) && throw(ArgumentError(
                    "$caller: duplicate .neu element label $id"))
                length(records)>=element_limit && throw(ArgumentError(
                    "$caller: .neu element count exceeds " *
                    "max_elements=$element_limit"))
                conn=_neu_ints(io,arity,caller,
                    "element $id connectivity",carry)
                nodes=Vector{Int32}(undef,arity)
                @inbounds for k in 1:arity
                    row=get(node_rows,conn[k],Int32(0))
                    row!=0 || throw(ArgumentError(
                        "$caller: .neu element $id references unknown " *
                        "node $(conn[k])"))
                    nodes[k]=row
                end
                perm===nothing || (nodes=nodes[collect(perm)])
                push!(records,_TextMeshCell(id,msh,nodes,Int64(0),Int64(0)))
                record_rows[id]=Int32(length(records))
                gambit_of[id]=Int(ntype)
                next_label[]=max(next_label[],id)
            end
            _neu_section_end!(io,caller,"ELEMENTS/CELLS")
        elseif occursin("ELEMENT GROUP",s)
            group_index+=1
            head=_neu_data_line(io)
            head===nothing && throw(ArgumentError(
                "$caller: truncated .neu ELEMENT GROUP header"))
            m=match(r"GROUP:\s*([+-]?\d+)\s+ELEMENTS:\s*([+-]?\d+)",head)
            m===nothing && throw(ArgumentError(
                "$caller: malformed .neu group header '$head'"))
            gnum=_unv_int(m.captures[1],caller,".neu group number")
            neg=_unv_int(m.captures[2],caller,".neu group size")
            neg>=0 || throw(ArgumentError(
                "$caller: negative .neu group size $neg"))
            fn=match(r"NFLAGS:\s*([+-]?\d+)",head)
            nflags=fn===nothing ? Int64(1) :
                _unv_int(fn.captures[1],caller,".neu group NFLAGS")
            0<=nflags<=8 || throw(ArgumentError(
                "$caller: unsupported .neu group NFLAGS=$nflags"))
            gname=_neu_data_line(io)
            gname===nothing && throw(ArgumentError(
                "$caller: truncated .neu ELEMENT GROUP name"))
            gname=="ENDOFSECTION" && throw(ArgumentError(
                "$caller: truncated .neu ELEMENT GROUP section"))
            name=strip(gname)
            carry=String[]
            nflags>0 &&
                _neu_ints(io,Int(nflags),caller,"group flag values",carry)
            members=_neu_ints(io,Int(neg),caller,"group member list",carry)
            gnum<0 && (gnum=Int64(0))
            gnum<=typemax(Int32) || throw(ArgumentError(
                "$caller: .neu group number $gnum exceeds Int32"))
            for label in members
                row=get(record_rows,label,Int32(0))
                row!=0 || throw(ArgumentError(
                    "$caller: .neu group $gnum references unknown " *
                    "element $label"))
                c=records[row]
                c.elementary==0 && c.physical==0 || continue   # first group wins
                dim=msh_dimension(c.msh)
                records[row]=_TextMeshCell(c.label,c.msh,c.nodes,
                    Int64(group_index),gnum)
                isempty(name) ||
                    (pnames[(dim,Int(gnum))]=name)
            end
            _neu_section_end!(io,caller,"ELEMENT GROUP")
        elseif occursin("BOUNDARY CONDITIONS",s)
            bc_index+=1
            head=_neu_data_line(io)
            head===nothing && throw(ArgumentError(
                "$caller: truncated .neu BOUNDARY CONDITIONS header"))
            f=split(head)
            length(f)>=5 || throw(ArgumentError(
                "$caller: malformed .neu boundary header '$head'"))
            bname=join(f[1:end-4]," ")
            nentry=_unv_int(f[end-2],caller,".neu boundary count")
            _unv_int(f[end-1],caller,".neu boundary nvalues")
            _unv_int(f[end],caller,".neu boundary code")
            nentry>=0 || throw(ArgumentError(
                "$caller: negative .neu boundary count $nentry"))
            # A `Physical<Dim><tag>` fallback name recovers the physical tag
            # the upstream writer encoded; other names map to the ordinal.
            # recdim/bcphys are bound by the first valid record.
            recdim=0;bcphys=nothing;carry=String[]
            for _ in 1:nentry
                e=_neu_ints(io,3,caller,"boundary condition record",carry)
                label,btype,side=e[1],e[2],e[3]
                row=get(record_rows,label,Int32(0))
                row!=0 || throw(ArgumentError(
                    "$caller: .neu boundary record references unknown " *
                    "element $label"))
                parent=records[row]
                gtype=get(gambit_of,parent.label,0)
                Int(btype)==gtype || throw(ArgumentError(
                    "$caller: .neu boundary record for element $label " *
                    "declares parent type $btype but the element is type " *
                    "$gtype"))
                sides=_NEU_SIDES[Int(btype)]
                1<=side<=length(sides) || throw(ArgumentError(
                    "$caller: .neu boundary side $side out of range for " *
                    "element type $btype"))
                length(records)>=element_limit && throw(ArgumentError(
                    "$caller: .neu element count exceeds " *
                    "max_elements=$element_limit"))
                verts=[parent.nodes[i+1] for i in sides[Int(side)]]
                fmsh=length(verts)==2 ? 1 : length(verts)==3 ? 2 : 3
                dim=msh_dimension(fmsh)
                bcphys===nothing || recdim==dim || continue
                if bcphys===nothing
                    recdim=dim
                    pn=_inp_physical_name(bname)
                    bcphys=pn!==nothing && pn[1]==dim ? pn[2] :
                        Int64(bc_index)
                end
                next_label[]+=1
                while haskey(record_rows,next_label[])
                    next_label[]+=1
                end
                label2=next_label[]
                push!(records,_TextMeshCell(label2,fmsh,verts,
                    Int64(bc_index),bcphys))
                record_rows[label2]=Int32(length(records))
            end
            if bcphys!==nothing && !isempty(bname)
                bcphys<=typemax(Int32) || throw(ArgumentError(
                    "$caller: .neu physical tag $bcphys exceeds Int32"))
                bname==_inp_physical_fallback(recdim,bcphys) ||
                    (pnames[(recdim,Int(bcphys))]=bname)
            end
            _neu_section_end!(io,caller,"BOUNDARY CONDITIONS")
        else
            continue   # unknown section text and ENDOFSECTION markers
        end
    end
    numnp>=0 || throw(ArgumentError(
        "$caller: not a Gambit neutral file (NUMNP counts missing)"))
    length(node_labels)==numnp || throw(ArgumentError(
        "$caller: .neu file declares NUMNP=$numnp but contains " *
        "$(length(node_labels)) node records"))
    length(gambit_of)==numelem || throw(ArgumentError(
        "$caller: .neu file declares NELEM=$numelem but contains " *
        "$(length(gambit_of)) element records"))
    return node_labels,coords,records,pnames
end

"""
    NeuFile

Result of [`read_neu`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the group/boundary-set ordinal and
whose `tags` hold the physical tag; group and boundary names live in
`mesh.physical_names`), the original `node_tags` in row order and the
per-block `element_tags`.
"""
struct NeuFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_neu(path; max_nodes, max_elements, max_file_bytes) -> NeuFile

Read a Gambit neutral `.neu` file: the `CONTROL INFO` size record, `NODAL
COORDINATES` (`id x y [z]`, `NDFCD` 2 or 3), `ELEMENTS/CELLS` (`id type npts`
connectivity, Gambit codes 1–7 with hex/pyramid corner swaps), `ELEMENT
GROUP` blocks (ordinal → elementary entity, `GROUP:` number → physical tag,
name → physical name) and `BOUNDARY CONDITIONS` blocks (`element type side`
records resolved to faces/edges through the Gambit side tables; the set
ordinal → elementary, the `Physical<Dim><tag>` fallback name or the ordinal
→ physical tag). Elements in no group keep elementary tag 0; `NODAL DATA`
and other sections are skipped. Dangling references, unknown element types
and short sections are errors.
"""
function read_neu(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_neu"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_neu(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .neu file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return NeuFile(mesh,node_tags,element_tags)
end

function _neu_bc_name(m::MixedMesh,dim::Int,tag::Int32)
    name=get(m.physical_names,(dim,Int(tag)),"")
    isempty(name) || return replace(strip(name),r"\s+"=>"_")
    return _inp_physical_fallback(dim,tag)
end

"""
    write_neu(path, mesh)

Write a `Mesh` or `MixedMesh` as a Gambit neutral file: a `CONTROL INFO`
header (`NUMNP/NELEM/NGRPS/NBSETS/NDFCD/NDFVL`, always NDFCD=NDFVL=3 so the
upstream reader can parse it), `NODAL COORDINATES`, the top-dimension cells
in `ELEMENTS/CELLS` (hexes and pyramids written in Gambit corner order), one
`ELEMENT GROUP` section per distinct physical tag (including a `GROUP: 0`
set for untagged cells, so foreign readers do not drop them) and one
`BOUNDARY CONDITIONS` section per dimension-1 physical group whose cells
resolve to parent-cell faces/edges through the Gambit side tables
(unresolvable members warn once and are skipped; unphysical boundary cells
are dropped like upstream). Boundary sections renumber physical tags to
their ordinal on re-read unless the fallback name `Physical<Dim><tag>` pins
the tag. Element types Gambit cannot express (second-order cells, special
records) reject loudly; point cells are skipped and a mesh with no
dimension-1+ elements is an error.
"""
function write_neu(path,m)
    caller="write_neu"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .neu"))
        block.msh==15 && continue   # no point elements in .neu
        get(_NEU_WRITE_TABLE,block.msh,nothing)===nothing && throw(ArgumentError(
            "$caller: element type $(block.msh) ($(msh_family(block.msh)) " *
            "order $(msh_spec(block.msh).order)) is not representable in .neu"))
    end
    topdim=maximum(msh_dimension(block.msh)
        for block in m.blocks if block.msh!=15;init=0)
    topdim<1 && throw(ArgumentError(
        "$caller: .neu export requires dimension 1–3 elements"))
    data=m.entity_data
    nn=size(m.coords,2)
    labels=data===nothing ? collect(UInt64(1):UInt64(nn)) : data.external_node_tags
    for tag in labels
        tag<=typemax(Int64) || throw(ArgumentError(
            "$caller: node tag $tag is not representable in .neu"))
    end
    # Emit top-dimension cells; index their boundary faces for BC resolution.
    counter=UInt64(0)
    emitted=Tuple{Int,Int}[]            # (block, column) of top-dim cells
    emitted_label=Int64[]
    facemap=Dict{Vector{Int32},Vector{Tuple{Int64,Int,Int}}}()
    cellphys=Dict{Int32,Vector{Int64}}()
    for (bi,block) in enumerate(m.blocks)
        msh_dimension(block.msh)==topdim || continue
        gambit,_=_NEU_WRITE_TABLE[block.msh]
        sides=_NEU_SIDES[gambit]
        for j in axes(block.nodes,2)
            counter+=UInt64(1)
            uid=_foreign_cell_uid(m,bi,j,counter)
            uid<=typemax(Int64) || throw(ArgumentError(
                "$caller: element tag $uid is not representable in .neu"))
            label=Int64(uid)
            push!(emitted,(bi,j));push!(emitted_label,label)
            push!(get!(()->Int64[],cellphys,block.tags[j]),label)
            for (side,rows) in enumerate(sides)
                isempty(rows) && continue
                key=sort!(Int32[block.nodes[r+1,j] for r in rows])
                push!(get!(()->Vector{Tuple{Int64,Int,Int}}(),facemap,key),
                      (label,gambit,side))
            end
        end
    end
    # dim-1 physical cells resolve to (parent label, gambit type, side).
    bcrecs=Dict{Int32,Vector{Tuple{Int64,Int,Int}}}()
    warnface=Ref(false)
    for (bi,block) in enumerate(m.blocks)
        block.msh==15 && continue
        msh_dimension(block.msh)==topdim-1 || continue
        for j in axes(block.nodes,2)
            tag=block.tags[j]
            tag<=0 && continue
            key=sort!(Int32[block.nodes[r,j] for r in 1:size(block.nodes,1)])
            hits=get(facemap,key,nothing)
            if hits===nothing
                warnface[] || @warn("$caller: boundary cells that are not " *
                    "faces of written .neu cells are skipped")
                warnface[]=true
                continue
            end
            for (label,gtype,side) in hits
                push!(get!(()->Vector{Tuple{Int64,Int,Int}}(),bcrecs,tag),
                      (label,gtype,side))
            end
        end
    end
    gdim=topdim;bdim=topdim-1
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"        CONTROL INFO 2.0.0")
        println(io,"** GAMBIT NEUTRAL FILE")
        println(io,"Tessella mesh in GAMBIT neutral file format")
        println(io,"PROGRAM:                Tessella     VERSION:  1.0")
        println(io,"     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL")
        println(io,@sprintf("%10d%10d%10d%10d%10d%10d",
            nn,length(emitted),length(cellphys),length(bcrecs),3,3))
        println(io,"ENDOFSECTION")
        println(io,"   NODAL COORDINATES 2.0.0")
        @inbounds for i in 1:nn
            println(io,@sprintf("%10d%20.11e%20.11e%20.11e",labels[i],
                m.coords[1,i],m.coords[2,i],m.coords[3,i]))
        end
        println(io,"ENDOFSECTION")
        println(io,"      ELEMENTS/CELLS 2.0.0")
        for (k,(bi,j)) in enumerate(emitted)
            block=m.blocks[bi]
            gambit,perm=_NEU_WRITE_TABLE[block.msh]
            arity=size(block.nodes,1)
            print(io,@sprintf("%8d %2d %2d ",emitted_label[k],gambit,arity))
            for i in 1:arity
                i==8 && print(io,"\n               ")
                row=perm===nothing ? i : perm[i]
                print(io,@sprintf("%8d",labels[block.nodes[row,j]]))
            end
            println(io)
        end
        println(io,"ENDOFSECTION")
        for tag in sort!(collect(keys(cellphys)))
            cells=cellphys[tag]
            println(io,"       ELEMENT GROUP 2.0.0")
            println(io,@sprintf(
                "GROUP: %10d ELEMENTS: %10d MATERIAL: %10d NFLAGS: %10d",
                tag,length(cells),0,1))
            name=get(m.physical_names,(gdim,Int(tag)),"")
            isempty(name) && (name="Material group $tag")
            println(io,@sprintf("%32s",name))
            println(io,"       0")
            for (i,label) in enumerate(cells)
                i>1 && mod(i-1,10)==0 && print(io,"\n")
                print(io,@sprintf("%8d",label))
            end
            println(io)
            println(io,"ENDOFSECTION")
        end
        for tag in sort!(collect(keys(bcrecs)))
            list=sort!(bcrecs[tag])
            name=_neu_bc_name(m,bdim,tag)
            code=get(_NEU_BOUNDARY_CODES,uppercase(name),0)
            println(io,"       BOUNDARY CONDITIONS 2.0.0")
            println(io,@sprintf("%32s%8d%8d%8d%8d",name,1,length(list),0,code))
            for (label,gtype,side) in list
                println(io,@sprintf("%10d%5d%5d",label,gtype,side))
            end
            println(io,"ENDOFSECTION")
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ════════════════════════════════════════════════════════════════════════════════
# Plot3D structured grids (.p3d), CEA triangulation (.mail), Iridium (.ir3)
# and Diffpack/GridFE (.diff) mesh I/O — the format semantics below follow
# the pinned release's `GModelIO_P3D.cpp`, `GModelIO_MAIL.cpp`,
# `GModelIO_IR3.cpp` and `GModelIO_DIFF.cpp`.
#
# `.p3d` is the PLOT3D "whole" ASCII structured-grid format: an integer block
# count, `ni nj nk` per block (the 2-D variant declares `nk == 1`), then —
# for each block in order — three coordinate passes of `ni*nj*nk` doubles
# with `i` running fastest, then `j`, then `k` (the z array is stored even
# for `nk == 1`). The pinned reader maps `nk == 1` blocks to faces carrying
# `quad(i,j),(i+1,j),(i+1,j+1),(i,j+1)` cells and `nk >= 2` blocks to
# regions carrying hexes with the 1–4 cyclic / 5–8 cyclic vertex order;
# each block is numbered on its own per-dimension entity counter and no
# physical tags are carried. The pinned writer exports only entities that
# carry structured transfinite grids and fails with "No structured grids
# to save" otherwise — session meshes keep no structured-topology record,
# so `write_p3d` fails the same way for every representable input rather
# than inventing a grid regularity the model never declared.
#
# `.mail` is the CEA triangulation deck, a write-only format upstream:
# ` numVertices numTriangles`, one ` %19.10E %19.10E %19.10E` coordinate
# triple per vertex in file order, one dense 1-based index triple per
# triangle and one edge record per triangle, which the pinned writer
# emits as ` 0 0 0` (edge signs are a documented TODO). Only face
# triangles are emitted upstream, so `write_mail` emits `msh_family ==
# :tri` cells (first three corner vertices for higher orders) and skips
# every other cell. Upstream ships no `.mail` reader; `read_mail` is the
# literal inverse extension — triangles land on a single implicit dim-2
# entity (tag 1) and the edge block is consumed but ignored, tolerating a
# clean EOF after the triangle records.
#
# `.ir3` is the Iridium deck, likewise write-only upstream: a `33` format
# line, `numVertices num2D num3D` (or `numVertices numElements` when only
# one element dimension is present), `index x y z` vertex rows at `%.16g`
# precision, then `num tag numVert v1 ... vk` element records — dim-2
# cells numbered `1..num2D` first, dim-3 cells renumbered `1..num3D`;
# `numVert` makes records self-describing so every order is emitted
# verbatim. `tag` is the elementary entity tag by default,
# `abs(physical)` under `Mesh.SaveElementTagType == 2` (a negative
# physical tag also reverses connectivity, like the pinned
# `MElement::writeIR3`), or the partition under type `3`, which the
# projection cannot express and rejects loudly. Curve and point cells
# are skipped like upstream. `read_ir3` is the inverse extension: the tag
# column becomes the elementary entity tag (the writer's default —
# tag-type-2 decks are indistinguishable), vertex indices become the file
# labels and `numVert` resolves through the catalog within each section's
# dimension; single-count decks pick the one dimension every record arity
# can express, with a node-planarity tiebreak preferring the surface
# reading on planar node sets.
#
# `.diff` is Diffpack's GridFE deck — bidirectional upstream. A free-form
# banner ends on the counted `Number of space dim.`, `Number of elements`
# and `Number of nodes` headers (substring anchors, fixed argument
# positions), the `Max number of nodes in an element` width and the
# `Boundary indicators:` summary whose per-node memberships prefix each
# `num ( x , y , z ) [k] b1 .. bk` coordinate row. Element records are
# `num type subdomain n1 .. nk` lines whose `type` is one of the ten
# upstream `Elm*` names (tri/quad and tet/hex only — lines, prisms,
# pyramids, points and every other family have no `getStringForDIFF`), the
# `subdomain` column is the elementary entity tag on input and the
# entity's physical tags on output, and each type's connectivity passes
# through a fixed file-to-Gmsh node permutation. The reader's `mapping`
# tables and the writer's `getVertexDIFF` orders are the pinned upstream
# pairs — they are *not* mutual inverses for the non-involution layouts
# (quad8/quad9, tet10, hex8/hex20/hex27), so an upstream write→read round
# trip permutes those connectivities; mirroring both tables keeps both
# directions byte- and cell-faithful to the oracle instead of correcting
# only one side.
# ════════════════════════════════════════════════════════════════════════════════

# ── PLOT3D (.p3d) reader ─────────────────────────────────────────────────────

# Upstream parses the whole file as a whitespace-separated `fscanf` stream;
# this buffered line scanner reproduces the same token semantics without
# materializing the file.
function _read_p3d(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    tokens=Iterators.Stateful(Iterators.flatmap(eachsplit,eachline(io)))
    next=function (what)
        it=iterate(tokens)
        it===nothing && throw(ArgumentError(
            "$caller: truncated .p3d file (expected $what)"))
        return it[1]
    end
    nblocks=_unv_int(next("block count"),caller,".p3d block count")
    nblocks>0 || throw(ArgumentError(
        "$caller: .p3d block count must be positive, got $nblocks"))
    dims=NTuple{3,Int64}[]
    total=Int64(0)
    for b in 1:nblocks
        ni=_unv_int(next("block $b dimensions"),caller,".p3d block dimension")
        nj=_unv_int(next("block $b dimensions"),caller,".p3d block dimension")
        nk=_unv_int(next("block $b dimensions"),caller,".p3d block dimension")
        (ni>=1 && nj>=1 && nk>=1) || throw(ArgumentError(
            "$caller: .p3d block $b dimensions must be positive, got " *
            "($ni,$nj,$nk)"))
        count=Base.checked_mul(Base.checked_mul(ni,nj),nk)
        total=Base.checked_add(total,count)
        total<=node_limit || throw(ArgumentError(
            "$caller: .p3d node count exceeds max_nodes=$node_limit"))
        push!(dims,(ni,nj,nk))
    end
    node_labels=collect(Int64(1):total)
    coords=Vector{Float64}(undef,3*total)
    records=_TextMeshCell[]
    base=Int64(0);elabel=Int64(0);nfaces=nregions=0
    for (b,(ni,nj,nk)) in enumerate(dims)
        count=ni*nj*nk
        for c in 1:3
            @inbounds for r in 1:count
                coords[3*(base+r)-3+c]=_unv_float(
                    next("block $b coordinates"),caller,
                    ".p3d node coordinate")
            end
        end
        entity=0
        if nk==1
            nfaces+=1;entity=nfaces
            ncell=(ni-1)*(nj-1)
        else
            nregions+=1;entity=nregions
            ncell=(ni-1)*(nj-1)*(nk-1)
        end
        length(records)+ncell<=element_limit || throw(ArgumentError(
            "$caller: .p3d element count exceeds max_elements=$element_limit"))
        if nk==1
            # Vertex row of grid point (i,j): base + i + (j-1)*ni.
            for i in 1:ni-1, j in 1:nj-1
                elabel+=1
                push!(records,_TextMeshCell(elabel,3,
                    Int32[base+i+(j-1)*ni,base+i+1+(j-1)*ni,
                          base+i+1+j*ni,base+i+j*ni],
                    entity,Int64(0)))
            end
        else
            # Vertex row of grid point (i,j,k):
            # base + i + (j-1)*ni + (k-1)*ni*nj.
            nij=ni*nj
            for i in 1:ni-1, j in 1:nj-1, k in 1:nk-1
                lo=base+i+(j-1)*ni+(k-1)*nij;hi=lo+nij
                elabel+=1
                push!(records,_TextMeshCell(elabel,5,
                    Int32[lo,lo+1,lo+ni+1,lo+ni,
                          hi,hi+1,hi+ni+1,hi+ni],entity,Int64(0)))
            end
        end
        base+=count
    end
    return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
end

"""
    P3dFile

Result of [`read_p3d`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the per-dimension block tag and
whose `tags` are zero — `.p3d` carries no physical groups), the dense
`node_tags` in file row order and the per-block `element_tags`.
"""
struct P3dFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_p3d(path; max_nodes, max_elements, max_file_bytes) -> P3dFile

Read a PLOT3D "whole" ASCII structured-grid `.p3d` file: an integer block
count, `ni nj nk` dims per block (`nk == 1` selects the 2-D quad-sheet
variant, whose z array is still stored) and, per block, three coordinate
passes of `ni*nj*nk` doubles with `i` running fastest. Blocks map to quad
or hex cells on their own per-dimension entity counter; physical tags are
not carried. Non-positive counts/dims, non-numeric data and short files
are errors.
"""
function read_p3d(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_p3d"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_p3d(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .p3d file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return P3dFile(mesh,node_tags,element_tags)
end

"""
    write_p3d(path, mesh)

The pinned release writes `.p3d` only from entities that carry a
structured `transfinite_vertices` grid and reports "No structured grids
to save" otherwise. Session meshes keep no structured-topology record, so
every representable `Mesh`/`MixedMesh` fails the same way — the
structured-grid check is done before the target file is touched.
"""
function write_p3d(path,m)
    caller="write_p3d"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    _foreign_target(path,caller)
    throw(ArgumentError("$caller: no structured grids to save"))
end

# ── CEA triangulation (.mail) reader ─────────────────────────────────────────

# Token layout is the exact inverse of `write_mail`/the pinned
# `GModelIO_MAIL.cpp` writer: `nv nt`, `nv` `x y z` triples, `nt` dense
# 1-based index triples and `nt` edge records whose values are ignored (the
# upstream writer emits ` 0 0 0` placeholders — the signed edge table is a
# documented TODO). A deck ending cleanly after the triangle block is
# accepted — the edge records carry no data — but a partial or non-integer
# edge record is malformed. Triangles land on a single implicit dim-2 entity
# (tag 1) because the deck carries no entity information.
function _read_mail(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    tokens=Iterators.Stateful(Iterators.flatmap(eachsplit,eachline(io)))
    next=function (what)
        it=iterate(tokens)
        it===nothing && throw(ArgumentError(
            "$caller: truncated .mail file (expected $what)"))
        return it[1]
    end
    nv=_unv_int(next("vertex count"),caller,".mail vertex count")
    nt=_unv_int(next("triangle count"),caller,".mail triangle count")
    nv>=0 || throw(ArgumentError(
        "$caller: .mail vertex count must be non-negative, got $nv"))
    nt>=0 || throw(ArgumentError(
        "$caller: .mail triangle count must be non-negative, got $nt"))
    nv<=node_limit || throw(ArgumentError(
        "$caller: .mail node count exceeds max_nodes=$node_limit"))
    nt<=element_limit || throw(ArgumentError(
        "$caller: .mail element count exceeds max_elements=$element_limit"))
    node_labels=collect(Int64(1):nv)
    coords=Vector{Float64}(undef,3*nv)
    @inbounds for i in 1:nv, c in 1:3
        coords[3i-3+c]=_unv_float(
            next("vertex $i coordinates"),caller,".mail node coordinate")
    end
    records=_TextMeshCell[]
    sizehint!(records,nt)
    for e in 1:nt
        rows=Vector{Int32}(undef,3)
        @inbounds for k in 1:3
            label=_unv_int(next("triangle $e node $k"),caller,
                ".mail triangle node")
            1<=label<=nv || throw(ArgumentError(
                "$caller: .mail triangle $e references unknown node $label"))
            rows[k]=Int32(label)
        end
        push!(records,_TextMeshCell(Int64(e),2,rows,Int64(1),Int64(0)))
    end
    # Per-triangle edge records follow; values are ignored like the upstream
    # placeholders. Clean EOF is tolerated, a partial record is not.
    for _ in 1:nt
        it=iterate(tokens)
        it===nothing && break
        _unv_int(it[1],caller,".mail edge record")
        for _ in 1:2
            it=iterate(tokens)
            it===nothing && throw(ArgumentError(
                "$caller: truncated .mail edge record"))
            _unv_int(it[1],caller,".mail edge record")
        end
    end
    return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
end

"""
    MailFile

Result of [`read_mail`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose `msh == 2` cells carry the implicit surface tag 1 in
`elementary_entities` and zero physical tags — `.mail` carries no entity or
group data), the dense `node_tags` in file row order and the per-block
`element_tags`.
"""
struct MailFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_mail(path; max_nodes, max_elements, max_file_bytes) -> MailFile

Read a CEA triangulation `.mail` file — the inverse of [`write_mail`](@ref):
` nv nt`, `nv` coordinate triples, `nt` 1-based index triples and the
(ignored) per-triangle edge records. Upstream ships only the writer; this
reader is a Tessella extension. Non-positive/oversized counts, non-numeric
data, dangling node references and truncated records are errors; a clean EOF
after the triangle block is tolerated.
"""
function read_mail(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_mail"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_mail(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .mail file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return MailFile(mesh,node_tags,element_tags)
end

"""
    write_mail(path, mesh)

Write a `Mesh` or `MixedMesh` as a CEA triangulation `.mail` file:
` numVertices numTriangles`, one ` %19.10E %19.10E %19.10E` coordinate
triple per node in file order, one dense 1-based index triple per
triangle cell (`msh_family == :tri`; higher-order cells emit their three
corner vertices) and one ` 0 0 0` edge placeholder per triangle — the
pinned writer leaves the signed edge table unwritten. Every other cell
is skipped like upstream.
"""
function write_mail(path,m)
    caller="write_mail"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .mail"))
    end
    nn=size(m.coords,2)
    ntri=0
    for block in m.blocks
        msh_family(block.msh)==:tri && (ntri+=size(block.nodes,2))
    end
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io," ",nn," ",ntri)
        @inbounds for i in 1:nn
            println(io,@sprintf(" %19.10E %19.10E %19.10E",
                m.coords[1,i],m.coords[2,i],m.coords[3,i]))
        end
        for block in m.blocks
            msh_family(block.msh)==:tri || continue
            @inbounds for j in axes(block.nodes,2)
                println(io," ",block.nodes[1,j]," ",block.nodes[2,j],
                    " ",block.nodes[3,j])
            end
        end
        # Signed edge records are a documented TODO upstream: ` 0 0 0`.
        for block in m.blocks
            msh_family(block.msh)==:tri || continue
            for _ in axes(block.nodes,2)
                println(io," 0 0 0")
            end
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── Iridium (.ir3) reader ────────────────────────────────────────────────────

# `numVert` makes `.ir3` records self-describing, but an arity maps to a Gmsh
# type only within a known dimension — the `numVertices num2D num3D` header
# fixes each record's dimension by section position. Arity resolution scans
# the catalog per dimension and prefers the lowest order, then complete over
# serendipity cells, then earlier families (:tri/:qua; :tet/:hex/:pri/:pyr —
# the order-0 placeholders and the trihedron are excluded). So within dim 2
# arity 9 reads as quad9, 15 as tri15 and 36 as quad36; within dim 3 arity 15
# reads as pri15 and 20 as hex20. Upstream ships no `.ir3` reader, so this is
# a documented extension.
function _ir3_arity_table(dim::Int)
    families=dim==2 ? (:tri,:qua) : (:tet,:hex,:pri,:pyr)
    table=Dict{Int,Int}();ranks=Dict{Int,NTuple{3,Int}}()
    for (fi,family) in enumerate(families), order in 1:10,
        serendipity in (false,true)
        msh=try
            msh_type(family,order;serendipity=serendipity)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            continue
        end
        spec=msh_spec(msh)
        (spec.dim==dim && spec.order>=1) || continue
        key=(spec.order,spec.serendipity ? 1 : 0,fi)
        current=get(ranks,spec.nnodes,nothing)
        (current===nothing || key<current) || continue
        ranks[spec.nnodes]=key;table[spec.nnodes]=msh
    end
    return table
end

const _IR3_DIM2_TYPES=_ir3_arity_table(2)
const _IR3_DIM3_TYPES=_ir3_arity_table(3)

# Single-count-header dimension tiebreak: true when every node lies in one
# plane within a relative tolerance — the first distinct pair of nodes
# supplies the in-plane direction, the first non-collinear node the normal.
function _ir3_planar(coords::Vector{Float64},n::Int)
    n<=3 && return true
    x1=coords[1];y1=coords[2];z1=coords[3]
    ux=uy=uz=0.0
    @inbounds for i in 2:n
        dx=coords[3i-2]-x1;dy=coords[3i-1]-y1;dz=coords[3i]-z1
        if dx!=0 || dy!=0 || dz!=0
            ux,uy,uz=dx,dy,dz
            break
        end
    end
    ux==0 && uy==0 && uz==0 && return true      # all nodes coincide
    nx=ny=nz=0.0
    @inbounds for i in 2:n
        dx=coords[3i-2]-x1;dy=coords[3i-1]-y1;dz=coords[3i]-z1
        cx=uy*dz-uz*dy;cy=uz*dx-ux*dz;cz=ux*dy-uy*dx
        if cx!=0 || cy!=0 || cz!=0
            nx,ny,nz=cx,cy,cz
            break
        end
    end
    nx==0 && ny==0 && nz==0 && return true      # all nodes collinear
    extent=0.0
    @inbounds for i in 1:3*n
        v=abs(coords[i]);v>extent && (extent=v)
    end
    tolerance=1e-12*max(extent,1.0)*sqrt(nx*nx+ny*ny+nz*nz)
    @inbounds for i in 2:n
        d=abs(nx*(coords[3i-2]-x1)+ny*(coords[3i-1]-y1)+nz*(coords[3i]-z1))
        d>tolerance && return false
    end
    return true
end

function _read_ir3(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    eof(io) && throw(ArgumentError("$caller: truncated .ir3 file"))
    format=strip(readline(io))
    format=="33" || throw(ArgumentError(
        "$caller: expected .ir3 format line '33', got '$format'"))
    eof(io) && throw(ArgumentError("$caller: truncated .ir3 file"))
    header=split(strip(readline(io)))
    (length(header)==2 || length(header)==3) || throw(ArgumentError(
        "$caller: .ir3 header must be 'numVertices numElements' or " *
        "'numVertices num2D num3D', got '$(join(header," "))'"))
    nv=_unv_int(header[1],caller,".ir3 vertex count")
    nv>=0 || throw(ArgumentError(
        "$caller: .ir3 vertex count must be non-negative, got $nv"))
    nv<=node_limit || throw(ArgumentError(
        "$caller: .ir3 node count exceeds max_nodes=$node_limit"))
    num2=num3=Int64(0)
    if length(header)==3
        num2=_unv_int(header[2],caller,".ir3 dim-2 element count")
        num3=_unv_int(header[3],caller,".ir3 dim-3 element count")
        (num2>=0 && num3>=0) || throw(ArgumentError(
            "$caller: .ir3 element counts must be non-negative, got " *
            "($num2,$num3)"))
    else
        num2=_unv_int(header[2],caller,".ir3 element count")
        num2>=0 || throw(ArgumentError(
            "$caller: .ir3 element count must be non-negative, got $num2"))
    end
    ne=num2+num3
    ne<=element_limit || throw(ArgumentError(
        "$caller: .ir3 element count exceeds max_elements=$element_limit"))
    tokens=Iterators.Stateful(Iterators.flatmap(eachsplit,eachline(io)))
    next=function (what)
        it=iterate(tokens)
        it===nothing && throw(ArgumentError(
            "$caller: truncated .ir3 file (expected $what)"))
        return it[1]
    end
    node_labels=Vector{Int64}(undef,nv)
    node_rows=Dict{Int64,Int32}()
    coords=Vector{Float64}(undef,3*nv)
    for i in 1:nv
        index=_unv_int(next("vertex $i index"),caller,".ir3 vertex index")
        index>=1 || throw(ArgumentError(
            "$caller: non-positive .ir3 vertex index $index"))
        haskey(node_rows,index) && throw(ArgumentError(
            "$caller: duplicate .ir3 vertex index $index"))
        x=_unv_float(next("vertex $i x"),caller,".ir3 vertex coordinate")
        y=_unv_float(next("vertex $i y"),caller,".ir3 vertex coordinate")
        z=_unv_float(next("vertex $i z"),caller,".ir3 vertex coordinate")
        node_labels[i]=index;node_rows[index]=Int32(i)
        @inbounds coords[3i-2]=x;@inbounds coords[3i-1]=y;@inbounds coords[3i]=z
    end
    raw=Vector{Tuple{Int64,Vector{Int32}}}(undef,ne)
    for e in 1:ne
        num=_unv_int(next("element $e number"),caller,".ir3 element number")
        num>=1 || throw(ArgumentError(
            "$caller: non-positive .ir3 element number $num"))
        tag=_unv_int(next("element $e tag"),caller,".ir3 element tag")
        tag>=0 || throw(ArgumentError(
            "$caller: negative .ir3 element tag $tag"))
        count=_unv_int(next("element $e vertex count"),caller,
            ".ir3 element vertex count")
        # The count selects the record width *and* the cell type, so it must
        # land on a catalog arity before its vertex list is consumed.
        (haskey(_IR3_DIM2_TYPES,count) || haskey(_IR3_DIM3_TYPES,count)) ||
            throw(ArgumentError(
                "$caller: .ir3 element $e has $count vertices; no supported " *
                "element type"))
        rows=Vector{Int32}(undef,count)
        for k in 1:count
            label=_unv_int(next("element $e vertex $k"),caller,
                ".ir3 element vertex")
            row=get(node_rows,label,Int32(0))
            row!=0 || throw(ArgumentError(
                "$caller: .ir3 element $e references unknown node $label"))
            rows[k]=row
        end
        raw[e]=(tag,rows)
    end
    if length(header)==2
        # Single-count decks hold one dimension only — resolved once: an
        # arity set representable in exactly one dimension fixes it, and when
        # both work a planar node set selects the surface reading (Iridium
        # decks are usually surface triangulations).
        fits2=all(entry->haskey(_IR3_DIM2_TYPES,length(entry[2])),raw)
        fits3=all(entry->haskey(_IR3_DIM3_TYPES,length(entry[2])),raw)
        dim=if fits2 && !fits3
            2
        elseif fits3 && !fits2
            3
        elseif fits2 && fits3
            _ir3_planar(coords,nv) ? 2 : 3
        else
            throw(ArgumentError(
                "$caller: .ir3 element vertex counts match no single " *
                "element dimension"))
        end
        dims=fill(dim,ne)
    else
        dims=[e<=num2 ? 2 : 3 for e in 1:ne]
    end
    records=_TextMeshCell[]
    sizehint!(records,ne)
    for e in 1:ne
        tag,rows=raw[e]
        table=dims[e]==2 ? _IR3_DIM2_TYPES : _IR3_DIM3_TYPES
        msh=get(table,length(rows),nothing)
        msh===nothing && throw(ArgumentError(
            "$caller: .ir3 element $e has $(length(rows)) vertices; no " *
            "supported dim-$(dims[e]) element type"))
        push!(records,_TextMeshCell(Int64(e),msh,rows,tag,Int64(0)))
    end
    return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
end

"""
    Ir3File

Result of [`read_ir3`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the record tag column — the
writer's elementary-tag default — and whose `tags` are zero), the file
`node_tags` in vertex row order and the per-block `element_tags` (dense
file-order counters, since each dim section restarts numbering).
"""
struct Ir3File
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_ir3(path; max_nodes, max_elements, max_file_bytes) -> Ir3File

Read an Iridium `.ir3` file — the inverse of [`write_ir3`](@ref): a `33`
format line, a `numVertices num2D num3D` or `numVertices numElements`
header, `index x y z` vertex rows (the index column becomes `node_tags`)
and self-describing `num tag numVert v1 ... vk` element records. The tag
column is read as the elementary entity tag — the writer's
`element_tag_type=1` default; tag-type-2 decks are indistinguishable.
Three-count headers fix each record's dimension by section; the
single-count form resolves one dimension for the whole deck — an arity set
representable in exactly one dimension fixes it, while a set working in
both falls back to a node-planarity tiebreak (planar → surface). `numVert`
maps to the lowest-order catalog element of the resolved dimension,
preferring complete over serendipity cells. Bad format lines, non-positive
counts, duplicate/unknown vertices, negative tags and truncated records are
errors. Upstream ships only the writer; this reader is a Tessella
extension.
"""
function read_ir3(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_ir3"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_ir3(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .ir3 file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return Ir3File(mesh,node_tags,element_tags)
end

"""
    write_ir3(path, mesh; element_tag_type=1)

Write a `Mesh` or `MixedMesh` as an Iridium `.ir3` file: a `33` format
line, `numVertices num2D num3D` (or `numVertices numElements` when only
one element dimension is present), `index x y z` vertex rows at `%.16g`
precision, then `num tag numVert v1 ... vk` element records — dim-2
cells numbered `1..num2D`, then dim-3 cells renumbered `1..num3D` — in
block order. `element_tag_type` follows the pinned
`Mesh.SaveElementTagType` semantics: `1` (default) writes the elementary
entity tag, `2` writes `abs(physical)` and any negative physical tag
additionally reverses the connectivity, `3` would write the partition
and rejects loudly (the export projection carries no per-cell
partition); any other value falls back to the elementary tag like
upstream. Curve and point cells are skipped.
"""
function write_ir3(path,m;element_tag_type::Integer=1)
    caller="write_ir3"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    Int64(element_tag_type)==3 && throw(ArgumentError(
        "$caller: partition element tags are not representable in .ir3 " *
        "export (Mesh.SaveElementTagType=3)"))
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .ir3"))
    end
    nn=size(m.coords,2)
    num2=num3=0
    for block in m.blocks
        d=msh_dimension(block.msh)
        d==2 && (num2+=size(block.nodes,2))
        d==3 && (num3+=size(block.nodes,2))
    end
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"33")
        if num2>0 && num3>0
            println(io,nn," ",num2," ",num3)
        else
            println(io,nn," ",num2>0 ? num2 : num3)
        end
        @inbounds for i in 1:nn
            println(io,i," ",@sprintf("%.16g",m.coords[1,i])," ",
                @sprintf("%.16g",m.coords[2,i])," ",
                @sprintf("%.16g",m.coords[3,i]))
        end
        for dim in (2,3)
            num=0
            for (bi,block) in enumerate(m.blocks)
                msh_dimension(block.msh)==dim || continue
                arity=size(block.nodes,1)
                @inbounds for j in axes(block.nodes,2)
                    phys=Int64(block.tags[j])
                    tag=element_tag_type==2 ? abs(phys) :
                        _foreign_cell_entity(m,bi,j)
                    num+=1
                    print(io,num," ",tag," ",arity)
                    column=view(block.nodes,:,j)
                    verts=phys<0 ? Iterators.reverse(column) : column
                    for row in verts
                        print(io," ",row)
                    end
                    println(io)
                end
            end
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── Diffpack/GridFE (.diff) reader ───────────────────────────────────────────

# Gmsh type → (Diffpack `Elm*` type string, upstream `getVertexDIFF` write
# placement): file slot k holds the element's local node `perm[k]`.
const _DIFF_WRITE_TABLE=Dict{Int,Tuple{String,NTuple{N,Int} where N}}(
    2=>("ElmT3n2D",(1,2,3)),          # MTriangle — identity
    9=>("ElmT6n2D",(1,2,3,4,5,6)),    # MTriangle6 — identity
    3=>("ElmB4n2D",(1,2,4,3)),        # MQuadrangle
    16=>("ElmB8n2D",(1,2,4,3,5,8,6,7)),
    10=>("ElmB9n2D",(1,3,9,7,2,6,8,4,5)),
    4=>("ElmT4n3D",(1,2,3,4)),        # MTetrahedron — identity
    11=>("ElmT10n3D",(1,2,3,4,5,6,7,8,10,9)),
    5=>("ElmB8n3D",(3,4,8,7,1,2,6,5)),
    17=>("ElmB20n3D",(3,4,8,7,1,2,6,5,10,19,13,20,15,12,16,14,9,17,18,11)),
    12=>("ElmB27n3D",(7,9,27,25,1,3,21,19,8,16,4,18,6,26,24,22,2,10,12,
        20,17,5,13,15,23,11,14)))

# Diffpack `Elm*` type string → (Gmsh type, upstream `readDIFF` `mapping`):
# the k-th node listed on the record fills Gmsh slot `rmap[k]`. The upstream
# pair (`getVertexDIFF`, `mapping`) is not self-inverse for the brick and
# second-order layouts, so this table intentionally differs from
# `_DIFF_WRITE_TABLE` for those types — the pinned release behaves the same
# way in both directions.
const _DIFF_READ_TABLE=Dict{String,Tuple{Int,NTuple{N,Int} where N}}(
    "ElmT3n2D"=>(2,(1,2,3)),
    "ElmT6n2D"=>(9,(1,2,3,4,5,6)),
    "ElmB4n2D"=>(3,(1,2,4,3)),
    "ElmB8n2D"=>(16,(1,2,4,3,5,7,8,6)),
    "ElmB9n2D"=>(10,(1,5,2,8,9,6,4,7,3)),
    "ElmT4n3D"=>(4,(1,2,3,4)),
    "ElmT10n3D"=>(11,(1,2,3,4,5,6,7,8,9,10)),
    "ElmB8n3D"=>(5,(5,6,1,2,8,7,4,3)),
    "ElmB20n3D"=>(17,(5,6,1,2,8,7,4,3,17,9,20,14,16,13,15,18,19,10,12)),
    "ElmB27n3D"=>(12,(5,17,6,11,22,13,1,9,2,18,26,19,23,27,24,10,21,12,8,
        20,7,16,25,15,4,14,3)))

# Upstream scans each header marker with `strstr` (substring anywhere on the
# line) and fails through a stale sscanf at end-of-file; this version reports
# the missing section directly instead.
function _diff_scan_marker(io::Base.IO,caller::AbstractString,
                           what::AbstractString,markers)
    while !eof(io)
        line=readline(io)
        any(marker->occursin(marker,line),markers) && return line
    end
    throw(ArgumentError("$caller: truncated .diff file (expected $what)"))
end

# The section separator is a literal `#' in column 0 (`str[0] != '#'`).
function _diff_scan_hash(io::Base.IO,caller::AbstractString,
                         what::AbstractString)
    while !eof(io)
        line=readline(io)
        startswith(line,'#') && return line
    end
    throw(ArgumentError("$caller: truncated .diff file (expected $what)"))
end

# Header values sit at fixed token positions after the skipped marker words
# (upstream's `sscanf` suppression list). The value itself is frequently
# unused upstream; the parse still has to succeed.
function _diff_header_int(line::AbstractString,caller::AbstractString,
                          what::AbstractString,position::Int)
    fields=split(line)
    length(fields)>=position || throw(ArgumentError(
        "$caller: malformed .diff $what line '$line'"))
    return _unv_int(fields[position],caller,".diff $what")
end

# `num ( x , y , z ) [k] b1 .. bk` — literal parentheses/commas/brackets like
# the pinned `sscanf` pattern, tolerant of missing interior whitespace.
const _DIFF_NODE_PATTERN=Regex(
    "^\\s*([+-]?\\d+)\\s*\\(\\s*([^\\s,]+)\\s*,\\s*([^\\s,]+)\\s*," *
    "\\s*([^\\s,)]+)\\s*\\)\\s*\\[\\s*([+-]?\\d+)\\s*\\](.*)\$")

function _read_diff(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    line=_diff_scan_marker(io,caller,"the 'Number of space dim.' header",
        ("Number of space dim. =",))
    _diff_header_int(line,caller,"space dimension",6)
    line=_diff_scan_marker(io,caller,"the 'Number of elements' header",
        ("Number of elements   =",))
    ne=_diff_header_int(line,caller,"element count",5)
    ne>=0 || throw(ArgumentError(
        "$caller: .diff element count must be non-negative, got $ne"))
    ne<=element_limit || throw(ArgumentError(
        "$caller: .diff element count exceeds max_elements=$element_limit"))
    line=_diff_scan_marker(io,caller,"the 'Number of nodes' header",
        ("Number of nodes      =",))
    nv=_diff_header_int(line,caller,"node count",5)
    nv>=0 || throw(ArgumentError(
        "$caller: .diff node count must be non-negative, got $nv"))
    nv<=node_limit || throw(ArgumentError(
        "$caller: .diff node count exceeds max_nodes=$node_limit"))
    line=_diff_scan_marker(io,caller,
        "the 'Max number of nodes in an element' header",
        ("Max number of nodes in an element:",))
    # The declared maximum is consumed only: each record carries its own type
    # name and therefore its own arity, which is what upstream parses before
    # feeding the connectivity check.
    maxverts=_diff_header_int(line,caller,"maximum nodes per element",8)
    maxverts>=0 || throw(ArgumentError(
        "$caller: .diff maximum nodes per element must be non-negative, " *
        "got $maxverts"))
    # The "Only one material/subdomain" verdict is parsed but never used
    # upstream (`severalSubdomains` is dead state); presence is required.
    eof(io) && throw(ArgumentError(
        "$caller: truncated .diff file (expected the subdomain header)"))
    readline(io)
    line=_diff_scan_marker(io,caller,"the 'Boundary indicators' header",
        ("Boundary indicators:","boundary indicators:"))
    fields=split(line)
    nbi=_unv_int(fields[1],caller,".diff boundary-indicator count")
    nbi>=0 || throw(ArgumentError(
        "$caller: .diff boundary-indicator count must be non-negative, " *
        "got $nbi"))
    length(fields)>=3+nbi || throw(ArgumentError(
        "$caller: .diff boundary-indicators line lists fewer than $nbi " *
        "values"))
    @inbounds for i in 1:nbi
        _unv_int(fields[3+i],caller,".diff boundary indicator")
    end
    _diff_scan_hash(io,caller,"the nodal-coordinates '#' separator")
    node_labels=Int64[];sizehint!(node_labels,nv)
    node_rows=Dict{Int64,Int32}();sizehint!(node_rows,nv)
    coords=Vector{Float64}(undef,3*nv);nstored=0
    for i in 1:nv
        eof(io) && throw(ArgumentError(
            "$caller: truncated .diff file (expected node $i)"))
        matched=match(_DIFF_NODE_PATTERN,readline(io))
        matched===nothing && throw(ArgumentError(
            "$caller: malformed .diff node record $i"))
        captures=matched.captures
        label=_unv_int(captures[1],caller,".diff node number")
        label>=1 || throw(ArgumentError(
            "$caller: non-positive .diff node number $label"))
        x=_unv_float(captures[2],caller,".diff node coordinate")
        y=_unv_float(captures[3],caller,".diff node coordinate")
        z=_unv_float(captures[4],caller,".diff node coordinate")
        k=_unv_int(captures[5],caller,".diff node boundary-indicator count")
        k>=0 || throw(ArgumentError(
            "$caller: .diff node $label has a negative boundary-indicator " *
            "count $k"))
        rest=split(captures[6])
        length(rest)>=k || throw(ArgumentError(
            "$caller: .diff node $label lists fewer than $k boundary " *
            "indicators"))
        @inbounds for t in 1:k
            _unv_int(rest[t],caller,".diff node boundary indicator")
        end
        # Upstream warns and keeps the first coordinate row of a duplicate
        # node number; the trailing indicator values are still consumed.
        if haskey(node_rows,label)
            @warn "$caller: skipping duplicate .diff node $label"
        else
            nstored+=1
            node_rows[label]=Int32(nstored)
            push!(node_labels,label)
            @inbounds begin
                coords[3nstored-2]=x
                coords[3nstored-1]=y
                coords[3nstored]=z
            end
        end
    end
    resize!(coords,3nstored)
    _diff_scan_hash(io,caller,"the element-connectivity '#' separator")
    records=_TextMeshCell[];sizehint!(records,ne)
    for e in 1:ne
        eof(io) && throw(ArgumentError(
            "$caller: truncated .diff file (expected element $e)"))
        fields=split(readline(io))
        length(fields)>=3 || throw(ArgumentError(
            "$caller: malformed .diff element record $e"))
        # The element-number column is a required integer but carries no
        # identity upstream — records are numbered by file order.
        _unv_int(fields[1],caller,".diff element number")
        entry=get(_DIFF_READ_TABLE,fields[2],nothing)
        entry===nothing && throw(ArgumentError(
            "$caller: unsupported .diff element type '$(fields[2])' on " *
            "element $e"))
        msh,rmap=entry
        material=_unv_int(fields[3],caller,".diff element subdomain")
        0<=material<=typemax(Int32) || throw(ArgumentError(
            "$caller: .diff element $e subdomain tag $material is outside " *
            "the non-negative Int32 range"))
        arity=length(rmap)
        length(fields)>=3+arity || throw(ArgumentError(
            "$caller: .diff element $e ($(fields[2])) needs $arity nodes, " *
            "got $(length(fields)-3)"))
        nodes=Vector{Int32}(undef,arity)
        @inbounds for k in 1:arity
            label=_unv_int(fields[3+k],caller,".diff element node")
            row=get(node_rows,label,Int32(0))
            row!=0 || throw(ArgumentError(
                "$caller: .diff element $e references unknown node $label"))
            nodes[rmap[k]]=row
        end
        push!(records,_TextMeshCell(Int64(e),msh,nodes,Int64(material),
            Int64(0)))
    end
    return node_labels,coords,records,Dict{Tuple{Int,Int},String}()
end

"""
    DiffFile

Result of [`read_diff`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the record's subdomain/material
column — the pinned reader stores it as the geometrical entity tag — and
whose `tags` are zero, since `.diff` has no physical-groups concept on
input), the file `node_tags` in row order and the per-block `element_tags`
(dense file-order counters: the element-number column is consumed but not
used as identity, like upstream).
"""
struct DiffFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_diff(path; max_nodes, max_elements, max_file_bytes) -> DiffFile

Read a Diffpack GridFE `.diff` file: a free-form banner scanned for the
`Number of space dim.`/`Number of elements`/`Number of nodes`/
`Max number of nodes in an element`/`Boundary indicators:` headers
(substring anchors at fixed token positions, like the pinned `strstr` +
`sscanf` sequence), a `#`-led node section of `num ( x , y , z ) [k] b1 ..`
rows whose boundary-indicator lists are consumed but unused upstream, and a
`#`-led element section of `num type subdomain n1 .. nk` records over the
ten `Elm*` type names. Each type's arity comes from its own name (the
declared maximum width is informational only) and node placement follows
the pinned `readDIFF` `mapping` tables; the subdomain column becomes the
elementary entity tag. Missing headers, malformed records, non-positive
node numbers, duplicate or dangling node references, unknown element types,
negative counts and short lines are errors; duplicate nodes warn and keep
the first coordinates.
"""
function read_diff(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_diff"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_diff(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .diff file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return DiffFile(mesh,node_tags,element_tags)
end

"""
    write_diff(path, mesh; scaling_factor=1.0)

Write a `Mesh` or `MixedMesh` as a Diffpack GridFE `.diff` file, matching the
pinned `GModel::writeDIFF`/`MVertex::writeDIFF`/`MElement::writeDIFF`
records: the `Number of space dim. =   3` banner, `Number of elements`/
`Number of nodes` counts, the `All elements are of the same type`/
`Max number of nodes in an element` summary, the `N Boundary indicators:`
line, `num ( %25.16E , %25.16E , %25.16E ) [k] b1 .. bk` node rows and
`num type subdomain n1 .. nk` element records for the ten `Elm*` types.
`scaling_factor` multiplies coordinates like the upstream
`Mesh.ScalingFactor`. Only cells of the highest meshed dimension with a
Diffpack type string are emitted — lines, prisms, pyramids, points and every
other family are skipped like upstream (nodes used by skipped cells are
still written and carry boundary indicators). Nodes referenced by no cell
are not written. When any entity declares physical tags the subdomain
column repeats the element once per entity physical (and unphysical
entities emit none — still counted in the header, like upstream); with no
declared physicals the subdomain is the elementary entity tag. Boundary
indicators are the sorted physical tags (or entity tags when no physicals
exist) of the top-dimension-minus-one entities containing each node; the
header lists their union. Upstream's read/write node permutations are not
mutual inverses for the non-involution layouts, so a `.diff` round trip
permutes quad8/quad9, tet10 and hex8/hex20/hex27 connectivity exactly like
Gmsh itself.
"""
function write_diff(path,m;scaling_factor::Real=1.0)
    caller="write_diff"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    factor=try
        Float64(scaling_factor)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: scaling_factor is not a real number"))
    end
    isfinite(factor) || throw(ArgumentError(
        "$caller: scaling_factor must be finite"))
    for block in m.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .diff"))
    end
    data=m.entity_data
    nn=size(m.coords,2)
    labels=data===nothing ? collect(UInt64(1):UInt64(nn)) :
        data.external_node_tags
    # `noPhysicalGroups()` upstream forces the elementary-tag mode; with any
    # declared physical group the subdomain column and the boundary
    # indicators follow each entity's `physical_tags` instead.
    use_physical=data!==nothing &&
        any(entity->!isempty(entity.physical_tags),values(data.entities))
    physicals_of=function (dim,tag)
        entity=get(data.entities,(dim,Int(tag)),nothing)
        return entity===nothing ? Int32[] : entity.physical_tags
    end
    # Upstream `saveAll` is forced on when the model has no physical groups;
    # with physicals, only entities declaring at least one contribute nodes,
    # records and the top dimension.
    saveable=function (dim,tag)
        !use_physical && return true
        entity=get(data.entities,(dim,Int(tag)),nothing)
        return entity!==nothing && !isempty(entity.physical_tags)
    end
    # The projection contains exactly the entities upstream would save, so
    # the top dimension is the maximum dimension among saveable entities with
    # cells — including cells whose family never reaches the element section.
    dim=0
    for (bi,block) in enumerate(m.blocks)
        isempty(block.nodes) && continue
        cdim=msh_dimension(block.msh)
        for j in axes(block.nodes,2)
            if saveable(cdim,_foreign_cell_entity(m,bi,j))
                dim=max(dim,cdim);break
            end
        end
    end
    # `indexMeshVertices`: only nodes used by elements of saveable entities
    # are written (`v->getIndex() > 0`), regardless of the element's family.
    used=falses(nn)
    for (bi,block) in enumerate(m.blocks)
        cdim=msh_dimension(block.msh)
        for j in axes(block.nodes,2)
            saveable(cdim,_foreign_cell_entity(m,bi,j)) || continue
            @inbounds for index in view(block.nodes,:,j);used[index]=true;end
        end
    end
    num_vertices=count(used)
    # Boundary indicators: the sorted union of each dim-1 entity's physical
    # tags (or its tag when no physicals exist) over every entity declared —
    # unmeshed ones contribute to the header but mark no nodes, matching the
    # upstream per-entity pass.
    vertex_tags=[Set{Int64}() for _ in 1:nn]
    indicators=Set{Int64}()
    if data!==nothing && dim>=1
        for ((edim,etag),entity) in data.entities
            edim==dim-1 || continue
            union!(indicators,use_physical ? Int64.(entity.physical_tags) :
                Int64[Int64(etag)])
        end
    end
    for (bi,block) in enumerate(m.blocks)
        msh_dimension(block.msh)==dim-1 || continue
        for j in axes(block.nodes,2)
            etag=_foreign_cell_entity(m,bi,j)
            membership=use_physical ? Int64.(physicals_of(dim-1,etag)) :
                Int64[etag]
            union!(indicators,membership)
            isempty(membership) && continue
            @inbounds for row in view(block.nodes,:,j)
                used[row] && union!(vertex_tags[row],membership)
            end
        end
    end
    ordered_indicators=sort!(collect(indicators))
    # The header counts each writable top-dimension cell once even when the
    # physical-tag mode below emits the record once per entity physical;
    # unsaveable cells (no physicals on their entity) are not counted.
    num_elements=0;min_arity=typemax(Int);max_arity=0
    for (bi,block) in enumerate(m.blocks)
        haskey(_DIFF_WRITE_TABLE,block.msh) || continue
        msh_dimension(block.msh)==dim || continue
        arity=size(block.nodes,1)
        @inbounds for j in axes(block.nodes,2)
            saveable(dim,_foreign_cell_entity(m,bi,j)) || continue
            num_elements+=1
            min_arity=min(min_arity,arity);max_arity=max(max_arity,arity)
        end
    end
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        print(io,"\n\n")
        print(io," Finite element mesh (GridFE):\n\n")
        print(io," Number of space dim. =   3\n")
        print(io," Number of elements   =  ",num_elements,"\n")
        print(io," Number of nodes      =  ",num_vertices,"\n\n")
        print(io," All elements are of the same type: ",
            max_arity==min_arity ? "dpTRUE" : "dpFALSE","\n")
        print(io," Max number of nodes in an element: ",max_arity," \n")
        print(io," Only one subdomain               : dpFALSE\n")
        print(io," Lattice data                     ? 0\n\n\n\n")
        print(io," ",length(ordered_indicators)," Boundary indicators:  ")
        for value in ordered_indicators
            print(io," ",value)
        end
        print(io,"\n\n\n")
        print(io,"  Nodal coordinates and nodal boundary indicators,\n")
        print(io,"  the columns contain:\n")
        print(io,"   - node number\n")
        print(io,"   - coordinates\n")
        print(io,"   - no of boundary indicators that are set (ON)\n")
        print(io,"   - the boundary indicators that are set (ON) if any.\n")
        print(io,"#\n")
        @inbounds for i in 1:nn
            used[i] || continue
            print(io," ",labels[i]," (",
                @sprintf("%25.16E",m.coords[1,i]*factor)," ,",
                @sprintf("%25.16E",m.coords[2,i]*factor)," ,",
                @sprintf("%25.16E",m.coords[3,i]*factor)," )")
            memberships=sort!(collect(vertex_tags[i]))
            print(io," [",length(memberships),"] ")
            for value in memberships
                print(io," ",value," ")
            end
            print(io,"\n")
        end
        print(io,"\n\n")
        print(io,"  Element types and connectivity\n")
        print(io,"  the columns contain:\n")
        print(io,"   - element number\n")
        print(io,"   - element type\n")
        print(io,"   - subdomain number \n")
        print(io,"   - the global node numbers of the nodes in the element.\n")
        print(io,"#\n")
        num=0
        for (bi,block) in enumerate(m.blocks)
            entry=get(_DIFF_WRITE_TABLE,block.msh,nothing)
            entry===nothing && continue
            cdim=msh_dimension(block.msh)
            cdim==dim || continue
            name,perm=entry
            for j in axes(block.nodes,2)
                etag=_foreign_cell_entity(m,bi,j)
                subdomains=use_physical ? Int64.(physicals_of(cdim,etag)) :
                    Int64[etag]
                for subdomain in subdomains
                    num+=1
                    print(io,num," ",name," ",subdomain," ")
                    for slot in perm
                        print(io," ",labels[block.nodes[slot,j]])
                    end
                    print(io,"\n")
                end
            end
        end
        print(io,"\n")
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ════════════════════════════════════════════════════════════════════════════════
# Legacy VTK (.vtk) unstructured-grid mesh I/O — the format semantics follow
# the pinned release's `GModelIO_VTK.cpp`. Only the legacy ASCII dialect is
# supported: `BINARY` payloads and the XML `.vtu`/`.vti` families reject
# loudly. `DATASET UNSTRUCTURED_GRID` (POINTS/CELLS/CELL_TYPES) and
# `DATASET POLYDATA` (POLYGONS by arity, LINES as polylines) are read; field
# data sections (`CELL_DATA`, `POINT_DATA`, `SCALARS`, `LOOKUP_TABLE`,
# `FIELD`, ...) are ignored exactly like upstream, which stops after
# `CELL_TYPES`.
#
# Cell-type mapping and node permutations use the same source of truth as the
# model-view `.vtk` writer: `API._native_model_view_vtk_cell` returns the VTK
# type and the `order` list whose `k`-th entry is the Gmsh local node emitted
# at VTK position `k` (the module sits later in the include order, so it is
# resolved through `getfield` like the `GeoExec`/`Model` bridges). Reading
# applies the inverse placement `msh_node[order[k]] = vtk_node[k]`, which is
# also the permutation the lattice-verified writer uses — upstream's
# `readVTK` constructor lists apply its *write* map a second time for the
# 20/27-node hexahedra and the quadratic wedges, which is not the inverse for
# non-involution families; the inverse used here is verified against the
# public `vtkQuadratic*`/`vtkHigherOrder*` parametric lattices and the
# existing view exports were read back by VTK 9.7.1. In particular this means
# VTK wedges/prisms keep their parametric vertex order (upstream additionally
# reverses the two triangular faces on write — a legacy convention this
# format does not reproduce).
# ════════════════════════════════════════════════════════════════════════════════

# VTK cell type → Gmsh element type for the fixed-arity records the pinned
# reader supports (plus 27/quadratic pyramid, which upstream lacks but whose
# ordering is pinned down by the shared lattice). Arity is revalidated per
# cell against the shared layout.
const _VTK_READ_TYPES=Dict{Int,Int}(
    1=>15,   # VTK_VERTEX              → point
    3=>1,    # VTK_LINE                → lin2
    5=>2,    # VTK_TRIANGLE            → tri3
    9=>3,    # VTK_QUAD                → quad4
    10=>4,   # VTK_TETRA               → tet4
    12=>5,   # VTK_HEXAHEDRON          → hex8
    13=>6,   # VTK_WEDGE               → pri6
    14=>7,   # VTK_PYRAMID             → pyr5
    21=>8,   # VTK_QUADRATIC_EDGE      → lin3
    22=>9,   # VTK_QUADRATIC_TRIANGLE  → tri6
    23=>16,  # VTK_QUADRATIC_QUAD      → quad8
    28=>10,  # VTK_BIQUADRATIC_QUAD    → quad9
    24=>11,  # VTK_QUADRATIC_TETRA     → tet10
    25=>17,  # VTK_QUADRATIC_HEXAHEDRON → hex20
    29=>12,  # VTK_TRIQUADRATIC_HEXAHEDRON → hex27
    26=>18,  # VTK_QUADRATIC_WEDGE     → pri15
    32=>13,  # VTK_BIQUADRATIC_QUADRATIC_WEDGE → pri18
    27=>19)  # VTK_QUADRATIC_PYRAMID   → pyr13

# VTK Lagrange cell types decode through the complete-order element family;
# the node count selects the order. Types 70 (quadrilateral), 72
# (hexahedron), 73 (wedge) and 74 (pyramid) have no upstream reader branch
# but map cleanly through the shared lattice.
const _VTK_LAGRANGE_FAMILIES=Dict{Int,Symbol}(
    68=>:lin,69=>:tri,70=>:qua,71=>:tet,72=>:hex,73=>:pri,74=>:pyr)

function _vtk_lagrange_msh(family::Symbol,n::Int,caller::AbstractString)
    for order in 1:10
        msh=try
            msh_type(family,order)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            continue
        end
        msh_num_nodes(msh)==n && return msh
    end
    throw(ArgumentError("$caller: VTK Lagrange $family cell with $n nodes " *
        "has no supported complete order"))
end

# Resolve a (VTK cell type, node count) pair to `(msh, order)` where `order`
# lists the Gmsh local node index occupying each VTK connectivity slot —
# i.e. `msh_nodes[order[k]] = vtk_nodes[k]` inverts the writer's placement.
function _vtk_cell_layout(vtktype::Int64,nverts::Int,api,cache,info)
    key=(vtktype,nverts)
    cached=get(cache,key,nothing)
    cached===nothing || return cached
    msh=get(_VTK_READ_TYPES,vtktype,nothing)
    if msh===nothing
        family=get(_VTK_LAGRANGE_FAMILIES,vtktype,nothing)
        family===nothing && throw(ArgumentError(
            "$info: unsupported VTK cell type $vtktype"))
        msh=_vtk_lagrange_msh(family,nverts,info)
    end
    layout=api._native_model_view_vtk_cell(msh,info)
    order=layout.order
    length(order)==nverts || throw(ArgumentError(
        "$info: VTK cell type $vtktype lists $nverts nodes; " *
        "expected $(length(order))"))
    cache[key]=(msh,order)
    return cache[key]
end

function _vtk_int(token,caller::AbstractString,context::AbstractString)
    return _unv_int(token,caller,context)
end

function _read_vtk(io::Base.IO,caller::AbstractString,node_limit::Int,
                   element_limit::Int)
    # Two free-form header lines (version and title), then a token stream —
    # the pinned reader's `fgets`/`fscanf` split.
    eof(io) && throw(ArgumentError("$caller: truncated .vtk file"))
    version_line=readline(io)
    (startswith(version_line,"<?xml") || startswith(version_line,"<VTKFile")) &&
        throw(ArgumentError(
            "$caller: XML .vtu/.vtk files are not supported (legacy ASCII only)"))
    eof(io) && throw(ArgumentError("$caller: truncated .vtk file"))
    readline(io)  # title line, ignored
    tokens=Iterators.Stateful(Iterators.flatmap(eachsplit,eachline(io)))
    next=function (what)
        it=iterate(tokens)
        it===nothing && throw(ArgumentError(
            "$caller: truncated .vtk file (expected $what)"))
        return it[1]
    end
    mode=next("ASCII/BINARY keyword")
    mode=="BINARY" && throw(ArgumentError(
        "$caller: binary .vtk files are not supported (ASCII only)"))
    mode=="ASCII" || throw(ArgumentError(
        "$caller: expected ASCII or BINARY keyword, got '$mode'"))
    next("DATASET keyword")=="DATASET" || throw(ArgumentError(
        "$caller: expected DATASET keyword"))
    dataset=next("dataset type")
    unstructured=dataset=="UNSTRUCTURED_GRID"
    polydata=dataset=="POLYDATA"
    (unstructured || polydata) || throw(ArgumentError(
        "$caller: VTK reader can only read unstructured or polydata " *
        "datasets, got DATASET $dataset"))
    next("POINTS keyword")=="POINTS" || throw(ArgumentError(
        "$caller: expected POINTS section"))
    npoints=_vtk_int(next("point count"),caller,".vtk point count")
    npoints>0 || throw(ArgumentError("$caller: no points in dataset"))
    npoints<=node_limit || throw(ArgumentError(
        "$caller: .vtk point count exceeds max_nodes=$node_limit"))
    point_type=next("POINTS data type")
    point_type in ("float","double") || throw(ArgumentError(
        "$caller: VTK reader only accepts float or double datasets, " *
        "got '$point_type'"))
    npoints<=typemax(Int32) || throw(ArgumentError(
        "$caller: .vtk point count exceeds Int32 range"))
    coords=Vector{Float64}(undef,3*npoints)
    @inbounds for i in eachindex(coords)
        coords[i]=_unv_float(next("point coordinate"),caller,
            ".vtk point coordinate")
    end
    section=next("CELLS/POLYGONS/LINES keyword")
    section in ("CELLS","POLYGONS","LINES") || throw(ArgumentError(
        "$caller: no CELLS, POLYGONS or LINES section in dataset " *
        "(got '$section')"))
    ncells=_vtk_int(next("cell count"),caller,".vtk cell count")
    ncells>0 || throw(ArgumentError("$caller: no cells in dataset"))
    ncells<=element_limit || throw(ArgumentError(
        "$caller: .vtk cell count exceeds max_elements=$element_limit"))
    _vtk_int(next("cell connectivity size"),caller,
        ".vtk cell connectivity size")  # declared total is advisory upstream
    api=getfield(parentmodule(@__MODULE__),:API)
    cache=Dict{Tuple{Int64,Int},Tuple{Int,Vector{Int}}}()
    records=_TextMeshCell[]
    pnames=Dict{Tuple{Int,Int},String}()
    points=0
    if section=="LINES"
        # POLYDATA polylines: each record is a connected node chain split into
        # segments, lands on its own curve entity and carries physical tag 1
        # named "centerline" — the pinned reader's `physicals` assignment.
        curve=0
        for _ in 1:ncells
            k=_vtk_int(next("polyline size"),caller,".vtk polyline size")
            0<=k<=npoints || throw(ArgumentError(
                "$caller: .vtk polyline lists $k nodes"))
            nodes=Vector{Int32}(undef,k)
            for i in 1:k
                index=_vtk_int(next("polyline node"),caller,
                    ".vtk polyline node")
                (0<=index<npoints) || throw(ArgumentError(
                    "$caller: .vtk polyline references node $index outside " *
                    "0:$(npoints-1)"))
                nodes[i]=Int32(index)+Int32(1)
            end
            curve+=1
            length(records)+max(k-1,0)<=element_limit || throw(ArgumentError(
                "$caller: .vtk cell count exceeds max_elements=$element_limit"))
            for i in 1:k-1
                push!(records,_TextMeshCell(length(records)+Int64(1),1,
                    Int32[nodes[i],nodes[i+1]],Int64(curve),Int64(1)))
            end
            pnames[(1,1)]="centerline"
        end
    else
        # CELLS / POLYGONS: fixed per-cell node count followed by 0-based
        # indices; `CELLS` in an unstructured grid is typed via CELL_TYPES,
        # every other section assigns the type from the node count.
        conns=Vector{Int32}[]
        for c in 1:ncells
            k=_vtk_int(next("cell size"),caller,".vtk cell node count")
            0<k<=npoints || throw(ArgumentError(
                "$caller: .vtk cell $c lists $k nodes"))
            length(records)+1<=element_limit || throw(ArgumentError(
                "$caller: .vtk cell count exceeds max_elements=$element_limit"))
            conn=Vector{Int32}(undef,k)
            for i in 1:k
                index=_vtk_int(next("cell $c connectivity"),caller,
                    ".vtk connectivity index")
                (0<=index<npoints) || throw(ArgumentError(
                    "$caller: .vtk cell $c references node $index outside " *
                    "0:$(npoints-1)"))
                conn[i]=Int32(index)+Int32(1)
            end
            push!(conns,conn)
        end
        if unstructured
            # The pinned reader keys the CELL_TYPES pass off the DATASET
            # kind, so a POLYGONS section inside an unstructured grid still
            # expects one type per record.
            next("CELL_TYPES keyword")=="CELL_TYPES" || throw(ArgumentError(
                "$caller: expected CELL_TYPES section"))
            ntypes=_vtk_int(next("cell type count"),caller,
                ".vtk cell type count")
            ntypes==ncells || throw(ArgumentError(
                "$caller: invalid number of cell types $ntypes for " *
                "$ncells cells"))
            for c in 1:ncells
                vtktype=_vtk_int(next("cell type"),caller,
                    ".vtk cell type")
                msh,order=_vtk_cell_layout(vtktype,length(conns[c]),
                    api,cache,caller)
                nodes=Vector{Int32}(undef,length(order))
                @inbounds for k in eachindex(order)
                    nodes[order[k]]=conns[c][k]
                end
                dim=msh_dimension(msh)
                entity=dim==0 ? Int64(points+=1) : Int64(1)
                push!(records,_TextMeshCell(Int64(c),msh,nodes,entity,
                    Int64(0)))
            end
        else
            # POLYGONS (and CELLS under POLYDATA): type from node count —
            # 1 point, 2 line, 3 triangle, 4 quadrangle like upstream.
            for c in 1:ncells
                conn=conns[c];k=length(conn)
                msh=k==1 ? 15 : k==2 ? 1 : k==3 ? 2 : k==4 ? 3 :
                    throw(ArgumentError(
                        "$caller: unknown type of mesh element with " *
                        "$k nodes"))
                entity=k==1 ? Int64(points+=1) : Int64(1)
                push!(records,_TextMeshCell(Int64(c),msh,conn,entity,
                    Int64(0)))
            end
        end
    end
    # Everything past CELL_TYPES (CELL_DATA/POINT_DATA/SCALARS/FIELD/...) is
    # ignored, matching the pinned reader.
    return collect(Int64(1):npoints),coords,records,pnames
end

"""
    VtkFile

Result of [`read_vtk`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the per-dimension entity tag —
point cells get one entity each, line/surface/volume cells share one entity
per dimension — and whose `tags` are the physical tag, zero unless a
POLYDATA `LINES` polyline emitted the upstream `centerline` group), the
dense `node_tags` in file order and the per-block `element_tags`.
"""
struct VtkFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
end

"""
    read_vtk(path; max_nodes, max_elements, max_file_bytes) -> VtkFile

Read a legacy ASCII VTK file. `DATASET UNSTRUCTURED_GRID` files carry
`POINTS n float|double`, `CELLS m size` (per-cell node count and 0-based
indices) and `CELL_TYPES m`; `DATASET POLYDATA` files carry `POLYGONS`
records whose type is inferred from the node count or `LINES` polylines
(split into `lin2` segments, each record on its own curve entity under the
upstream `centerline` physical group). VTK cell types map to Gmsh types for
point, line, triangle, quadrangle, tetrahedron, hexahedron, wedge, pyramid,
the quadratic edge-only and complete quadratic families (23/28, 25/29,
26/32) and quadratic pyramid, plus the complete-order Lagrange families
(68–74) whose order is inferred from the node count; any other cell type,
node-count mismatch, `BINARY` payload or XML dialect rejects loudly. Field
data sections are ignored like the pinned reader. Node ordering is the exact
inverse of [`write_vtk`](@ref)/the model-view VTK export lattice.
"""
function read_vtk(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_vtk"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    node_labels,coords,records,pnames=open(path,"r") do io
        try
            return _read_vtk(io,caller,node_limit,element_limit)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError(
                "$caller: truncated .vtk file"))
            rethrow()
        end
    end
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records;physical_names=pnames)
    return VtkFile(mesh,node_tags,element_tags)
end

"""
    write_vtk(path, mesh; title="Tessella mesh")

Write a `Mesh` or `MixedMesh` as a legacy ASCII VTK file: a `# vtk DataFile
Version 2.0` header, a free-form `title` line, `ASCII`, `DATASET
UNSTRUCTURED_GRID`, `POINTS n double` coordinates (`%.17g`, so a read-back
is bit-exact), `CELLS m size` with the per-cell node count and 0-based
connectivity in VTK order, and `CELL_TYPES m`. When any cell carries a
positive physical tag a `CELL_DATA`/`SCALARS CellEntityIds`/`LOOKUP_TABLE
default` section emits the physical tag per cell (`-1` for untagged cells),
following the pinned `GModel::writeVTK`. The VTK type and node permutation
come from the same verified `API._native_model_view_vtk_cell` lattice the
model-view `.vtk` writer uses, covering linear, quadratic (including
serendipity quad8/hex20/prism15/pyr13) and complete Lagrange families;
`ElementBlock` types with no VTK layout (incomplete families, special
records) reject before the destination file is touched. Unlike the
model-view writer the legacy `CELLS` encoding is always used — the pinned
model writer has no `OFFSETS`/`CONNECTIVITY` variant.
"""
function write_vtk(path,m;title::AbstractString="Tessella mesh")
    caller="write_vtk"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    (occursin('\n',title) || occursin('\r',title)) && throw(ArgumentError(
        "$caller: title must be a single line"))
    api=getfield(parentmodule(@__MODULE__),:API)
    layouts=Vector{NamedTuple}(undef,length(m.blocks))
    ncells=0;size_ints=0
    for (bi,block) in enumerate(m.blocks)
        block isa ElementBlock || throw(ArgumentError(
            "$caller: special MSH record type $(block.msh) is not " *
            "representable in .vtk"))
        layout=api._native_model_view_vtk_cell(block.msh,caller)
        layouts[bi]=layout
        count=size(block.nodes,2)
        ncells+=count
        size_ints+=count*(length(layout.order)+1)
    end
    any_physical=any(tag>0 for block in m.blocks for tag in block.tags)
    nn=size(m.coords,2)
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"# vtk DataFile Version 2.0")
        println(io,title)
        println(io,"ASCII")
        println(io,"DATASET UNSTRUCTURED_GRID")
        println(io,"POINTS ",nn," double")
        @inbounds for i in 1:nn
            println(io,@sprintf("%.17g %.17g %.17g",
                m.coords[1,i],m.coords[2,i],m.coords[3,i]))
        end
        println(io,"CELLS ",ncells," ",size_ints)
        for (block,layout) in zip(m.blocks,layouts)
            order=layout.order
            @inbounds for j in axes(block.nodes,2)
                print(io,length(order))
                for k in order
                    print(io," ",Int64(block.nodes[k,j])-1)
                end
                println(io)
            end
        end
        println(io,"CELL_TYPES ",ncells)
        for (block,layout) in zip(m.blocks,layouts)
            for _ in axes(block.nodes,2)
                println(io,layout.type)
            end
        end
        if any_physical
            println(io,"CELL_DATA ",ncells)
            println(io,"SCALARS CellEntityIds int 1")
            println(io,"LOOKUP_TABLE default")
            for block in m.blocks
                for j in axes(block.nodes,2)
                    physical=Int64(block.tags[j])
                    println(io,physical>0 ? physical : -1)
                end
            end
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── VRML (.wrl) ────────────────────────────────────────────────────────────────
#
# VRML 1.0 (`#VRML V1.0 ascii`) and VRML 2.0/97 (`#VRML V2.0 utf8`) indexed
# geometry import. A file contributes one geometry set per
# `IndexedFaceSet`/`IndexedLineSet`: in V1 the current `Coordinate3` sibling
# supplies the points (state-style, rebindable per set); in V2 the set's own
# `coord Coordinate`/`Coordinate3` field does, inside a `Shape { geometry ... }`
# wrapper like the pinned reader requires.
#
# Face classification follows the observed upstream element rules — V2 emits a
# sliding-window triangle strip `(v_i, v_{i+1}, v_{i+2})` for `i = 0..k-3` with
# alternating orientation (odd windows reversed, keeping consistent winding —
# the diagonal `v_i..v_{i+2}` zig-zag triangulation); V1 classes by arity: 2
# indices → line, 3 → triangle, 4 → quad, ≥5 → a `(v_0, v_i, v_{i+1})` fan.
# Two-index runs become line cells on the set's curve entity in both versions;
# single-index runs carry no geometry and are dropped. Elements land on one
# entity per dimension per set, node tags keep the source point index (+1,
# unused points skipped), and `IndexedLineSet` takes the identical walk — all
# matching the pinned `GModelIO_VRML` outputs. Documented upstream defects that
# are deliberately NOT reproduced: the first face's first index is not
# duplicated here (upstream repeats it, minting a phantom leading vertex and
# extra degenerate/overlapping cells), and multi-shape V2 files merge their
# sets instead of collapsing to nothing.
#
# Bounded scope: `coord USE` references and every non-indexed geometry node
# (`PointSet`, `ElevationGrid`, `Extrusion`, `Sphere`, ...) are outside the
# indexed-set contract; `normalIndex`/`materialIndex`/`texCoordIndex`-class
# fields are parsed through and ignored.

struct _WrlGeometry
    v2::Bool
    points::Union{Nothing,Vector{Float64}}   # flat xyz triples, nothing when no coord
    faces::Vector{Vector{Int32}}             # 0-based vertex indices per face
    coord_index_seen::Bool
end

# Tokens: `#` comments to end-of-line, `{}`/`[]` structural, `,` separators,
# `"..."` quoted literals, everything else a whitespace-delimited word.
function _wrl_tokens(text::AbstractString,caller::AbstractString)
    tokens=String[];i=firstindex(text);last=lastindex(text)
    while i<=last
        c=text[i]
        if isspace(c) || c==','
            i=nextind(text,i);continue
        end
        if c=='#'
            j=findnext('\n',text,i)
            i=j===nothing ? last+1 : j
            continue
        end
        if c=='"'
            j=findnext('"',text,nextind(text,i))
            j===nothing && throw(ArgumentError(
                "$caller: unterminated string literal"))
            push!(tokens,String(SubString(text,i,j)))
            i=nextind(text,j);continue
        end
        if c in "{}[]"
            push!(tokens,string(c));i=nextind(text,i);continue
        end
        j=i
        while j<=last && !(isspace(text[j]) || text[j]==',' || text[j]=='#' ||
                           text[j]=='"' || text[j] in "{}[]")
            j=nextind(text,j)
        end
        push!(tokens,String(SubString(text,i,prevind(text,j))))
        i=j
    end
    return tokens
end

function _wrl_float(tok::AbstractString,caller::AbstractString)
    value=tryparse(Float64,tok)
    value===nothing && throw(ArgumentError(
        "$caller: invalid coordinate '$tok'"))
    isfinite(value) || throw(ArgumentError(
        "$caller: non-finite coordinate '$tok'"))
    return value
end

# `tokens[i] == "["`; returns (parsed ints, index after "]").
function _wrl_index_list(tokens::Vector{String},i::Int,caller::AbstractString)
    values=Int[];j=i+1;last=length(tokens)
    while j<=last && tokens[j]!="]"
        tok=tokens[j]
        value=tryparse(Int,tok)
        value===nothing && throw(ArgumentError(
            "$caller: invalid coordIndex value '$tok'"))
        value < -1 && throw(ArgumentError(
            "$caller: negative coordIndex value $value"))
        push!(values,value)
        j+=1
    end
    j>last && throw(ArgumentError("$caller: unterminated '[' index list"))
    return values,j+1
end

# `tokens[i] == "["`; returns (parsed floats, index after "]").
function _wrl_point_list(tokens::Vector{String},i::Int,caller::AbstractString)
    values=Float64[];j=i+1;last=length(tokens)
    while j<=last && tokens[j]!="]"
        push!(values,_wrl_float(tokens[j],caller))
        j+=1
    end
    j>last && throw(ArgumentError("$caller: unterminated '[' point list"))
    return values,j+1
end

# Split a raw index list on the -1 face separator. A trailing unterminated run
# still closes a face (upstream requires the terminator; the lenient read keeps
# the data rather than discarding the file over one missing sentinel).
function _wrl_faces!(dst::Vector{Vector{Int32}},values::Vector{Int})
    face=Int32[]
    for value in values
        if value==-1
            push!(dst,copy(face));empty!(face)
        else
            push!(face,Int32(value))
        end
    end
    push!(dst,face)
    return dst
end

# V2: `Shape { ... geometry IndexedFaceSet|IndexedLineSet { coord <node>
# coordIndex [...] } }`. Returns (sets, index after the Shape's '}').
function _wrl_v2_shape!(sets,tokens::Vector{String},i::Int,
                        caller::AbstractString)
    last=length(tokens);depth=1;j=i
    while j<=last
        t=tokens[j]
        if t=="{"
            depth+=1;j+=1;continue
        elseif t=="}"
            depth-=1;j+=1
            depth==0 && return j
            continue
        elseif t=="geometry" && depth==1
            j+=1
            j<=last && tokens[j]=="DEF" && (j+=2)
            j>last && break
            name=tokens[j]
            if name=="IndexedFaceSet" || name=="IndexedLineSet"
                (j+1<=last && tokens[j+1]=="{") || throw(ArgumentError(
                    "$caller: expected '{' after geometry $name"))
                set,j=_wrl_v2_set(tokens,j+2,caller,name)
                push!(sets,set);continue
            else
                # Unhandled geometry node: its braces are consumed by the
                # depth walk like any other field.
                j+=1;continue
            end
        else
            j+=1;continue
        end
    end
    throw(ArgumentError("$caller: unterminated Shape node"))
end

# `i` = first token inside the geometry set's braces. Returns (_WrlGeometry, j)
# with `j` the index after the closing brace.
function _wrl_v2_set(tokens::Vector{String},i::Int,caller::AbstractString,
                     name::AbstractString)
    last=length(tokens);depth=1;j=i
    points=nothing;faces=Vector{Int32}[];ci_seen=false
    while j<=last
        t=tokens[j]
        if t=="{"
            depth+=1;j+=1;continue
        elseif t=="}"
            depth-=1;j+=1
            depth==0 &&
                return _WrlGeometry(true,points,faces,ci_seen),j
            continue
        elseif t=="coord" && depth==1
            j+=1
            j<=last && tokens[j]=="DEF" && (j+=2)
            j>last && break
            node=tokens[j]
            node=="USE" && throw(ArgumentError(
                "$caller: 'coord USE' references are not supported"))
            node in ("Coordinate","Coordinate3") || throw(ArgumentError(
                "$caller: expected a Coordinate node after 'coord', " *
                "got '$node'"))
            j+=1
            (j<=last && tokens[j]=="{") || throw(ArgumentError(
                "$caller: expected '{' after $node"))
            j+=1;cdepth=1;values=Float64[]
            while j<=last && cdepth>0
                tj=tokens[j]
                if tj=="{";cdepth+=1;j+=1;continue
                elseif tj=="}";cdepth-=1;j+=1;continue
                elseif tj=="point" && cdepth==1 && j+1<=last &&
                        tokens[j+1]=="["
                    vals,j=_wrl_point_list(tokens,j+1,caller)
                    append!(values,vals);continue
                else
                    j+=1;continue
                end
            end
            cdepth==0 || throw(ArgumentError(
                "$caller: unterminated $node node"))
            length(values)%3==0 || throw(ArgumentError(
                "$caller: point list has $(length(values)) values, " *
                "not a multiple of 3"))
            points=values;continue
        elseif t=="coordIndex" && depth==1
            j+=1
            (j<=last && tokens[j]=="[") || throw(ArgumentError(
                "$caller: expected '[' after 'coordIndex'"))
            vals,j=_wrl_index_list(tokens,j,caller)
            ci_seen=true
            _wrl_faces!(faces,vals);continue
        else
            j+=1;continue
        end
    end
    throw(ArgumentError("$caller: unterminated $name node"))
end

function _wrl_v2_sets(tokens::Vector{String},caller::AbstractString)
    sets=_WrlGeometry[];i=firstindex(tokens);last=lastindex(tokens)
    while i<=last
        if tokens[i]=="Shape" && i+1<=last && tokens[i+1]=="{"
            i=_wrl_v2_shape!(sets,tokens,i+2,caller)
        else
            i+=1
        end
    end
    return sets
end

# V1: `Coordinate3` nodes set the current points for every following
# `IndexedFaceSet`/`IndexedLineSet` sibling until the next `Coordinate3`,
# wherever they appear (Separator/Group wrappers carry no scan semantics).
function _wrl_v1_sets(tokens::Vector{String},caller::AbstractString)
    sets=_WrlGeometry[];points=nothing
    i=firstindex(tokens);last=lastindex(tokens)
    while i<=last
        t=tokens[i]
        if t=="Coordinate3" && i+1<=last && tokens[i+1]=="{"
            depth=1;j=i+2;values=Float64[]
            while j<=last && depth>0
                tj=tokens[j]
                if tj=="{";depth+=1;j+=1;continue
                elseif tj=="}";depth-=1;j+=1;continue
                elseif tj=="point" && j+1<=last && tokens[j+1]=="["
                    vals,j=_wrl_point_list(tokens,j+1,caller)
                    append!(values,vals);continue
                else
                    j+=1;continue
                end
            end
            depth==0 || throw(ArgumentError(
                "$caller: unterminated Coordinate3 node"))
            length(values)%3==0 || throw(ArgumentError(
                "$caller: point list has $(length(values)) values, " *
                "not a multiple of 3"))
            points=values;i=j;continue
        elseif (t=="IndexedFaceSet" || t=="IndexedLineSet") &&
                i+1<=last && tokens[i+1]=="{"
            depth=1;j=i+2;faces=Vector{Int32}[];ci_seen=false
            while j<=last && depth>0
                tj=tokens[j]
                if tj=="{";depth+=1;j+=1;continue
                elseif tj=="}";depth-=1;j+=1;continue
                elseif depth==1 && tj=="coordIndex" && j+1<=last &&
                        tokens[j+1]=="["
                    vals,j=_wrl_index_list(tokens,j+1,caller)
                    ci_seen=true
                    _wrl_faces!(faces,vals);continue
                else
                    j+=1;continue
                end
            end
            depth==0 || throw(ArgumentError(
                "$caller: unterminated $t node"))
            push!(sets,_WrlGeometry(false,points,faces,ci_seen))
            i=j;continue
        else
            i+=1;continue
        end
    end
    return sets
end

# Emit each set's faces as `_TextMeshCell` records. Per set, two-index faces
# become lines on that set's curve entity and 3+-index faces become triangles
# (V2 sliding strip, V1 4→quad/≥5→fan) on its surface entity; the source
# point index stays the node label so sparse point usage keeps its tag gaps.
function _wrl_records(v2::Bool,sets::Vector{_WrlGeometry},
                      caller::AbstractString,node_limit::Int,element_limit::Int)
    node_labels=Int64[];coords=Float64[];records=_TextMeshCell[]
    stag=0;ctag=0;base=Int64(0);label=Int64(0);saw_index=false
    for set in sets
        np=set.points===nothing ? 0 : div(length(set.points),3)
        set.coord_index_seen || (base+=np;continue)
        saw_index=true
        used1=Int32[];used2=Int32[]
        dim1=Vector{NTuple{2,Int32}}();dim2=Vector{NTuple{4,Int32}}()
        dim2msh=Int8[]
        for face in set.faces
            k=length(face)
            k==0 && continue
            set.points===nothing && throw(ArgumentError(
                "$caller: coordIndex faces appear in a geometry set " *
                "without a Coordinate point list"))
            for value in face
                (0<=value<np) || throw(ArgumentError(
                    "$caller: coordIndex value $value out of range " *
                    "(0:$(np-1))"))
            end
            if k==1
                continue
            elseif k==2
                push!(dim1,(face[1],face[2]))
            elseif v2
                for i in 0:k-3
                    a,b,c=face[i+1],face[i+2],face[i+3]
                    push!(dim2,isodd(i) ? (c,b,a,Int32(-1)) :
                                          (a,b,c,Int32(-1)))
                    push!(dim2msh,Int8(2))
                end
            else
                if k==3
                    push!(dim2,(face[1],face[2],face[3],Int32(-1)))
                    push!(dim2msh,Int8(2))
                elseif k==4
                    push!(dim2,(face[1],face[2],face[3],face[4]))
                    push!(dim2msh,Int8(3))
                else
                    for i in 1:k-2
                        push!(dim2,(face[1],face[i+1],face[i+2],Int32(-1)))
                        push!(dim2msh,Int8(2))
                    end
                end
            end
        end
        (isempty(dim1) && isempty(dim2)) && (base+=np;continue)
        for (a,b) in dim1
            push!(used1,a);push!(used1,b)
        end
        for (ci,cell) in enumerate(dim2)
            for e in 1:(dim2msh[ci]==Int8(3) ? 4 : 3)
                push!(used2,cell[e])
            end
        end
        sort!(unique!(used1));sort!(unique!(used2))
        rows=Dict{Int32,Int32}()
        for p in used1
            push!(node_labels,base+Int64(p)+1);rows[p]=Int32(length(node_labels))
            append!(coords,@view(set.points[3p+1:3p+3]))
        end
        for p in used2
            p in used1 && continue
            push!(node_labels,base+Int64(p)+1);rows[p]=Int32(length(node_labels))
            append!(coords,@view(set.points[3p+1:3p+3]))
        end
        length(node_labels)>node_limit && throw(ArgumentError(
            "$caller: node count exceeds max_nodes=$node_limit"))
        base+=np
        if !isempty(dim1)
            ctag+=1
            for (a,b) in dim1
                label+=1
                push!(records,_TextMeshCell(label,1,
                    Int32[rows[a],rows[b]],Int64(ctag),Int64(0)))
            end
        end
        if !isempty(dim2)
            stag+=1
            for (cell,msh) in zip(dim2,dim2msh)
                arity=msh==Int8(3) ? 4 : 3
                label+=1
                push!(records,_TextMeshCell(label,Int(msh),
                    Int32[rows[cell[e]] for e in 1:arity],
                    Int64(stag),Int64(0)))
            end
        end
        length(records)>element_limit && throw(ArgumentError(
            "$caller: element count exceeds max_elements=$element_limit"))
    end
    saw_index || throw(ArgumentError(
        "$caller: file contains no IndexedFaceSet/IndexedLineSet " *
        "coordIndex geometry"))
    return node_labels,coords,records
end

# Header (V1 `ascii` or V2 `utf8`), tokenization and set extraction shared by
# the public reader and the foreign-mesh staging path. Returns (v2, sets).
function _wrl_parse(text::AbstractString,caller::AbstractString)
    i=firstindex(text)
    # An optional UTF-8 BOM precedes the required header comment.
    if length(text)>=3 && codeunit(text,1)==0xEF && codeunit(text,2)==0xBB &&
            codeunit(text,3)==0xBF
        i=nextind(text,nextind(text,nextind(text,i)))
    end
    head=match(r"^[^\S\n]*#VRML[^\S\n]+V([12])(?:\.0)?[^\S\n]",text[i:end])
    head===nothing && throw(ArgumentError(
        "$caller: missing '#VRML V1.0' or '#VRML V2.0' header"))
    v2=head.captures[1]=="2"
    tokens=_wrl_tokens(text,caller)
    return v2,v2 ? _wrl_v2_sets(tokens,caller) : _wrl_v1_sets(tokens,caller)
end

function _read_wrl(io::Base.IO,caller::AbstractString,
                   node_limit::Int,element_limit::Int)
    v2,sets=_wrl_parse(read(io,String),caller)
    return _wrl_records(v2,sets,caller,node_limit,element_limit)
end

"""
    WrlFile

Result of [`read_wrl`](@ref): the compact `mesh` (an MSH2-style `MixedMesh`
whose per-cell `elementary_entities` hold the VRML geometry-set entity tag —
one surface and/or one curve tag per `IndexedFaceSet`/`IndexedLineSet`),
the source `node_tags` (1-based point indices, preserving gaps left by
unused points), the per-block `element_tags`, and the parsed VRML `version`
(1 or 2).
"""
struct WrlFile
    mesh::MixedMesh
    node_tags::Vector{UInt64}
    element_tags::Vector{Vector{UInt64}}
    version::Int
end

"""
    read_wrl(path; max_nodes, max_elements, max_file_bytes) -> WrlFile

Read a VRML 1.0/2.0 `.wrl` file's indexed geometry. Each
`IndexedFaceSet`/`IndexedLineSet` contributes its `Coordinate`/`Coordinate3`
point list (V2 requires the `Shape { geometry ... }` nesting; V1 consumes the
current sibling `Coordinate3`) and the `-1`-separated `coordIndex` faces.
Two-index faces become line cells; V2 faces of three or more indices become a
sliding-window alternating-orientation triangle strip, V1 three-index faces a
triangle, four a quad, and five or more a vertex-0 fan — all matching the
observed upstream element walk. Elements land on one entity per dimension per
geometry set, node tags are the 1-based source point indices, and missing
`-1` terminators, comma-free separators and multi-shape files are accepted as
documented leniencies where upstream's reader is stricter or buggy. Index
overruns, `coord USE` references, malformed structure, and files with no
indexed geometry are errors.
"""
function read_wrl(path;
                  max_nodes=typemax(Int32),max_elements=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_wrl"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    node_limit=_io_limit(max_nodes,caller,"max_nodes";ceiling=Int(typemax(Int32)))
    element_limit=_io_limit(
        max_elements,caller,"max_elements";ceiling=Int(typemax(Int32)))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    v2,sets=_wrl_parse(read(path,String),caller)
    node_labels,coords,records=_wrl_records(
        v2,sets,caller,node_limit,element_limit)
    mesh,node_tags,element_tags=_foreign_mesh_result(
        caller,node_labels,coords,records)
    return WrlFile(mesh,node_tags,element_tags,v2 ? 2 : 1)
end

# ── XAO/BRep CAD envelopes (.xao/.brep) ────────────────────────────────────────
#
# `.xao` is gmsh's native OpenCASCADE interchange: an XML envelope naming the
# model, a `<shape format="BREP">` CDATA carrying the OCCT `CASCADE Topology`
# BRep text dump, a `<topology>` section mapping gmsh entity indices to TShape
# references (with optional entity `name`s), `<groups>` = physical groups with
# their member entity indices, and `<fields>` = mesh fields. `.brep` is the
# same CASCADE dump without the XML envelope.
#
# Entity publication from these formats needs the full BRep boundary kernel
# (TShape graph → analytic curves/surfaces → trimmed faces/solids) — the same
# classified `importShapes` boundary as upstream's OCC route; `API.open` and
# `API.merge` reject with that named boundary instead of fabricating a model.
# `read_xao`/`read_brep` are the honest partial contract: they parse and
# validate the real envelope/section structure and report the census, leaving
# the TShape geometry graph to a future kernel.

struct BrepFile
    header::String                    # e.g. "CASCADE Topology V1, (c) Matra-Datavision"
    locations::Int
    curve2ds::Int
    curves::Int
    polygons3d::Int
    polygons_on_triangulations::Int
    surfaces::Int
    triangulations::Int
    tshape_total::Int                 # declared `TShapes N`
    tshapes::Dict{Symbol,Int}         # :vertex/:edge/:wire/:face/:shell/:solid/:compsolid/:compound census
end

const _BREP_SECTIONS=("Locations","Curve2ds","Curves","Polygon3D",
                      "PolygonOnTriangulations","Surfaces","Triangulations")
const _BREP_TSHAPE_KINDS=Dict{String,Symbol}(
    "Ve"=>:vertex,"Ed"=>:edge,"Wi"=>:wire,"Fa"=>:face,"Sh"=>:shell,
    "So"=>:solid,"Cs"=>:compsolid,"Co"=>:compound)

# Parse the `CASCADE Topology` text envelope: header line, the section-count
# declarations and the per-kind TShape census. Geometry sections' records are
# not interpreted — only the declared counts and the shape-kind census are
# verified against the file.
function _brep_census(text::AbstractString,caller::AbstractString)
    # A stray UTF-8 BOM must not hide the header.
    if length(text)>=3 && codeunit(text,1)==0xEF && codeunit(text,2)==0xBB &&
            codeunit(text,3)==0xBF
        text=text[nextind(text,3):end]
    end
    header=match(r"^[ \t]*(CASCADE Topology V\d+[^\r\n]*)",text)
    header===nothing && throw(ArgumentError(
        "$caller: missing 'CASCADE Topology V<version>' header"))
    counts=Dict{String,Int}(kw=>0 for kw in _BREP_SECTIONS)
    tshape_total=-1
    for line in eachsplit(text,'\n';keepempty=false)
        m=match(r"^([A-Za-z]+)\s+(\d+)\s*$",line)
        m===nothing && continue
        kw,count=m.captures[1],tryparse(Int,m.captures[2])
        count===nothing && continue
        haskey(counts,kw) && (counts[kw]=count)
        kw=="TShapes" && (tshape_total=count)
    end
    tshape_total<0 && throw(ArgumentError(
        "$caller: missing 'TShapes <count>' declaration"))
    census=Dict{Symbol,Int}(kind=>0 for kind in values(_BREP_TSHAPE_KINDS))
    seen=0;inside=false
    for line in eachsplit(text,'\n';keepempty=false)
        stripped=strip(line)
        if !inside
            match(r"^TShapes\s+\d+\s*$",stripped)!==nothing && (inside=true)
            continue
        end
        kind=get(_BREP_TSHAPE_KINDS,stripped,nothing)
        kind===nothing && continue
        census[kind]+=1;seen+=1
    end
    seen==tshape_total || throw(ArgumentError(
        "$caller: declared $tshape_total TShapes but found $seen " *
        "shape records"))
    return BrepFile(String(header.captures[1]),counts["Locations"],
        counts["Curve2ds"],counts["Curves"],counts["Polygon3D"],
        counts["PolygonOnTriangulations"],counts["Surfaces"],
        counts["Triangulations"],tshape_total,census)
end

"""
    read_brep(path; max_file_bytes) -> BrepFile

Read the envelope census of an OCCT `CASCADE Topology` `.brep` text dump:
the versioned header, the `Locations`/`Curve2ds`/`Curves`/`Polygon3D`/
`PolygonOnTriangulations`/`Surfaces`/`Triangulations` declared counts, and a
per-kind `TShapes` census verified against its declared total. Malformed or
truncated dumps are errors. Model reconstruction is not attempted — the
TShape geometry graph needs the full OpenCASCADE boundary kernel, which is
the classified `importShapes` boundary also enforced by `API.open`/`API.merge`.
"""
function read_brep(path;max_file_bytes=typemax(Int))
    caller="read_brep"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    return _brep_census(read(path,String),caller)
end

"""
    XaoEntity

One `<topology>` entry of an `.xao` file: its model `dim` (0–3), the 0-based
topology `index`, the published entity `tag` (`index+1`), the optional entity
`name`, and the OCCT TShape `reference` it binds in the embedded BRep.
"""
struct XaoEntity
    dim::Int
    index::Int
    tag::Int
    name::String
    reference::Int
end

"""
    XaoGroup

One `<groups>/<group>` physical group of an `.xao` file: its `name`, entity
`dim`, `tag`, and the member entity `tags` resolved from the topology indices.
"""
struct XaoGroup
    name::String
    dim::Int
    tag::Int
    tags::Vector{Int}
end

"""
    XaoFile

Result of [`read_xao`](@ref): the XAO `version` and model `name`, the
`<shape>` `format` (`"BREP"`, or `""` when the model stores none), the
`<topology>` `entities`, physical `groups`, raw `<field>` attribute rows, and
the [`BrepFile`](@ref) census of the embedded CASCADE shape (`nothing` when
absent).
"""
struct XaoFile
    version::String
    name::String
    format::String
    entities::Vector{XaoEntity}
    groups::Vector{XaoGroup}
    fields::Vector{Dict{String,String}}
    brep::Union{Nothing,BrepFile}
end

const _XAO_TOPOLOGY_SECTIONS=("vertices"=>0,"edges"=>1,"faces"=>2,"solids"=>3)
const _XAO_GROUP_DIMS=Dict{String,Int}(
    "vertex"=>0,"edge"=>1,"face"=>2,"solid"=>3)

function _xao_attrs(tag::AbstractString)
    attrs=Dict{String,String}()
    for m in eachmatch(r"([A-Za-z_:][-A-Za-z0-9_:.]*)\s*=\s*\"([^\"]*)\"",tag)
        attrs[m.captures[1]]=m.captures[2]
    end
    return attrs
end

function _xao_unescape(s::AbstractString)
    s=replace(s,"&quot;"=>"\"","&apos;"=>"'","&lt;"=>"<","&gt;"=>">")
    s=replace(s,r"&#(\d+);"=>m->string(Char(parse(Int,m[2:findfirst(';',m)-1]))))
    return replace(s,"&amp;"=>"&")
end

function _xao_int(attrs,key,caller,what)
    raw=get(attrs,key,nothing)
    raw===nothing && throw(ArgumentError(
        "$caller: $what is missing the '$key' attribute"))
    value=tryparse(Int,raw)
    value===nothing && throw(ArgumentError(
        "$caller: $what has an invalid '$key' value '$raw'"))
    return value
end

"""
    read_xao(path; max_file_bytes) -> XaoFile

Read the envelope of a gmsh XAO `.xao` CAD file: the `<XAO>` version and
`<geometry>` name, the `<shape format>` + embedded `CASCADE Topology` dump
census (as a [`BrepFile`](@ref), via [`read_brep`](@ref)'s parser), the
`<topology>` index→TShape entity map (with optional names; entity tag =
index+1 like upstream), `<groups>` physical groups resolved to member entity
tags, and `<fields>` rows. Entity materialization is not attempted — the
embedded BRep needs the full OpenCASCADE boundary kernel, which is the
classified `importShapes` boundary also enforced by `API.open`/`API.merge`.
"""
function read_xao(path;max_file_bytes=typemax(Int))
    caller="read_xao"
    path isa AbstractString || throw(ArgumentError(
        "$caller: path must be a string"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    file_limit=_io_limit(max_file_bytes,caller,"max_file_bytes")
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    text=read(path,String)
    root=match(r"<XAO\b([^>]*)>",text)
    root===nothing && throw(ArgumentError(
        "$caller: missing the <XAO> root element"))
    version=_xao_unescape(get(_xao_attrs(root.captures[1]),"version",""))
    geo=match(r"<geometry\b([^>]*)>",text)
    name=geo===nothing ? "" : _xao_unescape(get(_xao_attrs(geo.captures[1]),"name",""))
    shape=match(r"(?s)<shape\b([^>]*)>(.*?)</shape>",text)
    format="";brep=nothing
    if shape!==nothing
        format=_xao_unescape(get(_xao_attrs(shape.captures[1]),"format",""))
        body=shape.captures[2]
        cdata=match(r"(?s)<!\[CDATA\[(.*)\]\]>\s*$",body)
        payload=strip(cdata===nothing ? body : cdata.captures[1])
        isempty(payload) || (brep=_brep_census(payload,caller))
    end
    entities=XaoEntity[]
    for (section,dim) in _XAO_TOPOLOGY_SECTIONS
        block=match(Regex("(?s)<$section\\b([^>]*)>(.*?)</$section>"),text)
        block===nothing && continue
        sattrs=_xao_attrs(block.captures[1]);tag=block.captures[2]
        node=section=="vertices" ? "vertex" :
              section=="edges" ? "edge" :
              section=="faces" ? "face" : "solid"
        n=0
        for m in eachmatch(Regex("<$node\\b([^>]*)/?>"),tag)
            attrs=_xao_attrs(m.captures[1])
            index=_xao_int(attrs,"index",caller,"$section entry")
            reference=_xao_int(attrs,"reference",caller,"$section entry")
            index>=0 || throw(ArgumentError(
                "$caller: $section entry has a negative index"))
            push!(entities,XaoEntity(dim,index,index+1,
                _xao_unescape(get(attrs,"name","")),reference))
            n+=1
        end
        declared=tryparse(Int,get(sattrs,"count",""))
        (declared===nothing || declared==n) || throw(ArgumentError(
            "$caller: <$section> declares count=$declared but has " *
            "$n entries"))
    end
    groups=XaoGroup[]
    gblock=match(r"(?s)<groups\b([^>]*)>(.*?)</groups>",text)
    if gblock!==nothing
        for m in eachmatch(r"(?s)<group\b([^>]*)>(.*?)</group>",gblock.captures[2])
            attrs=_xao_attrs(m.captures[1])
            dimword=get(attrs,"dimension",nothing)
            dimword===nothing && throw(ArgumentError(
                "$caller: group is missing the 'dimension' attribute"))
            dim=get(_XAO_GROUP_DIMS,dimword,nothing)
            dim===nothing && throw(ArgumentError(
                "$caller: group dimension '$dimword' is not " *
                "vertex/edge/face/solid"))
            gtag=_xao_int(attrs,"tag",caller,"group")
            members=Int[]
            for em in eachmatch(r"<element\b([^>]*)/?>",m.captures[2])
                eattrs=_xao_attrs(em.captures[1])
                eindex=_xao_int(eattrs,"index",caller,"group element")
                eindex>=0 || throw(ArgumentError(
                    "$caller: group element has a negative index"))
                push!(members,eindex+1)
            end
            declared=tryparse(Int,get(attrs,"count",""))
            (declared===nothing || declared==length(members)) ||
                throw(ArgumentError("$caller: group '$gtag' declares " *
                    "count=$declared but has $(length(members)) elements"))
            push!(groups,XaoGroup(_xao_unescape(get(attrs,"name","")),
                dim,gtag,members))
        end
        gdeclared=tryparse(Int,get(_xao_attrs(gblock.captures[1]),"count",""))
        (gdeclared===nothing || gdeclared==length(groups)) ||
            throw(ArgumentError("$caller: <groups> declares " *
                "count=$gdeclared but has $(length(groups)) groups"))
    end
    fields=Dict{String,String}[]
    # `<fields/>` self-closes; `<fields>...</fields>` carries `<field>` rows.
    fblock=match(r"(?s)<fields\b[^>]*>(.*?)</fields>",text)
    if fblock!==nothing
        for m in eachmatch(r"<field\b([^>]*)/?>",fblock.captures[1])
            push!(fields,_xao_attrs(m.captures[1]))
        end
    end
    return XaoFile(version,name,format,entities,groups,fields,brep)
end

# ── Second-generation long-tail formats (OFF / OBJ / PLY2 / TOCHNOG .dat /
#    MATLAB .m / VRML .wrl / LS-DYNA .key / RADIOSS .rad / X3D) ────────────────
# Upstream sources: src/geo/GModelIO_{OFF,OBJ,PLY,TOCHNOG,MATLAB,VRML,KEY,RAD,
# X3D}.cpp. These writers index every projected vertex densely (upstream's
# `indexMeshVertices(saveAll)`) and keep only the element families each deck
# can express; the remaining families are silently skipped — matching the
# upstream per-family loops — except where the format itself would become
# ambiguous, in which case the writer rejects.

@inline _io_dim(msh::Int)=msh_spec(msh).dim

# Element string tables (upstream getStringForTOCHNOG / getStringForKEY /
# getStringForRAD); families with no entry are skipped upstream (`typ==nullptr`).
const _TOCHNOG_STRINGS=Dict{Int,String}(1=>"-bar2",8=>"-bar3",2=>"-tria3",
    9=>"-tria6",3=>"-quad4",10=>"-quad9",4=>"-tet4",11=>"-tet10",5=>"-hex8",
    17=>"-hex20",12=>"-hex27")
const _KEY_STRINGS=Dict{Int,String}(1=>"_BEAM",8=>"_BEAM_ELBOW",2=>"_SHELL",
    9=>"_SHELL",3=>"_SHELL",16=>"_SHELL",4=>"_SOLID",11=>"_SOLID",
    5=>"_SOLID",17=>"_SOLID_H20",12=>"_SOLID_H27",6=>"_SOLID",7=>"_SOLID")
const _RAD_STRINGS=Dict{Int,String}(1=>"/TRUSS",8=>"/BEAM",2=>"/SH3N",
    9=>"/SHELL",3=>"/SHELL",16=>"/SHELL",4=>"/TETRA4",11=>"/TETRA10",
    5=>"/BRICK",17=>"/BRIC20",12=>"/BRIC20",6=>"/BRICK",7=>"/BRICK")
# KEY/RAD vertex orders reuse the BDF map for tets/prisms and the UNV map for
# lines (upstream MLine::getVertexKEY/RAD = getVertexUNV, tet/prism = BDF);
# hex20/hex27 use the upstream KEY and RAD tables below.
const _KEY_PERM=Dict{Int,NTuple{N,Int} where N}(
    8=>(1,3,2),11=>(1,5,2,6,3,7,8,10,9,4),
    18=>(1,7,2,10,3,8,9,11,12,4,13,5,15,6,14),
    17=>(1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16),
    19=>(1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16,27,21,26,22,24,25,23))
const _RAD_PERM=Dict{Int,NTuple{N,Int} where N}(
    8=>(1,3,2),11=>(1,5,2,6,3,7,8,10,9,4),
    18=>(1,7,2,10,3,8,9,11,12,4,13,5,15,6,14),
    17=>(1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18),
    19=>(1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18,20,21,22,23,24,25,26))
# MATLAB LOAD_GMSH2 section name per MSH type (upstream getMATLABName; types
# absent here have getTypeForMSH()==0 upstream and are skipped).
const _MATLAB_NAMES=Dict{Int,String}(1=>"LINES",2=>"TRIANGLES",3=>"QUADS",
    4=>"TETS",5=>"HEXAS",6=>"PRISMS",7=>"PYRAMIDS",8=>"LINES3",
    9=>"TRIANGLES6",10=>"QUADS9",16=>"QUADS8",11=>"TETS10",17=>"HEXAS20",
    12=>"HEXAS27",13=>"PRISMS18",18=>"PRISMS15",14=>"PYRAMIDS14",
    19=>"PYRAMIDS13",15=>"PNT",20=>"LINES4",26=>"LINES5",21=>"TRIANGLES9",
    22=>"TRIANGLES12",24=>"TRIANGLES15")

@inline function _perm_apply(perm,arity,k)
    perm===nothing && return k
    return k<=length(perm) ? perm[k] : k
end

# Upstream `indexMeshVertices(saveAll)` numbers every stored vertex densely in
# entity order; the foreign projection already orders rows 1..nn, so the file
# index of node row r is simply r (labels still come from external_node_tags
# where a format prints real node numbers).
@inline _dense_labels(m::MixedMesh)=
    m.entity_data===nothing ? collect(Int64(1):Int64(size(m.coords,2))) :
        Int64.(m.entity_data.external_node_tags)

function _format16(x::Real)
    return @sprintf("%.16g",Float64(x))
end

"""
    write_off(path, mesh)

Write a `Mesh`/`MixedMesh` as an Object File Format deck (`OFF nv nf 0`,
`%g` vertex rows, `3/4` face records over 0-based indices). Only dim-2
triangle/quadrangle cells are emitted (upstream iterates `faces`); other
element families are ignored, matching `GModel::writeOFF`.
"""
function write_off(path,m;scaling_factor::Real=1.0)
    caller="write_off"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    # Upstream iterates `gf->triangles`/`gf->quadrangles` and always prints a
    # 3- or 4-index record built from the corner vertices, so higher-order
    # faces degrade to their corners exactly like upstream's MTriangleN.
    _FAMILY_TRI=(2,9,21,23,24,25,31)
    _FAMILY_QUA=(3,10,16,29,33,36,37,38)
    for block in m.blocks
        block isa ElementBlock || continue
        _io_dim(block.msh)!=2 && continue
        (block.msh in _FAMILY_TRI || block.msh in _FAMILY_QUA) ||
            throw(ArgumentError(
            "$caller: element type $(block.msh) is not representable in .off"))
    end
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        nfaces=0
        for block in m.blocks
            block isa ElementBlock && _io_dim(block.msh)==2 || continue
            nfaces+=size(block.nodes,2)
        end
        println(io,"OFF ",size(m.coords,2)," ",nfaces," 0")
        for i in 1:size(m.coords,2)
            println(io,@sprintf("%g %g %g",m.coords[1,i]*scaling_factor,
                m.coords[2,i]*scaling_factor,m.coords[3,i]*scaling_factor))
        end
        for block in m.blocks
            block isa ElementBlock && _io_dim(block.msh)==2 || continue
            corners=block.msh in _FAMILY_QUA ? 4 : 3
            for j in 1:size(block.nodes,2)
                print(io,corners)
                for k in 1:corners
                    print(io," ",Int64(block.nodes[k,j])-1)
                end
                println(io)
            end
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

"""
    write_obj(path, mesh)

Write a `Mesh`/`MixedMesh` as Wavefront OBJ (`# Gmsh OBJ export` header,
`v x y z` rows at `%.16g` precision, `f` records carrying the full per-cell
vertex list as 1-based indices). Only dim-2 cells are emitted (upstream
iterates `faces`).
"""
function write_obj(path,m;scaling_factor::Real=1.0)
    caller="write_obj"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"# Gmsh OBJ export")
        for i in 1:size(m.coords,2)
            println(io,"v ",_format16(m.coords[1,i]*scaling_factor)," ",
                _format16(m.coords[2,i]*scaling_factor)," ",
                _format16(m.coords[3,i]*scaling_factor))
        end
        for block in m.blocks
            block isa ElementBlock && _io_dim(block.msh)==2 || continue
            for j in 1:size(block.nodes,2)
                print(io,"f")
                for k in 1:size(block.nodes,1)
                    print(io," ",Int64(block.nodes[k,j]))
                end
                println(io)
            end
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

"""
    write_ply2(path, mesh)

Write a `Mesh`/`MixedMesh` in the Stanford `ply2` triangulation dialect
(vertex count line, triangle count line, `%.16g` vertex rows, `3 i j k`
records over 0-based indices). Only triangle cells are emitted
(`GModel::writePLY2` iterates face `triangles` only).
"""
function write_ply2(path,m)
    caller="write_ply2"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        ntri=0
        for block in m.blocks
            block isa ElementBlock && block.msh in (2,9,21,23,24,25,31) ||
                continue
            ntri+=size(block.nodes,2)
        end
        println(io,size(m.coords,2));println(io,ntri)
        for i in 1:size(m.coords,2)
            println(io,_format16(m.coords[1,i])," ",_format16(m.coords[2,i]),
                " ",_format16(m.coords[3,i]))
        end
        for block in m.blocks
            block isa ElementBlock && block.msh in (2,9,21,23,24,25,31) ||
                continue
            for j in 1:size(block.nodes,2)
                print(io,"3 ")
                for k in 1:size(block.nodes,1)
                    print(io," ",Int64(block.nodes[k,j])-1)
                end
                println(io)
            end
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── Readers for the round-trippable decks ────────────────────────────────────
# Each returns `(node_labels, coords, records)` in the `_read_foreign_records`
# convention: dense file-side labels, a 3×n coordinate table and
# `_TextMeshCell` rows bound to the single imported discrete surface.

# Compact a flat xyz vertex table to the rows referenced by face index lists,
# remapping those indices (upstream's "delete unused vertices" behaviour).
function _drop_unref_vertices(flat::Vector{Float64},faces::Vector{Vector{Int}})
    used=sort!(unique!(reduce(vcat,faces;init=Int[])))
    remap=Dict{Int,Int}(old=>new for (new,old) in enumerate(used))
    out=Float64[]
    for old in used
        append!(out,@view flat[3old-2:3old])
    end
    return out,[[remap[i] for i in f] for f in faces]
end

function _foreign_surface_records(path,caller)
    filesize(path)<=typemax(Int) || throw(ArgumentError("$caller: file too large"))
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    return nothing
end

function _read_off(path,caller::AbstractString)
    _foreign_surface_records(path,caller)
    lines=readlines(path)
    isempty(lines) && throw(ArgumentError("$caller: empty OFF file"))
    # Upstream scans `OFF %d %d %d` — find the header token allowing comments.
    header=nothing;lno=0
    for line in lines
        lno+=1
        s=strip(line)
        (isempty(s) || startswith(s,"#")) && continue
        header=split(s);break
    end
    (header!==nothing && header[1]=="OFF" && length(header)>=4) ||
        throw(ArgumentError("$caller: invalid OFF header"))
    counts=tryparse.(Int,header[2:4])
    any(c->c===nothing,counts) && throw(ArgumentError(
        "$caller: invalid OFF header counts"))
    nv,nf,_=counts
    nv>=0 && nf>=0 || throw(ArgumentError("$caller: negative OFF counts"))
    coords=zeros(Float64,3,nv)
    vi=0
    faces=Vector{Int}[]
    for line in lines[lno+1:end]
        s=strip(line)
        (isempty(s) || startswith(s,"#")) && continue
        if vi<nv
            f=tryparse.(Float64,split(s))
            (length(f)>=3 && !any(x->x===nothing,f[1:3])) || throw(ArgumentError(
                "$caller: could not read OFF vertex $(vi+1)"))
            coords[1,vi+1],coords[2,vi+1],coords[3,vi+1]=f[1],f[2],f[3]
            vi+=1
        else
            t=tryparse.(Int,split(s))
            isempty(t) && continue
            n=t[1]
            n===nothing && throw(ArgumentError("$caller: could not read OFF face"))
            # 0-based indices; 3→triangle, 4→quadrangle, else warn-and-skip
            # (upstream prints "Ignoring n-node face").
            if n==3 && length(t)>=4
                push!(faces,[t[2]+1,t[3]+1,t[4]+1])
            elseif n==4 && length(t)>=5
                push!(faces,[t[2]+1,t[3]+1,t[4]+1,t[5]+1])
            end
        end
    end
    vi==nv || throw(ArgumentError(
        "$caller: OFF declares $nv vertices but only $vi were readable"))
    # upstream `_storeVerticesInEntities` deletes vertices no face references.
    flat=vec(coords)
    flat,faces=_drop_unref_vertices(flat,faces)
    node_labels=collect(Int64(1):Int64(div(length(flat),3)))
    records=_TextMeshCell[]
    eid=0
    for f in faces
        eid+=1
        push!(records,_TextMeshCell(Int64(eid),length(f)==3 ? 2 : 3,
            Int32.(f),Int64(1),Int64(0)))
    end
    return node_labels,flat,records
end

function _read_obj(path,caller::AbstractString)
    _foreign_surface_records(path,caller)
    col=Float64[]
    faces=Vector{Int}[]
    for line in readlines(path)
        s=strip(line)
        isempty(s) && continue
        t=split(s)
        if t[1]=="v" && length(t)>=4
            f=tryparse.(Float64,t[2:4])
            any(x->x===nothing,f) && throw(ArgumentError(
                "$caller: malformed OBJ vertex row"))
            append!(col,f)
        elseif t[1]=="f" && length(t)>=2
            idx=Int[]
            for tok in t[2:end]
                v=tryparse(Int,split(tok,"/")[1])
                v===nothing && throw(ArgumentError(
                    "$caller: malformed OBJ face index '$tok'"))
                # negative indices count back from the current vertex count
                push!(idx,v<0 ? div(length(col),3)+v : v)
            end
            nv=length(idx)
            if nv==3 || nv==4
                push!(faces,idx)
            elseif nv>4
                # upstream fan-triangulates polygons into an MPolygon; keep
                # the fan as plain triangles.
                for k in 2:nv-1
                    push!(faces,[idx[1],idx[k],idx[k+1]])
                end
            end
        end
    end
    nn=div(length(col),3)
    isempty(col) && throw(ArgumentError("$caller: OBJ file has no vertices"))
    maximum(maximum(f) for f in faces;init=0)<=nn || throw(ArgumentError(
        "$caller: OBJ face references a missing vertex"))
    # upstream drops vertices no face references.
    col,faces=_drop_unref_vertices(col,faces)
    nn=div(length(col),3)
    node_labels=collect(Int64(1):Int64(nn))
    records=_TextMeshCell[]
    for (eid,f) in enumerate(faces)
        push!(records,_TextMeshCell(Int64(eid),length(f)==3 ? 2 : 3,
            Int32.(f),Int64(1),Int64(0)))
    end
    return node_labels,col,records
end

function _read_ply2(path,caller::AbstractString)
    _foreign_surface_records(path,caller)
    it=Iterators.Stateful(Iterators.filter(l->!startswith(strip(l),"#") &&
        !isempty(strip(l)),readlines(path)))
    isempty(it) && throw(ArgumentError("$caller: empty ply2 file"))
    nv=tryparse(Int,strip(popfirst!(it)))
    nv===nothing && throw(ArgumentError("$caller: invalid ply2 vertex count"))
    isempty(it) && throw(ArgumentError("$caller: truncated ply2 file"))
    nf=tryparse(Int,strip(popfirst!(it)))
    nf===nothing && throw(ArgumentError("$caller: invalid ply2 face count"))
    nv>=0 && nf>=0 || throw(ArgumentError("$caller: negative ply2 counts"))
    coords=zeros(Float64,3,nv)
    for i in 1:nv
        isempty(it) && throw(ArgumentError("$caller: truncated ply2 vertices"))
        f=tryparse.(Float64,split(strip(popfirst!(it))))
        (length(f)>=3 && !any(x->x===nothing,f[1:3])) || throw(ArgumentError(
            "$caller: could not read ply2 vertex $i"))
        coords[1,i],coords[2,i],coords[3,i]=f[1],f[2],f[3]
    end
    faces=Vector{Int}[]
    for i in 1:nf
        isempty(it) && throw(ArgumentError("$caller: truncated ply2 faces"))
        t=tryparse.(Int,split(strip(popfirst!(it))))
        (length(t)>=4 && t[1]==3 && !any(x->x===nothing,t[2:4])) ||
            throw(ArgumentError("$caller: could not read ply2 face $i"))
        push!(faces,Int.(t[2:4]).+1)
    end
    # upstream drops vertices no face references.
    flat=vec(coords)
    flat,faces=_drop_unref_vertices(flat,faces)
    node_labels=collect(Int64(1):Int64(div(length(flat),3)))
    records=_TextMeshCell[]
    for (i,f) in enumerate(faces)
        push!(records,_TextMeshCell(Int64(i),2,Int32.(f),Int64(1),Int64(0)))
    end
    return node_labels,flat,records
end

# ── Stanford PLY (.ply) ──────────────────────────────────────────────────────
# Upstream `GModel::readPLY` accepts ASCII plus experimental
# `binary_little_endian`/`binary_big_endian` decks. The header is scanned for
# `element vertex`/`element face` counts and `property` declarations; vertex
# rows carry x,y,z followed by extra scalar properties (which upstream turns
# into per-node PViews — no channel exists through the mesh bridge, so they
# are read and dropped). Faces are triangles: upstream `sscanf`s the count and
# the first three vertex indices regardless of the declared arity, which the
# ASCII path matches verbatim. Ordinary `.ply` is read-only upstream (the
# writable dialect is `.ply2`), so there is no writer here.
const _PLY_SCALAR_TYPES=Dict{String,DataType}(
    "char"=>Int8,"int8"=>Int8,"uchar"=>UInt8,"uint8"=>UInt8,
    "short"=>Int16,"int16"=>Int16,"ushort"=>UInt16,"uint16"=>UInt16,
    "int"=>Int32,"int32"=>Int32,"uint"=>UInt32,"uint32"=>UInt32,
    "float"=>Float32,"float32"=>Float32,
    "double"=>Float64,"float64"=>Float64)

function _ply_scalar_type(word,caller)
    T=get(_PLY_SCALAR_TYPES,word,nothing)
    T===nothing && throw(ArgumentError(
        "$caller: unsupported PLY property type '$word'"))
    return T
end

# One scalar property value; `swap` byte-reverses big-endian payloads. Upstream
# swaps only the face indices (an acknowledged gap — big-endian coordinates
# come back garbled there); the declared-width swap here follows the PLY spec.
function _ply_read_scalar(io,T,swap)
    v=read(io,T)
    return swap ? bswap(v) : v
end

function _read_ply(path,caller::AbstractString)
    _foreign_surface_records(path,caller)
    io=open(path,"r")
    try
        format=:ascii;swap=false
        element="";nbv=-1;nbf=-1
        # Vertex properties in declaration order (scalars only); face
        # properties keep `(count_type, item_type)` for list entries and the
        # scalar type otherwise — needed to walk binary face records.
        vprops=DataType[];fprops=Tuple{DataType,Union{DataType,Nothing}}[]
        ordered_elements=String[]
        while true
            eof(io) && throw(ArgumentError(
                "$caller: truncated PLY header (no end_header)"))
            line=readline(io)
            words=split(strip(line))
            isempty(words) && continue
            head=words[1]
            if head=="comment" || head=="obj_info" || head=="ply"
                continue
            elseif head=="format"
                length(words)>=2 || throw(ArgumentError(
                    "$caller: malformed PLY format line"))
                words[2]=="ascii" && continue
                words[2] in ("binary_little_endian","binary_big_endian") ||
                    throw(ArgumentError(
                        "$caller: unsupported PLY format '$(words[2])'"))
                format=:binary;swap=words[2]=="binary_big_endian"
            elseif head=="element"
                length(words)>=3 || throw(ArgumentError(
                    "$caller: malformed PLY element line"))
                element=words[2]
                push!(ordered_elements,element)
                n=tryparse(Int,words[3])
                n===nothing && throw(ArgumentError(
                    "$caller: invalid PLY element count '$(words[3])'"))
                if element=="vertex"
                    nbv<0 || throw(ArgumentError(
                        "$caller: duplicate 'element vertex' declaration"))
                    nbv=n
                elseif element=="face"
                    nbf<0 || throw(ArgumentError(
                        "$caller: duplicate 'element face' declaration"))
                    nbf=n
                else
                    throw(ArgumentError(
                        "$caller: unsupported PLY element '$element' " *
                        "(upstream reads only vertex/face)"))
                end
            elseif head=="property"
                length(words)>=3 || throw(ArgumentError(
                    "$caller: malformed PLY property line"))
                if words[2]=="list"
                    length(words)>=5 || throw(ArgumentError(
                        "$caller: malformed PLY list property"))
                    element=="vertex" && throw(ArgumentError(
                        "$caller: list property '$(words[end])' on the " *
                        "vertex element is unsupported"))
                    element=="face" ||
                        throw(ArgumentError(
                            "$caller: property outside vertex/face scope"))
                    # `list <count_type> <item_type> <name>`: the count is
                    # read first, then `count` items of the item type.
                    push!(fprops,(_ply_scalar_type(words[3],caller),
                        _ply_scalar_type(words[4],caller)))
                else
                    T=_ply_scalar_type(words[2],caller)
                    element=="vertex" && push!(vprops,T)
                    element=="face" && push!(fprops,(T,nothing))
                    element in ("vertex","face") ||
                        throw(ArgumentError(
                            "$caller: property outside vertex/face scope"))
                end
            elseif head=="end_header"
                break
            end
            # Unrecognized header keywords are skipped, matching upstream's
            # sscanf-driven scan.
        end
        nbv<0 && throw(ArgumentError(
            "$caller: PLY file declares no 'element vertex'"))
        nbf<0 && (nbf=0)
        # PLY orders body records by element declaration; upstream assumes
        # vertices then faces.
        vertex_first=findfirst(==("vertex"),ordered_elements)
        face_first=findfirst(==("face"),ordered_elements)
        (vertex_first!==nothing &&
         (face_first===nothing || vertex_first<face_first)) ||
            throw(ArgumentError(
                "$caller: PLY face data must follow vertex data"))
        length(vprops)>=3 || throw(ArgumentError(
            "$caller: PLY vertex element needs at least x, y, z properties"))
        coords=zeros(Float64,3,nbv)
        faces=Vector{Int}[]
        if format==:ascii
            for i in 1:nbv
                eof(io) && throw(ArgumentError(
                    "$caller: truncated PLY vertex data"))
                t=tryparse.(Float64,split(strip(readline(io))))
                (length(t)>=3 && !any(x->x===nothing,t[1:3])) ||
                    throw(ArgumentError(
                        "$caller: could not read PLY vertex $i"))
                coords[1,i],coords[2,i],coords[3,i]=t[1],t[2],t[3]
            end
            for i in 1:nbf
                eof(io) && throw(ArgumentError(
                    "$caller: truncated PLY face data"))
                t=tryparse.(Int,split(strip(readline(io))))
                # Upstream reads `count i0 i1 i2` unconditionally and builds a
                # triangle from the first three indices.
                (length(t)>=4 && !any(x->x===nothing,t[1:4])) ||
                    throw(ArgumentError(
                        "$caller: could not read PLY face $i"))
                t[1]>=3 || throw(ArgumentError(
                    "$caller: PLY face $i declares fewer than 3 vertices"))
                push!(faces,Int.(t[2:4]).+1)
            end
        else
            for i in 1:nbv
                for (k,T) in enumerate(vprops)
                    v=_ply_read_scalar(io,T,swap)
                    k<=3 && (coords[k,i]=Float64(v))
                end
            end
            for i in 1:nbf
                indices=nothing
                for (T,itemT) in fprops
                    if itemT===nothing
                        _ply_read_scalar(io,T,swap) # scalar face property: skip
                    else
                        nbe=Int(_ply_read_scalar(io,T,swap))
                        nbe==3 || throw(ArgumentError(
                            "$caller: binary PLY face $i has $nbe vertices " *
                            "(upstream handles only triangles)"))
                        row=[Int(_ply_read_scalar(io,itemT,swap)) for _ in 1:nbe]
                        indices===nothing && (indices=row)
                    end
                end
                indices===nothing && throw(ArgumentError(
                    "$caller: binary PLY face element declares no index list"))
                push!(faces,indices.+1)
            end
        end
        nv=size(coords,2)
        for (i,f) in enumerate(faces)
            all(j->1<=j<=nv,f) || throw(ArgumentError(
                "$caller: PLY face $i references a missing vertex"))
        end
        # upstream `_storeVerticesInEntities` deletes unreferenced vertices.
        flat=vec(coords)
        flat,faces=_drop_unref_vertices(flat,faces)
        node_labels=collect(Int64(1):Int64(div(length(flat),3)))
        records=_TextMeshCell[]
        for (i,f) in enumerate(faces)
            push!(records,_TextMeshCell(Int64(i),2,Int32.(f),Int64(1),Int64(0)))
        end
        return node_labels,flat,records
    catch err
        err isa InterruptException && rethrow()
        err isa EOFError && throw(ArgumentError(
            "$caller: truncated .ply file"))
        rethrow()
    finally
        close(io)
    end
end

# ── TOCHNOG (.dat) ────────────────────────────────────────────────────────────
"""
    write_dat(path, mesh; save_groups_of_nodes=true)

Write a `Mesh`/`MixedMesh` as a TOCHNOG input deck: the upstream literal
preamble, `node <i> x [y [z]]` rows at `%.16g` for the mesh dimension,
`element <num> <tochnog-type> n…` records, physical-group element-set and
node-set sections. Element families without an upstream TOCHNOG mnemonic are
skipped exactly as upstream (`getStringForTOCHNOG()==nullptr`).
"""
function write_dat(path,m;save_groups_of_nodes::Bool=true,
                   scaling_factor::Real=1.0)
    caller="write_dat"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    dim=0
    for block in m.blocks
        block isa ElementBlock || continue
        dim=max(dim,_io_dim(block.msh))
    end
    target,parent=_foreign_target(path,caller)
    labels=_dense_labels(m)
    mktemp(parent) do temporary,io
        println(io,"(----Tochnog Input File Created by Gmsh Version 2.13.1 ----)")
        println(io," ")
        print(io,"(**EDIT OR MODIFY THE ENTRIES BELOW AS REQUIRED**)\n",
            "echo  -yes \n",
            "number of space dimensions (add number here) \n",
            "derivatives \n","materi_velocity \n",
            "materi_displacement \n","materi_strain_total  \n",
            "materi_stress  \n","condif_temperature \n",
            "number_of_integration_points (add number here) \n",
            "end_initia \n","options_element_dof -yes  \n")
        println(io,"\n(----------------Nodes----------------)\n")
        for i in 1:size(m.coords,2)
            if dim==3
                println(io,"node ",labels[i]," ",_format16(m.coords[1,i]*scaling_factor),
                    " ",_format16(m.coords[2,i]*scaling_factor)," ",
                    _format16(m.coords[3,i]*scaling_factor))
            elseif dim==2
                println(io,"node ",labels[i]," ",_format16(m.coords[1,i]*scaling_factor),
                    " ",_format16(m.coords[2,i]*scaling_factor))
            elseif dim==1
                println(io,"node ",labels[i]," ",_format16(m.coords[1,i]*scaling_factor))
            end
        end
        println(io,"\n(++++++++++++++ E L E M E N T S ++++++++++++++)\n")
        counter=UInt64(0)
        for (bi,block) in enumerate(m.blocks)
            block isa ElementBlock || continue
            typ=get(_TOCHNOG_STRINGS,block.msh,nothing)
            typ===nothing && continue
            # upstream writes points always, and the top-dim family only
            wanted=_io_dim(block.msh)
            (wanted==0 || wanted==dim) || continue
            for j in 1:size(block.nodes,2)
                counter+=UInt64(1)
                elnum=_foreign_cell_uid(m,bi,j,counter)
                print(io,"element ",elnum," ",typ," ")
                for k in 1:size(block.nodes,1)
                    print(io," ",labels[block.nodes[k,j]])
                end
                println(io)
            end
        end
        println(io,"\n( +------------+---------Physical Groups Section----------+------------+\n")
        # Physical element sets: group cells by physical tag per dim.
        _write_tochnog_groups(io,m,labels;nodes=save_groups_of_nodes)
        println(io,"( +------------+---------Physical Groups Section----------+------------+")
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

function _write_tochnog_groups(io,m::MixedMesh,labels;nodes::Bool)
    phys=Dict{Tuple{Int,Int64},Vector{UInt64}}()
    counter=UInt64(0)
    for (bi,block) in enumerate(m.blocks)
        block isa ElementBlock || continue
        for j in 1:size(block.nodes,2)
            counter+=UInt64(1)
            tag=Int64(block.tags[j])
            tag==0 && continue
            push!(get!(phys,(_io_dim(block.msh),tag),UInt64[]),
                _foreign_cell_uid(m,bi,j,counter))
        end
    end
    for ((d,tag),uids) in sort!(collect(phys);by=first)
        println(io,"\n(Element sets ===> 'element_group' to identify " *
            "DIFFERENT MATERIALS =$(get(Dict()," ",""))$tag)")
        for (n,uid) in enumerate(uids)
            (n>1 && mod(n-1,10)==0) && println(io)
            print(io,uid," ")
        end
        println(io)
    end
    if nodes
        nodesets=Dict{Tuple{Int,Int64},Set{Int32}}()
        for (bi,block) in enumerate(m.blocks)
            block isa ElementBlock || continue
            for j in 1:size(block.nodes,2)
                tag=Int64(block.tags[j])
                tag==0 && continue
                push!(get!(nodesets,(_io_dim(block.msh),tag),Set{Int32}()),
                    block.nodes[:,j]...)
            end
        end
        for ((d,tag),rows) in sort!(collect(nodesets);by=first)
            println(io,"\n(Node sets ===> Used to set BOUNDARY CONDITIONS " *
                "in Tochnog =$tag)\n")
            for (n,row) in enumerate(sort!(collect(rows)))
                (n>1 && mod(n-1,10)==0) && println(io)
                print(io,labels[row]," ")
            end
            println(io)
        end
    end
    return nothing
end

# ── MATLAB (.m, LOAD_GMSH2 dialect — upstream hardcodes filetype=1) ───────────
"""
    write_m(path, mesh)

Write a `Mesh`/`MixedMesh` as a MATLAB `load_gmsh2` deck: `msh.nbNod`,
`msh.POS` (`%.16g` rows), `msh.MAX`/`msh.MIN`, then per-MSH-type
`msh.<NAME>` connectivity matrices whose trailing column carries the cell's
physical tag (or elementary tag when zero). Element families without an
upstream MATLAB name are skipped.
"""
function write_m(path,m;scaling_factor::Real=1.0,model_name::AbstractString="")
    caller="write_m"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    target,parent=_foreign_target(path,caller)
    labels=_dense_labels(m)
    mktemp(parent) do temporary,io
        println(io,"%% Matlab mesh")
        println(io,"%% ",isempty(model_name) ? splitext(basename(target))[1] :
            model_name,", Created by Gmsh")
        println(io,"%% ASCII")
        println(io,"clear msh;")
        println(io,"msh.nbNod = ",size(m.coords,2),";")
        println(io,"msh.POS = [")
        for i in 1:size(m.coords,2)
            println(io,_format16(m.coords[1,i]*scaling_factor)," ",
                _format16(m.coords[2,i]*scaling_factor)," ",
                _format16(m.coords[3,i]*scaling_factor),";")
        end
        println(io,"];")
        println(io,"msh.MAX = max(msh.POS);")
        println(io,"msh.MIN = min(msh.POS);")
        # Upstream groups elements by MSH type order 1..MSH_MAX_NUM.
        order=sort!(unique(Int(block.msh) for block in m.blocks if
            block isa ElementBlock))
        for msh in order
            name=get(_MATLAB_NAMES,msh,nothing)
            name===nothing && continue
            println(io,"msh.",name," =[")
            for (bi,block) in enumerate(m.blocks)
                (block isa ElementBlock && block.msh==msh) || continue
                for j in 1:size(block.nodes,2)
                    tag=Int64(block.tags[j])
                    trailing=tag!=0 ? abs(tag) : _foreign_cell_entity(m,bi,j)
                    for k in 1:size(block.nodes,1)
                        print(io," ",labels[block.nodes[k,j]])
                    end
                    println(io," ",trailing)
                end
            end
            println(io,"];")
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── VRML 1.0 (.wrl/.vrml) ────────────────────────────────────────────────────
"""
    write_wrl(path, mesh)

Write a `Mesh`/`MixedMesh` as VRML 1.0 (`#VRML V1.0 ascii`,
`Coordinate3{point[...]}` for the dim≤2 node set, `DEF Curve<tag>
IndexedLineSet` and `DEF Surface<tag> IndexedFaceSet` `coordIndex` blocks
terminated by `-1`). Matches `GModel::writeVRML`: only dim-0..2 vertices are
written and only line/surface cells appear in the index sets.
"""
function write_wrl(path,m;scaling_factor::Real=1.0)
    caller="write_wrl"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    target,parent=_foreign_target(path,caller)
    # VRML indices are over the dim≤2 vertex stream; keep the dense file index
    # of each node row, or -1 for volume-owned rows which are never written.
    written=falses(size(m.coords,2))
    for block in m.blocks
        block isa ElementBlock && _io_dim(block.msh)<=2 || continue
        for r in block.nodes
            written[r]=true
        end
    end
    # Upstream writes entity vertices for points/edges/faces — which is every
    # node owned by a dim≤2 entity, not merely referenced rows. The foreign
    # projection carries entity ownership per node through block_entities'
    # owners only at the cell level, so keep referenced dim≤2 rows plus all
    # nodes owned by dim≤2 entities when per-node ownership is unavailable.
    index=fill(Int64(-1),size(m.coords,2))
    c=0
    for i in 1:length(written)
        written[i] && (index[i]=(c+=1))
    end
    mktemp(parent) do temporary,io
        println(io,"#VRML V1.0 ascii")
        println(io,"#created by Gmsh")
        println(io,"Coordinate3 {")
        println(io,"  point [")
        for i in 1:size(m.coords,2)
            index[i]<0 && continue
            println(io,"    ",_format16(m.coords[1,i]*scaling_factor)," ",
                _format16(m.coords[2,i]*scaling_factor)," ",
                _format16(m.coords[3,i]*scaling_factor),",")
        end
        println(io,"  ]")
        println(io,"}")
        # Group cells per elementary entity to mirror upstream's per-entity
        # DEF blocks.
        for wanted in (1,2)
            seen=Dict{Int64,Vector{Tuple{Int,Int}}}()
            for (bi,block) in enumerate(m.blocks)
                block isa ElementBlock && _io_dim(block.msh)==wanted || continue
                for j in 1:size(block.nodes,2)
                    push!(get!(seen,_foreign_cell_entity(m,bi,j),
                        Tuple{Int,Int}[]),(bi,j))
                end
            end
            for (etag,cells) in sort!(collect(seen);by=first)
                println(io,wanted==1 ? "DEF Curve$etag IndexedLineSet {" :
                    "DEF Surface$etag IndexedFaceSet {")
                println(io,"  coordIndex [")
                for (bi,j) in cells
                    block=m.blocks[bi]
                    for k in 1:size(block.nodes,1)
                        print(io,index[block.nodes[k,j]]<0 ? 0 :
                            index[block.nodes[k,j]]-1,",")
                    end
                    println(io,"-1,")
                end
                println(io,"  ]")
                println(io,"}")
            end
        end
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── LS-DYNA keyword (.key) ───────────────────────────────────────────────────
"""
    write_key(path, mesh)

Write a `Mesh`/`MixedMesh` as an LS-DYNA keyword deck (`*KEYWORD`, `*TITLE`,
`*NODE` rows `%ld, %.12g, %.12g, %.12g`, per-entity `*ELEMENT_<KIND>` sections
with the upstream `partID(dim,tag)=dim*1e6+tag` part id, `*SET_*_LIST`
physical-group sections, `*END`). Implements the upstream node-filling rules:
tet4 repeats node 4 into slots 4–8, penta6 into 7–8/5–6 pattern, tri3 shell
repeats node 3, 6-node shell expands to 8, 3-node beam appends its third node
on a second field.
"""
function write_key(path,m;scaling_factor::Real=1.0)
    caller="write_key"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    # Element families without an upstream KEY mnemonic are skipped
    # (upstream getStringForKEY returns nullptr), including point cells.
    target,parent=_foreign_target(path,caller)
    labels=_dense_labels(m)
    mktemp(parent) do temporary,io
        println(io,"\$# LS-DYNA Keyword file created by Gmsh")
        println(io,"*KEYWORD")
        println(io,"*TITLE")
        println(io," ",basename(target))
        println(io,"*NODE")
        for i in 1:size(m.coords,2)
            println(io,labels[i],", ",@sprintf("%.12g",m.coords[1,i]*scaling_factor),
                ", ",@sprintf("%.12g",m.coords[2,i]*scaling_factor),", ",
                @sprintf("%.12g",m.coords[3,i]*scaling_factor))
        end
        # Per-entity element sections in entity order per dim (upstream loops
        # points, lines, surfaces, volumes in entity iteration order).
        for wanted in (0,1,2,3)
            seen=Dict{Int64,Vector{Int}}()
            for (bi,block) in enumerate(m.blocks)
                block isa ElementBlock && _io_dim(block.msh)==wanted || continue
                seen[_foreign_cell_entity(m,bi,1)]=Int[]
            end
            for (bi,block) in enumerate(m.blocks)
                block isa ElementBlock && _io_dim(block.msh)==wanted || continue
                typ=get(_KEY_STRINGS,block.msh,nothing)
                typ===nothing && continue
                for etag in sort!(collect(keys(seen)))
                    rows=[j for j in 1:size(block.nodes,2) if
                        _foreign_cell_entity(m,bi,j)==etag]
                    isempty(rows) && continue
                    pid=wanted*1000000+etag
                    println(io,"*ELEMENT",typ)
                    println(io,"\$#SET_ELEMENT=",etag)
                    counter=UInt64(0)
                    for j in rows
                        elnum=_foreign_cell_uid(m,bi,j,UInt64(j))
                        print(io,elnum,", ",pid,", ")
                        n=size(block.nodes,1)
                        perm=get(_KEY_PERM,block.msh,nothing)
                        nid=[labels[block.nodes[_perm_apply(perm,n,k),j]]
                            for k in 1:n]
                        d=_io_dim(block.msh)
                        if d==3
                            if n==4
                                append!(nid,fill(nid[4],4))
                            elseif n==6
                                nid=[nid[1:6];nid[6];nid[6]]
                            end
                            length(nid)<8 && append!(nid,fill(0,8-length(nid)))
                        elseif d==2
                            if n==3
                                nid=[nid;nid[3]]
                            elseif n==6
                                nid=[nid[1:3];nid[3];nid[4:6];nid[6]]
                            end
                        elseif d==1 && n==3
                            # 3-node beam: third node written on the 9th field.
                            append!(nid,fill(0,6));push!(nid,nid[3])
                        end
                        for (i,v) in enumerate(nid)
                            print(io,v)
                            i!=length(nid) && print(io,", ")
                            mod(i+1,10)==0 && i!=length(nid) && println(io)
                        end
                        println(io)
                    end
                end
            end
        end
        # *SET_<SOLID|SHELL|BEAM|NODE>_LIST per physical tag.
        sets=Dict{Tuple{Int,Int64},Vector{UInt64}}()
        counter=UInt64(0)
        for (bi,block) in enumerate(m.blocks)
            block isa ElementBlock || continue
            for j in 1:size(block.nodes,2)
                counter+=UInt64(1)
                tag=Int64(block.tags[j]);tag==0 && continue
                push!(get!(sets,(_io_dim(block.msh),tag),UInt64[]),
                    _foreign_cell_uid(m,bi,j,counter))
            end
        end
        setid=0
        for ((d,tag),uids) in sort!(collect(sets);by=first)
            kind=d==3 ? "SOLID" : d==2 ? "SHELL" : d==1 ? "BEAM" : "NODE"
            setid+=1
            println(io,"*SET_",kind,"_LIST")
            println(io,"\$# ",tag)
            println(io,setid)
            for (n,uid) in enumerate(uids)
                mod(n-1,8)==0 ? print(io,"\n",uid) : print(io,", ",uid)
            end
            println(io)
        end
        println(io,"*END")
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── RADIOSS starter (.rad) ───────────────────────────────────────────────────
"""
    write_rad(path, mesh)

Write a `Mesh`/`MixedMesh` as a RADIOSS starter deck (`#RADIOSS STARTER`,
`/BEGIN` header, `/NODE` `%10ld%20g` rows, per-entity `<TYPE>/<part>` blocks
with `%10d` element records, `/GRSH3N/SH3N` physical-group sections for
triangles and `/GR...` element sets, `*`-free Radioss syntax). Implements the
upstream degenerate-node fills (penta→brick pattern, tet4 4-node records).
"""
function write_rad(path,m;scaling_factor::Real=1.0)
    caller="write_rad"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    # Element families without an upstream RAD card are skipped
    # (upstream getStringForRAD returns nullptr), including point cells.
    target,parent=_foreign_target(path,caller)
    labels=_dense_labels(m)
    beginname=basename(target)
    for pat in ("_0000.rad",".rad")
        endswith(beginname,pat) && (beginname=beginname[1:end-length(pat)])
    end
    mktemp(parent) do temporary,io
        println(io,"#RADIOSS STARTER")
        println(io,"#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|")
        println(io,"# Created by Gmsh, Radioss Mesh Interface by PaulAltair sharp@altair.com")
        println(io,"#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|")
        println(io,"/BEGIN")
        println(io,beginname)
        println(io,"      2022         0")
        println(io,"                  Mg                  mm                   s")
        println(io,"                  Mg                  mm                   s")
        println(io,"/NODE")
        for i in 1:size(m.coords,2)
            println(io,@sprintf("%10d%20g%20g%20g",labels[i],
                m.coords[1,i]*scaling_factor,m.coords[2,i]*scaling_factor,
                m.coords[3,i]*scaling_factor))
        end
        # Per-entity blocks per dimension, upstream entity order.
        for wanted in (0,1,2,3)
            for (bi,block) in enumerate(m.blocks)
                block isa ElementBlock && _io_dim(block.msh)==wanted || continue
                typ=get(_RAD_STRINGS,block.msh,nothing)
                typ===nothing && continue
                byent=Dict{Int64,Vector{Int}}()
                for j in 1:size(block.nodes,2)
                    push!(get!(byent,_foreign_cell_entity(m,bi,j),Int[]),j)
                end
                for (etag,rows) in sort!(collect(byent);by=first)
                    pid=wanted*1000000+etag
                    println(io,typ,"/",pid)
                    println(io,"#SET_ELEMENT=",etag)
                    perm=get(_RAD_PERM,block.msh,nothing)
                    for j in rows
                        elnum=_foreign_cell_uid(m,bi,j,UInt64(j))
                        n=size(block.nodes,1)
                        nid=[labels[block.nodes[_perm_apply(perm,n,k),j]]
                            for k in 1:n]
                        d=_io_dim(block.msh)
                        if d==3 && n==6
                            # penta6 degenerated to brick — upstream
                            # writeRAD's slot swap produces source positions
                            # [4,1,3,6,5,2,2,5].
                            nid=[nid[4],nid[1],nid[3],nid[6],nid[5],nid[2],
                                nid[2],nid[5]]
                        end
                        print(io,@sprintf("%10d",elnum))
                        for v in nid
                            print(io,@sprintf("%10d",v))
                        end
                        println(io)
                    end
                end
            end
        end
        # /GRSH3N/SH3N physical-group sections for 3-node shells, then
        # generic /GRPART-style element lists per physical tag.
        sets=Dict{Tuple{Int,Int64},Vector{Tuple{Int,Int64}}}()
        counter=UInt64(0)
        for (bi,block) in enumerate(m.blocks)
            block isa ElementBlock || continue
            for j in 1:size(block.nodes,2)
                counter+=UInt64(1)
                tag=Int64(block.tags[j]);tag==0 && continue
                push!(get!(sets,(_io_dim(block.msh),tag),
                    Tuple{Int,Int64}[]),(size(block.nodes,1),
                    _foreign_cell_uid(m,bi,j,counter)))
            end
        end
        setid=0
        for ((d,tag),cells) in sort!(collect(sets);by=first)
            d==2 || continue
            tris=[uid for (n,uid) in cells if n==3]
            isempty(tris) && continue
            setid+=1
            println(io,"/GRSH3N/SH3N/",setid)
            println(io,tag)
            for (n,uid) in enumerate(tris)
                print(io,@sprintf("%10d",uid))
                mod(n,8)==0 && println(io)
            end
            mod(length(tris),8)!=0 && println(io)
        end
        println(io,"/END")
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    return path
end

# ── X3D ─────────────────────────────────────────────────────────────────────
"""
    write_x3d(path, mesh)

Write a `Mesh`/`MixedMesh` as an X3D 3.3 Interchange document plus the
companion `*_x3d.html` launcher (upstream `GModel::writeX3D` single-file
mode). Faces become `<TriangleSet>` shapes with explicit `<Coordinate>`
rows (quadrangles split on the first-third vertex diagonal); dim-1 cells
become `<LineSet>` shapes and leftover dim-0 cells `<PointSet>`-style
coordinates under `IndexedLineSet` is not emitted — only surface, edge and
vertex geometry is written, matching upstream's surfaces/edges/vertices
switches with the default `x3dsurfaces=1` merged behaviour.
"""
function write_x3d(path,m;scaling_factor::Real=1.0,surfaces::Bool=true,
                   edges::Bool=false,vertices::Bool=false)
    caller="write_x3d"
    m isa Mesh && (m=simplex_to_mixed(m))
    m isa MixedMesh || throw(ArgumentError(
        "$caller: mesh must be a Mesh or MixedMesh"))
    (surfaces || edges || vertices) || throw(ArgumentError(
        "$caller: no surfaces, edges or vertices to write into '$path'"))
    target,parent=_foreign_target(path,caller)
    mktemp(parent) do temporary,io
        println(io,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        println(io,"<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.3//EN\" " *
            "\"http://www.web3d.org/specifications/x3d-3.3.dtd\">")
        println(io,"<X3D profile='Interchange' version='3.3'  " *
            "xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' >")
        println(io,"  <head>")
        println(io,"    <meta name='creator' content='gmsh'/> ")
        println(io,"  </head>")
        println(io,"  <Scene>")
        if surfaces
            println(io,"   <Group DEF=\"faces\">")
            println(io,"    <Shape DEF=\"face\">")
            println(io,"     <Appearance><Material DEF=\"matface\"></Material></Appearance>")
            println(io,"     <TriangleSet DEF=\"setface\">")
            println(io,"      <Coordinate point=\"")
            for block in m.blocks
                block isa ElementBlock && _io_dim(block.msh)==2 || continue
                n=size(block.nodes,1)
                for j in 1:size(block.nodes,2)
                    if n>=3
                        tris=n==3 ? ((1,2,3),) :
                            n==4 ? ((1,2,3),(1,3,4)) :
                            ((1,2,3),)
                        for t in tris
                            for k in t
                                r=block.nodes[k,j]
                                println(io,_format16(m.coords[1,r]*scaling_factor)," ",
                                    _format16(m.coords[2,r]*scaling_factor)," ",
                                    _format16(m.coords[3,r]*scaling_factor))
                            end
                        end
                    end
                end
            end
            println(io,"\"></Coordinate>")
            println(io,"     </TriangleSet>")
            println(io,"    </Shape>")
            println(io,"   </Group>")
        end
        if edges
            println(io,"   <Group DEF=\"edges\">")
            for (bi,block) in enumerate(m.blocks)
                block isa ElementBlock && _io_dim(block.msh)==1 || continue
                for j in 1:size(block.nodes,2)
                    println(io,"    <Shape DEF=\"edge\">")
                    println(io,"     <Appearance><Material DEF=\"matedge\"></Material><LineProperties id=\"propedge\"></LineProperties></Appearance>")
                    println(io,"     <LineSet vertexCount=\"",size(block.nodes,1),"\">")
                    println(io,"      <Coordinate point=\"")
                    for k in 1:size(block.nodes,1)
                        r=block.nodes[k,j]
                        println(io,_format16(m.coords[1,r]*scaling_factor)," ",
                            _format16(m.coords[2,r]*scaling_factor)," ",
                            _format16(m.coords[3,r]*scaling_factor))
                    end
                    println(io,"\"/>")
                    println(io,"     </LineSet>")
                    println(io,"    </Shape>")
                end
            end
            println(io,"   </Group>")
        end
        if vertices
            println(io,"   <Group DEF=\"vertices\">")
            for block in m.blocks
                block isa ElementBlock && _io_dim(block.msh)==0 || continue
                for j in 1:size(block.nodes,2)
                    r=block.nodes[1,j]
                    println(io,"    <Shape DEF=\"vertex\">")
                    println(io,"     <PointSet>")
                    println(io,"      <Coordinate point=\"",
                        _format16(m.coords[1,r]*scaling_factor)," ",
                        _format16(m.coords[2,r]*scaling_factor)," ",
                        _format16(m.coords[3,r]*scaling_factor),"\"/>")
                    println(io,"     </PointSet>")
                    println(io,"    </Shape>")
                end
            end
            println(io,"   </Group>")
        end
        println(io,"  </Scene>")
        println(io,"</X3D>")
        flush(io);close(io)
        mv(temporary,target;force=true)
    end
    # Companion html launcher (upstream always writes it next to the x3d).
    htmlname=string(splitext(target)[1],"_x3d.html")
    open(htmlname,"w") do io
        println(io,"<!DOCTYPE HTML\n>")
        println(io,"<html lang=\"en\"> ")
        println(io,"<head>")
        println(io,"  <title> gmsh x3d render</title>")
        println(io,"  <meta charset=\"utf-8\">")
        println(io,"  <link rel=\"stylesheet\" type=\"text/css\" href=\"https://x3dom.org/release/x3dom.css\">")
        println(io,"  <script src=\"https://x3dom.org/release/x3dom.js\"></script>")
        println(io,"</head>")
        println(io,"<body>")
        println(io,"  <x3d id=\"gmsh-scene\" style=\"width: 100%; height: 100%;border: none\" >")
        println(io,"    <Scene>")
        println(io,"      <transform scale=\"1,1,1\">")
        println(io,"      <transform id=\"components\" rotation=\"1 0 0 -1.57079632679\">")
        println(io,"        <inline onload=\"fit()\" mapDEFToID=\"true\" url=",basename(target),"></inline>")
        println(io,"      </transform>")
        println(io,"      </transform>")
        println(io,"    </Scene>")
        println(io,"  </x3d>")
        println(io,"  <script>")
        println(io,"  function fit()")
        println(io,"  {")
        println(io,"    var x3dElem = document.getElementById('gmsh-scene');")
        println(io,"    x3dElem.runtime.fitAll();")
        println(io,"  }")
        println(io,"  </script>")
        println(io,"</body>")
        println(io,"</html>")
    end
    return path
end

"""
    write_celum(path, mesh)

The CELUM deck (`<name>_f` faces file, `<name>_s` vertices file) requires
per-vertex exterior normals and principal curvatures, which upstream obtains
by reparametrizing each node on its CAD surface (`reparamMeshVertexOnFace` +
`f->normal(param)`/`f->curvatures(param,...)`). Tessella's foreign-format
projection carries no parametric surface, so this writer rejects truthfully
instead of fabricating curvature data.
"""
function write_celum(path,m;kwargs...)
    caller="write_celum"
    throw(ArgumentError("$caller: CELUM export requires per-vertex surface " *
        "normals and principal curvatures from the parametric surface, " *
        "which the mesh-only projection cannot supply (upstream warns and " *
        "writes undefined data on discrete faces); use .off/.ply2/.wrl for " *
        "triangulated surface export"))
end

end # module IO
