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
using ..Elements: MSH_PHYSICAL_NAME_MAX_BYTES, _copy_physical_names
using ..GmshLibm: _gm_sin, _gm_cos, _gm_tan, _gm_asin, _gm_acos, _gm_atan,
                  _gm_atan2, _gm_sinh, _gm_cosh, _gm_tanh, _gm_exp, _gm_log,
                  _gm_log10, _gm_pow
using Printf: @printf, @sprintf, Format, format

export read_msh, write_msh, MshFile
export read_stl, read_geo_params, GeoParams, GeoFieldSpec

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
end

# Mutable accumulator used while parsing, converted to a `Mesh` at the end.
mutable struct _Accum
    node_x::Vector{Float64}
    node_y::Vector{Float64}
    node_z::Vector{Float64}
    tag2idx::Dict{Int,Int}          # gmsh node tag → compact 1-based index
    segs::Vector{NTuple{2,Int32}}
    tris::Vector{NTuple{3,Int32}}
    tets::Vector{NTuple{4,Int32}}
    seg_tag::Vector{Int32}
    tri_tag::Vector{Int32}
    tet_tag::Vector{Int32}
    pnames::Dict{Tuple{Int,Int},String}
    ep::Dict{Tuple{Int,Int},Int}    # (entityDim, entityTag) → physical tag (v4)
    element_tags::Set{Int}
    physical_name_records::Int
    element_blocks::Int
    _Accum() = new(Float64[], Float64[], Float64[], Dict{Int,Int}(),
                   NTuple{2,Int32}[], NTuple{3,Int32}[], NTuple{4,Int32}[],
                   Int32[], Int32[], Int32[], Dict{Tuple{Int,Int},String}(),
                   Dict{Tuple{Int,Int},Int}(),Set{Int}(),0,0)
end

@inline function _nodeidx!(acc::_Accum, tag::Int)
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
        elseif line == "\$Nodes"
            seen_nodes && throw(ArgumentError("IO: duplicate \$Nodes section"))
            seen_elements && throw(ArgumentError("IO: \$Nodes must precede \$Elements"))
            seen_nodes=true
            version != 0 || throw(ArgumentError("IO: \$Nodes appeared before \$MeshFormat"))
            version < 3 ? _read_nodes_v2!(acc, io, limits) :
                          _read_nodes_v4!(acc, io, limits)
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
            version == 4.1 || throw(ArgumentError("IO: \$Entities is only valid in MSH v4.1"))
            _read_entities_v4!(acc, io, limits)
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
    return MshFile(mesh, acc.pnames)
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
        tag = parse(Int, p[1])
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
        _record_element_tag!(acc,parse(Int,p[1]))
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
    _read_entity_block!(acc.ep, io, nP, 0)   # points
    _read_entity_block!(acc.ep, io, nC, 1)   # curves
    _read_entity_block!(acc.ep, io, nS, 2)   # surfaces
    _read_entity_block!(acc.ep, io, nV, 3)   # volumes
    _expect_end(io, "\$EndEntities")
end

function _read_entity_block!(ep, io, count, dim)
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
        for j in 1:nphys
            physj=parse(Int,p[idx+j])
            physj>0 || throw(ArgumentError("IO: physical tag on entity ($dim,$tag) must be positive"))
            _int32_tag(physj,"physical")
        end
        if nphys >= 1
            phys=parse(Int,p[idx+1])
            ep[(dim, tag)] = phys   # first physical group tag
        else
            ep[(dim,tag)]=0         # record existence, even without a physical group
        end
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
            end
        end
    end
end

function _read_nodes_v4!(acc, io, limits::_MshReadLimits)
    hdr = split(strip(readline(io)))
    length(hdr) == 4 || throw(ArgumentError("IO: malformed v4 node header"))
    numBlocks = parse(Int, hdr[1]); numNodes = parse(Int, hdr[2])
    declared_min=parse(Int,hdr[3]);declared_max=parse(Int,hdr[4])
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
    actual_min=typemax(Int);actual_max=typemin(Int)
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
        tags = Int[]
        blocktags=Set{Int}()
        for _ in 1:nInBlock
            tag=parse(Int, strip(readline(io)))
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
            if parametric==1
                all(isfinite(parse(Float64,c[j])) for j in 4:expected) ||
                    throw(ArgumentError("IO: node $(tags[i]) has non-finite parametric coordinates"))
            end
            push!(acc.node_x,x); push!(acc.node_y,y); push!(acc.node_z,z)
            acc.tag2idx[tags[i]] = length(acc.node_x)
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
    declared_min=parse(Int,hdr[3]);declared_max=parse(Int,hdr[4])
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
    actual_min=typemax(Int);actual_max=typemin(Int)
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
            elem_tag=parse(Int,p[1]);_record_element_tag!(acc,elem_tag)
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

@inline function _record_element_tag!(acc::_Accum,tag::Int)
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
        tag=parse(Int,t); i=get(acc.tag2idx,tag,0)
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
            _write_msh_v2(io, mesh, names, gmsh_compatible)
        else
            _write_msh_v4(io, mesh, names, gmsh_compatible)
        end
        flush(io); close(io)
        mv(tmp,target;force=true)
    end
    return path
end

function _validate_write_mesh(m::Mesh)
    diagnostic=validate(m)
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

function _write_msh_v2(io, m::Mesh, physical_names, gmsh_compatible::Bool)
    println(io, "\$MeshFormat"); println(io, "2.2 0 8"); println(io, "\$EndMeshFormat")
    _write_physical_names(io, physical_names, gmsh_compatible)
    # nodes
    println(io, "\$Nodes"); println(io, nnodes(m))
    @inbounds for i in 1:nnodes(m)
        p = node(m, i)
        @printf(io, "%d %.17g %.17g %.17g\n", i, p[1], p[2], p[3])
    end
    println(io, "\$EndNodes")
    # elements
    nel = nsegs(m) + ntris(m) + ntets(m)
    println(io, "\$Elements"); println(io, nel)
    eid = 0
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
    tris = NTuple{9,Float64}[]
    verts = NTuple{3,Float64}[]
    infacet=false;inloop=false
    for line in eachline(path)
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
            ratio=(Rational{BigInt}(value)-Rational{BigInt}(origin)) /
                  Rational{BigInt}(weld_tol)
            return floor(Int,ratio)
        end
        return floor(Int,(value-origin)*inv)
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
    mesh::Union{Nothing,Mesh}
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
    # Exec hook resolving `Point{t}`/`Physical X{t}` name reads in string
    # expressions: `(dim, tag, :entity|:physical) -> String`; `nothing` in pure
    # contexts.
    entity_name_lookup::Any
    # `yymsg` channels — level 0 diagnostics are errors that accumulate through
    # `gmsh_yyerrorstate` while the statement stream continues; level 1 are
    # warnings.
    exec_errors::Vector{String}
    exec_warnings::Vector{String}
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
    Dict{Tuple{String,Int,String},NTuple{4,Int}}(),0,nothing,nothing,
    String[],String[],0,false,:run,0,Dict{Int,Tuple{Int,Int}}(),
    Dict{Tuple{Int,Int},Vector{Int}}(),true,true,nothing,
    Dict{Int,Any}(),nothing,false)

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
_geo_syntax_abort(token::AbstractString,detail::AbstractString)=
    throw(_GeoSyntaxAbort(String(token),String(detail),false))

Base.showerror(io::Base.IO,err::_GeoSyntaxAbort)=
    print(io,isempty(err.detail) ? "syntax error ($(err.token))" : err.detail)

# `yymsg(0, ...)` — record a parse error; the stream continues and `execute_geo`
# throws the accumulated diagnostics at the end (or at `Exit`/`Abort`).
function _geo_yyerror!(context::_GeoNumericContext,message::AbstractString)
    context.suppress_diagnostics && return nothing
    context.stderr_diagnostics &&
        (print(stderr,"Error   : ",message,"\n");flush(stderr))
    push!(context.exec_errors,String(message))
    return nothing
end

# `yymsg(1, ...)` — a warning; printed like `Msg::Warning`, execution
# continues.
function _geo_yywarn!(context::_GeoNumericContext,message::AbstractString)
    context.suppress_diagnostics && return nothing
    context.stderr_diagnostics &&
        (print(stderr,"Warning : ",message,"\n");flush(stderr))
    push!(context.exec_warnings,String(message))
    return nothing
end

# `Msg::Info` — `Info    : <msg>` on stdout at `General.Verbosity ≥ 4` (the
# default is 5). Recorded nowhere upstream — purely observable output.
function _geo_yyinfo!(context::_GeoNumericContext,message::AbstractString)
    something(_geo_option_number(context,"General",0,"Verbosity"),5.0)<4 &&
        return nothing
    print(stdout,"Info    : ",message,"\n");flush(stdout)
    return nothing
end

# `Msg::Error` outside the parser diagnostics channel — option-table lookups,
# plugin and field actions. Prints immediately and marks the run failed, but
# does not bump `gmsh_yyerrorstate` (no >20 abort).
function _geo_msg_error!(context::_GeoNumericContext,message::AbstractString)
    context.suppress_diagnostics && return nothing
    context.stderr_diagnostics &&
        (print(stderr,"Error   : ",message,"\n");flush(stderr))
    context.msg_error_count+=1
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
            # floating-point fmod (Gmsh.y, FExpr).  Reject the C++ undefined
            # cases instead of depending on a platform-specific conversion or
            # integer trap.
            ta=trunc(a);tb=trunc(b)
            int_min=Float64(typemin(Int32));int_max=Float64(typemax(Int32))
            (int_min<=ta<=int_max && int_min<=tb<=int_max) ||
                _geo_expr_error(parser,
                    "modulo operands are outside Gmsh's signed 32-bit integer range",pos)
            ia=Int32(ta);ib=Int32(tb)
            iszero(ib) && _geo_expr_error(parser,
                "modulo divisor truncates to zero",pos)
            ia==typemin(Int32) && ib==Int32(-1) && _geo_expr_error(parser,
                "modulo is outside its signed 32-bit integer domain",pos)
            rem(ia,ib)
        else
            kind==:plus ? a+b : kind==:minus ? a-b : kind==:star ? a*b :
            kind==:slash ? a/b : _gm_pow(a,b)
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
    elseif name=="Cos"; _gm_cos
    elseif name=="Cosh"; _gm_cosh
    elseif name=="Exp"; _gm_exp
    elseif name=="Fabs" || name=="Abs"; abs
    elseif name=="Floor"; floor
    elseif name=="Log"; _gm_log
    elseif name=="Log10"; _gm_log10
    elseif name=="Round"; x->round(x,RoundNearestTiesUp)
    elseif name=="Sqrt"; sqrt
    elseif name=="Sin"; _gm_sin
    elseif name=="Sinh"; _gm_sinh
    elseif name=="Step"; x->x<0 ? 0.0 : 1.0
    elseif name=="Tan"; _gm_tan
    elseif name=="Tanh"; _gm_tanh
    else; nothing
    end
    binary=if name=="Atan2"; _gm_atan2
    elseif name=="Fmod" || name=="Modulo"; rem
    # Gmsh 4.15.2 spells this as sqrt(a*a + b*b); preserve its overflow and
    # underflow behavior (execution keeps the IEEE754 Inf/NaN like upstream).
    elseif name=="Hypot"; (a,b)->sqrt(a*a+b*b)
    elseif name=="Max"; max
    elseif name=="Min"; min
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

# Process CPU seconds — Gmsh's `Cpu()` via `getrusage`. The rusage head is two
# `timeval`s (user, sys): `{long sec; suseconds_t usec}` where suseconds_t is
# 32-bit on Darwin and 64-bit on Linux. Wall time is the fallback for systems
# without getrusage (finite, nondecreasing — the only portable contract).
function _geo_cpu_time()
    buf=Vector{UInt8}(undef,256)
    ret=ccall(:getrusage,Cint,(Cint,Ptr{UInt8}),0,buf)
    ret==0 || return time()
    @static if Sys.isbsd()
        secs=reinterpret(Int64,buf)
        usecs=reinterpret(Int32,buf)
        return Float64(secs[1])+Float64(usecs[3])/1e6+
               Float64(secs[3])+Float64(usecs[7])/1e6
    else
        vals=reinterpret(Int64,buf)
        return Float64(vals[1])+Float64(vals[2])/1e6+
               Float64(vals[3])+Float64(vals[4])/1e6
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

# `(int)value` — the C double→int cast is undefined outside the Int32 range
# and on non-finite input; fail explicitly there.
function _geo_int32_operand(parser::_GeoExprParser,value::Float64,
                            operation::AbstractString,pos::Int)
    t=trunc(value)
    (isfinite(value) && Float64(typemin(Int32))<=t<=Float64(typemax(Int32))) ||
        _geo_expr_error(parser,"$operation operand is outside Gmsh's signed " *
                               "32-bit integer range",pos)
    return Int32(t)
end

# `((int)a >> (int)b)` / `<<` — C int shifts; the count outside [0,31] is UB
# upstream, so it is rejected rather than emulated.
function _geo_apply_shift(parser::_GeoExprParser,kind::Symbol,a::Float64,
                          b::Float64,pos::Int)
    # The params scan keeps its documented arithmetic subset.
    parser.context.soft_unknown_reads || _geo_expr_error(parser,
        "operator $(kind==:shift_left ? "<<" : ">>") is outside the " *
        "supported arithmetic subset",pos)
    ia=_geo_int32_operand(parser,a,"shift",pos)
    ib=_geo_int32_operand(parser,b,"shift",pos)
    0<=ib<32 || _geo_expr_error(
        parser,"shift count is outside the 32-bit shift range",pos)
    return Float64(kind==:shift_left ? ia<<ib : ia>>ib)
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
     (ncodeunits(in)>3 && in[2]==':' && in[3] in ('/','\\'))) && return String(in)
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
        rhs_values=_geo_multi_term_values(rhs,context,caller,limit,depth+1)
        rhs_values===nothing && return nothing
        lhs=try
            _geo_eval_numeric(
                String(strip(source[firstindex(source):prevind(source,star)])),
                context,"$caller multiplier")
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
            [10.0^(a+(b-a)*Float64(i)/denominator) for i in 0:n-1]
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
_geo_unique_sorted(values::Vector{Float64})=unique(sort(values))

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
        return context.lists[name]
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
        if tag===:int_floor
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

# `Mesh.MeshSizeFactor` and `Mesh.CharacteristicLengthFactor` are two names
# for the single `CTX.mesh.lcFactor` (both DefaultOptions entries dispatch to
# `opt_mesh_lc_factor`), so a write to either name mirrors into the other key.
# The callback also ignores non-positive writes (`if(val > 0)`), which is
# replicated by skipping the store entirely.
function _geo_store_option_number!(context::_GeoNumericContext,
                                   family::String,index::Int,
                                   member::String,value::Float64,
                                   caller::AbstractString)
    v=_geo_option_store_value(family,member,value,caller)
    if family=="Mesh" && index==0 &&
       member in ("MeshSizeFactor","CharacteristicLengthFactor")
        v>0 || return nothing
        context.option_numbers[(family,index,member=="MeshSizeFactor" ?
            "CharacteristicLengthFactor" : "MeshSizeFactor")]=v
    end
    context.option_numbers[(family,index,member)]=v
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

# `NumberOption(GMSH_GET)` — `nothing` for unknown option names.
function _geo_option_number(context::_GeoNumericContext,family::String,
                            index::Int,member::String;warn::Bool=false)
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
        r"^Physical\s+(Point|Curve|Line|Surface|Volume)\s*\(\s*\"([^\"]*)\"\s*(?:,\s*(.*?))?\s*\)\s*(=|\+=|-=|\*=|/=)\s*(.*)$",
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
        r"^Physical\s+(Point|Curve|Line|Surface|Volume)\s*\(\s*(.*?)\s*\)\s*(=|\+=|-=|\*=|/=)\s*(.*)$",
        source)
    unnamed===nothing && return nothing
    tag_source=String(strip(unnamed.captures[2]))
    isempty(tag_source) && return nothing
    return (kind=String(unnamed.captures[1]),name=nothing,
            tag_source=tag_source,op=String(unnamed.captures[3]),
            membership=String(strip(unnamed.captures[4])))
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
    name in _GEO_SIDE_EFFECT_SYMBOLS && throw(ArgumentError(
        "read_geo_params: dynamic tag allocator $name is read-only"))
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
    name in _GEO_SIDE_EFFECT_SYMBOLS && throw(ArgumentError(
        "read_geo_params: dynamic tag allocator $name is read-only"))
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
    for raw in eachline(path)
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
            dim=_PHYS_DIM[physical.kind]
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
                word=="Embedded" || throw(ArgumentError(
                    "read_geo_params: unknown command 'Delete $word { ... }'"))
                return  # `Delete Embedded` clears embeddings — no params state
            end
            if (vm=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+([A-Za-z_][A-Za-z0-9_]*)\s*;?$",
                         tail))!==nothing
                vm.captures[1]=="Empty" && vm.captures[2]=="Views" && return
                throw(ArgumentError("read_geo_params: unknown command " *
                    "'Delete $(vm.captures[1]) $(vm.captures[2])'"))
            end
            if (fm=match(r"^([A-Za-z_][A-Za-z0-9_]*|\"[^\"]*\")\s*\[\s*(.*?)\s*\]\s*;?$",
                         tail))!==nothing
                base=String(strip(fm.captures[1],'"'))
                base=="Field" || throw(ArgumentError(
                    "read_geo_params: unknown command 'Delete $base[...]'"))
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

end # module IO
