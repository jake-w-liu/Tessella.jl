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
evaluator: global mesh sizing, physical groups, and the raw background-field graph.
"""
module IO

using ..MeshTypes: Mesh, nnodes, nsegs, ntris, ntets, node, validate
using Printf: @printf, @sprintf

export read_msh, write_msh, MshFile
export read_stl, read_geo_params, GeoParams, GeoFieldSpec

# gmsh element type codes we handle. (type => n_nodes)
const MSH_POINT = 15
const MSH_LINE  = 1
const MSH_TRI   = 2
const MSH_TET   = 4
const MSH_TET2  = 11    # 10-node (quadratic) tet — read as its 4 corner vertices
const _NN = Dict(MSH_POINT => 1, MSH_LINE => 2, MSH_TRI => 3, MSH_TET => 4, MSH_TET2 => 10)
const _EDIM = Dict(MSH_POINT => 0, MSH_LINE => 1, MSH_TRI => 2, MSH_TET => 3, MSH_TET2 => 3)

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
    _Accum() = new(Float64[], Float64[], Float64[], Dict{Int,Int}(),
                   NTuple{2,Int32}[], NTuple{3,Int32}[], NTuple{4,Int32}[],
                   Int32[], Int32[], Int32[], Dict{Tuple{Int,Int},String}(),
                   Dict{Tuple{Int,Int},Int}(),Set{Int}())
end

@inline function _nodeidx!(acc::_Accum, tag::Int)
    get(acc.tag2idx, tag, 0)
end

"""
    read_msh(path) -> MshFile

Read a gmsh ASCII `.msh` file (format 2.2 or 4.1). Node/element tags are
relabelled to a compact `1:N`; per-element physical tags are preserved.
"""
function read_msh(path::AbstractString)
    open(path, "r") do io
        try
            return _read_msh(io)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError("IO: truncated .msh file"))
            rethrow()
        end
    end
end

function _read_msh(io::Base.IO)
    acc = _Accum()
    version = 0.0
    seen_format=false;seen_names=false;seen_nodes=false;seen_elements=false;seen_entities=false
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
            seen_names && throw(ArgumentError("IO: duplicate \$PhysicalNames section"))
            seen_names=true
            _read_physical_names!(acc, io)
        elseif line == "\$Nodes"
            seen_nodes && throw(ArgumentError("IO: duplicate \$Nodes section"))
            seen_elements && throw(ArgumentError("IO: \$Nodes must precede \$Elements"))
            seen_nodes=true
            version != 0 || throw(ArgumentError("IO: \$Nodes appeared before \$MeshFormat"))
            version==4.1 && !seen_entities &&
                throw(ArgumentError("IO: v4 \$Entities must precede \$Nodes"))
            version < 3 ? _read_nodes_v2!(acc, io) : _read_nodes_v4!(acc, io)
        elseif line == "\$Elements"
            seen_elements && throw(ArgumentError("IO: duplicate \$Elements section"))
            seen_nodes || throw(ArgumentError("IO: \$Elements appeared before \$Nodes"))
            seen_elements=true
            version != 0 || throw(ArgumentError("IO: \$Elements appeared before \$MeshFormat"))
            version < 3 ? _read_elements_v2!(acc, io) : _read_elements_v4!(acc, io)
        elseif line == "\$Entities"
            seen_format || throw(ArgumentError("IO: \$Entities appeared before \$MeshFormat"))
            seen_entities && throw(ArgumentError("IO: duplicate \$Entities section"))
            (seen_nodes || seen_elements) &&
                throw(ArgumentError("IO: \$Entities must precede \$Nodes and \$Elements"))
            seen_entities=true
            version == 4.1 || throw(ArgumentError("IO: \$Entities is only valid in MSH v4.1"))
            _read_entities_v4!(acc, io)
        elseif startswith(line, "\$") && !startswith(line, "\$End")
            _skip_section(io, "\$End" * line[2:end])   # ignore unknown sections
        end
    end
    version != 0 || throw(ArgumentError("IO: missing \$MeshFormat section"))
    seen_nodes || throw(ArgumentError("IO: missing \$Nodes section"))
    return MshFile(_to_mesh(acc), acc.pnames)
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

function _read_physical_names!(acc, io)
    n = parse(Int, strip(readline(io)))
    n >= 0 || throw(ArgumentError("IO: negative physical-name count $n"))
    n<=typemax(Int32) || throw(ArgumentError("IO: physical-name count exceeds Int32"))
    for _ in 1:n
        line = strip(readline(io))
        # format: `dim tag "name"`. Take the name between the first and last quote VERBATIM
        # — splitting on whitespace and rejoining collapses runs of interior spaces/tabs and
        # breaks the name-preservation round-trip contract.
        m = match(r"^([+-]?\d+)\s+([+-]?\d+)\s+\"((?:\\.|[^\"])*)\"\s*$", line)
        m === nothing && throw(ArgumentError("IO: malformed physical-name record '$line'"))
        dim = parse(Int,m.captures[1]); tag = parse(Int,m.captures[2])
        0 <= dim <= 3 || throw(ArgumentError("IO: physical-name dimension $dim is outside 0:3"))
        tag > 0 || throw(ArgumentError("IO: physical-name tag must be positive (got $tag)"))
        _int32_tag(tag,"physical-name")
        haskey(acc.pnames,(dim,tag)) &&
            throw(ArgumentError("IO: duplicate physical name for dimension/tag ($dim,$tag)"))
        acc.pnames[(dim, tag)] = _unescape_name(m.captures[3])
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
function _read_nodes_v2!(acc, io)
    n = parse(Int, strip(readline(io)))
    n >= 0 || throw(ArgumentError("IO: negative v2 node count $n"))
    n <= typemax(Int32)-length(acc.node_x) ||
        throw(ArgumentError("IO: v2 node count exceeds Int32 indexing"))
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

function _read_elements_v2!(acc, io)
    n = parse(Int, strip(readline(io)))
    n >= 0 || throw(ArgumentError("IO: negative v2 element count $n"))
    n<=typemax(Int32) || throw(ArgumentError("IO: v2 element count exceeds Int32"))
    for _ in 1:n
        p = split(strip(readline(io)))
        length(p) >= 3 || throw(ArgumentError("IO: malformed v2 element record"))
        _record_element_tag!(acc,parse(Int,p[1]))
        etype = parse(Int, p[2])
        ntags = parse(Int, p[3])
        ntags >= 0 || throw(ArgumentError("IO: negative element tag count $ntags"))
        length(p) >= 3+ntags || throw(ArgumentError("IO: truncated v2 element tags"))
        phys = ntags >= 1 ? parse(Int, p[4]) : 0
        nodestart = 3 + ntags + 1
        needed = get(_NN,etype,0)
        needed == 0 || length(p) == nodestart+needed-1 ||
            throw(ArgumentError("IO: element type $etype expects $needed node tags"))
        _push_element!(acc, etype, phys, @view p[nodestart:end])
    end
    _expect_end(io, "\$EndElements")
end

# ── v4.1 ────────────────────────────────────────────────────────────────────────
function _read_entities_v4!(acc, io)
    hdr = split(strip(readline(io)))
    length(hdr)==4 || throw(ArgumentError("IO: malformed v4 entity header"))
    nP, nC, nS, nV = parse.(Int, hdr[1:4])
    all(>=(0),(nP,nC,nS,nV)) || throw(ArgumentError("IO: negative v4 entity count"))
    total=try Base.checked_add(Base.checked_add(nP,nC),Base.checked_add(nS,nV)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("IO: v4 entity count overflows Int"))
    end
    total<=typemax(Int32) || throw(ArgumentError("IO: v4 entity count exceeds Int32"))
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

function _read_nodes_v4!(acc, io)
    hdr = split(strip(readline(io)))
    length(hdr) == 4 || throw(ArgumentError("IO: malformed v4 node header"))
    numBlocks = parse(Int, hdr[1]); numNodes = parse(Int, hdr[2])
    declared_min=parse(Int,hdr[3]);declared_max=parse(Int,hdr[4])
    (numBlocks >= 0 && numNodes >= 0) || throw(ArgumentError("IO: negative v4 node count"))
    numNodes<=typemax(Int32) || throw(ArgumentError("IO: $numNodes nodes exceed Int32 indexing"))
    numBlocks<=numNodes || throw(ArgumentError("IO: v4 node-block count exceeds node count"))
    (numNodes==0)==(numBlocks==0) ||
        throw(ArgumentError("IO: v4 node and node-block counts must both be zero or both be positive"))
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
        haskey(acc.ep,(edim,etag)) ||
            throw(ArgumentError("IO: v4 node block references undeclared entity ($edim,$etag)"))
        parametric in (0,1) || throw(ArgumentError("IO: v4 node-block parametric flag must be 0 or 1"))
        nInBlock = parse(Int, bh[4])
        nInBlock > 0 || throw(ArgumentError("IO: v4 node blocks must be nonempty"))
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

function _read_elements_v4!(acc, io)
    ep = acc.ep
    hdr = split(strip(readline(io)))
    length(hdr) == 4 || throw(ArgumentError("IO: malformed v4 element header"))
    numBlocks = parse(Int, hdr[1])
    numElements = parse(Int,hdr[2])
    declared_min=parse(Int,hdr[3]);declared_max=parse(Int,hdr[4])
    (numBlocks >= 0 && numElements >= 0) || throw(ArgumentError("IO: negative v4 element count"))
    numElements<=typemax(Int32) || throw(ArgumentError("IO: v4 element count exceeds Int32"))
    numBlocks<=numElements || throw(ArgumentError("IO: v4 element-block count exceeds element count"))
    (numElements==0)==(numBlocks==0) || throw(ArgumentError(
        "IO: v4 element and element-block counts must both be zero or both be positive"))
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
        haskey(ep,(edim,etag)) ||
            throw(ArgumentError("IO: v4 element block references undeclared entity ($edim,$etag)"))
        haskey(_EDIM,etype) && _EDIM[etype]!=edim &&
            throw(ArgumentError("IO: element type $etype is incompatible with entity dimension $edim"))
        nInBlock > 0 || throw(ArgumentError("IO: v4 element blocks must be nonempty"))
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
    if etype == MSH_LINE
        a = idx(nodetoks[1]); b = idx(nodetoks[2])
        push!(acc.segs, (Int32(a), Int32(b))); push!(acc.seg_tag, ptag)
    elseif etype == MSH_TRI
        a = idx(nodetoks[1]); b = idx(nodetoks[2]); c = idx(nodetoks[3])
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
    write_msh(path, mesh; version=2.2, physical_names=Dict())

Write `mesh` to a gmsh ASCII `.msh` file. `version` is `2.2` or `4.1`.
`physical_names` maps `(dim, tag) => name`. Per-cell physical tags come from the
mesh's `*_tag` vectors; cells are grouped into entity blocks by `(dim, tag)`.
"""
function write_msh(path::AbstractString, mesh::Mesh; version::Real=2.2,
                   physical_names::AbstractDict=Dict{Tuple{Int,Int},String}())
    ver = try
        Float64(version)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("write_msh: version must be representable as Float64"))
    end
    ver in (2.2,4.1) ||
        throw(ArgumentError("write_msh: version must be 2.2 or 4.1 (got $version)"))
    _validate_write_mesh(mesh)
    target = abspath(path); parent = dirname(target)
    isdir(parent) || throw(ArgumentError("write_msh: parent directory does not exist: $parent"))
    # Build beside the destination and atomically rename only after a complete,
    # flushed write. A formatter/error cannot truncate a previously valid mesh.
    mktemp(parent) do tmp, io
        if ver == 2.2
            _write_msh_v2(io, mesh, physical_names)
        else
            _write_msh_v4(io, mesh, physical_names)
        end
        flush(io); close(io)
        mv(tmp,target;force=true)
    end
    return path
end

function _validate_write_mesh(m::Mesh)
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
                                           "\n"=>"\\n", "\r"=>"\\r")

function _write_physical_names(io, physical_names)
    isempty(physical_names) && return
    println(io, "\$PhysicalNames")
    println(io, length(physical_names))
    for ((dim, tag), name) in sort(collect(physical_names); by=x->x[1])
        (dim isa Integer && 0 <= dim <= 3 && tag isa Integer) ||
            throw(ArgumentError("write_msh: invalid physical-name key ($dim,$tag)"))
        tag > 0 || throw(ArgumentError("write_msh: physical-name tag must be positive (got $tag)"))
        _int32_tag(tag,"physical-name")
        name isa AbstractString ||
            throw(ArgumentError("write_msh: physical name for ($dim,$tag) must be a string"))
        println(io, dim, " ", tag, " \"", _escape_name(name), "\"")
    end
    println(io, "\$EndPhysicalNames")
end

function _write_msh_v2(io, m::Mesh, physical_names)
    println(io, "\$MeshFormat"); println(io, "2.2 0 8"); println(io, "\$EndMeshFormat")
    _write_physical_names(io, physical_names)
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

function _write_msh_v4(io, m::Mesh, physical_names)
    println(io, "\$MeshFormat"); println(io, "4.1 0 8"); println(io, "\$EndMeshFormat")
    _write_physical_names(io, physical_names)

    # Group cells into (dim, tag) entity blocks so physical tags survive.
    segblocks = _group_cells(m.segs, m.seg_tag, 2)
    triblocks = _group_cells(m.tris, m.tri_tag, 3)
    tetblocks = _group_cells(m.tets, m.tet_tag, 4)

    # Entity tags are positive geometric identifiers, distinct from physical
    # tags (where zero means "unclassified").  Number them independently per
    # dimension and reserve one volume entity to classify the node block.
    segblocks = [(phys,i,cols) for (i,(phys,cols)) in enumerate(segblocks)]
    triblocks = [(phys,i,cols) for (i,(phys,cols)) in enumerate(triblocks)]
    tetblocks = [(phys,i,cols) for (i,(phys,cols)) in enumerate(tetblocks)]
    node_entity = length(tetblocks)+1

    # Entities: one entity per (dim, physical tag) block, declaring its physical group.
    _write_entities_v4(io, m, segblocks, triblocks, tetblocks, node_entity)

    # Nodes: a single block on the synthetic volume entity when nonempty.
    println(io, "\$Nodes")
    println(io, nnodes(m)==0 ? 0 : 1, " ", nnodes(m), " ", nnodes(m) == 0 ? 0 : 1, " ", nnodes(m))
    if nnodes(m)>0
        println(io, 3, " ", node_entity, " ", 0, " ", nnodes(m))
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

# group cell columns by tag → Dict(tag => Vector{col index})
function _group_cells(cells::Matrix{Int32}, tags::Vector{Int32}, k::Int)
    g = Dict{Int32,Vector{Int}}()
    @inbounds for j in axes(cells, 2)
        push!(get!(g, tags[j], Int[]), j)
    end
    return sort(collect(g); by=x->x[1])
end

function _write_entities_v4(io, m, segblocks, triblocks, tetblocks, node_entity)
    println(io, "\$Entities")
    println(io, 0, " ", length(segblocks), " ", length(triblocks), " ", length(tetblocks)+1)
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
    read_stl(path; merge_tol=1e-9) -> Mesh

Read an ASCII or binary STL. Coincident vertices within `merge_tol` (relative to
the bounding-box diagonal) are merged so the result is a connected triangle
surface, not a triangle soup.
"""
function read_stl(path::AbstractString; merge_tol::Real=1e-9)
    mtol=try
        Float64(merge_tol)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_stl: merge_tol must be representable as Float64"))
    end
    (isfinite(mtol) && mtol >= 0) ||
        throw(ArgumentError("read_stl: merge_tol must be finite and non-negative (got $merge_tol)"))
    isbin = _stl_is_binary(path)
    tris_xyz = isbin ? _read_stl_binary(path) : _read_stl_ascii(path)
    out=_weld_triangles(tris_xyz, mtol)
    d=validate(out)
    d.ok || throw(ArgumentError("read_stl: welded surface is invalid — "*join(d.messages,"; ")))
    return out
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

function _read_stl_ascii(path)
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

function _read_stl_binary(path)
    open(path, "r") do io
        filesize(path) >= 84 || throw(ArgumentError("read_stl: truncated binary STL header"))
        skip(io, 80)
        ntri = Int(ltoh(read(io, UInt32)))
        ntri>0 || throw(ArgumentError("read_stl: binary STL contains no facets"))
        ntri<=typemax(Int32) || throw(ArgumentError("read_stl: binary facet count exceeds Int32"))
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

function _weld_triangles(tris_xyz::Vector{NTuple{9,Float64}}, reltol::Real)
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
    diag = hypot(hi[1]-lo[1],hi[2]-lo[2],hi[3]-lo[3])
    isfinite(diag) || throw(ArgumentError("read_stl: bounding-box diagonal is non-finite"))
    tol = diag * reltol
    isfinite(tol) || throw(ArgumentError("read_stl: absolute welding tolerance is non-finite"))
    # Quantize relative to the bbox min corner `lo`, not the absolute coordinate:
    # a constant per-axis bucket shift (does not change which vertices merge) that
    # bounds every key to [0, diag/tol] so far-from-origin coords can't overflow
    # Int64 in round(Int, ·) — an unshifted absolute coord ~1e11 does.
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    exactmap=Dict{NTuple{3,Float64},Int32}()
    buckets=Dict{NTuple{3,Int},Vector{Int32}}()
    inv = tol == 0 ? 0.0 : 1/tol
    if tol > 0
        cells = diag*inv
        (isfinite(cells) && cells <= typemax(Int)-2) ||
            throw(ArgumentError("read_stl: merge_tol is too small for Int bucket indices at this extent"))
    end
    keyof(p) = (floor(Int,(p[1]-lo[1])*inv),floor(Int,(p[2]-lo[2])*inv),floor(Int,(p[3]-lo[3])*inv))
    function getid(p)
        if tol == 0
            key=(p[1]==0 ? 0.0 : p[1],p[2]==0 ? 0.0 : p[2],p[3]==0 ? 0.0 : p[3])
            return get!(exactmap,key) do
                length(xs)<typemax(Int32) || throw(ArgumentError("read_stl: unique vertex count exceeds Int32"))
                push!(xs,p[1]);push!(ys,p[2]);push!(zs,p[3]);Int32(length(xs))
            end
        end
        key=keyof(p)
        for dz in -1:1, dy in -1:1, dx in -1:1
            ids=get(buckets,(key[1]+dx,key[2]+dy,key[3]+dz),nothing)
            ids === nothing && continue
            for id in ids
                hypot(p[1]-xs[id],p[2]-ys[id],p[3]-zs[id]) <= tol && return id
            end
        end
        length(xs)<typemax(Int32) || throw(ArgumentError("read_stl: unique vertex count exceeds Int32"))
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
# parser-populated five-argument form below preserves exact statement order for
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
name` map, dim ∈ {0:point,1:curve,2:surface,3:volume}), raw `fields`, and the
`background_field` tag. Missing numeric mesh options are represented by `NaN`;
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
end

# Preserve the former full positional constructor, adding the new global with
# its Gmsh default.
GeoParams(mesh_size_min::Real,mesh_size_max::Real,mesh_size_factor::Real,
          random_seed::Integer,physical_groups::Dict{Tuple{Int,Int},String},
          fields::Dict{Int,GeoFieldSpec},background_field::Integer,
          boundary_layer_fields::Vector{Int},geometry_tolerance::Real)=
    GeoParams(Float64(mesh_size_min),Float64(mesh_size_max),Float64(mesh_size_factor),
              Int(random_seed),physical_groups,fields,Int(background_field),
              boundary_layer_fields,Float64(geometry_tolerance),5)

# Preserve the eight-argument constructor introduced with boundary-layer
# declarations.
GeoParams(mesh_size_min::Real,mesh_size_max::Real,mesh_size_factor::Real,
          random_seed::Integer,physical_groups::Dict{Tuple{Int,Int},String},
          fields::Dict{Int,GeoFieldSpec},background_field::Integer,
          boundary_layer_fields::Vector{Int})=
    GeoParams(Float64(mesh_size_min),Float64(mesh_size_max),Float64(mesh_size_factor),
              Int(random_seed),physical_groups,fields,Int(background_field),
              boundary_layer_fields,NaN,5)

# Preserve the seven-argument field-aware constructor used before boundary-layer
# declarations became part of the parsed model.
GeoParams(mesh_size_min::Real,mesh_size_max::Real,mesh_size_factor::Real,
          random_seed::Integer,physical_groups::Dict{Tuple{Int,Int},String},
          fields::Dict{Int,GeoFieldSpec},background_field::Integer)=
    GeoParams(Float64(mesh_size_min),Float64(mesh_size_max),Float64(mesh_size_factor),
              Int(random_seed),physical_groups,fields,Int(background_field),Int[],NaN,5)

# Preserve the original public positional constructor.
GeoParams(mesh_size_min::Real, mesh_size_max::Real, random_seed::Integer,
          physical_groups::Dict{Tuple{Int,Int},String}) =
    GeoParams(Float64(mesh_size_min), Float64(mesh_size_max), 1.0, Int(random_seed),
              physical_groups, Dict{Int,GeoFieldSpec}(), 0, Int[],NaN,5)

const _PHYS_DIM = Dict("Point"=>0, "Curve"=>1, "Line"=>1, "Surface"=>2, "Volume"=>3)

const _MAX_GEO_STATEMENT_BYTES=1_000_000
const _MAX_GEO_EXPRESSION_BYTES=65_536
const _MAX_GEO_EXPRESSION_TOKENS=4_096
const _MAX_GEO_EXPRESSION_DEPTH=128
const _MAX_GEO_LIST_ITEMS=65_536

# This is deliberately a constant-expression evaluator, not a Julia evaluator
# and not a general Gmsh interpreter.  Keeping a small lexer/parser here makes
# unknown identifiers, assignments and stateful built-ins impossible to execute.
struct _GeoExprToken
    kind::Symbol
    text::String
    value::Float64
    pos::Int
end

mutable struct _GeoNumericContext
    values::Dict{String,Float64}
    unavailable::Dict{String,String}
end
_GeoNumericContext()=_GeoNumericContext(Dict{String,Float64}(),Dict{String,String}())

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
const _GEO_NONCONSTANT_FUNCTIONS=Set((
    "Rand","DefineNumber","GetNumber","GetValue","Exists","FileExists",
    "StringToName","S2N","Find","StrFind","StrCmp","StrLen","TextAttributes"))

@inline _geo_ascii_letter(c::Char)=('a'<=c<='z') || ('A'<=c<='Z') || c=='_'
@inline _geo_ascii_digit(c::Char)=('0'<=c<='9')
@inline _geo_ascii_ident(c::Char)=_geo_ascii_letter(c) || _geo_ascii_digit(c)

function _geo_expr_preview(source::AbstractString)
    ncodeunits(source)<=160 && return String(source)
    return String(first(source,120))*"…"*String(last(source,24))
end

function _geo_expr_error(parser::_GeoExprParser,message::AbstractString,
                         pos::Int=parser.token.pos)
    throw(ArgumentError("$(parser.caller): $message at byte $pos in expression `" *
                        _geo_expr_preview(parser.source)*"`"))
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
        value=tryparse(Float64,text)
        (value!==nothing && isfinite(value)) ||
            _geo_expr_error(parser,"numeric literal must be finite",start)
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
        elseif c in ('=','!','<','>','&','|','?',':')
            :unsupported_operator
        elseif c in ('"','\'')
            :quoted
        else
            :invalid
        end
        if kind==:side_effect
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
    isfinite(result) || _geo_expr_error(parser,"$operation produced a non-finite value",pos)
    return result
end

function _geo_apply_binary(parser::_GeoExprParser,kind::Symbol,a::Float64,
                           b::Float64,pos::Int)
    label=kind==:plus ? "addition" : kind==:minus ? "subtraction" :
          kind==:star ? "multiplication" : kind==:slash ? "division" :
          kind==:percent ? "modulo" : "exponentiation"
    value=try
        kind==:plus ? a+b : kind==:minus ? a-b : kind==:star ? a*b :
        kind==:slash ? a/b : kind==:percent ? rem(a,b) : a^b
    catch err
        err isa InterruptException && rethrow()
        (err isa DomainError || err isa OverflowError || err isa DivideError) || rethrow()
        _geo_expr_error(parser,"$label is outside its finite real domain",pos)
    end
    return _geo_finite_result(parser,value,label,pos)
end

function _geo_apply_function(parser::_GeoExprParser,name::String,
                             args::Vector{Float64},pos::Int)
    name in _GEO_NONCONSTANT_FUNCTIONS && _geo_expr_error(parser,
        "non-constant or externally stateful function $name is not supported",pos)
    unary=if name=="Acos"; acos
    elseif name=="Asin"; asin
    elseif name=="Atan"; atan
    elseif name=="Ceil"; ceil
    elseif name=="Cos"; cos
    elseif name=="Cosh"; cosh
    elseif name=="Exp"; exp
    elseif name=="Fabs" || name=="Abs"; abs
    elseif name=="Floor"; floor
    elseif name=="Log"; log
    elseif name=="Log10"; log10
    elseif name=="Round"; x->round(x,RoundNearestTiesUp)
    elseif name=="Sqrt"; sqrt
    elseif name=="Sin"; sin
    elseif name=="Sinh"; sinh
    elseif name=="Step"; x->x<0 ? 0.0 : 1.0
    elseif name=="Tan"; tan
    elseif name=="Tanh"; tanh
    else; nothing
    end
    binary=if name=="Atan2"; (y,x)->atan(y,x)
    elseif name=="Fmod" || name=="Modulo"; rem
    elseif name=="Hypot"; hypot
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
            _geo_expr_error(parser,"function $name is outside its finite real domain",pos)
        end
        return _geo_finite_result(parser,value,"function $name",pos)
    end
    _geo_expr_error(parser,"unknown numeric function $name",pos)
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
    value=_geo_parse_unary!(parser)
    while parser.token.kind in (:star,:slash,:percent)
        kind=parser.token.kind;pos=parser.token.pos;_geo_advance!(parser)
        value=_geo_apply_binary(parser,kind,value,_geo_parse_unary!(parser),pos)
    end
    return value
end

function _geo_parse_unary!(parser::_GeoExprParser)
    if parser.token.kind==:plus || parser.token.kind==:minus
        kind=parser.token.kind;pos=parser.token.pos;_geo_advance!(parser)
        _geo_enter!(parser)
        value=_geo_parse_unary!(parser)
        _geo_leave!(parser)
        kind==:minus && (value=_geo_finite_result(parser,-value,"unary minus",pos))
        return value
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
        value=_geo_parse_additive!(parser)
        parser.token.kind==:right_paren ||
            _geo_expr_error(parser,"expected closing parenthesis")
        _geo_advance!(parser);_geo_leave!(parser)
        return value
    elseif token.kind==:identifier
        name=token.text;_geo_advance!(parser)
        if parser.token.kind==:left_paren || parser.token.kind==:left_bracket
            opener=parser.token.kind
            closer=opener==:left_paren ? :right_paren : :right_bracket
            _geo_advance!(parser);_geo_enter!(parser)
            args=Float64[]
            if parser.token.kind!=closer
                while true
                    push!(args,_geo_parse_additive!(parser))
                    parser.token.kind==:comma || break
                    _geo_advance!(parser)
                end
            end
            parser.token.kind==closer || _geo_expr_error(parser,
                opener==:left_paren ? "expected closing parenthesis" :
                                      "expected closing bracket")
            _geo_advance!(parser);_geo_leave!(parser)
            return _geo_apply_function(parser,name,args,token.pos)
        elseif name=="Pi"
            return Float64(pi)
        elseif haskey(parser.context.values,name)
            return parser.context.values[name]
        elseif name in _GEO_SIDE_EFFECT_SYMBOLS
            _geo_expr_error(parser,"side-effecting Gmsh symbol $name is not supported",token.pos)
        elseif haskey(parser.context.unavailable,name)
            _geo_expr_error(parser,
                "scalar variable $name is unavailable ($(parser.context.unavailable[name]))",
                token.pos)
        else
            _geo_expr_error(parser,"unknown scalar identifier $name",token.pos)
        end
    elseif token.kind==:side_effect
        _geo_expr_error(parser,"increment and decrement operators are not supported",token.pos)
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

function _geo_eval_numeric(raw::AbstractString,context::_GeoNumericContext,
                           caller::AbstractString)
    source=String(strip(raw))
    isempty(source) && throw(ArgumentError("$caller: numeric expression must not be empty"))
    ncodeunits(source)<=_MAX_GEO_EXPRESSION_BYTES || throw(ArgumentError(
        "$caller: expression exceeds $_MAX_GEO_EXPRESSION_BYTES bytes"))
    parser=_GeoExprParser(source,firstindex(source),_GeoExprToken(:eof,"",0.0,1),
                          0,0,context,String(caller))
    _geo_advance!(parser)
    value=_geo_parse_additive!(parser)
    parser.token.kind==:eof || begin
        token=parser.token
        if token.kind==:side_effect
            _geo_expr_error(parser,"increment and decrement operators are not supported",
                            token.pos)
        elseif token.kind==:unsupported_operator
            _geo_expr_error(parser,
                "operator $(token.text) is outside the supported arithmetic subset",token.pos)
        end
        _geo_expr_error(parser,"unexpected token $(repr(token.text))",token.pos)
    end
    return _geo_finite_result(parser,value,"expression",1)
end

@inline _geo_number_source(value::Float64)=repr(value)

function _geo_int_value(value::Float64,caller::AbstractString)
    return try
        trunc(Int,value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError) || rethrow()
        throw(ArgumentError("$caller: value is outside the platform Int range"))
    end
end

const _GMSH_TAG_MAX=Int(typemax(Int32))

function _geo_positive_gmsh_tag(value::Float64,caller::AbstractString)
    tag=_geo_int_value(value,caller)
    tag>0 || throw(ArgumentError("$caller must evaluate to a positive tag"))
    tag<=_GMSH_TAG_MAX || throw(ArgumentError(
        "$caller exceeds Gmsh's signed 32-bit tag range"))
    return tag
end

function _geo_split_list(raw::AbstractString,caller::AbstractString)
    value=String(strip(raw))
    (startswith(value,"{") && endswith(value,"}")) ||
        throw(ArgumentError("$caller: expected a brace-delimited list"))
    body=String(strip(value[nextind(value,firstindex(value)):prevind(value,lastindex(value))]))
    isempty(body) && return String[]
    items=String[];start=firstindex(body);i=start;last=lastindex(body)
    parens=0;brackets=0
    while i<=last
        c=body[i]
        if c=='(';parens+=1
        elseif c==')'
            parens-=1;parens>=0 || throw(ArgumentError(
                "$caller: unmatched closing parenthesis in list"))
        elseif c=='[';brackets+=1
        elseif c==']'
            brackets-=1;brackets>=0 || throw(ArgumentError(
                "$caller: unmatched closing bracket in list"))
        elseif c=='{' || c=='}'
            throw(ArgumentError("$caller: nested brace lists are not supported"))
        elseif c==',' && parens==0 && brackets==0
            item=String(strip(body[start:prevind(body,i)]))
            isempty(item) && throw(ArgumentError("$caller: list contains an empty entry"))
            length(items)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
                "$caller: list exceeds $_MAX_GEO_LIST_ITEMS entries"))
            push!(items,item);start=nextind(body,i)
        end
        i=nextind(body,i)
    end
    parens==0 || throw(ArgumentError("$caller: unmatched opening parenthesis in list"))
    brackets==0 || throw(ArgumentError("$caller: unmatched opening bracket in list"))
    item=String(strip(body[start:last]))
    isempty(item) && throw(ArgumentError("$caller: list contains an empty entry"))
    length(items)<_MAX_GEO_LIST_ITEMS || throw(ArgumentError(
        "$caller: list exceeds $_MAX_GEO_LIST_ITEMS entries"))
    push!(items,item)
    return items
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

function _geo_normalize_field_option(raw::AbstractString,name_raw::AbstractString,
                                     context::_GeoNumericContext,caller::AbstractString)
    name=String(name_raw);caller_string=String(caller)
    value=String(strip(raw))
    name in _GEO_FIELD_RAW_OPTIONS && return value
    if name in _GEO_FIELD_FLOAT_LIST_OPTIONS || name in _GEO_FIELD_INTEGER_LIST_OPTIONS
        items=_geo_split_list(value,caller_string)
        normalized=String[]
        for item in items
            if name in _GEO_FIELD_ENTITY_LIST_OPTIONS &&
               occursin(r"^[+-]?[A-Za-z_][A-Za-z0-9_]*(?:\[\])?$",item)
                bare=replace(item,r"^[+-]"=>"");bare=endswith(bare,"[]") ? bare[1:end-2] : bare
                if !haskey(context.values,bare) && bare!="Pi"
                    push!(normalized,item)
                    continue
                end
            end
            number=_geo_eval_numeric(item,context,"$caller_string entry")
            if name in _GEO_FIELD_INTEGER_LIST_OPTIONS
                integer=name=="FieldsList" ?
                    _geo_positive_gmsh_tag(number,"$caller_string entry") :
                    _geo_int_value(number,"$caller_string entry")
                push!(normalized,string(integer))
            else
                push!(normalized,_geo_number_source(number))
            end
        end
        return "{"*join(normalized,", ")*"}"
    elseif name in _GEO_FIELD_NUMERIC_OPTIONS
        (startswith(value,"{") || endswith(value,"}")) && throw(ArgumentError(
            "$caller_string: numeric scalar option cannot use a brace list"))
        number=_geo_eval_numeric(value,context,caller_string)
        return name in _GEO_FIELD_INTEGER_OPTIONS ?
            string(_geo_int_value(number,caller_string)) : _geo_number_source(number)
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
    for name in keys(context.values)
        context.unavailable[name]=String(reason)
    end
    empty!(context.values)
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

function _geo_record_scalar!(context::_GeoNumericContext,name_raw::AbstractString,
                             raw::AbstractString)
    name=String(name_raw)
    # `Pi` is a lexical Gmsh constant, not a mutable scalar binding.
    name=="Pi" && return nothing
    try
        context.values[name]=_geo_eval_numeric(raw,context,
            "read_geo_params: scalar variable $name")
        delete!(context.unavailable,name)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        delete!(context.values,name)
        message=sprint(showerror,err)
        context.unavailable[name]=ncodeunits(message)<=240 ? message :
            String(first(message,220))*"…"
    end
    return nothing
end

function _geo_expression_field_tags(raw::AbstractString,context::_GeoNumericContext,
                                    caller::AbstractString)
    value=String(strip(raw))
    items=if startswith(value,"{") || endswith(value,"}")
        (startswith(value,"{") && endswith(value,"}")) ||
            throw(ArgumentError("$caller: malformed field-tag list $raw"))
        _geo_split_list(value,String(caller))
    else
        isempty(value) ? String[] : String[value]
    end
    isempty(items) && throw(ArgumentError("$caller: field-tag list must not be empty"))
    tags=Int[]
    for item in items
        numeric=_geo_eval_numeric(item,context,"$caller entry")
        tag=_geo_positive_gmsh_tag(numeric,"$caller entry")
        tag in tags || push!(tags,tag)
    end
    return tags
end

function _geo_positive_tag_value(raw::AbstractString,context::_GeoNumericContext,
                                 caller::AbstractString)
    numeric=_geo_eval_numeric(raw,context,caller)
    return _geo_positive_gmsh_tag(numeric,caller)
end

function _gmsh_random_seed(value::Float64)
    # Gmsh 4.15.2 stores this option as an unsigned 32-bit value: assignments
    # are truncated toward zero and clamped to the option range.
    value<=0 && return 0
    upper=min(Float64(typemax(Int)),Float64(typemax(UInt32)))
    value>=upper && return _geo_int_value(upper,"read_geo_params: Mesh.RandomSeed")
    return _geo_int_value(value,"read_geo_params: Mesh.RandomSeed")
end

function _scan_geo_statements(consume,path::AbstractString)
    buffer=IOBuffer();quote_char='\0';block_comment=false
    for raw in eachline(path)
        i=firstindex(raw);lastindex_raw=lastindex(raw)
        while i<=lastindex_raw
            c=raw[i];j=nextind(raw,i)
            nextc=j<=lastindex_raw ? raw[j] : '\0'
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
                block_comment=true;i=nextind(raw,j);continue
            elseif c=='"' || c=='\''
                quote_char=c;write(buffer,c)
            elseif c==';'
                write(buffer,c)
                statement=strip(String(take!(buffer)))
                isempty(statement) || consume(statement)
            else
                write(buffer,c)
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

"""
    read_geo_params(path) -> GeoParams

Scan a gmsh `.geo`/`.geo_unrolled` for mesh-sizing options, Physical group
declarations, and `Field[...]`/`Background Field` statements. Deterministic,
finite arithmetic expressions can use `Pi`, prior scalar variables and Gmsh's
pure numeric functions, including in explicit Physical-group tags and
`Field[...]` tags. Loops, macros, option reads, random/external functions, list
ranges, CSG and Boolean geometry are deliberately not evaluated. Numeric field
options are normalized to literals; geometric references and string or
point-dependent field expressions remain source strings.
"""
function read_geo_params(path::AbstractString)
    smin = NaN; smax = NaN; sfactor = 1.0; seed = 0; background = 0
    geometry_tolerance=NaN
    mesh_size_from_curvature=0;boundary_layer_fan_elements=5
    boundary_layers=Int[]
    groups = Dict{Tuple{Int,Int},String}()
    kinds = Dict{Int,String}()
    options = Dict{Int,Dict{String,String}}()
    option_order=Dict{Int,Vector{String}}()
    creation_curvature=Dict{Int,Int}()
    context=_GeoNumericContext()
    control_depth=0
    _scan_geo_statements(path) do line
        endswith(line,";") || throw(ArgumentError(
            "read_geo_params: internal statement scanner lost a semicolon"))
        raw_body=String(strip(line[firstindex(line):prevind(line,lastindex(line))]))
        raw_code=_geo_unquoted_code(raw_body)

        # Any control-flow/macro context could mutate prior scalar bindings.  We
        # do not interpret it, so invalidate those bindings instead of using a
        # stale value later.
        if occursin(r"\b(?:For|EndFor|If|ElseIf|Else|EndIf|Macro|Function|Return|Call|Include|DefineConstant|UndefineConstant)\b",raw_code)
            _geo_invalidate_context!(context,"an unsupported loop, conditional, macro or include may have changed it")
        end
        # Gmsh's EndIf/EndFor/Return do not carry semicolons; the streaming
        # scanner consequently receives them as a harmless prefix of the first
        # statement after the closed block.
        body,closed=_geo_strip_control_terminators(raw_body)
        control_depth=max(0,control_depth-closed)
        code=_geo_unquoted_code(body)
        opened=count(_ -> true,eachmatch(r"\b(?:For|If|Macro|Function)\b",code))
        control_depth+=opened
        if control_depth>0
            _geo_invalidate_context!(context,
                "an unsupported loop, conditional or macro may have changed it")
            _geo_relevant_code(code) && throw(ArgumentError(
                "read_geo_params: malformed relevant statement or unsupported control-flow/macro context: " *
                _geo_expr_preview(body)))
            return
        end

        mesh=match(r"^(Mesh\.(?:MeshSizeMin|MeshSizeMax|MeshSizeFactor|RandomSeed|MeshSizeFromCurvature|MinimumElementsPerTwoPi|BoundaryLayerFanElements|BoundaryLayerFanPoints)|Geometry\.Tolerance)\s*=\s*(.*)$",body)
        if mesh!==nothing
            key=mesh.captures[1];raw=String(strip(mesh.captures[2]))
            value=_geo_eval_numeric(raw,context,"read_geo_params: $key")
            if key=="Mesh.MeshSizeMin";smin=value
            elseif key=="Mesh.MeshSizeMax";smax=value
            elseif key=="Mesh.MeshSizeFactor";sfactor=value
            elseif key=="Mesh.RandomSeed";seed=_gmsh_random_seed(value)
            elseif key=="Geometry.Tolerance";geometry_tolerance=value
            elseif key=="Mesh.MeshSizeFromCurvature" ||
                   key=="Mesh.MinimumElementsPerTwoPi"
                mesh_size_from_curvature=_gmsh_int_option(value,key)
            else
                boundary_layer_fan_elements=_gmsh_int_option(value,key)
            end
            return
        end

        # Physical Volume("air", 1) = {...}; / Physical Surface("s", 3) = {...};
        pm=match(r"^Physical\s+(Point|Curve|Line|Surface|Volume)\s*\(\s*\"([^\"]*)\"\s*,\s*(.+?)\s*\)\s*=\s*\{.*\}$",body)
        if pm!==nothing
            dim = _PHYS_DIM[pm.captures[1]]
            name = pm.captures[2]
            tag = _geo_positive_tag_value(pm.captures[3],context,
                "read_geo_params: Physical $(pm.captures[1]) tag")
            groups[(dim, tag)] = name
            return
        end

        fm=match(r"^Field\s*\[\s*(.*)\s*\]\s*=\s*([A-Za-z][A-Za-z0-9_]*)$",body)
        if fm!==nothing
            tag=_geo_positive_tag_value(fm.captures[1],context,
                "read_geo_params: Field declaration tag")
            kind=fm.captures[2]
            haskey(kinds,tag) &&
                throw(ArgumentError("read_geo_params: duplicate declaration for Field[$tag]"))
            kinds[tag]=kind
            creation_curvature[tag]=mesh_size_from_curvature
            return
        end

        om=match(r"^Field\s*\[\s*(.*)\s*\]\.([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*)$",body)
        if om!==nothing
            tag=_geo_positive_tag_value(om.captures[1],context,
                "read_geo_params: Field option tag")
            name=om.captures[2];value=String(strip(om.captures[3]))
            isempty(value) &&
                throw(ArgumentError("read_geo_params: Field[$tag].$name has an empty value"))
            caller="read_geo_params: Field[$tag].$name"
            normalized=_geo_normalize_field_option(value,name,context,caller)
            get!(() -> Dict{String,String}(),options,tag)[name]=normalized
            push!(get!(() -> String[],option_order,tag),name)
            return
        end

        bm=match(r"^Background\s+Field\s*=\s*(.*)$",body)
        if bm!==nothing
            tags=_geo_expression_field_tags(bm.captures[1],context,
                                            "read_geo_params: Background Field")
            length(tags)==1 || throw(ArgumentError(
                "read_geo_params: Background Field requires exactly one field tag"))
            background=tags[1]
            return
        end

        blm=match(r"^BoundaryLayer\s+Field\s*=\s*(.*)$",body)
        if blm!==nothing
            for tag in _geo_expression_field_tags(blm.captures[1],context,
                                                   "read_geo_params: BoundaryLayer Field")
                tag in boundary_layers || push!(boundary_layers,tag)
            end
            return
        end

        scalar=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$",body)
        if scalar!==nothing
            _geo_record_scalar!(context,scalar.captures[1],scalar.captures[2])
            return
        end
        mutation=match(r"^(?:\+\+|--)?\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?:\+\+|--|\+=|-=|\*=|/=).*$",body)
        if mutation!==nothing
            name=mutation.captures[1]
            delete!(context.values,name)
            context.unavailable[name]="an unsupported increment or compound assignment changed it"
            return
        end
        array_assignment=match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*\[.*\]\s*(?:=|\+=|-=|\*=|/=).*$",body)
        if array_assignment!==nothing
            name=array_assignment.captures[1]
            delete!(context.values,name)
            context.unavailable[name]="an unsupported list assignment changed it"
            return
        end

        # Reject relevant assignments hidden in a loop/macro/prefix or malformed
        # on their left-hand side. Quoted occurrences were removed from `code`.
        if _geo_relevant_code(code)
            throw(ArgumentError(
                "read_geo_params: malformed relevant statement or unsupported control-flow/macro context: " *
                _geo_expr_preview(body)))
        end
    end
    undeclared=sort!(collect(setdiff(keys(options),keys(kinds))))
    isempty(undeclared) || throw(ArgumentError(
        "read_geo_params: options provided for undeclared field tag(s) $(join(undeclared, ", "))"))
    background==0 || haskey(kinds,background) || throw(ArgumentError(
        "read_geo_params: Background Field $background is not declared"))
    missing_boundary=filter(tag->!haskey(kinds,tag),boundary_layers)
    isempty(missing_boundary) || throw(ArgumentError(
        "read_geo_params: BoundaryLayer Field tag(s) $(join(missing_boundary, ", ")) are not declared"))
    fields=Dict{Int,GeoFieldSpec}()
    for (tag,kind) in kinds
        fields[tag]=GeoFieldSpec(tag,kind,get(options,tag,Dict{String,String}()),
                                 get(option_order,tag,String[]),
                                 get(creation_curvature,tag,0))
    end
    return GeoParams(smin,smax,sfactor,seed,groups,fields,background,boundary_layers,
                     geometry_tolerance,boundary_layer_fan_elements)
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
