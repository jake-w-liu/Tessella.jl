"""
    IO

Mesh file I/O (PLAN.md §3 "IO"): gmsh **MSH v2.2 and v4.1** (ASCII) read/write,
ASCII/binary **STL** ingest for boundary surface meshes, and a lightweight
**`.geo`** parameter/structure scanner.

The round-trip contract (DEVELOPMENT.md CRC gate for Stage 0) is *connectivity
preservation*: reading a mesh and writing it back — in either format version —
must reproduce the same topology, verified by `MeshTypes.mesh_crc`. gmsh allows
arbitrary node/element tags; we relabel to a compact `1:N` on read (connectivity
is invariant under a consistent relabel) while preserving physical-group tags.

Full evaluation of a `.geo` OpenCASCADE CSG script (Booleans → BREP → faces) is
the Stage-5 geometry kernel and is *not* attempted here; `read_geo_params` reads
only what is tractable without OCC (mesh sizing, physical-group declarations).
"""
module IO

using ..MeshTypes: Mesh, nnodes, nsegs, ntris, ntets, node
using Printf: @printf, @sprintf

export read_msh, write_msh, MshFile
export read_stl, read_geo_params, GeoParams

# gmsh element type codes we handle. (type => n_nodes)
const MSH_POINT = 15
const MSH_LINE  = 1
const MSH_TRI   = 2
const MSH_TET   = 4
const MSH_TET2  = 11    # 10-node (quadratic) tet — read as its 4 corner vertices
const _NN = Dict(MSH_POINT => 1, MSH_LINE => 2, MSH_TRI => 3, MSH_TET => 4, MSH_TET2 => 10)

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
    _Accum() = new(Float64[], Float64[], Float64[], Dict{Int,Int}(),
                   NTuple{2,Int32}[], NTuple{3,Int32}[], NTuple{4,Int32}[],
                   Int32[], Int32[], Int32[], Dict{Tuple{Int,Int},String}(),
                   Dict{Tuple{Int,Int},Int}())
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
        return _read_msh(io)
    end
end

function _read_msh(io::Base.IO)
    acc = _Accum()
    version = 0.0
    while !eof(io)
        line = strip(readline(io))
        isempty(line) && continue
        if line == "\$MeshFormat"
            fmt = split(strip(readline(io)))
            version = parse(Float64, fmt[1])
            filetype = length(fmt) >= 2 ? parse(Int, fmt[2]) : 0
            filetype == 0 || error("IO: binary .msh not supported (file-type $filetype); use ASCII")
            _expect_end(io, "\$EndMeshFormat")
        elseif line == "\$PhysicalNames"
            _read_physical_names!(acc, io)
        elseif line == "\$Nodes"
            version < 3 ? _read_nodes_v2!(acc, io) : _read_nodes_v4!(acc, io)
        elseif line == "\$Elements"
            version < 3 ? _read_elements_v2!(acc, io) : _read_elements_v4!(acc, io)
        elseif line == "\$Entities"
            _read_entities_v4!(acc, io)
        elseif startswith(line, "\$") && !startswith(line, "\$End")
            _skip_section(io, "\$End" * line[2:end])   # ignore unknown sections
        end
    end
    return MshFile(_to_mesh(acc), acc.pnames)
end

function _expect_end(io, tok)
    l = strip(readline(io))
    l == tok || error("IO: expected $tok, got '$l'")
end

function _skip_section(io, endtok)
    while !eof(io)
        strip(readline(io)) == endtok && return
    end
    error("IO: unterminated section (missing $endtok)")
end

function _read_physical_names!(acc, io)
    n = parse(Int, strip(readline(io)))
    for _ in 1:n
        parts = split(strip(readline(io)))
        dim = parse(Int, parts[1]); tag = parse(Int, parts[2])
        # name is quoted, possibly containing spaces — rejoin remaining tokens
        nm = join(parts[3:end], " ")
        nm = strip(nm, ['"'])
        acc.pnames[(dim, tag)] = nm
    end
    _expect_end(io, "\$EndPhysicalNames")
end

# ── v2 ──────────────────────────────────────────────────────────────────────────
function _read_nodes_v2!(acc, io)
    n = parse(Int, strip(readline(io)))
    sizehint!(acc.node_x, n); sizehint!(acc.node_y, n); sizehint!(acc.node_z, n)
    for _ in 1:n
        p = split(strip(readline(io)))
        tag = parse(Int, p[1])
        push!(acc.node_x, parse(Float64, p[2]))
        push!(acc.node_y, parse(Float64, p[3]))
        push!(acc.node_z, parse(Float64, p[4]))
        acc.tag2idx[tag] = length(acc.node_x)
    end
    _expect_end(io, "\$EndNodes")
end

function _read_elements_v2!(acc, io)
    n = parse(Int, strip(readline(io)))
    for _ in 1:n
        p = split(strip(readline(io)))
        etype = parse(Int, p[2])
        ntags = parse(Int, p[3])
        phys = ntags >= 1 ? parse(Int, p[4]) : 0
        nodestart = 3 + ntags + 1
        _push_element!(acc, etype, phys, @view p[nodestart:end])
    end
    _expect_end(io, "\$EndElements")
end

# ── v4.1 ────────────────────────────────────────────────────────────────────────
function _read_entities_v4!(acc, io)
    hdr = split(strip(readline(io)))
    nP, nC, nS, nV = parse.(Int, hdr[1:4])
    _read_entity_block!(acc.ep, io, nP, 0)   # points
    _read_entity_block!(acc.ep, io, nC, 1)   # curves
    _read_entity_block!(acc.ep, io, nS, 2)   # surfaces
    _read_entity_block!(acc.ep, io, nV, 3)   # volumes
    _expect_end(io, "\$EndEntities")
end

function _read_entity_block!(ep, io, count, dim)
    for _ in 1:count
        p = split(strip(readline(io)))
        tag = parse(Int, p[1])
        # points: tag x y z numPhysicalTags [physicalTag...] (no max bbox)
        # dim≥1: tag minX minY minZ maxX maxY maxZ numPhysicalTags [physicalTag...] ...
        idx = dim == 0 ? 5 : 8
        nphys = parse(Int, p[idx])
        if nphys >= 1
            ep[(dim, tag)] = parse(Int, p[idx+1])   # first physical group tag
        end
    end
end

function _read_nodes_v4!(acc, io)
    hdr = split(strip(readline(io)))
    numBlocks = parse(Int, hdr[1]); numNodes = parse(Int, hdr[2])
    sizehint!(acc.node_x, numNodes); sizehint!(acc.node_y, numNodes); sizehint!(acc.node_z, numNodes)
    for _ in 1:numBlocks
        bh = split(strip(readline(io)))
        nInBlock = parse(Int, bh[4])
        tags = Vector{Int}(undef, nInBlock)
        for i in 1:nInBlock
            tags[i] = parse(Int, strip(readline(io)))
        end
        for i in 1:nInBlock
            c = split(strip(readline(io)))
            push!(acc.node_x, parse(Float64, c[1]))
            push!(acc.node_y, parse(Float64, c[2]))
            push!(acc.node_z, parse(Float64, c[3]))
            acc.tag2idx[tags[i]] = length(acc.node_x)
        end
    end
    _expect_end(io, "\$EndNodes")
end

function _read_elements_v4!(acc, io)
    ep = acc.ep
    hdr = split(strip(readline(io)))
    numBlocks = parse(Int, hdr[1])
    for _ in 1:numBlocks
        bh = split(strip(readline(io)))
        edim = parse(Int, bh[1]); etag = parse(Int, bh[2])
        etype = parse(Int, bh[3]); nInBlock = parse(Int, bh[4])
        phys = get(ep, (edim, etag), 0)
        for _ in 1:nInBlock
            p = split(strip(readline(io)))
            _push_element!(acc, etype, phys, @view p[2:end])   # p[1] = element tag
        end
    end
    _expect_end(io, "\$EndElements")
end

# ── shared element push (relabel tags → compact indices) ────────────────────────
function _push_element!(acc, etype::Int, phys::Int, nodetoks)
    haskey(_NN, etype) || return   # ignore element types we don't model (higher order, quads…)
    idx(t) = acc.tag2idx[parse(Int, t)]
    if etype == MSH_LINE
        a = idx(nodetoks[1]); b = idx(nodetoks[2])
        push!(acc.segs, (Int32(a), Int32(b))); push!(acc.seg_tag, Int32(phys))
    elseif etype == MSH_TRI
        a = idx(nodetoks[1]); b = idx(nodetoks[2]); c = idx(nodetoks[3])
        push!(acc.tris, (Int32(a), Int32(b), Int32(c))); push!(acc.tri_tag, Int32(phys))
    elseif etype == MSH_TET || etype == MSH_TET2
        # linear tet, or quadratic tet read as its 4 corner vertices (first 4 nodes)
        a = idx(nodetoks[1]); b = idx(nodetoks[2]); c = idx(nodetoks[3]); d = idx(nodetoks[4])
        push!(acc.tets, (Int32(a), Int32(b), Int32(c), Int32(d))); push!(acc.tet_tag, Int32(phys))
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
    open(path, "w") do io
        if version < 3
            _write_msh_v2(io, mesh, physical_names)
        else
            _write_msh_v4(io, mesh, physical_names)
        end
    end
    return path
end

function _write_physical_names(io, physical_names)
    isempty(physical_names) && return
    println(io, "\$PhysicalNames")
    println(io, length(physical_names))
    for ((dim, tag), name) in sort(collect(physical_names); by=x->x[1])
        println(io, dim, " ", tag, " \"", name, "\"")
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

    # Entities: one entity per (dim, physical tag) block, declaring its physical group.
    _write_entities_v4(io, m, segblocks, triblocks, tetblocks)

    # Nodes: a single block (dim 3 entity tag 0) listing all nodes.
    println(io, "\$Nodes")
    println(io, 1, " ", nnodes(m), " ", nnodes(m) == 0 ? 0 : 1, " ", nnodes(m))
    println(io, 3, " ", 0, " ", 0, " ", nnodes(m))
    @inbounds for i in 1:nnodes(m); println(io, i); end
    @inbounds for i in 1:nnodes(m)
        p = node(m, i); @printf(io, "%.17g %.17g %.17g\n", p[1], p[2], p[3])
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

function _write_entities_v4(io, m, segblocks, triblocks, tetblocks)
    # entity tag == physical tag; each block becomes one entity with that phys group.
    println(io, "\$Entities")
    println(io, 0, " ", length(segblocks), " ", length(triblocks), " ", length(tetblocks))
    for (tag, _) in segblocks
        # curve: tag minx miny minz maxx maxy maxz numPhys phys... numBounding
        println(io, tag, " 0 0 0 0 0 0 ", tag == 0 ? 0 : 1, tag == 0 ? "" : " $tag", " 0")
    end
    for (tag, _) in triblocks
        println(io, tag, " 0 0 0 0 0 0 ", tag == 0 ? 0 : 1, tag == 0 ? "" : " $tag", " 0")
    end
    for (tag, _) in tetblocks
        println(io, tag, " 0 0 0 0 0 0 ", tag == 0 ? 0 : 1, tag == 0 ? "" : " $tag", " 0")
    end
    println(io, "\$EndEntities")
end

function _write_elem_blocks_v4(io, cells::Matrix{Int32}, blocks, dim::Int, etype::Int, eid)
    for (tag, cols) in blocks
        println(io, dim, " ", tag, " ", etype, " ", length(cols))
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
    isbin = _stl_is_binary(path)
    tris_xyz = isbin ? _read_stl_binary(path) : _read_stl_ascii(path)
    return _weld_triangles(tris_xyz, merge_tol)
end

function _stl_is_binary(path)
    open(path, "r") do io
        head = read(io, min(filesize(path), 84))
        # ASCII STL starts with "solid"; but some binary files do too, so also
        # check the 80-byte-header + uint32 triangle count against file size.
        if length(head) >= 84
            ntri = reinterpret(UInt32, head[81:84])[1]
            expected = 84 + 50 * Int(ntri)
            filesize(path) == expected && return true
        end
        return !(length(head) >= 5 && String(head[1:5]) == "solid" && _looks_ascii(head))
    end
end

_looks_ascii(bytes) = all(b -> b == 0x09 || b == 0x0a || b == 0x0d || (0x20 <= b <= 0x7e), bytes)

function _read_stl_ascii(path)
    tris = NTuple{9,Float64}[]
    verts = NTuple{3,Float64}[]
    for line in eachline(path)
        s = split(strip(line))
        isempty(s) && continue
        if s[1] == "vertex"
            push!(verts, (parse(Float64,s[2]), parse(Float64,s[3]), parse(Float64,s[4])))
            if length(verts) == 3
                push!(tris, (verts[1]..., verts[2]..., verts[3]...))
                empty!(verts)
            end
        end
    end
    return tris
end

function _read_stl_binary(path)
    open(path, "r") do io
        skip(io, 80)
        ntri = Int(read(io, UInt32))
        tris = Vector{NTuple{9,Float64}}(undef, ntri)
        for t in 1:ntri
            skip(io, 12)   # normal vector (3× Float32) — recomputed on demand
            v = ntuple(_ -> Float64(read(io, Float32)), 9)
            tris[t] = v
            skip(io, 2)    # attribute byte count
        end
        return tris
    end
end

function _weld_triangles(tris_xyz::Vector{NTuple{9,Float64}}, reltol::Real)
    isempty(tris_xyz) && return Mesh(Matrix{Float64}(undef, 3, 0))
    # bbox diagonal → absolute tolerance for quantization
    lo = (Inf,Inf,Inf); hi = (-Inf,-Inf,-Inf)
    for t in tris_xyz, k in (1,4,7)
        p = (t[k], t[k+1], t[k+2])
        lo = (min(lo[1],p[1]),min(lo[2],p[2]),min(lo[3],p[3]))
        hi = (max(hi[1],p[1]),max(hi[2],p[2]),max(hi[3],p[3]))
    end
    diag = sqrt((hi[1]-lo[1])^2 + (hi[2]-lo[2])^2 + (hi[3]-lo[3])^2)
    tol = max(diag * reltol, eps(Float64))
    inv = 1.0 / tol
    keyof(p) = (round(Int, p[1]*inv), round(Int, p[2]*inv), round(Int, p[3]*inv))
    idxmap = Dict{NTuple{3,Int},Int32}()
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    getid(p) = get!(idxmap, keyof(p)) do
        push!(xs,p[1]); push!(ys,p[2]); push!(zs,p[3]); Int32(length(xs))
    end
    tris = Matrix{Int32}(undef, 3, length(tris_xyz))
    ntri = 0
    for t in tris_xyz
        a = getid((t[1],t[2],t[3])); b = getid((t[4],t[5],t[6])); c = getid((t[7],t[8],t[9]))
        (a==b || b==c || a==c) && continue   # drop degenerate welded triangles
        ntri += 1
        tris[1,ntri]=a; tris[2,ntri]=b; tris[3,ntri]=c
    end
    coords = Matrix{Float64}(undef, 3, length(xs))
    @inbounds for i in eachindex(xs); coords[1,i]=xs[i]; coords[2,i]=ys[i]; coords[3,i]=zs[i]; end
    return Mesh(coords; tris=tris[:, 1:ntri])
end

# ════════════════════════════════════════════════════════════════════════════════
# .geo parameter / structure scan (no OCC evaluation — Stage 5)
# ════════════════════════════════════════════════════════════════════════════════

"""
    GeoParams

What `read_geo_params` can extract from a `.geo` without a geometry kernel:
`mesh_size_min/max`, `random_seed`, and `physical_groups` (a `(dim, tag) => name`
map, dim ∈ {0:point,1:curve,2:surface,3:volume}).
"""
struct GeoParams
    mesh_size_min::Float64
    mesh_size_max::Float64
    random_seed::Int
    physical_groups::Dict{Tuple{Int,Int},String}
end

const _PHYS_DIM = Dict("Point"=>0, "Curve"=>1, "Line"=>1, "Surface"=>2, "Volume"=>3)

"""
    read_geo_params(path) -> GeoParams

Scan a gmsh `.geo`/`.geo_unrolled` for mesh-sizing options and Physical group
declarations. Deliberately does **not** evaluate CSG/Boolean geometry.
"""
function read_geo_params(path::AbstractString)
    smin = NaN; smax = NaN; seed = 0
    groups = Dict{Tuple{Int,Int},String}()
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "//")) && continue
        _scan_opt!(line, "Mesh.MeshSizeMin", v -> (smin = v)) ||
        _scan_opt!(line, "Mesh.MeshSizeMax", v -> (smax = v))
        if occursin("Mesh.RandomSeed", line)
            m = match(r"Mesh\.RandomSeed\s*=\s*([0-9]+)", line)
            m !== nothing && (seed = parse(Int, m.captures[1]))
        end
        # Physical Volume("air", 1) = {...};  /  Physical Surface("s", 3) = {...};
        pm = match(r"Physical\s+(Point|Curve|Line|Surface|Volume)\s*\(\s*\"([^\"]*)\"\s*,\s*([0-9]+)\s*\)", line)
        if pm !== nothing
            dim = _PHYS_DIM[pm.captures[1]]
            name = pm.captures[2]
            tag = parse(Int, pm.captures[3])
            groups[(dim, tag)] = name
        end
    end
    return GeoParams(smin, smax, seed, groups)
end

function _scan_opt!(line, key, setter)
    occursin(key, line) || return false
    m = match(Regex(replace(key, "."=>"\\.") * raw"\s*=\s*([-+0-9.eE]+)"), line)
    m !== nothing && setter(parse(Float64, m.captures[1]))
    return true
end

end # module IO
