# SPDX-FileCopyrightText: 2016-2020 UCLouvain
# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: GPL-2.0-or-later
# Contains HXT adaptations from Gmsh 4.15.2 contrib/hxt/tetMesh/src (hxt_vertices.c,
# hxt_tetDelaunay.c, hxt_tetSync.c, hxt_tetUtils.c, hxt_tetNodalSize.c, hxt_tetRefine.c)
# and contrib/hxt/core (hxt_tools.h, hxt_sort.h, hxt_bbox.c). Upstream author: Célestin Marot.
# Julia adaptation/changes: Tessella contributors, 2026.
# See LICENSES/README.md and NOTICES.md#hxt.

"""
    Hxt

Port of Gmsh's HXT tetrahedral mesher (`Mesh.Algorithm3D = 10`), one thread,
step for step: the same vertex order (Moore/Hilbert curve, the same pseudo-random
shuffle), the same Delaunay kernel with its cell recycling, so that every array
is the one HXT builds. Indices follow HXT: vertices and cells are 0-based in the
stored values (`node`, `neigh`), `neigh[4t+f]` is `4·(neighbour)+facet`.
"""
module Hxt

import ..Predicates: orient3, insphere

const GHOST = typemax(UInt32)
const NO_ADJACENT = typemax(UInt64)
const COLOR_OUT = typemax(UInt32)

# HXTStatus values the kernel distinguishes
const ST_TRUE = Int32(0)
const ST_FALSE = Int32(1)
const ST_TRYAGAIN = Int32(-1026)
const ST_CONFLICT = Int32(-1027)
const ST_DOUBLE_PT = Int32(-1030)

const DELETED_BUFFER_SIZE = 8182
const SMALLEST_PASS = 2048

const DELETED_MASK = 0x0040
const ALL_EDGES_MASK = 0x003f
const ALL_FACETS_MASK = 0x0f00

include("HxtTables.jl")

# ── core tools (hxt_tools.h) ────────────────────────────────────────────────────
@inline lcg_u32(state::Base.RefValue{UInt32}) =
    (state[] = UInt32((UInt64(state[]) * UInt64(279470273)) % UInt64(0xfffffffb)))
@inline lcg_f64(state::Base.RefValue{UInt32}) = Float64(lcg_u32(state)) * (1.0 / 4294967296.0)

@inline function hash32(x::UInt32)
    x = ((x >> 16) ⊻ x) * UInt32(0x45d9f3b)
    x = ((x >> 16) ⊻ x) * UInt32(0x45d9f3b)
    return (x >> 16) ⊻ x
end

@inline function hash64(x::UInt64)
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

# ── the mesh (HXTMesh) ─────────────────────────────────────────────────────────
mutable struct NodalSizes
    array::Vector{Float64}
    min::Float64
    max::Float64
    factor::Float64
    enabled::Bool
end

mutable struct Mesh
    coord::Vector{Float64}          # 4 per vertex; the 4th is HXT's padding union
    nv::Int
    node::Vector{UInt32}
    neigh::Vector{UInt64}
    color::Union{Nothing,Vector{UInt32}}
    flag::Vector{UInt16}
    ntet::Int
    tri_node::Vector{UInt32}
    tri_color::Vector{UInt32}
    line_node::Vector{UInt32}
    line_color::Vector{UInt32}
    point_node::Vector{UInt32}
    point_color::Vector{UInt32}
    # brep
    surfaces_per_volume::Vector{Vector{UInt32}}
    curves_per_surface::Vector{Vector{UInt32}}
end

Mesh(coord::Vector{Float64}) = Mesh(coord, length(coord) ÷ 4, UInt32[], UInt64[], nothing, UInt16[], 0,
    UInt32[], UInt32[], UInt32[], UInt32[], UInt32[], UInt32[], Vector{UInt32}[], Vector{UInt32}[])

# padding of vertex `v` (0-based) as the union's members
@inline pad(m::Mesh, v) = @inbounds reinterpret(UInt64, m.coord[4Int(v)+4])
@inline setpad!(m::Mesh, v, x::UInt64) = (@inbounds m.coord[4Int(v)+4] = reinterpret(Float64, x); nothing)
@inline point(m::Mesh, v) = @inbounds (m.coord[4Int(v)+1], m.coord[4Int(v)+2], m.coord[4Int(v)+3])

@inline tnode(m::Mesh, i) = @inbounds m.node[Int(i)+1]                  # i = 4t+j
@inline settnode!(m::Mesh, i, v) = (@inbounds m.node[Int(i)+1] = UInt32(v); nothing)
@inline tneigh(m::Mesh, i) = @inbounds m.neigh[Int(i)+1]
@inline settneigh!(m::Mesh, i, v) = (@inbounds m.neigh[Int(i)+1] = UInt64(v); nothing)
@inline tflag(m::Mesh, t) = @inbounds m.flag[Int(t)+1]
@inline settflag!(m::Mesh, t, v) = (@inbounds m.flag[Int(t)+1] = UInt16(v); nothing)

@inline deleted(m::Mesh, t) = (tflag(m, t) & DELETED_MASK) != 0
@inline setdeleted!(m::Mesh, t) = settflag!(m, t, tflag(m, t) | DELETED_MASK)
@inline unsetdeleted!(m::Mesh, t) = settflag!(m, t, tflag(m, t) & ~DELETED_MASK)
@inline facet_constraint(m::Mesh, t, f) = tflag(m, t) & (UInt16(0x100) << f)
@inline edge_constraint(m::Mesh, t, e) = tflag(m, t) & (UInt16(1) << e)

@inline node0_of_facet(f) = (f + 1) & 3
@inline node1_of_facet(f) = (f & 2) ⊻ 3
@inline node2_of_facet(f) = (f + 3) & 2

function reserve!(m::Mesh, total::Integer)
    if total > length(m.flag)
        resize!(m.flag, total)
        m.color === nothing || resize!(m.color, total)
        resize!(m.node, 4total)
        resize!(m.neigh, 4total)
    end
    return m
end

# ── bounding box and the predicates' static filter (hxt_bbox.c, exactinit) ─────
struct Bbox
    min::NTuple{3,Float64}
    max::NTuple{3,Float64}
end

function bbox(m::Mesh)
    lo = [floatmax(Float64), floatmax(Float64), floatmax(Float64)]
    hi = [-floatmax(Float64), -floatmax(Float64), -floatmax(Float64)]
    @inbounds for i in 0:m.nv-1, k in 1:3
        x = m.coord[4i+k]
        x < lo[k] && (lo[k] = x)
        x > hi[k] && (hi[k] = x)
    end
    return Bbox((lo[1], lo[2], lo[3]), (hi[1], hi[2], hi[3]))
end

o3d_static_filter(b::Bbox) = 5.1107127829973299e-15 * (b.max[1] - b.min[1]) *
    (b.max[2] - b.min[2]) * (b.max[3] - b.min[3])

# exactinit's insphere filter, with its partial sort of the extents
function isp_static_filter(b::Bbox)
    maxx = b.max[1] - b.min[1]; maxy = b.max[2] - b.min[2]; maxz = b.max[3] - b.min[3]
    maxx > maxz && ((maxx, maxz) = (maxz, maxx))
    if maxy > maxz
        (maxy, maxz) = (maxz, maxy)
    elseif maxy < maxx
        (maxy, maxx) = (maxx, maxy)
    end
    return 1.2466136531027298e-13 * maxx * maxy * maxz * (maxz * maxz)
end

# orient3d and insphere behind exactinit's static filters: the floating-point
# determinant decides when it is beyond the filter, the exact predicate
# otherwise (the sign is exact either way)
@inline function orient3_static(a, b, c, d, filter::Float64)
    adx = a[1] - d[1]; ady = a[2] - d[2]; adz = a[3] - d[3]
    bdx = b[1] - d[1]; bdy = b[2] - d[2]; bdz = b[3] - d[3]
    cdx = c[1] - d[1]; cdy = c[2] - d[2]; cdz = c[3] - d[3]
    det = adz * (bdx * cdy - cdx * bdy) + bdz * (cdx * ady - adx * cdy) + cdz * (adx * bdy - bdx * ady)
    det > filter && return 1
    -det > filter && return -1
    return orient3(a, b, c, d)
end

@inline function insphere_static(a, b, c, d, e, filter::Float64)
    aex = a[1] - e[1]; bex = b[1] - e[1]; cex = c[1] - e[1]; dex = d[1] - e[1]
    aey = a[2] - e[2]; bey = b[2] - e[2]; cey = c[2] - e[2]; dey = d[2] - e[2]
    aez = a[3] - e[3]; bez = b[3] - e[3]; cez = c[3] - e[3]; dez = d[3] - e[3]
    ab = aex * bey - bex * aey; bc = bex * cey - cex * bey
    cd = cex * dey - dex * cey; da = dex * aey - aex * dey
    ac = aex * cey - cex * aey; bd = bex * dey - dex * bey
    abc = aez * bc - bez * ac + cez * ab
    bcd = bez * cd - cez * bd + dez * bc
    cda = cez * da + dez * ac + aez * cd
    dab = dez * ab + aez * bd + bez * da
    alift = aex * aex + aey * aey + aez * aez
    blift = bex * bex + bey * bey + bez * bez
    clift = cex * cex + cey * cey + cez * cez
    dlift = dex * dex + dey * dey + dez * dez
    det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd)
    det > filter && return 1
    -det > filter && return -1
    return insphere(a, b, c, d, e)
end

# ── Moore / Hilbert indices (hxt_vertices.c) ───────────────────────────────────
@inline function zorder(x::UInt32, y::UInt32, z::UInt32)
    M = _HXT_MORTON7
    @inbounds return M[(x & 0x7f)+1] | M[(y & 0x7f)+1] << 1 | M[(z & 0x7f)+1] << 2 |
        M[(x >> 7 & 0x7f)+1] << 21 | M[(y >> 7 & 0x7f)+1] << 22 | M[(z >> 7 & 0x7f)+1] << 23 |
        M[(x >> 14 & 0x7f)+1] << 42 | M[(y >> 14 & 0x7f)+1] << 43 | M[(z >> 14 & 0x7f)+1] << 44
end

# C's double -> uint32_t conversion of an in-range value
@inline c_u32(x::Float64) = unsafe_trunc(UInt32, unsafe_trunc(Int64, x))

"""
    moore!(m, bbox, n, mid01=(0.5,0.5,0.5))

`hxtMoore`: write the Moore-curve index of vertices `0:n-1` into their padding.
The vertices after the last whole block of 32 have their third integer
coordinate divided by four, exactly as HXT does.
"""
function moore!(m::Mesh, b::Bbox, n::Integer, mid01=(0.5, 0.5, 0.5))
    nmax = UInt32(0x200000)
    width_max = 0.0
    for i in 1:3
        w = b.max[i] - b.min[i]
        w > width_max && (width_max = w)
    end
    middle = zeros(3); f0 = zeros(3); f1 = zeros(3); sub1 = zeros(3)
    for i in 1:3
        if width_max > 1.5 * (b.max[i] - b.min[i])
            if mid01[i] >= 0.5
                xmin = b.min[i]; xmax = b.min[i] + width_max
            else
                xmax = b.max[i]; xmin = b.max[i] - width_max
            end
        else
            xmin = b.min[i]; xmax = b.max[i]
        end
        middle[i] = mid01[i] * (xmax - xmin) + xmin
        f0[i] = Float64(nmax ÷ 2) / (middle[i] - xmin)
        f1[i] = Float64(nmax ÷ 2) / (xmax - middle[i])
        sub1[i] = 2 * middle[i] - xmax
        while c_u32((xmax - sub1[i]) * f1[i]) >= nmax
            f1[i] = prevfloat(f1[i])
        end
    end
    full = (Int(n) ÷ 32) * 32
    @inbounds for v in 0:Int(n)-1
        xyz = ntuple(3) do k
            c = m.coord[4v+k]
            c < middle[k] ? c_u32((c - b.min[k]) * f0[k]) : c_u32((c - sub1[k]) * f1[k])
        end
        x, y, z = xyz
        v >= full && (z = z ÷ UInt32(4))
        zo = zorder(x, y, z)
        transform = UInt64(_HXT_MOORE_HILBERT2[Int(zo >> 54)+1])
        out = transform & 0x1ff
        for iter in 45:-9:0
            transform = (transform & ~UInt64(0x1ff)) | (zo >> iter & 0x1ff)
            transform = UInt64(_HXT_HILBERT3[Int(transform)+1])
            out = (out << 9) | (transform & 0x1ff)
        end
        setpad!(m, v, out)
    end
    return m
end

# ── the vertices to insert (HXTNodeInfo) ────────────────────────────────────────
struct NodeInfo
    hilbert::UInt64
    node::UInt32
    status::Int32
end

# hxtNodeInfoShuffle: a stable sort by the 32-bit hash of the position
function shuffle_info!(info::AbstractVector{NodeInfo})
    for i in eachindex(info)
        x = info[i]
        info[i] = NodeInfo(UInt64(hash32(UInt32(i - first(eachindex(info))))), x.node, x.status)
    end
    sort!(info; by=x -> x.hilbert % UInt32, alg=Base.Sort.DEFAULT_STABLE)
    return info
end

# hxtNodeInfoSort: a stable sort by the full Moore index (HXT's LSB radix sort)
sort_info!(info::AbstractVector{NodeInfo}) =
    sort!(info; by=x -> x.hilbert, alg=Base.Sort.DEFAULT_STABLE)

# hxtVerticesShuffle and hxtVerticesSort move the vertices themselves (and
# their padding) in `coord[4first+1 : 4(first+n)]`.
function permute_vertices!(m::Mesh, first::Integer, keys::Vector{UInt64})
    n = length(keys)
    order = sortperm(keys; alg=Base.Sort.DEFAULT_STABLE)
    block = m.coord[4Int(first)+1:4(Int(first)+n)]
    @inbounds for (to, from) in enumerate(order), k in 1:4
        m.coord[4(Int(first)+to-1)+k] = block[4(from-1)+k]
    end
    return m
end

function shuffle_vertices!(m::Mesh, first::Integer, n::Integer)
    keys = UInt64[UInt64(hash32(UInt32(i))) for i in 0:n-1]
    for i in 0:n-1
        setpad!(m, first + i, keys[i+1])
    end
    return permute_vertices!(m, first, keys)
end

sort_vertices!(m::Mesh, first::Integer, n::Integer) =
    permute_vertices!(m, first, UInt64[pad(m, first + i) for i in 0:n-1])

# ── Delaunay options (HXTDelaunayOptions) ───────────────────────────────────────
mutable struct DelaunayOptions
    bbox::Bbox
    nodal_sizes::Union{Nothing,NodalSizes}
    num_vertices_in_mesh::Int
    insertion_first::Int
    partitionability::Float64
    perfect_delaunay::Bool
    allow_outer_insertion::Bool
end

# ── per-thread work arrays (TetLocal) ──────────────────────────────────────────
struct BallFacet
    neigh::UInt64
    node::NTuple{3,UInt32}
    flag::UInt16
end

mutable struct Local
    ball::Vector{BallFacet}
    deleted::Vector{UInt64}
    # open-addressing table pairing the new cells' facets by edge, valid for
    # one fill (`epoch`)
    edge_key::Vector{UInt64}
    edge_facet::Vector{UInt64}
    edge_epoch::Vector{UInt32}
    epoch::UInt32
    # computeAdjacenciesFast: cavity-local vertex numbers and the 32x32 edge map
    vertex_id::Vector{UInt8}
    local_nodes::Vector{NTuple{3,UInt8}}
    edge_map::Vector{UInt64}
    # exactinit's static filters for this bounding box
    o3d::Float64
    isp::Float64
end
Local(b::Bbox) = Local(BallFacet[], UInt64[], zeros(UInt64, 256), zeros(UInt64, 256), zeros(UInt32, 256), 0,
    UInt8[], NTuple{3,UInt8}[], zeros(UInt64, 32 * 32), o3d_static_filter(b), isp_static_filter(b))

# ── initial tetrahedron and its four ghosts (hxtTetrahedraInit) ────────────────
function tetrahedra_init!(m::Mesh, info::AbstractVector{NodeInfo}, n::Integer, filter::Float64)
    n < 4 && throw(ArgumentError("HXT: cannot mesh less than four points"))
    reserve!(m, max(DELETED_BUFFER_SIZE + 8n, 10n))
    nodeof(i) = info[i+1].node
    orientation = 0
    i = 0; j = 1; k = 2; l = 3
    i = 0
    while orientation == 0 && i < n - 3
        d = point(m, nodeof(i))
        j = i + 1
        while orientation == 0 && j < n - 2
            c = point(m, nodeof(j))
            cdx = c[1] - d[1]; cdy = c[2] - d[2]; cdz = c[3] - d[3]
            if !(cdx == 0.0 && cdy == 0.0 && cdz == 0.0)
                k = j + 1
                while orientation == 0 && k < n - 1
                    b = point(m, nodeof(k))
                    bdx = b[1] - d[1]; bdy = b[2] - d[2]; bdz = b[3] - d[3]
                    crossx = bdy * cdz - bdz * cdy
                    crossy = bdz * cdx - bdx * cdz
                    crossz = bdx * cdy - bdy * cdx
                    if !(crossx == 0.0 && crossy == 0.0 && crossz == 0.0)
                        l = k + 1
                        while orientation == 0 && l < n
                            a = point(m, nodeof(l))
                            adx = a[1] - d[1]; ady = a[2] - d[2]; adz = a[3] - d[3]
                            det = adx * crossx + ady * crossy + adz * crossz
                            orientation = (det > filter) - (det < -filter)
                            l += 1
                        end
                    end
                    k += 1
                end
            end
            j += 1
        end
        i += 1
    end
    l -= 1; k -= 1; j -= 1; i -= 1
    orientation == 0 && throw(ArgumentError("HXT: all vertices are coplanar or nearly coplanar"))
    for (slot, at) in ((0, i), (1, j), (2, k), (3, l))
        tmp = info[at+1]
        info[at+1] = info[slot+1]
        info[slot+1] = NodeInfo(tmp.hilbert, tmp.node, ST_TRUE)
    end
    i = 0; j = 1; k = 2; l = 3
    orientation > 0 && ((i, j) = (j, i))
    N(x) = info[x+1].node
    neigh = (19, 15, 11, 7, 18, 10, 13, 3, 17, 14, 5, 2, 16, 6, 9, 1, 12, 8, 4, 0)
    nodes = (N(l), N(k), N(j), N(i), N(l), N(j), N(k), GHOST, N(l), N(k), N(i), GHOST,
        N(l), N(i), N(j), GHOST, N(k), N(j), N(i), GHOST)
    for s in 0:19
        settneigh!(m, s, neigh[s+1]); settnode!(m, s, nodes[s+1])
    end
    m.ntet = 5
    if m.color !== nothing
        for t in 0:4; m.color[t+1] = COLOR_OUT; end
    end
    for t in 0:4; settflag!(m, t, 0); end
    return m
end

function compute_passes(inserted::Integer, n::Integer, partitionability::Float64)
    alpha = 1.0                         # 1/(maxPartitions^2), one partition
    inserted = c_u32(Float64(UInt32(inserted)) * ((1.0 - alpha) * partitionability + alpha))
    passes = zeros(UInt32, 12)
    npasses = 0
    passes[1] = UInt32(n)
    for i in 0:9
        if passes[i+1] < SMALLEST_PASS || passes[i+1] ÷ UInt32(8) < inserted
            passes[i+2] = 0
            npasses = i + 1
            break
        end
        passes[i+2] = c_u32(passes[i+1] / 7.5)
    end
    for i in 0:npasses÷2
        passes[i+1], passes[npasses-i+1] = passes[npasses-i+1], passes[i+1]
    end
    return passes, npasses
end

# ── predicates (tetInsphere with the ghost vertex and symbolic perturbation) ───
function symbolic_perturbation(indices::NTuple{5,UInt32}, pts::NTuple{5,NTuple{3,Float64}})
    idx = collect(indices); pt = collect(pts)
    swaps = 0; n = 5
    while true
        count = 0
        n -= 1
        for it in 1:n
            if idx[it] > idx[it+1]
                pt[it], pt[it+1] = pt[it+1], pt[it]
                idx[it], idx[it+1] = idx[it+1], idx[it]
                count += 1
            end
        end
        swaps += count
        count > 0 || break
    end
    ori = orient3(pt[2], pt[3], pt[4], pt[5])
    if ori != 0
        isodd(swaps) && (ori = -ori)
        return ori
    end
    ori = -orient3(pt[1], pt[3], pt[4], pt[5])
    isodd(swaps) && (ori = -ori)
    return ori
end

# `base` is 4·(cell), the index of its first node
function tet_insphere(m::Mesh, L::Local, base::UInt64, vta::UInt32)
    n0 = tnode(m, base); n1 = tnode(m, base + 1); n2 = tnode(m, base + 2); n3 = tnode(m, base + 3)
    a = point(m, n0); b = point(m, n1); c = point(m, n2); e = point(m, vta)
    if n3 == GHOST
        det = orient3_static(a, b, c, e, L.o3d)
        det != 0 && return det
        opposite = tnode(m, tneigh(m, base + 3))
        o = point(m, opposite)
        sign = insphere_static(a, b, c, o, e, L.isp)
        sign == 0 && (sign = symbolic_perturbation((n0, n1, n2, opposite, vta), (a, b, c, o, e)))
        return -sign
    end
    d = point(m, n3)
    sign = insphere_static(a, b, c, d, e, L.isp)
    sign == 0 && (sign = symbolic_perturbation((n0, n1, n2, n3, vta), (a, b, c, d, e)))
    return sign
end

# ── walk (walking2Cavity, one partition) ───────────────────────────────────────
function walk!(m::Mesh, L::Local, cur::Base.RefValue{UInt64}, vta::UInt32)
    next = cur[]
    if tnode(m, 4next + 3) == GHOST
        next = tneigh(m, 4next + 3) ÷ 4
    end
    p = point(m, vta)
    entering = 4
    seed = Ref(UInt32(1))
    while true
        base = 4next
        index = 4
        outside = false
        random = lcg_u32(seed)
        for jj in 0:3
            f = Int((UInt32(jj) + random) % UInt32(4))
            f == entering && continue
            a = point(m, tnode(m, base + node0_of_facet(f)))
            b = point(m, tnode(m, base + node1_of_facet(f)))
            c = point(m, tnode(m, base + node2_of_facet(f)))
            if orient3_static(p, a, b, c, L.o3d) < 0
                outside = true
                nb = tneigh(m, base + f)
                if tnode(m, nb) == GHOST
                    cur[] = nb ÷ 4
                    return ST_TRUE
                end
                index = f
                break
            end
        end
        if index == 4
            outside && return ST_CONFLICT
            a = point(m, tnode(m, base)); b = point(m, tnode(m, base + 1))
            c = point(m, tnode(m, base + 2)); d = point(m, tnode(m, base + 3))
            if (orient3_static(a, b, c, p, L.o3d) >= 0) + (orient3_static(a, b, p, d, L.o3d) >= 0) +
               (orient3_static(a, p, c, d, L.o3d) >= 0) + (orient3_static(p, b, c, d, L.o3d) >= 0) > 2
                cur[] = next
                return ST_DOUBLE_PT
            end
            cur[] = next
            return ST_TRUE
        end
        nb = tneigh(m, base + index)
        entering = Int(nb & 3)
        next = nb ÷ 4
    end
end

# ── cavity (diggingACavity, perfect Delaunay) ──────────────────────────────────
@inline function bnd_push!(L::Local, flag, n1, n2, n3, neigh)
    push!(L.ball, BallFacet(UInt64(neigh), (UInt32(n1), UInt32(n2), UInt32(n3)), UInt16(flag)))
end

# the boundary flag a cavity facet carries to the new cell (facet f of `t` becomes
# facet 0 of the new cell, its edges the new cell's edges 0..2), as diggingACavity
@inline function dig_flag(m::Mesh, t, f)
    fc(k) = facet_constraint(m, t, k); ec(k) = edge_constraint(m, t, k)
    f == 0 && return tflag(m, t) & UInt16(0x107)
    f == 1 && return (fc(1) >> 1) | (ec(3) >> 3) | (ec(0) << 1) | (ec(4) >> 2)
    f == 2 && return (fc(2) >> 2) | (ec(1) >> 1) | (ec(3) >> 2) | (ec(5) >> 3)
    return (fc(3) >> 3) | (ec(4) >> 4) | (ec(2) >> 1) | (ec(5) >> 3)
end

# diggingACavity: breadth-first over the cells whose circumsphere holds vta.
# Returns whether a cell of the cavity has a constrained edge.
function dig!(m::Mesh, L::Local, first::UInt64, vta::UInt32, perfect::Bool)
    push!(L.deleted, first)
    setdeleted!(m, first)
    empty!(L.ball)
    edge_constraint = false
    start = length(L.deleted)
    while start <= length(L.deleted)
        cur = L.deleted[start]
        base = 4cur
        edge_constraint |= !perfect && (tflag(m, cur) & ALL_EDGES_MASK) != 0
        facet_constrained = !perfect && (tflag(m, cur) & ALL_FACETS_MASK) != 0
        nd = (tnode(m, base), tnode(m, base + 1), tnode(m, base + 2), tnode(m, base + 3))
        for (f, face) in ((0, (nd[2], nd[3], nd[4])), (1, (nd[3], nd[1], nd[4])),
                          (2, (nd[1], nd[2], nd[4])), (3, (nd[2], nd[1], nd[3])))
            nbi = tneigh(m, base + f)
            nb = nbi ÷ 4
            constrained = facet_constrained && facet_constraint(m, cur, f) != 0
            (!deleted(m, nb) || constrained) || continue
            if constrained || tet_insphere(m, L, 4nb, vta) >= 0
                bnd_push!(L, dig_flag(m, cur, f), face[1], face[2], face[3], nbi)
            else
                push!(L.deleted, nb)
                setdeleted!(m, nb)
            end
        end
        start += 1
    end
    return edge_constraint
end

const _FACETS2EDGE = ((-1, 0, 1, 2), (0, -1, 3, 4), (1, 3, -1, 5), (2, 4, 5, -1))
const _EDGE2FACET_MIN = (0, 0, 0, 1, 1, 2)
const _EDGE2FACET_MAX = (1, 2, 3, 2, 3, 3)
@inline edge_of_facets(f1, f2) = _FACETS2EDGE[f1+1][f2+1]
@inline facets_of_edge(e) = (_EDGE2FACET_MIN[e+1], _EDGE2FACET_MAX[e+1])
const UNDELETE_MASK = 0x1000
@inline to_undelete(m::Mesh, t) = (tflag(m, t) & UNDELETE_MASK) != 0
@inline set_to_undelete!(m::Mesh, t) = settflag!(m, t, tflag(m, t) | UNDELETE_MASK)
@inline unset_to_undelete!(m::Mesh, t) = settflag!(m, t, tflag(m, t) & ~UNDELETE_MASK)

# turn once around the edge entering `cur` by `in`, leaving by `out`
@inline function turn(m::Mesh, cur::UInt64, in::Int, out::Int)
    newv = tnode(m, 4cur + in)
    nb = tneigh(m, 4cur + out)
    cur = nb ÷ 4
    in = Int(nb % 4)
    out = 0
    while out < 3 && tnode(m, 4cur + out) != newv
        out += 1
    end
    return cur, in, out, newv
end

# respectEdgeConstraint (with HXT's own index arithmetic, `coord + newV`, kept)
function respect_edge_constraint!(m::Mesh, L::Local, vta::UInt32, color::UInt32, previous::Int)
    undelete = false
    colors = m.color
    for i in previous+1:length(L.deleted)
        colors[Int(L.deleted[i])+1] = 0
    end
    for i in previous+1:length(L.deleted)
        del = L.deleted[i]
        for edge in 0:5
            (edge_constraint(m, del, edge) != 0 && (colors[Int(del)+1] & (UInt32(1) << edge)) == 0) || continue
            in, out = facets_of_edge(edge)
            safe = false
            cur = del
            while true
                cur, in, out, _ = turn(m, cur, in, out)
                if deleted(m, cur) && !to_undelete(m, cur)
                    colors[Int(cur)+1] |= UInt32(1) << edge_of_facets(in, out)
                else
                    safe = true
                end
                cur == del && break
            end
            if !safe
                in, out = facets_of_edge(edge)
                cur = del
                containing = L.deleted[previous+1]
                undelete_tet = NO_ADJACENT
                dist_max = 0.0
                v = point(m, vta)
                while true
                    cur, in, out, newv = turn(m, cur, in, out)
                    if cur != containing
                        other = tnode(m, 4cur + in)
                        dist = 0.0
                        for l in 0:2
                            mean = (m.coord[Int(newv)+l+1] + m.coord[Int(other)+l+1]) * 0.5
                            diff = mean - v[l+1]
                            dist += diff * diff
                        end
                        dist > dist_max && (undelete_tet = cur)   # HXT never raises dist_max
                    end
                    cur == del && break
                end
                set_to_undelete!(m, undelete_tet)
                undelete = true
            end
        end
    end
    for i in previous+1:length(L.deleted)
        colors[Int(L.deleted[i])+1] = color
    end
    return undelete
end

function bnd_push_facet!(L::Local, m::Mesh, t::UInt64, f::Int, neigh::UInt64)
    n = (tnode(m, 4t), tnode(m, 4t + 1), tnode(m, 4t + 2), tnode(m, 4t + 3))
    fc(k) = facet_constraint(m, t, k); ec(k) = edge_constraint(m, t, k)
    if f == 0
        bnd_push!(L, fc(0) | (ec(1) >> 1) | (ec(0) << 1) | ec(2), n[3], n[2], n[4], neigh)
    elseif f == 1
        bnd_push!(L, (fc(1) >> 1) | ec(0) | (ec(3) >> 2) | (ec(4) >> 2), n[1], n[3], n[4], neigh)
    elseif f == 2
        bnd_push!(L, (fc(2) >> 2) | (ec(3) >> 3) | ec(1) | (ec(5) >> 3), n[2], n[1], n[4], neigh)
    else
        bnd_push!(L, (fc(3) >> 3) | (ec(2) >> 2) | (ec(4) >> 3) | (ec(5) >> 3), n[1], n[2], n[3], neigh)
    end
end

@inline set_ball_neigh!(L::Local, i, neigh) =
    (f = L.ball[i]; L.ball[i] = BallFacet(UInt64(neigh), f.node, f.flag))

# reshapeCavityIfNeeded: give back cells until vta sees every boundary facet
function reshape_cavity!(m::Mesh, L::Local, vta::UInt32, previous::Int, undelete::Bool)
    p = point(m, vta)
    blind = 0
    if !undelete
        blind = -1
        for i in 1:length(L.ball)
            f = L.ball[i]
            f.node[3] == GHOST && continue
            if orient3_static(p, point(m, f.node[1]), point(m, f.node[2]), point(m, f.node[3]), L.o3d) >= 0
                blind = i - 1; break
            end
        end
        blind < 0 && return nothing
    end
    ntet = length(L.deleted) - previous
    tets = view(L.deleted, previous+1:length(L.deleted))
    index = Dict{UInt64,Int}(tets[i] => i - 1 for i in 1:ntet)
    faces = Vector{UInt64}(undef, 4ntet)
    cur = 0
    for i in 0:ntet-1, f in 0:3
        t = tets[i+1]
        nb = tneigh(m, 4t + f)
        if facet_constraint(m, t, f) != 0 || !deleted(m, nb ÷ 4)
            set_ball_neigh!(L, cur + 1, 4i + f)
            faces[4i+f+1] = cur
            cur += 1
        else
            faces[4i+f+1] = 4index[nb÷4] + nb % 4
        end
    end
    if undelete
        for i in 0:ntet-1
            t = tets[i+1]
            to_undelete(m, t) || continue
            boundary = false
            for f in 0:3
                nb = tneigh(m, 4t + f)
                if facet_constraint(m, t, f) != 0 || !deleted(m, nb ÷ 4)
                    boundary = true; break
                end
            end
            boundary && continue
            unset_to_undelete!(m, t)
            unsetdeleted!(m, t)
            tets[i+1] = NO_ADJACENT
            for f in 0:3
                out = faces[4i+f+1]
                faces[Int(out)+1] = length(L.ball)
                bnd_push_facet!(L, m, t, f, out)
            end
        end
    end
    cur = blind
    while cur < length(L.ball)
        in = L.ball[cur+1].neigh
        t = tets[Int(in ÷ 4)+1]
        if !to_undelete(m, t)
            f = L.ball[cur+1]
            if f.node[3] == GHOST
                cur += 1; continue
            end
            if orient3_static(p, point(m, f.node[1]), point(m, f.node[2]), point(m, f.node[3]), L.o3d) < 0
                cur += 1; continue
            end
        else
            unset_to_undelete!(m, t)
        end
        unsetdeleted!(m, t)
        tets[Int(in ÷ 4)+1] = NO_ADJACENT
        for f in 0:3
            nb = tneigh(m, 4t + f)
            if facet_constraint(m, t, f) != 0 || !deleted(m, nb ÷ 4)
                face = Int(faces[4Int(in÷4)+f+1])
                if face < cur
                    cur -= 1
                    if face != cur
                        L.ball[face+1] = L.ball[cur+1]
                        faces[Int(L.ball[face+1].neigh)+1] = face
                    end
                    face = cur
                end
                last = length(L.ball) - 1
                if face != last
                    L.ball[face+1] = L.ball[last+1]
                    faces[Int(L.ball[face+1].neigh)+1] = face
                end
                pop!(L.ball)
            else
                out = faces[4Int(in÷4)+f+1]
                faces[Int(out)+1] = length(L.ball)
                bnd_push_facet!(L, m, t, f, out)
            end
        end
    end
    for i in 1:length(L.ball)
        in = L.ball[i].neigh
        set_ball_neigh!(L, i, tneigh(m, 4tets[Int(in÷4)+1] + in % 4))
    end
    shift = 0
    for i in 1:ntet
        if tets[i] == NO_ADJACENT
            shift += 1
        elseif shift != 0
            tets[i-shift] = tets[i]
        end
    end
    resize!(L.deleted, length(L.deleted) - shift)
    return nothing
end

function restore_deleted!(m::Mesh, L::Local, previous::Int)
    for i in previous+1:length(L.deleted)
        unsetdeleted!(m, L.deleted[i])
    end
    resize!(L.deleted, previous)
end

function mean_nodal_size!(L::Local, ns::NodalSizes, vta::UInt32)
    s = 0.0; denom = 0.0
    for facet in L.ball, j in 1:3
        node = facet.node[j]
        if node != GHOST && ns.array[Int(node)+1] > 0.0
            s += ns.array[Int(node)+1]
            denom += 1.0
        end
    end
    denom > 0.0 && (ns.array[Int(vta)+1] = s / denom)
    return nothing
end

# createNewDeleted: append fresh deleted cells (a block of DELETED_BUFFER_SIZE)
function create_new_deleted!(m::Mesh, L::Local, demand::Integer)
    needed = max(DELETED_BUFFER_SIZE, demand) - length(L.deleted)
    first = m.ntet
    m.ntet += needed
    m.ntet > length(m.flag) && reserve!(m, m.ntet + length(m.flag))
    for i in 0:needed-1
        push!(L.deleted, UInt64(first + i))
        settflag!(m, first + i, 0)
        setdeleted!(m, first + i)
        m.color === nothing || (m.color[first+i+1] = 0)
    end
    return nothing
end

# fillingACavity: the last `length(ball)` deleted cells become the star of vta
function fill_cavity!(m::Mesh, L::Local, cur::Base.RefValue{UInt64}, vta::UInt32, color::UInt32)
    blength = length(L.ball)
    start = length(L.deleted) - blength
    for i in 1:blength
        t = L.deleted[start+i]
        f = L.ball[i]
        m.color === nothing || (m.color[Int(t)+1] = color)
        base = 4t
        settnode!(m, base, vta)
        settnode!(m, base + 1, f.node[1]); settnode!(m, base + 2, f.node[2]); settnode!(m, base + 3, f.node[3])
        settneigh!(m, base, f.neigh)
        settflag!(m, t, f.flag)
        settneigh!(m, f.neigh, base)
    end
    # Adjacency inside the star: the facet (vta, x, y) of one new cell is
    # shared with exactly one other; pair them by their edge (x, y).
    # (Each edge appears exactly twice, so the pairing is unique however it
    # is found; HXT uses a small map or a hash table, as here.)
    if blength <= 58
        fill_adjacencies_fast!(m, L, start, blength)
    else
        fill_adjacencies_slow!(m, L, start, blength)
    end
    cur[] = L.deleted[start+1]
    resize!(L.deleted, start)
    return nothing
end

# computeAdjacenciesFast: at most 58 facets have at most 31 distinct vertices
# (the ghost is 0), so the directed edges index a 32x32 map
function fill_adjacencies_fast!(m::Mesh, L::Local, start::Int, blength::Int)
    id = L.vertex_id
    length(id) < m.nv && resize!(id, max(m.nv, 2length(id)))
    nodes = L.local_nodes
    resize!(nodes, blength)
    @inbounds begin
        for i in 1:blength
            f = L.ball[i].node
            id[Int(f[1])+1] = 0xff; id[Int(f[2])+1] = 0xff
            f[3] == GHOST || (id[Int(f[3])+1] = 0xff)
        end
        npts = 0x01
        for i in 1:blength
            f = L.ball[i].node
            a = id[Int(f[1])+1]; a == 0xff && (a = npts; id[Int(f[1])+1] = npts; npts += 0x01)
            b = id[Int(f[2])+1]; b == 0xff && (b = npts; id[Int(f[2])+1] = npts; npts += 0x01)
            c = 0x00
            if f[3] != GHOST
                c = id[Int(f[3])+1]; c == 0xff && (c = npts; id[Int(f[3])+1] = npts; npts += 0x01)
            end
            nodes[i] = (a, b, c)
        end
        map = L.edge_map
        for i in 1:blength
            base = 4L.deleted[start+i]
            a, b, c = nodes[i]
            map[32Int(a)+Int(b)+1] = base + 3
            map[32Int(b)+Int(c)+1] = base + 1
            map[32Int(c)+Int(a)+1] = base + 2
        end
        for i in 1:blength
            base = 4L.deleted[start+i]
            a, b, c = nodes[i]
            settneigh!(m, base + 1, map[32Int(c)+Int(b)+1])
            settneigh!(m, base + 2, map[32Int(a)+Int(c)+1])
            settneigh!(m, base + 3, map[32Int(b)+Int(a)+1])
        end
    end
    return nothing
end

function fill_adjacencies_slow!(m::Mesh, L::Local, start::Int, blength::Int)
    capacity = length(L.edge_key)
    if capacity < 8blength
        capacity = nextpow(2, 8blength)
        L.edge_key = zeros(UInt64, capacity); L.edge_facet = zeros(UInt64, capacity)
        L.edge_epoch = zeros(UInt32, capacity); L.epoch = 0
    elseif L.epoch == typemax(UInt32)
        fill!(L.edge_epoch, 0); L.epoch = 0
    end
    L.epoch += UInt32(1)
    epoch = L.epoch; mask = capacity - 1
    open = 0
    @inbounds for i in 1:blength
        t = L.deleted[start+i]; base = 4t
        n = (tnode(m, base + 1), tnode(m, base + 2), tnode(m, base + 3))
        # facet j+1 is opposite node j+1; its two non-vta nodes:
        for (slot, e) in ((1, (n[2], n[3])), (2, (n[1], n[3])), (3, (n[1], n[2])))
            key = e[1] < e[2] ? (UInt64(e[1]) << 32) | e[2] : (UInt64(e[2]) << 32) | e[1]
            h = Int(hash64(key) & mask)
            while L.edge_epoch[h+1] == epoch && L.edge_key[h+1] != key
                h = (h + 1) & mask
            end
            if L.edge_epoch[h+1] == epoch && L.edge_facet[h+1] != NO_ADJACENT
                other = L.edge_facet[h+1]
                settneigh!(m, base + slot, other)
                settneigh!(m, other, base + slot)
                L.edge_facet[h+1] = NO_ADJACENT
                open -= 1
            else
                L.edge_epoch[h+1] = epoch; L.edge_key[h+1] = key; L.edge_facet[h+1] = base + slot
                open += 1
            end
        end
    end
    open == 0 || error("HXT: cavity star is not closed")
    return nothing
end

@inline square_dist(a, b) = (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]) +
    (a[3] - b[3]) * (a[3] - b[3])

@inline function too_close(s0, s1, d2, ns::NodalSizes)
    (s0 > 0.0 && s1 > 0.0) || return false
    mean = min(ns.max, max(ns.min, 0.5 * (s0 + s1))) * ns.factor
    return d2 < mean * mean
end

function insertion!(m::Mesh, L::Local, opt::DelaunayOptions, cur::Base.RefValue{UInt64}, vta::UInt32)
    previous = length(L.deleted)
    status = walk!(m, L, cur, vta)
    status == ST_TRUE || return status
    color = m.color === nothing ? COLOR_OUT : m.color[Int(cur[])+1]
    m.color !== nothing && !opt.allow_outer_insertion && color == COLOR_OUT && return ST_FALSE
    ns = opt.nodal_sizes
    ns !== nothing && ns.enabled && filter_tet(m, ns, cur[], vta) && return ST_FALSE
    edge_constrained = dig!(m, L, cur[], vta, opt.perfect_delaunay)
    if !opt.perfect_delaunay
        undelete = edge_constrained && respect_edge_constraint!(m, L, vta, color, previous)
        reshape_cavity!(m, L, vta, previous, undelete)
    end
    if ns !== nothing
        ns.array[Int(vta)+1] <= 0.0 && mean_nodal_size!(L, ns, vta)
        if ns.enabled && filter_cavity(m, L, ns, vta)
            restore_deleted!(m, L, previous)
            return ST_FALSE
        end
    end
    length(L.ball) > length(L.deleted) && create_new_deleted!(m, L, length(L.ball))
    fill_cavity!(m, L, cur, vta, color)
    return ST_TRUE
end

function filter_tet(m::Mesh, ns::NodalSizes, t::UInt64, vta::UInt32)
    p = point(m, vta); s = ns.array[Int(vta)+1]
    for j in 0:3
        node = tnode(m, 4t + j)
        (j == 3 && node == GHOST) && continue
        other = ns.array[Int(node)+1]
        other <= 0.0 && (other = s)
        too_close(s, other, square_dist(p, point(m, node)), ns) && return true
    end
    return false
end

function filter_cavity(m::Mesh, L::Local, ns::NodalSizes, vta::UInt32)
    p = point(m, vta); s = ns.array[Int(vta)+1]
    for facet in L.ball, j in 1:3
        node = facet.node[j]
        node == GHOST && continue
        other = ns.array[Int(node)+1]
        other <= 0.0 && (other = s)
        too_close(s, other, square_dist(p, point(m, node)), ns) && return true
    end
    return false
end

# ── hxtRemoveDeleted ───────────────────────────────────────────────────────────
function remove_deleted!(m::Mesh)
    m.ntet == 0 && return m
    for t in 0:m.ntet-1
        deleted(m, t) || continue
        for j in 0:3
            index = 4t + j
            nb = tneigh(m, index)
            (nb != NO_ADJACENT && tneigh(m, nb) == index) && settneigh!(m, nb, NO_ADJACENT)
        end
    end
    right = m.ntet - 1; left = 0
    while true
        while left < right && deleted(m, right); right -= 1; end
        while left < right && !deleted(m, left); left += 1; end
        left >= right && break
        m.color === nothing || (m.color[left+1] = m.color[right+1])
        settflag!(m, left, tflag(m, right))
        for j in 0:3
            nr = tneigh(m, 4right + j)
            nr != NO_ADJACENT && settneigh!(m, nr, 4left + j)
            settnode!(m, 4left + j, tnode(m, 4right + j))
            settneigh!(m, 4left + j, nr)
        end
        left += 1; right -= 1
    end
    (left == right && !deleted(m, left)) && (left += 1)
    m.ntet = left
    return m
end

# ── parallelDelaunay3D, one thread ─────────────────────────────────────────────
function filter_on_moore!(m::Mesh, info::AbstractVector{NodeInfo}, ns::NodalSizes)
    skipped = 0
    p1 = nothing; p1size = 0.0
    for i in eachindex(info)
        x = info[i]
        x.status == ST_TRYAGAIN || continue
        p2 = point(m, x.node); p2size = ns.array[Int(x.node)+1]
        if p1 !== nothing && too_close(p1size, p2size, square_dist(p1, p2), ns)
            skipped += 1
            info[i] = NodeInfo(x.hilbert, x.node, ST_FALSE)
        else
            p1 = p2; p1size = p2size
        end
    end
    return skipped
end

function delaunay3d!(m::Mesh, opt::DelaunayOptions, info::Vector{NodeInfo}, n::Int, noreorder::Bool)
    passes, npasses = compute_passes(opt.num_vertices_in_mesh, n, opt.partitionability)
    first_pass_ever = opt.num_vertices_in_mesh == 0
    offset = m.nv - n                          # first vertex to insert
    ns = opt.nodal_sizes
    if ns === nothing && !noreorder
        (npasses > 1 || first_pass_ever) && shuffle_vertices!(m, offset, n)
        moore!(m, opt.bbox, m.nv)
        for i in Int(first_pass_ever):npasses-1
            sort_vertices!(m, offset + passes[i+1], passes[i+2] - passes[i+1])
        end
        for i in 1:n
            x = info[i]; info[i] = NodeInfo(pad(m, offset + i - 1), x.node, x.status)
        end
    else
        moore!(m, opt.bbox, m.nv)
        (npasses > 1 || first_pass_ever) && shuffle_info!(info)
        for i in 1:n
            x = info[i]; info[i] = NodeInfo(pad(m, x.node), x.node, x.status)
        end
        for i in Int(first_pass_ever):npasses-1
            sort_info!(view(info, passes[i+1]+1:passes[i+2]))
        end
        if !noreorder
            first = opt.insertion_first
            vcopy = Vector{Float64}(undef, 4n); scopy = Vector{Float64}(undef, n)
            for i in 1:n
                src = Int(info[i].node) - first
                for k in 1:4; vcopy[4(i-1)+k] = m.coord[4(offset+src)+k]; end
                scopy[i] = ns.array[first+src+1]
                info[i] = NodeInfo(info[i].hilbert, UInt32(first + i - 1), info[i].status)
            end
            m.coord[4offset+1:4(offset+n)] .= vcopy
            ns.array[first+1:first+n] .= scopy
        end
    end
    if m.ntet < 5
        tetrahedra_init!(m, info, n, o3d_static_filter(opt.bbox))
        opt.num_vertices_in_mesh = 4
        passes[1] = 4
    end
    L = Local(opt.bbox)
    total_inserted = 0; total_skipped = 0
    for ipass in 0:npasses-1
        while passes[ipass+2] - passes[ipass+1] != 0
            pstart = Int(passes[ipass+1]); pend = Int(passes[ipass+2]); plength = pend - pstart
            ns !== nothing && ns.enabled && filter_on_moore!(m, view(info, pstart+1:pend), ns)
            cur = Ref(NO_ADJACENT)
            for t in 0:m.ntet-1
                if !deleted(m, t)
                    cur[] = UInt64(t); break
                end
            end
            if cur[] != NO_ADJACENT
                for i in 0:plength-1
                    x = info[pstart+i+1]
                    if x.status == ST_TRYAGAIN
                        status = insertion!(m, L, opt, cur, x.node)
                        status == ST_DOUBLE_PT && (status = ST_FALSE)
                        status == ST_CONFLICT && (status = ST_TRYAGAIN)
                        info[pstart+i+1] = NodeInfo(x.hilbert, x.node, status)
                    else
                        info[pstart+i+1] = NodeInfo(x.hilbert, x.node, ST_FALSE)
                    end
                end
            end
            shift = 0; skipped = 0
            i = pend
            while i > pstart
                i -= 1
                if info[i+1].status != ST_TRYAGAIN
                    info[i+1].status == ST_FALSE && (skipped += 1)
                    shift += 1
                elseif shift != 0
                    info[i+1], info[i+shift+1] = info[i+shift+1], info[i+1]
                end
            end
            total_inserted += shift - skipped
            total_skipped += skipped
            passes[ipass+1] += UInt32(shift)
        end
    end
    opt.num_vertices_in_mesh += total_inserted
    for t in L.deleted, j in 0:3
        settneigh!(m, 4t + j, NO_ADJACENT)
    end
    remove_deleted!(m)
    if !noreorder && total_skipped != 0
        compact_skipped_vertices!(m, opt, info, n)
        m.nv -= total_skipped
    end
    return total_inserted, total_skipped
end

# "if reordering allowed, remove vertices we could not insert"
function compact_skipped_vertices!(m::Mesh, opt::DelaunayOptions, info::Vector{NodeInfo}, n::Int)
    first_shifted = m.nv - n
    count = n
    status = Dict{Int,Int32}()
    for i in 1:n
        setpad!(m, info[i].node, reinterpret(UInt64, Int64(info[i].status)) & 0x00000000ffffffff |
                (pad(m, info[i].node) & 0xffffffff00000000))
    end
    getstatus(v) = reinterpret(Int32, UInt32(pad(m, v) & 0xffffffff))
    i = 0
    while getstatus(first_shifted + i) == ST_TRUE; i += 1; end
    first_shifted += i + 1
    count -= i + 1
    start = first_shifted - 1
    for i in 0:count-1
        old = start
        getstatus(first_shifted + i) == ST_TRUE && (start += 1)
        setpad!(m, first_shifted + i, (pad(m, first_shifted + i) & 0xffffffff00000000) | UInt64(UInt32(old)))
    end
    index(v) = UInt32(pad(m, v) & 0xffffffff)
    for i in 0:4m.ntet-1
        v = tnode(m, i)
        (v >= first_shifted && v != GHOST) && settnode!(m, i, index(v))
    end
    ns = opt.nodal_sizes
    for v in first_shifted:m.nv-1
        to = Int(index(v))
        ns === nothing || (ns.array[to+1] = ns.array[v+1])
        for k in 1:4; m.coord[4to+k] = m.coord[4v+k]; end
    end
    return m
end

# ── nodal sizes (hxt_tetNodalSize.c) ───────────────────────────────────────────
function nodal_sizes_init(m::Mesh, min, max, factor)
    array = Vector{Float64}(undef, m.nv)
    big = floatmax(Float64)
    for i in 0:m.nv-1
        if m.coord[4i+4] <= 0.0
            m.coord[4i+4] = 0.0
            array[i+1] = 0.0
        else
            array[i+1] = m.coord[4i+4]
            m.coord[4i+4] = big
        end
    end
    function add!(n1, n2)
        X1 = point(m, n1); X2 = point(m, n2)
        (m.coord[4Int(n1)+4] == big && m.coord[4Int(n2)+4] == big) && return
        l = sqrt((X1[1] - X2[1]) * (X1[1] - X2[1]) + (X1[2] - X2[2]) * (X1[2] - X2[2]) +
                 (X1[3] - X2[3]) * (X1[3] - X2[3]))
        if m.coord[4Int(n1)+4] != big
            m.coord[4Int(n1)+4] += 1; array[Int(n1)+1] += l
        end
        if m.coord[4Int(n2)+4] != big
            m.coord[4Int(n2)+4] += 1; array[Int(n2)+1] += l
        end
    end
    for t in 0:length(m.tri_node)÷3-1, j in 0:1, k in j+1:2
        add!(m.tri_node[3t+j+1], m.tri_node[3t+k+1])
    end
    for t in 0:length(m.line_node)÷2-1
        add!(m.line_node[2t+1], m.line_node[2t+2])
    end
    for i in 0:m.nv-1
        m.coord[4i+4] == big && continue
        array[i+1] /= m.coord[4i+4] * factor
    end
    return NodalSizes(array, min, max, factor, false)
end

# ── boundary facets and lines in the tetrahedralization (hxt_tetFlag.c) ────────
@inline hash_less(a::UInt64, b::UInt64) = (((a ⊻ b) & 1) ⊻ UInt64(a < b)) != 0

@inline function sorted3(a::UInt32, b::UInt32, c::UInt32)
    a > b && ((a, b) = (b, a)); b > c && ((b, c) = (c, b)); a > b && ((a, b) = (b, a))
    return (a, b, c)
end

# hxtGetTri2TetMap: for each triangle the facet 4t+f that carries it (the side
# HXT registers), NO_ADJACENT when the triangle is missing.
function tri2tet_map(m::Mesh)
    # Each facet is registered from exactly one side (hash_less is a strict
    # order on distinct cells), so the map does not depend on the visit order.
    # The triangles are bucketed by their smallest vertex and each such facet
    # is looked up in its bucket.
    ntri = length(m.tri_node) ÷ 3
    first = zeros(Int, m.nv + 1)
    @inbounds for t in 0:ntri-1
        a, _, _ = sorted3(m.tri_node[3t+1], m.tri_node[3t+2], m.tri_node[3t+3])
        first[Int(a)+1] += 1
    end
    total = 0
    @inbounds for v in 1:m.nv
        c = first[v]; first[v] = total; total += c
    end
    first[m.nv+1] = total
    bucket = Vector{Int32}(undef, ntri)
    fill_at = copy(first)
    @inbounds for t in 0:ntri-1
        a, _, _ = sorted3(m.tri_node[3t+1], m.tri_node[3t+2], m.tri_node[3t+3])
        fill_at[Int(a)+1] += 1
        bucket[fill_at[Int(a)+1]] = Int32(t)
    end
    map = fill(NO_ADJACENT, ntri)
    @inbounds for i in 0:m.ntet-1
        tnode(m, 4i + 3) == GHOST && continue
        for j in 0:3
            nb = tneigh(m, 4i + j)
            if nb == NO_ADJACENT || hash_less(UInt64(i), nb ÷ 4) || tnode(m, nb) == GHOST
                a, b, c = sorted3(tnode(m, 4i + ((j + 1) & 3)), tnode(m, 4i + ((j + 2) & 3)),
                                  tnode(m, 4i + ((j + 3) & 3)))
                for k in first[Int(a)+1]+1:first[Int(a)+2]
                    t = Int(bucket[k])
                    x, y, z = sorted3(m.tri_node[3t+1], m.tri_node[3t+2], m.tri_node[3t+3])
                    (y == b && z == c) && (map[t+1] = UInt64(4i + j))
                end
            end
        end
    end
    missing = count(==(NO_ADJACENT), map)
    return map, missing
end

# hxtGetLines2TriMap: the lines that are not an edge of a triangle
function lines_not_in_triangles(m::Mesh)
    edges = Set{Tuple{UInt32,UInt32}}()
    for t in 0:length(m.tri_node)÷3-1
        a, b, c = sorted3(m.tri_node[3t+1], m.tri_node[3t+2], m.tri_node[3t+3])
        push!(edges, (a, b), (a, c), (b, c))
    end
    out = Int[]
    for l in 0:length(m.line_node)÷2-1
        a = m.line_node[2l+1]; b = m.line_node[2l+2]
        a == b && continue
        (min(a, b), max(a, b)) in edges || push!(out, l)
    end
    return out
end

function constrain_triangles!(m::Mesh, map::Vector{UInt64})
    face = falses(4m.ntet)
    for f in map
        f == NO_ADJACENT && continue
        face[Int(f)+1] = true
        face[Int(tneigh(m, f))+1] = true
    end
    for i in 0:m.ntet-1, j in 0:3
        face[4i+j+1] && settflag!(m, i, tflag(m, i) | (UInt16(0x100) << j))
    end
    return m
end

const _NODES_OF_EDGE = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

# hxtGetLines2TetMap + hxtConstrainLinesNotInTriangles: every cell around a
# constrained line gets the edge constraint (which cell registers it does not
# change the result). Returns the number of lines missing from the mesh.
function constrain_lines!(m::Mesh, lines::Vector{Int})
    isempty(lines) && return 0
    where = Dict{Tuple{UInt32,UInt32},Tuple{UInt64,Int}}()
    for t in 0:m.ntet-1, e in 0:5
        n0, n1 = _NODES_OF_EDGE[e+1]
        p0 = tnode(m, 4t + n0); p1 = tnode(m, 4t + n1)
        (p0 == GHOST || p1 == GHOST) && continue
        key = (min(p0, p1), max(p0, p1))
        haskey(where, key) || (where[key] = (UInt64(t), e))
    end
    missing = 0
    flags = falses(6m.ntet)
    for l in lines
        a = m.line_node[2l+1]; b = m.line_node[2l+2]
        found = get(where, (min(a, b), max(a, b)), nothing)
        if found === nothing
            missing += 1; continue
        end
        tet, edge = found
        in, out = facets_of_edge(edge)
        cur = tet
        while true
            flags[6Int(cur)+edge_of_facets(in, out)+1] = true
            cur, in, out, _ = turn(m, cur, in, out)
            cur == tet && break
        end
    end
    for t in 0:m.ntet-1, e in 0:5
        flags[6t+e+1] && settflag!(m, t, tflag(m, t) | (UInt16(1) << e))
    end
    return missing
end

# ── volumes (hxt_tetColor.c) ───────────────────────────────────────────────────
function color_mesh!(m::Mesh)
    colors = zeros(UInt32, max(length(m.flag), m.ntet))
    m.color = colors
    color = UInt32(1); color_out = UInt32(0)
    stack = UInt64[]
    while true
        first = findfirst(t -> colors[t] == 0, 1:m.ntet)
        first === nothing && break
        empty!(stack)
        push!(stack, UInt64(first - 1))
        colors[first] = color
        i = 1
        while i <= length(stack)
            t = stack[i]
            tnode(m, 4t + 3) == GHOST && (color_out = color)
            for j in 0:3
                nb = tneigh(m, 4t + j)
                if nb != NO_ADJACENT && facet_constraint(m, t, j) == 0
                    n = nb ÷ 4
                    if colors[Int(n)+1] == 0
                        push!(stack, n)
                        colors[Int(n)+1] = color
                    end
                end
            end
            i += 1
        end
        color += UInt32(1)
    end
    nvolumes = Int(color) - 2
    for i in 1:m.ntet
        c = colors[i]
        colors[i] = c == color_out ? COLOR_OUT : c > color_out ? c - UInt32(2) : c - UInt32(1)
    end
    return nvolumes
end

function map_colors_to_brep!(m::Mesh, map::Vector{UInt64})
    nvolumes = color_mesh!(m)
    nvolumes == 0 && throw(ArgumentError("HXT: no closed volume"))
    pairs = UInt64[]
    last0 = typemax(UInt64); last1 = typemax(UInt64)
    for i in eachindex(map)
        m.tri_color[i] == COLOR_OUT && continue
        tet = map[i] ÷ 4
        nb = tneigh(m, map[i])
        c1 = UInt64(m.color[Int(tet)+1])
        c2 = nb == NO_ADJACENT ? UInt64(COLOR_OUT) : UInt64(m.color[Int(nb÷4)+1])
        for c in (c1, c2)
            c == COLOR_OUT && continue
            pair = c << 32 | m.tri_color[i]
            if pair != last0 && pair != last1
                push!(pairs, pair)
                last1 = last0; last0 = pair
            end
        end
    end
    sort!(pairs); unique!(pairs)
    reconstructed = [UInt32[] for _ in 1:nvolumes]
    for p in pairs
        push!(reconstructed[Int(p >> 32)+1], UInt32(p & 0xffffffff))
    end
    isempty(m.surfaces_per_volume) && return nvolumes
    original = [sort(s) for s in m.surfaces_per_volume]
    length(original) > nvolumes && throw(ArgumentError("HXT: BRep contains more volumes than there really are"))
    volume_map = fill(COLOR_OUT, nvolumes)
    for (c, s) in enumerate(reconstructed)
        k = findall(==(s), original)
        length(k) > 1 && throw(ArgumentError("HXT: duplicate volume in the BRep"))
        length(findall(==(s), reconstructed)) > 1 && throw(ArgumentError("HXT: duplicate reconstructed volume"))
        isempty(k) || (volume_map[c] = UInt32(k[1] - 1))
    end
    count(!=(COLOR_OUT), volume_map) == length(original) ||
        throw(ArgumentError("HXT: some volumes of the BRep were not found"))
    next = UInt32(length(original))
    for c in eachindex(volume_map)
        volume_map[c] == COLOR_OUT && (volume_map[c] = next; next += UInt32(1))
    end
    for i in 1:m.ntet
        m.color[i] == COLOR_OUT || (m.color[i] = volume_map[Int(m.color[i])+1])
    end
    return nvolumes
end

const PROCESSED_MASK = 0x0080
@inline processed(m::Mesh, t) = (tflag(m, t) & PROCESSED_MASK) != 0
@inline setprocessed!(m::Mesh, t) = settflag!(m, t, tflag(m, t) | PROCESSED_MASK)

function process_only_brep_volumes!(m::Mesh)
    nvol = length(m.surfaces_per_volume)
    for t in 0:m.ntet-1
        f = tflag(m, t)
        settflag!(m, t, m.color[t+1] >= nvol ? f | PROCESSED_MASK : f & ~PROCESSED_MASK)
    end
    return m
end

# ── refinement (hxt_tetRefine.c) ───────────────────────────────────────────────
# C's sqrt: NaN for a negative argument, which the comparisons below reject.
@inline csqrt(x::Float64) = Base.Math.sqrt_llvm(x)

# getBestCenter: the circumcentre in the metric of the nodal sizes, or a point
# pulled towards the incentre when that is outside or too close to a corner.
function best_center(p::NTuple{4,NTuple{3,Float64}}, sizes::NTuple{4,Float64}, ns::NodalSizes)
    avg = 0.0; num = 0.0
    for s in sizes
        s > 0.0 && (avg += s; num += 1.0)
    end
    avg = num == 0 ? 1.0 : avg / num
    s0 = sizes[1] > 0.0 ? sizes[1] : avg
    s1 = sizes[2] > 0.0 ? sizes[2] : avg
    s2 = sizes[3] > 0.0 ? sizes[3] : avg
    s3 = sizes[4] > 0.0 ? sizes[4] : avg
    e0l2 = square_dist(p[1], p[2]) / (0.25 * (s0 + s1) * (s0 + s1))
    e1l2 = square_dist(p[1], p[3]) / (0.25 * (s0 + s2) * (s0 + s2))
    e2l2 = square_dist(p[1], p[4]) / (0.25 * (s0 + s3) * (s0 + s3))
    e3l2 = square_dist(p[2], p[3]) / (0.25 * (s1 + s2) * (s1 + s2))
    e4l2 = square_dist(p[2], p[4]) / (0.25 * (s1 + s3) * (s1 + s3))
    e5l2 = square_dist(p[3], p[4]) / (0.25 * (s2 + s3) * (s2 + s3))
    xa = csqrt(e0l2)
    invtwoxa = 1.0 / (2 * xa)
    xb = (e1l2 + e0l2 - e3l2) * invtwoxa
    yb = csqrt(e1l2 - xb * xb)
    xc = (e2l2 + e0l2 - e4l2) * invtwoxa
    yc = (e1l2 + e2l2 - e5l2 - 2 * xb * xc) / (2 * yb)
    zc = csqrt(e2l2 - xc * xc - yc * yc)
    xcrossbc = yb * zc
    ycrossbc = -zc * xb
    zcrossbc = xb * yc - xc * yb
    ycrossca = zc * xa
    zcrossca = -xa * yc
    zcrossab = xa * yb
    denominator = 0.5 / (xa * yb * zc)
    xcirca = 0.5 * xa
    ycirca = (e0l2 * ycrossbc + e1l2 * ycrossca) * denominator
    zcirca = (e0l2 * zcrossbc + e1l2 * zcrossca + e2l2 * zcrossab) * denominator
    bary0 = (xcirca * xcrossbc + ycirca * ycrossbc + zcirca * zcrossbc) * (2.0 * denominator)
    bary1 = (ycirca * ycrossca + zcirca * zcrossca) * (2.0 * denominator)
    bary2 = (zcirca * zcrossab) * (2.0 * denominator)
    bary3 = 1.0 - bary0 - bary1 - bary2
    combine(b0, b1, b2, b3) = (b0 * p[1][1] + b1 * p[2][1] + b2 * p[3][1] + b3 * p[4][1],
        b0 * p[1][2] + b1 * p[2][2] + b2 * p[3][2] + b3 * p[4][2],
        b0 * p[1][3] + b1 * p[2][3] + b2 * p[3][3] + b3 * p[4][3],
        b0 * s0 + b1 * s1 + b2 * s2 + b3 * s3)
    close(c) = too_close(s0, c[4], square_dist(p[1], c), ns) || too_close(s1, c[4], square_dist(p[2], c), ns) ||
        too_close(s2, c[4], square_dist(p[3], c), ns) || too_close(s3, c[4], square_dist(p[4], c), ns)
    outside = false; circum_close = false
    center = (0.0, 0.0, 0.0, 0.0)
    if bary0 > 0.0 && bary1 > 0.0 && bary2 > 0.0 && bary3 > 0.0
        center = combine(bary0, bary1, bary2, bary3)
        circum_close = close(center)
    else
        outside = true
    end
    other_close = false
    if outside || circum_close
        xsumcros = xcrossbc
        ysumcros = ycrossbc + ycrossca
        zsumcros = zcrossab + zcrossbc + zcrossca
        invA0x2 = 1.0 / csqrt(xsumcros * xsumcros + ysumcros * ysumcros + zsumcros * zsumcros)
        invA1x2 = 1.0 / csqrt(xcrossbc * xcrossbc + ycrossbc * ycrossbc + zcrossbc * zcrossbc)
        invA2x2 = 1.0 / csqrt(ycrossca * ycrossca + zcrossca * zcrossca)
        invA3x2 = 1.0 / zcrossab
        den = invA0x2 + invA1x2 + invA2x2 + invA3x2
        alpha = 0.5
        b0 = (1.0 - alpha) * 0.25 + alpha * invA0x2 / den
        b1 = (1.0 - alpha) * 0.25 + alpha * invA1x2 / den
        b2 = (1.0 - alpha) * 0.25 + alpha * invA2x2 / den
        b3 = (1.0 - alpha) * 0.25 + alpha * invA3x2 / den
        if !isfinite(b0) || !isfinite(b1) || !isfinite(b2) || !isfinite(b3)
            b0 = b1 = b2 = b3 = 0.25
        end
        other = combine(b0, b1, b2, b3)
        other_close = close(other)
        (outside || !other_close) && (center = other)
    end
    return center, other_close
end

"""
    refine!(m, opt, size_callback)

`hxtRefineTetrahedra`: put a point in every unprocessed cell, size it with
`size_callback(x, y, z, interpolated, color)` (Gmsh's `nodalSizesCallBack`),
drop those too close to the cell's corners and insert the others, until no
point is added (42 rounds at most).
"""
function refine!(m::Mesh, opt::DelaunayOptions, size_callback)
    ns = opt.nodal_sizes::NodalSizes
    newv = NTuple{4,Float64}[]
    pt2tet = UInt64[]
    for iter in 0:41
        opt.insertion_first = m.nv
        empty!(newv); empty!(pt2tet)
        for t in 0:m.ntet-1
            processed(m, t) && continue
            setprocessed!(m, t)
            nodes = ntuple(i -> tnode(m, 4t + i - 1), 4)
            p = ntuple(i -> point(m, nodes[i]), 4)
            s = ntuple(i -> ns.array[Int(nodes[i])+1], 4)
            center, close = best_center(p, s, ns)
            size_callback === nothing && close && (center = (center[1], center[2], center[3], -floatmax(Float64)))
            push!(newv, center); push!(pt2tet, UInt64(t))
        end
        if size_callback !== nothing
            for i in eachindex(newv)
                c = newv[i]
                newv[i] = (c[1], c[2], c[3], size_callback(c[1], c[2], c[3], c[4], m.color[Int(pt2tet[i])+1]))
            end
        end
        kept = 0
        for i in eachindex(newv)
            c = newv[i]
            c[4] <= 0.0 && continue
            t = pt2tet[i]
            nodes = ntuple(k -> tnode(m, 4t + k - 1), 4)
            if any(k -> too_close(ns.array[Int(nodes[k])+1], c[4], square_dist(point(m, nodes[k]), c), ns), 1:4)
                newv[i] = (c[1], c[2], c[3], -floatmax(Float64))
            else
                kept += 1
            end
        end
        first = m.nv
        resize!(m.coord, 4(m.nv + kept)); resize!(ns.array, m.nv + kept)
        v = first
        for c in newv
            c[4] <= 0.0 && continue
            m.coord[4v+1] = c[1]; m.coord[4v+2] = c[2]; m.coord[4v+3] = c[3]; m.coord[4v+4] = 0.0
            ns.array[v+1] = c[4]
            v += 1
        end
        m.nv += kept
        opt.insertion_first == m.nv && break
        opt.partitionability = 1.0 - 0.5^iter
        before = opt.num_vertices_in_mesh
        n = m.nv - opt.insertion_first
        info = NodeInfo[NodeInfo(0, UInt32(opt.insertion_first + i), ST_TRYAGAIN) for i in 0:n-1]
        delaunay3d!(m, opt, info, n, false)
        opt.num_vertices_in_mesh == before && break
    end
    return m
end

# ── hxtEmptyMesh ───────────────────────────────────────────────────────────────
function empty_mesh!(m::Mesh, opt::DelaunayOptions)
    m.ntet == 0 || throw(ArgumentError("HXT: the input mesh should only contain triangles"))
    isempty(m.tri_node) && throw(ArgumentError("HXT: the input mesh should contain triangles"))
    for i in 0:m.nv-1; m.coord[4i+4] = 0.0; end
    for v in m.tri_node; m.coord[4Int(v)+4] = 1.0; end
    for v in m.line_node; m.coord[4Int(v)+4] = 1.0; end
    for v in m.point_node; m.coord[4Int(v)+4] = 1.0; end
    info = NodeInfo[]
    for i in 0:m.nv-1
        m.coord[4i+4] == 1.0 && push!(info, NodeInfo(0, UInt32(i), ST_TRYAGAIN))
    end
    opt.perfect_delaunay = true
    delaunay3d!(m, opt, info, length(info), true)
    opt.perfect_delaunay = false
    any(x -> x.status != ST_TRUE, info) &&
        @warn "HXT: a point of the empty mesh was not inserted"
    return m
end

include("HxtOpti.jl")

"""
    tet_mesh!(m; min_size, max_size, size_callback, quality_min=0.3,
              optimize=true, refine=true) -> m

`hxtTetMesh` on one thread: empty mesh of the boundary vertices, boundary
facets and lines constrained, volumes coloured after the BRep, refinement with
`size_callback(x, y, z, interpolated, color)`, then optimisation of every cell
below `quality_min`. Cells of volume `k` (0-based, the order of
`m.surfaces_per_volume`) have colour `k`.
"""
function tet_mesh!(m::Mesh; min_size::Float64=0.0, max_size::Float64=floatmax(Float64),
                   factor::Float64=1.0, size_callback=nothing, quality_min::Float64=0.3,
                   optimize::Bool=true, refine::Bool=true, boundary_recovery=nothing)
    max_size <= 0.0 && (max_size = floatmax(Float64))
    factor <= 0.0 && (factor = 1.0)
    b = bbox(m)
    constrained = m.nv
    ns = refine ? nodal_sizes_init(m, min_size, max_size, factor) : nothing
    opt = DelaunayOptions(b, ns, 0, 0, 0.0, false, false)
    empty_mesh!(m, opt)
    map, missing = tri2tet_map(m)
    lines = lines_not_in_triangles(m)
    missing_lines = missing == 0 ? constrain_lines!(m, lines) : 0
    if missing != 0 || missing_lines != 0
        boundary_recovery === nothing && throw(ArgumentError(
            "HXT: $missing boundary facets and $missing_lines lines are missing from the " *
            "Delaunay tetrahedralization and no boundary recovery was given"))
        # hxt_boundary_recovery: TetGen/BR rebuilds the cells and their
        # adjacency (flags cleared); the map is taken again and must be whole.
        nv = m.nv
        node, neigh, ntet = boundary_recovery(m)
        reserve!(m, ntet)
        copyto!(m.node, 1, node, 1, 4ntet); copyto!(m.neigh, 1, neigh, 1, 4ntet)
        fill!(view(m.flag, 1:ntet), 0x0000)
        m.ntet = ntet
        opt.num_vertices_in_mesh += m.nv - nv
        map, missing = tri2tet_map(m)
        missing == 0 || throw(ErrorException(
            "HXT: $missing boundary faces still missing (after recovery step)"))
        missing_lines = constrain_lines!(m, lines)
        missing_lines == 0 || throw(ErrorException(
            "HXT: $missing_lines constrained edges still missing (after recovery step)"))
    end
    constrain_triangles!(m, map)
    map_colors_to_brep!(m, map)
    if refine
        process_only_brep_volumes!(m)
        ns.enabled = true
        refine!(m, opt, size_callback)
    end
    if optimize
        process_only_brep_volumes!(m)
        optimize!(m, b, quality_min, constrained)
    end
    return m
end

end # module
