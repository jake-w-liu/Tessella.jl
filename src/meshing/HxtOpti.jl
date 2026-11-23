# SPDX-FileCopyrightText: 2016-2020 UCLouvain
# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: GPL-2.0-or-later
# Contains HXT adaptations from Gmsh 4.15.2 contrib/hxt/tetMesh/src (hxt_tetOpti.c,
# hxt_edgeRemoval.c, hxt_smoothing.c, HXTSPR.c, hxt_tetQuality.c). Upstream author:
# Célestin Marot. Julia adaptation/changes: Tessella contributors, 2026.
# See LICENSES/README.md and NOTICES.md#hxt. Included by Hxt.jl.

# ── quality (hxt_tetQuality.c) ─────────────────────────────────────────────────
function aspect_fast_ratio(a, b, c, d)
    ab = (b[1] - a[1], b[2] - a[2], b[3] - a[3])
    ac = (c[1] - a[1], c[2] - a[2], c[3] - a[3])
    ad = (d[1] - a[1], d[2] - a[2], d[3] - a[3])
    acxad0 = ac[2] * ad[3] - ac[3] * ad[2]
    adxab0 = ad[2] * ab[3] - ad[3] * ab[2]
    abxac0 = ab[2] * ac[3] - ab[3] * ac[2]
    volume6 = ab[1] * acxad0 + ac[1] * adxab0 + ad[1] * abxac0
    volume6 <= 0.0 && return 0.0
    acxad1 = ac[3] * ad[1] - ac[1] * ad[3]
    acxad2 = ac[1] * ad[2] - ac[2] * ad[1]
    adxab1 = ad[3] * ab[1] - ad[1] * ab[3]
    adxab2 = ad[1] * ab[2] - ad[2] * ab[1]
    abxac1 = ab[3] * ac[1] - ab[1] * ac[3]
    abxac2 = ab[1] * ac[2] - ab[2] * ac[1]
    db = (b[1] - d[1], b[2] - d[2], b[3] - d[3])
    bc = (c[1] - b[1], c[2] - b[2], c[3] - b[3])
    cd = (d[1] - c[1], d[2] - c[2], d[3] - c[3])
    bcxcd0 = bc[2] * cd[3] - bc[3] * cd[2]
    bcxcd1 = bc[3] * cd[1] - bc[1] * cd[3]
    bcxcd2 = bc[1] * cd[2] - bc[2] * cd[1]
    area = sqrt(acxad0 * acxad0 + acxad1 * acxad1 + acxad2 * acxad2) +
        sqrt(adxab0 * adxab0 + adxab1 * adxab1 + adxab2 * adxab2) +
        sqrt(abxac0 * abxac0 + abxac1 * abxac1 + abxac2 * abxac2) +
        sqrt(bcxcd0 * bcxcd0 + bcxcd1 * bcxcd1 + bcxcd2 * bcxcd2)
    l = ab[1] * ab[1] + ab[2] * ab[2] + ab[3] * ab[3]
    l2 = ac[1] * ac[1] + ac[2] * ac[2] + ac[3] * ac[3]
    l3 = ad[1] * ad[1] + ad[2] * ad[2] + ad[3] * ad[3]
    l4 = bc[1] * bc[1] + bc[2] * bc[2] + bc[3] * bc[3]
    l5 = cd[1] * cd[1] + cd[2] * cd[2] + cd[3] * cd[3]
    l6 = db[1] * db[1] + db[2] * db[2] + db[3] * db[3]
    l2 > l && (l = l2); l3 > l && (l = l3); l4 > l && (l = l4)
    l5 > l && (l = l5); l6 > l && (l = l6)
    return volume6 * volume6 / (l * area * area)
end

@inline aspect_ratio(a, b, c, d) = sqrt(24.0 * aspect_fast_ratio(a, b, c, d))
@inline tet_quality(m::Mesh, p0, p1, p2, p3) =
    aspect_ratio(point(m, p0), point(m, p1), point(m, p2), point(m, p3))

# ── the optimizer's state (ThreadShared + ThreadLocal, one thread) ─────────────
mutable struct Opti
    m::Mesh
    quality::Vector{Float64}
    creation::Vector{UInt16}
    check::Vector{UInt16}
    date::UInt16
    threshold::Float64
    deleted::Vector{UInt64}
    constrained::UInt32
    last_check::UInt64
    last_creation::UInt64
    max_search_nodes::Int
    # scratch reused by every edge removal, smoothing and SPR attempt
    er_tets::Vector{UInt64}
    er_up::Vector{UInt64}
    er_down::Vector{UInt64}
    er_flag::Vector{UInt16}
    er_annulus::Vector{UInt32}
    qup::Vector{Float64}
    qdown::Vector{Float64}
    golden::Vector{Float64}
    golden_f::Vector{Float64}
    spr::Any
end

const ST_INTERNAL = Int32(-1024)
const ST_CONSTRAINT = Int32(-1028)
const ST_NOTBETTER = Int32(-1029)
const CONFLICT_MASK = 0x1000

@inline conflict_flag(m::Mesh, t) = (tflag(m, t) & CONFLICT_MASK) != 0
@inline set_conflict!(m::Mesh, t) = settflag!(m, t, tflag(m, t) | CONFLICT_MASK)
@inline unset_conflict!(m::Mesh, t) = settflag!(m, t, tflag(m, t) & ~CONFLICT_MASK)
@inline reset_unused!(m::Mesh, t) = settflag!(m, t, tflag(m, t) & ~0xf000)

function grow!(o::Opti, n::Integer)
    if n > length(o.quality)
        resize!(o.quality, n); resize!(o.creation, n); resize!(o.check, n)
    end
end

# createNewDeleted for the optimizer (its quality and dates arrays follow)
function new_deleted!(o::Opti, demand::Integer)
    m = o.m
    needed = max(DELETED_BUFFER_SIZE, demand) - length(o.deleted)
    first = m.ntet
    m.ntet += needed
    m.ntet > length(m.flag) && reserve!(m, m.ntet)
    grow!(o, m.ntet)
    for i in 0:needed-1
        push!(o.deleted, UInt64(first + i))
        settflag!(m, first + i, 0)
        setdeleted!(m, first + i)
    end
end

const _UP_FACET = ((-1, 2, 3, 1), (3, -1, 0, 2), (1, 3, -1, 0), (2, 0, 1, -1))
@inline up_facet(i, o) = _UP_FACET[i+1][o+1]
@inline down_facet(i, o) = _UP_FACET[o+1][i+1]

include("HxtOptiTables.jl")

# ── edge removal (hxt_edgeRemoval.c) ───────────────────────────────────────────
function edge_removal!(o::Opti, bad::UInt64, edge::Int)
    m = o.m
    in, out = facets_of_edge(edge)
    edge_constraint(m, bad, edge_of_facets(in, out)) != 0 && return ST_CONSTRAINT
    v_up = tnode(m, 4bad + down_facet(in, out))
    v_down = tnode(m, 4bad + up_facet(in, out))
    tets = o.er_tets; neigh_up = o.er_up; neigh_down = o.er_down; flag = o.er_flag; annulus = o.er_annulus
    empty!(tets); empty!(neigh_up); empty!(neigh_down); empty!(flag); empty!(annulus)
    cur = bad
    while true
        oldv = tnode(m, 4cur + out)
        newv = tnode(m, 4cur + in)
        push!(tets, cur)
        uf = up_facet(in, out); df = down_facet(in, out)
        push!(neigh_up, tneigh(m, 4cur + uf)); push!(neigh_down, tneigh(m, 4cur + df))
        ud = edge_of_facets(uf, df); uo = edge_of_facets(uf, out); ui = edge_of_facets(uf, in)
        dout = edge_of_facets(df, out); din = edge_of_facets(df, in)
        b(x) = UInt16(x != 0)
        push!(flag, b(facet_constraint(m, cur, uf)) + (b(edge_constraint(m, cur, uo)) << 1) +
            (b(edge_constraint(m, cur, ud)) << 2) + (b(edge_constraint(m, cur, ui)) << 3) +
            (b(facet_constraint(m, cur, df)) << 4) + (b(edge_constraint(m, cur, dout)) << 5) +
            (b(edge_constraint(m, cur, din)) << 6) + (b(edge_constraint(m, cur, ud)) << 7))
        push!(annulus, oldv)
        nb = tneigh(m, 4cur + out)
        (nb == NO_ADJACENT || facet_constraint(m, cur, out) != 0) && return ST_CONSTRAINT
        length(tets) >= 7 && return ST_INTERNAL
        cur = nb ÷ 4
        in = Int(nb % 4)
        out = 0
        while out < 3 && tnode(m, 4cur + out) != newv
            out += 1
        end
        cur == bad && break
    end
    n = length(tets)
    worst = floatmax(Float64)
    for t in tets
        q = o.quality[Int(t)+1]
        q < worst && (worst = q)
    end
    pat = _ER_PATTERNS[n]
    ntri = n - 2
    qup = o.qup; qdown = o.qdown
    mask = UInt64(0)
    for i in eachindex(pat.triangles)
        tr = pat.triangles[i]
        p0 = annulus[tr[1]+1]; p1 = annulus[tr[2]+1]; p2 = annulus[tr[3]+1]
        qup[i] = tet_quality(m, p0, p1, v_up, p2)
        qdown[i] = tet_quality(m, p0, p1, p2, v_down)
        (qup[i] <= worst || qdown[i] <= worst) && (mask |= pat.triangle_in_triangul[i])
    end
    best = -1; best_worst = 0.0
    for i in 0:length(pat.triangulations)-1
        mask & (UInt64(1) << i) == 0 || continue
        cur_worst = floatmax(Float64)
        for j in 1:ntri
            tri = pat.triangulations[i+1][j] + 1
            qu = qup[tri]; qd = qdown[tri]
            if qu < best_worst || qd < best_worst
                cur_worst = 0.0; break
            end
            qu < cur_worst && (cur_worst = qu)
            qd < cur_worst && (cur_worst = qd)
        end
        if cur_worst > best_worst
            best_worst = cur_worst; best = i
        end
    end
    best == -1 && return ST_NOTBETTER
    for t in tets
        setdeleted!(m, t)
        push!(o.deleted, t)
    end
    2ntri > length(o.deleted) && new_deleted!(o, 2ntri)
    start = length(o.deleted) - 2ntri
    color = m.color[Int(bad)+1]
    D(k) = o.deleted[start+k+1]
    for i in 0:ntri-1
        tri = pat.triangulations[best+1][i+1]
        tr = pat.triangles[tri+1]
        p0 = annulus[tr[1]+1]; p1 = annulus[tr[2]+1]; p2 = annulus[tr[3]+1]
        nn = pat.neigh[best+1]
        n0 = nn[4i+1]; n1 = nn[4i+2]; n2 = nn[4i+3]
        up = D(2i); down = D(2i + 1)
        # upper cell
        m.color[Int(up)+1] = color; settflag!(m, up, 0)
        o.quality[Int(up)+1] = qup[tri+1]; o.creation[Int(up)+1] = o.date
        settnode!(m, 4up, p0); settnode!(m, 4up + 1, p1); settnode!(m, 4up + 2, v_up); settnode!(m, 4up + 3, p2)
        newtet(nv) = 4D((nv ÷ 4) * 2) + (nv % 4) + (nv % 4) ÷ 2
        if n0 >= 0
            settneigh!(m, 4up, newtet(n0))
        else
            f = flag[-n0]
            settneigh!(m, 4up, neigh_up[-n0])
            settflag!(m, up, tflag(m, up) | (f & 1) << 8 | (f & 2) >> 1 | (f & 4) >> 1 | (f & 8) >> 1)
            neigh_up[-n0] != NO_ADJACENT && settneigh!(m, neigh_up[-n0], 4up)
        end
        if n1 >= 0
            settneigh!(m, 4up + 1, newtet(n1))
        else
            f = flag[-n1]
            settneigh!(m, 4up + 1, neigh_up[-n1])
            settflag!(m, up, tflag(m, up) | (f & 1) << 9 | (f & 2) << 3 | (f & 4) << 1 | (f & 8) >> 3)
            neigh_up[-n1] != NO_ADJACENT && settneigh!(m, neigh_up[-n1], 4up + 1)
        end
        settneigh!(m, 4up + 2, 4down + 3)
        if n2 >= 0
            settneigh!(m, 4up + 3, newtet(n2))
        else
            f = flag[-n2]
            settneigh!(m, 4up + 3, neigh_up[-n2])
            settflag!(m, up, tflag(m, up) | (f & 1) << 11 | (f & 2) << 1 | (f & 4) << 3 | (f & 8) << 1)
            neigh_up[-n2] != NO_ADJACENT && settneigh!(m, neigh_up[-n2], 4up + 3)
        end
        # lower cell
        m.color[Int(down)+1] = color; settflag!(m, down, 0)
        o.quality[Int(down)+1] = qdown[tri+1]; o.creation[Int(down)+1] = o.date
        settnode!(m, 4down, p0); settnode!(m, 4down + 1, p1); settnode!(m, 4down + 2, p2); settnode!(m, 4down + 3, v_down)
        newdown(nv) = 4D((nv ÷ 4) * 2 + 1) + (nv % 4)
        if n0 >= 0
            settneigh!(m, 4down, newdown(n0))
        else
            f = flag[-n0]
            settneigh!(m, 4down, neigh_down[-n0])
            settflag!(m, down, tflag(m, down) | (f & 16) << 4 | (f & 32) >> 5 | (f & 64) >> 5 | (f & 128) >> 5)
            neigh_down[-n0] != NO_ADJACENT && settneigh!(m, neigh_down[-n0], 4down)
        end
        if n1 >= 0
            settneigh!(m, 4down + 1, newdown(n1))
        else
            f = flag[-n1]
            settneigh!(m, 4down + 1, neigh_down[-n1])
            settflag!(m, down, tflag(m, down) | (f & 16) << 5 | (f & 32) >> 2 | (f & 64) >> 6 | (f & 128) >> 3)
            neigh_down[-n1] != NO_ADJACENT && settneigh!(m, neigh_down[-n1], 4down + 1)
        end
        if n2 >= 0
            settneigh!(m, 4down + 2, newdown(n2))
        else
            f = flag[-n2]
            settneigh!(m, 4down + 2, neigh_down[-n2])
            settflag!(m, down, tflag(m, down) | (f & 16) << 6 | (f & 32) >> 4 | (f & 64) >> 3 | (f & 128) >> 2)
            neigh_down[-n2] != NO_ADJACENT && settneigh!(m, neigh_down[-n2], 4down + 2)
        end
        settneigh!(m, 4down + 3, 4up + 2)
    end
    resize!(o.deleted, start)
    return ST_TRUE
end

# ── smoothing (hxt_smoothing.c) ────────────────────────────────────────────────
function smoothing!(o::Opti, face::UInt64)
    m = o.m
    previous = length(o.deleted)
    vertex = tnode(m, face)
    status = ST_TRUE
    if vertex < o.constrained
        status = ST_CONSTRAINT
    else
        push!(o.deleted, face)
        setdeleted!(m, face ÷ 4)
        start = length(o.deleted)
        while start <= length(o.deleted) && status == ST_TRUE
            cf = o.deleted[start]
            ct = cf ÷ 4; fc = Int(cf % 4)
            for j in 1:3
                f = (fc + j) % 4
                nb = tneigh(m, 4ct + f)
                if nb == NO_ADJACENT
                    status = ST_CONSTRAINT; break
                end
                nt = nb ÷ 4; nf = Int(nb % 4)
                deleted(m, nt) && continue
                if facet_constraint(m, ct, f) != 0
                    status = ST_CONSTRAINT; break
                end
                f1 = up_facet(f, fc); f2 = down_facet(f, fc)
                if edge_constraint(m, ct, edge_of_facets(f, f1)) != 0 || edge_constraint(m, ct, edge_of_facets(f, f2)) != 0
                    status = ST_CONSTRAINT; break
                end
                for k in 1:3
                    if tnode(m, 4nt + (nf + k) % 4) == vertex
                        push!(o.deleted, 4nt + (nf + k) % 4)
                        break
                    end
                end
                setdeleted!(m, nt)
            end
            start += 1
        end
    end
    for i in previous+1:length(o.deleted)
        unsetdeleted!(m, o.deleted[i] ÷ 4)
    end
    if status != ST_TRUE
        resize!(o.deleted, previous)
        return status
    end
    ncav = length(o.deleted) - previous
    old = point(m, vertex)
    cx = 0.0; cy = 0.0; cz = 0.0
    third = 1.0 / 3.0
    for i in previous+1:length(o.deleted)
        t = o.deleted[i] ÷ 4; fc = Int(o.deleted[i] % 4)
        a = point(m, tnode(m, 4t + (fc + 1) % 4)); b = point(m, tnode(m, 4t + (fc + 2) % 4))
        c = point(m, tnode(m, 4t + (fc + 3) % 4))
        cx += (a[1] + b[1] + c[1]) * third
        cy += (a[2] + b[2] + c[2]) * third
        cz += (a[3] + b[3] + c[3]) * third
    end
    center = (cx / Float64(ncav), cy / Float64(ncav), cz / Float64(ncav))
    status = golden_section!(o, previous, old, center, vertex)
    resize!(o.deleted, previous)
    return status
end

function worst_quality(o::Opti, previous::Int)
    worst = floatmax(Float64)
    for i in previous+1:length(o.deleted)
        t = o.deleted[i] ÷ 4
        q = tet_quality(o.m, tnode(o.m, 4t), tnode(o.m, 4t + 1), tnode(o.m, 4t + 2), tnode(o.m, 4t + 3))
        q < worst && (worst = q)
    end
    return worst
end

@inline function interpolate!(m::Mesh, v, start, stop, alpha)
    for i in 1:3
        m.coord[4Int(v)+i] = alpha * start[i] + (1 - alpha) * stop[i]
    end
end

function golden_section!(o::Opti, previous::Int, start, stop, vertex)
    m = o.m
    tol = 0.00001
    invphi = 0.6180339887498949
    v = o.golden; fv = o.golden_f
    v[1] = 1.0; v[2] = 0.0; v[3] = invphi; v[4] = 1 - invphi
    for i in 1:4
        interpolate!(m, vertex, start, stop, v[i])
        fv[i] = worst_quality(o, previous)
    end
    original = fv[1]
    while fv[3] > fv[1] && fv[3] > fv[2] && fv[4] > fv[1] && fv[4] > fv[2] && abs(fv[1] - fv[2]) > tol
        if fv[3] >= fv[4]
            fv[2] = fv[4]; v[2] = v[4]
            fv[4] = fv[3]; v[4] = v[3]
            v[3] = v[1] + invphi * (v[2] - v[1])
            ind = 3
        else
            fv[1] = fv[3]; v[1] = v[3]
            fv[3] = fv[4]; v[3] = v[4]
            v[4] = v[2] + invphi * (v[1] - v[2])
            ind = 4
        end
        interpolate!(m, vertex, start, stop, v[ind])
        fv[ind] = worst_quality(o, previous)
    end
    fv[2] > fv[1] && (v[1] = v[2])
    if v[1] >= 0.99 || fv[1] <= original
        m.coord[4Int(vertex)+1] = start[1]; m.coord[4Int(vertex)+2] = start[2]; m.coord[4Int(vertex)+3] = start[3]
        return ST_INTERNAL
    end
    interpolate!(m, vertex, start, stop, v[1])
    for i in previous+1:length(o.deleted)
        t = o.deleted[i] ÷ 4
        o.quality[Int(t)+1] = tet_quality(m, tnode(m, 4t), tnode(m, 4t + 1), tnode(m, 4t + 2), tnode(m, 4t + 3))
        reset_unused!(m, t)
        o.creation[Int(t)+1] = o.date
    end
    return ST_TRUE
end

include("HxtSPR.jl")

@inline function sort4(a, b, c, d)
    a > b && ((a, b) = (b, a)); c > d && ((c, d) = (d, c))
    a > c && ((a, c) = (c, a)); b > d && ((b, d) = (d, b))
    b > c && ((b, c) = (c, b))
    return (a, b, c, d)
end

# ── the optimization driver (hxtOptimizeTetrahedra, one thread) ────────────────
function optimize!(m::Mesh, b::Bbox, threshold::Float64, constrained::Integer)
    o = Opti(m, zeros(m.ntet), zeros(UInt16, m.ntet), zeros(UInt16, m.ntet), 0x0000, threshold,
        UInt64[], UInt32(constrained), 0, 0, 500, UInt64[], UInt64[], UInt64[], UInt16[], UInt32[],
        zeros(35), zeros(35), zeros(4), zeros(4), nothing)
    for t in 0:m.ntet-1
        processed(m, t) && continue
        o.quality[t+1] = tet_quality(m, tnode(m, 4t), tnode(m, 4t + 1), tnode(m, 4t + 2), tnode(m, 4t + 3))
    end
    moore!(m, b, m.nv)
    seed = Ref(UInt32(1))
    doSPR = true
    conflict_ratio = 0.0
    round_in = typemax(UInt64); round_out = typemax(UInt64)
    max_bad_in = UInt64(0); round_bad_in = typemax(UInt64)
    last_spr_ratio = 0.95
    bad = Tuple{UInt64,UInt64}[]
    while true
        # conflicts never happen on one thread: every round rebuilds the list
        empty!(bad)
        for t in 0:m.ntet-1
            (!processed(m, t) && !deleted(m, t) && o.quality[t+1] < threshold) && push!(bad, (UInt64(0), UInt64(t)))
        end
        nbad = UInt64(length(bad))
        nbad > max_bad_in && (max_bad_in = nbad)
        for (_, t) in bad; set_conflict!(m, t); end
        if round_out == 0 || nbad == 0 ||
           (max_bad_in - round_bad_in > 1000 * round_out && round_bad_in > 1000 * round_out)
            break
        end
        ratio = Float64(round_in) / Float64(round_bad_in)
        if doSPR
            last_spr_ratio = ratio
            doSPR = false
            round_out = UInt64(0)
        elseif 200.0 * ratio <= last_spr_ratio * last_spr_ratio * 0.75 + last_spr_ratio * 0.1 + 0.1
            doSPR = true
            round_out = UInt64(0)
        end
        round_in = UInt64(0)
        round_bad_in = nbad
        # threadLocals_update: order the bad cells along the Moore curve
        lcg_f64(seed)
        for i in eachindex(bad)
            t = bad[i][2]
            h = sort4(pad(m, tnode(m, 4t)), pad(m, tnode(m, 4t + 1)), pad(m, tnode(m, 4t + 2)), pad(m, tnode(m, 4t + 3)))
            mindist = typemax(UInt64); middle = UInt64(0)
            for j in 0:3
                d = h[(j+1)%4+1] - h[j+1]
                if d < mindist
                    mindist = d
                    middle = h[j+1] + mindist ÷ 2
                end
            end
            bad[i] = (middle, t)
        end
        sort!(bad; by=first, alg=Base.Sort.DEFAULT_STABLE)
        success = UInt64(0)
        for (_, t) in bad
            (!conflict_flag(m, t) || deleted(m, t)) && continue
            status = ST_TRUE
            if !doSPR
                for edge in 0:5
                    status = edge_removal!(o, t, edge)
                    (status == ST_TRUE || status > ST_INTERNAL) && break
                end
                if status <= ST_INTERNAL
                    for vertex in 0:3
                        status = smoothing!(o, 4t + vertex)
                        (status == ST_TRUE || status > ST_INTERNAL) && break
                    end
                end
            else
                status = spr_opti!(o, t)
            end
            if status == ST_TRUE
                success += 1
            elseif status <= ST_INTERNAL
                unset_conflict!(m, t)
            else
                error("HXT optimization failed with status $status")
            end
        end
        if doSPR
            last = UInt16(o.last_check % 65536)
            o.date = last
            o.last_check = last; o.last_creation = last
        end
        round_in += success
        round_out += success
        conflict_ratio = 0.0
    end
    for t in o.deleted, j in 0:3
        settneigh!(m, 4t + j, NO_ADJACENT)
    end
    remove_deleted!(m)
    return m
end
