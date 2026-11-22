# SPDX-FileCopyrightText: 2016-2020 UCLouvain
# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: GPL-2.0-or-later
# Contains HXT adaptations from Gmsh 4.15.2 contrib/hxt/tetMesh/src/HXTSPR.c
# (small polyhedron reconnection, growing cavity). Upstream author: Célestin Marot.
# Julia adaptation/changes: Tessella contributors, 2026.
# See LICENSES/README.md and NOTICES.md#hxt. Included by HxtOpti.jl.

const SPR_MAX_PTS = 32
const SPR_MAX_TETS = (SPR_MAX_PTS * SPR_MAX_PTS - 3SPR_MAX_PTS - 2) ÷ 2
const NO_FACE = 0xffff
const _NODE_FROM_FACET = ((1, 2, 3, 0), (3, 3, 1, 1), (2, 0, 0, 2))

mutable struct SPRStep
    candidate::Vector{UInt8}
    tet::NTuple{4,UInt8}
    num_candidates::Int
    quality::Float64
end

mutable struct SPR
    coord::Vector{NTuple{3,Float64}}
    user::Vector{UInt32}
    valence::Vector{Int}
    interior::Vector{Bool}
    npts::Int
    ciedges::Vector{NTuple{2,UInt8}}
    citris::Vector{NTuple{4,UInt8}}
    bnd::Vector{NTuple{4,UInt8}}
    tets::Vector{NTuple{4,UInt8}}
    quality::Float64
    steps::Vector{SPRStep}
    nsteps::Int
    faces::Vector{UInt16}
    qualities::Vector{Float64}
    max_search::Int
    search::Int
    adjacencies::Vector{UInt64}
    worst::Float64
end

function SPR()
    SPR(fill((0.0, 0.0, 0.0), SPR_MAX_PTS), zeros(UInt32, SPR_MAX_PTS), zeros(Int, SPR_MAX_PTS),
        falses(SPR_MAX_PTS), 0, NTuple{2,UInt8}[], NTuple{4,UInt8}[], NTuple{4,UInt8}[],
        NTuple{4,UInt8}[], 0.0, [SPRStep(zeros(UInt8, SPR_MAX_PTS - 3), (0, 0, 0, 0), -1, 0.0) for _ in 1:SPR_MAX_TETS+1],
        0, fill(NO_FACE, SPR_MAX_PTS^3), fill(NaN, SPR_MAX_PTS * (SPR_MAX_PTS - 1) * (SPR_MAX_PTS - 2) * (SPR_MAX_PTS - 3) ÷ 24),
        500, 0, UInt64[], 0.0)
end

@inline faceidx(v0, v1, v2) = Int(v0) * SPR_MAX_PTS^2 + Int(v1) * SPR_MAX_PTS + Int(v2) + 1
@inline function add_face_map!(S::SPR, v0, v1, v2, index)
    S.faces[faceidx(v0, v1, v2)] = index
    S.faces[faceidx(v1, v2, v0)] = index
    S.faces[faceidx(v2, v0, v1)] = index
end
@inline remove_face_map!(S::SPR, v0, v1, v2) = add_face_map!(S, v0, v1, v2, NO_FACE)
@inline get_face_map(S::SPR, v0, v1, v2) = S.faces[faceidx(v0, v1, v2)]

function compressed_index(v0, v1, v2, v3)
    v0 = UInt(v0); v1 = UInt(v1); v2 = UInt(v2); v3 = UInt(v3)
    v0 > v1 && ((v0, v1) = (v1, v0)); v2 > v3 && ((v2, v3) = (v3, v2))
    v0 > v2 && ((v0, v2) = (v2, v0)); v1 > v3 && ((v1, v3) = (v3, v1))
    v1 > v2 && ((v1, v2) = (v2, v1))
    aux = v3 * (v3 - 3)
    return Int(v0 + ((v1 * (v1 - 1)) >> 1) + v2 * (v2 - 1) * (v2 - 2) ÷ 6 + aux * (aux + 2) ÷ 24)
end

@inline function sorient(S::SPR, v0, v1, v2, v3)
    (v0 == v1 || v0 == v2 || v0 == v3 || v1 == v2 || v1 == v3 || v2 == v3) && return 0
    return orient3(S.coord[v0+1], S.coord[v1+1], S.coord[v2+1], S.coord[v3+1])
end

function add_face!(S::SPR, v0, v1, v2)
    push!(S.bnd, (UInt8(v0), UInt8(v1), UInt8(v2), UInt8(v0)))
    add_face_map!(S, v0, v1, v2, UInt16(length(S.bnd) - 1))
    S.valence[v0+1] += 1; S.valence[v1+1] += 1; S.valence[v2+1] += 1
end

function remove_face!(S::SPR, id::Int)
    d = S.bnd[id+1]
    remove_face_map!(S, d[1], d[2], d[3])
    S.valence[d[1]+1] -= 1; S.valence[d[2]+1] -= 1; S.valence[d[3]+1] -= 1
    if id != length(S.bnd) - 1
        l = S.bnd[end]
        add_face_map!(S, l[1], l[2], l[3], UInt16(id))
        S.bnd[id+1] = l
    end
    pop!(S.bnd)
end

function get_quality_map(S::SPR, v0, v1, v2, v3)
    sorient(S, v0, v1, v2, v3) >= 0 && return -floatmax(Float64), 0
    index = compressed_index(v0, v1, v2, v3)
    q = S.qualities[index+1]
    isnan(q) || return q, index
    q = aspect_ratio(S.coord[v0+1], S.coord[v1+1], S.coord[v2+1], S.coord[v3+1])
    q < 0.0 && (q = floatmin(Float64))
    S.qualities[index+1] = q
    return q, index
end

@inline tri_edge(t012v0, t012v1, t01v01, t12v01, t20v01) =
    t012v0 * t012v1 < 0 && t012v0 * t01v01 <= 0 && t012v0 * t12v01 <= 0 && t012v0 * t20v01 <= 0

function bbox4(S::SPR, n)
    lo = [S.coord[n[1]+1]...]; hi = copy(lo)
    for k in 2:4, i in 1:3
        c = S.coord[n[k]+1][i]
        c < lo[i] && (lo[i] = c)
        c > hi[i] && (hi[i] = c)
    end
    return lo, hi
end

function tet_tri_intersection(S::SPR, tb, tet, tri)
    lo, hi = tb
    tlo, thi = bbox4(S, tri)
    for i in 1:3
        (lo[i] > thi[i] || hi[i] < tlo[i]) && return false
    end
    o(a, b, c, d) = sorient(S, a, b, c, d)
    idx(x) = x == tet[1] ? 0 : x == tet[2] ? 1 : x == tet[3] ? 2 : x == tet[4] ? 3 : -1
    t0 = idx(tri[1]); t1 = idx(tri[2]); t2 = idx(tri[3])
    tt = [-1, -1, -1, -1]
    t0 != -1 && (tt[t0+1] = 0); t1 != -1 && (tt[t1+1] = 1); t2 != -1 && (tt[t2+1] = 2)
    # C assigns tetN_to_tri in order tri0, tri1, tri2; a later match overwrites
    (t0 != -1 && t1 != -1 && t2 != -1) && return false
    side(f0, f1, f2) = (f0 > 0 || (f0 == 0 && t0 != -1)) && (f1 > 0 || (f1 == 0 && t1 != -1)) && (f2 > 0 || (f2 == 0 && t2 != -1))
    f0v0 = o(tri[1], tet[2], tet[3], tet[4]); f0v1 = o(tri[2], tet[2], tet[3], tet[4]); f0v2 = o(tri[3], tet[2], tet[3], tet[4])
    side(f0v0, f0v1, f0v2) && return false
    f1v0 = o(tet[1], tri[1], tet[3], tet[4]); f1v1 = o(tet[1], tri[2], tet[3], tet[4]); f1v2 = o(tet[1], tri[3], tet[3], tet[4])
    side(f1v0, f1v1, f1v2) && return false
    f2v0 = o(tet[1], tet[2], tri[1], tet[4]); f2v1 = o(tet[1], tet[2], tri[2], tet[4]); f2v2 = o(tet[1], tet[2], tri[3], tet[4])
    side(f2v0, f2v1, f2v2) && return false
    f3v0 = o(tet[1], tet[2], tet[3], tri[1]); f3v1 = o(tet[1], tet[2], tet[3], tri[2]); f3v2 = o(tet[1], tet[2], tet[3], tri[3])
    side(f3v0, f3v1, f3v2) && return false
    ((t0 == -1 && f0v0 <= 0 && f1v0 <= 0 && f2v0 <= 0 && f3v0 <= 0) ||
     (t1 == -1 && f0v1 <= 0 && f1v1 <= 0 && f2v1 <= 0 && f3v1 <= 0) ||
     (t2 == -1 && f0v2 <= 0 && f1v2 <= 0 && f2v2 <= 0 && f3v2 <= 0)) && return true
    v0f = o(tet[1], tri[1], tri[2], tri[3]); v1f = o(tet[2], tri[1], tri[2], tri[3])
    v2f = o(tet[3], tri[1], tri[2], tri[3]); v3f = o(tet[4], tri[1], tri[2], tri[3])
    (v0f > 0 && v1f > 0 && v2f > 0 && v3f > 0) && return false
    (v0f < 0 && v1f < 0 && v2f < 0 && v3f < 0) && return false
    e(a, b) = (o(tet[1], tet[2], a, b), o(tet[1], tet[3], a, b), o(tet[1], tet[4], a, b),
               o(tet[2], tet[3], a, b), o(tet[2], tet[4], a, b), o(tet[3], tet[4], a, b))
    E0 = e(tri[1], tri[2])
    (tri_edge(f0v0, f0v1, -E0[4], -E0[6], E0[5]) || tri_edge(f1v0, f1v1, -E0[3], E0[6], E0[2]) ||
     tri_edge(f2v0, f2v1, -E0[1], -E0[5], E0[3]) || tri_edge(f3v0, f3v1, -E0[2], E0[4], E0[1])) && return true
    E1 = e(tri[2], tri[3])
    (tri_edge(f0v1, f0v2, -E1[4], -E1[6], E1[5]) || tri_edge(f1v1, f1v2, -E1[3], E1[6], E1[2]) ||
     tri_edge(f2v1, f2v2, -E1[1], -E1[5], E1[3]) || tri_edge(f3v1, f3v2, -E1[2], E1[4], E1[1])) && return true
    E2 = e(tri[3], tri[1])
    (tri_edge(f0v2, f0v0, -E2[4], -E2[6], E2[5]) || tri_edge(f1v2, f1v0, -E2[3], E2[6], E2[2]) ||
     tri_edge(f2v2, f2v0, -E2[1], -E2[5], E2[3]) || tri_edge(f3v2, f3v0, -E2[2], E2[4], E2[1])) && return true
    a0, a1, a2, a3 = tt
    ((a0 == -1 && a1 == -1 && tri_edge(v0f == 0 ? v1f : -v0f, -v1f, E0[1], E1[1], E2[1])) ||
     (a0 == -1 && a2 == -1 && tri_edge(v0f == 0 ? v2f : -v0f, -v2f, E0[2], E1[2], E2[2])) ||
     (a0 == -1 && a3 == -1 && tri_edge(v0f == 0 ? v3f : -v0f, -v3f, E0[3], E1[3], E2[3])) ||
     (a1 == -1 && a2 == -1 && tri_edge(v1f == 0 ? v2f : -v1f, -v2f, E0[4], E1[4], E2[4])) ||
     (a1 == -1 && a3 == -1 && tri_edge(v1f == 0 ? v3f : -v1f, -v3f, E0[5], E1[5], E2[5])) ||
     (a2 == -1 && a3 == -1 && tri_edge(v2f == 0 ? v3f : -v2f, -v3f, E0[6], E1[6], E2[6]))) && return true
    return false
end

function tet_edge_intersection(S::SPR, tb, tet, edge)
    lo, hi = tb
    l0 = S.coord[edge[1]+1]; l1 = S.coord[edge[2]+1]
    for i in 1:3
        (hi[i] < l0[i] && hi[i] < l1[i]) && return false
    end
    for i in 1:3
        (lo[i] > l0[i] && lo[i] > l1[i]) && return false
    end
    o(a, b, c, d) = sorient(S, a, b, c, d)
    idx(x) = x == tet[1] ? 0 : x == tet[2] ? 1 : x == tet[3] ? 2 : x == tet[4] ? 3 : -1
    e0 = idx(edge[1]); e1 = idx(edge[2])
    (e0 != -1 && e1 != -1) && return false
    side(f0, f1) = (f0 > 0 || (f0 == 0 && e0 != -1)) && (f1 > 0 || (f1 == 0 && e1 != -1))
    f0v0 = o(edge[1], tet[2], tet[3], tet[4]); f0v1 = o(edge[2], tet[2], tet[3], tet[4])
    side(f0v0, f0v1) && return false
    f1v0 = o(tet[1], edge[1], tet[3], tet[4]); f1v1 = o(tet[1], edge[2], tet[3], tet[4])
    side(f1v0, f1v1) && return false
    f2v0 = o(tet[1], tet[2], edge[1], tet[4]); f2v1 = o(tet[1], tet[2], edge[2], tet[4])
    side(f2v0, f2v1) && return false
    f3v0 = o(tet[1], tet[2], tet[3], edge[1]); f3v1 = o(tet[1], tet[2], tet[3], edge[2])
    side(f3v0, f3v1) && return false
    ((e0 == -1 && f0v0 <= 0 && f1v0 <= 0 && f2v0 <= 0 && f3v0 <= 0) ||
     (e1 == -1 && f0v1 <= 0 && f1v1 <= 0 && f2v1 <= 0 && f3v1 <= 0)) && return true
    ((f0v0 == 0 && f0v1 > 0) || (f0v1 == 0 && f0v0 > 0) || (f1v0 == 0 && f1v1 > 0) || (f1v1 == 0 && f1v0 > 0) ||
     (f2v0 == 0 && f2v1 > 0) || (f2v1 == 0 && f2v0 > 0) || (f3v0 == 0 && f3v1 > 0) || (f3v1 == 0 && f3v0 > 0)) && return false
    c0 = o(tet[1], tet[2], edge[1], edge[2]); c1 = o(tet[1], tet[3], edge[1], edge[2])
    c2 = o(tet[1], tet[4], edge[1], edge[2]); c3 = o(tet[2], tet[3], edge[1], edge[2])
    c4 = o(tet[2], tet[4], edge[1], edge[2]); c5 = o(tet[3], tet[4], edge[1], edge[2])
    (tri_edge(f0v0, f0v1, -c3, -c5, c4) || tri_edge(f1v0, f1v1, -c2, c5, c1) ||
     tri_edge(f2v0, f2v1, -c0, -c4, c2) || tri_edge(f3v0, f3v1, -c1, c3, c0)) && return true
    return false
end

function tet_contains_pt(S::SPR, tb, tet, pt)
    lo, hi = tb
    c = S.coord[pt+1]
    (pt == tet[1] || pt == tet[2] || pt == tet[3] || pt == tet[4] ||
     c[1] < lo[1] || c[2] < lo[2] || c[3] < lo[3] || c[1] > hi[1] || c[2] > hi[2] || c[3] > hi[3] ||
     sorient(S, tet[2], tet[4], tet[3], pt) > 0 || sorient(S, tet[3], tet[4], tet[1], pt) > 0 ||
     sorient(S, tet[4], tet[2], tet[1], pt) > 0 || sorient(S, tet[1], tet[2], tet[3], pt) > 0) && return false
    return true
end

function best_face_heuristic(S::SPR)
    min_valid = 255
    worst = -1
    minimax = floatmax(Float64)
    for id in length(S.bnd)-1:-1:0
        t = S.bnd[id+1]
        valid = 0
        mx = 0.0
        n3 = 0
        while n3 < S.npts && valid < min_valid
            if !(n3 == t[1] || n3 == t[2] || n3 == t[3] || S.valence[n3+1] == 0)
                q, _ = get_quality_map(S, t[1], t[2], t[3], n3)
                if q > S.quality
                    q > mx && (mx = q)
                    valid += 1
                end
            end
            n3 += 1
        end
        valid == 0 && return -1, minimax
        mx < minimax && (minimax = mx)
        if valid < min_valid
            min_valid = valid
            worst = id
        end
    end
    return worst, minimax
end

function compute_candidates!(S::SPR, step::SPRStep, qualmax)
    t = step.tet
    cand = UInt8[]; score = UInt8[]
    for n3 in 0:S.npts-1
        (n3 == t[1] || n3 == t[2] || n3 == t[3] || S.valence[n3+1] == 0) && continue
        q, _ = get_quality_map(S, t[1], t[2], t[3], n3)
        q <= S.quality && continue
        s = 0
        S.interior[n3+1] && (s += 4)
        s += q > step.quality ? 2 : q == step.quality ? 1 : 0
        s += q > qualmax ? 2 : q == qualmax ? 1 : 0
        tet = (t[1], t[2], t[3], UInt8(n3))
        for f in 0:2
            p0 = tet[_NODE_FROM_FACET[1][f+1]+1]; p1 = tet[_NODE_FROM_FACET[2][f+1]+1]; p2 = tet[_NODE_FROM_FACET[3][f+1]+1]
            get_face_map(S, p0, p1, p2) != NO_FACE && (s += 2)
        end
        push!(cand, UInt8(n3)); push!(score, UInt8(s))
    end
    count = zeros(Int, 16)
    for s in score; count[s+1] += 1; end
    sum = 0
    for i in 1:16
        c = count[i]; count[i] = sum; sum += c
    end
    for i in eachindex(cand)
        step.candidate[count[score[i]+1]+1] = cand[i]
        count[score[i]+1] += 1
    end
    step.num_candidates = length(cand)
    step.tet = (t[1], t[2], t[3], UInt8(S.npts))
end

function new_tetrahedralization!(S::SPR, step::SPRStep)
    S.quality = step.quality
    resize!(S.tets, S.nsteps)
    for i in 1:S.nsteps
        S.tets[i] = S.steps[i].tet
    end
end

function detect_interior!(S::SPR)
    for i in 1:S.npts
        if S.valence[i] == 0
            S.interior[i] = true
            S.valence[i] = 1
        end
    end
end

function init_step!(S::SPR)
    S.nsteps = 0
    S.steps[1].num_candidates = -1
    S.steps[1].quality = floatmax(Float64)
end

function spr_advanced!(S::SPR)
    step = S.steps[S.nsteps+1]
    while true
        S.search >= S.max_search && return ST_INTERNAL
        rewind = false
        if isempty(S.bnd)
            new_tetrahedralization!(S, step)
            rewind = true
        elseif step.num_candidates < 0
            id, qualmax = best_face_heuristic(S)
            if id == -1
                step.num_candidates = -1
                rewind = true
            else
                step.tet = S.bnd[id+1]
                compute_candidates!(S, step, qualmax)
                remove_face!(S, id)
            end
        end
        if step.num_candidates > 0
            step.num_candidates -= 1
            t = step.tet
            step.tet = (t[1], t[2], t[3], step.candidate[step.num_candidates+1])
            t = step.tet
            q, index = get_quality_map(S, t[1], t[2], t[3], t[4])
            (q <= S.quality || q == -floatmax(Float64)) && continue
            tb = bbox4(S, t)
            stop = false
            for pt in 0:S.npts-1
                if tet_contains_pt(S, tb, t, pt)
                    stop = true; S.qualities[index+1] = -floatmax(Float64); break
                end
            end
            stop && continue
            for e in S.ciedges
                if tet_edge_intersection(S, tb, t, e)
                    stop = true; S.qualities[index+1] = -floatmax(Float64); break
                end
            end
            stop && continue
            for tri in S.citris
                if tet_tri_intersection(S, tb, t, tri)
                    stop = true; S.qualities[index+1] = -floatmax(Float64); break
                end
            end
            stop && continue
            for tri in S.bnd
                if tet_tri_intersection(S, tb, t, tri)
                    stop = true; break
                end
            end
            stop && continue
            S.search += 1
            for f in 0:2
                p0 = t[_NODE_FROM_FACET[1][f+1]+1]; p1 = t[_NODE_FROM_FACET[2][f+1]+1]; p2 = t[_NODE_FROM_FACET[3][f+1]+1]
                index2 = get_face_map(S, p0, p1, p2)
                index2 == NO_FACE ? add_face!(S, p0, p2, p1) : remove_face!(S, Int(index2))
            end
            old = step.quality
            S.nsteps += 1
            step = S.steps[S.nsteps+1]
            step.quality = q < old ? q : old
            step.num_candidates = -1
        elseif step.num_candidates == 0
            add_face!(S, step.tet[1], step.tet[2], step.tet[3])
            rewind = true
        end
        if rewind
            S.nsteps == 0 && return ST_TRUE
            S.nsteps -= 1
            step = S.steps[S.nsteps+1]
            t = step.tet
            for f in 0:2
                p0 = t[_NODE_FROM_FACET[1][f+1]+1]; p1 = t[_NODE_FROM_FACET[2][f+1]+1]; p2 = t[_NODE_FROM_FACET[3][f+1]+1]
                index2 = get_face_map(S, p0, p2, p1)
                index2 == NO_FACE ? add_face!(S, p0, p1, p2) : remove_face!(S, Int(index2))
            end
            step.quality <= S.quality && (step.num_candidates = 0)
        end
    end
end

function spr_rewind!(S::SPR)
    step = S.steps[S.nsteps+1]
    while true
        step.num_candidates >= 0 && add_face!(S, step.tet[1], step.tet[2], step.tet[3])
        S.nsteps == 0 && return
        S.nsteps -= 1
        step = S.steps[S.nsteps+1]
        t = step.tet
        for f in 0:2
            p0 = t[_NODE_FROM_FACET[1][f+1]+1]; p1 = t[_NODE_FROM_FACET[2][f+1]+1]; p2 = t[_NODE_FROM_FACET[3][f+1]+1]
            index2 = get_face_map(S, p0, p2, p1)
            index2 == NO_FACE ? add_face!(S, p0, p1, p2) : remove_face!(S, Int(index2))
        end
    end
end

# ── SPR around a bad cell of the mesh ──────────────────────────────────────────
@inline known_quality!(S::SPR, v0, v1, v2, v3, q) = (S.qualities[compressed_index(v0, v1, v2, v3)+1] = q)

function spr_init!(S::SPR, o::Opti, bad::UInt64)
    m = o.m
    S.worst = o.quality[Int(bad)+1]
    fill!(S.faces, NO_FACE); fill!(S.qualities, NaN)
    S.npts = 4
    empty!(S.ciedges); empty!(S.citris); empty!(S.bnd)
    S.tets = [(0x00, 0x01, 0x02, 0x03)]
    S.quality = o.quality[Int(bad)+1]
    S.search = 0
    S.max_search = o.max_search_nodes
    setdeleted!(m, bad)
    push!(o.deleted, bad)
    for i in 0:3
        node = tnode(m, 4bad + i)
        S.coord[i+1] = point(m, node); S.user[i+1] = node; S.valence[i+1] = 0; S.interior[i+1] = false
    end
    known_quality!(S, 0, 1, 2, 3, o.quality[Int(bad)+1])
    S.adjacencies = UInt64[]
    for i in 0:3
        add_face!(S, node0_of_facet(i), node1_of_facet(i), node2_of_facet(i))
        push!(S.adjacencies, tneigh(m, 4bad + i))
    end
end

function find_best_point(S::SPR, o::Opti)
    m = o.m
    adj = S.adjacencies
    points = UInt32[]; connected = Int[]
    for i in 1:length(S.bnd)
        a = adj[i]
        (a == NO_ADJACENT || facet_constraint(m, a ÷ 4, a % 4) != 0) && continue
        node = tnode(m, a)
        p = findfirst(==(node), points)
        if p === nothing
            push!(points, node); push!(connected, 1)
        else
            connected[p] += 1
        end
    end
    most = UInt32(0); most_n = 0; equality = false
    for i in eachindex(points)
        if connected[i] > most_n
            most_n = connected[i]; most = points[i]; equality = false
        elseif connected[i] == most_n
            equality = true
        end
    end
    most_n == 0 && return nothing
    if equality
        best = -floatmax(Float64)
        for i in eachindex(points)
            connected[i] == most_n || continue
            node = points[i]; score = 0.0
            for t in 1:length(S.bnd)
                a = adj[t]
                (a == NO_ADJACENT || tnode(m, a) != node || facet_constraint(m, a ÷ 4, a % 4) != 0) && continue
                score -= o.quality[Int(a÷4)+1]
            end
            if score > best
                best = score; most = node
            end
        end
    end
    return most
end

function need_new_edge(m::Mesh, tet::UInt64, in::Int, out::Int)
    e = edge_of_facets(in, out)
    edge_constraint(m, tet, e) == 0 && return -1
    first = tet
    while true
        facet_constraint(m, tet, out) != 0 && return -1
        opp = tnode(m, 4tet + in)
        nb = tneigh(m, 4tet + out)
        (nb == NO_ADJACENT || !deleted(m, nb ÷ 4)) && return -1
        tet = nb ÷ 4; in = Int(nb % 4)
        out = 0
        while out < 3 && tnode(m, 4tet + out) != opp; out += 1; end
        tet == first && break
    end
    return e
end

function attach_best_point!(S::SPR, o::Opti)
    m = o.m
    newnode = find_best_point(S, o)
    newnode === nothing && return ST_INTERNAL
    nv = S.npts
    S.npts += 1
    S.coord[nv+1] = point(m, newnode); S.user[nv+1] = newnode; S.valence[nv+1] = 0; S.interior[nv+1] = false
    adj = S.adjacencies
    i = 0
    while i < length(S.bnd)
        a = adj[i+1]
        mt = a ÷ 4; facet = Int(a % 4)
        if a == NO_ADJACENT || facet_constraint(m, mt, facet) != 0
            i += 1; continue
        end
        opp = tnode(m, a)
        oppv = findfirst(==(opp), view(S.user, 1:S.npts))
        if oppv === nothing
            i += 1; continue
        end
        oppv -= 1
        setdeleted!(m, mt)
        push!(o.deleted, mt)
        tri = S.bnd[i+1]
        tv = [0x00, 0x00, 0x00, 0x00]
        tv[facet+1] = UInt8(oppv)
        for j in 0:2
            if S.user[tri[j+1]+1] == tnode(m, 4mt + node0_of_facet(facet))
                tv[node0_of_facet(facet)+1] = tri[j+1]
                tv[node1_of_facet(facet)+1] = tri[(j+2)%3+1]
                tv[node2_of_facet(facet)+1] = tri[(j+1)%3+1]
                break
            end
        end
        push!(S.tets, (tv[1], tv[2], tv[3], tv[4]))
        known_quality!(S, tv[1], tv[2], tv[3], tv[4], o.quality[Int(mt)+1])
        if o.quality[Int(mt)+1] < S.worst
            S.quality = o.quality[Int(mt)+1]
            S.worst = o.quality[Int(mt)+1]
        end
        adj[i+1] = adj[length(S.bnd)]
        remove_face!(S, i)
        for f in 0:3
            f == facet && continue
            p0 = tv[_NODE_FROM_FACET[1][f+1]+1]; p1 = tv[_NODE_FROM_FACET[2][f+1]+1]; p2 = tv[_NODE_FROM_FACET[3][f+1]+1]
            index = get_face_map(S, p0, p2, p1)
            if index == NO_FACE
                n = length(S.bnd)
                length(adj) <= n && resize!(adj, n + 1)
                adj[n+1] = tneigh(m, 4mt + f)
                add_face!(S, p0, p1, p2)
            else
                index = Int(index)
                if facet_constraint(m, mt, f) != 0
                    adj[index+1] = adj[length(S.bnd)]
                    remove_face!(S, index)
                    push!(S.citris, (p0, p1, p2, p0))
                else
                    for f2 in 0:3
                        f2 == f && continue
                        e = need_new_edge(m, mt, f, f2)
                        if e != -1
                            n0, n1 = _NODES_OF_EDGE[e+1]
                            push!(S.ciedges, (tv[n0+1], tv[n1+1]))
                        end
                    end
                    adj[index+1] = adj[length(S.bnd)]
                    remove_face!(S, index)
                    index < i && (i = index)
                end
            end
        end
    end
    return ST_TRUE
end

function constrained_edge_map(m::Mesh, S::SPR)
    map = falses(SPR_MAX_PTS, SPR_MAX_PTS)
    for (a, b) in S.ciedges
        map[a+1, b+1] = true; map[b+1, a+1] = true
    end
    for i in 1:length(S.bnd)
        a = S.adjacencies[i]
        a == NO_ADJACENT && continue
        mt = a ÷ 4; facet = Int(a % 4)
        ind = (node0_of_facet(facet), node2_of_facet(facet), node1_of_facet(facet))
        tri = S.bnd[i]
        rot = 0
        while rot < 3 && S.user[tri[rot+1]+1] != tnode(m, 4mt + ind[1])
            rot += 1
        end
        for j in 0:2
            n0 = ind[j+1]; n1 = ind[(j+1)%3+1]
            if edge_constraint(m, mt, 5 - edge_of_facets(n0, n1)) != 0
                # HXT indexes the second coordinate by the rotation itself
                map[tri[(rot+j)%3+1]+1, (rot+j+1)%3+1] = true
                map[(rot+j+1)%3+1, tri[(rot+j)%3+1]+1] = true
            end
        end
    end
    return map
end

function rebuild_mesh!(S::SPR, o::Opti, bad::UInt64)
    m = o.m
    color = m.color[Int(bad)+1]
    last = UInt16(o.last_check % 65536)
    o.last_creation = last
    ntets = length(S.tets)
    ntets > length(o.deleted) && new_deleted!(o, ntets)
    base = length(o.deleted) - ntets
    newtet(i) = o.deleted[base+i+1]
    cmap = constrained_edge_map(m, S)
    for i in 0:ntets-1
        mt = newtet(i)
        v = S.tets[i+1]
        for j in 0:3; settnode!(m, 4mt + j, S.user[v[j+1]+1]); end
        m.color[Int(mt)+1] = color
        settflag!(m, mt, 0)
        o.quality[Int(mt)+1] = S.qualities[compressed_index(v[1], v[2], v[3], v[4])+1]
        o.creation[Int(mt)+1] = last; o.check[Int(mt)+1] = last
        for j in 0:3
            a, b, c = v[node0_of_facet(j)+1], v[node1_of_facet(j)+1], v[node2_of_facet(j)+1]
            index = get_face_map(S, a, b, c)
            remove_face_map!(S, a, b, c)
            if index == NO_FACE
                settneigh!(m, 4mt + j, NO_ADJACENT)
            else
                adj = S.adjacencies[Int(index)+1]
                if adj == NO_ADJACENT
                    settneigh!(m, 4mt + j, NO_ADJACENT)
                    settflag!(m, mt, tflag(m, mt) | (UInt16(0x100) << j))
                else
                    settneigh!(m, 4mt + j, adj)
                    settneigh!(m, adj, 4mt + j)
                    facet_constraint(m, adj ÷ 4, adj % 4) != 0 && settflag!(m, mt, tflag(m, mt) | (UInt16(0x100) << j))
                end
            end
        end
        for e in 0:5
            n0, n1 = _NODES_OF_EDGE[e+1]
            cmap[v[n0+1]+1, v[n1+1]+1] && settflag!(m, mt, tflag(m, mt) | (UInt16(1) << e))
        end
    end
    for i in 0:ntets-1
        mt = newtet(i); v = S.tets[i+1]
        for j in 0:3
            tneigh(m, 4mt + j) == NO_ADJACENT &&
                add_face_map!(S, v[node0_of_facet(j)+1], v[node1_of_facet(j)+1], v[node2_of_facet(j)+1], UInt16(4i + j))
        end
    end
    for t in S.citris
        adj0 = get_face_map(S, t[1], t[2], t[3]); adj1 = get_face_map(S, t[1], t[3], t[2])
        tet0 = newtet(Int(adj0) ÷ 4); facet0 = Int(adj0) % 4
        tet1 = newtet(Int(adj1) ÷ 4); facet1 = Int(adj1) % 4
        settflag!(m, tet0, tflag(m, tet0) | (UInt16(0x100) << facet0))
        settflag!(m, tet1, tflag(m, tet1) | (UInt16(0x100) << facet1))
        settneigh!(m, 4tet0 + facet0, 4tet1 + facet1)
        settneigh!(m, 4tet1 + facet1, 4tet0 + facet0)
    end
    for i in 0:ntets-1
        mt = newtet(i); v = S.tets[i+1]
        for j in 0:3
            tneigh(m, 4mt + j) == NO_ADJACENT || continue
            adj = get_face_map(S, v[node0_of_facet(j)+1], v[node2_of_facet(j)+1], v[node1_of_facet(j)+1])
            adj != NO_FACE && settneigh!(m, 4mt + j, 4newtet(Int(adj) ÷ 4) + Int(adj) % 4)
        end
    end
    resize!(o.deleted, length(o.deleted) - ntets)
end

function deleted_reset!(o::Opti, previous::Int)
    for i in previous+1:length(o.deleted)
        unsetdeleted!(o.m, o.deleted[i])
    end
    resize!(o.deleted, previous)
end

function set_check_date!(o::Opti, bad::UInt64)
    o.check[Int(bad)+1] = UInt16((o.last_creation + 1) % 65536)
    o.last_check = o.last_creation + 1
end

const _SPR_BUFFER = Ref{Union{Nothing,SPR}}(nothing)

function spr_opti!(o::Opti, bad::UInt64)
    S = something(_SPR_BUFFER[], SPR())
    _SPR_BUFFER[] = S
    previous = length(o.deleted)
    lastcheck = o.check[Int(bad)+1]
    dates_ok = false
    spr_init!(S, o, bad)
    while true
        status = attach_best_point!(S, o)
        if status != ST_TRUE
            deleted_reset!(o, previous)
            return status
        end
        if !dates_ok
            for i in previous+1:length(o.deleted)
                if o.creation[Int(o.deleted[i])+1] >= lastcheck
                    dates_ok = true; break
                end
            end
        end
        if dates_ok
            detect_interior!(S)
            init_step!(S)
            original = copy(S.bnd)
            adjacencies = S.adjacencies[1:length(S.bnd)]
            status = spr_advanced!(S)
            if status == ST_INTERNAL
                if S.worst < S.quality
                    spr_rewind!(S)
                else
                    deleted_reset!(o, previous)
                    set_check_date!(o, bad)
                    return ST_INTERNAL
                end
            end
            for i in 1:length(S.bnd)
                t = original[i]
                index = get_face_map(S, t[1], t[2], t[3])
                S.adjacencies[Int(index)+1] = adjacencies[i]
            end
            S.worst < S.quality && break
        end
        if S.npts == SPR_MAX_PTS
            deleted_reset!(o, previous)
            set_check_date!(o, bad)
            return ST_INTERNAL
        end
    end
    rebuild_mesh!(S, o, bad)
    return ST_TRUE
end
