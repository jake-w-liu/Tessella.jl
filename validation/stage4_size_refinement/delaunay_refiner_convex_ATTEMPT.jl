# ════════════════════════════════════════════════════════════════════════════════
# Convex-domain size refiner (boundary-protected Delaunay refinement).
# Paste-ready for `module Mesh3D` (uses its internals + `orient3`, already imported).
#
# The triangulation is REBUILT with `delaunay3d(…; perturb=false)` from the growing
# point set each pass — NOT mutated incrementally. Reason (measured): incremental
# `insert_point3!` of the exact, highly-cospherical points that box refinement
# generates (face-diagonal midpoints, right-triangle circumcenters) leaves flat tets
# and vertex-on-facet T-junctions the exact kernel cannot repair, pinning maxedge at
# 2√2. A batch `delaunay3d` of the SAME points yields 0 flat tets and a clean
# empty-circumsphere Delaunay (verified). Real work stays in the exact kernel on
# exact coords; Float64 only CHOOSES points, and every chosen point is proven on-or-
# inside D (EXACT point-in-hull) so the volume is invariant.
# ════════════════════════════════════════════════════════════════════════════════

# squared distance between two coordinate tuples
@inline _rc_dist2(p, q) = (p[1]-q[1])^2 + (p[2]-q[2])^2 + (p[3]-q[3])^2

# tet circumcenter + R²  (Float64; `ok=false` for a flat/degenerate tet).
# Same determinant as MeshTypes.tet_circumradius:
#   A=b-a,B=c-a,C=d-a; den=2·A·(B×C); O=(‖A‖²(B×C)+‖B‖²(C×A)+‖C‖²(A×B))/den.
function _rc_tet_cc(pa, pb, pc, pd)
    A = _subn(pb, pa); B = _subn(pc, pa); C = _subn(pd, pa)
    BxC = _cross(B, C)
    den = 2.0 * _dot(A, BxC)
    den == 0.0 && return (false, (0.0,0.0,0.0), 0.0)
    la = _dot(A,A); lb = _dot(B,B); lc = _dot(C,C)
    CxA = _cross(C, A); AxB = _cross(A, B)
    Ox = (la*BxC[1] + lb*CxA[1] + lc*AxB[1]) / den
    Oy = (la*BxC[2] + lb*CxA[2] + lc*AxB[2]) / den
    Oz = (la*BxC[3] + lb*CxA[3] + lc*AxB[3]) / den
    return (true, (pa[1]+Ox, pa[2]+Oy, pa[3]+Oz), Ox*Ox+Oy*Oy+Oz*Oz)
end

# subface (in-plane) circumcenter + r²  (Float64; `ok=false` for a degenerate face).
#   A=b-a,B=c-a,n=A×B; O=((‖A‖²B−‖B‖²A)×n)/(2‖n‖²).
function _rc_subface_cc(pa, pb, pc)
    A = _subn(pb, pa); B = _subn(pc, pa)
    n = _cross(A, B); nl2 = _dot(n, n)
    nl2 == 0.0 && return (false, (0.0,0.0,0.0), 0.0)
    la = _dot(A,A); lb = _dot(B,B)
    t = (la*B[1]-lb*A[1], la*B[2]-lb*A[2], la*B[3]-lb*A[3])
    cr = _cross(t, n)
    inv = 1.0/(2.0*nl2)
    Ox = cr[1]*inv; Oy = cr[2]*inv; Oz = cr[3]*inv
    return (true, (pa[1]+Ox, pa[2]+Oy, pa[3]+Oz), Ox*Ox+Oy*Oy+Oz*Oz)
end

# Hull subfaces via the ghost shell: one (u,v,w,apex) per live ghost tet. `(u,v,w)`
# is the real hull face (oriented so orient3(u,v,w,interior)=+1); `apex` is the
# incident real tet's fourth vertex. Complete + closed (χ=2), independent of any
# flat interior tets.
function _rc_hull_subfaces(T::Triangulation3)
    subs = NTuple{4,Int32}[]
    @inbounds for g in eachindex(T.alive)
        (T.alive[g] && _is_ghost_tet(T, g)) || continue
        gs = _ghost_slot(T, g)
        f = _face(T, g, gs)                 # 3 real hull-face vertices
        u, v, w = f[1], f[2], f[3]
        r = _nbr(T, g, gs)                  # incident real tet across the hull face
        apex = Int32(0)
        for k in 1:4
            vv = _vert(T, r, k)
            (vv == u || vv == v || vv == w) && continue
            apex = vv; break
        end
        push!(subs, (u, v, w, apex))
    end
    return subs
end

# Precompute (ok, center, r²) for every hull subface (used by size + encroachment).
function _rc_subface_ccs(T::Triangulation3, subs::Vector{NTuple{4,Int32}})
    ccs = Vector{Tuple{Bool,NTuple{3,Float64},Float64}}(undef, length(subs))
    @inbounds for i in eachindex(subs)
        u, v, w, _ = subs[i]
        ccs[i] = _rc_subface_cc(_pt(T,u), _pt(T,v), _pt(T,w))
    end
    return ccs
end

# Subsegments (sharp hull creases): each hull edge is shared by exactly 2 subfaces;
# it is a subsegment iff the two facets are non-coplanar (EXACT orient3 ≠ 0). On
# exact box coords this is 0 on the coplanar face-diagonals and ≠0 on the 12 edges.
function _rc_subsegments(T::Triangulation3, subs::Vector{NTuple{4,Int32}})
    edgeapex = Dict{NTuple{2,Int32}, Tuple{Int32,Int32}}()   # edge → (apex1, apex2|0)
    @inbounds for (u, v, w, _) in subs
        for (a, b, opp) in ((u,v,w), (v,w,u), (w,u,v))
            key = _edgekey(a, b)
            cur = get(edgeapex, key, (Int32(0), Int32(0)))
            edgeapex[key] = cur[1] == 0 ? (opp, Int32(0)) : (cur[1], opp)
        end
    end
    segs = NTuple{2,Int32}[]
    @inbounds for (key, ap) in edgeapex
        (ap[1] != 0 && ap[2] != 0) || continue           # need both incident facets
        a, b = key
        if orient3(_pt(T,a), _pt(T,b), _pt(T,ap[1]), _pt(T,ap[2])) != 0
            push!(segs, key)
        end
    end
    return segs
end

# Relation of coord `p` to the hull (EXACT orient3 over all subfaces):
#   1  ⇒ strictly interior (orient3>0 on every subface)
#   0  ⇒ lies on ∂D (orient3==0 on some subface, none <0)
#  -1  ⇒ strictly outside (orient3<0 on some subface) — must never be inserted.
function _rc_hull_relation(T::Triangulation3, subs::Vector{NTuple{4,Int32}}, p)
    onb = false
    @inbounds for (u, v, w, _) in subs
        o = orient3(_pt(T,u), _pt(T,v), _pt(T,w), p)
        o < 0 && return -1
        o == 0 && (onb = true)
    end
    return onb ? 0 : 1
end

# The hull subface whose diametral sphere `p` is nearest to being inside
# (min ‖p−cc_f‖² − r²_f) — the boundary facet closest to a rejected circumcenter.
function _rc_nearest_subface(T::Triangulation3, subs::Vector{NTuple{4,Int32}},
                             subcc::Vector{Tuple{Bool,NTuple{3,Float64},Float64}}, p)
    bi = 0; bv = Inf
    @inbounds for fi in eachindex(subs)
        ok, cc, r2 = subcc[fi]
        ok || continue
        d = _rc_dist2(cc, p) - r2
        d < bv && (bv = d; bi = fi)
    end
    return bi
end

# subsegment (u,v) encroached by point p ⇔ p in the diametral ball ⇔ (u−p)·(v−p) < 0.
@inline function _rc_seg_encr(pu, pv, p)
    return (pu[1]-p[1])*(pv[1]-p[1]) + (pu[2]-p[2])*(pv[2]-p[2]) + (pu[3]-p[3])*(pv[3]-p[3]) < 0.0
end

# subsegment (a,b) encroached by ANY existing vertex (prototype scan of all vertices).
function _rc_seg_encr_vertex(T::Triangulation3, a::Int32, b::Int32)
    pa = _pt(T, a); pb = _pt(T, b)
    @inbounds for vv in 1:T.nreal
        (vv == a || vv == b) && continue
        _rc_seg_encr(pa, pb, _pt(T, vv)) && return true
    end
    return false
end

@inline function _rc_midpoint(pa, pb)
    return (0.5*(pa[1]+pb[1]), 0.5*(pa[2]+pb[2]), 0.5*(pa[3]+pb[3]))
end

# Midpoint of a subface's longest edge (on ∂D by convexity; a current edge ⇒ its
# midpoint is not an existing vertex). Fallback for obtuse facets whose cc escapes D.
function _rc_subface_longest_mid(T::Triangulation3, u::Int32, v::Int32, w::Int32)
    pu = _pt(T,u); pv = _pt(T,v); pw = _pt(T,w)
    duv = _rc_dist2(pu,pv); dvw = _rc_dist2(pv,pw); dwu = _rc_dist2(pw,pu)
    if duv >= dvw && duv >= dwu
        return _rc_midpoint(pu, pv)
    elseif dvw >= dwu
        return _rc_midpoint(pv, pw)
    else
        return _rc_midpoint(pw, pu)
    end
end

# Midpoint of tet t's longest edge (a current edge; both endpoints ∈ closed convex D
# ⇒ midpoint ∈ D). The always-new fallback for a tet whose circumcenter is unusable.
function _rc_tet_longest_mid(T::Triangulation3, t)
    vs = (_vert(T,t,1), _vert(T,t,2), _vert(T,t,3), _vert(T,t,4))
    best = -1.0; ba = vs[1]; bb = vs[2]
    @inbounds for i in 1:4, j in i+1:4
        d = _rc_dist2(_pt(T,vs[i]), _pt(T,vs[j]))
        d > best && (best = d; ba = vs[i]; bb = vs[j])
    end
    return _rc_midpoint(_pt(T, ba), _pt(T, bb))
end

# The point that SPLITS hull subface fi (returns a coord on ∂D):
#   • if the facet circumcenter encroaches a subsegment → that subsegment's midpoint;
#   • elseif the circumcenter lies in the hull → the circumcenter (in the facet plane);
#   • else (obtuse facet whose cc leaves the hull) → the facet's longest-edge midpoint.
function _rc_subface_split_point(T::Triangulation3, subs::Vector{NTuple{4,Int32}},
                                 segs::Vector{NTuple{2,Int32}},
                                 subcc::Vector{Tuple{Bool,NTuple{3,Float64},Float64}}, fi::Int)
    u, v, w, _ = subs[fi]
    ok, cc, _ = subcc[fi]
    ok || return _rc_subface_longest_mid(T, u, v, w)
    @inbounds for s in segs
        _rc_seg_encr(_pt(T,s[1]), _pt(T,s[2]), cc) && return _rc_midpoint(_pt(T,s[1]), _pt(T,s[2]))
    end
    return _rc_hull_relation(T, subs, cc) >= 0 ? cc : _rc_subface_longest_mid(T, u, v, w)
end

# ── point collection ────────────────────────────────────────────────────────────
# Collect every split point that any size/encroachment rule fires this pass; the
# driver adds them all and REBUILDS once. The domain is convex, so boundary
# refinement adds points ON the (fixed) hull and never changes its shape — a tet's
# interior/boundary classification is pass-invariant and the three priorities can be
# gathered together (still deferring a boundary-encroaching tet circumcenter to a
# boundary split), cutting the rebuild count from O(N) to O(log N). Returns
# (new points, `viol` = did any rule fire — separates genuine convergence from an
# unrepairable degenerate no-op).
function _rc_collect_points(T::Triangulation3, subs::Vector{NTuple{4,Int32}},
                            subcc::Vector{Tuple{Bool,NTuple{3,Float64},Float64}},
                            segs::Vector{NTuple{2,Int32}},
                            hmax::Float64, rho::Float64, rho2::Float64,
                            seen::Set{NTuple{3,Float64}})
    h2 = hmax*hmax; rho_2 = rho*rho
    out = NTuple{3,Float64}[]
    proposed = Set{NTuple{3,Float64}}()
    interior = NTuple{3,Float64}[]     # interior points accepted this pass (for spacing)
    viol = false                       # did any size/encroachment rule fire at all?
    # Only NEW points count. A proposal equal to an existing vertex is a no-op (an
    # exact-Delaunay degeneracy — e.g. a collinear hull-edge T-junction whose midpoint
    # already exists): skip it rather than stall; a later rebuild with more nearby
    # points typically heals it. Also dedups points proposed twice within a pass.
    add!(p) = (p in seen || p in proposed) || (push!(out, p); push!(proposed, p))
    # Interior Steiner points are spacing-filtered to ≥ρ apart PER PASS: batching every
    # bad tet's circumcenter otherwise over-refines ~3× (redundant nearby inserts one
    # rebuild makes moot), blowing up the BigInt-heavy rebuild cost. A tet still too big
    # next pass is refined then — spacing only caps per-pass additions, not convergence.
    function add_interior!(p)
        (p in seen || p in proposed) && return
        @inbounds for q in interior
            _rc_dist2(q, p) < rho_2 && return
        end
        push!(out, p); push!(proposed, p); push!(interior, p)
    end

    # P1 SUBSEGMENT: too long OR encroached by an existing vertex → midpoint.
    @inbounds for s in segs
        a, b = s[1], s[2]
        pa = _pt(T,a); pb = _pt(T,b)
        if _rc_dist2(pa, pb) > h2 || _rc_seg_encr_vertex(T, a, b)
            viol = true; add!(_rc_midpoint(pa, pb))
        end
    end

    # P2 SUBFACE: circumradius too big OR encroached by its incident apex.
    @inbounds for fi in eachindex(subs)
        ok, cc, r2 = subcc[fi]
        ok || continue
        apex = subs[fi][4]
        if r2 > rho_2 || _rc_dist2(_pt(T, apex), cc) < r2
            viol = true; add!(_rc_subface_split_point(T, subs, segs, subcc, fi))
        end
    end

    # P3 TET: positive real tet with R² > rho² (skip flats). A rejected circumcenter
    # DEFERS to a boundary split, never dropped.
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T, t)) || continue
        a=_vert(T,t,1); b=_vert(T,t,2); c=_vert(T,t,3); d=_vert(T,t,4)
        pa=_pt(T,a); pb=_pt(T,b); pc=_pt(T,c); pd=_pt(T,d)
        orient3(pa, pb, pc, pd) == 0 && continue          # flat ⇒ Inf circumradius, skip
        okc, cen, R2 = _rc_tet_cc(pa, pb, pc, pd)
        (okc && R2 > rho2) || continue
        viol = true
        # (a) c encroaches a subsegment → split that subsegment.
        hit = false
        for s in segs
            if _rc_seg_encr(_pt(T,s[1]), _pt(T,s[2]), cen)
                add!(_rc_midpoint(_pt(T,s[1]), _pt(T,s[2]))); hit = true; break
            end
        end
        hit && continue
        # (b) c strictly encroaches a subface → split that subface.
        for fi in eachindex(subs)
            ok, fcc, fr2 = subcc[fi]
            (ok && _rc_dist2(fcc, cen) < fr2) || continue
            add!(_rc_subface_split_point(T, subs, segs, subcc, fi)); hit = true; break
        end
        hit && continue
        # (c)/(d) classify c against ∂D. A strictly-outside c splits the nearest facet;
        # otherwise c (interior or on ∂D) is inserted — it lies in D so volume is kept.
        if _rc_hull_relation(T, subs, cen) == -1
            fi = _rc_nearest_subface(T, subs, subcc, cen)
            add!(fi == 0 ? _rc_tet_longest_mid(T, t) :
                           _rc_subface_split_point(T, subs, segs, subcc, fi))
        else
            add_interior!(cen)          # spacing-filtered genuine interior Steiner point
        end
    end
    return (out, viol)
end

"""
    refine_convex(surface::Mesh; hmax, maxiters=100_000)
        -> (; mesh, terminated, iters, nverts)

Boundary-protected Delaunay size refinement of the CONVEX domain `D` = convex hull
of `surface`'s vertices. Produces a tet [`Mesh`](@ref) with `maxedge ≤ hmax`,
`validate`d, boundary χ=2/closed-manifold, and volume exactly `vol(D)`.

The triangulation is rebuilt with `delaunay3d(…; perturb=false)` from the growing
point set each pass, because incremental insertion of the exact cospherical points
box refinement generates leaves flat / T-junction tets the kernel cannot repair,
whereas a batch build of the same points is a clean Delaunay (verified). Float64 is
used only to CHOOSE which points to insert; every chosen point is proven on-or-inside
`D` (EXACT point-in-hull) so no point escapes `D` and the volume is invariant.

Returns the exported mesh plus `terminated` (false ⇒ the iteration cap was hit, or a
round produced only already-present points — an explicit blocker, e.g. a sharp convex
feature with a small input angle, NOT a fabricated result), the number of rebuild
passes `iters`, and the vertex count. Assumes a convex domain with no small
dihedral/facet angles (a box is 90°).
"""
function refine_convex(surface::Mesh; hmax::Real, maxiters::Integer=100_000)
    hmax > 0 || throw(ArgumentError("refine_convex: hmax must be positive (got $hmax)"))
    nn = size(surface.coords, 2)
    nn >= 4 || throw(ArgumentError("refine_convex: need ≥4 surface vertices (got $nn)"))
    h = float(hmax); rho = h/2; rho2 = rho*rho

    # growing point set (deduplicated on exact coords)
    pts = NTuple{3,Float64}[]
    seen = Set{NTuple{3,Float64}}()
    @inbounds for i in 1:nn
        p = (surface.coords[1,i], surface.coords[2,i], surface.coords[3,i])
        (p in seen) || (push!(pts, p); push!(seen, p))
    end

    T = Triangulation3(Float64[], Float64[], Float64[])   # placeholder; rebuilt below
    it = 0; terminated = false
    while it < maxiters
        it += 1
        xs = Float64[p[1] for p in pts]; ys = Float64[p[2] for p in pts]; zs = Float64[p[3] for p in pts]
        T = delaunay3d(xs, ys, zs; perturb=false)
        subs  = _rc_hull_subfaces(T)
        subcc = _rc_subface_ccs(T, subs)
        segs  = _rc_subsegments(T, subs)
        newpts, viol = _rc_collect_points(T, subs, subcc, segs, h, rho, rho2, seen)
        if isempty(newpts)
            # No NEW points. Either every size/encroachment rule is satisfied
            # (viol=false ⇒ genuinely converged) or the only remaining violations are
            # degenerate no-ops that no new point can fix (viol=true ⇒ explicit blocker,
            # e.g. an unremovable collinear-hull T-junction).
            terminated = !viol; break
        end
        for p in newpts
            push!(pts, p); push!(seen, p)     # newpts are already dup-filtered
        end
    end

    # Export with the anti-fabrication guard: keep positive real tets, forbid any
    # negatively-oriented real tet, drop flats (exactly zero signed volume ⇒ Σvol over
    # kept tets = vol(D)).
    keep = falses(length(T.alive)); nneg = 0
    @inbounds for t in eachindex(T.alive)
        (T.alive[t] && !_is_ghost_tet(T, t)) || continue
        a=_vert(T,t,1); b=_vert(T,t,2); c=_vert(T,t,3); d=_vert(T,t,4)
        o = orient3(_pt(T,a), _pt(T,b), _pt(T,c), _pt(T,d))
        o > 0 && (keep[t] = true)
        o < 0 && (nneg += 1)
    end
    nneg == 0 || error("refine_convex: $nneg negatively-oriented real tet(s) — kernel bug, refusing to fabricate")
    m = to_mesh3(T; keep=keep)
    return (mesh = m, terminated = terminated, iters = it, nverts = T.nreal)
end