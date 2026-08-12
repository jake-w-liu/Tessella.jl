"""
    ExactMesh3D

Exact-coordinate 3-D Delaunay tetrahedralization on `Rational{BigInt}` coordinates —
the kernel that lets boundary-recovery Steiner points sit at non-Float64 rational
positions (e.g. the midpoint of two float vertices) and stay EXACTLY on-feature, which
the Float64 kernel cannot (rounding moves them off).

Construction: a single bounding **super-tetrahedron** provably containing every input
point is the (trivially Delaunay) start; points are inserted by exact Bowyer–Watson —
the cavity is the set of tets whose circumsphere contains the point (`insphere_rat`),
its boundary faces (each incident to one cavity tet) become the new tets, each made
positively oriented by `orient3_rat`. Every decision is an exact `Rational{BigInt}`
predicate with the shared +ε SoS, so cospherical/coplanar degeneracies (pervasive for
box/prism caps) are broken deterministically with NO coordinate jitter and NO ghost
tetrahedra. Super-tet-incident tets are dropped at the end, leaving the Delaunay
tetrahedralization of the input points (their convex hull, filled).
"""
module ExactMesh3D

using ..Predicates: orient3_rat, insphere_rat

export RB, delaunay3d_exact, is_delaunay_exact

const RB = Rational{BigInt}

@inline _sort3(a, b, c) = a <= b ? (a <= c ? (c <= b ? (a, c, b) : (a, b, c)) : (c, a, b)) :
                                   (b <= c ? (c <= a ? (b, c, a) : (b, a, c)) : (c, b, a))

"""
    delaunay3d_exact(pts::Vector{NTuple{3,RB}}) -> Vector{NTuple{4,Int}}

Delaunay tetrahedralization of the exact `Rational{BigInt}` points. Returns tets as
1-based index quadruples into `pts`, each positively oriented (`orient3_rat > 0`).
Requires ≥ 4 points, not all coplanar.
"""
function delaunay3d_exact(pts::Vector{NTuple{3,RB}})
    n = length(pts)
    n >= 4 || throw(ArgumentError("delaunay3d_exact: need ≥ 4 points (got $n)"))

    # bounding super-tetrahedron that provably contains every point. With m = min over
    # ALL coordinates and D = (max over all coords) − m, the box [m, m+D]³ ⊇ the point
    # cloud, and the tet s1=(m−D,m−D,m−D), s2=(m+Km·D,…), … with K=10 contains that box
    # (the plane through s2,s3,s4 is x+y+z = 3m+8D ≥ 3(m+D); the other three planes are
    # x=m−D, y=m−D, z=m−D, all below the box). K≥5 suffices; 10 is margin.
    m = pts[1][1]; Mx = pts[1][1]
    @inbounds for p in pts, k in 1:3
        v = p[k]; v < m && (m = v); v > Mx && (Mx = v)
    end
    D = Mx - m; D == 0 && (D = RB(1))
    s1 = (m - D,     m - D,     m - D)
    s2 = (m + 10D,   m - D,     m - D)
    s3 = (m - D,     m + 10D,   m - D)
    s4 = (m - D,     m - D,     m + 10D)
    X = Vector{NTuple{3,RB}}(undef, n + 4)
    @inbounds for i in 1:n; X[i] = pts[i]; end
    X[n+1] = s1; X[n+2] = s2; X[n+3] = s3; X[n+4] = s4

    tets = NTuple{4,Int}[]
    alive = Bool[]
    # initial super-tet, positively oriented
    if orient3_rat(X[n+1], X[n+2], X[n+3], X[n+4], n+1, n+2, n+3, n+4) > 0
        push!(tets, (n+1, n+2, n+3, n+4))
    else
        push!(tets, (n+1, n+2, n+4, n+3))
    end
    push!(alive, true)

    cav = Int[]
    for i in 1:n
        p = X[i]
        empty!(cav)
        @inbounds for ti in eachindex(tets)
            alive[ti] || continue
            t = tets[ti]
            insphere_rat(X[t[1]], X[t[2]], X[t[3]], X[t[4]], p, t[1], t[2], t[3], t[4], i) > 0 && push!(cav, ti)
        end
        isempty(cav) && error("delaunay3d_exact: empty cavity for point $i (super-tet containment bug)")
        # boundary faces of the cavity: a face incident to exactly ONE cavity tet.
        fcount = Dict{NTuple{3,Int},Int}()
        fverts = Dict{NTuple{3,Int},NTuple{3,Int}}()
        @inbounds for ti in cav
            t = tets[ti]
            for f in ((t[2],t[3],t[4]), (t[1],t[4],t[3]), (t[1],t[2],t[4]), (t[1],t[3],t[2]))
                key = _sort3(f[1], f[2], f[3])
                fcount[key] = get(fcount, key, 0) + 1
                fverts[key] = f
            end
        end
        for ti in cav; alive[ti] = false; end
        @inbounds for (key, c) in fcount
            c == 1 || continue
            f = fverts[key]
            nt = orient3_rat(X[f[1]], X[f[2]], X[f[3]], p, f[1], f[2], f[3], i) > 0 ?
                 (f[1], f[2], f[3], i) : (f[1], f[3], f[2], i)
            push!(tets, nt); push!(alive, true)
        end
    end

    # keep only all-real (super-tet-free) tets, flipped to the Mesh orientation
    # convention. Internally tets are `orient3_rat > 0` (so the insphere cavity test has
    # the right sign); the 4×4-homogeneous orient sign is the NEGATIVE of the standard
    # signed volume `det[b−a,c−a,d−a]/6` = `MeshTypes.tet_signed_volume`, so one vertex
    # swap flips each output tet to positive signed volume (what a valid Mesh needs).
    out = NTuple{4,Int}[]
    @inbounds for ti in eachindex(tets)
        alive[ti] || continue
        t = tets[ti]
        (t[1] <= n && t[2] <= n && t[3] <= n && t[4] <= n) || continue
        push!(out, (t[1], t[2], t[4], t[3]))
    end
    return out
end

"""
    is_delaunay_exact(pts, tets) -> (ok::Bool, nviol::Int)

Exact empty-circumsphere oracle: no input point lies strictly inside any tet's
circumsphere (`insphere_rat`). This DEFINES the Delaunay property, exactly.
"""
function is_delaunay_exact(pts::Vector{NTuple{3,RB}}, tets::Vector{NTuple{4,Int}})
    n = length(pts); nviol = 0
    @inbounds for t0 in tets
        a, b, c, d = t0[1], t0[2], t0[3], t0[4]
        # normalize to orient3_rat > 0 so insphere_rat > 0 means "strictly inside"
        orient3_rat(pts[a], pts[b], pts[c], pts[d], a, b, c, d) < 0 && ((c, d) = (d, c))
        for v in 1:n
            (v == a || v == b || v == c || v == d) && continue
            insphere_rat(pts[a], pts[b], pts[c], pts[d], pts[v], a, b, c, d, v) > 0 && (nviol += 1)
        end
    end
    return (nviol == 0, nviol)
end

end # module ExactMesh3D
