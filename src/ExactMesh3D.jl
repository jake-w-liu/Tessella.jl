"""
    ExactMesh3D

Exact-coordinate 3-D Delaunay tetrahedralization on `Rational{BigInt}` coordinates —
the kernel that lets boundary-recovery Steiner points sit at non-Float64 rational
positions (e.g. the midpoint of two float vertices) and stay EXACTLY on-feature, which
the Float64 kernel cannot (rounding moves them off).

Construction: a bounding **super-tetrahedron** provably containing every input point
is the (trivially Delaunay) start; points are inserted by exact Bowyer–Watson —
the cavity is the set of tets whose circumsphere contains the point (`insphere_rat`),
its boundary faces (each incident to one cavity tet) become the new tets, each made
positively oriented by `orient3_rat`. Every decision is an exact `Rational{BigInt}`
predicate with the shared +ε SoS, so cospherical/coplanar degeneracies (pervasive for
box/prism caps) are broken deterministically with NO coordinate jitter and NO ghost
tetrahedra. The super-tet scale is conditioned on the exact affine thickness of the
cloud and grows geometrically when necessary.  A topology/completeness gate rejects
any build that omits an input point, contains a flat/duplicate tet, or has a face with
incidence greater than two; super-tet-incident tets are dropped only after that gate
can be satisfied.  Thus a non-coplanar input never returns a silent partial/empty mesh.
"""
module ExactMesh3D

using ..Predicates: orient3_rat, insphere_rat
using ..MeshTypes: _tet_topology

export RB, delaunay3d_exact, is_delaunay_exact

const RB = Rational{BigInt}
const _MAX_SUPER_GROWTHS = 32

@inline _sort3(a, b, c) = a <= b ? (a <= c ? (c <= b ? (a, c, b) : (a, b, c)) : (c, a, b)) :
                                   (b <= c ? (c <= a ? (b, c, a) : (b, a, c)) : (c, b, a))

@inline function _sort4(a, b, c, d)
    a > b && ((a, b) = (b, a)); c > d && ((c, d) = (d, c))
    a > c && ((a, c) = (c, a)); b > d && ((b, d) = (d, b))
    b > c && ((b, c) = (c, b))
    return (a, b, c, d)
end

@inline _sub3(a, b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline _cross3(a, b) = (a[2]*b[3]-a[3]*b[2],
                         a[3]*b[1]-a[1]*b[3],
                         a[1]*b[2]-a[2]*b[1])
@inline _normsq3(a) = _dot3(a, a)

# Exact physical orientation determinant det[b-a,c-a,d-a].  Unlike `orient3_rat`,
# this deliberately has NO SoS tie-break: zero means a physically flat tetrahedron.
@inline _orient3_det(a, b, c, d) = _dot3(_sub3(b, a), _cross3(_sub3(c, a), _sub3(d, a)))

# Find a well-conditioned exact affine basis and use its dimensionless thickness to
# choose the first super-tet scale.  For a base of diameter D and altitude h,
# D^3/|det| is O(D/h), exactly the factor by which a sliver's circumsphere exceeds
# the cloud.  Maximising line/area/volume at each step avoids an accidentally tiny
# basis.  `nothing` means affine dimension < 3, for which there are no volume tets.
function _initial_super_scale(pts::Vector{NTuple{3,RB}}, D::RB)
    a = pts[1]
    ib = 0; bestd = zero(RB)
    @inbounds for i in 2:length(pts)
        q = _normsq3(_sub3(pts[i], a))
        q > bestd && (bestd = q; ib = i)
    end
    ib == 0 && return nothing
    ab = _sub3(pts[ib], a)

    ic = 0; besta = zero(RB)
    @inbounds for i in eachindex(pts)
        q = _normsq3(_cross3(ab, _sub3(pts[i], a)))
        q > besta && (besta = q; ic = i)
    end
    ic == 0 && return nothing

    bestv = zero(RB)
    @inbounds for i in eachindex(pts)
        q = abs(_orient3_det(a, pts[ib], pts[ic], pts[i]))
        q > bestv && (bestv = q)
    end
    bestv == 0 && return nothing

    condition = D*D*D / bestv
    kcondition = cld(numerator(condition), denominator(condition))
    # Preserve the historically sufficient, smaller K=1 construction for ordinary
    # clouds (condition ≤ 10); larger aspect ratios start at their measured affine
    # condition instead of paying for several doomed rebuilds.
    return kcondition <= 10 ? big(1) : kcondition
end

# Float64 circumsphere (centre, radius²) of a tet, for the cavity pre-filter ONLY.
# Returns r²=Inf for a (near-)degenerate tet so it is NEVER pruned — the exact
# `insphere_rat` then decides. The filter is used only to SKIP tets where the query is
# *safely* outside (see `_maybe_in_sphere`), so the exact result is unchanged. Takes
# Float64 coords ALREADY SHIFTED to a per-axis local origin (see `delaunay3d_exact`): the
# shift is done in exact rational, so these are small, well-conditioned values and the
# differences below carry no catastrophic cancellation (the far-from-origin failure mode).
@inline function _fcircum(a::NTuple{3,Float64}, b::NTuple{3,Float64}, c::NTuple{3,Float64}, d::NTuple{3,Float64})
    ax=a[1]; ay=a[2]; az=a[3]
    Ax=b[1]-ax; Ay=b[2]-ay; Az=b[3]-az
    Bx=c[1]-ax; By=c[2]-ay; Bz=c[3]-az
    Cx=d[1]-ax; Cy=d[2]-ay; Cz=d[3]-az
    bcx=By*Cz-Bz*Cy; bcy=Bz*Cx-Bx*Cz; bcz=Bx*Cy-By*Cx
    denom=2.0*(Ax*bcx+Ay*bcy+Az*bcz)
    # scale-relative degeneracy guard: a tet whose signed volume is a vanishing
    # fraction of its edge cube has an ill-conditioned circumcentre → don't prune it.
    e2=max(Ax*Ax+Ay*Ay+Az*Az, Bx*Bx+By*By+Bz*Bz, Cx*Cx+Cy*Cy+Cz*Cz)
    (isfinite(denom) && abs(denom) > 1e-9*(e2^1.5 + 1e-300)) || return (0.0,0.0,0.0,Inf)
    cax=Cy*Az-Cz*Ay; cay=Cz*Ax-Cx*Az; caz=Cx*Ay-Cy*Ax
    abx=Ay*Bz-Az*By; aby=Az*Bx-Ax*Bz; abz=Ax*By-Ay*Bx
    la=Ax*Ax+Ay*Ay+Az*Az; lb=Bx*Bx+By*By+Bz*Bz; lc=Cx*Cx+Cy*Cy+Cz*Cz
    ox=(la*bcx+lb*cax+lc*abx)/denom; oy=(la*bcy+lb*cay+lc*aby)/denom; oz=(la*bcz+lb*caz+lc*abz)/denom
    r2=ox*ox+oy*oy+oz*oz
    isfinite(r2) || return (0.0,0.0,0.0,Inf)
    return (ax+ox, ay+oy, az+oz, r2)
end

# Heuristic pre-filter: TRUE unless the query is well outside the tet's Float64
# circumsphere. Borderline and ill-conditioned cases keep the exact test. Because a
# fixed Float64 slack is not a formal roundoff proof, each completed build is checked
# by `is_delaunay_exact`; a failed check is rebuilt with this filter disabled.
@inline function _maybe_in_sphere(cs::NTuple{4,Float64}, px::Float64, py::Float64, pz::Float64)
    cs[4] == Inf && return true
    dx=px-cs[1]; dy=py-cs[2]; dz=pz-cs[3]
    (dx*dx+dy*dy+dz*dz) <= cs[4]*(1.0+1e-6) + 1e-12*(1.0+cs[4])
end

"""
    delaunay3d_exact(pts::Vector{NTuple{3,RB}}) -> Vector{NTuple{4,Int}}

Delaunay tetrahedralization of the exact `Rational{BigInt}` points. Returns tets as
1-based index quadruples into `pts`, each physically positively oriented. Requires
≥ 4 pairwise-distinct points. An affine-dimension < 3 cloud has no volume cells and
returns an empty vector. A non-coplanar cloud either returns a complete, gated tet
complex using every input point or throws an explicit error; it never silently returns
a partial/empty result.
"""
function delaunay3d_exact(pts::Vector{NTuple{3,RB}})
    n = length(pts)
    n >= 4 || throw(ArgumentError("delaunay3d_exact: need ≥ 4 points (got $n)"))
    n<=typemax(Int32) || throw(ArgumentError("delaunay3d_exact: point count exceeds Int32 indexing"))

    # Duplicate coordinates are not separate vertices of a geometric triangulation.
    # Reject them explicitly instead of retrying forever while requiring every index.
    length(Set(pts)) == n ||
        throw(ArgumentError("delaunay3d_exact: input points must be pairwise distinct"))

    # Global scalar extent, exact and translation-invariant.
    m = pts[1][1]; Mx = pts[1][1]
    @inbounds for p in pts, k in 1:3
        v = p[k]; v < m && (m = v); v > Mx && (Mx = v)
    end
    D = Mx - m
    # Keep the common path allocation-lean and output-identical: try the historical
    # K=1 construction first.  Only a failed completeness gate pays for the exact
    # affine-conditioning scan used to jump directly to an appropriate larger scale.
    K = big(1)

    lastreason = "no build attempted"
    for growth in 0:_MAX_SUPER_GROWTHS
        for allow_prefilter in (true,false)
            out = _delaunay3d_exact_build(pts,m,D,K,allow_prefilter)
            if out===nothing
                lastreason=allow_prefilter ? "a filtered insertion had an empty cavity" :
                                             "an exact insertion had an empty cavity"
                continue
            end
            ok, reason = _complete_exact_output(pts, out)
            if ok
                dok,nviol=_empty_sphere_exact(pts,out)
                dok && return out
                lastreason="exact empty-circumsphere certification found $nviol violation(s)"
            else
                lastreason=reason
            end
        end
        if growth == 0
            conditioned = _initial_super_scale(pts, D)
            conditioned === nothing && return NTuple{4,Int}[]
            K = max(big(4), conditioned)
        else
            K *= 4
        end
    end
    throw(ErrorException("delaunay3d_exact: failed to construct a complete valid " *
                         "tetrahedralization after $(_MAX_SUPER_GROWTHS + 1) exact " *
                         "super-tet scales; last gate failure: $lastreason"))
end

# One exact Bowyer–Watson build at a specified super-tet scale.  `nothing` requests
# a larger scale (an input point had no cavity in the current finite construction).
function _delaunay3d_exact_build(pts::Vector{NTuple{3,RB}}, m::RB, D::RB, K::BigInt,
                                 allow_prefilter::Bool=true)
    n = length(pts)

    # Bounding super-tetrahedron.  The point cloud lies in [m,m+D]³.  All four
    # vertices must recede as K grows: keeping the negative vertex fixed at m−D was
    # the actual high-aspect failure (a flat real tet's huge circumsphere kept
    # containing that vertex at every positive-only scale).  Here the coordinate
    # planes are x/y/z=m−KD and the fourth plane is x+y+z=3m+8KD, which contains
    # the box for every K≥1.  K=1 is the original well-tested construction for
    # ordinary clouds; larger K scales *all* four vertices without needlessly growing
    # BigInt operands on the common path.
    s1 = (m - K*D,     m - K*D,     m - K*D)
    s2 = (m + 10K*D,   m - K*D,     m - K*D)
    s3 = (m - K*D,     m + 10K*D,   m - K*D)
    s4 = (m - K*D,     m - K*D,     m + 10K*D)
    X = Vector{NTuple{3,RB}}(undef, n + 4)
    @inbounds for i in 1:n; X[i] = pts[i]; end
    X[n+1] = s1; X[n+2] = s2; X[n+3] = s3; X[n+4] = s4

    # Per-axis-shifted Float64 coordinates for the cavity PRE-FILTER only (the exact tests
    # below always use the unshifted `X`). Shifting by the per-axis minimum is exact
    # (rational), so `Xf` holds small, well-conditioned values — this removes the
    # catastrophic cancellation `Float64(b)−Float64(a)` suffers when the model sits far
    # from the origin. If the shifted extent is still enormous (a physically absurd
    # >~1e11-span model where Float64 cannot resolve tet sizes), the pre-filter is disabled
    # (`usePF=false`) and every tet gets the exact `insphere_rat` test — correct, just slower.
    mnx = X[1][1]; mny = X[1][2]; mnz = X[1][3]; mxx = mnx; mxy = mny; mxz = mnz
    @inbounds for i in 1:n+4
        p = X[i]
        p[1] < mnx && (mnx = p[1]); p[1] > mxx && (mxx = p[1])
        p[2] < mny && (mny = p[2]); p[2] > mxy && (mxy = p[2])
        p[3] < mnz && (mnz = p[3]); p[3] > mxz && (mxz = p[3])
    end
    extent = max(mxx-mnx, mxy-mny, mxz-mnz)
    fextent = Float64(extent)
    usePF = allow_prefilter && isfinite(fextent) && fextent < 1.0e11
    Xf = Vector{NTuple{3,Float64}}(undef, n + 4)
    if usePF
        @inbounds for i in 1:n+4
            Xf[i] = (Float64(X[i][1]-mnx), Float64(X[i][2]-mny), Float64(X[i][3]-mnz))
        end
    else
        fill!(Xf, (0.0, 0.0, 0.0))
    end

    tets = NTuple{4,Int}[]
    alive = Bool[]
    csph = NTuple{4,Float64}[]              # parallel: each tet's Float64 circumsphere (pre-filter)
    sizehint!(tets, max(16, 12n)); sizehint!(alive, max(16, 12n)); sizehint!(csph, max(16, 12n))
    function addtet!(t)
        length(tets)<typemax(Int32) ||
            throw(ErrorException("delaunay3d_exact: tetrahedron storage exceeds Int32"))
        push!(tets,t);push!(alive,true)
        push!(csph,usePF ? _fcircum(Xf[t[1]],Xf[t[2]],Xf[t[3]],Xf[t[4]]) :
                             (0.0,0.0,0.0,Inf))
    end
    # initial super-tet, positively oriented
    if orient3_rat(X[n+1], X[n+2], X[n+3], X[n+4], n+1, n+2, n+3, n+4) > 0
        addtet!((n+1, n+2, n+3, n+4))
    else
        addtet!((n+1, n+2, n+4, n+3))
    end

    cav = Int[]
    fcount = Dict{NTuple{3,Int},Int}()
    fverts = Dict{NTuple{3,Int},NTuple{3,Int}}()
    sizehint!(cav, 64); sizehint!(fcount, 128); sizehint!(fverts, 128)
    for i in 1:n
        p = X[i]; pfx=Xf[i][1]; pfy=Xf[i][2]; pfz=Xf[i][3]
        empty!(cav)
        @inbounds for ti in eachindex(tets)
            alive[ti] || continue
            # cheap Float64 circumsphere reject before the expensive exact insphere_rat;
            # only tets that might contain p (borderline/degenerate included) are tested.
            _maybe_in_sphere(csph[ti], pfx, pfy, pfz) || continue
            t = tets[ti]
            insphere_rat(X[t[1]], X[t[2]], X[t[3]], X[t[4]], p, t[1], t[2], t[3], t[4], i) > 0 && push!(cav, ti)
        end
        isempty(cav) && return nothing
        # boundary faces of the cavity: a face incident to exactly ONE cavity tet.
        empty!(fcount); empty!(fverts)
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
            addtet!(nt)
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

# Cheap independent completeness/topology gate for one build.  The Delaunay
# predicate is already exact inside the builder; this gate prevents vacuous success
# after removing super-incident cells and rejects physically unusable SoS-flat cells.
function _complete_exact_output(pts::Vector{NTuple{3,RB}}, out::Vector{NTuple{4,Int}})
    n = length(pts)
    isempty(out) && return (false, "all real tetrahedra were removed")
    length(out)<=typemax(Int32) || return (false,"tetrahedron count exceeds Int32")
    used = falses(n)
    tetkeys = Vector{NTuple{4,Int}}(undef, length(out))
    facerecords = Vector{NTuple{4,Int}}(undef, 4length(out)) # sorted face key + owner tet
    topology=Matrix{Int32}(undef,4,length(out))
    fi = 0
    @inbounds for (ti, t) in enumerate(out)
        for (k,v) in enumerate(t)
            1 <= v <= n || return (false, "tet $ti references vertex $v outside 1:$n")
            used[v] = true
            topology[k,ti]=Int32(v)
        end
        _orient3_det(pts[t[1]], pts[t[2]], pts[t[3]], pts[t[4]]) > 0 ||
            return (false, "tet $ti is physically flat or negatively oriented")
        tetkeys[ti] = _sort4(t[1], t[2], t[3], t[4])
        for f in (_sort3(t[2],t[3],t[4]), _sort3(t[1],t[3],t[4]),
                  _sort3(t[1],t[2],t[4]), _sort3(t[1],t[2],t[3]))
            fi += 1; facerecords[fi] = (f[1],f[2],f[3],ti)
        end
    end
    all(used) || return (false, "output omits $(count(!, used)) of $n input points")
    sort!(tetkeys)
    @inbounds for i in 2:length(tetkeys)
        tetkeys[i] == tetkeys[i-1] && return (false, "duplicate tetrahedron at output cell $i")
    end
    topok,topreason=_tet_topology(topology)
    topok || return (false,"tet complex is not a combinatorial manifold — $topreason")

    sort!(facerecords)
    parent=Int32[i for i in eachindex(out)]
    findroot(i::Int32)=begin
        while parent[i]!=i
            parent[i]=parent[parent[i]];i=parent[i]
        end
        i
    end
    function unite(a::Int32,b::Int32)
        ra=findroot(a);rb=findroot(b);ra!=rb&&(parent[ra]=rb);nothing
    end
    hasboundary = false; i = 1
    @inbounds while i <= length(facerecords)
        j = i + 1
        while j<=length(facerecords) &&
              facerecords[j][1]==facerecords[i][1] &&
              facerecords[j][2]==facerecords[i][2] &&
              facerecords[j][3]==facerecords[i][3]
            j+=1
        end
        incidence = j - i
        incidence <= 2 || return (false, "a triangular face has incidence greater than two")
        if incidence==2
            unite(Int32(facerecords[i][4]),Int32(facerecords[i+1][4]))
        else
            hasboundary=true
            a,b,c,owner=facerecords[i];tet=out[owner]
            opp=first(v for v in tet if v!=a&&v!=b&&v!=c)
            sopp=_orient3_det(pts[a],pts[b],pts[c],pts[opp])
            sopp!=0 || return (false,"a boundary face has a coplanar owner vertex")
            for v in 1:n
                sp=_orient3_det(pts[a],pts[b],pts[c],pts[v])
                sp==0&&continue
                (sp>0)==(sopp>0) || return (false,
                    "a boundary face is not a supporting face of the input convex hull")
            end
        end
        i = j
    end
    hasboundary || return (false, "tet complex has no boundary faces")
    root=findroot(Int32(1))
    @inbounds for t in 2:length(out)
        findroot(Int32(t))==root || return (false,"tet complex is disconnected across faces")
    end
    return (true, "complete")
end

"""
    is_delaunay_exact(pts, tets) -> (ok::Bool, nviol::Int)

Exact empty-circumsphere oracle for a complete tetrahedral complex. Malformed,
duplicate, flat, incomplete, or non-manifold input returns `(false,nviol)` rather
than being accepted vacuously. For a structurally valid complex, `nviol` is the
number of strict empty-circumsphere violations (`insphere_rat`).
"""
function is_delaunay_exact(pts::Vector{NTuple{3,RB}}, tets::Vector{NTuple{4,Int}})
    n = length(pts); nviol = 0
    n>=4 || return (false,1)
    length(Set(pts))==n || return (false,1)
    ok,_=_complete_exact_output(pts,tets)
    ok || return (false,1)
    return _empty_sphere_exact(pts,tets)
end

function _empty_sphere_exact(pts,tets)
    n=length(pts);nviol=0
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
