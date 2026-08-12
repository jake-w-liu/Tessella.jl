"""
    Predicates

Exact, adaptive geometric predicates — the correctness foundation of the whole
mesher (PLAN.md §3, DEVELOPMENT.md "Predicate policy"). Floating-point predicate
inconsistency is *the* root cause of the "overlapping facets / empty volume"
failures we are building Tessella to defeat, so every sign returned here is
provably the true sign of an exact determinant.

Scheme (Shewchuk-style staged precision, simplified to two stages):

1. **Fast filter.** Evaluate the determinant in `Float64` with a certified a
   priori error bound. If the magnitude exceeds the bound, the floating-point
   sign is guaranteed correct and is returned immediately.
2. **Exact fallback.** Otherwise recompute the determinant in exact
   `Rational{BigInt}` arithmetic (every finite `Float64` is a dyadic rational, so
   the conversion is loss-free) and return its exact sign.

The bare predicates (`orient2`, `orient3`, `incircle`, `insphere`) return
`-1`, `0`, or `+1` and are checked against an *independent* exact-rational oracle
on an exhaustive degenerate set (see `test/predicates_test.jl`).

On top of them, the `*_sos` variants take vertex **indices** and apply
Simulation of Simplicity (Edelsbrunner–Mücke): when the exact determinant is
zero, a deterministic symbolic perturbation keyed on the indices breaks the tie,
so a mesher never has to handle a genuine zero. The perturbation is consistent
(the same four points always resolve the same way), which is what keeps the
3-D kernel's decisions globally coherent.
"""
module Predicates

export orient2, orient3, incircle, insphere
export orient2_sos, orient3_sos, incircle_sos, insphere_sos
export orient2_rat, orient3_rat, incircle_rat, insphere_rat

# ── Shewchuk static error bounds for IEEE-754 double ────────────────────────────
# ε is the unit roundoff (half an ulp at 1.0). The A-bounds are the first-stage
# static filters from Shewchuk, "Adaptive Precision Floating-Point Arithmetic and
# Fast Robust Geometric Predicates" (1997), §4. We use the A-stage filter only,
# then jump straight to exact arithmetic — correctness is identical to the full
# adaptive scheme; only near-degenerate speed differs (a Stage-4 optimization).
const EPS = 2.0^-53
const CCWERRBOUND_A  = (3.0  + 16.0  * EPS) * EPS
const O3DERRBOUND_A  = (7.0  + 56.0  * EPS) * EPS
const ICCERRBOUND_A  = (10.0 + 96.0  * EPS) * EPS
const ISPERRBOUND_A  = (16.0 + 224.0 * EPS) * EPS

@inline _isign(x)::Int = x > 0 ? 1 : (x < 0 ? -1 : 0)

# Points are accepted either as tuples/vectors of coordinates. Small helpers keep
# the call sites readable while staying allocation-free and type-stable.
@inline _x(p) = @inbounds Float64(p[1])
@inline _y(p) = @inbounds Float64(p[2])
@inline _z(p) = @inbounds Float64(p[3])

# ════════════════════════════════════════════════════════════════════════════════
# orient2 — sign of the signed area of triangle (a, b, c).
#   > 0 : c is left of a→b (counter-clockwise)
#   < 0 : c is right (clockwise)
#   = 0 : a, b, c collinear
# ════════════════════════════════════════════════════════════════════════════════
"""
    orient2(a, b, c) -> Int

Exact sign of the signed area of the triangle `(a, b, c)`. `+1` if `c` lies to
the left of the directed line `a→b` (counter-clockwise), `-1` to the right,
`0` if the three points are exactly collinear. `a`,`b`,`c` are indexable 2-D
coordinates.
"""
function orient2(a, b, c)::Int
    ax = _x(a); ay = _y(a)
    bx = _x(b); by = _y(b)
    cx = _x(c); cy = _y(c)

    detleft  = (ax - cx) * (by - cy)
    detright = (ay - cy) * (bx - cx)
    det = detleft - detright

    # Certain-sign shortcuts + A-stage filter.
    if detleft > 0.0
        detright <= 0.0 && return 1
        detsum = detleft + detright
    elseif detleft < 0.0
        detright >= 0.0 && return -1
        detsum = -detleft - detright
    else
        return _isign(det)
    end
    errbound = CCWERRBOUND_A * detsum
    (det >= errbound || -det >= errbound) && return _isign(det)

    return _orient2_exact(ax, ay, bx, by, cx, cy)
end

# Every finite Float64 is dyadic (m·2^e), so a coordinate tuple shares a common
# denominator 2^K; scaling by that positive 2^K leaves every exact-predicate SIGN
# unchanged — orient2/orient3 are homogeneous in the coords, and incircle/insphere in
# the coords AND their paraboloid lift with matching total degree. Evaluating over
# these common-denominator BigInt INTEGERS skips Rational GCD/normalization on every
# op: ≈1.8× faster on the exact fallback path, sign-identical (verified against the
# Rational form and the module's exhaustive independent oracle).
@inline function _common_ints(vals::NTuple{N,Float64}) where {N}
    rs = ntuple(i -> Rational{BigInt}(vals[i]), N)
    ks = ntuple(i -> trailing_zeros(denominator(rs[i])), N)
    K  = maximum(ks)
    return ntuple(i -> numerator(rs[i]) << (K - ks[i]), N)
end

function _orient2_exact(ax, ay, bx, by, cx, cy)::Int
    Ax, Ay, Bx, By, Cx, Cy = _common_ints((ax, ay, bx, by, cx, cy))
    det = (Ax - Cx) * (By - Cy) - (Ay - Cy) * (Bx - Cx)
    return det > 0 ? 1 : (det < 0 ? -1 : 0)
end

# ════════════════════════════════════════════════════════════════════════════════
# orient3 — sign of the signed volume of tetrahedron (a, b, c, d).
#   > 0 : d is below the plane of a,b,c as seen so that a,b,c is CCW
#         (i.e. (a,b,c,d) is negatively/positively oriented per convention below)
# Convention (Shewchuk): returns > 0 iff d lies below the oriented plane through
# a,b,c — equivalently the determinant of [a-d; b-d; c-d] is positive.
# ════════════════════════════════════════════════════════════════════════════════
"""
    orient3(a, b, c, d) -> Int

Exact sign of the orientation determinant of the tetrahedron `(a, b, c, d)`:
`sign(det[a-d, b-d, c-d])`. `+1` if `d` is below the plane through `a,b,c`
oriented counter-clockwise (right-hand rule), `-1` above, `0` if the four points
are exactly coplanar.
"""
function orient3(a, b, c, d)::Int
    ax = _x(a); ay = _y(a); az = _z(a)
    bx = _x(b); by = _y(b); bz = _z(b)
    cx = _x(c); cy = _y(c); cz = _z(c)
    dx = _x(d); dy = _y(d); dz = _z(d)

    adx = ax - dx; ady = ay - dy; adz = az - dz
    bdx = bx - dx; bdy = by - dy; bdz = bz - dz
    cdx = cx - dx; cdy = cy - dy; cdz = cz - dz

    bdxcdy = bdx * cdy; cdxbdy = cdx * bdy
    cdxady = cdx * ady; adxcdy = adx * cdy
    adxbdy = adx * bdy; bdxady = bdx * ady

    det = adz * (bdxcdy - cdxbdy) +
          bdz * (cdxady - adxcdy) +
          cdz * (adxbdy - bdxady)

    permanent = (abs(bdxcdy) + abs(cdxbdy)) * abs(adz) +
                (abs(cdxady) + abs(adxcdy)) * abs(bdz) +
                (abs(adxbdy) + abs(bdxady)) * abs(cdz)
    errbound = O3DERRBOUND_A * permanent
    (det > errbound || -det > errbound) && return _isign(det)

    return _orient3_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz)
end

function _orient3_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz)::Int
    Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz = _common_ints((ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz))
    adx = Ax-Dx; ady = Ay-Dy; adz = Az-Dz
    bdx = Bx-Dx; bdy = By-Dy; bdz = Bz-Dz
    cdx = Cx-Dx; cdy = Cy-Dy; cdz = Cz-Dz
    det = adz * (bdx * cdy - cdx * bdy) +
          bdz * (cdx * ady - adx * cdy) +
          cdz * (adx * bdy - bdx * ady)
    return det > 0 ? 1 : (det < 0 ? -1 : 0)
end

# ════════════════════════════════════════════════════════════════════════════════
# incircle — is d inside the circle through a, b, c?
# Returns > 0 iff d is inside the circumcircle of the CCW triangle (a,b,c).
# ════════════════════════════════════════════════════════════════════════════════
"""
    incircle(a, b, c, d) -> Int

Exact in-circle test. With `(a, b, c)` in counter-clockwise order, returns `+1`
if `d` lies strictly inside their circumcircle, `-1` strictly outside, `0`
exactly on the circle. (If `(a,b,c)` is clockwise the sign is reversed, per the
determinant definition.)
"""
function incircle(a, b, c, d)::Int
    ax = _x(a); ay = _y(a)
    bx = _x(b); by = _y(b)
    cx = _x(c); cy = _y(c)
    dx = _x(d); dy = _y(d)

    adx = ax - dx; ady = ay - dy
    bdx = bx - dx; bdy = by - dy
    cdx = cx - dx; cdy = cy - dy

    bdxcdy = bdx * cdy; cdxbdy = cdx * bdy
    cdxady = cdx * ady; adxcdy = adx * cdy
    adxbdy = adx * bdy; bdxady = bdx * ady

    alift = adx * adx + ady * ady
    blift = bdx * bdx + bdy * bdy
    clift = cdx * cdx + cdy * cdy

    det = alift * (bdxcdy - cdxbdy) +
          blift * (cdxady - adxcdy) +
          clift * (adxbdy - bdxady)

    permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift +
                (abs(cdxady) + abs(adxcdy)) * blift +
                (abs(adxbdy) + abs(bdxady)) * clift
    errbound = ICCERRBOUND_A * permanent
    (det > errbound || -det > errbound) && return _isign(det)

    return _incircle_exact(ax, ay, bx, by, cx, cy, dx, dy)
end

function _incircle_exact(ax, ay, bx, by, cx, cy, dx, dy)::Int
    Ax,Ay,Bx,By,Cx,Cy,Dx,Dy = _common_ints((ax,ay,bx,by,cx,cy,dx,dy))
    adx = Ax-Dx; ady = Ay-Dy
    bdx = Bx-Dx; bdy = By-Dy
    cdx = Cx-Dx; cdy = Cy-Dy
    alift = adx * adx + ady * ady
    blift = bdx * bdx + bdy * bdy
    clift = cdx * cdx + cdy * cdy
    det = alift * (bdx * cdy - cdx * bdy) +
          blift * (cdx * ady - adx * cdy) +
          clift * (adx * bdy - bdx * ady)
    return det > 0 ? 1 : (det < 0 ? -1 : 0)
end

# ════════════════════════════════════════════════════════════════════════════════
# insphere — is e inside the sphere through a, b, c, d?
# Returns > 0 iff e is inside the circumsphere of the positively-oriented
# tetrahedron (a,b,c,d) (i.e. orient3(a,b,c,d) > 0).
# ════════════════════════════════════════════════════════════════════════════════
"""
    insphere(a, b, c, d, e) -> Int

Exact in-sphere test. If `orient3(a,b,c,d) > 0`, returns `+1` when `e` lies
strictly inside the circumsphere of `(a,b,c,d)`, `-1` outside, `0` on the sphere.
The raw determinant sign is relative to the orientation of `(a,b,c,d)`; callers
that need the geometric "inside" should combine with `orient3`.
"""
function insphere(a, b, c, d, e)::Int
    ax = _x(a); ay = _y(a); az = _z(a)
    bx = _x(b); by = _y(b); bz = _z(b)
    cx = _x(c); cy = _y(c); cz = _z(c)
    dx = _x(d); dy = _y(d); dz = _z(d)
    ex = _x(e); ey = _y(e); ez = _z(e)

    aex = ax - ex; aey = ay - ey; aez = az - ez
    bex = bx - ex; bey = by - ey; bez = bz - ez
    cex = cx - ex; cey = cy - ey; cez = cz - ez
    dex = dx - ex; dey = dy - ey; dez = dz - ez

    ab = aex * bey - bex * aey
    bc = bex * cey - cex * bey
    cd = cex * dey - dex * cey
    da = dex * aey - aex * dey
    ac = aex * cey - cex * aey
    bd = bex * dey - dex * bey

    abc = aez * bc - bez * ac + cez * ab
    bcd = bez * cd - cez * bd + dez * bc
    cda = cez * da + dez * ac + aez * cd
    dab = dez * ab + aez * bd + bez * da

    alift = aex * aex + aey * aey + aez * aez
    blift = bex * bex + bey * bey + bez * bez
    clift = cex * cex + cey * cey + cez * cez
    dlift = dex * dex + dey * dey + dez * dez

    det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd)

    aezplus = abs(aez); bezplus = abs(bez)
    cezplus = abs(cez); dezplus = abs(dez)
    aexbey = abs(aex * bey); bexaey = abs(bex * aey)
    bexcey = abs(bex * cey); cexbey = abs(cex * bey)
    cexdey = abs(cex * dey); dexcey = abs(dex * cey)
    dexaey = abs(dex * aey); aexdey = abs(aex * dey)
    aexcey = abs(aex * cey); cexaey = abs(cex * aey)
    bexdey = abs(bex * dey); dexbey = abs(dex * bey)
    permanent =
        ((cexdey + dexcey) * bezplus + (dexbey + bexdey) * cezplus +
         (bexcey + cexbey) * dezplus) * alift +
        ((dexaey + aexdey) * cezplus + (aexcey + cexaey) * dezplus +
         (cexdey + dexcey) * aezplus) * blift +
        ((aexbey + bexaey) * dezplus + (bexdey + dexbey) * aezplus +
         (dexaey + aexdey) * bezplus) * clift +
        ((bexcey + cexbey) * aezplus + (cexaey + aexcey) * bezplus +
         (aexbey + bexaey) * cezplus) * dlift
    errbound = ISPERRBOUND_A * permanent
    (det > errbound || -det > errbound) && return _isign(det)

    return _insphere_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez)
end

function _insphere_exact(ax, ay, az, bx, by, bz, cx, cy, cz,
                         dx, dy, dz, ex, ey, ez)::Int
    Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Ex,Ey,Ez =
        _common_ints((ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez))
    aex = Ax-Ex; aey = Ay-Ey; aez = Az-Ez
    bex = Bx-Ex; bey = By-Ey; bez = Bz-Ez
    cex = Cx-Ex; cey = Cy-Ey; cez = Cz-Ez
    dex = Dx-Ex; dey = Dy-Ey; dez = Dz-Ez

    ab = aex * bey - bex * aey
    bc = bex * cey - cex * bey
    cd = cex * dey - dex * cey
    da = dex * aey - aex * dey
    ac = aex * cey - cex * aey
    bd = bex * dey - dex * bey

    abc = aez * bc - bez * ac + cez * ab
    bcd = bez * cd - cez * bd + dez * bc
    cda = cez * da + dez * ac + aez * cd
    dab = dez * ab + aez * bd + bez * da

    alift = aex * aex + aey * aey + aez * aez
    blift = bex * bex + bey * bey + bez * bez
    clift = cex * cex + cey * cey + cez * cez
    dlift = dex * dex + dey * dey + dez * dez

    det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd)
    return _isign(det)
end

# ════════════════════════════════════════════════════════════════════════════════
# Simulation of Simplicity (SoS) wrappers.
#
# Edelsbrunner & Mücke, "Simulation of Simplicity" (1990). Each point i is
# perturbed by (ε^(2^(i·d − 1)), …, ε^(2^(i·d − d))) for ε → 0⁺, with a virtual
# "infinite" index ordering. The upshot for a determinant predicate that comes
# out exactly zero is: expand the perturbed determinant as a polynomial in ε and
# return the sign of the lowest-order nonzero term. For orient/incircle/insphere
# those leading terms are a fixed sequence of lower-degree minors evaluated in a
# canonical order determined by *sorting the point indices*.
#
# We realise SoS operationally: sort the input points by index (tracking the
# permutation parity), evaluate the exact primary determinant, and if it is zero
# fall to the documented minor sequence in sorted order. Sorting makes the
# tie-break independent of the argument order the caller happened to use, which
# is exactly the global consistency SoS guarantees. The returned sign is the
# geometric sign for the *original* argument order (parity reapplied).
# ════════════════════════════════════════════════════════════════════════════════

# Sort helper: returns (perm parity ∈ {+1,-1}, sorted point tuple) for indices.
@inline function _sort2(i, pi, j, pj)
    i < j ? (1, pi, pj) : (-1, pj, pi)
end

"""
    orient2_sos(pa, pb, pc, ia, ib, ic) -> Int (±1, never 0)

Simulation-of-Simplicity orient2: exact orient2, and on an exact collinearity a
deterministic symbolic tie-break keyed on the integer indices `ia,ib,ic`
(all distinct). Never returns 0.
"""
function orient2_sos(pa, pb, pc, ia::Integer, ib::Integer, ic::Integer)::Int
    s = orient2(pa, pb, pc)
    s != 0 && return s
    # Degenerate (collinear/coincident): resolve by the exact +ε Edelsbrunner–Mücke SoS
    # via the SAME canonical evaluator (`_orient_nd_sos_exact`) that orient3_sos/
    # incircle_sos/insphere_sos fall through to, so all four predicates break ties under
    # ONE consistent +ε scheme. (The earlier hand-coded first-order minor sequence agreed
    # on all distinct-point configs but disagreed with the canonical evaluator on
    # coincident-coordinate inputs — measured 144/4416 coincident vs 0/1584 distinct; the
    # canonical route makes the four-predicate consistency hold for every input.)
    return _orient2_sos_exact(pa, pb, pc, ia, ib, ic)
end

# Sort four (index, point) pairs ascending by index; return (parity, sorted).
@inline function _sort4_pts(t)
    v = [t[1], t[2], t[3], t[4]]
    perm = 1
    @inbounds for i in 1:3, j in 1:3-i+1
        if v[j][1] > v[j+1][1]
            v[j], v[j+1] = v[j+1], v[j]
            perm = -perm
        end
    end
    return perm, (v[1], v[2], v[3], v[4])
end

"""
    orient3_sos(pa, pb, pc, pd, ia, ib, ic, id) -> Int (±1, never 0)

Simulation-of-Simplicity orient3. Exact orient3; on an exact coplanarity, a
deterministic tie-break keyed on the distinct indices. Never returns 0.
"""
function orient3_sos(pa, pb, pc, pd,
                     ia::Integer, ib::Integer, ic::Integer, id::Integer)::Int
    s = orient3(pa, pb, pc, pd)
    s != 0 && return s
    # Exact coplanarity: resolve by the canonical +ε Edelsbrunner–Mücke SoS, the SAME
    # evaluator orient2_sos/incircle_sos/insphere_sos use — so all four predicates break
    # ties under ONE consistent scheme. (The earlier hand-coded 8-minor sequence agreed
    # with the canonical evaluator on most configs but MIScomputed the +ε leading term on
    # some coplanar-distinct configs — verified against a finite-ε exact-rational oracle:
    # e.g. points on x=2 gave −1 where the true +ε limit is +1. That path is reachable in
    # the perturb=false kernel, so it is replaced with the always-correct exact SoS.)
    return _orient3_sos_exact(pa, pb, pc, pd, ia, ib, ic, id)
end

"""
    incircle_sos(pa, pb, pc, pd, ia, ib, ic, id) -> Int (±1, never 0)

Simulation-of-Simplicity in-circle. Never returns 0.
"""
function incircle_sos(pa, pb, pc, pd,
                      ia::Integer, ib::Integer, ic::Integer, id::Integer)::Int
    s = incircle(pa, pb, pc, pd)
    s != 0 && return s
    perm, sp = _sort4_pts(((ia, pa), (ib, pb), (ic, pc), (id, pd)))
    b = sp[2][2]; c = sp[3][2]; d = sp[4][2]
    # Leading SoS term of a degenerate in-circle is the orient2 of the three
    # highest-sorted points (lifting the lowest away): the first nonzero of an
    # orient3-style minor sequence. Using orient2 of (b,c,d) resolves the
    # generic cocircular tie; deeper minors cover full symmetry.
    m = orient2(b, c, d); m != 0 && return perm * m
    a = sp[1][2]
    m = orient2(a, c, d); m != 0 && return perm * (-m)
    m = orient2(a, b, d); m != 0 && return perm * m
    m = orient2(a, b, c); m != 0 && return perm * (-m)
    # Fully degenerate (all orient2 minors vanish ⇒ the 4 points are collinear): the
    # fast minors cannot resolve it. Fall to the exact +ε SoS of the paraboloid-lifted
    # points, correct here where the constant `return perm` was only coincidentally so.
    return _incircle_sos_exact(pa, pb, pc, pd, ia, ib, ic, id)
end

"""
    insphere_sos(pa, pb, pc, pd, pe, ia, ib, ic, id, ie) -> Int (±1, never 0)

Simulation-of-Simplicity in-sphere. Never returns 0.
"""
function insphere_sos(pa, pb, pc, pd, pe,
                      ia::Integer, ib::Integer, ic::Integer,
                      id::Integer, ie::Integer)::Int
    s = insphere(pa, pb, pc, pd, pe)
    s != 0 && return s
    # Exact +ε SoS via the paraboloid lift — consistent with orient3/incircle. The old
    # hand-derived orient3-minor sequence evaluated the OPPOSITE (−ε) perturbation (a
    # sign inconsistency with orient3, which the 3-D kernel needs to agree — cavity by
    # insphere, orientation by orient3; verified 120/120 vs the +ε canonical oracle),
    # and its constant `return perm` fallthrough was also wrong for 5 coplanar points.
    # Reached only for exactly cospherical inputs (rare — the kernel pre-perturbs), so
    # the BigInt cost is irrelevant.
    return _insphere_sos_exact(pa, pb, pc, pd, pe, ia, ib, ic, id, ie)
end

# ════════════════════════════════════════════════════════════════════════════════
# Exact +ε Simulation-of-Simplicity leading-term sign (degenerate fallthrough).
#
# When the fast float-minor sequence leaves a *fully* degenerate tie unresolved (all
# minors vanish — e.g. 4 collinear points for orient3), fall to this exact
# evaluation: expand the determinant as a polynomial in ε under the Edelsbrunner–
# Mücke perturbation pᵢ[j] += ε^(2^(rankᵢ·d − j)), ε→0⁺, and return the sign of its
# lowest-order nonzero term. Each matrix entry is a sparse polynomial (ε-exponent ⇒
# Rational{BigInt} coeff), the determinant is the Leibniz sum, exact. This is the SoS
# definition itself — correct for every configuration and consistent with the +ε
# fast-path minors (which handle the common, non-fully-degenerate degeneracies).
# Reached rarely (exact full degeneracy), so the BigInt cost is irrelevant.
# ════════════════════════════════════════════════════════════════════════════════
const _SosPoly = Dict{BigInt, Rational{BigInt}}

@inline function _sos_pcoord(v::Real, r::Int, j::Int, d::Int)
    p = _SosPoly(); p[big(0)] = Rational{BigInt}(v)      # exact: Float64 is dyadic; lift is Rational
    e = big(2)^(r*d - j); p[e] = get(p, e, zero(Rational{BigInt})) + 1
    return p
end
_sos_one() = (p = _SosPoly(); p[big(0)] = one(Rational{BigInt}); p)

function _sos_mul(a::_SosPoly, b::_SosPoly)
    r = _SosPoly()
    for (ea, ca) in a, (eb, cb) in b
        e = ea + eb; r[e] = get(r, e, zero(Rational{BigInt})) + ca*cb
    end
    return r
end
@inline function _sos_add!(acc::_SosPoly, p::_SosPoly)
    for (e, c) in p; acc[e] = get(acc, e, zero(Rational{BigInt})) + c; end
    return acc
end
function _sos_leading_sign(acc::_SosPoly)
    best = nothing
    for (e, c) in acc
        c == 0 && continue
        (best === nothing || e < best) && (best = e)
    end
    best === nothing && return 1        # unreachable for distinct points
    return acc[best] > 0 ? 1 : -1
end

# all permutations of 1:n tagged with sign(σ), n small (≤5 here).
function _sos_perm_rec!(out::Vector{Tuple{Vector{Int},Int}}, p::Vector{Int}, k::Int)
    n = length(p)
    if k > n
        inv = 0
        @inbounds for i in 1:n, j in i+1:n; p[i] > p[j] && (inv += 1); end
        push!(out, (copy(p), isodd(inv) ? -1 : 1)); return
    end
    @inbounds for i in k:n
        p[k], p[i] = p[i], p[k]
        _sos_perm_rec!(out, p, k+1)
        p[k], p[i] = p[i], p[k]
    end
    return
end
function _sos_perms(n::Int)
    out = Tuple{Vector{Int},Int}[]
    _sos_perm_rec!(out, collect(1:n), 1)
    return out
end

# Leibniz determinant (sign of the leading ε-term) of an n×n matrix of `_SosPoly`.
function _sos_det_sign(rows::Vector{Vector{_SosPoly}})
    n = length(rows)
    acc = _SosPoly()
    for (perm, sgn) in _sos_perms(n)
        term = _SosPoly(); term[big(0)] = Rational{BigInt}(sgn)
        @inbounds for r in 1:n; term = _sos_mul(term, rows[r][perm[r]]); end
        _sos_add!(acc, term)
    end
    return _sos_leading_sign(acc)
end

# rank[k] = ascending-index rank (1-based) of the k-th argument point.
@inline _sos_ranks(idx::NTuple{N,Int}) where {N} = invperm(sortperm(collect(idx)))

# Exact +ε SoS sign of the orientation of M=(d+1) points in d-space. `coords[k]` is a
# length-d tuple of exact values (Float64 for real coords, Rational{BigInt} for a
# paraboloid-lift coordinate). Correct for every degenerate configuration.
function _orient_nd_sos_exact(coords, idx::NTuple{M,Int}, d::Int)::Int where {M}
    rk = _sos_ranks(idx)
    rows = Vector{Vector{_SosPoly}}(undef, M)
    @inbounds for k in 1:M
        row = Vector{_SosPoly}(undef, d + 1)
        for j in 1:d; row[j] = _sos_pcoord(coords[k][j], rk[k], j, d); end
        row[d + 1] = _sos_one()
        rows[k] = row
    end
    return _sos_det_sign(rows)
end

"""
    _orient2_sos_exact(pa,pb,pc, ia,ib,ic) -> Int (±1)

Exact +ε SoS sign of the 2-D orientation for the (degenerate) argument points, keyed
on the distinct indices — the 2-D analogue of `_orient3_sos_exact`, so `orient2_sos`
resolves ties by the same canonical evaluator as the higher-dimension predicates.
"""
function _orient2_sos_exact(pa, pb, pc, ia::Integer, ib::Integer, ic::Integer)::Int
    coords = ((_x(pa),_y(pa)), (_x(pb),_y(pb)), (_x(pc),_y(pc)))
    return _orient_nd_sos_exact(coords, (Int(ia),Int(ib),Int(ic)), 2)
end

"""
    _orient3_sos_exact(pa,pb,pc,pd, ia,ib,ic,id) -> Int (±1)

Exact +ε SoS sign of the 3-D orientation for the (degenerate) argument points, keyed
on the distinct indices. Correct for every configuration, collinear included.
"""
function _orient3_sos_exact(pa, pb, pc, pd,
                            ia::Integer, ib::Integer, ic::Integer, id::Integer)::Int
    coords = ((_x(pa),_y(pa),_z(pa)), (_x(pb),_y(pb),_z(pb)),
              (_x(pc),_y(pc),_z(pc)), (_x(pd),_y(pd),_z(pd)))
    return _orient_nd_sos_exact(coords, (Int(ia),Int(ib),Int(ic),Int(id)), 3)
end

# incircle SoS = orientation SoS of the points lifted to the paraboloid (x,y,x²+y²),
# with the lift an INDEPENDENT coordinate perturbed by the EM scheme (the algorithm-
# consistent formulation — verified against a canonical lifted-orientation oracle over
# every degenerate labeling). Reached only for the fully-degenerate (collinear) tie.
@inline _incircle_lift(p) = (xr = Rational{BigInt}(_x(p)); yr = Rational{BigInt}(_y(p)); (xr, yr, xr*xr + yr*yr))
function _incircle_sos_exact(pa, pb, pc, pd,
                             ia::Integer, ib::Integer, ic::Integer, id::Integer)::Int
    coords = (_incircle_lift(pa), _incircle_lift(pb), _incircle_lift(pc), _incircle_lift(pd))
    return _orient_nd_sos_exact(coords, (Int(ia),Int(ib),Int(ic),Int(id)), 3)
end

# insphere SoS = orientation SoS of the points lifted to (x,y,z,x²+y²+z²) in 4-D — the
# same paraboloid-lift construction as incircle, one dimension up. Correct-by-construction
# +ε for EVERY degenerate configuration (cospherical or coplanar).
@inline _insphere_lift(p) = (xr = Rational{BigInt}(_x(p)); yr = Rational{BigInt}(_y(p)); zr = Rational{BigInt}(_z(p));
                             (xr, yr, zr, xr*xr + yr*yr + zr*zr))
function _insphere_sos_exact(pa, pb, pc, pd, pe,
                             ia::Integer, ib::Integer, ic::Integer, id::Integer, ie::Integer)::Int
    coords = (_insphere_lift(pa), _insphere_lift(pb), _insphere_lift(pc), _insphere_lift(pd), _insphere_lift(pe))
    return _orient_nd_sos_exact(coords, (Int(ia),Int(ib),Int(ic),Int(id),Int(ie)), 4)
end

# ════════════════════════════════════════════════════════════════════════════════
# Exact-coordinate predicates (Rational{BigInt}) — for the exact recovery kernel.
#
# The predicates above take Float64 (dyadic) coordinates and exploit dyadic
# representability. The exact kernel places Steiner points at rational positions that
# are NOT Float64-representable (e.g. a crease midpoint of two float vertices), so it
# needs predicates that consume exact `Rational{BigInt}` coordinates directly. These
# evaluate the SAME homogeneous-determinant forms as `test/oracles.jl` (so on dyadic
# input they agree with orient2/3/incircle/insphere_sos), and on a zero (degenerate)
# determinant fall through to the SHARED +ε SoS `_orient_nd_sos_exact` — identical tie-
# break scheme as the Float64 predicates, keeping the two kernels mutually consistent.
# Points are `NTuple`s of `Rational{BigInt}` (or anything the arithmetic promotes).
# ════════════════════════════════════════════════════════════════════════════════
const _RAT = Rational{BigInt}

# exact determinant by recursive Laplace expansion (n ≤ 5 here; O(n!) but exact).
function _rat_det(M::AbstractMatrix{_RAT})::_RAT
    n = size(M, 1)
    n == 1 && return M[1, 1]
    n == 2 && return M[1, 1]*M[2, 2] - M[1, 2]*M[2, 1]
    acc = zero(_RAT); sgn = one(_RAT)
    @inbounds for j in 1:n
        keep = [c for c in 1:n if c != j]
        acc += sgn * M[1, j] * _rat_det(@view M[2:n, keep]); sgn = -sgn
    end
    return acc
end

# Fraction-free (Bareiss) exact integer determinant — O(n³) with exact integer
# division (each cross-multiply is exactly divisible by the previous pivot), no gcd,
# no fraction growth. Much faster than rational Laplace for the kernel's n≤5 matrices.
function _bigint_det(A::Matrix{BigInt})::BigInt
    n = size(A, 1)
    n == 1 && return A[1, 1]
    prev = one(BigInt); sgn = 1
    @inbounds for k in 1:n-1
        if A[k, k] == 0
            piv = 0
            for i in k+1:n; A[i, k] != 0 && (piv = i; break); end
            piv == 0 && return zero(BigInt)                 # singular
            for j in 1:n; A[k, j], A[piv, j] = A[piv, j], A[k, j]; end
            sgn = -sgn
        end
        for i in k+1:n, j in k+1:n
            A[i, j] = (A[i, j]*A[k, k] - A[i, k]*A[k, j]) ÷ prev
        end
        prev = A[k, k]
    end
    return sgn * A[n, n]
end

# Exact SIGN of a rational determinant, fast: clear each row's denominators (scaling a
# row by its positive lcm-of-denominators multiplies det by a POSITIVE number, so the
# sign is preserved), then take the fraction-free integer determinant sign.
function _rat_det_sign(M::AbstractMatrix{_RAT})::Int
    n = size(M, 1)
    A = Matrix{BigInt}(undef, n, n)
    @inbounds for i in 1:n
        d = one(BigInt)
        for j in 1:n; d = lcm(d, denominator(M[i, j])); end
        for j in 1:n; A[i, j] = numerator(M[i, j]) * (d ÷ denominator(M[i, j])); end
    end
    s = _bigint_det(A)
    return s > 0 ? 1 : (s < 0 ? -1 : 0)
end

"""
    orient3_rat(a,b,c,d, ia,ib,ic,id) -> Int (±1)

Exact 3-D orientation sign on `Rational{BigInt}` points (never 0: degeneracies break
by the shared +ε SoS). Matches `orient3`/`orient3_sos` on dyadic input.
"""
function orient3_rat(a, b, c, d, ia::Integer, ib::Integer, ic::Integer, id::Integer)::Int
    M = _RAT[a[1] a[2] a[3] one(_RAT); b[1] b[2] b[3] one(_RAT);
             c[1] c[2] c[3] one(_RAT); d[1] d[2] d[3] one(_RAT)]
    s = _rat_det_sign(M)
    s > 0 && return 1
    s < 0 && return -1
    return _orient_nd_sos_exact((a, b, c, d), (Int(ia), Int(ib), Int(ic), Int(id)), 3)
end

"""
    insphere_rat(a,b,c,d,e, ia,ib,ic,id,ie) -> Int (±1)

Exact in-sphere sign on `Rational{BigInt}` points: `>0` matches `insphere_sos` (for a
positively-oriented `a,b,c,d`, `e` strictly inside the circumsphere). Cospherical ties
break by the shared +ε SoS via the 4-D paraboloid lift.
"""
function insphere_rat(a, b, c, d, e, ia::Integer, ib::Integer, ic::Integer, id::Integer, ie::Integer)::Int
    lift(p) = p[1]*p[1] + p[2]*p[2] + p[3]*p[3]
    M = _RAT[a[1] a[2] a[3] lift(a) one(_RAT); b[1] b[2] b[3] lift(b) one(_RAT);
             c[1] c[2] c[3] lift(c) one(_RAT); d[1] d[2] d[3] lift(d) one(_RAT);
             e[1] e[2] e[3] lift(e) one(_RAT)]
    s = _rat_det_sign(M)
    s > 0 && return 1
    s < 0 && return -1
    lp(p) = (p[1], p[2], p[3], lift(p))
    return _orient_nd_sos_exact((lp(a), lp(b), lp(c), lp(d), lp(e)),
                                (Int(ia), Int(ib), Int(ic), Int(id), Int(ie)), 4)
end

"""
    orient2_rat(a,b,c, ia,ib,ic) -> Int (±1)

Exact 2-D orientation sign on `Rational{BigInt}` points; degeneracies break by SoS.
"""
function orient2_rat(a, b, c, ia::Integer, ib::Integer, ic::Integer)::Int
    M = _RAT[a[1] a[2] one(_RAT); b[1] b[2] one(_RAT); c[1] c[2] one(_RAT)]
    s = _rat_det_sign(M)
    s > 0 && return 1
    s < 0 && return -1
    return _orient_nd_sos_exact((a, b, c), (Int(ia), Int(ib), Int(ic)), 2)
end

"""
    incircle_rat(a,b,c,d, ia,ib,ic,id) -> Int (±1)

Exact in-circle sign on `Rational{BigInt}` 2-D points (paraboloid-lift SoS for ties).
"""
function incircle_rat(a, b, c, d, ia::Integer, ib::Integer, ic::Integer, id::Integer)::Int
    lift(p) = p[1]*p[1] + p[2]*p[2]
    M = _RAT[a[1] a[2] lift(a) one(_RAT); b[1] b[2] lift(b) one(_RAT);
             c[1] c[2] lift(c) one(_RAT); d[1] d[2] lift(d) one(_RAT)]
    s = _rat_det_sign(M)
    s > 0 && return 1
    s < 0 && return -1
    lp(p) = (p[1], p[2], lift(p))
    return _orient_nd_sos_exact((lp(a), lp(b), lp(c), lp(d)),
                                (Int(ia), Int(ib), Int(ic), Int(id)), 3)
end

end # module Predicates
