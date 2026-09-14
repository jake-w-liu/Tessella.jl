# ════════════════════════════════════════════════════════════════════════════════
# Adaptive expansion stages B–D for orient2/orient3 (Shewchuk 1997, §4.3–4.4).
#
# The A-stage filters in Predicates.jl decide every well-conditioned case.  The
# remaining near-degenerate and exactly degenerate cases (collinear or coplanar
# input is routine on faceted geometry) previously went straight to
# Rational{BigInt}, allocating on every call.  These stages evaluate the
# determinant exactly with error-free Float64 transformations and
# nonoverlapping expansions held in per-thread scratch, so the common exact
# decisions allocate nothing.  A magnitude guard routes inputs whose products
# could overflow or whose expansion components could underflow to the retained
# BigInt path, because the expansion theorems assume neither.
# ════════════════════════════════════════════════════════════════════════════════

const RESULTERRBOUND = (3.0 + 8.0 * EPS) * EPS
const CCWERRBOUND_B  = (2.0 + 12.0 * EPS) * EPS
const CCWERRBOUND_C  = (9.0 + 64.0 * EPS) * EPS * EPS
const O3DERRBOUND_B  = (2.0 + 12.0 * EPS) * EPS
const O3DERRBOUND_C  = (26.0 + 288.0 * EPS) * EPS * EPS

# Every nonzero difference and difference tail must lie in a magnitude band so
# that no product overflows and no expansion component underflows into the
# subnormal range.  orient3: degree-3 products stay below 2^900 and components,
# at most ~2^-320 relative to a product, stay above 2^-1010.  orient2: degree-2
# products stay below 2^960 and components (~2^-220 relative) above 2^-980.
const _ADAPTIVE3_MIN_MAGNITUDE = 2.0^-230
const _ADAPTIVE3_MAX_MAGNITUDE = 2.0^300
const _ADAPTIVE2_MIN_MAGNITUDE = 2.0^-380
const _ADAPTIVE2_MAX_MAGNITUDE = 2.0^480

@inline function _fast_two_sum(a::Float64, b::Float64)
    x = a + b
    return x, b - (x - a)
end
@inline function _two_sum_pair(a::Float64, b::Float64)
    x = a + b
    bvirt = x - a
    avirt = x - bvirt
    return x, (a - avirt) + (b - bvirt)
end
@inline function _two_diff_pair(a::Float64, b::Float64)
    x = a - b
    bvirt = a - x
    avirt = x + bvirt
    return x, (a - avirt) + (bvirt - b)
end
@inline function _two_diff_tail(a::Float64, b::Float64, x::Float64)
    bvirt = a - x
    avirt = x + bvirt
    return (a - avirt) + (bvirt - b)
end
@inline function _two_product_pair(a::Float64, b::Float64)
    x = a * b
    return x, fma(a, b, -x)
end
# (a1,a0) - b  ->  (x2,x1,x0)
@inline function _two_one_diff(a1::Float64, a0::Float64, b::Float64)
    i, x0 = _two_diff_pair(a0, b)
    x2, x1 = _two_sum_pair(a1, i)
    return x2, x1, x0
end
# (a1,a0) - (b1,b0)  ->  (x3,x2,x1,x0)
@inline function _two_two_diff(a1::Float64, a0::Float64, b1::Float64, b0::Float64)
    j, r0, x0 = _two_one_diff(a1, a0, b0)
    x3, x2, x1 = _two_one_diff(j, r0, b1)
    return x3, x2, x1, x0
end
# (a1,a0) * b  ->  (x3,x2,x1,x0)
@inline function _two_one_product(a1::Float64, a0::Float64, b::Float64)
    i, x0 = _two_product_pair(a0, b)
    j, r0 = _two_product_pair(a1, b)
    k, x1 = _two_sum_pair(i, r0)
    x3, x2 = _fast_two_sum(j, k)
    return x3, x2, x1, x0
end

@inline function _magnitude_ok3(v::Float64)
    a = abs(v)
    return a == 0.0 || (_ADAPTIVE3_MIN_MAGNITUDE <= a <= _ADAPTIVE3_MAX_MAGNITUDE)
end
@inline function _magnitude_ok2(v::Float64)
    a = abs(v)
    return a == 0.0 || (_ADAPTIVE2_MIN_MAGNITUDE <= a <= _ADAPTIVE2_MAX_MAGNITUDE)
end

# Per-thread scratch: expansions are written in place, never returned.
struct _AdaptiveScratch
    bc::Vector{Float64}; ca::Vector{Float64}; ab::Vector{Float64}
    adet::Vector{Float64}; bdet::Vector{Float64}; cdet::Vector{Float64}
    abdet::Vector{Float64}
    fin1::Vector{Float64}; fin2::Vector{Float64}
    at_b::Vector{Float64}; at_c::Vector{Float64}
    bt_c::Vector{Float64}; bt_a::Vector{Float64}
    ct_a::Vector{Float64}; ct_b::Vector{Float64}
    bct::Vector{Float64}; cat::Vector{Float64}; abt::Vector{Float64}
    u::Vector{Float64}; v::Vector{Float64}; w::Vector{Float64}
    B::Vector{Float64}; C1::Vector{Float64}; C2::Vector{Float64}; D::Vector{Float64}
end
_AdaptiveScratch() = _AdaptiveScratch(
    Vector{Float64}(undef, 4), Vector{Float64}(undef, 4), Vector{Float64}(undef, 4),
    Vector{Float64}(undef, 8), Vector{Float64}(undef, 8), Vector{Float64}(undef, 8),
    Vector{Float64}(undef, 16),
    Vector{Float64}(undef, 192), Vector{Float64}(undef, 192),
    Vector{Float64}(undef, 4), Vector{Float64}(undef, 4),
    Vector{Float64}(undef, 4), Vector{Float64}(undef, 4),
    Vector{Float64}(undef, 4), Vector{Float64}(undef, 4),
    Vector{Float64}(undef, 8), Vector{Float64}(undef, 8), Vector{Float64}(undef, 8),
    Vector{Float64}(undef, 4), Vector{Float64}(undef, 12), Vector{Float64}(undef, 16),
    Vector{Float64}(undef, 4), Vector{Float64}(undef, 8), Vector{Float64}(undef, 12),
    Vector{Float64}(undef, 16))

const _ADAPTIVE_SCRATCH = Vector{_AdaptiveScratch}()
const _ADAPTIVE_SCRATCH_LOCK = ReentrantLock()

@inline function _adaptive_scratch()
    tid = Threads.threadid()
    scratch = _ADAPTIVE_SCRATCH
    if tid > length(scratch)
        lock(_ADAPTIVE_SCRATCH_LOCK) do
            while length(scratch) < tid
                push!(scratch, _AdaptiveScratch())
            end
        end
    end
    @inbounds return scratch[tid]
end

# h = e * b with zero elimination; returns the component count.
function _scale_expansion_zeroelim!(elen::Int, e::Vector{Float64}, b::Float64,
                                    h::Vector{Float64})
    @inbounds begin
        Q, hh = _two_product_pair(e[1], b)
        hindex = 0
        if hh != 0.0
            hindex += 1; h[hindex] = hh
        end
        for eindex in 2:elen
            product1, product0 = _two_product_pair(e[eindex], b)
            sum, hh = _two_sum_pair(Q, product0)
            if hh != 0.0
                hindex += 1; h[hindex] = hh
            end
            Q, hh = _fast_two_sum(product1, sum)
            if hh != 0.0
                hindex += 1; h[hindex] = hh
            end
        end
        if Q != 0.0 || hindex == 0
            hindex += 1; h[hindex] = Q
        end
    end
    return hindex
end

# h = e + f with zero elimination; returns the component count.
function _fast_expansion_sum_zeroelim!(elen::Int, e::Vector{Float64},
                                       flen::Int, f::Vector{Float64},
                                       h::Vector{Float64})
    @inbounds begin
        enow = e[1]; fnow = f[1]
        eindex = 1; findex = 1
        if (fnow > enow) == (fnow > -enow)
            Q = enow; eindex += 1
            eindex <= elen && (enow = e[eindex])
        else
            Q = fnow; findex += 1
            findex <= flen && (fnow = f[findex])
        end
        hindex = 0
        if eindex <= elen && findex <= flen
            if (fnow > enow) == (fnow > -enow)
                Q, hh = _fast_two_sum(enow, Q)
                eindex += 1
                eindex <= elen && (enow = e[eindex])
            else
                Q, hh = _fast_two_sum(fnow, Q)
                findex += 1
                findex <= flen && (fnow = f[findex])
            end
            if hh != 0.0
                hindex += 1; h[hindex] = hh
            end
            while eindex <= elen && findex <= flen
                if (fnow > enow) == (fnow > -enow)
                    Q, hh = _two_sum_pair(Q, enow)
                    eindex += 1
                    eindex <= elen && (enow = e[eindex])
                else
                    Q, hh = _two_sum_pair(Q, fnow)
                    findex += 1
                    findex <= flen && (fnow = f[findex])
                end
                if hh != 0.0
                    hindex += 1; h[hindex] = hh
                end
            end
        end
        while eindex <= elen
            Q, hh = _two_sum_pair(Q, enow)
            eindex += 1
            eindex <= elen && (enow = e[eindex])
            if hh != 0.0
                hindex += 1; h[hindex] = hh
            end
        end
        while findex <= flen
            Q, hh = _two_sum_pair(Q, fnow)
            findex += 1
            findex <= flen && (fnow = f[findex])
            if hh != 0.0
                hindex += 1; h[hindex] = hh
            end
        end
        if Q != 0.0 || hindex == 0
            hindex += 1; h[hindex] = Q
        end
    end
    return hindex
end

@inline function _estimate(elen::Int, e::Vector{Float64})
    Q = 0.0
    @inbounds for i in 1:elen
        Q += e[i]
    end
    return Q
end

@inline function _store4!(h::Vector{Float64}, x3, x2, x1, x0)
    @inbounds h[1] = x0; h[2] = x1; h[3] = x2; h[4] = x3
    return nothing
end

# Sign of the exact determinant, or `nothing` when the magnitude guard sends
# the decision to the BigInt path.
function _orient2_adaptive(ax, ay, bx, by, cx, cy, detsum::Float64)
    acx = ax - cx; bcx = bx - cx; acy = ay - cy; bcy = by - cy
    (_magnitude_ok2(acx) && _magnitude_ok2(bcx) && _magnitude_ok2(acy) &&
     _magnitude_ok2(bcy)) || return nothing
    s = _adaptive_scratch()
    detleft, detlefttail = _two_product_pair(acx, bcy)
    detright, detrighttail = _two_product_pair(acy, bcx)
    B3, B2, B1, B0 = _two_two_diff(detleft, detlefttail, detright, detrighttail)
    _store4!(s.B, B3, B2, B1, B0)
    det = _estimate(4, s.B)
    errbound = CCWERRBOUND_B * detsum
    (det >= errbound || -det >= errbound) && return _isign(det)

    acxtail = _two_diff_tail(ax, cx, acx)
    bcxtail = _two_diff_tail(bx, cx, bcx)
    acytail = _two_diff_tail(ay, cy, acy)
    bcytail = _two_diff_tail(by, cy, bcy)
    if acxtail == 0.0 && acytail == 0.0 && bcxtail == 0.0 && bcytail == 0.0
        return _isign(det)
    end
    (_magnitude_ok2(acxtail) && _magnitude_ok2(bcxtail) && _magnitude_ok2(acytail) &&
     _magnitude_ok2(bcytail)) || return nothing

    errbound = CCWERRBOUND_C * detsum + RESULTERRBOUND * abs(det)
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail)
    (det >= errbound || -det >= errbound) && return _isign(det)

    s1, s0 = _two_product_pair(acxtail, bcy)
    t1, t0 = _two_product_pair(acytail, bcx)
    u3, u2, u1, u0 = _two_two_diff(s1, s0, t1, t0)
    _store4!(s.u, u3, u2, u1, u0)
    c1len = _fast_expansion_sum_zeroelim!(4, s.B, 4, s.u, s.C1)

    s1, s0 = _two_product_pair(acx, bcytail)
    t1, t0 = _two_product_pair(acy, bcxtail)
    u3, u2, u1, u0 = _two_two_diff(s1, s0, t1, t0)
    _store4!(s.u, u3, u2, u1, u0)
    c2len = _fast_expansion_sum_zeroelim!(c1len, s.C1, 4, s.u, s.C2)

    s1, s0 = _two_product_pair(acxtail, bcytail)
    t1, t0 = _two_product_pair(acytail, bcxtail)
    u3, u2, u1, u0 = _two_two_diff(s1, s0, t1, t0)
    _store4!(s.u, u3, u2, u1, u0)
    dlen = _fast_expansion_sum_zeroelim!(c2len, s.C2, 4, s.u, s.D)
    @inbounds return _isign(s.D[dlen])
end

function _orient3_adaptive(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz,
                           permanent::Float64)
    adx = ax - dx; bdx = bx - dx; cdx = cx - dx
    ady = ay - dy; bdy = by - dy; cdy = cy - dy
    adz = az - dz; bdz = bz - dz; cdz = cz - dz
    (_magnitude_ok3(adx) && _magnitude_ok3(bdx) && _magnitude_ok3(cdx) &&
     _magnitude_ok3(ady) && _magnitude_ok3(bdy) && _magnitude_ok3(cdy) &&
     _magnitude_ok3(adz) && _magnitude_ok3(bdz) && _magnitude_ok3(cdz)) ||
        return nothing
    s = _adaptive_scratch()

    bdxcdy1, bdxcdy0 = _two_product_pair(bdx, cdy)
    cdxbdy1, cdxbdy0 = _two_product_pair(cdx, bdy)
    x3, x2, x1, x0 = _two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0)
    _store4!(s.bc, x3, x2, x1, x0)
    alen = _scale_expansion_zeroelim!(4, s.bc, adz, s.adet)

    cdxady1, cdxady0 = _two_product_pair(cdx, ady)
    adxcdy1, adxcdy0 = _two_product_pair(adx, cdy)
    x3, x2, x1, x0 = _two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0)
    _store4!(s.ca, x3, x2, x1, x0)
    blen = _scale_expansion_zeroelim!(4, s.ca, bdz, s.bdet)

    adxbdy1, adxbdy0 = _two_product_pair(adx, bdy)
    bdxady1, bdxady0 = _two_product_pair(bdx, ady)
    x3, x2, x1, x0 = _two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0)
    _store4!(s.ab, x3, x2, x1, x0)
    clen = _scale_expansion_zeroelim!(4, s.ab, cdz, s.cdet)

    ablen = _fast_expansion_sum_zeroelim!(alen, s.adet, blen, s.bdet, s.abdet)
    finlength = _fast_expansion_sum_zeroelim!(ablen, s.abdet, clen, s.cdet, s.fin1)

    det = _estimate(finlength, s.fin1)
    errbound = O3DERRBOUND_B * permanent
    (det >= errbound || -det >= errbound) && return _isign(det)

    adxtail = _two_diff_tail(ax, dx, adx)
    bdxtail = _two_diff_tail(bx, dx, bdx)
    cdxtail = _two_diff_tail(cx, dx, cdx)
    adytail = _two_diff_tail(ay, dy, ady)
    bdytail = _two_diff_tail(by, dy, bdy)
    cdytail = _two_diff_tail(cy, dy, cdy)
    adztail = _two_diff_tail(az, dz, adz)
    bdztail = _two_diff_tail(bz, dz, bdz)
    cdztail = _two_diff_tail(cz, dz, cdz)

    if adxtail == 0.0 && bdxtail == 0.0 && cdxtail == 0.0 &&
       adytail == 0.0 && bdytail == 0.0 && cdytail == 0.0 &&
       adztail == 0.0 && bdztail == 0.0 && cdztail == 0.0
        return _isign(det)
    end
    (_magnitude_ok3(adxtail) && _magnitude_ok3(bdxtail) && _magnitude_ok3(cdxtail) &&
     _magnitude_ok3(adytail) && _magnitude_ok3(bdytail) && _magnitude_ok3(cdytail) &&
     _magnitude_ok3(adztail) && _magnitude_ok3(bdztail) && _magnitude_ok3(cdztail)) ||
        return nothing

    errbound = O3DERRBOUND_C * permanent + RESULTERRBOUND * abs(det)
    det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) +
            adztail * (bdx * cdy - cdx * bdy)) +
           (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) +
            bdztail * (cdx * ady - adx * cdy)) +
           (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) +
            cdztail * (adx * bdy - bdx * ady))
    (det >= errbound || -det >= errbound) && return _isign(det)

    finnow = s.fin1; finother = s.fin2

    if adxtail == 0.0
        if adytail == 0.0
            @inbounds s.at_b[1] = 0.0; at_blen = 1
            @inbounds s.at_c[1] = 0.0; at_clen = 1
        else
            negate = -adytail
            at_blarge, tail = _two_product_pair(negate, bdx)
            @inbounds s.at_b[1] = tail; s.at_b[2] = at_blarge; at_blen = 2
            at_clarge, tail = _two_product_pair(adytail, cdx)
            @inbounds s.at_c[1] = tail; s.at_c[2] = at_clarge; at_clen = 2
        end
    else
        if adytail == 0.0
            at_blarge, tail = _two_product_pair(adxtail, bdy)
            @inbounds s.at_b[1] = tail; s.at_b[2] = at_blarge; at_blen = 2
            negate = -adxtail
            at_clarge, tail = _two_product_pair(negate, cdy)
            @inbounds s.at_c[1] = tail; s.at_c[2] = at_clarge; at_clen = 2
        else
            adxt_bdy1, adxt_bdy0 = _two_product_pair(adxtail, bdy)
            adyt_bdx1, adyt_bdx0 = _two_product_pair(adytail, bdx)
            x3, x2, x1, x0 = _two_two_diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0)
            _store4!(s.at_b, x3, x2, x1, x0); at_blen = 4
            adyt_cdx1, adyt_cdx0 = _two_product_pair(adytail, cdx)
            adxt_cdy1, adxt_cdy0 = _two_product_pair(adxtail, cdy)
            x3, x2, x1, x0 = _two_two_diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0)
            _store4!(s.at_c, x3, x2, x1, x0); at_clen = 4
        end
    end
    if bdxtail == 0.0
        if bdytail == 0.0
            @inbounds s.bt_c[1] = 0.0; bt_clen = 1
            @inbounds s.bt_a[1] = 0.0; bt_alen = 1
        else
            negate = -bdytail
            bt_clarge, tail = _two_product_pair(negate, cdx)
            @inbounds s.bt_c[1] = tail; s.bt_c[2] = bt_clarge; bt_clen = 2
            bt_alarge, tail = _two_product_pair(bdytail, adx)
            @inbounds s.bt_a[1] = tail; s.bt_a[2] = bt_alarge; bt_alen = 2
        end
    else
        if bdytail == 0.0
            bt_clarge, tail = _two_product_pair(bdxtail, cdy)
            @inbounds s.bt_c[1] = tail; s.bt_c[2] = bt_clarge; bt_clen = 2
            negate = -bdxtail
            bt_alarge, tail = _two_product_pair(negate, ady)
            @inbounds s.bt_a[1] = tail; s.bt_a[2] = bt_alarge; bt_alen = 2
        else
            bdxt_cdy1, bdxt_cdy0 = _two_product_pair(bdxtail, cdy)
            bdyt_cdx1, bdyt_cdx0 = _two_product_pair(bdytail, cdx)
            x3, x2, x1, x0 = _two_two_diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0)
            _store4!(s.bt_c, x3, x2, x1, x0); bt_clen = 4
            bdyt_adx1, bdyt_adx0 = _two_product_pair(bdytail, adx)
            bdxt_ady1, bdxt_ady0 = _two_product_pair(bdxtail, ady)
            x3, x2, x1, x0 = _two_two_diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0)
            _store4!(s.bt_a, x3, x2, x1, x0); bt_alen = 4
        end
    end
    if cdxtail == 0.0
        if cdytail == 0.0
            @inbounds s.ct_a[1] = 0.0; ct_alen = 1
            @inbounds s.ct_b[1] = 0.0; ct_blen = 1
        else
            negate = -cdytail
            ct_alarge, tail = _two_product_pair(negate, adx)
            @inbounds s.ct_a[1] = tail; s.ct_a[2] = ct_alarge; ct_alen = 2
            ct_blarge, tail = _two_product_pair(cdytail, bdx)
            @inbounds s.ct_b[1] = tail; s.ct_b[2] = ct_blarge; ct_blen = 2
        end
    else
        if cdytail == 0.0
            ct_alarge, tail = _two_product_pair(cdxtail, ady)
            @inbounds s.ct_a[1] = tail; s.ct_a[2] = ct_alarge; ct_alen = 2
            negate = -cdxtail
            ct_blarge, tail = _two_product_pair(negate, bdy)
            @inbounds s.ct_b[1] = tail; s.ct_b[2] = ct_blarge; ct_blen = 2
        else
            cdxt_ady1, cdxt_ady0 = _two_product_pair(cdxtail, ady)
            cdyt_adx1, cdyt_adx0 = _two_product_pair(cdytail, adx)
            x3, x2, x1, x0 = _two_two_diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0)
            _store4!(s.ct_a, x3, x2, x1, x0); ct_alen = 4
            cdyt_bdx1, cdyt_bdx0 = _two_product_pair(cdytail, bdx)
            cdxt_bdy1, cdxt_bdy0 = _two_product_pair(cdxtail, bdy)
            x3, x2, x1, x0 = _two_two_diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0)
            _store4!(s.ct_b, x3, x2, x1, x0); ct_blen = 4
        end
    end

    bctlen = _fast_expansion_sum_zeroelim!(bt_clen, s.bt_c, ct_blen, s.ct_b, s.bct)
    wlength = _scale_expansion_zeroelim!(bctlen, s.bct, adz, s.w)
    finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, wlength, s.w, finother)
    finnow, finother = finother, finnow

    catlen = _fast_expansion_sum_zeroelim!(ct_alen, s.ct_a, at_clen, s.at_c, s.cat)
    wlength = _scale_expansion_zeroelim!(catlen, s.cat, bdz, s.w)
    finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, wlength, s.w, finother)
    finnow, finother = finother, finnow

    abtlen = _fast_expansion_sum_zeroelim!(at_blen, s.at_b, bt_alen, s.bt_a, s.abt)
    wlength = _scale_expansion_zeroelim!(abtlen, s.abt, cdz, s.w)
    finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, wlength, s.w, finother)
    finnow, finother = finother, finnow

    if adztail != 0.0
        vlength = _scale_expansion_zeroelim!(4, s.bc, adztail, s.v)
        finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, vlength, s.v, finother)
        finnow, finother = finother, finnow
    end
    if bdztail != 0.0
        vlength = _scale_expansion_zeroelim!(4, s.ca, bdztail, s.v)
        finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, vlength, s.v, finother)
        finnow, finother = finother, finnow
    end
    if cdztail != 0.0
        vlength = _scale_expansion_zeroelim!(4, s.ab, cdztail, s.v)
        finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, vlength, s.v, finother)
        finnow, finother = finother, finnow
    end

    if adxtail != 0.0
        if bdytail != 0.0
            adxt_bdyt1, adxt_bdyt0 = _two_product_pair(adxtail, bdytail)
            x3, x2, x1, x0 = _two_one_product(adxt_bdyt1, adxt_bdyt0, cdz)
            _store4!(s.u, x3, x2, x1, x0)
            finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
            finnow, finother = finother, finnow
            if cdztail != 0.0
                x3, x2, x1, x0 = _two_one_product(adxt_bdyt1, adxt_bdyt0, cdztail)
                _store4!(s.u, x3, x2, x1, x0)
                finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
                finnow, finother = finother, finnow
            end
        end
        if cdytail != 0.0
            negate = -adxtail
            adxt_cdyt1, adxt_cdyt0 = _two_product_pair(negate, cdytail)
            x3, x2, x1, x0 = _two_one_product(adxt_cdyt1, adxt_cdyt0, bdz)
            _store4!(s.u, x3, x2, x1, x0)
            finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
            finnow, finother = finother, finnow
            if bdztail != 0.0
                x3, x2, x1, x0 = _two_one_product(adxt_cdyt1, adxt_cdyt0, bdztail)
                _store4!(s.u, x3, x2, x1, x0)
                finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
                finnow, finother = finother, finnow
            end
        end
    end
    if bdxtail != 0.0
        if cdytail != 0.0
            bdxt_cdyt1, bdxt_cdyt0 = _two_product_pair(bdxtail, cdytail)
            x3, x2, x1, x0 = _two_one_product(bdxt_cdyt1, bdxt_cdyt0, adz)
            _store4!(s.u, x3, x2, x1, x0)
            finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
            finnow, finother = finother, finnow
            if adztail != 0.0
                x3, x2, x1, x0 = _two_one_product(bdxt_cdyt1, bdxt_cdyt0, adztail)
                _store4!(s.u, x3, x2, x1, x0)
                finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
                finnow, finother = finother, finnow
            end
        end
        if adytail != 0.0
            negate = -bdxtail
            bdxt_adyt1, bdxt_adyt0 = _two_product_pair(negate, adytail)
            x3, x2, x1, x0 = _two_one_product(bdxt_adyt1, bdxt_adyt0, cdz)
            _store4!(s.u, x3, x2, x1, x0)
            finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
            finnow, finother = finother, finnow
            if cdztail != 0.0
                x3, x2, x1, x0 = _two_one_product(bdxt_adyt1, bdxt_adyt0, cdztail)
                _store4!(s.u, x3, x2, x1, x0)
                finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
                finnow, finother = finother, finnow
            end
        end
    end
    if cdxtail != 0.0
        if adytail != 0.0
            cdxt_adyt1, cdxt_adyt0 = _two_product_pair(cdxtail, adytail)
            x3, x2, x1, x0 = _two_one_product(cdxt_adyt1, cdxt_adyt0, bdz)
            _store4!(s.u, x3, x2, x1, x0)
            finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
            finnow, finother = finother, finnow
            if bdztail != 0.0
                x3, x2, x1, x0 = _two_one_product(cdxt_adyt1, cdxt_adyt0, bdztail)
                _store4!(s.u, x3, x2, x1, x0)
                finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
                finnow, finother = finother, finnow
            end
        end
        if bdytail != 0.0
            negate = -cdxtail
            cdxt_bdyt1, cdxt_bdyt0 = _two_product_pair(negate, bdytail)
            x3, x2, x1, x0 = _two_one_product(cdxt_bdyt1, cdxt_bdyt0, adz)
            _store4!(s.u, x3, x2, x1, x0)
            finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
            finnow, finother = finother, finnow
            if adztail != 0.0
                x3, x2, x1, x0 = _two_one_product(cdxt_bdyt1, cdxt_bdyt0, adztail)
                _store4!(s.u, x3, x2, x1, x0)
                finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, 4, s.u, finother)
                finnow, finother = finother, finnow
            end
        end
    end

    if adztail != 0.0
        wlength = _scale_expansion_zeroelim!(bctlen, s.bct, adztail, s.w)
        finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, wlength, s.w, finother)
        finnow, finother = finother, finnow
    end
    if bdztail != 0.0
        wlength = _scale_expansion_zeroelim!(catlen, s.cat, bdztail, s.w)
        finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, wlength, s.w, finother)
        finnow, finother = finother, finnow
    end
    if cdztail != 0.0
        wlength = _scale_expansion_zeroelim!(abtlen, s.abt, cdztail, s.w)
        finlength = _fast_expansion_sum_zeroelim!(finlength, finnow, wlength, s.w, finother)
        finnow, finother = finother, finnow
    end

    @inbounds return _isign(finnow[finlength])
end
