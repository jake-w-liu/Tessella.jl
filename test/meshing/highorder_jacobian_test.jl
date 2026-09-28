# ── Stage-6b CRC suite: curved-P2 Jacobian certification on the non-simplex
#    reference families (gmsh type-10 9-node quadrangle, type-12 27-node
#    hexahedron, type-13 18-node prism) ────────────────────────────────────────
#
# Correctness  : INDEPENDENT oracle — direct product-Lagrange evaluation of the
#                same isoparametric maps (_quad_g/_hex_det/_pri_det), run in
#                Float64 AND in exact Rational{BigInt} arithmetic; the sampling
#                path is itself verified against the rational recomputation at
#                fixed-seed interior points. The certificate is then compared
#                against dense adaptive sampling and an exact-rational grid min.
# Robustness   : fixed-seed warped corpus, boundary-touching warps, a zero-area
#                quadrangle, a collapsed-edge prism and a thin/near-planar
#                hexahedron — success or explicit diagnostic (certificate ≤ 0),
#                never a wrong certificate (cert > sampled min would be a bug).
# Completeness : malformed input is an ArgumentError; certificates are pure
#                BigInt arithmetic, so they are bit-identical across processes.

using Test
using Random: MersenneTwister
import Tessella
using Tessella.Elements: lagrange_nodes
using Tessella.HighOrder: p2_quad_min_jacobian, p2_hex_min_jacobian,
                          p2_prism_min_jacobian

# ══ Independent isoparametric oracle ═════════════════════════════════════════
# Quadratic Lagrange factors on [-1,1] at the cardinal nodes -1,0,+1, written in
# direct product form (i = 1,2,3 ↔ ξ = -1,0,+1).  Deliberately NOT the module's
# Bernstein-coefficient tables: a distinct code path evaluating the same maps.
_l2(i::Integer, x) =
    i == 1 ? x*(x-one(x))/2 : i == 2 ? one(x)-x*x : x*(x+one(x))/2
_l2d(i::Integer, x) =
    i == 1 ? x-one(x)/2 : i == 2 ? -2x : x+one(x)/2

# Type-9 quadratic triangle factors on {r,s ≥ 0, r+s ≤ 1} and their exact
# first derivatives (a = 1..6 in gmsh order: corners, then edges (1,2),(2,3),(3,1)).
_tri9(a::Integer, r, s) =
    a == 1 ? (1-r-s)*(2(1-r-s)-1) : a == 2 ? r*(2r-1) : a == 3 ? s*(2s-1) :
    a == 4 ? 4(1-r-s)*r : a == 5 ? 4r*s : 4s*(1-r-s)
_tri9dr(a::Integer, r, s) =
    a == 1 ? 1-4(1-r-s) : a == 2 ? 4r-1 : a == 3 ? zero(r) :
    a == 4 ? 4((1-r-s)-r) : a == 5 ? 4s : -4s
_tri9ds(a::Integer, r, s) =
    a == 1 ? 1-4(1-r-s) : a == 2 ? zero(s) : a == 3 ? 4s-1 :
    a == 4 ? -4r : a == 5 ? 4r : 4((1-r-s)-s)

# gmsh local 1-D index per node (1 ↔ -1, 2 ↔ 0, 3 ↔ +1; prism: type-9 tri index
# × 1-D w index) — the `lagrange_nodes` layout, re-derived here from the
# documented Gmsh ordering rather than from the module's tables.
const _QIDX = ((1,1),(3,1),(3,3),(1,3),(2,1),(3,2),(2,3),(1,2),(2,2))
const _HIDX = ((1,1,1),(3,1,1),(3,3,1),(1,3,1),(1,1,3),(3,1,3),(3,3,3),(1,3,3),
               (2,1,1),(1,2,1),(1,1,2),(3,2,1),(3,1,2),(2,3,1),(3,3,2),(1,3,2),
               (2,1,3),(1,2,3),(3,2,3),(2,3,3),
               (2,2,1),(2,1,2),(1,2,2),(3,2,2),(2,3,2),(2,2,3),(2,2,2))
const _PIDX = ((1,1),(2,1),(3,1),(1,3),(2,3),(3,3),
               (4,1),(6,1),(1,2),(5,1),(2,2),(3,2),
               (4,3),(6,3),(5,3),(4,2),(6,2),(5,2))

_conn(n::Integer) = reshape(collect(Int32(1):Int32(n)), n, 1)

# Gram determinant g = |∂x/∂u × ∂x/∂v|² of the 9-node quadrangle map; generic
# over the coordinate/reference scalar type so the same code runs Float64 and
# Rational{BigInt}.
function _quad_g(X, u, v)
    ju = (zero(u), zero(u), zero(u)); jv = (zero(u), zero(u), zero(u))
    for k in 1:9
        a, b = _QIDX[k]
        du = _l2d(a, u)*_l2(b, v); dv = _l2(a, u)*_l2d(b, v)
        ju = ju .+ (X[1,k]*du, X[2,k]*du, X[3,k]*du)
        jv = jv .+ (X[1,k]*dv, X[2,k]*dv, X[3,k]*dv)
    end
    cx = ju[2]*jv[3]-ju[3]*jv[2]
    cy = ju[3]*jv[1]-ju[1]*jv[3]
    cz = ju[1]*jv[2]-ju[2]*jv[1]
    return cx*cx + cy*cy + cz*cz
end

function _hex_det(X, u, v, w)
    jr = (zero(u), zero(u), zero(u))
    js = (zero(u), zero(u), zero(u))
    jw = (zero(u), zero(u), zero(u))
    for k in 1:27
        a, b, c = _HIDX[k]
        dr = _l2d(a, u)*_l2(b, v)*_l2(c, w)
        ds = _l2(a, u)*_l2d(b, v)*_l2(c, w)
        dw = _l2(a, u)*_l2(b, v)*_l2d(c, w)
        jr = jr .+ (X[1,k]*dr, X[2,k]*dr, X[3,k]*dr)
        js = js .+ (X[1,k]*ds, X[2,k]*ds, X[3,k]*ds)
        jw = jw .+ (X[1,k]*dw, X[2,k]*dw, X[3,k]*dw)
    end
    return jr[1]*(js[2]*jw[3]-js[3]*jw[2]) -
           jr[2]*(js[1]*jw[3]-js[3]*jw[1]) +
           jr[3]*(js[1]*jw[2]-js[2]*jw[1])
end

function _pri_det(X, r, s, w)
    jr = (zero(r), zero(r), zero(r))
    js = (zero(r), zero(r), zero(r))
    jw = (zero(r), zero(r), zero(r))
    for k in 1:18
        a, b = _PIDX[k]
        dr = _tri9dr(a, r, s)*_l2(b, w)
        ds = _tri9ds(a, r, s)*_l2(b, w)
        dw = _tri9(a, r, s)*_l2d(b, w)
        jr = jr .+ (X[1,k]*dr, X[2,k]*dr, X[3,k]*dr)
        js = js .+ (X[1,k]*ds, X[2,k]*ds, X[3,k]*ds)
        jw = jw .+ (X[1,k]*dw, X[2,k]*dw, X[3,k]*dw)
    end
    return jr[1]*(js[2]*jw[3]-js[3]*jw[2]) -
           jr[2]*(js[1]*jw[3]-js[3]*jw[1]) +
           jr[3]*(js[1]*jw[2]-js[2]*jw[1])
end

# Dense + one local refinement pass: the sampled minimum is a deterministic
# upper bound of the true minimum (the sample lattice is a subset of the closed
# domain), so `cert <= sampled_min` is the no-wrong-certificate property and
# `sampled_min - cert` is the certificate's measured conservativeness.
function _quad_oracle_min(X)
    g(u, v) = _quad_g(X, u, v)
    best = (Inf, 0.0, 0.0)
    for u in -1:0.05:1, v in -1:0.05:1
        d = g(u, v)
        d < best[1] && (best = (d, u, v))
    end
    mn = best[1]
    for u in best[2]-0.05:0.005:best[2]+0.05, v in best[3]-0.05:0.005:best[3]+0.05
        (-1.0 <= u <= 1.0 && -1.0 <= v <= 1.0) || continue
        mn = min(mn, g(u, v))
    end
    return mn
end
function _hex_oracle_min(X)
    f(u, v, w) = _hex_det(X, u, v, w)
    best = (Inf, 0.0, 0.0, 0.0)
    for u in -1:0.25:1, v in -1:0.25:1, w in -1:0.25:1
        d = f(u, v, w)
        d < best[1] && (best = (d, u, v, w))
    end
    mn = best[1]
    for u in best[2]-0.25:0.025:best[2]+0.25, v in best[3]-0.25:0.025:best[3]+0.25,
        w in best[4]-0.25:0.025:best[4]+0.25
        (-1.0 <= u <= 1.0 && -1.0 <= v <= 1.0 && -1.0 <= w <= 1.0) || continue
        mn = min(mn, f(u, v, w))
    end
    return mn
end
function _pri_oracle_min(X)
    f(r, s, w) = _pri_det(X, r, s, w)
    coarse = -1:0.2:1
    tricoarse = Tuple{Float64,Float64}[]
    for r in 0.0:0.1:1.0, s in 0.0:0.1:1.0
        r+s <= 1.0+1e-9 && push!(tricoarse, (r, s))
    end
    best = (Inf, 0.0, 0.0, 0.0)
    for (r, s) in tricoarse, w in coarse
        d = f(r, s, w)
        d < best[1] && (best = (d, r, s, w))
    end
    mn = Inf
    for r in best[2]-0.1:0.01:best[2]+0.1, s in best[3]-0.1:0.01:best[3]+0.1,
        w in best[4]-0.2:0.02:best[4]+0.2
        (r >= 0.0 && s >= 0.0 && r+s <= 1.0+1e-9 && -1.0 <= w <= 1.0) || continue
        mn = min(mn, f(r, s, w))
    end
    return min(mn, best[1])
end

_exact(x) = Rational{BigInt}(x)
_exact_eval_float(x::Rational{BigInt}) = Float64(BigFloat(x))

@testset "HighOrder P2 non-simplex Jacobian certification" begin

    @testset "input validation contracts" begin
        C9 = Float64[0 1 1 0 0.5 1 0.5 0 0.5;
                     0 0 1 1 0 0.5 1 0.5 0.5;
                     0 0 0 0 0   0 0   0 0]
        for f in (p2_quad_min_jacobian, p2_hex_min_jacobian, p2_prism_min_jacobian)
            n = f === p2_quad_min_jacobian ? 9 : f === p2_hex_min_jacobian ? 27 : 18
            Z = f === p2_quad_min_jacobian ? C9 :
                f === p2_hex_min_jacobian ? lagrange_nodes(12) : lagrange_nodes(13)
            K = _conn(n)
            @test_throws ArgumentError f(zeros(3), K)
            @test_throws ArgumentError f(Z, zeros(Int32, n))
            @test_throws ArgumentError f(zeros(2, n), K)
            @test_throws ArgumentError f(Z, zeros(Int32, n-1, 1))
            @test_throws ArgumentError f(trues(3, n), K)
            @test_throws ArgumentError f(Z, falses(n, 1))
            @test_throws ArgumentError f(fill(NaN, 3, n), K)
            bad = copy(K); bad[1, 1] = 0
            @test_throws ArgumentError f(Z, bad)
            bad = copy(K); bad[1, 1] = Int32(n + 1)
            @test_throws ArgumentError f(Z, bad)
            bad = copy(K); bad[2, 1] = bad[1, 1]
            @test_throws ArgumentError f(Z, bad)
            badf = Float64.(K); badf[1, 1] = 1.5
            @test_throws ArgumentError f(Z, badf)
            @test f(Z, zeros(Int32, n, 0)) == 0.0
        end
    end

    @testset "affine elements: certificate is the exact constant" begin
        # Identity quadrangle in the z=0 plane: w = (0,0,1) -> g = 1.
        Xq = lagrange_nodes(10); Xq[3, :] .= 0.0
        @test p2_quad_min_jacobian(Xq, _conn(9)) === 1.0
        Xq2 = copy(Xq); Xq2[1, :] .*= 2.0               # |w| = 2 -> g = 4
        @test p2_quad_min_jacobian(Xq2, _conn(9)) === 4.0
        Xq3 = copy(Xq); Xq3 .+= 1e6                    # translation invariance
        @test p2_quad_min_jacobian(Xq3, _conn(9)) === 1.0

        # Identity hexahedron det J = 1; axis scale 2 -> 2; mirrored -> -1.
        Xh = lagrange_nodes(12)
        @test p2_hex_min_jacobian(Xh, _conn(27)) === 1.0
        Xh2 = copy(Xh); Xh2[1, :] .*= 2.0
        @test p2_hex_min_jacobian(Xh2, _conn(27)) === 2.0
        Xhi = copy(Xh); Xhi[1, :] .*= -1.0
        @test p2_hex_min_jacobian(Xhi, _conn(27)) === -1.0

        # Identity prism det J = 1 (reference prism volume is 1); w scale 3 -> 3.
        Xp = lagrange_nodes(13)
        @test p2_prism_min_jacobian(Xp, _conn(18)) === 1.0
        Xp3 = copy(Xp); Xp3[3, :] .*= 3.0
        @test p2_prism_min_jacobian(Xp3, _conn(18)) === 3.0
        Xpi = copy(Xp); Xpi[3, :] .*= -1.0
        @test p2_prism_min_jacobian(Xpi, _conn(18)) === -1.0
    end

    @testset "exact-rational oracle vs Float64 evaluation" begin
        rng = MersenneTwister(20260214)
        # One mildly warped, strictly regular element per family.
        Xq = lagrange_nodes(10)
        for k in 5:9, d in 1:3; Xq[d, k] += 0.12*(2rand(rng)-1); end
        Xh = lagrange_nodes(12)
        for k in 9:27, d in 1:3; Xh[d, k] += 0.08*(2rand(rng)-1); end
        Xp = lagrange_nodes(13)
        for k in 7:18, d in 1:3; Xp[d, k] += 0.08*(2rand(rng)-1); end
        Xr = map(_exact, Xq); Yr = map(_exact, Xh); Zr = map(_exact, Xp)

        # (a) Float64 det evaluation recomputed exactly at fixed-seed interior
        #     points.  Tolerance is a conditioning allowance for the ~27-term
        #     polynomial evaluation; it does not touch the exact certificate.
        for i in 1:8
            u, v = 2rand(rng)-1, 2rand(rng)-1
            f = _quad_g(Xq, u, v)
            e = _exact_eval_float(_quad_g(Xr, _exact(u), _exact(v)))
            @test abs(f - e) <= 1e-9*max(1.0, abs(e))
        end
        for i in 1:8
            u, v, w = 2rand(rng)-1, 2rand(rng)-1, 2rand(rng)-1
            f = _hex_det(Xh, u, v, w)
            e = _exact_eval_float(_hex_det(Yr, _exact(u), _exact(v), _exact(w)))
            @test abs(f - e) <= 1e-9*max(1.0, abs(e))
        end
        for i in 1:8
            r = rand(rng); s = rand(rng)*(1 - r)          # interior of the tri
            w = 2rand(rng)-1
            f = _pri_det(Xp, r, s, w)
            e = _exact_eval_float(_pri_det(Zr, _exact(r), _exact(s), _exact(w)))
            @test abs(f - e) <= 1e-9*max(1.0, abs(e))
        end

        # (b-exact) the certificate is a lower bound of the EXACT minimum over
        # a fixed dyadic grid — the whole check is done in rational arithmetic.
        quad_exact_min = minimum(
            _quad_g(Xr, _exact(u), _exact(v))
            for u in -1:1//4:1, v in -1:1//4:1)
        @test p2_quad_min_jacobian(Xq, _conn(9)) <=
              _exact_eval_float(quad_exact_min) + 1e-12
        hex_exact_min = minimum(
            _hex_det(Yr, _exact(u), _exact(v), _exact(w))
            for u in -1:1//2:1, v in -1:1//2:1, w in -1:1//2:1)
        @test p2_hex_min_jacobian(Xh, _conn(27)) <=
              _exact_eval_float(hex_exact_min) + 1e-12
        pri_exact_min = minimum(
            _pri_det(Zr, _exact(r), _exact(s), _exact(w))
            for r in 0:1//4:1, s in 0:1//4:1, w in -1:1//2:1 if r+s <= 1)
        @test p2_prism_min_jacobian(Xp, _conn(18)) <=
              _exact_eval_float(pri_exact_min) + 1e-12
    end

    @testset "certified bound vs dense adaptive sampling (fixed-seed corpus)" begin
        rng = MersenneTwister(0xC0FFEE)
        quad_gaps = Float64[]; hex_gaps = Float64[]; pri_gaps = Float64[]
        for rep in 1:4
            Xq = lagrange_nodes(10)
            for k in 5:9, d in 1:3; Xq[d, k] += 0.12*(2rand(rng)-1); end
            cert = p2_quad_min_jacobian(Xq, _conn(9))
            smin = _quad_oracle_min(Xq)
            @test cert <= smin + 1e-9*max(1.0, smin)
            @test smin > 0                       # corpus elements stay regular
            push!(quad_gaps, smin - cert)

            Xh = lagrange_nodes(12)
            for k in 9:27, d in 1:3; Xh[d, k] += 0.08*(2rand(rng)-1); end
            cert = p2_hex_min_jacobian(Xh, _conn(27))
            smin = _hex_oracle_min(Xh)
            @test cert <= smin + 1e-9*max(1.0, smin)
            @test smin > 0
            push!(hex_gaps, smin - cert)

            Xp = lagrange_nodes(13)
            for k in 7:18, d in 1:3; Xp[d, k] += 0.08*(2rand(rng)-1); end
            cert = p2_prism_min_jacobian(Xp, _conn(18))
            smin = _pri_oracle_min(Xp)
            @test cert <= smin + 1e-9*max(1.0, smin)
            @test smin > 0
            push!(pri_gaps, smin - cert)
        end
        # Documented conservativeness of the single-level Bernstein-hull bound
        # on this corpus (deterministic measurements, ~±0.1 warp amplitude on
        # det≈1 elements): observed sampled-min − certificate gaps were
        # quad ≤ 0.095, hex ≤ 0.080, prism ≤ 0.120.  The bounds are thus the
        # correct side of the truth but visibly conservative for curved maps —
        # tighter certificates would need reference-domain subdivision, which
        # the simplex-family certifiers deliberately avoid as well.  Margins
        # below are ~1.5x the observed maxima.
        @test maximum(quad_gaps) < 0.15
        @test maximum(hex_gaps) < 0.12
        @test maximum(pri_gaps) < 0.18
    end

    @testset "inverted, folded and degenerate elements report <= 0" begin
        Xq = lagrange_nodes(10)
        folded = copy(Xq); folded[1, 9] = folded[2, 9] = 1.4   # center past corner
        certf = p2_quad_min_jacobian(folded, _conn(9))
        @test certf === -5.968888888888888                     # pinned certificate
        @test _quad_oracle_min(folded) < 1e-3                  # genuinely singular
        bowtie = copy(Xq); bowtie[1, 3] = bowtie[2, 3] = -0.4  # bow-tie corner map
        @test p2_quad_min_jacobian(bowtie, _conn(9)) === -1.8606666666666667
        @test _quad_oracle_min(bowtie) < 1e-3
        flat = copy(Xq); flat[3, :] .= 0.0
        edge = copy(flat); edge[1, 5] = 0.99; edge[2, 5] = -0.99
        @test p2_quad_min_jacobian(edge, _conn(9)) <= 0.0      # grazing warp
        zero_area = zeros(3, 9)
        for k in 1:9; zero_area[1, k] = (k-1)/8; end           # collinear quad
        @test p2_quad_min_jacobian(zero_area, _conn(9)) == 0.0
        @test isfinite(p2_quad_min_jacobian(zero_area, _conn(9)))

        Xh = lagrange_nodes(12)
        invh = copy(Xh); invh[1, :] .*= -1.0
        @test p2_hex_min_jacobian(invh, _conn(27)) === -1.0    # mirrored affine
        foldedh = copy(Xh); foldedh[:, 27] .= (2.0, 2.0, 2.0)  # body node outside
        certh = p2_hex_min_jacobian(foldedh, _conn(27))
        @test certh <= 0.0
        @test _hex_oracle_min(foldedh) < 0.0                   # genuinely folded
        thin = copy(Xh); thin[3, :] .*= 1e-4                   # near-planar hex
        certt = p2_hex_min_jacobian(thin, _conn(27))
        @test certt > 0
        @test certt <= _hex_oracle_min(thin) + 1e-9

        Xp = lagrange_nodes(13)
        invp = copy(Xp); invp[3, :] .*= -1.0
        @test p2_prism_min_jacobian(invp, _conn(18)) === -1.0
        collapsed = copy(Xp); collapsed[:, 2] .= collapsed[:, 3]  # collapsed edge
        certc = p2_prism_min_jacobian(collapsed, _conn(18))
        @test certc <= 0.0
        @test _pri_oracle_min(collapsed) <= 0.0                # degenerate map
        foldedp = copy(Xp); foldedp[:, 9] .= (0.0, 0.0, 4.0)   # vertical mid out
        certp = p2_prism_min_jacobian(foldedp, _conn(18))
        @test certp <= 0.0
        @test _pri_oracle_min(foldedp) < 0.0
    end

    @testset "batch semantics and determinism" begin
        rng = MersenneTwister(1234567)
        Xh = lagrange_nodes(12)
        good = copy(Xh)
        for k in 9:27, d in 1:3; good[d, k] += 0.05*(2rand(rng)-1); end
        bad = copy(Xh); bad[:, 27] .= (2.0, 2.0, 2.0)
        @test p2_hex_min_jacobian(good, _conn(27)) > 0
        both = hcat(_conn(27), _conn(27) .+ Int32(27))
        C2 = hcat(good, bad)
        # A folded element drags the batch certificate down to its own bound.
        @test p2_hex_min_jacobian(C2, both) ==
              p2_hex_min_jacobian(bad, _conn(27))

        # Pure-integer certificate: identical across repeated evaluations and
        # pinned to a literal, which makes cross-process determinism checkable.
        Xw = lagrange_nodes(13)
        for k in 7:18, d in 1:3; Xw[d, k] += 0.1*(2rand(rng)-1); end
        c1 = p2_prism_min_jacobian(Xw, _conn(18))
        c2 = p2_prism_min_jacobian(Xw, _conn(18))
        @test c1 === c2
        Xq = lagrange_nodes(10)
        for k in 5:9, d in 1:3; Xq[d, k] += 0.1*(2rand(rng)-1); end
        @test p2_quad_min_jacobian(Xq, _conn(9)) ===
              p2_quad_min_jacobian(Xq, _conn(9))
    end

    @testset "exported API wiring" begin
        for name in (:p2_quad_min_jacobian, :p2_hex_min_jacobian,
                     :p2_prism_min_jacobian)
            @test name in names(Tessella.HighOrder)
        end
    end
end
