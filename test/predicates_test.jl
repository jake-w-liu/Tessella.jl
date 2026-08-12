# ── Stage-0 CRC suite: exact predicates vs an independent exact-rational oracle ──
#
# Correctness  : every predicate sign is cross-checked against Oracles.* (a
#                generic Laplace determinant over the full homogeneous matrix in
#                Rational{BigInt} — an independent computation).
# Robustness   : exhaustive degenerate configs (collinear / cocircular / coplanar
#                / cospherical) on small integer grids force the exact fallback,
#                plus fixed-seed random doubles that exercise the fast filter.
# Completeness : the SoS variants never return 0, agree with the base predicate
#                off-degeneracy, are deterministic and antisymmetric.
#
# Mutation sensitivity: flip any sign or swap any index in src/Predicates.jl and
# the oracle cross-check below fails on a degenerate config.

using Test
using Random
using Tessella.Predicates
include("oracles.jl")
using .Oracles

@testset "Predicates (Stage 0)" begin

    # ── orient2: exhaustive integer grid, every ordered triple ──────────────────
    @testset "orient2 vs oracle — integer grid" begin
        grid = [(x, y) for x in -3:3 for y in -3:3]
        mism = 0
        ncollinear = 0
        for a in grid, b in grid, c in grid
            (a == b || a == c || b == c) && continue
            got = orient2(a, b, c)
            want = oracle_orient2(a, b, c)
            got == want || (mism += 1)
            want == 0 && (ncollinear += 1)
        end
        @test mism == 0
        @test ncollinear > 0            # the grid really does contain collinear triples
    end

    # ── orient2: antisymmetry under argument transposition ──────────────────────
    @testset "orient2 antisymmetry" begin
        rng = 0
        ok = true
        for a in [(0,0),(1,0),(2,1),(3,3),(-2,1)], b in [(1,1),(0,2),(2,2)], c in [(2,0),(1,3),(-1,-1)]
            (a==b||a==c||b==c) && continue
            ok &= orient2(a,b,c) == -orient2(b,a,c)
            ok &= orient2(a,b,c) == -orient2(a,c,b)
            ok &= orient2(a,b,c) ==  orient2(b,c,a)   # cyclic → even
        end
        @test ok
    end

    # ── orient3: integer grid (smaller), every ordered 4-tuple sampled ──────────
    @testset "orient3 vs oracle — integer grid" begin
        grid = [(x, y, z) for x in -2:2 for y in -2:2 for z in -2:2]
        # full 4-tuple loop is 15625^... too big; sample deterministically.
        mism = 0; ncoplanar = 0; ntested = 0
        step = 7
        n = length(grid)
        i = 1
        for ia in 1:n
            a = grid[ia]
            for ib in ia+1:n
                b = grid[ib]
                for ic in ib+1:n
                    ((i += 1) % step == 0) || continue
                    c = grid[ic]
                    # pick a few d's deterministically
                    for id in (ic+3, ic+11, ic+29)
                        id > n && continue
                        d = grid[id]
                        (a==d||b==d||c==d) && continue
                        got = orient3(a,b,c,d); want = oracle_orient3(a,b,c,d)
                        got == want || (mism += 1)
                        want == 0 && (ncoplanar += 1)
                        ntested += 1
                    end
                end
            end
        end
        @test mism == 0
        @test ntested > 10_000
        @test ncoplanar > 0
    end

    # ── orient3: explicit coplanar degeneracies (force exact path) ──────────────
    @testset "orient3 coplanar degeneracies" begin
        # all points on z = 0 plane
        for a in [(0,0,0)], b in [(1,0,0),(2,1,0)], c in [(0,1,0),(3,3,0)], d in [(2,2,0),(-1,-1,0)]
            @test orient3(a,b,c,d) == 0
            @test oracle_orient3(a,b,c,d) == 0
        end
        # a general plane 2x + 3y - z = 5  → z = 2x+3y-5
        f(x,y) = (x, y, 2x+3y-5)
        pts = [f(0,0), f(1,0), f(0,1), f(2,3), f(-1,2), f(4,-2)]
        for i in 1:length(pts), j in i+1:length(pts), k in j+1:length(pts), l in k+1:length(pts)
            @test orient3(pts[i],pts[j],pts[k],pts[l]) == 0
        end
    end

    # ── incircle: integer grid triples + query points, vs oracle ────────────────
    @testset "incircle vs oracle — integer grid" begin
        grid = [(x, y) for x in -3:3 for y in -3:3]
        mism = 0; ncocirc = 0; ntested = 0
        n = length(grid)
        i = 0
        for ia in 1:n, ib in ia+1:n, ic in ib+1:n
            # ensure a,b,c not collinear
            oracle_orient2(grid[ia],grid[ib],grid[ic]) == 0 && continue
            i += 1
            (i % 13 == 0) || continue
            a=grid[ia]; b=grid[ib]; c=grid[ic]
            for id in 1:4:n
                d = grid[id]
                (d==a||d==b||d==c) && continue
                got = incircle(a,b,c,d); want = oracle_incircle(a,b,c,d)
                got == want || (mism += 1)
                want == 0 && (ncocirc += 1)
                ntested += 1
            end
        end
        @test mism == 0
        @test ntested > 5_000
        @test ncocirc > 0
    end

    # ── incircle: exact cocircular set (points on a circle radius² = 25) ─────────
    @testset "incircle cocircular degeneracies" begin
        circ = [(5,0),(4,3),(3,4),(0,5),(-3,4),(-4,3),(-5,0),(0,-5),(3,-4),(4,-3)]
        # any 4 points on this circle are cocircular ⇒ incircle == 0
        for i in 1:length(circ), j in i+1:length(circ), k in j+1:length(circ), l in k+1:length(circ)
            @test incircle(circ[i],circ[j],circ[k],circ[l]) == 0
        end
        # a point strictly inside / outside the radius-5 circle
        @test oracle_incircle((5,0),(0,5),(-5,0),(0,0)) == incircle((5,0),(0,5),(-5,0),(0,0))
        @test incircle((5,0),(0,5),(-5,0),(0,0)) != 0     # (0,0) is inside → nonzero
        @test incircle((5,0),(0,5),(-5,0),(9,9)) != 0     # far outside
    end

    # ── insphere: fixed-seed random + explicit cospherical, vs oracle ───────────
    @testset "insphere vs oracle" begin
        # deterministic pseudo-points on a small integer grid
        grid = [(x,y,z) for x in -2:2 for y in -2:2 for z in -2:2]
        n = length(grid)
        mism = 0; ntested = 0
        # sample tets + query
        for s in 1:1500
            ia = ((3s)   % n) + 1
            ib = ((5s+1) % n) + 1
            ic = ((7s+2) % n) + 1
            id = ((11s+3)% n) + 1
            ie = ((13s+4)% n) + 1
            idx = (ia,ib,ic,id,ie)
            length(unique(idx)) == 5 || continue
            a,b,c,d,e = grid[ia],grid[ib],grid[ic],grid[id],grid[ie]
            oracle_orient3(a,b,c,d) == 0 && continue   # need a real tetrahedron
            got = insphere(a,b,c,d,e); want = oracle_insphere(a,b,c,d,e)
            got == want || (mism += 1)
            ntested += 1
        end
        @test mism == 0
        @test ntested > 500
    end

    @testset "insphere cospherical degeneracies" begin
        # 8 vertices of a cube inscribed in a sphere: all cospherical
        cube = [(x,y,z) for x in (-1,1) for y in (-1,1) for z in (-1,1)]
        cnt = 0
        for i in 1:8, j in i+1:8, k in j+1:8, l in k+1:8, m in l+1:8
            a,b,c,d,e = cube[i],cube[j],cube[k],cube[l],cube[m]
            oracle_orient3(a,b,c,d) == 0 && continue
            @test insphere(a,b,c,d,e) == 0
            cnt += 1
        end
        @test cnt > 0
        # center of the cube is strictly inside the circumsphere
        a,b,c,d = (-1,-1,-1),(1,-1,-1),(-1,1,-1),(-1,-1,1)
        @test insphere(a,b,c,d,(0,0,0)) == oracle_insphere(a,b,c,d,(0,0,0))
        @test insphere(a,b,c,d,(0,0,0)) != 0
    end

    # ── fast-filter path exercised with random Float64 (agreement with oracle) ──
    @testset "float filter path vs oracle" begin
        # a tiny xorshift RNG for reproducibility without external deps
        state = UInt64(0x243F6A8885A308D3)
        nextf() = (state = state ⊻ (state << 13); state = state ⊻ (state >> 7);
                   state = state ⊻ (state << 17); (state >> 11) / Float64(2^53) - 0.5)
        for _ in 1:8_000
            a=(nextf(),nextf()); b=(nextf(),nextf()); c=(nextf(),nextf())
            @test orient2(a,b,c) == oracle_orient2(a,b,c)
        end
        for _ in 1:8_000
            a=(nextf(),nextf(),nextf()); b=(nextf(),nextf(),nextf())
            c=(nextf(),nextf(),nextf()); d=(nextf(),nextf(),nextf())
            @test orient3(a,b,c,d) == oracle_orient3(a,b,c,d)
        end
        for _ in 1:3_000
            a=(nextf(),nextf()); b=(nextf(),nextf()); c=(nextf(),nextf()); d=(nextf(),nextf())
            @test incircle(a,b,c,d) == oracle_incircle(a,b,c,d)
        end
        for _ in 1:2_500
            g() = (nextf(),nextf(),nextf())
            a,b,c,d,e = g(),g(),g(),g(),g()
            @test insphere(a,b,c,d,e) == oracle_insphere(a,b,c,d,e)
        end
    end

    # ── near-degenerate perturbations that defeat a naive Float64 predicate ─────
    @testset "near-degenerate robustness" begin
        # points almost collinear: c just off the line through a,b by 1 ulp
        a = (0.5, 0.5); b = (12.0, 12.0)
        for k in -5:5
            c = (1.0, nextfloat(1.0, k))
            @test orient2(a, b, c) == oracle_orient2(a, b, c)
        end
        # classic Shewchuk failure family: tiny grid near a shared point
        base = 0.5
        for i in 0:20, j in 0:20
            a = (base, base)
            b = (base + 1e-15*i, base + 1e-15*j)
            c = (base + 1e-15*j, base + 1e-15*i)
            @test orient2(a, b, c) == oracle_orient2(a, b, c)
        end
    end

    # ── SoS contract: never 0, agrees off-degeneracy, deterministic, antisym ────
    @testset "orient2_sos contract" begin
        grid = [(x,y) for x in -2:2 for y in -2:2]
        n = length(grid)
        for i in 1:n, j in 1:n, k in 1:n
            (i==j||i==k||j==k) && continue
            a,b,c = grid[i],grid[j],grid[k]
            s = orient2_sos(a,b,c, i,j,k)
            @test s != 0
            base = orient2(a,b,c)
            base != 0 && @test s == base
            @test s == orient2_sos(a,b,c, i,j,k)              # deterministic
            @test s == -orient2_sos(b,a,c, j,i,k)             # antisymmetric
            @test s == -orient2_sos(a,c,b, i,k,j)
        end
    end

    @testset "orient3_sos contract" begin
        grid = [(x,y,z) for x in -1:1 for y in -1:1 for z in -1:1]
        n = length(grid)
        cnt = 0
        for i in 1:n, j in i+1:n, k in j+1:n, l in k+1:n
            a,b,c,d = grid[i],grid[j],grid[k],grid[l]
            s = orient3_sos(a,b,c,d, i,j,k,l)
            @test s != 0
            base = orient3(a,b,c,d)
            base != 0 && @test s == base
            @test s == -orient3_sos(b,a,c,d, j,i,k,l)         # one transposition
            cnt += 1
        end
        @test cnt > 100
    end

    @testset "insphere_sos / incircle_sos never zero" begin
        cube = [(x,y,z) for x in (-1,1) for y in (-1,1) for z in (-1,1)]
        for i in 1:8, j in i+1:8, k in j+1:8, l in k+1:8, m in l+1:8
            a,b,c,d,e = cube[i],cube[j],cube[k],cube[l],cube[m]
            orient3(a,b,c,d) == 0 && continue
            @test insphere_sos(a,b,c,d,e, i,j,k,l,m) != 0
        end
        circ = [(5,0),(4,3),(3,4),(0,5),(-3,4),(-4,3),(-5,0),(0,-5)]
        for i in 1:8, j in i+1:8, k in j+1:8, l in k+1:8
            orient2(circ[i],circ[j],circ[k]) == 0 && continue
            @test incircle_sos(circ[i],circ[j],circ[k],circ[l], i,j,k,l) != 0
        end
    end

    # ── SoS +ε tie-break signs on FULLY-degenerate configs (regression pins) ─────
    # The four *_sos predicates must resolve a fully-degenerate determinant (4
    # collinear points / 5 coplanar points — where the fast-path minors all vanish)
    # by the SAME +ε Edelsbrunner–Mücke perturbation. The golden sign over every index
    # labeling below was captured from the fixed predicates and cross-verified against
    # an exact Rational{BigInt} leading-term SoS oracle AND a canonical paraboloid-
    # lifted-orientation oracle (STATUS.md "Audit findings"). They pin: orient2's
    # −ε→+ε sign flip; orient3's & incircle's exact collinear fallthrough (was a
    # constant `perm`); and insphere's −ε→+ε flip via the 4-D lift.
    @testset "SoS +ε signs on fully-degenerate configs (regression)" begin
        aperms(v::Vector{Int}) = length(v) <= 1 ? [v] :
            [ [v[i]; p] for i in eachindex(v) for p in aperms(v[[1:i-1; i+1:end]]) ]
        sigs(f, n) = [f(Tuple(t)) for t in aperms(collect(1:n))]

        p2 = ((0,0),(1,0),(2,0))                                   # collinear
        @test sigs(t -> orient2_sos(p2[1],p2[2],p2[3], t...), 3) == [1,1,-1,1,-1,1]

        p3 = ((-2,-2,0),(-1,-2,0),(0,-2,0),(1,-2,0))               # collinear (was wrong 12/24)
        @test sigs(t -> orient3_sos(p3[1],p3[2],p3[3],p3[4], t...), 4) ==
              [-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,1,-1,-1,1,-1,1,1,-1,-1,1]

        pic = ((0,0),(1,0),(2,0),(3,0))                            # collinear (was wrong 12/24)
        @test sigs(t -> incircle_sos(pic[1],pic[2],pic[3],pic[4], t...), 4) ==
              [-1,-1,-1,-1,-1,-1,1,1,-1,1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,1]

        pis = ((0,0,0),(1,0,0),(0,1,0),(1,1,0),(2,0,0))            # coplanar (was wrong 60/120)
        @test sigs(t -> insphere_sos(pis[1],pis[2],pis[3],pis[4],pis[5], t...), 5) ==
              [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,
               -1,-1,1,-1,1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,-1,-1,1,1,1,-1,
               -1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1,
               -1,1,1,1,-1,-1,-1,1,1,1,-1]

        # every one is a valid, antisymmetric, nonzero tie-break (no coincidental perm)
        @test all(!=(0), sigs(t -> insphere_sos(pis[1],pis[2],pis[3],pis[4],pis[5], t...), 5))
    end

    # ── allocation guard: fast path must not allocate ───────────────────────────
    @testset "fast path is allocation-free" begin
        a=(0.0,0.0); b=(1.0,0.0); c=(0.0,1.0)
        orient2(a,b,c)                        # warmup / compile
        @test (@allocated orient2(a,b,c)) == 0
        a3=(0.0,0.0,0.0); b3=(1.0,0.0,0.0); c3=(0.0,1.0,0.0); d3=(0.0,0.0,1.0)
        orient3(a3,b3,c3,d3)
        @test (@allocated orient3(a3,b3,c3,d3)) == 0
        incircle(a,b,c,(0.4,0.4)); insphere(a3,b3,c3,d3,(0.1,0.1,0.1))
        @test (@allocated incircle(a,b,c,(0.4,0.4))) == 0
        @test (@allocated insphere(a3,b3,c3,d3,(0.1,0.1,0.1))) == 0
    end

    @testset "orient2_sos degenerate ties match the canonical +ε SoS (four-predicate consistency)" begin
        # Regression: orient2_sos's degenerate branch now routes through the SAME
        # canonical evaluator (_orient_nd_sos_exact) as orient3_sos/incircle_sos/
        # insphere_sos, so all four agree under one +ε scheme on EVERY input. The old
        # hand-coded first-order minor sequence agreed on distinct points but disagreed
        # on coincident coordinates (measured 144/4416); it never returns 0.
        pts = [(Float64(x),Float64(y)) for x in 0:3 for y in 0:3]
        perms = [(1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1)]
        function scan()
            dis=0; nz=0; ntot=0
            for pa in pts, pb in pts, pc in pts
                orient2(pa,pb,pc) == 0 || continue
                for (ia,ib,ic) in perms
                    s1 = orient2_sos(pa,pb,pc, ia,ib,ic); ntot += 1
                    s1 == 0 && (nz += 1)
                    s2 = Predicates._orient_nd_sos_exact(((pa[1],pa[2]),(pb[1],pb[2]),(pc[1],pc[2])),(ia,ib,ic),2)
                    s1 == s2 || (dis += 1)
                end
            end
            (dis, nz, ntot)
        end
        d, z, n = scan()
        @test n > 0
        @test z == 0            # a *_sos predicate must never return 0
        @test d == 0            # full agreement with the canonical evaluator
        # the exact finding config (coincident): +1 canonical, was −1 under the old branch
        @test orient2_sos((0.0,0.0),(0.0,0.0),(0.0,1.0), 2,3,1) == 1
    end

    @testset "exact-rational predicates match production; orient3_sos coplanar consistency" begin
        # The exact-coordinate kernel's predicates (orient*/incircle/insphere_rat on
        # Rational{BigInt}) evaluate the SAME homogeneous determinants and share the
        # canonical +ε SoS, so on dyadic input they must agree with the Float64 SoS
        # predicates for EVERY configuration — distinct AND coincident, general AND
        # degenerate. This also pins the orient3_sos fix: its degenerate branch now
        # routes through the canonical evaluator (the old 8-minor path miscomputed the
        # +ε sign on some coplanar-distinct configs).
        Q(x) = Rational{BigInt}(Float64(x))
        r3(p) = (Q(p[1]), Q(p[2]), Q(p[3])); r2(p) = (Q(p[1]), Q(p[2]))
        function scan()
            G = Float64.(-2:2); rng = MersenneTwister(20260812)
            d = 0; n = 0
            for _ in 1:6000
                a=(rand(rng,G),rand(rng,G),rand(rng,G)); b=(rand(rng,G),rand(rng,G),rand(rng,G))
                c=(rand(rng,G),rand(rng,G),rand(rng,G)); dd=(rand(rng,G),rand(rng,G),rand(rng,G)); e=(rand(rng,G),rand(rng,G),rand(rng,G))
                ix=(2,4,6,8,10); n += 1
                orient3_sos(a,b,c,dd,ix[1:4]...)==orient3_rat(r3(a),r3(b),r3(c),r3(dd),ix[1:4]...) || (d+=1)
                insphere_sos(a,b,c,dd,e,ix...)==insphere_rat(r3(a),r3(b),r3(c),r3(dd),r3(e),ix...) || (d+=1)
                a2=(rand(rng,G),rand(rng,G));b2=(rand(rng,G),rand(rng,G));c2=(rand(rng,G),rand(rng,G));d2=(rand(rng,G),rand(rng,G))
                orient2_sos(a2,b2,c2,ix[1:3]...)==orient2_rat(r2(a2),r2(b2),r2(c2),ix[1:3]...) || (d+=1)
                incircle_sos(a2,b2,c2,d2,ix[1:4]...)==incircle_rat(r2(a2),r2(b2),r2(c2),r2(d2),ix[1:4]...) || (d+=1)
            end
            (d, n)
        end
        d, n = scan()
        @test n > 0
        @test d == 0
        # the confirmed coplanar-distinct config (all on x=2): true +ε sign is +1 (was −1)
        @test orient3_sos((2.,1.,0.),(2.,-2.,1.),(2.,-2.,2.),(2.,-2.,-2.), 2,4,6,8) == 1
        @test orient3_rat((Q(2),Q(1),Q(0)),(Q(2),Q(-2),Q(1)),(Q(2),Q(-2),Q(2)),(Q(2),Q(-2),Q(-2)), 2,4,6,8) == 1
    end
end
