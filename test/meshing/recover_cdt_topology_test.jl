# ── Float-assisted topology orientation (RecoverCDT._partition_tets_float) ─────
#
# Regression for the dense-curved-PLC recovery blocker (t19 pipe volume): the
# emitted working-tet orientation must be judged on the Float64 output
# coordinates — the same coordinates `validate` audits through the exact
# `orient3` predicate — not on the arbitrary-rational preimage.  Recovery
# Steiner vertices sit ~1 ulp off their Float64 images, so on a near-flat tet
# the two determinants can disagree in sign: orienting by the preimage swapped
# Float64-positive cells into negatively-oriented output, and every topology
# attempt then failed ("dyadic-exact topology was invalid — N negatively-
# oriented tets"), aborting boundary recovery on non-star PLCs.

using Test
using Tessella
using Tessella.MeshTypes
using Tessella.Geometry: box_surface
using Tessella.RecoverCDT: recover_boundary_cdt, recover_partition_cdt

const RC = Tessella.RecoverCDT
const RB = Rational{BigInt}
const RBv = NTuple{3,RB}

_rsub3(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
_det3(a,b,c) = a[1]*(b[2]*c[3]-b[3]*c[2]) -
               a[2]*(b[1]*c[3]-b[3]*c[1]) +
               a[3]*(b[1]*c[2]-b[2]*c[1])
_sdet(pts,q) = _det3(_rsub3(pts[q[2]],pts[q[1]]),
                     _rsub3(pts[q[3]],pts[q[1]]),
                     _rsub3(pts[q[4]],pts[q[1]]))

@testset "Float-assisted topology orients emitted tets on output coordinates" begin
    # Deterministic sign-flip fixture: point 4 sits exactly 1e-18 below the
    # (1,2,3) supporting plane — a sub-ulp offset whose Float64 image still
    # lands on the rounded side that keeps the dyadic determinant positive.
    # The dyadic Delaunay therefore emits a cell that is positive on the
    # Float64 output coordinates but negative on the exact preimage; the
    # preimage-oriented emitter swapped it into a negatively-oriented tet and
    # every attempt failed the mesh audit.
    epsr = RB(1)//big(10)^18
    pts = RBv[(0,0,0), (3,0,0), (0,3,-1), (1,1-epsr,-1//3), (1,1,3)]
    fpts = RBv[(RB(Float64(p[1])),RB(Float64(p[2])),RB(Float64(p[3]))) for p in pts]
    # The fixture must actually exercise the defect: some candidate tet has
    # Float64-positive but exact-negative determinant.
    candidates = ((1,2,3,4),(1,2,4,3),(1,3,2,4),(1,3,4,2),(1,4,2,3),(1,4,3,2))
    @test any(q -> _sdet(fpts,q) > 0 && _sdet(pts,q) < 0, candidates)

    limits = RC._recovery_limits("recover_partition_cdt")
    for exact_first in (true, false)
        etets = RC._partition_tets_float(pts, limits; exact_first=exact_first)
        @test !isempty(etets)
        used = falses(length(pts))
        for q in etets
            # Every emitted tet is positively oriented on the coordinates the
            # mesh audit certifies.
            @test _sdet(fpts,q) > 0
            for v in q; used[v] = true; end
        end
        @test all(used)
        tm = Matrix{Int32}(undef,4,length(etets))
        for (j,t) in enumerate(etets), i in 1:4; tm[i,j] = Int32(t[i]); end
        coords = Matrix{Float64}(undef,3,length(pts))
        for i in eachindex(pts), k in 1:3; coords[k,i] = Float64(pts[i][k]); end
        m = Mesh(coords; tets=tm)
        @test validate(m).ok
    end
end

@testset "recover_partition_cdt still conforms on a material partition" begin
    # A regular octahedron nested strictly inside a box: two material regions
    # joined by the non-prismatic octahedron interface.  The joint PLC routes
    # through the same float-assisted topology path, then certifies the
    # partition exactly (innermost surface first).
    r = 0.5
    oct_coords = Float64[ r 0 0;  -r 0 0;  0 r 0;  0 -r 0;  0 0 r;  0 0 -r ]'
    oct_tris = Int32[1 3 5; 3 2 5; 2 4 5; 4 1 5; 3 1 6; 2 3 6; 4 2 6; 1 4 6]'
    oct = Mesh(oct_coords; tris=oct_tris)
    outer = box_surface(-1,1,-1,1,-1,1)
    m = recover_partition_cdt(Any[oct, outer])
    @test validate(m).ok
    @test is_closed_manifold(m)
    tags = sort!(unique(collect(m.tet_tag)))
    @test tags == Int32[1,2]
    vol = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),
                         node(m,m.tets[3,t]),node(m,m.tets[4,t]))
              for t in 1:ntets(m); init=0.0)
    @test vol ≈ 8.0 rtol=1e-9
    # Inner octahedron region: |V| = 4r³/3 with r = 1/2.
    inner = [t for t in 1:ntets(m) if m.tet_tag[t] == Int32(1)]
    ivol = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),
                          node(m,m.tets[3,t]),node(m,m.tets[4,t]))
               for t in inner; init=0.0)
    @test ivol ≈ 4r^3/3 rtol=1e-9
end

# ── dense curved PLC: input-proportional maxpts ────────────────────────────
#
# Regression for the t13 verbatim arm (`Field[1].F = "4"`, ~2470-facet
# aneurysm wall): `recover_boundary_cdt` converges that PLC at 6238 points —
# ~5.1× its 1233 input vertices — but the flat default `maxpts=6000` rejected
# it while splitting creases / inserting pierce points.  A dense curved PLC's
# Steiner demand is proportional to its input density, so the default
# `maxpts=0` now resolves to `max(6000, 8·input vertices)` once the
# deduplicated PLC is known; explicit values remain literal hard caps.
#
# The fixture below is the real t19 helical-pipe `Volume[12]` PLC (dense
# curved wall, 182 vertices / 360 facets, non-star): recovery must insert 31
# exact boundary Steiner points (creases split, pierce points, coverage
# centroids) before the input facets are Delaunay-recoverable.  The exact
# certificate runs inside `recover_boundary_cdt`; returning without throwing
# is itself the conformance evidence.  Full-scale verification of the
# 1233-vertex case lives in the opt-in `validation/gmsh_parity/t13.jl`
# driver (~12 minutes of exact recovery).

function _load_plc_txt(path)
    lines = readlines(path)
    n, f = parse.(Int, split(lines[1]))
    coords = Matrix{Float64}(undef, 3, n)
    for i in 1:n
        coords[:, i] .= parse.(Float64, split(lines[1+i]))
    end
    tris = Matrix{Int32}(undef, 3, f)
    for c in 1:f
        tris[:, c] .= parse.(Int32, split(lines[1+n+c]))
    end
    return Mesh(coords; tris=tris)
end

_tri_area(a, b, c) = 0.5 * sqrt(
    ((b[2]-a[2])*(c[3]-a[3]) - (b[3]-a[3])*(c[2]-a[2]))^2 +
    ((b[3]-a[3])*(c[1]-a[1]) - (b[1]-a[1])*(c[3]-a[3]))^2 +
    ((b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]))^2)

function _surf_area(m)
    s = 0.0
    for t in 1:ntris(m)
        s += _tri_area(node(m,m.tris[1,t]), node(m,m.tris[2,t]),
                       node(m,m.tris[3,t]))
    end
    s
end

function _bnd_area(m)
    bf, _ = boundary_faces(m.tets)
    s = 0.0
    for f in bf
        s += _tri_area(node(m,f[1]), node(m,f[2]), node(m,f[3]))
    end
    s
end

@testset "Dense curved PLC converges under auto maxpts" begin
    # Sentinel resolution: 0 → max(6000, 8·deduplicated input vertices);
    # explicit values pass through unchanged.
    lims = RC._recovery_limits("recover_boundary_cdt")
    @test lims.points == 0
    @test RC._resolve_recovery_limits(lims, 8).points == 6000
    @test RC._resolve_recovery_limits(lims, 1233).points == 9864
    @test RC._resolve_recovery_limits(
        RC._recovery_limits("recover_boundary_cdt"; maxpts=999), 2000).points == 999

    plc = _load_plc_txt(joinpath(@__DIR__, "..", "fixtures",
                                 "plc_helical_pipe_wall_360.txt"))
    @test nnodes(plc) == 182
    @test ntris(plc) == 360
    @test validate(plc).ok

    m = recover_boundary_cdt(plc)
    @test validate(m).ok
    @test is_closed_manifold(m)
    # The recovery is deterministic: 182 input + 31 exact Steiner points.
    @test nnodes(m) == 213
    @test ntets(m) == 665
    # Boundary Steiner subdivision covers the input surface exactly.
    @test _bnd_area(m) ≈ _surf_area(plc) rtol=1e-12

    # Explicit caps remain literal hard bounds: the PLC genuinely needs the
    # 31 Steiner points above, so maxpts=200 must fail cleanly (not silently
    # truncate the boundary).
    @test_throws ArgumentError recover_boundary_cdt(plc; maxpts=200)
end
