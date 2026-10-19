# Mesh.HighOrderOptimize — curve-owned parametric boundary sliding (OPT01).
#
# Correctness : a Tri6 edge mid owned by a curve entity slides along its
#               stored [u] parametrization — candidates are re-evaluations of
#               the owning curve (`_P2ParametricSupport.lift`), never
#               free-space displacement, so a unit-circle mid keeps rho == 1
#               bit-for-bit.
# Acceptance  : a slide commits only when every incident element's exact
#               Bernstein positivity certificate survives AND the sampled
#               minimum scaled Jacobian strictly improves; the stored
#               parameter and the coordinates commit/rollback atomically.
#               Nodes also referenced by uncertified non-curve blocks, or
#               owned by dim-2 surface entities, still stay fixed.
# API         : HighOrderOptimize=1 slides a deliberately mis-parameterized
#               OCC circle-arc mid back toward the quality optimum on
#               set_order(2) and writes the accepted parameter back through
#               get_node; =0 leaves the whole cache bitwise identical.

using Test
using LinearAlgebra: norm, cross
import Tessella
using Tessella.HighOrder
const API = Tessella.API

# ── independent Tri6 scaled-Jacobian oracle (loop-form gradients on the
#    r,s >= 0, r+s <= 1 reference triangle, denser lattice than the engine's
#    interpolation samples) ────────────────────────────────────────────────
function _slide_tri_grads(r, s)
    L1, L2, L3 = 1.0 - r - s, r, s
    dr = (1 - 4L1, 4L2 - 1, 0.0, 4(L1 - L2), 4L3, -4L3)
    ds = (1 - 4L1, 0.0, 4L3 - 1, -4L2, 4L2, 4(L1 - L3))
    return dr, ds
end
function _slide_tri_scaled(coords, tri6, t, r, s)
    dr, ds = _slide_tri_grads(r, s)
    c1 = zeros(3); c2 = zeros(3)
    for k in 1:6
        v = tri6[k, t]
        for d in 1:3
            c1[d] += dr[k] * coords[d, v]
            c2[d] += ds[k] * coords[d, v]
        end
    end
    return norm(cross(c1, c2)) / (norm(c1) * norm(c2))
end
function _slide_tri_min(coords, tri6, ts = axes(tri6, 2))
    mn = Inf
    grid = range(0.0, 1.0; length = 5)
    for t in ts, r in grid, s in grid
        r + s > 1 && continue
        w = _slide_tri_scaled(coords, tri6, t, r, s)
        w < mn && (mn = w)
    end
    return mn
end

# Unit quarter-arc parametrization shared by the engine fixture and its
# lift: u ∈ [0,1] maps to angle θ = π(1-u)/2, i.e. u=0 → (1,0,0),
# u=1 → (0,1,0).  The arc midpoint is u = 0.5.
_slide_arc(u) = (cos(π * (1 - u) / 2), sin(π * (1 - u) / 2), 0.0)

# One Tri6 whose (2,3) edge mid is the curve-owned node: corners
# (0,0,0)–(1,0,0)–(0,1,0), straight mids on edges (1,2)/(3,1), arc mid at
# parameter `u0`.  Returns (coords, tri6, mid_node).
function _slide_arc_tri(u0)
    coords = zeros(3, 6)
    coords[:, 1] .= (0.0, 0.0, 0.0)
    coords[:, 2] .= (1.0, 0.0, 0.0)
    coords[:, 3] .= (0.0, 1.0, 0.0)
    coords[:, 4] .= (0.5, 0.0, 0.0)
    p = _slide_arc(u0)
    coords[:, 5] .= (p[1], p[2], p[3])
    coords[:, 6] .= (0.0, 0.5, 0.0)
    tri6 = reshape(Int32[1, 2, 3, 4, 5, 6], 6, 1)
    return coords, tri6, Int32(5)
end
_slide_support(u0; lo = 0.0, hi = 1.0) =
    HighOrder._P2ParametricSupport(Float64[u0], u -> _slide_arc(u[1]),
                                   Float64[lo], Float64[hi])

@testset "HighOrderOptimize curve parametric sliding (engine)" begin
    @testset "off-center arc mid slides on the curve and improves quality" begin
        u0 = 0.30
        coords, tri6, mid = _slide_arc_tri(u0)
        before = copy(coords)
        @test _slide_tri_min(coords, tri6) > 0          # certified start
        jac0 = _slide_tri_min(coords, tri6)
        support = _slide_support(u0)
        sliders = Dict{Int32,HighOrder._P2ParametricSupport}(mid => support)
        # Every node masked off — only the slider descriptor admits a move.
        moved = optimize_p2_interior!(coords, ((Int32(9), tri6),),
                                      falses(size(coords, 2));
                                      sliders = sliders,
                                      itermax = 60, passes = 10)
        @test moved[1] >= 1
        @test support.u[1] != u0                        # u state advanced
        # The committed coordinates are the curve evaluation at the new
        # parameter — the node stays exactly on the unit arc (rho == 1).
        p = _slide_arc(support.u[1])
        @test coords[:, mid] == collect(p)
        @test hypot(coords[1, mid], coords[2, mid]) ≈ 1.0 atol = 1e-14
        @test _slide_tri_min(coords, tri6) > jac0       # strict improvement
        @test HighOrder._p2opt_certify(Val(9), coords, tri6, 1)
        @test _slide_tri_min(coords, tri6) > 0
        for v in axes(coords, 2)
            v == Int(mid) && continue
            @test coords[:, v] == before[:, v]          # untouched bitwise
        end
    end

    @testset "no slider and full-fixed mask is a bitwise no-op" begin
        coords, tri6, _ = _slide_arc_tri(0.30)
        before = copy(coords)
        moved = optimize_p2_interior!(coords, ((Int32(9), tri6),),
                                      falses(size(coords, 2));
                                      itermax = 60, passes = 10)
        @test moved[1] == 0
        @test coords == before
    end

    @testset "degenerate parameter bounds pin the node bitwise" begin
        coords, tri6, mid = _slide_arc_tri(0.30)
        before = copy(coords)
        support = _slide_support(0.30; lo = 0.30, hi = 0.30)
        sliders = Dict{Int32,HighOrder._P2ParametricSupport}(mid => support)
        moved = optimize_p2_interior!(coords, ((Int32(9), tri6),),
                                      falses(size(coords, 2));
                                      sliders = sliders,
                                      itermax = 60, passes = 10)
        @test moved[1] == 0
        @test coords == before
        @test support.u == [0.30]
    end

    @testset "Seg3 sibling on the owning curve does not block sliding" begin
        # The same mid node also terminates the owning curve's own Seg3
        # discretization — an uncertified 1-D block shares the parameter
        # state rather than conflicting with it.
        u0 = 0.30
        coords, tri6, mid = _slide_arc_tri(u0)
        seg3 = reshape(Int32[2, 3, mid], 3, 1)
        before = copy(coords)
        support = _slide_support(u0)
        sliders = Dict{Int32,HighOrder._P2ParametricSupport}(mid => support)
        jac0 = _slide_tri_min(coords, tri6)
        moved = optimize_p2_interior!(coords,
                                      ((Int32(9), tri6), (Int32(8), seg3)),
                                      falses(size(coords, 2));
                                      sliders = sliders,
                                      itermax = 60, passes = 10)
        @test moved[1] >= 1
        @test _slide_tri_min(coords, tri6) > jac0
        @test hypot(coords[1, mid], coords[2, mid]) ≈ 1.0 atol = 1e-14
        for v in axes(coords, 2)
            v == Int(mid) && continue
            @test coords[:, v] == before[:, v]
        end
    end

    @testset "uncertified 2-D block sharing the node keeps it fixed" begin
        # A plain Tri3 (no Bernstein certificate) referencing the mid turns
        # it into a cross-support node — sliding it would be uncertifiable.
        coords, tri6, mid = _slide_arc_tri(0.30)
        tri3 = reshape(Int32[1, 2, mid], 3, 1)
        before = copy(coords)
        support = _slide_support(0.30)
        sliders = Dict{Int32,HighOrder._P2ParametricSupport}(mid => support)
        moved = optimize_p2_interior!(coords,
                                      ((Int32(9), tri6), (Int32(2), tri3)),
                                      falses(size(coords, 2));
                                      sliders = sliders,
                                      itermax = 60, passes = 10)
        @test moved[1] == 0
        @test coords == before
        @test support.u == [0.30]
    end

    @testset "rejected candidates leave u and coordinates bit-for-bit" begin
        # Bounds box clamped to the starting neighbourhood: every probe
        # either leaves the box or produces no strict improvement — whatever
        # the optimizer tries, a rejected pass restores both stores exactly.
        u0 = 0.45
        coords, tri6, mid = _slide_arc_tri(u0)
        before = copy(coords)
        support = _slide_support(u0; lo = u0, hi = u0)
        moved = optimize_p2_interior!(coords, ((Int32(9), tri6),),
                                      falses(size(coords, 2));
                                      sliders = Dict(mid => support),
                                      itermax = 60, passes = 10)
        @test moved[1] == 0
        @test coords == before
        @test support.u == [u0]
    end
end

@testset "HighOrderOptimize curve sliding (API)" begin
    @testset "OCC quarter-disk arc mid slides along the circle" begin
        api = API
        api.initialize()
        try
            api.option("Mesh.HighOrderOptimize", 1)
            api.option("Mesh.MeshSizeMax", 0.7)
            # Quarter disk: center + two axial endpoints, center-form circle
            # arc plus two radial lines closing the boundary loop.
            pc = api.model.occ.add_point(0.0, 0.0, 0.0)
            p1 = api.model.occ.add_point(1.0, 0.0, 0.0)
            p2 = api.model.occ.add_point(0.0, 1.0, 0.0)
            arc = api.model.occ.add_circle_arc(p1, pc, p2; center = true)
            l1 = api.model.occ.add_line(p2, pc)
            l2 = api.model.occ.add_line(pc, p1)
            loop = api.model.occ.add_curve_loop([arc, l1, l2])
            surf = api.model.occ.add_plane_surface([loop])
            api.model.occ.synchronize()
            api.mesh.generate(2)
            api.mesh.set_order(2)

            # The curve-owned mids shared with a Tri6 edge: the third slot
            # of each Seg3 on the arc is the mid node, and boundary Tri6
            # elements reference the same node index.
            types_c, _, enodes_c = api.mesh.get_elements(1, arc)
            seg_mids = Set{UInt64}()
            for (ti, ty) in enumerate(types_c)
                ty == 8 || continue
                flat = enodes_c[ti]
                for j in 1:length(flat)÷3
                    push!(seg_mids, flat[3j])
                end
            end
            @test !isempty(seg_mids)
            function curve_mid()
                for t in seg_mids
                    p, uv, dim, owner = api.mesh.get_node(t)
                    if dim == 1 && owner == arc && length(uv) == 1
                        return (t, p, uv)
                    end
                end
                return nothing
            end
            found = curve_mid()
            @test found !== nothing
            tag, p0, uv0 = found
            u0 = uv0[1]
            @test hypot(p0[1], p0[2]) ≈ 1.0 atol = 1e-12

            lo, hi = api.model.get_parametrization_bounds(1, arc)
            lo = lo[1]; hi = hi[1]
            span = hi - lo
            @test span > 0
            # Deliberately mis-parameterize: 4% of the range toward the
            # farther bound — enough to depress the incident min scaled
            # Jacobian while staying inside the [0.1, 2) optimization band —
            # coordinates re-evaluated through the public curve evaluator.
            u_bad = u0 + (u0 - lo < hi - u0 ? 0.04 : -0.04) * span
            q_bad = api.model.get_value(1, arc, [u_bad])
            api.mesh.set_node(tag, q_bad, [u_bad])
            p1v, uv1, dim1, own1 = api.mesh.get_node(tag)
            @test uv1 == [u_bad]
            @test p1v == q_bad

            # Re-running set_order(2) re-plans through the optimizer: the
            # mid slides back toward the quality optimum along the arc.
            api.mesh.set_order(2)
            p2v, uv2, dim2, own2 = api.mesh.get_node(tag)
            @test dim2 == 1 && own2 == arc
            @test length(uv2) == 1
            u_new = uv2[1]
            @test u_new != u_bad
            @test abs(u_new - u0) < abs(u_bad - u0) ||
                  abs(u_new - 0.5 * (lo + hi)) <= abs(u_bad - 0.5 * (lo + hi))
            @test lo <= u_new <= hi
            @test hypot(p2v[1], p2v[2]) ≈ 1.0 atol = 1e-12
            # Accepted parameter and coordinates stay consistent: the
            # stored coordinates are the curve evaluation at u_new.
            q_new = api.model.get_value(1, arc, [u_new])
            @test p2v ≈ q_new atol = 1e-12
        finally
            api.option("Mesh.HighOrderOptimize", 0)
            api.option("Mesh.MeshSizeMax", 1e22)
            api.finalize()
        end
    end

    @testset "OCC quarter-disk with =0 leaves nodes bitwise identical" begin
        api = API
        api.initialize()
        try
            @test api.option("Mesh.HighOrderOptimize") == 0
            api.option("Mesh.MeshSizeMax", 0.7)
            pc = api.model.occ.add_point(0.0, 0.0, 0.0)
            p1 = api.model.occ.add_point(1.0, 0.0, 0.0)
            p2 = api.model.occ.add_point(0.0, 1.0, 0.0)
            arc = api.model.occ.add_circle_arc(p1, pc, p2; center = true)
            l1 = api.model.occ.add_line(p2, pc)
            l2 = api.model.occ.add_line(pc, p1)
            loop = api.model.occ.add_curve_loop([arc, l1, l2])
            api.model.occ.add_plane_surface([loop])
            api.model.occ.synchronize()
            api.mesh.generate(2)
            api.mesh.set_order(2)
            tags, coords_f, _ = api.mesh.get_nodes()
            snapshot = copy(coords_f)
            api.mesh.set_order(2)
            _, coords_f2, _ = api.mesh.get_nodes()
            @test coords_f2 == snapshot
        finally
            api.option("Mesh.MeshSizeMax", 1e22)
            api.finalize()
        end
    end
end
