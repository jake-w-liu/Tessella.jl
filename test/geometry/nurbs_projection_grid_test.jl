# The suite includes every test file into Main; this file keeps its names to itself.
module NURBSProjectionGridTests
using Test
using Tessella
import Tessella.Model as M
import Tessella.NURBS as N
const A = Tessella.API

# Projection onto a tensor patch follows OpenCASCADE's Extrema_GenExtPS and
# evaluates the patch through a polynomial cache as BSplSLib does. Three things
# keep their results while costing less: the evaluator works in fixed space,
# the sample grid of a patch is built once, and an edge too far from a patch to
# cross it is not searched.

_lcg(seed) = let state = UInt64(seed)
    () -> (state = state * 6364136223846793005 + 1442695040888963407; (state >> 11) / 9007199254740992.0)
end

function _random_surface(next, pu, pv, nu, nv, mode)
    knots(n, p) = vcat(zeros(p + 1), sort!([next() for _ in 1:n-p-1]), ones(p + 1))
    controls = [(i + 0.3 * next(), j + 0.3 * next(), 0.5 * sin(0.7i) * cos(0.9j) + 0.2 * next()) for i in 1:nu, j in 1:nv]
    weights = mode === :plain ? ones(nu, nv) :
              mode === :rational ? [0.5 + next() for _ in 1:nu, _ in 1:nv] :
              [i <= 3 && j <= 3 ? 0.5 + next() : 1.0 for i in 1:nu, j in 1:nv]
    return N.NURBSSurface(pu, pv, knots(nu, pu), knots(nv, pv), controls, weights)
end

# The evaluator as it was written first: lists and dictionaries per call.
function _reference_evaluate(cache, u, v, order)
    M._model_nurbs_cache_prepare!(cache, u, v)
    inverse = (1.0 / cache.half_spans[1], 1.0 / cache.half_spans[2])
    parameter = order == 0 ?
        ((u - cache.centers[1]) / cache.half_spans[1], (v - cache.centers[2]) / cache.half_spans[2]) :
        ((u - cache.centers[1]) * inverse[1], (v - cache.centers[2]) * inverse[2])
    min_parameter, max_parameter = cache.swapped ? reverse(parameter) : parameter
    width = cache.global_rational ? 4 : 3
    max_values = M._model_nurbs_cache_polynomial(cache.coefficients, max_parameter, order)
    min_degree = cache.swapped ? cache.patch.surface.degree_v : cache.patch.surface.degree_u
    homogeneous = Dict{NTuple{2,Int},Vector{Float64}}()
    for b in 0:order
        pieces = Vector{Float64}[max_values[b+1][a*width+1:(a+1)*width] for a in 0:min_degree]
        min_values = M._model_nurbs_cache_polynomial(pieces, min_parameter, order - b)
        for a in 0:order-b
            homogeneous[a, b] = min_values[a+1]
        end
    end
    if order == 0
        point = homogeneous[0, 0]
        return cache.global_rational ? (point[1] / point[4], point[2] / point[4], point[3] / point[4]) :
               (point[1], point[2], point[3])
    end
    result = Dict{NTuple{2,Int},NTuple{3,Float64}}()
    if cache.global_rational
        inverse_weight = 1.0 / homogeneous[0, 0][4]
        zero_h = (0.0, 0.0, 0.0, 0.0)
        for b in 0:order, a in 0:order
            h = get(homogeneous, (a, b), zero_h)
            value = Float64[h[1], h[2], h[3]]
            for pp in 0:b-1, qq in 0:a
                coefficient = binomial(b, pp) * binomial(a, qq) * get(homogeneous, (a - qq, b - pp), zero_h)[4]
                prior = result[qq, pp]
                for component in 1:3
                    value[component] -= coefficient * prior[component]
                end
            end
            for qq in 0:a-1
                coefficient = binomial(a, qq) * homogeneous[a-qq, 0][4]
                prior = result[qq, b]
                for component in 1:3
                    value[component] -= coefficient * prior[component]
                end
            end
            result[a, b] = (value[1] * inverse_weight, value[2] * inverse_weight, value[3] * inverse_weight)
        end
    else
        for (key, value) in homogeneous
            result[key] = (value[1], value[2], value[3])
        end
    end
    u_slot = cache.swapped ? (0, 1) : (1, 0)
    v_slot = cache.swapped ? (1, 0) : (0, 1)
    point = result[0, 0]
    du = ntuple(i -> result[u_slot][i] * inverse[1], 3)
    dv = ntuple(i -> result[v_slot][i] * inverse[2], 3)
    order == 1 && return (point, du, dv)
    uu_slot = cache.swapped ? (0, 2) : (2, 0)
    vv_slot = cache.swapped ? (2, 0) : (0, 2)
    duu = ntuple(i -> result[uu_slot][i] * (inverse[1] * inverse[1]), 3)
    dvv = ntuple(i -> result[vv_slot][i] * (inverse[2] * inverse[2]), 3)
    duv = ntuple(i -> result[1, 1][i] * (inverse[1] * inverse[2]), 3)
    return (point, du, dv, duu, dvv, duv)
end

const _CASES = ((3, 3, 7, 6, :plain), (2, 3, 6, 7, :plain), (3, 2, 7, 5, :plain), (1, 3, 4, 6, :plain),
                (3, 3, 7, 6, :rational), (3, 2, 7, 5, :rational), (2, 2, 5, 5, :rational),
                (3, 3, 9, 8, :partly), (3, 2, 9, 7, :partly), (5, 4, 8, 7, :plain), (4, 5, 7, 9, :rational))

@testset "The patch evaluator in fixed space gives the first evaluator's numbers" begin
    next = _lcg(0x2545f4914f6cdd1d)
    seen = Set{Tuple{Bool,Bool}}()
    for (pu, pv, nu, nv, mode) in _CASES
        patch = M._ModelNURBSSurface(:bspline, _random_surface(next, pu, pv, nu, nv, mode))
        cache = M._ModelNURBSProjectionCache(patch); reference = M._ModelNURBSProjectionCache(patch)
        push!(seen, (cache.global_rational, pu > pv))
        same = true
        for _ in 1:60, order in 0:2
            u = next(); v = next()
            same &= M._model_nurbs_cache_evaluate(cache, u, v, order) === _reference_evaluate(reference, u, v, order)
        end
        @test same
        # Derivatives after a value, and a value after derivatives, in one span.
        u = 0.37; v = 0.41
        @test M._model_nurbs_cache_evaluate(cache, u, v, 2) === _reference_evaluate(reference, u, v, 2)
        @test M._model_nurbs_cache_evaluate(cache, u, v, 0) === _reference_evaluate(reference, u, v, 0)
        @test M._model_nurbs_cache_evaluate(cache, u, v, 1) === _reference_evaluate(reference, u, v, 1)
        if Base.JLOptions().compile_enabled == 1
            probe(c, order) = @allocated M._model_nurbs_cache_evaluate(c, 0.37, 0.41, order)
            probe(cache, 0); probe(cache, 2)
            @test probe(cache, 0) == 0
            @test probe(cache, 2) <= 256      # the six returned vectors, nothing per coefficient
        end
    end
    @test seen == Set([(false, false), (false, true), (true, false), (true, true)])
    @test_throws ArgumentError M._model_nurbs_cache_evaluate(
        M._ModelNURBSProjectionCache(M._ModelNURBSSurface(:bspline, _random_surface(next, 2, 2, 4, 4, :plain))), 0.5, 0.5, 3)
end

@testset "A patch's sample grid is built once and read for every point" begin
    next = _lcg(0x9e3779b97f4a7c15)
    surface = _random_surface(next, 3, 3, 7, 6, :rational)
    patch = M._ModelNURBSSurface(:bspline, surface)
    target(k) = (2.0 + 0.5k, 3.0 - 0.3k, 0.4 + 0.1k)
    delete!(M._MODEL_NURBS_PROJECTION_GRIDS, surface.controls)
    first_results = [M._model_nurbs_projection_cached(patch, target(k)) for k in 0:4]
    grid = M._MODEL_NURBS_PROJECTION_GRIDS[surface.controls]
    @test size(grid.points) == (length(grid.us), length(grid.vs))
    @test length(grid.us) >= 44 && length(grid.vs) >= 44
    @test all(grid.points[i, j] === M._model_nurbs_cache_evaluate(M._ModelNURBSProjectionCache(patch), grid.us[i], grid.vs[j], 0)
              for i in (1, 7, length(grid.us)), j in (1, 5, length(grid.vs)))
    # The same answers from the stored grid, and from a grid built again.
    again = [M._model_nurbs_projection_cached(patch, target(k)) for k in 0:4]
    @test M._MODEL_NURBS_PROJECTION_GRIDS[surface.controls] === grid
    @test again == first_results
    delete!(M._MODEL_NURBS_PROJECTION_GRIDS, surface.controls)
    @test [M._model_nurbs_projection_cached(patch, target(k)) for k in 0:4] == first_results
    # The projection is a foot: the residual is normal to the patch there.
    for (k, result) in enumerate(first_results)
        point, du, dv = M._model_nurbs_cache_evaluate(M._ModelNURBSProjectionCache(patch), result.uv[1], result.uv[2], 1)
        residual = target(k - 1) .- point
        inside = 1e-6 < result.uv[1] < 1 - 1e-6 && 1e-6 < result.uv[2] < 1 - 1e-6
        inside && @test abs(sum(residual .* du)) <= 1e-7 * sqrt(sum(abs2, du)) * (1 + sqrt(sum(abs2, residual)))
        inside && @test abs(sum(residual .* dv)) <= 1e-7 * sqrt(sum(abs2, dv)) * (1 + sqrt(sum(abs2, residual)))
        @test result.xyz === point
    end
    # A surface with the same control net and other knots does not read this grid.
    other = N.NURBSSurface(3, 3, vcat(zeros(4), [0.2, 0.5, 0.9], ones(4)), vcat(zeros(4), [0.3, 0.6], ones(4)),
                           copy(surface.controls), copy(surface.weights))
    other_patch = M._ModelNURBSSurface(:bspline, other)
    mixed = M._model_nurbs_projection_cached(other_patch, target(1))
    delete!(M._MODEL_NURBS_PROJECTION_GRIDS, other.controls)
    @test M._model_nurbs_projection_cached(other_patch, target(1)) == mixed
    @test mixed != first_results[2]
    # Candidates are kept in one work space per task: a smaller grid after a
    # larger one, and the larger one again, read only their own.
    small = M._ModelNURBSSurface(:bspline, _random_surface(next, 2, 2, 4, 4, :plain))
    small_result = M._model_nurbs_projection_cached(small, (2.0, 2.0, 1.0))
    @test M._model_nurbs_projection_cached(patch, target(2)) == first_results[3]
    @test M._model_nurbs_projection_cached(small, (2.0, 2.0, 1.0)) == small_result
end

@testset "An edge far from a patch is not searched, a crossing one is found" begin
    A.initialize()
    try
        A.option("General.Terminal", 0)
        occ = A.model.occ
        control = Int[]
        for j in 1:4, i in 1:4
            push!(control, occ.add_point(-1 + 2 * (i - 1) / 3, -1 + 2 * (j - 1) / 3, 0.2 * sin(1.3i) * cos(0.8j)))
        end
        face = occ.add_bspline_surface(control, 4; degree_u=3, degree_v=3)
        occ.synchronize()
        A._with_model() do m
            sg = m.surface_geometry[face]
            @test M._model_nurbs_surface(sg)
            line(o, d) = (kind=:line, o=o, d=d)
            # A vertical line through the patch crosses it once, at the patch's height there.
            crossing = line((0.1, -0.2, -1.0), (0.0, 0.0, 1.0))
            roots = M._brep_edge_pierce_tensor(sg, crossing, 0.0, 2.0, 1e-7)
            @test length(roots) == 1
            u, v, xyz = M._model_nurbs_surface_project(sg, M._brep_curve_eval(crossing, roots[1]))
            @test abs(xyz[1] - 0.1) < 1e-9 && abs(xyz[2] + 0.2) < 1e-9
            @test abs(M._brep_curve_eval(crossing, roots[1])[3] - xyz[3]) < 1e-9
            # The same line shifted far along z misses the box of the control points.
            far = line((0.1, -0.2, 5.0), (0.0, 0.0, 1.0))
            @test isempty(M._brep_edge_pierce_tensor(sg, far, 0.0, 2.0, 1e-7))
            if Base.JLOptions().compile_enabled == 1
                probe(g, c) = @allocated M._brep_edge_pierce_tensor(g, c, 0.0, 2.0, 1e-7)
                probe(sg, far)
                @test probe(sg, far) <= 256            # no projection was made
            end
            # A line beside the patch but inside the box's margin is still searched, and misses.
            beside = line((1.0 + 2e-5, 0.0, -1.0), (0.0, 0.0, 1.0))
            @test isempty(M._brep_edge_pierce_tensor(sg, beside, 0.0, 2.0, 1e-7))
        end
    finally
        A.finalize()
    end
    # Each kind of curve record stays inside its box.
    line = (kind=:line, o=(1.0, -2.0, 0.5), d=(0.3, 0.7, -0.2))
    circle = (kind=:circle, c=(0.5, 0.25, -1.0), x=(0.6, 0.8, 0.0), y=(-0.8, 0.6, 0.0), r=2.5)
    for (curve, lo, hi) in ((line, -1.5, 2.0), (circle, 0.0, 2pi))
        lower, upper = M._brep_curve_box(curve, lo, hi)
        @test all(all(lower .- 1e-12 .<= M._brep_curve_eval(curve, lo + k * (hi - lo) / 64) .<= upper .+ 1e-12) for k in 0:64)
    end
    @test M._brep_curve_box((kind=:ellipse,), 0.0, 1.0) === nothing
end
end
