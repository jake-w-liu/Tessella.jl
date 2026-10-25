# ══════════════════════════════════════════════════════════════════════════════
# Mesh.HighOrderOptimize — bounded interior high-order-node relocation
# (OPT01 slice).
#
# Upstream semantics (gmsh 4.15.2): 0 disables high-order optimization,
# -1 applies elastic smoothing only, 1 optimizes (relocates interior
# high-order nodes to improve the minimum scaled Jacobian while keeping
# boundary/geometry supports fixed), 2 applies both.  The optimization
# stage honors `Mesh.HighOrderThresholdMin/Max` — only nodes whose
# incident minimum lies inside the band are candidates.  The elastic
# stage runs the same certificate-gated bounded ascent without the band
# (`threshold_min=-Inf`, `threshold_max=Inf`), so elements below the
# optimize floor are smoothed too — the native realization of upstream's
# elasticity-analogy pass, restricted to the certified P2 families and
# still never able to create or repair inverted elements.
#
# The pass itself lives in HighOrder.optimize_p2_interior! and covers every
# full-P2 family with an exact Bernstein positivity bound — Tri6, Quad9,
# Tet10, Hex27, Prism18 and Pyr14 (the collapsed-map certificate of
# HighOrderPyramid.jl supplies its boolean bound) — while every uncertified
# block family keeps its nodes fixed.  Curve-owned mids additionally slide along their
# owning CAD curve in parameter space and face-owned mids slide in the
# owning surface's (u,v) parameters — both re-evaluated through the same
# model evaluators that produced the stored parameters.  This file
# translates option state, classification ownership, and
# geometric-support structure into the engine's movable mask, fixed
# edge/facet sets, and slider descriptors, then commits accepted parameter
# states back into the plan-local entity_data/classification atomically.
# It mutates the *unpublished* coordinate matrix of a freshly built
# plan/overlay, so a failure still leaves the previously published cache
# untouched.
# ══════════════════════════════════════════════════════════════════════════════

_high_order_optimize_level() = begin
    v = get(OPTIONS, "Mesh.HighOrderOptimize", 0.0)
    (v isa Real && isfinite(v)) || return 0
    return clamp(round(Int, v), typemin(Int32), typemax(Int32))
end

# Stage gates: the unbounded elastic-smoothing stage runs for upstream
# modes -1 (elastic only) and 2 (elastic then optimize); the banded
# optimization stage runs for modes 1 and 2.
_high_order_optimize_elastic() =
    (l = _high_order_optimize_level(); l == -1 || l == 2)
_high_order_optimize_optimize() =
    (l = _high_order_optimize_level(); l == 1 || l == 2)
_high_order_optimize_relocation() =
    _high_order_optimize_elastic() || _high_order_optimize_optimize()

function _high_order_optimize_int_option(name::AbstractString, default::Integer)
    v = get(OPTIONS, name, Float64(default))
    (v isa Real && isfinite(v)) || return Int(default)
    return clamp(round(Int, v), 0, typemax(Int32))
end

function _high_order_optimize_num_option(name::AbstractString, default::Real)
    v = get(OPTIONS, name, Float64(default))
    (v isa Real && isfinite(v)) || return Float64(default)
    return Float64(v)
end

# Candidate mask: a node may move freely only when it is owned by a surface
# or volume entity and carries no stored parametric coordinates (curve, point
# and parametric-surface supports keep their projected parameters).  Curve
# and surface owners instead enter through the `sliders` map built by
# `_high_order_optimize_sliders`, which relocates them in the owning entity's
# parameter space.  When the classification does not cover the mesh,
# ownership cannot be proven and the mask stays permissive — the engine's own
# mid-slot, single-family and boundary-support rules still keep corners and
# boundary supports fixed.
function _high_order_optimize_movable(mesh::MixedMesh, class)
    nn = size(mesh.coords, 2)
    movable = falses(nn)
    classified = class !== nothing &&
        length(class.node_entities) == nn &&
        length(class.node_parametric) == nn
    for v in 1:nn
        movable[v] = !classified ||
            (class.node_entities[v][1] in (2, 3) &&
             class.node_parametric[v] === nothing)
    end
    return movable
end

# Parametric sliders for `optimize_p2_interior!`: every node with owner
# dimension 1 and a stored finite `[u]` parameter gets a
# `_P2ParametricSupport` whose lift is the same
# `Model._model_high_order_curve_geometry`/`_model_high_order_curve_point`
# evaluator that produced the node's coordinates, clamped to the curve's
# native parameter bounds; every node with owner dimension 2 and a stored
# finite `[u,v]` pair gets a support whose lift is `model_value` — the
# same surface evaluation family `model_parametrization` inverted when the
# parameter was stored — clamped to `model_parametrization_bounds`.
# Returns `Dict{Int32,_P2ParametricSupport}` or nothing.
function _high_order_optimize_sliders(mesh::MixedMesh, class, model, caller)
    (class === nothing || model === nothing) && return nothing
    nn = size(mesh.coords, 2)
    entities = class.node_entities
    parameters = class.node_parametric
    (length(entities) == nn && length(parameters) == nn) || return nothing
    sliders = Dict{Int32,HighOrder._P2ParametricSupport}()
    for v in 1:nn
        stored = parameters[v]
        stored === nothing && continue
        owner = entities[v]
        if owner[1] == 1
            stored isa Vector{Float64} || continue
            length(stored) == 1 || continue
            isfinite(stored[1]) || continue
            tag = Int(owner[2])
            curve = try
                Model._model_high_order_curve_geometry(model, tag, caller)
            catch err
                (err isa InterruptException || err isa OutOfMemoryError) &&
                    rethrow()
                nothing
            end
            curve === nothing && continue
            (isfinite(curve.lower) && isfinite(curve.upper) &&
             curve.lower <= stored[1] <= curve.upper) || continue
            support = HighOrder._P2ParametricSupport(
                stored,                      # shared: replaced, never mutated
                let geo = curve
                    u -> begin
                        p = Model._model_high_order_curve_point(
                            geo, u[1], caller)
                        (p[1], p[2], p[3])
                    end
                end,
                Float64[curve.lower], Float64[curve.upper])
            sliders[Int32(v)] = support
        elseif owner[1] == 2
            support = _high_order_optimize_surface_support(
                model, Int(owner[2]), stored)
            support === nothing || (sliders[Int32(v)] = support)
        end
    end
    return isempty(sliders) ? nothing : sliders
end

# Surface-owned parametric slider: a node classified on a dim-2 entity with
# a stored finite `[u,v]` parameter inside the face's native parameter
# bounds gets a `_P2ParametricSupport` whose lift re-evaluates the face
# through `model_value` — the same evaluator family that produced the
# stored parameter via `model_parametrization`/`model_closest_point`
# (OCC analytic and NURBS patch records, built-in ruled surfaces honoring
# the `Geometry.OldRuledSurface` runtime selector, native plane frames and
# discrete surfaces).  A face that cannot supply finite bounds — or whose
# evaluator cannot lift the stored parameter to a finite point — leaves
# the node fixed; mid-evaluation failures degrade to rejected probes
# rather than aborting the pass.
function _high_order_optimize_surface_support(model, tag::Int, stored)
    stored isa Vector{Float64} || return nothing
    length(stored) == 2 || return nothing
    (isfinite(stored[1]) && isfinite(stored[2])) || return nothing
    bounds = try
        Model.model_parametrization_bounds(model, 2, tag)
    catch err
        (err isa InterruptException || err isa OutOfMemoryError) && rethrow()
        nothing
    end
    bounds === nothing && return nothing
    lower, upper = bounds
    (length(lower) == 2 && length(upper) == 2) || return nothing
    for i in 1:2
        (isfinite(lower[i]) && isfinite(upper[i]) &&
         lower[i] <= stored[i] <= upper[i]) || return nothing
    end
    lift = let m = model, face = tag
        function (u)
            values = try
                Model.model_value(m, 2, face, u;
                    old_ruled_surface =
                        !iszero(OPTIONS["Geometry.OldRuledSurface"]))
            catch err
                (err isa InterruptException || err isa OutOfMemoryError) &&
                    rethrow()
                return nothing
            end
            (length(values) == 3 && all(isfinite, values)) ||
                return nothing
            return (values[1], values[2], values[3])
        end
    end
    # The stored parameter must already lift to a finite surface point —
    # an evaluator that cannot produce one leaves the node fixed.
    lift(stored) === nothing && return nothing
    return HighOrder._P2ParametricSupport(
        stored,                            # shared: replaced, never mutated
        lift,
        Float64[lower[1], lower[2]], Float64[upper[1], upper[2]])
end

# Write accepted slider parameters back into the plan-local stores: the mesh
# `entity_data` (rebuilt — the MixedMesh constructor revalidates each stored
# parameter vector against its owner dimension) and the classification's
# `node_parametric` (a fresh vector so a shared rebind can never mutate the
# published class).  `slid` maps node index to the accepted `Vector{Float64}`.
function _high_order_optimize_commit_sliders(mesh, class, slid)
    nn = size(mesh.coords, 2)
    data = mesh.entity_data
    if data !== nothing && length(data.node_parametric) == nn
        updated = collect(data.node_parametric)
        for (v, u) in slid
            updated[v] = u
        end
        mesh = _mixed_rebuild_metadata(mesh, mesh.blocks;
            node_parametric=updated, keep_empty_blocks=true)
    end
    if class !== nothing && length(class.node_parametric) == nn
        updated = collect(class.node_parametric)
        for (v, u) in slid
            updated[v] = u
        end
        pclass = class
        class = _mixed_rebind_class(class, mesh; node_parametric=updated)
        if class !== nothing
            # The rebind rebuilds classification vectors only; the
            # edge-keyed mid-support tables survive the identical-index
            # copy — carry them over and refresh the entries of the nodes
            # whose parameter state was just accepted.
            merge!(class.mid_parametric, pclass.mid_parametric)
            merge!(class.mid_keys, pclass.mid_keys)
            for (v, u) in slid
                key = get(class.mid_keys, Int32(v), nothing)
                key === nothing || (class.mid_parametric[key] = u)
            end
        end
    end
    return mesh, class
end

# Per-owner-group boundary supports of every certified family: primary edges
# on a group's boundary facet whose mid nodes stay fixed, plus the sorted
# facet keys for the face-center test.  Grouping by element owner makes
# faces shared between two owners (internal interfaces of a conformal
# multi-entity mesh) count as each group's boundary even though they are
# interior to the combined complex — those mids carry interface supports.
# For surface families the facets are element edges, so the same loop marks
# the owner-group boundary wires.  Returns `(edges, facet_keys)`.
function _high_order_optimize_fixed_supports(mesh::MixedMesh, class)
    edges = Set{NTuple{2,Int32}}()
    facets = Set{Tuple{Vararg{Int32}}}()
    # Facet incidence per (family, owner) group across all blocks — blocks of
    # the same type and owner belong to the same entity group, so faces
    # between them are interior while faces between two owners stay fixed.
    # The sorted key is only the identity: facet edges come from the
    # connectivity's own cycle order.
    inc = Dict{Tuple{Int32,Int32,Tuple{Vararg{Int32}}},Int}()
    cycleof = Dict{Tuple{Vararg{Int32}},Tuple{Vararg{Int32}}}()
    for (msh, _, _, cells, owners) in _cache_catalog(mesh, class)
        fam = HighOrder._p2opt_family(msh)
        fam === nothing && continue
        @inbounds for t in axes(cells, 2)
            owner = t <= length(owners) ? owners[t] : Int32(0)
            for facet in fam.facets
                key = HighOrder._p2opt_facet_key(cells, t, facet)
                group = (msh, owner, key)
                inc[group] = get(inc, group, 0) + 1
                haskey(cycleof, key) ||
                    (cycleof[key] = Tuple(cells[facet[i], t]
                                          for i in eachindex(facet)))
            end
        end
    end
    for ((_, _, f), c) in inc
        c == 1 || continue
        push!(facets, f)
        cycle = cycleof[f]
        n = length(cycle)
        for i in 1:n
            a, b = cycle[i], cycle[i==n ? 1 : i+1]
            push!(edges, minmax(a, b))
        end
    end
    return edges, facets
end

# Relocate interior high-order nodes of an unpublished MixedMesh in place —
# every certified family present in the cache — and slide curve-owned
# parametric mids along their owning curve.  `elastic`/`banded` select the
# two stages explicitly: `Mesh.HighOrderOptimize` derives them from its level
# (-1/1/2) while `API.mesh.optimize("HighOrderElastic")`/(`"HighOrder"`)
# request exactly one.  Returns `(stats, slid)` where `slid` maps node index
# to the accepted parameter vector, or nothing when no stage is enabled or
# the mesh carries no certified P2 block.
function _high_order_optimize_mixed_stages!(mesh::MixedMesh, class, model,
        caller::AbstractString; elastic::Bool, banded::Bool)
    (elastic || banded) || return nothing
    catalog = _cache_catalog(mesh, class)
    blocks = Tuple{Int32,Matrix{Int32}}[
        (msh, cells) for (msh, _, _, cells, _) in catalog]
    any(blk -> HighOrder._p2opt_family(blk[1]) !== nothing, blocks) ||
        return nothing
    fixed_edges, fixed_facets = _high_order_optimize_fixed_supports(mesh, class)
    sliders = _high_order_optimize_sliders(mesh, class, model, caller)
    movable = _high_order_optimize_movable(mesh, class)
    itermax = _high_order_optimize_int_option("Mesh.HighOrderIterMax", 100)
    passes = _high_order_optimize_int_option("Mesh.HighOrderPassMax", 25)
    # Stage order follows upstream: elastic smoothing first (mode -1 and the
    # first half of mode 2, no threshold band), then the banded optimization
    # stage (modes 1 and 2 honoring Mesh.HighOrderThresholdMin/Max).
    moved = 0
    stats = (0, NaN, NaN, 0, 0)
    if elastic
        s = HighOrder.optimize_p2_interior!(mesh.coords, blocks, movable;
            itermax=itermax, passes=passes,
            fixed_edges=fixed_edges, fixed_facets=fixed_facets,
            sliders=sliders,
            threshold_min=-Inf, threshold_max=Inf)
        moved += s[1]
        stats = (moved, s[2], s[3], stats[4] + s[4], stats[5] + s[5])
        s[1] > 0 &&
            _native_diagnostic_emit_locked(:info,
                "High-order elastic smoothing moved $(s[1]) node(s) " *
                "(min scaled Jacobian $(round(s[2];sigdigits=4)) -> " *
                "$(round(s[3];sigdigits=4)))")
    end
    if banded
        s = HighOrder.optimize_p2_interior!(mesh.coords, blocks, movable;
            itermax=itermax, passes=passes,
            fixed_edges=fixed_edges, fixed_facets=fixed_facets,
            sliders=sliders,
            threshold_min=_high_order_optimize_num_option(
                "Mesh.HighOrderThresholdMin", 0.1),
            threshold_max=_high_order_optimize_num_option(
                "Mesh.HighOrderThresholdMax", 2.0))
        base = isnan(stats[2]) ? s[2] : stats[2]
        moved += s[1]
        stats = (moved, base, s[3], stats[4] + s[4], stats[5] + s[5])
        s[1] > 0 &&
            _native_diagnostic_emit_locked(:info,
                "High-order optimize moved $(s[1]) interior node(s) " *
                "(min scaled Jacobian $(round(s[2];sigdigits=4)) -> " *
                "$(round(s[3];sigdigits=4)))")
    end
    # Accepted slider parameter states — the engine replaced `support.u` only
    # where a move was committed, so a bitwise compare identifies them.
    slid = Dict{Int,Vector{Float64}}()
    if sliders !== nothing && stats[1] > 0 && class !== nothing
        for (v, support) in sliders
            isequal(support.u, class.node_parametric[Int(v)]) ||
                (slid[Int(v)] = support.u)
        end
    end
    return stats, slid
end

# Option-gated entry used by the generation and set_order commit sites.
function _high_order_optimize_mixed!(mesh::MixedMesh, class, model,
                                     caller::AbstractString)
    _high_order_optimize_relocation() || return nothing
    return _high_order_optimize_mixed_stages!(mesh, class, model, caller;
        elastic=_high_order_optimize_elastic(),
        banded=_high_order_optimize_optimize())
end

# Same pass for a legacy-cache P2 overlay. `mids[i]` is the owning entity of
# mid node `nnodes_linear + i`; only same-dimension-owned mids can move
# (volume mids for a Tet10 overlay, surface mids for a Tri6 overlay — the
# engine's own wire/interior rules still fix the actual boundary supports).
# The overlay has a single owner group, so the combined complex boundary
# already covers every geometric support — no extra fixed edges are needed.
function _high_order_optimize_overlay_stages!(p2::Union{P2Mesh,P2TriMesh},
        mids::AbstractVector{<:Tuple{Int,Int32}}; elastic::Bool, banded::Bool)
    (elastic || banded) || return nothing
    dim = p2 isa P2Mesh ? 3 : 2
    cells = p2 isa P2Mesh ? p2.tet10 : p2.tri6
    size(cells, 2) == 0 && return nothing
    nlin = nnodes(p2) - length(mids)
    movable = falses(nnodes(p2))
    for i in eachindex(mids)
        mids[i][1] == dim && (movable[nlin + i] = true)
    end
    msh = p2 isa P2Mesh ? Int32(11) : Int32(9)
    itermax = _high_order_optimize_int_option("Mesh.HighOrderIterMax", 100)
    passes = _high_order_optimize_int_option("Mesh.HighOrderPassMax", 25)
    moved = 0
    stats = (0, NaN, NaN, 0, 0)
    if elastic
        s = HighOrder.optimize_p2_interior!(p2.coords, ((msh, cells),),
            movable; itermax=itermax, passes=passes,
            threshold_min=-Inf, threshold_max=Inf)
        moved += s[1]
        stats = (moved, s[2], s[3], stats[4] + s[4], stats[5] + s[5])
        s[1] > 0 &&
            _native_diagnostic_emit_locked(:info,
                "High-order elastic smoothing moved $(s[1]) node(s) " *
                "(min scaled Jacobian $(round(s[2];sigdigits=4)) -> " *
                "$(round(s[3];sigdigits=4)))")
    end
    if banded
        s = HighOrder.optimize_p2_interior!(p2.coords, ((msh, cells),),
            movable; itermax=itermax, passes=passes,
            threshold_min=_high_order_optimize_num_option(
                "Mesh.HighOrderThresholdMin", 0.1),
            threshold_max=_high_order_optimize_num_option(
                "Mesh.HighOrderThresholdMax", 2.0))
        base = isnan(stats[2]) ? s[2] : stats[2]
        moved += s[1]
        stats = (moved, base, s[3], stats[4] + s[4], stats[5] + s[5])
        s[1] > 0 &&
            _native_diagnostic_emit_locked(:info,
                "High-order optimize moved $(s[1]) interior node(s) " *
                "(min scaled Jacobian $(round(s[2];sigdigits=4)) -> " *
                "$(round(s[3];sigdigits=4)))")
    end
    return stats
end

function _high_order_optimize_overlay!(p2::Union{P2Mesh,P2TriMesh},
                                       mids::AbstractVector{<:Tuple{Int,Int32}})
    _high_order_optimize_relocation() || return nothing
    return _high_order_optimize_overlay_stages!(p2, mids;
        elastic=_high_order_optimize_elastic(),
        banded=_high_order_optimize_optimize())
end

# Retained overlays are live published state: relocate on a deep copy so a
# failure cannot mutate the mesh the session still exposes.
function _high_order_optimize_overlay(p2::Union{P2Mesh,P2TriMesh},
                                      mids::AbstractVector{<:Tuple{Int,Int32}};
                                      fresh::Bool)
    if fresh
        _high_order_optimize_overlay!(p2, mids)
        return p2
    end
    _high_order_optimize_relocation() || return p2
    copy_p2 = p2 isa P2Mesh ?
        P2Mesh(copy(p2.coords), copy(p2.tet10); tet_tag=copy(p2.tet_tag),
               require_positive_tets=false) :
        P2TriMesh(copy(p2.coords), copy(p2.tri6); tri_tag=copy(p2.tri_tag),
                  _check_geometry=false)
    _high_order_optimize_overlay!(copy_p2, mids)
    return copy_p2
end

# Plan-level entry used by the generation and set_order commit sites. Runs the
# pass on an unpublished `plan.mesh`/`plan.preserved_overlay`. A plan may carry
# the *live* cache object unchanged (e.g. `set_order(2)` on an already
# quadratic mesh, or a retained overlay): those are detached into copies first
# so the relocation can never mutate published session state, and the plan is
# rebound so the commit stays atomic.
function _high_order_optimize_plan!(plan, caller::AbstractString)
    _high_order_optimize_relocation() || return plan
    if get(plan, :mesh, nothing) isa MixedMesh
        mesh = plan.mesh
        if mesh === LAST_MESH[]
            # `set_order(2)` on an already-quadratic cache carries the live
            # mesh object — detach into a copy and rebind the classification
            # so the published `class.mesh` still matches `plan.mesh`
            # (`_cached_classification_locked` rejects any other pairing).
            mesh = MixedMesh(mesh.coords, mesh.blocks;
                physical_names=mesh.physical_names,
                entity_data=mesh.entity_data,
                elementary_entities=mesh.elementary_entities,
                periodic_links=mesh.periodic_links,
                ancillary_sections=mesh.ancillary_sections,
                data_sections=mesh.data_sections,
                partition_data=mesh.partition_data)
            pclass = get(plan, :class, nothing)
            class = pclass === nothing ? pclass :
                _mixed_rebind_class(pclass, mesh)
            if pclass !== nothing && class !== nothing
                # The rebind rebuilds classification vectors only; mid-support
                # parametric tables are keyed by primary edge and survive the
                # identical-index copy — carry them over explicitly.
                merge!(class.mid_parametric, pclass.mid_parametric)
                merge!(class.mid_keys, pclass.mid_keys)
            end
            plan = merge(plan, (mesh=mesh, class=class))
        end
        result = _high_order_optimize_mixed!(mesh, get(plan, :class, nothing),
            get(plan, :model, nothing), caller)
        if result !== nothing
            slid = result[2]
            if !isempty(slid)
                # Commit the accepted curve parameters atomically with the
                # relocated coordinates: entity_data and class each get a
                # plan-local copy carrying the new values.
                class = get(plan, :class, nothing)
                mesh, class = _high_order_optimize_commit_sliders(
                    mesh, class, slid)
                updates = NamedTuple()
                mesh === get(plan, :mesh, nothing) ||
                    (updates = (mesh=mesh,))
                class === get(plan, :class, nothing) ||
                    (updates = merge(updates, (class=class,)))
                isempty(keys(updates)) || (plan = merge(plan, updates))
            end
        end
    end
    overlay = get(plan, :preserved_overlay, nothing)
    if overlay isa Union{P2Mesh,P2TriMesh}
        mids = get(plan, :preserved_mid_owners, Tuple{Int,Int32}[])
        optimized = _high_order_optimize_overlay(
            overlay, mids; fresh = overlay !== LAST_MESH_HIGH_ORDER[])
        optimized === overlay ||
            (plan = merge(plan, (preserved_overlay=optimized,)))
    end
    return plan
end

# ══════════════════════════════════════════════════════════════════════════════
# API.mesh.optimize — quadratic-simplex arms (OPT01).
#
# Upstream `OptimizeMesh` (src/mesh/Generator.cpp) never checks element order
# for the region kernels: the default/`Relocate3D` native entries already
# ingest the actual order-2 cells.  For the 2-D methods it iterates every
# model face — `laplaceSmoothing` refuses order > 1 faces (its Msg::Error
# precedes the niter early-return), while `RelocateVertices` has no order
# gate and relocates each face-owned vertex through the face's (u,v)
# parameters (`relocateVertex` + `Maximize_Quality_Golden_Section`,
# re-lifted by `GFace::point`).  The bounded arms below accept the curved
# quadratic-simplex class for `Relocate2D`: face-owned mid nodes slide in
# the owning surface's parameters through the certificate-gated
# `optimize_p2_interior!` pass — the same (u,v) sliding that
# Mesh.HighOrderOptimize commits — while corner vertices and curve-owned
# supports keep their positions (the engine relocates certified-family
# interior slots only).  `"HighOrder"`/`"HighOrderElastic"` are the dedicated
# upstream method names for the same two stages this module stages for
# `Mesh.HighOrderOptimize`; both are routed here — `"HighOrderElastic"` to
# the unbounded elastic stage, `"HighOrder"` to the
# Mesh.HighOrderThresholdMin/Max-banded stage.  Every arm works on detached
# coordinate copies and commits atomically through the standard cache
# replacement path, so a failure leaves the published mesh untouched.
# ══════════════════════════════════════════════════════════════════════════════

# Does the mesh carry an order>1 surface element anywhere? — the exact gate
# upstream's laplaceSmoothing applies per face (`getPolynomialOrder() > 1`,
# a pure order test — curvature is irrelevant).  For a MixedMesh that is any
# dim-2 block of order > 1, or any dim-3 block of order > 1: upstream keeps
# the order-2 boundary faces of a quadratic volume mesh on its GFaces, so a
# Tet10 volume still trips the gate even though the cache stores no
# materialized surface block (RelocateVertices reaches those same faces and
# relocates their face-owned mids).  For a legacy simplex cache the retained
# order-2 overlay implies the actual surface (Tri6 surface cells, or the
# order-2 boundary faces of a Tet10 elevation).
_optimize_order2_surface(cached::MixedMesh) = begin
    for block in cached.blocks
        block isa ElementBlock || continue
        dim = msh_dimension(block.msh)
        dim >= 2 && Elements.msh_order(block.msh) > 1 && return true
    end
    return false
end
_optimize_order2_surface(cached::Mesh) = _high_order_overlay(cached) !== nothing

# `RelocateVertices(GFace, niter)` on an order-2 surface inside a MixedMesh:
# face-owned mid nodes slide in their owning surface's parameters.  `pairs`
# keeps the same entity-scoping the linear kernel applies through
# `_optimize_movable_mask` (upstream ignores dimTags entirely, emitting only
# a warning — the session already warns through the native methods).
function _optimize_relocate2d_mixed!(cached::MixedMesh, class, model, pairs,
                                     iterations, caller::AbstractString)
    classified = class !== nothing &&
        length(class.node_entities) == nnodes(cached) &&
        length(class.node_parametric) == nnodes(cached)
    classified || return nothing
    selected = Set{Tuple{Int,Int32}}(
        (dim, Int32(tag)) for (dim, tag) in pairs)
    in_scope(owner) = owner[1] == 2 &&
        (isempty(selected) || (2, owner[2]) in selected)
    nn = nnodes(cached)
    movable = falses(nn)
    for v in 1:nn
        movable[v] = in_scope(class.node_entities[v])
    end
    sliders_all = _high_order_optimize_sliders(cached, class, model, caller)
    sliders = sliders_all === nothing ? nothing :
        Dict{Int32,HighOrder._P2ParametricSupport}(
            v => support for (v, support) in sliders_all
            if in_scope(class.node_entities[Int(v)]))
    sliders === nothing && !any(movable) && return nothing
    coords = copy(cached.coords)
    blocks = Tuple{Int32,Matrix{Int32}}[
        (block.msh, block.nodes) for block in cached.blocks
        if block isa ElementBlock]
    stats = HighOrder.optimize_p2_interior!(coords, blocks, movable;
        itermax=_high_order_optimize_int_option("Mesh.HighOrderIterMax", 100),
        passes=iterations, sliders=sliders,
        threshold_min=-Inf, threshold_max=Inf)
    stats[1] == 0 && return nothing
    slid = Dict{Int,Vector{Float64}}()
    if sliders !== nothing
        for (v, support) in sliders
            stored = class.node_parametric[Int(v)]
            isequal(support.u, stored) || (slid[Int(v)] = support.u)
        end
    end
    mesh = _mixed_rebuild_metadata(cached, cached.blocks; coords=coords)
    new_class = _mixed_rebind_class(class, mesh)
    if !isempty(slid)
        mesh, new_class = _high_order_optimize_commit_sliders(
            mesh, new_class, slid)
    end
    _replace_mesh_cache_with_overlay_locked!(mesh, new_class, caller;
        preserve_visibility=true)
    return nothing
end

# Same arm for a legacy simplex cache's retained curved overlay.  The linear
# skeleton and its corner vertices stay untouched; each face-owned mid with a
# stored surface parameter gets a slider like upstream's parametric
# relocateVertex.  Accepted parameter states are written into a plan-local
# copy of `class.mid_parametric` so the publish step keeps actual midpoint
# geometry and its parametric support together.
function _optimize_relocate2d_overlay!(cached::Mesh, class, overlay, pairs,
                                       iterations, caller::AbstractString)
    model = _model_locked()
    base = nnodes(cached)
    mids = LAST_MESH_HIGH_ORDER_MIDS[]
    cells = _p2_cells(overlay)
    selected = Set{Tuple{Int,Int32}}(
        (dim, Int32(tag)) for (dim, tag) in pairs)
    sliders = Dict{Int32,HighOrder._P2ParametricSupport}()
    for column in axes(cells, 2), (slot, i, j) in _p2_edge_slots(overlay)
        v = cells[slot, column]
        haskey(sliders, v) && continue
        index = v - base
        (1 <= index <= length(mids)) || continue
        owner = mids[index]
        owner[1] == 2 || continue
        (isempty(selected) || (2, owner[2]) in selected) || continue
        key = minmax(cells[i, column], cells[j, column])
        stored = class === nothing ? nothing :
            get(class.mid_parametric, key, nothing)
        stored === nothing && continue
        support = _high_order_optimize_surface_support(
            model, Int(owner[2]), stored)
        support === nothing || (sliders[v] = support)
    end
    isempty(sliders) && return nothing
    coords = copy(overlay.coords)
    stats = HighOrder.optimize_p2_interior!(coords,
        ((_p2_etype(overlay), cells),), falses(nnodes(overlay));
        itermax=_high_order_optimize_int_option("Mesh.HighOrderIterMax", 100),
        passes=iterations, sliders=sliders,
        threshold_min=-Inf, threshold_max=Inf)
    stats[1] == 0 && return nothing
    updated_mids = class === nothing ? nothing : copy(class.mid_parametric)
    if updated_mids !== nothing
        for (v, support) in sliders
            key = get(class.mid_keys, v, nothing)
            key === nothing && continue
            isequal(support.u, get(updated_mids, key, nothing)) ||
                (updated_mids[key] = copy(support.u))
        end
    end
    new_class = class === nothing || isequal(updated_mids, class.mid_parametric) ?
        class :
        _MeshClassification(cached, class.entity, class.entities,
            class.node_entities, class.boundaries, class.seg_entities,
            class.tri_entities, class.tet_entities, class.cell_entities,
            class.public_tags, class.edge_entities, class.face_entities,
            class.quad_entities, class.node_parametric, updated_mids,
            class.mid_keys)
    edited = overlay isa P2Mesh ?
        P2Mesh(coords, overlay.tet10; tet_tag=overlay.tet_tag,
               require_positive_tets=false, _check_geometry=false) :
        P2TriMesh(coords, overlay.tri6; tri_tag=overlay.tri_tag,
                  _check_geometry=false)
    _replace_mesh_cache_with_overlay_locked!(cached, new_class, caller;
        prepared_overlay=(edited, copy(mids)), preserve_visibility=true)
    return nothing
end

# `API.mesh.optimize("HighOrder")` and `("HighOrderElastic")` — the dedicated
# upstream method names for the staged interior relocation.  `niter` is not
# consulted upstream (Mesh.HighOrderIterMax/PassMax drive the engine) and
# `force`/`dimTags` are ignored like the native entries; a mesh without any
# certified quadratic block draws upstream's "linear elements" warning.
function _optimize_high_order_locked!(method, pairs, caller::AbstractString)
    model = _model_locked()
    cached = _cached_mesh_locked(caller)
    class = _cached_classification_locked(cached)
    isempty(pairs) ||
        @warn "Entity selection in the mesh optimizer is ignored"
    elastic = method == "HighOrderElastic"
    _native_diagnostic_emit_locked(:info, elastic ?
        "Applying elastic analogy to high-order mesh..." :
        "Optimizing high-order mesh...")
    if cached isa MixedMesh
        work = MixedMesh(copy(cached.coords), cached.blocks;
            physical_names=cached.physical_names,
            entity_data=cached.entity_data,
            elementary_entities=cached.elementary_entities,
            periodic_links=cached.periodic_links,
            ancillary_sections=cached.ancillary_sections,
            data_sections=cached.data_sections,
            partition_data=cached.partition_data)
        result = _high_order_optimize_mixed_stages!(work, class, model,
            caller; elastic=elastic, banded=!elastic)
        result === nothing && return _optimize_high_order_linear_notice()
        stats, slid = result
        stats[1] == 0 && return nothing
        mesh = work
        new_class = class === nothing ? nothing :
            _mixed_rebind_class(class, mesh)
        if !isempty(slid)
            mesh, new_class = _high_order_optimize_commit_sliders(
                mesh, new_class, slid)
        end
        _replace_mesh_cache_with_overlay_locked!(mesh, new_class, caller;
            preserve_visibility=true)
        return nothing
    end
    overlay = _high_order_overlay(cached)
    if overlay !== nothing
        edited = overlay isa P2Mesh ?
            P2Mesh(copy(overlay.coords), copy(overlay.tet10);
                tet_tag=copy(overlay.tet_tag), require_positive_tets=false,
                _check_geometry=false) :
            P2TriMesh(copy(overlay.coords), copy(overlay.tri6);
                tri_tag=copy(overlay.tri_tag), _check_geometry=false)
        stats = _high_order_optimize_overlay_stages!(edited,
            LAST_MESH_HIGH_ORDER_MIDS[]; elastic=elastic, banded=!elastic)
        stats === nothing && return _optimize_high_order_linear_notice()
        stats[1] == 0 && return nothing
        _replace_mesh_cache_with_overlay_locked!(cached, class, caller;
            prepared_overlay=(edited, copy(LAST_MESH_HIGH_ORDER_MIDS[])),
            preserve_visibility=true)
        return nothing
    end
    return _optimize_high_order_linear_notice()
end

function _optimize_high_order_linear_notice()
    # `HighOrderMeshOptimizer`/`HighOrderMeshElasticAnalogy` report this exact
    # warning and return the mesh untouched.
    _native_diagnostic_emit_locked(:warning,
        "Applying high-order mesh optimizer to mesh with linear elements")
    return nothing
end
