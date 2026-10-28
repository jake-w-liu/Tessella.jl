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

# ══════════════════════════════════════════════════════════════════════════════
# Linear-MixedMesh surface kernels (OPT01): upstream `laplaceSmoothing(GFace)`
# / `RelocateVertices(GFace)` — per-face parametric relocation.
#
# For each model face's classified linear `:tri`/`:qua` cells, face-owned
# vertices (owner dim 2, matching tag) are moved in the owning surface's
# (u,v) parameter space and re-lifted through `model_value` — the same
# evaluator family that produced their stored parameters
# (reparamMeshVertexOnFace maps every incident vertex onto the face; stored
# node_parametric covers face-owned vertices and `model_parametrization`
# projects edge/point-owned neighbours).  `"Laplace2D"` is upstream's
# `_relocate`: a metric-weighted centroid of edge neighbours (first
# fundamental form from `model_derivative`, mirroring `buildMetric`), then
# up to five damped trials accepted only when `_isItAGoodIdeaToMoveThatVertex`
# passes — parametric-space cell area may not grow and the minimum
# etaShapeMeasure may not drop.  `"Relocate2D"` is `relocateVertex`: the
# (u,v) centroid of every incident cell vertex gives the target, then
# `Maximize_Quality_Golden_Section` on [0,1] (tol 1e-6, the GFace default)
# against the worst incident shape measure — eta for 4-node cells,
# |gamma| otherwise.  A vertex without a stored parameter takes upstream's
# `relocateVertex2` real-space arm: golden search along the segment to the
# incident-vertex centroid with positions re-projected by
# `model_closest_point`.
#
# Bounded-scope deviations from upstream, kept truthful: Tessella does not
# track per-face boundary-layer column marks, so
# `getAllBoundaryLayerVertices` is an empty set and every eligible
# face-owned vertex participates; and periodic-seam vertices are not
# detected, so the `isSeam` bail-out is omitted.  Both affect only which
# vertices move, never the move semantics themselves.
# ══════════════════════════════════════════════════════════════════════════════

# (u,v) of `node` on face `tag`: stored parameters when face-owned, else a
# `model_parametrization` projection — `reparamMeshVertexOnFace` upstream.
function _surface_optimize_uv(model,class,tag,node,coords,params)
    stored=params[node]
    owner=class.node_entities[node]
    if owner[1]==2 && Int(owner[2])==tag && stored!==nothing &&
            length(stored)>=2 && isfinite(stored[1]) && isfinite(stored[2])
        return Float64[stored[1],stored[2]]
    end
    point=Float64[coords[1,node],coords[2,node],coords[3,node]]
    value=try
        Model.model_parametrization(model,2,tag,point)
    catch err
        (err isa InterruptException||err isa OutOfMemoryError)&&rethrow()
        return nothing
    end
    (value===nothing||length(value)<2||!isfinite(value[1])||!isfinite(value[2]))&&
        return nothing
    return Float64[value[1],value[2]]
end

# `gf->point(uv)`: lift through `model_value`; nothing on evaluator failure.
function _surface_optimize_lift(model,tag,uv)
    values=try
        Model.model_value(model,2,tag,uv;
            old_ruled_surface=!iszero(OPTIONS["Geometry.OldRuledSurface"]))
    catch err
        (err isa InterruptException||err isa OutOfMemoryError)&&rethrow()
        return nothing
    end
    (length(values)==3&&all(isfinite,values))||return nothing
    return (values[1],values[2],values[3])
end

# `buildMetric(gf,uv)`: first fundamental form E,F,G from `model_derivative`.
function _surface_optimize_metric(model,tag,uv)
    values=try
        Model.model_derivative(model,2,tag,uv;
            old_ruled_surface=!iszero(OPTIONS["Geometry.OldRuledSurface"]))
    catch err
        (err isa InterruptException||err isa OutOfMemoryError)&&rethrow()
        return nothing
    end
    length(values)>=6 || return nothing
    e=values[1]^2+values[2]^2+values[3]^2
    f=values[1]*values[4]+values[2]*values[5]+values[3]*values[6]
    g=values[4]^2+values[5]^2+values[6]^2
    (isfinite(e)&&isfinite(f)&&isfinite(g))||return nothing
    return (e,f,g)
end

# MTriangle::etaShapeMeasure — 1-|60°-a_min|/60 over interior angles.
function _surface_optimize_tri_eta(a,b,c)
    angle(p,q,r)=begin
        u=_native_cell_difference(p,q);v=_native_cell_difference(r,q)
        180.0*atan(_native_cell_norm(_native_cell_cross(u,v)),
            _native_cell_dot(u,v))/pi
    end
    amin=min(angle(a,b,c),angle(b,c,a),angle(c,a,b))
    return 1.0-abs(60.0-amin)*(1.0/60.0)
end

# The relocateVertex objective: etaShapeMeasure for 4-node cells,
# |gammaShapeMeasure| otherwise.
_surface_optimize_cell_objective(family,points)=
    length(points)==4 ? MeshElementQuality._quadrangle_native_eta(points) :
        abs(_native_cell_triangle_gamma(points[1],points[2],points[3]))

# The _isItAGoodIdeaToMoveThatVertex measure: etaShapeMeasure on every cell.
_surface_optimize_cell_eta(family,points)=
    length(points)==4 ? MeshElementQuality._quadrangle_native_eta(points) :
        _surface_optimize_tri_eta(points[1],points[2],points[3])

@inline _surface_uv_tri_area(u1,u2,u3)=0.5*abs(
    (u2[1]-u1[1])*(u3[2]-u1[2])-(u3[1]-u1[1])*(u2[2]-u1[2]))

# surfaceFaceUV(t,gf,false): tri → half-cross; quad → min over the two
# diagonal split sums.
function _surface_uv_cell_area(uvs)
    length(uvs)==3 && return _surface_uv_tri_area(uvs[1],uvs[2],uvs[3])
    a1=_surface_uv_tri_area(uvs[1],uvs[2],uvs[3])+
        _surface_uv_tri_area(uvs[3],uvs[4],uvs[1])
    a2=_surface_uv_tri_area(uvs[2],uvs[3],uvs[4])+
        _surface_uv_tri_area(uvs[4],uvs[1],uvs[2])
    return min(a1,a2)
end

# `_isItAGoodIdeaToMoveThatVertex`: probe `trial`, then restore.  Accepts
# only when the lifted point exists, parametric cell area does not grow
# (1+1e-10 slack) and the minimum eta does not drop.  The probe vertex's
# (u,v) is `trial` during the new-state pass and its committed (u,v)
# during the old-state pass — the shared `uvmap` cache is never mutated.
function _surface_optimize_good_move(coords,uv_for,node,incident,cells,
                                     model,tag,trial)
    pafter=_surface_optimize_lift(model,tag,trial)
    pafter===nothing && return false
    points(cell)=ntuple(i->(coords[1,cell[2][i]],coords[2,cell[2][i]],
        coords[3,cell[2][i]]),length(cell[2]))
    function sweep(probe)
        surface=0.0;minq=1.0
        for index in incident
            cell=cells[index]
            cell_uvs=Union{Nothing,Vector{Float64}}[
                slot==node ? probe : uv_for(slot) for slot in cell[2]]
            any(u->u===nothing,cell_uvs) && continue
            surface+=_surface_uv_cell_area(cell_uvs)
            minq=min(minq,_surface_optimize_cell_eta(cell[1],points(cell)))
        end
        return surface,minq
    end
    surface_old,minq=sweep(uv_for(node))
    px,py,pz=coords[1,node],coords[2,node],coords[3,node]
    coords[1,node],coords[2,node],coords[3,node]=pafter
    surface_new,minq_new=sweep(trial)
    coords[1,node],coords[2,node],coords[3,node]=px,py,pz
    (1.0+1e-10)*surface_old<surface_new && return false
    minq_new<minq && return false
    return true
end

# Upstream `_relocate` (laplaceSmoothing): metric-weighted neighbour
# centroid in (u,v), five damped trials gated by _isItAGoodIdea…  Returns
# the last accepted trial (u,v), or `nothing` when no move was taken.
function _surface_laplace_vertex!(coords,uv_for,node,incident,cells,
                                  edges_of,model,tag)
    before=uv_for(node)
    before===nothing && return nothing
    neighbours=Dict{Int,Vector{Float64}}()
    for index in incident
        family,nodes=cells[index]
        for (i,j) in edges_of[family]
            # Element edge tables are 0-based vertex slots.
            other=nodes[i+1]==node ? nodes[j+1] :
                  nodes[j+1]==node ? nodes[i+1] : 0
            (other==0||haskey(neighbours,other)) && continue
            uv=uv_for(other)
            uv===nothing || (neighbours[other]=uv)
        end
    end
    isempty(neighbours) && return nothing
    after=(0.0,0.0);weight=0.0
    for (other,adj) in neighbours
        du=adj[1]-before[1];dv=adj[2]-before[2]
        norm=sqrt(du*du+dv*dv)
        norm>0 || continue
        du/=norm;dv/=norm
        metric=_surface_optimize_metric(model,tag,adj)
        metric===nothing && continue
        f=sqrt(metric[1]*du*du+2*metric[2]*du*dv+metric[3]*dv*dv)
        isfinite(f) || continue
        after=(after[1]+adj[1]*f,after[2]+adj[2]*f);weight+=f
    end
    weight>0 || return nothing
    after=(after[1]/weight,after[2]/weight)
    factor=1.0;accepted=nothing
    for _ in 1:5
        trial=Float64[after[1]*factor+before[1]*(1.0-factor),
                      after[2]*factor+before[2]*(1.0-factor)]
        if _surface_optimize_good_move(coords,uv_for,node,incident,cells,
                                       model,tag,trial)
            pafter=_surface_optimize_lift(model,tag,trial)
            if pafter!==nothing
                coords[1,node],coords[2,node],coords[3,node]=pafter
                accepted=trial
            end
        end
        factor/=1.4
    end
    return accepted
end

# Parametric golden section on ξ∈[0,1] — mirrors
# Maximize_Quality_Golden_Section including its `worst` bookkeeping: the
# caller records f(0) first, accepts the full move when it already beats
# 0.5, and falls back to it when the refined best does not improve on it.
function _surface_golden_search(objective,tol,worst)
    lambda=0.5*(sqrt(5.0)-1.0);mu=0.5*(3.0-sqrt(5.0))
    a=0.0;b=1.0
    x1=b-lambda*(b-a);x2=a+lambda*(b-a)
    fx1=objective(x1);fx2=objective(x2)
    while abs(b-a)>=tol
        if fx1<fx2
            a=x1;abs(b-a)<tol&&break
            x1=x2;fx1=fx2;x2=b-mu*(b-a);fx2=objective(x2)
        else
            b=x2;abs(b-a)<tol&&break
            x2=x1;fx2=fx1;x1=a+mu*(b-a);fx1=objective(x1)
        end
    end
    final=objective(a)
    final<worst && return 0.0
    return a
end

# Upstream `relocateVertex(gf,ver,lt)`: parametric centroid + golden
# search; `relocateVertex2` real-space arm when the vertex stores no (u,v).
# Returns the accepted (u,v), `true` for the real-space arm (upstream
# leaves the vertex's stored parameters untouched there), or `nothing`.
function _surface_relocate_vertex!(coords,uv_for,node,incident,cells,
                                   model,tag)
    points(cell)=ntuple(i->(coords[1,cell[2][i]],coords[2,cell[2][i]],
        coords[3,cell[2][i]]),length(cell[2]))
    objective=function(point_at,xi)
        xyz=point_at(xi)
        xyz===nothing && return -1.0
        px,py,pz=coords[1,node],coords[2,node],coords[3,node]
        coords[1,node],coords[2,node],coords[3,node]=xyz
        worst=1.0
        for index in incident
            cell=cells[index]
            worst=min(worst,_surface_optimize_cell_objective(
                cell[1],points(cell)))
        end
        coords[1,node],coords[2,node],coords[3,node]=px,py,pz
        return worst
    end
    p2=uv_for(node)
    if p2===nothing
        # relocateVertex2 — real-space centroid + closestPoint projection.
        counter=0;cx=0.0;cy=0.0;cz=0.0
        for index in incident, slot in cells[index][2]
            cx+=coords[1,slot];cy+=coords[2,slot];cz+=coords[3,slot]
            counter+=1
        end
        counter==0 && return nothing
        cx/=counter;cy/=counter;cz/=counter
        point_at=function(xi)
            q=Float64[cx*(1.0-xi)+coords[1,node]*xi,
                      cy*(1.0-xi)+coords[2,node]*xi,
                      cz*(1.0-xi)+coords[3,node]*xi]
            closest,_=try
                Model.model_closest_point(model,2,tag,q)
            catch err
                (err isa InterruptException||err isa OutOfMemoryError)&&rethrow()
                return nothing
            end
            (length(closest)>=3&&all(isfinite,view(closest,1:3)))||return nothing
            return (closest[1],closest[2],closest[3])
        end
        worst=objective(point_at,0.0)
        xi_final=worst>0.5 ? 0.0 :
            _surface_golden_search(xi->objective(point_at,xi),1e-6,worst)
        xyz=point_at(xi_final)
        xyz===nothing && return nothing
        coords[1,node],coords[2,node],coords[3,node]=xyz
        return true
    end
    counter=0;cu=0.0;cv=0.0
    for index in incident, slot in cells[index][2]
        uv=uv_for(slot)
        uv===nothing || (cu+=uv[1];cv+=uv[2];counter+=1)
    end
    counter==0 && return nothing
    cu/=counter;cv/=counter
    point_at=function(xi)
        uv=Float64[cu*(1.0-xi)+p2[1]*xi,cv*(1.0-xi)+p2[2]*xi]
        return _surface_optimize_lift(model,tag,uv)
    end
    worst=objective(point_at,0.0)
    xi_final=worst>0.5 ? 0.0 :
        _surface_golden_search(xi->objective(point_at,xi),1e-6,worst)
    uv=Float64[cu*(1.0-xi_final)+p2[1]*xi_final,
               cv*(1.0-xi_final)+p2[2]*xi_final]
    xyz=_surface_optimize_lift(model,tag,uv)
    xyz===nothing && return nothing
    coords[1,node],coords[2,node],coords[3,node]=xyz
    return uv
end

# Linear MixedMesh driver for Laplace2D/Relocate2D: per-face vertex→cell
# adjacency over the classified dim-2 linear cells, per-face niter passes
# over the face-owned vertices in deterministic tag order, then the same
# atomic coordinate/parameter commit as the quadratic arm.
function _optimize_surface_mixed_linear!(cached::MixedMesh,class,model,pairs,
                                         iterations,caller::AbstractString,
                                         relocate::Bool)
    class===nothing && return nothing
    nn=nnodes(cached)
    (length(class.node_entities)==nn &&
     length(class.node_parametric)==nn) || return nothing
    selected=Set{Tuple{Int,Int32}}(
        (dim,Int32(tag)) for (dim,tag) in pairs)
    faces=Dict{Int32,Vector{Tuple{Symbol,Vector{Int32}}}}()
    for block in cached.blocks
        block isa ElementBlock || continue
        spec=msh_spec(block.msh)
        (spec.dim==2&&spec.order==1&&spec.family in (:tri,:qua))||continue
        owners=get(class.cell_entities,block.msh,Int32[])
        length(owners)==size(block.nodes,2) || return nothing
        for column in axes(block.nodes,2)
            tag=owners[column]
            (isempty(selected)||(2,tag) in selected)||continue
            push!(get!(Vector{Tuple{Symbol,Vector{Int32}}},faces,tag),
                (spec.family,collect(block.nodes[:,column])))
        end
    end
    isempty(faces) && return nothing
    edges_of=Dict(:tri=>Elements._TRI_EDGES,:qua=>Elements._QUA_EDGES)
    coords=copy(cached.coords)
    params=Union{Nothing,Vector{Float64}}[
        p===nothing ? nothing : copy(p) for p in class.node_parametric]
    moved=false
    for tag in sort!(collect(keys(faces)))
        cells=faces[tag]
        adjacency=Dict{Int,Vector{Int}}()
        for (index,(family,nodes)) in enumerate(cells),slot in nodes
            push!(get!(Vector{Int},adjacency,slot),index)
        end
        uvmap=Dict{Int,Union{Nothing,Vector{Float64}}}()
        # Per-face (u,v) lookup: stored parameters for face-owned nodes,
        # `model_parametrization` projections for edge/point-owned
        # neighbours — `reparamMeshVertexOnFace` upstream.  Accepted moves
        # are written back so later vertices see the updated parametric
        # position just like upstream's setParameter.
        uv_for=function(node)
            get!(uvmap,node) do
                _surface_optimize_uv(model,class,tag,node,coords,params)
            end
        end
        in_scope(owner)=owner[1]==2&&Int(owner[2])==tag&&
            (isempty(selected)||(2,Int32(tag)) in selected)
        for _ in 1:iterations
            for vertex in sort!(collect(keys(adjacency)))
                in_scope(class.node_entities[vertex]) || continue
                incident=adjacency[vertex]
                result=relocate ?
                    _surface_relocate_vertex!(coords,uv_for,vertex,incident,
                        cells,model,tag) :
                    _surface_laplace_vertex!(coords,uv_for,vertex,incident,
                        cells,edges_of,model,tag)
                result===nothing && continue
                moved=true
                if result isa Vector{Float64}
                    uvmap[vertex]=copy(result)
                    params[vertex]=copy(result)
                end
            end
        end
    end
    moved || return nothing
    slid=Dict{Int,Vector{Float64}}()
    for vertex in 1:nn
        stored=params[vertex]
        original=class.node_parametric[vertex]
        (stored===nothing&&original===nothing) && continue
        (stored===nothing||original===nothing||
         !isequal(stored,original)) &&
            stored!==nothing && (slid[vertex]=stored)
    end
    mesh=_mixed_rebuild_metadata(cached,cached.blocks;coords=coords)
    new_class=_mixed_rebind_class(class,mesh)
    if !isempty(slid)
        mesh,new_class=_high_order_optimize_commit_sliders(
            mesh,new_class,slid)
    end
    _replace_mesh_cache_with_overlay_locked!(mesh,new_class,caller;
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
