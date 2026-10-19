# ══════════════════════════════════════════════════════════════════════════════
# Mesh.HighOrderOptimize — bounded interior high-order-node relocation
# (OPT01 slice).
#
# Upstream semantics (gmsh 4.15.2): 0 disables high-order optimization,
# -1 applies elastic smoothing only, 1 optimizes (relocates interior
# high-order nodes to improve the minimum scaled Jacobian while keeping
# boundary/geometry supports fixed), 2 applies both.  This slice implements
# the certificate-gated relocation pass (the optimization stage of modes 1
# and 2); the dedicated elastic-smoothing stage (-1 and the second half of
# 2) remains an OPT01 gap — mode -1 currently performs no relocation.
#
# The pass itself lives in HighOrder.optimize_p2_interior! and covers every
# full-P2 family with an exact Bernstein positivity bound — Tri6, Quad9,
# Tet10, Hex27, Prism18 — while Pyr14 and every uncertified block family
# keep their nodes fixed.  Curve-owned mids additionally slide along their
# owning CAD curve in parameter space; face-owned (u,v) sliding remains a
# gap.  This file translates option state, classification ownership, and
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

# The relocation pass runs for upstream modes 1 (optimize) and 2 (both);
# mode -1 (elastic smoothing only) has no relocation stage.
_high_order_optimize_relocation() = _high_order_optimize_level() >= 1

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

# Curve-owned parametric sliders for `optimize_p2_interior!`: every node with
# owner dimension 1 and a stored finite `[u]` parameter gets a
# `_P2ParametricSupport` whose lift is the same
# `Model._model_high_order_curve_geometry`/`_model_high_order_curve_point`
# evaluator that produced the node's coordinates, clamped to the curve's
# native parameter bounds.  Face-owned `(dim 2, [u,v])` sliding is a separate
# OPT01 stage and not wired yet — those nodes stay fixed.  Returns
# `Dict{Int32,_P2ParametricSupport}` or nothing.
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
        owner[1] == 1 || continue
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
            stored,                          # shared: replaced, never mutated
            let geo = curve
                u -> begin
                    p = Model._model_high_order_curve_point(geo, u[1], caller)
                    (p[1], p[2], p[3])
                end
            end,
            Float64[curve.lower], Float64[curve.upper])
        sliders[Int32(v)] = support
    end
    return isempty(sliders) ? nothing : sliders
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
# parametric mids along their owning curve.  Returns `(stats, slid)` where
# `slid` maps node index to the accepted parameter vector, or nothing when
# the pass is disabled / the mesh carries no certified P2 block.
function _high_order_optimize_mixed!(mesh::MixedMesh, class, model,
                                     caller::AbstractString)
    _high_order_optimize_relocation() || return nothing
    catalog = _cache_catalog(mesh, class)
    blocks = Tuple{Int32,Matrix{Int32}}[
        (msh, cells) for (msh, _, _, cells, _) in catalog]
    any(blk -> HighOrder._p2opt_family(blk[1]) !== nothing, blocks) ||
        return nothing
    fixed_edges, fixed_facets = _high_order_optimize_fixed_supports(mesh, class)
    sliders = _high_order_optimize_sliders(mesh, class, model, caller)
    stats = HighOrder.optimize_p2_interior!(mesh.coords, blocks,
        _high_order_optimize_movable(mesh, class);
        itermax=_high_order_optimize_int_option("Mesh.HighOrderIterMax", 100),
        passes=_high_order_optimize_int_option("Mesh.HighOrderPassMax", 25),
        fixed_edges=fixed_edges, fixed_facets=fixed_facets,
        sliders=sliders,
        threshold_min=_high_order_optimize_num_option(
            "Mesh.HighOrderThresholdMin", 0.1),
        threshold_max=_high_order_optimize_num_option(
            "Mesh.HighOrderThresholdMax", 2.0))
    stats[1] > 0 &&
        _native_diagnostic_emit_locked(:info,
            "High-order optimize moved $(stats[1]) interior node(s) " *
            "(min scaled Jacobian $(round(stats[2];sigdigits=4)) -> " *
            "$(round(stats[3];sigdigits=4)))")
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

# Same pass for a legacy-cache P2 overlay. `mids[i]` is the owning entity of
# mid node `nnodes_linear + i`; only same-dimension-owned mids can move
# (volume mids for a Tet10 overlay, surface mids for a Tri6 overlay — the
# engine's own wire/interior rules still fix the actual boundary supports).
# The overlay has a single owner group, so the combined complex boundary
# already covers every geometric support — no extra fixed edges are needed.
function _high_order_optimize_overlay!(p2::Union{P2Mesh,P2TriMesh},
                                       mids::AbstractVector{<:Tuple{Int,Int32}})
    _high_order_optimize_relocation() || return nothing
    dim = p2 isa P2Mesh ? 3 : 2
    cells = p2 isa P2Mesh ? p2.tet10 : p2.tri6
    size(cells, 2) == 0 && return nothing
    nlin = nnodes(p2) - length(mids)
    movable = falses(nnodes(p2))
    for i in eachindex(mids)
        mids[i][1] == dim && (movable[nlin + i] = true)
    end
    msh = p2 isa P2Mesh ? Int32(11) : Int32(9)
    stats = HighOrder.optimize_p2_interior!(p2.coords, ((msh, cells),), movable;
        itermax=_high_order_optimize_int_option("Mesh.HighOrderIterMax", 100),
        passes=_high_order_optimize_int_option("Mesh.HighOrderPassMax", 25),
        threshold_min=_high_order_optimize_num_option(
            "Mesh.HighOrderThresholdMin", 0.1),
        threshold_max=_high_order_optimize_num_option(
            "Mesh.HighOrderThresholdMax", 2.0))
    stats[1] > 0 &&
        _native_diagnostic_emit_locked(:info,
            "High-order optimize moved $(stats[1]) interior node(s) " *
            "(min scaled Jacobian $(round(stats[2];sigdigits=4)) -> " *
            "$(round(stats[3];sigdigits=4)))")
    return stats
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
