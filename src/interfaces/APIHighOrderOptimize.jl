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
# keep their nodes fixed.  This file only translates option state,
# classification ownership, and geometric-support structure into the
# engine's movable mask and fixed edge/facet sets.  It mutates the
# *unpublished* coordinate matrix of a freshly built plan/overlay, so a
# failure still leaves the previously published cache untouched.
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

# Candidate mask: a node may move only when it is owned by a surface or
# volume entity and carries no stored parametric coordinates (curve, point
# and parametric-surface supports keep their projected parameters — sliding
# those mids along their actual geometry is the separate parametric-relocation
# gap).  When the classification does not cover the mesh, ownership cannot be
# proven and the mask stays permissive — the engine's own mid-slot,
# single-family and boundary-support rules still keep corners and boundary
# supports fixed.
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
# every certified family present in the cache.  Returns the engine
# statistics tuple, or nothing when the pass is disabled / the mesh carries
# no certified P2 block.
function _high_order_optimize_mixed!(mesh::MixedMesh, class, caller::AbstractString)
    _high_order_optimize_relocation() || return nothing
    catalog = _cache_catalog(mesh, class)
    blocks = Tuple{Int32,Matrix{Int32}}[
        (msh, cells) for (msh, _, _, cells, _) in catalog]
    any(blk -> HighOrder._p2opt_family(blk[1]) !== nothing, blocks) ||
        return nothing
    fixed_edges, fixed_facets = _high_order_optimize_fixed_supports(mesh, class)
    stats = HighOrder.optimize_p2_interior!(mesh.coords, blocks,
        _high_order_optimize_movable(mesh, class);
        itermax=_high_order_optimize_int_option("Mesh.HighOrderIterMax", 100),
        passes=_high_order_optimize_int_option("Mesh.HighOrderPassMax", 25),
        fixed_edges=fixed_edges, fixed_facets=fixed_facets,
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
        _high_order_optimize_mixed!(mesh, get(plan, :class, nothing), caller)
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
