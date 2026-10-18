# ══════════════════════════════════════════════════════════════════════════════
# Mesh.HighOrderOptimize — bounded interior mid-node relocation (OPT01 slice).
#
# Upstream semantics (gmsh 4.15.2): 0 disables high-order optimization,
# -1 applies elastic smoothing only, 1 optimizes (relocates interior
# high-order nodes to improve the minimum scaled Jacobian while keeping
# boundary/geometry supports fixed), 2 applies both.  This slice implements
# the certificate-gated relocation pass (the optimization stage of modes 1
# and 2); the dedicated elastic-smoothing stage (-1 and the second half of
# 2) remains an OPT01 gap — mode -1 currently performs no relocation.
#
# The pass itself lives in HighOrder.optimize_p2_tet_interior!; this file only
# translates option state, classification ownership, and geometric-support
# structure into the engine's movable mask and fixed-edge set.  It mutates the
# *unpublished* coordinate matrix of a freshly built plan/overlay, so a failure
# still leaves the previously published cache untouched.
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

# Candidate mask: a node may move only when it is owned by a volume entity and
# carries no stored parametric coordinates (boundary and interface supports
# keep their projected parameters).  When the classification does not cover
# the mesh, ownership cannot be proven and the mask stays permissive — the
# engine's own mid-slot/boundary-edge rules still keep corners and boundary
# supports fixed.
function _high_order_optimize_movable(mesh::MixedMesh, class)
    nn = size(mesh.coords, 2)
    movable = falses(nn)
    classified = class !== nothing &&
        length(class.node_entities) == nn &&
        length(class.node_parametric) == nn
    for v in 1:nn
        movable[v] = !classified ||
            (class.node_entities[v][1] == 3 && class.node_parametric[v] === nothing)
    end
    return movable
end

# Primary edges whose mid nodes are geometric supports. `boundary_faces` is
# applied per element-owner group, so faces shared between two volume owners
# (internal interfaces of a conformal multi-region mesh) count as each group's
# boundary even though they are interior to the combined tet complex.
function _high_order_optimize_fixed_edges(mesh::MixedMesh, class)
    groups = Dict{Int32,Vector{NTuple{4,Int32}}}()
    for (msh, dim, _, cells, owners) in _cache_catalog(mesh, class)
        (dim == 3 && msh == Int32(11)) || continue
        for t in axes(cells, 2)
            owner = t <= length(owners) ? owners[t] : Int32(0)
            push!(get!(() -> NTuple{4,Int32}[], groups, owner),
                  (cells[1,t], cells[2,t], cells[3,t], cells[4,t]))
        end
    end
    edges = Set{NTuple{2,Int32}}()
    for tets in values(groups)
        prim = Matrix{Int32}(undef, 4, length(tets))
        @inbounds for (t, tet) in pairs(tets), k in 1:4
            prim[k, t] = tet[k]
        end
        bnd, _ = boundary_faces(prim)
        for f in bnd
            a, b, c = f[1], f[2], f[3]
            push!(edges, (min(a,b),max(a,b)))
            push!(edges, (min(b,c),max(b,c)))
            push!(edges, (min(a,c),max(a,c)))
        end
    end
    return edges
end

# Relocate interior Tet10 mid nodes of an unpublished MixedMesh in place.
# Returns the engine statistics tuple, or nothing when the pass is disabled /
# the mesh carries no quadratic tetrahedra.
function _high_order_optimize_mixed!(mesh::MixedMesh, class, caller::AbstractString)
    _high_order_optimize_relocation() || return nothing
    any(blk -> blk.msh == 11, mesh.blocks) || return nothing
    blocks = [blk.nodes for blk in mesh.blocks if blk.msh == 11]
    stats = HighOrder.optimize_p2_tet_interior!(mesh.coords, blocks,
        _high_order_optimize_movable(mesh, class);
        itermax=_high_order_optimize_int_option("Mesh.HighOrderIterMax", 100),
        passes=_high_order_optimize_int_option("Mesh.HighOrderPassMax", 25),
        fixed_edges=_high_order_optimize_fixed_edges(mesh, class),
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

# Same pass for a legacy-cache P2Mesh overlay. `mids[i]` is the owning entity
# of mid node `nnodes_linear + i`; only volume-owned mids can move. The mesh
# has a single owner group, so the combined tet boundary already covers every
# geometric support — no extra fixed edges are needed.
function _high_order_optimize_overlay!(p2::P2Mesh,
                                       mids::AbstractVector{<:Tuple{Int,Int32}})
    _high_order_optimize_relocation() || return nothing
    ntets(p2) == 0 && return nothing
    nlin = nnodes(p2) - length(mids)
    movable = falses(nnodes(p2))
    for i in eachindex(mids)
        mids[i][1] == 3 && (movable[nlin + i] = true)
    end
    stats = HighOrder.optimize_p2_tet_interior!(p2.coords, p2.tet10, movable;
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
function _high_order_optimize_overlay(p2::P2Mesh,
                                      mids::AbstractVector{<:Tuple{Int,Int32}};
                                      fresh::Bool)
    if fresh
        _high_order_optimize_overlay!(p2, mids)
        return p2
    end
    _high_order_optimize_relocation() || return p2
    copy_p2 = P2Mesh(copy(p2.coords), copy(p2.tet10); tet_tag=copy(p2.tet_tag),
                     require_positive_tets=false)
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
    if overlay isa P2Mesh
        mids = get(plan, :preserved_mid_owners, Tuple{Int,Int32}[])
        optimized = _high_order_optimize_overlay(
            overlay, mids; fresh = overlay !== LAST_MESH_HIGH_ORDER[])
        optimized === overlay ||
            (plan = merge(plan, (preserved_overlay=optimized,)))
    end
    return plan
end
