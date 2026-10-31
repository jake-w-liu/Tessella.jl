# PAR01 — detached independent-Curve production under `General.NumThreads`.
#
# `_dim01_curve_phase!` still owns the serial `pending` traversal, the
# `_dim01_emit_curve!` publication merge, the keep gates, the deferred
# periodic-parameter copy and the stalled-cycle seed arm.  With
# `General.NumThreads >= 2` each pass additionally detects runs of
# independent singleton Curves and produces their `_model_mesh_curve!`
# receipts — the stored parameter vector, `:keep`, `:pending` or a
# `_ModelNativeCurveGrade` — on per-unit transaction copies
# (`Model._model_mesh_transaction_copy`) inside `Threads.@threads`, each task
# writing a disjoint result slot while the shared model and assembly stay
# untouched.  At merge time each slot is checked for staleness: an earlier
# merged Curve can invalidate a detached receipt only by changing
# `m.curve_params` for one of its footprint entities.  A stale or failed slot
# is recomputed by the identical live serial call; a fresh slot's receipt is
# fed through the same emit arm the serial call would have taken, so node and
# element UIDs are still minted by the live assembly in tag order and the
# produced mesh is bit-identical to the serial one.
#
# The eligibility contract is deliberately conservative, mirroring the
# Surface arm.  Only a Curve whose serial arm is the plain
# `_model_mesh_curve!` produce can detach: every keep-gate case (fully
# discrete, `Degenerated`, `MeshOnlyVisible`/`MeshOnlyEmpty`, retained
# `curve_params`), periodic slaves and masters (the deferred
# `_dim01_native_periodic_parameters` copy and `_model_periodic_cycle_seeds`
# cycle machinery), transfinite Curves, extrusion copies/sources and sweep
# gates (`_extrude_gate`/`extrude_sources`), embedded or embedding Curves,
# dim-1 compound members, Curves carrying a discrete/attached record and any
# stateful user size callback all keep the original serial arm — their
# producer contracts consume or publish shared state outside the recorded
# footprint.  A resolved background `options.inner_field` is *not* an
# exclusion: catalog fields are immutable read-only evaluators (the lone
# mutable `ExternalProcessField` serializes every query on its own lock), so
# the same object bound into `detached_options` is safe to grade against —
# the same contract the other arms already keep for `options.inner_field`.
#
# The RNG ledger is uniform across arms: each detached task runs a private
# `random_state` Ref (`_parallel_task_options`) entering at the batch's
# serial `Mesh.RandomSeed` position and reporting its exit value.
# `options.random_state` feeds only `ModelNativeInitialTet`'s perturbation
# kernel and the native hybrid batch/components paths — no 1-D grading
# reaches it — so detached Curves are expected to observe `seed_out==seed_in`
# every time; the ledger turns that contract into a check exactly like the
# Surface arm: an unexpected consumption is adopted only when the live serial
# position still equals the observed entry value, otherwise the slot is stale
# and re-runs serially.
#
# Non-contiguous scheduling mirrors the other arms: `_parallel_eager_indices`
# lets a fill also produce later eligible slots whose footprints are disjoint
# from the accumulated write-potential of every unit in between (a serial
# Curve between them cannot dirty their receipts).  Each eligible work
# records `gap_writes` — the write-set estimate of the serial slots since the
# previous eligible slot — and the merge-time staleness check remains the
# sole authority.
#
# Two determinism details the worker path must preserve, as in the Surface
# arm:
# * `GmshLibm._win_cli_precision()` is task-local (CLI x87 libm selector);
#   worker tasks re-run the produce under `_with_win_cli_precision` when the
#   flag is set on the orchestrating task, so trig/pow kernels cannot differ
#   between serial and parallel production.
# * Task exceptions are captured into the slot and the affected unit re-runs
#   serially, so a failing producer reports at the same merge position the
#   serial traversal would have reached.

import ..GmshLibm

struct _ParallelCurveWork
    tag::Int
    footprint::Set{Tuple{Int,Int32}}
    # Write-set estimate of the serial slots since the previous eligible
    # slot (`nothing` = an unbounded serial write precedes) — read only by
    # `_parallel_eager_indices`, not part of the staleness contract.
    gap_writes::Union{Nothing,Set{Tuple{Int,Int32}}}
end

# `General.NumThreads` requests the worker count; `Mesh.MaxNumThreads1D` caps
# the Curve stage when > 0 (upstream's per-dimension limit).
function _parallel_curve_workers()
    workers=_partition_workers()
    workers>=2 || return 1
    cap=get(OPTIONS,"Mesh.MaxNumThreads1D",0.0)
    cap>0 && (workers=min(workers,max(1,floor(Int,cap))))
    return workers
end

# The serial-read serial-write contract of one detached `_model_mesh_curve!`
# call: every entity whose stored discretization the produce can consult or
# rewrite.  An eligible Curve reads only its own `curve_params` slot (the
# `MeshOnlyEmpty` gate, dropped by `deMeshGEdge` before grading) plus
# immutable geometry — its `m.curves` endpoints and `m.points` coordinates,
# `curve_types`/`curve_geometry`, the `transfinite_curves`/`recombine`/
# `size_at_params`/`algorithm` attributes and `point_size` — none of which the
# merge mutates.  The endpoint Points still join the footprint so the dirty
# scan stays a faithful superset of the produce's read set.  A dirty mark on
# any member makes the detached receipt stale.
function _parallel_curve_footprint(m::GeoModel,tag::Int)
    footprint=Set{Tuple{Int,Int32}}()
    push!(footprint,(1,Int32(tag)))
    endpoints=get(m.curves,tag,nothing)
    if endpoints!==nothing
        push!(footprint,(0,Int32(endpoints[1])))
        push!(footprint,(0,Int32(endpoints[2])))
    end
    return footprint
end

# Per-Curve eligibility under the phase-level gates, evaluated for one pass's
# `ordered` pending list.  `nothing` means the whole phase keeps the serial
# traversal (zero added bookkeeping); an entry of `nothing` in the vector
# means that Curve runs serially in place.  The keep gates are folded into
# eligibility (their predicates are phase-stable — `have_cells` only gains a
# Curve's own owner when that same Curve merges) so a kept Curve never pays
# for a detached produce its serial arm would discard.
function _parallel_curve_prepare(m::GeoModel,ordered::Vector{Int},
        options::Model._ModelMesh1DOptions;
        have_cells::Set{Tuple{Int,Int32}},reuse_curve_params::Bool)
    _parallel_curve_workers()>=2 || return nothing
    # A stateful size callback is consulted from inside the grading kernel and
    # cannot cross a transaction copy safely — it keeps the serial arm.  The
    # resolved `inner_field` is different: catalog fields are immutable
    # read-only evaluators (the lone mutable `ExternalProcessField` serializes
    # every query on its own lock), and `_parallel_task_options` already binds
    # the same object into every detached produce — the exact contract the
    # Surface/Region arms keep for `options.inner_field`/`size_field`.
    (options.callback===nothing &&
        m.meshing.size_callback===nothing) || return nothing
    periodic_members=Set{Int}()
    for constraint in values(m.periodic)
        constraint.dim==1 || continue
        push!(periodic_members,Int(abs(constraint.slave_entity)))
        push!(periodic_members,Int(abs(constraint.master_entity)))
    end
    embedded_curves=Set{Int}()
    for (key,members) in m.embeds
        key[1] in (2,3) || continue
        for (dim,entity) in members
            dim==1 && push!(embedded_curves,entity)
        end
    end
    compound_members=Set{Int}()
    for pair in m.meshing.compounds,member in pair.second
        pair.first==1 && push!(compound_members,member)
    end
    # Curves serving as a dim-1 extrusion source: their emitted node chains
    # feed `_dim01_extruded_curve_points`/`_model_curve_extrude_params` copies
    # of other entities — shared state outside the footprint.
    extrude_source_curves=Set{Int}()
    for link in values(m.meshing.extrude_sources)
        link[1]==1 && push!(extrude_source_curves,abs(link[2]))
    end
    prepared=Vector{Any}(nothing,length(ordered))
    for (index,curve) in enumerate(ordered)
        owner=(1,Int32(curve))
        # The keep gates stay serial and are re-evaluated verbatim below.
        _dim01_fully_discrete(m,(1,curve)) && continue
        Model._model_curve_degenerate_kind(m,curve)!=-1 && continue
        (options.mesh_only_visible &&
            model_entity_visibility(m,1,curve)==0) && continue
        (options.mesh_only_empty && owner in have_cells) && continue
        (reuse_curve_params && haskey(m.curve_params,curve)) && continue
        haskey(m.curves,curve) || continue
        curve in periodic_members && continue
        haskey(m.meshing.transfinite_curves,curve) && continue
        Model._extrude_gate(m,1,curve)!==nothing && continue
        haskey(m.meshing.extrude_sources,(1,curve)) && continue
        curve in extrude_source_curves && continue
        _dim01_record(m,(1,curve))!==nothing && continue
        haskey(m.embeds,(1,curve)) && continue
        curve in embedded_curves && continue
        curve in compound_members && continue
        prepared[index]=_ParallelCurveWork(curve,
            _parallel_curve_footprint(m,curve),nothing)
    end
    any(!isnothing,prepared) || return nothing
    # Second pass: record each eligible slot's preceding serial write-set for
    # the eager scheduler.  Keep-gated and ineligible Curves alike take the
    # serial arm; their write-set estimate is their footprint (the same
    # contract the dirty scan relies on), so a `nothing` estimate — which
    # cannot arise from `_parallel_curve_footprint`'s get-only reads — would
    # stop the scan at the next eligible slot.
    serial_writes=Set{Tuple{Int,Int32}}()
    serial_unbounded=false
    for (index,curve) in enumerate(ordered)
        work=prepared[index]
        if work===nothing
            if !serial_unbounded
                estimate=try
                    _parallel_curve_footprint(m,curve)
                catch err
                    err isa InterruptException && rethrow()
                    nothing
                end
                estimate===nothing ? (serial_unbounded=true) :
                    union!(serial_writes,estimate)
            end
            continue
        end
        prepared[index]=_ParallelCurveWork(work.tag,work.footprint,
            serial_unbounded ? nothing : copy(serial_writes))
        empty!(serial_writes);serial_unbounded=false
    end
    return prepared
end

# `m.curve_params` slots an eligible `_model_mesh_curve!` may touch: its own
# tag only (the produce's `deMeshGEdge` drop and the caller's writeback).
# A write to any other slot — a dependency the footprint never claimed — can
# never be transported, so like the Surface/Region arms it flags the receipt
# `foreign` and the merge re-runs the identical serial call.
function _parallel_curve_foreign_params(working::GeoModel,m::GeoModel,tag::Int)
    for (curve,params) in working.curve_params
        curve==tag && continue
        get(m.curve_params,curve,nothing)!=params && return true
    end
    for curve in keys(m.curve_params)
        curve==tag && continue
        haskey(working.curve_params,curve) || return true
    end
    return false
end

# One detached Curve produce on a premade transaction copy.  The task's
# private `random_state` Ref enters at the batch seed and reports its exit
# value — the draw ledger the merge checks.  `foreign` is the transport
# certificate shared with the other arms: the produce must write only its
# own `curve_params` slot — a mutation anywhere else in the model stores
# (including `next_tag` entity births) could never reach the live model.
function _parallel_curve_attempt(working::GeoModel,work::_ParallelCurveWork,
        m::GeoModel,detached_options::Model._ModelMesh1DOptions,
        base_seed::UInt32,caller::AbstractString,precise::Bool)
    out=try
        body=()->begin
            result=Model._model_mesh_curve!(working,work.tag,
                detached_options,caller;_native_producer=true)
            foreign=_parallel_curve_foreign_params(working,m,work.tag) ||
                _parallel_model_foreign_writes(working,m)
            (result,base_seed,detached_options.random_state[],foreign)
        end
        (;ok=true,value=precise ?
            GmshLibm._with_win_cli_precision(body) : body())
    catch err
        err isa InterruptException && rethrow()
        (;ok=false,value=err)
    end
    return out
end

# Threaded detached production over one contiguous eligible run of a pass
# plus every predicted-independent later slot (`_parallel_eager_indices`).
# `range` shares the batch's serial merge position: `base_dirty` fixes the
# dirt the detached copies could already see, `base_seed` fixes the
# `Mesh.RandomSeed` ledger position each task's private Ref starts from, and
# every slot records the produce's outcome for the in-order merge.  A batch
# of one would pay the copy for no concurrency, so it falls back to the
# serial arm by leaving its slot empty.  The produce mints no public UIDs and
# writes only the working copy's own `curve_params` slot — the merge replays
# the identical writeback and emits on the live assembly, so no UID remap
# ledger is needed.
function _parallel_curve_fill!(results::Vector{Any},range,m::GeoModel,
        prepared,options::Model._ModelMesh1DOptions,
        dirty::Set{Tuple{Int,Int32}},precise::Bool,caller::AbstractString)
    base_dirty=copy(dirty)
    base_seed=options.random_state[]
    eager=_parallel_eager_indices(prepared,results,range)
    (length(range)<2 && isempty(eager)) && return nothing
    indices=Int[index for index in range]
    append!(indices,eager)
    Threads.@threads for index in indices
        work=prepared[index]
        sinks=Tuple{Symbol,String}[]
        detached_options=_parallel_task_options(options,base_seed;sinks)
        results[index]=(_parallel_curve_attempt(
            Model._model_mesh_transaction_copy(m),work,m,detached_options,
            base_seed,caller,precise),base_dirty,base_seed,sinks)
    end
    return nothing
end

# Staleness is measured against dirt accumulated after the batch's snapshot —
# dirt the detached copy could not have seen.  Entries already dirty at fill
# time were part of the produce's observed state.  The RNG ledger completes
# the contract as in the other arms: a consumed seed is admissible only when
# the live serial position still equals the recorded entry value.
function _parallel_curve_fresh(work::_ParallelCurveWork,out,
        dirty::Set{Tuple{Int,Int32}},base_dirty::Set{Tuple{Int,Int32}},
        options::Model._ModelMesh1DOptions,base_seed::UInt32)
    out.ok || return false
    for entity in work.footprint
        entity in base_dirty && continue
        entity in dirty && return false
    end
    _,seed_in,seed_out,foreign=out.value
    foreign && return false
    seed_in==base_seed || return false
    seed_out==seed_in || seed_in==options.random_state[] || return false
    return true
end

# Resolve a prepared slot into the produce result the serial arm would have
# returned at this merge position.  A fresh receipt replays the one live-side
# mutation a `:pending` arm owns — `_model_mesh_curve!` drops the stale stored
# discretization before grading (`deMeshGEdge`), which a pending Curve never
# writes back — restores the serial RNG ledger position the produce consumed,
# and returns the transported value for the unchanged emit arm.  Anything
# else re-runs the identical live serial call.
function _parallel_curve_result(m::GeoModel,work::_ParallelCurveWork,entry,
        dirty::Set{Tuple{Int,Int32}},options::Model._ModelMesh1DOptions,
        caller::AbstractString)
    head,base_dirty,base_seed,sinks=entry
    out=head isa Task ? fetch(head) : head
    if !_parallel_curve_fresh(work,out,dirty,base_dirty,options,base_seed)
        return Model._model_mesh_curve!(m,work.tag,options,caller;
            _native_producer=true)
    end
    result,seed_in,seed_out,foreign=out.value
    _PARALLEL_ADOPTED[]+=1
    _parallel_replay_sinks!(options,sinks)
    seed_out==seed_in || (options.random_state[]=seed_out)
    result===:pending && delete!(m.curve_params,work.tag)
    return result
end
