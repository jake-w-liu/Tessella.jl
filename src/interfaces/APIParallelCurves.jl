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
# dim-1 compound members and Curves carrying a discrete/attached record all
# keep the original serial arm — their producer contracts consume or publish
# shared state outside the recorded footprint.
#
# Unlike the Region arm there is no seed ledger: `options.random_state` feeds
# only `ModelNativeInitialTet`'s perturbation kernel and the native hybrid
# batch/components paths — no 1-D grading reach it — so the shared options
# bundle is read-only for detached Curves.  Two determinism details the
# worker path must preserve, as in the Surface arm:
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
    # A stateful size callback or an opaque background field is consulted from
    # inside the grading kernel; neither crosses a transaction copy safely.
    (options.callback===nothing && options.inner_field===nothing &&
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
            _parallel_curve_footprint(m,curve))
    end
    any(!isnothing,prepared) || return nothing
    return prepared
end

# Threaded detached production over one contiguous eligible run of a pass.
# `range` shares the batch's serial merge position: `base_dirty` fixes the
# dirt the detached copies could already see, and every slot records the
# produce's outcome for the in-order merge.  A run of one would pay the copy
# for no concurrency, so it falls back to the serial arm by leaving its slot
# empty.  The produce mints no public UIDs and writes only the working copy's
# own `curve_params` slot — the merge replays the identical writeback and
# emits on the live assembly, so no UID remap ledger is needed.
function _parallel_curve_fill!(results::Vector{Any},range,m::GeoModel,
        prepared,options::Model._ModelMesh1DOptions,
        dirty::Set{Tuple{Int,Int32}},precise::Bool,caller::AbstractString)
    length(range)<2 && return nothing
    base_dirty=copy(dirty)
    Threads.@threads for offset in eachindex(range)
        index=range[offset]
        work=prepared[index]
        out=try
            body=()->begin
                working=Model._model_mesh_transaction_copy(m)
                return Model._model_mesh_curve!(working,work.tag,options,caller;
                    _native_producer=true)
            end
            (;ok=true,value=precise ?
                GmshLibm._with_win_cli_precision(body) : body())
        catch err
            err isa InterruptException && rethrow()
            (;ok=false,value=err)
        end
        results[index]=(out,base_dirty)
    end
    return nothing
end

# Staleness is measured against dirt accumulated after the batch's snapshot —
# dirt the detached copy could not have seen.  Entries already dirty at fill
# time were part of the produce's observed state.
function _parallel_curve_fresh(work::_ParallelCurveWork,out,
        dirty::Set{Tuple{Int,Int32}},base_dirty::Set{Tuple{Int,Int32}})
    out.ok || return false
    for entity in work.footprint
        entity in base_dirty && continue
        entity in dirty && return false
    end
    return true
end

# Resolve a prepared slot into the produce result the serial arm would have
# returned at this merge position.  A fresh receipt replays the one live-side
# mutation a `:pending` arm owns — `_model_mesh_curve!` drops the stale stored
# discretization before grading (`deMeshGEdge`), which a pending Curve never
# writes back — and returns the transported value for the unchanged emit arm.
# Anything else re-runs the identical live serial call.
function _parallel_curve_result(m::GeoModel,work::_ParallelCurveWork,entry,
        dirty::Set{Tuple{Int,Int32}},options::Model._ModelMesh1DOptions,
        caller::AbstractString)
    out,base_dirty=entry
    if !_parallel_curve_fresh(work,out,dirty,base_dirty)
        return Model._model_mesh_curve!(m,work.tag,options,caller;
            _native_producer=true)
    end
    result=out.value
    result===:pending && delete!(m.curve_params,work.tag)
    return result
end
