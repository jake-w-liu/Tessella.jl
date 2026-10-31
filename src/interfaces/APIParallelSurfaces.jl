# PAR01 — detached independent-Surface production under `General.NumThreads`.
#
# `_higher_surface_phase!` still owns the serial `_generate_units` traversal,
# the `_higher_append_part!` publication merge and every authority/certificate
# update.  With `General.NumThreads >= 2` the phase additionally detects runs
# of independent singleton Surfaces and produces their
# `_ModelSurfaceProduction` receipts on per-unit transaction copies
# (`Model._model_mesh_transaction_copy`) inside `Threads.@threads`, each task
# writing a disjoint result slot while the shared model, authority and
# assembly stay untouched.  At merge time each slot is checked for staleness:
# an earlier merged unit can invalidate a detached receipt only by changing
# `m.curve_params` for one of its footprint entities (boundary writeback,
# lazy grading or a grown recovered Curve chain) or by moving one of its
# consumed source nodes (`production.source_node_updates`).  A stale or
# unsupported slot is recomputed by the identical live serial call, so the
# produced mesh is bit-identical to the serial one; a fresh slot has its
# detached `m.curve_params` writeback replayed and every freshly-minted UID
# shifted by the amount the serial allocator advanced since the batch ran.
#
# The eligibility contract is deliberately conservative.  Only a singleton
# unit whose serial arm is the plain `mesh_model_surface` call can detach:
# compounds (multi-tag units), periodic slaves and masters, extrusion
# sources and gates, no-new-vertices scopes, QuadTri hybrids, sheets embedded
# in a Region, user size callbacks and geometry-birth overlays all keep the
# original serial arm — their producer contracts consume or publish shared
# state outside the recorded footprint.  A resolved background `size_field`
# is *not* an exclusion: catalog fields are immutable read-only evaluators
# (the lone mutable `ExternalProcessField` serializes every query on its own
# lock), so the same object is threaded into the detached produce exactly
# like `options.inner_field` already is.
#
# Non-contiguous scheduling: eligible slots separated by serial units used to
# run serially (a one-element run never paid for a detached copy).  `_parallel_
# eager_indices` now lets a fill also produce *later* eligible slots whose
# footprints are disjoint from the accumulated write-potential of every unit
# in between — literally "distinct owners, disjoint products" — so a serial
# unit between them can never dirty their receipts.  Each eligible work
# records `gap_writes`: the write-set estimate of the serial slots since the
# previous eligible slot (`nothing` when an in-between serial unit's
# write-set cannot be bounded, which stops the eager scan entirely).  The
# merge-time staleness check remains the sole authority — the predictor only
# decides which slots are worth producing detached, so a wrong guess costs a
# discarded produce, never a wrong mesh.
#
# RNG ledger: every detached task runs under a private `random_state` Ref
# (`_parallel_task_options`) seeded at the batch's serial `Mesh.RandomSeed`
# position and reporting its exit value — the same draw ledger the Region arm
# keeps.  Surface and Curve produces are not expected to touch the seed (only
# `ModelNativeInitialTet`'s kernel and the native hybrid paths draw), but the
# ledger turns that contract from an assumption into a check: an unexpected
# consumption is adopted only when the live serial position still equals the
# observed entry value, otherwise the slot is stale and re-runs serially.
#
# Two determinism details the worker path must preserve:
# * `GmshLibm._win_cli_precision()` is task-local (CLI x87 libm selector);
#   worker tasks re-run the produce under `_with_win_cli_precision` when the
#   flag is set on the orchestrating task, so trig/pow kernels cannot differ
#   between serial and parallel production.
# * Task exceptions are captured into the slot and the affected unit re-runs
#   serially, so a failing producer reports at the same merge position the
#   serial traversal would have reached.
#
# Cross-phase scheduling entry points `_parallel_surface_ready` and
# `_parallel_surface_spawn!` additionally allow a caller to detach ONE
# Surface ahead of the phase barrier once every footprint dep exists — each
# dim-1 entry owning its stored `m.curve_params` discretization and every
# other dep an authority product in the `input` the produce will read ("a
# Surface cannot mesh before its Curves' products exist").  The spawned task
# touches only private state; the merge resolves the `Task` head at the
# unit's serial merge position under the unchanged staleness/ledger checks.

import ..GmshLibm

struct _ParallelSurfaceWork
    tag::Int
    footprint::Set{Tuple{Int,Int32}}
    # Write-set estimate of the serial units since the previous eligible
    # slot (`nothing` = an unbounded serial write precedes).  Read only by
    # `_parallel_eager_indices`; not part of the staleness contract.
    gap_writes::Union{Nothing,Set{Tuple{Int,Int32}}}
end

# A detached task must never touch the caller's shared `random_state` Ref —
# serial ordering reads the same position at merge.  Every other option field
# is shared read-only, so the task-local bundle copies the `Ref` slot only;
# `random_state` is the sole mutable field of `_ModelMesh1DOptions`.  The
# `warn`/`report_error` sinks are the second shared channel: they append to
# `context.exec_warnings`/`error_log` and bump `msg_error_count` on the
# calling context — a shared `push!`/`+=` the worker could race, and the
# messages would interleave at produce time rather than at the unit's serial
# merge position.  `sinks` — one task-private `Vector` per slot — captures
# them for `_parallel_replay_sinks!` at adoption instead.
function _parallel_task_options(options::Model._ModelMesh1DOptions,seed::UInt32;
        sinks=nothing)
    return Model._ModelMesh1DOptions(
        (name===:random_state ? Base.RefValue{UInt32}(seed) :
         name===:warn && sinks!==nothing ?
             (msg->push!(sinks,(:warn,String(msg)))) :
         name===:report_error && sinks!==nothing ?
             (msg->push!(sinks,(:err,String(msg)))) :
            getfield(options,name)
         for name in fieldnames(Model._ModelMesh1DOptions))...)
end

# Replay one adopted receipt's captured diagnostics through the live sinks at
# the unit's serial merge position — exactly where the serial traversal would
# have emitted them, in merge order.  A stale or failed slot drops its buffer
# with the receipt: the live rerun emits whatever the real state produces.
function _parallel_replay_sinks!(options::Model._ModelMesh1DOptions,sinks)
    for (kind,msg) in sinks
        kind===:warn ? Model._model_mesh_warn(options,msg) :
            Model._model_mesh_error(options,msg)
    end
    return nothing
end

# Non-contiguous eligible slots whose footprints are disjoint from the
# accumulated write-potential of every unit ahead of them in this phase may
# be produced in the same batch as the contiguous `range`: a serial unit
# between the batch and the slot writes only inside its own footprint, so a
# disjoint slot observes identical state at merge.  `work.gap_writes` carries
# each eligible slot's preceding serial write-estimate; `nothing` marks an
# unbounded serial writer (a footprint that could not be computed), which
# conservatively ends the scan.  Already-filled slots still contribute their
# footprints.  The merge-time staleness check stays authoritative — this is
# purely a "which copies are worth making" predictor.
function _parallel_eager_indices(prepared,results,range)
    seen=Set{Tuple{Int,Int32}}()
    for index in range
        work=prepared[index]
        work===nothing || union!(seen,work.footprint)
    end
    eager=Int[]
    for index in last(range)+1:length(prepared)
        work=prepared[index]
        work===nothing && continue
        gap=work.gap_writes
        gap===nothing && break
        union!(seen,gap)
        isdisjoint(work.footprint,seen) && results[index]===nothing &&
            push!(eager,index)
        union!(seen,work.footprint)
    end
    return eager
end

# Serial write-set estimate of one ineligible Surface unit — the union of its
# member tags' footprints (the same contract the dirty scan relies on).
# `nothing` marks an unbounded writer and stops the eager scan at the next
# eligible slot.
function _parallel_surface_serial_writes(m::GeoModel,tags)
    acc=Set{Tuple{Int,Int32}}()
    for tag in tags
        union!(acc,_parallel_surface_footprint(m,tag))
    end
    return acc
end

# `General.NumThreads` requests the worker count; `Mesh.MaxNumThreads2D` caps
# the Surface stage when > 0 (upstream's per-dimension limit).
function _parallel_surface_workers()
    workers=_partition_workers()
    workers>=2 || return 1
    cap=get(OPTIONS,"Mesh.MaxNumThreads2D",0.0)
    cap>0 && (workers=min(workers,max(1,floor(Int,cap))))
    return workers
end

# The serial-read serial-write contract of one detached Surface call: every
# entity whose stored discretization or node coordinates the produce can
# consult or rewrite.  Boundary curves come from the surface loops (or a
# discrete Surface record's boundary), embedded dim-0/dim-1 entities ride
# `m.embeds[(2,tag)]`, corner/discrete Curve endpoints ride `m.curves`/
# `m.discrete`, and a dim-1 periodic relation pulls the partner chain into
# `forced` — `_synchronize_periodic_parameters!` propagates partitions along
# the relation graph, so partners belong to the footprint too.  Curve
# writeback (`_model_surface_boundary_writeback!`), lazy grading
# (`_model_surface_mesh_curves!`) and recovered-chain remesh all stay inside
# this set; a dirty mark on any member makes the detached receipt stale.
function _parallel_surface_footprint(m::GeoModel,tag::Int)
    footprint=Set{Tuple{Int,Int32}}()
    curves=Int[]
    if haskey(m.surfaces,tag)
        for loop in m.surfaces[tag],signed in get(m.loops,loop,Int[])
            push!(curves,abs(signed))
        end
    else
        record=get(m.discrete,(2,tag),nothing)
        if record!==nothing
            for (dim,entity) in record.boundary
                dim==1 && push!(curves,abs(entity))
                dim==0 && push!(footprint,(0,Int32(abs(entity))))
            end
        end
    end
    for (dim,entity) in get(m.embeds,(2,tag),NTuple{2,Int}[])
        push!(footprint,(dim,Int32(entity)))
        dim==1 && push!(curves,entity)
    end
    for curve in curves
        push!(footprint,(1,Int32(curve)))
        endpoints=get(m.curves,curve,nothing)
        if endpoints!==nothing
            push!(footprint,(0,Int32(endpoints[1])))
            push!(footprint,(0,Int32(endpoints[2])))
        end
        record=get(m.discrete,(1,curve),nothing)
        if record!==nothing
            for (dim,entity) in record.boundary
                dim==0 && push!(footprint,(0,Int32(abs(entity))))
            end
        end
    end
    for constraint in Model.model_periodic_constraints(m)
        constraint.dim==1 || continue
        slave=(1,constraint.slave_entity);master=(1,constraint.master_entity)
        slave in footprint && push!(footprint,master)
        master in footprint && push!(footprint,slave)
    end
    return footprint
end

# Per-unit eligibility under the phase-level gates.  `nothing` means the
# whole phase keeps the serial traversal (zero added bookkeeping); an entry
# of `nothing` in the vector means that unit runs serially in place.
function _parallel_surface_prepare(m::GeoModel,units,caller::AbstractString;
        metric_surface,scope,periodic_masters,ordinary_sources,
        unmerged_nonew_sources,size_field,geometry_birth)
    _parallel_surface_workers()>=2 || return nothing
    (metric_surface===nothing && scope===nothing &&
        m.meshing.size_callback===nothing) || return nothing
    embedded_in_volumes=Set{Int}()
    for volume in keys(m.volumes),
        (dim,entity) in get(m.embeds,(3,volume),NTuple{2,Int}[])
        dim==2 && push!(embedded_in_volumes,entity)
    end
    compound_members=Set{Int}()
    for pair in m.meshing.compounds,member in pair.second
        push!(compound_members,member)
    end
    prepared=Vector{Any}(nothing,length(units))
    for (index,tags) in enumerate(units)
        length(tags)==1 || continue
        tag=tags[1]
        haskey(m.periodic,(2,tag)) && continue
        tag in periodic_masters && continue
        tag in ordinary_sources && continue
        tag in unmerged_nonew_sources && continue
        tag in embedded_in_volumes && continue
        tag in compound_members && continue
        tag in m.meshing.quad_tri && continue
        Model._extrude_gate(m,2,tag)!==nothing && continue
        prepared[index]=_ParallelSurfaceWork(tag,
            _parallel_surface_footprint(m,tag),nothing)
    end
    any(!isnothing,prepared) || return nothing
    # Second pass: record each eligible slot's preceding serial write-set for
    # the eager scheduler.  `serial_unbounded` marks an in-between unit whose
    # footprint could not be bounded; it stops the scan at the next slot.
    serial_writes=Set{Tuple{Int,Int32}}()
    serial_unbounded=false
    for (index,tags) in enumerate(units)
        work=prepared[index]
        if work===nothing
            if !serial_unbounded
                estimate=try
                    _parallel_surface_serial_writes(m,tags)
                catch err
                    err isa InterruptException && rethrow()
                    nothing
                end
                estimate===nothing ? (serial_unbounded=true) :
                    union!(serial_writes,estimate)
            end
            continue
        end
        prepared[index]=_ParallelSurfaceWork(work.tag,work.footprint,
            serial_unbounded ? nothing : copy(serial_writes))
        empty!(serial_writes);serial_unbounded=false
    end
    return prepared
end

# Content snapshot of the live stored discretizations.  Copies (not the
# stored vectors) guard against an in-place mutation that an identity
# comparison could never see.
function _parallel_surface_param_refs(m::GeoModel)
    return Dict{Int,Vector{Float64}}(curve=>copy(params)
        for (curve,params) in m.curve_params)
end

# `m.curve_params` writes a detached produce can own: boundary writeback,
# lazy grading of missing curves and the native lower-input `merge!` for the
# recorded footprint.  The write scan must cover the *whole* store, not just
# the footprint: a produce writing a slot outside its recorded read set means
# the transport contract was violated — replaying only footprint diffs would
# silently drop live state (and skip the dirty propagation a live write
# owns).  Such writes are flagged `foreign`; the slot re-runs the identical
# serial call instead of adopting an untransportable receipt.  A `nothing`
# value records a deleted entry for replay symmetry.
function _parallel_curve_param_diffs(params::Dict{Int,Vector{Float64}},
        base::Dict{Int,Vector{Float64}},footprint::Set{Tuple{Int,Int32}})
    diffs=Dict{Int,Union{Nothing,Vector{Float64}}}()
    foreign=false
    for (curve,current) in params
        before=get(base,curve,nothing)
        current===before && continue
        (before===nothing || current!=before) || continue
        if (1,Int32(curve)) in footprint
            diffs[curve]=current
        else
            foreign=true
        end
    end
    for curve in keys(base)
        haskey(params,curve) && continue
        (1,Int32(curve)) in footprint ? (diffs[Int(curve)]=nothing) :
            (foreign=true)
    end
    return diffs,foreign
end

# The second half of the transport certificate: a detached produce's
# contract is "read model state, write `curve_params`", so a mutation
# anywhere else in the model stores — entity births through `next_tag`,
# topology, names, constraints, discrete records, kernel records, meshing
# attributes — can never be transported and must send the slot back to the
# serial arm.  The parent's stores are safe to read from a worker: while a
# fill runs, the merge loop is suspended inside the same batch, and during a
# spawned attempt the only live `m` writes are `m.curve_params` (covered by
# the diff scan) and `m.meshing.mesh_authority` (rebound by merges —
# excluded below).  `size_callback` is the one shared field of the
# transaction copy and cannot diverge.  Every record type in these stores
# (`DiscreteEntity`, `_ModelNativeGeoKernel`, `_ModelNativeOCCKernel`,
# `ModelPeriodicConstraint`, ...) defines `==`, so a plain fieldwise
# comparison is the certificate.
function _parallel_model_foreign_writes(working::GeoModel,m::GeoModel)
    for name in (:points,:point_size,:point_w,:curves,:curve_control_points,
            :curve_types,:curve_geometry,:loops,:surfaces,:surface_loops,
            :surface_types,:surface_geometry,:volumes,:entity_names,
            :entity_visibility,:entity_colors,:attributes,:physical,
            :physical_names,:entity_physicals,:physical_tag_max,:box_extents,
            :cylinders,:spheres,:cones,:booleans,:boolean_operands,
            :boolean_components,:periodic,:embeds,:discrete,:next_tag,:name,
            :occ_kernel)
        getfield(working,name)==getfield(m,name) || return true
    end
    wa,ma=working.meshing,m.meshing
    for name in fieldnames(Model.ModelMeshingAttributes)
        name in (:mesh_authority,:size_callback) && continue
        getfield(wa,name)==getfield(ma,name) || return true
    end
    return false
end

function _parallel_replay_param_diffs!(m::GeoModel,diffs)
    for (curve,params) in diffs
        params===nothing ? delete!(m.curve_params,curve) :
            (m.curve_params[curve]=params)
    end
    return nothing
end

# Observability for the detach transport itself: the merge thread increments
# this once per adopted detached receipt (per entity), so tests can prove a
# run took the parallel path rather than silently collapsing to serial.
# Stores happen only inside the in-order merge — never inside workers.
const _PARALLEL_ADOPTED=Ref{Int}(0)
# Companion counter for the Region seed-ledger transport: incremented when
# an adopted detached produce actually advanced the shared `Mesh.RandomSeed`
# ledger (exit seed ≠ entry seed), so tests can prove an RNG-consuming
# Region detached and merged rather than re-running serially.
const _PARALLEL_SEED_ADOPTED=Ref{Int}(0)

# One detached Surface produce on a premade transaction copy.  The task's
# private `random_state` Ref enters at the batch seed and reports its exit
# value — the draw ledger the merge checks.  `foreign` is the transport
# certificate: the produce must write only `m.curve_params` slots inside its
# footprint — a write anywhere else could never reach the live model and
# sends the receipt back to the serial arm instead of being silently
# dropped.  Shared by the `@threads` batch arm and the `Threads.@spawn`
# cross-phase arm.
function _parallel_surface_attempt(working::GeoModel,work::_ParallelSurfaceWork,
        m::GeoModel,input,detached_options::Model._ModelMesh1DOptions,
        base_params::Dict{Int,Vector{Float64}},geom_points,
        hn::UInt64,hc::UInt64,base_seed::UInt32,precise::Bool;
        size_field=nothing)
    out=try
        body=()->begin
            produced=mesh_model_surface(working,work.tag;
                size_field,_extrude_scope=nothing,
                _mesh_input=input,_mesh1d_options=detached_options,
                _mesh_node_high_water=hn,_mesh_element_high_water=hc,
                _return_surface_production=true,_boundary_regrade=false,
                _geometry_points=geom_points)
            diffs,foreign=_parallel_curve_param_diffs(
                working.curve_params,base_params,work.footprint)
            foreign || (foreign=_parallel_model_foreign_writes(working,m))
            (produced,diffs,base_seed,detached_options.random_state[],foreign)
        end
        (;ok=true,value=precise ?
            GmshLibm._with_win_cli_precision(body) : body())
    catch err
        err isa InterruptException && rethrow()
        (;ok=false,value=err)
    end
    return out
end

# Threaded detached production over one contiguous eligible run plus every
# predicted-independent later slot (`_parallel_eager_indices`).  `range`
# shares the batch's serial merge position: `hn`/`hc` fix the allocator base,
# `base_params`/`base_dirty` fix the state the detached copies observed,
# `geom_points` fixes the geometry-birth overlay the produce can consult,
# `base_seed` fixes the `Mesh.RandomSeed` ledger position each task's private
# Ref starts from, and every slot records the produce's outcome for the
# in-order merge.  A batch of one would pay the copy for no concurrency, so
# it falls back to the serial arm by leaving its slot empty.
function _parallel_surface_fill!(results::Vector{Any},range,m::GeoModel,
        prepared,input,options,assembly,dirty::Set{Tuple{Int,Int32}},
        precise::Bool,geometry_birth;size_field=nothing)
    hn=assembly.node_max;hc=assembly.cell_max
    base_params=_parallel_surface_param_refs(m)
    base_dirty=copy(dirty)
    geom_points=_higher_geometry_source_points(assembly,geometry_birth)
    base_seed=options.random_state[]
    eager=_parallel_eager_indices(prepared,results,range)
    (length(range)<2 && isempty(eager)) && return nothing
    indices=Int[index for index in range]
    append!(indices,eager)
    Threads.@threads for index in indices
        work=prepared[index]
        sinks=Tuple{Symbol,String}[]
        detached_options=_parallel_task_options(options,base_seed;sinks)
        results[index]=(_parallel_surface_attempt(
            Model._model_mesh_transaction_copy(m),work,m,input,
            detached_options,base_params,geom_points,hn,hc,base_seed,precise;
            size_field),
            hn,hc,base_dirty,geom_points,base_seed,sinks)
    end
    return nothing
end

# A detached receipt stays usable only when its shape carries actual identity
# rows the merge can remap: a `_ModelSurfaceProduction` without a constructor
# mesh (periodic masters own constructors and are never eligible), a geometry
# production wrapping one, or a plain mesh for the write-back arm.  Anything
# else — an extrusion production, a constructor-carrying receipt — reverts to
# the live serial call rather than smuggling an unmapped product.
function _parallel_surface_remap_supported(produced)
    produced isa Model._ModelSurfaceGeometryProduction &&
        return _parallel_surface_remap_supported(produced.production)
    produced isa Model._ModelSurfaceProduction &&
        return produced.constructor_mesh===nothing
    return produced isa Union{Mesh,MixedMesh}
end

# Shift every UID the detached produce minted past the batch's allocator base
# by the amount the serial allocator advanced to this merge position.
# Consumed identities (`<= hn`/`<= hc`) stay untouched, matching the values a
# live call would have minted from the current high water.
function _parallel_surface_remap_production(produced,dn::UInt64,dc::UInt64,
        hn::UInt64,hc::UInt64)
    produced isa Model._ModelSurfaceGeometryProduction &&
        return Model._ModelSurfaceGeometryProduction(
            _parallel_surface_remap_production(produced.production,dn,dc,hn,hc),
            produced.geometry)
    produced isa Model._ModelSurfaceProduction || return produced
    dn==0 && dc==0 && return produced
    node_tags=dn==0 ? produced.node_tags :
        UInt64[tag>hn ? tag+dn : tag for tag in produced.node_tags]
    element_tags=dc==0 ? produced.element_tags :
        UInt64[tag>hc ? tag+dc : tag for tag in produced.element_tags]
    chains=dn==0 ? produced.curve_chains : Dict{Int32,Vector{UInt64}}(
        curve=>UInt64[tag>hn ? tag+dn : tag for tag in chain]
        for (curve,chain) in produced.curve_chains)
    updates=dn==0 ? produced.source_node_updates :
        Dict{UInt64,NTuple{3,Float64}}((tag>hn ? tag+dn : tag)=>point
            for (tag,point) in produced.source_node_updates)
    return Model._ModelSurfaceProduction(produced.mesh,node_tags,
        produced.node_owners,produced.node_parameters,produced.retained_nodes,
        produced.node_high_water+dn,element_tags,
        produced.element_high_water+dc,produced.birth_rows,
        produced.constructor_mesh,updates,
        produced.native_recombine_events,chains)
end

# Staleness is measured against dirt accumulated after the batch's snapshot —
# dirt the detached copy could not have seen.  Entries already dirty at fill
# time were part of the produce's observed state, as was the recorded
# geometry-birth overlay; a changed `_geometry_points` view is stale too.
# The RNG ledger completes the contract the same way the Region arm does:
# `seed_in==base_seed` certifies the task observed the recorded ledger
# position; a consumed seed is admissible only when the live serial position
# still equals that entry value (the exit state is then written back at
# merge), otherwise the slot is stale and re-runs serially.
function _parallel_surface_fresh(work::_ParallelSurfaceWork,out,
        dirty::Set{Tuple{Int,Int32}},base_dirty::Set{Tuple{Int,Int32}},
        assembly,geometry_birth,geom_points,options,base_seed::UInt32)
    out.ok || return false
    for entity in work.footprint
        entity in base_dirty && continue
        entity in dirty && return false
    end
    geom_points===nothing ||
        _higher_geometry_source_points(assembly,geometry_birth)==geom_points ||
        return false
    produced,_,seed_in,seed_out,foreign=out.value
    foreign && return false
    _parallel_surface_remap_supported(produced) || return false
    seed_in==base_seed || return false
    seed_out==seed_in || seed_in==options.random_state[] || return false
    return true
end

# Resolve a prepared slot into the unit's `products` list.  A fresh receipt
# replays its detached `m.curve_params` writeback onto the live model (the
# serial produce would have written it during the call), restores the serial
# RNG ledger position the produce consumed, and remaps its fresh UIDs;
# anything else re-runs the identical serial produce so the merge sees
# exactly the state it always owned.  A `Task` head (`Threads.@spawn`
# cross-phase arms) is fetched at the slot's merge position.
function _parallel_surface_products(m::GeoModel,work::_ParallelSurfaceWork,
        entry,dirty::Set{Tuple{Int,Int32}},assembly;
        input,scope,options,size_field,geometry_birth,caller)
    head,hn,hc,base_dirty,geom_points,base_seed,sinks=entry
    out=head isa Task ? fetch(head) : head
    if !_parallel_surface_fresh(work,out,dirty,base_dirty,assembly,
            geometry_birth,geom_points,options,base_seed)
        return [(work.tag,mesh_model_surface(m,work.tag;
            size_field=size_field,_extrude_scope=scope,_mesh_input=input,
            _mesh1d_options=options,_mesh_node_high_water=assembly.node_max,
            _mesh_element_high_water=assembly.cell_max,
            _return_surface_production=true,_boundary_regrade=false,
            _geometry_points=_higher_geometry_source_points(
                assembly,geometry_birth)))]
    end
    produced,diffs,seed_in,seed_out,foreign=out.value
    _PARALLEL_ADOPTED[]+=1
    _parallel_replay_sinks!(options,sinks)
    _parallel_replay_param_diffs!(m,diffs)
    seed_out==seed_in || (options.random_state[]=seed_out)
    return [(work.tag,_parallel_surface_remap_production(produced,
        assembly.node_max-hn,assembly.cell_max-hc,hn,hc))]
end

# Every mutation a merged unit leaves on the shared state a later detached
# produce reads: stored discretization changes (boundary writeback, lazy
# grading, remesh parameter sync) and actual source-node relocations through
# `source_node_updates` (whose owner entity carries the changed coordinates).
function _parallel_surface_mark_dirty!(dirty::Set{Tuple{Int,Int32}},
        m::GeoModel,assembly,products,refs::Dict{Int,Vector{Float64}})
    for (curve,params) in m.curve_params
        before=get(refs,curve,nothing)
        before!==nothing && before==params && continue
        push!(dirty,(1,Int32(curve)))
    end
    for (_,produced) in products
        binding=produced isa Model._ModelSurfaceGeometryProduction ?
            produced.production :
            produced isa Model._ModelSurfaceProduction ? produced : nothing
        binding===nothing && continue
        for uid in keys(binding.source_node_updates)
            row=get(assembly.node_positions,uid,Int32(0))
            row>0 && push!(dirty,assembly.nodes[row].owner)
        end
    end
    return nothing
end

# Cross-phase readiness, the Surface analogue of `_parallel_volume_ready`:
# the ordering constraint itself — a Surface cannot mesh before its Curve
# products exist.  Every dim-1 footprint entry (boundary loops, embedded
# Curves and their periodic partners, all pulled in by
# `_parallel_surface_footprint`) must already own a stored discretization on
# the live model — the produce would otherwise lazily grade it on the
# detached copy and the receipt could only fall back at merge once the live
# Curve merged anyway.  Every non-dim-1 dep (curve endpoints, discrete
# boundary Points, embedded Points) must already own an authority product in
# the exact `input` the detached produce will read.  Run on the caller's
# merge task — `m.curve_params` and `input` are merge-task state.
function _parallel_surface_ready(m::GeoModel,input,work::_ParallelSurfaceWork)
    input isa Model._ModelMeshAuthority || return false
    for dep in work.footprint
        if dep[1]==1
            haskey(m.curve_params,Int(dep[2])) || return false
        else
            haskey(input.products,dep) || return false
        end
    end
    return true
end

# Cross-phase scheduling entry point: spawn ONE detached Surface produce.
# The transaction copy and every snapshot are taken on the *calling* task —
# the spawned task then touches only private state, so it may run while the
# serial merge of an earlier phase (Curves, for a Surface scheduled ahead of
# the phase barrier) continues; a copy made inside the task would race the
# live merge.  The slot is shaped exactly like a `_parallel_surface_fill!`
# slot, so `_parallel_surface_products` resolves it at the unit's serial
# merge position under the unchanged freshness/ledger checks.  The caller
# owns the ordering contract: spawn only once every lower producer the unit
# can consult (its boundary Curves' `curve_params`, embedded-entity products)
# has been merged, and thread ONE `dirty` set through both phases so dirt
# accumulated between the spawn and the merge still invalidates the receipt.
function _parallel_surface_spawn!(results::Vector{Any},index::Int,m::GeoModel,
        work::_ParallelSurfaceWork,input,options,assembly,
        dirty::Set{Tuple{Int,Int32}},geometry_birth;precise::Bool,
        size_field=nothing)
    working=Model._model_mesh_transaction_copy(m)
    hn=assembly.node_max;hc=assembly.cell_max
    base_params=_parallel_surface_param_refs(m)
    geom_points=_higher_geometry_source_points(assembly,geometry_birth)
    base_seed=options.random_state[]
    sinks=Tuple{Symbol,String}[]
    detached_options=_parallel_task_options(options,base_seed;sinks)
    task=Threads.@spawn _parallel_surface_attempt(working,work,m,input,
        detached_options,base_params,geom_points,hn,hc,base_seed,precise;
        size_field)
    results[index]=(task,hn,hc,copy(dirty),geom_points,base_seed,sinks)
    return nothing
end
