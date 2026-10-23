# PAR01 — detached independent-Region production under `General.NumThreads`.
#
# `_higher_volume_phase!` still owns the serial `sort!(meshed)` traversal, the
# `_higher_append_part!` publication merge and every authority/certificate
# update.  With `General.NumThreads >= 2` the phase additionally detects runs
# of independent singleton Regions and produces their
# `_ModelVolumeGeometryProduction` receipts on per-unit transaction copies
# (`Model._model_mesh_transaction_copy`) inside `Threads.@threads`, each task
# writing a disjoint result slot while the shared model, authority and
# assembly stay untouched.  At merge time each slot is checked for staleness:
# an earlier merged unit can invalidate a detached receipt only by changing
# `m.curve_params` for one of its footprint entities (embedded-curve lazy
# grading or an embedded/missing-boundary Surface produce's writeback).  A
# stale or unsupported slot is recomputed by the identical live serial call,
# so the produced mesh is bit-identical to the serial one; a fresh slot has
# its detached `m.curve_params` writeback replayed and every freshly-minted
# UID shifted by the amount the serial allocator advanced since the batch ran.
#
# The eligibility contract is deliberately conservative, mirroring the
# Surface arm.  Only a singleton Region whose serial arm is the plain
# `mesh_model_volume` call can detach: extrusion gates and staged subdivision
# (`extruded`/`scope`), the native hybrid batch and per-Region hybrid
# constructor, transfinite Regions, dim-3 compound members, Regions whose
# boundary or embedded Surfaces sit in a dim-2 periodic relation, size fields
# and user size callbacks all keep the original serial arm — their producer
# contracts consume or publish shared state outside the recorded footprint.
#
# Non-contiguous scheduling mirrors the Surface arm: `_parallel_eager_indices`
# lets a fill also produce later eligible slots whose footprints are disjoint
# from the accumulated write-potential of every unit in between, so a serial
# Region between them can never dirty their receipts.  Each eligible work
# records `gap_writes` (the write-set estimate of the serial slots since the
# previous eligible slot; `nothing` when an in-between serial unit cannot be
# bounded, which stops the scan).  The merge-time staleness check remains the
# sole authority.
#
# Cross-phase scheduling entry points `_parallel_volume_ready` and
# `_parallel_volume_spawn!` detach a single Region ahead of the volume phase
# once every boundary/embedded product it can consult is present in the
# authority it will read — the volume analogue of "a Surface cannot mesh
# before its Curves' products exist".  The caller threads one `dirty` set
# through both phases so surface-phase merges still invalidate the receipt.
#
# A second serial-threading detail has no Surface analogue: the shared
# `_ModelMesh1DOptions.random_state` Ref seeds `ModelNativeInitialTet`'s
# perturbation kernel on the native lower-authority fill (`native_input`) and
# advances by a seed-dependent draw count, so later Regions can never replay
# the serial seed chain out of order.  `_parallel_volume_rng_free` predicts at
# prepare time which Regions provably avoid that kernel — implicit factory
# solids (analytic boundary + seeded-local `tetrahedralize`), non-native
# explicit shells (a boundary product absent from the authority kills
# `native_identity` and selects the `perturb=false` recovery path) and
# factory-encoded shells without embedded constraints.  Each task still runs
# under a private `options` copy whose `random_state` Ref starts at the batch
# seed and records its exit value; at merge a consumed seed is admissible only
# when the live serial state still equals the observed entry value, in which
# case the exit state is written back — the same interleaved RNG ledger the
# serial arm keeps.  A consumed-but-mis-seeded slot is stale and re-runs.
#
# Two determinism details the worker path must preserve:
# * `GmshLibm._win_cli_precision()` is task-local (CLI x87 libm selector);
#   worker tasks re-run the produce under `_with_win_cli_precision` when the
#   flag is set on the orchestrating task, so trig/pow kernels cannot differ
#   between serial and parallel production.
# * Task exceptions are captured into the slot and the affected unit re-runs
#   serially, so a failing producer reports at the same merge position the
#   serial traversal would have reached.

import ..GmshLibm

struct _ParallelVolumeWork
    tag::Int
    footprint::Set{Tuple{Int,Int32}}
    # Write-set estimate of the serial units since the previous eligible
    # slot (`nothing` = an unbounded serial write precedes) — read only by
    # `_parallel_eager_indices`, not part of the staleness contract.
    gap_writes::Union{Nothing,Set{Tuple{Int,Int32}}}
    # Authority products the produce can consult: its boundary Surfaces plus
    # every entity embedded in it.  `_parallel_volume_ready` requires each
    # one present before a cross-phase spawn is ordered — the detached
    # produce must never need a serial arm's just-emitted product.
    deps::Vector{Tuple{Int,Int32}}
end

# `General.NumThreads` requests the worker count; `Mesh.MaxNumThreads3D` caps
# the Region stage when > 0 (upstream's per-dimension limit).
function _parallel_volume_workers()
    workers=_partition_workers()
    workers>=2 || return 1
    cap=get(OPTIONS,"Mesh.MaxNumThreads3D",0.0)
    cap>0 && (workers=min(workers,max(1,floor(Int,cap))))
    return workers
end

# The serial-read serial-write contract of one detached Region call: every
# entity whose stored discretization or node coordinates the produce can
# consult or rewrite.  Boundary Surface meshes come from the shared authority
# (`input`) — never re-meshed when present — but a missing product falls back
# to a live `mesh_model_surface` call, so each boundary and embedded Surface
# contributes its full `_parallel_surface_footprint`.  Volume embeds add their
# own entities; an embedded Curve's endpoint Points, discrete boundary Points
# and — transitively, since `_volume_embedded_curve_params!` walks master
# chains — every dim-1 periodic partner join the footprint.  A dirty mark on
# any member makes the detached receipt stale.
function _parallel_volume_footprint(m::GeoModel,tag::Int,boundary)
    footprint=Set{Tuple{Int,Int32}}()
    push!(footprint,(3,Int32(tag)))
    surfaces=Set{Int}()
    for signed in boundary
        push!(surfaces,abs(signed))
    end
    embedded_curves=Int[]
    for (dim,entity) in get(m.embeds,(3,tag),NTuple{2,Int}[])
        push!(footprint,(dim,Int32(entity)))
        dim==2 && push!(surfaces,entity)
        dim==1 && push!(embedded_curves,entity)
    end
    for surface in surfaces
        union!(footprint,_parallel_surface_footprint(m,surface))
    end
    for curve in embedded_curves
        endpoints=get(m.curves,curve,nothing)
        endpoints===nothing || (push!(footprint,(0,Int32(endpoints[1])));
            push!(footprint,(0,Int32(endpoints[2]))))
        record=get(m.discrete,(1,curve),nothing)
        record===nothing || for (dim0,entity0) in record.boundary
            dim0==0 && push!(footprint,(0,Int32(abs(entity0))))
        end
    end
    changed=true
    while changed
        changed=false
        for constraint in Model.model_periodic_constraints(m)
            constraint.dim==1 || continue
            slave=(1,constraint.slave_entity);master=(1,constraint.master_entity)
            if slave in footprint && !(master in footprint)
                push!(footprint,master);changed=true
            end
            if master in footprint && !(slave in footprint)
                push!(footprint,slave);changed=true
            end
        end
    end
    return footprint
end

# Serial write-set estimate of one ineligible Region — its full footprint
# (the same contract the dirty scan relies on).  `nothing` marks an
# unbounded writer (the boundary could not even be listed), which stops the
# eager scan at the next eligible slot.
function _parallel_volume_serial_writes(m::GeoModel,tag::Int,caller)
    boundary=try
        Model._model_volume_boundary_surfaces(m,tag,caller)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
    return try
        _parallel_volume_footprint(m,tag,boundary)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
end

# Whether the serial produce is predicted to keep `options.random_state`
# untouched.  The shared seed is consumed only by `ModelNativeInitialTet`'s
# kernel: the native-initial fill (`native_input` — every boundary Surface
# resolving to an actual entity-data product through `_model_mesh_product`,
# modulo factory carve-outs) and the per-Region native hybrid path.  Every
# other fill — analytic factory surfaces through `tetrahedralize`, the
# `perturb=false` `delaunay3d`/`recover_boundary` fallbacks — uses local seeds
# and is order-independent.  A misprediction is still safe: the merge-time
# seed check on the task's private Ref catches it and re-runs serially.
# `false` therefore only costs the detached produce.
function _parallel_volume_rng_free(m::GeoModel,input,tag::Int,boundary,caller)
    # `explicit_geometry===nothing` — implicit factory solids and empty shell
    # lists never build the PLC-joined boundary, so the kernel cannot run.
    (isempty(m.volumes[tag]) || Model._implicit_volume_surface(m,tag)) &&
        return true
    factory=try
        Model._model_factory_volume_encoding(m,tag)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
    # The per-Region hybrid constructor is probed ahead of the explicit shell
    # and drives the same kernel.  It requires every boundary product present
    # with only Tri/Quad blocks and at least one nonempty Quad block.
    if !factory
        complete=true;quads=false
        for signed in boundary
            product=get(input.products,(2,Int32(abs(signed))),nothing)
            product===nothing && (complete=false;break)
            for block in product.blocks
                block.msh in (2,3) || (complete=false;break)
                quads|=block.msh==3 && !isempty(block.nodes)
            end
            complete || break
        end
        if complete && quads
            hybrid=try
                Model._model_native_hybrid_volume_input(m,tag,input,caller)
            catch err
                err isa InterruptException && rethrow()
                return false
            end
            hybrid===nothing || return false
        end
    end
    # `native_identity` holds when every boundary Surface supplies an actual
    # authority product (`_model_mesh_product` rebinds each to a MixedMesh
    # with entity data); any miss makes the produce fall back to a fresh
    # `mesh_model_surface`/OCC kernel result, a plain Mesh that disables the
    # native-initial fill — except an OCC curved Surface could still return
    # identity metadata, a misprediction the merge-time seed check absorbs.
    all_native=true
    for signed in boundary
        product=try
            Model._model_mesh_product(input,2,abs(signed);caller,
                producer_order=true)
        catch err
            err isa InterruptException && rethrow()
            return false
        end
        product===nothing && (all_native=false;break)
    end
    all_native || return true
    factory || return false
    # `native_input` needs `!factory_encoding || factory_constraints`: a
    # factory-encoded shell only takes the native path when the Region or a
    # declared Surface carries embedded entities.
    isempty(get(m.embeds,(3,tag),NTuple{2,Int}[])) || return false
    surfaces=Set{Int}(abs.(boundary))
    for (dim,entity) in get(m.embeds,(3,tag),NTuple{2,Int}[])
        dim==2 && push!(surfaces,entity)
    end
    for surface in surfaces
        nested=try
            Model._model_surface_embedding_tags(m,surface,caller)
        catch err
            err isa InterruptException && rethrow()
            return false
        end
        (isempty(nested[1]) && isempty(nested[2])) || return false
    end
    return true
end

# Per-Region eligibility under the phase-level gates.  `nothing` means the
# whole phase keeps the serial traversal (zero added bookkeeping); an entry
# of `nothing` in the vector means that Region runs serially in place.
function _parallel_volume_prepare(m::GeoModel,volume_tags,input,caller;
        size_field,scope,extruded,hybrid_batch,staged)
    _parallel_volume_workers()>=2 || return nothing
    (size_field===nothing && scope===nothing && extruded===nothing &&
        hybrid_batch===nothing && !staged &&
        m.meshing.size_callback===nothing) || return nothing
    input isa Model._ModelMeshAuthority || return nothing
    m.meshing.algorithm_3d in (1,3) || return nothing
    compound_members=Set{Int}()
    for pair in m.meshing.compounds
        pair.first==3 || continue
        append!(compound_members,pair.second)
    end
    prepared=Vector{Any}(nothing,length(volume_tags))
    for (index,tag) in enumerate(volume_tags)
        haskey(m.meshing.transfinite_volumes,tag) && continue
        Model._extrude_gate(m,3,tag)!==nothing && continue
        tag in compound_members && continue
        boundary=try
            Model._model_volume_boundary_surfaces(m,tag,caller)
        catch err
            err isa InterruptException && rethrow()
            continue
        end
        surfaces=Set{Int}(abs.(boundary))
        for (dim,entity) in get(m.embeds,(3,tag),NTuple{2,Int}[])
            dim==2 && push!(surfaces,entity)
        end
        # A dim-2 periodic relation on any of the Region's Surfaces routes the
        # produce through slave-facet synchronization and consumes the
        # operation's actual copy pairs — shared state outside the footprint.
        # A dim-3 relation on the Region itself is likewise excluded.
        any(constraint->constraint.dim==2 &&
            (Int(constraint.slave_entity) in surfaces ||
             Int(constraint.master_entity) in surfaces),
            Model.model_periodic_constraints(m)) && continue
        any(constraint->constraint.dim==3 &&
            (Int(constraint.slave_entity)==tag ||
             Int(constraint.master_entity)==tag),
            Model.model_periodic_constraints(m)) && continue
        footprint=try
            _parallel_volume_footprint(m,tag,boundary)
        catch err
            err isa InterruptException && rethrow()
            continue
        end
        _parallel_volume_rng_free(m,input,tag,boundary,caller) || continue
        deps=Tuple{Int,Int32}[(2,Int32(abs(signed))) for signed in boundary]
        for (dim,entity) in get(m.embeds,(3,tag),NTuple{2,Int}[])
            push!(deps,(dim,Int32(entity)))
        end
        prepared[index]=_ParallelVolumeWork(tag,footprint,nothing,deps)
    end
    any(!isnothing,prepared) || return nothing
    # Second pass: record each eligible slot's preceding serial write-set for
    # the eager scheduler.  `serial_unbounded` marks an in-between Region
    # whose footprint could not be bounded; it stops the scan at the next
    # eligible slot.
    serial_writes=Set{Tuple{Int,Int32}}()
    serial_unbounded=false
    for (index,tag) in enumerate(volume_tags)
        work=prepared[index]
        if work===nothing
            if !serial_unbounded
                estimate=_parallel_volume_serial_writes(m,tag,caller)
                estimate===nothing ? (serial_unbounded=true) :
                    union!(serial_writes,estimate)
            end
            continue
        end
        prepared[index]=_ParallelVolumeWork(work.tag,work.footprint,
            serial_unbounded ? nothing : copy(serial_writes),work.deps)
        empty!(serial_writes);serial_unbounded=false
    end
    return prepared
end

# One detached Region produce on a premade transaction copy.  The task's
# private `random_state` Ref enters at the batch seed and reports its exit
# value — the draw ledger the merge checks.  Shared by the `@threads` batch
# arm and the `Threads.@spawn` cross-phase arm.
function _parallel_volume_attempt(working::GeoModel,work::_ParallelVolumeWork,
        input,detached_options::Model._ModelMesh1DOptions,
        geom_points,base_params::Dict{Int,Vector{Float64}},
        base_seed::UInt32,hn::UInt64,hc::UInt64;
        periodic_surface_pairs,optimize_pyramids,precise::Bool)
    out=try
        body=()->begin
            produced=mesh_model_volume(working,work.tag;
                size_field=nothing,_extrude_scope=nothing,
                _extrude_working_model=true,_mesh_input=input,
                _mesh1d_options=detached_options,
                _return_volume_production=true,
                _periodic_surface_pairs=periodic_surface_pairs,
                _optimize_pyramids=optimize_pyramids,
                _mesh_node_high_water=hn,_mesh_element_high_water=hc,
                _geometry_points=geom_points)
            (produced,_parallel_curve_param_diffs(
                working.curve_params,base_params,work.footprint),
                base_seed,detached_options.random_state[])
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
# `base_seed` fixes the `Mesh.RandomSeed` ledger position the serial arm would
# consume from, `geom_points` fixes the geometry-birth overlay, and every slot
# records the produce's outcome plus its private RNG entry/exit state for the
# in-order merge.  A batch of one would pay the copy for no concurrency, so
# it falls back to the serial arm by leaving its slot empty.
function _parallel_volume_fill!(results::Vector{Any},range,m::GeoModel,
        prepared,input,options,assembly,dirty::Set{Tuple{Int,Int32}},
        precise::Bool,geometry_birth;periodic_surface_pairs,optimize_pyramids)
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
        detached_options=_parallel_task_options(options,base_seed)
        results[index]=(_parallel_volume_attempt(
            Model._model_mesh_transaction_copy(m),work,input,detached_options,
            geom_points,base_params,base_seed,hn,hc;
            periodic_surface_pairs,optimize_pyramids,precise),
            hn,hc,base_dirty,geom_points,base_seed)
    end
    return nothing
end

# A detached receipt stays usable only when it carries the ordinary geometry
# production the serial arm would have returned: transfinite, hybrid, swept
# and plain-mesh receipts revert to the live serial call rather than
# smuggling an unmapped product past the allocator ledger.
# A boundary write-back mutates shared Surface/Curve products and stored
# curve discretizations on the detached copy — that state cannot be remapped
# into the live merge, so such receipts revert to the serial arm like the
# transfinite/hybrid/swept cases.
_parallel_volume_remap_supported(produced)=
    produced isa Model._ModelVolumeGeometryProduction &&
    produced.boundary_writeback===nothing

# Shift every UID the detached produce minted past the batch's allocator base
# by the amount the serial allocator advanced to this merge position.
# Consumed boundary identities (`<= hn`) stay untouched, matching the values a
# live call would have minted from the current high water.
function _parallel_volume_remap_production(produced,dn::UInt64,dc::UInt64,
        hn::UInt64,hc::UInt64)
    produced isa Model._ModelVolumeGeometryProduction || return produced
    dn==0 && dc==0 && return produced
    boundary_tags=dn==0 ? produced.boundary_node_tags :
        UInt64[tag>hn ? tag+dn : tag for tag in produced.boundary_node_tags]
    element_tags=dc==0 ? produced.element_tags :
        UInt64[tag>hc ? tag+dc : tag for tag in produced.element_tags]
    return Model._ModelVolumeGeometryProduction(produced.mesh,boundary_tags,
        produced.node_high_water+dn,produced.temporary_elements,
        element_tags,produced.element_high_water+dc,produced.boundary_writeback)
end

# Staleness is measured against dirt accumulated after the batch's snapshot —
# dirt the detached copy could not have seen.  Entries already dirty at fill
# time were part of the produce's observed state, as was the recorded
# geometry-birth overlay.  Two volume-specific receipts checks complete the
# contract: the RNG ledger position (`seed_out==seed_in` certifies an
# untouched shared seed; a consumed seed must equal the live serial position
# and is written back at merge) and the node high water (`== hn + fresh`
# certifies the initial kernel's absolute `count+8` numbering never lifted the
# base — a lifted base would not commute with the +dn remap).
function _parallel_volume_fresh(work::_ParallelVolumeWork,out,
        dirty::Set{Tuple{Int,Int32}},base_dirty::Set{Tuple{Int,Int32}},
        assembly,geometry_birth,geom_points,options,base_seed,hn,hc)
    out.ok || return false
    for entity in work.footprint
        entity in base_dirty && continue
        entity in dirty && return false
    end
    geom_points===nothing ||
        _higher_geometry_source_points(assembly,geometry_birth)==geom_points ||
        return false
    produced,diffs,seed_in,seed_out=out.value
    _parallel_volume_remap_supported(produced) || return false
    seed_in==base_seed || return false
    seed_out==seed_in || seed_in==options.random_state[] || return false
    fresh_nodes=count(>(hn),produced.boundary_node_tags)
    produced.node_high_water==hn+fresh_nodes || return false
    return true
end

# Resolve a prepared slot into the unit's produce receipt.  A fresh receipt
# replays its detached `m.curve_params` writeback onto the live model (the
# serial produce would have written it during the call), restores the serial
# RNG ledger position the produce consumed, and remaps its fresh UIDs;
# anything else re-runs the identical serial produce so the merge sees exactly
# the state it always owned.  A `Task` head (`Threads.@spawn` cross-phase
# arms) is fetched at the slot's merge position.
function _parallel_volume_product(m::GeoModel,work::_ParallelVolumeWork,
        entry,dirty::Set{Tuple{Int,Int32}},assembly,options;
        input,scope,size_field,periodic_surface_pairs,optimize_pyramids,
        geometry_birth,caller)
    head,hn,hc,base_dirty,geom_points,base_seed=entry
    out=head isa Task ? fetch(head) : head
    if !_parallel_volume_fresh(work,out,dirty,base_dirty,assembly,
            geometry_birth,geom_points,options,base_seed,hn,hc)
        return mesh_model_volume(m,work.tag;size_field,
            _extrude_scope=scope,_extrude_working_model=true,_mesh_input=input,
            _mesh1d_options=options,_return_volume_production=true,
            _periodic_surface_pairs=periodic_surface_pairs,
            _optimize_pyramids=optimize_pyramids,
            _mesh_node_high_water=assembly.node_max,
            _mesh_element_high_water=assembly.cell_max,
            _geometry_points=_higher_geometry_source_points(
                assembly,geometry_birth))
    end
    produced,diffs,seed_in,seed_out=out.value
    _parallel_replay_param_diffs!(m,diffs)
    seed_out==seed_in || (options.random_state[]=seed_out)
    return _parallel_volume_remap_production(produced,
        assembly.node_max-hn,assembly.cell_max-hc,hn,hc)
end

# Every mutation a merged unit leaves on the shared state a later detached
# produce reads: stored discretization changes (embedded-curve grading and
# boundary writeback replayed from receipts, or written live by a serial arm).
# Region merges never relocate consumed source nodes — the volume receipts
# carry no `source_node_updates` — so only the stored Curve discretizations
# need the diff scan.
function _parallel_volume_mark_dirty!(dirty::Set{Tuple{Int,Int32}},
        m::GeoModel,refs::Dict{Int,Vector{Float64}})
    for (curve,params) in m.curve_params
        before=get(refs,curve,nothing)
        before!==nothing && before==params && continue
        push!(dirty,(1,Int32(curve)))
    end
    for curve in keys(refs)
        haskey(m.curve_params,curve) || push!(dirty,(1,Int32(curve)))
    end
    return nothing
end

# ── Cross-phase scheduling ──────────────────────────────────────────────────
# The ordering constraint lifted one dimension: a Region cannot detach before
# every product it can consult exists — its boundary Surface meshes plus the
# embedded entities (`work.deps`).  `_parallel_volume_ready` is that check:
# every dep must already own an authority product in the exact `input` the
# detached produce will read, so the produce never reaches the live
# `mesh_model_surface` fallback (which would mint consumed-side UIDs the
# merge remap cannot re-derive).  `_parallel_volume_spawn!` then spawns the
# produce: the transaction copy and all snapshots are taken on the *calling*
# task — a copy made inside the spawned task would race the live merge that
# keeps mutating `m` — and the slot is shaped exactly like a
# `_parallel_volume_fill!` slot, so `_parallel_volume_product` resolves it at
# the Region's serial merge position under the unchanged footprint/geometry/
# seed-ledger checks.
#
# The caller's obligations when wiring a surface↔volume overlap:
# * pass the authority (`input`) the serial arm would read — products for a
#   merged Surface are stable across rebuilds (they derive from that
#   Surface's own emitted cells), so a ready Region sees final data;
# * thread ONE `dirty` set through both phases (union the surface phase's
#   `parallel_dirty` into the volume phase's before its first merge) so
#   surface-phase merges after the spawn still invalidate the receipt;
# * spawn on the phase's merge task only — snapshots must be atomic with the
#   merge's own reads.
function _parallel_volume_ready(input,work::_ParallelVolumeWork)
    input isa Model._ModelMeshAuthority || return false
    for dep in work.deps
        haskey(input.products,dep) || return false
    end
    return true
end

function _parallel_volume_spawn!(results::Vector{Any},index::Int,m::GeoModel,
        work::_ParallelVolumeWork,input,options,assembly,
        dirty::Set{Tuple{Int,Int32}},geometry_birth;
        periodic_surface_pairs,optimize_pyramids,precise::Bool)
    working=Model._model_mesh_transaction_copy(m)
    hn=assembly.node_max;hc=assembly.cell_max
    base_params=_parallel_surface_param_refs(m)
    geom_points=_higher_geometry_source_points(assembly,geometry_birth)
    base_seed=options.random_state[]
    detached_options=_parallel_task_options(options,base_seed)
    task=Threads.@spawn _parallel_volume_attempt(working,work,input,
        detached_options,geom_points,base_params,base_seed,hn,hc;
        periodic_surface_pairs,optimize_pyramids,precise)
    results[index]=(task,hn,hc,copy(dirty),geom_points,base_seed)
    return nothing
end
