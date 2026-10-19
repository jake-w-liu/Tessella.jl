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
# in a Region, size fields, user size callbacks and geometry-birth overlays
# all keep the original serial arm — their producer contracts consume or
# publish shared state outside the recorded footprint.
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

struct _ParallelSurfaceWork
    tag::Int
    footprint::Set{Tuple{Int,Int32}}
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
    (metric_surface===nothing && size_field===nothing && scope===nothing &&
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
            _parallel_surface_footprint(m,tag))
    end
    any(!isnothing,prepared) || return nothing
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
# lazy grading of missing curves and the native lower-input `merge!`.  All
# land inside the footprint; anything outside is the serial arm's business.
# A `nothing` value records a deleted entry for replay symmetry.
function _parallel_curve_param_diffs(params::Dict{Int,Vector{Float64}},
        base::Dict{Int,Vector{Float64}},footprint::Set{Tuple{Int,Int32}})
    diffs=Dict{Int,Union{Nothing,Vector{Float64}}}()
    for (dim,curve) in footprint
        dim==1 || continue
        key=Int(curve)
        current=get(params,key,nothing);before=get(base,key,nothing)
        current===before && continue
        (current===nothing || before===nothing || current!=before) ||
            continue
        diffs[key]=current
    end
    return diffs
end

function _parallel_replay_param_diffs!(m::GeoModel,diffs)
    for (curve,params) in diffs
        params===nothing ? delete!(m.curve_params,curve) :
            (m.curve_params[curve]=params)
    end
    return nothing
end

# Threaded detached production over one contiguous eligible run.  `range`
# shares the batch's serial merge position: `hn`/`hc` fix the allocator base,
# `base_params`/`base_dirty` fix the state the detached copies observed,
# `geom_points` fixes the geometry-birth overlay the produce can consult, and
# every slot records the produce's outcome for the in-order merge.  A run of
# one would pay the copy for no concurrency, so it falls back to the serial
# arm by leaving its slot empty.
function _parallel_surface_fill!(results::Vector{Any},range,m::GeoModel,
        prepared,input,options,assembly,dirty::Set{Tuple{Int,Int32}},
        precise::Bool,geometry_birth)
    length(range)<2 && return nothing
    hn=assembly.node_max;hc=assembly.cell_max
    base_params=_parallel_surface_param_refs(m)
    base_dirty=copy(dirty)
    geom_points=_higher_geometry_source_points(assembly,geometry_birth)
    Threads.@threads for offset in eachindex(range)
        index=range[offset]
        work=prepared[index]
        out=try
            body=()->begin
                working=Model._model_mesh_transaction_copy(m)
                produced=mesh_model_surface(working,work.tag;
                    size_field=nothing,_extrude_scope=nothing,
                    _mesh_input=input,_mesh1d_options=options,
                    _mesh_node_high_water=hn,_mesh_element_high_water=hc,
                    _return_surface_production=true,
                    _geometry_points=geom_points)
                return (produced,_parallel_curve_param_diffs(
                    working.curve_params,base_params,work.footprint))
            end
            (;ok=true,value=precise ?
                GmshLibm._with_win_cli_precision(body) : body())
        catch err
            err isa InterruptException && rethrow()
            (;ok=false,value=err)
        end
        results[index]=(out,hn,hc,base_dirty,geom_points)
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
function _parallel_surface_fresh(work::_ParallelSurfaceWork,out,
        dirty::Set{Tuple{Int,Int32}},base_dirty::Set{Tuple{Int,Int32}},
        assembly,geometry_birth,geom_points)
    out.ok || return false
    for entity in work.footprint
        entity in base_dirty && continue
        entity in dirty && return false
    end
    geom_points===nothing ||
        _higher_geometry_source_points(assembly,geometry_birth)==geom_points ||
        return false
    return _parallel_surface_remap_supported(out.value[1])
end

# Resolve a prepared slot into the unit's `products` list.  A fresh receipt
# replays its detached `m.curve_params` writeback onto the live model (the
# serial produce would have written it during the call) and remaps its fresh
# UIDs; anything else re-runs the identical serial produce so the merge sees
# exactly the state it always owned.
function _parallel_surface_products(m::GeoModel,work::_ParallelSurfaceWork,
        entry,dirty::Set{Tuple{Int,Int32}},assembly;
        input,scope,options,size_field,geometry_birth,caller)
    out,hn,hc,base_dirty,geom_points=entry
    if !_parallel_surface_fresh(work,out,dirty,base_dirty,assembly,
            geometry_birth,geom_points)
        return [(work.tag,mesh_model_surface(m,work.tag;
            size_field=size_field,_extrude_scope=scope,_mesh_input=input,
            _mesh1d_options=options,_mesh_node_high_water=assembly.node_max,
            _mesh_element_high_water=assembly.cell_max,
            _return_surface_production=true,
            _geometry_points=_higher_geometry_source_points(
                assembly,geometry_birth)))]
    end
    produced,diffs=out.value
    _parallel_replay_param_diffs!(m,diffs)
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
