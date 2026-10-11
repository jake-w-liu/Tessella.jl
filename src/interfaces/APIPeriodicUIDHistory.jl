# The current-row mesh link certifies cell interpolation. The model record
# separately keeps native pair events whose vertex identities can be retired.
function _periodic_uid_history_source(m,mesh::MixedMesh,class,caller)
    histories=m.meshing.periodic_uid_history
    isempty(histories) && return mesh,class
    tags=_dim01_periodic_public_tags(mesh,class)
    tags===nothing && throw(ArgumentError(
        "$caller: periodic history requires the actual source UID table"))
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(tags))
    length(positions)==length(tags) || throw(ArgumentError(
        "$caller: periodic history source repeats a live node UID"))
    active=Set(key for key in keys(histories) if key[1]==2 && haskey(m.periodic,key))
    isempty(active) && return mesh,class
    links=Elements.MixedPeriodicLink[link for link in mesh.periodic_links
        if !((Int(link.dim),Int(link.slave_entity)) in active)]
    for key in sort!(collect(active))
        history=histories[key]
        constraint=m.periodic[key]
        constraint.master_entity==history.master && constraint.affine==history.affine ||
            throw(ArgumentError("$caller: periodic history changed its declared CAD relation"))
        isempty(history.primary) && continue
        slave=_periodic_surface_catalog(mesh,class,key[2],caller)
        master=_periodic_surface_catalog(mesh,class,history.master,caller)
        pairs=Pair{Int32,Int32}[]
        for event in history.primary
            event.slave_live && event.master_live || continue
            s=get(positions,event.slave,Int32(0));t=get(positions,event.master,Int32(0))
            s in slave.primary && t in master.primary || continue
            push!(pairs,s=>t)
        end
        # Cell admission remains a complete current primary bijection. Retired
        # or orphan events never manufacture a current interpolation vertex.
        roots=Dict(pairs)
        Set(keys(roots))==slave.primary && Set(values(roots))==master.primary &&
            length(roots)==length(pairs)==length(unique(values(roots))) ||
            throw(ArgumentError("$caller: periodic history lacks actual cell primary identities"))
        push!(links,Elements.MixedPeriodicLink(2,key[2],history.master,
            Int32[first(pair) for pair in pairs],Int32[last(pair) for pair in pairs];
            affine=history.affine))
    end
    rebuilt=_mixed_rebuild_metadata(mesh,mesh.blocks;periodic_links=links)
    return rebuilt,_mixed_rebind_class(class,rebuilt)
end

function _periodic_uid_history_live_uids(m,mesh,class)
    result=Set{UInt64}()
    for (_,_,record) in _discrete_mesh_records(m;include_authoritative=true)
        union!(result,record.node_tags)
    end
    if mesh!==nothing
        tags=_dim01_periodic_public_tags(mesh,class)
        tags===nothing || union!(result,tags)
    end
    return result
end

function _periodic_uid_history_clear!(staged,before_model,before_mesh,before_class,
        after_mesh,after_class,cleared)
    isempty(staged.meshing.periodic_uid_history) && return nothing
    before=_periodic_uid_history_live_uids(before_model,before_mesh,before_class)
    after=_periodic_uid_history_live_uids(staged,after_mesh,after_class)
    staged.meshing.periodic_uid_history=Model._model_periodic_uid_retire(
        staged.meshing.periodic_uid_history,setdiff(before,after))
    Model._model_periodic_uid_clear!(staged.meshing.periodic_uid_history,
        [(Int(dim),Int(tag)) for (dim,tag) in cleared])
    return nothing
end

function _periodic_uid_history_generate!(plan,bindings,caller)
    m=plan.model
    isempty(m.meshing.periodic_uid_history) && isempty(bindings) && return nothing
    replacement=Dict(key=>Model._model_periodic_uid_remap(history,plan.node_tag_map)
        for (key,history) in m.meshing.periodic_uid_history)
    tags=_dim01_periodic_public_tags(plan.mesh,plan.class)
    tags===nothing && throw(ArgumentError(
        "$caller: periodic generation publication requires actual target UIDs"))
    for (slave,(master,events)) in bindings
        constraint=get(m.periodic,(2,slave),nothing)
        constraint!==nothing && constraint.master_entity==master || throw(ArgumentError(
            "$caller: periodic generation changed its declared CAD master"))
        primary=Model._ModelPeriodicUIDPair[]
        for (s,t) in events
            haskey(plan.node_tag_map,s) && haskey(plan.node_tag_map,t) ||
                throw(ArgumentError("$caller: periodic copy loses an actual producer UID"))
            push!(primary,Model._ModelPeriodicUIDPair(
                plan.node_tag_map[s],plan.node_tag_map[t],true,true))
        end
        base=Set((event.slave,event.master) for event in primary)
        higher=Model._ModelPeriodicUIDPair[]
        links=[link for link in plan.mesh.periodic_links if
            link.dim==2 && link.slave_entity==slave]
        length(links)==1 || throw(ArgumentError(
            "$caller: periodic generation requires one completed current Surface link"))
        link=only(links)
        link.master_entity==master && link.affine==constraint.affine ||
            throw(ArgumentError("$caller: periodic generation link changed its CAD relation"))
        for (s,t) in zip(link.slave_nodes,link.master_nodes)
            pair=(tags[s],tags[t])
            pair in base && continue
            push!(higher,Model._ModelPeriodicUIDPair(pair[1],pair[2],true,true))
        end
        replacement[(2,slave)]=Model._ModelPeriodicUIDHistory(Int32(master),
            constraint.affine,primary,higher)
    end
    m.meshing.periodic_uid_history=replacement
    return nothing
end

function _periodic_uid_history_order!(m,mesh::MixedMesh,class,
        mapping::AbstractDict{UInt64,UInt64},caller)
    isempty(m.meshing.periodic_uid_history) && return nothing
    tags=_dim01_periodic_public_tags(mesh,class)
    tags===nothing && throw(ArgumentError(
        "$caller: periodic order publication requires actual target UIDs"))
    replacement=Dict{Tuple{Int,Int},Model._ModelPeriodicUIDHistory}()
    for (key,history) in m.meshing.periodic_uid_history
        higher=Model._ModelPeriodicUIDPair[]
        matching=[link for link in mesh.periodic_links if
            link.dim==key[1] && link.slave_entity==key[2]]
        length(matching)<=1 || throw(ArgumentError(
            "$caller: periodic order repeats a Surface link"))
        if !isempty(matching)
            link=only(matching)
            link.master_entity==history.master && link.affine==history.affine ||
                throw(ArgumentError("$caller: actual periodic target changed its CAD relation"))
            slave=_periodic_surface_catalog(mesh,class,key[2],caller)
            master=_periodic_surface_catalog(mesh,class,history.master,caller)
            for (s,t) in zip(link.slave_nodes,link.master_nodes)
                (s in slave.primary && t in master.primary) && continue
                haskey(slave.stencils,s) && haskey(master.stencils,t) ||
                    throw(ArgumentError("$caller: periodic support lacks an actual target cell"))
                push!(higher,Model._ModelPeriodicUIDPair(tags[s],tags[t],true,true))
            end
        end
        replacement[key]=Model._model_periodic_uid_set_order(history,mapping,higher)
    end
    m.meshing.periodic_uid_history=replacement
    return nothing
end

function _periodic_uid_history_refine!(m,source,mesh::MixedMesh,class,
        mapping::AbstractDict{UInt64,UInt64},element_mapping,receipt,caller)
    isempty(m.meshing.periodic_uid_history) && return nothing
    if receipt===nothing
        source isa MixedMesh &&
            all(link->link.dim!=2 || isempty(link.slave_nodes),source.periodic_links) &&
            all(link->link.dim!=2 || isempty(link.slave_nodes),mesh.periodic_links) &&
            all(history->isempty(history.higher),values(m.meshing.periodic_uid_history)) ||
            throw(ArgumentError("$caller: periodic refinement lacks its actual quadratic parent"))
        # An empty declaration stays empty through independent refinement.
        # No current Face correspondence means there is no periodic child
        # interpolation to certify or manufacture from the CAD declaration.
        replacement=Dict(key=>Model._model_periodic_uid_refine(history,mapping,
            Model._ModelPeriodicUIDPair[]) for (key,history) in m.meshing.periodic_uid_history)
        m.meshing.periodic_uid_history=replacement
        return nothing
    end
    receipt isa _SurfacePeriodicRefineReceipt || throw(ArgumentError(
        "$caller: periodic refinement publication requires its actual quadratic parent"))
    _surface_periodic_refine_validate(receipt,source,mesh,class,mapping,element_mapping,caller)
    source_tags=_dim01_periodic_public_tags(source,receipt.source_class)
    source_tags===nothing && throw(ArgumentError(
        "$caller: periodic quadratic parent has no actual UID table"))
    tags=_dim01_periodic_public_tags(mesh,class)
    tags===nothing && throw(ArgumentError(
        "$caller: periodic refinement publication requires actual child UIDs"))
    replacement=Dict{Tuple{Int,Int},Model._ModelPeriodicUIDHistory}()
    for (key,history) in m.meshing.periodic_uid_history
        # Refinement performs its own fresh SetOrder2. Its current supports,
        # rather than a prior public order operation's supports, are promoted.
        expected=Tuple{UInt64,UInt64}[]
        source_links=[link for link in source.periodic_links if
            link.dim==key[1] && link.slave_entity==key[2]]
        length(source_links)<=1 || throw(ArgumentError(
            "$caller: periodic quadratic parent repeats a Surface link"))
        if !isempty(source_links)
            parent=only(source_links)
            slave=_periodic_surface_catalog(source,receipt.source_class,key[2],caller)
            master=_periodic_surface_catalog(source,receipt.source_class,history.master,caller)
            for (s,t) in zip(parent.slave_nodes,parent.master_nodes)
                (s in slave.primary && t in master.primary) && continue
                push!(expected,(source_tags[s],source_tags[t]))
            end
        end
        all(event->event.slave_live && event.master_live,history.higher) &&
            sort([(event.slave,event.master) for event in history.higher])==sort(expected) ||
            throw(ArgumentError("$caller: periodic history is not bound to the fresh quadratic parent"))
        children=Model._ModelPeriodicUIDPair[]
        matching=[link for link in mesh.periodic_links if
            link.dim==key[1] && link.slave_entity==key[2]]
        length(matching)<=1 || throw(ArgumentError(
            "$caller: periodic refinement repeats a Surface link"))
        if !isempty(matching)
            link=only(matching)
            link.master_entity==history.master && link.affine==history.affine ||
                throw(ArgumentError("$caller: periodic children changed their declared CAD relation"))
            append!(children,[Model._ModelPeriodicUIDPair(tags[s],tags[t],true,true)
                for (s,t) in zip(link.slave_nodes,link.master_nodes)])
        end
        replacement[key]=Model._model_periodic_uid_refine(history,mapping,children)
    end
    m.meshing.periodic_uid_history=replacement
    return nothing
end
