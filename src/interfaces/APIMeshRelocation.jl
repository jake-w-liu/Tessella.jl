# Relocation uses retained native parameters. Current Cartesian coordinates
# cannot recover the original parameter after an explicit mesh edit.
function _relocation_legacy_parameters(model,mesh,owners,stored)
    parameters=deepcopy(stored)
    for index in eachindex(owners)
        dimension,tag=owners[index]
        dimension in (1,2) && parameters[index]===nothing || continue
        _relocation_native_entity(model,dimension,Int(tag)) || continue
        parameters[index]=_api_model_parametrization(model,dimension,Int(tag),
            collect(@view mesh.coords[:,index]))
    end
    return parameters
end

function _relocation_rebind_parameters(class,mesh,primary_map)
    primary_map===nothing && return class.node_parametric
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(mesh)]
    assigned=falses(nnodes(mesh))
    for old in eachindex(primary_map)
        target=primary_map[old]
        target==0 && continue
        assigned[target] && continue
        parameters[target]=class.node_parametric[old]
        assigned[target]=true
    end
    return parameters
end

function _relocation_rebind_mid_parameters(class,primary_map)
    primary_map===nothing && return class.mid_parametric
    parameters=Dict{NTuple{2,Int32},Union{Nothing,Vector{Float64}}}()
    overlay=_high_order_overlay(class.mesh)
    source_keys=overlay===nothing ? sort!(collect(keys(class.mid_parametric))) :
        [minmax(_p2_cells(overlay)[i,column],_p2_cells(overlay)[j,column])
         for column in axes(_p2_cells(overlay),2) for (_,i,j) in _p2_edge_slots(overlay)]
    for key in source_keys
        value=get(class.mid_parametric,key,nothing)
        a,b=primary_map[key[1]],primary_map[key[2]]
        (a==0 || b==0 || a==b) && continue
        mapped=minmax(Int32(a),Int32(b))
        get!(parameters,mapped,value)
    end
    return parameters
end

function _relocation_mid_key(class,node)
    key=get(class.mid_keys,Int32(node),nothing)
    key===nothing && throw(ArgumentError(
        "API mesh relocation: midpoint has no primary-edge support"))
    return key
end

function _relocation_prepare_mid_parameters(model,cache,class,overlay,owners)
    cells=_p2_cells(overlay)
    count=nnodes(cache)
    staged=copy(class.mid_parametric)
    keys=Dict{Int32,NTuple{2,Int32}}()
    for column in axes(cells,2),(slot,i,j) in _p2_edge_slots(overlay)
        key=minmax(cells[i,column],cells[j,column])
        node=cells[slot,column]
        keys[node]=key
        haskey(staged,key) && continue
        dimension,tag=owners[node-count]
        parameter=dimension in (1,2) && _relocation_native_entity(model,dimension,Int(tag)) ?
            _api_model_parametrization(model,dimension,Int(tag),collect(@view overlay.coords[:,node])) : nothing
        staged[key]=parameter
    end
    return staged,keys
end

function _relocation_publish_mid_parameters!(class,prepared)
    empty!(class.mid_parametric);merge!(class.mid_parametric,prepared[1])
    empty!(class.mid_keys);merge!(class.mid_keys,prepared[2])
    return nothing
end

function _relocation_overlay_coordinates(model,cache,class,overlay,owners,dimension,entity,caller)
    overlay===nothing && return nothing
    class===nothing && return nothing
    edited=deepcopy(overlay)
    changed=false
    visited=falses(nnodes(overlay)-nnodes(cache))
    cells=_p2_cells(overlay)
    for column in axes(cells,2),(slot,i,j) in _p2_edge_slots(overlay)
        node=cells[slot,column]
        index=node-nnodes(cache)
        visited[index] && continue
        visited[index]=true
        owner=owners[index]
        _relocation_selected(owner,dimension,entity) || continue
        key=minmax(cells[i,column],cells[j,column])
        point=_relocation_value(model,owner,get(class.mid_parametric,key,nothing),caller)
        point===nothing && continue
        for axis in 1:3
            changed |= edited.coords[axis,node]!=point[axis]
            edited.coords[axis,node]=point[axis]
        end
    end
    return changed ? edited : nothing
end

@inline _relocation_selected(owner,dimension,entity)=
    (dimension<0 || owner[1]==dimension) &&
    (dimension<0 || entity<0 || owner[2]==entity)

function _relocation_native_entity(model,dimension,tag)
    return dimension==0 ? haskey(model.points,tag) :
           dimension==1 ? haskey(model.curves,tag) :
           dimension==2 ? haskey(model.surfaces,tag) : false
end

function _relocation_value(model,owner,parameters,caller)
    dimension,tag=owner
    _relocation_native_entity(model,dimension,Int(tag)) || return nothing
    if dimension>0
        parameters===nothing && return nothing
        length(parameters)>=dimension || return nothing
    end
    coordinates=model_value(model,dimension,Int(tag),
        dimension==0 ? Float64[] : parameters[1:dimension];
        old_ruled_surface=!iszero(OPTIONS["Geometry.OldRuledSurface"]))
    length(coordinates)==3 && all(isfinite,coordinates) || throw(ArgumentError(
        "$caller: entity evaluation must produce three finite coordinates"))
    return coordinates
end

function _relocation_record_coordinates(model,record,owner,caller)
    _relocation_native_entity(model,owner[1],owner[2]) || return nothing
    isempty(record.node_tags) && return nothing
    count=length(record.node_tags)
    if owner[1]==0
        point=_relocation_value(model,owner,nothing,caller)
        coordinates=Matrix{Float64}(undef,3,count)
        for column in 1:count,axis in 1:3
            coordinates[axis,column]=point[axis]
        end
    else
        selected=Int[];parameters=Float64[]
        for position in eachindex(record.node_tags)
            stored=Model._model_mesh_node_parameters(record,position)
            stored!==nothing && length(stored)>=owner[1] || continue
            push!(selected,position);append!(parameters,@view stored[1:owner[1]])
        end
        isempty(selected) && return nothing
        values=model_value(model,owner[1],owner[2],parameters;
            old_ruled_surface=!iszero(OPTIONS["Geometry.OldRuledSurface"]))
        length(values)==3length(selected) && all(isfinite,values) || throw(ArgumentError(
            "$caller: entity evaluation must produce three finite coordinates per node"))
        coordinates=copy(record.node_coords)
        coordinates[:,selected]=reshape(values,3,length(selected))
    end
    return coordinates==record.node_coords ? nothing : coordinates
end

function _relocation_cache_coordinates(model,cached,class,dimension,entity,caller)
    class===nothing && return nothing
    parameters=class.node_parametric
    coordinates=copy(cached.coords)
    changed=false
    for column in eachindex(class.node_entities)
        owner=class.node_entities[column]
        _relocation_selected(owner,dimension,entity) || continue
        point=_relocation_value(model,owner,parameters[column],caller)
        point===nothing && continue
        for axis in 1:3
            changed |= coordinates[axis,column]!=point[axis]
            coordinates[axis,column]=point[axis]
        end
    end
    return changed ? coordinates : nothing
end

function _relocate_nodes(dim=-1,tag=-1)
    caller="API.mesh.relocate_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        dimension>=0 && entity>=0 && !_model_entity_known(model,dimension,entity) &&
            throw(ArgumentError("$caller: unknown entity ($dimension, $entity)"))
        cached=LAST_MESH[]
        class=cached===nothing ? nothing : _cached_classification_locked(cached)
        record_updates=Tuple{Tuple{Int,Int},Vector{UInt64},Matrix{Float64}}[]
        for (edim,etag,record) in _discrete_mesh_records(model)
            owner=(edim,etag)
            _relocation_selected(owner,dimension,entity) || continue
            coordinates=_relocation_record_coordinates(model,record,owner,caller)
            coordinates===nothing || push!(record_updates,(owner,copy(record.node_tags),coordinates))
        end
        coordinates=cached===nothing ? nothing :
            _relocation_cache_coordinates(model,cached,class,dimension,entity,caller)
        replacement=cached
        replacement_class=class
        overlay=cached===nothing ? nothing : _high_order_overlay(cached)
        edited_overlay=overlay===nothing ? nothing : _relocation_overlay_coordinates(
            model,cached,class,overlay,LAST_MESH_HIGH_ORDER_MIDS[],dimension,entity,caller)
        prepared=overlay===nothing ? nothing :
            (edited_overlay===nothing ? deepcopy(overlay) : edited_overlay,
             copy(LAST_MESH_HIGH_ORDER_MIDS[]))
        if coordinates!==nothing
            replacement=cached isa MixedMesh ?
                _mixed_rebuild_metadata(cached,cached.blocks;coords=coordinates) :
                Mesh(coordinates;segs=cached.segs,tris=cached.tris,tets=cached.tets,
                     seg_tag=cached.seg_tag,tri_tag=cached.tri_tag,tet_tag=cached.tet_tag)
            replacement_class=cached isa MixedMesh ?
                _mixed_rebind_class(class,replacement) : _api_rebind_class(class,replacement)
            if prepared!==nothing
                prepared[1].coords[:,1:nnodes(cached)]=coordinates
            end
        end
        # All CAD evaluation and allocation finishes before publication.
        if replacement!==cached || edited_overlay!==nothing
            _replace_mesh_cache_with_overlay_locked!(replacement,replacement_class,caller;
                prepared_overlay=prepared,preserve_visibility=true)
        end
        for (owner,tags,new_coordinates) in record_updates
            # Cache publication can replace mirrored records. Resolve the live
            # record by owner and UID after publication, never mutate a filtered
            # read view or retain a pointer to a replaced record.
            record=_raw_model_mesh_record(_model_locked(),owner...)
            positions=Dict(uid=>position for (position,uid) in enumerate(record.node_tags))
            for (column,uid) in enumerate(tags)
                record.node_coords[:,positions[uid]]=new_coordinates[:,column]
            end
        end
        return nothing
    end
end
