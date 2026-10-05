# Native generation retains GeoPoint mesh vertices. Their public identity is
# carried by the Point owner, even when a prior attachment was consumed into
# the cache; coordinates alone never establish a source-node identity.
function _generation_point_sources(model,previous,caller)
    sources=Dict{Int32,Tuple{UInt64,NTuple{3,Float64}}}()
    if previous!==nothing
        old_class=_cached_classification_locked(previous)
        if old_class!==nothing && old_class.public_tags!==nothing
            public=old_class.public_tags
            for (node,owner) in enumerate(old_class.node_entities)
                owner[1]==0 && haskey(model.points,Int(owner[2])) || continue
                point=Tuple(previous.coords[:,node])
                source=(UInt64(_cache_node_tag(public,node)),point)
                existing=get!(sources,owner[2],source)
                existing==source || throw(ArgumentError(
                    "$caller: native Point[$(owner[2])] has several published mesh vertices"))
            end
        end
    end
    for ((dim,tag),record) in Model._discrete_mesh_records_model(model)
        dim==0 && haskey(model.points,tag) && !isempty(record.node_tags) || continue
        length(record.node_tags)==1 || throw(ArgumentError(
            "$caller: preserving several mesh vertices on native Point[$tag] is not implemented"))
        source=(UInt64(only(record.node_tags)),Tuple(record.node_coords[:,1]))
        existing=get!(sources,Int32(tag),source)
        existing==source || throw(ArgumentError(
            "$caller: native Point[$tag] has several published mesh vertices"))
    end
    return sources
end

# A retained Point record can own its cell while its mesh vertex is cache-owned.
# Read-only Point query products use a detached geometric view; the published
# native cache and its support rows are never converted or replaced.
function _public_point_query_mesh(model,cached,entity,caller;selected_tags=nothing)
    # Cache-only Point queries already have a geometric adapter and require no
    # owner classification. This detached view adds retained raw Point cells.
    have_raw_points=false
    for (_,_,record) in _discrete_mesh_records(model)
        if 15 in record.element_types
            have_raw_points=true;break
        end
    end
    have_raw_points || return nothing
    tags,nodes=_get_elements_by_type(15,entity)
    selected=selected_tags===nothing ? collect(eachindex(tags)) :
        [index for index in eachindex(tags) if tags[index] in selected_tags]
    node_positions=Dict{UInt64,Int32}()
    external_nodes=UInt64[];coordinates=NTuple{3,Float64}[]
    node_owners=Tuple{Int,Int32}[];cell_owners=Int32[]
    cells=zeros(Int32,1,length(selected));entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    for (column,index) in enumerate(selected)
        external=nodes[index]
        coordinate=_mesh_public_node_coords(model,cached,Int32(external))
        coordinate===nothing && throw(ArgumentError("$caller: Point node tag $external has no coordinates"))
        _,_,_,point=_get_element(tags[index])
        position=get(node_positions,external,Int32(0))
        if position==0
            position=Int32(length(external_nodes)+1)
            node_positions[external]=position
            push!(external_nodes,external);push!(coordinates,coordinate)
            push!(node_owners,(0,Int32(point)))
        end
        cells[1,column]=position;push!(cell_owners,Int32(point))
        get!(entities,(0,point),Elements.MixedEntity(0,point,coordinate))
    end
    isempty(selected) && return nothing
    data=Elements.MixedEntityData(entities;node_entities=node_owners,
        node_parametric=fill(nothing,length(external_nodes)),external_node_tags=external_nodes,
        block_entities=[cell_owners],external_element_tags=[tags[selected]])
    coords=zeros(3,length(coordinates))
    for (index,coordinate) in enumerate(coordinates),axis in 1:3
        coords[axis,index]=coordinate[axis]
    end
    return MixedMesh(coords,[ElementBlock(15,cells)];entity_data=data)
end

function _public_point_record_mesh(model,cached,element_tag,caller)
    tag=_mesh_query_integer(element_tag,caller,"element_tag")
    0<tag<=typemax(Int32) || throw(ArgumentError("$caller: element_tag must be a positive Int32 tag"))
    for (_,entity,record) in _discrete_mesh_records(model)
        for index in eachindex(record.element_tags)
            record.element_tags[index]==tag && record.element_types[index]==15 || continue
            return _public_point_query_mesh(model,cached,entity,caller;selected_tags=Set([UInt64(tag)]))
        end
    end
    return nothing
end

function _public_point_locations(model,cached,p,strict,caller)
    mesh=_public_point_query_mesh(model,cached,-1,caller)
    mesh===nothing && return nothing,Int[]
    tags=_locate_elements(_MixedMeshLocator(mesh),p,0,strict,caller)
    return mesh,sort!(Int.(tags);by=tag->_cache_element_tag(mesh,tag))
end

function _generation_bind_point_records!(model,cache,class,sources,caller)
    bindings=Dict{Int,UInt64}()
    for (point,(tag,coordinate)) in sources
        positions=class===nothing ? Int[] : findall(==((0,point)),class.node_entities)
        length(positions)<=1 || throw(ArgumentError(
            "$caller: native Point[$point] maps to several generated node rows"))
        record=_raw_model_mesh_record(model,0,Int(point))
        if isempty(positions)
            # A formerly cached Point can be outside the new higher-dimensional
            # closure. Keep its existing vertex in the model record instead.
            if record===nothing || isempty(record.node_tags)
                Model.add_discrete_nodes!(model,0,Int(point),[tag],collect(coordinate))
            end
            continue
        end
        position=only(positions)
        Tuple(cache.coords[:,position])==coordinate || throw(ArgumentError(
            "$caller: retaining edited native Point[$point] during boundary regeneration is not implemented"))
        bindings[position]=tag
        if record!==nothing
            for (msh,nodes) in zip(record.element_types,record.element_nodes)
                msh==15 && nodes==Int32[tag] || throw(ArgumentError(
                    "$caller: native Point[$point] has an unsupported mesh element identity"))
            end
            empty!(record.node_tags)
            record.node_coords=zeros(3,0)
            record.node_params=zeros(0,0)
            empty!(record.aux_params)
        end
    end
    return bindings
end

function _generation_public_labels!(model,previous,cache,class,overlay,mids,
        bindings,activate,renumber,caller)
    (activate || !isempty(bindings)) || return cache,class
    class!==nothing || throw(ArgumentError(
        "$caller: public generation labels require classified node ownership"))
    actual_count=overlay===nothing ? nnodes(cache) : nnodes(overlay)
    owners=overlay===nothing ? class.node_entities : vcat(class.node_entities,mids)
    node_entries=Tuple{Tuple{Int,Int32},Int,Int,UInt64}[]
    for position in 1:actual_count
        push!(node_entries,(owners[position],0,position,get(bindings,position,UInt64(0))))
    end
    records=Model._discrete_mesh_records_model(model)
    for ((dim,entity),record) in records,(position,tag) in enumerate(record.node_tags)
        push!(node_entries,((dim,Int32(entity)),1,position,UInt64(tag)))
    end
    sort!(node_entries;by=entry->(entry[1],entry[2],entry[3]),alg=Base.Sort.MergeSort)
    node_tags=zeros(UInt64,actual_count)
    node_map=Dict{UInt64,UInt64}()
    next_node=renumber ? UInt64(0) : max(UInt64(NODE_TAG_MAX[]),
        maximum((entry[4] for entry in node_entries);init=UInt64(0)))
    used_nodes=Set{UInt64}()
    for (_,kind,position,old_tag) in node_entries
        tag=old_tag
        if renumber || old_tag==0
            next_node<typemax(Int32) || throw(ArgumentError(
                "$caller: generated public node tags exceed Int32"))
            next_node+=1;tag=next_node
        end
        tag in used_nodes && throw(ArgumentError(
            "$caller: retained public node tag $tag is assigned more than once"))
        push!(used_nodes,tag)
        kind==0 && (node_tags[position]=tag)
        old_tag==0 || (node_map[old_tag]=tag)
    end
    _validate_generated_cache_references(model,previous,caller;
        retained_nodes=Set(values(bindings)))
    _tagged_record_renumber!(model,node_map,true)
    _,_,cell_count=_mesh_element_offsets(cache)
    cell_entries=Tuple{Tuple{Int,Int32},Int,Int,Int,UInt64}[]
    for (msh,dim,offset,cells,cell_owners) in _cache_catalog(cache,class)
        overlay!==nothing && msh==_p2_skeleton_etype(overlay) && (msh=_p2_etype(overlay))
        for column in axes(cells,2)
            push!(cell_entries,((dim,cell_owners[column]),Int(msh),0,offset+column,UInt64(0)))
        end
    end
    for ((dim,entity),record) in records,(position,tag) in enumerate(record.element_tags)
        push!(cell_entries,((dim,Int32(entity)),Int(record.element_types[position]),1,position,UInt64(tag)))
    end
    sort!(cell_entries;by=entry->(entry[1],entry[2],entry[3],entry[4]),alg=Base.Sort.MergeSort)
    element_tags=zeros(UInt64,cell_count)
    element_map=Dict{UInt64,UInt64}()
    next_element=renumber ? UInt64(0) : max(UInt64(ELEMENT_TAG_MAX[]),
        maximum((entry[5] for entry in cell_entries);init=UInt64(0)))
    used_elements=Set{UInt64}()
    for (_,_,kind,position,old_tag) in cell_entries
        tag=old_tag
        if renumber || old_tag==0
            next_element<typemax(Int32) || throw(ArgumentError(
                "$caller: generated public element tags exceed Int32"))
            next_element+=1;tag=next_element
        end
        tag in used_elements && throw(ArgumentError(
            "$caller: retained public element tag $tag is assigned more than once"))
        push!(used_elements,tag)
        kind==0 && (element_tags[position]=tag)
        old_tag==0 || (element_map[old_tag]=tag)
    end
    _tagged_record_renumber!(model,element_map,false)
    if cache isa MixedMesh && cache.entity_data!==nothing
        per_block=Vector{UInt64}[];offset=0
        for block in cache.blocks
            width=length(block.tags)
            push!(per_block,element_tags[offset+1:offset+width]);offset+=width
        end
        cache=_mixed_rebuild_metadata(cache,cache.blocks;external_node_tags=node_tags,
            external_element_tags=per_block,keep_empty_blocks=true)
        class=_mixed_rebind_class(class,cache)
    end
    public=_cache_public_tags(cache;node_tags,element_tags,node_count=actual_count,caller)
    return cache,_classification_with_public_tags(class,public)
end
