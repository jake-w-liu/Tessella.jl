# Detached public-tag mutations. Dense coordinate rows and block cell columns
# keep their identity; only the public labels and mirrored record references
# change. A failed request never changes the live model or cache.
function _tagged_renumber_lists(old_tags,new_tags,caller,kind)
    (old_tags isa AbstractVector || old_tags isa Tuple) &&
        (new_tags isa AbstractVector || new_tags isa Tuple) || throw(ArgumentError(
            "$caller: old_tags and new_tags must be vectors or tuples"))
    length(old_tags)==length(new_tags) || throw(ArgumentError(
        "$caller: old_tags and new_tags must have the same length"))
    return map((old_tags,new_tags),("old_tags","new_tags")) do values,name
        parsed=UInt64[]
        for value in values
            tag=_cache_checked_public_tag(value,caller,kind)
            tag<=typemax(Int32) || throw(ArgumentError(
                "$caller: $name must contain positive Int32 tags"))
            push!(parsed,tag)
        end
        parsed
    end
end

function _tagged_renumber_order(model,cached,class,node::Bool,caller)
    public=class.public_tags
    entries=Tuple{Tuple{Int,Int32},Int,Int,UInt64}[]
    count=node ? length(public.node_tags) : length(public.element_tags)
    for ((dim,tag),record) in Model._discrete_mesh_records_model(model)
        _cache_covers_record(public,dim,tag) && continue
        count=Base.checked_add(count,length(node ? record.node_tags : record.element_tags))
    end
    sizehint!(entries,count)
    if node
        for (position,tag) in enumerate(public.node_tags)
            owner=position<=length(class.node_entities) ? class.node_entities[position] :
                LAST_MESH_HIGH_ORDER_MIDS[][position-length(class.node_entities)]
            push!(entries,(owner,0,position,tag))
        end
    else
        for (msh,dim,offset,_,owners) in _cache_catalog(cached,class)
            overlay=_high_order_overlay(cached)
            overlay!==nothing && msh==_p2_skeleton_etype(overlay) && (msh=_p2_etype(overlay))
            for (column,owner) in enumerate(owners)
                push!(entries,((dim,owner),Int(msh),column,
                               public.element_tags[offset+column]))
            end
        end
    end
    seen=Set(last(entry) for entry in entries)
    for ((dim,tag),record) in Model._discrete_mesh_records_model(model)
        _cache_covers_record(public,dim,tag) && continue
        tags=node ? record.node_tags : record.element_tags
        for (position,value) in enumerate(tags)
            external=UInt64(value)
            external in seen && throw(ArgumentError(
                "$caller: public $(node ? "node" : "element") tag $external " *
                "is assigned on more than one entity"))
            push!(seen,external)
            msh=node ? 0 : Int(record.element_types[position])
            push!(entries,((dim,Int32(tag)),msh,position,external))
        end
    end
    sort!(entries;by=entry->(entry[1],entry[2],entry[3]),alg=Base.Sort.MergeSort)
    return UInt64[last(entry) for entry in entries]
end

function _tagged_renumber_mapping(order,old_tags,new_tags,caller,kind)
    requested=Dict{UInt64,UInt64}()
    for (old,new) in zip(old_tags,new_tags)
        requested[old]=new # Gmsh uses the last pair when an old tag is repeated.
    end
    # Unmatched live tags are allocated after every requested new label,
    # including a label paired with an unknown old tag.
    next=maximum(values(requested);init=UInt64(0))
    mapping=Dict{UInt64,UInt64}()
    used=Set{UInt64}()
    sizehint!(mapping,length(order));sizehint!(used,length(order))
    for old in order
        new=get(requested,old,UInt64(0))
        if new==0
            next<typemax(Int32) || throw(ArgumentError(
                "$caller: renumbered $kind tags exceed Int32"))
            next+=UInt64(1)
            new=next
        end
        new in used && throw(ArgumentError(
            "$caller: renumbering produces duplicate $kind tag $new"))
        push!(used,new)
        mapping[old]=new
    end
    return mapping
end

function _tagged_record_renumber!(model,mapping,node::Bool)
    for (_,record) in Model._discrete_mesh_records_model(model)
        tags=node ? record.node_tags : record.element_tags
        for position in eachindex(tags)
            tags[position]=Int32(get(mapping,UInt64(tags[position]),UInt64(tags[position])))
        end
        node || continue
        for connectivity in record.element_nodes,position in eachindex(connectivity)
            connectivity[position]=Int32(get(mapping,UInt64(connectivity[position]),
                                            UInt64(connectivity[position])))
        end
        record.aux_params=Dict{Int32,Vector{Float64}}(
            Int32(get(mapping,UInt64(tag),UInt64(tag)))=>parameters
            for (tag,parameters) in record.aux_params)
    end
    return nothing
end

function _tagged_mesh_max_tags(model,mesh,public=nothing)
    public===nothing && (public=_cache_public_tags(mesh))
    node=maximum(public.node_tags;init=UInt64(0))
    element=maximum(public.element_tags;init=UInt64(0))
    for (_,record) in Model._discrete_mesh_records_model(model)
        node=max(node,UInt64(maximum(record.node_tags;init=Int32(0))))
        element=max(element,UInt64(maximum(record.element_tags;init=Int32(0))))
    end
    return node,element
end

function _tagged_renumber_plan(model,cached,class,old_list,new_list,kind,caller)
    node=kind===:node || kind=="node"
    node || kind===:element || kind=="element" || throw(ArgumentError(
        "$caller: renumber kind must be node or element"))
    cached isa Union{Mesh,MixedMesh} && class!==nothing && class.mesh===cached ||
        throw(ArgumentError("$caller: renumbering requires a classified mesh cache"))
    overlay=_high_order_overlay(cached)
    public=class.public_tags
    if public===nothing
        count=overlay===nothing ? nnodes(cached) : nnodes(overlay)
        public=_cache_public_tags(cached;node_count=count)
        class=_classification_with_public_tags(class,public)
    end
    label=node ? "node" : "element"
    old_tags,new_tags=_tagged_renumber_lists(old_list,new_list,caller,label)
    order=_tagged_renumber_order(model,cached,class,node,caller)
    all(tag->0<tag<=typemax(Int32),order) || throw(ArgumentError(
        "$caller: retained $label tags exceed Int32"))
    mapping=_tagged_renumber_mapping(order,old_tags,new_tags,caller,label)
    staged_model=deepcopy(model)
    _tagged_record_renumber!(staged_model,mapping,node)
    node_tags=node ? UInt64[mapping[tag] for tag in public.node_tags] : public.node_tags
    element_tags=node ? public.element_tags : UInt64[mapping[tag] for tag in public.element_tags]
    mesh=cached
    if cached isa MixedMesh && cached.entity_data!==nothing
        # Labels do not validate or reconstruct an existing mesh. In particular,
        # an edited empty block keeps its stored position and metadata arrays.
        mesh=deepcopy(cached)
        data=mesh.entity_data
        data.external_node_tags[:]=node_tags
        offset=0
        for (index,block) in enumerate(mesh.blocks)
            width=length(block.tags)
            data.external_element_tags[index][:]=element_tags[offset+1:offset+width]
            offset+=width
        end
    end
    table=_cache_public_tags(mesh;authority=public.authority,node_tags,element_tags,
        node_count=length(node_tags),caller)
    rebound=mesh===cached ? class : _mixed_rebind_class(class,mesh)
    new_class=_classification_with_public_tags(rebound,table)
    max_node_tag,max_element_tag=_tagged_mesh_max_tags(staged_model,mesh,table)
    # Renumbering does not lower Gmsh's historical maximum tag counters.
    previous_node,previous_element=_tagged_mesh_max_tags(model,cached,public)
    return (model=staged_model,mesh=mesh,class=new_class,
        element_tag_map=node ? nothing : mapping,authority=public.authority,
        preserved_overlay=overlay,preserved_mid_owners=LAST_MESH_HIGH_ORDER_MIDS[],
        max_node_tag=max(max_node_tag,previous_node,UInt64(NODE_TAG_MAX[])),
        max_element_tag=max(max_element_tag,previous_element,UInt64(ELEMENT_TAG_MAX[])))
end

function _mixed_remap_periodic_nodes(mesh,remap,caller;drop_removed::Bool=false)
    links=Elements.MixedPeriodicLink[]
    for link in mesh.periodic_links
        slaves=Int32[];masters=Int32[]
        by_slave=Dict{Int32,Int32}();by_master=Dict{Int32,Int32}()
        for (old_slave,old_master) in zip(link.slave_nodes,link.master_nodes)
            slave=remap[old_slave];master=remap[old_master]
            if slave==0 || master==0
                drop_removed && continue
                throw(ArgumentError(
                    "$caller: a periodic correspondence references a removed node"))
            end
            existing=get(by_slave,slave,Int32(0))
            existing==master && continue
            (existing==0 && !haskey(by_master,master)) || throw(ArgumentError(
                "$caller: merged periodic nodes have conflicting correspondences"))
            by_slave[slave]=master;by_master[master]=slave
            push!(slaves,slave);push!(masters,master)
        end
        push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,
            link.master_entity,slaves,masters;affine=link.affine))
    end
    return links
end

function _dim01_native_interpolation(cached,class,caller)
    interpolation=Dict{Tuple,Int32}()
    orders=Set{Int}()
    top=Set(class.entities)
    for (msh,dim,_,cells,owners) in _cache_catalog(cached,class)
        dim>=2 && any(owner->(dim,owner) in top,owners) || continue
        spec=msh_spec(msh)
        spec.order in (1,2) && !spec.serendipity || throw(ArgumentError(
            "$caller: native source completion supports complete linear and quadratic cells; " *
            "MSH $msh is unsupported"))
        push!(orders,spec.order)
        spec.order==1 && continue
        family=spec.family
        primary=Elements.msh_num_nodes(msh_type(family,1))
        reference=Elements.lagrange_nodes(msh)
        for slot in primary+1:size(cells,1)
            weights=MeshFunctionSpaces._first_order_values(Val(family),
                reference[1,slot],reference[2,slot],reference[3,slot],caller,slot)
            for column in axes(cells,2)
                (dim,owners[column]) in top || continue
                key=_dim01_support_key(@view(cells[1:primary,column]),weights)
                node=cells[slot,column]
                previous=get(interpolation,key,node)
                previous==node || throw(ArgumentError(
                    "$caller: native quadratic source has distinct interpolation " *
                    "nodes on the same shared reference entity"))
                interpolation[key]=node
            end
        end
    end
    length(orders)<=1 || throw(ArgumentError(
        "$caller: mixed linear/quadratic native source completion is unsupported"))
    return interpolation,2 in orders
end

function _dim01_lift_native_cell(primary_nodes,msh,quadratic,interpolation,caller)
    (!quadratic || msh==15) && return msh,Vector{Int32}(primary_nodes)
    family=msh_spec(msh).family
    target=msh_type(family,2)
    reference=Elements.lagrange_nodes(target)
    primary=length(primary_nodes)
    nodes=Vector{Int32}(undef,size(reference,2))
    nodes[1:primary]=primary_nodes
    for slot in primary+1:length(nodes)
        weights=MeshFunctionSpaces._first_order_values(Val(family),
            reference[1,slot],reference[2,slot],reference[3,slot],caller,slot)
        key=_dim01_support_key(primary_nodes,weights)
        node=get(interpolation,key,Int32(0))
        node!=0 || throw(ArgumentError(
            "$caller: native source lacks the actual MSH $target interpolation " *
            "node required by a projected boundary cell"))
        nodes[slot]=node
    end
    return target,nodes
end

function _dim01_native_source_model(m,cached,class,caller)
    source=deepcopy(m)
    primary=falses(nnodes(cached))
    edges=Set{NTuple{2,Int32}}()
    top=Set(class.entities)
    for (msh,dim,_,cells,owners) in _cache_catalog(cached,class)
        dim>=2 || continue
        linear=msh_type(msh_spec(msh).family,1)
        count=Elements.msh_num_nodes(linear)
        patterns=MeshEntityTopology._simplex_edge_patterns(linear)
        for column in axes(cells,2)
            (dim,owners[column]) in top || continue
            for row in 1:count;primary[cells[row,column]]=true;end
            for pattern in patterns
                a=cells[pattern[1],column];b=cells[pattern[2],column]
                push!(edges,a<b ? (a,b) : (b,a))
            end
        end
    end
    for ((dim,curve),_) in class.boundaries
        dim==1 && haskey(source.curves,Int(curve)) || continue
        a,b=source.curves[curve]
        eligible=BitVector([primary[node] &&
            (owner==(1,curve) || owner==(0,Int32(a)) || owner==(0,Int32(b)))
            for (node,owner) in enumerate(class.node_entities)])
        any(eligible) || continue
        entries,_=Model._curve_parameter_nodes(source,cached,Int(curve),
                                               eligible,edges,1e-12,caller)
        params=Float64[entry[1] for entry in entries]
        source.curve_params[Int(curve)]=params
        attribute=get(source.meshing.transfinite_curves,Int(curve),nothing)
        attribute===nothing || (source.meshing.transfinite_curves[Int(curve)]=
            merge(attribute,(num_nodes=length(params),)))
    end
    return source
end

function _dim01_native_projection_part(m,cached,class,entity,caller)
    dim,tag=entity
    catalog=_cache_catalog(cached,class)
    selected=[entry for entry in catalog if entry[2]==dim && any(==(tag),entry[5])]
    isempty(selected) && throw(ArgumentError(
        "$caller: native entity $entity has no source cells"))
    referenced=falses(nnodes(cached))
    for (msh,_,_,cells,owners) in selected
        primary=Elements.msh_num_nodes(msh_type(msh_spec(msh).family,1))
        for column in axes(cells,2)
            owners[column]==tag || continue
            for row in 1:primary
                referenced[cells[row,column]]=true
            end
        end
    end
    original=findall(referenced)
    remap=zeros(Int32,nnodes(cached))
    remap[original]=Int32.(1:length(original))
    blocks=ElementBlock[]
    for (msh,_,_,cells,owners) in selected
        linear=msh_type(msh_spec(msh).family,1)
        primary=Elements.msh_num_nodes(linear)
        columns=findall(==(tag),owners)
        push!(blocks,ElementBlock(linear,remap[cells[1:primary,columns]],
                                  zeros(Int32,length(columns))))
    end
    input=MixedMesh(cached.coords[:,original],blocks)
    projected=try
        model_to_mixed(m,input,dim,Int(tag))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: native source completion for $entity " *
                            "cannot certify its exact existing mesh — " *sprint(showerror,err)))
    end
    return projected,original
end

# Legacy native generation stores only its primary cell families. Complete a
# detached input for the 1-D planner through the certified projection of those
# exact cells. No model entity is remeshed, no coordinate is generated, and no
# unrelated coincident vertices are welded while combining entity projections.
function _dim01_complete_native_cache(m,cached,class,caller)
    (cached===nothing || class===nothing || class.public_tags!==nothing) &&
        return cached,class
    entities=sort!([entity for entity in class.entities if entity[1] in (2,3)])
    isempty(entities) && return cached,class
    any(record->!isempty(record.node_tags) || !isempty(record.element_tags),
        values(m.meshing.attached)) && throw(ArgumentError(
            "$caller: a legacy native higher-dimensional cache with attached " *
            "mesh records has untracked source tag allocations; completing " *
            "that allocation history is not supported"))
    interpolation,quadratic=_dim01_native_interpolation(cached,class,caller)
    source_model=_dim01_native_source_model(m,cached,class,caller)
    metadata=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    boundaries=copy(class.boundaries)
    node_entities=copy(class.node_entities)
    cells=Tuple{Int,Tuple{Int,Int32},Vector{Int32},Int32}[]
    seen=Set{Tuple{Int,Tuple{Int,Int32},Tuple}}()
    links=Elements.MixedPeriodicLink[]
    link_keys=Set{Tuple}()
    function append_cell!(msh,owner,nodes,physical)
        primary=msh==15 ? 1 : Elements.msh_num_nodes(msh_type(msh_spec(msh).family,1))
        key=(Int(msh),owner,Tuple(sort!(collect(@view nodes[1:primary]))))
        key in seen && return nothing
        push!(seen,key)
        push!(cells,(Int(msh),owner,Vector{Int32}(nodes),Int32(physical)))
        for node in nodes
            node_entities[node]=_merged_node_owner(node_entities[node],owner)
        end
        return nothing
    end
    for entity in entities
        projected,original=_dim01_native_projection_part(source_model,cached,class,entity,caller)
        data=projected.entity_data
        for (key,record) in data.entities
            previous=get(metadata,key,nothing)
            if previous!==nothing
                previous.boundaries==record.boundaries &&
                    previous.physical_tags==record.physical_tags || throw(ArgumentError(
                        "$caller: native projections disagree about shared entity $key"))
            else
                metadata[key]=record
            end
            boundaries[(key[1],Int32(key[2]))]=Int32[abs(tag) for tag in record.boundaries]
        end
        for (bi,block) in enumerate(projected.blocks),column in axes(block.nodes,2)
            owner=(msh_dimension(block.msh),data.block_entities[bi][column])
            primary_nodes=Int32.(original[block.nodes[:,column]])
            msh,nodes=_dim01_lift_native_cell(primary_nodes,Int(block.msh),
                                             quadratic,interpolation,caller)
            append_cell!(msh,owner,nodes,block.tags[column])
        end
        for link in projected.periodic_links
            slaves=Int32.(original[link.slave_nodes]);masters=Int32.(original[link.master_nodes])
            key=(link.dim,link.slave_entity,link.master_entity,Tuple(slaves),Tuple(masters))
            key in link_keys && continue
            push!(link_keys,key)
            push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,
                link.master_entity,slaves,masters;affine=link.affine))
        end
    end
    # Preserve existing lower-dimensional native cells that are not a projected
    # boundary duplicate. They already carry classified dense source identity.
    for (bi,(msh,dim,_,nodes,owners)) in enumerate(_cache_catalog(cached,class))
        dim<2 || continue
        physical=_cache_native_blocks(cached)[bi][3]
        for column in axes(nodes,2)
            spec=msh_spec(msh)
            primary=msh==15 ? 1 : Elements.msh_num_nodes(msh_type(spec.family,1))
            target,actual=_dim01_lift_native_cell(nodes[1:primary,column],
                msh==15 ? 15 : msh_type(spec.family,1),quadratic,interpolation,caller)
            append_cell!(target,(dim,owners[column]),actual,physical[column])
        end
    end
    sort!(cells;by=cell->(cell[2],cell[1]),alg=Base.Sort.MergeSort)
    families=Int[];positions=Dict{Int,Vector{Int}}()
    for (position,cell) in enumerate(cells)
        haskey(positions,cell[1]) || push!(families,cell[1])
        push!(get!(()->Int[],positions,cell[1]),position)
    end
    blocks=ElementBlock[];block_owners=Vector{Int32}[];external=Vector{UInt64}[]
    owners=Dict{Int,Vector{Int32}}()
    for msh in families
        selected=positions[msh]
        nodes=hcat((cells[position][3] for position in selected)...)
        physical=Int32[cells[position][4] for position in selected]
        classified=Int32[cells[position][2][2] for position in selected]
        push!(blocks,ElementBlock(msh,nodes,physical))
        push!(block_owners,classified);push!(external,UInt64.(selected))
        owners[msh]=classified
    end
    # Projection metadata only covers source entities with actual cells. Every
    # cached owner must have a certified geometric entity record.
    for owner in node_entities
        haskey(metadata,(owner[1],Int(owner[2]))) || throw(ArgumentError(
            "$caller: native source node owner $owner has no certified entity projection"))
    end
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in node_entities]
    for (node,(dim,tag)) in enumerate(node_entities)
        if dim==0 && haskey(source_model.points,Int(tag))
            parameters[node]=Float64[]
        elseif dim==1 && haskey(source_model.curve_params,Int(tag))
            record=_dim01_record(m,(1,Int(tag)))
            raw=record!==nothing && !haskey(m.curve_params,Int(tag)) &&
                size(record.node_params,1)!=1
            raw && continue
            parameters[node]=Model.model_parametrization(source_model,1,Int(tag),
                                                        collect(@view cached.coords[:,node]))
        end
    end
    data=Elements.MixedEntityData(metadata;node_entities,node_parametric=parameters,
        external_node_tags=UInt64.(1:nnodes(cached)),block_entities=block_owners,
        external_element_tags=external)
    names=cached isa MixedMesh ? cached.physical_names : m.physical_names
    mesh=MixedMesh(cached.coords,blocks;physical_names=names,
        entity_data=data,elementary_entities=block_owners,periodic_links=links)
    result=_mixed_classification(mesh,class.entity,class.entities,node_entities,boundaries,owners)
    return mesh,result
end
