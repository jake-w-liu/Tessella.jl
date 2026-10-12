# Actual partition products for bare, first-order tetrahedral discrete Regions.
# One authoritative cell owns each UID; ghost entities expose query aliases.
using ..GeoExec: _geo_import_entity_receipts!

function _partition_entity_birth!(model,dimension,caller)
    maximum_tag=max(model.next_tag[dimension+1],maximum(
        (tag for (dim,tag) in model_entities(model) if dim==dimension);init=0))
    maximum_tag<typemax(Int32) || throw(ArgumentError(
        "$caller: partition entity tags exceed Int32"))
    model.next_tag[dimension+1]=maximum_tag
    tag=Model._alloc_tag!(model,dimension,0,caller)
    return Model.add_discrete_entity!(model,dimension,tag)
end

function _partition_bounds(coordinates,rows)
    isempty(rows) && return ntuple(_->0.0,6)
    ntuple(axis->axis<=3 ? minimum(coordinates[axis,row] for row in rows) :
        maximum(coordinates[axis-3,row] for row in rows),6)
end

function _partition_tet_snapshot(model,cached,class,caller)
    class!==nothing && class.public_tags!==nothing || throw(ArgumentError(
        "$caller: topology and ghosts require actual classified UID authority"))
    _high_order_overlay(cached)===nothing || throw(ArgumentError(
        "$caller: partition topology and ghosts currently require first-order Tet4"))
    catalog=[row for row in _cache_catalog(cached,class) if !isempty(row[4])]
    all(row->row[1]==4,catalog) || throw(ArgumentError(
        "$caller: partition topology and ghosts currently support bare Tet4 Regions only; lower and mixed cells are unsupported"))
    cached isa MixedMesh && (cached.partition_data!==nothing || !isempty(cached.periodic_links)) &&
        throw(ArgumentError("$caller: existing partition or periodic topology is unsupported"))
    isempty(model.meshing.periodic_uid_history) && isempty(model.periodic) || throw(ArgumentError(
        "$caller: partition topology and ghosts do not yet support periodic constraints"))
    all(parameter->parameter===nothing || isempty(parameter),class.node_parametric) ||
        throw(ArgumentError("$caller: parametric partition node ownership is not implemented"))
    parents=sort!(unique(Int(owner) for (_,_,_,_,owners) in catalog for owner in owners))
    all(tag->tag>0 && haskey(model.discrete,(3,tag)) &&
        isempty(model.discrete[(3,tag)].boundary),parents) || throw(ArgumentError(
        "$caller: partition topology and ghosts require discrete Regions without pre-existing lower BRep"))
    public=class.public_tags
    for (_,record) in Model._discrete_mesh_records_model(model)
        all(uid->haskey(public.node_indices,uid),record.node_tags) &&
        all(uid->haskey(public.element_indices,uid),record.element_tags) || throw(ArgumentError(
            "$caller: partition topology and ghosts cannot omit caller-owned records outside the cache"))
    end
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),caller)
    _dim01_ingest!(assembly,model,cached,class;preserve_source_order=true)
    length(assembly.nodes)==nnodes(cached) && length(assembly.cells)==length(public.element_tags) ||
        throw(ArgumentError("$caller: actual cached and caller-owned partition rows do not agree"))
    entities,_=_dim01_entities(model,model_entities(model),assembly,eachindex(assembly.nodes))
    blocks=ElementBlock[];owners=Vector{Int32}[];uids=Vector{UInt64}[]
    for (msh,dim,offset,cells,classified) in catalog
        dim==3 && size(cells,1)==4 || throw(ArgumentError("$caller: incomplete Tet4 authority"))
        physical=Int32[_dim01_physical(model,(3,Int(owner))) for owner in classified]
        push!(blocks,ElementBlock(msh,cells,physical));push!(owners,copy(classified))
        push!(uids,copy(public.element_tags[offset+1:offset+size(cells,2)]))
    end
    cells=reduce(hcat,(block.nodes for block in blocks))
    all(row->class.node_entities[row][1]==3 &&
        Int(class.node_entities[row][2]) in parents,unique(vec(cells))) || throw(ArgumentError(
        "$caller: partition Tet vertices with pre-existing lower-dimensional classification are unsupported"))
    diagnostic=validate(Mesh(cached.coords;tets=cells))
    diagnostic.ok || throw(ArgumentError("$caller: partition input is not a positive manifold Tet4 mesh"))
    return (;entities,blocks,owners,uids,public,cells,parents)
end

function _partition_interface_components(faces)
    edges(face)=_simplex_edge_patterns(get(face,:msh,2))
    incidence=Dict{Tuple{UInt64,UInt64},Vector{Int}}()
    for (index,face) in enumerate(faces),(a,b) in edges(face)
        push!(get!(()->Int[],incidence,minmax(face.uids[a],face.uids[b])),index)
    end
    seen=falses(length(faces));components=Vector{Int}[]
    for seed in eachindex(faces)
        seen[seed] && continue
        component=Int[seed];seen[seed]=true;cursor=1
        while cursor<=length(component)
            face=faces[component[cursor]];cursor+=1
            for (a,b) in edges(face),neighbor in incidence[minmax(face.uids[a],face.uids[b])]
                seen[neighbor] && continue
                seen[neighbor]=true;push!(component,neighbor)
            end
        end
        sort!(component);push!(components,component)
    end
    return components
end

function _partition_tet_production_plan(model,cached,class,count,assignment,caller;
                                        topology::Bool,ghosts::Bool)
    if class!==nothing && any(row->row[1] in (1,2,3,8,9,10,15) && !isempty(row[4]),
                              _cache_catalog(cached,class))
        return _partition_brep_production_plan(model,cached,class,count,assignment,caller;topology,ghosts)
    end
    source=_partition_tet_snapshot(model,cached,class,caller)
    total=size(source.cells,2)
    count<=total && length(unique(assignment))==count || throw(ArgumentError(
        "$caller: topology and ghosts require every declared partition to own at least one Tet4"))
    staged=Model._model_mesh_transaction_copy(model)
    original=reduce(vcat,source.owners);pairs=sort!(unique(zip(original,assignment)))
    children=Dict{Tuple{Int32,Int32},Int32}();part_entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    cell_owners=Int32[]
    for (parent,part) in pairs
        tag=Int32(_partition_entity_birth!(staged,3,caller));children[(parent,part)]=tag
        groups=Model._model_projection_physical_tags(staged,3,Int(parent))
        staged.entity_physicals[(3,Int(tag))]=copy(groups)
        for group in groups
            members=staged.physical[(3,Int(group))]
            Int(tag) in members || push!(members,Int(tag))
        end
        selected=findall(index->original[index]==parent && assignment[index]==part,1:total)
        rows=unique(vec(source.cells[:,selected]))
        part_entities[(3,Int(tag))]=Elements.MixedEntity(3,tag,_partition_bounds(cached.coords,rows);
            physical_tags=groups,partitions=Int32[part],parent=(3,parent))
    end
    for index in 1:total;push!(cell_owners,children[(original[index],assignment[index])]);end
    node_owners=copy(class.node_entities)
    incident=fill(Int32(typemax(Int32)),nnodes(cached))
    for index in 1:total,row in source.cells[:,index]
        incident[row]=min(incident[row],cell_owners[index])
    end
    for row in eachindex(incident)
        incident[row]==typemax(Int32) || (node_owners[row]=(3,incident[row]))
    end
    boundaries=Dict(tag=>Int32[] for tag in values(children))
    blocks=copy(source.blocks);owners=Vector{Int32}[];offset=0
    for block in blocks
        n=size(block.nodes,2);push!(owners,cell_owners[offset+1:offset+n]);offset+=n
    end
    element_uids=copy(source.uids);maximum_uid=ELEMENT_TAG_MAX[]
    if topology
        incidence=Dict{NTuple{3,UInt64},Vector{NamedTuple}}();order=NTuple{3,UInt64}[]
        for index in 1:total,face in ((1,3,2),(1,2,4),(2,3,4),(3,1,4))
            rows=ntuple(slot->source.cells[face[slot],index],3)
            uids=ntuple(slot->source.public.node_tags[rows[slot]],3)
            key=Tuple(sort!(collect(uids)))
            haskey(incidence,key) || push!(order,key)
            push!(get!(()->NamedTuple[],incidence,key),(;rows,uids,index,child=cell_owners[index]))
        end
        grouped=Dict{Tuple{Int32,Int32},Vector{NamedTuple}}()
        for key in order
            hits=incidence[key]
            length(hits)<=2 || throw(ArgumentError("$caller: partition facet has more than two Tet carriers"))
            length(hits)==2 || continue
            left,right=hits
            reversed=any(left.uids==ntuple(slot->right.uids[mod1(4-slot+rotation,3)],3) for rotation in 0:2)
            reversed || throw(ArgumentError("$caller: shared Tet facet carriers have equal winding"))
            assignment[left.index]==assignment[right.index] && continue
            push!(get!(()->NamedTuple[],grouped,minmax(left.child,right.child)),
                (rows=left.rows,uids=left.uids,index=left.index,child=left.child,other=right.child))
        end
        for pair in sort!(collect(keys(grouped)))
            faces=grouped[pair]
            for component in _partition_interface_components(faces)
                primary=faces[first(component)].child
                parent=part_entities[(3,Int(primary))].parent
                tag=Int32(_partition_entity_birth!(staged,2,caller))
                rows=Matrix{Int32}(undef,3,length(component));tags=UInt64[]
                for (column,index) in enumerate(component)
                    face=faces[index]
                    rows[:,column]=face.child==primary ? collect(face.rows) :
                        Int32[face.rows[1],face.rows[3],face.rows[2]]
                    maximum_uid<typemax(UInt64) || throw(ArgumentError("$caller: interface cell UIDs exceed UInt64"))
                    maximum_uid+=UInt64(1);push!(tags,maximum_uid)
                end
                total+sum(size(block.nodes,2) for block in blocks if block.msh==2;init=0)+length(component)<=typemax(Int32) ||
                    throw(ArgumentError("$caller: partition cell rows exceed Int32"))
                membership=sort!(Int32[part_entities[(3,Int(child))].partitions[1] for child in pair])
                part_entities[(2,Int(tag))]=Elements.MixedEntity(2,tag,_partition_bounds(cached.coords,unique(vec(rows)));
                    partitions=membership,parent)
                push!(blocks,ElementBlock(2,rows));push!(owners,fill(tag,length(component)));push!(element_uids,tags)
                for child in pair;push!(boundaries[child],child==primary ? tag : -tag);end
                for row in rows
                    (node_owners[row][1]>2 || tag<node_owners[row][2]) && (node_owners[row]=(2,tag))
                end
            end
        end
    end
    for ((dim,tag),entity) in collect(part_entities)
        dim==3 || continue
        boundary=boundaries[Int32(tag)]
        staged.discrete[(3,tag)].boundary=Tuple{Int,Int}[(2,Int(value)) for value in boundary]
        part_entities[(dim,tag)]=Elements.MixedEntity(dim,tag,entity.bbox;
            physical_tags=entity.physical_tags,boundaries=boundary,
            partitions=entity.partitions,parent=entity.parent)
    end
    global_entities=copy(source.entities);ghost_entities=Tuple{Int32,Int32}[]
    ghost_elements=Elements.MixedGhostElement[]
    if ghosts
        node_parts=[Set{Int32}() for _ in 1:nnodes(cached)]
        for index in 1:total,row in source.cells[:,index];push!(node_parts[row],assignment[index]);end
        destinations=[Set{Int32}() for _ in 1:total]
        for index in 1:total,row in source.cells[:,index];union!(destinations[index],node_parts[row]);end
        for index in 1:total
            delete!(destinations[index],assignment[index]);isempty(destinations[index]) && continue
            push!(ghost_elements,Elements.MixedGhostElement(index,assignment[index],sort!(collect(destinations[index]))))
        end
        for part in 1:count
            tag=Int32(_partition_entity_birth!(staged,3,caller));push!(ghost_entities,(tag,Int32(part)))
            selected=findall(index->Int32(part) in destinations[index],1:total)
            rows=unique(vec(source.cells[:,selected]))
            global_entities[(3,Int(tag))]=Elements.MixedEntity(3,tag,_partition_bounds(cached.coords,rows))
        end
    end
    data=Elements._native_mixed_entity_data(global_entities;node_entities=node_owners,
        node_parametric=deepcopy(class.node_parametric),external_node_tags=copy(source.public.node_tags),
        block_entities=owners,external_element_tags=element_uids)
    partition=Elements.MixedPartitionData(count;entities=part_entities,ghost_entities,ghost_elements)
    target=MixedMesh(cached.coords,blocks;entity_data=data,elementary_entities=owners,
        physical_names=staged.physical_names,partition_data=partition,
        ancillary_sections=cached isa MixedMesh ? cached.ancillary_sections : [],
        data_sections=cached isa MixedMesh ? cached.data_sections : Elements.MshDataSection[])
    _geo_import_entity_receipts!(staged,target)
    identities=Dict(uid=>uid for uid in source.public.node_tags)
    declarations=union(Set(keys(global_entities)),Set(keys(part_entities)))
    _dim01_reconcile_records!(staged,declarations,target,identities;runtime_parameters=data.node_parametric)
    rebound=_higher_classification_from_data(target,caller)
    authority=Set{Tuple{Int,Int32}}((dim,Int32(tag)) for (dim,tag) in declarations)
    previous=model.meshing.mesh_authority
    ranks=previous!==nothing && previous.mesh===cached ? previous.private_rank : collect(UInt64,1:nnodes(target))
    staged.meshing.mesh_authority=Model._model_mesh_authority(target;caller,private_rank=ranks)
    cell_map=Dict(uid=>uid for uid in source.public.element_tags)
    return (;model=staged,mesh=target,class=rebound,authority,
        max_node_tag=NODE_TAG_MAX[],max_element_tag=maximum_uid,
        node_tag_map=identities,element_tag_map=cell_map,visibility_tag_map=cell_map)
end
