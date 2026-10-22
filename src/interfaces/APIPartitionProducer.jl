# Actual partition products for bare, first-order tetrahedral discrete Regions.
# One authoritative cell owns each UID; ghost entities expose query aliases.
using ..GeoExec: _geo_import_entity_receipts!
import ..MeshTypes

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

# `General.NumThreads` requests the parallel worker count for the O(cells)
# extraction stages below; `Threads.@threads` caps actual concurrency at the
# runtime thread pool. Every parallel stage writes disjoint per-index slots
# and merges serially in the original traversal order, so the produced mesh
# is identical under any requested worker count.
function _partition_workers()
    requested=get(OPTIONS,"General.NumThreads",1.0)
    (isfinite(requested) && requested>=2.0) || return 1
    return Int(requested)
end

# The O(cells) validation and element-quality scans in `MeshTypes`/
# `MeshElementQuality` honor the same `General.NumThreads` request through this
# registered source. Reading the live option lazily keeps session updates and
# `initialize`/`finalize` resets effective without a second store.
MeshTypes._set_parallel_worker_source!(_partition_workers)

# Declared boundary records for either entity representation: discrete
# entities carry them literally, while native `.geo` entities expose the same
# signed `(dim-1, tag)` pairs through their stored topology. `model_boundary`
# reports native Curve endpoints unsigned, so they are re-signed here to the
# discrete record convention (start positive, end negative); a closed curve
# has no boundary like a discrete record does.
#
# These are the *topological* signs — the direction the stored Curve Loop /
# Surface Loop / Curve endpoints declare. A produced mesh is not obliged to
# store cells along them: OCC faces carry shell-side `reversed` use flags,
# parametric frames can disagree with the wire traversal, `Reverse`/
# `ReorientMesh` flip stored cells, and the volume orientation pass may flip
# them again. Discrete records are user declarations and stay literal; the
# snapshot re-reads every native sign against the actual carrier cells (edge
# chain direction, surface border traversal, volume facet winding) so each
# emitted boundary describes the mesh it certifies.
function _partition_declared_boundary(model,dim,tag,caller)
    record=get(model.discrete,(dim,Int(tag)),nothing)
    record!==nothing && return record.boundary
    _model_entity_known(model,dim,Int(tag)) || throw(ArgumentError(
        "$caller: partition parent ($(dim),$(tag)) is not a model entity"))
    dim==0 && return Tuple{Int,Int}[]
    raw=Model._model_direct_boundary(model,dim,Int(tag),caller;canonical_orientation=true)
    dim!=1 && return raw
    if length(raw)==2 && all(entry->entry[1]==0,raw)
        first,last=raw[1][2],raw[2][2]
        first==last && return Tuple{Int,Int}[]
        return Tuple{Int,Int}[(0,first),(0,-last)]
    end
    isempty(raw) || throw(ArgumentError(
        "$caller: native Curve $tag boundary is malformed"))
    return raw
end

# Upstream partitions only the top-dimension element dual; every lower cell
# then inherits the partition of its first adjacent carrier at the nearest
# higher dimension, so a boundary Surface follows an incident Region instead
# of keeping an unrelated independent bucket. Apply the same correction to a
# dense catalog-order assignment — required for the declared-BRep facet
# matching in `_partition_brep_snapshot`, where each lower Surface cell's
# partition must match an actual volume facet carrier.
function _partition_lower_inherit(assignment::Vector{Int32},mesh)
    catalog=_cache_catalog(mesh)
    total=0;dims=Int[]
    for (_,dim,_,cells,_) in catalog
        n=size(cells,2)
        append!(dims,Iterators.repeated(Int(dim),n))
        total+=n
    end
    length(assignment)==total || return assignment
    top=isempty(dims) ? -1 : maximum(dims)
    top<1 && return assignment
    incidence=[Int[] for _ in 1:nnodes(mesh)]
    cell_nodes=Vector{Vector{Int32}}(undef,total)
    for (_,_,offset,cells,_) in catalog
        for column in axes(cells,2)
            index=offset+column
            rows=collect(Int32,Base.view(cells,:,column))
            cell_nodes[index]=rows
            for row in rows
                push!(incidence[row],index)
            end
        end
    end
    fixed=copy(assignment)
    for index in 1:total
        dim=dims[index]
        dim>=top && continue
        cell=cell_nodes[index]
        for step in 1:(top-dim)
            wanted=dim+step;hit=0
            for row in cell,other in incidence[row]
                dims[other]==wanted || continue
                # Only a cell actually carrying this lower cell (a Tet whose
                # facet is a Surface Triangle, a Triangle whose edge is a
                # Curve segment) may donate its partition — a neighbor that
                # merely shares a vertex or an edge belongs to a different
                # facet and can sit in another partition.
                all(node->node in cell_nodes[other],cell) || continue
                hit=Int(assignment[other]);break
            end
            hit!=0 && (fixed[index]=Int32(hit);break)
        end
    end
    return fixed
end

function _partition_tet_snapshot(model,cached,class,caller)
    class!==nothing && class.public_tags!==nothing || throw(ArgumentError(
        "$caller: topology and ghosts require actual classified UID authority"))
    _high_order_overlay(cached)===nothing || throw(ArgumentError(
        "$caller: partition topology and ghosts currently require first-order Tet4"))
    catalog=[row for row in _cache_catalog(cached,class) if !isempty(row[4])]
    all(row->row[1]==4,catalog) || throw(ArgumentError(
        "$caller: partition topology and ghosts currently support bare Tet4 Regions only; lower and mixed cells are unsupported"))
    cached isa MixedMesh && cached.partition_data!==nothing &&
        throw(ArgumentError("$caller: existing partition topology is unsupported"))
    parents=sort!(unique(Int(owner) for (_,_,_,_,owners) in catalog for owner in owners))
    all(tag->tag>0 && _model_entity_known(model,3,tag) &&
        isempty(_partition_declared_boundary(model,3,tag,caller)),parents) || throw(ArgumentError(
        "$caller: partition topology and ghosts require Regions without pre-existing lower BRep"))
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

# Connected components of the interface faces. With no `links` the default
# adjacency shares any edge — faces around a junction edge all land in one
# component. The lower-BRep caller instead passes explicit `(i,j)` links so
# each emitted component stays edge-manifold for the junction-fan
# certificate: manifold edges link their two incident faces and a
# nonmanifold edge only links the pairs its carriers share.
function _partition_interface_components(faces,links=nothing)
    edges(face)=_partition_edge_patterns(get(face,:msh,2))
    adjacency=[Int[] for _ in faces]
    if links===nothing
        incidence=Dict{Tuple{UInt64,UInt64},Vector{Int}}()
        for (index,face) in enumerate(faces),(a,b) in edges(face)
            push!(get!(()->Int[],incidence,minmax(face.uids[a],face.uids[b])),index)
        end
        for hits in values(incidence),i in hits,j in hits
            i<j && (push!(adjacency[i],j);push!(adjacency[j],i))
        end
    else
        for (i,j) in links
            i==j && continue
            push!(adjacency[i],j);push!(adjacency[j],i)
        end
    end
    seen=falses(length(faces));components=Vector{Int}[]
    for seed in eachindex(faces)
        seen[seed] && continue
        component=Int[seed];seen[seed]=true;cursor=1
        while cursor<=length(component)
            for neighbor in adjacency[component[cursor]]
                seen[neighbor] && continue
                seen[neighbor]=true;push!(component,neighbor)
            end
            cursor+=1
        end
        sort!(component);push!(components,component)
    end
    return components
end

function _partition_tet_production_plan(model,cached,class,count,assignment,caller;
                                        topology::Bool,ghosts::Bool)
    if class!==nothing && any(row->row[1] in _PARTITION_LOWER_TYPES && !isempty(row[4]),
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
    element_uids=copy(source.uids)
    maximum_uid=max(ELEMENT_TAG_MAX[],maximum(source.public.element_tags;init=UInt64(0)))
    if topology
        workers=_partition_workers()
        cell_facets=Vector{Vector{Tuple{NTuple{3,UInt64},NamedTuple}}}(undef,total)
        facet_keys=((1,3,2),(1,2,4),(2,3,4),(3,1,4))
        if workers>1
            Threads.@threads for index in 1:total
                entries=Tuple{NTuple{3,UInt64},NamedTuple}[]
                for face in facet_keys
                    rows=ntuple(slot->source.cells[face[slot],index],3)
                    uids=ntuple(slot->source.public.node_tags[rows[slot]],3)
                    key=Tuple(sort!(collect(uids)))
                    push!(entries,(key,(;rows,uids,index,child=cell_owners[index])))
                end
                cell_facets[index]=entries
            end
        else
            for index in 1:total
                entries=Tuple{NTuple{3,UInt64},NamedTuple}[]
                for face in facet_keys
                    rows=ntuple(slot->source.cells[face[slot],index],3)
                    uids=ntuple(slot->source.public.node_tags[rows[slot]],3)
                    key=Tuple(sort!(collect(uids)))
                    push!(entries,(key,(;rows,uids,index,child=cell_owners[index])))
                end
                cell_facets[index]=entries
            end
        end
        incidence=Dict{NTuple{3,UInt64},Vector{NamedTuple}}();order=NTuple{3,UInt64}[]
        for index in 1:total,(key,hit) in cell_facets[index]
            haskey(incidence,key) || push!(order,key)
            push!(get!(()->NamedTuple[],incidence,key),hit)
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
        workers=_partition_workers()
        destinations=Vector{Set{Int32}}(undef,total)
        if workers>1
            Threads.@threads for index in 1:total
                parts=Set{Int32}()
                for row in source.cells[:,index];union!(parts,node_parts[row]);end
                destinations[index]=parts
            end
        else
            for index in 1:total
                parts=Set{Int32}()
                for row in source.cells[:,index];union!(parts,node_parts[row]);end
                destinations[index]=parts
            end
        end
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
        periodic_links=cached isa MixedMesh ? cached.periodic_links : Elements.MixedPeriodicLink[],
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

# Assignment-only partition production (`Mesh.PartitionCreateTopology` and
# `Mesh.PartitionCreateGhostCells` both zero). Upstream still materializes the
# partition entity layer — `createNewEntities` runs before the topology gate —
# so one discrete child is born per (dimension, parent entity, partition)
# triple that owns at least one cell, carries the parent's physical
# memberships, and becomes the cell's new owner. No interface topology or
# ghost layer is produced; `$PartitionedEntities` therefore describes exactly
# the children below. This is the fallback `.geo` `PartitionMesh` uses for
# generated meshes whose declared lower BRep does not satisfy the production
# producers' prerequisites, and it accepts any classified first-order mesh —
# surfaces without volumes included.
function _partition_assign_production_plan(model,cached,class,count,assignment,caller)
    class!==nothing && class.public_tags!==nothing || throw(ArgumentError(
        "$caller: partitioning requires actual classified UID authority"))
    _high_order_overlay(cached)===nothing || throw(ArgumentError(
        "$caller: partitioning currently requires first-order cells"))
    cached isa MixedMesh && cached.partition_data!==nothing &&
        throw(ArgumentError("$caller: existing partition topology is unsupported"))
    catalog=[row for row in _cache_catalog(cached,class) if !isempty(row[4])]
    total=sum(size(row[4],2) for row in catalog;init=0)
    total>0 || throw(ArgumentError("$caller: the mesh has no elements"))
    length(assignment)==total || throw(ArgumentError(
        "$caller: partition assignment rows do not align with the cell catalog"))
    all(part->1<=part<=count,assignment) || throw(ArgumentError(
        "$caller: partition tags exceed num_part $count"))
    public=class.public_tags
    for (_,record) in Model._discrete_mesh_records_model(model)
        all(uid->haskey(public.node_indices,uid),record.node_tags) &&
        all(uid->haskey(public.element_indices,uid),record.element_tags) || throw(ArgumentError(
            "$caller: partitioning cannot omit caller-owned records outside the cache"))
    end
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),caller)
    _dim01_ingest!(assembly,model,cached,class;preserve_source_order=true)
    length(assembly.nodes)==nnodes(cached) && length(assembly.cells)==length(public.element_tags) ||
        throw(ArgumentError("$caller: actual cached and caller-owned partition rows do not agree"))
    global_entities,_=_dim01_entities(model,model_entities(model),assembly,eachindex(assembly.nodes))
    staged=Model._model_mesh_transaction_copy(model)
    triples=Tuple{Int,Int32,Int32}[]
    offset=0
    for (_,dim,_,cells,classified) in catalog
        for column in axes(cells,2)
            push!(triples,(dim,classified[column],assignment[offset+column]))
        end
        offset+=size(cells,2)
    end
    children=Dict{Tuple{Int,Int32,Int32},Int32}()
    child_info=Dict{Tuple{Int,Int},Tuple{Int,Int32,Int32}}()
    # Deterministic birth order — sorted (dimension, parent, partition) — so
    # child tags cannot depend on catalog storage order.
    for (dim,parent,part) in sort!(unique(triples))
        tag=Int32(_partition_entity_birth!(staged,dim,caller))
        children[(dim,parent,part)]=tag
        groups=Model._model_projection_physical_tags(staged,dim,Int(parent))
        staged.entity_physicals[(dim,Int(tag))]=copy(groups)
        for group in groups
            members=staged.physical[(dim,Int(group))]
            Int(tag) in members || push!(members,Int(tag))
        end
        child_info[(dim,Int(tag))]=(dim,parent,part)
    end
    rows_by_child=Dict{Tuple{Int,Int32},Vector{Int32}}()
    blocks=ElementBlock[];owners=Vector{Int32}[];uids=Vector{UInt64}[]
    offset=0
    for (msh,dim,_,cells,classified) in catalog
        n=size(cells,2)
        block_owner=Vector{Int32}(undef,n)
        for column in 1:n
            child=children[(dim,classified[column],assignment[offset+column])]
            block_owner[column]=child
            for row in cells[:,column]
                push!(get!(()->Int32[],rows_by_child,(dim,child)),row)
            end
        end
        physical=Int32[_dim01_physical(model,(dim,Int(owner))) for owner in classified]
        push!(blocks,ElementBlock(Int(msh),copy(cells),physical))
        push!(owners,block_owner)
        push!(uids,copy(public.element_tags[offset+1:offset+n]))
        offset+=n
    end
    part_entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    for (key,(dim,parent,part)) in child_info
        rows=unique!(sort!(get(rows_by_child,key,Int32[])))
        bounds=_partition_bounds(cached.coords,rows)
        part_entities[key]=Elements.MixedEntity(dim,Int32(key[2]),bounds;
            physical_tags=Model._model_projection_physical_tags(staged,dim,Int(parent)),
            partitions=Int32[part],parent=(dim,parent))
    end
    # Nodes keep their existing owner unless a partition child owns an
    # incident cell at a lower-or-equal dimension. Parameter-carrying nodes
    # prefer a child of their classified owner so the parameters remain valid
    # — the same `isless` merge the BRep producer applies.
    node_owners=copy(class.node_entities)
    sentinel=(4,Int32(typemax(Int32)))
    incident=fill(sentinel,nnodes(cached))
    incident_par=fill(sentinel,nnodes(cached))
    for (bi,(msh,dim,_,cells,_)) in enumerate(catalog)
        for column in axes(cells,2),row in cells[:,column]
            candidate=(dim,owners[bi][column])
            isless(candidate,incident[row]) && (incident[row]=candidate)
            parameters=class.node_parametric[row]
            if parameters!==nothing && dim==length(parameters) &&
                    isless(candidate,incident_par[row])
                entry=part_entities[(dim,Int(candidate[2]))]
                owner=(class.node_entities[row][1],Int(class.node_entities[row][2]))
                entry.parent==owner && (incident_par[row]=candidate)
            end
        end
    end
    for row in eachindex(incident)
        pick=class.node_parametric[row]===nothing ? incident[row] : incident_par[row]
        pick[1]<4 && (node_owners[row]=pick)
    end
    data=Elements._native_mixed_entity_data(global_entities;node_entities=node_owners,
        node_parametric=deepcopy(class.node_parametric),
        external_node_tags=copy(public.node_tags),
        block_entities=owners,external_element_tags=uids)
    partition=Elements.MixedPartitionData(count;entities=part_entities)
    target=MixedMesh(cached.coords,blocks;entity_data=data,elementary_entities=owners,
        physical_names=staged.physical_names,partition_data=partition,
        periodic_links=cached isa MixedMesh ? cached.periodic_links : Elements.MixedPeriodicLink[],
        ancillary_sections=cached isa MixedMesh ? cached.ancillary_sections : [],
        data_sections=cached isa MixedMesh ? cached.data_sections : Elements.MshDataSection[])
    validate(target).ok || throw(ArgumentError(
        "$caller: partition produced invalid or duplicate cells"))
    _geo_import_entity_receipts!(staged,target)
    identities=Dict(uid=>uid for uid in public.node_tags)
    declarations=union(Set(keys(global_entities)),Set(keys(part_entities)))
    _dim01_reconcile_records!(staged,declarations,target,identities;runtime_parameters=data.node_parametric)
    previous=model.meshing.mesh_authority
    ranks=previous!==nothing && previous.mesh===cached ? previous.private_rank : collect(UInt64,1:nnodes(target))
    staged.meshing.mesh_authority=Model._model_mesh_authority(target;caller,private_rank=ranks)
    cell_map=Dict(uid=>uid for uid in public.element_tags)
    return (;model=staged,mesh=target,class=_higher_classification_from_data(target,caller),
        authority=Set{Tuple{Int,Int32}}((dim,Int32(tag)) for (dim,tag) in declarations),
        max_node_tag=NODE_TAG_MAX[],max_element_tag=ELEMENT_TAG_MAX[],
        node_tag_map=identities,element_tag_map=cell_map,visibility_tag_map=cell_map)
end
